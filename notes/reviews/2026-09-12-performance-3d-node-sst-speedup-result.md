# codex レビュー: performance-3d-node-sst-speedup (result)

- **plan**: [`plans/active/performance-3d-node-sst-speedup.md`](../../plans/active/performance-3d-node-sst-speedup.md)
- **stage**: `result` (diff base `0512823d`)
- **date**: 2026-09-12
- **commit**: `3aaf33ef` (feature/perf-3d-speedup)
- **codex**: effort `high`, 6.8 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M5/m2

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**NO-GO**

高速化自体はローカル生ログで確認できました。しかし、WALE の float 化に overflow 回帰があり、温度反転の単体試験も本番経路を検証していません。さらに、場の差が採否基準を超え、収束・主要 A10G 実測の裏付けも不足しています。

1. **Major — WALE の float 化で、有限入力から `vis_turb=Inf` が発生する**

   根拠: [turbulent_viscosity_d.cu:98](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/turbulent_viscosity_d.cu:98)。`pow(SdijSdij, 3.0f/2.0f)` は、後段の `Ls²` を掛ける前に overflow します。

   コードの演算を NumPy の float32/float64 で再現すると、`∇U=diag(1e7,-1e7,0) [s⁻¹]`、`ρ=1`、`Ls=1e-6 m` に対し、旧式は `μt=8.6964e-7 Pa·s`、新式は **`Inf`**。分母は両方とも有限の `6.2593e35` です。これは GPU run の実測ではなく、変更式の数値再現ですが、型変更による演算範囲の欠陥を示しています。

   **対案:** 今回の変更から WALE の高次べき乗の float 化を戻す。スケーリングによる安定な float 評価は別途検証し、高勾配・ゼロ勾配・純せん断を回帰に含める。

2. **Major — ハイブリッド温度反転の「単体 PASS」は本番実装の保証になっていない**

   根拠: 本番は [dependentVariables_d.cu:157](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/dependentVariables_d.cu:157) で **12 反復**ですが、[test_thermo_float.cpp:72](/home/sano/work/forge-perf/solver_density_cuda/tools/test_thermo_float.cpp:72) は **20 反復**を渡しています。本番の float 組成正規化も試験では通らず、`cp/h` の Taylor 再構成・`roe` 再格納の反復ドリフトも採否対象外です。

   [thermo_d.cuh:608](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/thermo_d.cuh:608) は float 段の収束状態を確認せず、必ず double 1 段で終了します。同梱 H2O 係数・datum 298.15 K の式をホストで再現すると、目標 `5900 K`、初期推定 `50.1 K` で、12 反復後の研磨結果は約 `5900.0191 K`、相対誤差 **`3.24e-6`**。主張の `1e-8` を超えます。

   **対案:** 試験を本番の12反復・組成構築・再格納まで含むものにする。float 段未収束や研磨補正過大の場合は従来 double 経路へフォールバックし、その後に精度保証を更新する。

3. **Major — ローカル回帰の「全場差 ≤ ノイズ床の1.5倍」は実データと矛盾する**

   根拠: [plan:164](/home/sano/work/forge-perf/plans/active/performance-3d-node-sst-speedup.md:164) の主張に対し、`case/16.nozzle_wys/run_0450_perf_regress_node2d_tp/` の `base_r1`・`base_r2`・`new_tf` の各 `res_300.h5` を再集計しました。

   `max|base_r1|` で正規化した最大差は以下です。

   | 場 | base同士の差 | base_r1対new_tf | ノイズ比 |
   |---|---:|---:|---:|
   | `ro` | `4.34e-7` | `1.30e-6` | **3.00倍** |
   | `roY0` | `4.39e-7` | `1.32e-6` | **3.00倍** |
   | `Ux` | `1.05e-6` | `2.33e-6` | **2.21倍** |

   `new` でも `ro` は約2.60倍です。絶対的には小さい差ですが、[plan:107](/home/sano/work/forge-perf/plans/active/performance-3d-node-sst-speedup.md:107) の **2倍以内**という採否基準を満たしません。

   **対案:** 「全てノイズ床内」を撤回し、基準の反復実行を増やしてノイズ分布を確定する。そのうえで、超過する変数について変更項目別に原因を切り分ける。

4. **Major — `NOT CONVERGED` を根拠に「収束解不変」を確定している**

   根拠: [plan:110](/home/sano/work/forge-perf/plans/active/performance-3d-node-sst-speedup.md:110) は「未収束なら比較を確定しない」としながら、直後で `run_0234/0410/0411` を `NOT CONVERGED plateau` と記録しています。[recommended-settings.md:40](/home/sano/work/forge-perf/procedures/recommended-settings.md:40) には「収束解不変」が掲載されています。

   今回、ローカル回帰4ケースの `new` に `check_convergence.py` を実行した結果も、**全て `NOT CONVERGED (stalled/plateau)`** でした。`run_0450/new` と `run_0451/new` の準定常判定は、保存出力が2枚しかなく、双方とも **`TRANSIENT-UNSETTLED`** です。確認した最終場と残差CSVに NaN/Inf はありませんが、収束の証明にはなりません。

   また、`pmax/machmax` の `STEADY` は、報告対象である壁圧分布の定常性を保証しません。

   **対案:** 現状は「有限時間継続での差」と記述する。収束解不変を完了条件とするなら基準・変更版双方の `PASS` と、壁圧そのものの時系列に対する定常性判定を揃える。判定条件を「同じ未収束区分」へ緩めて完了扱いにしない。

5. **Major — 主要 A10G 結果と最終実装の回帰を独立検証できない**

   根拠: [case README:67](/home/sano/work/forge-perf/case/16.nozzle_wys/README.md:67) に記載された AWS の `run_0400/0401/0410–0414`、基準 `run_0234` はこの作業ツリーにありません。従って、`34.7 ms/step`、壁圧差、収束・準定常 VERDICT、入力同一性は README の記述以上には確認できませんでした。

   一方、`case/16.nozzle_wys/run_0455_perf_big_norcm_bench/` の生ログは確認でき、基準 `177.01/173.38`、変更版 `63.90/65.95 ms/step` でした。**ローカルの大幅な高速化は裏付けられますが、A10G の最終合否を代替しません。**

   回帰記録も [plan:163](/home/sano/work/forge-perf/plans/active/performance-3d-node-sst-speedup.md:163) に予定した bump・NACA node がなく、変更対象の WALE、周期、軸対称、凝縮、line-implicit を最終版で確認した証拠がありません。

   **対案:** 最終ソースとバイナリのハッシュ、入力ハッシュ、生ログ、比較出力、VERDICTをレビュー可能な場所に揃える。修正後の最終版で、変更された分岐を通る回帰を実施する。

6. **Minor — ベンチマークの実装と現在仕様が一致していない**

   根拠: [performance.md:13](/home/sano/work/forge-perf/methods/architecture/performance.md:13) は元 run の config 書換え・`FORGE_PROFILE=1`・旧ログ名を案内しますが、[bench_steps.sh:20](/home/sano/work/forge-perf/solver_density_cuda/tools/bench_steps.sh:20) は専用ディレクトリ・既定 `FORGE_PROFILE=0`・日時付きログです。

   また、計時は [main.cpp:1720](/home/sano/work/forge-perf/solver_density_cuda/main.cpp:1720) から全ステップを含み、plan が要求する warm-up 除外区間がありません。同じ label の再使用も `mkdir -p` で既存 run に上書きできます。

   **対案:** 同一プロセス内の warm-up 後を計時対象にし、既存計測ディレクトリへの再投入を拒否する。実装に合わせて手順を更新する。

7. **Minor — 現在仕様と残作業表に、撤回済み方針・未完了事項が混在する**

   根拠:

   - [thermophysics.md:73](/home/sano/work/forge-perf/methods/thermophysics.md:73) は現在も「反転と多項式評価は double」と説明しています。
   - [performance.md:56](/home/sano/work/forge-perf/methods/architecture/performance.md:56)、plan §5、`plans/README.md` に、撤回したセル前計算案が残っています。
   - [plan:152](/home/sano/work/forge-perf/plans/active/performance-3d-node-sst-speedup.md:152) は同じ行で CFL 8 を「試験済み」と「未試験」の両方にしています。
   - §2で対象外とした RCM が実装され、§5.1 #10はAWS試験待ちですが、完了条件は #1–#8 のままです。
   - ローカル回帰 `run_0450/0451` の run 索引がなく、NACA の対象ディレクトリには `README.md` 自体がありません。

   **対案:** 現在仕様を最終実装に揃え、§5.1を「完了・却下・未完了」に整理する。RCMの扱いと不足する回帰・証拠整備を残作業に明記し、run索引を補完する。`thermoFloat`既定1と推奨sweep4の決定自体は、planに記録されています。

推奨は、**`in_progress` を維持し、1・2の数値実装を修正した後、3〜5の検証を最終版でやり直すこと**です。その証拠と文書を揃えてから再度 result レビューを受けてください。高速化の成果は認められますが、現状を `accepted` に移すことは推奨しません。

ファイルは変更していません。本レビューの指摘は **plan 未反映**です。

指摘数: Critical 0 / Major 5 / Minor 2
