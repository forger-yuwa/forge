# codex レビュー: performance-3d-node-sst-speedup (result)

- **plan**: [`plans/active/performance-3d-node-sst-speedup.md`](../../plans/active/performance-3d-node-sst-speedup.md)
- **stage**: `result` (diff base `0512823d`)
- **date**: 2026-09-12
- **commit**: `9333df5c` (feature/perf-3d-speedup)
- **codex**: effort `high`, 5.9 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M4/m2

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**NO-GO**

高速化は生ログで確認でき、WALE の double 復帰・温度反転の研磨反復化も実装されています。  
ただし、改訂後の合否基準を実測が満たしておらず、最終版の検証証拠と完了条件にも不足があります。現状の `accepted` 移行は認められません。

1. **Major — 未収束のまま、収束保証を完了扱いにしている**

   [plan:112](/home/sano/work/forge-perf/plans/active/performance-3d-node-sst-speedup.md:112) は「未収束なら比較を確定しない」と規定しています。一方、[recommended-settings.md:48](/home/sano/work/forge-perf/procedures/recommended-settings.md:48) には、撤回したはずの「収束解不変」が残っています。

   `check_convergence.py` を再実行した結果は次のとおりです。

   - 基準 `run_0234_user_node3d_sst_dry_half_ext_ps_newdef`：**NOT CONVERGED**
   - 回収証拠の `run_0410_sweep5_thermofl`、`run_0413_sweep4_default`、`run_0414_cfl8_sweep4`：**NOT CONVERGED**

   基準は主に `roe` が plateau、変更版は `roUy/roK` が plateau です。同じ判定区分でも停滞している変数は異なります。

   壁圧 CSV の 10–95 mm 区間では、`run_0410` の輪郭壁相対差 `1.10e-6`、`run_0414` の `8.05e-5` を再現できました。ただし、これは最終出力間の差です。壁圧時系列の頭打ちを示す VERDICT は確認できませんでした。

   **対案:** 「収束解不変」を推奨設定から撤回し、有限時間継続で観測した差として記述する。収束保証と壁圧定常性の未達を §5.1 に残し、現在の完了条件を満たしたとは扱わない。

2. **Major — 改訂後の絶対基準にも不合格があり、「全場」の報告値も誤っている**

   [plan:108](/home/sano/work/forge-perf/plans/active/performance-3d-node-sst-speedup.md:108) は速度成分の正規化最大差を `≤1e-4` としています。しかし、[cmp_head3.txt:5](/home/sano/work/forge-perf/case/16.nozzle_wys/_aws_perf_evidence/run_0401_perf_verify/cmp_head3.txt:5) の `Uz` は **`4.84e-4`** です。基準同士も `4.78e-4` で、現行基準は基準バイナリ自身を不合格にします。

   また、[plan:178](/home/sano/work/forge-perf/plans/active/performance-3d-node-sst-speedup.md:178) の凝縮回帰「全場 ≤1.6e-5」は、`run_0456_perf_regress_node2d_cond/{base_r1,base_r2,final}/res_300.h5` の再集計と食い違います。

   | 変数 | 基準同士 | 基準対最終版 |
   |---|---:|---:|
   | `Uy` | `9.82e-6` | `1.60e-5` |
   | `h0` | `3.20e-5` | `2.09e-5` |
   | `vis_turb` | `9.23e-4` | **`1.08e-3`** |

   これは実装劣化の証明ではありません。しかし、採否基準を通過したという報告にはできません。

   **対案:** 速度の正規化尺度を成分別最大値か速度ベクトルの尺度か明記し、基準同士の再現性を踏まえて判定方法を確定する。その方法で全 `VALUE/*` を再集計し、除外項目・超過項目を明示する。「全場」の誤記は README も修正する。

3. **Major — 最終ソースと、主要性能・長時間検証の対応が確定していない**

   A10G の `small4` は **34.71 / 34.69 ms/step**、`h4` は **33.87 / 33.88 ms/step** で、速度改善自体は裏付けられます。

   ただし、後者のログが記録するのは `forge_head4`・SHA256 `ba31edbe…` です。[ハッシュ一覧:7](/home/sano/work/forge-perf/case/16.nozzle_wys/_aws_perf_evidence/binary_and_mesh_sha256.txt:7) に、そのバイナリが修正後の `ac9262e8` に対応するビルド記録はありません。場の比較ファイルも `cmp_head3.txt` までです。

   さらに、[case README:75](/home/sano/work/forge-perf/case/16.nozzle_wys/README.md:75) は `run_0413` のバイナリを **`308d799f`** と明記しています。これは研磨反復化前です。これを [plan:165](/home/sano/work/forge-perf/plans/active/performance-3d-node-sst-speedup.md:165) の「最終バイナリで12000 step」の代用にはできません。

   **対案:** 最終ソース commit・ビルド条件・バイナリ SHA256・入力ハッシュを結び付け、同じ最終バイナリで性能計測と長時間比較を実施する。準定常 VERDICT と壁圧時系列も保存し、§5.1 に未完了として記載する。

4. **Major — `lineImplicit` の回帰省略理由がコードと合わない**

   [plan:180](/home/sano/work/forge-perf/plans/active/performance-3d-node-sst-speedup.md:180) は、キャッシュ・パックが自動 off になることを理由に専用 run を省いています。

   しかし [timeIntegration_d.cu:858](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/timeIntegration_d.cu:858) の sweep 0 gather 省略と、共有カーネルの `const __restrict__` 化は `lineImplicit` でも適用されます。SLAU・粘性・スカラー輸送の変更も通ります。オプションの無効化は変更全体の無効化ではありません。

   **対案:** 最終版で `lineImplicit: 1`、特に `dualTime + lineKFreeze` の既存ケースを基準版と比較する。node 軸拘束も今回の追加回帰は cell 軸対称で代替されているため、対象分岐を通す検証を §5.1 に追加する。

5. **Minor — 温度反転試験は、まだ「float 再格納ドリフト」の試験になっていない**

   12反復・冷間開始の追加は確認しました。しかし [test_thermo_float.cpp:77](/home/sano/work/forge-perf/solver_density_cuda/tools/test_thermo_float.cpp:77) は再構成した `e2` を double のまま再反転しています。本番の [dependentVariables_d.cu:194](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/dependentVariables_d.cu:194) にある `ρ(e+ek)` の float 格納・再読込を再現していません。

   また、`maxdrift` は絶対温度差と相対温度差を混在させ、[試験:90](/home/sano/work/forge-perf/solver_density_cuda/tools/test_thermo_float.cpp:90) の合否には使っていません。

   **対案:** 温度反転精度と保存量再格納の試験を分離する。後者は本番同様の float 保存量・組成構築を通し、単位を揃えたドリフトと有限性を合否に含める。

6. **Minor — 現在仕様・索引・残作業表に撤回済み記述が残る**

   - [performance.md:61](/home/sano/work/forge-perf/methods/architecture/performance.md:61)：研磨を常に1段・double評価1回と説明。実装は最大3段。
   - [performance.md:66](/home/sano/work/forge-perf/methods/architecture/performance.md:66)：採否が依然「ノイズ床の2倍以内」。
   - [plan:135](/home/sano/work/forge-perf/plans/active/performance-3d-node-sst-speedup.md:135)、[plans/README.md:26](/home/sano/work/forge-perf/plans/README.md:26)：不採用の per-cell 前計算案が残存。
   - §5.1 #9/#10：見出しは「A/B待ち」、本文は結果確定。今回の未達検証は残作業として登録されていない。

   `methods/index.md` の登録と `check_plans.py` の **PASS** は確認しましたが、機械的 lint は上記の内容矛盾を検出しません。

   **対案:** §5.1 を現在の未完了作業の正本に整理し、過去の実験経過を変更ログへ移す。現在仕様と索引を同じ最終方針に揃える。

推奨は、**`in_progress` を維持し、合否基準と未達事項を確定したうえで、最終バイナリの同条件 A/B・不足分岐の回帰を揃えて再レビューすること**です。確認した差分では、面状態での熱力学評価、境界半割面の除外、拘束行、周期時のパック無効化は維持されています。今回の NO-GO は、新たな符号誤りを断定したものではなく、採用判断の根拠がまだ成立していないためです。

ファイルは変更していません。本レビューは **plan 未反映**です。

指摘数: Critical 0 / Major 4 / Minor 2
