# codex レビュー: condensation-float-speedup (result)

- **plan**: [`plans/active/condensation-float-speedup.md`](../../plans/active/condensation-float-speedup.md)
- **stage**: `result` (diff base `3b77cc4b`)
- **date**: 2026-09-13
- **commit**: `de095c4b` (feature/perf-3d-speedup)
- **codex**: effort `high`, 7.7 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M4/m2
- **focus**: 凝縮 float 化の result レビュー 2 回目。1 回目 (NO-GO M6/m2) の採用内容 §5.1 #8–#14 (表範囲外の double 退避 + SASS 同一, 反転成功条件の統一 + 強制失敗試験 run_0461, cond_axis_h0.py の cell 対応, 単体の有限性/蒸発端, 周期保存 run_0053, cell H2O run_0457, 発達場速度 run_0463, 文書同期) を検証してください。A10G はインスタンスが起動できず未計測 (plan は達成判定に使わないと明記)。この状態で accepted に移してよいか、A10G 未計測をどう扱うべきかも判定してください。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**

反転成功条件の統一、`roe` 保護、cell 用 `h0` 評価、cell H2O 回帰は確認できました。  
ただし、表範囲外への退避に実装漏れがあり、周期保存・絶対 `h0` 保存・A10G 性能の完了条件も満たしていません。現状の `accepted` 移動には反対します。

1. **Major — 表範囲外でも、double へ退避する前に核生成を消してしまう。**

   [condensationSourceKernels_d.cuh:322](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/condensationSourceKernels_d.cuh:322) で端クランプした表から `lnS` を作り、350 行で dry 判定して return します。範囲外の double 退避は、その後の353行です。したがって、**実際には過飽和でも、表の端値では未飽和と判定される低温セル**が退避できません。

   コードの式を再評価した反例は、pure N2、`T=15 K`、`g=Q0=0`、`p_v=0.9·p_sat(20 K)=7.6068e-9 Pa`。密度は `1.7086e-12 kg/m³` で kernel の密度ガードを通ります。旧式では `S≈1.979e6`、上限適用前の `lnJ≈149.8` ですが、float kernel は `S≈0.9` と判定してソースをゼロにします。これは実 run の観測ではなく、[既存物性式:130](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/condensationProperties_d.cuh:130) に対する再評価です。

   **対案:** 範囲外判定を早期退出より前に移す。範囲外の dry 最適化を残すなら、旧 double 飽和圧で未飽和を確認する。device 試験には現在の高温・液相あり状態に加え、**低温・液相なし・実際には過飽和**の状態を追加する。

2. **Major — 周期試験は「融合による差が小さい」証拠であり、要求した保存性試験の合格ではない。**

   [plan:132](/home/sano/work/forge-perf/plans/active/condensation-float-speedup.md:132) は周期20 step の保存誤差 `≤1e-6` を要求しています。`case/09.Taylor-Green/run_0053_cond_periodic_conservation/` の HDF5 から、`ΣρQ0 V` の相対変化を再計算すると次の結果です。

   | ラベル | step | 相対変化 |
   |---|---:|---:|
   | `merged1` | 1 | `5.59965e-6` |
   | `condf1` | 1 | `5.60020e-6` |
   | `merged` | 20 | `1.14315e-3` |
   | `condf` | 20 | `1.11121e-3` |

   `Q1/Q2` も同程度です。既存版にも誤差がある点は確認できますが、1 step の前後差だけで保存性ゲートを完了にはできません。`ρg` は平衡射影を含むため、移流単独の検証にもなっていません。

   また、[plan:192](/home/sano/work/forge-perf/plans/active/condensation-float-speedup.md:192) の「`condensation-followups` に記録」は、[転記先の残作業表:40](/home/sano/work/forge-perf/plans/active/condensation-followups.md:40) に存在しません。

   **対案:** 同じ凍結状態で、移流だけの4モーメント残差と周期 seam の収支を直接比較する。20 step の保存誤差は未解決として残し、後続 plan に実際に転記する。既存誤差の原因を「dual-time と clamp」と確定するには、処理別の収支が必要です。

3. **Major — `h0` の PASS は、当初の必須絶対基準を緩和した結果である。**

   cell 座標対応と参照側の定常性検査は修正されています。再実行結果も §5.1 #10 の数値を裏付けます。

   | run／比較ラベル | 新版／基準の偏差 [J/kg] | 現ツール |
   |---|---:|---|
   | case/34 `run_0100_merge_regress_cell_air`、`condf/merged` | 428.542／434.125 | `PASS` |
   | case/34 `run_0103_condf_regress_n2_cell`、`condf/merged` | 482.370／481.952 | `PASS` |
   | case/34 `run_0104_condf_regress_n2_node`、`condf/merged_r1` | 128.914／128.935 | `PASS` |
   | case/16 `run_0460_condf_phys_node2d_cond`、`condf/merged` | 271.919／271.850 | `PASS` |

   全比較で末尾は `STEADY`。ただし cell の2件は、[plan:150](/home/sano/work/forge-perf/plans/active/condensation-float-speedup.md:150) の必須絶対上限 `300 J/kg` を超えています。[cond_axis_h0.py:108](/home/sano/work/forge-perf/solver_density_cuda/tools/cond_axis_h0.py:108) が `max(300, ref+10)` を許しているため PASS です。

   **対案:** 「基準からの非劣化 PASS」と「絶対保存基準 FAIL」を別判定で残す。§4.3・§8を満たした扱いにせず、cell の既存保存誤差を残作業として明示する。高速化による悪化が小さいという結論は支持できます。

4. **Major — A10G の性能目標は未検証で、現行バイナリの起動区間検証も不足している。**

   `case/16.nozzle_wys/run_0463_perf_cond3d_local_dev_bench/` のログから、発達場の **220.55／221.62／220.81 → 114.38／112.10／114.61 ms/step**、約 **1.94倍**は確認できました。新版ログに反転失敗警告はありません。

   ただし、この新版の SHA は `e009473e06f9323f…`。起動区間 `run_0459_perf_cond3d_local_bench/` の **197→84 ms/step** は、それ以前の `0c82a4d3bf5e5869…` です。double 退避・蒸発 Jacobian の変更後の起動速度を裏付けません。

   [plan:159](/home/sano/work/forge-perf/plans/active/condensation-float-speedup.md:159) と [完了条件:228](/home/sano/work/forge-perf/plans/active/condensation-float-speedup.md:228) は、A10G の発達場 `≤50`、起動 `≤45 ms/step` を依然要求しています。起動不能は計測延期の理由であり、完了条件の免除にはなりません。

   **対案:** A10G は「未検証」のまま残す。修正後の同一バイナリで dry・起動・発達場を計測し、SHA・入力・ログを対応させる。帯域比換算を達成判定に使わない方針は正しいです。

5. **Minor — #14 の文書同期は未完了で、現在の分岐・許容値と食い違う。**

   根拠は以下です。

   - [methods/condensation.md:874](/home/sano/work/forge-perf/methods/condensation.md:874) は範囲外を依然「端クランプ」「物理的に凝縮しない領域」と説明。
   - [同:908](/home/sano/work/forge-perf/methods/condensation.md:908) は蒸発端の許容を `|Δsj|dt≤1e-6` としていますが、[device 試験:150](/home/sano/work/forge-perf/solver_density_cuda/tests/unit/test_cond_float_device.cu:150) は `1e-4`。
   - [plan:191](/home/sano/work/forge-perf/plans/active/condensation-float-speedup.md:191) は蒸発 `λ≤5e-4` を「§4.3 に明記」としていますが、153行は `1e-5` のまま。
   - FP64ゼロの監査結果は、現在の double 退避・double 蒸発 Jacobian を含む実装全体には適用できません。
   - 「平衡形の EOS は常に double」も不正確です。[dependentVariables_d.cu:140](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/dependentVariables_d.cu:140) が別扱いするのは `condEquilibrium==2` で、`==1` の湿潤 TP セルはハイブリッド反転へ進みます。

   **対案:** 現在仕様、§4、§5.1を同じ分岐表・許容値に統一し、旧監査値には対象コミットを付ける。double 退避を含む最新監査結果を別に残す。

6. **Minor — 「VERDICT 保存済み」と実ファイルが一致せず、訂正前の説明も残っている。**

   [plan:194](/home/sano/work/forge-perf/plans/active/condensation-float-speedup.md:194) に保存済みとありますが、確認した `run_0456–0458`、`run_0460`、case/34 `run_0100–0104` に `CONVERGENCE_VERDICT.txt` 等はありません。

   `run_0460` の再実行結果は、新旧とも **`NOT CONVERGED (stalled/plateau)`**、`pmax/machmax` は **`OVERALL: ALL STEADY`**。場比較は **32/32 PASS**ですが、`rog_0` 差 `1.31e-5` はノイズ `5.71e-6` の2倍を超え、絶対基準での合格です。[README:75](/home/sano/work/forge-perf/case/16.nozzle_wys/README.md:75) と plan の変更ログには「場ノイズ床内」が残っています。

   **対案:** ラベルごとに実行コマンドと VERDICT を保存し、「未収束状態の回帰 PASS」と記載する。README・変更ログも合格理由を訂正する。

確認できた修正も明確です。`run_0461_frozen_roe_protect/{protect_cf0,protect_cf1}` は双方 **20/20セルで `roe` 不変、T=6000 K**。`run_0458` の `merged/new_cf0` は T・P・凝縮診断6量がビット一致し、double source 本体と double clamp 本体も旧コードとの同一性を確認しました。`run_0457` の cell H2O 回帰は **30/30 PASS**を再現しています。

**推奨は、`in_progress` を維持することです。** 優先順は、退避漏れの修正 → 移流単独の保存性検証と未解決事項の記録 → 判定基準・文書・証拠の統一 → 修正後バイナリの回帰とA10G実測です。

ファイルは変更しておらず、指摘は **plan 未反映**です。単体試験の再ビルドは行っていません。SASS の再取得は `cuobjdump` が一時ファイル作成を要求し、read-only 制約で失敗したため、命令列同一性までは独立確認できていません。

指摘数: Critical 0 / Major 4 / Minor 2
