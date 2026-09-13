# codex レビュー: condensation-vapour-cp-source (result)

- **plan**: [`plans/active/condensation-vapour-cp-source.md`](../../plans/active/condensation-vapour-cp-source.md)
- **stage**: `result` (diff base `HEAD~2`)
- **date**: 2026-09-14
- **commit**: `e7ff51d3` (feature/perf-3d-speedup)
- **codex**: effort `high`, 5.1 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M2/m2

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**

`kirchhoffCpv` の分離は設計どおりで、今回の diff に新たな符号・単位・境界処理の誤りは見つかりませんでした。
主要な A/B 数値とノイズ床判定は再現できましたが、検証の表現・判定対象・液相比熱の説明を修正してから移すべきです。

1. **Major — H2O の「一致」は、300 step の回帰比較に限定すべきです。**

   根拠: [plan §6](/home/sano/work/forge-perf/plans/active/condensation-vapour-cp-source.md:82) は「ノイズ床以内で一致」を合格条件にしています。しかし、`case/16.nozzle_wys/run_0466_h2ocp_node2d_{pre,cur}` と `run_0467_h2ocp_cell2d_{pre,cur}` の再判定は、すべて **`NOT CONVERGED`**。保存場は各２枚で、`check_quasisteady.py --quantity pmax` も **`TRANSIENT-UNSETTLED`** でした。cell では乱流・凝縮残差も下降途中です。

   一方、`diff_res.py` の node/cell **exit 0 は再現**しました。cell の合成床も、旧・現行それぞれの内部反復差の最大値であり、A/B 差の混入はありません。

   **対案:** §6・§8.1 を「未収束の同一反復数で、旧実装との差が反復ノイズ許容内」と修正し、定常解の一致とは区別してください。この回帰目的だけなら、長時間計算を追加する必要はありません。

2. **Major — onset の合格条件と、自動判定の対象が揃っていません。**

   根拠: [onset_analysis.py:69](/home/sano/work/forge-perf/case/34.arthur_n2_nozzle/onset_analysis.py:69) の返却値と同ファイルの `TOL` は、圧力基準 onset だけを判定します。plan が合格条件に含めた **`g>10^-4` onset は `--series` の判定対象外**です。

   また、`case/34.arthur_n2_nozzle/run_0107_gascp_air_cell_base` の圧力 onset は、8000–12000 step で **2.20680–2.22626 in** と揺れています。最終場の「2.21 in が同じ」は再現しますが、[methods の「onset は動かない」](/home/sano/work/forge-perf/methods/condensation.md:598) は精度を言い過ぎています。

   **対案:** 両 onset を時系列判定に登録し、閾値・窓・VERDICT を保存してください。「変化なし」は「格子判定と時間変動の分解能内」に限定すべきです。今回、抽出系列をメモリ上で `check_quasisteady.classify` に渡した追加照合では、空気 node/cell の両 onset・出口湿り度・壁圧比は **`STEADY`** でした。結果を覆す指摘ではなく、合格根拠の欠落です。

3. **Minor — 残作業表の実効液相比熱が、現行実装と違います。**

   根拠: [plan §5.1 #1](/home/sano/work/forge-perf/plans/active/condensation-vapour-cp-source.md:72) は `1014.7`、実効液比熱 `1976` としています。しかし、[dependentVariables_d.cu:245](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/dependentVariables_d.cu:245) と [EOS:103](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/condensationEOS_d.cuh:103) に空気 config を代入すると、70 K 未満では以下です。

   - `cv = 1008.7 / 1.4 = 720.5`
   - `cv + Rw = 1017.3`
   - 現行の実効液比熱 = `1017.3 − (1008.7 − 2000) = 2008.6 J/(kg·K)`

   **対案:** §4・§5.1・`methods/condensation.md` を訂正してください。残る不整合は約 **+0.43%** で、「自己整合」は厳密な熱力学的整合を意味しない、と明示すべきです。

4. **Minor — コメントに、今回排除した誤実装の説明が残っています。**

   根拠: [condensationProperties_d.cuh:63](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/condensationProperties_d.cuh:63) は依然として「Kantrowitz に使う」「`cp` を上書きする」と記載。[test_cond_air.cpp:172](/home/sano/work/forge-perf/solver_density_cuda/tests/unit/test_cond_air.cpp:172) も「carrier では渡さない」としています。

   **対案:** 「CPG carrier を含め、`kirchhoffCpv` のみ更新」に統一してください。H2O の `kirchhoffCpv` は実際には未使用なので、その点も明記してください。

推奨は、**現実装を維持し、上記順に検証記録と文書を修正して accepted へ移すこと**です。空気の `ρQ0` 差約2.43%／2.92%、純 N2 の床内判定は再現し、対象25 run・193保存場の `VALUE/*` に NaN/Inf はありませんでした。run 索引は [case/34 README](/home/sano/work/forge-perf/case/34.arthur_n2_nozzle/README.md:41) と [case/16 README](/home/sano/work/forge-perf/case/16.nozzle_wys/README.md:54) にあります。

単体「８本全 PASS」は完全には再確認できませんでした。既存バイナリのCPU試験５本は通りましたが、GPU試験は CUDA driver/runtime 不整合で再現不能でした。ファイルは変更しておらず、本レビューの指摘は **plan 未反映**です。

指摘数: Critical 0 / Major 2 / Minor 2
