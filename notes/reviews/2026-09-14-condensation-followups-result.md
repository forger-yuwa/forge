# codex レビュー: condensation-followups (result)

- **plan**: [`plans/active/condensation-followups.md`](../../plans/active/condensation-followups.md)
- **stage**: `result` (diff base `HEAD~1`)
- **date**: 2026-09-14
- **commit**: `e0ffca34` (feature/perf-3d-speedup)
- **codex**: effort `high`, 5.6 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M3/m2
- **focus**: F-cf9 (蒸気 c_p,v を CPG では physProp.cp に統一) の実装 diff と case/34 の A/B 結果に集中。空気キャリアで onset 不変・液滴モーメント 2.4 % 変化という判断の妥当性、accepted な condensation-air.md の検証を壊していないか、判定手順 (ノイズ床 + diff_res --factor 2, onset_analysis --series) の使い方が正しいかを見てほしい

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**  
case/34 の「格子点で判定した onset は同じ」「凝縮モーメント差はノイズ床を超える」は再現できました。  
ただし、比熱比の不整合、TP/H2O を含む検証漏れ、plan の未完了状態があり、accepted への移行は認められません。

1. **Major — `cp` の上書きで、Kantrowitz が使う `cp/cv` を壊しています。**

   [condensationTransport_d.cuh:38](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/condensationTransport_d.cuh:38) は空気 CPG に `cp=1008.7` を渡しますが、[condensationProperties_d.cuh:420](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/condensationProperties_d.cuh:420) は `cp` だけを更新し、N2 の `cv=742`、`R=296.8` を保持します。結果は `cp−cv=266.7≠R`、`gamma_v=1.35943` です。

   この比を核生成の double/float 両経路が使用します（[condensationSourceKernels_d.cuh:66](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/condensationSourceKernels_d.cuh:66)、[同:319](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/condensationSourceKernels_d.cuh:319)）。Kantrowitz の係数 `2(γ−1)/(γ+1)` は **0.333333→0.304678** に変わります。これは潜熱の傾き変更とは別の作用です。

   今回の空気 A/B は実効 `condKantrowitz=0` で、この問題を通りません（[forge_run.log:175](/home/sano/work/forge-perf/case/34.arthur_n2_nozzle/run_0105_gascp_air_node_new/forge_run.log:175)）。

   **対案:** Kirchhoff の傾き用の有効熱容量を核生成用の `cp/cv` から分離してください。F-cf9 の config 値は低温 `L/p_sat` に限定し、核生成には種固有の整合した `cp/cv/R` を渡す設計を推奨します。config 経由で `condKantrowitz=1` を通す double/float 試験も必要です。

2. **Major — 実際の変更範囲を、空気・純 N2 の node A/B だけでは検証できていません。**

   取得した `HEAD~1...HEAD` の diff では、TP の `gasCp` も **1038.8→0** に変わります。`condProps_make` の上書きが解除されるため、TP/H2O の実効 `cp` は **1038.8→1855**、`cp/cv` は **0.74546→1.33118** になります。

   [condensationSource_d.cuh:36](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/condensationSource_d.cuh:36) は負の θ をゼロにするので、直前版では mode 1 の非等温抑制が消えていました。今回の変更はそれを復活させます。**修正方向は正しいものの、TP を不変の対照として扱える変更ではありません。**

   また、共有物性を使う cell の F-cf9 回帰がありません。accepted な [condensation-air.md:184](/home/sano/work/forge-perf/plans/accepted/condensation-air.md:184) の node/cell 比較ゲートを、今回も満たす根拠は未取得です。

   **対案:** TP/H2O mode 1 の θ・核生成率を実カーネル経由で照合し、case/16 の場への影響を確認してください。case/34 は空気 cell A/B を追加し、旧計画の node/cell ゲートを再判定してください。これらを §5.1 に明示する必要があります。

3. **Major — `condensation-followups.md` 全体を完了扱いにできません。**

   実ファイルでも [§4:33](/home/sano/work/forge-perf/plans/active/condensation-followups.md:33) は「未着手」、§7 は未記述、[§6.1:74](/home/sano/work/forge-perf/plans/active/condensation-followups.md:74) は起票時の免除行だけです。F-cf9 の設計・検証計画を記した plan レビューはありません。

   [§8:82](/home/sano/work/forge-perf/plans/active/condensation-followups.md:82) は全項目の完了または別 plan への分離を要求していますが、多数の未着手項目が残っています。F-cf9 自身の「CPG carrier EOS の液相顕熱」という残件も、完了行の文章中に埋まっています。`check_plans.py` は PASS でしたが、これは `draft` の構造検査であり、完了条件の充足ではありません。

   **対案:** この集約 plan は active に残し、F-cf9 を独立 plan に分離してください。設計・影響範囲・追加検証・レビュー採否を記録し、EOS の残件には独立した未完了行を設けてください。

4. **Minor — 「液滴モーメント `Q_n` が最大 2.4 %」は測定量の表記が違います。**

   [methods/condensation.md:593](/home/sano/work/forge-perf/methods/condensation.md:593) は `Q_n` と書いていますが、[diff_vs_base.txt:11](/home/sano/work/forge-perf/case/34.arthur_n2_nozzle/run_0105_gascp_air_node_new/diff_vs_base.txt:11) の **2.429 % は `roQ0_0`** です。`Q0_0` は **2.11 %** です。どちらも尺度は `max|Δfield|/max|base field|` であり、局所相対差や出口値ではありません。

   3 反復の全ペアからノイズ床を再計算すると、保存済み JSON と全キーで一致しました。`--factor 2` の再実行も空気 **exit 1**、純 N2 **exit 0** を再現しています。「ノイズ床超」という判断は正当です。

   **対案:** 「保存モーメント `ρQ0` の最大差／基準場最大値が2.429 %、`Q0` は2.11 %」と明記してください。

5. **Minor — onset と収束について、判定の適用範囲を明記する必要があります。**

   再実行結果は次のとおりです。パスの基点は `case/34.arthur_n2_nozzle/`、索引は [README.md:41](/home/sano/work/forge-perf/case/34.arthur_n2_nozzle/README.md:41) です。

   | run | `check_convergence` | `onset_analysis --series` |
   |---|---|---|
   | `run_0105_gascp_air_node_base` | **NOT CONVERGED — plateau** | **STEADY** |
   | `run_0105_gascp_air_node_new` | **NOT CONVERGED — plateau** | **STEADY** |
   | `run_0106_gascp_n2_node_{base,new}` | 両方 **NOT CONVERGED — plateau** | 今回の再実行対象外 |

   空気の圧力 onset は両方 **2.197 in**、壁圧比 @3/4/5 in は旧 **1.2180/1.3569/1.4590**、新 **1.2185/1.3573/1.4595** でした。全対象 run の CSV の残差列に NaN/Inf はなく、凝縮 A/B・反復の全保存場も有限でした。

   一方、[onset_analysis.py:67](/home/sano/work/forge-perf/case/34.arthur_n2_nozzle/onset_analysis.py:67) の系列判定には、副指標の `g>10⁻⁴` onset とモーメントが含まれません。また onset は補間なしの格子点判定です。参考に隣接点間を線形補間すると圧力 onset は約 **2.1866→2.1813 in** と動きます。計画の **0.02 in** 許容幅内ですが、連続位置まで不変とはいえません。

   **対案:** 「未収束の準定常比較で、定義済みの格子点 onset ゲートを維持」と記述し、副 onset も系列判定に追加してください。現データは旧空気 plan の node onset 結論を覆していませんが、完全な回帰合格を意味しません。

**推奨:** 集約 plan は active に残し、F-cf9 を独立させて、まず熱容量の用途を分離し、次に TP/H2O・空気 cell の検証を補って再レビューしてください。今回は読み取りと再解析のみで、ファイル変更・plan 反映はしていません。

指摘数: Critical 0 / Major 3 / Minor 2
