# codex レビュー: boundary-cht-conjugate-benchmarks (result)

- **plan**: [`plans/active/boundary-cht-conjugate-benchmarks.md`](../../plans/active/boundary-cht-conjugate-benchmarks.md)
- **stage**: `result` (diff base `feature/cht-axisym-graetz`)
- **date**: 2026-09-30
- **commit**: `609ef2ca` (feature/cht-conjugate-benchmarks)
- **codex**: effort `high`, 4.6 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M5/m1
- **extra**: `case/64.conjugate_pipe_wall/README.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**NO-GO**  
A に限定しても、積分量と不確かさの評価に修正が必要です。準定常判定には NaN を PASS にする欠陥があり、保存済み証拠にも版の不整合があります。  
A/C の分割は妥当ですが、現状の「登録許容内」を確定結果として accepted に移すことは承認できません。

指定の `git diff feature/cht-axisym-graetz...HEAD -- solver_density_cuda design methods procedures` は空でした。ソルバ・既定値の変更による回帰はありません。追加された参照解の自己検査は再実行で **`VERDICT: PASS`**。製造解の誤差列、軸上散逸、熱収支は plan の数値を再現しました。ただし、これは以下の評価器・証拠の問題を解消しません。

1. **Major — `Q_up/Q_tot` が登録した比になっておらず、分母の不確かさも欠落している。**

   **根拠:** [eval_conj.py:212](/home/sano/work/forge-cht/case/64.conjugate_pipe_wall/eval_conj.py:212) で参照解の `Qtot` を取得し、221–222 行では forge・参照の両方の上流熱量をその同じ値で割っています。算出している差は `(Q_up,forge − Q_up,ref)/Q_tot,ref` です。`U` も上流熱流束の積分だけで、各感度計算による総入熱の変化を含みません。

   `case/64.conjugate_pipe_wall/run_0006_a1_r32/` では、保存された界面履歴の総入熱は約 **0.016123878 W/rad**、参照は約 **0.0161472 W/rad**。同じ forge 上流熱量でも、現行評価の **0.153254** と自身の総入熱で割った **0.153476** は異なります。

   **対案:** forge と参照でそれぞれ外面 Robin 入熱を計算し、各々の比を比較する。参照格子・写像・微分・領域変更ごとに**比そのもの**を再計算して `U` を求め、A 全 run を再評価してください。

2. **Major — 領域切断の試験が事前登録と異なり、未評価領域の不確かさをゼロにしている。**

   **根拠:** [plan:110](/home/sano/work/forge-cht/plans/active/boundary-cht-conjugate-benchmarks.md:110) は A の上流長を **1.5 倍**にする登録です。しかし [eval_conj.py:91](/home/sano/work/forge-cht/case/64.conjugate_pipe_wall/eval_conj.py:91) は入口を **−80R から約 −60R へ短縮**し、新入口には forge の温度を与えます。さらに [同:203](/home/sano/work/forge-cht/case/64.conjugate_pipe_wall/eval_conj.py:203) は切り落とした節点の感度を `nan_to_num` でゼロにします。それで「全長」の `U` を満たしたとは言えません。

   §5.1 #6b に短縮試験の記録はありますが、§4.6 との矛盾と代替試験の妥当性は未解決です。A2 r32 の積分量は **U=0.0016515、上限=0.0016667** と余裕が約 **1.5e−5**しかなく、この問題は合否に直結します。

   **対案:** 登録した領域延長試験を実施し、延長部分の固定流れ場の定義も記録する。評価できない領域をゼロ扱いせず、必要な `U` が得られるまで判定不能としてください。

3. **Major — 準定常評価が NaN を PASS にし、指定の正式判定経路も使っていない。**

   **根拠:** [series_conj.py:83](/home/sano/work/forge-cht/case/64.conjugate_pipe_wall/series_conj.py:83) の `judge()` は有限性を確認せず、`drift > D` と `fl > D` の件数だけで判定します。この関数をそのままメモリ内で実行し、末尾に NaN を一つ入れると、**「最大 nan、閾値超え 0、合格 True」**を再現しました。

   `check_quasisteady.py` はパスを定義するだけで呼ばれていません。また [eval_conj.py:173](/home/sano/work/forge-cht/case/64.conjugate_pipe_wall/eval_conj.py:173) の検査は存在する壁節点の `iface_ok` が中心で、§4.4 が要求する期待節点集合との照合がありません。欠落節点を評価対象から黙って落とせます。

   **対案:** 全系列・座標・規格化尺度の有限性と、期待節点集合の完全一致を必須にする。登録尺度による判定を正式ツールへ統合し、NaN/Inf・節点欠落・重複・ゼロ尺度の負例を追加してください。現在の実測に NaN があったという指摘ではなく、ゲートの実装欠陥です。

4. **Major — 本番結果を独立再検証できず、ログと入力 CSV の対応も崩れている。**

   **根拠:** A の採用候補を含む本番6本すべてで、ローカルに `residual_history.csv`、流体・壁・固体の保存場、`CONVERGENCE_VERDICT.txt` がありません。`check_convergence.py` を実行した結果は、全件 **`NO residual_history.csv`** でした。保存された `CONVERGENCE_CHECK.txt` の PASS は読めますが、残差全系列・保存場の NaN/Inf を独立確認できません。

   さらに `case/64.conjugate_pipe_wall/run_0006_a1_r32/SERIES_CONJ.txt` は温度・熱流束とも **255 節点**と記録していますが、対応する `series_conj.csv` は温度 **164 列**・熱流束 **127 列**です。A 本番6本すべてで同様の不一致があります。旧 CSV の積分量は正式ツールで **`OVERALL: ALL STEADY`** を再現しましたが、全節点の証拠にはなりません。

   A1 r16 も、`EVAL_CONJ.txt` の加熱区間温度差 **1.1905%**に対し、同居 CSV からの再集計は **1.24856%**で、4水準 A/B の値です。[eval_conj.py:243](/home/sano/work/forge-cht/case/64.conjugate_pipe_wall/eval_conj.py:243) は水準数によらず同じ CSV 名へ出力します。

   **対案:** 原出力をレビュー可能な場所へ揃え、評価条件・コード版・入力ハッシュ別に成果物を保存する。3水準評価、4水準 A/B、全節点準定常評価を分離し、ログと CSV を同時に再生成してください。

5. **Major — 「固体が効いている」の登録ゲートが未完了。**

   **根拠:** [solid_effect.py:33](/home/sano/work/forge-cht/case/64.conjugate_pipe_wall/solid_effect.py:33) は最終スナップショットだけを読み、中間参照格子との値を表示して終了します。`ΔT_s`・`Q_ax` 自体の不確かさ、時系列の定常性、[plan:152](/home/sano/work/forge-cht/plans/active/boundary-cht-conjugate-benchmarks.md:152) の「効果 ≥ 5U」の判定はありません。

   また A2 r32 の `SOLID_EFFECT.txt` は **ΔT_s=0.2588 K**。温度許容幅は **0.01×7.7719=0.077719 K**なので約 **3.33 倍**です。§5.1 #6c の「許容より1〜2桁大きい」は、この指標について誤りです。

   **対案:** 両指標の時系列と不確かさを評価し、`効果/U`・`効果/許容幅` を数値で判定する。未実施分を §5.1 の残作業として明示してください。

6. **Minor — run 索引・現在仕様・残作業表が最新結果に同期していない。**

   **根拠:** [case/64 README:15](/home/sano/work/forge-cht/case/64.conjugate_pipe_wall/README.md:15) は r16 を単に判定不能、A2 r64 を PASS としています。保存ログの最新判定は次のとおりです。

   | `case/64.conjugate_pipe_wall/` 配下の run | `EVAL_CONJ.txt` の判定 |
   |---|---|
   | `run_0005_a1_r16` | FAIL（一部判定不能） |
   | `run_0006_a1_r32`、`run_0007_a1_r64`、`run_0009_a2_r32` | PASS |
   | `run_0010_a2_r64` | 判定不能 |
   | `run_0011_a2_r16_df20` | FAIL（一部判定不能） |

   [methods/boundary.md:350](/home/sano/work/forge-cht/methods/boundary.md:350) は固体抵抗が効く共役を未検証とする旧記述のまま、[plans/README.md:28](/home/sano/work/forge-cht/plans/README.md:28) も分割前の説明です。§5.1 には終了済み作業と残作業が混在しています。C の未解決事項を後継 plan に残した点は適切です。

   **対案:** 再評価後の判定へ索引を更新し、現在仕様には FP64・node・流れ場固定という検証範囲を明記する。残作業表を更新し、移動時に後継 plan・case README のリンクも同期してください。

**推奨は、A を active に残し、評価器を修正して原出力から再評価した後に result レビューをやり直すことです。** 優先順は、原証拠の整備 → 指摘1〜3の修正 → 固体効果のゲート → 文書同期です。C は後継 plan で未完了を維持してください。ファイルは変更しておらず、本レビューの指摘は **plan 未反映**です。

指摘数: Critical 0 / Major 5 / Minor 1
