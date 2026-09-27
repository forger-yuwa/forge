# codex 諮問 (diagnose): graded-close

- **brief**: [`notes/reviews/briefs/2026-09-27-graded-close.md`](../../notes/reviews/briefs/2026-09-27-graded-close.md)
- **plan**: [`plans/active/axisymmetric-graded-grid-static-gas.md`](../../plans/active/axisymmetric-graded-grid-static-gas.md)
- **date**: 2026-09-27
- **commit**: `1b113c86` (feature/cht-axisym-fem2d)
- **codex**: effort `high`, 2.2 min, rc=0
- **結論**: **追加の forge run は投入せず、40000 step と CHT 起動手順の過大な主張を訂正した閉鎖方針を、まず `plan` レビューに回す。**
- **extra**: `case/62.conjugate_disk/README.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（Critical 0 / Major 3 / Minor 1）

| 重大度 | 対象・採否 | 根拠と対案 |
|---|---|---|
| **Major** | 「40000 step で既に静止、以後不変」は**却下** | `run_0020_hold_A_ext100k/static_hold_series.csv` を `check_quasisteady.py` の `classify_series` で再判定すると、累積 **30100–40000** は **67 STEADY / 4 DRIFTING / 1 OSCILLATING**。`qerr_hot`・`qerr_cj` の drift はそれぞれ **256.1 / 252.3 %/tail**。対案：**「40000 時点では振幅閾値内。登録窓 90000–100000 で全条件合格」**と書く。40000 を連成開始の推奨値にしない。 |
| **Major** | 非連成結果を根拠に「warmup 延長で CHT を解決」は**要再検証** | 非連成の合格は、連成開始後の安定性を保証しない。[発注元 plan:217](/home/sano/work/forge-cht/plans/accepted/boundary-cht-axisymmetric-fem2d.md:217) は純伝導 IC に限定している。対案：保証範囲の据え置きは**採用**し、実証済みの起動手順は純伝導 IC とする。一様 IC＋長い warmup は未検証候補として区別する。 |
| **Major** | `result` レビュー 1 回だけで閉じる案は**却下** | [レビュー手順:19](/home/sano/work/forge-cht/procedures/codex-review.md:19) は `diagnose` と `plan/result` を別扱いとし、[AGENTS.md:314](/home/sano/work/forge-cht/AGENTS.md:314) も両レビューの維持を要求する。対案：修正した閉鎖方針の **`plan` レビューを補完し、その後 `result` レビュー**。事前レビューの欠落は履歴に残し、事後レビューを事前実施扱いにしない。 |
| **Minor** | 限定した結論で閉じ、機構追究を保留する案は**採用** | 登録窓の合格は確認できる。一方、[plan:20](/home/sano/work/forge-cht/plans/active/axisymmetric-graded-grid-static-gas.md:20) の目的、§4、§8 には旧判断が残る。対案：完了対象を「登録構成での非連成静止保持の確認と、検証済み起動条件の文書化」に揃える。原因未特定・既知 slip 欠陥未解決を残し、コード変更がないなら §6 の回帰再実行を不要とする理由も記録する。 |

結論: **追加の forge run は投入せず、40000 step と CHT 起動手順の過大な主張を訂正した閉鎖方針を、まず `plan` レビューに回す。**

第 1 仮説: **登録された FP64・非連成構成では、一様 IC から生じた大振幅状態は減衰し、静止保持の許容状態に到達する。** 確度: **高**（減衰機構の確度ではない）  
  根拠: `case/62.conjugate_disk/run_0020_hold_A_ext100k/` の累積 **90000–100000** を再集計すると、max|U| ≤ **3.63756e-4 m/s**、熱流束誤差 hot ≤ **0.00245464 %** / cj ≤ **0.00278025 %**、相対圧力差 ≤ **2.77943e-7**。CSV 全800行に NaN/Inf はない。準定常はツール関数による再判定で **`OVERALL: ALL STEADY`（72系列）**。残差は保存記録 [CONVERGENCE_VERDICT_concat.txt:3](/home/sano/work/forge-cht/case/62.conjugate_disk/run_0020_hold_A_ext100k/CONVERGENCE_VERDICT_concat.txt:3) が累積 0–99999 に **`PASS (converged)` / `OVERALL: ALL PASS`**、有効な全保存量で **13.3–14.3 桁低下**。  
  反証条件: 同一設定・同一節点の元データから再構成した登録窓で、閾値超過または非定常が検出されること。連結区間の設定不一致が判明した場合も、この根拠は成立しない。

第 2・第 3 仮説: **追加しない。** 非一様格子で減衰が遅い機構も、既知 slip 欠陥の寄与も、今回の資料では特定できない。

判別 A/B: **追加計算 0 step。既存 CSV の評価窓だけを変える再解析**とする。A＝累積 20000–40000、B＝80000–100000 を入力し、双方 `tail=0.5, drift=0.001, osc=0.001`、同じ72系列を判定する。事前の読み方は、**A・B とも全系列 STEADY なら40000到達説を支持、A 不合格・B 合格なら後半の合格だけを採用**。今回の再解析は後者だった。これは40000到達説の検査であり、登録済み最終判定窓の変更ではない。

やらない方がよいこと: **一様 IC・warmup 60000 の CHT run を、閉鎖の必須条件として追加すること。** 現在の保証範囲を広げる必要がないため、今は実施価値が低い。また「起動手順で欠陥を修正した」「FP32 でも成立する」と書かない。

呼び出し側の前提への異議: **「振幅が小さい」「準定常」「連成開始後も安全」は別の判定である。** 40000時点は最初の条件しか裏付けていない。`run_0014` の失敗が大きな過渡と整合することは認めるが、それだけで長い warmup の CHT 成功を確定できない。

不足情報: 指定された閲覧範囲では、残差連結の元データ・継続前後の実効設定差分・集計元の全保存場は独立照合していない。`result` レビューではその来歴を確認する。本回答は **plan 未反映**であり、反映先は当該 plan の **§4・§5.1 #4・§6・§8**。
