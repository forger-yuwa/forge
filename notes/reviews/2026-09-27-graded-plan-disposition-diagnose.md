# codex 諮問 (diagnose): graded-plan-disposition

- **brief**: [`notes/reviews/briefs/2026-09-27-graded-plan-disposition.md`](../../notes/reviews/briefs/2026-09-27-graded-plan-disposition.md)
- **plan**: [`plans/active/axisymmetric-graded-grid-static-gas.md`](../../plans/active/axisymmetric-graded-grid-static-gas.md)
- **date**: 2026-09-27
- **commit**: `77e270c8` (feature/cht-axisym-fem2d)
- **codex**: effort `high`, 1.8 min, rc=0
- **結論**: **採否と修正方針は妥当なので、評価器テストの実行証跡の確認を含めて `result` レビューへ進んでよい。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表:

| 指摘 | 採否 | 根拠・対応評価 |
|---|---|---|
| M1 **Major**：逆向き熱流を合格させる | **採用。修正は妥当** | `9c0444f2` の `eval_static_hold.py:82` は `abs(q−q_ref)` に修正。判定式を独立検算し、正常／両壁反転／片側反転が新式では **PASS／FAIL／FAIL**、旧式では全部 PASS。反転時の誤差は **200 %**。 |
| m2 **Minor**：撤回記述の残存 | **採用** | 同コミットの `methods/boundary.md:500`、case README:114、plan 変更ログで訂正を確認。`plans/README.md` は閲覧禁止のため同期を独立確認していない。 |
| m3 **Minor**：連結判定の区間不整合 | **採用。記録の修正は妥当** | [連結判定記録](/home/sano/work/forge-cht/case/62.conjugate_disk/run_0018_hold_B_ext40k/CONVERGENCE_VERDICT_concat.txt:2) に入力2本、後段への **+20000**、**120000行・step 0–39999** を記載。保存された判定は **`OVERALL: ALL PASS`、8.0–8.3桁低下**。今回は入力CSVからの再判定までは行っていない。 |
| m4 **Minor**：同一質量の平衡解との混同 | **採用** | [plan:23](/home/sano/work/forge-cht/plans/active/axisymmetric-graded-grid-static-gas.md:23) に質量差 **30.1 %** と保証範囲を明記。静止保持の合格から質量保存まで主張しない処置でよい。 |

結論: **採否と修正方針は妥当なので、評価器テストの実行証跡の確認を含めて `result` レビューへ進んでよい。**

第 1 仮説: M1 の判定漏れは修正され、限定した閉鎖方針を妨げない。確度: 高  
  根拠: 上記の新旧判定式の検算と、[テスト:29–35](/home/sano/work/forge-cht/case/62.conjugate_disk/test_eval_static_hold.py:29) の評価器本体・終了コードを検査する構成。  
  反証条件: 評価器本体を実行すると符号反転ケースが PASS、または正常ケースが FAIL になる。

第 2・第 3 仮説: 無し。

判別 A/B: **同じ合成入力で旧／新の熱流束判定式だけを変更**し、既存テストの3時点を評価する（数秒、`forge` 不要）。→ 旧で反転2件を見逃し、新で正常 PASS・反転2件 FAIL なら修正を支持／新でも見逃すなら M1 は未解消。

やらない方がよいこと: この確認のための追加流体計算、原因特定済みという扱い、一様 IC からの CHT 起動への保証拡大。

呼び出し側の前提への異議: 「負例3件」は正確には**正常対照1件＋負例2件**。また、本レビューで確認した **3/3 は判定式の検算**であり、評価器本体の実行成功とは区別する。

不足情報: **Minor — 本体テストの実行証跡が今回の確認では不足**。実行は `test_eval_static_hold.py:27` の一時ディレクトリ作成で停止した。対案は、書き込み可能な呼び出し側環境で既存テストの出力・終了コードを `result` レビューに添付すること。ファイル変更なし、本判断は **plan 未反映**。
