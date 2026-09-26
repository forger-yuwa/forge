# codex 諮問 (diagnose): cht-result8-disposition

- **brief**: [`notes/reviews/briefs/2026-09-26-cht-result8-disposition.md`](../../notes/reviews/briefs/2026-09-26-cht-result8-disposition.md)
- **plan**: [`plans/active/boundary-conjugate-heat-transfer.md`](../../plans/active/boundary-conjugate-heat-transfer.md)
- **date**: 2026-09-26
- **commit**: `fe8142c8` (feature/cht-phase2-fem2d)
- **codex**: effort `high`, 2.0 min, rc=0
- **結論**: **負例試験の判定を厳密化して再実行し、期待結果を確認後、限定付き成果として `status: done` → `accepted/` へ移す。全面的な result 9 巡目は不要。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表:

| 指摘 | 採否 | 確認結果 |
|---|---|---|
| M1・Major | 採用 | 前壁 `q_eff`・後壁 `-q_eff` に変更され、符号付き収支も検査する。修正は妥当。 |
| M2・Major | 採用 | `load_node_log()` の共通化は妥当。変更コードをメモリ上で実行し、正常入力の受理、座標 NaN・時刻不一致・ID 重複の拒否を確認した。 |
| M3・Major | 採用 | §4.11(a) の比較値訂正と、他ケースへの一律「1 K 未満」保証の撤回を確認した。 |
| m4・Minor | 採用 | `HF[0]`／`HF[1]` の区別、固体 CFL、温度抽出元の訂正を確認した。 |
| m5・Minor | 採用 | `methods` の追記を確認。case/48 の同期は #104 の記録に依拠する。索引同期は移動時に必要。 |

**修正が不十分な点は1つ：Minor — 負例試験が異常終了を成功扱いする。**

[test_v6p_negative.py:71](/home/sano/work/forge-cht/case/58.conjugate_slot/test_v6p_negative.py:71) の `("PASS" in v) == expect_pass` は、負例で `(no VERDICT)` が返っても成功になる。`run()` も終了コードを検査していない。この判定式の誤受理を確認した。

対案は、終了コードと期待する判定を両方照合すること。正例は `PASS`、符号反転は `FAIL`、節点ログ破損は `REFUSED` を要求し、判定欠落・例外終了は試験失敗にする。

結論: **負例試験の判定を厳密化して再実行し、期待結果を確認後、限定付き成果として `status: done` → `accepted/` へ移す。全面的な result 9 巡目は不要。**

第 1 仮説: 残る問題は負例試験の判定強度であり、M1〜M3 の修正方針を再設計する必要はない。確度: 高  
  根拠: commit `a73ad11f` の差分、M2 のメモリ内試験、上記の誤受理。  
  反証条件: 厳密化した試験で、正常入力が拒否される、または破損入力が所定の判定で拒否されない。

第 2・第 3 仮説: 無し。

判別 A/B: 負例入力を固定し、評価器の終了だけを「正常な拒否／VERDICT なしの異常終了」に切り替えて各1回実行する。→ A の正常拒否だけ試験成功なら修正成立／B の異常終了も試験成功なら穴が残る。CFD の追加計算は不要。

やらない方がよいこと: V2 の `FAIL` を `PASS` に変更すること、V6′ の固定深部帯の承認を全域収束へ拡張すること。

呼び出し側の前提への異議: 「6/6」だけでは意図した拒否を保証しない。閉じる際は、V2 のユーザによる例外受け入れ、V6′ の範囲限定、放射の保留を維持する。

不足情報: 実 run を使う負例試験は、read-only 制約で一時ディレクトリを作れず開始できなかった。修正後4 run の再評価結果は #104 の記録に依拠し、今回は独立再実行していない。

ファイル変更なし。plan 未反映。
