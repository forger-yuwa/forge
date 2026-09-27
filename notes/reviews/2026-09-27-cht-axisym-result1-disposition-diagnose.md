# codex 諮問 (diagnose): cht-axisym-result1-disposition

- **brief**: [`notes/reviews/briefs/2026-09-27-cht-axisym-result1-disposition.md`](../../notes/reviews/briefs/2026-09-27-cht-axisym-result1-disposition.md)
- **plan**: [`plans/active/boundary-cht-axisymmetric-fem2d.md`](../../plans/active/boundary-cht-axisymmetric-fem2d.md)
- **date**: 2026-09-27
- **commit**: `0f1bab8b` (feature/cht-axisym-fem2d)
- **codex**: effort `high`, 1.6 min, rc=0
- **結論**: **残った手順表の保留記述を修正し、§5.1 #12 の限定・例外を維持して `done`・`accepted` へ移すことを推奨する；全面的な result 2 巡目は不要。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表:

| 指摘 | 採否・完了判断 | 根拠・対案 |
|---|---|---|
| **M1 / Major** | **採用・閉じてよい** | [conjugateWall.cpp:617](/home/sano/work/forge-cht/solver_density_cuda/conjugateWall.cpp:617) は部分欠落を状態復元前に拒否。[check_cht_balance.py:96](/home/sano/work/forge-cht/solver_density_cuda/tools/check_cht_balance.py:96) も同じ条件。評価器をメモリ上の HDF5 で独立実行し、属性有無の全組合せ×平面／軸対称 **16 件で期待との不一致 0**。C++ の実行結果は[保存記録](/home/sano/work/forge-cht/case/61.conjugate_annulus/vax4_contract_evidence/README.md:5)と整合する。追加修正は不要。 |
| **m1 / Minor** | **採用・一部未完了** | [case/61 README:59](/home/sano/work/forge-cht/case/61.conjugate_annulus/README.md:59) に「系列の組み方が決まるまで作らない」が残り、同文書103行の「決着・実行済み」と矛盾する。**この行を登録済みの連成格子感度・交差条件と run 一覧への参照に置き換える**。 |

結論: **残った手順表の保留記述を修正し、§5.1 #12 の限定・例外を維持して `done`・`accepted` へ移すことを推奨する；全面的な result 2 巡目は不要。**

第 1 仮説: M1 の契約検査漏れは解消し、今回の反映に残る不足は文書同期である。確度: 高  
　根拠: 上表のコード、16 件の独立確認、README 内の明示的な矛盾。  
　反証条件: 部分欠落した状態が現在の C++ で復元される、または3属性なしの旧平面状態が契約検査で拒否されること。

第 2・第 3 仮説: 無し。

判別 A/B: **評価器について実施済み、追加 run 不要**。平面 run に対し、A＝軸対称状態の3属性を保持、B＝`geometry` だけ削除。各1回、契約判定のみ確認。**A 拒否・B 受理なら穴が残存、両方拒否なら元の穴は解消**。旧版は前者、修正版は後者を再現した。

やらない方がよいこと: 文書修正の確認だけのために長時間の連成計算を再実行すること、既存の例外受け入れを全面 PASS に書き換えること。

呼び出し側の前提への異議: 「m1 は同期済み」は上記の残記述があるため受け入れない。平面2件の安全停止は契約通過の否定にはならないが、**停止原因が流体場未継承であること自体は今回独立検証していない**。

不足情報: `plans/README.md` は閲覧禁止のため同期内容を未確認。C++ の9件は保存記録とコードの照合であり、再実行していない。今回の推奨は **plan 未反映**。呼び出し側で §5.1 #12 に判断を記録し、移動時に索引・リンクを同期すること。
