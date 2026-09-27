# codex レビュー: axisymmetric-graded-grid-static-gas (plan)

- **plan**: [`plans/active/axisymmetric-graded-grid-static-gas.md`](../../plans/active/axisymmetric-graded-grid-static-gas.md)
- **stage**: `plan`
- **date**: 2026-09-27
- **commit**: `809fee38` (feature/cht-axisym-fem2d)
- **codex**: effort `high`, 4.1 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M1/m3
- **focus**: 事後の補完レビュー: §4/§6 は diagnose で決めて plan 段レビューを回していなかった。現時点の閉鎖方針 (§1 完了の範囲、§8 決着) が妥当か、閉鎖前に足りないものを指摘すること
- **extra**: `case/62.conjugate_disk/README.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
登録された FP64・node 構成に限定し、原因未特定のまま閉じる方針は妥当です。  
追加の `forge` 計算は不要ですが、熱流束の判定漏れと記録の不整合を修正してから `result` レビューへ進むべきです。

1. **Major — 熱流束の合格判定が、逆向きの熱流を通してしまう。**

   根拠: [`eval_static_hold.py:79`](/home/sano/work/forge-cht/case/62.conjugate_disk/eval_static_hold.py:79) は `abs(abs(q)-Q_REF)` を評価しています。符号付き系列の `STEADY` 判定を追加しても、向きの正しさは検査されません。[符号規約](/home/sano/work/forge-cht/methods/boundary.md:354) に従えば期待値は hot が **−120.5**、cj が **+120.5 W/m²** です。両者を逆転させた定数系列で、現行の誤差は **0 %**、準定常ツールは **`STEADY`**、正しい符号付き誤差は **200 %** になることを確認しました。

   対案: §6 と評価器を `max|q−q_ref|/120.5 ≤0.005` に揃え、符号反転の負例を追加してください。今回の保存済み66節点系列は正しい符号で、再集計した最大誤差も B **0.005653 %**、A **0.002780 %**。**既存の合格結果は覆りません。再計算ではなく評価器の修正と再判定で足ります。**

2. **Minor — 撤回済みの結論が、現在の結論と併存している。**

   根拠: [plan:110](/home/sano/work/forge-cht/plans/active/axisymmetric-graded-grid-static-gas.md:110) に「累積40000で到達」が残っています。[`methods/boundary.md:500`](/home/sano/work/forge-cht/methods/boundary.md:500) は「最終到達状態は未確認」と「累積100000で静止保持に届く」を同時に記載しています。[case README:114](/home/sano/work/forge-cht/case/62.conjugate_disk/README.md:114) の「加熱なので圧力は上がるはず」も、§6で棄却した判断です。[計画一覧:27](/home/sano/work/forge-cht/plans/README.md:27) の「閉包欠損は棄却済み」も、正確には「今回の大振幅過渡の支配要因として支持されなかった」です。

   対案: 旧記述には撤回・訂正先を付け、現行結論を「登録窓90000–100000で合格」「減衰機構は未特定」「CHTの実証済み起動は純伝導IC」に統一してください。§5.1には、評価器修正・文書同期・レビュー記録補完・`result` レビューを閉鎖前の残作業として明記してください。

3. **Minor — Bの保存された連結判定は、宣言した終了区間と一致しない。**

   根拠: [`CONVERGENCE_VERDICT_concat.txt:3`](/home/sano/work/forge-cht/case/62.conjugate_disk/run_0018_hold_B_ext40k/CONVERGENCE_VERDICT_concat.txt:3) は最終 step **39979**、入力は一時パスです。保存済みの前段・継続段を全件連結すると **0–39999、40000行** になり、末尾残差も記録と異なります。

   対案: 完全な区間で判定記録を更新し、入力2本・stepオフセット・使用コマンドを残してください。今回、書き込みを避けて同じ `check_convergence.py` の判定関数に連結データを渡した結果は **`PASS`、8.0–8.3桁低下**でした。これは合否の誤りではなく、証跡の不整合です。

4. **Minor — 「静止保持の合格」と「同一質量の平衡解」の違いを、実測付きで残すべき。**

   根拠: 最終保存場を双対体積×双対重心半径で積分すると、Bの `run_0018_hold_B_ext40k/res_20000.h5` は平均圧力 **1051.736 Pa**・質量 **1.01841e−8 kg/rad**、Aの `run_0020_hold_A_ext100k/res_80000.h5` は **735.618 Pa**・**7.12309e−9 kg/rad**。終端スナップショットの質量は約 **30.1 %** 異なります。

   対案: §1・§8に「合格は速度・熱流束・相対圧力差の登録基準についてであり、同一質量の平衡解への到達を保証しない」と追記してください。[局所擬似時間の実装](/home/sano/work/forge-cht/solver_density_cuda/cuda_forge/setDT_d.cu:362) と整合する限界であり、この差だけから境界漏れを断定したり、今回の閉鎖条件に保存機構の修正を追加したりする必要はありません。

**推奨は、保証範囲を据え置き、上記を修正して閉鎖することです。** 優先順は **1 → 2 → 3 → 4 → `result` レビュー**です。

独立再確認では、`case/62.conjugate_disk/` 配下の以下を確認しました。

| run | 判定区間 | 再確認結果 |
|---|---|---|
| `run_0018_hold_B_ext40k` | 残差0–40000、準定常30000–40000 | 残差 **PASS**、72系列 **ALL STEADY** |
| `run_0020_hold_A_ext100k` | 残差0–100000、準定常90000–100000 | 残差 **PASS**、72系列 **ALL STEADY** |
| `run_0019_disk_r32g_condic` | 保存履歴全体、準定常は登録設定 | 残差 **PASS**、66系列 **ALL STEADY** |

非一様メッシュも **`VERDICT: PASS`**（AR最大4.6、skew最大0）。成果物の索引は [case READMEの計算run一覧](/home/sano/work/forge-cht/case/62.conjugate_disk/README.md:102) です。中間場の一部はローカルにないため、全期間の場からCSVを再生成する照合までは行っていません。

`hoop` 修正・既知の `slip` 欠陥修正との重複実装は不要です。ソルバを変更しない今回の限定閉鎖に、cell・FP32・周期・軸上の追加回帰を要求する理由もありません。ファイル変更は行っておらず、本指摘は **plan未反映**です。

指摘数: Critical 0 / Major 1 / Minor 3
