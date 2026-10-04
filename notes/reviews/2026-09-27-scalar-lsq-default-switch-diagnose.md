# codex 諮問 (diagnose): scalar-lsq-default-switch

- **brief**: [`notes/reviews/briefs/2026-09-27-scalar-lsq-default-switch.md`](../../notes/reviews/briefs/2026-09-27-scalar-lsq-default-switch.md)
- **plan**: [`plans/active/gradient-scalar-lsq-unification.md`](../../plans/active/gradient-scalar-lsq-unification.md)
- **date**: 2026-09-27
- **commit**: `0e40da0a` (feature/sern-design)
- **codex**: effort `high`, 1.4 min, rc=0
- **結論**: **現案のままの既定切り替えは保留し、SERN の学習採否を実効 `scalarGradient` でも判定するゲートを先に追加する。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論: **現案のままの既定切り替えは保留し、SERN の学習採否を実効 `scalarGradient` でも判定するゲートを先に追加する。**

採否表:

| 対象 | 採否 | 理由・対案 |
|---|---|---|
| node 省略時 `lsq`、cell は `gg`、明示指定を維持 | 採用 | 現行の分岐に沿う変更。提示された起動エコー・1 step 比較で確認する。 |
| ユーザ決定「B」 | 採用 | [plan:120](/home/sano/work/forge-sern-design/plans/active/gradient-scalar-lsq-unification.md:120) に記録あり。#2g 完了を再び前提にはしない。 |
| `FLAG_POLICY` 更新だけで DB を分離 | **却下・Major** | 実効設定を検査していない。各作動点の実効値を保存し、node の新方針では `lsq` と確認できる評価だけを学習対象にする。 |

第 1 仮説: **新しい `FLAG_POLICY` を持つ旧 `gg` 評価が、新方針の学習対象へ混入する。** 確度: **高**（受け入れるコード経路を確認。実キャンペーンでの発生は未確認）。

  根拠:
  - [runner_sern.py:940](/home/sano/work/forge-sern-design/design/forge_design/evaluate/runner_sern.py:940) は Python の定数を無条件に付与する。バイナリの実効値から決めていない。
  - [driver_sern.py:190](/home/sano/work/forge-sern-design/design/forge_design/opt/driver_sern.py:190) も評価行へ定数を付与し、[同:236](/home/sano/work/forge-sern-design/design/forge_design/opt/driver_sern.py:236) の学習選別は `PASS` と定数一致だけを見る。
  - 実際の `_XF` を AST から取り出して読み取り専用で実行した。新 policy の合成行に実効値 `gg`・`lsq`・不明を与えると **3/3 行を採用**、旧 policy なら **0/3 行**。これは選別処理の試験であり、CFD の実測ではない。
  - 実効値の DB 保存自体も [plan:93](/home/sano/work/forge-sern-design/plans/active/gradient-scalar-lsq-unification.md:93) の #2h に含まれるが、今回の実装案では抜けている。

  反証条件: 新 policy でも実効 `gg`・不明の node 評価が学習前に拒否され、確認済み `lsq` だけが採用されること。既存の稼働中キャンペーンにもその判定が適用される必要がある。

第 2・第 3 仮説: 追加しない。

判別 A/B: **forge は回さず、評価結果取り込みから `_XF` までの試験を各 1 回行う。** node・新 policy・数値ゲート `PASS` を固定し、起動記録の実効 `scalarGradient` だけを `gg`／`lsq` に変える。見る量は保存された実効値と学習採用件数。
→ **A: `gg` は 0 件、`lsq` は 1 件**なら、この経路での混入仮説を棄却。
→ **B: 両方 1 件**なら、日付識別だけでは分離できないことが確定。現行 `_XF` 単体は B。

やらない方がよいこと: `FLAG_POLICY` の更新と警告だけで、走行中キャンペーンまで切り替わったと扱うこと。既存 Python プロセスは import 済み定数を保持するため、切り替え時にはキャンペーンを停止し、更新後に再起動して新ゲートを適用する。

呼び出し側の前提への異議: 「前日までの評価行を外す」は日付判定ではなく、**policy 文字列の一致判定**である。「DB は分けられる」は未成立。ユーザ決定「B」と Phase 1 の合格判断を覆す指摘ではない。

不足情報: 稼働中キャンペーンの停止・再起動予定と使用バイナリの対応。Phase 1 の run 数値・VERDICT は今回再検証しておらず、合格記録を独立に追認したものではない。**plan 未反映**（読み取り専用。呼び出し側が §4.4・§5.1 #2h・§6 に反映）。
