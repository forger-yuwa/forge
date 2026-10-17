# codex 諮問 (diagnose): orchestrator-approval-design

- **brief**: [`notes/reviews/briefs/2026-10-17-orchestrator-approval-design.md`](../../notes/reviews/briefs/2026-10-17-orchestrator-approval-design.md)
- **plan**: [`plans/active/tooling-design-orchestrator.md`](../../plans/active/tooling-design-orchestrator.md)
- **date**: 2026-10-17
- **commit**: `9c52fdd4` (feature/card-local-mesh)
- **分類**: `ok`
- **backend**: `codex` — 選んだ理由: 既定 (auto: codex を試し、unavailable のときだけ claude)
- **要求したモデル**: `既定 (~/.codex/config.toml)`
- **実効のモデル**: `gpt-6-astra` (codex のログの見出し)
- **CLI の版**: `codex-cli 0.162.1`
- **依頼側との関係**: 不明 (`--requester-model` 未指定)
- **permission_denials**: 対象外 (codex は read-only サンドボックス)
- **実行**: effort `high`, 5.3 min, rc=0
- **試行のログ**: `notes/reviews/2026-10-17-orchestrator-approval-design-diagnose.codex.log` (codex: ok)
- **結論**: §4.7.1 を、承認する実効条件と変更可能な試行入力を分離し、判定環境と全必須中間判定まで固定する契約へ修正する。
- **extra**: `design/forge_design/feedback/wt_qualify.py`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

**方向は採用しますが、§4.7.1 は修正してから実装すべきです。** 主な穴は、承認する条件の不足、カード全体の固定と試行変更の矛盾、実行環境の固定範囲、寸法以外の中間判定です。

確認対象は `9c52fdd41112f4f1aa1d6e712a257524665daf8b` と未 commit の §4.7.1。ファイル変更・forge 起動はしていません。以下は **plan 未反映**です。

採否表：

| 対象 | 採否・重大度 | 根拠と対案 |
|---|---|---|
| ユーザ欄と判定条件の範囲 | **要再検証・Major** | [plan:308](/home/sano/work/forge-cardmesh/plans/active/tooling-design-orchestrator.md:308) の `qualification` には評価点数と定常性の条件がありません。現行では `line_n` を run の基準から採ります（[wt_qualify.py:203](/home/sano/work/forge-cardmesh/design/forge_design/feedback/wt_qualify.py:203)）。また `_OVERRIDE_FORBIDDEN` は限定的で、`evaluate.settle`・`gas.transport` を拒否しません（[wt_card.py:115](/home/sano/work/forge-cardmesh/design/forge_design/intake/wt_card.py:115)）。**固定対象を YAML の所属節で決めず、生成後の実効条件で定義してください。** 評価位置・点数・量の定義・定常性の閾値と時間窓・物理モデル・物性を承認対象に含めます。`recipe`・`overrides` は、その固定条件を変えない範囲で選択可能にします。 |
| カード全体の sha256 一致 | **却下・Major：試行カードに対する必須条件として** | [plan:324](/home/sano/work/forge-cardmesh/plans/active/tooling-design-orchestrator.md:324) は変更のたびに再承認を要求します。一方、既存手順は較正係数を加えた `card_p1.yaml` を生成します（[quickstart:614](/home/sano/work/forge-cardmesh/procedures/design-chain-quickstart.md:614)）。このままでは許可した較正・格子変更でも止まります。**承認元カードの原本と sha256 は保存し、試行入力を別に持ってください。** 試行入力には承認記録への参照と自身の sha256 を記録し、解決後の保護対象が承認値と一致することを検査します。 |
| 別 checkout・依存環境の固定 | **採用・Major：現在の検査条件では不足** | [plan:329](/home/sano/work/forge-cardmesh/plans/active/tooling-design-orchestrator.md:329) の HEAD・追跡ファイル・Python・`pip freeze` だけでは、実際に使うコードと道具を限定できません。現行は `FORGE_MESH_PY` で別の Python を選びます（[deliverables.py:48](/home/sano/work/forge-cardmesh/design/forge_design/export/deliverables.py:48)）。STEP 読み取りも FreeCAD に依存します（[dimension_gate.py:50](/home/sano/work/forge-cardmesh/design/forge_design/export/dimension_gate.py:50)）。**判定専用環境を作り、import 元・外部 Python・FreeCAD/Gmsh・必要な変換器まで固定**してください。作業側の `PYTHONPATH`・実行ファイル選択用環境変数を引き継がず、コード探索先の未追跡・ignored ファイルも無視しない契約が必要です。 |
| 固定版での寸法測り直し | **採用・Major：寸法だけでは不足** | CAD には転送・回転面・独立読み直しの検査、mesh には再構築と品質の検査があります（[deliverables.py:246](/home/sano/work/forge-cardmesh/design/forge_design/export/deliverables.py:246)、[同:429](/home/sano/work/forge-cardmesh/design/forge_design/export/deliverables.py:429)）。現行の入口は総合 `verdict` と `failed_parts` を読み、必須部品について `ok=false` を見るだけです（[wt_qualify.py:773](/home/sano/work/forge-cardmesh/design/forge_design/feedback/wt_qualify.py:773)）。**合否に使う中間判定を列挙し、固定版での実行・入力との結び付きを確認できないものは再検査**してください。必須ファイル一覧も MANIFEST の自己申告から独立させます。停止・整定の判定も、この対象から外さないでください。 |
| `sources` の節・葉の指定 | **採用・Minor：正規化規則を補足** | [plan:312](/home/sano/work/forge-cardmesh/plans/active/tooling-design-orchestrator.md:312) の「最も細かく一致」は入力の簡略化として妥当です。承認時には実在する葉へ展開し、未知パスを拒否してください。`recipe_default` はラベルだけで認めず、固定した既定の提供元・キー・値と照合します。CSV から追加される表・ハッシュなどは導出情報なので（[wt_card.py:243](/home/sano/work/forge-cardmesh/design/forge_design/intake/wt_card.py:243)）、元のユーザ指定と結び付け、生成された葉を「省略された既定」と誤分類しないようにします。 |
| AI が承認記録を書く運用 | **採用・Minor：規約としての限定は妥当** | [plan:127](/home/sano/work/forge-cardmesh/plans/active/tooling-design-orchestrator.md:127)、[同:323](/home/sano/work/forge-cardmesh/plans/active/tooling-design-orchestrator.md:323) の強制と規約の区別は維持してください。ただし **sha256 と `user_ref` は承認の真正性を証明しません**。表示した表を固定する `preview` と、その同じ内容を書き出す操作を分け、表示後の再解決による値の変更を防ぎます。承認記録自体を書き換えられる構成では、保証は「指定された承認記録との整合性」に限ります。 |

`pip freeze` の sha256 は**補助記録としては採用、依存環境の固定の証拠として単独採用は却下**です。別 Python・外部実行ファイルを捉えない一方、無関係なパッケージの変更でも不一致になります。判定専用環境について、固定した配布物のハッシュを持つ依存一覧と、実際に起動する道具の識別を記録する方が目的に合います。editable install や作業ツリーへの参照を残すなら、その内容まで固定対象です。

評価位置は必須の固定対象です。反復・較正回数の上限は品質公差とは別の実行方針として開始時に解決して固定し、実施済み回数は台帳から取得してください。現行は `iteration`・`max_iterations` の欠落を許し、上限を判定しません（[wt_qualify.py:898](/home/sano/work/forge-cardmesh/design/forge_design/feedback/wt_qualify.py:898)）。

読み取り専用の関数単体確認では、次を確認しました。外部入力はメモリ上の fixture で代用しており、実 run の誤合格の再現ではありません。

- `_flatten_overrides` は `evaluate.settle`・`gas.transport` をエラーなしで受理。
- 同じカードの公差・η・評価窓でも、基準側の `n` により `line_n` が **2001 → 2** に変化。
- `FORGE_MESH_PY` の変更で、選ばれる外部 Python が変更。
- 来歴照合を `MATCH` と仮定した `deliverables_check` は、総合 `PASS`・`failed_parts=[]` の下で **`cad.pass=false` を渡しても `PASS`**。列挙から落ちた必須ファイルも独立には検出しません。

§5.1 #2 は、まず承認対象の解決と照合を作り、次に固定環境での判定・中間判定の再検査へ進める順を推奨します。既存の負例に、次を追加してください。

- 承認後の `line_n`・定常性の窓／閾値・物性変更を、**起動前に拒否**する。
- 原本を保った格子・較正係数だけの試行変更は、保護対象が同じなら**再承認なしで通る**。
- 同じ HEAD・追跡ファイルでも、外部 Python や import 元が違えば**判定前に止まる**。
- 総合 `PASS` と個別検査の FAIL／欠落が矛盾する、必須ファイルが列挙から消える場合は**accepted を出さない**。
- 出所の未知パス・既定値との不一致・正式入口で必要な予算の欠落を拒否する。`requirements_without_gate` の削除も承認条件の変更として拒否する。

正常例を併置し、「全部止める実装」を合格にしないことが必要です。

結論: §4.7.1 を、承認する実効条件と変更可能な試行入力を分離し、判定環境と全必須中間判定まで固定する契約へ修正する。

第 1 仮説: 最大の穴は、承認する実効条件が不足している一方、変更可能な試行カードまで丸ごと固定していることである。 確度: 高  
　根拠: `tooling-design-orchestrator.md:308・:324`、`wt_qualify.py:203`、`design-chain-quickstart.md:614`。  
　反証条件: 評価点数・定常性・物理条件の変更を起動前に拒否し、承認条件を保つ較正・格子変更は再承認なしで通る契約と検査が示されること。

第 2 仮説: 別 checkout が無変更でも、作業側の環境指定が判定処理へ入る。確度: 高。根拠: `deliverables.py:48・:229`。  
第 3 仮説: 寸法だけを再評価しても、他の必須検査の欠落・不合格を中間記録の総合 PASS が隠す。確度: 高。根拠: `wt_qualify.py:773` と上記の関数単体確認。

判別 A/B: `line_n=2001` を承認した同一契約に対し、生成された基準の `line_mach.n` だけを A＝2001、B＝2 とする。各1回、GPU 0 step、正式入口の起動前照合まで実行する。古い停止記録とのハッシュ不一致に頼らず、承認値との比較を確認する。**A が通り B が承認条件不一致で止まれば、この経路の固定漏れ仮説を棄却。B も通れば支持。A も止まれば試験不成立。**

やらない方がよいこと: CLI 公差の拒否だけで固定完了とすること、カードの sha256 不一致を回避するため承認記録を自動更新すること、寸法の再評価だけで成果物全体を信用すること。

呼び出し側の前提への異議: 「現行 `_OVERRIDE_FORBIDDEN` がユーザの意味に触る上書きを拒否する」「別 checkout と `pip freeze` の照合で作業側の変更が効かなくなる」は受け入れない。前者は検査範囲が不足し、後者は実行環境の選択経路が残る。

不足情報: 実装検証時には正式入口の起動方法と判定専用環境の構成が必要。今回の設計修正判断を妨げる不足ではない。

証拠: 足りる
