# codex 諮問 (diagnose): orchestrator-executor

- **brief**: [`notes/reviews/briefs/2026-10-17-orchestrator-executor.md`](../../notes/reviews/briefs/2026-10-17-orchestrator-executor.md)
- **plan**: [`plans/active/tooling-design-orchestrator.md`](../../plans/active/tooling-design-orchestrator.md)
- **date**: 2026-10-17
- **commit**: `7cc17b9f` (feature/card-local-mesh)
- **分類**: `ok`
- **backend**: `codex` — 選んだ理由: 既定 (auto: codex を試し、unavailable のときだけ claude)
- **要求したモデル**: `既定 (~/.codex/config.toml)`
- **実効のモデル**: `gpt-6-astra` (codex のログの見出し)
- **CLI の版**: `codex-cli 0.162.1`
- **依頼側との関係**: 不明 (`--requester-model` 未指定)
- **permission_denials**: 対象外 (codex は read-only サンドボックス)
- **実行**: effort `high`, 4.2 min, rc=0
- **試行のログ**: `notes/reviews/2026-10-17-orchestrator-executor-diagnose.codex.log` (codex: ok)
- **結論**: §4.5.1 を、上記の recipe 分岐・復旧・予算・入力固定の契約を含む形に修正する。

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

2 層の service は採用します。ただし、§4.5.1 は修正してから実装すべきです。特に、現在の recipe の分岐では正常な run も成果物生成へ進めません。

ファイル変更・forge 起動は行っていません。以下は作業ツリーのコード・手順書・既存の試験記録に基づく設計レビューです。**Critical 0 / Major 5 / Minor 1。plan 未反映です。**

採否表:

| 重大度・対象 | 採否 | 根拠・理由・対案 |
|---|---|---|
| **Major 1：recipe の分岐** | **却下して修正** | [plan:328](/home/sano/work/forge-cardmesh/plans/active/tooling-design-orchestrator.md:328) は `accepted → deliverables`、その他を停止としている。しかし現行手順は **`incomplete / make_deliverables → deliverables`** で、`accepted` は成果物がある場合だけ出る（[手順書:667](/home/sano/work/forge-cardmesh/procedures/design-chain-quickstart.md:667)）。`gate_env` も判定器の終了コードをそのまま返す（[gate_env.py:326](/home/sano/work/forge-cardmesh/design/forge_design/feedback/gate_env.py:326)）。また、初回較正は p0 の NS の後に行い、再較正の条件判定とは別である（[手順書:579](/home/sano/work/forge-cardmesh/procedures/design-chain-quickstart.md:579)、[同:653](/home/sano/work/forge-cardmesh/procedures/design-chain-quickstart.md:653)）。**p0 は関門通過後に初回較正、p1 以降は既存の「判定＋次の手」で分岐**させる。初回較正を再較正回数に混ぜない。 |
| **Major 2：driver の復旧・二重起動・識別** | **2 層構成は採用、復旧契約は要修正** | [plan:309](/home/sano/work/forge-cardmesh/plans/active/tooling-design-orchestrator.md:309) の単一書き手と原子的更新だけでは、service 起動と台帳更新の間の障害を扱えない。`is-active` と完了証拠の照合にも、起動途中・停止途中・照会失敗・unit 消失の扱いが必要。**試行 ID・unit 名・出力先・予約を先に永続化 → service 起動 → 結果を試行 ID ごとに一度だけ取り込み・精算**とする。driver は生存期間全体でロックを保持し、競合する `--resume` は既存 driver を止めず状態を返す。unit 名の `<name>` は表示名に限定し、永続的な campaign ID と試行 ID を使う。現状の campaign 名は空でない文字列という検査だけで、case をまたぐ一意性はない（[campaign.py:53](/home/sano/work/forge-cardmesh/design/forge_design/campaign/campaign.py:53)）。run 番号の割当ても README 更新と同じ case ロックの対象にする。 |
| **Major 3：予算停止と停止記録の順序** | **現在の記述は却下** | [plan:340](/home/sano/work/forge-cardmesh/plans/active/tooling-design-orchestrator.md:340) は「service を止めて停止記録を書く」という順序。自分も停止対象なので、後続の記録を保証できない。既存試験も、実行担当が死ねば子孫を止めることを確認したもの（[試験記録:14](/home/sano/work/forge-cardmesh/notes/investigations/2026-10-17-attempt-service-ab/README.md:14)）。**予算到達・計測値・停止要求を先に原子的に記録し、その後に停止**する。停止要求と「全子孫の終了確認」は分け、後者を確認してから driver が精算する。担当への `SIGKILL` で記録が欠けた場合は、復旧した driver が根拠付きで補完し、原因は確認できなければ「不明」とする。 |
| **Major 4：予算の測定範囲・超過量** | **単位と残額予約は採用、監視仕様は要修正** | [plan:338](/home/sano/work/forge-cardmesh/plans/active/tooling-design-orchestrator.md:338) の run サイズだけでは、run 外の格子・カード・ログ等が対象外になる。同じ NS run に後段が成果物を追加するため、段ごとのサイズの単純合算も二重計上になる。**依頼が新たに占有した出力を、重複しない所有範囲で集計**する。GPU 時間は前段を含む各 forge 起動の経過時間を永続記録し、担当死亡時の未確定区間をゼロ扱いしない。数分おきのディスク監視は超過検出であり厳密な上限ではない。最大出力サイズ・監視間隔・停止遅延を踏まえた余裕と、停止記録用の空き容量を確保する。AWS のインスタンス稼働時間とは別の指標であることも明記する。 |
| **Major 5：`fixed` の照合** | **再生成して全体の sha256 を比較する案は却下** | [plan:345](/home/sano/work/forge-cardmesh/plans/active/tooling-design-orchestrator.md:345) は生成器の実体の固定を規定していない。同じ変更済み生成器で再生成すれば、照合が一致しても「検証済み recipe の入力」とは言えない。現在の recipe も作業側の runner を import して構成される（[plan.py:51](/home/sano/work/forge-cardmesh/design/forge_design/campaign/plan.py:51)）。**生成器の版・依存する台本と外部入力を固定し、一度生成した入力のハッシュと生成元を記録して、使用直前に照合する**形を推奨する。再生成比較は V1a の検証に残す。`role` は campaign に固定し、`--resume` で黙って変更させない。これは #14 完了前には運用上の検査であり、権限による強制ではない。 |
| **Minor 1：ログと台帳の書き手** | **補足して採用** | [plan:313](/home/sano/work/forge-cardmesh/plans/active/tooling-design-orchestrator.md:313) は状態と停止記録しか置き場所を定めていない。driver の起動ごと・試行ごとに独立したログを残し、既存ログを上書きしない。journal だけを復旧の根拠にしない。また `forge-design note` が直接台帳を書くと「driver が唯一の書き手」と衝突する（[plan:348](/home/sano/work/forge-cardmesh/plans/active/tooling-design-orchestrator.md:348)）。CLI は変更要求を渡し、driver が入力のハッシュ・理由とともに取り込む形にする。 |

問い 2 の「見積もりの 1.5 倍で人に聞く」は、**一律の必須条件にはしません**。依頼ごとの承認という決定（[plan:89](/home/sano/work/forge-cardmesh/plans/active/tooling-design-orchestrator.md:89)）に対して、固定倍率を新たな承認境界にする根拠がありません。単一試行だけを動かす今回の範囲では、残額全部を上限として予約して構いません。ただし、起動前には段の見積もりが残額内か確認し、見積もりとの差は記録・警告に使います。承認総額の超過が見込まれる場合は次段を始めず返します。

問い 4 の順序は、**#6a → #6b → #6c → #6d を維持し、最小 driver と復旧処理だけ #6b に前倒し**します。試行 service 単体では driver 障害を検証できません。合格条件には、起動要求と台帳記録の間の障害、driver 不在中の試行完了、結果取り込み・精算途中の障害、同名の別 campaign、停止途中の再開、ログ保存先の容量不足を追加してください。既存試験は予算・台帳を検証していないと明記されています（[試験記録:42](/home/sano/work/forge-cardmesh/notes/investigations/2026-10-17-attempt-service-ab/README.md:42)）。

問い 5 の手元 GPU による V4(a) 前倒しは、**却下**します。#7 後に GPU の通しを行う計画を維持します（[plan:539](/home/sano/work/forge-cardmesh/plans/active/tooling-design-orchestrator.md:539)）。ブリーフにある手元障害の原因は今回確認していませんが、AWS 優先を覆す材料はありません。GPU なしの #6 を先に完了できます。#6 の先行実装を、未完了の #5・#8 の関門を省いて通しを合格にする理由にはしないでください。

結論: §4.5.1 を、上記の recipe 分岐・復旧・予算・入力固定の契約を含む形に修正する。

第 1 仮説: 現案をそのまま実装すると、流れの要求を満たした正常な run でも成果物生成前に停止する。確度: 高  
根拠: `plans/active/tooling-design-orchestrator.md:328` と `procedures/design-chain-quickstart.md:667` の分岐が不一致。  
反証条件: 成果物のない正常な run が、承認条件を変えず、現案の分岐のまま成果物生成と最終 `accepted` まで到達すること。

第 2 仮説: driver 障害をまたぐ起動・結果取り込み・精算は、現在の記述だけでは一度だけの処理を保証できない。確度: 高（契約の不足。実装上の発生は未確認）。

第 3 仮説: 再生成ハッシュの一致だけでは、`fixed` が変更済みの生成器や台本を使うことを排除できない。確度: 高（生成器を別途固定する実装は未確認）。

判別 A/B: 既存の正式な入口の正例で使った run の写しを用い、変えるのは必須成果物の有無だけとする。A は成果物あり、B は成果物なし。同じ承認・場・判定条件で `gate_env judge` を各 1 回実行し、判定と recipe の遷移先を見る。GPU 不要。A が `accepted`、B が `incomplete / make_deliverables` なのに停止すれば第 1 仮説を支持する。B が成果物生成へ進めば第 1 仮説は反証される。A も失敗すれば試験不成立。

やらない方がよいこと: 成果物前の停止を避けるために必須成果物を任意化すること、復旧時に同じ runner を再実行すること、停止記録を自分の service 停止後に書くこと。

呼び出し側の前提への異議: 「今の手順書の流れ」とする recipe 表は現行手順と一致していない。また、既存 service 試験の PASS は子孫停止の証拠であり、2 層構成の復旧・予算精算の証拠ではない。

不足情報: 実装時には、予算計測の単位・対象パス・監視間隔・停止遅延・障害時の未確定消費量の扱いを確定する必要がある。今回の設計判断は妨げない。

証拠: 足りる
