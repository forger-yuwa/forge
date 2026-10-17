# codex 諮問 (diagnose): orchestrator-direction

- **brief**: [`notes/reviews/briefs/2026-10-17-orchestrator-direction.md`](../../notes/reviews/briefs/2026-10-17-orchestrator-direction.md)
- **plan**: [`plans/active/tooling-design-orchestrator.md`](../../plans/active/tooling-design-orchestrator.md)
- **date**: 2026-10-17
- **commit**: `4039c301` (feature/card-local-mesh)
- **分類**: `ok`
- **backend**: `codex` — 選んだ理由: 既定 (auto: codex を試し、unavailable のときだけ claude)
- **要求したモデル**: `既定 (~/.codex/config.toml)`
- **実効のモデル**: `gpt-6-astra` (codex のログの見出し)
- **CLI の版**: `codex-cli 0.162.1`
- **依頼側との関係**: 不明 (`--requester-model` 未指定)
- **permission_denials**: 対象外 (codex は read-only サンドボックス)
- **実行**: effort `high`, 4.9 min, rc=0
- **試行のログ**: `notes/reviews/2026-10-17-orchestrator-direction-diagnose.codex.log` (codex: ok)
- **結論**: 次の一手は、§4・§6を「承認した要求と判定条件を固定し、予算・既存データは実行時に守り、途中の修正はその範囲内でAIに任せる」方針へ書き直すこと。
- **extra**: `procedures/design-chain-quickstart.md`, `design/forge_design/feedback/wt_qualify.py`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

改訂の方向は採用します。ただし、**A1〜A4をそのまま追加すれば「4つをプログラムで守れる」という判断は採用できません。** 承認した要求・判定条件との照合と、実行中の予算・書き込み制御が必要です。専用の判断役を省く方針には賛成です。

確認対象は `4039c301` と調査案 `85712d7d` の§3・§4です。依頼ファイルに埋め込まれたplan・手順書・`wt_qualify.py`は作業ツリーと一致しました。ファイル変更・forge起動はしていません。以下は改訂案であり、**plan未反映**です。

採否表：

| 対象 | 採否・重大度 | 根拠と対案 |
|---|---|---|
| A1 判定器の保護 | **採用・Major：固定commitの別checkout方式に修正** | 現行はCLIがカードの公差より優先される（[wt_qualify.py:179](/home/sano/work/forge-cardmesh/design/forge_design/feedback/wt_qualify.py:179)）。コードのshaだけでは防げない。動く`main`ではなく、依頼開始時に承認したcommit・依存環境・判定条件を固定し、作業側から変更できない場所で判定する。抽出位置・窓・必須成果物も固定対象。 |
| A2 予算の保護 | **採用・Major：終了後加算と次回起動拒否だけでは不足** | 現行はバイナリを同期実行するだけで、時間・ディスクの予算処理がない（[run_case.sh:55](/home/sano/work/forge-cardmesh/solver_density_cuda/tools/run_case.sh:55)）。起動前の予約、実行中の監視・停止、終了時の精算を同じ台帳に置く。失敗・計時run・前段も計上し、再開で消費量をリセットしない。 |
| A3 要求と物理設定の照合 | **採用・Major：最優先** | `qualify`は`M_exit`を照合する一方、`r_throat`不一致は警告で、入口条件・壁温等の要求との意味的な比較を行っていない（[wt_qualify.py:949](/home/sano/work/forge-cardmesh/design/forge_design/feedback/wt_qualify.py:949)）。承認した要求を、生成後の実効設定・境界条件・外部の物性／分布ファイルと照合する。 |
| A4 accepted成果物のchmod | **却下・Major：主たる保護として不足** | 保護対象には入力参照や失敗した試行も含まれる。現行の再実行は来歴・ログ・判定を上書きする（[run_case.sh:52](/home/sano/work/forge-cardmesh/solver_density_cuda/tools/run_case.sh:52)、[wt_qualify.py:1118](/home/sano/work/forge-cardmesh/design/forge_design/feedback/wt_qualify.py:1118)）。試行ごとの新規出力先と、確定済みデータへの書込権限の分離を使う。所有者が戻せるchmodは補助にとどめる。 |
| A5 判断役の口を省く | **採用・Major：§4.3も同時に修正** | §4.3には「判断役は番号を選ぶだけ」「一覧外は人に返す」「問いの組へ追加」が残っている（[orchestrator:129](/home/sano/work/forge-cardmesh/plans/active/tooling-design-orchestrator.md:129)）。そのままでは改訂と矛盾する。既知の復旧手順は残すが、網羅的な選択肢一覧と専用の判断役実装は延期する。 |
| A6 依頼の課題集 | **採用・Major：6件は初期試走として扱う** | GPUなしでは、流れの公差達成や成果物までの到達を判定できない。現planも入力検証V1aとGPU検証V1bを分離している（[orchestrator:240](/home/sano/work/forge-cardmesh/plans/active/tooling-design-orchestrator.md:240)）。課題集に、下記の安価な負例検査を組み合わせる。 |

**Major：出口で固定すべきなのは、判定器のファイルだけでなく「何を合格とするか」です。**

読み取り専用で現行の`resolve_tolerances`を呼び、同じ入力について次を確認しました。

```text
カードの qualification.exit_core_M_tol_rel = 0.001
CLI指定なし                         → value: 0.001, source: カードの qualification
CLIで exit_core_M_tol_rel = 0.1     → value: 0.1,   source: CLI
```

コード変更なしで許容値が100倍になります。これは関数単体の確認であり、実runの`accepted`を再現したものではありません。また、生成器が受け付けるカードには`qualification`をまだ書けないと明記されています（[wt_qualify.py:54](/home/sano/work/forge-cardmesh/design/forge_design/feedback/wt_qualify.py:54)）。「承認値をカードに書けば守れる」は未実装です。

最小の対案は、**承認済みカードと解決済みの判定条件を、既存のcampaign／`execution_plan.json`に結び付けること**です。別の汎用台帳は不要です。固定する内容は、要求する物理・寸法、評価量の定義と対象領域、公差、定常性の条件、必須成果物、判定器の版です。CLIや問題YAMLによる相違は、出所を記録して採用するのではなく拒否します。

固定版の`wt_qualify`を呼ぶだけでも十分ではありません。例えば寸法判定は既存の`DIMENSIONS.json`の`verdict`を読みます（[wt_qualify.py:760](/home/sano/work/forge-cardmesh/design/forge_design/feedback/wt_qualify.py:760)）。その結果を作った判定処理も固定版に結び付け、確認できなければ再評価する必要があります。ハッシュは内容の同一性を示しても、承認や判定の正しさを単独では示しません。

**Major：要求の一致、解析の成立、要求の達成を区別してください。**

現行の成果物照合には壁・格子・設定・場との結び付けが既にあります（[wt_qualify.py:678](/home/sano/work/forge-cardmesh/design/forge_design/feedback/wt_qualify.py:678)）。これを捨てて作り直す必要はありません。足りないのは、それらが**承認した問題を解いているか**の照合です。

出口で必須にするものと、人の報告確認に残せるものの境界は次です。

| 扱い | 内容 |
|---|---|
| 機械的な照合が必要 | 入口条件・組成とその基準、壁熱条件、採用した物理モデル、軸対称／平面の別、対象形状・解析領域、評価位置、公差、必須成果物。要求された判定が欠ければ`incomplete`。 |
| 起動前・実行中に必要 | 承認の欠落、予算、既存データへの書き込み、禁止設定、メッシュ品質。出口で発見しても消費や上書きを取り消せない。 |
| 人の報告確認でよい | 試行理由、原因についての仮説、設計上の考察、要求範囲外への適用可能性。ただし、これらを必須検査の代わりにしてはならない。 |

A3の照合は単純なYAML文字列比較では足りません。組成の換算、物理スロート径と設計側の`r_throat`、壁温分布などの導出関係を扱う必要があります。親planの実効設定と全段照合の契約を再利用するのが適切です（[campaign-recipe:89](/home/sano/work/forge-cardmesh/plans/active/tooling-design-problem-campaign-recipe.md:89)）。

また、**数値設定を調整してよいことと、数値的な妥当性を検査しなくてよいことは別です。** 初期化用の1次精度と本段の2次精度は区別されています（[solver-settings.md:9](/home/sano/work/forge-cardmesh/procedures/solver-settings.md:9)）。既存の`launch_gate`と最終評価の要件を維持し、変更した離散化・格子に必要な検証を通してください。全CFL・リミッタ値を固定する必要はありません。

**Major：新しい機種で判定器がない場合、承認済みの基準文だけでは自動合格にできません。**

最小の仕組みは、必須の要求ごとに「評価方法・公差・証拠・判定担当」を事前に決め、**未評価が1つでもあれば`accepted`を出さない**ことです。これは親planの「要求された評価の未実行・判定ファイルの欠落は不合格」と整合します（[campaign-recipe:104](/home/sano/work/forge-cardmesh/plans/active/tooling-design-problem-campaign-recipe.md:104)）。

判定器がない探索は、予算内で進めて結果を出してよい。ただし結果は`incomplete`として未評価項目を列挙します。人が採用する場合も、人の判断を別途記録し、数値VERDICTを書き換えません。新しい判定器をAIが実装した場合は、既存の開発レビューを通して固定してから使います。専用の昇格審査は不要です。

**Major：エスカレーション条件4の免除は、作業の範囲で書き分けてください。**

「設計依頼の途中なら適用しない」という一括免除は広すぎます。次の限定を推奨します。

> 承認済みの要求・判定条件・予算を保ち、既存データを変更せず、新しい試行で行う設定調整・格子再生成・case内の探索台本の変更は、条件4の事前諮問を不要とする。正式なソルバ・判定器・共通ツールの変更、承認条件の変更、原因や結果解釈の確定には従来の条件を適用する。

条件2の「2回対処しても解けない」、条件3の事前登録した検証のFAIL、条件6の数値カーネル変更、条件7の解釈確定は残します（[AGENTS.md:339](/home/sano/work/forge-cardmesh/AGENTS.md:339)）。通常の設計候補の公差未達と、開発検証のFAILも区別してください。

同期が必要な箇所は以下です。

- `design-chain-quickstart.md:14`の手順厳守・修正禁止を、決められた流れを回す役に限定する。権限はモデルの自己申告でなく、起動時の役割で決める。
- 同手順書の「寸法公差は聞いてなければこのまま」を、承認済みの既定値に限定する（[quickstart:277](/home/sano/work/forge-cardmesh/procedures/design-chain-quickstart.md:277)）。
- `design-intake`の「許容値変更はコメントに理由」を、承認した版との結び付けに置き換える（[SKILL.md:138](/home/sano/work/forge-cardmesh/.claude/skills/design-intake/SKILL.md:138)）。
- 同skillの未対応メニュー分岐では、既存ツールの新しい組合せによる探索と、未実装の物理機能の開発を分ける（[SKILL.md:128](/home/sano/work/forge-cardmesh/.claude/skills/design-intake/SKILL.md:128)）。
- orchestratorの§1・§2・§4.3・§4.8・§5・§7も更新する。§4.4だけの差し替えでは旧方針が残る。

**Major：`run_case.sh`を「唯一の入口」と見なす根拠は不足しています。**

フック自身がシェル変数経由の呼び出しを検知できないと明記しています（[hook_forge_guard.py:12](/home/sano/work/forge-cardmesh/solver_density_cuda/tools/hook_forge_guard.py:12)）。Claude以外にも効く保護を主張するなら、AIから見える実行権限と書込権限で担保する必要があります。

最初は、既知の流れと探索台本の両方が同じ小さな実行入口を使い、そこで承認予算・attempt ID・新規出力先を必須にしてください。GPUを直接起動できる権限や保護先を書き換える権限がAIに残る構成なら、その構成は規約による運用であって、プログラムによる強制ではありません。

§4.4の停止記録は、試行ID、段、判定器・VERDICT、対象量・値・閾値・区間、証拠パスと必要なハッシュ、消費予算、プロセスの状態で十分です。原因不明を4分類へ無理に押し込めないでください。SSH切断だけで一過性の再試行を選ばず、§4.5のプロセス照合を先に行います。

V3は、正常3件・停止3件を残し、次を事前登録してください。

- **GPUなしの到達点**はカード→YAML→格子まで。これを設計全体の`accepted`に数えない。流れの公差未達は、固定した既存の不合格データの再判定か、後続のGPU試験で扱う。
- **安価な負例**として、根拠のないユーザ欄、CLIでの公差・評価窓・必須成果物の変更、物理設定の不一致、判定器なし、予算の実行中超過、再開による二重起動、既存出力先への書込みを入れる。
- **誤合格**は機械結果とAI報告で別々に数える。JSONが`incomplete`でもAIが「要求達成」と報告したら報告側の誤合格。分母と途中停止も記録する。
- モデル版・権限・入力・予算・人の介入条件を固定する。初回は各依頼を各モデルで3回程度繰り返し、到達率・正しい返却率・誤合格件数・費用を出す。観測された誤合格0を一般的な保証とは扱わない。

§5.1の優先順は次を推奨します。

1. **承認要求と判定条件の固定**。現#9のカード出所を先頭へ移し、A1・A3と判定器なしの扱いを定義する。
2. **親planの最小契約と出口の統合**。入力参照・実効設定・`execution_plan`・`judge_run`・qualification。
3. **予算・新規出力先・台帳・再開を一体で実装**。停止理由の記録もここへ含める。
4. **既知の線形フローと格子の包み、最小CLI**を作り、V1a・V2・GPUなしV3を通す。
5. **AWSのidle停止を解決してからGPUの通し検証**を行う。

専用判断役、`local`バックエンド、引用照合、判断役用の問いの組、`forge-design answer`は現スコープから外してよいです。手の一覧は実際に使う既知の復旧手順だけ残します。格子係数の一般化等の承認済み課題は削除せず別マイルストーンに保ち、直接のNS補正経路も今回の整理を理由に削除しません（[orchestrator:233](/home/sano/work/forge-cardmesh/plans/active/tooling-design-orchestrator.md:233)）。

結論: 次の一手は、§4・§6を「承認した要求と判定条件を固定し、予算・既存データは実行時に守り、途中の修正はその範囲内でAIに任せる」方針へ書き直すこと。

第 1 仮説: 最大の穴は、判定器の版と来歴が揃っても、実際に使う判定条件・解析条件が承認要求に固定されていないことである。 確度: 高  
　根拠: `wt_qualify.py:179`のCLI優先と関数単体の実測`0.001 → 0.1`、`:949`以降の要求照合範囲、`:749`の`optional_parts`による必須成果物の変更。  
　反証条件: 正式な入口が承認済みの条件を固定し、公差・評価対象・必須成果物・物理条件の未承認変更を判定前にすべて拒否することが確認される。

第 2・第 3 仮説:  
　第 2 仮説: A2・A4の現案では、終了前の予算超過と既存データの変更を防げない。確度: 高。根拠は`run_case.sh:52`・`:55`と、A2の次回起動時だけの検査、A4のaccepted限定。  
　第 3 仮説: §4.3の全手列挙を維持すると、専用判断役を削っても改訂の自由度と簡素化が実現しない。確度: 中。根拠は`tooling-design-orchestrator.md:129`。実際の運用負荷は未測定。

判別 A/B: 承認済み許容値`exit_core_M_tol_rel=0.001`を固定した同一入力で、Aは追加指定なし、BはCLIの同項目だけ`0.1`にする。各1回、GPU 0 step、正式入口の条件解決までを実行し、拒否の有無と解決後の値を見る。Bが`0.1`を採用すれば第1仮説を支持し、承認条件との不一致として拒否すれば当該経路について棄却する。現行の関数単体では前者を確認済み。

やらない方がよいこと: 判定器のshaを記録しただけで保護完了にすること、§4.3を据え置いて§4.4だけ削ること、GPUなし6件の試走で設計全体の自律性を合格にすること。

呼び出し側の前提への異議: 「出口だけで4つを守れる」「`run_case.sh`は全エージェントに対する強制された唯一の入口」「数値設定は照合不要」「基準を文章で承認すれば判定器なしでもAIが合格を出せる」は受け入れない。上記のコード・契約では成立しない。

不足情報: 実装時には、forge-gate統合後の正式入口と権限構成、依頼予算の計数単位、V3の固定入力が必要。今回の設計判断を妨げる不足ではない。

証拠: 足りる
