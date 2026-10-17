# 設計の実行役: campaign 実行器と入口 (`tooling-design-problem-campaign-recipe` §5 ステップ 5 の実装)

## メタ

- **area**: `その他` (設計チェーンの道具)
- **status**: `draft`
- **related_docs**:
  - `procedures/nozzle-design-workflow.md` (今の設計の流れ)
  - `procedures/nozzle-design-outputs.md` (成果物)
  - `procedures/solver-settings.md` (本段の監視用の出力と協調停止)
- **related_plans**:
  - **親**: [`tooling-design-problem-campaign-recipe.md`](tooling-design-problem-campaign-recipe.md) (problem / campaign / recipe、`resolve_evaluation` → `execution_plan.json`、3 値判定 `accepted / rejected / incomplete`、campaign 実行器 + ledger。この plan はその §5 ステップ 5 と、要求カードからの入口の実装。語彙はそちらに合わせる)
  - `plans/accepted/tooling-settle-convergence-and-monitor.md`: 判定・停止規則 v2・受け渡しの関門 #7。この plan の前提。
  - `plans/accepted/tooling-nozzle-card-local-mesh.md`: case/71 の要求カード・手組みの格子・3 係数の壁。
  - `tooling-generator-defaults-and-launch-gate` (作業ツリー `forge-gate`): 推奨設定の機械可読の正本と、起動前の関門。統合する。
  - 推力ノズル (`tooling-nozzle-bell-two-mode`・`tooling-nozzle-dual-bell`)、SERN (`tooling-nozzle-sern-chain`): 後の機種。それぞれのセッションの担当。
- **created**: `2026-10-16`
- **owner**: M6ノズルの設計 (主セッション)

## 1. 目的

ユーザ (2026-10-16):

- 「何も知らない、あほなAIでも設計できる状況になってると思う？」
- 「このやり方を発展させて、いろんな機種、いろんな設計要求に対応させられるようにできるかな？」
- 「実行役が止まったときにだけ、上位の判断役 これいいやん。ただ必ずしもfableが使えるとも限らず、ローカルLLMかもしれん」
- 「そのやり方でよいか、diagnosticianにも相談した上で、進めてほしい」

ユーザ (2026-10-17):

- 「ツール群（forgeや、CEAや、熱評価ツール、gmshなど）がたくさん備わった中から、ユーザのニーズに基づいて設計プロセスを組んでくれるような状況を最終的に想定している」
- 「AIにある程度じゆうにうごいてもらうかたちのほうが良いのでは」
- 調査のセッションへの指摘 (調査ファイル §3): 「ややこしさが増しそうで心配」「ルールが増えすぎるのも心配」「ある程度 AI の実力に期待する形で」

この plan の目的は次の 3 つ。

- 既知の機種 (最初は風洞ノズル) は、要求カードから成果物までを 1 つの入口 (CLI) で通す。段を回すのは決定的なプログラム (実行役) で、止まったら理由を 1 つのファイルに書いて終わる。
- 賢い AI には途中 (数値設定・台本・格子の作り直し・手順の並べ替え) を直させ、新しい流れを台本で組ませる。守りは決まりの文章でなく、
  (a) 承認した要求と判定条件の固定、(b) 予算と既存データの実行時の保護、(c) 出口の判定、の 3 つをプログラムに置く (§4.0)。
- 知識の少ないモデル (ローカル LLM など) には、検証済みの決まった流れを回させるだけにする (直させない)。

## 2. スコープ

- **やる (最初の範囲)**: 風洞ノズル (axismach) で、search の無い campaign を要求カードから回す。
  - 承認した要求と判定条件の固定 (カードの出所・判定条件の節・判定器の版・判定器の無い要求の扱い)。
  - 要求の一致の照合 (承認した問題を解いているか)。
  - 親 plan の未実装の契約: 最小の campaign、入力の参照の解決、pass ごとの `execution_plan.json`、全段の設定の照合。
  - 実行の入口・予算・試行ごとの新しい出力先・台帳・再開・停止の記録。
  - 格子の段: 包むだけ。
  - 道具の `--help`・skill・手順書・AGENTS.md の同期。
- **やる**: case/71 (手組みの格子) と case/70 (雛形の格子) で検証する。答えの分かっている依頼の課題集 (§6 V3) を、賢い AI と知識の少ないモデルで回して測る。
- **やらない**:
  - 汎用 DAG・新しい状態機械の枠組み・新しい VERDICT の語彙 (親 plan §2)。
  - 判断役の専用の仕組みなど、§4.9 の「今は作らないもの」。
  - ベル・デュアルベル・SERN のレシピ (各セッションの成果を後から載せる)。
  - 機種ごとの物理の公差を決めること (ユーザと決める)。
  - ソルバの演算の変更。
- **別のマイルストーン** (この plan の §5.1 に置くが、完了条件 §8 には入れない): 格子の係数の一般化 (#10)、別の寸法での収束の判定の確認 (#11)、
  生の排除厚さの経路を生産から外す (#12)。

## 3. 関連 docs と前提

- 前提: settle plan の #7 (受け渡しの関門の新方式) と V6 (設計チェーンの通しの 1 pass) が通っていること (どちらも accepted)。この plan は、その流れを部品にして束ねる。
- 改訂案の出典: `notes/investigations/2026-10-17-ai-design-survey.md` §3・§4 (セッション `AI設計の調査`、commit 85712d7d)。
- 今ある部品:

| 部品 | 役割 |
| --- | --- |
| `design/forge_design/intake/wt_card.py` | 要求カード → 問題 YAML |
| `design/forge_design/intake/recipes/axismach_*.yaml` | 雛形 |
| `case/71.m55_combustion_wt/_mesh/hand_mesh.py` | 手組みの格子。case/71 専用 |
| `design/forge_design/evaluate/runner_axismach.py` | Euler・NS の準備と起動 |
| `design/forge_design/evaluate/stop_rule.py` | 停止規則 v2 |
| `design/forge_design/feedback/deltastar_loop.py` | 境界層の補正の繰り返しと受け渡しの関門 |
| `design/forge_design/feedback/wt_qualify.py` | 出口の判定 (qualification)。公差の出所は今は CLI が最優先 (`resolve_tolerances`、179 行) |
| `design/forge_design/export/deliverables.py` | 成果物 |
| `design/forge_design/report/nozzle_report.py` | 設計報告 |
| `solver_density_cuda/tools/run_case.sh` | forge の起動 (同期の実行だけ。予算は扱わない) |
| `solver_density_cuda/tools/hook_forge_guard.py` | forge の直接実行を拒む Claude Code のフック (シェル変数経由の呼び出しは検知できない、12 行) |
| `solver_density_cuda/tools/codex_review.py` | 外部レビューと諮問 (この plan の実行役からは呼ばない) |
| `solver_density_cuda/tools/aws_instance.sh` | AWS の起動と状態 |

## 4. 設計方針

**ユーザ決定 (2026-10-17)**:

- 方向: 設計方針の改訂案 (`notes/investigations/2026-10-17-ai-design-survey.md` §3) で進める (ユーザ「おｋ」)。ローカル LLM について「出口の判定で守る、というのがちゃんとわかんない」と
  質問があり、「知識の少ないモデルは決まった流れを回すだけで直さない、賢い AI が直す」と説明して了承された。LLM が流れ (台本) を組む形は残る。
- **予算**: 一律の上限は置かない (「ケースバイケース」: パレート解を取るとき・メニューしだいで大きく変わる)。依頼ごとに AI がメニューと本数 × その PC で測った
  1 本の時間から GPU 時間・ディスクを見積もり (`design-intake` §8「コスト見積もりと承認」)、ユーザが承認した値をその依頼の予算とする (§4.5・§4.7)。
- 依頼の原文の置き場所: 設計を行うリポジトリの case のディレクトリ (§4.7)。

経緯: 2026-10-16 に Fable の `diagnostician` の諮問と codex の plan 段のレビューを全件採用した前の案 (判断役の口・問いの組など) を、
2026-10-17 のユーザ決定と codex (diagnose) の諮問 (全件採用、§6.1) で書き直した。前の案の本文は commit 4039c301 にある。

### 4.0 原則・役割・AI が変えてはいけない 4 つ

**原則** (調査ファイル §3.1 を codex の諮問で修正):

1. 守りは「承認した要求と判定条件の固定」「予算と既存データの実行時の保護」「出口の判定」の 3 つに置き、途中の修正はその範囲の中で AI に任せる。
   出口の判定だけでは守れない (判定器が測っていない量・出口で見つけても取り消せない消費や上書きがある; codex 2026-10-17)。
2. 仕組みは、同じことを繰り返すようになってから作る。実績が少ないうちは人に聞く。
3. 決まりを足すときは、文章ではなくプログラムの検査にする。検査を足したら、置き換わる文章の決まりを消す。

**役割**:

| 担い手 | 持つもの |
| --- | --- |
| ユーザ | 依頼 (原文)。要求の解釈の確認 (カードのうちユーザが決める欄、判定条件)。予算の承認。新しい手順を正式な道具にするかの判断。判定器の無い要求の結果を採用するかの判断 |
| 直す役の AI (主セッションの賢いモデル) | カードを作る。道具を選んで回す。失敗したら原因を調べて直す・別の手を試す (台本を書いてよい; §4.4 の範囲)。結果を判定の出力つきで報告する |
| 決まった流れを回す役 (知識の少ないモデル) | 実行役を起動し、止まったら停止の記録を人に見せる。直さない |
| プログラム | 承認した条件の固定と照合。合否の判定。予算の予約・監視・精算と停止。試行ごとの新しい出力先。記録を自動で残す |

- 役割は**起動のときに決める** (`forge-design run --role fixed|fixer`、台帳に記録)。モデルの自己申告では決めない。
  `fixed` では、実行役はレシピが生成した入力と違う入力 (設定・格子・台本) を拒む。

**AI が変えてはいけない 4 つと、守る仕組み**:

| 変えてはいけないもの | 守る仕組み | 節 |
| --- | --- | --- |
| 1. 判定器と判定条件 (公差・評価位置と窓・定常性の条件・必須成果物・判定器の版) | 承認の時点で固定し、作業側から変更できない場所 (固定した commit の別の checkout と依存環境) で判定する。CLI・問題 YAML による相違は出所を記録して採用するのでなく**拒否**する | §4.2・§4.7 |
| 2. ユーザが決めた欄 (入口条件・組成・壁温・寸法・成果物の選択など) | カードの出所と承認の版 (§4.7)、実効の設定・境界条件・外部の物性と分布のファイルとの照合 (§4.2 の要求の一致) | §4.2・§4.7 |
| 3. 予算 (GPU 時間・ディスク) | 実行の入口が起動前に予約し、実行中に監視して止め、終了時に精算する | §4.5 |
| 4. 既存の run と成果物 (入力の参照・失敗した試行を含む) | 試行ごとの新しい出力先 (既にあれば止める)、確定済みのデータへの書き込みの権限の分離。所有者が戻せる chmod は補助にとどめる | §4.5 |

- 4 つのどれかを変える必要があると AI が判断したら、ユーザに聞く。
- **強制と規約の区別** (codex 2026-10-17): AI が GPU を直接起動できる、または保護先を書き換えられる権限の構成のままなら、上の仕組みは**規約による運用と、出口での検出**であって、
  プログラムによる強制ではない。最初の範囲 (手元と今の AWS) はこの状態であることを報告と記録に書く。強制にするには、実行の権限と書き込みの権限を AI から分ける構成が要る (§5.1 #14)。
- `hook_forge_guard.py` は Claude Code にしか効かず、シェル変数経由の呼び出しも検知できない。`run_case.sh` を「すべてのエージェントに強制された唯一の入口」とはみなさない。

### 4.1 位置づけと語彙

- この plan は、親 plan `tooling-design-problem-campaign-recipe` の §5 ステップ 5 (campaign 実行器・ledger) と、要求カードからの入口の実装である。
- **汎用 DAG は作らない**。段は線形の列で、δ* の繰り返しは recipe の内部に置く。
- 語彙は親 plan に合わせる:
  - レシピ名は `axismach.*/vN`
  - 段と campaign の判定は `accepted / rejected / incomplete`
  - 来歴は `execution_plan.json` (+ `StageManifest`・`MAIN_STAGE_STARTED`・`HANDOFF_GATE.json`)
- **前提 (2026-10-16 codex M1)**: 親 plan の未実装の契約を、実行器より先に作る。
  - search の無い最小の campaign。
  - 入力の参照の解決 (snapshot・格子・物性・バイナリの識別の固定)。
  - pass ごとの `execution_plan.json`。
  - 全段の設定の照合。
  - 例: `prepare_ns` は `ic_from` のディレクトリの最大番号の場を選ぶ (`runner_axismach.py:2511`)。計画を確定した後に参照の run の出力が増えても、入力が変わらないようにする (負例の試験)。
- 承認した要求と判定条件 (§4.7) は、別の台帳を作らず、この campaign と `execution_plan.json` に結び付ける。

#### 4.1.1 親 plan の契約のうち今回作るもの (§5.1 #4 の設計、2026-10-17 codex 諮問で修正)

最初の案 (準備の前に入力の一覧を書く・走った後に段の列を照合する) は、codex (diagnose) の諮問
([記録](../../notes/reviews/2026-10-17-orchestrator-minimal-campaign-diagnose.md)) で次のように直した。
理由: 今の `run_staged_ns` は起動の時に `solverConfig.yaml` を読み直し、本段の `nStepOuter` を停止規則の上限に書き換える (`runner_axismach.py:2946`・`:3002`)。
準備の前の入力の一覧だけでは、実際に起動した設定と計画が結び付かない。段の記録の `key` (18 キー) は軸対称・輸送などを含まず、全段の照合の証拠にならない。

**計画の 3 段 (1 つの run に 1 つの `execution_plan.json`、段の記録は上書きしない)**:

1. **準備の前 (入力の固定)**: カード・問題 YAML・承認の記録の sha256、格子の msh の sha256、初期場の**ファイル** (パス・step・sha256)、Euler の参照の場、δ の表、
   forge の実行ファイルと変換器の sha256、数値に効く環境変数。**初期場の元の run から継承する値** (`limiter_ref_inherit_from_ic` のリミッタの基準値、
   `runner_axismach.py:1777`) はここで解決して計画に書き、runner には明示の値として渡す (準備の中で元の run を読み直さない)。
   初期場の適合性 (元の run の輸送種の配置・エネルギー基準 `thermoHrefTemp`・形状と格子の識別) を検査し、合わなければ止める。
2. **準備の後・起動の前 (実効の計画の確定)**: 準備が作った `solverConfig.yaml` (と `PRESTAGE_CONFIG`) から `ns_stage_configs` (forge を回さずに全段の実効の設定を返す、
   `runner_axismach.py:2676`) で**全段の起動直前の設定の全文**を求め、本段の `nStepOuter` は停止規則の上限を入れた形で確定する。境界条件・停止規則の基準・
   節点の表・準備の後の初期場 (元の snapshot とは別の生成物として識別) の sha256 も書く。段の列と、段の間で許す設定の差 (CFL・空間の次数・内反復・step 数・
   出力の間隔・B0 のパス・本段の出力のパス) は**版付きの recipe** (`axismach.wt_card_chain/v1`) が持つ。物理 (入口・壁・軸対称・熱力学・輸送・乱流) の段の間の差は許さない
   (今の NS の段の生成は物理を保って次数・CFL を変える、`runner_axismach.py:2668`)。
3. **実行と照合**: `run_staged_ns` は確定した計画を受け取り、各段の起動の直前に、起動しようとする設定の sha256 が計画と同じかを確かめ、違えば起動しない。
   終わった後、段の記録・起動の記録 (forge の sha256)・準備の記録 (格子・初期場) が計画と同じかを確かめる。承認との照合 (#3) は別に残す。
   起動の後に、走った設定を計画に書き足して「一致」とすることはしない。数値の設定を直すのは、試行の計画を確定する前だけ (確定の後の記録の無い変更は許さない)。

- **search の無い campaign** (`campaign.yaml`、`forge.campaign/v1`): 名前・承認の記録 (パスと sha256)・recipe の名前と版・case のディレクトリ・台帳の置き場所。
  親 plan の search の無い campaign (`evaluations` を持つ形) とは「evaluation 1 つ = 承認した依頼 1 つ」で対応させる。未対応の `search`・知らない recipe・知らないキーは拒否する。
- **親 plan との分担** (Minor の採用): 親 plan §5.1 #6・#7 を、今回の限定の実装 (#6a・#7a: 上の 3 段と search の無い campaign、未知のキーの拒否、入力の参照、全段の照合) と、
  探索の需要が出てからの分 (#6b・#7b: 探索のスキーマ、adapter の `variables()`・`recipes()`・`feasible_range()`、副作用の無い `resolve_evaluation`) に分ける。
  親 plan §4.3 に、今回の限定の実装では準備を先に行い、生成した入力から各段の計画を確定する例外と、その照合の条件を書く。
- **判別の A/B** (実装の後、GPU・forge なし; 起動を記録するダミーの forge): 計画を確定した後、A は変更なし、B は `S1_soft` の設定の CFL だけを計画の値から変える。
  A は通り、B は起動の回数 0 で止まれば、計画と起動が結び付いている。A も止まれば試験の不成立。
- §6 V1a は、この修正だけでは終わらない (基準・節点の表・`RECIPE.json` などの比較と、旧の台本の検査から実行役の関門への対応表が要る; #6 で)。

### 4.2 判定: 要求の一致・解析の成立・要求の達成を分ける

判定は次の 3 つに分ける (codex 2026-10-17)。

| 判定 | 問い | 中身 |
| --- | --- | --- |
| 要求の一致 | 承認した問題を解いているか | 下の照合 (機械的に) |
| 解析の成立 | 計算が成り立っているか | forge-gate の `judge_run.py` (`RUN_VERDICT.json`)・停止規則 v2 (`SETTLE_VERDICT.json`)・壁解像・メッシュ品質・`launch_gate`・受け渡しの関門 |
| 要求の達成 | 要求の公差を満たすか | qualification (風洞ノズルは V5b・V5c・寸法・必須の成果物) |

**要求の一致の照合** (承認した要求 ↔ 生成後の実効の設定・境界条件・外部の物性と分布のファイル):

- 機械的に照合するもの: 入口条件・組成とその基準、壁の熱条件、採用した物理モデル、軸対称と平面の別、対象の形状と解析領域、評価位置、公差、必須成果物。
- 単純な YAML の文字列の比較にしない。導出の関係を扱う (組成の換算、物理スロート径と設計側の `r_throat`、壁温の分布など)。
  親 plan §4.3 の実効の設定と全段の照合の契約を使う。
- 今の `qualify` は `M_exit` だけを照合し、`r_throat` の不一致は警告にとどまる (`wt_qualify.py:949` 以降)。これを照合に置き換える。
- 起動の前に照合し (`execution_plan.json` を確定するとき)、出口で同じ照合をもう一度行う。出口で不一致なら `incomplete` (その run は承認した問題の評価になっていない)。
- 既存の成果物の来歴の照合 (`provenance_check`、`wt_qualify.py:678`: 壁・格子・msh・設定・場・カードの sha256) は残す。足すのは「承認した問題を解いているか」の照合である。

**数値の設定**: CFL・リミッタの値などは照合しない (AI が直してよい)。ただし数値の妥当性の検査は省かない。
既存の `launch_gate` と最終の評価の要件 (例: 本段は 2 次精度、初期化の 1 次精度と区別する; `procedures/solver-settings.md`) は残し、
離散化や格子を変えたら、その変更に要る検証を通す。

**判定器と判定条件の固定**:

- 判定は、承認の時点で固定した版 (commit・依存環境) の別の checkout で行う。作業側の作業ツリーで書き換えても判定に効かないようにする。
  判定器のファイルの sha256 を記録するだけでは保護にならない (内容が同じことを示しても、承認や判定の正しさは示さない)。
- 判定器が読む中間の判定 (例: 寸法の関門は `DIMENSIONS.json` の `verdict` を読む、`wt_qualify.py:760`) も、それを作った処理を固定した版に結び付け、確認できなければ固定した版で再評価する。
- 正式な入口では、CLI・問題 YAML の公差・評価窓・必須成果物 (`optional_parts`) が承認した値と違えば拒否する。今の `resolve_tolerances` は CLI を最優先で採用する
  (例: カード 0.001 に CLI 0.1 を渡すと 0.1 を採用する; codex の関数単体の確認)。開発のための CLI の上書きは、正式な入口の外でだけ使う。

**判定表** (要求の達成; 要求の一致と解析の成立が通ったときだけ評価する):

| 判定 | 条件 |
| --- | --- |
| `accepted` | 要求の公差 (例: 風洞ノズルの V5b・V5c)・寸法・壁解像・必要な成果物のファイルをすべて満たす |
| 補正へ戻す | 評価済みで公差が未達 |
| `rejected` | 反復の上限に達した |
| `incomplete` | 必要な判定・成果物が欠けた、要求の一致が不一致、評価していない要求がある |

- 関門 (`handoff_gate_v2`)・`judge_run` は計算の成立と受け渡しを扱う。公差は別の判定にする。CENSORED の扱い (延長 1 回 / 人に返す) はレシピが持つ。

**判定器の無い要求** (新しい機種・新しい解析領域):

- 必須の要求ごとに「評価方法・公差・証拠・判定の担当」を事前に決めてユーザの承認を得る (§4.7 のユーザ欄に入れる)。
- 評価していない要求が 1 つでもあれば `accepted` を出さない。結果は `incomplete` として未評価の項目を列挙する (親 plan §4.4「要求された評価の未実行は不合格」と同じ)。
- 判定器の無い探索も予算の内で進めて結果を出してよい。人が採用するときは人の判断を別に記録し、数値の VERDICT は書き換えない。
- AI が新しい判定器を実装したら、今の開発フロー (plan・codex レビュー・試験) を通して固定してから使う。専用の昇格の審査は作らない。

### 4.3 失敗の返り値と停止の記録

- **返り値**: 新しい列挙は作らない。既存の語彙の和にする。
  - `RUN_VERDICT.checks`
  - 停止規則 v2 の status + failures
  - `HANDOFF_GATE` の blocked の条件
  - `LAUNCH_GATE` の FAIL
  - `MESH_QUALITY`
  - 要求の一致の不一致、qualification の未達
  - 予算の超過、出力先の衝突
- **停止の記録** (実行役が止まったときに 1 つのファイルに書いて終わる; 名前は実装で決める):
  - 試行 ID・段
  - 判定器と VERDICT、対象量・値・閾値・判定区間
  - 証拠のパスと必要なハッシュ
  - 消費した予算 (その試行・依頼の累計)
  - プロセスの状態
  - 原因の分類は付けられるときだけ付ける (一過性・資源・数値・照合の不一致・**不明**)。原因が分からないものを無理に分類しない。
- **一過性の再試行**: ssh の切断だけで「もう 1 回」を選ばない。§4.5 のプロセスの照合 (生きていれば監視に戻る) を先に行う。
- **既知の復旧の手順**: 実際に使っているものだけをレシピのデータに残す (例: プロセスの照合の後の再開)。網羅的な手の一覧と、手ごとのコストの表は作らない。
- **止まった後**: 直す役の AI は、§4.4 の範囲で原因を調べて直してよい (新しい試行で、記録を残す)。決まった流れを回す役は直さず、停止の記録を人に見せる。

### 4.4 途中の修正の範囲とエスカレーション

- **AI が事前の諮問なしに直してよいもの**: 承認済みの要求・判定条件・予算を保ち、既存のデータを変えず、新しい試行で行う
  設定の調整・格子の作り直し・case の中の探索の台本の変更。直したことは報告に書く (何を試し、何を直し、どの判定で通ったか、判定の出力つき)。
- **AGENTS.md のエスカレーション条件 4 の書き分け** (codex 2026-10-17 の文案、§5.1 #9 で AGENTS.md に入れる):

  > 承認済みの要求・判定条件・予算を保ち、既存データを変更せず、新しい試行で行う設定調整・格子再生成・case内の探索台本の変更は、条件4の事前諮問を不要とする。正式なソルバ・判定器・共通ツールの変更、承認条件の変更、原因や結果解釈の確定には従来の条件を適用する。

- 条件 2 (手順で 2 回対処しても解けない)・条件 3 (事前登録した検証の FAIL)・条件 6 (数値カーネルの変更)・条件 7 (解釈の確定) は残す。
  設計の候補が公差に届かないこと (補正へ戻す) は通常の流れで、開発の検証の FAIL (条件 3) とは区別する。
- **探索** (新しい機種・道具の新しい組み合わせ): case の下の探索用のディレクトリで台本を書いて試す。forge は既知の流れと同じ実行の入口 (§4.5) から起動し、
  予算・試行 ID・新しい出力先を必須にする。判定器の無い要求は §4.2 のとおり `incomplete`。同じ台本を 2〜3 回使ったら今の開発フローで道具にする。
- `design-intake` の未対応のメニューの分岐 (§6) では、「既存の道具の新しい組み合わせによる探索」と「未実装の物理の機能の開発 (plan を起こす)」を分ける。

### 4.5 実行の入口・予算・台帳・再開 (2026-10-16 codex M3 + 2026-10-17 codex)

- **実行の入口**: 既知の流れと探索の台本の両方が使う小さな入口を 1 つ置き、`run_case.sh` を包む。承認された予算の ID・試行 ID・新しい出力先を必須にする。
- **予算** (ユーザ決定: 依頼ごとに見積もって承認):
  - 起動の前に予約する (その試行の見積もり)。台帳の残りが足りなければ起動しない。
  - 実行中に監視し、承認値を超えたら止める。終了時に精算する。
  - 失敗した試行・時間を測る run・前の段も計上する。再開で消費量を戻さない。
  - 単位は GPU 時間 (GPU 台数 × 経過時間) とディスクの増分を案とし、実装の時に決める。
  - 超えそうなときは次の計算を始めずに止めて人に聞く (追加の承認は台帳に記録する)。
- **実行の場所**: 実行器は run と同じホスト (AWS) で、試行ごとの service (下) として動き、台帳もそこに置く。手元は rsync の写しを読むだけ (`forge-design status`)。
- **台帳**:
  - 単一の書き手のロック、原子的な更新 (一時ファイル → rename)。
  - 段の試行ごとの attempt ID と、起動したプロセスとの結び付け (PID・開始時刻・ホスト)。
  - 再試行の回数と、予算の予約・消費。再開の後も保持する。
- **試行ごとの実行担当** (2026-10-17 codex plan 段 Major 1): 試行 1 つに実行担当のプロセスを 1 つ置き、runner・停止規則の監視・予算の監視をまとめて持たせる。
  - 今の停止規則 v2 は forge を別のセッション (`start_new_session=True`) で起動し、監視の状態と `Popen` を呼び出し元が持つ (`stop_rule.py:1901`)。
    監視の側が死んでも forge は残りうる。同じ関数の再呼び出しは既存の出力を見て拒否する (`stop_rule.py:1879`)。したがって PID の生存の確認だけでは監視を戻せず、
    再開までの予算の超過も防げない。
  - CLI (`forge-design run --resume`・`status`) は実行担当に接続するだけにする。CLI が死んでも実行担当は監視と予算の監視を続ける。
  - **停止の責任は実行担当の外に置く** (2026-10-17 codex 諮問 [plan2-disposition](../../notes/reviews/2026-10-17-orchestrator-plan2-disposition-diagnose.md)):
    今の停止規則 v2 は例外のときの forge の停止も同じ監視のプロセスが行う (`stop_rule.py:1980`) ので、監視そのものへの `SIGKILL` では停止の処理が動かない
    (codex が同じ `start_new_session=True` の構成をダミーの子で確かめ、親の `SIGKILL` の後も子は生きていた)。死んだ実行担当に後始末を期待しない。
  - **契約**: 試行ごとに systemd の service を 1 つ作り、実行担当とその全子孫 (`run_case.sh`・forge) を収める。実行担当を main process とし、
    `ExitType=main`・`KillMode=control-group`・有限の停止期限 (`TimeoutStopSec`)・`Restart=no` にする。main process が終わると service が止まり、残った子孫ごと止める。
    これで停止規則 v2 の `start_new_session=True` はそのまま使える。`nohup` での起動はやめる。
  - 停止期限の内に全子孫の終了を確かめられなければ、その試行の新しい試行を起動しない。確かめられるまで予算の予約を解放しない (精算は確認の後)。
  - 実装の最初に、対象の環境 (AWS・手元) で service の管理 (systemd の user service など) が使えるかを検査する。使えなければ `nohup` に黙って戻さず、起動を止める。
  - **Linger が要る** (2026-10-17 実測、[記録](../../notes/investigations/2026-10-17-attempt-service-ab/README.md)): `Linger=no` (AWS B・手元の既定) では、ssh を切ると約 12 秒後に
    systemd が user service を止めた。`Linger=yes` では 25 秒後も生きていた。実行役の事前検査は `Linger=yes` を要り、無ければ起動を止める。
    インスタンスで `sudo loginctl enable-linger ubuntu` を 1 回行う (ユーザに確認してから; 試験の後は `Linger=no` に戻してある)。
  - 判別の A/B (GPU・forge なし; 実装の先頭、§5.1 #15): 同じ service の中で別のセッションのダミーの `run_case` とその子のダミーの forge を動かし、100 ms ごとに生存を観測する。
    A は介入なしで 10 秒、B は 2 秒の時点で実行担当だけを `SIGKILL`。停止期限は 5 秒に固定する。A が続き B の全子孫が期限の内に消えれば契約の成立、B が残れば不成立。
    A も止まったら試験の不成立。
  - 実行担当が死んだ試行は再接続しない。完了の証拠 (`STOPPED_AT`・`MAIN_STAGE_STARTED`・関門の記録) を確かめ、足りなければ新しい試行 (新しい run 名) に引き継ぐ。
    既存の runner を再び呼び出さない (開始・停止の記録や途中の出力を消す `runner_axismach.py:2946`・`:2981`、来歴とログを上書きする `run_case.sh:52`)。
  - 不明なら人に返す。
- **既存のデータの保護**: 出力先が既にあれば止める (今の道具の振る舞いを入口でも検査する)。確定済みのデータ (accepted の run・成果物・入力の参照) への書き込みの権限を分ける。chmod は補助。
- **巨大な出力の扱い**: sha256 で照合しない。
  - 格子は不変部分のハッシュ (forge-gate `mesh_provenance`)。
  - 全場は停止の場だけ sha256、ほかは大きさ + 完了の印。
  - 既存の記録は sha で参照し、中身を複製しない。
- **idle 自動停止との関係**: 実行器の準備の隙間が「forge なし」と数えられる。GPU の通しの検証の前に決める (§5.1 #7)。

#### 4.5.1 実行役の作り (§5.1 #6 の設計、2026-10-17 codex 諮問で修正)

codex (diagnose) の諮問 ([記録](../../notes/reviews/2026-10-17-orchestrator-executor-diagnose.md)、C0/M5/m1) で、最初の案の recipe の分岐・driver の復旧・
予算の停止の順序と測る範囲・`fixed` の照合を直した。

**2 層の service**:

- **driver** (campaign 1 つに 1 つ、systemd の user service): 台帳の**唯一の書き手**。生きている間ずっと台帳のロック (`<ledger>/ledger.lock` の flock) を持ち、
  更新は一時ファイル → rename。二つ目の `--resume` はロックを取れないので、今の driver を止めずに状態を返す。
- **試行** (段の試行 1 つに 1 つ、別の user service、§5.1 #15 の契約): 実行担当が段の処理と予算の監視を持ち、自分の状態 (`attempts/<試行 ID>/status.json`) と
  停止の記録 (`STOP_REASON.json`) とログだけを書く。driver とは別の unit なので、driver が死んでも試行は続く。
- **識別**: campaign の ID (初めて起動したときに driver が作る UUID、台帳に保存) と試行 ID で unit 名・出力先を決める。campaign の名前は表示だけに使う
  (同じ名前の別の case の campaign と衝突させない)。
- **起動の順序** (一度だけの処理): (1) 試行 ID・unit 名・run の出力先・予算の予約を台帳に書く → (2) service を起動する → (3) 終わった試行の結果を、試行 ID ごとに
  **一度だけ**取り込んで精算する (取り込み済みの印を台帳に書く)。run の番号の割り当ては、case の README の更新と同じ case のロックの下で行う。
- **再開** (`--resume`、driver を起動し直す): 台帳の各試行を状態で照合する。unit が動いている → 待つ。止まっていて試行の状態が終わり → 取り込む。起動の途中・停止の途中 →
  期限つきで待つ。unit の照会に失敗した・unit が無くて完了の証拠も無い → 人に返す。同じ試行の段の処理は呼び直さない (新しい試行・新しい run)。
- **ログ**: driver の起動ごと・試行ごとに別のファイルに書き、上書きしない。journal だけを復旧の根拠にしない。
- **変更の要求**: CLI (`note` など) は要求のファイルを置くだけで、driver が入力の sha256 と理由と一緒に台帳に取り込む (台帳の書き手を 1 つに保つ)。
- 事前検査: systemd の user service が使えること・`Linger=yes`。満たさなければ起動しない。

**recipe `axismach.wt_card_chain/v1` の段と分岐** (今の手順書 §2〜§9 と同じ; codex Major 1 で直した):

| 段 | 処理 (今の道具) | 次の段 |
| --- | --- | --- |
| `problem_p{k}` | `wt_card`、`approval check` (k ≥ 1 は試行の入力) | 終了 0・寸法 PASS・照合 OK → `mesh_p{k}` |
| `mesh_p{k}` | 格子の段 (§4.6; 今は case/71 の `hand_mesh.py build/check`) | `VERDICT: PASS` (NS・Euler) → `euler_p{k}` |
| `euler_p{k}` | `runner_axismach` (Euler) を新しい run で | forge 終了 0・`VERDICT: PASS`・非有限 0 → `ns_p{k}` |
| `ns_p{k}` | `deltastar_loop --init-integral` (Euler の参照・初期場を固定、§4.1.1 の計画を確定して `run_staged_ns(plan=…)`) | 停止規則 v2 REACHED・受け渡しの関門 進む・`check_settled` SETTLED・非有限 0 → k = 0 は `calibrate_p0`、k ≥ 1 は `qualify_p{k}` |
| `calibrate_p0` | **最初の較正** (`contur_calibrate`、手順書 §6) → 係数だけを変えた試行の入力 `card_p1` | 終了 0・`approval check` OK → `problem_p1` (較正し直しの回数に数えない) |
| `qualify_p{k}` (k ≥ 1) | `gate_env judge` (承認のモード) | `incomplete`・`make_deliverables` → `deliverables`。`back_to_correction`・`recalibrate` → `recalibrate_p{k}`。`accepted` (成果物が既にある) → 終わり。ほか (`escalate`・`rejected`・ほかの `incomplete`) → 止まる |
| `recalibrate_p{k}` | 判定の記録の係数 (`wt_qualify/contur_overrides.yaml`) → `card_p{k+1}` | 照合 OK → `problem_p{k+1}` (較正し直しの回数 +1; 上限は承認の `run_policy`) |
| `deliverables` | `gate_env deliverables` → `gate_env judge` (最後の判定) | `accepted` → 終わり。ほか → 止まる |

- 反復の上限は承認の `run_policy.max_iterations` (pass の数)、較正し直しの上限は `max_recalibrations`。上限で止まる。
- 成果物の前で止まるのを避けるために、必須の成果物を任意にしない。
- case の README の「計算 run 一覧」に、その campaign の行を印の間に台帳から作って書く (ほかの行は残す、case のロック)。

**予算** (承認の記録の `budget`; 一律の上限は置かない (ユーザ決定)):

- **GPU 時間**: GPU の段の forge の起動ごとの経過時間 (GPU 1 台) を、起動の時刻と終わりの時刻で台帳に残す (前の段の起動・時間を測る run・失敗した試行も数える)。
  実行担当が死んで終わりの時刻が無い区間は 0 にしない (止まったと分かった時刻か出力の最後の更新の時刻までで数え、根拠を書く)。AWS のインスタンスの稼働時間とは別の指標。
- **ディスク**: この依頼が新たに占めた出力 (campaign が作った run・台帳・生成したカード・問題 YAML・格子など、台帳に登録した**重ならない範囲**) の大きさ。
  後の段が同じ run に成果物を足しても二重に数えない。数分おきの監視は超過の検出で、厳密な上限ではない (最大の出力の大きさ・監視の間隔・止めるまでの遅れぶんの余裕と、
  停止の記録を書く空きを先に確保する)。
- **予約**: 起動の前に、その段の見積もり (`design-intake` §8 の見積もりを台帳に写したもの) が残りの予算の内かを確かめ、残り全部をその試行の上限として予約する。
  見積もりとの差は記録して警告に使う (固定の倍率で止めることはしない)。承認の総額を超えると見込まれるときは次の段を始めずに返す。
- **予算で止める順序** (codex Major 3): 予算に届いたこと・測った値・停止の要求を**先に**原子的に記録し、その後で試行を止める。停止の要求と「全子孫の終了の確認」を分け、
  driver が終了を確かめてから精算する。実行担当が `SIGKILL` で記録を書けなかったときは、再開した driver が根拠つきで補い、原因が確かめられなければ「不明」とする。

**役割** (campaign に固定、`--resume` で変えない):

- `fixed` (知識の少ないモデル): 段の処理 (生成器・台本) は**承認の commit の checkout** (§4.7.1 (5) の判定専用の環境と同じ checkout) から回し、生成器の版・依存する台本・
  外部の入力を固定する。一度生成した入力の sha256 と生成元 (commit・台本の sha256) を台帳に書き、使う直前に照合する。作り直して比べる照合は §6 V1a の検証に残す。
  止まったら停止の記録を残して終わり、人に見せる。強制は #14 まで運用の検査であって権限の強制ではない。
- `fixer` (賢い AI): 段の処理は作業ツリーから回してよい。止まった後に AI が §4.4 の範囲で直し、変更の要求 (入力の sha256 と理由) を出して `--resume` で新しい試行から回す。
  保護する条件・判定・予算は #2・#3・#4 の照合が守る (判定は判定専用の環境)。

**#6 の分け方と順** (codex: #6a → #6b → #6c → #6d、最小の driver と復旧は #6b に前倒し):

- #6a: 台帳 (ロック・原子的な更新・campaign の ID・試行の記録・取り込みの一度だけ)・停止の記録・予算の勘定 (GPU なし)。
- #6b: 試行の service・予算の監視・最小の driver と復旧 (ダミーの forge で障害の試験)。
- #6c: recipe の段と分岐・CLI・README の同期・変更の要求。
- #6d: §6 V1a (case/71 の 1 pass の決定的な入力を V6 の run と比べる、旧の台本の検査と関門の対応表)。
- GPU を使う通し (§6 V1b・V4 (a)) は #7 の後 (手元の GPU で前倒ししない)。#6 を先に作ったことを、未完の #5 (forge-gate の関門)・#8 (格子の上限の照合) を省いて
  通しを合格にする理由にしない。

### 4.6 格子の段

- 契約: 「問題 YAML の mesh 節 (表を含む) + 変換後の格子の不変部分のハッシュ + 品質の VERDICT + **依頼ごとの承認済みの節点数の上限と承認元**」。
- **節点数の上限** (2026-10-17 codex plan 段 Major 2): 生成器が算出する必要な節点数と、ユーザが許した上限を分け、格子を作る前に照合する。
  - 今の `hand_mesh.py` は点数から `max_nodes` を自動で広げ、理由に固定の文「ユーザ決定『a』2026-10-15」を書く (`hand_mesh.py:354`)。
    `runner_axismach.node_budget` は上限が正の整数・理由が空でないことしか見ない (`runner_axismach.py:790`)。そのまま包むと、別の依頼でも過去の承認の文で上限を通せる。
  - 上限の承認はカードのユーザ欄 (§4.7) に入れ、承認元 (ユーザの指示の記録・カードの版) に結び付ける。case/71 の承認は case/71 のカードだけに効く。Euler と NS の格子の両方を照合する。
  - 規則は `plans/accepted/tooling-nozzle-mesh-size-guard.md` §4.3 (上限を超えるときはユーザの了承を取ってから書く、AI が了承なしに書かない) と同じ。
- 実装は 2 つ:
  - 雛形: `wt_card.estimate_mesh` → mesh2d (case/70 など)。
  - 手組み: `hand_mesh.py` の格子の生成の数値の処理を再利用する (case/71)。承認と上限を作る部分 (固定の承認の文・点数からの上限の拡大) は #8 で、カードの承認済みの上限との照合に置き換える。
- 係数の一般化と生成器への組み込みは別のマイルストーン (#10)。

### 4.7 要求カードと承認の固定 (2026-10-16 codex M5 + 2026-10-17 codex)

- カードの値に出所を付けて保存する。出所は 4 種類: ユーザの入力の参照 (発言・指示の記録)、レシピの既定、導出、未確認の提案。
- ユーザの指定が要る欄 (公差・壁温・寸法・成果物の選択・格子の節点数の上限 (§4.6)・判定器の無い要求の評価方法など、欄ごとの表で決める) は、ユーザの指示か確認済みのカードの版に結び付ける。
  未確認なら起動の前に止める。レシピの既定の値も、ユーザに示して承認した版にだけ使う。
- **判定条件の節**: カードに `qualification` (公差・評価窓・必須成果物) と予算の承認値を書けるようにする。今の書式 `forge.wt_card/v1` は `qualification` を受け付けない
  (`wt_qualify.py:54`) ので、書式の版を上げる。承認したカードの版と解決済みの判定条件は campaign・`execution_plan.json` に結び付ける (§4.1)。
- 文章の決まりの置き換え:
  - `design-intake` §7「既定から変える場合は理由を YAML コメントに残す」→ 承認した版への結び付け。
  - 手順書 `design-chain-quickstart.md` のカードの「寸法の公差 (聞いて無ければこのまま)」(277 行) → 承認済みの既定値に限る。
- 負例: 欠けたカードだけでなく、「値はあるが根拠が無い」カード。
- **依頼の原文の置き場所 (2026-10-17 ユーザ決定「そもそも機微な情報はそもそも入れないので、気にしなくていい」「このレポジトリを取り込んでから、個別のgitlabにuploadする想定であり、全く問題ない。世の中に公開のレポジトリには載せていい情報しか載せない」)**:
  依頼の原文・やり取り・AI の報告は、**設計を行うリポジトリの中**に、カード・判定・成果物と同じ場所 (case のディレクトリ) に置く。forge には載せてよい情報だけを置き、
  forge を取り込んだ個別のリポジトリでも同じ置き方で使える。

#### 4.7.1 承認する実効条件・試行の入力・判定専用の環境・中間の判定 (§5.1 #2 の設計、2026-10-17 codex 諮問で修正)

最初の案 (カードの節で承認の範囲を決め、カード全体の sha256 で固定し、別 checkout と `pip freeze` で判定の版を固定し、寸法だけを測り直す) は、
codex (diagnose) の諮問 ([記録](../../notes/reviews/2026-10-17-orchestrator-approval-design-diagnose.md)) で次のように直した。

**(1) 承認の対象は「生成した後の実効の条件」(保護する条件) で決める** (カードの節で決めない):

| 種類 | 保護する条件 |
| --- | --- |
| 要求 | 組成とその基準・入口の Pt・Tt・M_exit・物理スロート径・入口の形 (半径・直管・縮流部)・拡大部の長さ・壁の熱条件 (壁温の値と分布)・背圧・寸法の公差 |
| 物理モデル・物性 | ガスのモデル・熱力学の出典 (種の DB と sha256)・輸送 (粘性・熱伝導)・乱流モデル・軸対称の別・壁の境界条件の種類 |
| 評価 | 評価の位置 (x_E・線の η・x の窓)・点数 (`line_n`)・量の定義 (V5b のコアの範囲・V5c)・公差・定常性の条件 (停止規則のプリセットと解決した閾値・時間窓)・必須の成果物 |
| 実行の方針 | 反復の上限・較正し直しの上限 (開始のときに解決して固定; 済んだ回数は台帳から)・予算・格子の節点数の上限 |
| 判定器の無い要求 | `requirements_without_gate` (削除も承認の条件の変更として拒否) |

- `recipe`・`overrides` は、保護する条件を変えない範囲で選んでよい。今の `_OVERRIDE_FORBIDDEN` は範囲が足りない (`evaluate.settle`・`gas.transport` を通す; codex の関数単体の確認) ので、
  上書きの可否で守らず、生成した後の保護する条件を承認の値と照合して守る。
- 今の `wt_qualify` は `line_n` を run の基準から取る (`wt_qualify.py:203`)。同じカードでも基準しだいで 2001 → 2 に変わる (codex の関数単体の確認)。`line_n` も承認の値で照合する。
- `wt_qualify` は `iteration`・`max_iterations` の欠落を許し上限を判定しない (`wt_qualify.py:898`)。正式な入口では必須にする。

**(2) 承認の原本と試行の入力を分ける**:

- 承認の記録 (`card_approval.json`) は、承認したカードの原本 (写しと sha256) と、解決した保護する条件の値を持つ。
- 試行の入力 (例: 較正の係数を足した `card_p1.yaml`、手順書 §7.2) は、承認の記録への参照と自分の sha256 を持つ。検査は「試行の入力から解決した保護する条件 = 承認の値」で行う。
  格子・較正の係数だけの変更は再承認なしで通り、保護する条件の変更は起動の前に止まる。承認の記録を自動で更新して不一致を避けることはしない。

**(3) カードの版 2 と出所**: v1 の欄に、省略可の節 `qualification` (公差・評価の位置と点数・必須成果物)・`budget`・`mesh_limits`・`requirements_without_gate`・`sources` を足す。

- `sources` は節 (`conditions`) か葉 (`conditions.Pt`) で書いてよいが、承認のときに実在する葉へ展開し、未知のパスは拒否する。
- `recipe_default` はラベルだけで認めず、固定した既定の提供元 (雛形のファイル・キー・値) と照合する。
- CSV から作る表・ハッシュなどの導出の葉は、元のユーザの指定 (CSV のパスと sha256) に結び付け、「書かれなかった既定」と取り違えない。
- 正式な入口では予算 (`budget`) が必須。

**(4) 承認の記録の書き方**: `approval preview` が表示する表 (保護する条件の値と出所) をファイルに固定し、`approval write --preview <そのファイル>` が同じ内容を書き出す
(表示の後に解決し直して値が変わるのを防ぐ)。AI がユーザに表を見せて承認を得てから write する。sha256 と `user_ref` は承認の真正性を証明しない。
承認の記録そのものを書き換えられる構成では、保証は「指定した承認の記録との整合」に限る (強制は #14)。

**(5) 判定専用の環境**: 別 checkout の HEAD・追跡ファイル・`pip freeze` だけでは、実際に使うコードと道具を限定できない
(`FORGE_MESH_PY` で別の Python を選ぶ `deliverables.py:48`、STEP の読み取りは FreeCAD `dimension_gate.py:50`)。

- 判定専用の環境を作る: 固定した commit の checkout と、ハッシュつきの依存の一覧から作った venv。外部の Python・FreeCAD・Gmsh・変換器の実体 (パスと sha256・版) を記録して照合する。
- 判定は作業側の環境変数 (`PYTHONPATH`・`FORGE_*` の道具の選択) を引き継がずに起動する (隔離のモードの Python、許す変数だけを渡す)。import の元が checkout の中であることを確かめる。
- コードの探索先 (`design/forge_design` など) の未追跡・ignored のファイルも見る (`.py` が足されていたら止める)。
- `pip freeze` の sha256 は補助の記録にとどめる (別の Python・外部の実行ファイルを捉えず、無関係のパッケージの変更でも不一致になる)。

**(6) 合否に使う中間の判定を全部再検査する** (寸法だけでは足りない):

- 対象を列挙する: 成果物の各部品 (CAD の転送・回転面・独立の読み直し、mesh の再構築と品質、Fluent、Excel、壁の比較、寸法)、停止・整定の判定 (`SETTLE_VERDICT`)、壁解像、メッシュ品質、起動の関門。
- 各々を固定した環境で実行し直すか、固定した環境で作られ入力に結び付いていることを確かめる。どちらもできないものは accepted にしない。
- 必須のファイルの一覧は MANIFEST の自己申告でなく、固定した版のコードが決める。総合の PASS と個別の検査の FAIL・欠落が矛盾したら accepted にしない
  (今の `deliverables_check` は総合の `verdict` と `failed_parts` を読み、`cad.pass=false` でも PASS になる; codex の関数単体の確認)。

### 4.8 入口・道具の説明・skill・手順書・AGENTS.md

- **入口**: CLI は次の 2 つ。
  - `forge-design run <card> --role fixed|fixer [--resume] [--executor local|aws]`
  - `forge-design status <ledger>`
- **道具の説明**: 各道具の `--help` と冒頭の説明に、何をするか、入力と単位、出力、終了コード、所要時間の目安を書く。索引は `design/CAPABILITIES.md`。目録の専用ファイルは作らない。
- **skill**: 薄くする。中身は「4 つの変えてはいけないもの」「道具の索引」「カードを作る (出所つき)」「回す」「報告する」。
- **手順書**: `design-chain-quickstart.md` の「この手順書の使い方」3 (自分で直さない、16〜17 行) は決まった流れを回す役 (`--role fixed`) にだけ適用すると書き分ける。
- **AGENTS.md**:
  - 規則は「レシピのある機種 (最初は風洞ノズル) の設計の run は `forge-design` 経由」に限る。
  - エスカレーション条件 4 の書き分け (§4.4 の文案)。
  - 同じ節を書き換える plan が 3 つある。書き換えは 1 セッションが次の順でまとめて行い、調整板に予告する: settle (収束確認、済み) → gate (起動の関門) → この plan (入口)。
- **集計**: 依頼 1 件あたりの、試行の数・人に返した回数・直した回数 (直した種類ごと)・消費した予算を台帳から出す。
- **case の run 一覧** (2026-10-17 codex plan 段 Minor 5・諮問で補足): 実行担当は run を作ったとき・終えたとき・破棄したときに、case の `README.md` の「計算 run 一覧」に、その campaign の行を台帳から生成して書く
  (他の run の行は残す)。自動の再試行・時間を測る run も台帳だけに残さない (AGENTS.md「case の run 索引」)。実行担当が死んで書けなかった分は、再開のときに台帳から直す。
  複数の campaign が同じ README を書くので、case ごとのロックと原子的な更新で他の行を保つ。

### 4.9 今は作らないもの

どれも、必要が測れたら作る。

- 判断役の専用の呼び出し (digest・引用の照合・安全集合・問いの組・`local` の backend・`forge-design answer`)。再開の条件: 知識の少ないモデルに止まった後の判断をさせる要求が出たとき。
  そのときも選べる手を「もう一度 (一過性のときだけ)」と「人に返す」に絞る制限 1 つから始める。
- 網羅的な手の一覧と手ごとのコストの表 (§4.3)。
- ツールごとの専門家の AI、MCP のサーバ。
- 検証済みの範囲の一覧と、その判定の規則 (実績が数件たまったら、accepted の結果を集めて表にするスクリプトを 1 本置く)。
- 目録の専用ファイル、新しい台帳の形式、探索から recipe への専用の審査。

## 5. 実装ステップ

codex (2026-10-17) の優先順に従う。ただし先頭に、試行ごとの service の契約の事前検査と障害の A/B (GPU なし) を置く (2026-10-17 codex 諮問)。

0. service の管理の利用可否の検査と、`SIGKILL` の判別の A/B (§4.5)。
1. 承認した要求と判定条件の固定: 保護する条件の解決と承認の記録 (#2a) → 判定専用の環境 (#2b) → 中間の判定の再検査 (#2c) の順 (§4.2・§4.7.1)。
2. 要求の一致の照合と、親 plan の最小の契約・出口の統合: 入力の参照・実効の設定・`execution_plan.json`・forge-gate の `judge_run`・qualification (§4.1・§4.2)。
3. 実行の入口・試行ごとの実行担当 (service)・予算・新しい出力先・台帳・再開・停止の記録を一体で (§4.3・§4.5)。
4. 既知の線形の流れと格子の包み、最小の CLI。V1a・V2・GPU なしの V3 を通す (§4.6・§4.8)。
5. AWS の idle 自動停止を解決してから、GPU の通しの検証 (V1b・V4)。
6. 道具の `--help`・skill・手順書・AGENTS.md の同期 (§4.4・§4.8)。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | ~~設計方針の諮問 (diagnostician、ユーザ指定)~~ (2026-10-16 済み) | 判断: 2026-10-16・方向は妥当、新しい枠組みにせず campaign-recipe の実行器 + ledger と forge-gate の判定・関門の上に載せる (全件採用) | F |
| 1b | ~~codex の plan 段レビュー~~ (2026-10-16 済み) | GO-with-changes、C0/M5/m2、全件採用 (§6.1) | F |
| 1c | ~~改訂の方向の諮問 (codex diagnose)~~ (2026-10-17 済み) | 判断: 2026-10-17・方向は採用、ただし出口の判定だけでは 4 つを守れない。承認した要求と判定条件の固定・予算と既存データの実行時の保護を足し、途中の修正はその範囲で AI に任せる (全件採用、§4 を書き直し、§6.1) | F |
| 1d | ~~書き直した §4・§6 の codex の plan 段レビュー~~ (2026-10-17 済み) | 判断: 2026-10-17・GO-with-changes、C0/M4/m1。採否は codex (diagnose) に諮って全件採用、M1 は停止の責任を試行ごとの service に置く契約まで plan で決めた (§4.5・§6.1) | F |
| 15 | ~~試行ごとの service の契約の事前検査と障害の判別の A/B~~ (2026-10-17 済み、PASS) (§4.5、GPU・forge なし; 実装の先頭) | 対象の環境 (AWS B・手元) で systemd の user service (`ExitType=main`・`KillMode=control-group`・`TimeoutStopSec`・`Restart=no`) が使えるかを調べ、§4.5 の A/B を回す。合格: A は 10 秒続き、B は実行担当の `SIGKILL` の後 5 秒の内に全子孫が消える。使えない環境は「起動を止める」とし、別の仕組みは諮ってから。**結果 (2026-10-17)**: 手元 (WSL2)・AWS B とも systemd 255 で A PASS・B PASS (全子孫が 0.10 秒で消えた)、対照 (service なし) は子孫が残った。`Linger=no` では ssh を切ると service が止まる (§4.5)。台本 `design/tests/run_attempt_service_ab.py`、記録 `notes/investigations/2026-10-17-attempt-service-ab/` | O |
| 2a | ~~保護する条件の解決と承認の記録~~ (2026-10-17 実装・試験 ALL PASS) (§4.7.1 (1)〜(4)) | カードの版 2 (`qualification`・`budget`・`mesh_limits`・`requirements_without_gate`・`sources`)、保護する条件の抽出 (カード + 生成した問題 YAML + 雛形の既定)、`approval preview / write / check`、試行の入力 (`card_p1.yaml` など) から承認の記録への参照、正式な入口での CLI・問題 YAML の相違の拒否、`requirements_without_gate` は `incomplete`。合格 (正例と負例を並べる; 全部止める実装を合格にしない): **正例** 承認したカードと、格子・較正の係数だけを変えた試行の入力が再承認なしで通る。**負例** 根拠の無いユーザ欄・出所の未知のパス・既定との不一致・予算の欠落・`requirements_without_gate` の削除・CLI での公差・窓・必須成果物の指定 (値が同じでも)・承認の後の `line_n`・定常性の窓と閾値・物性の変更、で起動の前に止まる。codex の判別 A/B 2 つ: (i) 承認値 0.001 に CLI 0.1 → 拒否、(ii) `line_n` 2001 を承認し、生成した基準の `line_mach.n` だけ A 2001 / B 2 → A 通過・B 停止 (A も止まれば試験の不成立) **結果 (2026-10-17)**: `design/forge_design/intake/approval.py` (preview・write・check・check_run)、`wt_card` の版 2 の検査、`wt_qualify --approval` (`approval_mode`・`approved_tolerances`)。試験 `design/tests/run_approval_tests.py` (固定の入力 `design/tests/fixtures/approval_trial_a/`) ALL PASS: 正例 4 (承認・原本・較正の係数だけの試行の入力・run の基準)、負例 N1〜N15、判別 A/B (i) CLI 0.1 → 拒否 (承認なしの今の経路は 0.1 を採用)、(ii) A (2001) OK・B (2) MISMATCH。実物の run_0003 の基準は承認より古い生成器のもの (`residual_finite.required` が無い) で MISMATCH になる (N15)。既存の試験 9 本 (wt_qualify 176 件・wt_card・mknot_v2_card・deliverables ほか) は通る。物理の照合のうち run の設定 (solverConfig・bcond・物性の DB) は #3 | O |
| 2b | ~~判定専用の環境~~ (2026-10-17 実装・試験 ALL PASS) (§4.7.1 (5)) | 固定した commit の checkout とハッシュつきの依存の一覧の venv、外部の Python・FreeCAD・Gmsh・変換器の実体の記録と照合、環境変数を引き継がない起動、import の元と未追跡・ignored のファイルの検査。合格: 正例 同じ環境で判定できる。負例 同じ HEAD・追跡ファイルでも外部の Python・import の元が違う、探索先に `.py` が足された、で判定の前に止まる **結果 (2026-10-17)**: `design/forge_design/feedback/gate_env.py` (承認の commit の sparse checkout `design`・`solver_density_cuda/tools`・`solver_density_cuda/data` 約 8 MB、判定の Python の指紋 (実体・prefix・配布物ごとの RECORD の sha256・editable)、外部の道具 (gmsh 用の Python・FreeCAD・変換器) の実体と sha256、隔離のモードと許した環境変数だけでの起動、`wt_qualify --approval` は印・import の元・checkout の HEAD を確かめる)。承認の記録の `gate` は sha256 で守る。試験 `design/tests/run_gate_env_tests.py` ALL PASS: 正例 3 (check・probe・judge が門を通り場の読み込みで止まる)、負例 8 (別の Python・印なし・import の元が外・探索先の .py・道具の中身・作業ツリーから直接・gate の書き換え・未 commit の preview)。最初の judge は物性の DB が checkout に無くて止まった (試験の条件が甘く通っていた) → data を足し、条件を「場が無いことで止まる」に締めた | O |
| 2c | ~~中間の判定の再検査~~ (2026-10-17 実装・試験 ALL PASS) (§4.7.1 (6)) | 合否に使う中間の判定の列挙と、固定した環境での実行し直し (または結び付けの確認)、必須のファイルの一覧を固定した版のコードから。正式な入口では `iteration`・`max_iterations` を必須に。合格: 正例 `run_0014_chain_v2` が accepted のまま。負例 総合 PASS と個別の FAIL・欠落の矛盾 (例 `cad.pass=false`)、必須のファイルが列挙から消える、で accepted を出さない **結果 (2026-10-17)**: (1) 両方のモード: `deliverables_check` が必須の部品の `pass` (cad・mesh・fluent; `deliverables.PASS_REQUIRED`)、MANIFEST の総合と部品・寸法の矛盾、必須のファイル (`deliverables.REQUIRED_FILES`、MANIFEST の列挙と別) を見る。(2) 承認のモード: 停止の判定 (`check_settled --end-step <停止> --no-write`) と計算の格子の品質 (承認で固定した AR・skew) を判定専用の環境で回し直す (`approval_rechecks`)。成果物は MANIFEST の `gate` (判定専用の環境で書き出した印; `gate_env deliverables`) が承認と合うものだけを使う。壁解像・V5b・V5c・線の定常は元から判定の側で計算している。試験: `run_wt_qualify_tests.py` ok 184 件・FAIL 0 (4c の負例 7: cad.pass=false・mesh の pass なし・必須のファイルの欠け・gate なし/違う)、`run_approval_rechecks_tests.py` (run_0014: 正例・停止 20000・AR 10・停止なし)、`run_gate_deliverables_tests.py` (run_0014 の写しから判定専用の環境で成果物一式が PASS・gate が承認と合う・判定の側が PASS)。実物の run_0014 は承認なしの経路で accepted のまま。**残り**: 起動の関門 (`launch_gate`) は forge-gate との統合 (#5) の後。承認のモードで accepted まで通る正例は、承認の後に作った run が要る (§6 V1b・V4 (a)) | O |
| 3 | ~~要求の一致の照合~~ (2026-10-17 実装・試験 ALL PASS; 起動の前の呼び出しは #6) (§4.2) | 承認した要求 ↔ 実効の設定・境界条件・物性と分布のファイル。導出の関係 (組成の換算・物理と設計の `r_throat`・壁温の分布) を扱う。起動の前と出口の 2 回。合格: case/71 `run_0014_chain_v2` と試走 A (`case/90.qs_trial_a`) で一致、入口の Tt・壁温・組成・軸対称の別をずらした負例で止まる **結果 (2026-10-17)**: `approval.check_run_physics` (本段の solverConfig・bcondConfig・準備の記録 ↔ 承認と照合済みの問題 YAML): 入口の Pt・Tt・組成 (輸送種の配置の質量分率; `cfd_initial_line` の規則を使い回す)・熱力学 (thermalMethod・species・thermoHrefTemp)・種ごとの輸送と粘性のモデル・NS の乱流 (承認の時点の `runner_wt._config_sst_node` の雛形の値を承認に固定)・軸対称・壁 (一様の等温・分布 (準備の記録)・断熱)・出口の静圧。`check_run` に入れ、`wt_qualify --approval` が問題 YAML を渡す。試験: `run_approval_tests.py` (試走 A の正例、物理の変更 14 種で MISMATCH、数値の設定 (CFL・リミッタ) だけの変更は通る)、`run_approval_rechecks_tests.py` (case/71 run_0014 の正例 2 例目、停止規則の基準 wt_nozzle_v2 も一致)。**正式な入口の通しの正例** (`run_gate_deliverables_tests.py` (P4)): run_0014 の写しで、判定専用の環境で成果物を作り `wt_qualify --approval` が accepted (出口コア M 5.502842、V5c 0.0623 %・0.0573 %、保存済みの判定と同じ)。途中で FreeCAD の中の Python が checkout に `__pycache__` を書いて判定が止まったので、照合の前に `__pycache__` だけを消すようにした (3d840047)。入口の分布 (`inletProfile`) の run は照合の対象外で止める (今の風洞ノズルは一様入口) | O |
| 4a | ~~初期場の固定・継承する値の解決・適合性~~ (2026-10-17 実装・試験 ALL PASS) (§4.1.1 の 1) | `prepare`・`prepare_ns` に初期場のファイルの引数、準備の記録にパス・step・sha256。継承するリミッタの基準値を計画の時点で解決して明示で渡す。初期場の適合性 (輸送種の配置・エネルギー基準・形状と格子の識別) の検査。合格: 計画を確定した後に参照の run に新しい場を足しても使う入力が変わらない (負例)、適合しない初期場で止まる **結果 (2026-10-17)**: `runner_axismach.select_ic` (固定したファイル `ic_file` と計画の sha256 `ic_sha256` の照合、`ic_from` だけのときも選んだファイルを記録)、`prepare`・`prepare_ns` の引数 `ic_file`・`ic_sha256` と準備の記録の `ic_snapshot` (パス・step・sha256)、`design/forge_design/campaign/ic.py` の `resolve_ic` (元の run の種類・継承するリミッタの基準値を明示の値 {ro, p, a, L} に・適合性: 熱力学のモデル・エネルギー基準・輸送種の配置 (`physProp.species`、まとめた種の中身まで)・設計マッハ数)。試験 `design/tests/run_campaign_ic_tests.py` ALL PASS (試走 A: 解決した基準値が run_0005 が実際に使った値と同じ、新しい場を足しても固定したファイルは変わらない、中身が変われば止まる、適合しない 3 種で止まる)。準備の記録の NS 側は `species` を持たないので、配置は元の run の solverConfig から読む。既存の試験 (リミッタの継承・格子の上限・停止規則・本段の B0・ほか 9 本) は通る。`run_wall_temperature_tests` の 4 件の FAIL は変更前の runner でも同じ (手元の変換器の置き場所と固定の入力の化学種の属性; 今回の変更と無関係) | O |
| 4b | ~~実効の計画の確定と、計画からの実行・照合~~ (2026-10-17 実装・試験 ALL PASS) (§4.1.1 の 2・3) | `execution_plan.json` (準備の前・準備の後の確定、上書きしない)、版付きの recipe の段の列と許す差、`run_staged_ns` が計画を受け取り起動の直前に照合、終わった後の照合。合格: 判別の A/B (B は起動 0 回)、段の間の物理の差・計画に無い段で止まる、run_0014 の段の設定から計画を作り直すと一致 (GPU なし) **結果 (2026-10-17)**: `design/forge_design/campaign/plan.py` (`fixed_inputs`・`confirm`・`load_confirmed`・`assert_launch_set`・`assert_launch`・`verify_after`、recipe `axismach.wt_card_chain/v1` の段の列と許す差・禁止する差 (physProp・turbulence・mesh))、`runner_axismach.staged_ns_launch_configs` (run_staged_ns と同じ関数で全段の起動直前の設定、本段は停止規則の上限)、`run_staged_ns(plan=…)` (最初の段の前に全段と境界条件、各段の起動の直前にその段を照合)。確定した計画は sha256 で守り、上書きしない。試験 `design/tests/run_execution_plan_tests.py` ALL PASS (run_0014 の写し、ダミーの forge): **計画の 3 段の設定の FNV が run_0014 の実際の段の記録と一致**、判別の A/B (A 起動 1 回・B (S1 の CFL だけ 0.5 → 0.6) 起動 0 回)、計画の書き換え・二重の確定・前の段の乱流の差・計画の後の境界条件の変更で止まる、走った後の照合。既存の試験 (停止規則 v1・v2・整定・本段の B0・実データの整定 70・リミッタの継承) は通る | O |
| 4c | ~~search の無い campaign と親 plan の分担~~ (2026-10-17 実装・試験 ALL PASS) (§4.1.1) | `forge.campaign/v1` (search なし、未知のキー・search・知らない recipe の拒否)。親 plan §5.1 #6・#7 を #6a/#6b・#7a/#7b に分け、§4.3 に例外を書く。合格: 正例 1・負例 3 以上 **結果 (2026-10-17)**: `design/forge_design/campaign/campaign.py` (`forge.campaign/v1`: 名前・承認の記録 (sha256 の照合)・recipe `axismach.wt_card_chain/v1`・case・台帳。`search`・`evaluations`・`qualification`・知らないキー・知らない recipe と版を拒否)。試験 `design/tests/run_campaign_schema_tests.py` ALL PASS (正例 1・負例 7)。親 plan の §4.3 の例外と #6・#7 の分割は d912df91 で書いた | O |
| 5 | forge-gate と feature/card-local-mesh の統合、`judge_run` の v2 対応、qualification (§4.2) | 合格: run_0011 で RUN_VERDICT が SETTLE_VERDICT と整合、`SETTLED` でも壁解像 FAIL なら総合 FAIL、公差未達でも整定している例・必須の CAD の欠けの例の試験 | O |
| 6a | ~~台帳・停止の記録・予算の勘定~~ (2026-10-17 実装・試験 ALL PASS) (§4.5.1) | 台帳 (flock を生きている間持つ・一時ファイル → rename・campaign の ID・試行の記録・取り込みの一度だけ)、停止の記録 (§4.3 の欄)、予算の予約・消費・精算 (GPU 時間は起動ごとの時刻、ディスクは重ならない所有の範囲)。campaign に `role` を足す。合格: 二つ目の書き手が取れない、取り込みを 2 回しても精算は 1 回、終わりの時刻の無い区間を 0 にしない、二重に数えない、同じ名前の別の campaign で ID が違う **結果 (2026-10-17)**: `design/forge_design/campaign/ledger.py` (flock を持ち続ける唯一の書き手・一時ファイル → fsync → rename・campaign の ID (UUID)・起動の前の試行の記録・試行 ID ごとに一度だけの取り込みと精算・変更の要求のファイル)、`budget.py` (起動ごとの GPU 時間・終わりの時刻の無い起動は根拠の時刻まで・重ならない所有の範囲のディスク・予約の検査)、`stop.py` (停止の記録の欄・上書きしない・driver の補いは根拠が無ければ unknown)、campaign に `role` (fixed / fixer)。試験 `design/tests/run_ledger_tests.py` ALL PASS (13 件)、`run_campaign_schema_tests.py` ALL PASS (負例 8) | O |
| 6b | ~~試行の service・予算の監視・最小の driver と復旧~~ (2026-10-17 実装・試験 ALL PASS) (§4.5.1) | 合格 (ダミーの forge、GPU なし): 起動の要求と台帳の記録の間の障害、driver のいない間の試行の完了、取り込み・精算の途中の障害、二重の `--resume`、停止の途中の再開、CLI が死んでも監視が続く、実行担当が死んだ後に forge だけが残らない、予算の超過で先に記録してから止まる、ログの置き場所の容量不足、既存の出力先に書かない **結果 (2026-10-17)**: `design/forge_design/campaign/step.py` (試行の実行担当: 段の処理を子として回し、予算を監視し、超えたら先に停止の記録を書いてから子のグループを止める)、`driver.py` (最小の driver: 台帳の記録 → systemd の user service の起動 → 起動の印 → 待ち → 一度だけの取り込み、再開の照合、起動の前の検査 (空き容量・出力先・予算)、試験用の障害の仕掛け `FORGE_DRIVER_CRASH_AT`)。試験 `design/tests/run_driver_tests.py` ALL PASS (13 件、実際の user service とダミーの段): T1 ふつう、T2 記録の後に落ちる、T3 起動の後に落ちる、T4 driver のいない間の完了、T5 取り込みの前・後に落ちる、T6 二重の driver、T7 予算の超過 (先に記録、2.6 s で止まる、子が残らない)、T8 driver が死んでも予算の監視が続く、T9 実行担当の SIGKILL (子が残らない、driver が unknown で補う、GPU 時間を 0 にしない)、T10 停止の途中の再開、T11 空き容量、T12 既存の出力先。試験の後に試行の unit は残っていない。実行中のディスクの容量不足 (書き込みの失敗) は試していない | O |
| 6c | ~~recipe の段と分岐・CLI・README の同期・変更の要求~~ (2026-10-17 実装・試験 ALL PASS; GPU の段の本物の実行は V1b・V4) (§4.5.1) | 合格: 判別の A/B (codex: 同じ run の写しで成果物あり → accepted で終わる、成果物なし → `make_deliverables` で成果物の段へ進む; 止まれば第 1 仮説の支持)、p0 の後の最初の較正が較正し直しの回数に入らない、上限で止まる、README の行が出てほかの行が残る、`fixed` で生成した入力の sha256 が使う直前に違えば止まる **結果 (2026-10-17)**: `design/forge_design/campaign/recipes/wt_card_chain.py` (段の列と分岐、各段の処理 `exec` が結果のファイル `result.json` を書き、`next_step` は結果のファイルから位置を決める、前の段が作った入力の sha256 を使う直前に照合 (変更の要求 note があれば通す)、NS の段は Euler の参照と初期場を固定し `deltastar_loop --plan-inputs` で実行の計画を確定)、`deltastar_loop` の `--ic-file`・`--ic-sha256`・`--plan-inputs`、driver の走らせる環境の固定 (初回に写し、再開で照合)・解釈の記録・README の同期 (`campaign/readme.py`)、`campaign/cli.py` (`run` (systemd・Linger の検査、役割 fixed は承認の commit の checkout から)・`status`・`note`)。試験: `run_recipe_tests.py` ALL PASS (分岐 15 件・README 2 件)、`run_recipe_gate_tests.py` ALL PASS (**判別の A/B: 成果物の無い run の判定は make_deliverables で成果物の段へ進み、判定専用の環境で成果物を作った後の最後の判定は accepted で終わる**、MANIFEST の gate は承認の commit)、`run_cli_tests.py` ALL PASS (この PC は Linger=no なので run は止まる)。判定専用の checkout に recipe の格子の台本の場所 (`case/71…/_mesh`) を足した (12.8 MB)。**残り**: 格子・Euler・NS・較正の段の本物の実行 (GPU) は §6 V1b・V4 (#7 の後)、役割 fixed の通し・run の起動の service での driver は GPU の通しで確かめる | O |
| 6d | §6 V1a (§4.5.1) | case/71 の 1 pass の決定的な入力 (設定・境界条件・基準・節点の表・初期場・バイナリ・変換器・数値に効く環境変数) を V6 の run と比べる (許した正規化以外の差 0)、旧の台本 (`p5_ns.py`・`p5_all.sh`) の検査と実行役の関門の対応表 (未対応 0)、登録した負例の通過 0 | O |
| 7 | idle 自動停止と実行器の準備の隙間 (§4.5) | GPU の通しの検証の前提。決め方を諮る | F |
| 8 | 格子の段の契約と 2 つの実装、節点数の上限の承認 (§4.6) | `hand_mesh.py` の固定の承認の文と自動の上限の拡大を、カードの承認済みの上限との照合に置き換える (#10 を待たない)。合格: §6 V2a・V2b・V2c | O |
| 9 | CLI・`--role`・道具の `--help`・skill・手順書・AGENTS.md (§4.4・§4.8) | 手順書 16〜17 行と 277 行、`design-intake` §6・§7、skill `nozzle-design` の書き分け。AGENTS.md の条件 4 の書き分けは settle → gate → この plan の順で、調整板に予告してから | O |
| 13 | §6 V3 依頼の課題集 (GPU なし) | 事前登録 (依頼の文・期待する結果・固定する条件) を commit してから回す。#2・#3・#6・#8 の後 | O |
| 14 | 規約から強制へ: 実行の権限と書き込みの権限を AI から分ける構成 (§4.0) | 例: forge の起動と確定済みのデータへの書き込みを別のユーザ (または別のサービス) に限る。手元と AWS で何ができるかを調べて決め方を諮る。V3 の負例の結果を見てから | F |
| 10 | **格子の係数の一般化と生成器への組み込み** (2026-10-17 ユーザ決定「推奨でいい」・「係数は普通試行錯誤して決めるよね」で `tooling-nozzle-card-local-mesh` から移した; 別のマイルストーン) | 今は case/71 の台本 `case/71.m55_combustion_wt/_mesh/hand_mesh.py` が `wt_card` の問題 YAML の格子の節を書き換える (M5.5 系の 1 周目の係数 配管 1.0・縮流部 0.6・スロート 0.5・拡大部 0.49)。(1) 生成器に組み込む: `hand_mesh.py` の中身を `design/forge_design` に移し、雛形の `mesh.method: local_v1` (係数は `mesh.first_layer_calibration` の表) を既定にし、手順書の §3 を不要にする。(2) **係数の試行錯誤のループ**: 今の係数を初めの値として格子 → NS (停止規則 v2) → `check_wall_resolution.py` (壁ごとの y₁⁺ > 1 の面積 ≤ 2 %・超えた位置) → 外れた区間だけ係数を下げて格子と NS をやり直す、回数の上限 (例 2 回) を超えたら人に返す。直した係数と根拠 (どの run のどの区間) を記録し、機種の系列で継承する。(3) 格子の感度の A/B (同じ係数の網を細かくして目的の量 (出口コア M・V5c・δ_E) が許容の内か) を事前登録してから。合格: 寸法が大きく違うノズル 1 本 (#11 と同じ run を使ってよい) でループが止まる・人に返る、M5.5 系は今の格子とビット同一 (組み込みの退行なし) | F |
| 12 | **NS の生の排除厚さを壁に直接入れる経路を生産から外す** (`tooling-nozzle-card-local-mesh` §5.1 #3v から移した、2026-10-17; 別のマイルストーン) | ユーザ 2026-10-15「ns計算で、生の排除厚さ入れるやり方、もう取る意味ないよね」。要求カードの流れの既定を「較正した CONTUR (係数を NS の δ_E から決める; `contur_calibrate`)」にし、`wt_card --ns-correction` (直接の表) は既定から外して使うと警告。コードの削除は、直接の経路を使う他のセッションの plan (`tooling-nozzle-isothermal-wall-chain` §4.9 の受け渡しの関門・③ベル など) と範囲を確かめてから (相談する)。この plan の方針の改訂を理由には削除しない (codex 2026-10-17) | F |
| 11 | **収束の判定 (停止規則 v2・プリセット `wt_nozzle_v3`) の別の寸法のノズルでの確認** (2026-10-16 に plan `tooling-settle-convergence-and-monitor` の完了条件から移した; ユーザ決定「後続の plan へ移す」; 別のマイルストーン) | 寸法が大きく違うノズル (小さいノズル) を 1 本、新しい出力を有効にして停止規則 v2 で回し、オンラインと再生の判定の一致・停止の時の出力の揃い・停止の後の区間 (V3(c) と同じ形の自由継続の再生) を確かめる。プリセットの閾値は暫定のまま使い、結果を見て較正しない。合格の形は settle plan §6 V3(b)・V3(c)・V4 に倣う | F |

## 6. 検証 (事前登録は各項目の実装の前に確定する)

- **前提**: settle plan の V6 (`run_0014_chain_v2`) が通っていること (accepted)。
- **V1a (GPU 不要; 2026-10-16 codex M4)**: 次の 2 つを分けて確かめる。
  - **決定的な入力の比較**: 実行器が case/71 の 1 pass のために作る起動の入力 (設定・境界条件・基準・節点の表・`RECIPE.json`・初期場・バイナリ・変換器・数値に効く環境変数) を、V6 の run のものと比べる。
    - 生成時刻・配置のパスは、項目を限って正規化する。
    - 前の段の成果物は V6 のものを固定して使う。
  - **来歴の照合**: 参照先の内容とハッシュの結び付けを別に検査する。
  - あわせて、台本 (`p5_ns.py`・`p5_all.sh`) が行っていた検査がすべて実行器の関門として存在し、負例で止まること。
- **V1a の事前登録 (2026-10-17、§5.1 #6d、結果を見る前)**:
  - 対象: case/71 の 1 pass = V6 (`case/71.m55_combustion_wt/run_0014_chain_v2/`、台本 `_runs/run_0014_chain_v2/p6_chain.sh`) と同じ入力
    (`problem_m55_p5_settle_v2.yaml`、Euler の参照 `run_0009_euler/res_12000.h5` sha256 2b2d02ce…、初期場 `run_0008_ns_b0/res_90000.h5` sha256 c64a73f1… (AWS B から写した))。
  - 実行役の側: recipe の NS の段と同じ CLI (`deltastar_loop --init-integral --ic-file … --ic-sha256 … --euler-ref … --euler-ref-sha256 …`) を
    `--prepare-only` で一時の run に準備し、`campaign.plan.confirm` で実行の計画を確定する (forge は起動しない、手元)。
  - 比べるものと許す正規化: (1) 全段の起動直前の設定: 計画の 3 段の FNV = run_0014 の `stage_manifest.json` の `cfg_fnv` (正規化なし)。
    (2) `bcondConfig.yaml`・(3) `solverConfig_prestage.yaml`: バイト一致。(4) `settle_criteria.yaml`: 中身が一致、ただし `source.generator_sha256` (生成器のコードの版) は除く。
    (5) `settle_nodes.txt`: コメントの行を除いて一致。(6) 格子: `nozzle.h5` の座標とトポロジのハッシュ・`nozzle.msh` の sha256 が一致。
    (7) 準備の記録の主な欄 (`initializer`・`Md`・`x_E`・`scale_m`・`physical_wall`・`stop_rule` の cap・`main_line`・継承したリミッタの基準値) が一致 (パス・時刻を除く)。
    (8) 初期場: 実行役は `ic_snapshot` (res_90000.h5 の sha256) を固定して記録する。run_0014 は `ic_from` (同じディレクトリの最大番号の場)。
    (9) バイナリ・変換器: 実行役が sha256 を計画に書く。値の一致は求めない (run_0014 は AWS B の FP64 のビルド、手元は別のビルド; 手元で回す V1a の制約)。
  - 合格: (1)〜(7) で許した正規化以外の差 0、(8)・(9) は記録の有無。旧の台本 (`p5_ns.py`・`p5_all.sh`・`p6_chain.sh`) の検査と実行役の関門の対応表で未対応 0。
    登録した負例 (#2〜#6 の試験) の通過 0。
- **V1a の結果 (2026-10-17)**: **FAIL** (事前登録の文言どおり)。(1) 全段の設定・(2) 境界条件・(3) 前の段の設定・(8) 初期場の固定・(9) 記録は一致・あり。
  (4) 停止の基準は `source.written` (書いた時刻; 事前登録の正規化に書き漏らした) が違い、(5) 節点の表 (3,830 行 vs 3,831 行)・(6) 格子 (座標と msh が違う、トポロジは一致)・
  (7) 初期壁の数値 (x_E 40.633257 vs 40.630797) が違う。記録 `case/71.m55_combustion_wt/_runs/v1a_executor/`。codex (diagnose) に諮った
  (`notes/reviews/2026-10-17-orchestrator-v1a-fail-diagnose.md`): FAIL は残し、呼び出し経路の差を別の検証 V1a′ で切り分ける。比較の台本 `compare.py` は部分文字列で除外していて
  緩い (計測時間・配置のパスを除いた、停止規則を事前登録より広く比べた) ことも指摘された (記録として残し、V1a′ では使わない)。
- **V1a′ の事前登録 (2026-10-17、結果を見る前; codex 諮問 v1a-fail の判別の A/B)**:
  - 因子は呼び出し経路だけ。A = 旧の台本 `p6_chain.sh` と同じ呼び方 (`run_pass0_integral(prob, eref, run, ic_from=<run_0008 の写しのディレクトリ>, stop_rule=True, prepare_only=True)`)、
    B = 実行役の CLI (`deltastar_loop --init-integral --ic-file … --ic-sha256 … --prepare-only`)。同じコード (commit)・Python と依存・変換器・数値に効く環境変数・問題 YAML・
    Euler の参照・初期場 (写しのディレクトリには res_90000.h5 だけ置き、A が選ぶ最大番号の場が sha256 c64a73f1… であることを確かめる)。各 1 回準備 (forge 0 step)。
  - **揺らぎの欄** (A と B を比べる前に決める): B を 2 回 (B1・B2) 準備し、B1 と B2 で違うキー (完全なキーパス) を記録する。A と B1 の比較では、この欄だけを除く。
  - **比べるもの**: 全段の起動直前の設定 (`staged_ns_launch_configs` の全文)、`bcondConfig.yaml`・`solverConfig_prestage.yaml`・`solverConfig.yaml` (継承したリミッタの基準値を含む全キー)、
    `settle_criteria.yaml` (全キー)、`settle_nodes.txt` (コメントの行を除く)、`nozzle.h5` の座標・トポロジと**補間した後の初期場** (VALUE の全データセットの値)、`nozzle.msh`、
    `prepare_info.json` の全キー。
  - **別に記すもの** (経路の違いそのもの; 合否に入れない): `prepare_info` の `stop_rule.explicit`・`stop_rule.requested_source` (A は明示、B は既定)、`ic_snapshot.how`・`ic_from`。
  - **合格**: 比べるもので、揺らぎの欄と別に記すもの以外の差 0 → この入力について呼び出し経路の差を除外する (V1a とは別の検証として記録)。差が残れば、版の差だけの説明を棄却し、
    最初に差が出た生成の段を調べる。**解釈の限定**: A = B でも、V6 との差の唯一の原因が #6e とは確定しない (環境・変換器も違う)。A を「V6 の再現」と呼ばない。
- **V1b (GPU 1 pass、記録)**: 次の種類が V6 と同じであること。
  - REACHED・関門の pass・`delta_E.status ok`・成果物の実在・非有限 0。
  - qualification の判定は `accepted` (V6 の `case/71.m55_combustion_wt/run_0014_chain_v2/` は 3 係数で補正した壁で、保存した `WT_QUALIFICATION.json` は accepted:
    出口コア M 5.502842、V5c の最大の偏差 0.0623 %・片振幅 0.0573 %; 2026-10-17 codex plan 段 Major 3 で前の記述「補正へ戻す」を直した)。成果物までそろえる。
  - 「補正へ戻す」の正例は別の fixture に分け、**試走 A の `case/90.qs_trial_a/run_0003_ns_p0/` に固定する** (2026-10-17 codex 諮問)。
    解析は成立している (`check_settled.py` 32,000〜42,000・停止の観測 42,060 で `VERDICT: SETTLED`、メッシュ品質 PASS、壁解像 PASS、非有限 0) が、
    出口コア M 5.519563 (要求 5.5 から +0.356 %)・V5c の最大の偏差 0.404 % で、公差 0.1 % では V5b と V5c が FAIL する。
    固定するもの: `card_p0.yaml`・`res_42060.h5` (sha256 を事前登録に書く; 場は git 追跡外なので置き場所も書く)・公差・評価窓・反復の上限。
    期待の結果は `back_to_correction`。次の手は較正の規則しだいで、規則の `apply=true` が無ければ `escalate` (`wt_qualify.py:877`・`:908`;
    「補正へ戻す」と「自動の較正を続ける (`recalibrate`)」は今のコードでは別の判定)。
    case/71 の 1 回目の壁は未整定・ディスク満杯・延長も CENSORED (case/71 README) なので、公差未達の正例に使わない。同じ参照に相反する期待の結果を置かない。
  - 段の間の引き継ぎと完了の記録。
  - 数値の差は表にして記録する。数値の許容は、同じホスト・同じ入力の反復 2 本で床を取るまで事前登録しない。
- **V2a**: 手組みの格子を包んだ段が、case/71 の節点の座標と一致する。
- **V2b**: 雛形の格子で、別のカード (case/70 M5 の空気) が品質の関門を通る。
- **V2c (節点数の上限の負例)**: 別のカード + 過去の承認の文 (case/71 の「ユーザ決定『a』2026-10-15」) で上限を広げようとすると、格子を作る前に止まる。承認した上限を超える作り直しも止まる。case/71 のカードと承認では V2a と同じ格子ができる。
- **V3 依頼の課題集 (GPU なし; 2026-10-17 codex で前の案の「問いの組で判断役の成績」を置き換えた)**:
  - **依頼**: 正常 3 件 (既知の機種で条件違い。例: case/71 のカード・試走 A のカード・case/70 M5 の空気) と、止まるべき 3 件 (公差が届かない依頼・ユーザ欄が欠けた依頼・予算が足りない依頼)。
    依頼の文・期待する結果・固定する条件は、回す前に commit する。
  - **GPU なしの到達点**はカード → 問題 YAML → 格子 (品質の VERDICT) まで。これを設計全体の `accepted` に数えない。流れの公差の未達は、固定した既存の不合格のデータ
    (例: case/71 の 1 回目の壁の NS) の再判定で扱う。
  - **安価な負例** (AI を使わず、プログラムの検査だけで): 根拠の無いユーザ欄、出所の未知のパス・既定との不一致・予算の欠落、`requirements_without_gate` の削除、CLI での公差・評価窓・必須成果物の指定、承認の後の `line_n`・定常性の窓と閾値・物性の変更、物理の設定の不一致、判定器なし、
    外部の Python・import の元の違い、総合 PASS と個別の FAIL・欠落の矛盾、予算の実行中の超過、再開による二重の起動、既存の出力先への書き込み。全部が止まる・拒否される・`incomplete` になること。
  - **正例を並べる** (全部止める実装を合格にしない): 承認したカード、格子・較正の係数だけを変えた試行の入力は止まらずに通る。
  - **回し方**: 賢い AI (`--role fixer`) と知識の少ないモデル (Haiku、使えればローカル LLM; `--role fixed`)。モデルの版・役割 (権限)・入力・予算・人の介入の条件を固定し、
    各依頼を各モデルで 3 回程度繰り返す。
  - **数えるもの**: 到達率、正しく返した率、誤合格の件数、費用。分母と途中の停止も記録する。**誤合格は機械の結果と AI の報告で別に数える**
    (JSON が `incomplete` なのに AI が「要求を達成した」と報告したら、報告の側の誤合格)。
  - 合格: 機械の結果の誤合格 0、安価な負例が全部止まる。AI の報告の誤合格は件数を記録する (0 でなければ skill・手順書に戻して直す)。
    観測された誤合格 0 を一般の保証とは扱わない。
- **V4 GPU の通し (#7 の後; 2026-10-17 codex plan 段 Major 4 で成功の条件とモデルの測定を分けた)**:
  - **(a) 決定的な recipe の通し (成功の条件)**: 既知の正常のカード 1 件 (例: case/71 のカード) について、`forge-design run` が入口から、カードの生成・1 回目の NS・較正・次の壁への引き継ぎ・
    最終の判定を通って、必須の成果物と最終の `accepted` まで届く。人もモデルも途中で直さない。**未較正のカード (係数なし) と入力の参照を固定し、既存の補正済みの壁を途中で入れて
    合格する余地をなくす** (2026-10-17 codex 諮問)。
  - **(b) モデルの測定 (成功の条件と分ける)**: 同じ依頼を知識の少ないモデル (`--role fixed`) と賢い AI (`--role fixer`) で回し、止まった所・直したこと (fixer)・人に返した所・
    カードのどの欄を自分で決めたか (出所の記録)・消費した予算を記録する。両モデルの全試行の成功は求めない。
- 判定不能は合格にしない。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| diagnose (Fable、ユーザ指定) | 2026-10-16 | diagnostician の返答 (ブリーフ [briefs/2026-10-16-design-orchestrator.md](../../notes/reviews/briefs/2026-10-16-design-orchestrator.md); 本文は §4 に統合) | C1/M10/m3 | 全件**採用**: 親 plan の実行器 + ledger と forge-gate の判定・関門の上に載せる (二重化をやめる)、判定の 1 本化、返り値は既存の語彙の和 + 原因の分類、答えは引用で照合 (確信は記録だけ)、digest、問いの組の作り方、実行の場所、格子は包むだけ、台帳の巨大な出力、AGENTS.md の範囲と順序、カードの欄の表、集計。**2026-10-17 の改訂で、digest・引用の照合・問いの組は「今は作らない」(§4.9) に移した** |
| plan | 2026-10-16 | [2026-10-16-tooling-design-orchestrator-plan.md](../../notes/reviews/2026-10-16-tooling-design-orchestrator-plan.md) | GO-with-changes, C0/M5/m2 | 全件**採用** (§2・§4・§5・§6・§7 を書き直し): M1 親 plan の未実装の契約 (最小の campaign・入力の固定・execution_plan・設定の照合) を実行器より先に (§4.1、§5.1 #4)。M2 計算の成立と設計の要求の達成を分ける qualification (§4.2、#5)。M3 再開の契約 (ロック・原子的な更新・attempt ID・プロセスとの結び付け・新しい run での再試行・障害の注入の試験; idle 自動停止を GPU の検証の前提へ) (§4.5、#6・#7)。M4 V1a を決定的な入力の比較と来歴の照合に分ける (§6)。M5 カードの値の出所 (§4.7、#2)。m6 行動の選択の呼び出しを外部レビュー用と分ける (2026-10-17 の改訂で、行動の選択の呼び出しそのものを §4.9 へ)。m7 旧い記述を残さず書き直した |
| diagnose (codex) | 2026-10-17 | [2026-10-17-orchestrator-direction-diagnose.md](../../notes/reviews/2026-10-17-orchestrator-direction-diagnose.md) (ブリーフ [briefs/2026-10-17-orchestrator-direction.md](../../notes/reviews/briefs/2026-10-17-orchestrator-direction.md)) | 方向は採用、Major 多数 (採否表 6 件 + 本文 6 件) | 全件**採用** (主セッションが根拠の行 `wt_qualify.py:54`・`:179`・`run_case.sh:52`・`:55`・`hook_forge_guard.py:12`・`design-chain-quickstart.md:277`・`design-intake` §7 を確かめた): A1 判定器の保護は sha の記録でなく固定した版の別 checkout + 判定条件の固定、CLI・問題 YAML の相違は拒否 (§4.0・§4.2、#2)。A2 予算は起動前の予約・実行中の監視と停止・終了時の精算、失敗・計時 run・前段も計上 (§4.5、#6)。A3 要求の一致の照合を最優先、導出の関係を扱う (§4.2、#3)。**A4 (accepted の chmod) は却下に従う** → 試行ごとの新しい出力先と書き込みの権限の分離、chmod は補助 (§4.5)。A5 §4.4 を停止の記録にし §4.3 も直す (§4.3・§4.9)。A6 V3 は初期の試走として扱い、安価な負例・誤合格の 2 系統の数え方・3 回の反復を足す (§6 V3)。要求の一致・解析の成立・要求の達成の 3 分け (§4.2)。数値の設定は値を照合しないが妥当性の検査は残す (§4.2)。判定器の無い要求は `incomplete` (§4.2)。条件 4 は範囲を限った書き分け (§4.4、#9)。`run_case.sh` は強制された唯一の入口とみなさない・規約と強制を区別 (§4.0、#14)。§5.1 の優先順を codex の順に (§5)。#10・#12 は削除しない |
| plan | 2026-10-17 | [2026-10-17-tooling-design-orchestrator-plan-2.md](../../notes/reviews/2026-10-17-tooling-design-orchestrator-plan-2.md) | GO-with-changes, C0/M4/m1 | 全件**採用** (採否は下の diagnose で諮った): M1 実行担当の障害 → 試行ごとの実行担当 + 停止の責任を外の service に (§4.5、#15・#6)。M2 格子の節点数の承認の流用 → 依頼ごとの承認済みの上限と承認元、Euler・NS とも照合 (§4.6・§4.7、#8、V2c)。M3 V1b の期待 → `accepted` + 成果物、「補正へ戻す」は試走 A の `run_0003_ns_p0` に固定 (§6 V1b)。M4 完了条件 → V4 (a) 決定的な recipe の通しを必須に (§6・§8)。m5 case の run 一覧の同期 (§4.8、#6) |
| diagnose (codex) | 2026-10-17 | [2026-10-17-orchestrator-plan2-disposition-diagnose.md](../../notes/reviews/2026-10-17-orchestrator-plan2-disposition-diagnose.md) (ブリーフ [briefs/2026-10-17-orchestrator-plan2-disposition.md](../../notes/reviews/briefs/2026-10-17-orchestrator-plan2-disposition.md)) | 5 件とも採用、M1 は反映が未完了 | 全件**採用**: M1 停止の責任を実行担当の外 (試行ごとの systemd の service、`ExitType=main`・`KillMode=control-group`・停止期限・`Restart=no`) に置き、全子孫の終了を確かめるまで新しい試行と予算の解放をしない、service が使えなければ起動を止める、GPU なしの `SIGKILL` の A/B を実装の先頭に (§4.5・§5・#15)。M2 「そのまま包む」→ 数値の処理を再利用し承認と上限の部分を置き換え、Euler・NS とも照合 (§4.6)。M3 fixture を試走 A の `run_0003_ns_p0` に固定、期待は `back_to_correction` (次の手は較正の規則しだい) (§6 V1b)。M4 未較正のカードと入力の参照を固定 (§6 V4 (a))。m5 破棄のときも同期・再開で台帳から修復・case ごとのロック (§4.8) |
| diagnose (codex) | 2026-10-17 | [2026-10-17-orchestrator-approval-design-diagnose.md](../../notes/reviews/2026-10-17-orchestrator-approval-design-diagnose.md) (ブリーフ [briefs/2026-10-17-orchestrator-approval-design.md](../../notes/reviews/briefs/2026-10-17-orchestrator-approval-design.md)) | 方向は採用、§4.7.1 は修正してから実装 (Major 4・Minor 2) | 全件**採用** (主セッションが根拠の行 `wt_qualify.py:203`・`:773`・`:898`・`deliverables.py:48`・手順書 §7.2 を確かめた): 承認の対象を生成した後の実効の条件 (評価の点数・定常性の条件・物性・実行の方針を含む) で決める、**カード全体の sha256 の一致を試行の入力の必須条件にする案は却下に従う** → 承認の原本と試行の入力を分ける、判定専用の環境 (外部の道具の実体・環境変数を引き継がない・import の元・未追跡のファイル)、`pip freeze` は補助の記録だけ、合否に使う中間の判定を全部再検査、`sources` の正規化と既定の照合、preview と write を分ける (§4.7.1、§5.1 #2a〜#2c、§6 V3) |
| diagnose (codex) | 2026-10-17 | [2026-10-17-orchestrator-minimal-campaign-diagnose.md](../../notes/reviews/2026-10-17-orchestrator-minimal-campaign-diagnose.md) (ブリーフ [briefs/2026-10-17-orchestrator-minimal-campaign.md](../../notes/reviews/briefs/2026-10-17-orchestrator-minimal-campaign.md)) | C0/M3/m1 | 全件**採用** (主セッションが根拠の行 `runner_axismach.py:2946`・`:3002`・`:2676`・`:2523`・`:1777` を確かめた): **(b)・(d) の元の案は却下に従う** → 準備の前の固定・準備の後の確定・計画からの実行と前後の照合の 3 段、版付きの recipe が段の列と許す差を持つ、段の記録の key を全設定の証拠にしない (§4.1.1、#4b)。(c) は採用し、継承する値の解決と初期場の適合性を足す (#4a)。親 plan の #6・#7 を限定の実装と探索待ちに分ける (#4c) |
| diagnose (codex) | 2026-10-17 | [2026-10-17-orchestrator-executor-diagnose.md](../../notes/reviews/2026-10-17-orchestrator-executor-diagnose.md) (ブリーフ [briefs/2026-10-17-orchestrator-executor.md](../../notes/reviews/briefs/2026-10-17-orchestrator-executor.md)) | C0/M5/m1 | 全件**採用** (主セッションが根拠の行 手順書 §6・§8 の表 (667 行)・`gate_env.py:326` を確かめた): M1 recipe の分岐を手順書に合わせる (`make_deliverables` → 成果物、p0 の後の最初の較正は判定と別)、M2 driver の復旧 (先に記録してから起動・取り込みは一度・ロックを持ち続ける・二重の再開・ID)、M3 予算で止めるときは先に記録、M4 予算の測る範囲 (重ならない所有の範囲・終わりの無い区間を 0 にしない)、**M5 作り直して比べる `fixed` の照合は却下に従う** → 承認の commit の checkout から回して生成器を固定し、生成した入力の sha256 を使う直前に照合、m1 ログと変更の要求。固定の倍率 (1.5 倍) で止める案は入れない。手元の GPU で V4 (a) を前倒ししない (§4.5.1、§5.1 #6a〜#6d) |

## 7. 影響範囲

- 新しい道具を足す: campaign の最小の実装・実行の入口 (予算・試行 ID・新しい出力先)・実行器・台帳・停止の記録・要求の一致の照合・CLI。
- 変えるもの:
  - 要求カードの書式 (`forge.wt_card` の版を上げる: `qualification`・予算・承認の版)。
  - `wt_qualify.py` (正式な入口での相違の拒否、要求の一致の照合、判定器の版の固定、中間の判定の結び付け)。
  - forge-gate の `judge_run.py` (v2 の run の出典)。
  - AGENTS.md の 2 項目 (レシピのある機種の設計の run は `forge-design` 経由、エスカレーション条件 4 の書き分け)。
  - 手順書 `design-chain-quickstart.md`、skill `design-intake`・`nozzle-design`、各道具の `--help`。
- 既存の部品 (runner・停止規則・関門・成果物) は呼び出すだけにする。振る舞いを変えるときは、それぞれの plan に書く。

## 8. 完了条件

V1a・V1b・V2a・V2b・V2c・V3・V4 (a)。V4 (a) は既知の正常のカード 1 件が要求カードから必須の成果物と最終の `accepted` まで届くこと (正しく止まるだけでは完了にしない)。V4 (b) は「詰まる所の一覧と、その対処」で完了とする (全部が通ることは求めない)。§5.1 #10・#11・#12 は別のマイルストーンで、完了条件に入れない。

## 9. 変更ログ

- `2026-10-17` — §5.1 #6c を実装: recipe `axismach.wt_card_chain/v1` (`campaign/recipes/wt_card_chain.py`)、CLI、README の同期、走らせる環境の固定、`deltastar_loop` の計画のつなぎ。判別の A/B (成果物なし → 成果物の段へ、成果物の後 → accepted で終わる) を本物の判定で確かめた。
- `2026-10-17` — §5.1 #6b を実装: 試行の実行担当 (`campaign/step.py`) と最小の driver (`campaign/driver.py`)。障害の試験 `run_driver_tests.py` ALL PASS (13 件、実際の user service)。
- `2026-10-17` — §5.1 #6a を実装: 台帳・予算の勘定・停止の記録 (`campaign/ledger.py`・`budget.py`・`stop.py`)、campaign の `role`。試験 `run_ledger_tests.py` ALL PASS。
- `2026-10-17` — §5.1 #6 の設計 (§4.5.1) を codex (diagnose) に諮り、全件採用して書き直した: recipe の分岐を手順書に合わせる (成果物の前で止めない)、driver の復旧、予算の停止の順序と測る範囲、`fixed` は承認の commit から回す。#6 を #6a〜#6d に分けた。
- `2026-10-17` — §5.1 #4c を実装: search の無い campaign (`campaign/campaign.py`)。試験 `run_campaign_schema_tests.py` ALL PASS。#4 (#4a〜#4c) はこれで終わり。
- `2026-10-17` — §5.1 #4b を実装: 実行の計画 (`campaign/plan.py`) と `run_staged_ns(plan=…)`。run_0014 の写しで、確定した計画の段の設定が実際の段の記録と一致、判別の A/B (B は起動 0 回) を確かめた。試験 `run_execution_plan_tests.py` ALL PASS。
- `2026-10-17` — §5.1 #4a を実装: 初期場をファイルで固定 (`select_ic`・`ic_file`・`ic_sha256`・準備の記録の `ic_snapshot`)、計画の時点の解決 (`campaign/ic.py` の `resolve_ic`: 継承するリミッタの基準値の明示化・適合性)。試験 `run_campaign_ic_tests.py` ALL PASS。
- `2026-10-17` — §5.1 #4 の設計 (§4.1.1) を codex (diagnose) に諮り、全件採用して書き直した: 計画の 3 段 (準備の前の固定・準備の後の確定・計画からの実行と照合)、継承する値の解決と初期場の適合性、親 plan の #6・#7 の分割。#4 を #4a・#4b・#4c に分けた。
- `2026-10-17` — §5.1 #3 を実装: run の本段の物理の設定と承認の照合 (`approval.check_run_physics`)。正式な入口の通しの正例 (判定専用の環境で成果物 → `wt_qualify --approval` → accepted) を case/71 run_0014 の写しで確かめた。FreeCAD が checkout に書く `__pycache__` を照合の前に消す (3d840047)。
- `2026-10-17` — §5.1 #2c を実装: 成果物の部品の `pass`・総合との矛盾・必須のファイル (両方のモード)、承認のモードで停止の判定と格子の品質の回し直し、成果物の gate の照合、`gate_env deliverables`。判定専用の環境で run_0014 の写しから成果物一式を作れた (FreeCAD・gmsh・変換器は承認のときの記録だけ)。#2 (#2a〜#2c) はこれで終わり。
- `2026-10-17` — §5.1 #2b を実装: 判定専用の環境 (`forge_design/feedback/gate_env.py`)。試験 `design/tests/run_gate_env_tests.py` ALL PASS (正例 3・負例 8)。最初の版は sparse checkout に物性の DB (`solver_density_cuda/data`) が無く、試験の条件が甘くて通っていたのを直した (c4f69474)。
- `2026-10-17` — §5.1 #2a を実装: 承認の記録 (`forge_design/intake/approval.py`)、要求カードの版 2 (`forge.wt_card/v2`)、判定の正式な入口 `wt_qualify --approval`。試験 `design/tests/run_approval_tests.py` ALL PASS (正例 4・負例 15・判別 A/B 2)。現在仕様は `methods/design/overview.md`「カードの版 2 と承認の記録」。
- `2026-10-17` — §5.1 #2 の設計 (§4.7.1) を codex (diagnose) に諮り、全件採用して書き直した: 承認の対象は生成した後の実効の条件、承認の原本と試行の入力を分ける、判定専用の環境、中間の判定の全部の再検査。#2 を #2a (保護する条件と承認の記録)・#2b (判定専用の環境)・#2c (中間の判定の再検査) に分け、§6 V3 に負例と正例を足した。
- `2026-10-17` — §5.1 #15 を実施: 試行ごとの systemd の user service (`ExitType=main`・`KillMode=control-group`・`TimeoutStopSec=5`・`Restart=no`) で、実行担当への `SIGKILL` の後に全子孫が 0.10 秒で消えた (手元 WSL2・AWS B とも A PASS・B PASS、対照は子孫が残った)。`Linger=no` では ssh を切ると約 12 秒で service が止まるので、実行役の事前検査に `Linger=yes` を足した (§4.5)。台本 `design/tests/run_attempt_service_ab.py`、記録 `notes/investigations/2026-10-17-attempt-service-ab/`。
- `2026-10-17` — codex の plan 段レビュー (2 回目、GO-with-changes、C0/M4/m1) を、codex (diagnose) に採否を諮って全件採用した (§6.1)。§4.5 に試行ごとの実行担当と、停止の責任を外の service に置く契約、§4.6 に依頼ごとの節点数の上限の承認、§4.8 に case の run 一覧の同期、§6 V1b の期待を `accepted` に直し「補正へ戻す」を試走 A の `run_0003_ns_p0` に固定、V2c を足し、V4 を (a) 決定的な通し (必須) と (b) モデルの測定に分けた。§5.1 に #15 (service の事前検査と障害の A/B、実装の先頭)。
- `2026-10-17` — codex (diagnose) の諮問 ([記録](../../notes/reviews/2026-10-17-orchestrator-direction-diagnose.md)) を全件採用し、§1〜§8 を書き直した (前の案は commit 4039c301)。守りを「承認した要求と判定条件の固定」「予算と既存データの実行時の保護」「出口の判定」の 3 つにし、途中の修正はその範囲で AI に任せる (§4.0)。判定を要求の一致・解析の成立・要求の達成に分けた (§4.2)。判断役の口・問いの組・`local` の backend・`forge-design answer` は今は作らない (§4.9)。§4.4 に途中の修正の範囲と AGENTS.md の条件 4 の書き分けの文案、§4.5 に予算の予約・監視・精算、§6 V3 を依頼の課題集に置き換えた。§5.1 は #2〜#9 を新しい中身にし、#13 (V3)・#14 (規約から強制へ) を足した。#10〜#12 は番号を保ち別のマイルストーンとした。
- `2026-10-17` — ユーザ決定: 設計方針の改訂案の方向で進める・予算は依頼ごとに見積もって承認する (§4 の冒頭)。本文の書き直しは codex 諮問の後。
- `2026-10-17` — ユーザ決定: 依頼の原文は設計を行うリポジトリの中 (case のディレクトリ) に置く (forge には載せてよい情報だけ; §4.7)。
- `2026-10-16` — codex の plan 段レビュー (GO-with-changes、C0/M5/m2) を全件採用し、§2・§4・§5・§5.1・§6・§7 を採用済みの方針に直接書き直した (旧い記述を残さない)。
- `2026-10-16` — diagnostician (Fable) の諮問を全件採用: 親 plan (campaign-recipe) の §5 ステップ 5 の実装に位置づけ直し、forge-gate の判定・関門の上に載せる、判断役の答えは引用で縛る、格子は包むだけ、V1 を V1a (入力の一致) と V1b (種類の一致) に分けた (§4.7)。

- `2026-10-16` — 起票 (draft)。ユーザ「そのやり方でよいか、diagnosticianにも相談した上で、進めてほしい」。
