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
  - `plans/active/tooling-settle-convergence-and-monitor.md`: 判定・停止規則 v2・受け渡しの関門 #7。この plan の前提。
  - `plans/active/tooling-nozzle-card-local-mesh.md`: case/71 の要求カード・手組みの格子・3 係数の壁。
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

この plan の目的は次の 2 つ。

- 要求カードから成果物までを 1 つの入口 (CLI) で通す。
  - 段を回すのは決定的なプログラム (実行役)。
  - AI が要るのは入口 (カードを作る)・出口 (報告)・実行役が止まったときの判断だけにする。
- 機種の違いは、データ (レシピ) とそれが呼ぶ段の部品で足す。判断役は Fable・codex・Claude・ローカル LLM のどれでもよいようにする。

## 2. スコープ

- **やる (最初の範囲)**: 風洞ノズル (axismach) の 1 機種で、次の 4 つを作る。
  - 実行役
  - 格子の段の一般化
  - 判断役の口 (形の決まった問いと答え、判断役の切り替え)
  - 失敗の返り値ごとの次の手の表
- **やる**: case/71 の今の台本を段の部品に置き換え、同じ結果が出ることを確かめる。
  - 置き換える台本: `_runs/run_*/p*_ns.py`・`p*_all.sh`・`judge_*.sh`
- **やる**: 知識のない AI での再現 (Haiku・ローカル LLM) で、詰まる所を潰す。
- **やらない (別の plan・後の範囲)**:
  - ベル・デュアルベル・SERN のレシピ (各セッションの成果を載せる形で後から)。
  - 判断役そのものの能力の改善。
  - 機種ごとの物理の公差を決めること (ユーザと決める)。

## 3. 関連 docs と前提

- 前提: settle plan の #7 (受け渡しの関門の新方式) と V6 (設計チェーンの通しの 1 pass) が通っていること。この plan は、その流れを部品にして束ねる。
- 今ある部品:

| 部品 | 役割 |
| --- | --- |
| `design/forge_design/intake/wt_card.py` | 要求カード → 問題 YAML |
| `design/forge_design/intake/recipes/axismach_*.yaml` | 雛形 |
| `case/71.m55_combustion_wt/_mesh/hand_mesh.py` | 手組みの格子。case/71 専用 |
| `design/forge_design/evaluate/runner_axismach.py` | Euler・NS の準備と起動 |
| `design/forge_design/evaluate/stop_rule.py` | 停止規則 v2 |
| `design/forge_design/feedback/deltastar_loop.py` | 境界層の補正の繰り返しと受け渡しの関門 |
| `design/forge_design/export/deliverables.py` | 成果物 |
| `design/forge_design/report/nozzle_report.py` | 設計報告 |
| `solver_density_cuda/tools/codex_review.py` | 判断役の切り替え (codex / claude) |
| `solver_density_cuda/tools/aws_instance.sh` | AWS の起動と状態 |

## 4. 設計方針

### 4.1 3 つの層

1. **要求カード (データ)**: 何を作るかを書く。
   - 機種、運転条件、要求と許容 (例: 試験部の一様性 V5b・V5c の公差)、欲しい成果物、使ってよい計算資源 (手元か AWS か、上限の時間)。
   - 今の `wt_card` (forge.wt_card/v1) を、機種の欄を持つ形に広げる。
2. **機種ごとのレシピ (データ、版つき)**: どう作るかを書く。
   - 段の並び、各段の既定の設定 (推奨設定の正本を参照)、判定のプリセット。
   - 各段の失敗の返り値ごとの次の手 (選べる手の一覧とその前提)。
   - 設計を終える条件と、繰り返しの上限。
3. **段の部品と実行役 (コード)**: 段の部品は関数で、契約を持つ。
   - 契約の中身: 入力、出力、合格の判定 (判定ツールの VERDICT)、失敗の返り値 (列挙)。
   - 実行役は決定的な状態機械。カードとレシピを読み、段を順に回し、台帳に記録して、止まっても続きから再開できる。

### 4.2 段 (風洞ノズル)

| # | 段 | 内容 |
| --- | --- | --- |
| 1 | カード → 寸法の解決 | M_knot・r_t。今の `wt_card` |
| 2 | 格子 | case/71 の `hand_mesh.py` を一般化した格子の段。品質の関門と AR の緩和は今の道具 |
| 3 | Euler の参照 | |
| 4 | NS | 本段 B0・停止規則 v2・`wt_nozzle_v2` |
| 5 | 受け渡しの関門 | settle plan §4.9 の新方式 |
| 6 | 壁の補正 | δ_E → 次の壁。規則つきの繰り返し |
| 7 | 成果物 | 手元か AWS か (FreeCAD の有無) |
| 8 | 設計報告 | |

- 各段は手元か AWS かで動き、計算資源の差 (起動・待ち合わせ・容量) は実行役の「実行先」の部品が隠す。
- AWS の部品は skill `forge-aws-run` の規則を実装したもの。

### 4.3 失敗の返り値と次の手

- 各段の失敗は列挙の返り値で返す。例:
  - 停止規則の `CENSORED`・`FAILED:diverged`・`FAILED:forge_exit`
  - 関門の `blocked:<条件>`
  - 格子の品質の `FAIL`
- レシピがその返り値に「選べる手の一覧」を対応づける。例:
  - 同じ設定でもう 1 回 (上限 1 回)
  - 前の段からやり直す
  - CFL を下げる段階起動 (`procedures/divergence-and-startup.md` の手順)
  - 壁の補正をもう 1 回 (上限つき)
  - 打ち切って人に返す
- 一覧に 1 つしか無い・前提で 1 つに決まる場合は、実行役がそのまま選ぶ。
- 2 つ以上が残る場合と、一覧に無い失敗の場合は、判断役に問う (§4.4)。

### 4.4 判断役の口 (止まったときだけ)

- **問い** (実行役が書く、JSON と同じ内容の Markdown):
  - どの段で、何の判定が、どの値で止まったか
  - 根拠のファイル (パスと sha256)
  - 選べる手の一覧 (各手の前提と、選んだら実行役が何をするか)
  - 過去の同じ問いの記録
- **答え** (判断役が返す): 一覧からちょうど 1 つの手、理由 (数行)、確信の度合い (high / medium / low)。
  - 自由に台本を書かせない。閾値・規則を変えさせない。
- **実行役の検査**: 次のどれかに当たれば、その答えは採らず「人に返す」で止まる。
  - 形式が壊れている
  - 一覧に無い手を選んでいる
  - 手の前提を満たさない (上限を超える等)
  - 確信が low
  - 判断役につながらない、時間切れ
- **判断役の切り替え**: 設定ファイル 1 つで選ぶ (今の `~/.config/forge/diagnose-backend` を広げる)。
  - 選べる判断役: `fable` (Claude Code の `diagnostician`)・`codex`・`claude` (別インスタンス)・`local` (OpenAI 互換の API。llama.cpp・vLLM・Ollama など、URL とモデル名を設定)。
  - `codex_review.py` の切り替えの仕組みに `local` を足す。
- **記録**: 問いと答えと検査の結果を台帳に残す。後から別の判断役で同じ問いを再生して比べられる。
- **判断役ごとの試験**: 過去の止まり方を問いの形にした組を作る。
  - 例: run_0012 の未成立、run_0010 の V5b FAIL、run_0002 の未整定、ディスクの満杯、発散。
  - 期待する手 (または「人に返す」) と比べて、判断役ごとの成績を出す。
  - ローカル LLM に任せてよい場面の線引きは、この成績で決める (決めるのはユーザ)。

### 4.5 入口と skill

- CLI: `forge-design run <card> [--resume] [--executor local|aws]`、`forge-design status <ledger>`、`forge-design answer <ledger> <手>` (人が答えるとき)。
- skill は薄くする。中身は次の 3 つだけで、判断の中身はコードとレシピが持つ。
  1. ユーザの要求から要求カードを作る (カードの欄と許容の決め方)
  2. `forge-design run`
  3. 結果か「人に返す」の報告をユーザに伝える
- 知識のない AI が台本を書いて先へ進む余地を作らないため、AGENTS.md にも次を書く: 「設計の run は `forge-design` 経由。止まったら問いの記録を読んでユーザに伝える」。

### 4.6 台帳 (再開と追跡)

- run の外に置く 1 つの JSON。名前は設計の単位で 1 つ (例: `design/<card>/ledger.json`)。
- 中身:
  - カードとレシピの sha256
  - 段ごとの状態 (未着手・実行中・完了・失敗)
  - 入出力のパスと sha256
  - VERDICT の要約
  - 失敗の返り値、選んだ手、判断役の問いと答え
- 再開すると、完了した段は入出力の sha256 を照らして飛ばす。照合が合わなければ止める。

### 4.7 諮問 (diagnostician、2026-10-16) を受けた改訂 — 上の §4.1〜§4.6 より優先する

Fable の `diagnostician` に諮った (ユーザ指定; 結論: 「方向は妥当だが新しい枠組みとして作らず、承認済みの campaign-recipe の実行器 + ledger と forge-gate の `judge_run.py`・`launch_gate.py`・`RECIPE.json` の上に載せる。判断役の口は確信の自己申告でなく検証できる引用と手ごとの最悪コストで縛る。格子は一般化せず包むだけ。V1 は GPU 不要の入力の一致を合格に」)。全件採用:

1. **二重化をやめる (Critical)**: この plan は親 plan の §5 ステップ 5 (campaign 実行器・ledger) と入口 (要求カード → search の無い campaign) の実装。**汎用 DAG・新しい状態機械の枠組みは作らない** (親 plan §2 で codex が却下: δ* は壁を変えて再評価する反復で段の直列に乗らない)。段は線形の列、δ* の繰り返しは recipe の内部。語彙を 1 つに: レシピ名 `axismach.*/vN`、段の判定 `accepted / rejected / incomplete`、来歴 `execution_plan.json` (+ `StageManifest`・`MAIN_STAGE_STARTED`・`HANDOFF_GATE.json`)。§4.1 の「決定的な状態機械」は「親 plan の実行器」と読み替える。
2. **判定の正本の 1 本化 (Major)**: forge-gate の `judge_run.py` (`RUN_VERDICT.json`、PASS / FAIL / INDETERMINATE / NOT_APPLICABLE) を run の総合の判定の正本とし、v2 の run では収束・停止の出典に `SETTLE_VERDICT.json`・停止規則 v2 の状態を使う (今は互いに知らない)。forge-gate のブランチと feature/card-local-mesh を先に統合してから実行器を作る (同じ主セッションの持ち物)。`judge_run` は CENSORED を FAIL に写すが、実行器にとって CENSORED は分岐 (延長 1 回 / 人に返す) — 写像はレシピが持つ。
3. **失敗の返り値は既存の語彙の和 (Major)**: `RUN_VERDICT.checks`、停止規則 v2 の status + failures、`HANDOFF_GATE` の blocked の条件、`LAUNCH_GATE` の FAIL、`MESH_QUALITY`。新しい列挙を作らない。これに**原因の分類** (一過性 (ssh 切断・インスタンス停止) / 資源 (ディスク満杯 `Write failed`) / 数値 (detectNaN・発散) / 照合の不一致) を足す。「同じ設定でもう 1 回」は原因が一過性のときだけ。一覧に無い失敗の手は「人に返す」だけ (判断役に問わない)。一覧への追加は「ツールで検出できる印 + 許す手 + 問いの組への 1 件」の 3 点で。手ごとに最悪コスト (GPU 時間・ディスク・やり直す段) と前提を書く。数値を変える手 (CFL を下げる等) は `RECIPE.json` の理由つき逸脱 + `launch_gate` を通す形でだけ (判断役は番号を選ぶだけ)。
4. **判断役の答え (Major)**: 一覧から 1 つ + 理由 + **問いの digest からの引用 1 つ以上** (量の名前・値・出典)。実行器は引用が digest に実在するかを照合し、無い・違う → 人に返す。**確信の度合いは記録だけ** (自己申告は較正されていない — 推測、一般論)。「証拠不足: <要るもの>」を正式な選択肢に (= 人に返す + 要るもの)。判断役の種類ごとに選べる手を**安全集合** (上限つきで安い手) に制限できる (レシピで指定)。行動に結びつかない「所見」の問い (記録だけ) を別の型で。
5. **問いは自己完結の digest (Minor)**: 実行器が抽出した判定の要約値 (各 VERDICT の種類・末尾の窓の統計・blocked の条件・原因の分類) を本文に入れ、ファイルはパスと sha256 だけ。ローカル LLM (道具なし) でも答えられる。判断役の比較は全員 digest だけで。
6. **問いの組 (Major)**: 過去の止まり方を値を摂動した変種で増やした正例、正解が「人に返す」の負例 (一覧に無い失敗・照合の不一致・情報の欠けた問い) を同数、自明な方針の対照 (常に先頭・常に人に返す・多数決)、正解と digest を判断役を回す前に commit (事前登録)、成績はコストの重み (誤って高い手を選んだ回数 × その手のコスト)。ローカル LLM に任せる範囲は「誤っても損失が上限 1 回の安全集合の手だけ」から (決めるのはユーザ)。
7. **実行器の動く場所 (Major)**: run と同じホスト (AWS) で nohup、台帳もそこ。手元は rsync の写しを読むだけ (`forge-design status`)。再開は同じホストで `--resume`。段の完了は完了の印 (`STOPPED_AT`・`MAIN_STAGE_STARTED` と同型) を確かめてから次へ。idle 自動停止との関係 (実行器の準備の隙間が「forge なし」と数えられる) は未決 (§5.1)。
8. **格子は包むだけ (Major)**: 段の契約は「問題 YAML の mesh 節 (表を含む) + 変換後の格子の不変部分のハッシュ (forge-gate `mesh_provenance`) + 品質 VERDICT」。実装は 2 つ: 雛形 (`wt_card.estimate_mesh` → mesh2d; case/70 など) と手組み (`hand_mesh.py` をそのまま包む; case/71)。係数の一般化は別項目 (F、格子の感度の A/B つき) で後。
9. **台帳 (Minor)**: 巨大な出力は sha256 でなく、格子は不変ハッシュ、全場は停止の場だけ sha256・ほかは大きさ + 完了の印。既存の記録 (`HANDOFF_GATE.json`・`STOP_RULE_STATE.json`・`RUN_VERDICT.json`) は sha で参照して中身を複製しない。
10. **AGENTS.md (Major)**: 規則は「**レシピのある機種** (最初は風洞ノズル) の設計の run は `forge-design` 経由」に限定。settle (収束確認) → gate (起動の関門) → この plan (入口) の 3 つが同じ節を触るので、書き換えは 1 セッションがこの順でまとめて行い、調整板に予告。`run_case.sh` に触る merge は run の走っていないときだけ。
11. **要求カード (Major)**: カードの欄ごとに「ユーザから聞く / レシピの既定 / 導出」を表にし、聞くべき欄を聞かずに埋めたら止まる (既定の無い欄は必須)。最大の漏れ口。
12. **集計 (Minor)**: 設計 1 件あたりの「判断役に問うた回数」「人に返した回数」「決定表で自動に決まった回数」を台帳から出す。2 つ以上が残る分岐はできるだけ決定表 (レシピのデータ) で潰す。

## 5. 実装ステップ

1. 段の契約と台帳の形式、失敗の返り値の列挙 (`design/forge_design/orchestrator/` 新規)。
2. 風洞ノズルの段の部品。今の部品を包み、case/71 の台本の中身を移す。
3. 格子の段の一般化 (`hand_mesh.py` → レシピから呼べる形)。
4. 判断役の口と切り替え (`local` を足す)。判断役ごとの試験の組。
5. CLI・skill・AGENTS.md。
6. 検証 (§6)。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | ~~設計方針の諮問 (diagnostician、ユーザ指定)~~ (2026-10-16 済み) | 判断: 2026-10-16・方向は妥当、新しい枠組みにせず campaign-recipe の実行器 + ledger と forge-gate の判定・関門の上に載せる (§4.7、全件採用)。次は codex の plan 段レビュー | F |
| 1b | codex の plan 段レビュー | §4.7 を反映した plan で `codex_review.py --stage plan` | F |
| 2 | forge-gate と feature/card-local-mesh の統合、`judge_run.py` が v2 の `SETTLE_VERDICT.json`・停止規則 v2 を受ける 1 本化 | §4.7 の 2。合格: judge_run の試験 + v2 の run (run_0011) で RUN_VERDICT が SETTLE_VERDICT と矛盾しない | O |
| 3 | 失敗の返り値の和・原因の分類・手の一覧とコスト (レシピのデータ) | §4.7 の 3。合格: 単体試験 (既存の記録から分類が決まる、一覧に無い失敗は人に返す) | O |
| 4 | 実行器 (親 plan の §5 ステップ 5) と台帳・入口の CLI | §4.7 の 1・7・9。合格: §6 V1a | O |
| 5 | 格子の段の契約と 2 つの実装 (雛形 / 手組みを包む) | §4.7 の 8。合格: §6 V2a・V2b | O |
| 6 | 判断役の口 (digest・引用の照合・安全集合)・`local` の backend・問いの組 | §4.7 の 4〜6。合格: §6 V3 | O |
| 7 | 要求カードの欄の表と止まり方、skill、AGENTS.md (レシピのある機種に限定) | §4.7 の 10・11。合格: §6 V4 の準備 | O |
| 8 | idle 自動停止と実行器の準備の隙間 | §4.7 の 7 の未決 | F |
| 9 | 格子の係数の一般化 (後) | 格子の感度の A/B つき | F |

## 6. 検証 (案、codex の plan 段の後に確定)

- **前提**: settle plan の V6 (`run_0014_chain_v2`、設計チェーンの通しの 1 pass) が通っていること。
- **V1a (GPU 不要、合格の条件)**: 実行器が case/71 の 1 pass (3 係数の壁、V6 と同じ入力) のために作る入力一式 (`solverConfig.yaml`・`solverConfig_prestage.yaml`・`bcondConfig.yaml`・`settle_criteria.yaml`・`settle_nodes.txt`・`RECIPE.json`・格子の不変ハッシュ・`MAIN_STAGE_STARTED` の sha) が V6 の run のものとバイト同一、かつ台本 (`p5_ns.py`・`p5_all.sh`) が行っていた検査がすべて実行器の関門として存在し、負例で止まる。
- **V1b (GPU 1 pass、記録)**: VERDICT の種類が V6 と同じ (REACHED・関門 pass・`delta_E.status ok`・成果物の実在・非有限 0)。停止の step・δ_E・`delta_r_next.csv` の差は表にして記録 (数値の許容は、同じホスト・同じ入力の反復 2 本で床を取るまで事前登録しない)。
- **V2a**: 手組みの格子を包んだ段が case/71 の節点の座標と一致。**V2b**: 雛形の格子で別のカード (例: case/70 M5 の空気) が品質の関門を通る。
- **V3 判断役**: 摂動した正例・負例 (同数)・自明な対照・コストの重み・事前登録 (正解と digest を先に commit)・全判断役を digest だけで。形の壊れた答え・一覧に無い手・前提を満たさない手・引用の不一致を実行器が全部拒む (偽の判断役で)。
- **V4 知識のない AI** (V1a の後): Haiku (と、使えればローカル LLM) に要求の文と skill だけを渡し、カードから成果物まで。止まった所、台本を書こうとしたか、カードのどの欄を自分で決めたか、を記録。
- **集計**: 設計 1 件あたりの判断役に問うた回数・人に返した回数・決定表で決まった回数。
- 判定不能は合格にしない。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| diagnose (Fable、ユーザ指定) | 2026-10-16 | diagnostician の返答 (ブリーフ [briefs/2026-10-16-design-orchestrator.md](../../notes/reviews/briefs/2026-10-16-design-orchestrator.md); 本文は §4.7 に要約) | C1/M10/m3 | 全件**採用** (§4.7・§5.1・§6 を書き直し) |

## 7. 影響範囲

新しい道具を足す。既存の部品は包むだけで振る舞いを変えない。変えるのは AGENTS.md の 1 項目 (設計の run は `forge-design` 経由) と、skill の追加。

## 8. 完了条件

V1〜V4。ただし V4 は「詰まる所の一覧と、その対処」で完了とする (全部が通ることは求めない)。

## 9. 変更ログ

- `2026-10-16` — diagnostician (Fable) の諮問を全件採用: 親 plan (campaign-recipe) の §5 ステップ 5 の実装に位置づけ直し、forge-gate の判定・関門の上に載せる、判断役の答えは引用で縛る、格子は包むだけ、V1 を V1a (入力の一致) と V1b (種類の一致) に分けた (§4.7)。

- `2026-10-16` — 起票 (draft)。ユーザ「そのやり方でよいか、diagnosticianにも相談した上で、進めてほしい」。
