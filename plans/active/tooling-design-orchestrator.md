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

- **やる (最初の範囲)**: 風洞ノズル (axismach) で、search の無い campaign を要求カードから回す。
  - 親 plan の未実装の契約: 最小の campaign、入力の参照の解決、pass ごとの `execution_plan.json`、全段の設定の照合。
  - 実行器と台帳: 失敗の返り値と原因の分類、手の一覧 (レシピのデータ)。
  - 格子の段: 包むだけ。
  - 判断役の口: 行動の選択の専用の呼び出し。
  - 要求カードの値の出所の記録。
- **やる**: case/71 (手組みの格子) と case/70 (雛形の格子) で検証する。知識のない AI (Haiku・ローカル LLM) で再現する。
- **やらない**:
  - 汎用 DAG・新しい状態機械の枠組み・新しい VERDICT の語彙 (親 plan §2)。
  - 格子の係数の一般化 (別項目、F)。
  - ベル・デュアルベル・SERN のレシピ (各セッションの成果を後から載せる)。
  - 機種ごとの物理の公差を決めること (ユーザと決める)。
  - ソルバの演算の変更。

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

経緯: Fable の `diagnostician` に諮り (2026-10-16、ユーザ指定)、codex の plan 段のレビュー (GO-with-changes、C0/M5/m2) を受けた。どちらも全件採用した現行の仕様を、以下に 1 つにまとめる (§6.1)。

### 4.1 位置づけと語彙

- この plan は、親 plan `tooling-design-problem-campaign-recipe` の §5 ステップ 5 (campaign 実行器・ledger) と、要求カードからの入口の実装である。
- **汎用 DAG は作らない**。段は線形の列で、δ* の繰り返しは recipe の内部に置く。
- 語彙は親 plan に合わせる:
  - レシピ名は `axismach.*/vN`
  - 段と campaign の判定は `accepted / rejected / incomplete`
  - 来歴は `execution_plan.json` (+ `StageManifest`・`MAIN_STAGE_STARTED`・`HANDOFF_GATE.json`)
- **前提 (codex M1)**: 親 plan の未実装の契約を、実行器より先に作る。
  - search の無い最小の campaign。
  - 入力の参照の解決 (snapshot・格子・物性・バイナリの識別の固定)。
  - pass ごとの `execution_plan.json`。
  - 全段の設定の照合。
  - 例: `prepare_ns` は `ic_from` のディレクトリの最大番号の場を選ぶ (`runner_axismach.py:2511`)。計画を確定した後に参照の run の出力が増えても、入力が変わらないようにする (負例の試験)。

### 4.2 判定の正本と段の判定

- **run の総合の判定**: forge-gate の `judge_run.py` (`RUN_VERDICT.json`)。
  - v2 の run では、収束・停止の出典に `SETTLE_VERDICT.json` と停止規則 v2 の状態を使う (1 本化; forge-gate のブランチと統合してから)。
  - `SETTLED` でも壁解像が FAIL なら総合は FAIL。
  - CENSORED の扱い (延長 1 回 / 人に返す) はレシピが持つ。
- **計算の成立と設計の要求の達成を分ける (codex M2)**: recipe に最終の qualification を書く。

| 判定 | 条件 |
| --- | --- |
| `accepted` | 要求の公差 (例: 風洞ノズルの V5b・V5c)・寸法・壁解像・必要な成果物のファイルをすべて満たす |
| 補正へ戻す | 評価済みで公差が未達 |
| `rejected` | 反復の上限に達した |
| `incomplete` | 必要な判定・成果物が欠けた |

  - 関門 (`handoff_gate_v2`)・`judge_run` は計算の成立と受け渡しを扱う。公差は別の判定にする。

### 4.3 失敗の返り値・原因の分類・手

- **返り値**: 新しい列挙は作らない。既存の語彙の和にする。
  - `RUN_VERDICT.checks`
  - 停止規則 v2 の status + failures
  - `HANDOFF_GATE` の blocked の条件
  - `LAUNCH_GATE` の FAIL
  - `MESH_QUALITY`
  - qualification の未達
- **原因の分類**:
  - 一過性 (ssh の切断・インスタンスの停止)
  - 資源 (ディスクの満杯 `Write failed`)
  - 数値 (detectNaN・発散)
  - 照合の不一致
- **手の一覧 (レシピのデータ)**: 返り値と原因に対応づける。
  - 手ごとに最悪のコスト (GPU 時間・ディスク・やり直す段) と前提を書く。
  - 「同じ設定でもう 1 回」は原因が一過性のときだけ。
  - 数値を変える手 (CFL を下げる段階起動など) は、`RECIPE.json` の理由つきの逸脱 + `launch_gate` を通す形でだけ。判断役は番号を選ぶだけ。
  - 一覧に無い失敗の手は「人に返す」だけ。判断役に問わない。
  - 一覧に足すときは「ツールで検出できる印 + 許す手 + 問いの組への 1 件」の 3 点をそろえる。
- 2 つ以上が残る分岐は、できるだけ決定表 (レシピのデータ) で潰す。

### 4.4 判断役の口 (止まって 2 つ以上の手が残るときだけ)

- **問い**: 自己完結の digest にする。
  - 判定の要約値: 各 VERDICT の種類・末尾の窓の統計・blocked の条件・原因の分類。
  - ファイル: パスと sha256。
  - 選べる手の一覧: 各手の前提・コスト・選んだら実行器が何をするか。
- **答え**: 一覧から 1 つ、理由、**digest からの引用 1 つ以上** (量の名前・値・出典)。
  - 「証拠不足: <要るもの>」も正式な選択肢 (= 人に返す)。
  - 確信の度合いは記録だけ。拒否の条件にはしない。
- **実行器の検査**: 次のどれかなら人に返す。
  - 形式が壊れている
  - 一覧に無い手
  - 前提を満たさない
  - 引用が digest に無い・値が違う
  - 判断役につながらない、時間切れ
- **安全集合**: 判断役の種類ごとに、選べる手を「上限つきで安い手」に制限できる (レシピで指定)。
  - ローカル LLM は、誤っても損失が上限 1 回の手だけから始める (広げるのはユーザの決定)。
- **呼び出しの分離 (codex m6)**: 行動の選択の呼び出しは、外部レビュー用 (`codex_review.py`) と分ける。
  - 全 backend (`fable`・`codex`・`claude`・`local` (OpenAI 互換の API)) で、道具なし・同じ digest・同じ答えの形に限定する。
  - 共用するのは、接続・timeout・試行の記録の部分だけ。
- **所見の問い**: 行動に結びつかない所見の問いは別の型にする (記録だけ・人が読む)。
- **問いの組**: 判断役ごとの成績を出す。
  - 中身: 過去の止まり方を値を摂動して増やした正例、正解が「人に返す」の負例 (同数)、自明な方針の対照 (常に先頭・常に人に返す・多数決)。
  - 正解と digest は、判断役を回す前に commit する (事前登録)。
  - 成績はコストの重みで出す。

### 4.5 実行器・台帳・再開 (codex M3)

- **実行の場所**: 実行器は run と同じホスト (AWS) で nohup で動き、台帳もそこに置く。手元は rsync の写しを読むだけ (`forge-design status`)。
- **台帳**:
  - 単一の書き手のロック、原子的な更新 (一時ファイル → rename)。
  - 段の試行ごとの attempt ID と、起動したプロセスとの結び付け (PID・開始時刻・ホスト)。
  - 再試行の回数と累積の予算 (GPU 時間・ディスク)。再開の後も保持する。
- **`--resume`**: 実行中の段を次の手順で照合する。
  - 起動したプロセスが生きていれば、再起動せず監視に戻る。
  - 止まっていれば、完了の証拠 (`STOPPED_AT`・`MAIN_STAGE_STARTED`・関門の記録) を確かめる。
  - 不明なら人に返す。
  - 再試行は新しい run (新しい run 名) に分ける。既存の runner の再呼び出しが開始・停止の記録や途中の出力を消す (`runner_axismach.py:2946`・`:2981`) のを避けるため。
- **巨大な出力の扱い**: sha256 で照合しない。
  - 格子は不変部分のハッシュ (forge-gate `mesh_provenance`)。
  - 全場は停止の場だけ sha256、ほかは大きさ + 完了の印。
  - 既存の記録は sha で参照し、中身を複製しない。
- **idle 自動停止との関係**: 実行器の準備の隙間が「forge なし」と数えられる。AWS 対応の前提として先に決める (§5.1)。

### 4.6 格子の段

- 契約: 「問題 YAML の mesh 節 (表を含む) + 変換後の格子の不変部分のハッシュ + 品質の VERDICT」。
- 実装は 2 つ:
  - 雛形: `wt_card.estimate_mesh` → mesh2d (case/70 など)。
  - 手組み: `hand_mesh.py` をそのまま包む (case/71)。
- 係数の一般化は後 (F、格子の感度の A/B つき)。

### 4.7 要求カード (codex M5)

- カードの値に出所を付けて保存する。出所は 4 種類: ユーザの入力の参照 (発言・指示の記録)、レシピの既定、導出、未確認の提案。
- ユーザの指定が要る欄 (公差・壁温・成果物の選択など、欄ごとの表で決める) は、ユーザの指示か確認済みのカードの版に結び付ける。未確認なら起動の前に止める。
- 負例: 欠けたカードだけでなく、「値はあるが根拠が無い」カード。

### 4.8 入口・skill・AGENTS.md

- **入口**: CLI は次の 3 つ。
  - `forge-design run <card> [--resume] [--executor local|aws]`
  - `forge-design status <ledger>`
  - `forge-design answer <ledger> <手>` (人が答えるとき)
- **skill**: 薄くする。中身は「カードを作る (出所つき)」「回す」「報告する」だけ。
- **AGENTS.md**:
  - 規則は「レシピのある機種 (最初は風洞ノズル) の設計の run は `forge-design` 経由」に限る。
  - 同じ節を書き換える plan が 3 つある。書き換えは 1 セッションが次の順でまとめて行い、調整板に予告する: settle (収束確認、済み) → gate (起動の関門) → この plan (入口)。
- **集計**: 設計 1 件あたりの判断役に問うた回数・人に返した回数・決定表で決まった回数を、台帳から出す。

## 5. 実装ステップ

1. 親 plan の未実装の契約: 最小の campaign・入力の参照の解決・`execution_plan.json`・設定の照合 (`design/forge_design/` の campaign の置き場所; 親 plan §4.1〜§4.4)。
2. forge-gate との統合と `judge_run` の v2 対応、qualification の判定。
3. 失敗の返り値の和・原因の分類・手の一覧 (レシピのデータ)。
4. 実行器と台帳・再開 (親 plan §5 ステップ 5)。
5. 格子の段 (包む)。
6. 判断役の口 (行動の選択の専用の呼び出し)・問いの組。
7. 要求カードの出所・CLI・skill・AGENTS.md。
8. 検証 (§6)。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | ~~設計方針の諮問 (diagnostician、ユーザ指定)~~ (2026-10-16 済み) | 判断: 2026-10-16・方向は妥当、新しい枠組みにせず campaign-recipe の実行器 + ledger と forge-gate の判定・関門の上に載せる (全件採用、§4) | F |
| 1b | ~~codex の plan 段レビュー~~ (2026-10-16 済み) | GO-with-changes、C0/M5/m2、全件採用 (§4 を書き直し、§6.1) | F |
| 2 | 親 plan の未実装の契約 (§4.1 の前提) | search の無い最小の campaign・入力の参照の解決 (snapshot・格子・物性・バイナリの固定)・`execution_plan.json`・設定の照合。合格: 計画を確定した後に参照の run の出力が増えても入力が変わらない負例・照合の不一致で止まる負例 | O |
| 3 | forge-gate と feature/card-local-mesh の統合、`judge_run` の v2 対応、qualification (§4.2) | 合格: run_0011 で RUN_VERDICT が SETTLE_VERDICT と整合、`SETTLED` でも壁解像 FAIL なら総合 FAIL、公差未達でも整定している例・必須の CAD の欠けの例の試験 | O |
| 4 | idle 自動停止と実行器の準備の隙間 (§4.5) | AWS 対応の前提。決め方を諮る | F |
| 5 | 失敗の返り値の和・原因の分類・手の一覧とコスト (§4.3) | 合格: 既存の記録から分類が決まる、一覧に無い失敗は人に返す、一過性でない原因で「もう 1 回」を選ばない | O |
| 6 | 実行器・台帳・再開 (§4.5) | 合格: §6 V1a・障害の注入の試験 (起動の直後・結果の公開の直後の強制終了、同時の `--resume`、古い問いへの答え) | O |
| 7 | 格子の段の契約と 2 つの実装 (§4.6) | 合格: §6 V2a・V2b | O |
| 8 | 判断役の口・`local` の backend・問いの組 (§4.4) | 合格: §6 V3 | O |
| 9 | 要求カードの出所・CLI・skill・AGENTS.md (§4.7・§4.8) | 合格: 「値はあるが根拠が無い」カードで止まる、§6 V4 の準備 | O |
| 10 | 格子の係数の一般化 (後) | 格子の感度の A/B つき | F |
| 11 | **収束の判定 (停止規則 v2・プリセット `wt_nozzle_v3`) の別の寸法のノズルでの確認** (2026-10-16 に plan `tooling-settle-convergence-and-monitor` の完了条件から移した; ユーザ決定「後続の plan へ移す」) | 寸法が大きく違うノズル (小さいノズル) を 1 本、新しい出力を有効にして停止規則 v2 で回し、オンラインと再生の判定の一致・停止の時の出力の揃い・停止の後の区間 (V3(c) と同じ形の自由継続の再生) を確かめる。プリセットの閾値は暫定のまま使い、結果を見て較正しない。合格の形は settle plan §6 V3(b)・V3(c)・V4 に倣う | F |

## 6. 検証 (事前登録は各項目の実装の前に確定する)

- **前提**: settle plan の V6 (`run_0014_chain_v2`) が通っていること。
- **V1a (GPU 不要; codex M4 で書き直し)**: 次の 2 つを分けて確かめる。
  - **決定的な入力の比較**: 実行器が case/71 の 1 pass のために作る起動の入力 (設定・境界条件・基準・節点の表・`RECIPE.json`・初期場・バイナリ・変換器・数値に効く環境変数) を、V6 の run のものと比べる。
    - 生成時刻・配置のパスは、項目を限って正規化する。
    - 前の段の成果物は V6 のものを固定して使う。
  - **来歴の照合**: 参照先の内容とハッシュの結び付けを別に検査する。
  - あわせて、台本 (`p5_ns.py`・`p5_all.sh`) が行っていた検査がすべて実行器の関門として存在し、負例で止まること。
- **V1b (GPU 1 pass、記録)**: 次の種類が V6 と同じであること。
  - REACHED・関門の pass・`delta_E.status ok`・成果物の実在・非有限 0。
  - qualification の判定 (V5b・V5c は今の壁では未達の見込み → 「補正へ戻す」が正しい結果)。
  - 段の間の引き継ぎと完了の記録。
  - 数値の差は表にして記録する。数値の許容は、同じホスト・同じ入力の反復 2 本で床を取るまで事前登録しない。
- **V2a**: 手組みの格子を包んだ段が、case/71 の節点の座標と一致する。
- **V2b**: 雛形の格子で、別のカード (case/70 M5 の空気) が品質の関門を通る。
- **V3 判断役 (§4.4)**: 問いの組で判断役ごとの成績を出す。偽の判断役で、形の壊れた答え・一覧に無い手・前提を満たさない手・引用の不一致を、実行器が全部拒むこと。
- **V4 知識のない AI (V1a の後)**: Haiku (と、使えればローカル LLM) に要求の文と skill だけを渡し、カードから成果物まで回させる。次を記録する。
  - 止まった所
  - 台本を書こうとしたか
  - カードのどの欄を自分で決めたか (出所の記録で)
- **集計**: 設計 1 件あたりの、判断役に問うた回数・人に返した回数・決定表で決まった回数。
- 判定不能は合格にしない。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| diagnose (Fable、ユーザ指定) | 2026-10-16 | diagnostician の返答 (ブリーフ [briefs/2026-10-16-design-orchestrator.md](../../notes/reviews/briefs/2026-10-16-design-orchestrator.md); 本文は §4 に統合) | C1/M10/m3 | 全件**採用**: 親 plan の実行器 + ledger と forge-gate の判定・関門の上に載せる (二重化をやめる)、判定の 1 本化、返り値は既存の語彙の和 + 原因の分類、答えは引用で照合 (確信は記録だけ)、digest、問いの組の作り方、実行の場所、格子は包むだけ、台帳の巨大な出力、AGENTS.md の範囲と順序、カードの欄の表、集計 |
| plan | 2026-10-16 | [2026-10-16-tooling-design-orchestrator-plan.md](../../notes/reviews/2026-10-16-tooling-design-orchestrator-plan.md) | GO-with-changes, C0/M5/m2 | 全件**採用** (§2・§4・§5・§6・§7 を書き直し): M1 親 plan の未実装の契約 (最小の campaign・入力の固定・execution_plan・設定の照合) を実行器より先に (§4.1、§5.1 #2)。M2 計算の成立と設計の要求の達成を分ける qualification (§4.2、#3)。M3 再開の契約 (ロック・原子的な更新・attempt ID・プロセスとの結び付け・新しい run での再試行・障害の注入の試験; idle 自動停止を AWS 対応の前提へ) (§4.5、#4・#6)。M4 V1a を決定的な入力の比較と来歴の照合に分ける (§6)。M5 カードの値の出所 (§4.7、#9)。m6 行動の選択の呼び出しを外部レビュー用と分け、全 backend で道具なし・同じ digest (§4.4)。m7 旧い記述を残さず書き直した |

## 7. 影響範囲

- 新しい道具を足す: campaign の最小の実装・実行器・台帳・判断役の口・CLI・skill。
- 変えるもの:
  - forge-gate の `judge_run.py` (v2 の run の出典)。
  - 判断役の切り替え (`local` の backend; 行動の選択の呼び出しは別)。
  - AGENTS.md の 1 項目 (レシピのある機種の設計の run は `forge-design` 経由)。
- 既存の部品 (runner・停止規則・関門・成果物) は呼び出すだけにする。振る舞いを変えるときは、それぞれの plan に書く。

## 8. 完了条件

V1〜V4。ただし V4 は「詰まる所の一覧と、その対処」で完了とする (全部が通ることは求めない)。

## 9. 変更ログ

- `2026-10-16` — codex の plan 段レビュー (GO-with-changes、C0/M5/m2) を全件採用し、§2・§4・§5・§5.1・§6・§7 を採用済みの方針に直接書き直した (旧い記述を残さない)。
- `2026-10-16` — diagnostician (Fable) の諮問を全件採用: 親 plan (campaign-recipe) の §5 ステップ 5 の実装に位置づけ直し、forge-gate の判定・関門の上に載せる、判断役の答えは引用で縛る、格子は包むだけ、V1 を V1a (入力の一致) と V1b (種類の一致) に分けた (§4.7)。

- `2026-10-16` — 起票 (draft)。ユーザ「そのやり方でよいか、diagnosticianにも相談した上で、進めてほしい」。
