# Claude 代替レビュー: tooling-review-fallback-backend (plan)

> **Claude の別インスタンスによる代替** (独立性は文脈の分離だけ。依頼側と同じ系統のモデルでありうる; codex・Fable と別系統の判定ではない)

- **plan**: [`plans/active/tooling-review-fallback-backend.md`](../../plans/active/tooling-review-fallback-backend.md)
- **stage**: `plan`
- **date**: 2026-10-14
- **commit**: `88675245` (feature/nozzle-cooled-wall-settling)
- **分類**: `ok`
- **backend**: `claude` — 選んだ理由: 明示 (`--backend claude`)
- **要求したモデル**: `haiku`
- **実効のモデル**: `claude-haiku-5-5` (modelUsage: `claude-haiku-5-5` out 50066)
- **CLI の版**: `2.1.296 (Claude Code)`
- **依頼側との関係**: 別 (依頼側 `claude-opus-5-5`、実効 `claude-haiku-5-5`) (同じ Claude 系統)
- **permission_denials**: 1 件 — Grep {"pattern": "--safe-mode|--no-session-persistence|--strict-mcp-config", "path": "/home/sano/.local/bin/claude", "output_mode": "count"}
- **道具**: `Read,Grep,Glob` (`--safe-mode`、Bash なし) / 管理ポリシー: 無し (`/etc/claude-code/managed-settings.json`, `/etc/claude-code/managed-settings.d`, `/etc/claude-code/managed-mcp.json` が無い)
- **実行**: effort `high`, 4.1 min, rc=0, total_cost_usd 0.1002
- **試行のログ**: `notes/reviews/2026-10-14-tooling-review-fallback-backend-plan-claude.claude.log` (claude: ok)
- **証拠一式**: `notes/reviews/2026-10-14-tooling-review-fallback-backend-plan-claude.evidence/` (git 追跡外)
- **判定**: **GO-with-changes**, 指摘 C1/M9/m8

本文は Claude の別インスタンスの最終メッセージ (JSON の `result`) をそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**

計画の方向 (判定役が使えないときに文脈を分けた Claude で代える) は妥当で、課題の同定も正しい。ただし、受入れの試験 V1 は答え合わせが現行コードと合わず、記録ゲートが未成立の記録を通す。規則 (§4.3) を採用する前に、この 2 点とゲートを直す必要がある。実装は作業ツリーに既にあるので、「着手してよいか」より「規則を採用してよいか」が実際の判断点になる。

確認範囲: 計画・`codex_review.py`・`check_plans.py`・AGENTS.md・証拠一式を読んだ。`claude` の CLI フラグ (`--safe-mode` など) は、バイナリの読取りが拒否されたため確認できていない。計画 §3 の実測は計画の自己申告のままで、再現していない。子のインスタンスの既存ログ (`notes/reviews/…-plan-claude.claude.log`) は独立性のため読んでいない。

## 指摘一覧

### Critical

**C1. V1 の正答が現行コードと一致しない (受入れ試験の答え合わせが逆になる)**
- 根拠: 計画 §6 V1 は「`codex_review.py:295` 付近で『rc 非 0 でも result 非空なら `if not result` の分岐に入らない』が正答」とする。作業ツリーの 295 行は `read_switch` の docstring で、`if not result` の文字列は現行版に無い (Grep 0 件)。実際の分岐は 429 行 `elif rc != 0 or not at.result:` で、rc≠0 かつ result 非空でもこの分岐に入り、436–440 行で `failed` になる。したがって B の前提「非 0 なら必ずこの分岐に入る」は**現行版では正しく**、正答が逆になる。
- 対案: V1 の問いを、現行版 (commit を固定し、その版を証拠一式に入れる) で答えが一意に決まるものに差し替える。例: 「rc=0 で最終メッセージが『Error: rate limit reached』の文のとき、codex の unavailable として claude に切り替わるか」→ 正答は「切り替わらず `incomplete`」(441–442 行、rc=0 では `classify_body` だけを通り、`UNAVAILABLE_PATTERNS` は見ない)。この正答は V2(d) の偽 CLI 試験で機械的に先に確かめ、そのうえで LLM の試験に使う。

### Major

**M1. V1 は「誤った前提に騙されない」しか測らず、「欠陥を見つける」を測らない**
- 根拠: V1 の合格条件は、誤った前提を退けることだけ。どんな問いにも「GO」「問題なし」と返す判定役は、前提に同調しない答え方ならこの試験を通りうる。計画 §6 は V1 を「最低限の棄却試験」と書いているが、受入れ条件として「判定役が機能する」の根拠にはならない。
- 対案: 既知の欠陥を仕込んだ検出試験を追加する。リポジトリ外の写しに、既知の欠陥を 2〜3 個入れる (例: 指摘数の要件を外す、rc=0 の誤り文を ok にする、切り替えの条件を反転させる)。欠陥ごとに 3 回中 2 回以上で指摘されること、欠陥を入れていない対照で誤指摘が 1 件以下であることを合格にする。

**M2. 子が `notes/reviews/` を読めるため、独立性が保証されない**
- 根拠: `run_claude` は `--tools Read,Grep,Glob` で、cwd をリポジトリ全体 (`cwd=ROOT`、`codex_review.py:592–593, 620`) にしている。`notes/reviews/` には先行の子の `*.claude.log` (結果の JSON を含む) や `*.claude.prompt.md` が残り、同じ名前の `…-plan-claude.claude.log` が既にある (Glob で確認)。V1 の「各 3 回 (独立の子)」は、先行の子の答えを読めれば成り立たない。§4.2 の「親の信念・記憶を持ち込まない」は自動読み込みを止めるだけで、Read による参照は止めていない。
- 対案: 子を `notes/reviews/` を含まないスナップショットで動かす (`git archive HEAD` + 作業ツリーの差分パッチ + 証拠一式 + 対象 plan を一時ディレクトリに置き、そこを cwd にする)。パス拒否ルール (`Read(notes/reviews/**)` 相当) は代替になりうるが、構文は未検証なので使うなら V3 で拒否が効くことを確かめる。

**M3. 記録ゲート `check_plans.py` は記録の実在だけを見る**
- 根拠: `check_plans.py:79–84` は `notes/reviews/*.md` の実在しか確かめない。`**未成立 (…)**` の記録も、判定が `?` の記録も、`[Claude 代替]` の記録も同じく PASS になる。計画 §3 が挙げる欠陥 (誤りの文が rc 0 で記録になる) は生成側を直しても、§4.3 の Copilot 経路のように手で表に貼る場合はゲートで止まらない。
- 対案: `check_plans` が参照先の中身を読む。先頭が `**未成立` なら FAIL。`判定: **GO|GO-with-changes|NO-GO**` と `指摘数:` の行が無ければ FAIL (diagnose は `結論:`)。`[Claude 代替]` は PASS にしつつ一覧に表示する。V5 の試験項目に入れる。

**M4. 記録の commit と実際に読まれた版がずれる (作業ツリーに未コミットの変更がある)**
- 根拠: 証拠一式の `head.txt` は `git rev-parse HEAD` だけを書く (`build_evidence`、`codex_review.py:455–460`)。子は作業ツリーを読むが、セッション開始時の gitStatus では `codex_review.py` と plan 類が変更中。計画の行番号も現在の版と合わない (計画 §3「352 行」「:233」、現行は 909 行、PATH 確認は 372 行)。記録の `commit: HEAD` は「どの版を読んだか」を保証しない。
- 対案: 証拠一式に `git status --porcelain` と対象ファイルの `git diff HEAD` を入れ、記録の見出しに `dirty: yes/no` と差分のハッシュを書く。前書きに「引用する行番号は証拠一式の版で確かめる」を足す。

**M5. 同一系統の代替だけで result 段を通せる**
- 根拠: §4.1 は result 段も「codex → Claude」とし、§4.3 は「弱い証拠」と書くだけ。check_plans は区別しないため (M3)、AGENTS.md:230 の accepted への移動が同一系統だけで通る。
- 対案 (ユーザの判断事項): result 段の代替は、ユーザの明示の承認行 (`result 代替承認: 日付・理由`) が無い限り status `done` の条件を満たさないものとする。plan 段と諮問は代替を許す。

**M6. 実装が計画より先に進み、順序も逆になっている**
- 根拠: 作業ツリーの `codex_review.py` には §4.2 の実装 (`--backend`、`build_evidence`、`classify_body`) が既にあり 909 行。計画は「352 行」「実装 O」のまま。V2〜V4 の偽 CLI 試験は見当たらない (テストの Glob は該当なし)。§5.1 は 2 (実装) → 2b (V1) の順だが、V1 が失敗すれば実装の大半は要らなくなる。
- 対案: §5.1 と status を実態に合わせて直す。順序は V1 (C1 の作り直し後) → V2 (偽 CLI) → 実装の残り → V3・V4。

**M7. 証拠一式が子の文脈に入り切らない恐れ**
- 根拠: `EVIDENCE_DIFF_MAX_BYTES = 512 * 1024` (`codex_review.py:71`)。日本語とコードの混在で、1 KiB を数百トークンと見ると 512 KiB は 100k トークン超になる (推定。未計測)。子の Read は既定で 2000 行に切られ、計画 §4.2 は「読み切ったか」を記録させていない。
- 対案: 上限を 128 KiB 程度に下げる (推定で数十k トークン)。子は判定の前に読んだファイルの一覧を記録に書き、読めていないファイルがあれば `証拠不足` を出す。

**M8. diagnose の代替は数値の判定を構造的にできない**
- 根拠: STAGE_DIAGNOSE は `res_*.h5` の統計量を python で出すよう求める (`codex_review.py:170`)。子は Bash を持たず、計画自身が「HDF5・判定ツールを走らせることはできない」と書く。前書きは「実行しないと得られない数値は既存のファイルにあるものだけ」と指示する (`claude_preface`、543–568 行)。よって数値が主題の諮問は `証拠不足` になりやすい。自動切り替えは unavailable しか見ないので、`incomplete` では代替が回らない。
- 対案: 親が統計量を事前に計算して証拠一式に入れる (`check_convergence.py` の VERDICT と残差の要約など)。計画と AGENTS.md:250–256 に「数値の諮問は代替では incomplete になりやすい」を明記する。

**M9. 依頼側との関係が既定では記録されず、完了条件 2 を満たさない**
- 根拠: 計画 §1 (2) は「同一系統であること」が記録に残ることを完了条件にしている。しかし `--requester-model` の既定は空 (`codex_review.py:798–799`) で、空なら記録は「不明」になる (`requester_relation`、676–682 行)。
- 対案: `--backend claude` のときは `--requester-model` を必須にする (無ければ終了)。または完了条件を「`--requester-model` を渡した実行の記録がある」に改める。

### Minor

- **m1.** `.prompt.md` と `.claude.prompt.md` が gitignore に無い (`.gitignore:16` は `*.log`、`:18` は `notes/reviews/*.evidence/` のみ)。96 kB 超のプロンプトで未追跡ファイルが増え、`git add -A` で入る。対案: `notes/reviews/*.prompt.md` を追加する。
- **m2.** claude の JSON を「`{` で始まる行」だけで探している (`codex_review.py:635–641`)。整形された複数行の JSON だと `data` が None になり、全試行が `failed` になる。対案: まず stdout 全体を `json.loads` し、だめなら行を探す。V4 の前に実出力の形を一度確かめる。
- **m3.** 判定の正規表現は本文の最初の一致を取る (`VERDICT_RE`、101 行)。「NO-GO の場合は…」が先に出ると誤読する。`判定:` の行に限る。また `INSUFFICIENT_NONE_RE` (104 行) は「該当なし」を含まないので、「証拠不足: 該当なし」が `incomplete` になる。
- **m4.** unavailable の型が広い (`forbidden`・`authenticat`、`codex_review.py:86–87`)。codex が rc≠0 で終わったときの誤判定で、報告されない切り替えが起きうる。`stream disconnected` (90 行) は一時的な通信断でも即座に claude へ切り替える。対案: 一時的な型 (stream・timeout) は codex を 1 回だけ再試行してから切り替え、型の語は狭める。
- **m5.** 予算が書かれていない。auto で codex が 1800 s でタイムアウトすると、claude も 1800 s で、前景実行は計 3600 s になる (`codex_review.py:802`)。claude 側には `--max-turns` も費用の上限も無い。AGENTS.md:287 の `run_in_background` の規則が切り替え後の経路にも効くことを書き、上限の有無は CLI で確かめて入れる。
- **m6.** V2(b) の unavailable の型に `auth` (401/403、`codex_review.py:86–87・96`) が無い。V3 は子の道具を「指示した」だけでは弱い。`--output-format stream-json` の init イベントで実際の道具の一覧を記録し、そこで確かめる。
- **m7.** rc=0 で本文が rate limit・usage 上限の文だと `incomplete` になり、切り替えない (`codex_review.py:441–442`)。計画 §4.1 は「usage 上限・rate limit」を unavailable に入れている。対案: 計画の文言を「rc≠0 のとき」に限るか、本文の誤りの文も切り替えの候補に加える。どちらにするかを決める。
- **m8.** 計画の引用行がずれている: §3 の `codex_review.py:233`・`check_plans.py:73` (現行は 372 行・79–84 行)。文書の同期のときに直す。

## 推奨

実装を進めたり規則 (§4.3) を採用したりする前に、**V1 を作り直す**。現行版に固定し、答えが一意に決まる問い (C1 の例) にして、既知の欠陥を仕込んだ検出試験 (M1) を加える。同時に子の読み取り範囲を `notes/reviews/` を外したスナップショットに閉じる (M2)。これが通るまで AGENTS.md・手順書・skill への反映は止める。

指摘数: Critical 1 / Major 9 / Minor 8
