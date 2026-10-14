# codex レビュー: tooling-review-fallback-backend (plan)

- **plan**: [`plans/active/tooling-review-fallback-backend.md`](../../plans/active/tooling-review-fallback-backend.md)
- **stage**: `plan`
- **date**: 2026-10-14
- **commit**: `88675245` (feature/nozzle-cooled-wall-settling)
- **codex**: effort `high`, 5.0 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M4/m1

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
代替の判定役を用意する目的は妥当で、既存実装では未解決です。  
ただし、文脈の分離、証拠へのアクセス、レビュー成立条件が未定義のため、現状のまま実装には進めません。

目的・重複については、現行コードが codex 不在で停止することを確認しました（[codex_review.py:233](/home/sano/work/forge-deliv/solver_density_cuda/tools/codex_review.py:233)）。`plans/README.md` と `accepted/` の照合でも同じ代替機能は確認できず、関連する [tooling-rerun-conditions.md:25](/home/sano/work/forge-deliv/plans/accepted/tooling-rerun-conditions.md:25) は計算入力の準備が対象です。

数値カーネルや解析設定を変更しないため、node/cell、float32、block-DPLUR、周期・軸対称への直接の影響はありません。本件の検証に新しい CFD 計算や段階起動を追加する必要はありません。ただし、**代替レビューが既存の数値検証を点検できること**は必要です。

1. **Major — 読み取り専用ツールだけでは、result・数値診断の要求を満たせない。**

   **根拠:** [plan:51](/home/sano/work/forge-deliv/plans/active/tooling-review-fallback-backend.md:51) は `Read,Grep,Glob` のみに制限し、同じプロンプトを渡します。一方、既存プロンプトは result で `git diff` の取得を要求し（[codex_review.py:73](/home/sano/work/forge-deliv/solver_density_cuda/tools/codex_review.py:73)）、diagnose では Python による HDF5 の統計抽出を想定しています（[同:99](/home/sano/work/forge-deliv/solver_density_cuda/tools/codex_review.py:99)）。`sed` を `Read` に置き換えるだけでは解消しません。[plan:62](/home/sano/work/forge-deliv/plans/active/tooling-review-fallback-backend.md:62) の新規チャットにも、リポジトリへアクセスできる保証がありません。

   **対案:** 親のラッパーが機械的に作る**証拠一式**を入力契約にしてください。対象版・diff・必要なコード・既存判定ツールの出力を、run パス、実行引数、判定区間、入力ファイルのハッシュとともに渡します。既存の `--extra` を再利用できます。追加証拠が必要なら要求を返し、証拠不足のままレビュー完了にしない仕様にします。Claude 向けプロンプトも、その能力に合わせて変更が必要です。

2. **Major — 提案の起動引数では「文脈を共有しない」「コマンドを実行できない」を保証できない。**

   **根拠:** [plan:20](/home/sano/work/forge-deliv/plans/active/tooling-review-fallback-backend.md:20) の独立性の前提に対し、起動仕様はセッション保存と MCP の制限だけです。実際のリポジトリには、コマンドを起動する `Stop` フックがあります（[.claude/settings.json:14](/home/sano/work/forge-deliv/.claude/settings.json:14)）。これはモデルに公開するツールの制限とは別経路です。

   また、`--no-session-persistence` は会話の保存を止める機能です。外部シェルや Copilot 経由の起動も対象にするなら、別セッションで共有される auto memory の扱いも明示する必要があります。[公式 CLI 仕様](https://code.claude.com/docs/en/cli-reference)、[公式メモリ仕様](https://code.claude.com/docs/en/memory)。

   **対案:** 自動メモリ、フック、プラグイン、ユーザー由来の追加指示を分離する起動仕様を先に決めてください。対応版の `--safe-mode` 等を検証し、必要なリポジトリ規則は明示的に入力します。管理ポリシーのフックなど隔離できないものが残る場合も検出・記録します。「技術的に書き込み禁止」とするなら、対象リポジトリの OS レベルの読み取り専用化まで条件に含めるべきです。

3. **Major — 「本文が取れた」と「レビューが成立した」を区別していない。**

   **根拠:** [plan:49](/home/sano/work/forge-deliv/plans/active/tooling-review-fallback-backend.md:49) は最終メッセージの有無を切り替え基準にしています。現行実装は本文が非空なら記録を作り、判定・指摘数を抽出できなくても `?` を埋め、`rc=0` なら成功を返します（[codex_review.py:295](/home/sano/work/forge-deliv/solver_density_cuda/tools/codex_review.py:295)、[同:323](/home/sano/work/forge-deliv/solver_density_cuda/tools/codex_review.py:323)、[同:348](/home/sano/work/forge-deliv/solver_density_cuda/tools/codex_review.py:348)）。

   **実行確認:** ファイル出力と子プロセスをメモリ上のモックに置き換え、現行 `main()` に以下の本文を与えました。

   - `Error: rate limit reached`
   - `レビューに必要なファイルを読み取れません。`

   どちらも **終了コード 0、レビュー md の書き込み経路に到達、`判定=? 指摘=C?/M?/m?`** でした。さらに、実在記録を指す `?` のレビュー行を `check_reviews()` に渡すと合格しました。これは [check_plans.py:73](/home/sano/work/forge-deliv/solver_density_cuda/tools/check_plans.py:73) が内容を検査しない仕様どおりです。

   **対案:** ラッパー側で「正常なレビュー」「利用不能」「実行失敗」「不完全な回答」を区別してください。正常終了に加え、stage ごとの必須出力を満たした場合だけ完了行を発行します。自動切り替えは定義した利用不能に限定し、**有効な `NO-GO` は切り替え理由にしない**こと。失敗した各試行のログも保持します。`check_plans.py` の構造 lint をレビュー成立の保証に使ってはいけません。

4. **Major — V1〜V4 では、書き込み禁止と主要な失敗経路を検証できない。**

   **根拠:** [plan:83](/home/sano/work/forge-deliv/plans/active/tooling-review-fallback-backend.md:83) の `git status` 前後一致では、既に変更済みのファイルへの追加変更や、無視対象ファイルの変更を検出できません。実際に [ .gitignore:65](/home/sano/work/forge-deliv/.gitignore:65) は `case/**/run_*/` を無視しています。反対に、期待するレビュー md の新規作成自体が status を変える場合もあります。

   また、V3 は plan 段だけで、result 経路、タイムアウト、両 backend の失敗、非空のエラー本文、長大プロンプトを検証しません。

   **対案:** 小さな隔離 fixture を使い、親が生成する許可済み成果物を除いた**ファイル集合と内容ハッシュの一致**で確認してください。§6 には、上記の失敗条件に加え、codex 成功時の Claude 呼び出し回数 0、切り替え対象時 1、有効な `NO-GO` 時 0、96 KiB 超の入力、result の証拠不足を事前登録します。ライブモデルの試験は疎通確認に使い、分岐と拒否条件は偽 CLI で再現可能にします。

5. **Minor — 記録の定型文が、選択理由とモデルの関係を事実以上に断定する。**

   **根拠:** [plan:50](/home/sano/work/forge-deliv/plans/active/tooling-review-fallback-backend.md:50) は明示的な `--backend claude` を許しますが、[plan:55](/home/sano/work/forge-deliv/plans/active/tooling-review-fallback-backend.md:55) は一律に「codex・Fable が使えないため」「同一系統」と記録します。明示選択時に利用不能を確認したとは限らず、依頼側が別モデルなら同一系統とも限りません。

   **対案:** backend、要求モデル、取得できれば実効モデル、CLI 版、選択理由を分けて記録してください。選択理由は「明示指定／スイッチ／自動切り替え」、モデルの関係は「同一／別／不明」とし、不明を推測で埋めないことです。

**推奨は、証拠一式を入力とする、環境を分離した読み取り専用 Claude 経路の追加です。** 実装前に #1 → #2 → #3 の入力・起動・完了条件を §4 に確定し、#4 を §6、#5 を記録仕様に反映してください。その後にラッパー実装、規則・手順の同期、試験、result レビューへ進む順序を推奨します。

確認はコード閲覧、CLI の版・help、書き込みを伴わないモック実行で行いました。ファイル変更・CFD 計算・Claude への実レビュー投入は行っていません。

指摘数: Critical 0 / Major 4 / Minor 1
