あなたは forge (CUDA 圧縮性 CFD ソルバ) リポジトリの開発ワークフローの外部レビュアである。
コミット df82b46e で入れた「モデル分担とエスカレーション」の仕組みをレビューしてほしい。日本語で答えること。

## 背景
- 開発は Claude Code (モデル階層: 上位 Fable > 中位 Opus > 下位 Sonnet) が実装・run・後処理を行い、
  codex (あなた) が plan 段 / result 段の外部レビューを行う体制。
- 上位モデル Fable が高価なので「判断の場面だけ」に使いたい。これまでは既定モデルが Fable (1M コンテキスト) で
  全セッション・全ターンが上位モデルだった。費用は「コンテキスト長 × ターン数」で効くので、run 監視・ログ読み・
  後処理・README 同期といった長時間の機械的作業が主因という見立て (未計測。見立てである)。
- 分担の軸は「難しさ」ではなく「誤りをゲートが検知できるか」とした。収束 (check_convergence.py)・準定常
  (check_quasisteady.py)・メッシュ品質・plan 構造 (check_plans.py) はツールと Stop/PostToolUse フックが拾う。
  拾えないのは「もっともらしいが誤った真因」と「設計判断」。
- 仕組み: (1) `.claude/agents/` にモデル固定のサブエージェント 3 本 (run-watcher=sonnet / implementer=opus /
  diagnostician=fable・読み取り専用)。主セッションの Claude が description を見て自動で呼ぶか、ユーザが明示で呼ぶ。
  サブエージェントは主セッションの会話文脈を引き継がない (渡されたプロンプトだけ)。
  (2) AGENTS.md にエスカレーション条件 6 つを「外から観測できる行為」で明文化。
  (3) plan テンプレート §5.1 残作業表に担当列 (F/O)。
- サブエージェントを実際に呼んだ試験はまだ 1 回もしていない。

## 当方が自分で疑っている弱点 (これに縛られなくてよい)
(a) エスカレーション条件 3 (「予想と食い違った」) と 4 (「真因を書く前」) は、結局は中位モデルの自己認識に依存していて
    発火しないのではないか。より機械的に検知できる条件に置き換えられるか。
(b) run-watcher (サブエージェント) 経由で forge を起動したとき、PreToolUse フック hook_forge_guard.py の
    session_id による case の claim と、Stop フック hook_convergence_gate.py の検査対象が整合するか。
    (サブエージェントのツール呼び出しのフック入力が親と同じ session_id を持つかは未確認。持たない場合に何が壊れるか)
(c) diagnostician は文脈を引き継がないので、呼び出し側 (中位モデル) が書くブリーフの質が律速になる。
    中位モデルが既に誤った枠組みで症状を要約して渡すと、上位モデルもその枠内で答えてしまう。防ぐ形式はあるか。
(d) 既存ルール (run ディレクトリ複製・case README 同期・codex レビュー 2 回・並行セッションの plan 不可侵・
    commit 運用) との矛盾や二重規定、抜け穴 (例: サブエージェントが README を編集することと並行セッションの衝突)。
(e) 長時間 run (数十分〜数時間) を sonnet サブエージェントに持たせる設計は妥当か。サブエージェントが run 完了まで
    待てずに返る/判定を省く、といった失敗形は無いか。
(f) 費用対効果: この構成で本当に上位モデルの消費が減るか。逆に増える経路 (diagnostician が毎回リポジトリを
    読み直す、AGENTS.md+メモリ索引の固定コストが各サブエージェントに乗る等) は無いか。

## 禁止事項 (厳守)
- ファイルを変更しない (read-only サンドボックスで動いている)。
- **`*.log`, `residual_history.csv`, `res_*.h5`, `*.vtu`, `plans/README.md` を読まない** (巨大)。
- 下に列挙した `sed -n 'A,Bp' <file>` 以外のファイル読みをしない。grep は可。
- 両論併記で逃げず、**推奨は 1 つに絞る**。根拠は `ファイル:行` で示す。
- Claude Code のサブエージェント/フックの仕様について確信が無い点は「仕様未確認」と明記し、推測で断定しない。

## 読んでよいもの
- sed -n '229,271p' AGENTS.md                 (新設した節)
- sed -n '44,64p' AGENTS.md                   (計算・実行ルールの冒頭: run 複製・索引)
- sed -n '76,98p' AGENTS.md                   (収束確認: 段階起動の区間・Stop フックのセッション別スコープ)
- sed -n '182,228p' AGENTS.md                 (決定トリガ・codex レビュー節)
- sed -n '1,52p' .claude/agents/run-watcher.md
- sed -n '1,46p' .claude/agents/implementer.md
- sed -n '1,55p' .claude/agents/diagnostician.md
- sed -n '44,62p' plans/_template.md
- sed -n '1,38p' .claude/settings.json
- sed -n '1,45p' solver_density_cuda/tools/hook_forge_guard.py
- sed -n '1,91p' solver_density_cuda/tools/forge_gate_claims.py
- sed -n '1,63p' solver_density_cuda/tools/hook_convergence_gate.py

## 出力してほしいもの
1. 判定: GO / GO-with-changes / NO-GO
2. 指摘を Critical / Major / Minor に分けて列挙。各指摘に根拠 (`ファイル:行`) と、具体的な修正案 (どのファイルのどの文をどう変えるか)。
3. 上の (a)〜(f) それぞれに 1〜3 行で答える。
4. 最後に「最初に直すべき 1 点」を 1 つだけ挙げる。
