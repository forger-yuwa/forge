# 諮問ブリーフ: plan 段レビューの採否の確認 (plan tooling-review-fallback-backend の改訂後 §4.2・§6)

AGENTS.md 条件 5 (codex の Major の採否)。判断役のスイッチが codex。採否は「全件採用」とし、改訂を plan に書いた (§3 の実測、§4.2、§6 V1〜V5、§6.1)。

## 観測事実
- plan 段レビュー `notes/reviews/2026-10-14-tooling-review-fallback-backend-plan.md` (M1 証拠一式、M2 文脈の分離、M3 結果の分類、M4 偽の CLI とハッシュの試験、m5 記録)。
- 諮問 `notes/reviews/2026-10-14-review-fallback-backend-diagnose.md` (同調の対照試験)。
- 親の実測 (plan §3): `claude -p --safe-mode` の子は CLAUDE.md・MEMORY.md を読まない (入力 約 4.7k トークン)、`--allowedTools "Bash(git log:*)"` 等の外の Bash は拒否されて
  `permission_denials` に残る、Write・Edit は無い、`modelUsage` に実効のモデル名。`--bare` は API キー認証だけ。

## 問い
**改訂後の §4.2・§6 に、採用した指摘の反映として残る最大の穴を 1 つ挙げ、防ぐ変更を 1 つ推奨してほしい** (無ければ「無し」と言ってよい)。
候補 (縛られなくてよい): (a) `--permission-mode default` で -p のとき、許可外の Bash が「拒否」になることの保証 (版で変わる)、(b) `git show`/`git diff` で親の未 commit の作業ツリーが見えることは
入力の契約と矛盾しないか、(c) 結果の分類の正規表現が codex の既存の記録の書式で誤って incomplete になる (後方互換)、(d) V1 の合格条件 6/6 の厳しさ。

## 禁止事項 (厳守)
- ファイルを変更しない (read-only サンドボックスで動いている)。
- **`*.log`, `residual_history*.csv`, `res_*.h5`, `*.vtu`, `*.msh`, `plans/README.md` を読まない** (巨大)。
- 下に列挙した `sed -n 'A,Bp' <file>` 以外のファイル読みをしない。grep は可。
- 両論併記で逃げず、**推奨は 1 つに絞る**。根拠は `ファイル:行` か提示した数値で示す。

## 読んでよいもの
- `sed -n '1,140p' plans/active/tooling-review-fallback-backend.md`
- `sed -n '280,352p' solver_density_cuda/tools/codex_review.py`
- `grep -n "結論\|指摘数\|GO-with-changes" notes/reviews/2026-10-1*-*.md | head -30`
