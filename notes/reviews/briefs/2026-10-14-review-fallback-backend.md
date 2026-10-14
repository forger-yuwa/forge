# 諮問ブリーフ: 判定役の 3 段目 (Claude の別インスタンス) の設計 (plan tooling-review-fallback-backend §4・§6)

AGENTS.md 条件 1 (plan §4・§6 を新規に書く)。判断役のスイッチが codex。

## 観測事実
- ユーザ (2026-10-14): 別の環境では codex も Fable も使えず、諮問・レビューができない (規則上、諮るまで先へ進めない)。「同じ言語モデルになるかもしれん」が、測る相手を作ってほしい。
- `solver_density_cuda/tools/codex_review.py` は codex exec を read-only で起動して記録を残す。codex が無ければ終了。
- `claude` CLI: `-p`・`--model`・`--effort`・`--tools Read,Grep,Glob` (組み込みの道具を読むだけに絞る)・`--no-session-persistence`・`--strict-mcp-config` がある。
- `check_plans.py` は記録ファイルの実在だけを見る。

## 方針 (plan §4)
- 順序: 諮問 Fable → codex → Claude の別インスタンス、レビュー codex → Claude。`--backend auto` は codex の失敗 (無い・最終メッセージが取れない) で claude に回し直す。スイッチ値 `claude` を追加。
- claude は読むだけの道具で、codex と同じプロンプト + 前書き (道具の読み替え、同調せず反証を優先)。記録に「同一系統の代替」。
- 主セッションで Fable が使えないときに `Agent` の `model` 上書きで `diagnostician` を呼ぶのは禁止し、`codex_review.py --backend claude` に統一 (記録が残るため)。

## 問い
**この設計のうち、実装後に「代替の判定が判定として機能しない」(依頼側と同じ誤りに同調する、記録が誤って本物の外部レビューに見える、読むだけのはずが書ける、など) 形で失敗する可能性が最も高い点を 1 つ挙げ、防ぐ変更を 1 つ推奨してほしい。**
候補 (縛られなくてよい): (a) 同一モデルの同調を前書きだけで抑えられるか (プロンプトの順序・依頼側の仮説を隠すなど)、(b) `--tools` の制限は本当に書き込みを防ぐか (CLAUDE.md・hook の読み込み、`--strict-mcp-config`)、
(c) auto の切り替えが codex の一時的な失敗でも黙って弱い判定に落ちる、(d) Agent の model 上書きを禁止する理由の妥当性。

## 禁止事項 (厳守)
- ファイルを変更しない (read-only サンドボックスで動いている)。
- **`*.log`, `residual_history*.csv`, `res_*.h5`, `*.vtu`, `*.msh`, `plans/README.md` を読まない** (巨大)。
- 下に列挙した `sed -n 'A,Bp' <file>` 以外のファイル読みをしない。grep は可。
- 両論併記で逃げず、**推奨は 1 つに絞る**。根拠は `ファイル:行` か提示した数値で示す。

## 読んでよいもの
- `sed -n '1,110p' plans/active/tooling-review-fallback-backend.md`
- `sed -n '1,60p' solver_density_cuda/tools/codex_review.py` / `sed -n '150,300p' solver_density_cuda/tools/codex_review.py`
- `sed -n '244,300p' AGENTS.md`
