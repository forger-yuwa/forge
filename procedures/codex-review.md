# codex 外部レビューの手順 (計画立案時・検証結果時)

ルールの正本は [`AGENTS.md`](../AGENTS.md) 「codex レビュー (計画立案時と検証結果時の 2 回)」。本書はその手順。

## 1. いつ回すか

| 段階 | タイミング | 何を見てもらうか |
| --- | --- | --- |
| `plan` | `plans/active/<plan>.md` の §4 設計方針と §6 検証計画が書けた時点、**実装着手前** | 課題の同定・方針の健全性・既存構造との整合・検証計画の定量性・見落としリスク。判定 GO / GO-with-changes / NO-GO |
| `result` | 検証 run の VERDICT (`check_convergence.py` / `check_quasisteady.py`) が出そろい、`status: done` にして `accepted/` へ**移す前** | 実装 diff (`git diff <base>...HEAD`) と plan の主張の突き合わせ、run 実測との整合、回帰、文書整合。判定 GO / GO-with-changes / NO-GO |

対象は [`AGENTS.md`](../AGENTS.md) 「開発フロー」を踏む変更 (新規機能・スキーム/設計方針の変更)。例外 (typo、1 ファイル内バグ修正、
振る舞い同一のリファクタ、docs のみ) は不要。判断がつかなければ回す側を選ぶ。

## 2. 回し方

```bash
# plan 段 (実装前)
python3 solver_density_cuda/tools/codex_review.py plans/active/<plan>.md --stage plan

# result 段 (done にする前)。--base は diff の基準 ref (既定 main。feature ブランチを積んでいるなら分岐元)
python3 solver_density_cuda/tools/codex_review.py plans/active/<plan>.md --stage result --base main

# 補助: 重点と参考ファイル (case README 等) を渡す / プロンプトだけ確認
python3 solver_density_cuda/tools/codex_review.py PLAN --stage plan --focus "§4.2 に集中" --extra case/46.sern_design/README.md
python3 solver_density_cuda/tools/codex_review.py PLAN --stage plan --dry-run
```

- 所要 5〜15 分。Claude Code からは **`run_in_background` + timeout 1200 s 以上**で呼ぶ。進捗は
  `notes/reviews/<名前>.log` のサイズ、または `~/.codex/sessions/YYYY/MM/DD/rollout-*.jsonl` の末尾 timestamp で見る。
  フリーズ疑いでも kill する前に `ps -o pid,etime,time` の CPU 時間で進行中か確かめる。
- ツールが守っている作法 (手で `codex exec` を打つときも同じ): `--sandbox read-only`、**stdin は `/dev/null`**
  (開いたままだと EOF 待ちで永久ブロック)、**stdout をパイプに通さない** (tail/head を挟むと終了まで何も見えない)、
  `model_reasoning_effort` を `high` に上書き (`~/.codex/config.toml` は `low`)。
- 出力: `notes/reviews/<日付>-<plan stem>-<stage>.md` (メタ情報ヘッダ + codex の最終メッセージ)。同名があれば `-2` を付ける。
  生ログ `*.log` は git 追跡外。

## 3. 記録と採否

1. ツール末尾に表示される行を plan の **§6.1 レビュー記録**の表に貼る (列: 段階 / 日付 / 記録 / 判定・指摘数 / 対応)。
2. **Critical / Major は必ず採否を決める**: 採用なら §5.1 残作業表に入れる (優先順つき)、却下なら §6.1 の対応欄に理由を書く。
   Minor は対応欄にまとめて可。
3. 指摘を鵜呑みにしない。根拠 (`ファイル:行` / run の数値) が示されていない指摘、または示された根拠が実際と違う指摘は
   再検証してから採否を決め、その結果も対応欄に書く (codex の誤指摘も記録として価値がある)。
4. レビューを回した/反映した**応答には記録ファイルのパスと採否を書く** (run パス明示・VERDICT 貼付と同じ扱い)。
5. NO-GO / GO-with-changes で方針を変えたら、その場で plan §4 と残作業表を直す (「決定」がトリガ。実装を待たない)。

## 4. 強制の仕組み

- `solver_density_cuda/tools/check_plans.py` が各 `plans/active/*.md` の「レビュー記録」節を見る:
  `status` が `in_progress` 以上なら `plan` 行、`done` なら `result` 行が必要。行の記録欄の `notes/reviews/*.md` は実在すること。
- 免除は段階欄を `plan 免除` / `result 免除` にし、対応欄に理由を書く。**2026-09-09 以前に起票済みで実装が進んでいる plan は
  `plan 免除` で通す** (以後の `result` 段は免除しない)。
- plan を編集すると PostToolUse フック (`hook_plan_todo_gate.py`) が編集したファイルだけ lint し、不足を返す
  (既存 plan を一斉に FAIL にはしない「触ったら直す」方式)。
- 別セッションが同じツリーで並行作業しているときは、相手の plan にレビュー行を書き込まない。

## 5. 過去のレビュー (ツール化前)

ツール化 (2026-09-09) 以前は `notes/sessions/*-prompt-codex.md` にプロンプトを置き、結果は plan の変更ログや
`notes/investigations/` に直接書いていた (例: [`axisym-node-review-prompt-codex.md`](../notes/sessions/axisym-node-review-prompt-codex.md)、
[`nozzle-deltastar-throat-review.md`](../notes/investigations/nozzle-deltastar-throat-review.md))。以後は本手順に統一する。
