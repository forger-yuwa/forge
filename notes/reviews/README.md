# codex レビュー出力 (`notes/reviews/`)

`solver_density_cuda/tools/codex_review.py` の出力置き場。ファイル名は `<日付>-<plan stem>-<stage>.md` (`stage` = `plan` / `result`)。
本文は codex の最終メッセージの転記で、**採否と対応は各 plan の §6.1 レビュー記録に書く** (ここは編集しない)。
手順は [`procedures/codex-review.md`](../../procedures/codex-review.md)。`*.log` (生ログ) は git 追跡外。

| レビュー | plan | stage | 判定 |
| --- | --- | --- | --- |
| [2026-09-09-tooling-nozzle-sern-chain-plan.md](2026-09-09-tooling-nozzle-sern-chain-plan.md) | `tooling-nozzle-sern-chain` | plan | NO-GO, C2/M7/m0 |
| [2026-09-19-tooling-divergence-triage-minimal-plan.md](2026-09-19-tooling-divergence-triage-minimal-plan.md) | `tooling-divergence-triage-minimal` (草案 v1) | plan | NO-GO, C0/M11/m1 |
| [2026-09-19-tooling-divergence-triage-minimal-plan-2.md](2026-09-19-tooling-divergence-triage-minimal-plan-2.md) | `tooling-divergence-triage-minimal` (草案 v2) | plan | NO-GO, C0/M8/m2 |
| [2026-09-22-codex-model-tiering.md](2026-09-22-codex-model-tiering.md) | 無し (`AGENTS.md`「モデル分担とエスカレーション」。採否は記録本文) | 自由形式 | GO-with-changes, C0/M6/m2 |
