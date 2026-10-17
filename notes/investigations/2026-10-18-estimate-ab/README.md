# 段の見積もりの判別の A/B (2026-10-18)

plan [`tooling-design-orchestrator`](../../../plans/active/tooling-design-orchestrator.md) §4.5.1「段の見積もりの契約」・§6「V3 の 13a の事前登録」の判別の A/B
(codex 諮問 [v3](../../reviews/2026-10-18-orchestrator-v3-diagnose-v3.md))。GPU・forge なし。

- 台本 [`estimate_ab.py`](estimate_ab.py)、結果 [`result.json`](result.json)。
- 回した版: cc5b9d18 の作業ツリーに、見積もりの受け渡し (campaign の `estimates` → 台帳に固定 → recipe の段の `estimate`) だけを足した状態。
  **「GPU の段の見積もりが無い campaign は読み込みで止める」を足す前**に回した (足した後は A の腕が読み込みで止まるので、この台本は記録として残す)。
- 変えたのは Euler の段の見積もりの `gpu_hours` だけ (A: 無し、B: 4/60)。承認の予算 GPU 0.01 時間、格子の段 (mesh_p0) が済んだ台帳、本物の driver・recipe、
  起動は到達の記録だけの代わり。

| 腕 | Euler の試行の追加 | 起動への到達 | 止まり方 |
| --- | --- | --- | --- |
| A (見積もり無し) | 1 | 1 | 起動に到達 (代わりで止めた) |
| B (4/60 GPU 時間) | 0 | 0 | `stopped`: 見積もり 0.0667 が残りの予算 0.01 を超える — 次の段を始めない |

事前登録どおり (A が到達・B が予算で止まる)。見積もりの未接続が「予算が足りない依頼でも Euler を起動しようとする」原因だったことを支持する。
B の値は試験の入力で、生産の見積もりの妥当性は示さない。
