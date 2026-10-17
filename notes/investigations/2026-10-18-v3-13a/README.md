# §6 V3 の 13a (依頼の課題集の機械の側) の 1 回目の結果 (2026-10-18)

plan [`tooling-design-orchestrator`](../../../plans/active/tooling-design-orchestrator.md) §6「V3 の 13a の事前登録」(cc5b9d18)。
試験 `design/tests/run_v3_requests_tests.py` (b1a02511)、依頼の課題集 `design/tests/fixtures/v3_requests/`。手元、GPU・forge なし。
出力の全文 [`run_v3_requests_tests.log`](run_v3_requests_tests.log)。

| 依頼・項目 | 結果 |
| --- | --- |
| (K) 安価な負例 15 (11 種 + 13a で足した見積もりの欠け) | 全部実行されて ok (skip・欠け 0) |
| (K) 正例 4 | 全部 ok |
| N1 case/71 | 承認の前の生成 → 承認 → 生成し直して `approval.check` OK → `hand_mesh.py build --approval` → `check` で NS・Euler の品質 PASS |
| N2 試走 A | **FAIL**: `hand_mesh.py build --approval` が「実効の節点数 412160 (ni 4480 × nj 92) が上限 412000 を超える」で問題 YAML を書かずに止まった (当時の格子は ni 4470 × 92 = 411,240) |
| N3 case/70 | `mesh_stage.check_template` で NS・Euler の品質 PASS |
| S1 公差が届かない | iteration 1/3 → `back_to_correction` (V5b FAIL M 5.519563・V5c FAIL 0.404 %、解析の成立 PASS)、iteration 3/3 → `rejected` |
| S2 ユーザ欄の欠け | preview に `wall.thermal: 出所が無い`、承認の記録を書かない |
| S3 予算が足りない | Euler の前に予算で止まる (試行の追加 0・起動 0・格子の試行と止まった記録は残る)、見積もりの無い campaign は読み込みで止まる |

13a の合格は V3 の完了ではない (13b の後)。N2 の扱いは plan §6 に記録する。
