# 諮問: 事前登録した V2c (b)・(c)・(g) の「生成の関数の呼び出し 0 回」が FAIL (2026-10-18)

plan `plans/active/tooling-design-orchestrator.md` §6 の「V2c の事前登録 (2026-10-18、#8 の実装の前)」と §4.6「照合の場所と順」。
前の諮問: `notes/reviews/2026-10-18-orchestrator-mesh-limits-diagnose-mesh-limits.md`。

聞きたいこと: この FAIL の扱い。(A) 実装を直して事前登録のまま全部を回し直す、(B) 他の扱い。基準は緩めない前提で、(A) の直し方
(下の案) に穴が無いか。

## 1. 観測事実

試験 `design/tests/run_mesh_limits_tests.py` (未 commit、事前登録 af510413 の後に書いた)。`hand_mesh.build` を同じプロセスで呼び、
`wt_card.estimate_mesh` と `mesh2d.generate_axisym_mesh` を数える (モジュールの属性を差し替え)。

| 項目 | 結果 |
| --- | --- |
| V2a (変更の前後の A/B) | 全部 PASS: 問題 YAML の mesh・mesh_euler (max_nodes_reason 以外)・Euler の写し・`check` の NS と Euler の nozzle.msh (バイト)・NS の nozzle.h5 の MESH/*・BCONDS/* (91 データセット) が一致。参考: case/71 の problem_m55_p3.yaml (10-15) との差は `mesh.wall_first_frac_table`・`mesh.x_density_table`・理由 |
| (a) overrides の max_nodes | PASS (カードの検査で誤り) |
| (b) mesh_limits なし (既定 30 万) の承認 | 問題 YAML を書かずに止まった (`check_node_budget`: mesh 442520 > 300000)。**ただし `generate_axisym_mesh` 1 回** → 事前登録の文言で FAIL |
| (c) mesh_limits {400000, 467000} | 同じ形 (442520 > 400000、`generate_axisym_mesh` 1 回) → FAIL |
| (d) approval.check | PASS (MISMATCH) |
| (e) case/70 の承認を case/71 のカードに | PASS (`estimate_mesh` 0 回、ApprovalError) |
| (f) 判別の A/B (Pt だけ) | PASS: A は estimate_mesh に到達 (1 回)、B は 0 回・ApprovalError (`requirements.conditions.Pt が承認の値と違う`) |
| 正例 (係数だけ変えた試行) | PASS (`check_card` OK、kind trial) |
| (g) Euler だけ超過 (axis_gap_frac 0.005) | 書かずに止まった (mesh_euler 966810 > 467000)。**`generate_axisym_mesh` 1 回** → FAIL |
| (h) | 試験の途中で recipe の import の誤り (`from ..intake` が forge_design.campaign.intake を指す、#6c からの潜在の誤り; 同じ形 7 か所を直した) で落ちた。直して回し直し中 |

`generate_axisym_mesh` の 1 回の出所: `hand_mesh.build` の先頭の「式の写しの照合」が呼ぶ `wt_card.estimate_mesh` の中の `_mesh_ar_full`
(`design/forge_design/intake/wt_card.py:1107`)。これは**入力の問題 YAML (生成器の出力、mesh.ni 3000 × nj 97 = 291,000 点) の雛形の格子**で、
`estimate_mesh` 自身が作る前に入力の上限 (既定 30 万) で `check_node_budget` を通している (`wt_card.py:1096` の `guard_nodes`)。上限を超えた
手組みの格子 (44.3 万点・Euler 96.7 万点) は作っていない。`hand_mesh.build` が `estimate_mesh` から使うのは断面の値 (`sections`、1 次元の式) だけで、
全格子の結果 (AR・ni の繰り上げ) は使わない (`hand_mesh.py:212〜221`)。手組みの格子の節点数は、第一層の表と密度の表から ni を決めるまで分からない。

## 2. 期待値と出典

- 事前登録 (plan §6 V2c、af510413): 「負例は『生成の関数の呼び出し 0 回・出力の問題 YAML 0 件』で止まることを見る (最後に例外が出ただけを合格にしない)」。
  (b)・(c)・(g) にも適用すると書いた。
- 前の諮問 (mesh-limits) の要求: 承認の照合 (M1) は最初の全格子の生成の前。節点数 (M2) は共通の `check_node_budget` で YAML を書く前。
  「生成の関数の呼び出し 0 回」は前の諮問の判別の A/B ((f)) の文言で、(b)・(c)・(g) への一般化は事前登録を書いた私の拡張。

## 3. 案

(A) 実装を直す: `wt_card.estimate_mesh` に `full_grid: bool = True` を足し、`hand_mesh.build` は `full_grid=False` で呼ぶ (断面の値と
式の写しの照合だけ。全格子の AR の測定と ni の繰り上げをしない)。`hand_mesh` の出力は変わらない見込み (全格子の結果を使わないので)。
事前登録は変えずに V2a〜V2c を全部回し直す (V2a の B も新しいコードで)。生成器の既定 (`full_grid=True`) は変えない。

(B) FAIL として記録し、(b)・(c)・(g) の条件を「上限を超える格子の生成 0 回」と読み直す — これは結果を見た後の基準の読み替えなので採らない予定。
