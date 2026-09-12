# case/20 naca_ml 001.test

既定の検証ケース (procedures/verification/20-naca-ml.md)。run ディレクトリは `run_*` (2026-09-12 以降は本表で索引)。

## 計算 run 一覧

| run | 目的・主要設定差分 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `run_perf_regress_cell_slau_explicit` / `run_perf_regress_cell_tpair_explicit` (ローカル, 2026-09-12) | 高速化ブランチ `feature/perf-3d-speedup` の回帰: run_slau (cell, CPG 層流, 陽解法 RK3, res_4000) / run_tp_air (cell, TP 空気, 陽解法, res_2000) から 300 step 継続、基準 ×2 vs 最終 (`solver_density_cuda/tools/perf_regress.py`) | 絶対基準 or 2×ノイズ床で全場 PASS (TP は基準 4 本、ρ がノイズの 2.0 倍 = 2.2e-5)。2D 平面で面外速度 Uz が 4e-10→3e-6 m/s (主流の 1e-8, float 丸め) | active (回帰) |
