# case/23 axi_nozzle

## 計算 run 一覧

| run | 目的・主要設定差分 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `run_0100_perf_regress_cell_axisym_tp_implicit` (ローカル, 2026-09-12) | 高速化ブランチの回帰: run_tp_air_axisym_implicit (cell, 軸対称, TP 空気, 陰解法) の res_600 から 300 step 継続、基準 ×2 vs 最終 | 全場 (P/T/ρ/U/ρe/μt) の差 ≤3.4e-6 で base×base (≤3.2e-6) と同等。config の旧キー `LESorRANS/LESmodel` は `model: "wale"` に置換して実行 (WALE double 復帰後の最終版) | active (回帰) |
