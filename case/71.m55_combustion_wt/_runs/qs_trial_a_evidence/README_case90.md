# 90.qs_trial_a — case/71 (M5.5 燃焼ガス風洞ノズル) の入口条件と寸法を変えたノズル (手順書の試走 A)

手順: procedures/design-chain-quickstart.md。要求は card_p0.yaml (係数を入れたカードは card_p1.yaml)。AWS B (g5.xlarge・A10G・CUDA 13.2・float のビルド) で、HOME を
`/home/ubuntu/qs_trial` にした素の clone から打った (2026-10-16)。記録の写しは forge リポジトリの `case/71.m55_combustion_wt/_runs/qs_trial_a_evidence/`。

## 計算 run 一覧

| run | 目的・主要設定差分 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `run_0001_euler_p0` | Euler (壁 p0 = CONTUR の式のまま、problem_p0_euler.yaml) | 出口の軸の M 5.50017、`dM_max_rel_Md` 0.0006、非有限 0、121 s。`res_12000.h5` が NS p0 の Euler の参照 | ref |
| `run_0002_ns_timing` | NS の時間の試し (壁 p0、本段 2,000 step、停止規則なし) | 本段 15.4 ms/step、前の段 + 本段 2,000 step で 1.6 分 | 破棄予定 |
| `run_0003_ns_p0` | NS 1 回目 (壁 p0、停止規則 v2、プリセット wt_nozzle_v3) | 停止規則 v2 REACHED (本段 42,060 step)、δ_E の抽出 ok、関門「進む」、`check_settled` SETTLED (判定の区間 32,000〜42,000)、非有限 0、NS 13 分。`contur_calibrate/` の 3 係数 c₀ 0.688・c₁ 0.0215 /m・δ₀ −0.31 mm | ref |
| `run_0004_euler_p1` | Euler (壁 p1 = run_0003 で較正した 3 係数) | 出口の軸の M 5.50039、非有限 0、116 s | ref |
| `run_0005_ns_p1` | NS 2 回目 (壁 p1、初期場 run_0003、M_knot 2.70) | 停止規則 v2 REACHED (本段 30,088 step)、関門「進む」、SETTLED (判定の区間 20,000〜30,000)、非有限 0、NS 10 分。`wt_qualify_tol01/`: accepted (V5b 5.503596 = +0.065 %・V5c 片振幅 0.061 %・最大偏差 0.065 %、公差 0.1 %)。成果物 `deliverables/` PASS (寸法の関門 PASS・`thermal_use` PASS)。旧の公差 (0.05 %) の判定は `wt_qualify/`・`wt_qualify_final/` (escalate) | ref |
