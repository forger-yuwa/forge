# case/70 — M5 空気の風洞ノズルのデモ (要求カード → 問題 → 成果物一式)

plan [`tooling-nozzle-wt-deliverables`](../../plans/active/tooling-nozzle-wt-deliverables.md) の検証用 (§6 V5・V7・V8)。
**条件は試験用** (ユーザの実条件ではない): 空気 (semiperfect)・**Tt 1250 K (ユーザ指示 2026-10-13)**・Pt 3.0 MPa (出口静圧が下がりすぎないよう選定: 出口 221.8 K・5.1 kPa)・出口マッハ 5・物理スロート径 50 mm・
入口半径 150 mm・直管 0.20 m・縮流部 0.30 m・設計スロート→出口 1.20 m・等温壁 350 K。
要求カード: `card_m5_air_demo.yaml` (`design/forge_design/intake/examples/m5_air_demo.yaml` の写し)。

- CFD 前の段 (カード → YAML → 寸法 → 物理壁 → メッシュ → 成果物) は GPU を使わない。
- 熱量の単位検査 (V5) のために、粗い格子の短い 2D の Euler → NS を手元で回す (数分)。設計の評価には使わない。

## 計算 run 一覧

| run | 目的・主要設定差分 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `run_0001_precfd_final` | **CFD 前**: 要求カード `card_m5_air_demo.yaml` → `problem_m5_air_final.yaml` (雛形 axismach_final_v1、r_t 0.024968985 m、ni 1600×97、壁第 1 層 1.3e-5 / スロート 3.1e-6) → `deltastar_loop --init-integral --prepare-only` (CONTUR 積分法の δ_r の物理壁、`--euler-ref none_precfd`) → `forge_design.export.deliverables` | 品質 PASS (AR 最大 4841.8・skew 0.419)。成果物 `deliverables/` (20 ファイル): 寸法ゲート PASS (拡大部長 +0.066 mm、他は ≤2e-6 mm)・STEP→gmsh 再構築 PASS・Fluent 検査 PASS・`thermal_use` 未評価 (CFD 前)。出口半径 140.075 mm | active (V7) |
| `run_0002_euler` | Euler の参照 (`problem_m5_air_final_euler.yaml` = 最終 YAML から physical_wall_repr を除いた写し)、soft 起動 → 本段 cfl 6・implicitRelax 0.7・12000 step、手元 (RTX 3060) | 完走 (1 分 16 秒)・NaN 0・**NOT CONVERGED** (全列が下降中、rms_ro 1.4 桁)。metrics: 軸出口 M 4.9971・dM_max 0.06 % M_d・出口コア M 4.9973 (未収束の値) | 破棄予定 (run_0003 で置換) |
| `run_0003_euler_36k` | run_0002 と同じで本段 36000 step | 完走 (約 4 分)・NaN 0・**NOT CONVERGED** (2.6 桁、rms_roUy 横ばい・他は下降; case/45 の Euler と同じ warm 床の型)。参考値 (未収束): 軸出口 M 4.9996・出口コア M 4.9993・オーバーシュート +0.002 % | active (Euler 参照) |
| (`run_0004_ns_pass0_iface`) | NS pass 0 を Euler の場から始めようとした | **prepare で停止**: interp_field が化学種の照合で拒否 (Euler に輸送の記録が無い)。計算はしていない。ディレクトリは削除 | 削除済 |
| `run_0005_ns_pass0_iface` | NS pass 0 (CONTUR 積分法の物理壁、等温壁 350 K、`output: {interfaceDiag: 1}`)、IC = prepare_ns の 1 次元等エントロピー場、段階起動 full + 本段 cfl 1・12000 step、`run_ns_ifacediag.sh` | 完走 (1 分 50 秒)・NaN 0・**NOT CONVERGED** (約 1 桁、全列下降中)。壁ダンプに `iface_Qf_eff` (valid)。成果物の CFD 表: `thermal_use` FAIL (未収束)。予備の V5: 表の q と −qwall・q_compact の差 1.7〜4.9 % | 破棄予定 (run_0006 で置換) |
| `run_0006_ns_pass0_iface_60k` | run_0005 と同じで本段 60000 step (cfl 1) | 完走 (約 8 分)・NaN 0・**NOT CONVERGED** (`check_convergence --segment` main: ρ 1.5 桁・ρE 2.0 桁が下降中、他 2.9〜3.1 桁)。総熱量は 12k の 1.51 MW → 60k で 0.73 MW と**未整定**。成果物 `deliverables/` は PASS (寸法ゲート PASS)、`thermal_use` **FAIL** (未収束; 熱量を熱設計に使わない)。**V5 合格**: 表の q (2π iface_Qf_eff / A_node) と −qwall・iface_q_compact の差 1.3〜2.6 % (r = 25〜150 mm の 9 点、許容 5 %)、Σ Q = Σ q·A (相対 0)、符号一致、iface_q_eff / q = r (旧定義なら 7〜40 倍のずれ) | active (V5 の正本) |
| `run_0007_precfd_ltol005` | **V7 の対照** (codex result レビュー m6): run_0001 と同じカード・同じ手順で、from_length の内部許容差だけを既定 0.05 r_t にした写し `problem_m5_air_final_ltol005.yaml` (CFD 前) | 寸法ゲート **FAIL**: 拡大部長 −1.907 mm (公差 1 mm)。カード由来の 0.02 r_t の run_0001 は +0.066 mm で PASS → 内部許容差をカードから決める理由の裏付け | active (V7 の対照) |
