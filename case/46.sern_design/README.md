# case/46.sern_design — ⑤ SERN 設計チェーンの検証ケース

計画: [`plans/active/tooling-nozzle-sern-chain.md`](../../plans/active/tooling-nozzle-sern-chain.md)。
問題定義 YAML (`problem_*.yaml`) → `design/forge_design/evaluate/runner_sern.py` で
逆設計 → 2 バンド構造メッシュ (スリットカウル) → forge 平面 2D → 力係数 (`metrics.json`)。

実行例 (design/ で):

```bash
PYTHONPATH=. .venv-opt/bin/python -m forge_design.evaluate.runner_sern \
    ../case/46.sern_design/problem_smoke_euler.yaml ../case/46.sern_design/run_NNNN_<slug>
```

境界タグ: inlet_nozzle=1 / inlet_ext=2 / outlet=3 / ramp=4 / cowl_in=5 / cowl_out=6 / bottom=7 / top_out=8。
壁は `outputHDFflg: 1` で `res_wall_<id>_<step>.h5` を吐き、力係数はそこから積分する
(`metrics/sern_forces.py`)。規約: 推力 = 壁力の −x、揚力 = +y、モーメント = 頭上げ正。

**受理ゲート (2026-09-13, plan §4.7/§4.13 R1)**: 各 run の `metrics.json` の `gates` に verdict (PASS/FAIL) と fail_class
(`DIVERGED` / `RESIDUAL_RISING` / `NOT_CONVERGED` / `UNSTEADY` / `NO_FORCES`) を残す。力係数履歴は `force_history.csv` に書き、
`python3 solver_density_cuda/tools/check_quasisteady.py --series-csv <run>/force_history.csv --series-cols C_T_with_shear,C_L,C_M --drift 0.02 --osc 0.05`
で正式ツールの VERDICT を再取得できる。既存キャンペーンの再判定は `driver_sern <problem> <campaign> --rejudge <out_dir>`。

## 計算 run 一覧

| run | 目的・主要設定差分 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `run_0001_smoke_euler_cell` | S4(a) スモーク初回。**TE ノードを重複させたメッシュに幅 0 の隙間 (未タグ境界辺 2 本)** があった版。段階起動 soft 3000 + 本段 6000 (cfl 4) | C_T 0.9660 / C_L 0.1426 / C_M −0.938 (run_0002 と同値; 隙間の影響なし)。残差 1.4 桁プラトー | 破棄予定 (run_0002 で置換) |
| `run_0002_smoke_euler_cell` | S4(a) スモーク: S1 テスト 4 の設計 (M_in 2.5, ランプ 15°, カウル 5°/1H, M_c 3.9, f 0.45, p_ext/p_in 0.05, M∞ 6) を cell Euler (slip) で評価し MOC と照合。修正メッシュ (TE 共有ノード) | **C_T 0.9660 (MOC 0.9666, −0.05 %) / C_L 0.1427 (0.1392) / C_M −0.939 (−0.946)**、力係数は STEADY (500 step 毎の出力で 1e-5 以内)。残差は `NOT CONVERGED` (1.4 桁プラトー、せん断層・カウル衝撃の cell 床)。`metrics.json`, `mach_field.png`, `wall_pressure_cfd_vs_moc.png`, `residual_history.png`, `MESH_QUALITY.txt` PASS (AR 120, skew 0.29) | active (ref) |
| `run_0003_nasa_Lc20_euler` | S4(b) NASA TM X-71972 傾向照合: 平板ランプ 20°/18.54H, カウル 6°, **カウル長 2.0H**, M∞ 10 (q 71850 Pa), γ 1.3, `geometry.mode: straight` | C_T 0.9605 / C_L +0.035 / C_M −0.601 (内面のみ C_T 0.9640 vs MOC 0.9628)。力係数 STEADY、残差プラトー (`NOT CONVERGED`)。`metrics.json`, `mach_field.png`, `wall_pressure_cfd_vs_moc.png` | active (ref) |
| `run_0004_nasa_Lc312_euler` | 同上、カウル長 **3.12H (NASA 基準形)** | C_T 0.9704 / C_L −0.000 / C_M −0.685 (内面のみ 0.9760 vs MOC 0.9754)。力係数 STEADY、残差プラトー (`NOT CONVERGED`)。`metrics.json`, `mach_field.png`, `wall_pressure_cfd_vs_moc.png` | active (ref) |
| `run_0005_nasa_Lc45_euler` | 同上、カウル長 **4.5H** | C_T 0.9720 / C_L −0.014 / C_M −0.641 (内面のみ 0.9800 vs MOC 0.9795)。力係数 STEADY、残差プラトー (`NOT CONVERGED`)。`metrics.json`, `mach_field.png`, `wall_pressure_cfd_vs_moc.png` | active (ref) |
| `run_0006_nasa_cowl3deg_euler` | 同上、カウル長 3.12H、**カウル角 3°** | C_T 0.9786 / C_L −0.031 / C_M −0.904 (内面のみ 0.9797 vs MOC 0.9791)。力係数 STEADY、残差プラトー (`NOT CONVERGED`)。`metrics.json`, `mach_field.png`, `wall_pressure_cfd_vs_moc.png` | active (ref) |
| `run_0007_nasa_cowl12deg_euler` | 同上、カウル長 3.12H、**カウル角 12°** | C_T 0.9314 / C_L +0.120 / C_M −0.287 (内面のみ 0.9662 vs MOC 0.9656; カウル外面の衝撃圧が推力 −0.035)。力係数 STEADY、残差プラトー (`NOT CONVERGED`)。`metrics.json`, `mach_field.png`, `wall_pressure_cfd_vs_moc.png` | active (ref) |

| `run_0008_smoke_sst_node` | RANS 化 (S3): node + SST (壁関数), 本段 cfl 4 | **発散** (本段 step 3, ω が inf) | 破棄予定 |
| `run_0009_smoke_sst_node_cfl1` | 同上、本段 cfl 1 | **発散** (step 7, カウル角部直下流 x≈0.13H の内面側で roOmega inf。node+SST 固有、未解決) | 破棄予定 |
| `run_0011_smoke_sst_cell_cfl2` | RANS 化: **cell + SST (wallTreatmentSST 1, 第一セル 4e-3 H, y+~20–40)**, 本段 cfl 2, 8000 step | **C_T(p) 0.9685 / C_T(p+τ) 0.9602 (摩擦 −0.0085) / C_L 0.154 / C_M −0.980**、力係数 STEADY。残差 rms_roK/roOmega は増加後プラトー (`NOT CONVERGED`)。ランプ壁圧に 3% の偶奇鋸歯。`deltastar_vs_euler.png`, `wall_offset_deltastar_vs_euler.json` | active (ref) |
| `run_0012_smoke_sst_cell_dstar` | S5 δ* 補正 (単独抽出版 δ*: 膨張扇の ρu 変化を欠損と誤認し δ*/H 0.24 の過大値) | C_T 0.958 — δ* が誤りで無効 | 破棄予定 |
| `run_0013_smoke_sst_cell_dstar_euler` | **S5 δ* 一発補正**: run_0011 と run_0002 (Euler) の質量流束欠損差から δ*(x) を抽出 (ランプ 0.009→0.11 H, カウル 0.006→0.013 H) し壁を法線オフセット、cell SST で再評価 | **C_T(p) 0.9685 (不変) / C_T(p+τ) 0.9602 / C_L 0.146 (無補正 0.154 → Euler 0.143) / C_M −0.930 (−0.980 → Euler −0.939)**。壁圧が MOC の破線に全域で乗る (`wall_pressure_cfd_vs_moc.png`) | active (ref) |
| `run_0010_moo_euler_2op/` | **S6 多作動点 MOO** (cell Euler; **注意: 作動点ごとに kernel を再設計していたバグ入り** — 形状が作動点間で一致しない。傾向確認用) (cell Euler, 作動点 cruise M∞6/p 0.05 p_in [w 0.6] + accel M∞4/0.15 [w 0.4], dv 5 個, LHS 12 + EHVI infill 3×2) | 18 評価 / 17 PASS / 1 INFEASIBLE、HV 1.179。パレート 8 点 (L 3.9–12.5 H, C_T,w 0.960–0.977)。`pareto.json`, `pareto.png`, `ledger.jsonl`, 各点 `doe_NNN_<op>/`, `inf_II_J_<op>/` | active (ref) |

| `run_0014_smoke_sst_node_thick` | node+SST 再挑戦 1: カウル板厚 0.2 % H (入口側で 0 に絞る台形) + interp 移植 | **soft 段 step 4 で発散** — 入口角 (x=−L_up, y=0) の 1 ノードが inlet×2 + wall×2 の 4 境界を持つ形になったため | 破棄予定 |
| `run_0015_smoke_sst_node_thick_idx` | node+SST 再挑戦 2: 板厚 0.2 % H (入口から一定、TE 手前で 0) + **stage 間の場移植を index コピーに変更** | **完走**: C_T(p) 0.9691 / C_T(p+τ) 0.9626 (摩擦 −0.0065) / C_L 0.155 / C_M −0.978、力 STEADY。rms_roOmega は本段開始直後から 3e18 一定 (壁ノード ω ピン留めの診断値、場は健全) | active (ref) |
| `run_0016_smoke_sst_node_t0_idx` | node+SST 再挑戦 3: **板厚 0 (元のスリット) + index コピー移植** — 真因の切り分け | **完走**、run_0015 と同値 (C_T(p) 0.9691 / C_L 0.155 / C_M −0.981) → **真因は interp_field の最近傍移植が座標一致の双子壁ノードを同じ元ノードに写していたこと** (板厚は不要) | active (ref) |
| `run_0017_moo_sst_node_2op/` | **S6 多作動点 MOO (RANS 版; 作動点ごとの再設計バグ入り、run_0019 で置換)**: node + SST、作動点 cruise (M∞6, p_ext/p_in 0.05, w 0.6) + accel (M∞3.5, 0.20 過膨張, w 0.4)、目的 = 摩擦込み C_T の重み付き平均と L_ramp、**制約 C_M,w ≥ −2.5** (x_ref −20H)、LHS 10 + EHVI 2×2 | 14 評価 / 10 PASS / 3 は C_M 制約で除外 / 1 は key point 不成立。HV 0.953、パレート 5 点 (L 5.5–9.7 H で C_T,w 0.960–0.967)。剥離割合はランプ全点で 0 (accel の p_ext 0.2 p_in では剥離せず)。1 点 ≈ 76 s (2 作動点)。`pareto.json`, `pareto.png`, `ledger.jsonl` | active (ref) |
| `run_0018_smoke_sst_node_lownpr` | 低 NPR 作動点の単点評価: 設計点 (p_ext 0.05) の形状を **M∞1.5, p_ext/p_in 0.6** (遷音速加速相当の過膨張) で node SST 評価 | C_T(p) 0.928 / C_T(p+τ) 0.920 / **C_L −0.316 / C_M +1.82 (頭上げ)**、剥離ゼロ (ランプ後半で 0.2→0.6 p_in へ滑らかに再圧縮)、STEADY | active (ref) |
| `run_0019_moo_sst_node_3op/` | **S6 RANS 版 MOO (3 作動点・設計点固定版)**: node SST、cruise (0.05, w 0.5) / accel (0.20, w 0.3) / lownpr (0.60, M∞1.5, w 0.2)、C_M,w ≥ −2.5、LHS 10 + EHVI 2×2。**run_0010/0017 は作動点ごとに kernel (自由境界圧) を再設計していたバグがあり、作動点間で形状が一致しない → 本 run で置換** | 14 評価 / 10 PASS / 3 は C_M 制約 / 1 は key point 不成立、HV 0.983。**低 NPR 点が支配的**: 長いランプは lownpr で C_T 0.83–0.90 に落ち、前線は最短ランプ (L 3.8H, C_T,w 0.9607) の 1 点に退化。lownpr でランプ剥離が出始める (sep_frac 最大 0.10 = doe_004)。`pareto.json`, `pareto.png` | active (ref) |
| `run_0020_smoke_sst_node_3d` | **S7 3D 確認 (外側空間あり)**: スモーク設計を有限スパン W=2H + 側壁 (x ≤ L_cowl) + 外側空間 1.5H の 3D hex (52.5 万セル, 品質 PASS) で node SST 評価。z=0 対称面 | **未解決**: 4 回投入して全て soft 段 step 3–8 で NaN。判明した原因と対処: ① 入口面の幅外側が inlet_nozzle になっていた (修正)、② ランプ∩側壁の共有ノードが 2 種の入口に属した (ランプ線も 2 重化)、③ index IC がカウル外面ノードを排気にしていた (修正)。それでも入口面直下・カウル横端 (z≈W/2) の外部流ノードから発散 → 横端トポロジ (カウル横端の凸角 + 側壁∩入口線) が残課題 | active (未解決) |
| `run_0021_diag3d_euler_sw` | 3D 切り分け (a): 側壁あり・Euler slip (最終版: index IC・ランプ線 2 重化) | **soft 段 (1 次) 完走、本段 (2 次, cfl 1) step 5 で NaN**: ノズル幅外 (z>W/2) のランプ角部 (x≈0) 直下、スパン全域。M∞6 の外部流が 15° 凸角で p/10 以下に膨張する場所 (ノズル内は M2.5 で無害) = 2 次再構成の負圧 | 診断 (ref) |
| `run_0022_diag3d_sst_nosw` | 3D 切り分け (b): 側壁なし (排気が横方向に開放)・SST、重心ベース IC | soft 段 step 3 で ω NaN (ノズル幅端 z=W/2 の 20 倍圧力不連続、非物理な構成) | 破棄予定 |
| `run_0023_smoke_sst_node_3d_noouter` | **S7 3D 基準 (外側空間なし)**: z=W/2 を側壁の境界壁にした 2D 押し出し + 側壁境界層 (33 万セル)。node SST | **完走・2D と一致**: C_T(p) 0.9699 / C_T(p+τ) 0.9619 / C_L 0.156 / C_M −0.969 (2D run_0016: 0.9691 / 0.9626 / 0.155 / −0.981)。3D node パイプライン (メッシュ・IC・段階起動・quad 力積分) は健全 | active (ref) |
| `run_0024_diag3d_euler_accel` | 3D (外側空間あり) を加速作動点 (M∞3.5, p_ext/p_in 0.2) で Euler slip: 外側ランプ角部の膨張を緩めて横端トポロジを検証 | **soft 段 (1 次) 2000 step 完走** (横端トポロジは Euler で成立)、本段 (2 次) step 425 でランプ後縁下流の top_out (静圧出口、流れが境界に平行) から NaN | 診断 (ref) |
| `run_0025_smoke_sst_node_3d_accel` | 同上を node SST (冷間) | soft 段 step 5 でカウル後縁∩側壁後縁直後の 3 重点 (下側せん断層 × 横せん断層, x 1.14H) で ω 発散 | 破棄予定 |
| `run_0026_smoke_sst_node_3d_accel_warm` | 同上、Euler (run_0024 soft 段末尾) から暖機 | 3 重点は通過 (step 237 まで) したが top_out 付近で ω 発散。**暖機コピーが Euler の wall_dist (番兵値) まで上書きしていたバグ** → 壁距離全滅が真因の可能性大 | 破棄予定 |
| `run_0027_diag3d_euler_accel_slip` | run_0024 の top_out を slip 壁 (ランプ延長) にした Euler | **完走** (soft 1500 + 本段 3000, 力 STEADY): C_T 0.932 / C_L −0.204 / C_M +1.07。ランプ幅内: 揚力 −0.019・推力 −0.005、幅外 (機体下面): 揚力 −0.104・推力 −0.025 (外部流の膨張による負圧)。**外側空間ありの 3D は Euler で成立** | active (ref) |
| `run_0028_smoke_sst_node_3d_accel_warm_slip` | node SST、top_out slip、Euler 暖機 (状態量のみコピー、soft CFL 0.25) | soft 段 step 686 で **カウル後縁∩側壁後縁の 3 重点** (x 1.0H, y −0.1, z=W/2, 壁距離 0.4–3 mm) で ω → inf。SST の 2 せん断層交差の頑健性問題 (Euler は同格子で成立) | 破棄予定 |
| `run_0029_ref2d_euler_node_accel_slip` | run_0027 の 2D 参照: node Euler slip、加速点、top_out slip | C_T 0.976 / C_L +0.004 / C_M −0.001 (ランプ T +0.014 L +0.085, cowl_in T +0.005 L −0.091) | active (ref) |
| `run_0030_cycle3op_prep_{m6_on,m10_on,m4_off}` | **作動点サイクル値化の検算** (plan §4.10, `problem_moo_sst_node_cycle3op.yaml`)。CFD 無しの prepare-only ×3。設計点固定 (作動点が inflow/gas を上書きしても形状不変) と作動点ごとの CFD 条件を確認 | 3 点で**輪郭 sha 一致** (L_ramp 11.357 H, key point (5.767, 0.923), θ_e −2.659°)。CFD 条件は NPR 35.4 / 56.5 / 2.8、γ 1.183 / 1.226 / 1.399、音速 968.1 / 1032.2 / 368.4 m/s (CEA と一致)。`MESH_QUALITY.txt` **PASS** (AR max 443, skew 0.464) | active (ref) |
| `run_0031_cycle_euler_m6on` | 同上の設計点 m6_on を node Euler で回して MOC と突き合わせ (`problem_cycle_euler_node_m6on.yaml`, cfl 4, soft 3000 + 本段 6000) | **C_T 0.9059 (MOC 0.9073, −0.15 %) / C_L 0.3135 (0.3105) / C_M −8.227 (−8.168, −0.7 %)**。力係数は 12 スナップショットで **STEADY** (C_T trend 1.2e-6, fluct 8.8e-6)。`check_convergence.py` は **NOT CONVERGED** (1.0–1.7 桁のプラトー: run_0002 と同性格のせん断層・カウル衝撃の床)。`metrics.json` (`on_design_point: true`), `MESH_QUALITY.txt` PASS | active (ref) |
| `run_0032_cycle_euler_m4off` | 最も過膨張が強い作動点 m4_off (NPR 2.8, φ=0 空気 γ 1.399) の成立確認。形状・メッシュは run_0031 と同一 (設計点固定)。**ランプ側に外気が無い旧 2 バンド** (§4.11 前) | C_T 0.9736 / C_L −0.044 / C_M **+0.959**、力係数 STEADY、残差 2.3–2.8 桁プラトー (NOT CONVERGED)。Euler なので `sep_frac` は None。**後縁に外圧が無いため剥離評価には使えない** (run_0034 と比較用) | active (ref) |
| `run_0033_cycle_euler_m10on` | 3 点目 m10_on (NPR 56.5, φ=1.5, γ 1.226)。旧 2 バンド | C_T 0.9226 / C_L 0.224 / C_M −6.105、STEADY、残差 1.6–1.7 桁プラトー (NOT CONVERGED) | active (ref) |
| `run_0034_cycle_euler_m4off_exttop` | **ランプ側外部流ブロック (plan §4.11, `mesh.ext_top: 1`) の初回**: m4_off を node Euler で。鉛直 base + wake ブロック版 | **soft 段 step 5 で発散** (NaN はランプ後縁 x 11.07–11.9, y 2.68–2.85 = base 角に集中)。node の 90° 二重 slip 角 (ランプ∩base、上面∩base) が原因。`MESH_QUALITY.txt` PASS (AR 443, skew 0.464, 60.2k cells)、`mesh_blocks.png` | 破棄予定 (run_0036 で置換) |
| `run_0035_cycle_sst_m4off_exttop` | 同メッシュで node SST | **soft 段 step 702 で rms_roOmega → inf**。NaN は**カウル TE** (x 1.20, y −0.11) で ext_top とは無関係 (カウル板厚 0 の node 双子ノード問題の疑い → run_0037 で切り分け) | 破棄予定 |
| `run_0036_cycle_euler_m4off_taper` | ext_top を**機体後端テーパ** (`vehicle_taper: 2.0`: 後縁手前 2H で厚さ 0 に絞り TE 共有、base/wake 無し = カウルと同じ構造) に変更した node Euler m4_off | **成立**。C_T 0.9736 / C_L −0.044 / C_M +0.973 (run_0032: 0.9736 / −0.044 / +0.959)、力係数 STEADY、`check_quasisteady` ALL STEADY、残差 0.8 桁 still converging (NOT CONVERGED)。**後縁圧 p/p∞ 1.317 (run_0032, 領域の産物) → 0.932 (機体上面 boat-tail 膨張後の外気)**: 外気は届いた。ただし影響は後縁 0.1H のみで、ランプ圧 (0.5–7H で 0.84 p∞、以降カウル側プルーム境界からの圧縮で 1.3 p∞ へ) はカウル側外気が支配。`ramp_te_pressure_vs_run_0032.png`, `mach_field.png` | active (ref) |
| `run_0037_cycle_sst_m4off_2band_thick` | node SST m4_off を**旧 2 バンド** + カウル板厚 2e-3 (`problem_cycle_sst_node_2band_thick.yaml`) で。run_0035 のカウル TE ω 発散の切り分け | **完走 (NaN なし)** → 板厚 2e-3 で解決。C_T 0.9759 (摩擦込み 0.9701) / C_L −0.036 / C_M +0.796、STEADY、残差 1.1–1.4 桁プラトー。`sep_frac_ramp` = 0.0 | active (ref) |
| `run_0038_cycle_sst_m4off_taper` | node SST m4_off を ext_top テーパ版 + カウル板厚 2e-3 で (`problem_moo_sst_node_cycle3op.yaml`)。後縁の外気込みで剥離が出るか | 完走 (NaN なし)。C_T 0.9758 (摩擦込み **0.9701** = run_0037 と同値) / C_L −0.035 / C_M +0.75、STEADY、`check_quasisteady` ALL STEADY、残差 1.1–1.4 桁プラトー。**`sep_frac_ramp` = 0.0、τ_w>0 が 325/325 点**: 外気を入れても L_cowl 1.2 の形状は m4_off で剥離しない (ランプは 0.84 p∞ の軽い過膨張のみ) → 領域の問題ではなく形状の性質 | active (ref) |
| `run_0039_cycle_sst_m4off_Lc25` | 同条件で **L_cowl を dv 上限 2.5** に (`problem_cycle_sst_node_Lc25.yaml`)。m4_off が MOO の剥離制約として効く作動点かの感度確認 | 完走 (NaN なし)。L_ramp 13.37、C_T 摩擦込み 0.9625 / C_M +2.27、STEADY、残差プラトー。**`sep_frac_ramp` = 0.0、τ_w>0 が 327/327**。**cowl TE 圧は 1.35 (Lc 1.2) / 1.59 (Lc 2.5) p∞ で不足膨張** → NPR 2.8 ではカウル側衝撃が存在せず、dv 箱のどの形状も剥離しない。ランプ最小圧 0.79 / 0.57 p∞。`wall_pressure_m4off_Lc12_vs_Lc25.png` | active (ref) |
| `run_0040_moo_cycle3op` | S6 MOO 再取得の**第 1 回 (失敗)**: 板厚 2e-3・soft cfl 0.5・本段 cfl 1.0 の構成 | `doe_000` の m6_on は PASS (C_T 0.930 / C_M −6.81) だが **m10_on が soft 段 step 3 で `roOmega` 発散** → 全評価が同じ経路で落ちるので即停止。起動レシピの確立 (run_0041–0053) に移行 | 破棄予定 |
| `run_0041_m10on_softcfl05` / `run_0042_m10on_softcfl02` | m10_on の発散切り分け: soft 段 CFL 0.5 / 0.2 | 0.5 = step 3 発散、NaN は **x −0.45, y −0.0095 = 入口面 ∩ カウル外面の角**。0.2 = step 1719 (本段) 発散、NaN は **x 1.196, y −0.108 = カウル後縁**。**CFL では治らず場所が移るだけ** | 破棄予定 |
| `run_0043_m10on_thick005` / `run_0044_m10on_thick010` | カウル板厚 5e-3 / 1e-2 (soft cfl 0.2) | 両方完走、C_T 0.9246 / 摩擦込み 0.9175 / C_M −6.22 と −6.20 で**板厚は力係数を 0.3 % 以内でしか動かさない** (数値措置であって物理ではない) | active (ref) |
| `run_0045_prod_check_*` | 板厚 5e-3 + soft cfl 0.2 で 3 作動点確認 | m6_on・m10_on は PASS だが **m4_off が step 283 でカウル後縁発散** → 板厚は作動点で要求が逆転し ( m4_off は 2e-3、m10_on は 5e-3)、MOO では使えないと確定 | 破棄予定 |
| `run_0046_warmlam_m4off_th005` / `run_0047_warmlam_m10on_th005` | **層流暖機段** (`turbulence: none` + 粘性、1 次) を soft の前に入れる案。暖機 cfl 0.5 | m4_off **PASS** (C_T 0.9759 / 摩擦込み 0.9701 = run_0038 と一致)。m10_on は暖機自体が入口角で step 3 発散 → 暖機も cfl 0.2 が要る | 破棄予定 |
| `run_0048_recipe_*` | 3 段レシピ (層流暖機 cfl 0.2 → SST soft cfl 0.2 → SST 本段 cfl 1.0) を 3 作動点で | **3 点とも完走・全て STEADY**。m6_on 0.9089 (0.9018) / m10_on 0.9246 (0.9175) / m4_off 0.9759 (0.9701)、sep 0 | active (ref) |
| `run_0049_wired_m10on` | 同レシピを `run_staged` に組み込んだ版で m10_on | **本段 step 171 で発散**。run_0048 とメッシュ・config がバイト一致なのに結果が割れる = **本段 cfl 1.0 は限界的** | 破棄予定 |
| `run_0050_cfl05_m10on_a` / `run_0051_cfl05_m10on_b` / `run_0052_cfl05_m6on` / `run_0053_cfl05_m4off` | 本段 cfl 0.5 / 6000 step に下げて再現性確認 (m10_on は 2 回) | **4 本とも完走・全て STEADY**。m10_on は 2 回とも C_T 0.9246 / 摩擦込み 0.9175 / C_M −6.216 で**完全一致**。m6_on 0.9089 (0.9018) / m4_off 0.9759 (0.9701)、cfl 1.0 の値と一致 → **確定レシピ** (plan §4.12) | active (ref) |
| `run_0054_moo_cycle3op` | **S6 MOO 再取得の第 2 回**: 確定レシピ (plan §4.12) で `problem_moo_sst_node_cycle3op.yaml`。3 作動点 × node SST × ext_top テーパ × 板厚 5e-3 × 3 段起動、dv 5 変数、`cm_min` −7.0、剥離制約なし。DOE 40 + infill 8×2 = 56 評価、HV ref (−0.75, 20.5)。1 評価 ≈ 8 分 (1 run ≈ 156 s × 3) → 約 7.5 時間 | 実行中。`pareto.json` / `ledger.jsonl` | active |
| `run_0060`–`run_0067` | 起動レシピの切り分け: ランプ角部の丸め (`mesh.ramp_fillet`) 導入、dv 箱の隅と既定 dv での確認、緩レシピ試験 | 丸めで角部発散は解消 (θ_r0 22° でも通る)。既定 dv の力係数は丸め前と一致 (0.9091/0.9246/0.9758)。残る m4_off×短カウルは緩レシピでも救えず → 梯子方式へ (plan §4.13) | active (ref) / 一部破棄予定 |
| `run_0069_moo_cycle3op` | MOO 第 3 回 (丸め + 梯子, L_cowl 上限 2.5)。**22 評価 / PASS 14 で打ち切り** | 失敗は L_cowl ≲ 0.7 に集中。最良 doe_004 C_T_w 0.9434 @ Lc 2.30。**C_T が L_cowl 上限に張り付いていたため箱を広げて再投入** (plan §8-8) | 破棄予定 |
| `run_0070_aws_smoke_m6on` | **AWS g5 (A10G, CUDA 13.2)** での動作確認。設計点 m6_on | C_T 0.9091 / 摩擦込み 0.9019 / C_M −8.259、STEADY。ローカル run_0063 と一致 (C_L のみ 2e-4 差) | active (ref) |
| `run_0071_moo_cycle3op_wide` | **S6 MOO 本番 (AWS)**: L_cowl 上限 4.5 H の広い箱。3 作動点 × node SST × 丸め × 3 段起動 + 梯子、DOE 40 + infill 8×2 = 56 評価 | 実行中 (AWS `~/forge/case/46.sern_design/`) | active |
| `run_0072_warm_*` / `run_0073_warm_*` | **作動点間 warm start の検証** (codex レビュー A′): m6_on の収束場を熱力学整合リマップ (ρ,u,P を入口比でスケール、目標 γ で roe 再構成、k~u², ω~u) → 適応段 500 step → mid → 本段。2 形状 (既定 dv / L_cowl 3.5) で m10_on を cold と warm 両方 | **cold と一致**: A C_T 0.92462 = 0.92462 (差 0.000 %)、摩擦込み 0.91749 = 0.91749、C_M 0.042 % 差。B C_T 0.95902 = 0.95902、C_M 0.020 % 差。**短縮 27–28 %** (148→109 s)。飛ばしたのは暖機+soft 4000 step、適応段 500 step 追加で正味 3500/12000 = 29 % — 計算どおり | active (ref) |
| `run_0074_tapercheck_m10on` | テーパ smoothstep 化後の確認 (落ちた dv inf_01_1) | **FAIL** (mid 段 step 87)。NaN は (3.07, 2.25) = 機体上面テーパ区間 | 破棄予定 |
| `run_0075_diverge_watch` | **発散の過程を捉える診断 run**: mid 段を **5 step 刻み**で出力 (元は 2000 刻みで NaN ダンプしか無く過程が見えなかった) | `CONVERGENCE_VERDICT.txt` = **DIVERGED (NaN/Inf)** (意図どおり)。順序を確定: **P 床 (1 Pa) 着地 → 負密度 (step 65, 2 ノード) → 圧力 1.5e7 Pa 暴走 → ω 発散**。膨張の中央値 0.235 p∞ は PM 予測 0.080 より緩く**物理的**、床に落ちるのは 18 % のノードのみ。h5 (step 0/65/70/75/80/85/90) と図 `diverge_taper_vacuum.png` / `vehicle_surface_pressure.png` | **active (ref・診断の正本)** |
| `run_0076_relax_*` | `implicitRelax` / `pMin` の切り分け (5 通り) | relax 0.7 単独・relax 0.7 + cfl 2.0 とも **mid 段で FAIL** → **implicitRelax は無効**。pMin 50 Pa の 2 例は**私のログ出力バグ**で未測定 | 一部 active (ref) / pMin は再試行 |
### S4(b) NASA TM X-71972 傾向照合のまとめ (2026-09-04, run_0003–0007, 図 `nasa_trends.png`, 表 `nasa_trend_table.py`)

- **内面 (ランプ + カウル内面) の力は forge Euler と MOC が全 5 形状で C_T +0.0006〜+0.0012、C_M 0.01〜0.05 以内で一致**。差の残りは
  カウル外面 (M∞10 の外部流がカウル角部で圧縮される衝撃圧) で、MOC は扱わない。カウル角 12° では外面だけで C_T −0.035。
- **カウル長** (θ_c 6°): 2.0 → 3.12 → 4.5H で C_T 0.9605 → 0.9704 → 0.9720。短縮で推力が大きく落ち、延長の追加利得は小
  — NASA の結論と同じ。ピッチモーメントは基準点依存: 入口下端基準では −0.60 / −0.69 / −0.64 だが、機体 CG 相当の前方基準
  (x_ref = −20H) では **−1.31 / −0.68 / −0.36** となり、NASA の「短いカウルで大きな頭下げ、延長で頭上げ増分」を再現する。
- **カウル角** (L_cowl 3.12H): 3° → 6° → 12° で C_T 0.9786 → 0.9704 → 0.9314 (単調減少、NASA の「角を開くと推力劣化」と同じ)。
  C_M は −0.90 → −0.69 → −0.29 (入口基準) で角を開くほど頭下げが緩む。NASA が 6° を最良としたのは飛行域全体のトリム舵角で、
  1 作動点の本比較では順位付けできない (S6 の多作動点束ねで評価する量)。
- 残差は全 run で 1.3〜1.5 桁プラトー (せん断層・衝撃の cell 床)。力係数は 500 step 出力で 1e-5 以内 (STEADY)。

### S5 / S6 のまとめ (2026-09-04)

- **RANS 化**: node + SST は本段でカウル角部直下流 (内面側 x≈0.13H) から ω が発散 (cfl 4/1 とも)。cell + SST (壁関数) は
  cfl 2 で安定し、圧力推力は Euler 比 +0.25 %、摩擦で −0.85 % (twall は forge の「流体に働く traction」規約、
  `viscousFlux_d.cu` L527 で符号確認)。node 側は未解決の課題として plan に残す。
- **δ* 一発補正 (S5)**: 単独の 0.99 max 縁判定は SERN の非一様コア (膨張扇) を欠損と誤認して 20 倍過大になる。
  **Euler 場を基準にした質量流束欠損** (`deltastar_sern_vs_euler`) で妥当な δ* (ランプ末端 0.11 H) が取れ、
  オフセット壁の RANS は壁圧・揚力・モーメントを非粘性設計値に戻した。推力は δ* に鈍感 (0.9685 不変)。
- **MOO (S6)**: 1 点 ≈ 135 s (2 作動点、GPU 共有時)。パレートは「長いランプほど C_T,w 高い」の単調前線
  (3.9H で 0.960 → 12.5H で 0.977) で、短い側は θ_r0 14–16°・M_c 3.2 付近、長い側は M_c 3.9・f 0.47。
  C_M (−20H 基準) は前線上で −0.19〜−2.6 と大きく変わり、カウル角 2° 台が頭下げ最小 — C_M を制約に入れると前線の選択が変わる。

### node + SST 発散の真因と解決 (2026-09-04, run_0014–0016)

- run_0009 の本段初期場 (soft 段末尾を `interp_field.py` で移植したもの) を調べると、カウル内面と外面の壁ノード (座標が一致する
  スリットの双子) が全 station で**同一の圧力**を持っていた (排気側 15 kPa のはずが外部流側の 2 kPa など)。合成場で
  `interp_field.py` を試すと双子 134 station 全部が同じ元ノードに写った = **最近傍補間が座標一致ノードを区別できない**。
  壁ノードと隣接内部ノードの 6 倍の圧力段差が 2 次で step 7 に爆発した。solver 側の欠陥ではない。
- 対策: stage 間の同一メッシュ移植を `runner_sern.restart_by_index` (VALUE の index コピー) に変更。板厚 0 のスリットのまま
  node+SST が完走 (run_0016)。板厚オプション (`mesh.cowl_thickness`) は残すが必須ではない (run_0015 と同値)。
- node の twall は「壁ノードに働く力」(cell は「流体に働く力」) で符号が逆。`sern_forces` は離散化で規約を切替 (摩擦 node −0.0065、
  cell −0.0085)。
- 残る node 固有の観察: (i) rms_roOmega が本段で 3.3e18 一定 (壁ノードの ω ピン留め残差が混入する診断値。場は STEADY で ω max 1.5e7 は
  外部流側の壁ノード)、(ii) 入口角 (inlet+wall) と後縁・ランプ後縁の**単ノードの圧力が外れる** (cowl_in 入口角 2.2 p_in、cowl_out 入口角
  0.36 p_in)。力積分への影響は小 (node/cell の C_L 差 0.001) だが、壁圧分布を読むときは端点を除く。

### S6 RANS 版 MOO のまとめ (2026-09-04, run_0017)

- node + SST を評価器にした 2 作動点 MOO が 1 点 76 s (GPU 専有時) で回った。摩擦込み C_T は Euler 版より 0.5〜1 % 低く、
  cruise (高 NPR) より accel (p_ext 0.2 p_in) の方が C_T が高い (外部圧が高いぶん (p−p_a) の積分が有利)。
- C_M ≥ −2.5 の制約で 3 点 (カウル角 6〜7°・長いランプ) が除外され、パレートは L 5.5〜9.7 H の短い側に寄った。
  前線上はカウル角 2〜8°、M_c 3.5〜3.8、f 0.35〜0.44。C_M 最良 (−1.47) は L 9.66H・カウル角 2.7° (doe_004)。
- **剥離は全点でゼロ**: accel の p_ext/p_in 0.2 ではランプ末端 (0.13 p_in) がやや過膨張でも超音速の順圧力勾配で付着したまま。
  剥離 (RSS/FSS) を評価するには遷音速加速に相当する p_ext/p_in ≳ 0.5 の作動点が要る (次の課題)。
- Euler 版 (run_0010) との差: 同程度の L で C_T,w が 0.5〜1 % 低い (摩擦)。前線の形 (長いほど高推力) は同じ。
| `run_0078_vehsurf_pmin1` | 機体上面を**後縁 θ_e 接線 Hermite + くさび 3°** に作り替えた形状で、`pMin` を既定 1.0 Pa に戻して落ちた dv (θ_r0 22°) を再現 | **DIVERGED**: soft 段 step 61 で NaN。ただし NaN 位置は **x 0.42, y 0.223 = ノズル内部**で、機体上面 (y≈2.29) ではない。上面の凸角は消えたが `pMin` 1 Pa は別の場所でも割れる = 発散の根因が圧力床であることの追加証拠 | 診断 (ref) |
| `run_0079_vehsurf_pmin20` | 同形状 + `p_min: 20.0` (生産設定)。新上面の妥当性確認 | **完走**。C_T **0.9184** (STEADY) / C_M −0.240 (STEADY) / C_L 0.0022 (OSCILLATING, 振幅 7.7e-5) / 剥離 0。旧上面の run_0077 (C_T 0.9187, C_M −0.266) と C_T が **0.03 % 一致** = 上面の作り替えは物理を動かさない。残差は plateau (NOT CONVERGED) だが §4.13 の受入方針どおり力係数の定常性で判定 | active (ref) |
| `run_0081_moo_cycle3op_vehsurf` | **S6 MOO 本番 (再取得)**: 新機体上面 + `p_min 20` + L_cowl 上限 4.5 H。3 作動点 × node SST、DOE 40 + infill 8×2 = 56 評価 | 実行中 (ローカル、log `/tmp/moo_0081.log`) | active |
| `run_0080_3d_sst_ladder` | **3D SST 再挑戦**: 2D の段階起動 (層流暖機 2000 + soft cfl 0.2 + mid 2000 + 本段 cfl 0.5) と `p_min 20` を 3D に移植 | **暖機段は完走** (run_0028 が落ちた三重点を通過) したが、SST 投入直後の soft 段 **step 3** で同じ三重点 (x 1.081 H, y −0.140 H, z 1.004 H) に ω 発散。soft 段 step 0 の **rms_roK = 1.1e7**。段階起動だけでは足りない | 診断 (ref) |
| `run_0082_3d_sst_lsw08` | **3D node SST 初成功**: 上に加えて `mesh3d.L_sw: 0.8` でカウル後縁と側壁後縁を x 方向に分離 | **全段完走** (524628 セル, メッシュ PASS AR 535.8 / skew 0.365, 本段 4000 step 262 s)。**C_T 0.9385 (圧力) / 0.9262 (摩擦込み) / C_L −0.1787 / C_M +0.9892**、力係数 3 つとも **STEADY** (8 スナップショット)。同条件 3D Euler run_0027 (0.932 / −0.204 / +1.07) と整合。揚力の主役はランプ幅外 (C_L −0.0903 vs 幅内 −0.0014) | **active (3D の基準 run)** |
| `run_0083_3d_sst_cycle_m6on` | 3D SST を**生産作動点 m6_on** (dv = 設計点, L_cowl 1.2, `L_sw` 0.8) で評価 | **DIVERGED**: mid 段 (2 次) step 48。位置 x/H 1.187–1.200, z/H 0.972–0.986 = **カウル後縁の刃先** (スパン端の 0.014–0.028 H 内側)。`L_sw` を離しても残る別の特異点 | 診断 (ref) |
| `run_0084_3d_sst_cycle_pmin60` | 同上 + `p_min` 20 → 60 Pa | **同じ step 48・同じ 22 ノード・同じ位置** → 圧力床は原因でない (床に着地するのは結果) | 診断 (ref) |
| `run_0085_3d_sst_cycle_barth` | 同上 + リミッタ Barth | step 34 に**悪化** → リミッタでもない | 診断 (ref) |
| `run_0086_3d_sst_cycle_thick` | 3D メッシャにカウル板厚 0.005 を実装 (z 一様) | **暖機段 step 4** で ro NaN。NaN は側壁後縁 (x/H 0.785–0.824, z/H 1.000–1.010) = 厚さ 0 の側壁スリットの内外で中間線がずれた | 破棄予定 |
| `run_0087_3d_sst_cycle_thick_ztaper` | 板厚を側壁の 2 セル手前から 0 に絞る | 暖機段は完走 (残差 1.4–2.3 桁低下)、mid 段 **step 57** (48 → 57 と改善するが未解決)。位置は run_0083 と同じ刃先 | 診断 (ref) |
| `run_0088_3d_sst_cycle_cfl01` | 同上 + mid 段 CFL 0.2 → 0.1 | mid 段 **step 109** = 57 のほぼ倍 = **同じ物理時刻** → CFL はつまみでない。次はカウル後縁の鈍頭化 | 診断 (ref) |
| `run_0089_inletprof_supersonic_tt` | **超音速入口 (`inlet_uniformVelocity`) の入口分布検証**: run_0029 (node CPG Euler 2D) の config/IC + `inletProfile: 1`, Tt(y)=2025+200·exp(−((y−0.05)/0.02)²) K, M 2.5, Ps 20 kPa を `tools/gen_inlet_profile.py gen --Tt --M --Ps` で ρ/Ux/Ps に換算 (往復誤差 1e-12 K), 600 step。手順 [procedures/inlet-profile.md](../../procedures/inlet-profile.md) | `verify`: 入口ノード ρ・U が目標と 1e-5、h0 由来 T0 が換算元 Tt と 6e-4 K で一致。NaN 0 (600 step の過渡, 収束評価対象外) | active (機能検証 ref; 入力 config/CSV を git 追跡) |
| `run_0100_species_regress_euler_lumped` / `run_0101_species_regress_euler_full` (+ `run_0102`/`0103` の 24000 step 版) | **species 統一スキーマの回帰** (plan thermophysics-cea-mole-fraction-species §6 SERN): m6_on, node Euler, 6000 step, cfl 1.0。`0100` = `evaluate.tp_species` 省略 (旧 `[EXH, AIR]` 別名 = 流れごとの lump), `0101` = `{mode: full}` (排気 11 種 ∪ 外気 = 11 種 TP + 排気トレーサ `roXi`)。問題 YAML `problem_species_regress_euler_{lumped,full}.yaml` | 両 run NaN 0、力の時系列 C_T/C_L/C_M は 0.1 % で ALL STEADY (既存 gate の 2 %/5 % より厳しい)、**C_T 0.922735 / C_L 0.274206 / C_M −7.1161 が full と lumped で相対 1e-8 / 5e-7 / 4e-7 一致** (ゲート 1e-3)、場 P/T 2.5e-5・ro 1.5e-4、排気率 ξ (full: `Xi`, lumped: `Y_EXH`) は平均 4.7e-7・最大 1.7e-4 (カウル後縁 200 ノードに局在)、ΣY=1±2.4e-7、`rms_roY{s}` 列あり。`check_convergence` は両方 NOT CONVERGED (plateau 2.3–2.7 桁)。**24000 step 版 `run_0102`/`0103`** も plateau (2.4 桁, rms_ro 8e-7; node Euler SERN の限界サイクルで PASS は取れない): 力は 0.1 % で ALL STEADY, C_T/C_L/C_M 相対 3e-8 / 1.2e-6 / 1.0e-6, 場 P/T 3.5e-4 / 5.7e-4 (限界サイクル上のスナップショット位相差), ξ 平均 4.7e-7・最大 1.7e-4。step 時間 2.33 → 4.77 ms (11 種)。GATES PASS。**`run_0104` (lumped + `tracer: exhaust`, IC roXi=roY0)**: 同一 run 内の `\|Xi − Y_EXH\|` が最大 1.8e-4 / 平均 4.6e-7 (>1e-6 は 4380 ノード) → ξ の差は配置差でなく**カーネル差** (汎用スカラ輸送 vs 化学種輸送の再構成・質量整合補正); lumped+tracer の Xi と full の Xi は最大 9.6e-5。**`run_0105` (m10_on lumped, `warm_from_run` で `run_0102` [m6_on 24k] から作動点変更)**: 種順序・datum・トレーサ有無を照合し lump の MW 変化は許容、ΣρY/ρ−1 1.2e-7、roe は目標作動点の輸送種ガスで T'=P'/(ρ'R_mix) から再構成 (反転 T との差 1.8e-3 K)、適応段 (cfl 0.5, 1 次, 2000 step) で rms_ro 3.1e-2 → 9.0e-7・rms_roY0 8.6e-4 → 8.8e-7、NaN 0、ΣY 1.2e-7 | active |
| `run_0090_rejudge_r1_gates/` | **R1 評価ゲート (plan §4.7/§4.13, 2026-09-13) で既存キャンペーンを CFD 無しで再判定** (`driver_sern --rejudge`; 元 run は読むだけ)。`run_0069_moo_cycle3op/` (rejudge_summary.md, ledger_rejudged.jsonl, pareto_rejudged.json, 各 run の metrics.json + force_history.csv)、`run_0054_moo_cycle3op/`、`run_0027_3d_bookkeeping/` + `run_0029_2d_ref/` (R2 集計分離) | run_0069: 旧 PASS 14 → **新 10** (HV 2.6835 → 2.6450)。除外 = RESIDUAL_RISING 3 (doe_004/015 m6_on: 本段末尾で ρ,ρu,ρe 残差 2 倍リバウンド; doe_014 m4_off: step 5461 から rms_roOmega 1e1→5e18 の本物の発散、力係数は STEADY のまま) + UNSTEADY 1 (doe_016 m4_off の C_M TRANSIENT-UNSETTLED, 1.3 %/tail)。run_0054: 1/1 DIVERGED 不変。run_0027 3D 再集計: ノズル C_T 0.9578 / 機体 −0.0253 / 総計 0.9325 (codex 検算と一致), 2D run_0029 0.9762 比 −1.9 % (両 run NOT CONVERGED) | active (判定記録) |
| `run_0091_smoke_r1_driver/` | **R1 ゲート入り driver の実機スモーク**: run_0069 doe_001 (旧 PASS, C_T_w 0.9333, L 11.89) を現行 `driver_sern.evaluate` で再評価 (`replay.py`; 3 作動点 m6_on / m10_on [warm from m6_on] / m4_off, 生産 YAML `problem_moo_sst_node_cycle3op.yaml`)。バイナリ = 主 checkout `solver_density_cuda/build/forge` (2026-09-10 21:16 ビルド ≈ 5db69d93 世代、worktree から symlink) | **PASS** (degraded なし, 172 s / 3 作動点): C_T_w **0.93308** (旧 run_0069 doe_001 0.93331, −0.02 % = バイナリ差・run 間ノイズ), C_M_w −3.43, L 11.89。作動点別 C_T(摩擦込)/C_M: m6_on 0.9330/−5.93, m10_on 0.9304/−4.01, m4_off 0.9373/+3.69。全 op で rc 0・場有限・残差 NaN/rising 無し (verdict は NOT CONVERGED stalled/plateau = 本ケースの性格)・C_T_with_shear,C_T,C_L,C_M 全て STEADY (`check_quasisteady --series-csv` ALL STEADY)。`pareto.json` に tag/degraded/ops ゲート要約あり。`residual_history.png` 各 run | active |
| `run_0092_r2_3d_euler_accel_vehicle/` | **R2 (集計分離) の 3D 再計算**: run_0027 と同じ加速点 Euler (node, top_out slip) を、メッシャの `vehicle` タグ (W/2 < z ≤ W_vehicle/2 = 2.0 H, その外は top_out) 入りで本段 6000 step (`problem_r2_euler_node_3d_accel_vehicle.yaml`) | **GATES PASS**, 残差 2.4–3.2 桁でプラトー (NOT CONVERGED stalled), C_T/C_L/C_M STEADY。**ノズル C_T 0.95782 (run_0027 の再集計 0.95782 と一致 = タグ分離は帳簿だけを変える)**, 機体 C_T −0.0190 (W_vehicle 4H; 遠方まで取ると −0.0253), C_L −0.0998 / 機体 −0.0836, C_M +0.724。運動量収支 (`momentum_balance.json`, `metrics/sern_momentum.py`): 閉じ残差 x 1.6 % / y 0.4 % / z 1.7 % of F_ideal, 質量 0.000 % | active |
| `run_0093_r2_2d_euler_accel_ref/` | 同上の 2D 参照 (run_0029 相当, 8000 step) | GATES PASS, 残差 2.7–3.1 桁プラトー。C_T 0.97621 / C_L 0.0043 / C_M −0.0010。**3D ノズル C_T / 2D = −1.88 %** (旧「−4.5 %」は幅外機体面の混入)。運動量収支 閉じ残差 x 0.09 % / y 0.04 % | active |
| `run_0094_frozen_tp_doe001/` | **R3 (frozen_tp) の実機確認**: run_0069 doe_001 の dv を `problem_moo_frozen_tp_cycle3op.yaml` (排気 = CEA 凍結組成 EXH, 外気 = AIR, thermalMethod 2, thermoHrefTemp 298.15) で driver 評価 (`replay.py`, 3 作動点, m10_on は m6_on から warm)。prepare 段の整合確認 (ΣY=1, IC の T 反転 2328.0/220.8 K 厳密, q∞ 71.9 kPa) は済 | **m6_on 標準レシピは暖機段 step 568 で ro NaN** (ランプ後縁ノード x/H 10.4–11.0, y/H 3.75–3.9 = 排気 EXH と機体側 AIR の接触で T が Newton 床 50 K に落ちる; CPG 同レシピは完走)、**緩レシピ (cfl 半分・step 倍) は完走 GATES PASS**: C_T(p) 0.9520 / 摩擦込み 0.9447 / C_L / C_M −5.11, L_ramp 10.64 (凍結 γ 1.247・M_in 1.631 で形状が CPG 版 11.89 から変わる)。m10_on は warm_from_run の frozen 経路のバグ (作動点未適用の YAML から組成を読む) で ERROR → 修正、run_0095 で再実行 | 診断 (ref) |
| `run_0095_frozen_tp_doe001_warm01/` | run_0094 の再実行: `opt.warm_lam_cfl` 0.2 → **0.1** (steps 4000) を frozen_tp の標準レシピにし、warm_from_run の frozen 経路 (P を相似スケール、Y は場のまま、T' = P'/(ρ'R_mix)、目標擬似種で roe 再構成) を初使用 | **m6_on: 標準レシピで GATES PASS** (retry 不要; C_T(p) 0.9520 / 摩擦込み 0.9447 / C_M −5.11 = run_0094 retry と 6 桁一致)。**m10_on: warm/cold とも FAIL/RESIDUAL_RISING** (rc 0, 力係数は STEADY で 0.93528 が両者 5 桁一致; 残差は step ~2300 の最小から末尾で 3–4 倍に増加)。原因箇所 = 上流区間 x/H −0.48 の**カウル板下面 (外気側, M∞10)**: 壁 T 4570 K (M10 断熱回復温度) の 4–5 ノード外側で T が Newton 床 50 K・P 160 Pa (p∞ 1026) に落ちる冷点。CPG (run_0091) も同じ場所に T 116 K / P 397 Pa の冷点があり (残差はプラトーで PASS)、TP でより極端になった = R3 固有でなく M10 外気×板下面の近壁病理 (排気側の板上面壁 T も 4600 K で回復温度 3000 K を超える)。m4_off は未実行 (m10_on で評価が止まる) | active (診断) |
| `run_0096_frozen_m10on_long/` | frozen m10_on の残差上昇が発散か過渡かを見る 12000 step run (cold 標準レシピ)。**注意: `runner_sern` 直叩きなので dv は YAML 既定値 (L_ramp 10.09) = doe_001 (10.64) とは別形状** | 残差は step ~4000 の最小 2.3e-5 から 6000 で 5.4e-5 に上がった後**飽和** (12000 まで 4.7–5.6e-5 のプラトー) → GATES PASS (stalled/plateau)。力係数 STEADY (摩擦込み C_T 0.92586)。= 6000 step 時点の RESIDUAL_RISING は過渡リバウンドで、本段を伸ばせば通る | 診断 (ref) |
| `run_0097_frozen_m10on_doe001_long/` | 同上を **doe_001 形状** (`run_0095/doe_001_frozen.yaml`, cold 標準レシピ, 本段 12000) で | **GATES PASS**: 残差は 6000 以降 5–7e-5 で飽和 (stalled/plateau)、摩擦込み C_T **0.93528** = run_0095 の warm/cold-relaxed と 5 桁一致 (状態は一意)。→ frozen_tp 生産 YAML の本段を 12000 step に | active |
| `run_0098_frozen_tp_doe001_full/` | frozen_tp 3 作動点評価 (暖機 cfl 0.1 × 4000, soft/mid cfl 0.2, 本段 12000): doe_001 dv を driver で (m10_on は m6_on から warm) | **PASS (degraded)**: C_T_w 0.9415 / C_M_w −3.00 / L 10.64。m6_on 0.9447 (p: 0.9520, C_M −5.11)、m10_on 0.9353 (p: 0.9418, C_M −3.47)、m4_off 0.9430 (p: 0.9483, C_M +2.96)、全 op で残差 stalled/plateau・力係数 STEADY。**m10_on の warm 標準は mid 段 (2 次, cfl 0.2) step 1386 でランプ後縁 (x/H 10.2–10.9, y/H 3.74–3.9 = EXH∩AIR 接触) NaN** → cold 緩レシピ (cfl 半分) で完走。run_0095 では同レシピが通ったので限界状態 | active |
| `run_0099_frozen_tp_doe001_recipe/` | **frozen_tp 生産レシピ確定 run**: soft/mid/暖機 cfl 0.1 (steps 4000), 本段 12000 (`problem_moo_frozen_tp_cycle3op.yaml` の最終値) で retry 無しに 3 作動点が通るかの確認 | **PASS (retry 無し, 463 s / 3 作動点)**: C_T_w 0.94152 / C_M_w −3.003 / L 10.64。m6_on 0.94470 (p 0.95196, C_M −5.109) / m10_on 0.93528 (p 0.94181, C_M −3.468) / m4_off 0.94295 (p 0.94834, C_M +2.958) = run_0098 と 5–6 桁一致。全 op 残差 stalled/plateau・力係数 STEADY。**これが frozen_tp の生産レシピ** (1 評価 ≈ 8 分) | active (ref) |
| `run_0100_r4b_corner_m6on/` | **R4b(i) A/B**: doe_001 frozen m6_on を `mesh.node_inlet_corner_wall: 1` (変換時 `nodeInletCornerWall`) で | GATES PASS。**却下**: 入口角の壁圧は 1.385 p_in (無し) → 0.909 p_in になるが、上流ダクト壁 x∈[−0.5,−0.15] 全体が 0.88–0.96 p_in の**過小** (無しは 1.005 = 正解に近い)。C_T −0.0023 / 摩擦込み −0.0033 / C_M +0.057 | 診断 (ref) |
| `run_0101_r4b_corner_m10on/` | 同上 m10_on (cold 標準レシピ) | **暖機段 step 843 で NaN** (x/H −0.48〜−0.11, y/H −0.03〜0 = カウル板下面の外気側, P 床 20 Pa)。角所有の変更で M10 板下面がさらに不安定 → 却下 | 診断 (ref) |
| `run_0102_r4_domain/` | **R4 領域独立性**: smoke 設計・加速点・3D Euler node・**3D ext_top** (機体上面テーパ + 自由流バンド, `r4_domain_study.py`): base / Z_ext 3 / x_out 4 / bot_depth 1.5 / top_depth 4 / 格子 1.25 倍 | (実行中 — `domain_study.md`) | active |
| `run_0103_r4b_isothermal_m10on/` | **R4b(ii) A/B**: doe_001 frozen m10_on を等温壁 1000 K (`spec.wall_thermal`) で。M10 板下面の断熱回復温度 4600 K と冷点 (T 床 50 K) が消えるか | **GATES PASS・残差 3.3–3.5 桁低下** (rms_ro 4.4e-7; 断熱は 5–7e-5 のプラトー = 100 倍改善)。T max 2320 K (= 排気 T, 回復温度の極値が消える)、冷点は 143 K / P 504 Pa に緩和 (残る冷点はカウル板 LE が入口面にある上流延長の産物)。C_T(p) −0.0017 / 摩擦込み −0.0037 (摩擦 +30 %) / C_M +0.11 vs 断熱 (run_0099)。→ **生産 YAML の既定を等温壁 1000 K に** (plan §8-11) | active (ref) |
| `run_0105_r7_doe8_frozen/` | **R7 小規模探索**: frozen_tp 生産 YAML (等温壁 1000 K, cfl 0.1 レシピ, 本段 12000) で LHS 8 点 (seed 1, infill 無し) の判別能力確認 (C_T_w の広がり vs 再現ノイズ 1e-6, ゲート通過率, 1 評価の時間) | (実行中 — `ledger.jsonl` / `pareto.json`) | active |

### S7 3D の現状 (2026-09-05)

- **外側空間なし** (側壁 = 境界壁、横方向膨張なし) は node SST で完走し、2D と力係数が 1 % 以内で一致 (run_0023)。
  3D の hex 生成・スリット・index IC・段階起動 (index コピー)・quad 力積分は動く。壁面出力 CONNE は 1 面 5 整数 [5, n0..n3]。
- **外側空間あり** (カウルと側壁の横端が外部流に露出) は 4 回とも soft 段 step 3–8 で NaN。入口タグ・共有ノードの 2 種入口・
  IC の 3 つの実バグを潰した後も、入口面直下・カウル横端 (z≈W/2) の外部流から発散する。横端 (カウル外面 ⟂ 側壁外面の凸角線、
  側壁∩入口線の壁+入口ノード) の node 境界処理が疑わしい。Euler (slip) で同じ IC の切り分けを実行中。
- **外側空間あり・Euler は成立** (run_0027, 加速点 M∞3.5, top_out slip): 3 つの実バグ (入口タグ / 2 種入口の共有ノード /
  IC) と M∞6 の幅外ランプ角部の真空膨張 (2 次で負圧) を避ければ回る。**横方向膨張の効果**: 側壁がカウル後縁で終わるため、
  その先で排気が幅外 (0.05 p_in 相当の低圧) へ横に逃げ、幅内ランプの圧力が 2D の ~0.3 p_in から ~0.1 p_in に落ちる。
  面別 (半スパン正規化): ランプ幅内 T −0.005 / L −0.019 (2D: +0.014 / +0.085)、幅外 T −0.025 / L −0.104、カウル内外は 2D と同じ。
  合計 C_T 0.932 (2D 0.976)、C_L −0.204 (2D +0.004)、C_M +1.07 (2D −0.001)。→ **短い側壁の SERN は 2D 設計値から推力 −4.5 %、
  揚力は符号反転** (この加速点で)。側壁長 L_sw を dv/仕様に持つ意味が大きい。
- **SST 3D は未解決**: 暖機 + 正しい壁距離 + top_out slip でも soft 段 step 686 でカウル後縁∩側壁後縁の 3 重点 (2 せん断層の交差、
  壁距離 0.4–3 mm) で ω → inf (run_0028)。候補: 後縁に板厚を持たせて 2 つの剥離線を分離、`FORGE_FREEZE_TURB=1` で soft 段を凍結乱流で
  通してから解放、3 重点近傍の ω 上限。
