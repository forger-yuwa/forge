# case/42.isobutane_wt — イソブタン燃焼ガス風洞 (semi-perfect, axis-Mach チェーン)

案件 (2026-08-17): $P_t$=5 MPa, $T_t$=1000 K, C₄H₁₀/空気 φ=0.9 生成物
(CO₂ 16.8 / H₂O 8.6 / O₂ 2.2 / N₂ 72.4 wt%), 出口径 1 m, 入口配管径 1 m, **$M_d$=4.2** (当初 4.0、同日訂正)。
$r_t$=0.12575 m (semi-perfect の $A/A^*(M4.2)$=15.809 から)。
計画: [plans/active/tooling-nozzle-semiperfect-gas.md](../../plans/active/tooling-nozzle-semiperfect-gas.md)。

**熱力学**: 設計 (MOC) は NASA-9/CEA の semi-perfect (frozen 組成、γ 1.31→1.38)。
CFD は forge の TP (単一擬似種 `species_db.yaml`, **`thermoHrefTemp: 298.15` 必須** — 絶対基準では
陰解法 Jacobian の χ_eos が桁違いになり軸で発散 [run_0021–0025 で特定]) で回る (run_0026)。
R/L スタディ (run_0012–0019) は特定前に CPG(γ\*) で実施した相対比較。

**M6 拡張 (2026-08-16)**: 同ガス組成・$P_t$=5 MPa で **$T_t$=1550 K・$M_d$=6・出口/入口 1 m** ($r_t$=53.75 mm、$A/A^*$=86.55)。
単一 quintic 軸 M 則は M6 で全点 REJECT (θ_w>μ_w の逆 MOC fold) → **knot 則 (A6)** で成立。全長 ≥ $L_c$+55 $r_t$ (≈5–7 m) は
終端特性線 ($\mu_d$=9.6°) の物理的必然。計画: [tooling-nozzle-axismach-knot-law.md](../../plans/accepted/tooling-nozzle-axismach-knot-law.md)。

## M4.2 ベスト結果の参照・引き継ぎ (他マシン向け, 2026-08-17)

**ベスト = `run_0035_ib_R2_LU6_Lc8_tp` (R2 / L_U6 / L_c8, ‖ΔM‖∞ 0.084 % $M_d$)**、次点
`run_0036_ib_R3_LU6_Lc8_tp` (R3, 0.091 %, ε_M 0.009 % で出口一様性は最良)。差は入口 BC 床と
同程度で **R2/R3 は同等**と読む。CFD 成果物 (`res_*.h5`, メッシュ h5/msh, log) は git 対象外だが、
この 2 run だけ**入力・設計成果物をリファレンスとして git に残してある**:

| ファイル | 内容 |
| --- | --- |
| `problem_ib_R2_LU6_Lc8_tp.yaml` / `problem_ib_R3_LU6_Lc8_tp.yaml` | 問題定義 (設計変数・ガス・評価設定)。**これだけで全て再生成できる** |
| `run_0035_…/wall_design.csv` | 逆 MOC 壁テーブル [x_m, r_m, θ_rad, M_wall] (m 単位)。CAD/別ソルバへ渡す形状はこれ |
| `run_0035_…/target_axis_M.csv` | 軸 Mach 目標 M(x) |
| `run_0035_…/solverConfig.yaml` / `bcondConfig.yaml` / `probe.yaml` / `species_db.yaml` | forge 入力 (単一擬似種 TP, `thermoHrefTemp: 298.15`, 本段 cfl2) |
| `run_0035_…/prepare_info.json` / `metrics.json` / `achieved_vs_target.csv` | 設計メタ (x_A, x_E, アンカー, メッシュ ni/nj) と CFD 結果指標 |
| `run_0035_…/MESH_QUALITY.txt` / `CONVERGENCE_VERDICT.txt` | 品質・収束の記録 |
| `study_cfd_m42_tp.json`, `run_study_m42_tp.py` | 8 点スタディの表と投入ドライバ |

**再生成 (別 PC)**: forge をビルド (`solver_density_cuda/`, [procedures/development-environment.md](../../procedures/development-environment.md)) し、
`design/` に venv を作ってから
```
cd design && ./.venv-opt/bin/python -m forge_design.evaluate.runner_axismach \
    ../case/42.isobutane_wt/problem_ib_R2_LU6_Lc8_tp.yaml ../case/42.isobutane_wt/run_00NN_ib_R2_LU6_Lc8_tp [--prepare-only]
```
`--prepare-only` で設計 + メッシュ (`nozzle.msh`/`nozzle.h5`) + config だけ作れる (メッシュは TFI 構造化を
`design` 側で直接書くので Gmsh 不要、変換に forge の `convertGmshToForge` を使う)。CFD は 12000 step ≈ 20 s
(RTX 3060)。**現行 binary は TP 入口 BC 修正済み** (`run_0050` で残差床 3e-2→9e-7)。run_0031–0038 は修正前の
binary の結果なので、再実行すると残差は深く落ちるが軸 M/出口指標は ~0.01 % 以内で同等 (run_0033 vs run_0050 実測)。

**凝縮計算へ引き継ぐときの注意**:
- **凝縮は `evaluate.tp_species: split_h2o` で回せる (2026-08-17 実装・検証済)**: H₂O 以外を擬似種 `MIXDRY`、H₂O を
  独立種にした 2 種 TP + `evaluate.condensation: {condensation: 1, nCondSpecies: 1, condModel: 1, condKantrowitz: 1, condGrowthModel: 0}`
  (`condGasSpecies` は自動)。実例 `problem_ib_R2_LU6_Lc8_h2o5_split_cond.yaml` → `run_0065`。dry では pseudo と同一
  (run_0063 vs 0064)。Wyslouzil (case/16) で既存 CPG 凝縮結果を再現済 ([計画](../../plans/accepted/tooling-nozzle-tp-split-h2o-condensation.md))。
  多成分 TP × 陰解法は cfl ≤2 (段階起動、`thermoHrefTemp` 298.15) で安定。
- 形状・メッシュ・IC は `problem_ib_R2_LU6_Lc8_tp.yaml` から `--prepare-only` で再生成し、`solverConfig.yaml` の
  `physProp`/`condensation` ブロックだけ差し替える。IC は `paste_isentropic_ic` (設計 1D 等エントロピー) が既定。
- 出口静温 ~247 K・静圧 ~20–28 kPa (M4.2, Tt 1000 K) — H₂O 分圧はここから。

## M4.2 ベスト形状の数式・点列・gmsh (2026-08-17)

**Euler 版** (`run_0035` が実際に解いた非粘性設計壁) と **NS 版** (物理壁 = 非粘性壁 + 排除厚 δ\* の
法線オフセット、A13 `PhysicalNozzleWall`) の 2 つを出力する。実機の壁 (境界層込みで設計コアを実現する形状)
は **NS 版**、CFD/理論の基準形状は Euler 版。

再生成: `design/.venv-opt/bin/python case/42.isobutane_wt/export_wall_m42.py` →
`points_m42_best_{euler,ns}.csv` (全輪郭点列 [m] + 無次元)、`nozzle_m42_best_{euler,ns}.geo` (gmsh)。
gmsh: `gmsh -2 nozzle_m42_best_euler.geo -format msh41 -o nozzle.msh` (軸対称半平面・構造化 quad 365×64、
physID は forge 規約 inlet 1/outlet 2/wall 3/axis 4/fluid 5。両版とも `check_mesh_quality.py` PASS)。

### 主要寸法 [m]

| | Euler 版 | NS 版 (δ\* オフセット) |
| --- | --- | --- |
| スロート | x=0, r=0.125750 (φ0.25150) | x=−0.001023, r=0.126818 (φ0.253636) |
| 出口 | x=3.094893, r=0.499027 (φ0.998054) | x=3.094893, r=0.511326 (φ1.022652) |
| 入口 (直管) | x=−0.817375, r=0.499982 (φ1.0) | 同左 |
| 収縮開始 | x=−0.754500 (=−6 r_t) | 同左 |
| 全長 / 拡大部 | 3.9123 / 3.0949 | 同左 |
| 最大壁角 | 21.24° @ x=0.150 | ほぼ同じ |
| δ\* (壁オフセット) | — | スロート 1.07 mm → 出口 12.3 mm |

$r_t=0.12575$ m、$R$ (スロート曲率半径) $=2\,r_t=0.25150$ m、$L_U=6\,r_t$、$L_c=8\,r_t$、$M_d=4.2$。

### 数式 (Euler 版、無次元 $\hat x=x/r_t$, $\hat r=r/r_t$)

1. **入口直管** $-6.5\le\hat x\le-6$: $\hat r=3.976$ (φ1 m)。
2. **収縮 U→T (5 次 Hermite)** $-6\le\hat x\le0$、$\xi=(\hat x+6)/6$:
   $$\hat r(\xi)=3.976-20.76\,\xi^3+26.64\,\xi^4-8.856\,\xi^5$$
   端条件 $\hat r(0)=3.976,\ \hat r'=\hat r''=0$ / $\hat r(1)=1,\ \hat r'=0,\ \hat r''=1/R=0.5$ を厳密に満たす
   (単調収縮の判定量 $\mu=L_U^2/(R\,\Delta r)=6.05\le20$)。
3. **スロート下流 (逆 MOC 設計壁)** $0\le\hat x\le24.6115$: 逆 MOC の壁テーブル
   (`run_0035_ib_R2_LU6_Lc8_tp/wall_design.csv`, 129 点 [x, r, θ, M_wall]) を**端条件クランプ付き
   補間 5 次 B-spline** で表現 (スロート端で $\hat r'=0,\ \hat r''=1/R$、$C^2$ 接続)。閉形式は無いので
   点列 (`points_m42_best_euler.csv`, 1028 点) か上記スプラインを使う。

NS 版は 3. の壁を δ\*(s) だけ法線オフセットし、真の幾何スロート ($r'=0$ の点) を再探索して上流 Hermite を
その $(x_t, r_t, \kappa_t)$ で作り直したもの (数式は同形、係数が変わる)。

### 代表点 (Euler 版, [m])

| x | −0.8174 | −0.6288 | −0.3773 | −0.1258 | 0 | 0.1258 | 0.2515 | 0.5030 | 1.0060 | 1.7605 | 2.2635 | 3.0949 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| r | 0.49998 | 0.49034 | 0.34823 | 0.15723 | 0.12575 | 0.15686 | 0.20494 | 0.28674 | 0.39481 | 0.47286 | 0.49292 | 0.49903 |

## M5 推奨形状 (R3/L_U9/L_c14) の点列・gmsh (2026-08-30)

設計経緯 8 フェーズの総括と NS 軸 M 診断 (δ* 面積誤差 / 抽出信頼域 / 軸集束 / law 側帰還欠如):
[M5 設計報告ページ](https://claude.ai/code/artifact/22b65472-545e-4a05-aa7c-ac064af2e424)。

M5 スタディ (③) のパレートエルボ **R3/L_U9/L_c14** (`run_0096`, dM 0.043 %$M_d$) の壁輪郭を
`export_wall_m42.py --problem problem_ib_m5_R3_LU9_Lc14.yaml --tag m5_best` で出力
(スクリプトは `--problem` で任意の axismach 点に使える汎用形。ヘッダの M 表記も M_d 連動に修正済)。

| ファイル | 内容 |
| --- | --- |
| `points_m5_best_euler.csv` / `nozzle_m5_best_euler.geo` | **非粘性設計壁** (= run_0096 が解いた形状)。CFD/理論の基準 |
| `points_m5_best_ns.csv` / `nozzle_m5_best_ns.geo` | **物理壁 v3** (`PhysicalNozzleWall` A13 = 非粘性壁 + **CFD 抽出 δ\* (run_0100 v1 場から抽出・平滑化)** の法線オフセット, `--dstar-csv`)。**実機の壁はこちら。run_0102 で NS/SST 直接検証済み** (壁一致 4 μm) |

物理壁 v3 の要点: 幾何スロート半径 **81.30 mm** (設計 80.505 mm、0.45 mm 上流)、出口半径 499.5 → **522.7 mm**
(境界層排除厚を吸収してコア径 ~1 m を維持)。**NS/SST y+~1 直接検証 (run_0102)**: 出口コア (r≤0.44 m)
mean M **4.985 (−0.29 % vs M_d)**・ε_M_rms 0.10 %・δ99≈50 mm・有効コア径 ≈0.94 m — 0.5 % ゲート内。
δ\* は 2 巡目抽出で固定点 (new/prev 1.001±0.005)。相関 v1 壁は出口コア −0.71 % でゲート外だった
([δ\* plan](../../plans/accepted/tooling-nozzle-axismach-viscous-deltastar.md) の確立フローどおり v3 で締めた)。
残る軸 M うねり ±0.5 % は δ\* 非表現 (law 側帰還が残件、同 plan)。

## 計算 run 一覧

| run | 目的・設定差分 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `run_0001_ib_base` | **TP (semi-perfect) CFD の初回投入** (`problem_ib_m4.yaml`, 5 回試行: interp_field CPG 再構成バグ → 背圧 1000→28000 Pa → cfl4→mid 段+cfl2 → axisymMethod 1 → method 0 + axisRFloor 2.8 mm) | 全て発散。最終形は soft 段 step 271 で膨張部中域 (x≈5–6) の軸近傍から。**同一条件の CPG (run_0002) は完走** → **TP × node 軸対称は forge 側の穴** (nodeAxisDirichlet が TP 不可、代替なし)。別セッション (軸対称高度化) へ申し送り | 削除済み(2026-08-31) (**forge 側問題の記録**) |
| `run_0021_ib_tp_constcp` / `run_0022_ib_tp_lowTfreeze` / `run_0023_ib_tp_explicit` / `run_0024_ib_tp_hoopjacfix` / `run_0025_ib_tp_href` | **TP 発散の切り分け系列** (run_0020 の soft 段設定 [1次 cfl0.5 3000 step] を 1 要素ずつ変える): 0021 = 一定 cp 擬似種 (γ 1.309)、0022 = 実 NASA-9 で 200 K 未満 cp 凍結、0023 = 実 NASA-9 陽解法 (timeIntegration 3)、0024 = 軸ホップ Jacobian を一般 EOS の ∂P/∂Q (χ_eos 項) に修正した binary、0025 = **`thermoHrefTemp: 298.15` (sensible-enthalpy datum) + IC を同 datum** | 0021 完走 / 0022 step 191 発散 / 0023 完走 / 0024 step 190 発散 / **0025 完走**。→ **原因 = 陰解法 Jacobian の χ_eos = c² − κh に絶対基準 h (生成エンタルピー込み, −1.8e6 J/kg) が入り桁違いになる**。低温外挿・TP 経路自体・ホップ Jacobian は無罪 (ホップ修正は正しいので残置)。打ち手 = runner の TP 経路に thermoHrefTemp 298.15 を既定化 + IC を同 datum で構成 | 削除済み(2026-08-31) (**原因特定の記録**) |
| `run_0028_ib_R3_tp_inletdir` / `run_0029_ib_R3_tp_pipe2` | **TP 残差床 (rms_ro 3e-2 プラトー) の切り分け**: 場は凍っているが入口面で P/Pt=1.03 (全圧超過) が振動 = TP `inlet_Pressure` (静圧/速度参照 rf=0.5 緩和) の病理。0028 = `inlet_Pressure_dir` に変更、0029 = 配管長 0.5→2 r_t | **0028: step 9059 で発散** (更に悪い)。**0029: 配管 2 r_t でも P/Pt 1.03 振動・残差床 3.6e-2 不変** (内部 0.114%) → **TP の入口 BC 2 種とも大径低速入口 (M≈0.03) で問題**。CPG 入口 BC (同 case 完走・rms_ro 1e-9) にはない → **forge 側 TP 入口 BC の申し送り**。ノズル内部 (軸 M・出口) は凍っており設計評価には使えるが「収束」ではない | 削除済み(2026-08-31) (**申し送りの根拠**) |
| `run_0030_ib_tp_lowTfreeze_href` | datum 修正済みで **200 K 未満 cp 凍結**の効果を再検証 (run_0027 と同 IC・6000 step) | 発散なし、軸 M 差 3e-3、出口 M 同一 (4.2011/4.2010)。収束場の T min 247 K で外挿域に触れない。run_0022 の「凍結でも発散」は datum 未修正下の結果で凍結の判定ではなかった (訂正) | 削除済み(2026-08-31) (記録) |
| **`run_0031`〜`run_0038_ib_*_tp`** | **R/L スタディ TP 正本** (semi-perfect 設計 + TP CFD, thermoHrefTemp 298.15, `problem_ib_*_tp.yaml`; CPG 代替 run_0012–0019 と同 8 点) | **全点 0.5% ゲート内、L_U≥6 は 0.08–0.15%**: R1.5/2/3 Lc=max 0.108/0.117/0.115%、**Lc=8 が最良級 R2 0.084% / R3 0.091%** (CPG の「Lc=8 は R≥2 でゲート外」は熱力学不整合の産物 → 撤回)、L_U4 は 0.386% と TP でも明確に悪い (出口 M 4.216 過大)、L_U9 0.102%。ε_M 0.009–0.028%、ε_θ 0.004–0.022°。全 run 品質 PASS・NaN 0。残差は TP 入口 BC の床 (rms_ro 3e-2、軸 M 変動 3–4e-3) で **VERDICT は NOT CONVERGED** — ~0.1% 以下の設計点間差は床と同程度で有意でない。`study_cfd_m42_tp.json`、[結果ページ](https://claude.ai/code/artifact/07b52b7f-54b1-4cc8-97a5-4ee9af380b2d) §3b | active (**TP スタディの根拠・推奨の正本**) |
| **`run_0026_ib_R3_tp_full`** / `run_0027_ib_R3_tp_cont` | **TP (semi-perfect 設計 + semi-perfect CFD) フル run** (`problem_ib_R3_LU6_Lcmax_tp.yaml`, thermoHrefTemp 298.15, soft→mid→本段 cfl2 12000; 0027 = 継続 +24000) | **完走。‖ΔM‖∞ 0.116% Md・出口 ε_M 0.014%・ε_θ 0.0036°・overshoot +0.02%** — CPG(γ*) 検証 (0.256%) の半分以下 (熱力学が設計と完全整合)。出口軸 M 4.2010、コア M 4.2005、**出口径 0.999 m** (目標 1.0)。品質 PASS・NaN 0。残差は rms_ro 3e-2 プラトー・軸 M 8k→12k で 2.9e-3 とまだ動く (CPG 系列の 8e-6 より深さ不足) → 0027 で継続確認 | active (**TP CFD の基準**) |
| `run_0020_ib_R3_tp_newbin` | **TP (semi-perfect) CFD の再挑戦 — 新 binary** (`ea08bcbe`: nodeAxisDirichlet 撤去、軸ノードが真の DOF)。`problem_ib_R3_LU6_Lcmax_tp.yaml`, cfd_gas 無し・axisRFloor 無し | **再び発散** (soft 段 step 190、x≈7.0–7.5 の軸上で P が床 1.46 Pa まで落下)。IC は健全 (T 247–1000 K)。旧 binary の method 0 + Dirichlet 0 (step 271, x≈5–6 軸) と同じ性質 → **TP × node 軸対称の発散は nodeAxisDirichlet と無関係、forge の TP 経路 (EOS/軸ソース項) 固有**。同一メッシュ・IC の CPG は新旧 binary とも完走。申し送り継続 | 削除済み(2026-08-31) (**forge 側問題の記録・第 2 弾**) |
| `run_0002_cpgcheck` | 同一壁 (semi-perfect 設計) を **CPG(γ\*) で CFD** (`problem_ib_cpgcheck.yaml`) — TP×軸対称の切り分け | **完走**: ‖ΔM‖∞ 0.284% $M_d$・出口 ε_M 0.013%/ε_θ 0.0096°。品質 PASS、軸 M は 4000 step 間 8.8e-6 で凍結 (残差は warm 床プラトー)。CFD 経路の基準 | active (**CPG 経路の基準**) |
| `run_0003`〜`run_0011_ib_R*_LU*_Lc*` (初回, 破棄) | 設計 semi-perfect × CFD CPG(γ\*) の**熱力学不整合** run | 全 9 点で ΔM≈4.5%・ε_M 3.5%・出口コア M 3.86 — 壁 (A/A*=13.1) と CFD (CPG 1.309 なら 15.3) の食い違いで M4 に届かない。**設計と CFD の熱力学は必ず一致させる**教訓。同名で再投入 (下) | 破棄済 (記録のみ) |
| `run_0003`〜`run_0011_ib_R*_LU*_Lc*` (2 回目) | 整合 CPG だが **M4 と M4.2 が混在** (途中でユーザ訂正 4.0→4.2 を YAML に反映したため 0003–0006 が M4、0007–0011 が M4.2) | 個々の値は健全 (0.30–0.68%) だが系列として使わない | 削除済み(2026-08-31) (混在の記録) |
| **`run_0012`〜`run_0019_ib_R*_LU*_Lc*`** | **R/L スタディ CFD (M4.2 統一, 設計・CFD とも CPG γ\*=1.309, `cfd_gas: cpg`)**: R∈{1.5,2,3}×L_c∈{max,8} @L_U=6、L_U∈{4,9} @R2 Lc=max。CPG では出口径 1.09 m (semi-perfect 正本 1.0 m) — 相対比較用 | **R3/L_U6/L_c=max が最良 (‖ΔM‖∞ 0.256%, ε_θ 0.003°)**。L_c=8 は R≥2 でゲート外 (0.52/0.68%)。L_U は軸 M に効かず ε_M に効く (L_U4 で 0.061% vs L_U9 で 0.010%)。全 run 品質 PASS・NaN 0・軸 M 凍結 ~1e-5。`study_cfd_m42.json`、[plan §3.5](../../plans/accepted/tooling-nozzle-semiperfect-gas.md)、[結果ページ](https://claude.ai/code/artifact/07b52b7f-54b1-4cc8-97a5-4ee9af380b2d) | active (**スタディの根拠**) |
| **`run_0051_ib_m6_R3_Lc45_MK2.5_tp_inletfix`** / **`run_0061_ib_m6_R3_Lc50_MK2.5_tp`** / `run_0052`〜`run_0060_ib_m6_*_tp` | **M6 スタディ (Tt 1550 K, Pt 5 MPa, 出口/入口 1 m, r_t 53.75 mm; `problem_ib_m6*.yaml`; knot 軸 M 則 A6 + `axis_dx0` 0.03 + `inlet_Pressure` TP 修正 binary)**: R3 × L_c∈{30,40,45,50,60} @M_K 2.5; R∈{2,5} @L_c50; M_K∈{2,3} @R3/L_c50; L_U∈{8,16} @R3/L_c50 (基底 12)。0051 のみ 36000 step (定常性確認)、他 12000 | **全 11 点 0.03–0.08% $M_d$ (ゲート 0.5% の 1/6 以下)**、ε_M 0.002–0.017%、ε_θ 0.002–0.007°、M_exit 6.001–6.004、overshoot ≤0.07%。傾向: L_c 30→60 で θ_max 19.7→11.2°・dM 0.061→0.035%・全長 5.3→6.9 m; R2 が最良 (0.029%)、R5 も可 (0.040%); M_K 2/2.5/3 は同等 (0.040%、M_K2 は ε_M 最小 0.002%); **L_U8 は悪化 (0.078%, ε_M 0.017%)**、L_U16 は L_U12 と同等 (0.049% だが軸 M がまだ 1.8e-3 動く)。全 run 品質 PASS・NaN 0。残差 rms_ro 1e-6〜3e-5 (soft/mid 段で落ち切った warm 床、`check_convergence` は NOT CONVERGED 判定)、設計区間の軸 M は 8k→12k で 2–5e-5 に凍結 (0051 は 36k まで 3e-5)。`study_cfd_m6.json` / `study_design_m6.json` (設計 36 点全合格)、[knot plan](../../plans/accepted/tooling-nozzle-axismach-knot-law.md) | active (**M6 スタディの正本**) |
| `run_0039_ib_m6_R3_Lc45_MK2.5_tp` / `run_0040_..._cont` / **`run_0050_ib_R3_LU6_Lcmax_tp_inletfix`** | **TP `inlet_Pressure` BC の発散切り分けと修正**: 0039 = M6 初回 (旧 binary, 12000 step で 0.041% だが設計区間の軸 M が 2e-3 動く)、0040 = 継続 +24000 (cfl2) で**入口 P/Pt が 1.04→1.14→1.87→5.1→6.2 と指数成長し M_exit 6.21 まで場が壊れる** = TP 入口 BC (静圧参照ブレンド rf 0.5) の不安定が M_in≈0.011 で顕在化。**修正 (rf→0, CPG 分岐と同じ速度参照)** 後の M4.2 再実行 0050: 残差床 3e-2 → **9e-7**、入口 P/Pt 0.9986–0.9993、‖ΔM‖∞ 0.138%・ε_M 0.025% (0033 と同等) | 0039/0040 削除済み(2026-08-31) (**発散の記録**)、0050 active (**修正の根拠**) |
| **`run_0062_ib_m6_R3_Lc45_bsplineM_tp`** | **軸則比較 B (monotone B-spline) の CFD** (`problem_ib_m6_R3_Lc45_bsplineM.yaml`, R3/L_c45/M_d6, A [`run_0051`] と同一条件で `axis_law: bspline_M` のみ差し替え; [軸則比較 plan](../../plans/accepted/tooling-nozzle-axislaw-smoothness.md)) | **完走 (12000 step, PASS/NaN 0)**。‖ΔM‖∞ **0.183%** (A 12k 相当 0.042% の ~4.4 倍)、ε_M 0.023%、ε_θ 0.013°、M_exit 6.0035。B は $M'''$ ジャンプを消す (A の 0.30→B の 0.06) が、その代償として θ_max が 23.5° (A 14.1°) まで増え、`min(μ_w−θ_w)` 余裕が 0.09° (A 1.78°) とほぼ限界まで縮む — **軸則の見た目の滑らかさは壁品質の改善に直結しない**、が結論。C (非負 dν/dx B-spline) は前倒れ対策後も内部特性線網に恒常的な向き反転 (322–496 箇所、解像度で悪化) と壁 spline リンギング (~1°) が残り hard gate 不合格 (詳細は比較 plan)。詳細比較図・表: `compare_axislaw_ABC.json`/`.py`、[結果ページ](https://claude.ai/code/artifact/cf8a4c74-8d1f-44f4-a256-652cc122d00d) | active (**軸則比較 CFD の根拠**) |
| **`run_0063_ib_R2_LU6_Lc8_h2o5_pseudo`** / **`run_0064_…_split`** / **`run_0065_…_split_cond`** | **M4.2 ベスト形状の H₂O 5 % 版 + 凝縮引き継ぎ** (`problem_ib_R2_LU6_Lc8_h2o5_{pseudo,split,split_cond}.yaml`, 組成 N2 0.7528/CO2 0.1743/O2 0.0229/H2O 0.0501; node Euler TP, 段階起動 12000 step; [plan](../../plans/accepted/tooling-nozzle-tp-split-h2o-condensation.md)) — 0063 = 単一擬似種 MIX、0064 = `tp_species: split_h2o` (`[MIXDRY,H2O]`, dry)、0065 = split + H₂O 凝縮 Kw+HK (`condGasSpecies 1`) | **pseudo と split dry は同一** (‖ΔM‖∞ 0.062 %/ε_M 0.0117 %/ε_θ 0.0194° が完全一致、軸 M 差 max 1.7e-5)。**凝縮 ON**: onset x≈6.7 r_t (T≈250 K)、軸 g_max **0.0126** (H₂O の 25 %)、潜熱で出口軸 T 245→278 K・P +4.5 %・**出口軸 M 4.20→3.94** (‖ΔM‖∞ 2.2 %、ε_M 7 % — 凝縮域が半径方向に不均一)。全 run 品質 PASS・NaN 0。残差は warm 床 (rms_ro ~1e-6, NOT CONVERGED 判定; 0065 は 2.7 桁降下)。比較図 `compare_h2o5_split.py` → `compare_h2o5_split.png`/`.json` | active (**凝縮引き継ぎの正本**) |
| **`run_0066_ib_R2_LU6_Lc8_h2o5_split_cond_evap`** / `run_0067_…_split_cond_noevap_ctrl` | run_0065 の `res_12000.h5` から同一メッシュ restart 12000 step ([蒸発 plan](../../plans/accepted/condensation-evaporation.md)) — 0066 = `condEvaporation: 1`、0067 = 蒸発 OFF (restart 対照)。`interp_field.py` は凝縮モーメント dataset を新規作成するよう修正 | **蒸発が想定どおり作動**: 0065 では試験部圧縮帯 (x≈2.6–3.07 m, y≈0.27–0.40 m) に S 0.71–0.95 かつ g≈1.2 % が凍結 (514 ノード, 液相の 9.2 %) していたが、0066 では S_min 0.67→0.87、帯内 T 最大 293.2→289.1 K (−4.4 K)、g −10 % (帯平均 1.170→1.143 %)、S→1 へ平衡蒸発。上流 (x<2.5) は不変 (|ΔT| 4e-3 K)。診断 `condDrdt_0` は HK 式と 0.7 % 一致、λ 律速非発火 (λ≥0.86)。**step 4000 で定常・4000/8000/12000 同一 (STEADY)**。0067 は 0065 と同一。両 run NaN 0、残差 warm 床 (rms_ro 9e-7, NOT CONVERGED 判定 = 0065 と同じ) | active (**蒸発検証の正本 (Euler)**) |
| `run_0068_ib_R2_LU6_Lc8_h2o5_split_ns_dry` / **`run_0069_…_ns_cond_noevap`** / **`run_0070_…_ns_cond_evap`** / `run_0071_…_ns_cond_evap_default` | **H₂O 5 % split の NS 版** (`problem_ib_R2_LU6_Lc8_h2o5_split_ns.yaml`: 物理壁 δ\* + no-slip 断熱壁, node y+~1 SST, 561×97 wall_first 4.5e-5, `prepare_ns`+`run_staged_ns` soft/mid/本段 24000, IC=run_0064) — 0068 = dry NS (品質 PASS AR 752/skew 0.45)。0069/0070/0071 = 0068 res_24000 から凝縮 ON restart 12000 (0069 蒸発 OFF, 0070 `condEvaporation: 1` 明示, 0071 既定 [=1]) — [蒸発 plan](../../plans/accepted/condensation-evaporation.md) §6-3 高温壁デモ | 壁 T≈920 K (回復温度)。**蒸発 OFF (0069)**: T>400 K のノードに液相 43 個 (g≤2.5e-6)・液相質量の 46 % が S<1 (S 中央値 0.89, 最小 0.007) に凍結。**蒸発 ON (0070)**: T>400 K の液相 **0**、S<1 の液相質量 35 % に減り S 中央値 0.997 (平衡蒸発)、最大 ΔT −11.2 K (出口 x=3.08 m, y=0.34 m: g 1.28→0.97 %)。壁列 (wall_dist<0.1 mm) は OFF でも液相 1.9e-24 = 液滴は境界層に入らない (モーメント無拡散; 境界層流体は上流でも高温)。0071 は 0070 と同一 (|ΔT| 1.5e-2 K)、既定化後は g=0 モーメント塵も掃除 (Q1²>Q0Q2 違反 1247→63, S<1 で 0)。4000/8000/12000 で S<1 液相 35.0/35.0/35.1 % (STEADY)。全 run NaN 0、残差 2.6 桁 warm 床 (NOT CONVERGED 判定, 0068 と同じ) | active (**蒸発検証 (NS 高温壁)**; 0069 は対照として保持) |
| **`run_0072_ib_m5_R2_LU6_Lc17_tp`** | **M5 スタディ (③) の smoke/センター点** (`problem_ib_m5_R2_LU6_Lc17.yaml`: M_d5, Tt 1550 K, r_t 80.505 mm, knot MK2.5 + axis_dx0, R2/L_U6/L_c17, split_h2o dry TP; [design plan](../../plans/accepted/design-isobutane-wt-m5-sweep.md)) — M5 は quintic 全点不成立 (リンギング/許容窓)・knot L_c≥14 が成立域 (prepare-only プローブ 20 点) | **完走**: dM_max **0.083 %M_d**・出口 ε_M_rms 0.026 %・M_core 5.0013・overshoot 0.072 %。品質 PASS (AR 10.4/skew 0.64)・NaN 0・check_quasisteady ALL STEADY。残差 rms_ro 1.5e-6 の warm 床 (NOT CONVERGED 判定 = M4.2/M6 系列と同じ既知パターン)。Tt 1550 K で出口 T 308 K = H₂O 飽和線 +30 K (dry 妥当) | active (**M5 センター点・smoke の根拠**) |
| **`run_0073`〜`run_0098_ib_m5_R*_LU*_Lc*_tp`** | **M5 スタディ本体 (③, 2026-08-30)**: R {1.5,2,3} × L_U {4,6,9} × L_c {14,17,20} の 27 点 (センターは run_0072 流用; `run_study_m5_tp.py`, 各 `problem_ib_m5_<tag>.yaml`; [design plan](../../plans/accepted/design-isobutane-wt-m5-sweep.md)) | **27/27 完走・NaN 0・全品質 PASS・残差床 8.5e-7〜3.2e-5 (warm 床)**。パレート (全長 [亜音速込] vs dM) 9 点: 短尺 `run_0076` (R1.5/LU6/Lc14) 4.133 m/0.123%、**エルボ推奨 `run_0096` (R3/LU9/Lc14) 4.366 m/0.043%**、品質側 `run_0088` (R2/LU9/Lc17) 4.612 m/0.035%/ε_M 0.006%。**L_U4 は M5 でも悪い** (0.6–0.9%, R3_LU4_Lc14 は 0.5% ゲート超)。勝者 3 点 ALL STEADY・**Tsat_post 余裕 +29 K/S_max 0.16 = 凝縮なし・dry 確定**。集計 `study_cfd_m5_tp.json` / `study_cfd_m5_tp_pareto.json`、[結果ページ](https://claude.ai/code/artifact/f1b68d46-b8de-4d04-aa75-cab7bd342bb6) | active (**M5 スタディの正本**) |
| `run_0099_ib_m5_R3_LU9_Lc14_ns_coarse` / **`run_0100_ib_m5_R3_LU9_Lc14_ns`** / `run_0101_..._ns_v2` / **`run_0102_..._ns_v3`** | **M5 推奨点の物理壁 NS/SST 直接検証** (`problem_ib_m5_R3_LU9_Lc14_ns{,_coarse}.yaml`, `run_ns_m5.py`; A13 物理壁 + no-slip 断熱壁 + y+~1 SST, coarse y+~50 中継→本計算 800×97, 24000 step): 0099 = coarse 中継 (IC=run_0096)、0100 = **v1 相関 δ\***、0101 = 生抽出 CSV 全域採用 (失敗)、0102 = **v3 平滑化 CFD 抽出 δ\*** (`build_dstar_v3_m5.py --x-lo 12 --x-hi-trust 42`, δ\* plan の確立フロー) | 全 run 品質 PASS・NaN 0・コア量凍結 ≤3e-5。**v1 (0100): 出口コア M 4.965 = −0.71 % でゲート外** (相関が中流 1.29 倍過大→出口 0.94 倍過小)。0101 は壁波打ちで軸 M −5 % ディップ ([δ\* plan](../../plans/accepted/tooling-nozzle-axismach-viscous-deltastar.md) の「初回 v2 の失敗記録」の再演 = 生 CSV を使うな)。**v3 (0102): 出口コア M 4.985 (−0.29 %)・ε_M_rms 0.10 %・δ99≈50 mm、2 巡目抽出で固定点 1.001 → 検証完了・実機壁の正本** (points_m5_best_ns.csv は 0102 壁と 4 μm 一致)。M5 は抽出信頼下限 x_lo=12 (x<12 はコア一様性不足で抽出破綻、case/44 の x≥3 とは case 依存) | 0099/0100 active (中継・v1 記録)、0101 削除済み(2026-08-31) (失敗記録)、**0102 active (NS 検証の正本)** |
| `run_0103_ib_m6_R2_Lc50_ns_coarse` / `run_0104_ib_m6_R2_Lc50_ns` / **`run_0105_ib_m6_R2_Lc50_ns_v3`** | **M6 推奨点 (R2/L_U12/L_c50/MK2.5) の物理壁 NS/SST 直接検証** — M5 (run_0099–0102) と同フロー (`problem_ib_m6_R2_Lc50_ns{,_coarse}.yaml`; `run_ns_m5.py`/`build_dstar_v3_m5.py` は `--problem` 汎用化で流用)。メッシュ確立に 2 回失敗: ni1100/wf4.5e-5 は AR1538 FAIL、**ni1700/wf4.5e-5 (第一セル 2.4 µm) は median-dual 閉性破綻 (`dual faces not closed` 0.099, 壁面積 3 % 欠損) = converter 側の極薄セル問題 (要申し送り)** → ni1250/wf6.5e-5 (第一セル ~3.5 µm, y+~1.4) で成立 | 全 run 品質 PASS・NaN 0・コア量凍結 ≤1.5e-4。**v1 相関 δ\*: 出口コア M 5.966 (−0.57 %) でゲート外級** (M5/case44 と同パターン)。δ\* 抽出信頼下限 **x_lo=14** (x=10 は比 1.77 の汚染)。**v3: 出口コア M 5.974 (−0.43 %)・ε_M_rms 0.04 %、2 巡目抽出で固定点 0.998 [0.976, 1.004] → M6 の到達点** (0.5 % ゲート内ぎりぎり)。設計区間 dM は v1 0.38 %→v3 0.89 % (うねり化、M5 と同構造 = law 側帰還の領域) | 0103/0104 active (中継・v1)、**0105 active (M6 NS 検証の正本)** |
| `run_0106_m5_R3_LU9_Lc14_stepscan` | **評価 step 数の過剰確認 (2026-09-04)**: run_0096 (M5 エルボ) を同設定で再実行し、本段の出力間隔を 500 に詰めて軸 M 指標の step 依存を測定 (soft 3000 + 本段 12000) | 設計+メッシュ 3.8 s / soft 11.3 s / 本段 40.7 s。**dM_max は本段 3500 step 以降 0.043±0.005 %Md のリミットサイクル**、max\|M−M_12000\| は 5000 step で 0.003 %Md → 本段は 5000–6000 step で十分 (1 評価 65 s → ~35 s の見込み) | 破棄予定 (一時 run) |
| `run_0106_ib_m5_R3_LU9_Lc14_tp_geomfix` / **`run_0107_ib_m6_R2_Lc50_ns_yp1`** | **2D 双対幾何+CW 判定桁落ち修正の検証** ([plan](../../plans/accepted/discretization-median-dual-2d-facevect-precision.md)): 0106 = M5 Euler センター点を修正後 converter の新変換で再実行 (回帰)、0107 = 修正で初めて変換可能になった **y+~1 真値** (第一セル 2.4 µm, ni1700/wf4.5e-5, `problem_ib_m6_R2_Lc50_ns_yp1.yaml`) の M6 NS v1 | **回帰: 完全一致** (dM 0.000430 vs run_0096 0.000429, M_axis_exit 差 7e-6, ε_M 同一)。**payoff: 出口コア M 5.9640 (−0.60 %)** = y+1.4 版 run_0104 の 5.9660 (−0.57 %) と 0.03 % 整合 — y+1.4 退避は結果を歪めていなかったと確認、以後 y+1 真値が使える。閉性は全メッシュ normalized ~1e-7 (修正前 1e-4〜0.099)・壁面積 1e-5 一致 | 0106 削除済み(2026-08-31) (回帰記録)、0107 active (**y+1 真値の基準**) |
| **`run_0108_ib_m5_R3_LU9_Lc14_ns_ib_pass0`** | **排除厚さ更新計画 V4 (移植性)**: M5 推奨点に**積分法 (CONTUR) 初期壁 pass 0** を case/45 と同一設定で適用 (`feedback/deltastar_loop.py --init-integral`, 固定 Euler 参照 `run_0096`, 半径方向オフセット, IC=run_0102 cross-mesh, 24000 step; [plan](../../plans/accepted/tooling-nozzle-deltastar-core-matched-euler.md)) | 完走・NaN 0・品質 PASS・ALL STEADY・残差 warm 床。**ṁ_NS/ṁ_E = 0.9990 (v3 run_0102 は 1.0164)、出口面コア M 5.0008〜5.0027 (+0.02〜+0.05 %; v3 は −0.32〜−0.52 %)、軸 M は試験部で Euler +0.04〜0.13 % で平坦 (v3 は −1 % の波)**。帯局所抽出 δ_r との差 ≤2.5 % (固定点)。スロート補正 0.0012 r_t (v3 は 0.0098 = 8 倍過大) | active (**M5 の新チェーン正本**) |
| — | **軸則 D 案 (内部補間点 1 点 C⁴) の設計比較** (CFD なし): `compare_axislaw_D.py` → `compare_axislaw_D.json` / `compare_axislaw_D_summary.csv`。L_c 60→30 continuation、n_axis 600/1200/2400。D 不採用 ([plan](../../plans/accepted/tooling-nozzle-axislaw-onepoint.md), [結果ページ](https://claude.ai/code/artifact/b3142a8a-0c48-47d5-b7bc-f4653cc939c0)) | | ref |
| — | **A knot の $M_K$ 詳細感度** (CFDなし): `sweep_axislaw_A_MK.py` → `.json` / `_summary.csv` / `.png`。R3、L_c 30–60、M_K 1.5–4.0/0.1、n_axis 600 + 代表1200/2400。45は2.6、50は2.7で壁角・marginが小改善するが topology margin は6/13%低下。既定2.5を維持 ([knot plan](../../plans/accepted/tooling-nozzle-axismach-knot-law.md)) | | ref |
| — | **A knot の最短ロバスト軸長探索** (CFDなし): `optimize_axislaw_A_shortest.py` → `.json` / `_summary.csv` / [`report_axislaw_A_shortest.html`](report_axislaw_A_shortest.html)。hard gate・単峰・μ−θ≥1°・最小sin角≥0.02を制約に$x_F$最小化。600点で35–40を粗→0.01刻み、境界のみ1200点。**最短合格 $L_c=39.05,M_K=2.76$** ($x_F=94.34865r_t$, θ_max 16.1281°、margin 1.00095°、flip/交差0)。図: [`axis Mach`](best_axislaw_A_shortest_axis_mach.png) / [`Mach contour`](best_axislaw_A_shortest_mach_contour.png)、図示場 `best_axislaw_A_shortest_moc_field.npz`。2400点・CFDなし ([plan](../../plans/accepted/tooling-nozzle-shortest-robust-axis.md)) | | ref |
| — | 設計のみスイープ 36 点 → `study_design.json` (R=5 全滅、L_c=6 は R=1.5 のみ、L_U=9 は R=1.5 で μ>20) | | ref |

関連 case/41: `run_0076_cpg_axisymmethod1` = CPG でも `axisymMethod:1` は出口軸コーナーで発散
(TP 発散の切り分け副産物、別セッションへの申し送り)。
