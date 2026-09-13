# 09. Taylor-Green Vortex (3D)

3 次元 Taylor-Green 渦 (TGV)。全方向周期境界の立方体 $[0,2\pi]^3$、構造化 hex $32^3$ (節点 $33^3$)。
低散逸スキームの **運動エネルギー (KE)・エントロピー総和の保存性** を見るための検証ケース。

- メッシュ: [`mesh/Taylor-Green.geo`](mesh/Taylor-Green.geo) (`nx=ny=33`, z 方向 32 層 extrude)。
- 初期条件: `initial: "Taylor-Green"` (無次元, $M_0=0.4$, $\rho_0=1$, $P_0=1/\gamma\approx0.714$)。
  $u=M_0\sin x\cos y\cos z$, $v=-M_0\cos x\sin y\cos z$, $w=0$。
- 数値: pure KEEP (`solver: KEEP`, `keepDissipation: 0` = Roe 散逸無し・非散逸中心流束のみ), 陽解法 RK4 (`timeIntegration: 4`), 乱流モデル無し。

## 重要な設定上の注意

- **無次元・低圧ケースなので EOS フロアを下げること**: 既定の圧力フロア `pMin=1.0` Pa は本ケースの圧力場
  ($P\in[0.65,0.78]$) を全域クランプし、初期場を破壊する (step0 で $P\equiv1.0$)。`physProp.pMin/roMin/tMin`
  を小さく (例 `1e-6`) 指定する。floor の config 化は [`plans/accepted/thermophysics-eos-positivity-floor-config.md`](../../plans/accepted/thermophysics-eos-positivity-floor-config.md) 参照。
- **node モードのメッシュは node 設定で変換すること**: `discretization: "node"` のとき `convertGmshToForge` は
  median-dual を構築して primal を置換する。cell 設定で変換した h5 を node run に流用すると primal hex 上で
  解いてしまい不正 (双対 nCV=35937 ≠ primal hex 32768)。各 node run dir で自身の node config を使って変換する。

## 計算 run 一覧

全 run 共通: pure KEEP は `keepDissipation: 0`、陽解法 RK4、乱流なし、`pMin: 1e-6`、出力間隔 50 step。

全 run で全運動量 $\sum\rho\mathbf{u}\,V$ (理論=0) と全質量を保存量チェックに用いる。**node の res は周期マージ双対
セルを重複格納する** (VALUE 長 35937 = ノード数、一意 CV は 32³=32768) ため、後処理は各境界 CV を多重度
$m=2^{(境界方向数)}$ で割って一意 1 個分に補正する (`plot_ke_entropy_history.py` の `cv_weight`)。cell は補正不要。

以下は **cell 周期 partnerCellID device 転送バグ修正後** ([`plans/accepted/boundary-cell-periodic-conservation.md`](../../plans/accepted/boundary-cell-periodic-conservation.md)) の結果。
保存量チェックは全運動量 $\sum\rho\mathbf{u}\,V$ (理論=0) と全質量。**node の res は周期マージ双対セルを重複格納する**
(VALUE 長 35937 = ノード数、一意 CV は 32³=32768) ため、後処理は各境界 CV を多重度 $m=2^{(境界方向数)}$ で割る
(`plot_ke_entropy_history.py` の `cv_weight`)。cell は補正不要。

| run_* | 離散化 | 設定差分 | 主要結果 (多重度補正後) | 状態 |
| --- | --- | --- | --- | --- |
| `run_0009_node_keep_pure_eul` | node (median-dual) | pure KEEP, **非粘性** | **質量厳密・運動量 ~1e-7・KE 0.3%・エントロピー ~1e-5 保存**。教科書通りの KEEP 保存性 | active ✅ |
| `run_0010_cell_pure_eul` | cell (primal hex 32768) | pure KEEP, **非粘性** | **修正後: 質量厳密・運動量 ~1e-7・KE 0.4%・エントロピー ~1e-5 保存**。node と一致 (修正前は step358 発散) | active ✅ |
| `run_0007_node_keep_pure_visc` | node (median-dual) | pure KEEP, 粘性 (Re≈160) | 質量厳密・運動量 ~1e-7 保存。**KE 物理減衰 K/K0→0.644、エントロピー +1.3e-2 (第二法則)** | active ✅ |
| `run_0008_cell_pure_visc` | cell (primal hex) | pure KEEP, 粘性 | **修正後: 運動量 ~1e-7 保存・KE 物理減衰 K/K0→0.661** (node とほぼ一致; 残差は primal/dual メッシュ差)。修正前は KE×8 スプリアス増殖 | active ✅ |
| `run_keep`, `run_keep_M0.1` | cell | 旧参照入力 (旧スキーマ config) | — | ref |

### KEEP 陰解法 (block-DPLUR dual-time) 物理CFL掃引 (2026-06-29)

KEEP を陰解法化するときの flux Jacobian 方針 ([`plans/active/convection-keep-revive-node.md`](../../plans/active/convection-keep-revive-node.md) §7)
の検証。**新規コードゼロ** (`solver: KEEP` + `dualTime: 1` + `blockDPLUR: 1`, `timeIntegration: 11`) で、LHS は
既存の Roe 分割ヤコビアン (`accumulate_split_jacobian_cf`) を流用。全 cell・非粘性・共通 T≈3.5。explicit は
物理CFL≈0.05 (dt=0.007)、implicit は dt を上げて物理CFLを掃引。掃引解析は
[`analyze_implicit_sweep.py`](analyze_implicit_sweep.py) → [`implicit_sweep_ke_entropy.png`](implicit_sweep_ke_entropy.png)。

| run_* | 物理CFL | 主要設定差分 | K/K0 (T末, 理想=1) | dS/\|S0\| | 状態 |
| --- | --- | --- | --- | --- | --- |
| `run_0011_cell_keep_expl_ref` | 0.05 | explicit RK4 (基準) | 1.0033 | -1.1e-6 | active ✅ |
| `run_0012_cell_keep_impl_cfl005` | 0.05 | implicit, dt=explicitと一致 | **1.0033 (explicitと一致)** | -5.7e-7 | active ✅ |
| `run_0013_cell_keep_impl_cfl02` | 0.2 | implicit dt=0.028 | 1.0032 | -8.0e-7 | active ✅ |
| `run_0014_cell_keep_impl_cfl05` | 0.5 | implicit dt=0.07 | 1.0032 | 1.8e-6 | active ✅ |
| `run_0015_cell_keep_impl_cfl1` | 1.0 | implicit dt=0.14 | 1.0031 | 1.1e-5 | active ✅ |
| `run_0016_cell_keep_impl_cfl2` | 2.0 | implicit dt=0.28 | 1.0020 | 4.7e-5 | active ✅ |
| `run_0017_cell_keep_impl_cfl4` | 4.0 | implicit dt=0.56 | 0.9983 | 1.8e-4 | active ✅ |
| `run_0018_cell_keep_impl_cfl8` | 8.0 | implicit dt=1.12 | 0.9865 | 6.1e-4 | active ✅ |
| `run_0019_cell_keep_impl_cfl16` | 16 | implicit dt=2.24 | 0.9542 | 1.8e-3 | active ✅ |
| `run_0020_cell_keep_impl_cfl2_sub50` | 2.0 | CFL2 で nSubIterDualTime 20→50 | 1.0020 (sub20と完全一致) | 4.7e-5 | active ✅ |
| `run_0021_cell_keep_impl_cfl16_sub80` | 16 | CFL16 で nSubIterDualTime 20→80 | 0.9536 (sub20とほぼ同一) | 1.8e-3 | active ✅ |

### keepDissType (ES 散逸レイヤ) L2 較正 run ([plan](../../plans/accepted/convection-keep-es-dissipation.md))

| `run_*` | 目的・設定差分 | KE drop (500step, M0.4, RK4) | 状態 |
| --- | --- | --- | --- |
| `run_0022_cell_keep_dissoff_l2` | keepDissType=0 (新バイナリで無散逸経路の不変確認) | **−0.33%** (run_0011 参照と同挙動, 差は atomicAdd ノイズ級) | ref (不変確認) |
| `run_0023_cell_keep_diss015_l2` | keepDissType=1, σ=0.15, lowMachPrecond=1 | **8.37%** (市松6桁減衰の対価) | ref (σ 較正上限) |
| `run_0024_cell_keep_diss005_l2` | 同 σ=0.05 | **2.71%** (市松~4桁減衰で十分 → **既定採用**) | ref (σ 較正・既定根拠) |
| `run_0025_cell_keep_dissmat005_l2` | **keepDissType=2 (matrix ES)**, σ=0.05 | **1.36%** (scalar の半分。市松減衰は同等以上=case/35 run_0018) — せん断/エントロピー波に \|Un\| のみの選択性が効く | ref (Step 2 matrix L2) |
| `run_0026_cell_keep_dissmat005fullc_l2` | matrix σ=0.05 **フル c** (c' 無し) | **4.35%** (同一 σ 比較では c' 比 3.2 倍) — ただし下記 σ 掃引で「σ 弱化すれば同等 Pareto 点に到達」と判明 (同一 σ 比較は設計比較としては誤導) | ref (c' 切り分け・同一σ) |
| `run_0027_cell_keep_dissmat_fullc_s0015_l2` | フル c **σ=0.015** (σ 弱化で c' 代替) | **1.10%** (市松 3.9e-8 = case/35 run_0021)。**単一マッハ領域では c' 不要・σ 弱化で十分**。c' の残存価値はマッハ混在流の面ごと自動調整のみ | ref (σ掃引・c'代替検証) |
| `run_0028_cell_keep_dissmat_fullc_s004_l2` | フル c σ=0.04 | 3.43% (σ にほぼ線形) | ref (σ掃引) |

**所見**: (1) KEEP+block-DPLUR は **CFL≈16 (explicit の320倍) まで NaN なしで完走** — Roe-Jacobian LHS の流用で
陰解法が安定動作。(2) 物理CFL≤2 では explicit と KE/エントロピー保存が一致 (matched CFL.05 は K/K0=1.0033 で
explicit と4桁一致) → **LHS の Roe 散逸は収束解を汚染しない**。(3) CFL を上げると KE が単調減衰 (CFL16 で
-4.6%)・エントロピー増。**ただしこの劣化は内部反復不足 (LHS散逸漏れ) ではなく BDF2 dual-time の時間離散誤差**が
支配: CFL2 で内部残差は 6.8e-5→3e-9 と完全収束し sub20≡sub50、CFL16 でも sub20≈sub80 で KE が回復しない
(サブ反復を増やしても改善しない=温度離散誤差は irreducible)。→ **LHS の Jacobian 選択は解の質にほぼ無関係**で、
大 dt での精度限界は外側の時間積分 (BDF2) が決める。LES では物理 dt を時間精度が許す範囲に保つこと。

成果物: KE・エントロピー時間履歴 [`ke_entropy_history.png`](ke_entropy_history.png)、ポストスクリプト
[`plot_ke_entropy_history.py`](plot_ke_entropy_history.py) (旧 `plot_taylorGreen.py` の現行版=res ベース・多重度補正付き)。
各 run の VERDICT は `CONVERGENCE_VERDICT.txt` (非定常 TGV のため定常収束ツールは `NOT CONVERGED`=正常)。

### 結論

- **pure KEEP は cell・node とも KE/エントロピー/運動量を保存する** (修正後)。非粘性で KE 0.3–0.4%・エントロピー
  ~1e-5・運動量 ~1e-7、粘性で KE 物理減衰・エントロピー物理増加。cell と node が一致 (差は primal hex vs median-dual)。
- **根本原因 (修正済)**: cell 全周期で運動量が線形注入され seam が非周期化していたのは、**`bint_d["partnerCellID"]`
  (device) が host から未転送**だったため (`setPeriodicPartner` は host のみ充填、device コピーは yaml uniform 経路だけ)。
  periodic_d が未初期化 device partnerCellID を読み ghost が誤値に。`setPeriodicPartner` 直後に H2D コピーで根治。
  node は host の partnerCellID を直接使う (`buildPeriodicNodeGroups`) ため無傷だった。float/double・ghost 更新
  タイミングは無関係 (切り分け済)。
- 切り分けの決め手: ghost セル値 vs 解析 TGV(ghost 重心) が修正前 6112/6144 不一致 (誤差 0.43) → 修正後 0 (1.5e-7)。
- 補足: 過去 cell で「KEEP」と思っていた流束は legacy で `if(false)` 無効化されており実体は **MUSCL+Roe**。本物の
  KEEP 中心流束は `Revive KEEP` (79b4e67) で初めて有効化された (この periodic バグが顕在化したのもそのため)。
- 旧 `run_keep` が現バイナリで回らないのは solver 退行ではなく①config スキーマ進化②圧力フロア (`pMin` 既定 1.0 Pa
  が無次元低圧場をクランプ) であり、modern config + `pMin` 引き下げで解消する。

### node L2 / L3 (WALE 込み粘性 TGV Re=1600) — keepDiss 検証続き

| `run_*` | 目的・設定差分 | 主要結果 | 状態 |
| --- | --- | --- | --- |
| `run_0029_node_keep_dissmat005_l2` | **node L2**: node TGV M0.4 + matrix σ=0.05 (run_0009 σ=0 と比較) | KE cost **+1.371%**/500step = cell (+1.356%) と一致。**散逸レイヤは node 内部で cell 同等** | ref (node L2 PASS) |
| `run_0030_cell_keep_wale_re1600_s0` | **L3 基準**: 粘性 Re=1600 (visc 2.5e-4, Pr0.71) + WALE, keepDissType=0, 3600step (t*≈10) | ε* ピーク 0.066 @ t*=10.1、K/K0(10)=0.71。層流期はほぼ無散逸 (正)。32³ ゆえピーク過小・遅め | ref (L3 基準) |
| `run_0031_cell_keep_wale_re1600_diss` | L3 + matrix **σ=0.05** | ピーク 0.083 @ t*=8.1、K/K0(10)=0.49。**層流期 (t*<4) に 5-6% 食う = 滑らか場への連続ドレインが強すぎ**。LES 用には過大 | ref (L3 σ過大の記録) |
| `run_0032_cell_keep_wale_re1600_diss0015` | L3 + matrix **σ=0.015** | ピーク 0.080 @ t*=8.4、K/K0(10)=0.60、層流期損失 2.6% | ref (L3 σ掛引) |
| `run_0033_cell_keep_wale_re1600_diss002` | L3 + matrix **σ=0.02** | ピーク 0.082 @ t*=8.4、K/K0(10)=**0.594**、層流期 3.1% | ref (L3 σ掛引) |
| `run_0034_cell_keep_wale_re1600_diss003` | L3 + matrix **σ=0.03** | ピーク 0.083 @ t*=8.4、K/K0(10)=**0.557 (DNS帯 0.50-0.57 内)**、層流期 4.1% | ref (L3 σ掛引) |

L3 の判定は定性ゲート (ピークが t*≈9±1・カーブ滑らか・NaN なし): 全 run 通過。32³ では層流期と終値を同時に満たす σ は無い (解像度律速)。~~DNS 参考値 (K/K0(10)≈0.50-0.57...) はアンカー2点+積分近似で ±0.05 幅~~ **← この近似帯は誤り** (下の 64³ 節の実データ参照: 真値 K/K0(10)=0.596。32³ σ=0.03 を「帯内」とした判定は近似帯の誤りによるもので、実データでは σ=0.02 の 0.594 がほぼ厳密一致だった)。

### 64³ L3 定量化 (2026-07-19, 実 DNS データで σ 確定)

DNS 参照を実データ化: [`ref_dns/TGV_Re1600.dat`](ref_dns/README.md) (**Dairay+2017 JCP, Incompact3d 512³**, t=0..20 の Ek/ε_t 表; 真のアンカー **K/K0(10)=0.596, ε* ピーク=0.1029 @ t*=8.98**)。`plot_dissipation_rate.py --dns` で重ね描き可。

| run_* | σ | K/K0(4) 誤差 | K/K0(10) 誤差 | ε* ピーク (DNS 0.103) | ピーク t* (DNS 8.98) | 状態 |
| --- | --- | --- | --- | --- | --- | --- |
| `run_0035_cell_keep_wale64_re1600_s0` | 0 | **+0.8%** | +1.7% | 0.085 | **8.96** | ref (64³ L3 基準) |
| `run_0036_cell_keep_wale64_re1600_diss002` | 0.02 | −0.6% | **−1.5%** | 0.081 | 8.68 | ref (64³ σ掛引・**推奨確定**) |
| `run_0037_cell_keep_wale64_re1600_diss003` | 0.03 | −1.3% | −2.8% | 0.079 | 8.68 | ref (64³ σ掛引) |
| `run_0038_node_keep_wale64_re1600_diss002` | 0.02 (**node**) | −0.7% | −2.1% | 0.082 | **8.96** | ref (64³ node L3) |
| `run_0039_node_keep_wale64_re1600_s0` | 0 (**node**) | +0.7% | −0.4% | 0.092 | 9.24 | ref (64³ node σ=0 基準) |
| `run_0040_node_keep_wale64_re1600_diss002_jump` | 0.02 (**node, keepDissJump=1**) | **+0.7%** | **−0.8%** | **0.088** | **8.96** | ref (**recon-jump L3**) |
| `run_0041_node_keep_wale64_re1600_diss002_jump2` | 0.02 (**node, keepDissJump=2**) | +0.7% | −1.4% | 0.089 | **8.96** | ref (**sign-property L3**) |
| `run_0042_node_keep_wale64fix_s0` | 0 (**修正後の本物 WALE**) | −0.3% | −3.4% | 0.087 | **7.84 (早すぎ)** | ref (real-WALE 基準) |
| `run_0043_node_keep_wale64fix_diss002_jump2` | 0.02 jump2 (**本物 WALE 併用**) | −0.4% | −3.6% | 0.085 | 7.84 | ref (real-WALE+ES) |
| `run_0044_node_keep_sigma64_s0` | 0 (**σ-model** LESmodel:2) | −0.5% | −5.5% | 0.082 | 7.84 | ref (σ-model 基準) |
| `run_0045_node_keep_sigma64_diss002_jump2` | 0.02 jump2 (**σ-model 併用**) | −0.5% | −5.6% | 0.082 | 7.84 | ref (σ-model+ES) |
| `run_0046_sst_hdt_uniform_ek1` / `run_0047_sst_hdt_uniform_ek0` / `run_0048_sst_hdt_sine_ek1` | **一様乱流減衰 (SST の E_t 保存試験)**: run_0007 の node 周期メッシュ流用、SLAU + SST、u=0・ρ=1.4・P=1・ω=1、k=0.1 一様 (uniform) / k=0.1(1+0.5 sin x) (sine)、RK4 explicit dt 0.005 × 400。ek1=`sstEnergyIncludesK: 1`, ek0=0 (plan [turbulence-sst-energy-includes-k](../../plans/active/turbulence-sst-energy-includes-k.md) §6.1) | **ek1: ΣV(E_m+ρk) 変化 ≤1.3e-7、max\|u\| 1.9e-7 (一様 p* で偽の力なし)、k 0.100→0.0847 の減少が c_v ΔT に 5 桁一致 (Δe=+0.01534=−Δk)。ek0: E_t −0.81 % (k の散逸が消える)、T 不変。sine: p* 勾配で u≈0.02 が立つが E_t 変化 ≤1.1e-7・質量 4e-8** | active (保存試験) |
| `run_0049_sst_hdt_uniform_ek1_dual` / `run_0050_sst_hdt_sine_ek1_dual` | 同上を **dual-time (timeIntegration 11 + dualTime 1, block-DPLUR)** で。エネルギー行の BDF に ρk を含める修正後 | E_t 変化 uniform +3.5e-7 / sine −3.1e-6 (subiter 収束で決まる)、c_p ΔT = +0.02148 (explicit と同一)。修正前は E_t −0.81 % で T 不変だった | active (保存試験) |
| `run_0051_monitor_dualtime` | run_0050 の複製 (100 step)。**console モニタ行** ([plans/accepted/architecture-runtime-monitor-line.md](../../plans/accepted/architecture-runtime-monitor-line.md)) の dual-time 陰解法 (node) 検証: `t`/`dt`/`maxCFL` (物理 CFL) + ms/step + 残差要約。`_baseline/` は旧バイナリ比較用 | 完走・NaN 0。`t/dt/maxCFL 0.13` 表示、log 4330→230 行、CSV 旧×2/新×2 列別比較 (`csv_compare.txt`) 合格 (1.2e-4 vs 床 1.3e-4) | 破棄予定 (ログ形式検証) |

**結論 (σ 確定)**: 64³ では層流期・終値とも全 σ≤0.03 で ±3% 内 = 32³ の「両立不可」は解像度律速と確定。σ=0 (+1.7%) と σ=0.02 (−1.5%) が DNS をほぼ対称に挟み、**解像 LES の推奨は σ=0.02** (市松頑健性 L1 を持ちつつ KE 追従 ≤1.6%)。σ=0.03 は 64³ では一律に劣る (旧「0.03=帯内」判定は撤回)。TGV 追従だけなら最適 σ≈0.01 だが、σ は物理較正ノブではなく市松ロバスト性の床 — 解像度を上げるほど最適値は下がる (32³: ~0.02 → 64³: ~0.01) ので、「L1 を満たす最小 σ」で選ぶ。ε* ピーク高さは 64³+WALE の解像限界で全 run 8 割どまり (LES として正常)。成果物 `dissipation_rate_L3_64.png`。

**node L3 (run_0038)**: 同条件 node (median-dual 65³=274625 CV, 周期 DOF 合併 12481=厳密) は
cell とカーブがほぼ重なり (K/K0(10) で node −2.1% vs cell −1.5%、ピーク時刻は node 8.96 が DNS
8.98 により近い)、**node×σ=0.02 も解像 LES として cell 同等に検証済**。KE 集計は周期 slave CV の
ミラー重複を除くため `plot_dissipation_rate.py --node-mesh` (半開区間 dedupe) を使うこと
(素朴な総和は周期面を 2-8 倍過重み)。成果物 `dissipation_rate_L3_64_node.png`。

**再構成ジャンプ散逸 (keepDissJump=1, run_0039/0040)**: matrix 散逸の Δw を再構成後ジャンプで
組むと ([recon-jump plan](../../plans/accepted/convection-keep-diss-recon-jump.md))、σ=0.02 の
KE コストが**終値 −2.1%→−0.8%・層流期 −0.7%→+0.7% (=σ=0 と同一・実質ゼロ)** に縮み、遷移期の
ε\* の早すぎる立ち上がりも解消、ピーク時刻 8.96 (σ=0 単独は 9.24 と遅い) は維持 — つまり
**市松頑健性と σ=0 並みの KE 追従を両立**。市松減衰は無傷 (case/35 run_0039/0040)。
**node LES の新推奨: `keepDissType: 2, keepDissCoeff: 0.02, keepDissJump: 2`**。
成果物 `dissipation_rate_L3_64_node_jump.png`。

**sign-property クリップ (keepDissJump=2, run_0041)**: 再構成ジャンプの特性射影を生ジャンプ射影と
minmod し、各波で積 ≥0 を構造保証 = **エントロピー散逸性が証明付きで復活** (TeCNO 型)。KE コストは
jump=1 (−0.8%) と生 (−2.7%) の間の **−1.4%**、ピーク時刻 8.96・層流期 +0.7% は jump=1 と同一。
市松減衰も無傷 (case/35 run_0042)。**証明と性能の両立点として jump=2 を推奨既定とする**
(jump=1 は証明なしの最軽量オプションとして残置)。成果物 `dissipation_rate_L3_64_node_jump2.png`。

**重大訂正 (2026-07-19, WALE 不活性バグ発覚)**: 本節までの全「WALE 併用」run (32³ run_0030-0034・
64³ run_0035-0041, cell/node とも) は、**WALE が不活性 (vis_turb≡0) の ILES だった**。真因は
① 壁なしメッシュで wall_dist≡0 → Ls=min(κd,CwΔ)=0、② WALE の Sd テンソルが成分 2 乗の誤式
(行列 2 乗が正; 純せん断で Sd=0 になるべき性質を破っていた)。両方修正済
([turbulence-wale-fix](../../plans/accepted/turbulence-wale-fix.md))。**数値結果は「ILES
(KEEP+分子粘性+ES散逸)」として全て有効** — ラベルのみ訂正。修正後の本物 WALE で再検証した結果
(run_0042/0043): **WALE は遷移期に先回りして散逸し (ピーク t*=7.84 ≪ DNS 8.98)、終値も −3.4% と
ILES より悪化**。64³ TGV 級の解像/遷移流では **WALE off の ILES + matrix ES (σ=0.02, jump=2) が
最良 (−1.4%, ピーク 8.96)** = **解像 LES の推奨は LESorRANS: 0 + ES 散逸**。WALE は「本当に
未解像の高 Re 乱流」用オプションとし、その適用検証は今後の課題。成果物
`dissipation_rate_L3_64_walefix.png`。

**σ-model (run_0044/0045)**: 静的モデル改良版 (LESmodel:2) も検証したが、層流期 ν_t は設計どおり
WALE 比 −25% になるものの圧縮性 TG の早期 3 次元化でゼロ性質の効く時間が短く、発達後は C_σ=1.35 の
散逸が勝って **−5.5% (WALE より悪い)・ピーク 7.84**。**64³ TGV 級では静的 SGS はどちらも逆効果、
ILES+ES (σ=0.02, jump=2) が最良のまま**。σ-model の本領 (壁乱流・回転流・未解像高 Re) の検証は今後
([plan](../../plans/accepted/turbulence-sigma-model.md))。成果物 `dissipation_rate_L3_64_sigma.png`。
| `run_0052_perf_regress_node_periodic_dualtime` (ローカル, 2026-09-12) | 高速化ブランチ `feature/perf-3d-speedup` の回帰: run_0049 (node 周期, dual-time 陰解法, SST, sstEnergyIncludesK=1) の res_400 から 100 step 継続、基準バイナリ ×2 vs 最終 (`tools/perf_regress.py`) | P/T/k/ω/ρK/ρΩ/μt の差 ≤4.3e-7 で base×base と同等。速度 Ux/Uy は base×base 自体が 0.32 (減衰乱流の run 間ばらつき) で比較不能 | active (回帰) |
