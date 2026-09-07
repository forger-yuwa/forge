# 境界条件

forge は密度ベース有限体積で **ゴーストセル方式** を採用し、各境界条件は
内部セル状態と境界条件パラメータからゴーストセル値を構成する。
内部面ループでは境界面もそのまま処理されるため、ゴースト値が境界条件を
そのまま表現する形になる。

本ドキュメントは理論(係数・方程式)と実装(ソース対応)をまとめる。

## 理論

### 共通モデル

非周期境界面 $f$ について内部セルを $L$、ゴーストセルを $R$ と書く。
境界条件は次の二系統を返す。

1. **ゴーストセル状態** $(\rho, \rho \mathbf u, \rho e, U, P, T, H_t, c)_R$ —
   対流フラックス計算で R 状態として使用される。
2. **境界面値 `bvar_d[*]`** — 粘性束計算・post 処理 (壁面摩擦・$y^+$) で参照される。

非周期境界では対流再構成が `scheme = -1` (1 次風上) に強制されるため、
ゴーストセル値はそのまま面値として使われる ([`methods/convection/`](convection) 参照)。

### 提供される境界種別

| `kind` | 用途 | ゴースト構成の概要 |
| --- | --- | --- |
| `slip` | 滑り壁 (オイラー壁) | 法線速度を反転 ($\mathbf u_R = \mathbf u_L - 2 U_n \hat{\mathbf n}$)。圧力・密度は同値 |
| `wall` | 非滑り断熱壁 | $\mathbf u_R = -\mathbf u_L$、$P_R = P_L$、$T_R = T_L$ (断熱) |
| `wall_isothermal` | 非滑り等温壁 (移動壁可) | $\mathbf u_R = 2\mathbf u_{\text{wall}} -\mathbf u_L$ (壁速度は floats `Ux/Uy/Uz`、既定 0 で従来の $-\mathbf u_L$ とビット同一。粘性仕事 $\tau\cdot\mathbf u_{\text{wall}}$ は壁カーネルが bvar 面速度で加算。Couette 厳密解で検証: `case/24` `run_couette_cell`, q_w ±0.2%)、$T_R = 2 T_{\text{wall}} - T_L$ (負温度ガード $T_R \ge 0.2\,T_{\text{wall}}$)。旧実装は $T_R = T_{\text{wall}}$ 直置きで、壁熱流束 $(T_R-T_L)/d_{cc}$ ($d_{cc}=2y_1$) が正しい値の 1/2 になるバグがあった (2026-07-20 修正。純伝導厳密解 `case/24` `run_isoT_cond*` で検証: 修正後は線形場が機械精度の不動点・$q_w$ 誤差 +0.02%) |
| `inlet_uniformVelocity` | 均一速度流入 | $\mathbf u_R, \rho_R$ を指定値に固定、$P_R = P_L$ |
| `inlet_fluctVelocity` | 速度変動つき流入 | uniformVelocity に変動成分を加算 (`fluct_variables`) |
| `outlet_statPress` | 静圧固定流出 | $P_R = P_{\text{back}}$ を課し、$\rho$・速度は内部エントロピー＋外向き Riemann 不変量で構成 (亜音速)。逆流時も同じ静圧アンカー |
| `inlet_Pressure` | 全圧・全温固定流入 | 全条件 ($P_t, T_t$) から内部マッハで $P, T$ を再構成 |
| `inlet_Pressure_dir` | 方向指定全圧流入 | inlet_Pressure に流入方向ベクトルを併用 |
| `outflow` | サブソニック流出 | リーマン不変量に基づく Non-reflecting 流出 |
| `periodic` | 周期境界 | 対応するペア面のセル値をコピー (`scheme` 強制なし) |

### 例: 滑り壁

法線方向速度 $U_n = \mathbf u_L \cdot \hat{\mathbf n}$ を用いて

$$
\mathbf u_R = \mathbf u_L - 2 U_n \hat{\mathbf n},
\quad P_R = P_L, \quad \rho_R = \rho_L,
\quad \rho e_R = \frac{P_L}{\gamma - 1} + \tfrac{1}{2}\rho_L |\mathbf u_L|^2.
$$

これは法線フラックスから運動量の壁面垂直成分が消え、圧力束のみ残る効果を持つ。

### 例: 非滑り壁

$\mathbf u_R = -\mathbf u_L$ とすることで、面値 $\mathbf u_f = \tfrac{1}{2}(\mathbf u_L + \mathbf u_R) = 0$
(無滑り)。粘性束は壁面用カーネル ([`methods/diffusion/`](diffusion)) で別途加算する。

### 例: 全圧固定流入

総温・総圧 $T_t, P_t$ と局所マッハ $M_L$ を用い、

$$
T_R = \frac{T_t}{1 + \tfrac{\gamma - 1}{2} M_L^2},\quad
P_R = P_t \left(\frac{T_R}{T_t}\right)^{\gamma/(\gamma-1)},\quad
\rho_R = \frac{P_R}{R T_R},
$$

速度方向は外挿 (`inlet_Pressure`) または指定方向 (`inlet_Pressure_dir`)。

### 例: 静圧固定流出 (特性ベース・逆流統一)

亜音速流出では SU2 `CEulerSolver::BC_Outlet` と同様、指定静圧 $P_{\text{exit}}$ のみを境界条件として課し、
残りは内部状態から自己整合に構成する (1-incoming-characteristic)。内部エントロピー $s = P_L/\rho_L^\gamma$ と
外向き Riemann 不変量 $R^+ = U_n + 2c_L/(\gamma-1)$ を保存量として、

$$
\rho_R = \left(\frac{P_{\text{exit}}}{s}\right)^{1/\gamma},\quad
c_R = \sqrt{\gamma P_{\text{exit}}/\rho_R},\quad
V_n = R^+ - \frac{2 c_R}{\gamma-1},\quad
\mathbf u_R = \mathbf u_L + (V_n - U_n)\hat{\mathbf n}
$$

(接線速度は内部外挿、法線のみ $V_n$ へ補正)。超音速流出 ($M_L\ge 1$ かつ $U_n>0$) では境界条件不要で全量外挿。

**逆流 (局所流入, $U_n<0$) も同じ静圧アンカーで扱う**: 上式で $V_n<0$ となるのを許容し (クランプ無)、
incoming/outgoing 特性の捌きは upwind フラックスに委ねる。出口逆流で全圧 ($P_t,T_t$) の stagnation
流入へ切替えると、剥離域 (壁∩出口コーナー) へ高 stagnation エンタルピを注入し過加圧→発散させる
(これは `inlet_Pressure` の構成であり出口に流用すべきでない)。乱流スカラー $k,\omega$ は出口で
ゼロ勾配 (Neumann) であり、逆流時も内部値を再循環させる (固定値注入はしない)。

### 周期境界

周期境界は他境界と異なり、対応するペアセル値を直接コピーする。
対流再構成は内部面と同じ MUSCL 経路で処理されるため、`scheme` の強制 1 次降格は無い。

### 物理 ID と YAML 設定

各境界面はメッシュ生成時に物理 ID を持ち、`bcondConfig.yaml` で
ID → `kind` → 数値パラメータ (`ints`, `floats`) を紐付ける。
パラメータの種類 (定数 / 配列入力) は `bcondConfFormat::valueTypesOfBC` で定義する。

### 参考

- 適用順序: [`methods/time_integration/`](time_integration) の "ループ全体" を参照。
- 粘性壁面寄与: [`methods/diffusion.md`](diffusion.md#実装) の
  `viscousFlux_wall_d` を参照。
- 対流側の境界処理 (ゴースト → 1 次風上): [`methods/convection/`](convection) を参照。

## 実装

### ソースファイル

| ファイル | 役割 |
| --- | --- |
| [`solver_density_cuda/boundaryCond.hpp`](../solver_density_cuda/boundaryCond.hpp) | `bcond`, `bcondConfFormat` 構造体、各 BC kind の値型テーブル |
| [`solver_density_cuda/boundaryCond.cpp`](../solver_density_cuda/boundaryCond.cpp) | YAML 読み込みと `applyBconds` ディスパッチ |
| [`solver_density_cuda/cuda_forge/boundaryCond_d.cuh`](../solver_density_cuda/cuda_forge/boundaryCond_d.cuh) | GPU カーネル / ラッパ宣言 |
| [`solver_density_cuda/cuda_forge/boundaryCond_d.cu`](../solver_density_cuda/cuda_forge/boundaryCond_d.cu) | 全 BC kind のカーネル本体 (約 1500 行) |

CPU 経路 (`cfg.gpu == 0`) は未対応 (`applyBconds` 冒頭で `exit`)。

### YAML 設定

ケース直下の `bcondConfig.yaml` で各物理 ID の BC を定義する。

```yaml
1:
  physID: 1
  kind: wall
  outputHDFflg: 1
  ints: {}
  floats: { Ts: 300.0 }
```

`boundaryCond.cpp` の [`readBcondConfig`](../solver_density_cuda/boundaryCond.cpp#L11) が
YAML を読み、各 `bcond` (メッシュ生成側で物理 ID 付き) に `kind`, パラメータ,
値型テーブル (`valueTypes`) を紐付け、`bcondInitVariables` で
GPU/CPU 双方の境界変数バッファを確保する。

#### node 等温壁の壁ノード温度ピン (2026-07-20)

node (median-dual) の `wall_isothermal` は、壁ノードの温度状態を BC 壁温に**強制ピン**する
(運動量の `nodeWallDirichlet` の熱版。`nodeWallDirichlet: 1` [既定] のとき有効)。3 点セット:

1. **状態ピン** (`applyNodeIsothermalWallPin`): applyBconds 位相で壁ノードの T/roe/P/sonic を
   Tw に整合させる (WMLES 等温壁と同一カーネル `pin_wall_node_temperature_d` を共用)。
2. **残差ゼロ化** (`zeroNodeIsothermalEnergyResidual`): 壁ノードの `res_roe` を 0 に射影
   (Dirichlet ノードの残差は BC 強制であり物理不均衡でない)。
3. **陰解法 Jacobian 整合** (`iso_wall_flag_d`): block-DPLUR の 5×5 対角ブロックで壁ノードの
   エネルギー行 (row 4) を単位行化 + rhs=0 (SU2 `DeleteValsRowi` 相当、壁運動量 3 行 decouple と
   同型)。**これ無しで状態ピンだけ行うと implicit は数 step で発散する** (2026-07-20 実測)。

実装の所在: **node 壁強境界の状態層は [`cuda_forge/nodeWallDirichlet_d.cu`](../solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu) に集約**
(u=0 状態ピン `enforceWallNoSlip` / 残差射影 `zeroWallDirichletResiduals` [運動量+SST ω+等温 roe] /
温度ピン)。壁関数 (`ransWallFunction_d`)・壁モデル (`wmlesWallModel_d`) は「モデル層」として分離され、
流束層 (`viscousFlux_d` の Tau_Wall/Qw_Wall マーカ) 経由でのみ結合する — 状態層は
`wallTreatmentSST`/`wallModelLES` を参照しない。所与温度の熱力学状態は `thermo_d.cuh` の
`thermo_state_at_T` / `thermo_cell_Y` (共通 helper) を使う。

検証 (case/24 純伝導・厳密解シード): 壁ノード T = BC 値厳密 (350.0000/300.0000)、中央部熱流束
厳密 (+0.00%)。既知の残課題: 第 1 スペーシングの勾配に −15% 程度の局所バイアス (W-I 弱形式閉包の
離散化誤差・O(Δ))、implicit の実用擬似 CFL 上限 ~5 (エネルギー行 decouple により境界が実効陽扱いの
ため。cfl_pseudo 20 は発散)。ピン導入前は壁ノード T が壁 CV 平均に緩み (~0.1 K オフセット)、
第 1 スペーシング勾配 −24% だった。

#### WMLES 壁モデルの指定 (`wallModelLES`)

`wall` / `wall_isothermal` の `ints:` に `wallModelLES: 1` を書くと、その壁の粘性流束が
代数壁応力モデル (Reichardt + Kader, [`methods/turbulence/theory.md`](turbulence/theory.md) §10)
で置き換わる (既定 0 = 従来の解像壁)。LES/ILES (`turbulence.model` が `sst*` 以外) でのみ有効で、
SST automatic wall treatment (`wallTreatmentSST`) とはコードパスが分離されている。

```yaml
2:
  physID: 2
  kind: wall            # 断熱 WMLES 壁 (等温は kind: wall_isothermal + floats: {Ts: ...})
  outputHDFflg: 1
  ints: { wallModelLES: 1 }
  floats: {}
```

### 入口分布プロファイル (inlet profile)

入口 (`inlet_*`) の境界値を**一様でなく分布として**与える機能。inlet カーネルは per-face の
境界値 `bvar` (`Uxb[ib]`,`Uyb[ib]`,…) をそのまま読むため、**face 重心座標に応じて `bvar` を
非一様にセットすれば kernel 無改修で分布入口が実現**できる。これを
[`applyInletProfiles`](../solver_density_cuda/boundaryCond.cpp) が担う。

- **有効化**: 対象 inlet の YAML に `ints: {inletProfile: 1}` を付ける。
- **入力**: run dir の `inlet_profile_<physID>.csv`。**1 行目ヘッダで補間方向と量を指定**する。
  - 先頭の連続する `x`/`y`/`z` 列 = **補間座標**。1 列 (例 `y`) → **1D 線形補間** (その軸でテーブルを
    昇順ソートし線形補間、範囲外は端値クランプ)。3 列 (`x y z`) → **3D 最近傍**。
  - 残り列 = **bvar 量名** (`Ux Uy Uz Tt Pt Ts Ps ro k omega`、多成分 TP は `Y0..Y{n-1}` 等)。その inlet が
    持つ bvar のみ反映 (無い列名は黙って無視。`inlet_uniformVelocity` に `Tt` を書いても効かない)。
  - 例 (1D-y): ヘッダ `y Ux Uy Uz` に続けて行 `1.0 0 0 0` / `1.02 30 0 0` …。
- **種別ごとに分布化できる量** (2026-09-08 検証): `inlet_Pressure` (亜音速) は `Tt Pt Y{s} k omega`、
  `inlet_uniformVelocity` (超音速・全量固定) は `ro Ux Uy Uz Ps Y{s} k omega`。全温分布を超音速入口に与えるには
  静的状態 (ρ, U, Ps) に換算する。`k/omega` は `rans_dirichlet_scalar_boundary_d`、`Y{s}` は
  `species_dirichlet_boundary_d` がそれぞれ per-face `bvar` を読む。凝縮モーメントは入口で常に 0 (分布不可)。
- **node の化学種入口ピン** (2026-09-08 修正): node は境界半割面の化学種流束に ghost を使わず境界ノード自身の
  組成を使うため、ghost だけの Dirichlet では入口ノードの `Y{s}` が初期値のまま凍結していた (一様値の変更も
  分布も入らない)。k/ω と同じく `species_dirichlet_boundary_d` (node 分岐) が境界ノードを `Y_s^in` にピンし
  `scalarDirichletPin` を立て、`speciesPinResidual_d_wrapper` (化学ソース集計後) が `res_roY{s}`/`src_jac_Y{s}`
  を 0 化する。連成陰解法 (`speciesImplicitCoupling: 1/2`) の `species_dplur_sweep_d` はピン行を δ(ρY)=0 に
  拘束する (残差 0 だけでは隣接補正が Jacobi sweep で漏れる)。cell は不変 (ghost Dirichlet のまま)。
  `applyInletProfiles` は反映列 (`applied:`) と無視列 (`IGNORED`) をログに出し、多成分では各 face の
  0≤Y_s≤1・ΣY_s=1 を検査する。`inlet_Pressure_dir` は組成 `Yb` を受け取らない単成分熱物性なので多成分非対応。検証: case/16 run_0317 (node) / run_0318 (cell) /
  run_0319 (node 一様入口の回帰, 旧バイナリと 1e-6 以下)。
- **生成・照合ツール**: [`tools/gen_inlet_profile.py`](../solver_density_cuda/tools/gen_inlet_profile.py)
  (`gen`: 座標の式または測定表から CSV。TP は NASA-9 で Tt/M/Ps→ρ,U,Ps 換算。`verify`: 境界ノード/第 1 セルの
  T0 (h0 由来)・Y・k・ω を目標と比較)。運用手順は [`procedures/inlet-profile.md`](../procedures/inlet-profile.md)。
- **適用箇所**: `main.cpp` で `readBcondConfig` (bvar セット) の直後・最初の `applyBconds` より前に
  `applyInletProfiles(cfg, msh)` を呼ぶ。face 重心 `msh.planes[ip].centCoords` で補間し host `bvar` に
  セット → `bvar_d` に再アップロード。`inletProfile` 未指定の inlet は一様のまま (挙動不変)。
- **壁法則 helper**: [`tools/gen_inlet_walllaw.py`](../solver_density_cuda/tools/gen_inlet_walllaw.py) が
  Reichardt 合成則でチャネル/片側 BL の $u(y)$ を生成し CSV を書く ($u_\tau$ は中央/外縁速度 = `Uc` に
  なるよう二分法)。設計判断は
  [`boundary-inlet-profile.md`](../plans/accepted/boundary-inlet-profile.md)。

### ディスパッチ

[`applyBconds`](../solver_density_cuda/boundaryCond.cpp#L116) が
`msh.bconds` をループし、各 `bc.bcondKind` で次のラッパを呼ぶ。

| `bcondKind` | ラッパ | カーネル位置 |
| --- | --- | --- |
| `slip` | `slip_d_wrapper` | [L86](../solver_density_cuda/cuda_forge/boundaryCond_d.cu#L86) |
| `wall` | `wall_d_wrapper` | [L214](../solver_density_cuda/cuda_forge/boundaryCond_d.cu#L214) |
| `wall_isothermal` | `wall_isothermal_d_wrapper` | [L363](../solver_density_cuda/cuda_forge/boundaryCond_d.cu#L363) |
| `outlet_statPress` | `outlet_statPress_d_wrapper` | [L531](../solver_density_cuda/cuda_forge/boundaryCond_d.cu#L531) |
| `inlet_uniformVelocity` | `inlet_uniformVelocity_d_wrapper` | [L651](../solver_density_cuda/cuda_forge/boundaryCond_d.cu#L651) |
| `inlet_fluctVelocity` | (`fluct_variables_d` 経由) | — |
| `inlet_Pressure` | `inlet_Pressure_d_wrapper` | [L863](../solver_density_cuda/cuda_forge/boundaryCond_d.cu#L863) |
| `inlet_Pressure_dir` | `inlet_Pressure_dir_d_wrapper` | [L1033](../solver_density_cuda/cuda_forge/boundaryCond_d.cu#L1033) |
| `outflow` | `outflow_d_wrapper` | [L1184](../solver_density_cuda/cuda_forge/boundaryCond_d.cu#L1184) |
| `periodic` | `periodic_d_wrapper` | [L1310](../solver_density_cuda/cuda_forge/boundaryCond_d.cu#L1310) |

呼び出し後に `cudaPeekAtLastError + cudaDeviceSynchronize` で同期。

#### `outlet_statPress` の特性ベース構成と逆流統一

`outlet_statPress_d` は亜音速流出で指定静圧 `Psb[ib]` のみを境界値とし、$\rho$・速度は内部セル `ic` の
エントロピーと外向き Riemann 不変量から構成する (理論は 本ドキュメントの「理論」節 参照)。**逆流 (`Un<0`) も
同じ静圧アンカーで扱う** (`Vn_exit<0` 許容)。値の取得元に注意: 静圧のみ境界 `Psb[ib]`、エントロピー・
Riemann・接線速度はすべて内部 `ic`、面法線は plane の `sx/sy/sz`。乱流スカラーは `rans_neumann_scalar_boundary_d`
でゼロ勾配 (`k[ig]=k[ic]`)。

旧実装は逆流時に全圧 `Ptb/Ttb` の stagnation 流入 (=`inlet_Pressure` 構成) へ切替え、かつ速度を
`-Ux[ic]*nx` で構成していたため、剥離 BL が出口に達するケース (擬似衝撃ダクト, 壁∩出口コーナー) で
過加圧・forward↔backflow のばたつきにより発散していた。検証は
[`.github/plans/boundary-outlet-characteristic-backflow.md`](../plans/accepted/boundary-outlet-characteristic-backflow.md)。

### 入力データ構造

各 `bcond` (境界グループ) は次の GPU データを持つ。

- `map_bplane_plane_d[ib]` — グループ内インデックス `ib` → メッシュ面 ID `ip`
- `map_bplane_cell_d[ib]` — `ip` の内部セル ID
- `map_bplane_cell_ghst_d[ib]` — 対応するゴーストセル ID (`>= nCells`)
- `bvar_d[*]` — 境界面値 (`ro`, `roUx`, …, `Ts`, `Psb` ほか)
- `inputInts`, `inputFloats` — YAML から渡された定数

### カーネル構造 (例: `slip_d`)

[`slip_d`](../solver_density_cuda/cuda_forge/boundaryCond_d.cu#L3) は境界面 1 並列。

```cpp
ip = bplane_plane[ib];
ic = bplane_cell[ib];
ig = bplane_cell_ghst[ib];

Un = (sx[ip]*Ux[ic] + sy[ip]*Uy[ic] + sz[ip]*Uz[ic]) / ss[ip];

// ghost cell  ← 法線速度を反転
ro[ig]    = ro[ic];
Ux[ig]    = Ux[ic] - 2*Un*sx[ip]/ss[ip];
P[ig]     = P[ic];
roe[ig]   = P[ic]/(ga-1) + 0.5*ro[ic]*|U_L|^2;
sonic[ig] = sqrt(ga*P[ig]/ro[ig]);

// boundary value  ← 接線速度のみ (Un を 1 回引く)
Uxb[ib] = Ux[ic] - Un*sx[ip]/ss[ip];
…
```

他の BC kind も同じパターンで、`ig` (内部面で R 状態として参照) と `ib`
(境界面値で粘性束等から参照) の双方を書き込む。

### copyBcondsGradient (現状無効)

[`copyBcondsGradient_d_wrapper`](../solver_density_cuda/cuda_forge/boundaryCond_d.cu#L1439)
は境界面の勾配をゴーストセルにコピーする補助カーネルだが、現状の `boundaryCond.cpp`
からは呼び出されていない (コメントアウト)。境界勾配を厳密に扱う実装に
切り替える場合に有効化する。

### 既知の TODO / 注意点

- CPU 経路 (`cfg.gpu == 0`) は未対応。
- `inlet_fluctVelocity` は `fluct_variables_d` モジュールの変動生成と組み合わせて使う。
- 周期境界は対流再構成で 1 次風上強制が外れる (`ic1 < nCells` を維持)。
  ペアリングはメッシュ生成側で行う。
- `bcondConfig.yaml` の `kind` を新設する場合は、`boundaryCond.hpp` の
  `valueTypesOfBC` に値型テーブルを追加し、`applyBconds` にディスパッチ分岐を増やす。
