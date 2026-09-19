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

### 共役熱伝達 (CHT) — 壁温を固体と連立して解く

> **状態 (2026-09-19)**: 契約は下記のとおり確定。**実装済み = 界面診断の出力 (`output.interfaceDiag`)・
> 壁温分布の入力 (`wallProfile`)・共有 CV の壁温競合の起動時拒否**。
> **未実装 = 固体モデル・連成反復・`conjugate` 属性** (外部ループとソルバ内連成)。
> 設計判断と検証計画は [`plans/active/boundary-conjugate-heat-transfer.md`](../plans/active/boundary-conjugate-heat-transfer.md)
> (codex plan レビュー 3 巡: NO-GO → NO-GO → GO-with-changes、全件採用)。実装時は本節と実装の整合を確認する。

等温壁は $T_w$ を**与件**とするが、すきま・深いキャビティ・冷却壁では $T_w(x)$ は解の一部である。
CHT は流体の壁熱流束と固体の伝導を連立して $T_w(x)$ を決める。

#### 対応範囲 (初版)

| 項目 | 対応 |
| --- | --- |
| 離散化 | **node のみ**。cell は `wallProfile` ($T_w$ 分布の入力) までで、連成は起動時に拒否する (cell は `vizBfaceNodes` が空でシェルと 1 対 1 に対応づけられない) |
| 時間積分 | **定常陰解法** (`advanceImplicitSteady` → `implicitNonlinearUpdate`)。**dual-time 連成は対象外** (物理時間ステップ境界でのみ更新する別契約として後続) |
| 固体 | 薄肉シェル (`local1d` / `shell2d`) と一般 2D 領域 (`fem2d`) |
| 壁種別 | `wall_isothermal` + `ints: {conjugate: 1}`。**新種別を作らない** (種別名は `iso_wall_flag`・温度ピン・粘性壁・壁距離・block-DPLUR のエネルギー行切離しの 5 経路で直書き判定されており、新種別はそこから漏れる) |
| 対象外 | 表面間放射、非定常 (thin-skin 過渡)、軸対称の面内伝導、壁関数 (`wallTreatmentSST: 1` / `sstEnergyWallFunction: 1`) 併用、接触熱抵抗の同定 |

#### 界面量の定義と符号

面流束 $F^E$ は**流体 CV から外向きを正**、拘束反力 $C$ は**流体への供給を正**とする
([`plans/active/tooling-energy-balance-diagnostics.md`](../plans/active/tooling-energy-balance-diagnostics.md) と同一規約)。
壁 CV $i$ について**流体から固体へ入る熱量**は

$$Q_{f,i} \;=\; \sum_{f\in\partial_w} F^{E}_{if} \;-\; C_i \qquad [\mathrm{W}]\;(\text{平面 2D は } \mathrm{W/m})$$

- 定常の Dirichlet 行では $C_i=-R_i^{raw}$。**過渡では $C_i = D_t(V_iE_i) - R_i^{raw}$** で、定常式を瞬時入熱に使えない。
- 検算: $\sum F^E=80$, $C=-20$ なら $Q_f=100$。
- **これが連成の正本**である。node 等温壁は温度ピン後に `res_roe` を 0 化するため、
  壁 CV に実際に入った熱は「壁面の物理境界流束 + 拘束反力」であり、
  次の 2 つは**診断**として併記するだけで界面には渡さない。
  - **コンパクト差分形** $k_{\rm eff}(T_1-T_w)/d_1$ (固体向き正) — $D_f$ の推定と精度診断に使う。SU2 CHT の界面転送と同じ形。
  - **再構成勾配形** $k_{\rm eff}\nabla T\cdot\mathbf S$ — `viscousFlux_d.cu` が `qwall` に保存している値。
  - 実測差の例: case/48 `run_0011` の $x\approx0.5$ m で コンパクト 96.184 / 2 次片側 98.820 kW/m² (2.67 %)。
- 界面の積分 (面積重み、軸対称の $r$ 重み) は **host・double** で行う。
- 幾何は **primal facet 単位** (`bc.vizBfaceNodes`) を正本にする。node の合成半割面ベクトルは
  $|\sum_f\mathbf S_f|\ne\sum_f|\mathbf S_f|$ なので**面積として使わない**。

#### 固体モデル

未知数は**ガス側表面温度** $T_w$。背面環境 $T_b$ までの**全抵抗**は

$$R_{\rm tot}=\frac{t}{k_s}+R_{\rm back},\qquad
R_{\rm back}=\begin{cases}0&\text{背面等温}\\ 1/h_c&\text{冷却剤}\\ \sum_i t_i/k_i&\text{多層}\\ \infty&\text{断熱}\end{cases}$$

で、**背面等温でも $t/k_s$ を落とさない** (落とすと $T_w=T_b$ に退化する)。シェル方程式は

$$\nabla_{\!s}\!\cdot\!\left(k_s t\,\nabla_{\!s}T_w\right)+q_{\rm gas}-\frac{T_w-T_b}{R_{\rm tot}}=0 .$$

- **断熱・孤立系** ($R_{\rm tot}=\infty$ かつ端部断熱) は定数零空間を持ち、正味入熱が非零なら定常解が無い →
  起動時に拒否するか、適合条件と零空間の固定を明示する。**「SPD なので CG」は Robin 項がある構成に限る**。
- `fem2d` では未知数が固体全節点 $u$ になる。界面抽出を $E$、固体剛性を $K_s$ として

  $$\left(K_s+E^{\mathsf T}D_fE\right)u^{k+1}=b_s+E^{\mathsf T}\!\left[Q_f(Eu^{k})+D_f\,Eu^{k}\right]$$

  とし、流体と固体の外周節点を一致させる (補間を挟まない)。$Q_f$ は**積分済み節点荷重**なので
  $E^{\mathsf T}$ で載せるときに**面積を再乗算しない**。共有角の反力は**一度だけ**計上する。

#### 反復と受理判定

共役定常解は $A_sT^{*}=b_s+Q_f(T^{*})$。反復は**固定点を保存する**形で書く:

$$\left(A_s+D_f\right)T^{k+1}=b_s+Q_f(T^{k})+D_f\,T^{k}$$

- $D_f$ は収束速度だけを決め、**固定点は $D_f$ に依らない**。初期推定は $D_f^{(0)}=k_{\rm eff}A/d_1$ だが、
  これは**上界ではない** (実効応答 $H=-\partial Q_f/\partial T_w$ は非対角を持ち、発散する反例がある)。
- **受理はメリット関数の降下で判定する**。未緩和残差 $r^k=A_sT^k-b_s-Q_f(T^k)$ に対し、
  比較の間は重みを固定した $\Phi(r)=r^{\mathsf T}(A_s+D_f)^{-1}r$ を使い、降下しなければ line search → $D_f$ 増加 →
  再試行上限で失敗を報告する。**残差最大ノルムの単調減少を受理条件にしない** (収束する反復を棄却する反例がある)。
  **$\Delta\Phi$ が丸め以下の停滞を合格にしない**。局所最大ノルムは最終ゲート (下記 G-if) に使う。

#### 壁温分布の入力 (`wallProfile`)

`Ts` は `valueTypes==1` の **per-face bvar** で、起動時に YAML の一様値で 1 度埋めた後は
カーネルが書き換えない。したがって面ごとに違う $T_w$ を入れればそのまま効く
(cell ゴースト `wall_isothermal_d`、node 温度ピン `pin_wall_node_temperature_d` とも `Tsb[ib]` を読む)。
`ints: {wallProfile: 1}` で `wall_profile_<physID>.csv` から埋める (入口分布 `applyInletProfiles` の一般化。
CSV 書式は入口分布と同じで、先頭の連続 x/y/z 列が補間座標、残りが量名 = 壁では `Ts`)。

- **補間位置**: node は**ノード座標**、cell は面重心。
  face 重心をそのまま node に使うと位置がずれる (case/48 `run_0011` で実測 0.679 mm)。
- **verify は `bvar` の再出力では不十分**。壁ダンプの `Ts` は入力の再表示なので、
  場に入ったことは `VALUE/T` と EOS 整合まで見て確認する。
- CHT 内部の転送は座標補間でなく**安定なノード ID** を正本にする。
- **実測 (case/48 `run_0015_wallprofile`)**: $T_w = 300+100x$ を与えると、**場の `VALUE/T` が壁ノードで
  同じ分布になる (最大差 0.068 K)**。`bvar` の `Ts` を見るだけでは「入力が入ったこと」しか分からない。
- **入口分布 (`inletProfile`) は従来どおり face 重心で補間する** (既存 run のビット不変を守るため)。
  node の入口にも同じ位置ずれの問題はあるので、必要になったら別途 opt-in で直す (本節の契約の範囲外)。

#### 起動時に拒否する構成

契約を散文で守らず、次はいずれも**起動時にエラーで落とす**。

1. `cell` 離散化で `conjugate: 1`。
2. **温度を拘束する壁どうしが CV を共有**していて、**そこに与える壁温が食い違う**場合
   (温度ピンは bcond 順に適用され、角ノードは複数 bcond に重複するので**後勝ち**になる)。
   連成の有無を問わず拒否する。**実装済み** (`conjugateWall::checkWallTemperatureSharing`):
   壁温を陽に扱っている run (`wallProfile` か `interfaceDiag` が有効) でのみ走り、
   競合 CV の数・physID・それぞれの $T_w$ を出して起動時に落とす。それ以外の run は挙動不変。
   検証: case/48 の平板で `sym` (前縁上流の対称面) を $T_s$=500 K の等温壁にすると、
   前縁で共有する 1 CV を検出して exit 1 する。
3. 断熱・孤立固体で正味入熱が非零 (定常解が無い)。
4. 未定義の構成: 周期同一視・軸対称の面内伝導・壁関数併用・dual-time 連成。

#### 診断出力とゲート

**実装済み**: `output: {interfaceDiag: 1}` (既定 0) で、壁 (`wall` / `wall_isothermal`) のダンプ
`res_wall_<physID>_*.h5` に次を追加する (host 側で作るので device の `bvar` を汚さず、既定 run の出力は不変)。

| データセット | 中身 |
| --- | --- |
| `iface_T1` / `iface_d1` | 第一内部点の温度と法線距離 (定義は上記。`tools/check_wall_resolution.py` と同一規則) |
| `iface_keff` | $k_{\rm eff}=k_{\rm lam}+c_p\mu_t/Pr_t$ (viscousFlux の壁経路と同じ) |
| `iface_q_compact` | $k_{\rm eff}(T_1-T_w)/d_1$ (**固体向き正**) |
| `iface_q_recon` | $-$`qwall` = viscousFlux が残差に入れた再構成勾配形 (**固体向き正**に反転済み) |
| `iface_q_2nd` | 3 点非等間隔の 2 次片側差分による $k_{\rm eff}\,dT/dn$ |
| `iface_ok` / `iface_align` | 第一内部点が定まったか / 整列度 $|d\cdot\hat n|/|d|$ |

**実測 (case/48 `run_0014_iface_diag`, 壁法線に整列した node メッシュ, $y_1$=3 µm)**:

- `iface_d1` = 3.0001 µm = **第一層厚と一致**、`align` = 1.000、1001/1001 点が評価可。
- **`q_recon` は `q_compact` と 1.6e-7 相対で一致する**。この配置では再構成勾配がコンパクト差分に帰着するため。
- **`q_2nd` は `q_compact` と中央値 2.79 % 違う** (最大 3.1 %)。つまり過去に見えた ~2.7 % の食い違いは
  **後処理の差分形式の選択**であって、ソルバ内部の不整合ではない。滑らかな分布では 2 次片側の方が正確なので、
  **カーネルの壁熱流束はこの解像度で ~3 % の 1 次打ち切り誤差を持つ** (誤差予算に入れる)。
- 斜交・非整列メッシュでは 3 つが分かれるはずなので、**どれを見ているかを列名で明示する**。

$q_{\rm eff}$ (拘束反力込み、= 連成の正本) は依存診断の完成後に追加する。

- **G-cons (熱収支)**: $\varepsilon=\big|\sum_i Q_{f,i}-(\text{固体正味入熱})\big|$ を
  $\max(\sum_i|Q_{f,i}|,\,Q_{\rm floor})$ で規格化する (正味量で割ると符号相殺で分母が消える)。
- **G-if (界面収束)**: ① 局所面積で規格化した残差 $\max_i|r_i|/A_i$ [W/m²] の絶対条件、
  ② 相対条件 $\max_i|r_i|/\max_i|Q_{f,i}|$、③ 必要連続反復数、の 3 つを**別々に**満たすこと。
  欠損・非有限・量子化停滞は不合格 (float32 では 1000 K 付近の 1 ULP が $6.1\times10^{-5}$ K)。
- 派生量の準定常判定は **drift と振動幅の両方**を比較許容の 1/5 以下にする
  (`check_quasisteady.py` の既定 5 % / 10 % では 0.2–0.5 % の比較を支えられない)。

拘束反力 $C_i$ の採取は [`tooling-energy-balance-diagnostics`](../plans/active/tooling-energy-balance-diagnostics.md) が提供する。
同診断は初版で dual-time・周期・軸対称を対象外としているため、**dual-time は 1 次元検証の前、周期は翼列検証の前**に
解除するマイルストーンを同 plan の残作業に登録済み。

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
