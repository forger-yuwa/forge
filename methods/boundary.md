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

**速度の大きさは新しい音速に整合させる** (2026-09-20 修正、`boundaryCond_d.cu` CPG 分岐)。
$M_L$ は**内点の音速**で測った値なので、上式で $T_R$ を作ると音速が変わる。ここで速度を内点値のまま
残すと境界状態の実マッハが $|u_L|/a_R \ne M_L$ になり、**指定した $T_t$・$P_t$ を再現しない**。

$$|u_R| = |M_L|\,a_R,\qquad a_R = \sqrt{(\gamma-1)c_p T_R},\qquad \mathbf{u}_R = -|u_R|\,\mathbf{n}.$$

反例 ($\gamma$ 1.4, $c_p$ 1005, $T_L$ 100 K, $|u_{n,L}|$ 100 m/s、指定 $T_t$ 293.15 K / $P_t$ 100 kPa):
修正前は $T_t$ **284.232 K** / $P_t$ **89751 Pa** を返していた (3.0 % / 10.2 % 低い)。修正後は厳密に指定値を返す。
TP 分岐 (`thermalMethod: 2`) は `thermo_isentropic_from_total_*` の $|u|$ で速度を作り直しており、
**CPG 分岐だけが取り残されていた**。

**⚠ 定常の設計点では現れない**: 内点が指定全条件と等エントロピーで整合していると $a_R = a_L$ になり
不整合が消える (`case/08.bump` は修正前でも境界 $T_t$/$P_t$ の誤差 0.000 %)。
**過渡・オフデザインでのみ出る**ので、検査は過渡で行うこと。

**⚠ 超音速流入に使ってはいけない**: この閉包は内点から境界を決めるので、全特性が流入する超音速入口では
**$\rho\downarrow \to u\uparrow \to P_s\downarrow \to \rho\downarrow$ の正帰還**になり、
**入口の 1 節点だけが十数 step で枯れて発散する**。超音速入口は `inlet_uniformVelocity` を使う
(詳細と指紋は [`procedures/divergence-and-startup.md`](../procedures/divergence-and-startup.md))。

**`inlet_Pressure_dir` の退避 2 点** (2026-09-20 修正、`boundaryCond_d.cu`)。どちらも過渡で
入口面の一部の節点だけが非有限になる形で出る (実測: case/53 翼列で step 214、49 節点中 24 節点)。

1. **$P_c > P_t$ の根号**: 境界マッハは内点静圧から
   $M_b=\sqrt{2\{(P_c/P_t)^{-(\gamma-1)/\gamma}-1\}/(\gamma-1)}$ で作るが、起動過渡で内点静圧が
   指定全圧を超えると根号内が負になり NaN。亜音速全圧入口としては $M=0$ が正しい極限なので
   0 に落とし、逆に $P_c$ が落ち込んだときの暴走を避けるため $M_b\le1$ に制限する。
2. **方向ベクトルの消失**: `bvar` の `Ux/Uy/Uz` は config 指定の**方向**ベクトル (`valueTypes==1`) だが、
   カーネルは毎 step これを**次元付き速度で上書き**する。$M\to0$ で長さが 0 になると次 step の
   正規化が 0/0 になり、以後ずっと NaN。長さが $10^{-10}$ 未満なら面法線 (内向き) を方向に使う。

**用途上の注意**: `inlet_uniformVelocity` は「速度 3 成分 + config エントロピー ($P_s/\rho^\gamma$)」を
課し $R^+$ を内点から取る正しく posed な亜音速入口だが、**全温を固定しない**。内部静圧が config の
アンカー状態からずれるとその分だけ $T_t$ がずれる (case/53 で入口 $T_t$ が 786→738 K と 48 K 低下)。
$T_t$ 基準の熱伝達を比較する run では `inlet_Pressure_dir` (全圧・全温指定) を使う。

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

#### 等温壁のエネルギー境界: 強制 と 弱形式 (`mesh.nodeIsothermalEnergyBC`)

連続系の条件は同じ $T|_{\rm wall}=T_w$ だが、**離散化に 2 通りある**。既定は従来どおり強制 (0)。

| | 強制 `0` (既定) | 弱形式 `1` (opt-in, SU2 型) |
| --- | --- | --- |
| 壁ノード $T$ | BC 値へ上書き (`pin_wall_node_temperature_d`) | 上書きしない (方程式を解く) |
| 壁エネルギー残差 | ゼロ化 (`zero_res_roe_bplane_d`) | 残す |
| 陰解法エネルギー行 | 単位行 | 対角ブロックに線形化を加算 |
| 壁半割面の伝導 | $\nabla T\cdot S$ (ghostless) | $k_{\rm eff}(T_I-T_w)/d_1\cdot A_{\rm half}$ で**置換** |
| 運動量 no-slip | 強制 | 強制 (変えない) |

弱形式の残差寄与は壁半割面ごとに

$$R^{\rm roe}_W \mathrel{-}= k_{\rm eff}\,\frac{T_I-T_w}{d_1}\,A_{\rm half},\qquad
k_{\rm eff}=k_{\rm lam}+\frac{c_p\mu_t}{{\rm Pr}_t}$$

で、$I$ は第一内部点 (壁法線との alignment 最大の非壁隣接)、$d_1$ は**法線投影距離**
$\lvert(\mathbf x_I-\mathbf x_W)\cdot\hat n\rvert$。**点間距離ではない** — 法線から 30° 傾いた辺では
点間距離版の熱流束が 13.4 % 小さくなる。**$T_W$ 自身は使わない**。

陰解法には**近似対角項**を足す。$g=k_{\rm eff}A_{\rm half}/d_1$ と置くと壁寄与は
$R_W^{\rm wall}=-g(T_I-T_w)$ なので、**$(\rho e)_W$ による厳密微分は 0** であり温度微分は
内部点 $I$ の列にある。SU2 も同じく内部点温度による残差に対して壁点の対角項を加える
**近似線形化**である。forge は `res_roe` を右辺に取り行列へは符号を反転して組むので、
CPG でのエネルギー対角への追加は

$$\Delta A_{WW}^{\rm energy}=+\frac{g}{\rho\,c_v}$$

**初版は `thermalMethod: 0` (CPG) 限定**。TP は $e(T)$ を使うので密度微分の CPG 式を一般化できない。

**壁半割面が複数ある壁ノード**では、半割面ごとに残差と対角の両方を `atomicAdd` で積む (回数が構造的に一致する)。
**内部側は 1 点**で、複数の内部隣接に分配しない (SU2 と同型。診断 `iface_q_compact` と CHT の
$D_f=k_{\rm eff}A/d_1$ が同じ 1 点を使うので、BC・診断・連成が同じ幾何を見る)。

**初版に over-relaxed 非直交補正は入っていない** — 内部面と違い純粋な 2 点差分である。
非直交補正には壁面の接線勾配が要り、弱形式では隣の壁ノードの $T$ が自由 DOF なので**隣接壁ノード同士が
再結合する**。市松を減衰させる方向かもしれないし、旧弱形式の「CV 平均に緩む」問題を呼び戻すかもしれないので、
**第 2 の変数**として分離し最初の A/B には混ぜない (plan §5.1 #6b)。

**2026-07-20 に棄却した旧弱形式とは別物**である。旧実装は壁ノード $T$ を浮かせ、**その緩んだ
$T$ で壁流束を作った**ので壁 CV 平均へ落ち (~0.1 K オフセット)、第 1 スペーシング勾配が $-24\,\%$
だった。本節の弱形式は流束を**指定 $T_w$** から作るのでこの経路が無い。SU2 の実装は
`SU2_CFD/src/solvers/CNSSolver.cpp` の `BC_Isothermal_Wall_Generic` で、コメントに
"Apply a weak boundary condition for the energy equation" と明記され、
運動量だけ `Jacobian.DeleteValsRowi` で強制している。実測で SU2 の壁ノード保存温度は
566.011–566.886 K (指定 566 K)。

**動機**: 強制の既知の副作用 3 つ — (a) 壁熱流束の節点交番が同一メッシュの SU2 の **20 倍**
(C3X 負圧面 0.919 % vs 0.046 %)、(b) 第 1 スペーシング勾配の $-15\,\%$ バイアス、
(c) エネルギー行 decouple による擬似 CFL 上限 $\sim5$ — がこの閉包に帰属するかを切り分ける。
対流スキーム (ROE 0.901 % / HLLE 2.10 %)・低マッハ前処理 (破綻)・float32 の丸め
(観測 0.083 K は 1360 ULP) は**いずれも棄却済み**。
計画は [`plans/active/boundary-weak-isothermal-wall.md`](../plans/active/boundary-weak-isothermal-wall.md)。

**CHT との整合**: 保存的界面熱量 $Q_f=\sum F^E-C$ の拘束反力 $C=-R^{\rm raw}$ は、弱形式では
**拘束が無いので $C=0$**。$Q_f=\sum F^E$ がそのまま収支に一致する。`iface_q_eff` の式は変えない。

**併用不可 (起動時に拒否)**: `wallTreatmentSST: 1`、`wallModelLES: 1`、cell 方式、`thermalMethod != 0`、
軸対称、壁ノードが周期で同一視される構成、移動壁、内部点なし / $d_1$ 退化 / 選ばれた $I$ が別の壁ノード。
`nodeWallDirichlet: 1` は必須条件。**`nodeWallDirichlet: 0` で代用してはならない** (運動量まで弱くなる)。

**置換するのは `viscousFlux_wall_d` の壁半割面の伝導だけ**で、内部双対面には触らない。
`Qw_Wall` は W–I **内部双対面**を置換する別機構で、`> -0.5` を有効判定に使いマーカと符号付き値を
兼用しているため、冷却壁 (負の流束) に流用できない。

**診断との整合**: `iface_q_compact` は壁ノードの**保存温度**を使うので、強制側では $T_W=T_w$ で隠れるが
弱形式では境界流束と別量になる。比較用のコンパクト熱流束は**両枝とも指定 $T_w$** で定義し、
保存温度からの勾配は別名の診断量にする。`Tw_bc` と `T_W` を別々に壁ダンプへ記録する。

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

> **状態 (2026-09-23)**: 契約は下記のとおり確定。
> **実装済み**: 界面診断の出力 (`output.interfaceDiag`)、壁温分布の入力 (`wallProfile`)、
> 共有 CV の壁温競合の起動時拒否、**固体側モデル** ([`tools/solid_shell.py`](../solver_density_cuda/tools/solid_shell.py))、
> **外部弱連成ループ** ([`tools/cht_loop.py`](../solver_density_cuda/tools/cht_loop.py))、
> **ソルバ内連成** (`conjugate:` ブロック + bcond `ints: {conjugate: 1}`): `mode: local1d` (点ごとの 1 次元抵抗) と
> **`mode: fem2d`** (一般 2D 固体を全節点系で解く。下の節)。界面熱量は既定で**保存形** $q_{\rm eff}$。
> **未実装**: ソルバ内の `shell2d` (帯メッシュを `fem2d` に食わせる方針)、dual-time 連成。
> 検証: 1 次元純伝導の共役解を解析解と照合 (`case/52.conjugate_slab`) — $T_w$ 誤差 0.025 %、
> 両側 $q$ の不一致 0.0053 % で **PASS**。
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
**離散化 (`shell2d`, 実装は [`tools/solid_shell.py`](../solver_density_cuda/tools/solid_shell.py))**:
境界面 (primal facet) の線形 FE + 集中質量。**四角形は 2 通りの対角線分割を 1/2 ずつ使う** —
片方だけで割ると集中面積が非対称になり、**境界節点に O(1) の荷重不均衡**が残って
**収束次数が 2 次 → 1 次**に落ちる (帯フィン問題の実測: rate 1.00 → 対称化で 2.00)。
**断熱・孤立系は定数零空間を持つ**ので、正味入熱が非零なら例外にする (解が無い)。

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
  **比較の間は重みを固定した** $\Phi(r)=r^{\mathsf T}(A_s+D_f)^{-1}r$ を使い、降下しなければ line search → $D_f$ 増加 →
  再試行上限で失敗を報告する。**残差最大ノルムの単調減少を受理条件にしない** (収束する反復を棄却する反例がある)。
  **$D_f$ を動かしながら $\Phi$ を比べてはいけない**: 重み変更による見かけの降下で受理が通り、
  $D_f$ が発散的に増えて停滞する (実測)。**収束判定は物理量** (max$|\Delta T|$ [K] と max$|r|$/スケール) で行う。
- **加速器 (Anderson) を既定に含める**。スカラー $D_f$ は非対角な流体応答を表せないので、
  素の固定点反復は収束しないことがある (実測: 三重対角 SPD の応答で 200 反復未収束 → Anderson 深さ 5 で 57 反復)。
  加速候補は $\Phi$ が降下しなければ棄却し、$D_f$ を変えたら履歴を捨てる。
  **$\Delta\Phi$ が丸め以下の停滞を合格にしない**。局所最大ノルムは最終ゲート (下記 G-if) に使う。

#### ソルバ内連成 (`conjugate:`, Phase 2a)

```yaml
conjugate: {mode: local1d, flux: q_eff, thickness: 1.0e-3, k_solid: 0.217, back: isothermal,
            T_b: 300.0, interval: 50, warmup: 500, relax: 1.0}
```

と書き、対象壁の bcond に `ints: {conjugate: 1}` を付ける (種別は `wall_isothermal` のまま)。
`interval` step ごとに、**ステップ完了後** (次の残差組立ての前) に壁温を更新する。

**`flux: q_eff` (既定、保存形)** — 上の「界面熱量」で定義した $q_{\rm eff}$ を使い、
[plan §4.2](../plans/active/boundary-conjugate-heat-transfer.md) の**固定点を保存する形**で解く:

$$(g_s + D_f)\,T_w^{k+1} = g_s T_b + q_{\rm eff}(T_w^k) + D_f\,T_w^k,\qquad
  g_f = \frac{k_{\rm eff}}{d_1},\quad g_s = \frac{1}{R_{\rm tot}},\quad D_f = g_f$$

収束点は $g_s(T_w - T_b) = q_{\rm eff}$、すなわち**保存形の界面熱量**と固体の 1 次元法則の釣り合いで、
外部ループ (`cht_loop.py --flux q_eff`) と**同じ不動点**である。$D_f$ は界面抵抗の初期推定で、
収束速度だけを決める (固定点は $D_f$ に依らない)。**`output: {interfaceDiag: 1}` が要る** (無ければ起動時に拒否)。

**`flux: q_compact` (旧実装、A/B 専用)** — 抵抗加重平均 $T_w^{new} = (g_f T_1 + g_s T_b)/(g_f+g_s)$
(SU2 の `AVERAGED_TEMPERATURE` と同型)。これは $q_{\rm compact}$ の固定点であり、**保存形とは一致しない**。
差は壁半 CV 内の粘性加熱 $\tau\cdot u$ と流動仕事で、第一層厚 $d_1$ に比例する。
**実測** (`case/48.flat_plate_cooled_m4`, $d_1$=3.0 µm、[plan §5.1 #66](../plans/active/boundary-conjugate-heat-transfer.md)):
同一状態の G-cons が $q_{\rm eff}$ 形の更新では **0.0028 % (PASS)**、$q_{\rm compact}$ 形では **1.77 % (FAIL)**。
壁温は平均 599.96 → **604.17 K**、前縁の最大 987.1 → **1148.5 K** と動く。

- **起動時に拒否**: `node` 以外、`unsteady: 1` (dual-time)、`mode != local1d`、背面断熱、
  対象壁が `wall_isothermal` でない、第一内部点が定まらない壁 CV が 1 つでもある場合、
  `flux: q_eff` なのに `interfaceDiag != 1`、$q_{\rm eff}$ が 1 節点でも非有限 (適用範囲外の構成)。
- **界面の収束判定 (G-if)**: 更新ごとに run 直下の `conjugate_history.csv` に
  `step, physID, n, Tw_mean, Tw_min, Tw_max, dTw_max, res_abs_Wm2, res_max_W, res_rel, q_total` を追記する。
  $r_i = Q_{f,i} - g_s(T_{w,i}-T_b)A_i$ は**更新前 (未緩和)** の界面残差で、`res_abs_Wm2` $=\max_i|r_i|/A_i$、
  `res_rel` $=\max_i|r_i|/\max_i|Q_{f,i}|$。判定 (許容と連続回数) は判定ツール側で行う。
- **面内伝導が要るなら `mode: fem2d`** (下の節)。`local1d` は点ごとの 1 次元抵抗で、面内伝導を落とした極限。
  `shell2d` は外部ループ ([`tools/cht_loop.py`](../solver_density_cuda/tools/cht_loop.py) + `solid_shell.py`) の担当のまま。
##### `mode: fem2d` — 一般 2D 固体をソルバ内で解く (2026-09-23)

```yaml
conjugate:
  mode: fem2d
  solid: solid.h5        # tools/solid_mesh_to_h5.py が作る (メッシュ・孔 Robin・k_s(T) の正本)
  flux: q_eff            # 既定。保存形の界面熱量
  interval: 50           # K step ごとに更新
  flux_avg: 42           # 界面熱量を N 更新の後方移動平均にする (既定 1 = 平均しない)
  Df_scale: 1.0          # D_f の倍率 (発散したとき手で上げる。自動調整はしない)
  refactorDT: 1.0        # [K] 固体温度がこれ以上動いたら分解し直す
  gate: {eps_rel: 1.0e-3, eps_abs_Wm2: 150.0, dT_K: 1.0e-2, tol_solid: 1.0e-2, n_consec: 80}
```

毎更新、**固体の全節点系を 1 回解く** (Schur 補元は作らない)。解き方は**残差補正形**:

$$r = K_s(u^k)u^k - b_s - E^{\mathsf T}\bar Q_f,\qquad
  \Delta = -\bigl(K_s(u^{\rm fact}) + E^{\mathsf T}D_fE\bigr)^{-1} r,\qquad u^{k+1} = u^k + \Delta$$

残差 $r$ は**常に現在の $k_s(u^k)$** で組むので、固定点は $K_s(u)u = b_s + E^{\mathsf T}\bar Q_f$ になり、
分解を再利用しても物性が凍らない (左辺は前処理としてしか効かない)。
$r$ は **$D_f$ を含まない形で直接組む** — 相殺に頼ると界面温度の精度差 $D_f(Eu-T_s)$ が残り、
$D_f$ が大きいとき量子化で止まった状態を合格にできる。

- 未知数は**全節点温度** $u$ なので内部温度が状態になり、$k_s(T)$ の自己整合に「復元してから組み直す」操作が要らない。
- $D_f$ は**界面対角のみ** $g_fA_i$ ($g_f=k_{\rm eff}/d_1$、$A_i$ は**固体側の集中辺長**)。非対角性は左辺の $K_s$ が持つ。
  **$Q_f$ は流体側の積分済み荷重 `iface_Qf_eff` $=R^{raw}-F_w-e_wR_\rho$ をそのまま渡す**。
  面積で割って固体側の集中辺長を掛け直すと、両者が違う角で荷重が歪む (Mark II 後縁で +29.9 %)。
  集中辺長は**熱流束に換算する段でだけ**使う。`surfArea` を荷重に使わない (押し出し疑似 2D で奥行きが乗る)。
  初版は $z\equiv0$ の平面 2D 以外を拒否する。
- 行列は SPD なので**下三角バンド Cholesky** (RCM 並べ替えは変換時に済ませる)。
  **分解は `refactorDT` 以内なら再利用する** (上の残差補正形なので固定点は動かない)。
  **直接求解 $Au^{k+1}=b$ のまま再利用してはいけない**: 行列側の $D_f$ や $k_s$ が古いまま右辺だけ新しくなり、
  固定点がずれる (実測: `res_rel` 1.8e-4 → 1.5e-2、窓内の壁温の振れ 2.34 K)。
- **受理判定・line search・Anderson は持ち込まない**。ソルバ内の $Q_f(T)$ は step ごとに動く写像なので、
  更新間のメリット比較は同じ関数の比較にならない (外部ループが凍った機構)。
- **`flux_avg`**: 流体側に局所振動があると瞬時の $Q_f$ では界面ゲートが床に当たる。
  C3X は吸込面の $k$ オンセット前線が**周期 1012 step** で揺れ、最悪節点の $Q_f$ が中央値の 78 倍ばらつく。
  窓は「$F_N\le\epsilon_{\rm abs}L_i/2$ を全節点で満たす最小 N」で選ぶ (C3X では 21 = 1 周期。周期揺らぎ対策で 42 を採用)。
- **出力**: 更新ごとに `conjugate_history.csv` (G-if の素材)、流体の出力間隔で `res_solid_<physID>_<step>.h5`+`.xmf`
  (`T` / `k_s` / `q_iface` / `q_hole`) と再開用の `conjugate_state_<physID>.h5`。
- **起動時に拒否**: 平面 2D 以外、界面節点が流体の壁節点と 1 対 1 でない (**内挿しない**)、Robin 辺が 1 本も無い
  (定数零空間)、$q_{\rm eff}$ が 1 節点でも非有限、`flux: q_eff` なのに `interfaceDiag != 1`。
- **判定**: 界面の収束は [`tools/check_cht_interface.py`](../solver_density_cuda/tools/check_cht_interface.py) が
  `conjugate_gate.json` の**事前登録値**で行う (`check_convergence.py` は流体の保存量しか見ない)。
  **帯を限った判定**: `conjugate: {node_log: 1}` にすると、更新ごとに界面全節点の
  更新前物理残差 $r_i$ [W/m]・流体荷重 $Q_{f,i}$・壁温の変化 $\Delta T_i$・$T_{w,i}$ を `conjugate_iface_log_<physID>.csv` に、
  節点の $y_i$・集中辺長 $A_i$ を `conjugate_iface_nodes_<physID>.csv` に書く (host 側の出力のみで更新式は不変)。
  `check_cht_interface.py --band-y YTOP YBOT` がこれを読み、①②③を**帯内節点だけ**で、④固体内部残差は全域で判定する
  (結果は `CHT_INTERFACE_BAND_VERDICT.txt`)。節点表の有限性・正の面積、ID の連番・一意性・完全性、更新内の時刻一致、
  履歴との更新数の一致を検査し、崩れていれば `REFUSED`。全域 max だけを持つ `conjugate_history.csv` では、
  lip 行のような局所のリップルが深部の判定を支配するため (CHT plan §5.1 #101、`case/58`)。
  **固体が受け持つ熱の保存検査**は固体の物理作用素から $Q_{\rm sol}=(K_su-b_s)_{\rm iface}$ を組んで $Q_f$ と比べる。
  固体ダンプの `q_iface` は**流体荷重 $Q_f$ のコピー**なので、これを使うと渡した荷重を渡した荷重と比べることになる
  (CHT plan §5.1 #100)。

- **再開**: 出力ステップごとに `conjugate_Tw_<physID>.csv` を書く。続きを回すときは
  これを `wall_profile_<physID>.csv` にコピーして `ints: {conjugate: 1, wallProfile: 1}` にすると、
  収束した壁温から再開できる (`wallProfile` が初期値、`conjugate` がその後の更新)。
  **`mode: fem2d` では CSV だけでは足りない** — 固体の内部温度・平均バッファ・更新位相・累積 step は
  `conjugate_state_<physID>.h5` にあるので、**これも次の run ディレクトリへコピーする**。
  起動時に `content_sha1` を照合し、違えば拒否する。収支ゲートも
  **このチェックポイントを必須**とし、無ければ `REFUSED` にする (初期壁温から復元した固体で代替しない)。
- **実測 (case/52, V1 の 1 次元共役解)**: 解析解 $T_w$=316.2618 K に対し
  **ソルバ内 316.2371 K (温度上昇の −0.152 %)**、外部ループ (shell2d) 316.2659 K (+0.025 %)。
  両側 $q$ の不一致はそれぞれ **0.0002 % / 0.0053 %**。両者の差 0.029 K は
  **流体側の離散解の差** (同じ壁温での $q$ が ±0.2 % 動く) の範囲。

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


#### 界面熱量の 4 つの定義 — **保存的な `iface_q_eff` が正本**

`output.interfaceDiag: 1` の壁ダンプ (`res_wall_<physID>_*.h5`) は 4 つ出す。**どれを見ているかを
必ず明記する** (実測で 2〜19 % 違う)。符号はすべて**固体向きが正**。

| 列 | 定義 | 用途 |
| --- | --- | --- |
| **`iface_q_eff`** | $Q_f=\sum F^E_{\partial w}-C$ を面積で割ったもの。**正本** | 連成の熱量 (`cht_loop --flux q_eff`)、収支ゲート |
| `iface_q_compact` | $k_{\rm eff}(T_1-T_w)/d_1$ | $D_f$ の初期推定・精度診断 |
| `iface_q_recon` | 再構成勾配 (`qwall` の符号反転) | 診断 |
| `iface_q_2nd` | 3 点非等間隔の 2 次片側差分 | 診断 |

**SU2 と比べるときだけはコンパクト差分を使う** (2026-09-21): SU2 の表面出力 `Heat_Flux` は
等温境界では `HeatFlux = thermal_conductivity * (There - Twall) / dist_ij`
(`SU2_CFD/include/solvers/CFVMFlowSolverBase.inl:2638`) で、`iface_q_compact` と**同一構成**である。
SU2 には保存形の対応物が無いので、`iface_q_eff` と `Heat_Flux` を並べるのは別量の比較になる。
実測 ($h$ は金属の熱収支から出る) との比較は `iface_q_eff`、SU2 との比較は `iface_q_compact` 同士、
と使い分けること。C3X run 108 で両者は平均 2.2 % 違う (保存形が大きい)。

**なぜコンパクト差分ではいけないか**: 壁 CV に実際に入った熱は `viscousFlux` の再構成勾配と
内部面の離散の和であって、$k_{\rm eff}(T_1-T_w)/d_1$ ではない。これを連成に使うと**熱量が
閉じない解を合格させてしまう**。実測 (case/53 C3X 翼列, 480 節点): `q_eff` は `q_compact` に対し
**bias +2.26 % / rms 3.28 % / 局所最大 19.1 %** 違う。

**$Q_f$ の作り方** (符号規約はここが正本。散文で保証せず V1 で検算する):
面流束 $F^E$ は**流体 CV から外向きを正**、拘束反力 $C$ は**流体への供給を正**とすると、壁 1 節点で

$$Q_{f,i}=\sum_{f\in\partial_w}F^{E}_{if}-C_i,\qquad \text{定常の Dirichlet 行では } C_i=-R_i^{raw}$$

実装は残差の言葉で次の 2 つを採る。

- `ifaceFw[ib]`: **壁半割面が `res_roe` に入れた寄与そのもの** ([`viscousFlux_d.cu`](../solver_density_cuda/cuda_forge/viscousFlux_d.cu) の
  `res_roe_temp`)。残差規約が $R=-\sum F$ なので $\sum F^E_{\partial w}=-\,$`ifaceFw`。
- `ifaceRraw[ib]`: **壁残差射影の直前**の `res_roe` ([`nodeWallDirichlet_d.cu`](../solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu) の
  `zero_res_roe_bplane_d` が 0 化する前に退避)。**後から `res_roe` を読んでも 0 しか出ない**ので、
  この位置で採るしかない。

よって `iface_q_eff = (ifaceRraw - ifaceFw) / 面積`。これは恒等的に
$-\sum_{f\in \text{内部面}}F^E$ (= 流体内部から壁 CV へ入る正味熱) に等しく、
**壁面流束の離散化に依らず保存する**。

**未収束の擬似時間での蓄積項 (2026-09-21)**。定常解では壁 CV の質量残差 $R_\rho$ は 0 だが、擬似時間が落ちきらない場
(翼の衝撃足のように解が動き続ける所) では 0 でない。壁 CV は固定体積・$u=0$・$T=T_w$ なので、そのエネルギーは
$d(V\rho e_w)/d\tau=e_wR_\rho$ だけ変わる。これは壁 CV の**蓄積**であって固体へ渡る熱ではない (下の dual-time の $D_t(VE)$ と同じ位置づけ)。そこで

$$\texttt{iface\_q\_eff}=\frac{R^{raw}-F_w-e_wR_\rho}{A},\qquad e_w=\frac{\rho E}{\rho}\Big|_w$$

**係数は $H_w$ ではない**。対流エネルギー残差は $H_wR_\rho$ として現れる (実測: `case/53.c3x_vane_cht/run_0124_fxhalf_rro` の衝撃足で
相関 1.000、傾き 571 kJ/kg = $c_pT_w$) が、差 $(P/\rho)R_\rho$ は等温のまま質量を押し込む流動仕事で、壁が実際に受け取る熱である
(codex 2026-09-21)。block-DPLUR の実際の更新量 $V\Delta\rho/\Delta\tau$ は $R_\rho$ と一般に一致しないので、これは**半離散式に基づく推定**。
とし、引く前の値を `iface_q_eff_raw` に残す。$R_\rho$ は壁ピンの直前に `ifaceRro` へ退避する。定常では両者一致
(1 次元スラブ `case/52.conjugate_slab/run_0006_ctrl_newbin`: 上壁 −81.51316、`q_compact` −81.51314)。
壁エネルギー残差の内訳は `ifaceRconv` (対流の直後)・`ifaceRpre` (粘性の直前) で対流・ソース・粘性に分けられ、
`FORGE_WI_FORCE_DIAG=1` の `wi_eheat`/`wi_ework` で内部面の熱伝導と粘性仕事に分けられる
(`case/53.c3x_vane_cht/tools/resid_split.py`)。$q_{eff}$ と壁面勾配流束の差 (C3X で +2.5 %) の平均は
壁半 CV 内の粘性仕事 $\tau_w U_1(1-f)$ で、第一層厚に比例する (2 µm で 2.47 %、1 µm で 1.31 %)。

**適用範囲** (超えたら `NaN` を出す。ここに無い構成では使わない):

- **定常 (`unsteady: 0`) のみ**。dual-time は $C_i=D_t(V_iE_i)-R_i^{raw}$ で式が違う。
- **node の等温壁 (`wall_isothermal`) で Dirichlet ピンがあるもののみ**。断熱壁には $R^{raw}$ が無い。
- **周期境界に属する壁ノードを除く**。root 単位の 1 回集計が要る (依存 plan の解除待ち)。

**検証** (V1, `case/52.conjugate_slab/run_0005_qeff`): 1 次元純伝導の解析解 $q$=81.3090 W/m² に対し
`iface_q_eff` **81.1837 W/m² (−0.154 %)** で、`q_compact` (−0.155 %)・`q_recon` (−0.155 %) と同等。
**符号と絶対値を再現している**。


#### 合否ゲート — **G-if (界面反復) と G-cons (収支) を別々に満たす**

**G-if** (`cht_loop`、判定は `solid_shell.FixedPointDriver.advance`): 次を**独立に**満たし、
かつ `--n-consec` 回連続したときだけ収束とする。1 つでも欠けると、$D_f$ を上げて更新が
小さくなっただけの状態を収束と認めてしまう。

| 量 | 意味 | 渡し方 |
| --- | --- | --- |
| `dTw` | 温度更新の絶対値 [K] | `--tol-K` |
| `res_abs` | 界面残差の**絶対値** [W] (単位奥行きなら W/m) | `--tol-abs-W` (**事前登録**) |
| `res_rel` | $\max\lvert r\rvert/\max\lvert Q_f\rvert$。**規格化は $Q_f$ のみ** | `--tol-rel` |
| `res_solid` | 固体**内部**の残差 (Schur 縮約と内部復元の整合) | `--tol-solid` |
| 退避していないこと | $D_f$ を上げた反復は収束と認めない | — |

**規格化に $b$ を混ぜない**のが要点。旧実装は $\max(\lvert Q_f\rvert,\lvert b\rvert)$ で割っており、
背面温度で $b$ が大きいと $A_s=1000$ W/K, $b=3\times10^5$ W, $Q_f=1$ W, $D=10^9$ W/K のような構成で
`res_rel` が **3.33e−6** に見え、物理的な不釣合いが **100 %** でも合格した
(回帰試験 `test_solid_shell.py` T7)。

**G-cons** (`tools/check_cht_balance.py`): 同一状態で

$$\varepsilon=\Big|\sum_i Q_{f,i}-Q_{\rm solid}\Big|,\qquad
  \text{分母}=\max\Big(\sum_i\lvert Q_{f,i}\rvert,\ Q_{\rm floor}\Big)$$

とし、$\varepsilon/\text{分母}\le$ `--tol-rel` (既定 0.5 %) **かつ** $\varepsilon\le$ `--tol-abs` で合格。
**正味量で割らない** (正負が相殺する構成で分母が消える)。$Q_{\rm floor}$ は
**ケースごとに計算前に登録**する絶対床。`iface_q_eff` が `NaN` の節点が 1 つでもあれば**不合格**
(適用範囲外の構成をそのまま通さない)。

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
