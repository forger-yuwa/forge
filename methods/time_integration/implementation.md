# 時間積分 — 実装

forge の時間積分・更新カーネルの実装とソース対応をまとめる。
理論的背景は [theory.md](theory.md) を参照。

## ソースファイル

| ファイル | 役割 |
| --- | --- |
| [`solver_density_cuda/cuda_forge/timeIntegration_d.cuh`](../../solver_density_cuda/cuda_forge/timeIntegration_d.cuh) | カーネル / ラッパ宣言 |
| [`solver_density_cuda/cuda_forge/timeIntegration_d.cu`](../../solver_density_cuda/cuda_forge/timeIntegration_d.cu) | RK / 陰解法カーネル本体 (約 950 行) |
| [`solver_density_cuda/cuda_forge/update_d.cu`](../../solver_density_cuda/cuda_forge/update_d.cu) | 外側・内側ループ前後の Q コピー、陰解法補正反映 |
| [`solver_density_cuda/cuda_forge/implicitCorrection_d.cu`](../../solver_density_cuda/cuda_forge/implicitCorrection_d.cu) | dual-time 陽 (簡易) スキームの補正 |
| [`solver_density_cuda/cuda_forge/setDT_d.cu`](../../solver_density_cuda/cuda_forge/setDT_d.cu) | 局所 $\Delta t_{\text{loc}}$ 計算 |
| [`solver_density_cuda/update.cpp`](../../solver_density_cuda/update.cpp) | 旧 CPU 経路 (現状未使用) |

## エントリポイント

[`timeIntegration_d_wrapper`](../../solver_density_cuda/cuda_forge/timeIntegration_d.cu#L780)
が共通入口。`cfg.timeIntegration` の値に応じて次を呼ぶ。

| 値 | 呼び出し |
| --- | --- |
| `4` | `runge_kutta_exp_4th_d` |
| `1`, `3` | `runge_kutta_exp_d` |
| `11`, `blockDPLUR == 1` | `implicit_defect_correction_block_d`（5×5 block DPLUR、LU-SGS $A^\pm$。**推奨既定**） |
| `11`, `blockDPLUR == 0` | `implicit_defect_correction_d`（scalar 対角版＝スペクトル半径。軽量だが擬似 CFL が低く収束が遅い） |

> `solverConfig::initTimeIntegrationScheme` の `case 11` は `blockDPLUR ∈ {0,1}` を受理（それ以外を throw）。両者とも古典 DPLUR 制御フロー（残差固定 + `nStepInner` sweep + 単一 commit）に対応。`unsteady == 1`（dual-time）は dispatcher 側で throw する（本体未実装）。

> **軸対称ソース項の Jacobian**: scalar 版 `implicit_defect_correction_d` も block と整合させ、軸対称フープ源 `res_roUy += (P − τ_θθ)·A_planar` の Jacobian 対角成分 `A_pl·((γ−1)u_y + 2μ/(ρ r_eff))`（非負側、per-cell γ で TP 整合）を roUy 方程式の対角に陰化する。ただし scalar は式間連成 (源 Jacobian の非対角) を表現できないため、軸近傍以外が律速の強膨張ケース（例 `case/29` 出口コーナーの 2 次 MUSCL オーバーシュート）は救えない。**2 次精度の陰解法は block DPLUR 推奨**、scalar DPLUR は 1 次（起動・ロバスト用）に限るのが実用指針（切り分けは [`time_integration-scalar-dplur-axisym-source.md`](../../plans/accepted/time_integration-scalar-dplur-axisym-source.md) / `case/29.bell_vs_conical/README.md`）。

## ループ全体 ([`main.cpp`](../../solver_density_cuda/main.cpp))

`advanceOneStep` は巨大 lambda を廃し、`StepContext`（cfg, cuda_cfg, msh, mat_ns, var, fluct,
pprobes, profiler, residual_logger, implicit_diag_logger, iStep を束ねた参照集約構造体）を
受け取る自由関数群に分解する。スキームは 3 階層構造で、dual-time が後付けできる形にする。

```
advanceOneStep(ctx):                       // dispatcher
  if (isImplicit):
     if (unsteady) advanceImplicitDualTime(ctx)   // 後続フェーズ: 今は throw
     else          advanceImplicitSteady(ctx)
  else            advanceExplicitRK(ctx)

assembleResidual(ctx, stage):              // 残差組み立ての単一情報源（旧 assembleCurrentState）
  updateInner → dependentVariables → gasProperties → applyBconds → applyRansScalarBoundaries
  → calcGradient → axisymmetricGeomTerms → limiter → ducrosSensor → turbulent_viscosity
  → convectiveFlux → ransTransport → ransGradient + ransSource
  → axisymmetricSource → viscousFlux
  → [addUnsteadyTimeTerm(ctx)]            // dual-time の BDF 物理時間項フック（定常は no-op）

advanceExplicitRK(ctx):                    // tI 1/3/4（挙動不変）
  for iloop in perStepIterationCount():
     updateVariablesInner; assembleResidual(ctx,iloop+1); logResidualSnapshot
     timeIntegration_d_wrapper(iloop); ransTimeIntegration_d_wrapper(iloop)
  updateVariablesOuter; writeStepOutputs; setDT; logOuterEnd

implicitNonlinearUpdate(ctx):              // 定常・dual-time 共有の核
  assembleResidual(ctx, 1)                 // 残差・フラックスは 1 回（ransSource が src_jac_k/ω も出力）
  setDT_d_wrapper                          // 局所擬似時間 dτ（diag の V/Δτ）
  blockDPLURSolve(ctx)                     // 下記（平均流 5 式）
  applyBlockImplicitCorrection(ctx)        // Q = Q_baseline + dq を 1 回 commit
  if scalarResidualEnabled:                // RANS(SST) のとき
     applySSTPointImplicit(ctx)            // k/ω を segregated point-implicit で更新（凍結解除）

blockDPLURSolve(ctx):                      // 古典 DPLUR 線形ソルバ（res・Q 固定）
  for iSweep in nStepInner:
     implicit_defect_correction_block_d(...)   // 固定 res + lagged dq_old → dq_new
     swapBlockImplicitCorrectionBuffers(var)

advanceImplicitSteady(ctx):                // 定常: 擬似時間=メインループ、1 更新/step の縮退形
  updateVariablesInner; logOuterBegin
  implicitNonlinearUpdate(ctx); logResidualSnapshot
  updateVariablesOuter; writeStepOutputs; setDT; logOuterEnd
```

`updateVariablesOuter_d` は外側ループ開始時に `Q_N`, `Q_M` の両方を現在値に揃え、
`updateVariablesInner_d` は内側ステージ後の `Q_M` のみを更新する。

## RK カーネル詳細

### `runge_kutta_exp_d` (Jameson 多段)

ステージ係数 `coef_N, coef_M, coef_Res` を内部ループ index `loop` で参照し、

```cpp
Q[ic] = coef_N * Q_N[ic] + coef_M * Q_M[ic]
      + coef_Res * res[ic] * dt_local[ic] / vol[ic];
```

を全成分に適用。`dt_local` は `setDT_d_wrapper` で事前に書き込まれている。

### `runge_kutta_exp_4th_d`

低 storage 4 段 4 次 RK。`loop == 0` で残差累積バッファ `res_*_m` をゼロクリアし、
各段で `res_*_m += coef_Res * res * dt_local / vol`。
`loop < 3` の中間段では `Q = Q_N + coef_DT * res * dt_local / vol`、
最終段 (`loop == 3`) で `Q = Q_N + res_*_m` として確定。

### `runge_kutta_exp_scalar_d`（スカラー k/ω の陽解法 RK ＋ point-implicit 源項）

[`scalarTransport_d.cu`](../../solver_density_cuda/cuda_forge/scalarTransport_d.cu) の
`ransTimeIntegration_d_wrapper`（[`ransTransport_d.cu`](../../solver_density_cuda/cuda_forge/ransTransport_d.cu)）が平均流 RK と同じ段で k/ω を別カーネルで積分する
（`timeIntegration==1/3` は `runge_kutta_exp_scalar_d`、`==4` は `runge_kutta_exp_scalar_4th_d`）。

RANS (SST) の消散項・輸送項は stiff なため、`timeIntegration==1/3` の更新は**残差増分のみ源項+輸送ヤコビアンで減衰**する:

```cpp
const flow_float fac = 1.0 + coef_Res * dt_l * (src_jac[ic] + transport_diag[ic] / v); // ≥ 1
rho_phi[ic] = coef_N * rho_phi_N[ic] + coef_M * rho_phi_M[ic]
            + (coef_Res * res_rho_phi[ic] * dt_l / v) / fac;
```

- `src_jac`（消散 $\beta^\*\omega,\,2\beta\omega$）は `ScalarTransportDesc` 経由で k=`src_jac_k`、ω=`src_jac_omega`
  （[`ransSource_d.cu`](../../solver_density_cuda/cuda_forge/ransSource_d.cu) が毎 `assembleResidual` で出力）。
- `transport_diag`（移流+拡散の対角 $\Lambda^{T}_\phi$ [m³/s]）は `scalar_advection`/`scalar_diffusion` カーネルが
  面ループで集計（k=`transport_diag_k`、ω=`transport_diag_omega`）。$V$ で割って $[1/s]$ 化して `fac` に入る。
- 非 RANS (`model` が `sst*` 以外: 乱流なし/LES) では両者 0 で `fac=1`、従来の純陽的更新と一致（LES は無影響）。
- 減衰係数は `applySSTPointImplicit_d` の対角 $D_\phi=V/\Delta\tau+V\cdot\text{src\_jac}+\text{transport\_diag}$ に $\Delta\tau/V$ を掛けた形と整合。
- `runge_kutta_exp_scalar_4th_d`（4 次）は本減衰未適用＝RANS 非対応のまま（平均流 4 次自体が RANS 未検証）。
- 理論は [theory.md](theory.md) §"陽解法 RK での point-implicit 源項"。

## 陰解法カーネル詳細

### `implicit_defect_correction_block_d`（block DPLUR、LU-SGS $A^\pm$ 分割）

5×5 ブロック版。`block_dplur` 名前空間に補助 device 関数群を持つ。カーネルは線形 solve の内部精度
`ST` (float/double) で `template<typename ST>` 化され、状態/残差 (float) を `ST` へキャストして取り込む
(混合精度。下記「閉形式 FVS と混合精度」参照)。

- `accumulate_split_jacobian_cf<T>`（**既定の閉形式 FVS**）— $R,\Lambda,L$ の 5×5×5 三重積を作らず、
  **音響右/左固有ベクトルのみ**で $M(g)=g_2 I+(g_1-g_2)\,r_1\!\otimes l_1+(g_5-g_2)\,r_5\!\otimes l_5$
  ($=R\,\mathrm{diag}(g)\,L$) を構成し、$a^{+}=M(\Lambda^{+})$ を対角へ・$k_{\rm off}=M((-\Lambda)^{+})$ を近傍へ
  **直接畳み込む**（$a^{+}/k_{\rm off}$ を materialize しない）。`build_jacobian_split` (下記) と軸対称 (nz=0) で
  数値等価かつ ~10% 高速 (レジスタ・スピル削減)。
- `build_jacobian_split` — （旧経路・precond 版が使用）固有分解 $R,\Lambda,L$ から **対角用 $A^{+}=R\,\Lambda^{+}L$** と
  **RHS 近傍用 $K=-A^{-}=R\,(-\Lambda^{-})L=\tfrac12(|\widetilde A|-\widetilde A)$** を同時に返す
  （$\Lambda^{+}=\max(\Lambda,0)$、$-\Lambda^{-}=\max(-\Lambda,0)$、共に非負）。
- `add_identity_scaled`, `add_scaled_5x5` — $V/\Delta\tau\,I$・粘性対角・面寄与の加算。
- `solve_5x5` — 部分ピボット付き Gauss 消去（`diag` を破壊して in-place）。`|pivot| < 1e-20` でゼロ解にフォールバック。
- `multiply_add_5x5_vec` — 行列ベクトル積。

各セルで近傍寄与を集約し、対角 $D_i = V/\Delta\tau\,I + \sum_f A^{+}_f S_f + \sum_f \Lambda^{\nu}_f\,I$、
RHS = $-\mathbf R + \sum_f K_f S_f \cdot \Delta\mathbf Q_{\text{nbr}}^{\text{old}}$（$K_f=-A^{-}_f$）を構築し、
$\Delta\mathbf Q_{\text{new}} = D_i^{-1}\,\text{RHS}$ を解く。`cfg.implicitRelax` で $\Delta\mathbf Q$ を緩和。

ここで粘性対角は $\Lambda^{\nu}_f = 2\nu_f\,\dfrac{|S_f|^2}{\Delta\mathbf{cc}_f\cdot S_f}$（$\nu_f=(\mu_{\rm lam}+\mu_t)/\rho$）。
これは粘性流束 residual ([`viscousFlux_d.cu`](../../solver_density_cuda/cuda_forge/viscousFlux_d.cu)) の
法線拡散項 $\mu_f(\Delta U/|\Delta\mathbf{cc}|)\,\delta$（$\delta=|\Delta\mathbf{cc}|\,|S_f|^2/(\Delta\mathbf{cc}_f\cdot S_f)$）の
**Jacobian の大きさと整合**する（$\Lambda^{\nu}_f = 2\nu_f\,\delta/|\Delta\mathbf{cc}|$）。

> **2026-06-14 修正（粘性対角の幾何是正）**: 旧コードは粘性対角を $|S_f|\cdot(2\nu_f/\delta)$ と書いていたが、
> $\delta$ は**面積次元** ($\approx|S_f|$) なので $|S_f|$ が約分されて $\approx 2\nu_f$ に潰れ、(1) 軸対称近軸で
> 本来 $\propto r$ で消えるべき内側面の寄与を過大評価し、(2) residual に無いゼロ面積(対称/軸)面にも
> スプリアス項を載せていた（residual 側は `ip<nNormalPlanes` で境界面を除外）。これが **float block-DPLUR が
> 軸対称近軸第一セルの $U_r$ を収束させきれず固着する真因**だった。上記の residual 整合形
> $2\nu_f\,|S_f|^2/(\Delta\mathbf{cc}_f\cdot S_f)$ に是正すると $\propto r$ でゼロ面積面では消え、**float のまま固着が解消**
> （case 29 laminar conical 第一セル $U_r$: $+1.4\to+17.9$ で double solve と一致、1 次では未収束→収束）。
> 修正は `timeIntegration_d.cu` の scalar (`implicit_defect_correction_d`) / block (`implicit_defect_correction_block_d`) /
> precond (`implicit_defect_correction_block_precond_d`) の 3 箇所。**LHS のみの変更**で defect-correction の
> 定常解は不変（planar 回帰 bump で base/fix 場が $L2\sim10^{-5}$ 一致・RANS で残差レベル同一を確認）。

> **2026-06 修正**: 旧コードは対角に $A^{+}$ ではなく $|\widetilde A|$ を、近傍に $-A^{-}$ ではなく $+|\widetilde A|$ を
> 使っていた（符号付き分割でなく絶対値の誤用）。対角が upwind 自己 Jacobian と不一致・近傍結合が逆符号となり、
> block DPLUR は収束せず発散していた。`build_jacobian_split` による $A^{+}/{-}A^{-}$ 分割でこれを修正。
> `res_*` の符号は $-\mathbf R$（陽解法 `runge_kutta_exp_d` の `Q=Q_N+\mathrm{res}\,\Delta t/V` と整合）なので、
> カーネルでは `rhs` を `res_*` で初期化し近傍寄与 $K_f S_f\,\Delta\mathbf Q_{\text{nbr}}$ を**加算**する。

**古典 DPLUR の構造（重要）**: カーネルは `dq_block_new` の生成のみを行い、
`Q`（ro..roe）を**インライン更新しない**。`blockDPLURSolve` が `nStepInner` 回 sweep を回し、
各 sweep 後にドライバ側で [`swapBlockImplicitCorrectionBuffers`](../../solver_density_cuda/cuda_forge/timeIntegration_d.cuh)
（block 専用 swap、`.cuh` に公開）で `dq_block_old <-> dq_block_new` を入れ替える。
全 sweep 後に [`applyBlockImplicitCorrection`](../../solver_density_cuda/cuda_forge/update_d.cu)
（`applyScalarImplicitCorrection` と対称、`update_d.cu` に新設）で `Q = Q_baseline + dq_block` を
**1 度だけ** commit する。残差 `res_*` と $|\widetilde A_f|$ は sweep 中固定（matrix-free のため固定 Q から毎 sweep 再構築してよい）。

#### 閉形式 FVS と混合精度 (`implicitSolvePrecision`)

`accumulate_split_jacobian_cf<T>` は固有ベクトル行列 $R,L$ を陽に作らず、$\mathrm{diag}(g)-g_2 I$ が
shear/entropy の 3 モードを消すことを使って **音響右/左固有ベクトル $r_1,r_5,l_1,l_5$ のみ**で
$M(g)=R\,\mathrm{diag}(g_1,g_2,g_2,g_2,g_5)\,L = g_2 I+(g_1-g_2)\,r_1 l_1^\top+(g_5-g_2)\,r_5 l_5^\top$
を構成し、$a^{+}=M(\Lambda^{+})$ を対角へ・$k_{\rm off}=M((-\Lambda)^{+})$ を近傍へ直接畳み込む。
$R\,\Lambda\,L$ の 5×5×5 三重積と $a^{+}/k_{\rm off}/\text{solve\_mat}$ 保持を排除し、float 陰解法で ~10% 高速
(レジスタ・スピル削減)。**軸対称・平面 (nz=0) では legacy `build_jacobian_split` と数値厳密一致**
(一般 3D は legacy の $R,L$ が厳密逆行列でない=$RL\neq I$ ため僅差だが、固有値を厳密に $\max(\lambda,0)$ とする
valid な FVS で defect-correction の定常解は不変)。

カーネルは線形 solve の内部精度 `ST` で `template<typename ST>` 化。`solverConfig` の
**`implicitSolvePrecision`** (`time.deltaT`、既定 `0`) で切替える:

- `0` (float): 状態/残差 (float) のまま float で組立・solve（既定・高速）。
- `1` (double): 状態/残差を **double へキャスト**して Jacobian 構築・5×5 solve・近傍 sweep を **double** で行い、
  補正 $\Delta\mathbf Q$ を float `dq_new` へ書戻す（混合精度 iterative refinement）。

**動機と位置づけ（重要・2026-06-14 更新）**: 当初 float32 の block-DPLUR が軸対称 近軸第一セルで平均速度 `Uy` を
収束させきれず偽固着する (laminar conical で `Uy` が物理値 $-15$ でなく $-0.6$) 問題に対し、`implicitSolvePrecision=1`
(線形 solve を double) を root-fix とした。**その後、真因は精度ではなく上記「粘性対角の幾何不整合」であることが判明**
(粘性が無い Euler は float でも固着せず、固着は粘性 LHS 由来。詳細は本節冒頭「粘性対角の幾何是正」)。
粘性対角を residual 整合形に直せば **float のまま固着が解消**し double solve は不要。
したがって `implicitSolvePrecision=1` は**根治ではなく、悪条件 LHS を倍精度で押し切る検証/保険用の手段**として残す
(幾何是正後は通常 `0` で良い)。RTX 3060 では FP64=FP32 の 1/32 ゆえ ~×2.8 遅い。
切り分け・速度の詳細は [`.github/plans/precision-mixed-axisym.md`](../../plans/archived/precision-mixed-axisym.md)
と [`.github/plans/architecture-axisym-axis-singularity.md`](../../plans/accepted/architecture-axisym-axis-singularity.md)。
現状 `blockDPLUR=1`・`lowMachPrecond` 0/1 経路のみ対応 (precond>=2 / scalar 版は float のまま)。

**低マッハ前処理 (LHS 固有値) は不採用** (2026-06 検証・実装後 revert)。`build_jacobian_split` の
固有値 `lambda[5]={U+c,U,U,U,U-c}` を前処理固有値 `U±c'` に差し替える案を実装・検証したが、
block DPLUR では**対角優位性の源である大きい音響固有値を縮めてしまい有害**（フラックス散逸前処理
単独で安定だった `eps=0.15` すら発散させ、安定 `eps` 範囲を狭めた。収束加速も根治もなし）。よって
LHS は従来の $A^\pm$（物理音速 `sonic`）のまま。低マッハ前処理は `SLAU_d` の散逸スケールにのみ適用する。
根拠・データは [`theory.md`](theory.md#低マッハ前処理固有値-weisssmith--試行したが不採用) と計画
[`time_integration-lowmach-preconditioning.md`](../../plans/accepted/time_integration-lowmach-preconditioning.md) §9。

#### 一般EOS固有系 (TP gas, `thermalMethod==2`) — 実装・検証完了 (閉形式)

閉形式 `accumulate_split_jacobian_cf` は `inv_chi=sonic/(γ-1)` で全エンタルピーを $K+c^2/(\gamma-1)$ と再構成し、
接触波エネルギー成分に $K$ を使う。これは **CPG 専用**で TP では真の $\partial\mathbf F/\partial\mathbf Q$ と一致しない
(理論は [theory.md](theory.md) 「一般EOS固有系」、FD 検証で TP 誤差 469%)。修正方針 (plan
[`time_integration-general-eos-jacobian.md`](../../plans/accepted/time_integration-general-eos-jacobian.md)):

- **分岐保持**: `thermalMethod==0` は現行閉形式のまま (ビット不変・回帰基準)。`==2` のみ一般EOS固有系へ。
  数値 LU は閉形式と演算順序が違うため CPG ビット一致は望めない → CPG 経路を残すのが回帰構成。
- **thermo helper**: `thermo_d.cuh` に `ThermoDerivatives{p,h,cp,cv,R,kappa,chi,a2}` を 1 セル状態から返す
  `thermo_derivatives_mix` を追加。$c^2=\gamma RT,\ \kappa=\gamma-1,\ h=e+RT$ を**同一 $T$・同一組成・同一 NASA
  data・同一 datum**で評価 (`sonic[ic]`/`Ht[ic]`/`gamma[ic]` の別経路混在を避ける)。$\chi=c^2-\kappa h$。
- **固有系+LU**: 一般EOS右固有ベクトル ($\mathbf r_\mp,\mathbf r_c,\mathbf r_{sk}$, 実 $H_t$・接触 $H_t-c^2/\kappa$) を
  構築し、**部分ピボット付き double 5×5 LU** で $L=R^{-1}$ を解く ($R$ 構築・LU・$R\Lambda L$ 積算は試作段階で
  double; $H_t-c^2/\kappa=K-\chi/\kappa$ が大きな二数の差で float だと条件数誤差)。$A^\pm=R\Lambda^\pm L$。
- **セル/面の不混在**: DPLUR の対角・非対角ブロックは同一セル $i$ の $(\rho_i,\mathbf u_i,H_{t,i},c_i,\kappa_i,\chi_i)$ で
  構築 (RHS 数値流束が面/Roe 平均でも可、LHS は近似ヤコビアン)。$H_t$ だけ面・$c,\kappa$ セルの混成は不可。
- **検証 (デバイス単体 → ノズル)**: Level1 `‖LR−I‖`/`‖RΛL−A_FD‖` (CPG/TP 250・1000・高温/M=0/亜音速/音速近傍/
  超音速/一般方向)、Level2 split `A⁺+A⁻=A`・法線反転、Level3 DPLUR 行列作用、Level4 ノズル cfl 2/5/20/50/100
  (発散 step・収束率・残差床・出口M・massflux・全エンタルピー・Newton 反転回数)。**残差床は別トラック**で EOS
  反転誤差 ($\max_i|e(T_i)-e_{\text{target},i}|/\max(|e_{\text{target}}|,e_{\rm ref})$) と切り分け、機械ゼロは保証しない。
- **実装 (閉形式・確定)**: `accumulate_split_jacobian_cf` (共有ヘッダ `block_dplur_jacobian_d.cuh`) の音響右固有ベクトル
  エネルギーを実 `Ht`、左密度成分に `χ_eos=c²−κh` を加える 3 項改変 (`thermallyPerfect` 分岐、CPG ビット不変)。
  数値 L=R⁻¹ は検証参照のみ。**結果**: TP cfl 上限 2→≥100・残差 9.6e-8→4e-11、Level1/2/3 PASS (`tools/test_eos_jacobian.cpp`)。
- **precond=2 経路 (一般EOS統一・CPG 回帰確認済・TP end-to-end 未検証)**: `implicit_defect_correction_block_precond_d` は
  TP で 3 箇所 CPG 仮定 (`build_jacobian_split` の R[4]・Γ_c の g ベクトル `Htot`・∂p/∂Q の `rvec[0]` の χ 欠落) だった。
  **`build_jacobian_split` を検証済み `eos_split_jacobian_general_closed` を呼ぶ薄いラッパに統一** (CPG 専用べた打ち R/L=近似を廃止)、
  Γ_c の g/r も実 Ht・χ_eos=c²−κh に統一し `thermallyPerfect` 分岐を撤去 (CPG は χ_eos≈0 で簡約)。precond=2 が標準経路と同じ
  検証済み固有系を共有。**CPG precond=2 回帰確認** (`case/23` inviscid precond=2 eps0.15: step4000 rms_ro 1.60e-5 vs 旧 1.38e-5,
  NaN 無, 同一収束域)。ただし **TP precond=2 の収束を正で示す低マッハ TP ケースが無く** (超音速 wys は CPG でも precond=2 発散)、
  TP の end-to-end 検証は未。**TP は標準経路 (precond=0/1) が検証済・推奨**。

### `implicit_defect_correction_block_precond_d`（Phase 4: 完全 $\Gamma^{-1}A$ 前処理・`lowMachPrecond>=2`）

上の LHS 固有値前処理 (不採用) と異なり、**前処理を一貫させた別カーネル**。`blockDPLUR==1 && lowMachPrecond>=2`
のとき wrapper がこちらを起動する (既存 `implicit_defect_correction_block_d` は 0/1 専用・ビット/レジスタ不変)。

- **`lowMachPrecond=2`**: RHS のフラックス散逸も `c'` に是正 (`SLAU_d`) しつつ本 LHS 前処理を併用 → 低マッハ域の
  **収束解そのものを変える** (自励振動フロアの根治)。
- **`lowMachPrecond=3` (LHS-only)**: 本 LHS 前処理だけを行い RHS 散逸は `c_hat` のまま (`SLAU_d` 側で `==1||==2`
  のみ `c'` を使う)。本カーネルは保存形 $(\Gamma_c V/\Delta\tau' + A_c)\Delta Q=-R$ で**前処理は時間項のみ・$A_c$ と
  $R$ は非前処理**なので、**収束解は `lowMachPrecond=0` とビット一致 (保存性・解ともに不変)**。純粋に低マッハ
  剛性を除去して収束を加速する LHS 操作として使う。`setDT` の $\Delta\tau'$ 拡大も `>=2` で両モード共通に効く。

- **Sherman-Morrison 解法 (FP64 回避・高速)**: $\Gamma_c=I+\alpha g r^\top$ がランク 1 ゆえ対角ブロックは
  $D=D_0+\gamma g r^\top$ ($D_0$=V/Δτ'·I+物理 FVS+粘性+軸対称=既存 block と同形・良条件、$\gamma=(V/\Delta\tau')\alpha$)。
  $D_0$ を **float** で 2 RHS 同時 ([`solve_5x5_2rhs`](../../solver_density_cuda/cuda_forge/timeIntegration_d.cu)) に解き
  $y=D_0^{-1}b,\ z=D_0^{-1}g$、$x=y-[\gamma(r^\top y)/(1+\gamma(r^\top z))]z$。悪条件 $\sim1/\beta$ は分母スカラーのみ double。
  consumer GPU (FP64=FP32/64) でも物理 block とほぼ同速 (per-step 16.9 vs 17.8ms)。$g,r,\alpha$ はカーネル内インライン。
- **時間項**: $\Gamma_c\,V/\Delta\tau'$ (上の $\gamma g r^\top$ 寄与)。$\Delta\tau'$ は `setDT_d` の
  `setDTlocal_precond_scale_d` が `dt_local *= (|u|+c)/ρ'` で拡大済 ($\rho'$=前処理スペクトル半径)。
- **フラックス分割**: **物理の厳密 FVS** をそのまま使う (対角に $a^{+}=A_c^{+}$、近傍に $k_{\rm off}=-A_c^{-}$。既存 block と同一)。
  保存形 $(\Gamma_c V/\Delta\tau' + A_c)\Delta Q=-R$ より前処理は**時間項のみ**で、$A_c$ は残差の真のヤコビアンゆえ非前処理が正しい
  ([`theory.md`](theory.md) 参照)。収束は $(\Delta\tau'/V)\Gamma_c^{-1}A_c$ の固有値 $\lambda'$ で一様に前処理される。
  ※当初フラックスも $\hat A_\Gamma=\Gamma_c^{-1}A_c$ のスペクトル半径分割で前処理したが、別系で過散逸ゆえ撤回 (上記が正)。
- **その他**: 粘性スペクトル半径・軸対称ソースヤコビアン・dual-time 物理 BDF 項 (非前処理) も倍精度で踏襲。
- **`SLAU_d`** は `lowMachPrecond==1||==2` で `c'` 散逸を使う (==2 は散逸是正を併用、==3 は RHS 散逸を触らず `c_hat`)。
  $\beta=1$ で $\Gamma_c=I$・$\Delta\tau'=\Delta\tau$・フラックス同一ゆえ現行カーネルと同一組み立て (倍精度の丸め差 ~1e-7、解一致)。

> **検証結果 (2026-06-09・採用)**: `case/23.axi_nozzle` で**低マッハ自励振動を根治**。Phase1 (物理 LHS) が発散した
> $\epsilon=0.05$ を前処理 LHS が安定化し、chamber 圧振幅 (M<0.08, 4k–20k) を 0.882%→**0.087% (定常収束・振動消滅)**。
> 安定 `cfl_pseudo` も m1~1→m2~5-7 と拡大 (収束加速は per-step 2.54× で等 wall-clock 互角。価値は根治)。データは計画 §9。

### `implicit_defect_correction_d`（scalar 対角版＝スペクトル半径）

平均流 5 式を、5×5 ブロックの代わりに**スカラー対角**で解く軽量版（`blockDPLUR == 0`）。
各セルで対角 $D = V/\Delta\tau + \sum_f(|U_n|+c+\rho^\nu)S_f$、off-diagonal $\tfrac12(|U_n|+c+\rho^\nu)S_f$ で
`dq_*_new = relax·(res_* + Σ offdiag·dq_*_nbr^{old}) / D`（5 変数とも同じスカラー $D$）を作り、
[`blockDPLURSolve`](../../solver_density_cuda/main.cpp) が `nStepInner` 回 sweep（各 sweep 後に
`swapScalarImplicitCorrectionBuffers`）、最後に [`applyScalarImplicitCorrection`](../../solver_density_cuda/update.cpp)
で `Q = Q_N + dq_*_old` を 1 回 commit する（block 版と同じ古典 DPLUR 制御フロー）。
スペクトル半径は符号不変なので block の $A^\pm$ 法線符号問題は持たない。

**block 版との比較（2026-06, `case/20.naca_ml`）**: 収束先は同一（収束場から再開すると同じ解を保持、
壁面静圧 平均 0.02% 一致）。ただし近似ヤコビアンが粗いため**安定 `cfl_pseudo` が大幅に低い**
（scalar ≲ 1〜2、block は 20〜50）。supercritical 始動では block が cfl_pseudo=20 で 4000 step / 25s で
roe→0.5 に収束する一方、scalar は cfl_pseudo=1 で 12000 step でも収束せず大きな過渡オーバーシュート
（roe ピーク ~186 vs block/explicit ~85）を示す。よって **block DPLUR が既定、scalar 対角版は
5×5 を避けたい軽量用途・低レジスタ用途向けのフォールバック**と位置づける。

### `applySSTPointImplicit_d`（SST k-ω の segregated point-implicit）

平均流 commit 後に呼ぶスカラー陰解法。源項+輸送項のヤコビアン対角を陰化して k/ω の凍結を解き、
壁近傍の陽的輸送 stiff 性を緩和して安定 `cfl_pseudo` を一桁以上引き上げる。

- 消散ヤコビアン `src_jac_k=β*ω`・`src_jac_omega=2βω` は
  [`ransSource_d.cu`](../../solver_density_cuda/cuda_forge/ransSource_d.cu) `rans_sst_source_d` が出力。
- 輸送ヤコビアン `transport_diag_k`/`transport_diag_omega` [m³/s] は
  [`scalarTransport_d.cu`](../../solver_density_cuda/cuda_forge/scalarTransport_d.cu) の `scalar_advection_first_order_d`
  （1次風上 $\sum_f\max(\pm\dot m,0)/\rho$）と `scalar_diffusion_first_order_d`（$\sum_f(\mu_{\text{face}}/\rho)|\delta|/dcc$）が
  面ループで atomicAdd 集計する（`ransTransport_d_wrapper` 冒頭で毎 `assembleResidual` ゼロ初期化）。
- [`update_d.cu`](../../solver_density_cuda/cuda_forge/update_d.cu) の `applySSTPointImplicit_d` が各セルで
  $D_\phi = V/\Delta\tau + V\cdot\text{src\_jac}_\phi + \text{transport\_diag}_\phi$、
  $\delta(\rho\phi)=\text{relax}\cdot\text{res}_{\rho\phi}/D_\phi$、$\rho\phi=\max(\rho\phi^{N}+\delta,\ \text{floor})$ を適用
  （$\rho k\ge0$, $\rho\omega>0$）。`dt_local` は平均流と共用。
- [`main.cpp`](../../solver_density_cuda/main.cpp) `implicitNonlinearUpdate` で `scalarResidualEnabled` のとき
  `applyBlockImplicitCorrection` 直後に `applySSTPointImplicit` を呼ぶ。生産・近傍 ΔQ は `res` に含む lagged。
- 近傍 ΔQ 結合なしの純 point-implicit（消散+輸送の対角のみ陰化）。defect-correction のため定常解は不変。
  理論・検証は [theory.md](theory.md) §"輸送項 (移流+拡散) の point-implicit 対角"。

## 局所時間刻み

`setDT_d_wrapper` ([`setDT_d.cu`](../../solver_density_cuda/cuda_forge/setDT_d.cu)) が
セル中心スペクトル半径を集計して `dt_local[ic] = CFL * V / λ_max` を書き込む。
`cfg.dt`, `cfg.cfl` で挙動を制御。

**setDT の低マッハ前処理は不採用** (2026-06 検証)。`setCFL_pln_d` のスペクトル半径の音速 `sonic` を
前処理音速 `c'` に置換する案は、低マッハ域で `dt_local` を増大させ陰解法対角 `V/Δτ` を縮め、block DPLUR の
対角優位性を崩して発散させた。よって `setCFL_pln_d` は従来の `sonic` のまま。低マッハ前処理は対流フラックス
([`convectiveFlux_d.cu`](../../solver_density_cuda/cuda_forge/convectiveFlux_d.cu) `SLAU_d`) の散逸スケールにのみ
適用する。詳細は計画 [`time_integration-lowmach-preconditioning.md`](../../plans/accepted/time_integration-lowmach-preconditioning.md) §9。

## 入出力

入力: 残差 `res_ro, res_roUx, …`、$\Delta t_{\text{loc}}$ (`dt_local`)、
過去ステージの `Q_N, Q_M`、対角ヤコビアン構築用の `Ux, Uy, Uz, sonic, Ht, vis_turb`、
ステージ係数 (`cfg.coef_*`)。

出力: 更新後の `Q = (ro, roUx, roUy, roUz, roe)`。
陰解法では補助バッファ `dq_*_old/new`, `dq_block_old/new_k`, `diag_block_*`, `rhs_block_k`。

## 非定常 dual-time 陰解法（実装済み 2026-06）

`advanceImplicitDualTime` ([`main.cpp`](../../solver_density_cuda/main.cpp)) が 1 物理ステップを担当する。
使用条件 `unsteady=1, dualTime=1, timeIntegration=11, blockDPLUR=1, time.deltaT.control=0`（それ以外は throw）。

1 物理ステップの流れ:
1. **時間レベルシフト** `shiftDualTimeLevels_d_wrapper`: `roNN←roN`, `roN←ro`（$\mathbf Q^{n-1}\!\leftarrow\!\mathbf Q^n$, $\mathbf Q^n\!\leftarrow$現在）。
2. **BDF 係数**: 初回ステップ or `bdfOrder==1` は BDF1 $(a,b,c)=(1,1,0)$、以降 BDF2 $(\tfrac32,2,\tfrac12)$。
   `cfg.unsteadyDiagCoef = a/\Delta t` を設定（陰解法カーネルが対角に $V\cdot$この係数を加える）。
3. **擬似時間サブ反復** `nSubIterDualTime` 回:
   - `assembleResidual` → `addUnsteadyTimeTerm_d_wrapper(a,b,c)` で `res_* -= (V/\Delta t)(a\mathbf Q - b\mathbf Q^n + c\mathbf Q^{n-1})`
     （`include_scalar` で k/ω も）。
   - `setDT`（擬似 $\Delta\tau$）→ `blockDPLURSolve`（対角に $V\,a/\Delta t$ 込み）。
   - **in-place commit** `applyBlockImplicitCorrectionInPlace_d_wrapper`（$\mathbf Q\mathrel{+}=\delta\mathbf Q$。`roN`=$\mathbf Q^n$ 固定のため）。
     RANS のとき `applySSTPointImplicit`（こちらも in-place）。
4. `cfg.totalTime += \Delta t`、`unsteadyDiagCoef=0` リセット、出力。

実装した CUDA（[`update_d.cu`](../../solver_density_cuda/cuda_forge/update_d.cu)）: `addUnsteadyTimeTerm_d`,
`applyBlockImplicitCorrectionInPlace_d`, `shiftDualTimeLevels_d`。block/scalar/SST 各カーネルに対角の物理時間係数
`unsteady_diag` (`cfg.unsteadyDiagCoef`) を追加。第 2 時間レベル `roNN`/`roKNN` 系は [`variables.hpp`](../../solver_density_cuda/variables.hpp) に登録済。

## 一様体積力と質量流量一定制御（bodyForce / bodyForceCtrl）

周期境界系（チャネル・周期丘）を駆動する空間一様体積力。理論的には周期分解
$p_\mathrm{total} = -\beta(t)\,x + p'$ の非周期線形成分 $\beta$ を運動量ソース項に移項したものであり、
「平均圧力勾配駆動」と数学的に同一（局所の圧力勾配変動は解かれる $p'$ が担う）。
実装は [`bodyForce_d.cu`](../../solver_density_cuda/cuda_forge/bodyForce_d.cu):
運動量に $f_i V$、エネルギーに $(\mathbf f\cdot\mathbf u)V$ を residual へ加算（仕事項を落とすとエネルギー収支が破れる）。
閉じた周期系では体積力の仕事で系が加熱し続けるため、等温壁で排熱して温度を定常化させる。

### 固定値駆動（`bodyForce: [fx, fy, fz]`）

チャネルのように目標 $u_\tau$ から $f_x=\rho u_\tau^2/\delta$ を先験的に決められる場合はこれで足りる。

### 質量流量一定制御（`bodyForceCtrl: 1`）

周期丘のように断面が $x$ で変わる系では目標がバルク速度（=質量流量）で与えられ、抗力（壁摩擦+圧力抗力）が
先験的に分からないため、$\beta(t)$ を毎物理ステップ調整する（Benocci & Pinelli 1990 の圧縮性版）。
制御量は**体積平均 streamwise 運動量密度** $\langle\rho u_x\rangle_V = \frac1V\int \rho u_x\,dV$
（質量保存の下で統計定常では任意 $x$ 断面の質量流束と等価。目標値はケース側で
$\langle\rho u_x\rangle_V^{\,t} = \rho_0 U_b A_\mathrm{crest} L_x / V$ と換算して与える）。

更新則は 2 時刻の運動量収支から抗力を推定する deadbeat 型:

$$
f_x^{\,n} = f_x^{\,n-1} + \gamma\,\frac{M_t - 2M^n + M^{n-1}}{V\,\Delta t}
$$

（$M=\int\rho u_x\,dV$、$M_t$ は目標、$\gamma$=`bodyForceCtrlRelax` 既定 1.0。
初回は $M^{n-1}:=M^n$ とし P 項のみで立ち上がる。）
導出: $M^{n+1}=M^n+\Delta t(f_x V - D)$ で $M^{n+1}=M_t$ を課し、抗力 $D$ を前ステップの実収支
$D^{n-1}=f_x^{n-1}V-(M^n-M^{n-1})/\Delta t$ で推定して代入したもの。

実装 ([`bodyForce_d.cu`](../../solver_density_cuda/cuda_forge/bodyForce_d.cu) の `bodyForceCtrlUpdate`):

- 呼び出しは `advanceOneStep` 冒頭（物理ステップ境界）で 1 回。dual-time ではサブ反復を通じて $f_x$ 固定。
- $M^n$ は `thrust::inner_product(roUx, volume)`（once-per-step の D2H 同期 1 スカラー、コスト無視可）。
- 制御は $x$ 成分のみ。`cfg.bodyForceX` を書き換える（`bodyForceY/Z` は固定値のまま）。
- 全 CV 体積 $V$ は初回に `thrust::reduce` でキャッシュ。
- 履歴を `bodyforce_history.csv`（step, time, fx, M, M_target）へ 10 step ごとに追記。
- **`unsteady: 1` 必須**（物理 $\Delta t$ を使うため。定常局所 dt では throw）。定常 RANS の段階起動では
  固定 `bodyForce` を使い、非定常へ引き継ぐ際に `bodyForceCtrl: 1` へ切り替える運用。
- 陰解法での扱いは固定値版と同じ明示ソース（ステップ内定数のため Jacobian 不要。$\mathbf f\cdot\mathbf u$
  仕事項の対角寄与は微小につき省略）。

## 陰的更新の正値性ガード (`updateGuardAlpha`, 2026-09-02)

計画: [`plans/active/time_integration-update-positivity-guard.md`](../../plans/active/time_integration-update-positivity-guard.md)。

block-DPLUR の commit (`applyBlockImplicitCorrection_d` / 同 InPlace) で、セルごとに
縮小率 $s\in\{1,1/2,\dots,1/32\}$ (半減列) を選び $q\leftarrow q_N+s\Delta q$ とする:

$$\rho(q_N+s\Delta q)\ge\alpha\rho(q_N),\qquad e_i(q_N+s\Delta q)\ge\alpha e_i(q_N),\qquad e_i=\rho e-\tfrac12|\rho\mathbf u|^2/\rho$$

$\alpha$ = `time.deltaT.updateGuardAlpha` (既定 0.0 = OFF・ビット同一迂回、推奨 0.5)。
P でなく $(\rho, e_i)$ を使うのは TP で EOS 反転を避けつつ P 崩落と等価な検知をするため
(CPG では $P=(\gamma-1)e_i$ で厳密に等価)。5 回半減しても満たせないセルは $s=1/32$ で
commit し、EOS 床は最後の防波堤として残す。既に $\rho\le0$/$e_i\le0$ のセルは対象外。
目的は「P アンダーシュート → pMin 床洗浄 → NaN」(case/45 run_0015_cflsweep で特定した
cfl_pseudo 上限の律速) を、全域 `implicitRelax` なしにセル局所で抑えること。

## line-implicit (`lineImplicit`, 2026-09-02)

計画: [`plans/active/time_integration-line-implicit.md`](../../plans/active/time_integration-line-implicit.md)。

壁法線ライン (積層方向の CV 鎖、壁 CV 種の greedy 構築) 上の隣接結合を、point-DPLUR の
lag から **block 三重対角の直接解 (block-Thomas, 1 ライン 1 スレッド・内部 double)** に
昇格する。sweep カーネルはライン CV について点解せず diag/rhs/近傍行列 K (単位ベクトル列
抽出 = 点経路と同一 Jacobian) を保存し、各 sweep 直後の `lineThomas_d` が dq_new を上書き
する。反復に残る lag はライン外 (流れ方向) のみになる。**ただし case/45 M6 ノズルの実測では cfl 上限は不変** (律速は壁法線でなく streamwise lag と判明) で、効果は同一設定の収束 −14 %/step に留まる — 詳細と罠 (lu5 ピボット/printf 引数上限) は plan 参照。`lineImplicit: 1`
(既定 0 = 挙動不変)。blockDPLUR==1 専用、lowMachPrecond>=2 と併用不可。

### v2: factor/solve 分離・K 凍結・粘性結合・粘性 dt 割引 (2026-09-02)

計画: [`plans/accepted/time_integration-line-implicit-viscous-v2.md`](../../plans/accepted/time_integration-line-implicit-viscous-v2.md)。

- **factor/solve 分離 (常時, 厳密)**: D̃ の LU 分解・W=D̃⁻¹Knext・Kprev·W は rhs に依存しない
  ため `lineThomasFactor_d` (storeLU 時 1 回) と `lineThomasSolve_d` (毎 sweep, 保存因子で代入
  のみ) に分離。v1 モノリシックの毎 sweep 再分解が DDES で 2.44× だったコストの主犯で、分離
  だけで 1.59× に落ちる。`FORGE_LINE_MONO=1` で旧動作。因子は `line_LU_d`/`line_piv_d`。
- **`lineKFreeze: 1`**: dual-time サブ反復間で K 抽出・LU 分解を凍結 (subiter 0 のみ構築)。
  LHS 近似の強化で収束経路のみ変化 (実測で収束軌道は非凍結と一致)。コスト 1.32× まで低減。
- **`lineViscCoupling: 1`**: line 面にスカラー粘性結合 K += α·I (α=ν_eff·δ/dcc)、対角は
  2α→α で line 内に真の拡散行 [−α, 2α, −α] を完成。**圧縮性 pseudo-dt の壁法線律速は音響
  (λ_visc/λ_ac = 2ν/(Δn·c) ≪ 1) なので効果は僅差** — 意味を持つのは Δn < 2ν/c の超極薄セルのみ。
- **`lineViscousDtRelief: θ`**: on-line セルの擬似 dt 粘性スペクトル半径を (1−θ) 倍 (`setDT_d`
  で面ごとに割引、対流+音響分は残す)。θ=1 でも安定 (上と同じ理由で利得も僅差)。
- **`lineDtDirectional: 1`**: 方向別 dt — line 面 (Thomas が厳密に解く結合) の λ を音響込みで
  CFL の max から除外し、Δτ を off-line 面 (lag 側) の λ だけで決める。**注意: 除外は
  `line_prev/next` に一致する内部面のみで、壁ノードの境界半割面は残る** — 壁 CV 自身の Δτ は
  境界面の音響制約のまま、Δτ が streamwise 基準 (×AR) に伸びるのは壁の 1 つ内側以降。
  case/39 DDES で cp4 安定・ωバースト最大 9.5→6.15 (directional による低減, 帰属機構は未分離)。
- **実測の価値は pseudo-CFL 引き上げ** (case/39 ny160 DDES): point は cfl_pseudo 2 で発散、
  line は 8 まで安定。**cfl_pseudo 4 + nSub 13 で point (cp1+nSub20) の 0.88 倍時間・同品質**、
  同時間ならより深い収束・ωバースト低減。定常 (M6) と dual-time DDES で律速モードが違う点に注意。
  サブ反復収縮の残る律速候補は off-line lag / segregated SST / 2次 KEEP RHS×1次 FVS LHS の
  defect-correction 不整合の 3 者。

## 既知の TODO / 注意点

- 非定常 dual-time 陰解法（`tI==11 && unsteady==1 && dualTime==1`）は実装済（2026-06、`blockDPLUR==1` のみ、物理 $\Delta t$ 固定 `control=0`）。`implicitCorrection_d.cu` の `dualtime_explicit_d` は SLAU/Roe 用の別系統補助で本流とは独立（未使用）。
- scalar 対角陰解法（`tI==11 && blockDPLUR==0`）は有効（2026-06）。block より低 `cfl_pseudo` で収束も遅いため既定は block（上記比較参照）。
- block 陰解法でも SST(k/ω) は `applySSTPointImplicit` で segregated point-implicit 更新され、**凍結しない**（2026-06）。消散+輸送(移流+拡散)の対角を陰化、生産・近傍 ΔQ は lagged。
- `matrix mat_ns` は陰解法では未使用だが非陰解法のシグネチャに残るため `StepContext` に保持（除去は別途）。
- 旧 CPU の [`update.cpp`](../../solver_density_cuda/update.cpp) は使用されていない。
- 局所時間刻みの粘性スペクトル半径寄与は `setDT_d` 側で個別に実装されている。詳細は同ファイルを参照。
