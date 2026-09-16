# 時間積分 — 理論

forge の時間積分は次の 3 系統を提供する。

| `cfg.timeIntegration` | 概要 |
| --- | --- |
| `1` / `3` | Jameson 型多段陽解法 (擬時間ステップ) |
| `4` | 4 段 4 次 Runge–Kutta (低 storage、係数は `coef_DT_4thRunge`, `coef_Res_4thRunge`) |
| `11` | 陰解法 (block DPLUR、LU-SGS $A^\pm$)。`blockDPLUR`=1(5×5 ブロック, 既定)/0(スカラー対角)。定常・非定常 dual-time とも実装済 |

> **実装状態 (2026-06)**: `timeIntegration == 11` は **block DPLUR (`blockDPLUR == 1`)・定常 (`unsteady == 0`)** のみを正式サポートする。
> スカラー対角版 (`blockDPLUR == 0`) も有効（block より低 `cfl_pseudo`・低速）。非定常 dual-time (`unsteady=1 && dualTime=1`) も実装済（`blockDPLUR=1` のみ、物理 $\Delta t$ 固定）。

## 共通枠組み

時間 1 ステップは内部ループ (`loop = 0 .. nInner-1`) に分かれ、
1 ループあたり次の手順を踏む。

```
applyBconds → gradient → limiter → ducrosSensor
   → convectiveFlux + viscousFlux  (res_* に累積)
   → implicitCorrection (任意)
   → timeIntegration_d_wrapper (Q を更新)
   → dependentVariables (派生量)
```

`res_*` は対流・粘性両方からの符号付き面寄与の合計。これを各セルで体積で
正規化したものが $-\partial \mathbf{Q}/\partial t \cdot V$ に相当する。

## 陽解法 (Jameson 多段)

低 storage 形式: 1 ステージあたり

$$
\mathbf{Q}^{(k+1)} = c_N \mathbf{Q}^N + c_M \mathbf{Q}^M
+ c_R \frac{\Delta t_{\text{loc}}}{V}\,\mathbf{R}^{(k)},
$$

ここで $\mathbf{Q}^N$ は前回外部ステップ、$\mathbf{Q}^M$ は前回内部ステージの値。
係数 $c_N, c_M, c_R$ は `solverConfig` の `coef_N`, `coef_M`, `coef_Res` に格納する
(段数 = `cfg.nInner`)。

定常計算では局所時間刻み $\Delta t_{\text{loc}}$ を使うことで収束を加速する。

## 4 段 4 次 Runge–Kutta (low-storage)

参考: Allampalli et al. (2009)。`coef_DT_4thRunge`, `coef_Res_4thRunge` を用い、
内部 4 段で

$$
\mathbf{Q}^{(k+1)} = \mathbf{Q}^N + c_{DT}^{(k)} \frac{\Delta t_{\text{loc}}}{V}\,\mathbf{R}^{(k)},
$$

累積残差 $\mathbf{R}_m = \sum_k c_R^{(k)} \mathbf{R}^{(k)} \,\Delta t_{\text{loc}}/V$ を保持し、
最終段で $\mathbf{Q}^{N+1} = \mathbf{Q}^N + \mathbf{R}_m$ とする。LES など時間精度の必要な解析向け。

## 陰解法 block DPLUR (LU-SGS 型 $A^\pm$ 分割)

擬似時間 $\tau$ を導入し、$\Delta \mathbf{Q} = \mathbf{Q}^{n+1} - \mathbf{Q}^n$ を未知に取る。
1 次精度 upwind フラックス $\mathbf{F}_f = \tfrac12(\mathbf{F}_i+\mathbf{F}_j) - \tfrac12|\widetilde{A}_f|(\mathbf{Q}_j-\mathbf{Q}_i)$ の
線形化から、残差 $\mathbf{R}_i=\sum_f \mathbf{F}_f S_f$ の Jacobian は**通量分割行列**

$$
A^{\pm} = \tfrac12\big(\widetilde{A} \pm |\widetilde{A}|\big) = R\,\Lambda^{\pm} L,\qquad
\Lambda^{+}=\max(\Lambda,0),\ \ \Lambda^{-}=\min(\Lambda,0)
$$

で表される。これを使うと **対角ブロックは upwind スキームの厳密な自己 Jacobian** $\partial\mathbf{R}_i/\partial\mathbf{Q}_i=\sum_f A^{+}_f S_f$、
近傍ブロックは $\partial\mathbf{R}_i/\partial\mathbf{Q}_j = A^{-}_f S_f$ となり、解くべき線形系は

$$
\underbrace{\left[\frac{V}{\Delta \tau} I + \sum_f A^{+}_f S_f + \sum_f \rho^{\nu}_f S_f\, I\right]}_{D_i\ (\text{対角ブロック})} \Delta \mathbf{Q}_i
= -\mathbf{R}_i(\mathbf{Q}^n) - \sum_f A^{-}_f S_f\, \Delta \mathbf{Q}_{\text{nbr}}
$$

となる（$\rho^{\nu}_f$ は粘性スペクトル半径による対角強化）。ここで $\widetilde{A}_f=R\Lambda L$ は
セル中心状態から作る面法線フラックス Jacobian（固有値 $\Lambda=\{U{+}c,\,U,\,U,\,U,\,U{-}c\}$）で、$R\Lambda L$ は
真の流束 Jacobian を厳密に再現し、$|\widetilde{A}|=R|\Lambda|L$ の固有値は全て非負（matrix-free・近似 BC）。

> **2026-06 修正（重要）**: 旧実装は対角・近傍ともに**絶対値 Jacobian $|\widetilde{A}|$ を符号付き $A^\pm$ の代わりに誤用**しており、
> 対角が upwind 自己 Jacobian と一致せず（対流中心項を欠く）、近傍結合の符号も誤っていた。このため block DPLUR は
> 新旧コードとも **一度も収束せず発散**していた（対角のみ＝頭打ち、近傍 sweep 有効＝NaN）。
> $A^{+}$ 対角・$A^{-}$ 近傍の正しい LU-SGS 分割に置換することで反復行列のスペクトル半径が
> $\rho \approx \dfrac{V/\Delta\tau}{V/\Delta\tau+\lambda^{+}} < 1$ となり収束する。詳細は [implementation.md](implementation.md)。

### 一般EOS固有系 (TP gas) — 凍結 γ の閉形式は CPG 専用

$\widetilde A = R\Lambda L$ の固有ベクトルは EOS に依存する。閉形式 FVS は EOS を**単一の $\gamma$** で表し、
音響右固有ベクトルのエネルギー成分に $H_t = K + c^2/(\gamma-1)$（$K=\tfrac12|\mathbf u|^2$）、接触波に $K$ を使う。
これは **CPG (定数 $c_p$) でのみ厳密**な 2 つの恒等式

$$ \kappa \equiv \left.\frac{\partial p}{\partial(\rho e)}\right|_\rho = \gamma-1,\qquad h = \frac{c^2}{\gamma-1} $$

に依存している。**TP (thermally-perfect, $c_p(T)$) では②が破れる**: $h(T)=h(T_\text{ref})+\int_{T_\text{ref}}^T c_p\,d\theta \ne c_p(T)\,T = c^2/(\gamma(T)-1)$。
等価に、CPG で恒等的に 0 の項 $\chi \equiv \left.\partial p/\partial\rho\right|_{\rho e}$ が TP では非ゼロになる。一般EOSでは

$$ c^2 = \chi + \kappa\,h,\qquad \chi = c^2 - \kappa\,h,\qquad \kappa = \gamma(T)-1 = \frac{R}{c_v(T)} $$

（固定組成 TP 理想気体）。閉形式は暗黙に $\chi=0$ を仮定するため、TP では $\widetilde A=R\Lambda L$ が
**真の $\partial\mathbf F/\partial\mathbf Q$ と一致しない**（FD 検証で相対誤差 $\approx 4.7$、すなわち別方程式を線形化している。
`case/16.nozzle_wys` の CPG/TP 対照で TP のみ $c_{\rm fl}\approx2$ 頭打ち・残差プラトーになる真因）。

**一般EOS右固有ベクトル**（法線 $\mathbf n$、接線 $\mathbf t_1,\mathbf t_2$、$U_n=\mathbf u\cdot\mathbf n$、固有値 $U_n-c,\,U_n,\,U_n,\,U_n,\,U_n+c$）:

$$
\mathbf r_\mp=\begin{bmatrix}1\\ \mathbf u\mp c\mathbf n\\ H_t\mp cU_n\end{bmatrix},\quad
\mathbf r_c=\begin{bmatrix}1\\ \mathbf u\\ H_t-\dfrac{c^2}{\kappa}\end{bmatrix},\quad
\mathbf r_{s k}=\begin{bmatrix}0\\ \mathbf t_k\\ \mathbf u\cdot\mathbf t_k\end{bmatrix}
$$

ここで $H_t=h(T)+K$ は**実 NASA 全エンタルピー**。接触波のエネルギー成分は $H_t-c^2/\kappa = K-\chi/\kappa$ で、
CPG では $\chi=0$ ゆえ $K$ に簡約され従来式に戻る。$\Lambda^\pm$ 分割は固有値のみに作用するので不変。

数値検証（`/tmp/eos_eig_verify.py` → plan `time_integration-general-eos-jacobian.md` で恒久化）: 上記固有系は
**CPG・TP（NASA N2）双方で $\|R\Lambda L - A_\text{FD}\|/\|A_\text{FD}\|\approx 3\times10^{-8}$**（FD 機械精度）に一致。
実装は $L=R^{-1}$ を**部分ピボット付き double 5×5 LU** で数値的に解く（Method A）。$H_t,c^2,\kappa,\chi$ は
**同一セル状態・同一 thermo 評価**から作る（$c^2=\gamma RT,\ \kappa=\gamma-1,\ h=e+RT$ を同じ $T$・組成・datum で）。
CPG 経路は閉形式のまま分岐保持（回帰基準、ビット不変）。実装条件・検証レベルは [implementation.md](implementation.md)。

### 古典 DPLUR 反復 (残差固定 + 複数 sweep)

非線形残差 $\mathbf{R}(\mathbf{Q}^n)$ と対角ブロック $D_i$ を**1 回の非線形更新あたり 1 度だけ構築**し、
$\mathbf{Q}^n$ を固定したまま $\Delta\mathbf{Q}$ を Jacobi 型に複数回緩和する（lagged neighbor）。

$$
\Delta\mathbf{Q}_i^{(s+1)} = D_i^{-1}\Big(-\mathbf{R}_i - \sum_f A^{-}_f S_f\, \Delta\mathbf{Q}_{\text{nbr}}^{(s)}\Big),
\quad s = 0 \ldots n_{\text{sweep}}-1
$$

各 sweep で隣接セルの $\Delta\mathbf{Q}$ は前 sweep の値（`dq_old`）を参照し（Gauss–Seidel ではなく Jacobi）、
sweep 間で `dq_old`/`dq_new` バッファを入れ替える。全 sweep 後に **1 度だけ** $\mathbf{Q} = \mathbf{Q}^n + \Delta\mathbf{Q}$ を commit する
（sweep 中は $\mathbf{Q}$ を固定するのが古典 DPLUR の要件）。`implicitRelax` で $\Delta\mathbf{Q}$ を緩和する。

`build_jacobian_split` がセル中心状態から $A^{+}$ と RHS 用の $-A^{-}=\tfrac12(|\widetilde{A}|-\widetilde{A})$ を同時に構築、
`solve_5x5` が部分ピボット付き Gauss 消去で $D_i^{-1}(\cdots)$ を解く（ピボット `< 1e-20` でゼロ解にフェイルセーフ）。
擬似 CFL（`cfl_pseudo`）を上げるほど $V/\Delta\tau$ が小さくなり $\rho$ が下がる（収束加速）。

### 反復回数の意味付け

| config | 意味 |
| --- | --- |
| `nStepInner` | **DPLUR 緩和 sweep 回数** $n_{\text{sweep}}$（1 非線形更新内の線形反復） |
| `nStepOuter` | 定常では擬似時間ステップ数（メインループ＝擬似時間行進） |

定常計算はメインループ（擬似時間）を回し、1 ステップあたり 1 回の非線形更新
（残差構築 → DPLUR sweep ×`nStepInner` → commit）を行う。

### 低マッハ前処理固有値 (Weiss–Smith) — 試行したが不採用

> **結論: block DPLUR では有害のため不採用**（2026-06 検証。実装後 revert）。理論メモとして残す。

着想は、$A^\pm$ 分割に使う固有値を基準速度 $U_r=\min(c,\max(|\mathbf u|,\epsilon c))$、$\beta=(U_r/c)^2$ で
前処理した

$$
\lambda_{1,2,3}' = U_n,\qquad
\lambda_{4,5}' = \tfrac12(1+\beta)U_n \pm \tfrac12\sqrt{(1-\beta)^2 U_n^2 + 4U_r^2}
$$

に差し替え、低マッハの音響 stiffness を除く、というもの。**しかし block DPLUR (LU-SGS $A^\pm$) では
これが逆効果**だった。理由: 対角 $D_i=V/\Delta\tau\,I+\sum_f A^{+}_f S_f$ の**対角優位性の源は大きい音響
固有値 $U_n\pm c$ を含む $A^{+}$** であり（2026-06 の収束化はこれに依拠）、$U_n\pm c'$ (小) に落とすと
$A^{+}$ が縮んで対角優位性が下がる。検証では、フラックス散逸前処理単独 ($\epsilon=0.15$) では安定だった
ものが、固有値前処理を加えると発散し、**安定 $\epsilon$ 範囲をむしろ狭めた**（収束加速も根治もなし）。
あわせて検討した setDT のスペクトル半径前処理 ($c\to c'$ による $\Delta t_{\rm loc}$ 増大) も、
$V/\Delta\tau$ を縮めて対角を二重に潰し、最速で発散した。

根因は、`build_jacobian_split` が Euler 固有ヤコビアンで**SLAU 低マッハ散逸項を表さない**こと。固有値だけ
前処理しても aggressive なフラックス前処理を陰的に安定化できない。真の安定化には固有ベクトルも前処理する
**完全 preconditioned-flux Jacobian**（大改修）が要る。よって現状はフラックス散逸前処理のみ採用し、LHS は
従来の $A^\pm$ のまま。詳細・データは計画
[`time_integration-lowmach-preconditioning.md`](../../plans/accepted/time_integration-lowmach-preconditioning.md) §9。

### 低マッハ前処理固有系 (Weiss–Smith・完全 $\Gamma^{-1}A$) — Phase 4 (採用・振動根治)

固有値だけの前処理が失敗したので、**固有ベクトルまで前処理した完全形**を `lowMachPrecond=2` の新モードとして
実装する (Phase 1 の flux 散逸 `c'` とは排他)。要点は、擬似時間項とフラックス分割を**同一の前処理計量**で
揃えること。

**保存変数前処理行列 $\Gamma_c$ (ランク 1 閉形)。** 原始変数 $Q_p=(p,u,v,w,T)$ の Weiss–Smith 前処理
$\Gamma_p$ は変換 $M=\partial Q_c/\partial Q_p$ に対し $\rho_p=\partial\rho/\partial p|_T$ のみを
$\Theta=1/U_r^2+(\gamma-1)/c^2$ に置換したものである ($U_r=c$ で $\Theta=\rho_p$)。保存形の擬似時間項に必要な
$\Gamma_c=\Gamma_p M^{-1}$ は、$\Theta-\rho_p=(1-\beta)/(\beta c^2)$ を使って**恒等行列のランク 1 更新**に閉じる:

$$
\Gamma_c = I + \frac{1-\beta}{\beta c^2}\, g\, r^\top,\qquad
g=(1,u,v,w,H)^\top,\quad
r^\top=(\gamma-1)\,(e_k,-u,-v,-w,1)=\frac{\partial p}{\partial Q_c},
$$

ここで $H=c^2/(\gamma-1)+e_k$ (全エンタルピ)、$e_k=\tfrac12|\mathbf u|^2$、$Q_c=(\rho,\rho u,\rho v,\rho w,\rho E)$。
$r^\top g = c^2$ が成り立ち、Sherman–Morrison で逆行列も閉じる:

$$
\Gamma_c^{-1} = I - \frac{1-\beta}{c^2}\, g\, r^\top.
$$

$\beta=1$ (超音速、$U_r=c$) で両者とも係数が 0 となり $\Gamma_c=\Gamma_c^{-1}=I$ に**厳密復帰**する。
これらは `lowMachPrecond_d.cuh` の `lowMachGammaC` / `lowMachGammaCinv` が `double[5][5]` に組む
(条件数 $\sim1/\beta$ のため倍精度。ランク 1 形・$r^\top g=c^2$・$\beta{=}1$ 復帰は記号的に検証済)。

**前処理ヤコビアンと固有系。** 安定性・CFL を律する行列は $P=\Gamma_c^{-1}A_c$ で、その固有値は前処理済み

$$
\lambda'_{1,2,3}=U_n,\qquad
\lambda'_{4,5}=\tfrac12(1+\beta)U_n\pm c',\quad
c'=\tfrac12\sqrt{(1-\beta)^2U_n^2+4U_r^2},
$$

**保存形の厳密な陰解法系は前処理を時間項のみに掛ける。** 原始系から保存形へ移すと
$$
\big(\Gamma_c\,V/\Delta\tau' + A_c\big)\,\Delta Q = -\mathbf R
$$
となる。$A_c$ は**残差の真のヤコビアン**なので**物理 (非前処理)** でなければならず、前処理は擬似時間項 $\Gamma_c$ に入る。
収束は積 $(\Gamma_c V/\Delta\tau')^{-1}A_c=(\Delta\tau'/V)\,\Gamma_c^{-1}A_c$ の固有値 $(\Delta\tau'/V)\lambda'$ で**一様に前処理される**
(フラックスを物理のままにしても収束は前処理される)。したがって LHS は既存と同じ物理 FVS 分割をそのまま使う。

**LHS 構成 (block DPLUR)** は
$$
D = \Gamma_c\,\frac{V}{\Delta\tau'} \;+\; a\frac{V}{\Delta t}I \;+\; \sum_f A^{+}_{c,f}\,S_f \;+\;(\text{粘性}),\qquad
\text{近傍} = -\sum_f A^{-}_{c,f}S_f
$$
で、$A^{+}_c=a^{+}$・$-A^{-}_c=k_{\rm off}$ は既存 `build_jacobian_split` の物理分割そのもの。$\Delta\tau'$ は前処理スペクトル半径
$\rho'=\tfrac12(1+\beta)|U_n|+c'$ 基準に拡大する (`setDT_d` で $\Delta\tau'=\Delta\tau\cdot(|\mathbf u|+c)/\rho'$、低マッハで $\sim1/\epsilon$ 倍・有界)。
dual-time の物理 BDF 項 $aV/\Delta t\,I$ は**非前処理**。悪条件 ($\sim1/\beta$) は $\Gamma_c$ 時間項だけが持ち込むが、
$\Gamma_c=I+\alpha g r^\top$ がランク 1 ゆえ $D=D_0+\gamma g r^\top$ ($D_0$=物理ブロック・良条件) と分解でき、**Sherman-Morrison**
$x=y-[\gamma(r^\top y)/(1+\gamma(r^\top z))]z$ ($y=D_0^{-1}b,\ z=D_0^{-1}g$) で**$D_0$ 解は float**、悪条件は分母スカラーのみ
double に隔離する (FP64 を排し物理 block とほぼ同速)。$\beta=1$ (超音速) では $\Gamma_c=I,\ \Delta\tau'=\Delta\tau$・フラックス同一ゆえ
既存カーネルと解一致。

**前処理の2面と $U_r$ フロア $\epsilon$ (要点)。** Weiss–Smith 前処理は同じ $\Gamma$(同じ基準速度 $U_r$)が**2つの面**で効く:

1. **散逸面 (精度・振動・解を変える)**: 数値フラックスの散逸を前処理音速 $c'$ でスケールし、低マッハで $O(1/M)$ に
   過大化する upwind 散逸 (Guillard–Viozat) を是正する。**収束解そのものを変える**面で、forge では `lowMachPrecond=1`
   ([`convection/theory.md`](../convection/theory.md) の `c'` 散逸) が担う。**低マッハ自励振動を消すのはこの面**。
2. **時間項面 (収束・安定・解を変えない)**: 擬似時間微分項に $\Gamma_c$ を掛ける本節の前処理。**擬似時間項なので収束解は
   不変**で、役割は剛性除去による収束加速と、強い散逸 (小 $\epsilon$) を発散させずに回す安定化。forge では `lowMachPrecond=2`。

標準 Weiss–Smith ではこの2面が**束ねて**効く (同じ $\Gamma$)。**振動を直すのは 1.(散逸)であって 2.(時間項)単独では直らない** —
時間項は解を変えないため。今回 $\epsilon=0.05$ で振動が消えたのは、消したのは 1. の `c'` 散逸で、2. の $\Gamma_c$ 時間項は
その強い散逸を安定収束させる役だった (物理 LHS=`=1` は $\epsilon=0.05$ で発散、前処理 LHS=`=2` で初めて回せる)。

$\epsilon$ (`precondEps`) は $U_r=\min(c,\max(|\mathbf u|,\epsilon c))$ の**停留点フロア**で、$|\mathbf u|\to0$ で
$U_r\to0\Rightarrow\beta\to0$ となり $\Gamma_c$ が特異化 ($1/\beta\to\infty$) するのを防ぐ。これは**標準 Weiss–Smith に必須の
要素** (文献により $K|V_\infty|$・粘性/非定常速度スケール等、$\epsilon c$ は最も素朴な形) で、追加の細工ではない。
副作用は $\epsilon$ を下げるほど: 散逸増 (振動は消えるが**真の低マッハ圧力変動も削る**)、条件数 $\sim1/\beta=1/\epsilon^2$ 悪化
(倍精度で吸収するが下限あり)、$\Delta\tau'\sim1/\epsilon$ 拡大で**安定限界に接近**。1.と 2.が同じ `precondEps` 一本で
動くため独立調整できない点も含め、ドメイン最小マッハに応じたスイートスポット (本ケースで $\epsilon=0.05$–$0.15$) がある。

> **当初の誤り (記録)**: 一旦フラックスも $\hat A_\Gamma=\Gamma_c^{-1}A_c$ のスペクトル半径分割 (Rusanov 型) で前処理したが、
> これは解くべき系の $A_c$ を $\Gamma_c^{-1}A_c$ に置換した**別系**で過散逸だった。保存形では前処理は時間項のみが正しい (上式)。

> **結果 (2026-06-09・採用)**: `implicit_defect_correction_block_precond_d` (物理 FVS + $\Gamma_c$ 時間項 + 前処理 $\Delta\tau'$、倍精度)、
> `setDT_d` の $\Delta\tau'$ スケーリング、SLAU の `c'` 散逸を `lowMachPrecond>=1` に拡張、wrapper の `==2` dispatch を実装。
> **`case/23.axi_nozzle` で低マッハ自励振動を根治**: Phase1 (物理 LHS) が発散した $\epsilon=0.05$ を前処理 LHS が安定化し、
> chamber 圧振幅 (M<0.08, 4k–20k) を 0.882%→**0.087% (定常収束・振動消滅)**、$\epsilon=0.15$ でも 0.603%→0.333% (超音速域不変)。
> 安定 `cfl_pseudo` は m1~1→m2~7-9 と広がる。当初は倍精度で per-step 2.54× で収束加速は互角だったが、
> Sherman-Morrison 化で per-step を物理 block と同速 (15.7 vs 15.6ms) にした結果、cfl 余裕がそのまま wall-clock
> 加速に転化し **rms_roe<1.0 到達で約 2.2× の net 加速** (深い収束ほど有利)。計画
> [`time_integration-lowmach-preconditioning.md`](../../plans/accepted/time_integration-lowmach-preconditioning.md) §9 `2026-06-09`。

## 非定常 dual-time 陰解法 (実装済み 2026-06)

非定常は物理時間 $t$ に対し BDF（後退差分）で物理時間微分項を残差に加え、各物理ステップ内で
擬似時間サブ反復により残差を収束させる二重時間刻み法を用いる。残差は

$$
\mathbf{R}^*_i = \mathbf{R}_i(\mathbf{Q}) + \frac{V}{\Delta t}\big(a\mathbf{Q} - b\mathbf{Q}^{n} + c\mathbf{Q}^{n-1}\big)
$$

の形（BDF1: $a{=}1,b{=}1,c{=}0$／BDF2: $a{=}\tfrac32,b{=}2,c{=}\tfrac12$）で、物理時間項は対角 $D_i$ にも $aV/\Delta t$ として加わる。
第 2 時間レベル $\mathbf{Q}^{n-1}$ は `roNN` 系に保持する。擬似時間サブ反復・DPLUR sweep の構造は定常と共有する
（3 階層: 物理ステップ × 擬似時間サブ反復 × DPLUR sweep。定常は前 2 つが縮退）。

| config | 意味 |
| --- | --- |
| `nStepOuter` | 物理時間ステップ数 |
| `nSubIterDualTime` | 物理ステップ内の擬似時間サブ反復数（既定 20） |
| `bdfOrder` | 物理時間 BDF 次数 (1 or 2、初回ステップは自動的に BDF1) |
| `nStepInner` | 各サブ反復内の DPLUR sweep 回数（定常と同義） |

使用条件: `unsteady=1, dualTime=1, timeIntegration=11, blockDPLUR=1, time.deltaT.control=0`（物理 $\Delta t$ 固定）。

> **化学種・受動種の物理時間項 (2026-09-17 実装中, plan [species-passive-scalar-unification](../../plans/active/species-passive-scalar-unification.md) §4.4)**: 2026-09-16 までの
> main は dual-time で化学種を更新しておらず (ρ が動くと ΣY=ρ⁰/ρ になる; [[dualtime-species-frozen-bug]])、凝縮モーメント・トレーサは物理時間項なしの擬似時間更新だった。
> chem ブランチ (e296f0d0) の修正を移植し、化学種 $\rho Y_s$ と受動種 ($\rho\xi$, モーメント) に時間レベル `*P/*PP` と BDF 残差 $-(V/\Delta t)(a\,\rho\phi - b\,\rho\phi^P + c\,\rho\phi^{PP})$・
> 対角 $Va/\Delta t$ を付ける。履歴有効数と BDF 係数は流れ・化学種・受動種で共有する 1 つの状態 (最初の 1 物理 step は BDF1)。サブ反復の処理順: 空間残差 (+ピン除去) →
> 周期 gather → BDF 項 (合併体積で一度) → ピン再除去 → 流れ block 更新 → 化学種更新 (coupling 0/1/2) → 再正規化・primitive → 入口再適用 → 受動種更新 → 状態 mirror。
> checkpoint は流れ・化学種・受動種の履歴・時刻・刻み・履歴有効数を一括で書き、1 つでも欠ければ全系を BDF1 から再開する。
>
> **受動種の物理 step 末尾の保存的 FCT 補正 (`passiveFct`, plan §4.7 v5, 2026-09-17)**: BDF2 は $b=a+c$ から
> $\tfrac{V}{\Delta t}(q^{n+1}-q^n)=\tfrac1a[\sum_f s_fF_f(q^{n+1})+SV]+\tfrac ca\tfrac{V}{\Delta t}(q^n-q^{n-1})$ と書け、前 step の増分率を面流束 $G^n$ と局所残り $H^n$ に
> 厳密分解 (確定状態から $G^{n+1}_f=F_{L,f}(q_L)+\alpha_fA^{pre}_f$, $H^{n+1}=\tfrac{V}{\Delta t}(q^{n+1}-q^n)-\sum s_fG^{n+1}$) すると「BE + 有効面流束 $F^{eff}=\tfrac1aF+\tfrac caG^n$」になる。
> sub-iter 終了後に (1) 終了状態で残差を再評価、(2) 同じ履歴・凍結ソース・拡散で移流だけ 1 次風上、質量流束を流れの BE 形連続式と整合する $\dot m^{eff}=\tfrac1a\dot m+\tfrac ca\dot m^{eff,n}$
> にした低次 BE 陰解 $q_L$ (M 行列; node 逆流境界は $\phi^n_{own}$ の定数 RHS) を Jacobi で解き、(3) 生の反拡散流束 $A^{raw}_f=F^{eff}_{H,f}-F_{L,f}(q_L)$ (全輸送面)、
> 基点 $q_B=q_H-M^{-1}BA^{raw}$、(4) Zalesak ($P^\pm$, 局所極値 $\{\phi^n,\phi_B\}$ ∩ 物理限界, $R^\pm$, 面共有 $\alpha_f$; 境界は外部側 $R=1$)、
> (5) $q_C=q_H-M^{-1}B(A^{raw}-\alpha A^{pre})$ (保存; $\alpha=1$ の面は $q_H$ のまま = 時間 2 次の完全陰的スキーム不変)。後処理は入口再適用 → floor → $g$ クランプ →
> 二相 EOS 更新 → その $T$ でモーメント実現可能性射影 → $G/H/\dot m^{eff}$ (checkpoint 必須履歴)。定常 (局所 Δτ) と RK には掛けない。有界性の穴 (負の有効ソース、
> 前 step の局所残り、反復残差) は逸脱量として積算し、`tools/check_passive_budget.py` が総量比 1e-6 で FAIL にする。
擬似時間 $\Delta\tau$ は `cfl_pseudo` から決まり（物理 $\Delta t$ とは独立）、物理 $\Delta t$ は固定。
commit は in-place（$\mathbf Q \leftarrow \mathbf Q + \delta\mathbf Q$。`roN`=$\mathbf Q^n$ は BDF 基準で固定のため）。

> **検証 (2026-06, `case/20.naca_ml`)**: 各物理ステップ内で擬似サブ反復が $\mathbf R^*$ を ~4 桁低減
> （例: rms_roe 22.9→0.0035, サブ反復 ~10 で収束）。陽解法 CFL 上限を超える $\Delta t$ でも安定
> （$\Delta t{=}2\times10^{-6}$, 陽解法安定限界 ~$9\times10^{-7}$）。同一 $\Delta t{=}5\times10^{-7}$ では陽解法
> 時間精度解と壁面静圧が一致（平均 0.006%、最大 0.2% は RK3 と BDF2 の時間スキーム差 $O(\Delta t^2)$）。

## スカラー (k/ω) の陰解法 (segregated point-implicit)

平均流 5 式とスカラー (k/ω) は常に別カーネルで分離 (segregated) して積分する。
block 陰解法 (`timeIntegration==11`) では、平均流 block DPLUR を解いて commit した後、
**同じ擬似時間ステップで k/ω を segregated point-implicit で更新**する（旧実装はここで k/ω を凍結していた）。

SST ソースの消散項は $\mathbf{Q}$ に比例する負ソースで、時間刻みに依らず stiff:

$$
D_k = \beta^\* \rho k \omega = \beta^\*(\rho k)\,\omega,\qquad
D_\omega = \beta \rho \omega^2 = \beta \frac{(\rho\omega)^2}{\rho}
$$

$$
\frac{\partial D_k}{\partial(\rho k)} = \beta^\*\omega,\qquad
\frac{\partial D_\omega}{\partial(\rho\omega)} = 2\beta\omega
$$

これを対角に陰に取り込む point-implicit 更新を行う:

$$
D_\phi = \frac{V}{\Delta\tau} + V\,\frac{\partial D}{\partial(\rho\phi)} + \Lambda^{T}_\phi,\qquad
\delta(\rho\phi) = \mathrm{relax}\cdot\frac{\text{res}_{\rho\phi}}{D_\phi},\qquad
\rho\phi \leftarrow \max(\rho\phi^{N} + \delta(\rho\phi),\ \text{floor})
$$

ここで $\Lambda^{T}_\phi$ は次節の**輸送項 (移流+拡散) の対角**である。
擬似時間 $\Delta\tau$ は平均流と同じ `dt_local` を流用する。realizability のため $\rho k \ge 0$、$\rho\omega > 0$ を課す。
近傍 $\Delta Q_{nbr}$ 結合・生産項は持たない純 point-implicit（消散+輸送の対角のみを陰化する最小構成）。
これにより block 陰解法でも k/ω が凍結せず収束し、乱流ケースを陰解法で回せる。

### 輸送項 (移流+拡散) の point-implicit 対角 $\Lambda^{T}_\phi$ (2026-06)

消散項のみを陰化した初期実装では、壁法則検証ケース `case/26.flat_plate_sst`（乱流平板, M=0.2,
$y^+<1$）で安定 `cfl_pseudo` が発達場 restart でも **5〜6** に制限されていた。発散の素因を切り分けると、
平均流 5 式は block DPLUR で安定（残差フロアに静止）したまま **k/ω だけが先に指数発散**し、
`scalarDiffusion=0`（拡散を切る）と発散が消えることから、律速は **k/ω 輸送項の陽的（lagged）扱い**
だと判明した。壁第一セル（高さ $\sim4\times10^{-6}$ m, 高アスペクト比）で陽的拡散の安定限界
$\Delta t < \delta^2/(2\nu_{\text{eff}})$ を、また外層で陽的移流の CFL 限界を、`cfl_pseudo`$\gtrsim6$ で超える。

平均流 block DPLUR が対流ヤコビアン $A^\pm$・粘性スペクトル半径を対角に持つのと同様に、k/ω の
**輸送スペクトル半径を point-implicit 対角に加える**:

$$
\Lambda^{T}_\phi = \underbrace{\sum_f \frac{\max(\dot m_f,0)+\max(-\dot m_f,0)}{\rho}}_{\text{移流 (1次風上)}}
  \;+\; \underbrace{\sum_f \frac{\mu_{\text{face}}}{\rho}\frac{|\delta_f|}{|\Delta\mathbf{r}_f|}}_{\text{拡散}}
  \;\;[\mathrm{m^3/s}]
$$

各項は 1 次風上移流／2 点拡散の流束を $\rho\phi$ で微分した対角寄与 $-\partial\text{res}/\partial(\rho\phi)$
で、流出側セルにのみ正値を生む（$\dot m_f$ は面質量流束、$\mu_{\text{face}}=\mu+\sigma_\phi\mu_t$、
$\delta_f$ は拡散の幾何係数）。$\Lambda^{T}_\phi\ge 0$ は単位 $[\mathrm{m^3/s}]$ で $V/\Delta\tau$ と同じ
ため、$D_\phi$ にそのまま加える（src_jac のように $V$ を掛けない）。

これは block DPLUR と同じ **defect-correction** であり、**定常解を変えない**。$\Lambda^{T}_\phi$ は更新量
$\delta(\rho\phi)=\mathrm{relax}\cdot\text{res}_{\rho\phi}/D_\phi$ の係数にのみ入り、収束時 $\text{res}_{\rho\phi}\to0$
では $\delta\to0$ となるため、不動点（壁法則・$C_f$）は不変で緩和経路と安定 `cfl_pseudo` の上限だけが変わる。
`case/26.flat_plate_sst` 発達場 restart で安定 `cfl_pseudo` は **5〜6 → 120**（約 20 倍）に向上し、
`cfl_pseudo=50` の壁法則・$C_f$ は baseline（`cfl_pseudo=5`, 120k step）と $<0.1\%$ で一致した（§検証）。

> **生産項の陰化は不採用**。生産は $\phi$ を増やすため真の寄与 $-\partial P/\partial(\rho\phi)\le0$ は
> 対角を弱め発散させ（Patankar も正の生産勾配は陰化しない）、逆に振幅 $|\partial P/\partial(\rho\phi)|$ を
> 正対角に足す under-relaxation も試したが、上記の通り律速は輸送項のため安定 `cfl_pseudo` を一切
> 上げず（120 で不変）、ω に至っては $P_\omega=\alpha\rho S^2$ が $\omega$ に依らず勾配ゼロ。よって採用しない。

### 陽解法 RK での point-implicit 源項 (2026-06)

陽解法 RK (`timeIntegration==1/3`) で RANS (SST) を回す場合、平均流は陽的でも
**スカラー (k/ω) の消散項だけは上と同じ消散ヤコビアンで point-implicit 化**する必要がある。
純陽的に積分すると、壁近傍の大きな $\omega$ と $D_\omega=\beta\rho\omega^2$ の stiff 性により
CFL を下げても（$\mathrm{CFL}=0.05$ でも）初回サブステップで $\omega$ が発散する（陽的安定限界が
源項の時間スケール $1/(2\beta\omega)$ で律速されるため）。

陽解法 RK のスカラー更新は本来

$$
\rho\phi \leftarrow c_N\,\rho\phi^{N} + c_M\,\rho\phi^{M}
  + c_{\text{res}}\,\frac{\Delta t_{\text{loc}}}{V}\,\text{res}_{\rho\phi}
$$

だが、残差増分のみを源項+輸送ヤコビアンで割って減衰させる:

$$
\rho\phi \leftarrow c_N\,\rho\phi^{N} + c_M\,\rho\phi^{M}
  + \frac{1}{1 + c_{\text{res}}\,\Delta t_{\text{loc}}\,(\partial D/\partial(\rho\phi) + \Lambda^{T}_\phi/V)}
    \cdot c_{\text{res}}\,\frac{\Delta t_{\text{loc}}}{V}\,\text{res}_{\rho\phi}
$$

ここで `src_jac_k`/`src_jac_omega`（消散項 $\beta^\*\omega,\ 2\beta\omega$）と `transport_diag_k`/`omega`
（輸送項 $\Lambda^{T}_\phi$）は block 陰解法と同じ対角ヤコビアンを共有する。減衰係数
$1 + c_{\text{res}}\Delta t_{\text{loc}}(\text{src\_jac}+\text{transport\_diag}/V)\ge 1$ は陰解法対角
$D_\phi = V/\Delta\tau + V\,\text{src\_jac} + \text{transport\_diag}$ に $\Delta\tau/V$ を掛けた形と整合する
（消散+輸送を陰化し、生産・近傍結合は lagged のまま）。これにより陽解法 RK でも RANS が安定に回り、
平均流の時間積分法 (陽 RK / block 陰) を揃えた比較ができる。realizability の floor は陰解法と共通。
4 段 4 次 RK (`timeIntegration==4`) のスカラー積分は本減衰を未適用（RANS 非対応のまま）。

## 局所時間刻み

`setDT_d_wrapper` (`cuda_forge/setDT_d.cu`) が各セルの $\Delta t_{\text{loc}}$ を
CFL 条件から決定する。スペクトル半径

$$
\lambda_{\max} = \sum_f \big( |U_n| + c \big)_f A_f + \text{粘性項}
$$

から $\Delta t_{\text{loc}} = \mathrm{CFL}\,V / \lambda_{\max}$ を求める。
`cfg.dt` が指定された場合はそれを上限とする。

> **setDT の低マッハ前処理は不採用** (2026-06 検証)。スペクトル半径の $c\to c'$ 置換は低マッハ域で
> $\Delta t_{\text{loc}}$ を増大させ陰解法対角 $V/\Delta\tau$ を縮める。block DPLUR では対角優位性が崩れて
> 発散するため採用しない（LHS 固有値前処理と併用しても二重に対角が潰れさらに悪化。上記「低マッハ前処理
> 固有値」節を参照）。低マッハ前処理は対流フラックスの散逸スケールにのみ適用する。

## 参考

- 各カーネル分岐の詳細: [implementation.md](implementation.md)。
- 境界条件適用 (`applyBconds`): [`methods/boundary/`](../boundary/)。
- 残差を作る対流・粘性: [`methods/convection/`](../convection/), [`methods/diffusion/`](../diffusion/)。
