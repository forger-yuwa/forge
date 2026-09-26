# 勾配再構成

forge では各セル中心におけるセル平均量 $\phi_c$ から、勾配 $\nabla \phi_c$ を
**Green-Gauss 法 (cell-centered, 線形補間)** で求める。
対流フラックス再構成・粘性フラックス・速度発散 $\nabla \cdot \mathbf{u}$ の評価に用いる。

本ドキュメントは理論(係数・方程式)と実装(ソース対応)をまとめる。

## 理論

### 支配関係

セル $c$ (体積 $V_c$、面集合 $\partial c$) について発散定理から

$$
\int_{V_c} \nabla \phi \, dV = \oint_{\partial c} \phi \, \mathbf{n}\, dS .
$$

セル中心勾配を体積平均で近似し、面ごとに代表値 $\phi_f$ を与えると

$$
\nabla \phi_c \;\approx\; \frac{1}{V_c} \sum_{f \in \partial c} \phi_f \, \mathbf{S}_f ,
\qquad
\mathbf{S}_f = \mathbf{n}_f \, A_f ,
$$

ここで $\mathbf{S}_f$ は外向き面ベクトル (面積込み)。

### 面値の補間

forge では面 $f$ を共有する 2 セル $c_1, c_2$ に対し、面ごとの幾何重み $f_x \in [0,1]$
(メッシュ生成時に算出) を用いた線形補間を採用する。

$$
\phi_f = f_x \, \phi_{c_1} + (1 - f_x) \, \phi_{c_2}.
$$

$f_x$ は面中心からセル $c_2$ 中心までの距離を $c_1$–$c_2$ 中心間距離で割った値で、
等間隔メッシュでは $f_x = 1/2$ となる。
歪んだメッシュにおける線形精度を保つための重みであり、本実装では
スキュー補正・最小二乗による高次補正は行っていない。

### 面ループによる加算則

両側のセルへ符号反転して累積する離散版を採用する。
内部面 (両側に実セルが存在) について

$$
\nabla \phi_{c_1} \mathrel{+}= \phi_f \, \mathbf{S}_f, \qquad
\nabla \phi_{c_2} \mathrel{+}= -\phi_f \, \mathbf{S}_f,
$$

最後にセルごとに $V_c$ で除する。境界面については、ghost セルを $c_2$ 側として
扱い、ghost セル値は境界条件モジュール (`methods/boundary/`) が事前に与える。

### 計算する量

各セル中心で次の勾配と発散を同時に計算する。

| 量 | 用途 |
| --- | --- |
| $\nabla U_x, \nabla U_y, \nabla U_z$ | 粘性応力テンソル, 対流再構成 |
| $\nabla T$ | 熱伝導フラックス |
| $\nabla P$ | 対流フラックス再構成 (KEEP, SLAU, AUSM 系) |
| $\nabla \rho$ (GPU 経路のみ) | 密度ベース再構成 |
| $\nabla \cdot \mathbf{u}$ | 衝撃波・人工粘性指標, 圧縮性指標 |

### 既知の制約

- 制限関数 (limiter) は適用していない。スカラー量も含め純粋な Green-Gauss。
- 高次再構成 (二次以上) は未実装。最小二乗 (LSQ) 勾配は **node-centered モード限定**で
  `gradLSQ=1` のとき利用可 (近壁 checkerboard 対策、既定は Green-Gauss)。詳細は
  [discretization.md §7.3](discretization.md#実装)。cell モードは GG のみ。
  `gradLSQ=1` (毎ステップ solve) は退化ノードを持つメッシュ (例 case/29 軸対称ノズル近壁) で発散する
  既知の負結果があり回帰対照のみ。**推奨は `gradLSQ=2` (係数事前計算 + スペクトル打ち切りフォールバック**、
  [discretization.md §7.3.1](discretization.md#実装))。
  ([plans/active/discretization-lsq-gradient.md](../plans/active/discretization-lsq-gradient.md) §0/§9)。
- **node のスカラー勾配** ($k,\omega$・化学種 $Y_s$・受動種 $\xi$・凝縮モーメント): 既定は GG
  (境界半割面を owner 値で積算、軸対称は `A_planar`)。NS の LSQ と同じ事前計算係数に揃える opt-in 経路
  `mesh.scalarGradient: lsq` を実装中 ([plans/active/gradient-scalar-lsq-unification.md](../plans/active/gradient-scalar-lsq-unification.md))。
  LSQ 経路は差分形 $\sum_j c_{ij}(\phi_j-\phi_i)$ で定数場が厳密に 0 (GG は壁節点で float32 の閉包誤差約 $5\,\varepsilon|\phi|/h$)、
  境界は NS と同じ内部隣接のみ、並進周期は NS の合併係数を共有、軸対称×周期・回転周期は片側 (gather しない)。
  既定の生産設定で生きているスカラー勾配は $k,\omega$ だけ (`speciesFaceReconstruction` の既定 0)。
- ジッタ格子での LSQ 勾配の最大誤差の次数は 0.8–1.0 (格子ごとに評価点集合と stencil が変わるため)。
  形状固定の縮小では 1 次整合で、継ぎ目の誤差定数は内部以下 (plan boundary-node-periodic-gradient-fix §6.2 #8b)。欠陥でなく作用素の性質。
- **GG は非一様メッシュで線形場非厳密**: `fx` 射影補間の GG 勾配は一様直交では線形場を機械精度で
  再現するが、ノードジッタ/三角形/高 AR メッシュでは O(1) の相対勾配誤差が残る
  (30% ジッタ quad で最大 66%、`tools/verify_linear_recon.py` で定量・全 PASS)。
  勾配を使う面再構成 (MUSCL, `keepDissJump`) の高次性はこの範囲でしか成立しない。
  LSQ (double) は任意メッシュで線形場厳密 (同ツールで確認)。
- **境界面の寄与**: cell-centered ではゴーストセル値を `c_2` として扱うことで内部面ループへ吸収する
  (境界条件モジュールが事前にゴースト値を書く)。node-centered (median-dual) は別の閉じ方をする —
  境界半割面は内部面ループから除外し、専用カーネル `calcGradient_b_d` が **owner ノードの状態値を境界面値**
  として加算する (ゴースト・bvar のどちらも参照しない)。詳細・理由は
  [discretization.md §6.2](discretization.md#62-弱形式境界-weak-form-boundary) / §7.2.2 を参照。

### 参考

- 実装詳細とソース対応は 本ドキュメントの「実装」節 を参照。
- アーキテクチャ全体の中での位置付けは
  [`methods/architecture/overview.md`](architecture/overview.md) を参照。

## 実装

### ソースファイル

| ファイル | 役割 |
| --- | --- |
| [`solver_density_cuda/gradient.hpp`](../solver_density_cuda/gradient.hpp) | CPU エントリ `gradientGauss` の宣言 |
| [`solver_density_cuda/gradient.cpp`](../solver_density_cuda/gradient.cpp) | CPU 実装 + GPU ラッパへの委譲 |
| [`solver_density_cuda/cuda_forge/calcGradient_d.cuh`](../solver_density_cuda/cuda_forge/calcGradient_d.cuh) | GPU カーネル / ラッパ宣言 |
| [`solver_density_cuda/cuda_forge/calcGradient_d.cu`](../solver_density_cuda/cuda_forge/calcGradient_d.cu) | GPU カーネル実装 |

### エントリポイントとディスパッチ

`gradientGauss(cfg, cuda_cfg, msh, v)` が CPU/GPU 共通の入口。
`cfg.gpu == 1` の場合は GPU ラッパ `calcGradient_d_wrapper` に委譲し、
CPU 経路では同関数内で全処理を行う。実運用 (`main.cpp`) では GPU ラッパが
直接呼ばれている (CPU 経路は CPU デバッグ用)。

呼び出し箇所:

- 初期化直後の勾配計算 ([`main.cpp` L624](../solver_density_cuda/main.cpp#L624))
- 時間ステップ内の各サブステップ ([`main.cpp` L697](../solver_density_cuda/main.cpp#L697), [L727](../solver_density_cuda/main.cpp#L727))

### CPU 実装の構造

[`gradient.cpp`](../solver_density_cuda/gradient.cpp) の処理は次の 3 段。

1. **ゼロクリア**: 全セル (`nCells_all`、ゴースト含む) について
   $\nabla U_{x,y,z}, \nabla P, \nabla T$ を 0 に初期化。
2. **面ループ加算**: `for ip in [0, nPlanes)`
   - 両側セル `ic1 = iCells[0]`, `ic2 = iCells[1]` と面ベクトル `surfVect` を取得。
   - 重み `f = fxp[ip]` で線形補間: `Uxf = f*Ux[ic1] + (1-f)*Ux[ic2]` ほか。
   - `dXdx[ic1] += sv[0]*Xf`、`dXdx[ic2] -= sv[0]*Xf` (符号反転で対称加算)。
3. **体積除算**: `for ic in [0, nCells)` でセル体積 `msh.cells[ic].volume` で除する。

CPU 経路では密度勾配 $\nabla \rho$ は計算していない (上位で必要としていないため)。

### GPU 実装の構造

[`calcGradient_d_wrapper`](../solver_density_cuda/cuda_forge/calcGradient_d.cu#L384) が次のカーネルを順に発行する。

1. **`cudaMemset` 群**: $\nabla U_*, \nabla \rho, \nabla P, \nabla T,
   \nabla \cdot \mathbf{u}$ をデバイスバッファ上で 0 クリア。
2. **`calcGradient_1_d`**: 面 1 並列。`atomicAdd` で
   両側セルの勾配ベクトル成分と $\nabla \cdot \mathbf{u}$ に寄与を累積。
   `plane_cells[2*ip+0/1]` でセル ID を取得し、`fx[ip]` で線形補間。
3. **`calcGradient_2_d`**: セル 1 並列。`vol[ic]` で除算して勾配を確定。

CPU 版と異なり、密度 $\rho$ の勾配 $\nabla \rho$ も同時に計算する。
$\rho U_*, \rho e, H_t$ の勾配計算用コードは保留 (コメントアウト済み)。

#### 境界寄与

**cell-centered**: 境界面の寄与はゴーストセルを内部面ループの `ic1` として処理する設計
(ゴーストセル値は境界条件モジュールが事前に書き込む)。

**node-centered**: `calcGradient_cellgather_d` が境界半割面 (`ip>=nNormalPlanes`) を skip し、
[`calcGradient_d.cu`](../solver_density_cuda/cuda_forge/calcGradient_d.cu) の `calcGradient_b_d` が
非 periodic の全 bcond について **owner ノードの状態値**を境界面値として加算する (bvar は参照しない)。
periodic は DOF 同一視・gradient gather (§discretization.md §4.5) に委ね寄与を加えない。

**node × 並進周期の継ぎ目** (plan [`boundary-node-periodic-gradient-fix.md`](../plans/accepted/boundary-node-periodic-gradient-fix.md))。
適用条件は `periodicSeamMergeActive` (node ∧ 周期 group あり ∧ 非軸対称 ∧ 周期 bcond がすべて並進 `type: 0`)。以下はすべてこの条件で有効になる。

- **LSQ (`gradLSQ: 2`、NS の原始量)**: 事前計算で継ぎ目越しの**合併 stencil** を組む (`calcGradient_d.cu` の `lsqPre_mergePeriodic`)。
  group の全部分 CV の incidence を集め、同じ物理隣接 (隣接の `periodicRoot` が一致し、変位が $10^{-4}h_{min}$ 以内) を同値類にまとめて
  配分係数 $\alpha=1/\text{重複数}$ を決める。行列と係数は**各 incidence の実変位** $d_{mj}$ で組む:
  $M_r=\sum\alpha_{mj}w_{mj}d_{mj}d_{mj}^{\mathsf T}$ ($w=1/|d|^2$)、スペクトル打ち切りは $M_r$ に 1 回、係数 $c_{mj}=M_{r,\tau}^{+}\alpha_{mj}w_{mj}d_{mj}$。
  毎 step の `periodicGradientGather` は和のまま (部分和が合併 LSQ になる)。
- **Green–Gauss (SST の $k,\omega$、化学種、受動種・凝縮モーメント)**: 各部分 CV は**周期半割面を積算しない** (面フラグ `mesh.planePeriodic_d`)。
  合併体積で割った部分寄与を和で合併する。$k,\omega$ は `ransGradient` の直後 (`ransBlendF1` の前) に専用の gather、
  化学種・受動種は `periodicGradientGather` に登録された gather を使う。化学種・受動種は `species_gradient_d` の同じ呼び出し経路。
- **既知の制約**: 壁の CV は壁半割面を φ[ic0] で積算するので、float32 の格納面ベクトルでは定数場の GG が閉じず、壁節点で約 $5\,\varepsilon|\phi|/h$ の偽勾配が出る (2026-09-26 channel 実測、継ぎ目に依らない)。LSQ 化 (後続 plan) で解消する見込み。
- **この条件の外** (軸対称×周期、回転周期): スカラー勾配は片側 GG + 半割面込みのまま (既存の未修正挙動で、本修正では不変)。
  回転周期は plan `boundary-node-rotational-periodic` で扱う。

**SST の F1 の初期値** (周期に依らず全 SST run): `sstF1` は配列確保時 (`variables.cpp` の `allocVariables`) に 1 で初期化し、`buildScalarDescs` は副作用を持たない。
2026-09-26 までは `buildScalarDescs` の初回呼び出しが計算済みの F1 を 1 で上書きしていたので、`sstSigmaBlend: 1` (既定) の SST run は**非周期でも初回 step の k/ω が変わる**
(case/48 で step 1 の roK/roOmega が約 5 万点変化、`sstSigmaBlend: 0` で消えることを確認)。restart 直後の 1 step だけ k の残差が跳ねる現象 (case/39 旧バイナリで床の 408 倍) もこれが原因
(`1266aba1` + この修正だけの版で 0.99 倍に消えた)。

**履歴 (2026-09-26 まで)**: LSQ は各部分 CV が片側の隣接だけで完全な勾配を解き、和が線形場で正確に 2 倍になっていた (case/09 TGV 実測)。
SST の $k,\omega$ 勾配は gather の後に `ransGradient` が作り直して合算されていなかった。GG の周期半割面の除外条件 `ic1 < nCells` は、
周期 bcond にもゴーストが付く (`mesh.cpp`) ため一度も成立せず、継ぎ目に $\phi(S_a+S_b)/V$ の誤差があった。

理論・設計判断は [discretization.md §6.2/§7.2.2](discretization.md#62-弱形式境界-weak-form-boundary) を参照。

### 並列化メモ

- 面 1 並列の累積は `atomicAdd` を使用。`flow_float` が `double` の場合は
  Compute Capability 6.0 以降が必要。
- `__syncthreads()` 呼び出しは入退出の安全策 (ループ内同期はない)。
- `dimGrid_plane`, `dimGrid_cell` は `cudaConfig` で `nPlanes`, `nCells` から計算。

### 入出力データ

入力 (セル中心配列): `ro`, `Ux`, `Uy`, `Uz`, `P`, `T`, `roe`, `Ht`, セル体積、
面ベクトル (`sx, sy, sz, ss`)、補間重み `fx`、`plane_cells` 接続、
GPU 側のミラー `var.c_d[*]`, `var.p_d[*]`。

出力 (セル中心配列):

```
dUxdx, dUxdy, dUxdz
dUydx, dUydy, dUydz
dUzdx, dUzdy, dUzdz
dPdx,  dPdy,  dPdz
dTdx,  dTdy,  dTdz
drodx, drody, drodz   (GPU only)
divU
```

### 既知の TODO / 注意点

- limiter / 高次再構成は未実装。コード内では `gradientGauss` の単一経路のみ。
- $\rho U_*, \rho e, H_t$ の勾配計算用変数とカーネル分岐は `cu` 内に
  コメントとして残されている。必要時に有効化する。
- 境界カーネル `calcGradient_b_d` は無効化中。境界条件設計と合わせて再評価する。
