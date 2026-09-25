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

**node × 周期の継ぎ目 (既知の欠陥、修正中: plan [`boundary-node-periodic-gradient-fix.md`](../plans/active/boundary-node-periodic-gradient-fix.md))**:
継ぎ目で割れた CV の部分勾配は `periodicGradientGather` で**和**を取る。Green–Gauss (合併体積で割った部分寄与) なら和が合併勾配になるが、
node の既定の **LSQ (`gradLSQ: 2`) では各部分 CV が片側の隣接だけで完全な勾配を解くので、和は線形場で正確に 2 倍**になる
(2026-09-26 実測、case/09 TGV)。修正は、事前計算で継ぎ目越しの**合併 stencil** の LSQ を組むこと:
重み $w=1/|\Delta\mathbf x|^2$、同じ物理隣接が両側に現れる継ぎ目接線エッジは配分係数 $\alpha$ (同一隣接で総和 1) で重複を除き、
$M_r=\sum\alpha w\,\Delta\mathbf x\Delta\mathbf x^{\mathsf T}$、スペクトル打ち切りは $M_r$ に 1 回、各部分 CV の係数 $c=M_r^{+}\alpha w\Delta\mathbf x$ を焼き込む。
毎 step の gather は和のまま (部分和が合併 LSQ になる)。化学種・受動種・凝縮モーメントの勾配は Green–Gauss なので現行の和が正しい。
SST の $k,\omega$ 勾配は `ransGradient` が合算の後に作り直していて合算されていない (同 plan で修正)。
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
