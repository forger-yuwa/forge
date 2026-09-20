# リミッタ

MUSCL 再構成 (2 / 3 次) は不連続近傍で振動を生む。forge ではセル中心勾配にスカラ係数
$\psi \in [0, 1]$ を掛けて再構成を抑制する。

$$
\phi_L = \phi_C + \psi_C \,\nabla \phi_C \cdot (\mathbf{r}_f - \mathbf{r}_C).
$$

本ドキュメントは理論(係数の定義とスキーム)と実装(ソース対応)をまとめる。

## 理論

### 評価対象

セルごとに $\rho, U_x, U_y, U_z, P$ の 5 成分について独立にリミッタ係数を持つ。
温度 $T$ と全エンタルピ $H_t$ は派生量として扱い、独立リミッタは持たない。

### $\delta^+, \delta^-, \delta_m$ の定義

セル $C$ について隣接セル値の振れ幅を

$$
\delta^+_C = \phi_C^{\max} - \phi_C ,\qquad
\delta^-_C = \phi_C^{\min} - \phi_C
$$

($\phi_C^{\max/\min}$ は近傍セル値の最大・最小)。
面 $f$ への線形再構成による予測差分を

$$
\delta_m = \tilde\phi_f - \phi_C
= \nabla \phi_C \cdot (\mathbf r_f - \mathbf r_C)
$$

とする。$\delta_m > 0$ なら $\delta^+$ を、$\delta_m < 0$ なら $\delta^-$ を比較対象に取る。

### Barth–Jespersen リミッタ

$$
\psi_f = \min\!\left(1,\, \frac{\delta^\pm}{\delta_m}\right),
\qquad \psi_C = \min_{f \in \partial C} \psi_f .
$$

単純で TVD 性が強く、滑らかな領域でも 2 次精度の頭打ちを招きやすい。

### Venkatakrishnan リミッタ

Barth–Jespersen を滑らか化したもの。$\epsilon^2 = (K |V_C|^{1/3})^3$ (`K=1.0`) として

$$
\psi(\delta^+, \delta_m) =
\frac{1}{\delta_m}\,
\frac{(\delta^{+\,2} + \epsilon^2)\,\delta_m + 2\delta_m^2 \delta^+}
     {\delta^{+\,2} + 2\delta_m^2 + \delta^+ \delta_m + \epsilon^2}.
$$

$\delta_m \to 0$ で $\psi \to 1$ (連続)。$|\delta_m|$ が体積スケール以下では実質
無リミットになり、滑らかな領域での精度低下を避ける。

#### ⚠ 既定の $\epsilon^2$ は**次元が合っていない**

$\epsilon^2 = (K|V_C|^{1/3})^3 = K^3 |V_C|$ は**長さの 3 乗の次元**を持つのに、比べる相手の
$\delta$ は変数そのものの次元 ($\rho$ なら kg/m³、$P$ なら Pa) である。結果:

- **メッシュを拡大すると実質 OFF になる**。座標を ×1024 すると $\epsilon^2$ が面積とともに ×2²⁰ 増え、
  分子分母を支配して常に $\psi \approx 1$ を返す (実測: $\bar\psi_\rho$ 0.9718 → 0.99998、
  62069 ノード中 58761 で $\psi$ が変化)。
- **変数ごとに効き方が桁違いになる**。同じ $\epsilon^2$ を $\Delta\rho \sim O(0.1)$ と
  $\Delta P \sim O(10^5)$ に当てるので、片方では無リミット・片方では通常動作になる。

SU2 が同じ形の式で壊れないのは、**解を無次元化して解いている**ため $\delta$ が $O(1)$ だからである
(`Common/src/CConfig.cpp:5021` `RefElemLength = 1.0`, `VENKAT_LIMITER_COEFF = 0.05` → $\epsilon^2 = 1.25\times10^{-4}$)。
forge は SI 次元のまま解くので同じ式が成立しない。

#### 無次元化 Venkatakrishnan (`space.limiterScaled: 1`, opt-in)

$\delta$ を**変数ごとの固定参照** $q_\mathrm{ref}$ で割ってから Venkatakrishnan を当てる。

$$
\hat\delta = \delta / q_\mathrm{ref}, \qquad
\hat\epsilon^2 = \left(\frac{K\,h_i}{L_\mathrm{ref}}\right)^3 ,
$$

- $q_\mathrm{ref}$: $\rho$ は `limiterRoRef`、$P$ は `limiterPRef`、速度 3 成分は `limiterARef` (音速)。
  **明示指定 (>0) が最優先**で、0 (既定) なら起動時に初期場の体積加重平均から決める。
  自動だと restart のたびに値が変わるのでログに `(auto)` と警告を出す。貼り付け用の行も印字する。
- $h_i$: **軸対称は `A_planar`、平面 2D は `volume` (奥行 1 の面積) の平方根**、3D は $V^{1/3}$。
  `A_planar` は軸対称のときしか device へ転送されないので、**平面 2D で `A_planar` を読んではいけない**
  (読むと $\epsilon^2=0$ になり `venkatK` が効かなくなる)。
- $L_\mathrm{ref}$: `limiterRefLength`。0 ならメッシュ境界箱の対角。

**`space.venkatK` の推奨は 0.05** (SU2 既定と同値)。既定の 1.0 は Sod で全変数を悪化させる
(密度の近傍逸脱が K=1.0 で 86674 面側、K=0.05 で **0**)。SERN でも K=1.0 は $\rho$ 797 / $U_y$ 1777 に対し
K=0.05 は $\rho$/$U_x$/$P$ が 0。

**`space.limiterScaled > 0` は `space.limiterMatchRecon: 1` を要求する** (通常経路は `matchRecon==0` で
`limiterScaled` を無視するのに周期経路は適用するため、契約が割れる)。

#### 比の形 (`space.limiterScaled: 2`) — **棄却**

$\psi$ を $y=\delta^+/\delta_m$ だけの関数にする形。基準値も長さも要らないが、
**滑らか域でリミッタが切れない**ので定常残差の床が 2〜5 倍上がる (case/44 で実測)。使わないこと。

### リミッタの評価点 (`space.limiterMatchRecon`)

リミッタが $\delta_m$ を評価する点は、**流束が再構成する点と一致していなければならない**。

| 設定 | 評価点 | 増分の形 |
| --- | --- | --- |
| `0` (既定) | 双対面重心 | $g\cdot d$ のみ |
| `1` | **流束と同じ点** (node は常にエッジ中点) | `convMethod` と同じ (2 なら隣接値差の項も含む) |

node の流束は `matchRecon` に関係なく**常にエッジ中点**で再構成する (`convectiveFlux_d.cu` の
`g_reconEdgeMid` は node なら無条件 1) ので、既定の `0` では**リミッタと流束が別の点を見ている**。
実測ではこれが $U_x$ と $P$ の近傍逸脱の原因で、`1` にすると 0 になる
(密度と $U_y$ の逸脱は次元不整合が原因なので `limiterScaled` 側で直す)。

### 有界性の診断 (`space.limiterDiag: N`, 既定 0)

$\psi$ 確定後に**流束と同じ点・同じ増分関数**で face 値を作り直し、そのノードの近傍 min/max を
外れた face-side を変数ごとに数える。N 回の flux 呼び出しごとに印字し、**窓と累計 (`CUM`) の両方**を出す。

- 許容幅は $\max(\text{近傍レンジ},\ |q|_\max) \times 10^{-5}$ に**絶対床 $10^{-6} q_\mathrm{ref}$** を敷く
  (恒等 0 の変数では近傍レンジも 0 になり float ノイズを逸脱と数えてしまうため)。
- 非有限は別枠で数える (NaN は大小比較が両方 false になるので「逸脱なし」に化ける)。
- 検査件数 0 は `VERDICT=NO-CHECKS(FAIL)`。幾何尺度 $h_i \le 0$ も警告する。
- **後処理で再構成を再現して判定してはいけない** — 出力ファイルは状態が更新後・勾配が更新前で 1 step ずれる。

### 面単位の非物理フォールバック (`space.badReconDiag` / `badReconFallback`, 既定 0)

SU2 の `bad_recon` 相当。流束が消費する L/R 状態が $\rho \le$ `roMin` / $P \le$ `pMin` / 非有限なら、
**その面だけ**両側をセル値に戻し (`conv_scheme = -1`)、N 回の訪問だけ 1 次に保つ。
**⚠ カウンタはフォールバック前に数える**ので、発火数から介入の効果は読めない。
SERN のベース発散には効かなかった (opt-in で残置)。

### Nishikawa R1 リミッタ (未有効)

`nishikawa_r1_limiter` の実装は残されているがコメントアウト済み。
$\delta'_C$ と幾何比 $r_{ik}$ から不連続適応する形式。

### Ducros センサとの併用

[`methods/convection/`](convection/) で記述したように、
Roe カーネルでは Ducros センサ値 $\mathrm{duc}$ で

$$
\psi \leftarrow \begin{cases}
\psi & (\mathrm{duc} \le 0.8) \\
\max(0,\, (1 - \mathrm{duc})\,\psi) & (\mathrm{duc} > 0.8)
\end{cases}
$$

と二段がけする。

### 全リミッタ 1 の場合 (`cfg.limiter == 0`)

リミッタを評価せず $\psi = 1$ を一律で配るモード。低次再構成 (`convMethod = 0`)
や smooth な計算ですべての面で線形再構成を使いたい場合に選ぶ。

## 実装

forge のスカラ・リミッタ計算は GPU でのみ動作する (CPU 経路は無い)。

### ソースファイル

| ファイル | 役割 |
| --- | --- |
| [`solver_density_cuda/cuda_forge/limiter_d.cuh`](../solver_density_cuda/cuda_forge/limiter_d.cuh) | カーネル / ラッパ宣言、`deltas` 構造体 |
| [`solver_density_cuda/cuda_forge/limiter_d.cu`](../solver_density_cuda/cuda_forge/limiter_d.cu) | リミッタ device 関数とセル並列カーネル (約 325 行) |

### エントリポイント

[`limiter_d_wrapper`](../solver_density_cuda/cuda_forge/limiter_d.cu#L198)
が共通入口。次の順で動く。

1. `fill_limiter_d` で `limiter_ro, limiter_Ux, limiter_Uy, limiter_Uz, limiter_P` を
   `1.0` に初期化 (全セル `nCells_all`)。
2. `cfg.limiter == 0` なら早期 return (全リミッタ 1)。
3. 各成分について `limiter_r1_d` をセル並列で呼び、近傍探索を行って $\psi$ を確定。

呼び出しは時間積分ループの `gradient → limiter → convectiveFlux` の中間で発生する
([`main.cpp` L700](../solver_density_cuda/main.cpp#L700), [L730](../solver_density_cuda/main.cpp#L730))。

### `limiter_r1_d` 構造 ([L75](../solver_density_cuda/cuda_forge/limiter_d.cu#L75))

セル 1 並列。`ic0 < nCells` で動作する。

1. `limiter_scheme == 0` なら `limiter_Q[ic0] = 1.0` で抜ける。
2. `cell_planes_index[ic0]`, `cell_planes[*]` でセルが接する面を列挙。
3. 各内部面 (`ip < nNormalPlanes`) について反対側セル `ic1` を取得 (`plane_cells[2*ip+0] + plane_cells[2*ip+1] - ic0`)、
   `Q_max, Q_min` を更新。
4. もう一度ループして面ごとに $\delta_m$ を求め、`limiter_function` で $\psi_f$ を計算、
   最小値を取って `limiter_Q[ic0]` に書く。
5. 結果は `clamp(ψ, 0, 1)` で安全側に丸める。

#### `limiter_function` の選択

`limiter_scheme` (= `cfg.limiter`) で device 関数ポインタを差し替える。

| `limiter_scheme` | 関数 | 種別 |
| --- | --- | --- |
| `0` | — (早期 return) | 無リミッタ ($\psi = 1$) |
| `1` | `barth_Jespersen_limiter` ([L39](../solver_density_cuda/cuda_forge/limiter_d.cu#L39)) | TVD 強 |
| `2`, `-1` | `venkata_limiter` ([L17](../solver_density_cuda/cuda_forge/limiter_d.cu#L17)) | 滑らか版 |
| `3` (TODO) | `nishikawa_r1_limiter` | 実装あり・コメントアウト |

`venkata_limiter` の $\epsilon^2$ は `K = 1.0` 固定 (`K^3 * volume`)。
`|δ_m| < 1e-20` のとき $\psi = 1$ にフォールバック。

### 適用先 (面ごとの参照)

各スキームカーネル ([`convectiveFlux_d.cu`](../solver_density_cuda/cuda_forge/convectiveFlux_d.cu))
は L 側に `limiter_ro[ic0]` を、R 側に `limiter_ro[ic1]` を `interp_dispatch` に渡す
(他成分も同様)。Roe カーネルでは追加で Ducros 補正
[`apply_ducros_limiter`](../solver_density_cuda/cuda_forge/convectiveFlux_d.cu#L189)
を経由する (SLAU, HLLE 経路は Ducros 補正を直接通さず素のリミッタを使う)。

### 入出力

入力: セル中心保存量・原始量のうち $\rho, U_x, U_y, U_z, P$、対応する勾配
(`drod*, dUxd*, dUyd*, dUzd*, dPd*`)、`cell_planes_index`, `cell_planes`,
`plane_cells`, `pcx/y/z`, `ccx/y/z`, `vol`, `fx`。

出力: `limiter_ro, limiter_Ux, limiter_Uy, limiter_Uz, limiter_P` (各セル 1 値, `[0, 1]`)。

### 既知の TODO / 注意点

- `Ht` 用のリミッタ呼び出しはコメントアウト済み (現状 $H_t$ は再構成しない)。
- Nishikawa R1 は実装ありで未有効。$\delta'_C$ と $r_{ik}$ の計算を別関数 `calcDeltaIJ_dash`
  ([L19-L31](../solver_density_cuda/cuda_forge/limiter_d.cuh#L19)) で行う。
- リミッタ計算は近傍ループ 2 周回 (max/min 取得 → ψ 評価) でレジスタを多く使う。
  `dimGrid_normalcell_small`, `dimBlock_small` のブロック構成は `cudaConfig` で
  個別に調整されている。
- 全境界周期 / 滑らかな計算 (Taylor–Green 等) では `cfg.limiter = 0` (無リミッタ) が有効。

## 関連

- リミッタの参照側 (`interp_dispatch`): [`methods/convection/`](convection/)。
- 勾配の生成: [`methods/gradient/`](gradient/)。
