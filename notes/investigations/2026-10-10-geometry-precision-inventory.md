# 幾何の座標の使われ方と精度の洗い出し (2026-10-10)

plan [`architecture-float-state-double-geometry`](../../plans/active/architecture-float-state-double-geometry.md) の前提調査。調査用エージェント (Explore) が読んだ結果を、呼び出し側でまとめ直した。行番号は `fd2d4c3c` 時点。根は `solver_density_cuda/`。

## 0. 前提を縛る 2 つの事実

1. **`geom_float = double` だけにはできない。**
   - 幾何の配列はすべて `flow_float` の表 `variables::c/p/c_d/p_d` に入っている (`variables.hpp:22-25`)。
     - セルの量: `volume, ccx..ccz, A_planar, A_closure_*, delta_les` (`variables.hpp:94-99, 213`)
     - 面の量: `sx..ss, *_planar, pcx..pcz, fx, dcc` (`variables.hpp:227-232`)
     - 確保は `variables.cpp:357-397`。
   - カーネルは `geom_float*` として受け、`flow_float*` を渡されている (例 `calcStructualVariables_d.cu:76-79`)。
   - ホスト側の写しは `sizeof(geom_float)` で `flow_float` のバッファへコピーしている (`variables.cpp:676-696, 718`)。
   - 流束のカーネルの多くは、座標を `flow_float` のローカル変数に写してから差を取る (`convectiveFlux_roe_d.inc.cuh:101-111`、`viscousFlux_d.cu:113-119`)。
2. **既定のビルドでは、幾何は double で作られていない。**
   - 変換器も同じ typedef を使い、座標は `geom_float` が double でなければ `stof` で読む (`gmshReader.hpp:505-514`)。
   - `vector<geom_float>` で書くので、HDF5 も float32 になる。
   - FP64 の変換器で書いたファイルでも、ソルバは `std::vector<geom_float>` に読むので float に丸まる (`mesh.cpp:236-260, 304-307`、`variables.cpp:752-753`)。

略号:
- **cc-cc**: 節点 (node 離散化では中心 = 節点座標、`mesh.cpp:356-375`) の間のベクトル。
- **pc-cc**: 面重心 − 中心。
- **f32**: float で差を取る。
- **f64(q)**: double で差を取るが、元の座標は float に丸めた後。

## 1. 使う側ごとの表

### 1a. 読み込み時の構造量

| 場所 | 量 | 精度 |
| --- | --- | --- |
| `variables.cpp:565-587` | sx..ss・pcx..pcz・ccx..ccz・volume の写し | 差なし |
| `variables.cpp:597-629` | 軸対称の r 重み (面は pcy、セルは rEff/ccy)、`A_planar` | 絶対値の積 (問題なし) |
| `variables.cpp:634-651` | `A_closure_x/y = Σ±S_f·r_f` | **打ち消し**: 薄い半径方向のセルで S·r_top − S·r_bot ≈ Δx·Δr。r は float の pcy |
| `variables.cpp:705-720` | `delta_les` = max\|cc_j − cc_i\| | cc-cc・f32 (DES のみ) |
| `calcStructualVariables_d.cu:28-61` | `dcc` (どのカーネルも読まない)・`fx` | cc-cc・pc-cc・f32。node の内部面は fx = 0.5 |
| `mesh.cpp:467-493` | ゴーストの中心 = xc + 2((xp−xc)·n)n (cell 離散化) | pc-cc・f32 |

### 1b. 勾配

- Green-Gauss (`calcGradient_d.cu:189-450`) は座標を使わない (fx・sx・vol だけ)。
- LSQ (`gradLSQ=1`、`calcGradient_d.cu:509-552`): d = cc_j − cc_i、w = 1/\|d\|² は **f64(q)**。M・b は `flow_float` に格納し、解くのは double。
- LSQ の前処理 (`gradLSQ=2`、`calcGradient_d.cu:699-834`): 係数 c_ij = M⁺ w d は **f64(q)** で作り、`flow_float` の `cInt` に入れる。
- 周期の継ぎ目の合わせ込み (`calcGradient_d.cu:733-807`): **f64(q)**。許容差は 1e-4·h_min。

### 1c. MUSCL / 辺中点の再構成

node の方式では常に辺中点 (`convectiveFlux_d.cu:218-219`) で、dc0p = ±0.5·dcc。

- Roe・HLLE・SLAU (`convectiveFlux_{roe,hlle,slau}_d.inc.cuh`) は、dcc = cc1 − cc0 と pc − cc を f32 のローカル変数で作る。SLAU は化学種・受動種の面の再構成でも同じ。
- KEEP (`convectiveFlux_keep_d.inc.cuh:235-236, 417-418`) と、`#include` されている legacy の 2 本も同じ。

### 1d. リミタ

- `limiter_d.cu:220-222, 314-416`、`limiterPeriodic_d.cuh:97-165`、`passiveLimiter_d.cuh:54-56`: pc − cc0、または辺中点の 0.5(cc1 − cc0)。f32。

### 1e. 粘性の流束 (壁際の dcc に最も敏感)

- `viscousFlux_d.cu:113-124, 151-158`: dcc = cc1 − cc0 (f32)、D = dcc·S、δ = dcc·\|S\|²/D。
  - SU2 の補正・τ の法線項・熱流束・k のエネルギー項で使う。
  - 壁節点と内部節点の辺もここを通る。
- `viscousFlux_wall_d` (`:502-670`): cell 離散化はゴーストとの dcc (f32)。node は ∇U·S。
- 弱形式の等温壁の d1 (`conjugateWall.cpp:46-130`): **f64(q)** で作って `flow_float` に入れる。CHT も同じ d1/d2 を使う。

### 1f. スカラー・化学種・受動種の拡散

すべて cc-cc・f32・δ/dcc。

- `scalarTransport_d.cu:105-111, 202-208`: k/ω・トレーサ・遷移・凝縮モーメント。`transport_diag` の `geo = |δ|/dcc` もここで作る。
- `speciesTransport_d.cu:283-289, 2191-2199, 2808-2814`
- `passiveKernels_d.cuh:184-188`、`passiveFct_d.cuh:28-32`

### 1g. 局所時間刻み

- `setDT_d.cu`: dx = vol/ss。**座標の差は使わない**。

### 1h. 陰解法

- `timeIntegration_d.cu:587-602` (スカラー DPLUR): cc-cc・f32 で δ、粘性の対角。
- `timeIntegration_d.cu:952-963` (block DPLUR、`<ST>`): `static_cast<ST>(ccx[j]) - static_cast<ST>(ccx[i])`。`implicitSolvePrecision=1` で double だが f64(q)。
- `timeIntegration_d.cu:1327-1335` (前処理版): f32。

### 1i. ラインの構築

- `mesh.cpp:1045-1123`: 節点座標を double に写して dist2 と整列の余弦を計算する (f64(q))。トポロジーだけを決める。

### 1j. 壁距離と壁関数

- `wall_dist` は変換器だけで作る (`input/calcWallDistance_kdtree.cpp:133-174`、geom_float)。ソルバは絶対値として読む (`variables.cpp:752-782`)。
- node の壁関数の代表点 (`ransWallFunction_d.cu:146-166`、`wmlesWallModel_d.cu:158-175`): d = cc_I − cc_W、y = −d·n (f32)。

### 1k. 軸対称のソース

- r の絶対値だけを使う (float で問題なし)。打ち消しは 1a の closure だけ。

### 1l・1m. 境界条件・周期

- 境界条件のカーネルは sx/ss/fx だけ。
- 周期の相手の合わせ込み (`mesh.cpp:549-660`) は面重心で、geom_float。

## 2. ファイルから読む量と作る量

| 量 | 出どころ | 読み込み | ソルバでの型 |
| --- | --- | --- | --- |
| 節点座標 | `/MESH/COORD` | `mesh.cpp:236-249` | geom_float |
| 面ベクトル・面積・面重心 | `/PLANES/*` | `mesh.cpp:254-294` | geom_float → flow_float (device) |
| 体積・セル中心 | `/CELLS/*` | `mesh.cpp:304-351` | 同 |
| wall_dist | `/VALUE/wall_dist` (変換器) | `variables.cpp:752-782` | flow_float |
| node の中心・rEff | 読み込み時に作る | `mesh.cpp:356-375` | geom_float |
| r 重みの S/vol・A_planar・A_closure | 読み込み時に作る | `variables.cpp:589-658` | flow_float |
| LSQ の係数 | 読み込み時のカーネル | `calcGradient_d.cu:1053-1088` | flow_float |
| 弱形式の等温壁の d1/d2 | 読み込み時 (double) | `conjugateWall.cpp:76-156` | flow_float |
| dcc・pc − cc・δ・壁関数の y | **毎 step カーネルで作る** | §1 | f32 |

- ソルバは `/DUAL/*` を読まない。変換器が双対の幾何を `/PLANES`・`/CELLS` に書いている (`gmshReader.hpp:2222-2266`)。
- 2D の node の双対は、変換器の中では辺の中点を原点にした double で作る。ただし、体積・重心の和は geom_float に貯め、頂点に使う一次の重心も丸めた後のもの。

## 3. `qAccumulatorFP64` と軸対称

- 実装: `variables.hpp:27-35`、`cuda_forge/qAccumulator_d.cuh`、`update_d.cu:36-56, 173-186, 272-298`。設定 `input/solverConfig.hpp:101-106`。
- 軸対称を拒否するのは `main.cpp:3371-3374`。
  - 理由: `enforceAxisSymmetry` が commit の基準 `roeN`/`roUyN` を射影するが (`axisymmetricSource_d.cu:312-323`)、Qacc は Q_N を読まないので、その射影が Qacc に効かない (S1b-① 未対応)。
  - 幾何の精度とは関係ない。対象は流れの保存量 5 つだけ。

## 4. 既存の精度の切り替え

- typedef (全部か無しか): FP64 のビルドは約 2 倍遅い。幾何だけを double にする設定は無い。
- `implicitSolvePrecision`: block DPLUR の行列と 5×5 の解を double にする。dcc は f64(q) のまま。
- ライン陰解法の Thomas、`gradLSQ=2` の前処理、`qAccumulatorFP64` は、どれも幾何の精度を上げない。

## 5. 座標の丸めに最も弱いところ (直す順)

1. 粘性・拡散の辺の δ/dcc (特に壁節点と内部節点の辺)
2. LSQ の係数
3. node の壁関数の y
4. 変換器の `wall_dist` と、cell 離散化の一次の幾何
5. 弱形式の等温壁の d1
6. 軸対称の closure Σ S·r
7. 影響が小さいもの: MUSCL の辺ベクトルとリミタ (勾配の増分を掛けるだけ)、DPLUR の粘性の対角 (収束にだけ効き、収束解には効かない)
