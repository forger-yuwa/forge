# node の対流再構成: リミッタ評価点の不一致と非物理再構成のフォールバック

## メタ

- **area**: `convection / boundary`
- **status**: `draft`  <!-- 2026-09-19 起票。case/46 R4e の原因調査 (§4.15.6/§4.15.7) から分離。実装着手前 -->
- **related_docs**:
  - [`methods/boundary.md`](../../methods/boundary.md) 「node 等温壁の壁ノード温度ピン」節 (運動量 Dirichlet の 3 点セット)
  - [`methods/convection/`](../../methods/convection/) (MUSCL 再構成)
  - [`methods/limiter.md`](../../methods/limiter.md)
- **related_plans**:
  - 発見元: [`tooling-nozzle-sern-3d.md`](tooling-nozzle-sern-3d.md) §4.15.6 / §4.15.7 (機体ベースの発散)
  - 隣接: [`discretization-node-boundary-ghostless.md`](discretization-node-boundary-ghostless.md) (壁半割面を消す設計。本件はその**内部面**側)
  - 隣接: [`architecture-node-centroid-value-position.md`](architecture-node-centroid-value-position.md) (値の位置と重心のズレ)
- **created**: `2026-09-19`
- **owner**: `sano`

## 1. 目的

node (median-dual) の**対流再構成まわりで SU2 と食い違っている点**を直す。
起票のきっかけは case/46 の機体ベース発散だったが、**調査の結果それは本 plan の対象ではなくなった** (§4.5/§4.6)。
現在の主題は (a) **リミッタの評価点が再構成点とずれている**こと、(b) **再構成・更新の汎用セーフガードが無い**ことの 2 つ。

## 2. スコープ

**対象**: リミッタの評価点、再構成値の物理性検査、更新率の上限、非物理状態のロールバック。
**対象外**: 粘性流束・勾配計算そのもの・壁関数 (`wallTreatmentSST`)・cell モード。
`mesh.bndFirstOrder` は**使用禁止**なので代替案に含めない (粘性応力まで壊し、疑似 2D で全域に効く)。

## 3. 関連 docs と前提

- 壁ノードの Dirichlet は 3 点セット (状態ピン / 残差射影 / 陰解法 Jacobian 行) で
  [`methods/boundary.md`](../../methods/boundary.md) に記載済み。**本件はその 3 点とは別の層 (対流流束の再構成)**。
- 壁半割面の質量流束が 0 であることは [`discretization-node-boundary-ghostless.md`](discretization-node-boundary-ghostless.md) §4 の前提そのもの。

## 4. 設計方針

### 4.1 現行の実装 (実測で確認済み)

壁ノードまわりは 4 層に分かれている。

| 層 | 実装 | 壁ノードでの挙動 |
| --- | --- | --- |
| 状態 | `cuda_forge/nodeWallDirichlet_d.cu` `enforceWallNoSlip_d` | `roU*` と `U*` を 0 に、KE を `roe` から除去 |
| 残差 | 同 `zeroWallDirichletResiduals_d` | `res_roUx/roUy/roUz` (+SST `res_roOmega`、等温 `res_roe`) を 0 に。**`res_ro` は射影しない** |
| 境界半割面 | `convection/convectiveFlux_boundary_d.inc.cuh` | `mdot = sss*(ro_R*Vn_m)` で bvar 状態のみ・**再構成なし**。no-slip なら `Vn_m = 0` で**質量流束は厳密 0** |
| 内部面 | `convection/convectiveFlux_common_d.cuh` `interp_MUSCL_2nd` | `phif = phiC + psi*(grad . cpd)`。node では `cpd` = エッジ中点 (`convectiveFlux_d.cu` の `rem = (discretization=="node") ? 1 : 0`、**config で切れない**) |

**噛み合っていない**: 状態層が `phiC = 0` にしても、内部面の再構成は勾配を足す。
壁のすぐ外が高速なら勾配は大きく、面に有限の速度が乗る。
境界半割面は塞がっているので、内部面の正味が流出ならその CV は補充されずに空になる。

**1 次ではこの問題が起きない**: `phif = phiC = 0` なので、風上が壁ノード側のとき `rho u = 0` となり、
**壁 CV は対流で質量を失えない** (流れが入る向きなら隣接側が風上になるので、入ることはある)。
2 次はこの性質を壊す。

**forge には再構成値の物理性検査もクランプも無い** (確認済み: `convectiveFlux_slau_d.inc.cuh` で
`ro_L`/`P_L` は `interp_dispatch` の結果をそのまま使う。`contactBlend` は組成勾配向けの別機構で既定 off)。

### 4.2 症状と証拠 (case/46, 2 次元平面, node, SST, frozen_tp)

正本は [`tooling-nozzle-sern-3d.md`](tooling-nozzle-sern-3d.md) §4.15.6 / §4.15.7。要点:

- 機体ベース (厚さ 0.02 H、流れに正対する no-slip 壁) の壁ノードで、再構成が面に **457 m/s** を置く
  (節点値は 0)。実装式をそのまま再現した厳密値。
- その 9 ノードの密度が 27 step で **0.00358 → 0.00003 (120 倍減)**、圧力 1159 → 32 Pa。
  周囲のベース帯 41 ノードの平均圧力は 7 % しか落ちない = **壁 CV だけが空になる**。
- **1 次は 40 step 完全に不動**。`nodeWallDirichlet: 0` でも 40 step 完走 (ただし壁ノードが 1500 m/s で動く)。
- **CFL 律速でない** (1/10 の CFL で 6.6 倍の step)。**`implicitRelax` は 1.4 倍の延命のみ**。
- **格子を細かくすると悪化する**: 壁 CV 体積 3.3 倍減 → 落ちる step 2.8 倍早い (ほぼ体積比例)。

**測定器**: `solver_density_cuda/tools/check_wall_cv_drain.py` (本 plan で追加済み)。
**ふるいであって予測器ではない**: 風上化を含まないので偽陽性がある
(本ツールが最大と出したカウル後縁ノードは実際には密度が上がる。判定は `ro` の時系列で行う)。

### 4.3 他ケースとの比較

`case/18.backstep` (後ろ向きステップ = 同じく流れに正対する壁) の node・2 次・`nodeWallDirichlet: 1` の
**収束解**では、速度ピン壁ノード 474 個の漏れは最大 4.1 /s で、1 ms 以内に空になるノードは 0。
SERN のベースは 1.2e5 /s。**同じ機構が 5 桁弱い**。
ただし**収束解と過渡の比較なので公平ではない** (収束解では正味流束は定義上ほぼ 0)。
言えるのは「backstep はその状態に到達できている、SERN のベースは到達前に空になる」まで。

### 4.4 SU2 の実装 (`.external/su2-src`, v8.5.0 を実際に読んだ結果)

**結論: 壁ノードの扱いは forge とほぼ同じで、ここに差は無い。** 期待していた「SU2 は壁ノードを特別扱いしている」は**外れ**だった。

| 論点 | SU2 | forge |
| --- | --- | --- |
| no-slip の課し方 | **強制 Dirichlet のみ** (`CNSSolver.cpp:508-525` `SetVelocity_Old` + `LinSysRes(iPoint,iDim+1)=0`、`:577-581` `Jacobian.DeleteValsRowi`) | 同型 (状態ピン + 残差射影 + 陰解法の行単位化) |
| 壁面からの対流質量流束 | 静止格子では **0** (`Res_Conv = 0`、書き込みはエネルギー行のみ `:543`) | 0 (bvar の `Vn_m = 0`) |
| **壁ノードの連続の式** | **固定しない** (`LinSysRes(iPoint,0)` に触る処理は存在しない) | **固定しない** |
| **壁ノードでの MUSCL** | **特別扱いも除外も存在しない** (`Upwind_Residual` に境界判定なし。`GetPhysicalBoundary` は JST のセンサーにしか使われない) | 同じく特別扱いなし |
| 再構成の目標点 | **辺の中点** (`CEulerSolver.cpp:1925` `V_i + 0.5*lim_i*Project_Grad_i`) | **辺の中点** (node は `g_reconEdgeMid` 常時 1) |
| **リミッタの評価点** | **辺の中点** (`computeLimiters_impl.hpp:175-180` `dist_ij = 0.5*(coord_j-coord_i)`) = **再構成点と厳密一致** | **双対面重心** = **不一致** ← 唯一の実質的な差 |

**forge に無い SU2 の汎用セーフガード (壁専用ではない)**:

1. **再構成の物理性検査 → そのエッジを 20 反復 1 次に固定** (`CEulerSolver.cpp:1940-1962` + `CFlowVariable.hpp:80-87` `UpdateNonPhysicalEdgeCounter`)。負の p/rho/音速のみ。
2. **密度・エネルギー更新率の自動 under-relaxation** (`CFVMFlowSolverBase.inl:599-630`、既定 `MAX_UPDATE_FLOW = 0.2`) = 1 非線形反復で 20 % 以上変われない。
3. **非物理状態になった点の解を前反復値へ巻き戻す** (`CNSVariable.cpp:151-178`)。

### 4.5 診断の座標を間違えていた (2026-09-19, codex plan レビュー Critical 1)

**§4.5/§4.6 の当初の数値と結論は、測定器がソルバと違う点で再構成していたため無効**だった。

`CELLS/centCoords` は**実行時にノード座標へ置換される** (`main.cpp` の `nodeValueAtNode` → `mesh.cpp`。
`run_0208/forge_run.log` に `centCoords <- node coords for 62317 CVs (max centroid shift 0.03593721241)` と出ている)。
メッシュ HDF5 の `centCoords` をそのまま使うと別の点で再構成したことになる。訂正後:

| 診断量 | 誤 (centCoords) | 正 (ノード座標) |
| --- | --- | --- |
| 62100→62107 の `Ux_L` | 457.42 m/s | **651.77 m/s** |
| 62104→62111 の `Ux_L` | 365.77 m/s | **521.17 m/s** |
| `Ux_L / (0.5 x Ux_隣接)` | 0.70 | **0.999996** |
| 中点とのずれが辺長 1 % 超の内部面 | 58899 (47.7 %) | **48328 (39.1 %)** |

**なお成立するもの**: 再構成は片側勾配の理論値 `0.5 x u_隣接` に**厳密に一致**しており、
SU2 が出す値と同じ。**「壁ノードからの再構成が異常」という当初仮説の棄却は維持する**
(ただし「SU2 より保守的 (0.70 倍)」は誤りで、**一致**が正しい)。

**崩れたもの**: 「非物理な再構成はベース無しでも同じように出るから原因ではない」は**誤り**。
`tools/check_face_reconstruction.py` (本 plan で作り直し) でノード座標で数え直すと:

| step | ベース有り `run_0208` | ベース無し `run_0209` |
| --- | --- | --- |
| 1 | 0 面 (最小 1056.56 Pa) | 0 面 (最小 1375.90 Pa) |
| 5 | 0 面 (429.45 Pa) | 0 面 |
| **6** | **1 面 (−66.78 Pa)** | 0 面 |
| 27 / 40 | **3 面 (−296.25 Pa)** | **0 面 (40 step 通して)** |

**対照は完全にクリーンで、ベース run だけが非物理な再構成を作る。**
最初に出るのは x = −0.05 m (入口面) の 47703→47702 で、ベース自身の面も step 27 には負になる。
**SU2 の `bad_recon` 検査はベース run で発火し、対照では発火しない** → §4.4 の (1) は
「本件に無関係」ではなく**直接の候補**。

### 4.6 「1 次は壁 CV が質量を失えない」も誤り (codex Major 2)

SLAU の質量流束は `convectiveFlux_slau_d.inc.cuh:493`:

```
mdot = sss*0.5f*((ro_L*(Vn_p+Vn_hat_p_abs)+ro_R*(Vn_m-Vn_hat_m_abs)) - chi/(c_diss)*P_del)
```

**左右の速度がともに 0 でも、圧力差 `P_del` だけで質量が動く**。
したがって「1 次では壁ノード側が風上のとき `rho u = 0` なので壁 CV は質量を失えない (一方弁)」という
説明は成り立たない。**1 次が 40 step 不動だったのは観測事実として残すが、その理由づけは撤回する。**
原因は SLAU 面流束の符号付き総和・密度残差・陰的補正量で論じ直すこと。

同じ理由で、旧 `check_wall_cv_drain.py` の質量流束見積り (`rho_L u_L . S` を積むだけ) は
**座標と欠落項の二重に不正確**だったので削除し、`check_face_reconstruction.py`
(再構成状態の物理性のみを、実装と同じ式・同じ座標で検査) に置き換えた。

### 4.7 SU2 との比較で言えることの範囲 (codex Major 3)

§4.4 の比較表で引いた `CNSSolver.cpp:508` は**熱流束壁**の処理で、対象 run は**等温壁**。
SU2 の `BC_Isothermal_Wall_Generic` は温度差から熱流束を作りエネルギー残差と Jacobian に加える
(行削除は運動量のみ)。forge は `nodeWallDirichlet_d.cu` で `T/roe/P` をピンし等温壁のエネルギー残差を 0 にする。
**熱的閉包は同型ではない。**

また SU2 の Venkatakrishnan の平滑化係数は参照長さと設定係数依存 (`CLimiterDetails.hpp`)、
forge は局所 `volume` 依存で、勾配の作り方も未比較。
**「同条件なら SU2 のリミッタは 1」「SU2 より小さい」はソースの形だけからは言えない。**

言えるのは**構造上の一致に限る**: 速度 Dirichlet / 連続の式を解く / 内部辺を再構成する / 再構成の目標点は辺中点。
**同一メッシュ・対応 BC で実際に走らせるまで SU2 の成否は予測しない。**

### 4.8 W1 の設計: 再構成増分の列挙と共通化 (2026-09-19, codex Major 5 対応)

**ずれは 2 成分ある。** 「座標を揃える」だけでは足りない。

**(1) 流束側が実際に適用する増分** (`convectiveFlux_common_d.cuh` の `interp_dispatch`)。
`phi_face = phiC + psi * Δ` の Δ は `convMethod` で違う:

| `convMethod` | 関数 | 増分 Δ |
| --- | --- | --- |
| 0 / −1 | `interp_1stUp` | **0** (`return phiC`。ψ は無関係) |
| 1 | `interp_MUSCL_2nd` | `g · cpd` |
| **2** | `interp_MUSCL_3rd` | `0.5·k·(phiD − phiC) + (1−k)·(g · cpd)`, k = 1/3 ← **隣接値差を含む** |
| その他 | `interp_MINMOD` | 別式 (`f`・`dcc` を使う) |

`cpd` は目標点までのオフセットで、**node では常に `±0.5·dcc` (エッジ中点)**
(`convectiveFlux_d.cu` の `rem = (discretization=="node") ? 1 : 0`)、cell では `pc[ip] − cc[ic]` (双対面重心)。

**(2) リミッタ側が評価する増分**:

| 経路 | 評価点 | 増分 |
| --- | --- | --- |
| 通常 (`limiter_d.cu` の `limiter_r1_fused5_d` pass2) | `pc[ip] − cc[ic0]` (**双対面重心**) | `g · dcp` のみ |
| 周期 node (`limiterPeriodic_d.cuh`) | 同じく `pc[ip] − cc[ic0]` | `g · dcp` のみ (`SCALED` 版は `× inv_ref`) |

**したがって node では、点が違う (重心 vs 中点) うえに、`convMethod: 2` では形も違う**
(リミッタは隣接値差の項を**一度も見ない**)。SU2 はリミッタ側でも `umusclProjection` を適用して形まで揃えている
(`computeLimiters_impl.hpp`)。

**設計**: 増分を作る関数を 1 つにし、流束とリミッタが**同じ関数・同じ引数**で呼ぶ。

```
__device__ flow_float recon_increment(int scheme, flow_float phiC, flow_float phiD,
                                      flow_float gx, flow_float gy, flow_float gz,
                                      flow_float cpdx, flow_float cpdy, flow_float cpdz, ...);
// 流束:     phi_face = phiC + psi * recon_increment(...)
// リミッタ: delta_m  = recon_increment(...)   ← psi を求める材料
```

**必要な入力はすべて揃っている**: リミッタ pass2 で partner は pass1 と同じ一行
`ic1 = plane_cells[2*ip+0] + plane_cells[2*ip+1] - ic0` で取れる。`phiD = Q[k][ic1]` も pass1 で既読。
エッジ中点オフセットは `0.5*(cc[ic1] − cc[ic0])`。

**float32 の桁落ちも同時に直す** (codex Major 5): 現行は `Qt = qc + g·dcp` を作ってから `Qt − qc` を渡しており、
**足して引く往復で桁を落としている**。共通関数は Δ を直接返すのでこの往復が消える。
変更前後で `Δ` の相対差を float32/float64 で測って記録する。

**周期経路も同じ関数を通す** (`limiterPeriodic_d.cuh`)。`SCALED` 版のスケーリングは Δ を作った後に掛ける。
**cell は対象外** (cell の目標点は双対面重心のままで整合しているので、分岐を維持する)。

### 4.9 W1 の合否ゲート (実装前に固定, codex Major 6/7 対応)

**Barth と Venkatakrishnan で要求を分ける。**

| ゲート | 対象 | 要求 |
| --- | --- | --- |
| G1 厳密有界性 | `limiter: 1` (Barth) | 再構成した面値が**近傍の min/max を外れる面が 0** |
| G2 平滑化込みの逸脱 | `limiter: 2` (Venkatakrishnan) | 逸脱は許すが**逸脱量が改修前より減る**こと。`eps2 = volume` が正なので厳密有界は保証されない (例: `δmin=0, δm=−1, volume=1` で ψ=1/3) |
| G3 正値性 | 両方 | `check_face_reconstruction.py` の VERDICT。**W1 単独で 0 になることは要求しない** (それは W2 の役割) |
| G4 解の非退行 | 両方 | `check_field_regress.py --boundary` をノイズ床基準で。**リミッタが強く効くので解は動きうる** — 動いたら SU2 クロスチェック (`procedures/su2-cross-check.md`) で是非を決める |

**共通 IC と設定差分を固定する**: `run_0205/res_0.h5` を共通初期場に、`output.level: 2`・`outStepInterval: 1`。
差分は `space.limiter` と新キーのみ。メッシュ品質 VERDICT・全残差・`check_quasisteady.py` を毎回貼る。

**標準ケース回帰** (`procedures/verification/README.md` から具体名で指定):
`case/05.sod_shock_tube` (衝撃波)、`case/36.passive_pseudoshock_control` (擬似衝撃波・node 壁)、
`case/09.Taylor-Green` (周期 = W1 の周期経路)、`case/23.axi_nozzle` (軸対称)、
`case/44.vitiated_air_wt` (TP・凝縮)。変更範囲に応じて増やす。

### 4.10 W1 の実装と A/B の実測 (2026-09-19)

**実装**: `space.limiterMatchRecon` (既定 0)。`limiter_d.cu` の `limiter_r1_fused5_d` と
`limiterPeriodic_d.cuh` の `limiter_psi_merged_d` に分岐を足し、`=1` のとき

- 目標点を **node ならエッジ中点** (`0.5*(cc[ic1]-cc[ic0])`) にする
- `convMethod: 2` なら **隣接値差の項も含めた増分**にする (`interp_MUSCL_3rd` と同形)
- 増分 `delta` を**直接**作る (`Qt` を経由しないので float32 の足して引く桁落ちが無い)

周期経路には `plane_cells` を引数追加して partner を取れるようにした。cell は分岐を維持。

**既定パスの非退行**: 変更を外した HEAD の worktree でバイナリを作って比較した。
**ビット同一ではないが、差は run 間非決定性の水準**:

| | 同一バイナリの反復 (ノイズ床) | 変更前バイナリ vs `mr:0` |
| --- | --- | --- |
| step 1 `ro` | 1.96e-07 | 4.51e-07 |
| step 10 `ro` | 8.08e-07 | 6.52e-07 |
| step 27 `ro` | 7.90e-07 | 1.34e-06 |

(この診断 run は step 28 で発散するのでノイズがカオス的に増幅する。`limiter_P` はノイズ床自体が 0.96。)
**カーネル引数を足すとレジスタ割り当てと命令スケジューリングが変わる**ので、式が同一でもビット同一にはならない。
コード中の「ビット不変」というコメントは実測に合わせて「式は変更前と同一」に訂正した。

**A/B (共通初期場 `run_0205/res_0.h5`, cfl 1.0, 2 次)**:

| 構成 | 結果 |
| --- | --- |
| Venkatakrishnan + `mr:0` (既定) | step 28 で NaN |
| Venkatakrishnan + `mr:1` | step 31 (**3 step の延命のみ**) |
| Barth + `mr:0` | step 27 で NaN (ベース壁ノードの `ro` が **−0.001003** = 負に) |
| **Barth + `mr:1`** | **400 step 完走・NaN 無し** |

**ゲート判定**:

- **G1/G3**: Barth + `mr:1` は step 1〜400 すべてで**非物理な再構成ゼロ**
  (`check_face_reconstruction.py` VERDICT OK、再構成 P の最小は 338〜1065 Pa)。
- ベース壁ノードの密度: `mr:0` は 26 step で負になるが、`mr:1` は 0.003581 → 0.001039 で**頭打ち**
  (step 100 → 400 で 0.001156 → 0.001039 = 300 step で −10 %)。**排出が止まる**。
- **G2**: Venkatakrishnan は 3 step しか延びない。`eps2 = volume` の平滑化が残るので厳密有界にならない
  = codex Major 6 の分離がそのまま出た。

**ただし根治ではない**。全レシピ (`problem_r4e_2d_base_barth_mr1.yaml`, `run_0220_r4e_barth_mr1`) を回すと
**mid 段 step 129 で NaN** (既定は 28。4.6 倍)。NaN は同じベース近傍 (x/H 10.66–11.19, y/H 2.70–2.77 の 55 節点)。

**読み**: W1 はリミッタの正しさの問題として**独立に価値がある** (Barth の厳密有界性が実際に効くようになる)。
ベース発散も大きく緩和するが**消さない**ので、W2 (面単位の非物理フォールバック) と
SERN 側の解像度・形状の検討はどちらも残る。

## 5. 実装ステップ

1. 設計確定 (§4.5 の未決を codex レビューで決める)
2. 実装 (1 案。opt-in の config キーで既定 off から始め、A/B が取れる形にする)
3. §6 の検証
4. 既定化の可否を判断

### 5.1 残作業 (優先順)

| # | 項目 | 内容 |
| --- | --- | --- |
| W0 | ~~測定器の座標~~ **完了 (2026-09-19)**: `check_wall_cv_drain.py` を削除し `check_face_reconstruction.py` に置換 (ノード座標・再構成状態の物理性のみ)。§4.5/§4.6 の結論を訂正 |
| W1 | ~~リミッタ評価点を再構成点に揃える~~ **実装 + A/B 完了 (2026-09-19) → §4.8 設計 / §4.9 ゲート / §4.10 実測**。`space.limiterMatchRecon` (既定 0)。**Barth と組むと G1/G3 を満たし** (400 step 非物理ゼロ、ベース壁ノードの `ro` が頭打ち)、Venkatakrishnan では 3 step しか延びない。**既定パスの差は run 間非決定性の水準** (ビット同一ではない)。**ベース発散の根治ではない** (全レシピは mid 段 step 28 → 129)。残 = 標準ケース回帰 (§4.9 の 5 ケース) と既定化の判断 |
| W2 | **面単位の非物理フォールバック** (SU2 の `bad_recon` 相当) を**導入候補として独立検証**。§4.5 のとおりベース run でのみ発火し対照では発火しないので**本件の直接候補**。組成・エンタルピーも整合して再計算すること |
| W3 | **更新率制限と点ロールバックは別 plan へ**。forge には既に `update_d.cu` の `updateGuardScale` があり、[`accepted/time_integration-update-positivity-geard.md`](../accepted/time_integration-update-positivity-guard.md) に「試験した CFL 上限を改善しなかった」実測がある。SU2 の `MAX_UPDATE_FLOW` をそのまま移植する設計では、流れ 5 変数の後に別途更新される SST・化学種・凝縮との整合 (ΣρY=ρ、組成依存 EOS、周期共有 DOF、TP のエネルギー基準) を失う。**3 件一括の既定 on は推奨されない** |
| W4 | ~~壁ノードからの速度再構成を制限する~~ **撤回 (§4.5)**。再構成は片側勾配の理論値に厳密一致しており異常でない |
| W5 | **タイトル・索引・§7 の同期** (codex m1)。`plans/README.md` の説明が旧方針のまま。`update_d.cu` が §7 に無い |

**ベース発散の原因は未確定のまま残す** (codex Major 7)。「未解像が主因」と断定するのは早い —
§4.5 で SU2 の検査が発火することが分かった以上、W2 の A/B を局所格子試験より**先に**行う。
[`tooling-nozzle-sern-3d.md`](tooling-nozzle-sern-3d.md) R4e にもそう書き戻す。

## 6. 検証

**W1 (リミッタ評価点)** — 主眼は「有界性が実際に効くようになること」と「解が壊れないこと」:

1. `tools/check_wall_cv_drain.py` 系の面再構成診断で、**辺中点での再構成が近傍 min/max を外れる面が 0** になること
   (現行はカウル後縁で毎 step 負圧。`case/46` の 2D で測る)。
2. 既存 run の非退行: `check_field_regress.py --boundary` をノイズ床基準で。
   **リミッタが強く効くようになるので解は動きうる** — 動いた場合は「どちらが正しいか」を
   SU2 クロスチェック (`procedures/su2-cross-check.md`) で決める。
3. 標準ケース回帰 (`procedures/verification/README.md`)。特に衝撃波を持つケース (`case/05` Sod、`case/36`)。

**W2 (セーフギヤード)** — 主眼は「起動が頑健になること」と「収束解を動かさないこと」:

4. 発動回数を出力する (SU2 の `counter_local` 相当)。**収束後は 0 に落ちること**を確認する
   — 収束解で発動し続けるなら解を歪めている。
5. 既存の起動困難ケースで段数・step 数が減るか (`case/46` の梯子、`case/16`)。

許容差は実装案を決めてから数値で宣言する。

### 6.1 レビュー記録 (codex)

| stage | 日付 | 記録 | 判定 | 採否 |
| --- | --- | --- | --- | --- |
| plan | 2026-09-19 | [2026-09-19-convection-node-wall-reconstruction-plan.md](../../notes/reviews/2026-09-19-convection-node-wall-reconstruction-plan.md) | NO-GO, C1/M6/m1 | **全件採用**。C1 → 測定器の座標を修正し §4.5/§4.6 の結論を訂正 (対照はクリーン・ベースのみ非物理再構成 = SU2 の検査は発火する)、M2 → 「1 次は一方弁」を撤回 (SLAU の質量流束に圧力差項)、M3 → SU2 比較を構造一致に限定 (等温壁・勾配・リミッタ係数は未比較)、M4 → W3 を別 plan へ (既存 `updateGuardScale` と連成状態)、M5 → W1 に周期経路と 3 次を含める、M6 → 有界性ゲートを Barth と Venkat で分ける、M7 → ベース原因は未確定のまま残し W2 を格子試験より先に、m1 → 索引と §7 を同期 |

## 7. 影響範囲

`cuda_forge/convection/convectiveFlux_*.inc.cuh` / `convectiveFlux_common_d.cuh`、
`cuda_forge/limiter_d.cu` + `limiterPeriodic_d.cuh` (W1 は周期経路も対象)、
`cuda_forge/update_d.cu` (W3 を別 plan に出すまでの参照先: 既存 `updateGuardScale`)、
`input/solverConfig.cpp` (新キー)、`solver_density_cuda/tools/check_face_reconstruction.py`。
**node の 2 次以上の全計算に影響する** — SERN に閉じない。

## 8. 完了条件

- W1: リミッタの評価点が再構成点と一致し (通常・周期・次数すべて)、解の変化が宣言した許容内
- W2: 面単位フォールバックの発火回数が収束後 0 に落ち、ベース最小再現への効果が測れている
- 標準ケース回帰が宣言した許容内 (対象ケースは実装前に具体名で固定する)
- 既定化するかどうかの判断が plan に書かれている

## 9. 変更ログ

- `2026-09-19` — 起票。case/46 R4e の原因調査から分離。測定器 `check_wall_cv_drain.py` を追加。
- `2026-09-19` — **codex plan レビュー NO-GO (C1/M6/m1) を全件採用**。**測定器が `CELLS/centCoords` を使っていたが
  実行時はノード座標に置換される**ため、§4.5/§4.6 の数値と結論が無効だった。ツールを
  `check_face_reconstruction.py` に置換して訂正。**対照 (ベース無し) は 40 step 非物理な再構成ゼロ、
  ベース run は step 6 から発生** = SU2 の `bad_recon` はベースで発火する → W2 を格子試験より先に。
  「1 次は一方弁」も SLAU の圧力差項で成り立たないので撤回。
- `2026-09-19` — SU2 v8.5.0 のソースを読み、当初の仮説を棄却 (§4.4)。壁ノードの扱い
  (強制 Dirichlet・連続の式は自由・MUSCL の特別扱い無し) は SU2 も forge と同型で差が無い。
  forge の再構成値は隣接の半分の 0.70 倍 = **SU2 が出す値より小さく異常でない**。
  SU2 の 3 つのセーフガードもこの排出を止めない (20 % 上限に掛かるのは 27 step 中 5 step だけ)。
  **確認された唯一の実質差はリミッタの評価点**。plan の主題をそちらへ移し、ベース発散は
  `tooling-nozzle-sern-3d.md` R4e (未解像) に戻した。
