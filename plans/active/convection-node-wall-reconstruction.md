# node 速度ピン壁ノードからの対流再構成 (壁 CV の質量漏れ)

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

### 4.5 **当初の仮説 (壁ノードの再構成が異常) は棄却** (2026-09-19)

SU2 を読んだうえで forge の再構成値を測り直した結果、**再構成は異常ではなかった**。

| 壁ノード | 節点値 `Ux` | 隣接ノードの `Ux` | 再構成 `uL` | `uL / (0.5 x Ux_隣接)` |
| --- | --- | --- | --- | --- |
| 62100 | 0 | 1303.54 | 457.42 | **0.70** |
| 62104 | 0 | 1042.35 | 365.77 | **0.70** |

壁ノードから辺中点への再構成は、片側勾配なら `0.5 x u_隣接` になる。forge はその **0.70 倍**で、
**SU2 が出す値より小さい**。`psi_Ux` も 1.0 で、SU2 の Venkatakrishnan も同条件で 1 を返す
(壁ノードは近傍の最小値なので `fieldMin - field_i = 0`、`limMax` 側は `delta = u_max - 0` で有界)。
**したがって「壁法線成分を落とす」「壁ノードの再構成を制限する」という当初案は前提が成立しない。撤回する。**

**SU2 のセーフガードも救わない**。ベース壁ノード 62104 の 1 step あたり密度変化率は
4.2 % → 6.4 % → 11.9 % → 21.8 % → 60.4 % と育ち、**27 step のうち 20 % を超えるのは最後の 5 step だけ**。
`MAX_UPDATE_FLOW = 0.2` は終盤を遅らせるだけで、最初の 22 step の単調な排出は素通りする。
再構成も非物理にならない (p も rho も正) ので (1) も (3) も発火しない。
**SU2 で回しても同じように抜けた可能性が高い。**

### 4.6 残る仮説: ベース近傍が未解像で、離散的な定常解が無い

いま最も確からしい読み:

- ベース (厚さ 0.02 H) を横切るノードは 9 個、その直後の流れ方向間隔は**ベース厚の 4.1 倍**。
  背面の近傍後流 (ベース圧を決める再循環) は**まったく解像されていない**。
- 初期値のベース圧 1159 Pa は外気 2851 Pa の 0.41 倍で、超音速外流の base pressure として**物理的に妥当な値**。
  そこから 32 Pa (真空) へ落ちるのは数値的破綻であって物理ではない。
- **1 次が安定なのは解けているからではない**: 1 次では壁ノード側が風上のとき `rho u = 0` なので、
  壁 CV は**対流で質量を失えない** (増えることはある = 一方弁)。これも人工物であり「正しい定常解」ではない。
- 私が回した格子試験 (§4.15.7 の `split_plume_at_te`) は **Δx を縮めたが壁 CV も 3.3 倍小さくした**ので、
  後流を解像する試験になっていない。**ベースを横切る点数と後流の解像を分けて振る必要がある**。

**したがって本 plan の主題は「壁ノードの再構成の直し方」から次へ移す**:

1. **リミッタ評価点を再構成点に揃える** (SU2 と同じにする)。これは**確認された唯一の実質差**で、
   ベース発散の原因ではないが直すべき欠陥 ([`tooling-nozzle-sern-3d.md`](tooling-nozzle-sern-3d.md) §4.15.5:
   内部面の 47.7 % でずれ、カウル後縁で毎 step 負圧が流束に入る)。
2. **SU2 流の汎用セーフガードを入れる** (非物理再構成のエッジ単位 1 次化 / 更新率の上限 / 点ロールバック)。
   本件は救わないが、起動時の頑健性として独立に価値がある。
3. **ベース近傍の解像は SERN 側の課題に戻す** ([`tooling-nozzle-sern-3d.md`](tooling-nozzle-sern-3d.md) の R4e)。
   codex が指定した順序 (④ 局所格子 → ⑤ 肩の丸めと厚み) は正しかった。**⑤ を「筋が悪い」とした私の判断は撤回する。**

## 5. 実装ステップ

1. 設計確定 (§4.5 の未決を codex レビューで決める)
2. 実装 (1 案。opt-in の config キーで既定 off から始め、A/B が取れる形にする)
3. §6 の検証
4. 既定化の可否を判断

### 5.1 残作業 (優先順)

| # | 項目 | 内容 |
| --- | --- | --- |
| W1 | **リミッタ評価点を再構成点に揃える** (最優先)。`limiter_d.cu` の pass2 は試行値を双対面重心で評価するが、node の流束は辺中点へ再構成する。SU2 は両方とも辺中点で厳密一致 (`computeLimiters_impl.hpp:175-180` vs `CEulerSolver.cpp:1925`)。内部面の 47.7 % でずれ、カウル後縁では毎 step 負圧が流束に入る。**opt-in キーで入れて A/B を取る** |
| W2 | **SU2 流の汎用セーフガード**を入れるか判断する。(a) 再構成の物理性検査 → そのエッジを N 反復 1 次化 (`UpdateNonPhysicalEdgeCounter` 相当)、(b) 密度・エネルギー更新率の上限 (`MAX_UPDATE_FLOW` 相当)、(c) 非物理状態の点ロールバック。**本件のベース発散は救わない**が起動時の頑健性として独立に価値がある。既定 on にするかは要判断 (解を動かしうる) |
| W3 | **検証** (§6)。W1 は「解が変わらないこと」が主、W2 は「起動が頑健になること」が主で、判定条件が違う |
| W4 | ~~壁ノードからの速度再構成を制限する~~ **撤回 (§4.5)**。再構成は SU2 と同等かそれより小さく、異常ではなかった |

**本 plan から外したもの**: ベース近傍の発散は [`tooling-nozzle-sern-3d.md`](tooling-nozzle-sern-3d.md) R4e に戻す
(§4.6 のとおり未解像が主因の可能性が高く、ソルバ側の課題ではない)。

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
| plan | (未実施) | — | — | W1 のため実施予定 |

## 7. 影響範囲

`solver_density_cuda/cuda_forge/convection/convectiveFlux_*.inc.cuh` / `convectiveFlux_common_d.cuh`、
`cuda_forge/limiter_d.cu` (W5 を取り込む場合)、`input/solverConfig.cpp` (新キー)。
**node の 2 次以上の全計算に影響する** — SERN に閉じない。

## 8. 完了条件

- 最小再現が完走し、壁 CV の `ro` が保たれる
- 標準ケースの Cf・壁熱流束・再付着位置が宣言した許容内
- 既定化するかどうかの判断が plan に書かれている

## 9. 変更ログ

- `2026-09-19` — 起票。case/46 R4e の原因調査から分離。測定器 `check_wall_cv_drain.py` を追加。
- `2026-09-19` — **SU2 v8.5.0 のソースを読み、当初の仮説を棄却** (§4.4/§4.5)。壁ノードの扱い
  (強制 Dirichlet・連続の式は自由・MUSCL の特別扱い無し) は SU2 も forge と同型で差が無い。
  forge の再構成値は隣接の半分の 0.70 倍 = **SU2 が出す値より小さく異常でない**。
  SU2 の 3 つのセーフガードもこの排出を止めない (20 % 上限に掛かるのは 27 step 中 5 step だけ)。
  **確認された唯一の実質差はリミッタの評価点**。plan の主題をそちらへ移し、ベース発散は
  `tooling-nozzle-sern-3d.md` R4e (未解像) に戻した。
