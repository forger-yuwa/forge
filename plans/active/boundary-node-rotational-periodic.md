# node 周期 (seam 合算経路) の回転周期対応

## メタ

- **area**: `boundary`
- **status**: `draft`
- **related_docs**:
  - [`methods/boundary.md`](../../methods/boundary.md) (「周期境界」節、node の回転周期の表)
  - [`methods/discretization.md`](../../methods/discretization.md) §2.5.5 (periodic 双対面対応)
- **related_plans**:
  - [`discretization-median-dual-3d.md`](discretization-median-dual-3d.md) §4.5 (node 周期の骨格)、§4.5.8 (回転周期の既存設計 — 本 plan で実装)
  - [`convection-slau-wall-normal-chi-default.md`](convection-slau-wall-normal-chi-default.md) (周期域の既定化。本 plan の V4 が回転周期での検証)
  - [`convection-slau-wall-normal-chi.md`](../accepted/convection-slau-wall-normal-chi.md) (V6 P6-per: seam 両メンバの `wall_flag` 一致)
- **created**: `2026-09-26`
- **owner**: `CFD Dev`

## 1. 目的

**ユーザ指示 (2026-09-26)**: node の周期は seam 合算の別経路 (`periodicNode_d.cu`) で、回転 (`type: 1`, `dtheta`) の処理が無い
(回転は cell のゴースト経路 `boundaryCond_d.cu:1728` のみ)。**seam 経路に回転を実装し、各種スカラーも含めて対応したうえで、SLAU
(`slauWallNormalChi` 込み) を回転周期で検証する**。完了時、node で軸まわりのセクタ (例: 90°) を回転周期で解けるようになる。

## 2. スコープ

- **やる**: 角の割当 (重み付き union-find、閉ループ検査)、node でのパートナー探索 (回転座標一致)、ベクトル量の gather / broadcast / ミラー / dq /
  勾配 / 陰解法の対角ブロックへの回転の挿入、リミタの seam 方針、検証ケース (環状セクタと回転コピーの全周)。
- **やらない**: cell 経路の変更。軸 ($r=0$) を含むセクタ。軸対称との併用 (起動エラーにする)。複数の回転軸 (x 軸のみ)。
  リミタの速度成分の frame 整合 bound (円筒成分での bound は将来)。

## 3. 関連 docs と前提

- node 周期の骨格: 継ぎ目で割れた CV を union-find (root = 最小 index) で group にし、合併体積のもとで残差を slave→root に atomicAdd →
  broadcast。保存量・k/ω・遷移・化学種・凝縮・受動種・dq をミラー、勾配を gather (`periodicNode_d.cu`、`mesh.cpp:buildPeriodicNodeGroups`)。
- 既存設計 `discretization-median-dual-3d.md` §4.5.8: 運動量の授受に $T(d\theta)$、円筒一様流で free-stream、軸を含まない環状で先に、corner 連鎖で回転合成の閉ループ恒等。
- cell の回転周期: `setPeriodicPartner` (`mesh.cpp:531-579`) が座標を回してペアリング。

## 4. 設計方針 (2026-09-26 `diagnostician`)

### 4.1 角の割当 (`mesh.cpp` `buildPeriodicNodeGroups`)

- 重み付き union-find。partner 対 $(i,j)$ (type 1、bcond の `dtheta`、$i$ が $\theta=0$ 側) に対し $\theta_j=\theta_i+d\theta$。merge 時に offset を伝播し、
  root は $\theta=0$。
- **閉ループ検査**: 全 partner 対で $(\theta_j-\theta_i-d\theta) \bmod 2\pi < 10^{-9}$。違反は起動エラー (タイル張り不整合)。
- 角は double で持ち、node 配列 `perCos_d` / `perSin_d` (flow_float、double の角から 1 回だけ作る) と、group に非零角があるかの `perRotFlag`。
- パートナー探索は `setPeriodicPartner` type 1 (座標回転) を node でも使い、回転後の座標一致 ≤ $10^{-9}H$ を変換時に検査。
- **軸対称との併用は起動エラー**。

### 4.2 変換の一覧 ($R_m$ = root 座標系 → member $m$ 座標系、x 軸まわりの 2×2、$(\cdot)_y,(\cdot)_z$ に作用)

| 量 | 経路 | 変換 |
| --- | --- | --- |
| 運動量残差 | `periodicGatherToRoot` / `Broadcast` | gather: $\mathrm{res}_{root} \mathrel{+}= R_m^{\mathsf T}\,\mathrm{res}_m$、broadcast: $\mathrm{res}_m = R_m\,\mathrm{res}_{root}$ |
| 保存量 state | `periodicMirrorNSState` | $(\rho u_y,\rho u_z)_m = R_m(\cdot)_{root}$。$\rho,\rho u_x,\rho e$・化学種・凝縮・$\rho\xi$ は恒等 |
| dq | `periodicMirrorDq` | 運動量 2 成分に $R_m$ |
| スカラー勾配 ($\nabla\rho,\nabla P,\nabla T,\nabla k,\nabla\omega,\nabla Y,\nabla\xi$) | `periodicGradientGather` | gather $R_m^{\mathsf T}\mathbf g_m$、broadcast $R_m\mathbf g$ |
| 速度勾配テンソル $G=\partial u_i/\partial x_j$ | 同 | gather $R_m^{\mathsf T} G_m R_m$、broadcast $R_m G R_m^{\mathsf T}$ |
| $\nabla\cdot\mathbf u$・スカラー量 | 同 | 恒等 |
| block-DPLUR 対角 | 陰解法の fold | $D_{root}=D_L+\mathcal R_m^{\mathsf T} D_m \mathcal R_m$ ($\mathcal R_m$ = 5×5、運動量ブロックに $R_m$)。broadcast は $\mathcal R_m D \mathcal R_m^{\mathsf T}$ |
| リミタ min/max | `periodicGatherMax/Min` | $\rho,P,T$・全スカラー・$u_x$: 従来どおり gather。**$u_y,u_z$: 回転 group では seam 越し gather をしない** (自側の隣接だけで bound。reduced な max/min は回転できない) |
| `wall_flag` | mask | 回転でも group 内で全員一致を要求 (前 plan V6 P6-per の検査、違反は起動エラー) |

### 4.3 ビット不変

`perRotFlag == 0` の mesh (並進のみ) は wrapper が現行カーネルを呼ぶ (新カーネルは回転 group のみ)。

## 5. 実装ステップ

1. codex plan 段。
2. 角の割当・閉ループ・`wall_flag` 一致・座標一致検査 (`mesh.cpp`・`mesh.hpp`)。
3. 回転の gather / broadcast / mirror / dq / 勾配 (`periodicNode_d.{cu,cuh}`)。
4. リミタ方針 (`limiter_d.cu`・`limiterPeriodic_d.cuh`)。
5. 陰解法の対角 fold の相似変換。
6. 検証ケース `case/57.annular_sector` (gmsh の 90° 環状セクタと、回転コピーの全周)。
7. §6 の run。docs。codex result。

### 5.1 残作業 (優先順)

**計算資源**: 単体・小規模はローカル、全周など重いものは AWS (ユーザ指示: ローカルで大きな計算をかけない。AWS は他セッションと共有、手動 stop しない)。

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | codex plan 段 | Critical/Major は `diagnostician` に諮る | F |
| 2 | 角の割当・検査 | `solver_density_cuda/mesh/mesh.cpp` (`setPeriodicPartner` の node 対応、`buildPeriodicNodeGroups`)、`mesh.hpp`。合格: §6 U0・U1 | O |
| 3 | 回転の授受 (**cuda_forge の数値変更、条件 6: 実装前に `diagnostician` レビュー**) | `solver_density_cuda/cuda_forge/periodicNode_d.{cu,cuh}`。合格: §6 V1・V2 | O (F レビュー) |
| 4 | リミタ方針 | `limiter_d.cu`、`limiterPeriodic_d.cuh`。合格: §6 V6 | O |
| 5 | 陰解法の対角 fold (F) | 陰解法の fold 箇所。合格: §6 V7 | O (F レビュー) |
| 6 | メッシュ | `case/57.annular_sector/mesh/*.geo`、README。合格: §6 U1 | O |
| 7 | 検証 run | case/57、case/39、case/09。README の run 表 | O (結論 F) |
| 8 | docs | `methods/boundary.md` の「実装中」を外す、`discretization-median-dual-3d.md` §4.5.8 に「実装済 → 本 plan」 | O |
| 9 | codex result → accepted | | F |

## 6. 検証

**判定基準はすべて測る前に固定する。許容は既存のものを流用する。**

| # | 試験 | 合格 | 資源 |
| --- | --- | --- | --- |
| U0 | 単体: $R$ の直交性、$R(\theta)R(-\theta)=I$、3 群 corner 連鎖の閉ループ $=I$ (合成誤差 ≤ 1e-12)、閉ループ違反の合成入力で起動エラー | 全 PASS | ローカル |
| U1 | メッシュ: gmsh 環状セクタ (x 軸、内径/外径、90°、ヘキサ ≈ 5 万節点) と、同セクタを 4 回回転コピーした 360° 全周 (節点が一致)。node 変換、`check_mesh_quality` PASS、全周期ノードに partner あり、回転後座標一致 ≤ 1e-9 H、`wall_flag` 群一致 | 全 PASS | ローカル |
| V1 | free-stream: 純軸流 $u_x=U$ ($u_y=u_z=0$) を Euler・slip で 100 step | seam 残差が並進周期の free-stream 試験と同水準 (既存の閾値を流用) | ローカル |
| V2 | 同値性 (Euler、slip、旋回入口 $u_\theta=\Omega r$): 全周 vs セクタ、同 step・同 cfl | セクタ領域の節点で場の差 ≤ 同一バイナリ反復 3+3 のノイズ床 ×2 (`check_field_regress`) | ローカル (重ければ AWS) |
| V3 | スカラー保存: 2 成分 (径方向プロファイル入口) + k-ω + 受動種 | seam 越しの gather 後の残差総和 = 境界流束 (閉包誤差 ≤ 1e-12 相対)、全周との差 ≤ ノイズ床 ×2 | ローカル |
| V4 | **SLAU + `slauWallNormalChi` (auto 1)**: 滑りなし内外筒、SST 低 Re、`nodeWallDirichlet: 1` | 成立条件 (`wall_flag` 面 > 0、初回 `FORGE_DUMP_MASSFLUX` で flag 差 ≥ 1 面、seam 群一致)。全周 vs セクタで壁 $p$ の $L^\infty$ ≤ 0.5 % (case/16 V3)、内外筒 $C_f$ 相対 L2 ≤ 1 % (case/48 V3)、床到達 0・NaN 0 | ローカル / AWS |
| V5 | 並進回帰: case/39 `run_0007` 終端から旧/新バイナリ 200 step、case/09 TGV 同様 | 初回面流束ビット同一 + 場はノイズ床 ×2 以内 | ローカル |
| V6 | リミタ seam: V4 の 2 次 + limiter 2 | seam 節点の $u_y,u_z$ が全周解の [min, max] ± 0.5 % 内 (新極値なし) | (V4 と同時) |
| V7 | 陰解法: V4 のセクタ / 全周で `check_convergence` | セクタ側 rms 末尾平均 ≤ 2× 全周、RISING なし | (同時) |

- **不合格時**: V1/V2 → 変換の符号・経路のバグとして設計へ戻る。V6 のみ不合格 → リミタ方針を「円筒成分 bound」へ差し替える §5.1 項目を起こす。
  V4 不合格 → 回転周期を `slauWallNormalChi` の auto → 0 の条件に加える (既定化 plan の除外表と整合)。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |

### 6.2 結果

(未実施)

## 7. 影響範囲

- `solver_density_cuda/mesh/mesh.{cpp,hpp}`、`solver_density_cuda/cuda_forge/periodicNode_d.{cu,cuh}`、リミタ、陰解法の fold。
- 並進周期のみの mesh はビット不変 (V5)。
- docs: `methods/boundary.md`、`discretization-median-dual-3d.md` §4.5.8。

## 8. 完了条件

- [ ] `methods/boundary.md` の「実装中」を外す
- [ ] 実装・検証完了 (§6)
- [ ] codex レビュー 2 回 (`plan` / `result`) を §6.1 に記録し、Critical / Major の採否を残作業表に反映済み
- [ ] `status` を `done` に変更し、§9 に変更ログを記載
- [ ] `plans/active/` → `plans/accepted/` へ移動、[`plans/README.md`](../README.md) を同期

## 9. 変更ログ

- `2026-09-26` — 初稿。ユーザ指示で起票。§4/§6 は `diagnostician` の設計 (既存 §4.5.8 を具体化)。
