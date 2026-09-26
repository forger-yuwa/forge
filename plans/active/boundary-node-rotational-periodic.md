# node 周期 (seam 合算経路) の回転周期対応

## メタ

- **area**: `boundary`
- **status**: `in_progress`
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

## 4. 設計方針 (2026-09-26 `diagnostician`、codex plan で改訂)

### 4.0 最優先: LSQ 勾配の seam 合算 (codex plan M4 の訂正、並進周期にも関わる)

node の勾配は `gradLSQ: 2` 固定 (`solverConfig.cpp:240-246`)。現行の `periodicGradientGather` は**片側 stencil で解いた LSQ 勾配を和にしている**
(`periodicNode_d.cu:166-212`)。Green–Gauss なら合併体積で割った部分勾配の和が合併勾配になるが、**LSQ では各側が完全な勾配を出すので、和は線形場で 2 倍**に
なる疑いがある。**まず測る** (線形場 $\phi=ax+by+cz$、seam 節点 / 内部節点の比。期待: 現行 2.0 なら欠陥、1.0 なら修正不要)。
修正するなら: 事前計算で**正規方程式 $M=\Sigma\Delta\mathbf x\Delta\mathbf x^{\mathsf T}$ を seam 越しに合算してから逆行列**を取り、各側の隣接に合併 $M^{-1}$ から重みを配る
(毎 step の gather は現行の和のまま、これで厳密に合併勾配)。回転では $M_m \to R_m^{\mathsf T} M_m R_m$、$\Delta\mathbf x$ を root 系に回す。
**並進周期の seam 勾配は意図して変わる** (修正) ので、V5 はそれを前提に書き換える。

### 4.1 角の割当 (`mesh.cpp` `buildPeriodicNodeGroups`)

- 重み付き union-find、partner 対 $(i,j)$ に $\theta_j=\theta_i+d\theta$、root は $\theta=0$。
- **`dtheta` は YAML から double で保持** (codex M3: 現行は `boundaryCond.cpp:60` で float32 に丸まり、π/2 × 4 の閉ループで 1.7e-7 rad ずれる)。
- 閉ループ検査は `remainder(Δθ, 2π)`。**double の簿記 ≤ 1e-12 rad**、**float32 格納後の直交性 ≤ 1e-6**、**座標一致 ≤ 1e-6 H** (H = セクタの径方向幅 $r_{out}-r_{in}$、ノルム = 節点の max)。
- 試験: 逆向き BC、root 交換、3 群の連鎖、60° (非 90°)。

### 4.2 変換の一覧 ($R_m$ = root 系 → member 系、x 軸まわり)

| 量 | 経路 | 変換 |
| --- | --- | --- |
| 運動量残差 | gather / broadcast | $\mathrm{res}_{root} \mathrel{+}= R_m^{\mathsf T}\mathrm{res}_m$、$\mathrm{res}_m=R_m\mathrm{res}_{root}$ |
| 保存量 state・dq | ミラー | 運動量 2 成分に $R_m$。他は恒等 |
| スカラー勾配 | gather / broadcast | $R_m^{\mathsf T}\mathbf g_m$ / $R_m\mathbf g$ (§4.0 の合併 $M$ と整合させる) |
| 速度勾配テンソル | 同 | $R_m^{\mathsf T}G_mR_m$ / $R_mGR_m^{\mathsf T}$ |
| $\nabla\cdot\mathbf u$・スカラー量 (化学種・k/ω・遷移・凝縮・受動種) | — | 恒等 (回転不変) |

### 4.3 陰解法 (codex M1 を事実訂正として採用)

node 周期の block-DPLUR は、**各部分 CV が自分側の面だけで対角と近傍和を組み、右辺は合算済み残差、各 sweep 後に root の dq を member にミラー**している
(`main.cpp:1591`、`periodicMirrorDq`)。厳密な行の fold は未実装 (`discretization-median-dual-3d.md` §4.5.7)。
**本 plan は現行構造を維持し、回転は右辺の gather / broadcast と dq ミラーの 2 か所にだけ入れる** (並進と同じ近似で、回転が新しい近似を足さない)。
行縮約 (空間対角の相似合算 + Σ $R_m^{\mathsf T}$ neighbor_accum + 時間項 1 回) は**並進の収束率改善でもある**ので別 plan に分ける (§7)。

### 4.4 リミタ — 円筒成分 bound (codex M2、今回入れる)

現行は extrema の gather に加えて `limiter_Q` 自体を成分別に min-gather している (`limiter_d.cu:70-90`)。90° 回転では $\psi_y$ と $\psi_z$ が入れ替わるので、
成分別の制限は回転と可換でない。回転周期 mesh では**横 2 成分のリミタを $(u_r, u_\theta)$ で評価**する ($u_r, u_\theta$ は周期スカラーなので extrema と
`limiter_Q` の gather がそのまま成立)。
$\nabla u_r=\cos\theta\nabla u_y+\sin\theta\nabla u_z+u_\theta\nabla\theta$、$\nabla u_\theta=-\sin\theta\nabla u_y+\cos\theta\nabla u_z-u_r\nabla\theta$、$\nabla\theta=(0,-z,y)/r^2$。
**$\psi_\perp=\min(\psi_r,\psi_\theta)$ をスカラーとして横ベクトル $(u_y,u_z)$ の再構成に掛ける** (再構成が回転と可換)。$u_x$・スカラーは現行。
キー `limiterCylindricalVelocity` (回転周期 mesh で自動 ON、他は 0)。**全周の対照も同じモードで回す**。

### 4.5 受付範囲 (起動エラー + 負例試験)

回転周期は `gradLSQ: 2` のみ。`lineImplicit`・`lowMachPrecond ≥ 1`・軸対称・軸を含むセクタは起動エラー。

### 4.6 並進の保証

「type 0 の演算経路を維持」と定義し、root 配列・合併体積・面対応を旧版と一致させる (codex m8)。ただし §4.0 の修正が入れば seam 勾配は意図して変わる。

## 5. 実装ステップ

1. codex plan 段。
2. 角の割当・閉ループ・`wall_flag` 一致・座標一致検査 (`mesh.cpp`・`mesh.hpp`)。
3. 回転の gather / broadcast / mirror / dq / 勾配 (`periodicNode_d.{cu,cuh}`)。
4. リミタ方針 (`limiter_d.cu`・`limiterPeriodic_d.cuh`)。
5. 陰解法の対角 fold の相似変換。
6. 検証ケース `case/58.annular_sector` (gmsh の 90° 環状セクタと、回転コピーの全周)。
7. §6 の run。docs。codex result。

### 5.1 残作業 (優先順)

**計算資源**: 単体・小規模はローカル、全周など重いものは AWS (ユーザ指示: ローカルで大きな計算をかけない。AWS は他セッションと共有、手動 stop しない)。

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| ~~1~~ (**済 2026-09-26**: GO-with-changes C0/M7/m2、全件採用) | codex plan 段 | | F |
| **0a** (2026-09-26 登録、plan gradient-scalar-lsq-unification §5.1 #7 から) | **node + type 1 周期は回転実装まで起動エラーにする** (最優先) | 現状は速度も勾配も回さずに走ってしまう (gradient-scalar-lsq-unification で起動**警告**のみ入れる)。あわせて NS の `periodicGradientGather_d_wrapper` (`periodicNode_d.cu:171-206`) は並進かどうかを判定しないので、回転周期では合併係数の無い片側 LSQ の和 (2 倍の疑い) になっている — 回転の実装で直す | O |
| **0** (**測定済み 2026-09-26: 2 倍を確認**、§6.2。修正は codex に諮ってから) | **LSQ seam 合算の測定と修正** (§4.0、最優先) | (1) 並進周期 node mesh に線形場を置き 1 step の勾配、seam / 内部の比を記録 (期待: 現行 2.0、修正後 1.0 ± 1e-5)。2 倍が出なければ「確認済み・修正不要」。(2) 修正: 事前計算で $M$ を seam 越しに合算 (`calcGradient_d.cu`、`mesh.cpp`)。(3) 線形場試験 PASS。**区切りで codex** | O (F レビュー) |
| 2 | 角の簿記・pairing・閉ループ・受付表 + 負例試験 | `mesh.cpp`・`mesh.hpp`・`boundaryCond.cpp` (dtheta を double で)。合格: §6 U0・U1 | O |
| 3 | 回転の gather / broadcast / ミラー / dq / 勾配 ($M$ の回転を含む) | `periodicNode_d.{cu,cuh}`。合格: V1 (純軸流 free-stream) + V2-i (固定状態の作用素比較)。**区切りで codex** | O (F レビュー) |
| 4 | 円筒成分リミタ | `limiter_d.cu`・`limiterPeriodic_d.cuh`、キー `limiterCylindricalVelocity`。全周で単独検証 (§6 V6)。**区切りで codex** | O (F レビュー) |
| 5 | 検証 run | `case/58.annular_sector` (新設)、case/39、case/09。V2-ii・V3・V4・V7 | O (結論 F) |
| 6 | docs + codex result | `methods/boundary.md`、`discretization-median-dual-3d.md` §4.5.8 | F |

## 6. 検証

**判定基準はすべて測る前に固定する。許容は既存のものを流用する。**

| # | 試験 | 合格 (測る前に固定) |
| --- | --- | --- |
| G0 | §4.0 線形場の勾配 (並進周期、回転周期) | seam 節点の勾配 / 真値 = 1 ± 1e-5 |
| U0 | 角の簿記: 逆向き BC、root 交換、3 群連鎖、60° | double ≤ 1e-12 rad、float32 直交性 ≤ 1e-6、閉ループ違反で起動エラー |
| U1 | メッシュ: `case/58` の環状セクタ (90° と 60°) と回転コピーの全周。node 変換、`check_mesh_quality` PASS、全周期ノードに partner、座標一致 ≤ 1e-6 H、`wall_flag` 群一致 | 全 PASS |
| U2 | 受付範囲の負例 (`lineImplicit`・`lowMachPrecond`・軸対称・軸を含むセクタ) | 起動エラー |
| V1 | free-stream: 純軸流 $u_x=U$ を Euler・slip で 100 step | seam 残差が並進周期の free-stream と同水準 |
| V2-i | **固定状態の空間作用素**: 旋回 (非零横運動量) を含む同じ状態で、全周とセクタの残差・勾配・面流束を比較 (新ツール `sector_vs_full.py`: 回転コピー節点の対応 ID、ベクトルは $R$ で回して比較、合併体積・壁距離・面対応の一致を先に検査) | 差 ≤ Σ\|流束\| × 1e-6 (float32 丸め) |
| V2-ii | **共通 dt** (unsteady 固定刻み) で N step 更新後の場 | 差 ≤ 反復ノイズ床 ×2 + V2-i から見積もった系統差の上限 (併記) |
| V3 | スカラー: 輸送のみ (ソース OFF) とソース込みを分け、合算後は root のみ集計。登録配列 (化学種・k/ω・凝縮・受動・遷移) の転送試験 + 非零ソースの小規模試験 | 収支誤差 ≤ Σ\|面流束\| × 1e-6 |
| V4 | SLAU + `slauWallNormalChi` (auto 1): 滑りなし内外筒、SST 低 Re。成立条件 (`wall_flag` 面 > 0、flag 差のある面 ≥ 1、seam 群一致) | 両 run 全残差 PASS (同一区間)、報告量の準定常 VERDICT、局所 $y_1^+$、全周 vs セクタで壁 $p$ の $L^\infty$ ≤ 0.5 %・$C_f$ 相対 L2 ≤ 1 %。不合格時は 2×2 (セクタ/全周 × flag 0/1) と面流束照合で切り分け (auto→0 と決め打ちしない) |
| V5 | 並進回帰: case/39 (再現入力を確保。`run_0007` は残差履歴が無いので 6000 step 再生成)、case/09 | 面流束ビット同一 (勾配以外の経路)、§4.0 修正による seam 勾配の変化は意図どおり (G0) |
| V6 | 円筒成分リミタ: 全周 (周期なし) で Cartesian モード vs 円筒モード | 差 ≤ ノイズ床 ×2、滑らかな旋回場で面値 L2 誤差が 2 次、`limiterDiag` の面再構成検査で seam 面の逸脱 0 |
| V7 | 陰解法: セクタ / 全周とも `check_convergence` PASS (同一区間) | セクタの低下桁数 ≥ 全周 − 0.5 桁 |

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| plan | `2026-09-26` | [2026-09-26-boundary-node-rotational-periodic-plan.md](../../notes/reviews/2026-09-26-boundary-node-rotational-periodic-plan.md) | **GO-with-changes**, C0/M7/m2 | **全件採用** (2026-09-26 `diagnostician` 判断)。M1 → 事実訂正 (fold は未実装)。回転は右辺 gather と dq ミラーにだけ入れ、行縮約は並進にも効く別改善として分離 (§4.3)。M2 → 円筒成分リミタを今回入れる (§4.4)。M3 → dtheta を double で、閾値を精度経路で分離 (§4.1)。**M4 → codex 案の「Green–Gauss 限定」は node が LSQ 固定なので不成立。代わりに並進周期でも LSQ seam 勾配が 2 倍になっている疑いを最優先で測定・修正 (§4.0、§5.1 #0)**。M5 → 同値性を固定状態の作用素比較と共通 dt の更新比較に分割 (V2-i/ii)。M6 → 収支を輸送のみ/ソース込みに分け root のみ集計 (V3)。M7 → 受入に全残差 PASS・準定常・y1+、V4 不合格の決め打ちをやめる。m8 → 「演算経路維持」に。m9 → case/58 |

### 6.2 結果

**G0 並進周期 (2026-09-26) — LSQ の seam 勾配は 2 倍 (欠陥を確認)**。case/09 TGV メッシュ (32³、35,937 節点、node、SLAU、三重周期) に
線形場 $U_x=10+y$ を置き 1 step。y・z の継ぎ目から 2.5 格子以上離れた節点で $\partial U_x/\partial y$: **x の継ぎ目 1458 点 平均 2.000000 [1.999990, 2.000012]**、
内部 19683 点 平均 1.000000 [0.999995, 1.000006] → **比 2.000000**。証拠 `case/09.Taylor-Green/_g0_lsq_seam/G0_translational.txt`。
**並進周期の node 計算すべてに効く** (継ぎ目の節点で LSQ 勾配が 2 倍: 粘性応力・2 次再構成・リミタ・SST 生成など)。修正方針 (§4.0) は実装前に codex に諮る。

## 7. 影響範囲

- **別 plan に分けたもの**: 陰解法の周期の行縮約 (空間対角の相似合算 + 近傍寄与の合算 + 時間項 1 回)。並進の収束率改善でもある。本 plan では未実装で、収束率は V7 で全周と比べて受け入れる。

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

- `2026-09-26` — codex plan 段 GO-with-changes (C0/M7/m2) を全件採用して改訂。陰解法の fold は未実装だった事実を訂正、円筒成分リミタを今回に、LSQ seam 勾配の 2 倍疑いを最優先で測る。ユーザ指示: 実装前・区切りごとに codex に諮る。
- `2026-09-26` — 初稿。ユーザ指示で起票。§4/§6 は `diagnostician` の設計 (既存 §4.5.8 を具体化)。
