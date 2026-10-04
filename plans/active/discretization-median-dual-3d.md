# 3D median-dual (M4) — 3次元双対メッシュ生成と periodic 双対面対応

## メタ

- **area**: `architecture` (discretization 横断)
- **status**: `in-progress`
- **related_docs**:
  - [`methods/discretization.md`](../../methods/discretization.md) (median-dual 理論・2D 実装解説)
  - [`methods/architecture/overview.md`](../../methods/architecture/overview.md)
- **related_plans**:
  - [`discretization-median-dual.md`](discretization-median-dual.md) (親計画。M1–M3 実装済、本計画が M4 を独立詳細化)
  - [`discretization-node-boundary-ghostless.md`](discretization-node-boundary-ghostless.md) (node 境界弱形式)
  - [`turbulence-iddes-sst.md`](turbulence-iddes-sst.md) (3D node DDES の最終目的)
- **created**: `2026-06-28`
- **owner**: `CFD Dev`

## 1. 目的

node-centered (median-dual) モードを **3D 混合要素 (tet/prism/hex) + spanwise periodic** に対応させ、
3D node DDES/LES を可能にする。現状 [`buildMedianDual()`](../../solver_density_cuda/mesh/gmshReader.hpp#L1305)
は `if (is3D)` で即 bail し (`dualBuilt=false`)、node モードは 2D 平面のみ。periodic は 2D dual builder にも未実装
(双対面を周期境界で対応付ける処理が無い)。本計画で 3D 双対生成と periodic 双対面対応を実装する。

**ブロッカー確認 (2026-06-28)**: SU2 三者比較で「市松はスキーム (RHS 散逸) で決まる」と判明し
([notes/investigations/backstep-lowmach-checkerboard.md](../../notes/investigations/backstep-lowmach-checkerboard.md))、
3D node では roe (内在散逸) を検証軸にする方針。その前提として 3D node 自体が必要。

## 2. スコープ

- **やる**:
  - 3D primal (tet/prism/hex/pyramid) から **一意エッジ列**を構築 (現状 3D の `planes` は primal *面*でエッジ列が無い)。
  - **3D 双対面** (エッジ {A,B} ごとのポリゴン: 各 incident cell で `edge-mid M / 面重心 F1,F2 / cell 重心 G` の
    四角パッチ集約) の面ベクトル・面積・重心。
  - **3D 双対体積** (各ノード A の sub-polyhedron: A・edge-mid・面重心・cell 重心で構成する多面体、符号付き四面体分割)。
  - **3D 境界半割面** (境界面 tri/quad を median 分割しノードへ分配、マルチマーカ corner 処理は 2D と同型)。
  - **periodic 双対面対応**: 周期境界をまたぐエッジの双対面を partner 側ノード CV と接続 (spanwise 解像乱流に必須)。
  - 検証: 双対体積総和 = primal 体積、ノード closure `Σ dS + 半割面 = 0`、**free-stream 保存** (3D 一定流で残差 0)。
  - I/O: `/DUAL` の 3D 書き出し拡張、solver 接続マップ、`output.cpp` の `Center='Node'` 3D 化、`setInitial`/probe の 3D node 化。
- **やらない**:
  - 非平面四角面の厳密分割 (Newell 法で近似、歪み大は警告)。
  - polyhedral (ieleType=0) primal の双対化 (tet/prism/hex/pyramid のみ)。
  - rotation periodic の 3D 双対 (Cartesian translation periodic を先行。rotation は後続)。
  - 移動メッシュ・AMR。

## 3. 関連 docs と前提

- 2D median-dual の現在仕様: [`methods/discretization.md`](../../methods/discretization.md)。本計画着手前に同 doc へ
  「3D 双対生成」節 (エッジ列・双対面ポリゴン・多面体体積・periodic 対応) を追記する (AGENTS 開発フロー)。
- 既存 2D 実装 [`gmshReader.hpp:1305-1590`](../../solver_density_cuda/mesh/gmshReader.hpp#L1305) のデータ構造
  (`dualVolume[nN]`, `dualCentroid[nN*3]`, `dualFaceCells[nDF*2]`, `dualFaceVect[nDF*3]`, `dualFaceArea[nDF]`,
  `dualFaceCent[nDF*3]`, 境界 `dualBcond*`/`dualBnode*`) を 3D に踏襲。2D は `nDF=nPlanes` だったが 3D は
  `nDF=一意エッジ数`。
- periodic: [`mesh::setPeriodicPartner()`](../../solver_density_cuda/mesh/mesh.cpp#L488) の node 対応付け
  (Cartesian dx/dy/dz) を双対面構築に接続する。

## 4. 設計方針

### 4.1 一意エッジ列の構築
3D 要素 (tet/prism/hex/pyramid) の局所エッジ定義テーブルから、全セルのエッジを `(min(n0,n1),max(n0,n1))` で
重複除去し `edges[]` を作る。各エッジに **incident cells** と、各 cell 内で**エッジを含む 2 面**を記録
(`edge -> {cell, faceA, faceB}` リスト)。primal 面 (`planes`) は面重心算出に流用。

### 4.2 3D 双対面 (エッジ {A,B})
各 incident cell の寄与 = 四角パッチ `(M, F1, G, F2)` (M=edge midpoint, F1/F2=エッジを含む 2 面の重心, G=cell 重心)。
面ベクトルは Newell 法 (パッチ頂点の外積総和の 1/2)、A→B 向きに統一して全 incident cell 分を加算 = 双対面ベクトル。
面積=ノルム、重心=パッチ面積加重。2D の「midpoint→centroid 区間」を 3D ポリゴンに一般化したもの。

### 4.3 3D 双対体積 (ノード A)
各 incident cell で、A 周りの sub-polyhedron を `A・(A 接続エッジの midpoint)・(A を含む面の重心)・G` を頂点とする
四面体群に分割し符号付き体積を加算。全 cell 分の総和が `dualVolume[A]`。重心 `dualCentroid` は四面体重心の体積加重。

### 4.4 3D 境界半割面
境界面 (tri/quad) を median 分割: 各構成ノード N に、`N・(N 接続エッジ midpoint)・(面重心)` のサブポリゴンの
外向き面ベクトルを分配。マルチマーカ corner (壁∩出口∩periodic) は 2D と同型に各 incident bcond へ 1 枚ずつ
(`halfByOwner`)。閉性集計 `bnodeAccum` は所有非依存で全半割面。

### 4.5 periodic 双対面対応 (核心) — 設計 (2026-06-28 詳細化)

> **方針確定 (user)**: `setPeriodicPartner` 拡張で実装し、partner ノードは**別 CV index のまま維持**する
> (CV index == primal node index の可視化 1:1 対応を壊さないため)。`nNormalPlanes` 等のループ境界変数・
> 移流/粘性カーネルのループへの影響に注意する。

#### 4.5.1 幾何的真実 (なぜ素朴な「partner 面接続」では駄目か)
周期 partner 境界ノード $N_L$ (例 $x=0$) と $N_R$ ($x=L_x$) は **period シフト $(dx,dy,dz)$ で同一視される別ノード**。
median-dual で物理的に正しい CV は**継ぎ目をまたいで両側の部分 CV を合併したもの** ($N_L$ 内側の部分 CV +
$N_R$ 内側の部分 CV、両者は周期で隣接)。継ぎ目の 2 つの境界半割面 $S_L,S_R$ は外向きが逆で**大きさが等しく相殺**
($S_L+S_R=0$) するため、**合併 CV は継ぎ目に面を持たず内部双対面だけで閉じる**。
→ 当初案「周期エッジの双対面で $N_L$ の CV ↔ partner 側 $N_R$ の CV を繋ぐ内部面を作る」は、$N_L,N_R$ が周期で
同一視されるため**退化面 (dcc≈0)** になり不正。正しいのは「**面を作る**」でなく「**CV を合併 (DOF 同一視)**」。

#### 4.5.2 現状 (ghost 半割面) が発散する理由 (run_0003 step4, データフロー調査で確定)
- 周期半割面が**境界 plane** として書かれ、`setMeshMap_d` で convective 主ループの halo 帯 (`normal_halo_planes`
  の periodic 区間) に入り、`periodic_d` が partner 実セル状態を ghost にコピーして「ghost 結合の境界面」として処理される。
- だが (a) **node の境界ノードは周期面上に乗る**ので ghost 鏡像の `dcc≈0` で退化 (`calcStructualVariables_d` は
  シフト未適用)、(b) **粘性ループは `ip<nNormalPlanes` で periodic を完全除外** (継ぎ目に粘性が入らない)。
  → 幾何退化 + 粘性欠落で発散。

#### 4.5.3 採用設計: node DOF 同一視 (master/slave, `setPeriodicPartner` 拡張) — 最小形
partner ペアの一方を master $N_m$、他方を slave $N_s$ とし、**別 index のまま「同一 DOF」として扱う**。
本質は「**残差も体積も足すだけ + 更新後ミラー**」で、毎ステップ状態 scatter や面の幾何合成は不要 (下記)。

- **(a) 周期半割面を emit しない**: `replacePrimalWithDual` で periodic マーカの bcond は半割面 plane を作らず、
  partner ノード対応リストだけ残す。→ 周期 plane は移流/粘性/halo のどのループにも現れない。
- **(b) 残差 gather (足すだけ)**: assemble 後 $\mathrm{res}_m \mathrel{+}= \mathrm{res}_s$。
- **(c) 体積も足す**: $V_m \leftarrow V(N_m)+V(N_s)$ を time integration の master に使う (これだけは省けない。
  さもないと $V_m\neq V_s$ でペアがドリフトし「同一 DOF」が崩れる)。slave は独立更新しない。
- **(d) 更新後ミラー**: $Q(N_s)\leftarrow Q(N_m)$ を更新の最後に 1 回コピー。
- **scatter は不要**: 周期 IC で $Q_m=Q_s$ から始まり、毎ステップ「和を同じだけ当てる ((b)–(d))」を守れば
  両者は自動同期し続ける ((c) の合併体積で更新が一致するため)。→ 旧 §4.5.3 の `periodicScatter_d` は廃止。
- **free-stream 保存の証明**: closure より $\mathrm{res}_L=F\cdot(-S_{L,\text{半割}})$, $\mathrm{res}_R=F\cdot(-S_{R,\text{半割}})$、
  継ぎ目で $S_L=-S_R$。よって $\mathrm{res}_L+\mathrm{res}_R = -F\cdot(S_L+S_R)=0$ → 一様流で残差和=0 が自動成立。
- **多重周期の corner/edge** (例 三重周期 box の角): ノードが複数周期境界に属する。**union-find** で全 partner を
  単一 master に縮約 (slave チェーン $N_s\to N_m\to N_{master}$ を 1 つに)。残差・体積は全メンバを master へ集約。
- 継ぎ目の移流・**粘性は両側の内部双対面 (既に `nNormalPlanes` 帯) が担う**ので、periodic を粘性ループへ
  追加する必要なし (4.5.2(b) のギャップも自動解消)。

#### 4.5.4 ループ境界変数の扱い (user の注意点)
- `nNormalPlanes`: **不変** (新規内部面を作らない。両側の内部双対面は既存)。
- 周期半割面を emit しないので、`nBPlanes` / `setMeshMap_d` の periodic 区間 / `nBoundaryHaloPlanes` から周期分が消える。
  `convPlaneBound`(node) $= n\_normal\_halo - nBoundaryHaloPlanes$ は周期 plane を含まなくなる (= 純 `nNormalPlanes`)。
  **stale な周期 plane 参照が残らないこと**を確認する (normal_halo 構築・`periodic_d` 呼び出しを node では無効化)。
- 新規追加カーネル: `periodicGather_d` (slave→master res 加算) と更新後ミラー (master→slave Q コピー) の 2 つだけ。
  time integration の体積を master で合併値にする (slave は更新スキップ or master と同値に保つ)。`periodicScatter_d` は不要 (§4.5.3)。

#### 4.5.5 代替案 (build 時 union-find 合併) — 参考
`buildMedianDual3D` で周期境界ノードを shift マッチングし **union-find で 1 CV に物理合併**して書き出せば、
出力メッシュは継ぎ目で完全内部化し **solver 側 periodic コードは一切不要**。ただし CV 数が減り
**CV↔primal node の 1:1 可視化対応が壊れる** (output で 1 CV→複数 primal node の展開が要る)。user 方針 (別 CV 維持)
では非採用だが、solver 簡潔さ最優先なら有力。

#### 4.5.6 未確定・実装時に詰める点
- gather/ミラーの実行位置 (`applyBconds` 内 periodic ブランチを node 用に差し替えるか、専用フェーズか)。
- master/slave 選定規約 (physID 小を master 等) と union-find の corner 縮約順序の決定性。

#### 4.5.7 陰解法 (block-DPLUR) の partner 行 fold
陽解法の「残差を足す」は陰解法では「**行を足す (slave 行を master 行へ fold)**」に一般化される。$Q_m=Q_s$ より
block-DPLUR の各成分が partner で加算:

| 成分 | 合併ルール |
| --- | --- |
| 対角ブロック $D$ | $D_m = D_L + D_R$ (体積項 $V/dt$ を内包 → 合併体積が自動) |
| RHS 残差 | $\mathrm{res}_m = \mathrm{res}_L + \mathrm{res}_R$ (= 陽解法と同じ) |
| 非対角 (LU sweep 和) | $LU_m = LU_L + LU_R$ (両側隣接面の寄与を union で加算) |

- partner 間の**直接結合ブロック $J(N_L,N_R)$ は存在しない** (継ぎ目に面が無い)。両者は「同じ 1 行」で、変な結合項は不要。
- **sweep ごとにミラー** $\Delta Q_s\leftarrow\Delta Q_m$ (outer sweep 後に 1 回): slave を隣接に持つセルが正しい補正を
  読めるようにする。これで DPLUR の固定点が合併方程式に一致する。
- forge の既存**壁/軸の行 decouple 機構** ([discretization-node-wall-implicit-dirichlet](discretization-node-wall-implicit-dirichlet.md),
  `axis_flag`/wall Dirichlet) の隣に periodic の行 fold を足す形。LES は陽 RK3 / dual-time 主体なので**陽解法を先に成立させ
  陰解法 fold は後続**。

> **実装状況 (2026-06-28, 実装・検証済)**: forge の block-DPLUR は対角/LU を kernel 内で cell の plane 走査で構築し、
> master/slave は別 plane を持つため別 dq が出て drift→発散する (run_0008 step64 NaN)。これを **`blockDPLURSolve` の各
> sweep 後の dq ミラー** ([`periodicMirrorDq_d_wrapper`](../../solver_density_cuda/cuda_forge/periodicNode_d.cu)、root→member の
> `dq_*_old` broadcast) で解消。**ミラーだけで安定化を達成**: implicit periodic TG が step100 完走・NaN なし・場健全
> (ρ 0.83–1.02, P/T/roe 有界、`case/35` run_0010)、非 periodic slip (run_0009) と同等のプラトー挙動になり periodic 固有の
> 発散が消えた。対角 $D_m=D_L+D_R$・LU 和の厳密 fold は**収束「率」の改善**であり安定化には不要だったため未実装 (将来精緻化)。
> ミラーは implicit 専用パスなので explicit free-stream (run_0004 機械ゼロ) は不変。

#### 4.5.8 回転 (円周) 周期 — 運動量回転 $T(d\theta)$ の挿入 (将来対応)

> **2026-09-26: 実装は [`boundary-node-rotational-periodic.md`](boundary-node-rotational-periodic.md) で行う** (ユーザ指示。スカラー・勾配・リミタ・陰解法の対角まで含めて具体化)。
回転周期では partner $N_L$ ($\theta=0$ 側) と $N_R$ ($\theta=\Delta\theta$ 側) は**別の物理位置**で、流れは軸まわり回転
$T(d\theta)$ で関係する (スカラー $\rho,\rho e,P,T,$ species は等しく、運動量だけ回転)。→ Cartesian の「同一 DOF」が
「**回転で結ばれた DOF**」になる。骨格 (§4.5.3/4.5.7) は不変で、ベクトル (運動量3成分) の授受に $T$ が挟まるだけ:

| 操作 | Cartesian | 回転周期 |
| --- | --- | --- |
| 残差 gather | $\mathrm{res}_m=\mathrm{res}_L+\mathrm{res}_R$ | $\mathrm{res}_m=\mathrm{res}_L+T^{-1}(\mathrm{res}_R)$ |
| state ミラー | $Q_s=Q_m$ | $Q_s=T(d\theta)\,Q_m$ |
| Jacobian 対角 | $D_m=D_L+D_R$ | $D_m=D_L+T^{-1}D_R\,T$ (運動量ブロックの相似変換) |

- $T$ は軸 ($x$ 想定) まわりの 2D 回転で $(\rho u_y,\rho u_z)$ に作用、$\rho u_x$・スカラーは恒等。回転は直交=内積保存
  なので継ぎ目相殺 $\mathrm{res}_L+T^{-1}(\mathrm{res}_R)=0$ は frame 整合後に成立 (法線が逆向き)。散逸も入らない。
- **既存資産流用**: cell モード `periodic_d` が既に `dtheta` で運動量を回してゴーストにコピー済 / `setPeriodicPartner`
  type=1 が回転マッチング済。この $T$ を node の gather/ミラー/Jacobian-fold に流用する。
- **検証状態が変わる**: 回転周期では Cartesian 一様流は周期的でない。free-stream 検証は**円筒成分一様** (一定軸流 +
  剛体旋回) で $\mathrm{res}_L+T^{-1}(\mathrm{res}_R)=0$ を確認する。
- **軸 ($r=0$) 干渉**: セクタが回転軸を含むと軸特異性
  ([architecture-axisym-axis-singularity](architecture-axisym-axis-singularity.md)) と二重に絡む。環状 (軸を含まない)
  passage で先に成立させる。
- **回転合成整合**: corner で複数周期/回転が連鎖するとき union-find チェーンで合成した回転が**閉ループで恒等**になること
  (正しいタイル張りの条件) をサニティチェック。単一 passage なら問題にならない。
- **順序**: Cartesian 並進周期を先に成立 → 回転周期を後続で追加。

#### 4.5.9 保存量 state ミラー (実装済・seam 圧力欠陥の真因と解消)
§4.5.3 は「残差 gather + 合併体積で master/slave は自動同期し続ける」と設計したが、これは**両者が初期から一致し
かつ更新が完全一致する場合のみ**。実際には次で desync する:
- **初期 desync**: 非周期な seed 摂動 (例: 剪断層に `Uz±3` を z 非対称に与える)、cross-mesh restart、丸め。
- 残差 gather は「**同じ res を異なる state に足す**」だけなので、初期差を**永続させる** (補正しない)。

desync すると、継ぎ目隣接面 (master の z 内面 / slave の z 内面) が **master/slave で別 state** を読み、フラックスが
不整合になる。低散逸 KEEP + 非定常で顕在化し、**継ぎ目に ~90–150 Pa の定常 P 欠陥**を生む (backstep `run_node3d_keep_wale`)。

**切り分け** (`run_node3d_keep_{slip,pure,1storder,euler}`): slip では消える=periodic 固有 / 純粋KEEPは全域 checkerboard=
Roe散逸は必須 / **1 次・Euler でも seam 残存**=2 次再構成も粘性勾配も無関係。実測で master/slave の **Uz が最大 14.9 m/s・
Uy 4.8 m/s 食い違い** (ro/roe/Ux は同期) → seed の Uz 摂動由来の desync と確定。合併体積は比=1.0 で無罪。

**解消**: NS 保存量 `(ρ,ρu_x,ρu_y,ρu_z,ρe)` を root→member ミラー (`periodicMirrorNSState_d_wrapper`、§4.5.8 表の
「state ミラー $Q_s=Q_m$」の Cartesian 版を実体化)。**(a) 初期化時** (`readValueHDF5` 後・`dependentVariables` 前)、
**(b) 各 RK stage の保存量更新直後** に呼び slave=master を強制 → DOF を真に 1 個にする。結果: **継ぎ目 P 欠陥 150→0.2 Pa、
master/slave 差→機械ゼロ、NaN なし** (`run_node3d_keep_mirror`)。陰解法は init ミラー + §4.5.7 dq ミラーで同期維持。

> 教訓: 「残差を足すだけ」では DOF を 1 個にしたことにならない。**state の同一視が必須**。回転周期 (§4.5.8) では
> ミラーを $Q_s=T(d\theta)Q_m$ にする。

### 4.6 検証量
- 双対体積総和 = primal 体積総和 (relErr < 1e-6)。
- ノード closure: `Σ_f (符号付き双対面ベクトル) + Σ 境界半割面 = 0` (periodic 内部面は両側で相殺)。
- **free-stream 保存**: 3D 一定流 (ρ,u,p 一様) で全ノード残差が機械精度。非直交・周期で特に重要。
- 負体積 0。

## 5. 実装ステップ

1. **docs 先行**: [`methods/discretization.md`](../../methods/discretization.md) に「3D 双対生成」節を追記 (§4 の設計)。
2. **エッジ列**: `gmshReader.hpp` に要素別ローカルエッジテーブル + 一意エッジ構築 + `edge→{cell,face,face}` マップ (§4.1)。
3. **3D 双対面/体積/重心**: `buildMedianDual()` の `is3D` 分岐を実装 (§4.2–4.3)。Newell 法・四面体体積ヘルパ追加。
4. **3D 境界半割面**: §4.4 を実装 (2D の `halfByOwner`/`hcentByOwner` を 3D 面サブポリゴンへ一般化)。
5. **periodic (node DOF 同一視)**: §4.5.3。`setPeriodicPartner` を拡張し partner ノード対応 (union-find) を作る。
   `replacePrimalWithDual` で周期半割面を emit しない。`periodicGather_d` (残差和) + 合併体積 + 更新後ミラーを実装
   (陽 RK3)。陰解法は §4.5.7 の partner 行 fold、回転周期は §4.5.8 の $T(d\theta)$ を後続で追加。
   `mesh.cpp` の periodic 接続と整合。
6. **検証チェック**: §4.6 を `buildMedianDual()` 末尾に追加 (体積・closure・負体積)、free-stream は別途 run。
7. **I/O・solver**: `/DUAL` 3D 書き出し、`output.cpp` Center=Node 3D、`setInitial`/`point_probes` の 3D node 化、
   `replacePrimalWithDual` の 3D 対応確認。
8. **検証 run** (§6)。

## 6. 検証

- **単体/ビルド**: native ビルド。`convertGmshToForge` を 3D tet/hex メッシュへ適用し体積総和・closure・負体積を確認。
- **検証ケース** (易→難):
  1. **3D 直方体一定流** (tet/hex): free-stream 保存 (残差機械精度) — 双対面ベクトル整合の最小確認。
  2. **3D periodic チャネル/box**: spanwise periodic で一定流保存 + 周期方向連続性。
  3. `case/05.sod_shock_tube` の 3D hex 押し出し: shock 位置・一定流が cell/解析解と整合。
  4. **`case/18.backstep` 3D (`backstep.h5`, spanwise periodic)**: roe+RANS で cell/2D node と整合 → roe+LES →
     KEEP 3D。DDES (`turbulence-iddes-sst`) は本計画完了後。
- **判定基準**: cell モード完全不変。3D node の free-stream 保存 (機械精度)・体積/closure 合格・periodic 連続。
  `verify-node-and-cell-both`・AGENTS 収束/準定常手順に従う。

## 7. 影響範囲

- ファイル: `mesh/gmshReader.hpp` (中核), `mesh/mesh.cpp` (periodic 接続), `output/output.cpp`,
  `input/setInitial.hpp`, `probe/point_probes.cpp`, `methods/discretization.md`。
- 既存: cell モード・2D node は不変 (`is3D` 分岐の新規追加なので 2D 経路ビット不変)。
- docs: `methods/discretization.md` + `methods/index.md` 整合。

## 8. 完了条件

- [ ] `methods/discretization.md` に 3D 双対生成節を追記
- [ ] 実装・検証完了 (§6: free-stream 保存・体積/closure・periodic 連続・cell 不変)
- [ ] `status` を `done`、§9 に変更ログ
- [ ] `plans/active/` → `plans/accepted/` へ移動、親 plan の M4 を done 反映、[`plans/README.md`](../README.md) 同期

## 9. 変更ログ

- `2026-06-28` — 初稿 (計画)。3D node DDES 移行のため M4 を親 plan から独立詳細化。2D 実装の構造
  ([gmshReader.hpp:1305-1590](../../solver_density_cuda/mesh/gmshReader.hpp#L1305)) を踏襲し 3D 双対面/体積/
  境界半割面/periodic 対応を設計。SU2 検証で「3D node は roe を散逸軸にする」方針を前提化。
- `2026-06-28` — **ステップ 1–4・6 (検証チェック) を実装・検証** (ブランチ `feature/median-dual-3d`)。
  - docs: [`methods/discretization.md`](../../methods/discretization.md) §2.5「3D 双対生成 (M4)」を追記
    (エッジ列・双対面・多面体体積・境界半割面・periodic・検証量)。
  - 実装: [`gmshReader.hpp`](../../solver_density_cuda/mesh/gmshReader.hpp) に `buildMedianDual3D()` を追加し
    `buildMedianDual()` の `is3D` 分岐から呼ぶ。§4.1 一意エッジ列 (要素局所面テーブル `nodesOrderPlanes` の
    周回ノード対から `edge→{cell,2面}` を導出; 別エッジテーブル不要)、§4.2 Newell 法による双対面パッチ
    `(M,F1,G,F2)` を A→B 統一集約、§4.3 四面体群 `(node,M,F,G)` による双対体積・面積加重重心、§4.4 境界面
    (tri/quad) の median 分割半割面、§4.6 体積総和/closure/負体積チェックを 3D 化。コードは要素タイプ非依存
    (tet/prism/hex/pyramid 共通パス)。
  - **付随バグ修正**: primal 体積計算が tet を `eleType.name=="tet"` で判定していたが登録名は `"tetra"` で
    不一致 → **全 tet セルの primal 体積が 0** だった既存バグを修正 (cell モードの tet メッシュも影響)。
  - 検証 (native `convertGmshToForge`, node モード):
    - **hex** (`case/35` tg.msh, 32³ box): 体積 dual=248.05=primal=248.05 (relErr 3.5e-8)、closure
      normalized=2.5e-6、境界半割面積=primal 面積、負体積 0。
    - **tet** (生成した単位 box, 197 tetra): 体積 dual=1=primal=1 (relErr 3.1e-7)、closure 9.7e-8。
    - prism/pyramid は同一コードパス (tri/quad 面プリミティブ・table 駆動) で未実行 run のみ。
  - **スコープ結果 (solver/IO 3D node = ステップ 7)**: dual-as-primary-mesh 設計のおかげで **solver 側は
    ほぼ無変更で 3D node が動く**ことを確認。非 periodic の tet box (slip 箱, Taylor-Green IC) で forge が
    10 step 完走し全保存量残差が NaN なしで低下 (`/tmp.../scratchpad`)。mesh read・dual 接続・流束・残差・dt は
    3D node で機能する。→ ステップ 7 の大半は既存機構で充足。
  - **periodic は要 §4.5 (確認済)**: 全面 periodic box (`case/35` run_0003) は **step 4 で発散**。ログ
    `conv main planes = internal 104544 + periodic 6534` の通り `setPeriodicPartner` は周期半割面を主面接続するが、
    **半割面は外向き境界パッチで CV 間を繋ぐ面でない**ため幾何が不整合 → 発散。§4.5 は周期エッジの双対面を
    「N の CV ↔ partner 側隣接ノード CV」を繋ぐ内部双対面として構築し、周期半割面から除外する必要がある
    (現状の半割面化では不可)。
  - **未着手**: §4.5 periodic 双対面 (核心・次の主タスク)、free-stream 保存 run (検証ケース 1: closure
    2.5e-6 で幾何的には既に保証、inlet/outlet box の solver run で確認予定)、output.cpp Center=Node 3D 可視化の
    確認、検証ケース 2–4。
- `2026-06-28` — **§4.5 periodic を設計詳細化 (実装は次タスク)**。cell mode periodic のデータフロー調査
  (`setPeriodicPartner`→`bint[partnerCellID]`→`periodic_d` が partner 実セル状態を毎ステップ ghost コピー、
  `calcStructualVariables_d` はシフト未適用、`viscousFlux_d` は `ip<nNormalPlanes` で **periodic を除外**) を踏まえ、
  当初案「partner 側ノード CV を繋ぐ内部双対面」は **$N_L\equiv N_R$ 同一視で退化面 (dcc≈0) になり誤り**と判明。
  正しい設計は **node DOF 同一視 (master/slave)**: 周期半割面を emit せず、master→slave 状態 scatter・slave→master
  残差 gather・合併体積で「両側部分 CV を 1 つの CV」として扱う (継ぎ目面は $S_L+S_R=0$ で消える)。継ぎ目の
  粘性も両側内部双対面が担うので粘性ループ改変不要。user 方針で partner ノードは別 CV index 維持 (可視化 1:1)。
  §4.5.1–4.5.6 に記載。実装は陽 RK3 で先行成立 → 陰解法は後続。
- `2026-06-28` — **§4.5 を最小形に整理 + 陰解法/回転を追記** (user 議論を反映)。
  - §4.5.3 最小形: 本質は「**残差も体積も足すだけ + 更新後ミラー**」。周期 IC で $Q_m=Q_s$ から始まれば
    毎ステップの合併体積更新で両者が自動同期するため、**毎ステップ状態 scatter (`periodicScatter_d`) は不要**と判明し廃止。
    残すは `periodicGather_d` (残差和) + 合併体積 + ミラー 1 回。free-stream 保存を closure と $S_L=-S_R$ から証明。
  - §4.5.7 追加 (陰解法 block-DPLUR): 陽の「残差を足す」は「**slave 行を master 行へ fold**」に一般化。対角 $D$・RHS・
    非対角 LU 和がすべて partner で加算され、partner 間直接結合ブロックは無し。sweep ごとに $\Delta Q$ をミラー。
    既存の壁/軸 行 decouple 機構の隣に置く。
  - §4.5.8 追加 (回転/円周周期, 将来対応): 骨格不変で**運動量3成分の授受に回転 $T(d\theta)$ を挟むだけ**
    (スカラーは恒等)。gather $\mathrm{res}_m=\mathrm{res}_L+T^{-1}\mathrm{res}_R$、ミラー $Q_s=T\,Q_m$、Jacobian $D_m=D_L+T^{-1}D_R T$。
    既存 `periodic_d` の `dtheta` 回転と `setPeriodicPartner` type=1 を流用。検証は円筒成分一様流。Cartesian を先に成立 → 回転は後続。
- `2026-06-28` — **§4.5 Cartesian periodic を陽 RK3 で実装・検証** (ブランチ `feature/median-dual-3d`)。
  - 実装: (1) [`mesh.cpp`](../../solver_density_cuda/mesh/mesh.cpp) `buildPeriodicNodeGroups()` — `setPeriodicPartner` の
    `partnerCellID` から周期ノードを **union-find** でグループ化 (`periodicRoot[]`)、各 group の合併体積を
    `var.c_d["volume"]` へ書き戻す。(2) [`periodicNode_d.cu`](../../solver_density_cuda/cuda_forge/periodicNode_d.cu)
    `periodicNodeGather_d_wrapper` — slave 残差を root へ atomicAdd 集約 → group 全員へ broadcast (gather+broadcast)。
    `assembleResidual` 末尾 (壁/軸射影の後) で毎ステージ呼ぶ。(3) [`convectiveFlux_d.cu`](../../solver_density_cuda/cuda_forge/convection/convectiveFlux_d.cu)
    node の `convPlaneBound` を `nNormalPlanes` にし **periodic を移流ループから除外** (継ぎ目に双対面を作らない)。
    `periodic_d` は setDT が ghost を読むため維持 (ghost を valid に保つ、flux では未使用)。`periodicScatter_d` は §4.5.3 の
    通り不要 (合併体積更新で自動同期)。
  - 検証 (native, 構造化 periodic 立方体 8³ hex 全方向 periodic, `case/35` run_0004):
    **free-stream (一様 +x 流) で全 50 step 残差が機械ゼロ** = periodic 境界 CV が gather で正しく閉じる。
    union-find `217 merged`=729−512 (triperiodic 一意 CV) で corner 8-way/edge 4-way まで完全一致。
    `conv main planes periodic 0` で除外確認。
  - **非自明流れ (Taylor-Green M0.1) は別要因で未達**: node 全面 periodic TG は step4 発散だが、**periodic を slip に
    置換しても同様に発散** (`case/35` run_0006 切り分け) → 発散は periodic 実装でなく **node+低マッハ+explicit RK3 の
    既知不安定** ([[node-mode-periodic-and-backstep-status]] backstep step18 と同類)。lowMachPrecond=2 explicit も step1 発散。
    → 非自明 periodic の安定走行は **implicit (§4.5.7 の Jacobian 行 fold) が前提**。RHS gather は実装済なので陰解法対応で
    完結する見込み。cell モード・非 periodic node はガードで不変。
  - **implicit スコープ (§4.5.7 の必要性を実証)**: node periodic TG を implicit block-DPLUR で実行 (`case/35` run_0008) すると
    step64 まで粘って NaN 発散。一方 **implicit + slip (periodic 無し, run_0009) は step99 まで NaN なし**。explicit では
    periodic/slip が同一発散 (純 node-lowMach) だったのに対し、**implicit では periodic だけ NaN** → RHS gather だけでは
    陰解法 LHS が master/slave 不整合 (合併 RHS・合併体積だが対角/非対角は未合併) で、**§4.5.7 の Jacobian 行 fold が
    implicit periodic 安定化に必要**と確定。これが非自明 periodic 実走の次の主タスク。
  - **残**: 勾配/RANS/species の periodic gather (現状 NS 保存量5のみ; 2次再構成・粘性で boundary node 勾配が片側になる)、
    §4.5.8 回転、検証ケース 3-4 (sod 3D・backstep 3D)。
- `2026-06-28` — **§4.5.7 陰解法 periodic 安定化を実装・検証** (dq ミラー)。block-DPLUR の `blockDPLURSolve` 各 sweep 後に
  [`periodicMirrorDq_d_wrapper`](../../solver_density_cuda/cuda_forge/periodicNode_d.cu) で root の補正 dq を member へ broadcast
  (`dq_block_old_*` / scalar `dq_*_old`)。master/slave が別 plane で別 dq になる drift→発散を解消。**ミラーだけで安定化達成**:
  implicit periodic TG が step100 完走・NaN なし・場健全 (`case/35` run_0010)、非 periodic slip と同等プラトー (periodic 固有発散消失,
  従来 run_0008 step64 NaN)。対角/LU の厳密 fold は収束率改善で安定化には不要のため未実装。explicit free-stream は不変 (run_0004 再確認)。
  これで **Cartesian periodic は陽・陰とも node で安定**。残は勾配/RANS/species gather・§4.5.8 回転・検証ケース 3-4。
- `2026-06-28` — **勾配の periodic gather を実装** (2次再構成・粘性の seam 精度)。Green-Gauss 勾配
  $\nabla\phi[ic]=(1/V[ic])\sum\phi_f S_f$ は合併体積 ($V_{master}=V_{slave}=V_{合併}$) のおかげで
  $\nabla\phi_{master}+\nabla\phi_{slave}=$ 真の合併勾配になる → 残差と同じ「和→broadcast」で boundary periodic node の
  片側勾配を厳密合併に直す。実装: [`calcGradient_d.cu`](../../solver_density_cuda/cuda_forge/calcGradient_d.cu) の node
  `calcGradient_b_d` ループで **periodic 半割面を除外** (勾配=内部のみ片側) + [`periodicGradientGather_d_wrapper`](../../solver_density_cuda/cuda_forge/periodicNode_d.cu)
  で 18 勾配配列 (dUxd*,dUyd*,dUzd*,drod*,dPd*,dTd*,divU) を gather。`calcGradient` 直後に呼ぶ。非軸対称限定
  (grad_volume=volume; 軸対称 A_planar は §4.5.8)。検証: free-stream 機械ゼロ維持・implicit TG 安定不変 (seam viscous が
  片側→両側へ微修正, run_0011)・非 periodic node/cell は guard で byte 不変。残: RANS(k,ω)/species 勾配の gather・§4.5.8 回転・検証ケース 3-4。
- `2026-06-28` — **勾配 gather の定量検証 + 初期 setup gather 追加** (user 指摘で判明)。各 periodic partner の per-CV 勾配は
  「**半割面を片側無視 + 体積2倍 (合併体積) → 0.5·∇φ**」で、gather (和) が 0.5+0.5=∇φ に合成する。Taylor-Green の解析勾配
  $\partial u_x/\partial x=M_0\cos x\cos y\cos z$ と照合: **境界 (合併) 勾配 mean|err|=5.3e-4 (|g|max 0.398 vs 解析 0.4)、内部 5.8e-4 と同等精度**で
  master==slave 一致 (diff 2.9e-7) — periodic 境界の勾配が内部と同等に正しいと確認 (`case/35` run_0012)。当初 res_0 が 0.5×解析に見えたのは
  **初期 calcGradient (時間ループ前) に gather を入れていなかった**ため (res_0 が片側出力)。[`main.cpp`](../../solver_density_cuda/main.cpp)
  の初期 setup calcGradient 直後にも `periodicGradientGather` を追加し res_0・初期診断も合併勾配にした。
- `2026-06-28` — **RANS k/ω periodic gather を実装 + backstep 3D 検証着手 (SST 壁関数 blocker)**。
  - 実装: [`periodicNode_d.cu`](../../solver_density_cuda/cuda_forge/periodicNode_d.cu) に k/ω 対応を追加 — 残差 gather
    (res_roK/res_roOmega を 1配列版で merged CV へ合算)・勾配 gather (dKd*/dOmegad* を keys に追加)・状態ミラー
    `periodicMirrorScalarState` (roK/roOmega を root→member、`applySSTPointImplicit` 直後; point-implicit は per-cell 対角で
    drift するため)。全 SST-guard (LESorRANS==2 && RANSmodel==1)。
  - 検証ケース backstep 3D (`case/18` run_0060): 2D planar backstep を spanwise 4H×8層 hex 押し出し、side1↔side2 Cartesian
    periodic、run_0059 収束場を 2D→3D 最近傍移植 (interp_field は node 双対 CONNE 非対応のため CELLS/MESH COORD 直読の自作)
    で restart。**periodic 機構は動作** (union-find 10997 spanwise merge、periodic を移流から除外、NS 残差 step0 有限)。
  - **blocker (periodic 非依存)**: `wallTreatmentSST=1` で rms_roK=inf at step0。**slip 置換でも同 inf** → periodic でなく
    **3D node SST 壁関数バグ** (wf_pk が 3D 壁で inf)。`wallTreatmentSST=0` では inf 消失・NaN なし完走だが coarse mesh で
    低-Re 不適。→ **次タスク = 3D node SST 壁関数修正** ([[turbulence-node-wall-function-coverage]] 系)。periodic + k/ω gather
    自体は完成・有効化済。残: SST 壁関数修正後の backstep 3D RANS x_R 照合 (2D ~7.6 と一致確認)、§4.5.8 回転、sod 3D。
- `2026-06-28` — **backstep 3D spanwise periodic RANS で periodic を完全検証 (3D node SST コーナーバグ修正後)**。
  blocker だった omega inf の真因: SST omega 壁 BC $\omega=60\nu/(\beta_1 y^2)$ で、壁ノードの実効 y を内部隣接の
  wall_dist に置換する `compute_wall_y_eff_d` が、段差凹コーナー (2壁交線) では全エッジ隣接が壁ノードで `cnt==0` に
  なり kSmall fallback → 3D は wall_dist=0 で $y=$kSmall → omega float overflow (2D は偶々 wall_dist≈0.015 で救われていた)。
  **修正**: `cnt==0` のとき 2-ring 探索 (隣接壁ノードの隣接=対角内部ノード) で第一オフ壁距離を拾う
  ([`ransBoundary_d.cu`](../../solver_density_cuda/cuda_forge/ransBoundary_d.cu) commit `da22836`)。flat 壁 (cnt>0)・cell は不変。
  **検証 (`case/18` run_0060, 修正後)**: 3D backstep RANS が安定収束 (rms_ro 1.57e-5)。**periodic 完全検証**:
  (1) spanwise で x_R=5.163–5.170 (mean 5.167, std 0.003) = 解が z 非依存 = periodic が継ぎ目を完全接続、
  (2) **3D periodic x_R 5.16 == 2D 同コード x_R 5.18 (run_0062) に 0.4% 一致**。
  → **Cartesian periodic (NS+RANS, 陽・陰) が実 RANS ケースで完全検証完了**。残: §4.5.8 回転、sod 3D、x_R の完全 settle と
  wall_y_eff 修正の SST 較正影響 (旧 7.6→5.18、exp 6.26 との比較は別軸の SST 課題)。
- `2026-06-28` — **sod 3D 検証 (検証ケース 3): periodic setup は健全、node explicit は shock 不安定で未走**。
  2D/1D sod を 3D box (200×8×8 hex, transverse y,z periodic) に押出し node periodic で検証 (`case/05` run_node3d_periodic_sod)。
  **cell periodic は完走** (ro 0.125..1.0 正常・x=0.6 で ro std=0=完全 spanwise 均一, run_sod3d_cell_ref) → periodic setup・
  メッシュ・幾何 gather は正しい。**node explicit は step16-24 で発散** だが slip/z-only/yz periodic 全て同 step16-18 発散
  (cell は完走) → periodic 非依存の **node-explicit-shock 一般不安定** (TG/backstep/entropy wave と同類)。user 指摘の
  「multi-way overlap での多重加算」は **無し** (free-stream 3 方向 periodic=8-way コーナーが機械ゼロ=幾何 gather 厳密)。
  4-way エッジの微小異常 (Ux 7e-7=内部の10倍) は float closure 誤差の蓄積 (8 半割面相殺の項数増) で論理バグでない。
  → node sod は implicit dual-time 等の安定化が要 (node-explicit-shock 限界、periodic とは別軸)。
- `2026-06-28` — **【訂正】sod 3D の node 発散は誤診だった → node sod 3D periodic は動く (検証完結)**。
  前項で「node-explicit-shock 一般不安定」とした sod node 発散は **壊れた IC が真因** (user 指摘)。`initial:"sod"` は
  非次元 (P=1) で forge CPG thermo (cp1038.8,R296.8) と非互換 → T≈0.003K で roe 破壊 → **P 均一の異常状態** (停止接触)
  に node checkerboard が乗り発散、加えて私の cfl 0.4 が過大だった。**run_0002_slau_cpg の物性・IC をそのまま使用**
  (cpg 物理単位: 左 P1e6/T2000, 右 P1e5/T400, roe=P/(γ-1), cfl0.2) すると **node sod 3D periodic が安定完走** (200 step
  NaN なし)・物理 sod (ro0.84..1.92, 衝撃右伝播)・**spanwise std=0〜1e-6 完全均一** (`case/05` run_node3d_periodic_sod)。
  node は cell より衝撃前面がクリーン (決定的 vs atomicAdd ノイズ)。→ **periodic は transient 強衝撃でも node で健全**。
  双対面ベクトルは closure 3.6e-6 で整合 (user の面ベクトル確認要望も OK)。残: §4.5.8 回転、x_R settle、SST 較正影響。
