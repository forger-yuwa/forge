# 離散化レイアウト (cell-centered / node-centered)

forge は密度ベース圧縮性 Navier–Stokes を有限体積法 (FVM) で解く。本ドキュメントは、
解変数を **どの幾何要素に乗せるか** = 制御体積 (control volume, CV) の取り方に関する理論を扱う。

- **cell-centered FVM** (現行): primal メッシュのセルを CV とし、解をセル重心に代表させる。
- **node-centered FVM** (本整備の追加対象): primal メッシュの**ノードを CV 中心**とし、各ノードの
  周囲に **中点双対 (median dual)** で CV を構成する。vertex-centered FVM とも呼ぶ。

両者は同じ保存則の積分形を異なる CV 上で離散化したものであり、**数値スキーム
(対流・勾配・粘性・陰解法) の定式化は CV と面 (CV 間の界面) の抽象だけで書ける**。
forge ではこの事実を利用して両対応を実現する (→ 本ドキュメントの「実装」節)。

本ドキュメントは理論(係数・方程式)と実装(ソース対応)をまとめる。

## 理論

### 1. 有限体積法の積分形 (CV 非依存)

保存量 $\mathbf{U}$ について、任意の CV $\Omega_i$ 上で

$$
\frac{d}{dt}\int_{\Omega_i}\mathbf{U}\,dV
+ \oint_{\partial\Omega_i}\big(\mathbf{F}_c - \mathbf{F}_v\big)\cdot\mathbf{n}\,dS
= \int_{\Omega_i}\mathbf{S}\,dV
$$

が成り立つ。$\mathbf{F}_c$ は対流フラックス、$\mathbf{F}_v$ は粘性フラックス、$\mathbf{S}$ はソース項。
離散化では

$$
V_i\frac{d\mathbf{U}_i}{dt}
= -\sum_{f\in\partial\Omega_i}\big(\mathbf{F}_c - \mathbf{F}_v\big)_f\cdot\mathbf{n}_f\,S_f
+ V_i\mathbf{S}_i
$$

となる。ここで必要な幾何量は **CV 体積 $V_i$**、**CV 代表点 $\mathbf{x}_i$**、各 CV を囲む
**面 $f$ の法線 $\mathbf{n}_f S_f$ と面重心 $\mathbf{x}_f$**、そして **面が繋ぐ 2 つの CV** だけである。
この抽象に cell / node の区別は現れない。差異はすべて「CV と面をどう作るか」に押し込められる。

### 2. cell-centered の CV と面

- CV = primal セル。$V_i$ = セル体積、$\mathbf{x}_i$ = セル重心。
- 面 = primal セル間の界面 (内部面) と境界面。面は 2 つのセルを繋ぐ。
- 境界はゴーストセル法で扱う (境界面の外側に鏡像 CV を置く)。

### 3. node-centered の CV と面 (中点双対)

primal メッシュは変えずに、その上に**双対メッシュ**を構成する。

#### 3.1 双対 CV (median dual)
ノード $p$ の CV $\Omega_p$ は、$p$ に接続する各 primal セルを
**「セル重心 — 面重心 — エッジ中点 — ノード」** の小片に細分し、$p$ に帰属する小片を集めたもの。

- 2D: 各セル内で、$p$ から出る 2 本のエッジ中点とセル重心を結ぶ四辺形が $p$ の取り分。
- 3D: セル重心・面重心・エッジ中点・ノードから成る四面体片の和。
- 性質: $\sum_p V_p = \sum_{\text{cell}} V_{\text{cell}}$ (双対体積の総和は領域体積に一致)。

#### 3.2 双対面 (edge ↔ dual face)
primal の各**エッジ** $(p,q)$ にちょうど 1 枚の双対面 $f_{pq}$ が対応する。
$f_{pq}$ はエッジ $(p,q)$ を囲むセル重心・面重心を結ぶパッチ群であり、
$\Omega_p$ と $\Omega_q$ の界面をなす。法線ベクトルはパッチの面積ベクトル和

$$
\mathbf{n}_{pq}S_{pq} = \sum_{\text{patch}}\mathbf{S}_{\text{patch}}, \qquad
\mathbf{x}_{pq} = \frac{\sum_{\text{patch}} \mathbf{S}_{\text{patch}}\,\mathbf{x}_{\text{patch}}}{\sum_{\text{patch}} \mathbf{S}_{\text{patch}}}
$$

で与える。よって **1 エッジ = 1 面 = 2 CV を繋ぐ** という cell-centered と同型の接続が得られる。

#### 3.3 閉性条件 (consistency)
各内部ノード $p$ について、$p$ を囲む全双対面の符号付き面積ベクトル和が消える:

$$
\sum_{q\in N(p)} \mathbf{n}_{pq}S_{pq} = \mathbf{0}.
$$

これは「一定場で勾配がゼロ」「一定流で対流フラックスが保存」を保証する離散整合条件で、
median dual は構成上これを厳密に満たす。実装では前処理で数値的に検証する (→ implementation.md)。

#### 3.4 境界 CV と半割双対面
境界ノードは CV 中心が**境界面上に乗る**ため、ゴースト鏡像が定義できない。
代わりに各 primal 境界面を構成ノードへ細分した **半割双対面 (half median-dual face)** を境界ノードに
割り当て、そこで境界フラックスを弱形式で課す。境界面積の分配は保存的でなければならない:

$$
\sum_{p\in f_b} S_{b,p} = S_{f_b}.
$$

### 4. cell-centered と node-centered の比較

| 観点 | cell-centered | node-centered (median dual) |
| --- | --- | --- |
| CV | primal セル | ノード周りの双対セル |
| 自由度数 (tet メッシュ) | $\approx 5.5\,N_{\text{node}}$ | $N_{\text{node}}$ (約 1/5.5) |
| 面ループ単位 | primal 面 | primal エッジ |
| 境界 | ゴーストセル | 半割双対面 + 弱形式 |
| tet 非構造での効率 | 自由度過多 | 自由度少で有利 |
| 軸上特異性 (軸対称) | セル CV が軸から離れる | 境界ノード CV を半割で扱える |

tet 非構造メッシュでは node-centered の自由度が大幅に少なく、計算コスト削減が見込める
(本整備の主動機)。一方 hex メッシュでは自由度比が逆転し得るため、両対応として残す。

### 5. 精度・安定性に関する注意

- **勾配**: Green–Gauss を双対 CV 上で適用する。§3.3 の閉性が一次精度の前提。
- **対流の高次再構成 (MUSCL)**: CV 重心の勾配と CV 重心→面重心ベクトルで外挿する形式は
  双対 CV 上でもそのまま意味を持つ。limiter も CV 値基準で不変。
- **粘性 over-relaxed 補正**: 双対では面法線とエッジ方向が概ね揃い、非直交補正項が小さく安定しやすい。
- **混合要素双対**: セル重心がノード平均近似のとき双対体積も近似を継ぐ。高歪み・非凸要素では
  負体積に注意 (前処理でチェック)。
- **★ 均質方向に 2 ノードしかないメッシュを使ってはならない (2 次精度が破綻する)**:
  下記 §5.1。

#### 5.1 2 ノード方向で MUSCL 散逸が厳密に消える (押し出し疑似 2D メッシュの禁止)

cell-centered で 2D 問題を「1 層押し出し + spanwise slip」で解くのは常套手段だが、
**node-centered では同じメッシュが spanwise に *ノード 2 枚* を作り、2 次精度が構造的に破綻する**。

z=0 のノード $i$、z=dz のノード $j$、その間の双対面 (z=dz/2) を考える。両端は slip 境界なので
GG 勾配は境界半割面 (owner-state) と内部 z 面の寄与だけで決まり、両ノードとも

$$\left(\frac{\partial\phi}{\partial z}\right)_i=\left(\frac{\partial\phi}{\partial z}\right)_j=\frac{\phi_j-\phi_i}{dz}$$

となる。これを面へ MUSCL 外挿すると

$$\phi_L=\phi_i+\frac{\partial\phi}{\partial z}\Big|_i\frac{dz}{2}=\frac{\phi_i+\phi_j}{2},\qquad
\phi_R=\phi_j-\frac{\partial\phi}{\partial z}\Big|_j\frac{dz}{2}=\frac{\phi_i+\phi_j}{2}$$

で **$\phi_L=\phi_R$ が任意の場について厳密に成立する**。face jump が恒等的に 0 なので
**Riemann/上流差分の散逸がこの方向で完全に消え、純中心差分になる**。2 ノード方向が表現できる
唯一のモードは spanwise 市松 ($\phi_i\neq\phi_j$) なので、**そのモードだけが無減衰**となり
丸め誤差から指数成長する。

重要なのは、**リミッタでは止まらない**ことである。再構成値が隣接 2 点の**算術平均**なので
常に min/max 内に収まり、Barth も Venkatakrishnan も作動しない (limiter = 1)。CFL を下げても
散逸の欠落は空間離散の性質なので効かない。LSQ 勾配も 2 点ステンシルでは同じ相殺を起こし、
GG の $f_x$ 平均によるわずかな平滑化すら失うため**さらに悪化する**。

- **運用ルール**: node で 2D 問題を解くときは**押し出しをせず平面 2D メッシュ**
  (`Physical Curve` タグ, z レベル 1) を作る。3D で均質方向を持つ場合は**最低 3 ノード
  (2 層以上)** 確保する。
- **実測**: case/26 平板で 1 層押し出し (z ノード 2 枚) は 2 次精度が step 1000 で
  spanwise 非対称 120 m/s に暴走、4 層 (z ノード 5 枚) では同 step 数で 0.065 m/s
  (**350 分の 1**)、平面 2D では 40000 step 完走。詳細と全観測の整合は
  [`case/26.flat_plate_sst/README.md`](../case/26.flat_plate_sst/README.md) の該当節。
- **将来課題**: 変換時 (`convertGmshToForge`) に「node モード × ある方向のノード数 2」を
  検出して警告/エラーにする。

### 6. 境界条件: ゴースト法 vs 弱形式 (node-centered)

#### 6.1 なぜゴースト法が node-centered で問題になるか
cell-centered では境界面の外側に鏡像ゴースト CV を置き、対流フラックスも Green–Gauss 勾配も
「ゴーストを隣接 CV とする内部面」として一括処理できる。node-centered では境界ノードが
**境界面上に乗る**ため鏡像ゴーストが退化し (§3.4)、ゴースト位置を人為的にずらす便宜が要る。
この便宜のもとでは:

- **対流フラックス**は境界面でゴースト経由 1 次精度になり、低マッハでは妥当だが、
- **Green–Gauss 勾配**は境界面値を `0.5(φ_node + φ_ghost)` で評価し、境界ノードの半割 CV と
  あいまいなゴースト幾何が混ざる。低マッハでは緩く問題ないが、**高マッハ (強膨張・衝撃) の壁近傍で
  境界ノード勾配が不正確になり 2 次 MUSCL 再構成が過大 → 発散**する (M2 で bump Mach1.65 の
  下壁後縁 x≈1.98 で実証)。

#### 6.2 弱形式境界 (weak-form boundary)
ゴーストを介さず、境界 CV の開いた面 (半割双対面) に対し**境界フラックス・境界面値を直接課す**:

- **対流**: 境界半割面 $f_b$ を通る数値フラックスを境界状態 $Q_b$ と内部状態から評価し残差へ直接加える
  $$ \mathbf{R}_p \mathrel{+}= \hat{\mathbf{F}}(\mathbf{Q}_i, \mathbf{Q}_b, \mathbf{n}_b)\, S_b. $$
  $Q_b$ (bvar) は BC 種別で定まる境界フラックス状態 (slip 鏡像・inlet 規定・outlet 特性構成) で、
  対流フラックスと壁モデル・境界面出力にのみ使う。

- **勾配 (Green–Gauss の閉曲面積分の境界寄与)**: **owner ノードの状態値をそのまま境界面値とする**
  $$ \nabla\phi_p = \frac{1}{V_p}\Big( \sum_{\text{内部面}} \phi_f \mathbf{S}_f + \sum_{\text{境界半割面}} \phi_p\, \mathbf{S}_b \Big), \qquad \phi_b := \phi_p. $$
  対流と異なり、**勾配の境界寄与に bvar を使わない** — density/velocity/pressure/temperature の全 primitive、
  wall/wall_isothermal/inlet/outlet/slip の全 non-periodic bcond で共通の原則。理由は 2 つ:
  1. **node の境界点は境界面上にある** (§3.4) — cell-centered の「セル中心は境界の内側にあり、境界値は
     セル値と異なる」という前提が node には成り立たず、owner の状態値そのものが「その空間位置での場の値」
     である。境界固有の情報 (壁温・入口全圧等) を勾配に混ぜる根拠がない。
  2. **bvar は状態と独立に動く量を含む** (例: 規定背圧・回復温度モデル $T_{aw}$) — これらを勾配へ入れると、
     状態にピンしていない量が勾配経由で内部場を駆動する意図しない結合が生まれる
     (実例: node 出口列の背圧汚染 [§7.7.1]、SST 断熱壁 $T_{aw}$ の弱閉包 [`methods/turbulence/theory.md` §6.5(f)])。
  $\phi_b=\phi_p$ は「勾配をゼロにする」のではなく、**owner の実状態を境界面代表値として発散定理の面積分を閉じる**
  という意味 (SU2 v8.5.0 の Green–Gauss も境界寄与に `field(iPoint)` を使う — 同じ設計)。
  no-slip ($u_p=0$)・等温壁 ($T_p=T_w$)・axis (半径ミラー) のように**状態が既に BC 値へ Dirichlet ピンされている**
  境界では、これは「bvar と同じ値を状態から使う」ことと数学的に同値になる (§6.3, §7.2.1)。境界条件を追加で
  課したい場合は**状態層で明示的に行う**のが原則で、勾配のジャンプで代用しない
  (symmetry/slip の接線射影のような固有の勾配補正が要る場合は別途専用の射影ロジックを設計する)。

境界に固有の情報は、勾配ではなく **状態ピン (§6.3) / 境界フラックス (bvar) / 壁モデル出力**のいずれかで表現する
のが本設計の基本原則である。実装は 本ドキュメントの「実装」節 §7 (Green-Gauss) と §7.3 (LSQ) を参照。

#### 6.3 壁の強形式 (Dirichlet) と壁ゴースト撤廃

粘性壁 (no-slip) は弱形式フラックスより**強形式 (Dirichlet) の方が node-centered では自然かつ堅牢**である。
理由は §3.4 の通り**境界ノードが壁面そのものの上に乗る**ため、その点の速度は物理的に厳密にゼロでなければ
ならない:

$$ \mathbf{u}_p = \mathbf{0} \qquad (p \in \text{wall}). $$

cell-centered ではセル中心が壁から半セル内側にあり、壁面での $\mathbf{u}=0$ は鏡像ゴースト
($\mathbf{u}_{\text{ghost}}=-\mathbf{u}_{\text{interior}}$) で**面平均**を 0 にして表す。これはセル中心速度が非ゼロ
(近壁速度) であることと整合する。しかし node-centered で同じ鏡像を使うと、面平均は 0 になるが
**壁ノード自身の速度は 0 に固定されない** (残差で更新される)。これは壁ノードが壁上にあるという幾何と矛盾する。

**壁ゴーストは不要**:
- **断熱壁を通る対流フラックスは恒等的に 0**。$\mathbf{u}\cdot\mathbf{n}=0$ より質量・エネルギー流束 0、運動量は
  Dirichlet で固定 (壁反力が圧力を吸収)、断熱で熱流束 0。よって壁半割面の寄与は丸ごと不要で、CV を
  「壁側で閉じない」まま扱ってよい (落とした圧力起因の偽運動量は運動量 Dirichlet が吸収する)。
- **壁せん断の物理は内部双対面が担う**。近壁の内部ノード $q$ が感じる粘性せん断 $\mu\,\partial u_t/\partial n$ は、
  $q$ と壁ノード $p$ ($\mathbf{u}_p=0$) を繋ぐ**通常の内部双対面** $f_{pq}$ の粘性フラックスとして自然に現れる。
  したがって壁ゴーストを撤廃しても近壁せん断は失われない (むしろ鏡像で実効距離が 2 倍になる歪みが消える)。

> **実装上の補足 (検証で更新)**: 残差射影 Dirichlet は壁の**圧力寄与も落とすため不正**で発散した。また壁ゴーストを
> 完全に消すと壁ノードの Green-Gauss 勾配が閉じず悪化した。よって実際の採用形は: **壁ゴーストは勾配/粘性の
> ミラー用に残し、対流のみ壁境界値 (bvar 速度=0) で pressure-only に弱形式化**する (no-slip の本質=対流の運動量は
> 圧力のみ・せん断は内部双対面が担う、を満たす)。さらに**近壁極小 CV の粘性は陰解法では十分 implicit 化されず
> explicit (cfl≤0.1) が安定**。詳細・検証は 本ドキュメントの「実装」節 §7.2。

##### 6.3.1 陰解法での Dirichlet 整合 (Jacobian 行 decouple)

壁 $\mathbf{u}_p=\mathbf{0}$ を**残差射影だけ** ($\mathbf{res}_{\rho\mathbf{u}}=0$) で課すと、explicit では $\Delta(\rho\mathbf{u})=0$ になり厳密に効くが、**implicit (block-DPLUR) では不十分**である。block-DPLUR の対角 5×5 ブロックは壁ノードの運動量を連続・エネルギーと連成したままなので、線形 solve は残差 0 でも $\Delta(\rho\mathbf{u})\neq0$ を返し、**壁速度が drift する** (実測: 壁 $|\rho u_x|$ が時間発展で増大、状態射影 `enforceWallNoSlip` が毎ステージ余分な KE を剥ぐ持続的エネルギーシンク化)。

SU2 (同じ vertex-centered median-dual) はこれを **`Jacobian.DeleteValsRowi`** で解決する: 壁ノードの**運動量行を単位行に置換** ($\Delta(\rho u)=\Delta(\rho v)=\Delta(\rho w)=0$ を solve から直接得る)。連続・エネルギー行は保存式のまま残し、$u=0$ より $\rho e = \rho e_\text{int}$、圧力は EOS が復元する (CPG: $P=(\gamma-1)\rho e$、TP: $\rho e_\text{int}\to T\to P=\rho R_\text{mix}T$)。decouple は運動量行のみを触り thermo 行に依らないので **CPG/TP 共通** に動く。forge では既存の軸 `axis_flag` 行 decouple 機構を**壁の運動量3行へ拡張**して実装する。計画: [`discretization-node-wall-implicit-dirichlet.md`](../plans/accepted/discretization-node-wall-implicit-dirichlet.md)。

> 軸 (§7.1) の `roUy` 単行 decouple は r=0 の極小 CV が多方程式的に ill-posed なため単独では発散したが、**壁は r=0 特異点を持たない**ので、運動量3行 decouple が drift 機構を直接除去できる。

#### 6.4 コーナー所有優先 (boundary∩boundary)

2 つの境界に同時に乗るノード (例: 壁∩出口の出口リップ) は、ゴースト法では**境界ごとに別々のゴースト**が
立ち、同一 CV の残差に矛盾する条件 (壁の $\mathbf{u}=0$ と出口の超音速外挿) が同時に課されて発散する。

これを避けるため、各境界ノードを**優先度に従い 1 つの境界条件にのみ所有させる**:

$$ \text{wall} > \text{inlet} > \text{outlet} > \text{slip/axis}. $$

物理的に、壁とその他境界の交線上の点は**壁点**である (リップは固体壁)。よって壁が最優先で所有し、当該ノードは
壁 Dirichlet として扱い、他境界 (出口等) の条件・半割面からは除外する。これによりコーナーの矛盾が解消する。

> 段階導入: まず壁を強形式化＋壁ゴースト撤廃 (Phase 1)、続いて流入出/slip/axis のゴーストも半割面弱形式へ
> 置き換え node モードを完全 ghost-free にする (Phase 2)。

### 参考

- median-dual / edge-based FVM は非構造ソルバ (SU2, Fluent の node-based など) で広く使われる定式化。
- 実装上の対応とデータ構造は 本ドキュメントの「実装」節、
  実装計画は [`.github/plans/discretization-median-dual.md`](../plans/active/discretization-median-dual.md) を参照。

## 実装

### 1. 設計の核: discretization mesh 抽象

forge の GPU カーネル群 (対流・勾配・粘性・block-DPLUR) は、消費する幾何が次の 4 群に限られる:

1. **CV 中心量**: `volume[ic]`, `ccx/ccy/ccz[ic]` ([variables.cpp](../solver_density_cuda/variables.cpp))
2. **面量**: `sx/sy/sz[ip]`, `ss[ip]`, `pcx/pcy/pcz[ip]`, `fx[ip]`, `dcc[ip]`
3. **面 → CV**: `map_plane_cells_d[2*ip]` ([mesh.cpp](../solver_density_cuda/mesh/mesh.cpp))
4. **CV → 面 CSR**: `map_cell_planes_index_d` / `map_cell_planes_d`

これらは「CV と面」の抽象だけで、primal セル形状を参照しない。`fx`/`dcc` も CV 重心と
面重心からランタイム再計算される ([calcStructualVariables_d.cu](../solver_density_cuda/cuda_forge/calcStructualVariables_d.cu))。

→ **median-dual を前処理で構築し、同じ `nCells/nPlanes/map_plane_cells/CSR` レイアウトに
「CV=ノード, 面=エッジ双対面」として詰めれば、内部カーネルは無変更で動く。**

`mesh` 構造体 ([mesh.hpp](../solver_density_cuda/mesh/mesh.hpp)) は既にこの汎用抽象になっている。
よって solver はこの抽象だけを消費し、**供給側 (mesh を作る前処理) を 2 実装に分ける**:

- **cell モード**: 現行 gmshReader 経路。`/CELLS` を CV とする。**現状維持**。
- **node モード**: 新 `dualMeshBuilder` が双対量を HDF5 `/DUAL` グループに書き、
  `readMesh` がそこから `nCells := nNodes`, `nCells_ghst := 0` として読む。

切替は `solverConfig` の `discretization: cell | node` (既定 `cell`)。
**`nCells` を `nCV` にリネームしない** (差分肥大回避)。node モードでは `nCells == nNodes` と読み替える
(その旨を [mesh.hpp](../solver_density_cuda/mesh/mesh.hpp) の doc コメントに明記)。

### 2. 双対メッシュ生成 (前処理)

配置: `convertGmshToForge` 前処理 ([gmshReader.hpp](../solver_density_cuda/mesh/gmshReader.hpp))。
既存の primal 幾何構築 (セル重心・面重心・向き合わせ) と、**既に構築済みだが未使用の
`node.iCells` / `node.iPlanes`** ([mesh.hpp:19-20](../solver_density_cuda/mesh/mesh.hpp)) を起点に再利用する。

生成物:
- **双対体積** $V_p$: 各 primal セルを「セル重心–面重心–エッジ中点–ノード」で細分し各ノードへ集計。
  `Σ双対体積 == Σprimal体積` を assert。
- **双対面** (primal エッジ ↔ 1 枚): `surfVect = Σ パッチ面積ベクトル`, `surfArea = |Σ|`,
  `centCoords = 面積加重重心`。`plane_cells = {n0, n1}`、法線は n0→n1 向きに統一。
- **境界半割双対面**: 各 primal 境界面を構成ノードへ分配し `bnode_face_{sx,sy,sz,ss,pc}` を格納。
  `Σ半割面積 == primal face 面積` を assert。
- **閉性チェック**: 各内部ノードで `Σ(双対面ベクトル, 符号付き)` の max ノルムを print。
- 出力: HDF5 `/DUAL` グループ (可視化用 primal は `/MESH /CELLS` に残す二層構造)。

### 2.5. 3D 双対生成 (M4)

2D では primal **エッジ = `plane`** だったため双対面はそのまま `nPlanes` 枚だが、3D では `plane` は
primal **面 (tri/quad)** であり、双対面の単位となる**エッジ列が陽に存在しない**。よって 3D 双対生成は
**一意エッジ列の構築**から始める。配置は 2D と同じ [`buildMedianDual()`](../solver_density_cuda/mesh/gmshReader.hpp)
の `is3D` 分岐。計画: [`discretization-median-dual-3d.md`](../plans/active/discretization-median-dual-3d.md)。

#### 2.5.1 一意エッジ列と `edge → {cell, face, face}`
要素 (tet/prism/hex/pyramid) の各局所面 (`elementType.hpp` の `nodesOrderPlanes`) は、その周回ノード対が
そのままエッジである。全セル・全面の周回エッジを正規化キー $(\min(a,b),\max(a,b))$ で重複除去し `edges[]` を作る。
**1 エッジは各 incident cell 内でちょうど 2 つの面に属する** (多面体エッジの定義) ので、各 (edge, cell) に対し
その 2 面 (`F1,F2`) を `edge → {cell, faceA, faceB}` として記録する。これは別途エッジテーブルを持たず、
既存 primal 面 (`planes` の重心 `centCoords`) をそのまま面重心 $F$ に流用できる利点がある。3D では
$n_{DF} = $ 一意エッジ数 (2D の $n_{DF}=n_{Planes}$ に相当)。

#### 2.5.2 3D 双対面 (エッジ $\{A,B\}$)
各 incident cell の寄与は四角パッチ $(M, F_1, G, F_2)$ — $M$=エッジ中点, $F_1,F_2$=エッジを含む 2 面の重心,
$G$=セル重心。面ベクトルは **Newell 法** (パッチ頂点の外積総和の $1/2$) で非平面四角でも一貫に求め、$A\to B$ 向きに
統一して全 incident cell 分を加算したものが双対面ベクトル `dualFaceVect`。面積 `dualFaceArea` $=$ ノルム、
重心 `dualFaceCent` $=$ パッチ面積加重。これは 2D の「midpoint→centroid 区間 (単位厚みで rotate(-90))」を
3D ポリゴンに一般化したものに相当する。

#### 2.5.3 3D 双対体積 (ノード $A$)
各 incident cell で $A$ 周りの sub-polyhedron を、$A$・($A$ に接続するエッジの中点)・($A$ を含む面の重心)・$G$ を
頂点とする**四面体群**へ分割し、符号付き体積 $\frac{1}{6}(\mathbf{r}_1-\mathbf{r}_0)\cdot[(\mathbf{r}_2-\mathbf{r}_0)\times(\mathbf{r}_3-\mathbf{r}_0)]$
を加算する。全 cell 分の総和が `dualVolume[A]`、重心 `dualCentroid` は四面体重心の体積加重。
$\sum_A V_A == \sum_{\text{cell}} V_{\text{cell}}$ を検証する。

#### 2.5.4 3D 境界半割面
境界面 (tri/quad) を median 分割し、各構成ノード $N$ へ $N$・($N$ 接続エッジ中点)・(面重心) のサブポリゴンの
外向き面ベクトルを分配する。マルチマーカ corner (壁∩出口∩periodic 等) は 2D と同型に各 incident bcond へ
1 枚ずつ (`halfByOwner`)、閉性集計 `bnodeAccum` は所有非依存で全半割面を集計する。

#### 2.5.6 双対幾何の数値精度 (ローカル原点 + double)

節点座標は `geom_float` (float32) 格納であり丸め ~1e-9 m 級を持つ。Newell 法は
$(a-b)(a+b)$ 形で**差分 (パッチ寸法) に絶対座標の和 (領域寸法) が掛かる**ため、絶対座標のまま
評価すると float32 丸めが領域寸法で増幅され、薄壁スリバーパッチ (例: 100 μm × 1 μm 級) では
面ベクトルの数十 % 誤差・向き判定の誤反転になり得る (2D shoelace 桁落ちで実害が出たのと同根 —
壁 CV 重心が CV 外へ出て `wall_dist` を破壊、plan boundary-node-nozzle-wall-outlet-stability §2.9)。
このため **内部双対面パッチはエッジ中点 $M$ 相対、境界半割面サブ四角はノード $N$ 相対の
ローカル原点で Newell を評価**し (平行移動不変なので厳密演算では同値)、境界半割面の蓄積
(`bnodeAccum`/`halfByOwner`/`hcentByOwner`) も double で行って最後に `geom_float` へ cast する。
`tetVol`・双対重心 (正重み加重和) は元から局所差分 double で桁落ちしない。
計画: [`architecture-median-dual-3d-double-geometry.md`](../plans/accepted/architecture-median-dual-3d-double-geometry.md)。

**本規約は 2D 側の全幾何量にも適用する** (2026-08-31 完了)。2D で先行修正されたのは
sub-CV shoelace (体積・重心) のみで、**2D 双対面ベクトル** (`Mx`/`sx=G−M`/rotate/向き整合内積) と
**2D 境界半割面** (`planes[].surfVect` 半分の float 蓄積) が絶対座標 float32 のまま残っていた。
実害: 第一セル ~2.4 μm (M6 NS y+1 級, 全長 6.3 m) で**向き整合の内積 $\vec n\cdot\vec e_{AB}$ が
float で誤反転**し面 1 枚が逆符号で入り、閉性 normalized 0.099・壁 bcond 面積 3 % 欠損で
変換拒否 (`dual faces not closed`)。修正 = ① makeMesh の CW 判定を先頭ノード相対 + double 化 (真犯人)、② 面ベクトルを
エッジ中点 $M$ 相対 + double で評価 (向き内積も double)、③ 境界半割面は節点座標から
double で再構成 (向きは整向済み `surfVect` の符号に合わせる) し蓄積も double。
さらに **makeMesh の整向判定も同規約**: 2D CW 判定 (shoelace) は先頭ノード相対 + double
(float 絶対座標では積の丸めがスリバー符号付き面積を 8 桁上回り境界 surfVect が誤反転する —
これが閉性破綻の真犯人だった)、3D の整向内積 `Db`/`D` も面先頭ノード相対 + double で
面重心・セル重心を再構成して評価する (2026-08-31、健全メッシュではビット同一)。
計画: [`discretization-median-dual-2d-facevect-precision.md`](../plans/accepted/discretization-median-dual-2d-facevect-precision.md)。

#### 2.5.5 periodic 双対面対応
周期境界面上のエッジは partner 面側に同形のエッジが存在する。`setPeriodicPartner` のノード対応 (Cartesian
$dx,dy,dz$ 平行移動) で $A\leftrightarrow A'$, $B\leftrightarrow B'$ を対応付け、**周期エッジの双対面を境界半割面に
せず、両側 incident cell を直接繋ぐ内部双対面**として扱う (周期 ghost ではなく直接接続、cell モードの
`periodicPartner` と同思想)。これにより spanwise の解像乱流が周期方向に連続する。閉性は周期内部面が両側で
相殺するため通常の内部双対面と同様に成立する。

#### 2.5.6 検証量 (3D)
$\sum_A V_A == \sum V_{\text{cell}}$ (relErr $<10^{-6}$)、ノード closure $\sum_f \mathbf{S}_f + \sum$ 境界半割面 $= 0$、
**3D 一定流 (free-stream) で全ノード残差が機械精度**、負体積 0。非平面四角面は Newell 法で近似 (歪み大は警告)、
polyhedral primal の双対化は対象外。

### 3. カーネル別の対応状況

| カテゴリ | ソース | node 対応 |
| --- | --- | --- |
| 対流 (Roe/SLAU/AUSM) MUSCL | [convectiveFlux_d.cu](../solver_density_cuda/cuda_forge/convectiveFlux_d.cu) | 無変更で成立 |
| Green–Gauss 勾配 | [calcGradient_d.cu](../solver_density_cuda/cuda_forge/calcGradient_d.cu) | 無変更 (閉性が前提) |
| 粘性 over-relaxed | [viscousFlux_d.cu](../solver_density_cuda/cuda_forge/viscousFlux_d.cu) | 無変更 |
| block-DPLUR 陰解法 | [timeIntegration_d.cu](../solver_density_cuda/cuda_forge/timeIntegration_d.cu) | 無変更 |
| 構造量 (fx/dcc) | [calcStructualVariables_d.cu](../solver_density_cuda/cuda_forge/calcStructualVariables_d.cu) | 無変更 (CV 重心と面重心から再計算) |

### 3.5. M1 実装方針 (確定): 双対を primary mesh として書き出す

設計検討の結果、**M1 は solver 側を一切変更せずに実現できる**ことが判明した (当初計画の
`readMesh` dual 分岐や新 BC カーネルは不要だった)。具体的には:

- `convertGmshToForge` が node モードで `buildMedianDual()` → `replacePrimalWithDual()` を実行し、
  双対メッシュ (CV=ノード、内部面=双対面、境界面=半割双対面) を **そのまま `/MESH /PLANES /CELLS
  /BCONDS /VALUE` に primary mesh として書き出す** ([gmshReader.hpp](../solver_density_cuda/mesh/gmshReader.hpp))。
- solver の `readMesh` はそれを通常メッシュとして読み、**境界半割面 (1 セルの境界 plane) から既存ロジックで
  自動的にゴースト CV を生成する**。よって `setMeshMap_d`・構造量カーネル・対流フラックス・既存 BC
  カーネル (`wall_d`/`inlet_*`/`outlet_*`/`slip_d`) はすべて無変更で node-centered を扱える。
- **境界はゴースト経由の弱形式**: 各境界ノードの半割面に対しゴースト CV を置き、既存 BC カーネルが
  ゴースト状態を設定 → 対流フラックスがゴースト経由で境界フラックスを与える (cell-centered と同型)。
  ゴースト位置は `pc = node + h·n_out` (h=0.5√双対体積) で非退化にする (M1 は非粘性 1 次のため位置は
  フラックスに影響しない)。**より厳密な弱形式境界カーネルは M2+ の精度改善として残す**。
- **可視化 (対応済み)**: node モードでは `replacePrimalWithDual()` が primal セル接続を `/VIZMESH/CONNE`
  に退避し、`writeInputH5`・`output.cpp` がそれを使って **primal セルトポロジ + `Center='Node'`** で
  XDMF を書く (CV index == primal node index)。cell モードは従来どおり `Center='Cell'`。
  両モードとも自己整合な XDMF/HDF5 を出力し ParaView で読める。

### 4. 書き換えが必要な箇所 (M2 以降)

- **境界条件** ([boundaryCond_d.cu](../solver_density_cuda/cuda_forge/boundaryCond_d.cu)):
  現状はゴーストセル法 (`bplane_cell` / `bplane_cell_ghst` に対称/反射値)。node では境界ノードが
  境界上に乗りゴーストが置けない。新 `boundaryFlux_node_d.cu` で半割双対面に弱形式フラックスを
  `atomicAdd(res_*[node])`。ディスパッチは [boundaryCond.cpp](../solver_density_cuda/boundaryCond.cpp) でモード分岐。
- **readMesh のゴースト生成** ([mesh.cpp:401-455](../solver_density_cuda/mesh/mesh.cpp)):
  node モードでスキップし、`nCells_ghst = 0`、境界ノード + 半割面を構築。
- **I/O** ([output.cpp](../solver_density_cuda/output/output.cpp)): node モードで XDMF
  `Center='Cell'`→`'Node'`、`Dimensions=nCells`→`nNodes`。Topology/Geometry は可視化用 primal を流用。
- **初期条件** ([setInitial.hpp](../solver_density_cuda/input/setInitial.hpp)): 領域判定ループを node 座標へ。
- **プローブ** ([point_probes.cpp](../solver_density_cuda/probe/point_probes.cpp)): KD-tree を CV 中心=ノード座標で構築。
- **implicit/AMG の edge 分類** ([mesh.cpp:742-780](../solver_density_cuda/mesh/mesh.cpp)):
  `nNormalPlanes` ベースの CV 隣接で内部 edge と境界 edge を厳密に分類。
- **軸対称の軸 (R=0) 注意点 (対応済み)**: solver は回転体積を `volume = A_planar × ccy`(ccy=セル重心 R)で
  実行時計算する。node-centered で CV 重心にノード座標 (軸上で R=0) を使うと**軸ノードの回転体積が ≈0 → 時間積分が
  除算爆発**する (case/29 conical で 272 軸 CV が体積比 1.5e20、step3 で NaN)。対策: `replacePrimalWithDual` が
  dual セル重心に**双対 CV の面積加重重心** (`dualCentroid`, 軸上ノードでも R>0) を使う。これで回転体積が正しくなり
  (体積比 2.4e4, cell の 1.2e4 と同等)、平面ケースは無回帰 (内部ノードでは重心≈ノード)。FV 的にも CV 重心が正しい中心。
  **ただし残課題**: これは軸上のゼロ体積爆発を解消するだけで、node-centered の near-axis 精度問題は別途残る。
  case/29 conical Euler の収束場で、node は軸-inlet 角に非物理 P オーバーシュート (P>chamber Pt)・軸-exit に Mach
  スパイクが出る (bulk は cell と一致)。軸 BC・inlet 角・r 重みの near-axis 処理は今後の課題 ([architecture-axisym-axis-singularity] と関連)。
  **試行 (いずれも発散→不採用)**: (1) 軸ノードで `roUy=0` 毎ステージ強制 (roe から半径 KE 除去) → 過渡で roe<0→全域 NaN。
  (2) 軸ノードで半径方向圧力ソース抑制 → ソースは軸 CV の釣り合いに load-bearing で step~51 破綻。
  **結論: 半径ソースは軸 CV にも必要で安易な対称強制は不可**。baseline (ソース維持) は収束し corner オーバーシュートのみ残る。
  infra (`mesh.axis_flag_d`, `enforceAxisSymmetry_d`) は残置・既定 off。corner の正攻法は別途要検討。
  **SU2 流の確認 (ユーザ提案)**: SU2 はノードベース (中点双対) で、軸を `MARKER_SYM=(axis)` の**対称面**として扱い、
  対称ノードで**法線 (=半径) 運動量成分を残差から射影**する。forge で同方式を試したが、forge の implicit は
  **block-DPLUR (5×5 連成)** のため:
  - (3) `res_roUy=0` 射影 (全 flux+source 後) → 単独では連成 solve が補正を漏らし**無効** (Uy 164→136 でほぼ不変)。
  - (4) implicit commit で `roUy=0` 強制 → 線形 solve と不整合で**発散級に悪化 (Mach~1000)**。
  → **SU2 流対称面を効かせるには block-DPLUR の Jacobian row を軸ノードで整合修正する必要**があり (SU2 は
  Jacobian を修正している)、commit/残差だけの hack では不可。これは delicate な implicit-solver 改修なので
  open issue とし、infra (`zeroAxisRadialResidual_d`、commit の `axis_flag` 引数) は残置・既定無効。
  baseline (無介入) は収束し corner 1 セルのオーバーシュートのみ。
  **explicit 検証 (ユーザ要請)**: conical Euler を explicit で回すと **node は step1 で即発散** (cfl=0.1 でも)。
  原因は **near-axis の explicit 剛性**: 軸 CV は r 重みで面積→0 のため spectral radius→0、explicit 更新
  `res·CFL/spectral` が発散 (cell は step3483 まで持つが非粘性ノズル startup でやはり発散)。よって本ノズルは
  **explicit では node が回らず implicit (block-DPLUR 無条件安定) が前提**。explicit 検証には near-axis の
  dt/spectral フロアが別途要る。射影 (roUy 残差=0) は今 node+axisym で既定 ON だが、explicit が回らない・implicit は
  連成で漏れるため、単独では corner を直せない。
  **入口境界の確認 (ユーザ要請)**: 収束場の入口面 (x=inlet) を r 方向に見ると **r≳0.003 で入口 BC は正常**
  (P≈4MPa=chamber Pt・軸方向流入正)。異常は **軸-inlet 角 (r<0.002) のみ**: 軸方向**逆流 Ux≈-798**・Uy≈+164・P 6.82MPa。
  → **入口 BC 自体は問題なし。軸特異性が角を多方程式的に汚染** (radial だけでなく axial 逆流・P も)。roUy 単独射影では
  直らず、near-axis は全方程式整合 (対称面 + dt/spectral フロア + 連成 Jacobian) の包括対策が要る。
  **SU2 処方の全試行 (ユーザ提案: y<EPS で軸ソース OFF + roUy/Uy=0)**: 両処理を実装 (`axisymmetricSource` の
  `axis_flag` でソース OFF、`enforceAxisSymmetry_d` で roUy=0/Uy=0/roe 整合) し、**併用かつ全配置** (applyBconds
  =更新前、block-DPLUR commit 内、commit 後 post-update) で試したが、**implicit は全て Mach~1000 に発散**。
  理由: **block-DPLUR (defect-correction) は solve の外で roUy/roe を手術されると非整合**になり増幅する。SU2 が安定なのは
  **対称条件を Jacobian (行) の中で課す**から。explicit は recipe 併用でも軸 CV が **step1 で発散** (別系統の near-axis
  剛性)。**結論: SU2 処方は原理的に正しいが、forge では「外部状態手術」では不可。block-DPLUR の軸ノード 5×5 Jacobian で
  roUy 行を decouple する実装 (= SU2 と同じ内挿) が必須**。全フック (axis_flag/enforceAxisSymmetry/zeroAxisRadialResidual)
  は実装済み・既定無効で残置。baseline (無介入・ソース ON) は収束し corner 1 セルのみオーバーシュート。
  **in-Jacobian decouple も試行 (SU2 と同じ「Jacobian 内で対称化」)**: block-DPLUR の軸ノード 5×5 で roUy 行 (index2) を
  単位行に置換し rhs[2]=0 (`timeIntegration_d.cu` の `axis_flag`) → solve の中で一貫して dq_roUy=0。結果:
  **bulk は正しい (mean Mach 1.883) が、corner 2 セル (inlet-axis, exit-axis) が Mach~1007 に発散**(P 22MPa)。
  **重要: corner の発散は roUy ではなく軸方向運動量・密度**(Ux 逆流 −798 と同系)。**= corner CV (r=0 で両境界半割面が
  r 重み→0 の極小 CV) は多方程式的に ill-posed で、roUy の対称化だけ (外部手術でも in-Jacobian でも) では直らない**。
  → corner は構造的処置が要る可能性 (例: 軸-境界角ノードを隣接へ merge/除外する、専用の角クロージャ)。
  baseline (無介入) が現状ベスト (収束・bulk 正しい・corner は P≤6.8MPa の局所オーバーシュートに留まる)。全フックは既定無効。
  **SU2 メッシュ/角処理の調査結果 (ユーザ要請, `run_su2cmp_su2_euler`)**:
  - SU2 メッシュ = forge node メッシュと**同一** (NPOIN=27472、軸-inlet 角は**同じ共有ノード id 2**, x=-0.0587/y=0、
    marker は inlet101/axis272/wall272/outlet101 で forge bcond と一致)。**トポロジは差が無く、差は BC 処理のみ**。
  - **SU2 の共有角は clean**: 角で Ux=+54.78 (正常流入, 逆流でない)・**Uy=0 厳密**・P=3.99MPa (=chamber, overshoot 無し)。
    SU2 は**軸全体で Uy=0 を厳密強制** (max|Uy|=0)。
  - forge の roUy 対称化 (外部手術 / in-Jacobian 単行 decouple) は**全て発散**。SU2 が clean なのは、SU2 の対称面が
    **残差 + Jacobian + 解を一体で扱う coherent な対称面演算子 (CSymmetryPlane 相当)** だから。forge の piecemeal な
    1 行 decouple は、軸 CV の roUy 連成が持つ陰的減衰を壊し、軸方向・密度を不安定化させる (Mach~1000)。
  - **結論**: corner を SU2 並に clean にするには、**SU2 の対称面演算子を一体移植**する必要がある:
    (i) 対称面フラックスを圧力のみ (質量/エネルギー流束 0)、(ii) 残差から法線運動量成分を除去、(iii) 同じ射影を
    block-DPLUR の 5×5 Jacobian へ整合適用 (1 行 hack でなく射影行列 P=I-nn^T で両側変換)、(iv) 解にも v_n=0。
    SU2 ソースは本リポジトリにバイナリのみで未読。bulk は完成済みのため残作業は corner=この対称面演算子の移植。
- **軸対称** (副次, [variables.cpp](../solver_density_cuda/variables.cpp) /
  [axisymmetricSource_d.cu](../solver_density_cuda/cuda_forge/axisymmetricSource_d.cu)):
  双対体積に `r_node` 重み、`r_eff = volume/A_planar` を双対で再定義、軸上半割マスクをノード版に移植。

### 7. M3: node-centered 弱形式境界 (実装方針)

理論は 本ドキュメントの「理論」節 §6。**cell-centered のゴースト経路 ([boundaryCond_d.cu](../solver_density_cuda/cuda_forge/boundaryCond_d.cu))
は一切変更せず**、node 専用の弱形式を**別ファイル** `cuda_forge/boundaryNode_d.cu` (+ `.cuh`) に実装する。

**設計判断 (既存ロジックの流用)**: 各 BC 種別の境界状態 $Q_b$ は、既存 BC カーネルが**従来どおり
ゴースト CV に書き込む値**を流用する (slip 鏡像・inlet 規定・outlet 外挿の物理を再実装しない)。
弱形式カーネルはこの $Q_b$ (= ゴースト値) と半割面幾何を使い、勾配・フラックスへの寄与だけを置き換える。
ゴースト CV は残すが、**勾配・対流フラックスの境界寄与をゴースト平均から弱形式へ差し替える**のが要点。

**node モードのパイプライン差し替え** (`cfg.discretization=="node"` で分岐、cell は不変):
1. **勾配** ([calcGradient_d.cu](../solver_density_cuda/cuda_forge/calcGradient_d.cu)):
   `calcGradient_cellgather_d` は CSR 中の**内部面のみ** (`ip < nNormalPlanes`) を集約 (node モードで境界面をスキップ)。
   その後 `boundaryGradientNode_d` が各境界半割面で $\phi_b \mathbf{S}_b$ (φ_b=ゴースト値) を raw 勾配へ atomicAdd。
   最後に `calcGradient_2_d` が体積正規化 (閉曲面が内部+境界で閉じる → §6.2)。
2. **対流フラックス** ([convectiveFlux_d.cu](../solver_density_cuda/cuda_forge/convectiveFlux_d.cu) wrapper):
   node モードは内部面 `[0,nNormalPlanes)` のみ処理し、境界は `boundaryFluxNode_d` が
   1 次 Riemann/SLAU フラックス $\hat F(Q_i,Q_b,n_b)S_b$ を残差へ atomicAdd。
3. **ロバスト化 (任意)**: 境界ノードに隣接する内部面の MUSCL 再構成を 1 次に落とすオプション
   (高マッハ壁近傍の過大再構成対策)。まず弱形式勾配で改善するか確認し、不足なら導入。

**検証**: bump hiM (Mach1.65) の node 2 次 MUSCL が**発散しなくなる**ことを `check_convergence.py` で確認
(M2 で step~400 発散 → M3 で安定収束が目標)。cell モードの全既存ケースが**従来とビット一致**(無変更) を回帰確認。

#### 7.0 [解決] node-centered 軸対称の near-axis corner

長く open だった「軸-境界 (inlet/outlet) 角 CV の P オーバーシュート/逆流」を **2 修正の併用で解決** (cell/SU2 一致):
1. **軸 (R≤eps) で軸対称ソース OFF** (`axisymmetricSource_d` の `axis_flag`、node+axisym のみ)。ソース
   `res_roUy += P·A_planar` は**唯一 r 重みされない項**で `source/volume = P/r → ∞` (r→0) が発散源。対流項は
   r_face/r_cell がキャンセルし有限。軸でソースを切るのが本質 (SU2 の y<EPS 相当、ユーザ提案)。
2. **境界半割面の重心に真の面積加重重心 (R≥0) を使う** (`dualBnodeCent`)。旧 `node+h·n_out` 便宜は pcy≈0/<0 で
   入口/出口 BC の r 重み実効面積を ~0 にし **BC が軸近傍 corner CV に届いていなかった**。真の重心 (R>0) で入口が
   コーナーに届き P=chamber に。(`axisCentroidShift=1` で CV セル中心に面積加重重心を使うのも維持: r_cell と r_face を
   一致させ対流項の r をキャンセルさせるため。)

**検証 (case/29 conical & bell, 軸対称 Euler implicit)**: PASS。corner P=3.99MPa(=chamber, overshoot 0 セル、旧 6.82e6)、
Ux=+51.7(流入, 旧 -798 逆流)、Uy≈0 — **SU2 (P=3.99,Ux=+54.8,Uy=0) と一致**。near-axis<3% mean|ΔM|=0.016/max0.049
(旧 max0.95) と中/外帯並み。平面 (bump loM) 無回帰。**explicit は startup の exit-wall 角で依然脆く発散** (cell explicit も
同様; 軸でなく supersonic-startup ロバスト性) → implicit が実用。

#### 7.1 診断結果 (重要): hiM 発散は「境界値」でなく「壁近傍 2 次再構成」の問題
弱形式実装の前に切り分け診断を実施 (`bndFirstOrder`: 境界隣接 CV の再構成勾配を 0 にし境界側を 1 次化、
[calcGradient_d.cu](../solver_density_cuda/cuda_forge/calcGradient_d.cu) の `zeroBndNodeGradient_d`、
フラグは `mesh.bnode_flag_d`)。結果:
- **`bndFirstOrder=1` で bump hiM node の NaN 発散 (step~400) が解消** (explicit は ~1e-3 のリミットサイクル、
  **implicit + bndFirstOrder で完全収束 PASS**)。cell モードは flag 既定 0 で無影響 (ビット不変)。
- **→ 真因は近壁 2 次 MUSCL 再構成のロバスト性であって、境界値 (φ_b) ではない**。slip 壁では弱形式の φ_b は
  現状の `0.5(node+ghost)` と一致するため、**弱形式境界だけでは hiM は直らなかったはず**。`bndFirstOrder` が
  実効的な最小修正。
- **弱形式境界 (§7) は引き続き粘性 (壁せん断・熱流の正しい評価) で必要**だが、hiM 安定化とは別件として整理する。
  よって `boundaryNode_d.cu` の新規作成は粘性着手時に回す (既存 BC カーネルの Q_b 流用方針は不変)。

#### 7.2 node モード: 壁の弱形式 pressure-only 対流 ＋ 壁優先コーナー所有 ＋ explicit (viscous)

理論は 本ドキュメントの「理論」節 §6.3/§6.4。専用計画:
[`discretization-node-boundary-ghostless.md`](../plans/active/discretization-node-boundary-ghostless.md)。
ねらいは case/29 viscous node の壁発散の解消。試行錯誤の結論を実装に反映:

**(A) 壁優先コーナー所有** ([gmshReader.hpp](../solver_density_cuda/mesh/gmshReader.hpp), commit 93ef041):
各境界ノードを優先度 (wall>inlet>outlet>slip/axis) で 1 bcond に所有。壁∩出口リップは壁所有となり矛盾する
2 重ゴーストが消える。`wall_flag_d` ([mesh.cpp](../solver_density_cuda/mesh/mesh.cpp)) も同時に構築。

**(B) 壁の弱形式 pressure-only 対流** (本節の要):
- 壁境界値 (bvar) を no-slip に: `wall_d` ([boundaryCond_d.cu](../solver_density_cuda/cuda_forge/boundaryCond_d.cu))
  が壁 bvar 速度=0, `roUb=0`, `roeb=内部エネルギーのみ` を明示。
- `normal_halo_planes_d` で**壁境界 plane を末尾に並べ** (`nWallHaloPlanes`)、node モードの主対流ループ
  ([convectiveFlux_d.cu](../solver_density_cuda/cuda_forge/convectiveFlux_d.cu)) は壁を除外
  (`convPlaneBound = nNormal_halo_Planes - nWallHaloPlanes`)。壁だけ `convectiveFlux_boundary_d` (bvar) で扱う:
  u=0 → `mdot=0` → **運動量=圧力のみ p·n·S** (ユーザ指摘の「壁残差に圧力寄与」を保持)。
  これで壁ノード (壁上に乗りミラーゴースト幾何が退化=fx≠0.5 で質量貫通) の発散源を断つ。
- **inlet/outlet/axis は従来どおり主ループ+ゴースト**で処理 → 既存の near-axis corner 修正 (§7.0) を保つ。
  最初に試した「全境界 弱形式 (convectiveFlux_boundary_d)」は inflow flux が SLAU と非等価で inlet-axis corner を
  6.1MPa に悪化させたため、**壁のみ**に限定した (対流のみの話。勾配境界は §6.2/§7.5 の owner-state 方式が全境界共通)。
- 勾配は cellgather が境界半割面 (`ip>=nNormalPlanes`) を skip し、`calcGradient_b_d` (§7.5) が owner ノードの
  状態値 $\phi_p$ を境界面値として加算して閉じる (§6.2)。壁は $u_p=0$ が Dirichlet ピン済みのため、これは事実上
  「壁面値 0 の GG 境界寄与」と同値になり、壁せん断は内部双対面が担う (theory §6.3)。

**(D) 入口∩壁コーナーの半割面所有 (`mesh.nodeInletCornerWall: 1`, 変換時, 2026-09-04)**:
マルチマーカ方式では入口∩壁の角ノード CV は入口側半割面 (入口速度で質量流入) と壁側半割面を両方持つが、
壁ノードは `nodeWallDirichlet` で $\mathbf u=0$ に固定されるため、入ってきた質量を内部双対面から流し出せず
角 CV に質量が溜まり $P$ が数十 bar に暴走する (case/47 Burrows–Kurkov: 1 次風上では延命、2 次で NaN。入口 BL
プロファイル + Crocco–Busemann ρ でも解消しない)。対策として **`inlet_*` 境界エッジの壁ノード側半割面を壁 bcond
に帰属**させる (`gmshReader::buildMedianDual`, `inletCornerWall`)。角 CV の境界は壁半割面 (pressure-only, 流入なし) だけに
なり $\mathbf u=0$ ピンと整合する。閉性は変わらない (半割面ベクトルは壁側に合算)。落とす流入質量は第 1 スペーシングの
半分 × 壁近傍速度で、BL プロファイル併用なら無視できる。既定 0 (従来のマルチマーカ・ビット不変)。出口∩壁は流入が
無いので対象外 (`inlet_` プレフィクスのみ)。軸∩入口 (§7.0) は壁でないので不変。

**(C) viscous は explicit (cfl≤0.1)**: 近壁の極小双対 CV (vol~1.8e-9) に対する粘性が **block-DPLUR の近似対角では
十分 implicit 化されず** (cfl 2→0.1 で発散 step13→610=構造的)、陰解法は viscous node で発散する。一方 **explicit
(timeIntegration=3 RK3) は `setDT` の粘性スペクトル半径 `2(visc+vt)/(ro·dx_min)` で局所 dt が近壁で縮むため安定**。
cfl=0.5 は発散、**cfl=0.1 で安定収束** (検証: case/29 conical, 30k step 完走・NaN なし・残差 ~3 桁低下、場は物理的
P≤Pt/ro>0/T>0; 高 Re 層流で BL は微小厚=未解像のため場は概ね Euler 相当=正しい挙動)。陰解法での viscous node は
近壁粘性 Jacobian (スペクトル半径) 強化が要 (gpu-implicit-plan の後続)。

**棄却した方式**: 残差射影 Dirichlet (`zeroWallMomentumResidual_d`, `nodeWallDirichlet`) は壁圧力寄与も落とすため
不正で壁全域発散 (ユーザ指摘どおり)。既定 OFF で残置。壁 plane 不出力 (ghost 完全撤廃) は勾配閉性を壊し悪化したため
不採用 (壁 ghost は勾配/粘性 mirror 用に残す)。

##### 7.2.1 陰解法での壁 Dirichlet 整合: Jacobian 行 decouple (SU2 `DeleteValsRowi` 相当)

専用計画: [`discretization-node-wall-implicit-dirichlet.md`](../plans/accepted/discretization-node-wall-implicit-dirichlet.md)。

**症状**: `nodeWallDirichlet=1` + implicit (block-DPLUR) で壁ノード速度が drift する。残差射影 (`zeroWallMomentumResidual_d`)
は RHS を 0 にするが、block-DPLUR の対角 5×5 は壁運動量を連続・エネルギーと連成したまま (`wall_flag` が Jacobian に
出てこない) なので、solve が $\Delta(\rho\mathbf{u})\neq0$ を返す。状態射影 `enforceWallNoSlip` が毎ステージ再射影して
KE を剥ぐが、これは持続的なエネルギーシンクになり、近壁圧力場を蹴る (実測 case/36 node 層流 implicit: 壁
$|\rho u_x|$ mean 5.6→34 に増大、KE 除去 mean ~4000/stage)。**explicit (RK3) は残差 0 → 更新で $\Delta(\rho u)=0$ なので
drift しない**ので、本症状は implicit 特有。

**処方 (SU2 と同じ in-Jacobian)**: `implicit_defect_correction_block_d` ([timeIntegration_d.cu](../solver_density_cuda/cuda_forge/timeIntegration_d.cu))
の既存 `axis_flag` 行 decouple 機構を**壁の運動量3行 (index 1=roUx, 2=roUy, 3=roUz) へ拡張**する:

```cpp
if (wall_flag != nullptr && wall_flag[ic] == 1) {
    for (int row = 1; row <= 3; ++row) {
        for (int jj = 0; jj < 5; ++jj) diag_block[row][jj] = 0;
        diag_block[row][row] = 1; rhs[row] = 0;   // Δ(ρu_i)=0
    }
}
```

連続 (行0)・エネルギー (行4) は不変なので $\rho,\rho e$ は保存式で発展し、圧力復元 thermo (CPG/TP) も無変更。
起動側で `(discretization=="node" && nodeWallDirichlet==1) ? msh.wall_flag_d : nullptr` を渡す (現状 `nullptr`)。

> **軸 (§7.1) との違い**: 軸 `roUy` 単行 decouple は r=0 の極小 CV が多方程式的に ill-posed (Mach~1000 発散) で
> 単独では効かなかったが、**壁は r=0 特異点を持たず**、運動量3行 decouple が drift 機構を直接除去する。decouple 後は
> `enforceWallNoSlip`/`zeroWallMomentumResidual` は実質 no-op (剥ぐ KE→0) になり、保存性も回復する見込み。検証は
> case/36 node (壁 $|\rho u_x|$ 非増大) + case/29 回帰 (デグレ無し) + cell ビット不変。

**Phase 2 (後続) — 残り (inlet/outlet/slip/axis) のゴースト撤廃**: `convectiveFlux_boundary_d` を SLAU 等価にした上で
半割面弱形式へ、スカラ輸送 (rans/species) も境界カーネル化、その後 node モードで ghost 生成停止。

##### 7.2.2 node モード: Green-Gauss 境界寄与 (`calcGradient_b_d`) — owner-state 方式 (2026-08-11 全面改訂)

理論は §6.2 (弱形式境界 — 勾配の境界寄与は owner ノードの状態値を使う)。実装は
[calcGradient_d.cu](../solver_density_cuda/cuda_forge/calcGradient_d.cu) の `calcGradient_b_d`。

- **対象**: 全 non-periodic bcond (`wall`/`wall_isothermal`/`inlet_*`/`outlet_statPress`/`slip`/`axis` 等)、
  全 primitive (`ro, Ux, Uy, Uz, P, T`)。periodic は §4.5 の DOF 同一視・gradient gather に委ね、境界寄与を加えない
  (`calcGradient_d_wrapper` のループが `bc.bcondKind=="periodic"` を skip)。
- **境界面値**: $\phi_b := \phi[\text{owner node}]$。bvar (`bc.bvar_d[...]`) は**一切参照しない**
  — 旧実装は `ro/Ux/Uy/Uz` を常に bvar から、`P/T` は bcond ごとの特例 (`interiorPT` フラグ, §7.7.1 の
  node outlet 修正で導入) から取っていたが、これを全変数・全 bcond で owner state 参照に統一した
  (`interiorPT` フラグと `P_c`/`T_c` 引数付き特例は撤去)。
- **理由**: bvar には状態と独立に動く量が混ざる (規定背圧、SST 断熱壁の回復温度モデル $T_{aw}$ 等) ため、
  勾配へ入れると内部場が意図せず駆動される (§6.2 参照)。no-slip/等温壁/axis のように状態が既に
  BC 値へピンされている境界では owner state = bvar 相当なので**数値は変わらない** — 実際に効果が出るのは
  状態と bvar が乖離し得る境界 (outlet, SST 断熱壁の $T_{aw}$) のみ。
- **旧実装からの移行**: node outlet_statPress の P/T だけを owner state にする特例 (commit ddbe9ce5,
  §7.7.1 [`boundary-node-nozzle-wall-outlet-stability.md`](../plans/active/boundary-node-nozzle-wall-outlet-stability.md) §2.11)
  は、この一般化によって「outlet だけの特例」ではなく「全境界の既定」に格上げされ、コード上の特例分岐は不要になった。

#### 7.3 node モード: 最小二乗 (LSQ) 勾配 (`gradLSQ`, 既定 OFF)

専用計画: [`discretization-lsq-gradient.md`](../plans/active/discretization-lsq-gradient.md)。

**動機**: node-centered (median-dual) の近壁では Green-Gauss 面勾配 (面値 $\phi_f$ を両側中心の線形補間で作る、
[gradient.md](gradient.md#理論)) が **checkerboard (奇偶デカップリング) モード**を持ち、薄い粘性境界層で
速度勾配 $\to$ せん断応力を振動させる。GG は隣接「面」を介した広い stencil の平均で、median-dual の歪んだ近壁双対面では
中心点の値変動を平滑化しきれない。**重み付き最小二乗 (LSQ) 勾配**は近傍 CV「中心」の差分を直接フィットするため、
この checkerboard を持たず近壁で素直な勾配を返す。

**手法**: セル $i$ の勾配 $\mathbf{g}$ を、近傍 $j$ について逆距離二乗重み $w_{ij}=1/|\mathbf{d}_{ij}|^2$
($\mathbf{d}_{ij}=\mathbf{x}_j-\mathbf{x}_i$) を付けた正規方程式

$$
\underbrace{\Big(\textstyle\sum_j w_{ij}\,\mathbf{d}_{ij}\mathbf{d}_{ij}^{\mathsf T}\Big)}_{M\ (\text{対称}3\times3)} \mathbf{g}
= \underbrace{\textstyle\sum_j w_{ij}\,\mathbf{d}_{ij}\,(\phi_j-\phi_i)}_{\mathbf b}
$$

で解く。**境界に疑似点を追加しない** (2026-08-11 改訂。§7.2.2 の GG owner-state 方式と同じ原則 — 理由は §6.2)。
境界ノードの近傍は実在する内部隣接ノードのみで構成し、φ_b=φ_i (ジャンプ 0) の疑似点であっても正規行列 $M$ に
入れると「その方向の勾配をゼロへ暗黙に正則化する」効果を持つため、**疑似点を作らずそもそも $M$ に加えない**
(ゼロジャンプの $\mathbf b$ だけ足しても $M$ が入っていれば正則化は残るので、$M$ からの除外が必須)。
境界隣接ノードで近傍配置が退化する場合は §7.3.1 のスペクトル打ち切り擬似逆で扱う。

**実装** ([calcGradient_d.cu](../solver_density_cuda/cuda_forge/calcGradient_d.cu)): 2 カーネルで構成
(旧 `lsqGrad_accumBoundary_d` は境界疑似点ごと撤去)。
1. `lsqGrad_accumInternal_d` — per-cell gather で**内部面のみ** (`ip < nNormalPlanes`) を走査し $M$ と $\mathbf b$ を蓄積。
   $\mathbf b$ は既存の勾配配列を流用して書く (専用バッファ不要)。$M$ は static な scratch (`Mxx..Mzz`、`nCells` 変化時のみ再確保)。
   境界半割面は最初から走査対象外 (`if (ip>=nNormalPlanes) continue`) なので追加の除外処理は不要。
2. `lsqGrad_solve_d` — 対称 $3\times3$ を変数ごとに解いて勾配を**その場で上書き** (`lsq_solve_sym3`; $m_{zz}\le\text{tol}$ を
   2D と判定し $xy$ の $2\times2$ を解き $g_z=0$)。最後に `divU` を更新。$M$ の組み立て・solve は倍精度。

`gradLSQ=1` かつ `discretization=="node"` のときのみ有効で、LSQ 経路は wrapper 内で**早期 `return`** し GG の正規化
(`calcGradient_2_d`) と弱形式境界加算 (§7.2 の `calcGradient_b_d` ループ) をスキップする (LSQ は勾配を直接解くため正規化不要)。
**cell モード、および `gradLSQ=0` (既定) の node モードは従来の GG 経路のまま**で、数値はビット不変。

**設定**: `mesh.gradLSQ` (0:GG, 1:LSQ 毎ステップ solve, 2:LSQ 係数事前計算 — 推奨)。
[solverConfig.hpp](../solver_density_cuda/input/solverConfig.hpp)。

**検証状況**: `gradLSQ=1` は case/29 (軸対称ノズル) の GPU 検証で**近壁発散の負結果** (計画 §9)。実メッシュ診断
([tools/check_lsq_gradient.py](../solver_density_cuda/tools/check_lsq_gradient.py)) により、真因は (a) 正規方程式を
float32 格納で解く条件数 2 乗の増幅と、(b) 一部メッシュに実在する**真の退化ノード** (近傍方向が共線/共面で
勾配を決める情報が無い; case/29 は近壁 2.6%) の複合と判明。これを解決するのが mode 2。

##### 7.3.1 `gradLSQ=2`: 係数事前計算 + スペクトル打ち切りフォールバック

メッシュが静的なら正規方程式の解は**幾何のみの線形演算子**に畳める:

$$
\mathbf g_i = M_i^{-1}\mathbf b_i = \sum_j \underbrace{\bigl(M_i^{+}\,w_{ij}\,\mathbf d_{ij}\bigr)}_{\mathbf c_{ij}\ (幾何のみ)}\,(\phi_j-\phi_i)
$$

- **setup (初回 1 回, device, 倍精度)**: ノードごとに $M$ を**内部双対面の隣接ノードのみ**で組む
  (2026-08-11 改訂: 非 periodic 境界半割面重心を LSQ 点として $M$ に加える旧処理は撤去した — §7.3 の理由により、
  ジャンプが常に 0 の疑似点でも $M$ に入れば暗黙の正則化になるため。periodic は従来どおり DOF 合併機構に委ねスキップ)。
  対称 3×3 の**解析固有分解**から
  $\lambda_k < \max(\texttt{gradLSQDegenThresh}\cdot\lambda_{\max},\ \varepsilon)$ のモードを落とした
  **スペクトル打ち切り擬似逆行列** $M^{+}=\sum_{\text{kept}}\lambda_k^{-1}\mathbf v_k\mathbf v_k^{\mathsf T}$ を取り、
  内部 incidence の係数 $\mathbf c_{ij}$ を float32 テーブル (`cell_planes` CSR 対応) に焼き込む
  (2026-08-11 改訂: 境界疑似点の係数テーブル `cBnd` は撤去 — 境界は $M$ に寄与しないため焼き込む係数が無い)。
  旧 `lsq_solve_sym3` の 2D 分岐 ($m_{zz}\le\text{tol}$) はこの打ち切りが一般化して包含する
  (1 セル厚 2D は $\lambda_z\approx0$ が自動で落ちる)。退化ノード数は起動ログに出す。
- **runtime (毎ステップ, 全 float32)**: $\mathbf g_i=\sum_j\mathbf c_{ij}(\phi_j-\phi_i)$ の内部隣接 gather のみ
  (per-node、境界寄与なし)。$M$ の組立・solve は消滅し、係数は 6 変数で共有 → mode 1 より速い。
  悪条件 solve を setup 側 double に追い出したので float32 場の丸め以外の増幅は無い。
- **フォールバックの設計判断**: 退化方向の勾配成分は**ゼロ化** (その方向 1 次精度化)。GG 係数への差し替えでなく
  打ち切りを選んだのは、追加幾何 (面積ベクトル・体積・軸対称 planar 変種) への依存が無く全メッシュ形態で一様に成立し、
  実績ある 2D $m_{zz}$ 分岐の自然な一般化だから。壁法線方向の粘性は面法線コンパクト差分が主担なので、
  少数 (case/29 で 2.6%) の退化ノードの再構成勾配 1 次化は安全側。閾値は `mesh.gradLSQDegenThresh` (既定 1e-2)。

#### 7.4 node モード: 内部双対面の面補間係数を中点 (fx=0.5) に固定

node-centered (median-dual) では、双対面を挟む 2 つの CV 中心は**ノード**であり、面値は
ノード–ノード中点での平均 $\phi_f=\tfrac12(\phi_A+\phi_B)$ を取るのが標準的なエッジベース補間 (2 次精度)。
しかし forge の `dualFaceCent` は各セルの (エッジ中点 $M$ → セル重心 $G$) 区間の**面積加重重心**で、
$M$ そのものではない ([gmshReader.hpp](../solver_density_cuda/mesh/gmshReader.hpp))。このため面補間係数
`fx = dc1p/(dc0p+dc1p)` ([calcStructualVariables_d.cu](../solver_density_cuda/cuda_forge/calcStructualVariables_d.cu)、
面重心を法線方向に射影した距離比) は対称メッシュでは 0.5 だが**歪み面では 0.5 からずれ**、ノード間の中点で
値を取らない。

**修正 (2026-06-14, フラグは 2026-08-16 に撤去 — 値=ノード座標では幾何 fx が中点相当)**: フラグ `nodeMidpointFx` (既定 0=OFF) を追加。`discretization=="node"` かつ
`nodeMidpointFx==1` のとき**内部双対面** (`ip < nNormalPlanes`) を `fx[ip]=0.5` に固定し
$\phi_f=\tfrac12(\phi_A+\phi_B)$ の中点補間にする ([calcStructualVariables_d.cu](../solver_density_cuda/cuda_forge/calcStructualVariables_d.cu))。
cell モード・境界半割面 (`ip>=nNormalPlanes`) は対象外。`fx` は対流・粘性の全面補間で使われる
(gradient の over-relaxed 法線項が使う `dcc` は CV 中心間距離のまま)。

**検証 (case/29 conical, node laminar viscous 40k)**: `fx=0.5` は近壁 `dUxdy` の checkerboard roughness を
低減 (99pct 12.84→8.23, −36%)。SU2 (axisym laminar 同条件) との**壁圧比較で fx ON/OFF は区別不能** (<0.5% 差、
両者とも SU2 に平均 3.2% 一致、最大は未収束の超音速出口)。局所 Ux は出口リップ近傍で大きく変わる
(1376→430) が圧力場に伝播しない近傍速度の細部。**SU2/forge とも本ケースは未収束のため near-wall 速度細部の
厳密な是非は未確定**だが、fx=0.5 は SU2 一致を悪化させず checkerboard を下げる。既定 OFF・opt-in で残置。
計画は [`discretization-median-dual.md`](../plans/active/discretization-median-dual.md)。

### 5. 設定

[solverConfig.hpp](../solver_density_cuda/input/solverConfig.hpp) に `discretization` (`cell`/`node`, 既定 `cell`)。
設定の意味は [procedures/solver-settings.md](../procedures/solver-settings.md) を参照。

### 6. 検証

[.github/plans/discretization-median-dual.md](../plans/active/discretization-median-dual.md) §6 を参照。
要点: cell モードで既存ケースが従来と完全一致すること (デグレ無し) を回帰確認した上で、
node モードを同一メッシュで実行し sod の shock 位置、tet ケースの DOF/コスト/精度を比較する。
