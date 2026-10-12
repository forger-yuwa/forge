# 風洞ノズルの実寸の要求カードと成果物一式 (STEP・gmsh・Fluent・Excel)、skill の整備

## メタ

- **area**: `tooling / design export`
- **status**: `in_progress`  <!-- 2026-10-12 起票。2026-10-13 plan 段の諮問・レビュー済み、実装中 -->
- **related_docs**:
  - [`methods/design/overview.md`](../../methods/design/overview.md) 「実寸の要求カードと成果物一式」節 (本計画と同時に起草)
  - [`methods/boundary.md`](../../methods/boundary.md) 「界面熱量の 4 つの定義」(壁への熱流束の正本 `iface_q_eff`)
  - [`procedures/nozzle-design-workflow.md`](../../procedures/nozzle-design-workflow.md) (axismach の作業手順の正本)
  - [`procedures/nozzle-design-outputs.md`](../../procedures/nozzle-design-outputs.md) (標準報告 pptx)
- **related_plans**:
  - [`../accepted/tooling-nozzle-wall-single-bspline.md`](../accepted/tooling-nozzle-wall-single-bspline.md) (物理壁 1 本の B-spline と STEP 書き出し `export/wall_step.py`、W5 の転送誤差)
  - [`tooling-design-problem-campaign-recipe.md`](tooling-design-problem-campaign-recipe.md) (problem/campaign/recipe の入力構造の統一。§5.1 #15 寸法固定の campaign。本計画の要求カードは**現行の problem YAML を作る薄い入力層**で、#6 の新スキーマができたらそちらへ写す)
  - [`tooling-nozzle-isothermal-wall-chain.md`](tooling-nozzle-isothermal-wall-chain.md) (等温壁 `spec.wall_thermal`、CONTUR の壁温の危険 §5.1 #10〜#12)
- **created**: `2026-10-12`
- **owner**: `sano`

## 1. 目的

ユーザは今後、いろいろなマッハ数・形状の風洞ノズルを **Claude に指示して** 作らせる (2026-10-12 ユーザ)。入力は実寸
(全圧・全温・ガス・出口マッハ数・物理スロート径・入口半径・入口直管長・縮流部長・スロート→出口長・壁温) で、
欲しいのは (1) 設計した形状の STEP、(2) Fluent で同じ解析メッシュを使うための gmsh 入力と Fluent 用メッシュ、
(3) STEP から gmsh で**同じ構造のメッシュ**を作り直す手段、(4) Excel で扱える形状表と、1 次元の壁熱伝導の評価に使える壁面量の表。

完了時の状態:

1. 実寸の**要求カード**を書けば、問題 YAML (`wind_tunnel_axisym_axismach`) が生成され、$r_t$ を解いたあとも実寸の長さが保たれる。
2. 最終 NS の run から **1 コマンドで成果物一式** (STEP 2D/3D・STL・gmsh msh・STEP→gmsh 再構築スクリプト・Fluent `.msh`・CGNS・境界条件表・Excel 表・README) が出て、各成果物に独立の読み直し検査の結果が付く。
3. skill `nozzle-design` (と `design-intake`) が、ユーザの依頼文からカード作成 → チェーン → 成果物までを Claude が迷わず進められる記述になっている。

## 2. スコープ

- **やる**
  - 要求カード `forge.wt_card/v1` と生成器 (カード → problem YAML、$r_t$ の解き直しに伴う長さの割り直し、壁解像の第 1 層と `ni` の初期見積もり、空気の液化の事前警告)。
  - 成果物の書き出し `forge_design.export.deliverables` (下表の全部) と各成果物の読み直し検査。
  - STEP→gmsh の再構築スクリプト (gmsh Python API) と構造ファイル `mesh_structure.json`。
  - Fluent 形式 (ASCII `.msh`、2D 軸対称) の書き出し器。gmsh 経由の CGNS も併せて出す。
  - Excel 向けの表 (CSV、ASCII の列名) と、`.venv-pptx` に `xlsxwriter` があれば同じ内容の `.xlsx`。
  - skill `nozzle-design`・`design-intake` と `procedures/nozzle-design-workflow.md`・`nozzle-design-outputs.md` の同期。
  - デモ: M5 の空気・壁温 350 K のノズルを**CFD 前まで** (カード → MOC → CONTUR 初期の物理壁 → メッシュ → 成果物) 通す。
- **やらない**
  - 設計チェーン本体 (MOC・δ\* ループ・C2 較正) の変更。出口の後ろの直管 (ユーザ回答 2026-10-12: 出口長さは「スロート→出口」)。
  - 3D の解析メッシュ (Fluent の 3D 計算用)。2D 軸対称だけ。
  - 格子線が縦でない格子 (壁法線の近壁層、別ブランチ `feature/nozzle-wall-fit-and-pipeline` の `wall_normal_layer`) の再構築。v1 は検出して拒否する (§5.1 で後続)。
  - Fluent での読み込み試験 (手元に Fluent が無い)。ユーザに 1 回確認してもらう (§6 V4c)。
  - M5 デモの CFD (Euler・NS) 本体。GPU を使う段は、ユーザが実際の条件で指示したときに回す。
  - 壁と固体の連成計算。1 次元の壁伝導は既存の `cht_loop` + `solid_shell` があることを skill に書くだけ。

## 3. 関連 docs と前提 (観測事実, 2026-10-12)

- 設計チェーンは `main` で skill `nozzle-design` + `procedures/nozzle-design-workflow.md` にまとまっている。入力キーは $r_t$ 単位
  (`geometry.r_inlet`・`L_U`・`L_pipe`、`dv.L_total` + `Lc_mode: from_length`)。手順書 §0 に「r_t を変えても r_inlet・L_total は
  r_t 単位のまま (実寸がずれる)」と落とし穴として書かれている。
- 物理スロート径から $r_t$ を決める `deltastar_loop --solve-rt-throat R` (CFD 前: 許容差スロート 1e-7 m) と、CFD なしで NS の準備だけする
  `deltastar_loop --init-integral --prepare-only` がある。
- 壁温は `spec.wall_thermal: {mode: isothermal, Tw}` (一様のみ)。CONTUR の熱閉包は温度形で、冷却壁で δ_r が最大 2.8 % 変わりうる (isothermal plan §5.1 #10)。
- STEP: `design/forge_design/export/wall_step.py` が物理壁 (表現 `single_bspline`) の 2D 曲線を FreeCAD で書き、読み直して転送誤差
  (位置 1e-6 mm 等) を検査する。3D の回転面・STL・座標表は case/45 の `cad_deliverables.py` (case 専用のパス) にしかない。
- メッシュ: `meshing/mesh2d.py` (main) は格子線 $x=x_i$ の縦線、$r_{ij}=r_w(x_i)\,s_{ij}$、msh4.1 を 10 有効桁で書く。物理タグは
  inlet 1 / outlet 2 / wall 3 / axis 4 / fluid 5。
- gmsh 4.15.2 (`.venv-mesh`) は msh・CGNS・UNV・BDF・STEP を書ける (2026-10-12 確認)。`.venv-opt` に openpyxl は無く、`.venv-pptx` に xlsxwriter 3.2.9 がある。
- Fluent の面の向き (2026-10-13 測定):
  - **2D の実物** (TGrid 2D、OpenFOAM-10 `tutorials/incompressible/icoFoam/elbow/elbow.msh`、537 節点・918 セル・1,454 面): **全 1,454 面で
    $c_0$ が $n_0\to n_1$ の進行方向の左側** (内部 1,300・境界 154)。ファイル冒頭の注記は `(node-0 node-1 right-cell left-cell ...)`。
    取得元 `https://raw.githubusercontent.com/OpenFOAM/OpenFOAM-10/master/tutorials/incompressible/icoFoam/elbow/elbow.msh`、sha256 `a7a3c88110d19f5534b5a225eb99102f1416d5113331c61ee47de1a416c299e3`。
    OpenFOAM の配布物 (GPL) なのでリポジトリには入れない (公開リポジトリ)。試験は手元の写しがあれば走らせ、無ければ飛ばしたことを出力する。
  - 3D の実物 `case/32.fan_duct_fluent/mesh/fan.msh.h5` (内部面 24,542・境界面 1,484) と `case/30.mixing_elbow_fluent/mesh/mixing_elbow.msh.h5`
    (86,782・4,799) では、節点順の右手系の法線 (Newell) が全面で $c_0$ 側を向いた。
  - 起票時は 3D の規則を $z$ に押し出して「2D は右側が $c_0$」と推定したが、2D の実物で**逆**だった。推定で書き出し器の規則を決めない。
  - ParaView の `vtkFLUENTReader` は $c_0$ のセルで面の節点順をそのまま使ってセルを組む。右側 $c_0$ で書いた試作は 15 セル全部が時計回り
    (面積の総和は正しく 0.65) で読めた — 読み込み器は向きの誤りを止めない (向きの検証には使えない)。
- 壁面量: `res_wall_<physID>_*.h5` に `Ps`・`Ts`・`qwall` (流体向き正、診断)・`twall_{x,y,z}`・`utau`・`ypls`。`output.interfaceDiag: 1` で
  `iface_q_eff`・`iface_Qf_eff` (正本、固体向き正) ほかが加わる。`ypls` は壁解像の根拠にしない (AGENTS.md「壁解像確認」)。
- **軸対称 (`axisymMethod: 0`) の `iface_q_eff` は $q$ でなく $r\,q$** (codex plan レビュー M1 を当方でコードで確認, 2026-10-13):
  device の面ベクトルは $r$ を掛けている (`solver_density_cuda/variables.cpp:599-604` の `ss[ip] *= r_face`) が、`iface_q_eff` は
  $r$ を掛けていない host の `msh.planes[ip].surfArea` で割る (`conjugateWall.cpp:229`)。積分済みの `iface_Qf_eff` は 1 rad あたりの
  節点熱量 [W/rad] で、全周は $2\pi$ 倍。`qwall`・`twall` は $r$ 重みの面積で割るので物理量 [W/m²]・[Pa] のまま。
  `iface_q_eff` 自体の定義は本計画では変えない。**`feature/cht-conjugate-benchmarks` の `conjugateWall.cpp:223-250` は既に $r$ 重みの面積で割って
  `iface_q_eff` = $q$ になっている** (main は $r q$ のまま; 2026-10-13 implementer の確認、case/63 run_0013 で本計画の $q$ と同ブランチの `iface_q_eff` がビット一致)。
  `iface_Qf_eff` から作る方法はどちらの版でも正しい。
- `iface_q_eff` 系は等温壁 (`wall_isothermal`)・定常 (`unsteady: 0`)・`interfaceDiag: 1` のときだけ値が入り、それ以外も**データセットは
  出て中身が NaN** (`conjugateWall.cpp:192-215`)。「データセットがある」は「使える」ではない。

## 4. 設計方針

### 4.1 要求カード `forge.wt_card/v1`

実寸 [m, Pa, K] で書く YAML。未知のキーは拒否、必須キーの欠けも拒否。

```yaml
format: forge.wt_card/v1
name: m5_air_demo
gas: {preset: air}                 # air = N2 0.7552 / O2 0.2314 / AR 0.0129 / CO2 0.0005 (質量分率) の semiperfect。species: {...} で直接も可
conditions: {Pt: 1.0e6, Tt: 400.0}
M_exit: 5.0                         # = spec.M_design
throat: {diameter_physical_m: 0.050}   # 物理スロート径 (δ 込み)。diameter_design_m (設計 = 非粘性) でも可、どちらか 1 つ
inlet: {radius_m: 0.150, pipe_length_m: 0.20, contraction_length_m: 0.40}
divergent_length_m: 1.20            # 設計スロート (x=0) → 物理出口 x_F
wall: {thermal: isothermal, Tw: 350.0}   # adiabatic も可
p_ambient: design                  # 既定: 設計出口静圧 (等エントロピー、semiperfect)
tolerances: {length_m: 1.0e-3, diameter_m: 2.0e-5}   # 出荷ゲート (§4.2a) の公差。既定値は仮 (ユーザの実用公差が決まったら置き換える)
recipe: axismach_explore_v1        # 形状の詰め (joint・CFD ピン・C2) は skill の手順で後から足す
overrides: {}                      # 生成した YAML のキーを上書き (geometry.R、M_knot など)。上書きしたキーは YAML のコメントに残す
```

生成器 `python -m forge_design.intake.wt_card CARD --out problem.yaml [--r-throat R]`:

1. **長さの写し**: $r_t$ (初回は設計スロート径/2、または前回値) で割る。`geometry.r_inlet = R_in/r_t`、`L_pipe`・`L_U` 同様、
   `Lc_mode: from_length`・`dv.L_total = L_div/r_t`。**`geometry.L_total_tol` はまずカードの長さ公差の半分 / $r_t$** にする (既定の 0.05 $r_t$ は
   $r_t$ = 25 mm で 1.25 mm。codex 諮問 2026-10-13)。届かなければ `n_axis_inv` を 1 回 2 倍にし、それでも届かなければ**内部許容差を既定 0.05 $r_t$ に
   緩めて作り直し、緩めたことを YAML・`_card` 節・card_log に記録する** (`L_total_tol_relaxed`; codex 諮問 2026-10-13 の 2 回目
   `notes/reviews/2026-10-13-wt-length-accuracy-diagnose.md`)。**カードの公差は変えず、合否は §4.2a の出荷ゲートだけが決める** (未達は FAIL のまま返す)。
   止まった理由は「from_length の固定点反復が改善しなくなった」(`runner_axismach.py:606-616`) と書き、離散化のノイズ床とは書かない (区別できない)。
   各反復の $(L_c, x_F, x_e, x_F - x_e)$ を card_log に残す (§5.1 #15 の切り分け材料)。
2. **$r_t$ の固定点** (物理スロート径を与えたとき): `solve_rt_throat` (CFD 前の δ 経路) で $r_t$ を解き、1. をやり直す。
   $|\Delta r_t|/r_t \le 10^{-6}$ で止め、10 回で止まらなければ例外。CFD 後に解き直したときも同じ割り直しを生成器に `--r-throat` で渡して行う。
   設計値を固定する場合 (`diameter_design_m`) は 1 回で済む。
3. **縮流部の検査**: $\mu = L_U^2/(R\,\Delta r)$ (手順書の単調収縮の条件 ≤ 20) と、$R_{in} > r_t$ を検査し、外れたら例外。
4. **壁解像の初期見積もり** (`mesh` の `wall_first_frac`・`_throat`・`ni`): 等エントロピーの 1 次元分布と Eckert の参照温度・平板の $C_f$ 相関で
   $y_1^+ \le 0.8$ となる第 1 層厚を断面ごとに見積もり、スロート近傍と下流の最小値を取る。`ni` は壁際の AR ≤ `ar_max` (5000、境界層層の例外) から下限を出す。
   **見積もりは初期値であって合格の根拠ではない** — NS 後に `check_wall_resolution.py` の VERDICT で判定し、不合格なら細分する (手順書 §2 ⑥)。
5. **空気の液化の事前警告**: 出口の静温・静圧 (等エントロピー) を N₂ の飽和曲線と比べ、過冷却なら警告を YAML のコメントと標準出力に出す
   (止めない。判断はユーザ)。
6. **recipe**: `axismach_explore_v1` は探索段階の既定 (Hall 初期線・補間壁・`physical_wall_repr: single_bspline` は付けない)。
   STEP を出すには物理壁の表現 `single_bspline` が要るので、**成果物を出す最終の NS 用 YAML には `axismach_final_v1`** (joint 壁・`pw_upstream: poly`・
   `single_bspline`・MOC analytic/converge) を使う。どちらも雛形は `design/forge_design/intake/recipes/*.yaml` に置き、雛形の各キーに出典 (case/45 の生産 YAML) を書く。
   形状パラメータ (R・M_knot) の既定は雛形の値で、成立しない (law の窓・壁 QA で例外) ときは skill の手順で振る (自動探索は campaign plan #15)。

### 4.2 成果物の書き出し `python -m forge_design.export.deliverables RUN [--out DIR] [--res N] [--skip ...]`

入力は最終 NS の run (`nozzle.msh`・`wall_repr.json` (表現 `single_bspline`)・`prepare_info.json`・`bcondConfig.yaml`・`res_wall_<wallID>_*.h5`)。
足りないものがあれば、その成果物だけを飛ばして理由を `deliverables/MANIFEST.json` に書く (黙って欠けさせない)。出力:

```
deliverables/
  README.md                 単位・原点・座標系・各ファイルの出典 (run・res・sha256)・検査結果の要約
  MANIFEST.json             全ファイルの sha256・検査の数値・飛ばしたものと理由
  cad/wall_physical.step    2D 曲線 (wall_step.py)          cad/wall_physical_step.json (転送誤差)
  cad/wall_surface.step     3D 回転内面 (FreeCAD)           cad/wall_surface.stl
  mesh/nozzle.msh           run の解析メッシュ (msh4.1、そのまま)   mesh/nozzle_v22.msh (msh2.2 ASCII、gmsh で変換)
  mesh/mesh_structure.json  構造ファイル (§4.3)              mesh/rebuild_from_step.py (§4.3)
  fluent/nozzle_2d_axi.msh  Fluent 形式 (§4.4)               fluent/nozzle_2d_axi.cgns (gmsh 経由)
  fluent/boundary_conditions.md  Fluent の設定表 (§4.4)
  excel/wall_shape_mm.csv   形状表 (§4.5)                    excel/wall_cfd_nodes.csv (§4.5)
  excel/conditions.csv      条件と寸法                        excel/nozzle_tables.xlsx (xlsxwriter があれば)
```

### 4.2a 寸法の出荷ゲート (codex 諮問 2026-10-13 採用)

YAML の値 × $r_t$ の往復は単位換算の検査であって、できた形状の検査ではない。成果物の書き出しの最後に、**配布する STEP とメッシュから寸法を測り直し**、
カードの値と比べる (`deliverables/DIMENSIONS.json`)。1 つでも公差外なら MANIFEST の判定を `FAIL` にし、README の冒頭に書く (書き出しは止めない — 何が外れたかを見せるため)。

| 寸法 | 測り方 (STEP の曲線 $r(x)$、$x$ は設計スロートが原点) | 公差 |
| --- | --- | --- |
| 入口半径 | 曲線の始点の $r$ | `diameter_m`/2 |
| 直管長 `L_pipe` | **STEP の B-spline の制御点とノットから**定半径区間の終端を同定する: 始点から、影響する制御点の $r$ がすべて $r_{in}$ (差 ≤ 1e-9 mm) であるノット区間の最後のノット $u^*$ ($x(u^*)$ − 始点)。半径の閾値で探さない (接合点は $r'=r''=0$ なので閾値は下流にずれる; codex plan レビュー M2 の反例で +15〜18 mm) | `length_m` |
| 縮流部長 `L_U` | $x(u^*)$ から設計スロート $x=0$ まで (物理スロート $r'=0$ の位置も併記) | `length_m` |
| スロート径 | `diameter_physical_m` のとき: $\min r$ × 2 (物理)。`diameter_design_m` のとき: 設計壁 (`wall_repr.json` の設計壁) の最小半径 × 2 と比べ、物理の最小径は参考値として併記 | `diameter_m` |
| スロート→出口長 | **物理出口 = 最終 STEP の下流端** (設計スロートが原点)。MOC の $x_F$ で代用しない。メッシュの `outlet` の $x$ とは端点の転送・丸めの範囲 (1e-3 mm) で一致することを別に検査し、外れたら FAIL (途中の弦誤差で許容しない) | `length_m` |
| 出口マッハ数 | 測らない (CFD の判定。報告の出口コア M を併記するだけ) | — |

**誤差の量を分ける** (codex plan レビュー m6): (1) CAD の転送誤差 (`wall_step.py` の W5、STEP と保存した B-spline の差)、
(2) 再構築の節点誤差 (§4.3、元メッシュと作り直したメッシュの節点の差)、(3) 弦の形状誤差 (メッシュの壁の多角形と曲線の差の最大と、
最小半径の取り逃し = 曲線の $\min r$ と壁節点の $\min r$ の差; `wall_step.chord_deviation` を使う)。メッシュから測った寸法は
(3) の範囲で STEP と一致すればよい。直管の終端はメッシュでは区間 (接合点を挟む 2 節点の $x$) として返し、測定の不確かさを併記する。
NS 後に $r_t$・k_f を解き直して壁を作り直した最終 run にも同じゲートをかける (カード生成時の合格は最終形の合格ではない)。

### 4.3 STEP→gmsh の「同じ構造」のメッシュ

- **構造ファイル** `mesh_structure.json` は run の `nozzle.msh` から作る (メッシャを再実装しない。メッシャはブランチごとに違い、今後も変わるため)。
  中身: `ni`・`nj`・station $x_i$ [mm]・半径方向の比 $s_{ij}=r_{ij}/r_{i,n_j-1}$ (全 station で共通なら 1 行)・四角形の接続の規則
  (`mesh2d` の $(i,j)\to i\,n_j+j$、反時計回り)・物理名とタグ・単位。
- **縦線の検査**: 各 station で $\max_j |x_{ij}-x_{i0}| \le 10^{-9}\,L$ ($L$ = 全長)。外れたら構造ファイルを作らず理由を書く (壁法線の近壁層など)。
  接続が規則と違う (節点番号の並びが $(i,j)$ 規則でない) ときも拒否する。
- **再構築スクリプト** `rebuild_from_step.py` (gmsh Python API、`.venv-mesh`) は単体で動く (forge に依存しない):
  1. STEP を `occ.importShapes` で読み、曲線が 1 本で助変数 $u$ の範囲が構造ファイルの $[x_0, x_{n_i-1}]$ を含むことを確かめる。
  2. 軸・入口・出口の直線と面を OCC で作る (壁は STEP の曲線そのもの)。
  3. 節点を作る: **$x(u)=x_i$ を $u$ について解く** (`wall_step.py` が書いた曲線は $x(u)=u$ だが、CAD で編集した曲線では成り立たない;
     codex plan レビュー M4 の反例で 0.25 mm ずれる)。事前に曲線を密に評価して、平面 ($|z| \le 10^{-9}$ mm)・$r>0$・$x(u)$ の狭義単調
     (多価なら拒否)・全 station を含む、を検査する。解いた $u_i$ で $r_w = r(u_i)$、$r_{ij}=s_{ij}\,r_w(x_i)$。
     節点・要素は `addNodes`/`addElements` で各 CAD 要素 (曲線・面) に載せる (壁の節点は解いた助変数 $u_i$ つき)。
  4. 物理グループ (inlet/outlet/wall/axis/fluid、タグは元と同じ) を付け、msh4.1・msh2.2・CGNS を書く。
  5. 元の `nozzle.msh` があれば節点ごとの距離と接続の一致を検査して JSON に書く (§6 V3)。
- STEP の曲線を CAD で修正した場合も、同じ構造ファイルで作り直せる (station と比は保ち、壁だけ新しい曲線に沿う)。
  ただし修正後の壁で格子品質が保たれるかは `check_mesh_quality.py` で確かめる (README に書く)。

### 4.4 Fluent 形式のメッシュと境界条件表

- **書き出し器** (`export/fluent_msh.py`): Fluent の ASCII メッシュ (節 0 コメント・2 次元・10 節点・12 セル・13 面・45 ゾーン名)。2D (`(2 2)`)。
  - ゾーン: 1 = fluid (セル、型 1)、面ゾーンは interior (型 2)・wall (3)・pressure-inlet (4)・pressure-outlet (5)・axis (37)。名前は msh の物理名
    (`inlet`→`inlet` ほか、内部面は `interior-fluid`)。
  - セル型 3 (四角形)、面型 2 (線分)。座標 [m] を 17 有効桁。
  - **向き**: 面 $(n_0,n_1,c_0,c_1)$ は $n_0\to n_1$ の進行方向の**左側が $c_0$** (§3 の 2D の実物の規則)。
    境界面は $c_1=0$ で、境界を内部が左に来る向き (領域を反時計回り) にたどる。
- **gmsh 経由の CGNS** も出す (Fluent の `File > Import > CGNS`)。境界の型は Fluent 側で付け直しが要ることを README に書く。
- **境界条件表** `boundary_conditions.md`: `bcondConfig.yaml` と問題 YAML から、Fluent の 2D 軸対称・$x$ 軸が対称軸・単位 m・
  **`Operating Pressure = 0 Pa` を明示し、各圧力を絶対値と Fluent の入力欄の値 (= 絶対値 − 動作圧 = 絶対値) で併記** (Fluent の圧力入力は
  動作圧に対するゲージ値; codex plan レビュー M5)・入口 = pressure-inlet (全圧・全温、初期ゲージ静圧の目安 = 入口の等エントロピー静圧)・
  出口 = pressure-outlet (背圧)・壁 = 等温 $T_w$ (または断熱)・
  軸 = axis、ガスの組成と forge 側の物性モデル (NASA-9・輸送物性) を書く。Fluent 側の物性の選び方は書かない (forge と同じにする方法を断定できないため、差の要因として注記)。
- forge は節点中心 (節点 = 未知数)、Fluent はセル中心。同じ一次メッシュでも第 1 未知数の壁距離は Fluent が約半分になることを README に書く。

### 4.5 Excel の表

- **形状表** `wall_shape_mm.csv` (CAD・加工・1 次元熱伝導の幾何): 1 mm ごと + スロート付近 ($|x| \le 50$ mm と $|x|\le 2r_t$ の広い方) は 0.1 mm ごと。
  列 `x_mm, r_mm, wall_angle_deg, curvature_radius_mm, s_mm (弧長、入口端から), perimeter_mm (2πr), lateral_area_cum_mm2 (∫2πr ds)`。
  曲線は `wall_repr.json` の物理壁 (STEP と同じ B-spline)。
- **CFD 壁面量表** `wall_cfd_nodes.csv` (壁節点ごと、最終の `res_wall`): `x_m, r_m, s_m, ring_area_m2 (2π r ds の節点持ち分), p_w_Pa, T_w_K,
  tau_w_Pa (接線成分), q_into_wall_W_m2, Q_node_W (= q × ring_area), y1plus` 。
  - **熱量は積分済みの `iface_Qf_eff` から作る** (codex plan レビュー M1): 全周の節点熱量 $Q_{node} = 2\pi\,$`iface_Qf_eff` [W]、
    熱流束 $q = Q_{node} / A_{node}$、$A_{node} = 2\pi\sum_f r_f |S_f|$ (その壁節点の壁面の半割面の $r$ 重み面積、メッシュ h5 の
    `PLANES/surfArea` と面重心から)。`iface_q_eff` (= $r\,q$) は使わない。
  - **状態を区別する** (codex plan レビュー M3): 各熱量の列に対し `status` を `valid` (正本、有限) / `not_applicable` (データセットはあるが
    断熱壁・非定常などで全 NaN) / `missing` (データセットなし) / `diagnostic` (`-qwall` を使った) / `nonfinite` (一部 NaN) で返す。
    **正本が無効なとき黙って診断値に置き換えない** — 診断値を出すのは列名 `q_into_wall_diag_W_m2` と状態 `diagnostic` を付けたときだけ。
  - **熱設計に使える判定** (`thermal_use`): 対象 snapshot と同じ区間の `check_convergence.py` (段の manifest があれば `--segment`)、
    壁の総熱量の時系列 (全 `res_wall` から) に `check_quasisteady.py --series-csv`、`check_wall_resolution.py` の 3 つが合格し、熱量の状態が
    `valid` のときだけ `PASS`。それ以外は `未評価` か `FAIL` を MANIFEST・README・xlsx の注記・CSV と同じ名前の `.status.json` に書く
    (CFD 前・未収束の書き出しも許すが状態を必ず付ける)。使う `res_wall` の step は 1 つに揃え、異なる step の混入は例外。
  - $y_1^+$ は `check_wall_resolution.py` と同じ定義 (第 1 内部節点の距離・接線応力)。ソルバの `ypls` は載せない。
- **列名は ASCII** (日本語の Excel が UTF-8 の CSV を Shift-JIS と読んで化けるのを避ける)。説明は README と `.xlsx` の別シートに置く。
- `.xlsx` (xlsxwriter があれば): シート「形状」「壁面量」「条件」と、1 次元の壁伝導の**記入用シート** (厚さ $t$・熱伝導率 $k_s$ を入れると
  $T_{back} = T_w - q\,t/k_s$ を列で出す式)。式は等温壁 ($T_w$ 固定) の単純な抵抗の式だけで、連成が要るなら `cht_loop` + `solid_shell` を使うと書く。

### 4.6 skill と手順書

- `nozzle-design`: 「ユーザの言葉 → カードのキー」の対応表 (入口長さ = 直管 `pipe_length_m` と縮流部 `contraction_length_m` の両方、
  出口長さ = 設計スロート→物理出口 `divergent_length_m`、スロート径は物理か設計かを必ず確認)、カード → YAML → 既存の §1・§2 の流れ、
  最後に `deliverables` を回して README のパスを報告する、を足す。「未実装」の節を更新する。
- `design-intake`: 成果物の要否 (STEP・Fluent・Excel) と Fluent 側の用途 (同じメッシュで比較するか) を聞く項目を足す。
- `procedures/nozzle-design-workflow.md` §1 ⑤「形を見る・出す」を成果物コマンドに置き換え、`nozzle-design-outputs.md` に成果物の節を足す。

## 5. 実装ステップ

1. 要求カードの生成器と雛形 (`design/forge_design/intake/`)、単体試験 (`design/tests/test_wt_card.py`)。
2. 成果物の書き出し: 形状表・3D STEP/STL (FreeCAD ジョブ)・構造ファイル・再構築スクリプト・Fluent 書き出し・CGNS・壁面量表・README/MANIFEST
   (`design/forge_design/export/`)、試験 (`design/tests/test_deliverables.py`、Fluent の向きの測定スクリプト)。
3. 検証 §6 (case/45 の生産の壁 + M5 デモ)。
4. skill・手順書・methods の同期、case README の記録。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | plan 段の諮問とレビュー | §4・§6 を codex (diagnose; 判断役のスイッチが codex) に諮り、続けて `codex_review.py --stage plan`。採否を §6.1 に書く。判断: 2026-10-13 諮問済み (寸法の出荷ゲートを採用) と plan 段レビュー済み (GO-with-changes、全件採用)。§6.1 | F |
| 2 | ~~要求カードの生成器~~ **完了 2026-10-13** | §4.1 の 1〜6。触る: `design/forge_design/intake/{__init__,wt_card}.py`・`intake/recipes/*.yaml`・`design/tests/test_wt_card.py`。合格: §6 V1 | O |
| 3 | ~~形状表と 3D STEP/STL~~ **完了** | §4.5 の形状表、§4.2 の cad/。FreeCAD ジョブを `export/` に一般化 (case/45 の `cad_deliverables.py` の surface を移す)。合格: §6 V2 | O |
| 4 | ~~構造ファイルと STEP→gmsh 再構築~~ **完了** | §4.3。触る: `export/mesh_structure.py`・`export/rebuild_from_step.py`。合格: §6 V3 | O |
| 5 | ~~Fluent 書き出し器と CGNS・境界条件表~~ **完了** (Fluent での読み込みは #10) | §4.4。触る: `export/fluent_msh.py`。合格: §6 V4a・V4b | O |
| 6 | ~~壁面量表と xlsx~~ **完了** | §4.5。合格: §6 V5 | O |
| 7 | ~~一括コマンドと README/MANIFEST・寸法ゲート~~ **完了** | §4.2・§4.2a。合格: §6 V6・V7 のゲート | O |
| 8 | ~~M5 デモ (CFD 前)~~ **完了** (case/70 run_0001; V7 の対照 run_0007) | カード (**Tt 1250 K (ユーザ指示 2026-10-13)、Pt 3.0 MPa (出口静圧が下がりすぎないよう親が選定: semiperfect 空気で出口 221.8 K・5139 Pa)**、実条件ではない) → `--solve-rt-throat` → `--init-integral --prepare-only` → 成果物。壁 350 K は断熱壁温 ~1150 K に対し強い冷却 (CONTUR の温度形の熱閉包の危険、isothermal plan §5.1 #10) なので、CFD 前の壁の精度は M6 の断熱の実績からは言えない。合格: §6 V7 | O |
| 9 | ~~skill・手順書・methods の同期~~ **完了** (V8 は未実施のまま記録) | §4.6。合格: `check_plans.py`、skill を読んだだけの別エージェントが V7 を再現できる (§6 V8) | O |
| 10 | Fluent での読み込み確認 (ユーザ) | `fluent/nozzle_2d_axi.msh` を Fluent で `Mesh > Check` (負の体積なし・ゾーン名と型)。結果を §9 に記録。**2026-10-13 ユーザ: 手元に Fluent が無いので先に進める** → 未確認のまま保留 (Fluent が使えるようになったら確認する。読めなければ CGNS が代わり) | ユーザ |
| 11 | result 段の解釈とレビュー | V1〜V8 の結果を codex (diagnose) に諮ってから `--stage result` | F |
| 13 | ~~変換器のタグ順の検査~~ **完了** | 作り直したメッシュ (gmsh が CAD 要素ごとに節点を書き、元の番号を保つ) を `convertGmshToForge` に通すと、変換器 (`mesh/gmshReader.hpp`) が節点のタグを読まずにファイル順で番号を振るため接続が黙って壊れた (skew 1.0 が 16.5 %、元は最大 0.44)。**変換器にタグ順の検査を追加 (2026-10-13、1 ファイルのバグ修正: 正しい msh の結果は不変、壊れた入力は止まる)**。再構築は forge 形 (節点ブロック 1 つ・元と同じ番号) と CAD 形 (`_cad.msh`) を分けて書く。合格: 元と作り直しの両方が `check_mesh_quality.py --ar-max 5000` で PASS (V3(d)) | O |
| 14 | 壁の端の節点の熱量 | case/63 run_0013 の加熱壁の両端 (断熱壁と接する節点) で q / (−qwall) が 3.81・1.93 (内部は 1.8〜3.2 % の差)。原因は未確認。ノズルでは入口・出口の角に当たる。表は値を変えず `.status.json` の `walls.*.end_nodes` に隣の節点との比を記録し、README に「壁の両端の節点の熱量は信頼できない (未確認)」と書く。原因の切り分けは CHT 側の論点として別に | F |
| 15 | **(別 plan に切り出す) 出口の切断と from_length の解法** | codex 諮問 2026-10-13 (2 回目) の採用。(1) **[2026-10-13 implementer の測定で訂正]** 今の雛形の閉包 (`wall_mode: cplus`) では $x_F$ は終端 C⁺ 自身の壁点の上に来る (`moc_inverse.py:609`) が、丸めで $x_F$ がその点より 1e-13 程度下になると切る規則 `wall[wall[:,0] <= x_F]` (`moc_inverse.py:613`) がその点を落とし、壁 = メッシュ = STEP の終点 $x_e$ が**ちょうど 1 点間隔**短くなる (起票時に引いた `:903-904` は streamline の経路)。M5 デモの from_length の 8 回の評価のうち 5 回で落ち (0.845 mm)、case/45 の実寸では拡大部長の誤差 −6.06 mm のうち 5.35 mm (0.0698 $r_t$) がこれだった。`explicit` も同じ経路 (`runner_axismach.py:554`) なので**生産の壁 (case/45) にも当たる** (実害は未確認)。交点 $(x_F, r_F)$ を有効な壁区間で補間して追加するのが方向 (最後の点の $x$ の書き換えは却下、$x_e$ を逆問題の目標にするだけの修正も却下)。(2) from_length は最初の非改善で打ち切る。判別 A/B (同じ MOC 場で切断 vs 交点追加、出口移動が $r_t(x_F-x_e)$ と 0.05 mm 以内で一致するか) を含めて別 plan で。候補: `tooling-design-problem-campaign-recipe` §5.1 #15 (寸法固定) の前提作業として新規 plan。判断: 2026-10-13 codex 諮問 — 本 plan では扱わない | F |
| 16 | **格子の雛形は M6 系の最新の格子の結論に追随する** | 要求カードの雛形の `mesh` (と生成器の第 1 層の見積もり) は 2026-10-07 の生産 G1 (2000×97、壁に寄せた等比) を元にしている。冷却壁の格子はその後、別セッションで ni 4719 × nj 121 の壁法線の層つきの格子 (`cold_pair_mesh.py`、isothermal plan §5.1 #13〜#15) と、節点を減らして軸付近を細かくする core grid (`feature/nozzle-core-grid` の plan `tooling-nozzle-core-grid`、draft) へ進んでいる。**本 plan は格子の仕様を決めない** (ユーザ 2026-10-13: 本セッションはフローを流す役、細かい検証は別セッション)。core grid の plan が accepted になったら雛形と見積もりをその結論に差し替える。それまで生成器の格子は「初期値」と明記したまま | O |
| 17 | **Euler → NS の場の移植が化学種の照合で拒否される (フローの穴)** | 2026-10-13 case/70: 手順書 §2 ③ の `deltastar_loop --init-integral --ic-from <Euler run>` は、`interp_field.py` が species_hash の不一致で拒否する (Euler は輸送の記録が無く、NS は `transport: schema=forge_transport_v1` を持つ; 熱力学は同じ)。生産 (case/45) は NS の場から NS を始めるので当たらない。新しい設計の最初の NS は `--ic-from` なし (prepare_ns の 1 次元等エントロピー場 + 段階起動) で始める (case/70 run_0005 で起動確認)。直すなら化学種の照合で「輸送の記録だけの差」を移植可にする判断が要る (`thermophysics-solver-owned-species-db` の領分)。本 plan では手順書と skill に回避を書くだけ | F |
| 12 | (後続) 壁法線の近壁層の格子の構造ファイル | `wall_normal_layer` が main に入ったら、station の壁上の点と法線方向の距離で構造を表す版 (v2)。それまでは拒否 | F |

## 6. 検証

合格条件は実装前に固定する (結果を見て動かさない)。

- **V1 カード**: 単体試験 — (a) 正例 2 (物理スロート径・設計スロート径)・負例 5 以上 (未知キー・欠けたキー・$R_{in}\le r_t$・$\mu>20$・スロート径の二重指定)、
  (b) 単位換算の往復: 生成した YAML の値 × $r_t$ が実寸と相対 1e-12 で一致 (**単位換算の試験であって形状の合格根拠ではない**; 形状は V7・§4.2a)、
  (c) 固定点: case/45 の生産 YAML の実寸 (入口半径 0.5 m 相当・`L_U`・`L_pipe`・物理スロート径) をカードにして、`solve_rt_throat` の固定点が
  10 回以内に $|\Delta r_t|/r_t\le10^{-6}$ で止まる。
- **V2 STEP 3D**: FreeCAD の読み直しで面が有効・1 面、gmsh (OCC、FreeCAD と別経路) の読み直しで回転面の面積が $\int 2\pi r\,ds$ (wall_repr の曲線を
  20001 点で数値積分) と相対 $10^{-6}$ 以内、**面積は FreeCAD ではノット区間ごとに切った面の和、gmsh では面の助変数空間の求積で測る**
  (OCC の `Shape.Area`/`getMass` は長い B-spline の 1 枚の回転面で −4.12e-4 外れた (M6) ので使わない; 2026-10-13 実装で判明し測り方を明記、
  許容は変えていない。既存の `wall_step` の添え書きの `revolve.area_mm2` も同じ粗い値なので面積として引用しない)、bbox の $x$ 範囲が曲線の定義域と 1e-6 mm 以内。STL の三角形の弦誤差の設定値を README に書く。
- **V3 STEP→gmsh**: case/45 の生産の壁 (`_band_ab/prod_confirm/prep`、wall_repr + 生産の格子パラメータで `mesh2d` から作った msh) と M5 デモで、
  再構築したメッシュの (a) 四角形の接続が元と同一 (全要素で節点番号の並びが一致)、(b) 節点の距離 $\max|\Delta| \le 10^{-5}$ mm かつ
  各 station の第 1 層厚の 1 % 以下、(c) 物理グループの名前・タグ・要素数が一致、(d) 元と再構築の**両方**が `check_mesh_quality.py` の明示した閾値 (`--ar-max 5000`) で PASS。
  負例: 縦でない格子 (人工的に 1 節点を $x$ 方向にずらした msh) は拒否。**助変数を付け直した同じ形状の STEP** (gmsh で点列から作り直した
  曲線、$x(u)\ne u$) は正例として同じ節点になる (≤ 1e-5 mm)。$x$ が折り返す曲線は拒否。
- **V4 Fluent**: (a) 自前の検査 — **書き出した面について、実物の 2D の Fluent メッシュ (§3 の elbow.msh) と同じ測定スクリプトで「$c_0$ が左側」の割合が
  100 % であること** (書き出し器の規則を使わず幾何から測る。elbow.msh 自体も同じスクリプトで 100 % 左になることを試験に入れる)、
  面積の総和が壁と軸で囲む面積と相対 1e-12、境界ゾーンの面数が msh の物理グループの辺数と一致、境界がつながった 1 本の閉路。(b) 独立の読み直し — ParaView (`pvpython`) の Fluent 読み込みで節点数・セル数・ゾーン名が一致、節点座標が msh の座標の float32 丸めと一致
  (`vtkFLUENTReader` は点を float32 で持つ; 起票時の「相対 1e-15」はこの読み込み器では原理的に満たせないので 2026-10-13 に訂正。17 桁の照合は自前の
  読み込みで行うが、書き出し器と同じモジュールなので独立ではない)。CGNS は gmsh で読み直して節点数・要素数・物理名が一致 (物理タグの番号は gmsh が振り直すので名前で対応付ける)。
  (c) Fluent での確認はユーザ (§5.1 #10)。結果が出るまで「Fluent で読めることは未確認」と README と応答に書く。
  (d) 境界条件表: 動作圧 0 Pa の明記、bcond の各圧力 (絶対) と表の入力値が一致 (自動検査)。
- **V5 Excel**: 形状表のリング面積の総和 (台形) と V2 の面積が相対 1e-6、CFD 表の `ring_area` の総和と壁の多角形の側面積が相対 1e-12。CSV の列名がすべて ASCII。
  **熱量の単位と符号** (codex plan レビュー M1): 軸対称で `interfaceDiag: 1` の等温壁の run (手元の case/63 系など; 無ければ小さい
  軸対称の管の run を作る) で、表の $q$ (= $2\pi$`iface_Qf_eff`$/A_{node}$) と物理量の `qwall` (反転) と `iface_q_compact` が、半径の違う
  2 箇所以上で 5 % 以内 (3 つは別の離散式なので一致はしない; 桁と符号の検査)、総熱量 $\sum Q_{node}$ と $\int q\,2\pi r\,ds$ が相対 1e-9
  **[2026-10-13 result レビュー M5: 不合格 → 再定義]** 実測 2.880e-6 (case/70 run_0006) で 1e-9 を満たさない。原因は面積の定義の差で、
  ソルバの面積 $A_{node}$ (半割面の 1/4 点の $r$ で重み) と幾何のリング面積 (多角形の台形) の節点ごとの比 ring/A − 1 が最大 2.1e-5・中央値 7e-8、
  その Q 重みの平均が 2.880e-6 で積分差を説明しきる (単位・符号の誤りではない)。1e-9 は「同じ面積なら一致する」自己一致を想定した値で、独立な面積では
  原理的に届かない。**再定義**: 正本は保存形の $Q_{node}$ と $q = Q/A_{node}$ (q × A_node は Q に厳密に戻る)。幾何の面積との差は
  |Σq·ring − ΣQ| / ΣQ ≤ 1e-4 (熱流束の形式間の差 1〜3 % より 2 桁以上小さい; 工学的に無視できる) を合格とし、実測値と節点ごとの比の最大を記録する。
  結果を見てから許容を変えたことをここに残す (codex result レビューの推奨どおり、自己一致でない検査に置き換えた)、
  符号は流体→壁が正。**$r$ 倍のずれ (旧 `iface_q_eff` をそのまま使う) を負例として検出する**。
  **状態** (M3): 断熱壁の run (データセットはあるが NaN) で `not_applicable`、一部 NaN を入れた写しで `nonfinite`、異なる step の混入で例外、
  未収束の run で `thermal_use` が `PASS` にならない。
- **V6 一括**: M5 デモと case/45 で 1 コマンドが成功し、MANIFEST の sha256 が実ファイルと一致、欠けた入力 (wall_repr の無い run) で該当成果物だけ
  飛ばして理由が書かれる (負例)。
- **V7 M5 デモ (CFD 前)**: カード → YAML → 寸法 → 物理壁 → メッシュ → 成果物が通り、`check_mesh_quality.py` の VERDICT を記録
  (FAIL ならメッシュの見積もりを直す; CFD をしないので壁解像の判定は「未判定」と書く)。**§4.2a の出荷ゲートが STEP とメッシュから測った寸法で合格**
  (カードの既定公差: 長さ 1 mm・直径 0.02 mm)。あわせて codex 諮問の A/B として、`L_total_tol` を既定 0.05 $r_t$ にした YAML で同じ手順を回し、
  測ったスロート→出口長の残差を記録する (既定のままでは公差を外れることがあるかの確認。外れなくても §4.1 の写し方は変えない)。
  負例: カードの長さをわざと公差外の値にした成果物で、ゲートが FAIL を出す。
- **V8 skill**: カードの例と skill だけを渡した別エージェント (Opus) が、ユーザ風の依頼文 1 つから V7 と同じ成果物まで到達する (途中で手順書の外の判断をしたら記録)。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| plan | 2026-10-13 | [2026-10-13-tooling-nozzle-wt-deliverables-plan.md](../../notes/reviews/2026-10-13-tooling-nozzle-wt-deliverables-plan.md) | GO-with-changes, C0/M5/m1 | **全件採用** (根拠を当方で確認してから): M1 軸対称の `iface_q_eff` は $r q$ (`variables.cpp:599-604` と `conjugateWall.cpp:229` を確認) → §3・§4.5 を積分済み `iface_Qf_eff` からに、V5 に単位と符号の検査。M2 直管長の閾値測定が +15〜18 mm ずれる (接合点が $r'=r''=0$) → §4.2a を制御点とノットからの同定に、設計スロート径の比較を設計壁に。M3 NaN のデータセット・収束との結び付き → §4.5 の状態と `thermal_use`、V5 の負例。M4 CAD 編集後に $x(u)=u$ が崩れる → §4.3 で $x(u)=x_i$ を解く、V3 に正例・負例。M5 Fluent の動作圧 → §4.4・V4(d)。m6 誤差の量の混同・品質判定 → §4.2a・V3(d)。判断役 (codex) 自身の指摘で却下が無いため、採否の別途諮問は省略 |
| result | 2026-10-13 | [2026-10-13-tooling-nozzle-wt-deliverables-result.md](../../notes/reviews/2026-10-13-tooling-nozzle-wt-deliverables-result.md) | **NO-GO**, C0/M5/m2 | **全件採用** (codex がメモリ上の負例で再現したものを当方でも試験に): M1 寸法ゲートがメッシュ側の差を合否に使っていない → 入口半径・最小半径・直管の終端の区間・出口端・弦の形状誤差を必須の判定に (`dimension_gate.py`; 負例 5 本を `run_deliverables_tests.py` に)。M2 3D STEP の検査の失敗・欠落が合格になる → 全検査を必須に (`deliverables.part_cad`; FreeCAD の読み直し無効・面 2 枚・gmsh 失敗の負例)。M3 `thermal_use` が表の JSON と xlsx に伝わらない → `.status.json` の `thermal_use`・`use_for_thermal_design`、conditions、xlsx の notes と wall_1d に専用の欄 (case/70 run_0006 で FAIL が載ることを確認)。M4 公差を省いたカードでゲートが止まる → 生成器と同じ既定の公差で補う。M5 V5 の積分検査が自己一致で、独立な面積では 2.88e-6 > 1e-9 → 不合格を記録し、原因 (面積の定義の差) を測って説明、許容を根拠つきで再定義 (§6 V5)。m6 V7 の A/B を M5 で実施 (run_0007: 内部許容差 0.05 r_t で拡大部長 −1.907 mm → ゲート FAIL、カード由来 0.02 r_t の run_0001 は +0.066 mm で PASS)。m7 文書の状態の同期 (methods・手順書・plans/README・§7・§5.1) |
| diagnose (result 前) | 2026-10-13 | [2026-10-13-wt-deliverables-result-diagnose.md](../../notes/reviews/2026-10-13-wt-deliverables-result-diagnose.md) | 結論: accepted は保留し、総熱量が一定のまま局所の熱流束だけが変わる負例で `thermal_use` が誤って PASS しないかを確かめる。Major 1 | **採用**: `thermal_use` は `Q_total_W` しか準定常を見ていなかった (`wall_tables.py` の `--series-cols Q_total_W`)。用途 (1 次元の壁伝導) は局所の q を使うので、壁ごとに x の順で 16 区間に分けた熱量 `Qband*` を時系列に足し、全部の準定常を要件にした。判別 A/B (合成、`run_wall_tables_tests.py`): 腕 A (局所一定) PASS、腕 B (総熱量は相対差 0 で同じ・2 節点の熱量を ±20 % で交互) は旧判定なら総熱量が STEADY で PASS になるところ、新判定で FAIL (区間 0 DRIFTING 40 %・区間 6 OSCILLATING)。V5 の単位検査は未収束の snapshot でも有効 (桁・符号の検査) とする codex の判断、−qwall と q_compact は別の式 (`viscousFlux_d.cu:615-617,654`・`conjugateWall.cpp:264`) との指摘も記録 |
| diagnose | 2026-10-13 | [2026-10-13-wt-length-accuracy-diagnose.md](../../notes/reviews/2026-10-13-wt-length-accuracy-diagnose.md) | 結論: 寸法固定の本体修正は別 plan、今回は内部許容差の緩和を明示・記録しカードの公差で出荷ゲートが判定 (未達は FAIL)。Major 3 / Minor 1 | Major 3 件は採用 (§4.1 1.・§4.2a・§5.1 #15)。Minor (「YAML は必ず作れる」「explicit は影響しない」) は当方の前提の誤りとして撤回 — 緩和しても形状の成立は保証せず、explicit も同じ切断を通る (§5.1 #15)。判断役 (codex) 自身の結論で却下なし |
| diagnose | 2026-10-13 | [2026-10-13-wt-deliverables-diagnose.md](../../notes/reviews/2026-10-13-wt-deliverables-diagnose.md) | 結論: 最終の STEP とメッシュから寸法を測る出荷ゲートを入れる。Major 1 / Minor 1 | Major (V1(b)・V7 が自己一致で、`from_length` の許容差 0.05 $r_t$ が寸法を保証しない; `runner_axismach.py:571,591-611` を当方で確認) → 採用: §4.1 1. の `L_total_tol`、§4.2a、V1(b) の注記、V7 のゲートと A/B。Minor (V4(a) が書き出し器と同じ規則で自己検証) → 採用、かつ当方の測定で**規則そのものが誤り**と判明 (2D の実物で $c_0$ は左側、§3) → §4.4 を訂正し V4(a) を実物と同じ測定スクリプトでの判定に変更。V4(c) (ユーザの Fluent 確認) は必須のまま |

## 7. 影響範囲

- 新規: `design/forge_design/intake/`・`design/forge_design/export/{deliverables,mesh_structure,rebuild_from_step,fluent_msh,wall_tables,freecad_surface_job}.py`・`design/tests/`。
- 新規 (実際): `design/forge_design/intake/{wt_card.py, recipes/*.yaml, examples/m5_air_demo.yaml}`・`export/{deliverables, dimension_gate, mesh_structure, rebuild_from_step, fluent_msh, wall_tables, cad_check_gmsh}.py`・試験 6 本・`case/70.m5_air_wt_demo`。
- 既存の変更: `export/wall_step.py` (`write_surface` の追加のみ)・`export/freecad_wall_step_job.py` (mode surface の追加のみ)、skill 2 本・手順書 2 本・`methods/design/overview.md`・`design/CAPABILITIES.md`、
  **`solver_density_cuda/mesh/gmshReader.hpp`** (変換器: 節点のタグがファイル順に 1..N でない msh を止める検査。正しい msh の結果は不変; §5.1 #13)。
- ソルバの計算・設計チェーン本体・既存 run の結果は変わらない。

## 8. 完了条件

- §6 V1〜V7 が合格 (V8 は記録)。§5.1 #10 (Fluent での読み込み) は 2026-10-13 のユーザ判断で**未確認のまま進める**: accepted にしてよいが、
  README・skill・plan に「Fluent で読めることは未確認 (実物の 2D の Fluent メッシュと書式・面の向きを照合済み、ParaView で読み直し済み)」と残し、
  確認できたら §9 に記録する。
- skill と手順書から成果物コマンドに到達でき、`check_plans.py` が通る。plan 段・result 段の codex レビューの行がある。

## 9. 変更ログ

- 2026-10-12: 起票。ユーザの要求 (実寸入力・STEP・gmsh・Fluent・Excel・1 次元熱伝導、いろいろなマッハ数を Claude に作らせる) と回答
  (出口長さ = スロート→出口、入口長さ = 直管長と縮流部長の両方) を反映。Fluent の面の向きを case/30・32 の実物で測定 (§3)。
- 2026-10-13: codex 諮問 (§6.1) を反映 — 寸法の出荷ゲート §4.2a、`L_total_tol` をカードの公差から決める (§4.1)、V1/V4/V7 の改訂。
  2D の Fluent の面の向きを実物 (TGrid 2D の elbow.msh) で測り、起票時の推定 (右側 $c_0$) が逆だったので §3・§4.4 と methods を訂正。
- 2026-10-13: plan 段の codex レビュー (GO-with-changes, C0/M5/m1) を全件採用 — 軸対称の熱量の単位 (`iface_Qf_eff` から作る)・直管長の同定・
  熱量の状態と `thermal_use`・STEP の $x(u)=x_i$ の逆解き・Fluent の動作圧・誤差の量の分離 (§6.1)。status を in_progress に。
- 2026-10-13: M5 デモの条件をユーザ指示で Tt 1250 K に変更、Pt は「温度・圧力が下がりすぎないように」の指示から 3.0 MPa を選定 (出口 221.8 K・5.1 kPa、スロートの第 1 層の粗い見積もり ~80 nm = r_t の 3.2e-6、float32 で ~27 ulp)。
- 2026-10-13: 結合で見つけた変換器の欠陥 (節点のタグを読まずファイル順で番号を振る → gmsh が CAD 要素ごとに書いた msh で接続が黙って壊れる) を
  `solver_density_cuda/mesh/gmshReader.hpp` の検査で止めるようにした (§5.1 #13)。元の M6 メッシュは変換・品質 PASS のまま、gmsh が普通に作る msh も変換できることを確認。
  V4(b) の座標の許容を ParaView の float32 に合わせて訂正 (結果を見て緩めたのではなく、読み込み器の精度の限界)。
- 2026-10-13: #3・#6 の実装の返却を反映 — V2 の面積の測り方 (OCC の面積は −4.12e-4 外れるので使わない) を明記、壁の端の節点の熱量の異常を §5.1 #14 に、
  CHT ブランチでは `iface_q_eff` が修正済みであることを §3 に。
- 2026-10-13: codex 諮問 (拡大部長の精度) を反映 — 内部許容差の明示的な緩和と記録 (§4.1)、物理出口の定義と出口端の一致の検査 (§4.2a、`dimension_gate.py`)、
  出口の切断と from_length の解法は別 plan へ (§5.1 #15; 生産の壁にも同じ切断が当たる)。
- 2026-10-13: ユーザ指摘 (M6 の格子は他セッションで数十万節点に進んでいる; 本セッションはフローを流す役) を受け、格子の仕様は決めず最新の結論に追随すると §5.1 #16 に記録。
- 2026-10-13: ユーザ「手元に fluent ないから、先に進めてください」→ §5.1 #10 は未確認のまま保留、§8 の完了条件をそれに合わせて記録。
- 2026-10-13: case/70 で Euler → NS の場の移植が化学種の照合 (輸送の記録の差) で拒否されることを確認 (§5.1 #17)。最初の NS は等エントロピー場から段階起動で始める。最終の雛形の YAML は Euler の prepare が拒否する (physical_wall_repr・pw_upstream) ので、生成器が Euler 用の写し `<stem>_euler.yaml` も書くようにした (設計壁がビット同一、試験あり)。生成器の ni の見積もりに AR の余裕 3 % を入れた (case/70 で設計壁 4981 → 物理壁 5002.4)。
- 2026-10-13: **検証の結果** (§6; 試験は design/.venv-opt、cwd = 作業ツリー):
  - V1 カード: `design/tests/run_wt_card_tests.py` 80 ok / FAIL 0 (`--slow` で 81 ok: case/45 の実寸で from_length の緩和が記録され YAML ができる)。
    Euler 用の写し `<stem>_euler.yaml` の設計壁がビット同一 (差 0.0 r_t)。
  - V2 3D STEP: `run_wall_tables_tests.py` 58 ok。M6 の内面の面積は FreeCAD (ノット区間の和)・gmsh (助変数の求積) とも ∫2πr ds と相対 1.2〜1.5e-14
    (OCC の面積は −4.12e-4 なので使わない)。M5 デモで相対 2.2e-12。
  - V3 STEP→gmsh: `run_mesh_structure_tests.py` 49 ok。M6: 最大距離 1.391e-7 mm・第 1 層厚比 1.442e-4・四角形 191,904 一致・元と作り直し (forge 形) の
    両方が `check_mesh_quality --ar-max 5000` PASS。助変数を付け直した STEP で作り直せる・折り返す曲線は拒否。
  - V4 Fluent: `run_fluent_msh_tests.py` 47 ok。実物の 2D の Fluent メッシュ (elbow.msh) で c0 左 1,454/1,454、M6 で書いた 385,903 面が左 100 %、
    ParaView の読み直しで全セル正の面積、CGNS の読み直し一致、境界条件表の動作圧 0 Pa と圧力の一致。**Fluent での読み込みは未確認** (§5.1 #10)。
  - V5 Excel・熱量: 合成データで q が相対 2.2e-16 で戻り、旧定義は r 倍ずれる (負例)。**実 run (case/70 run_0006、軸対称・等温壁 350 K・r 25〜150 mm)**:
    表の q と −qwall・iface_q_compact の差 1.3〜2.6 % (許容 5 %)・Σ Q = Σ q·A・符号一致・iface_q_eff/q = r。状態の負例 (not_applicable・nonfinite・step 混入)。
    `thermal_use` は未収束の run で FAIL (run_0005・0006)。
  - V6 一括: `run_deliverables_tests.py` 17 ok (M6 で 1 分 48 秒、20 ファイル、MANIFEST の sha256 一致、2.2 mm ずらしたカードで FAIL、wall_repr の無い run は理由つきで飛ばす)。
  - V7 M5 デモ (CFD 前): case/70 run_0001 — カード → YAML (r_t 3 回で収束、物理スロート径の差 2e-9 m) → 物理壁・メッシュ (品質 PASS、AR 4841.8) →
    成果物 PASS (寸法ゲート: 拡大部長 +0.066 mm、他 ≤2e-6 mm、出口端のメッシュと STEP の差 3e-7 mm)。所要: カード 56 s、準備 21 s、成果物 23 s。
    from_length の内部許容差の緩和の A/B は case/45 の実寸で記録 (緩和 0.00652 → 0.05 r_t、寸法ゲートは拡大部長 −6.06 mm で FAIL を正しく出す; 内訳は §5.1 #15)。
  - V8 (skill だけで別エージェントが再現): 未実施。
  - 結合で見つけて直したもの: 変換器のタグ順 (§5.1 #13)、再構築の出力を forge 形と CAD 形に分離、生成器の AR の余裕、Euler 用の写し、Euler → NS の移植の拒否 (§5.1 #17、回避を手順書に)。
- 2026-10-13: result 前の codex 諮問を採用 — `thermal_use` に局所の熱量 (壁ごと 16 区間) の準定常を追加、判別 A/B の試験を追加 (§6.1)。
- 2026-10-13: result 段の codex レビュー (NO-GO, C0/M5/m2) を全件採用して修正・再検証 (§6.1)。V5 の積分検査は事前登録の 1e-9 に対し 2.880e-6 で不合格、
  原因を面積の定義の差と特定して許容を再定義 (§6 V5、実測 2.880e-6 ≤ 1e-4)。V7 の A/B を M5 で実施 (case/70 run_0007)。試験: deliverables 23 ok (負例 7 本を追加)・
  wall_tables 61 ok (判別 A/B を含む)。成果物は case/70 run_0001・run_0006 を新しいコードで作り直した (run_0006 の `thermal_use` FAIL が表の JSON・xlsx・conditions に載る)。
