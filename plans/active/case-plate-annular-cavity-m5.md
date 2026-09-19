# 超音速平板の環状深キャビティ内部温度場 (等温壁)

## メタ

- **area**: `その他 (解析ケース: boundary / turbulence / heat transfer)`
- **status**: `draft`
- **related_docs**:
  - [`methods/boundary.md`](../../methods/boundary.md) — 壁 / 等温壁 / slip / 入口分布プロファイルの現在仕様
  - [`methods/turbulence/theory.md`](../../methods/turbulence/theory.md) / [`implementation.md`](../../methods/turbulence/implementation.md) — SST (低 Re 壁解像) の現在仕様
  - [`procedures/recommended-settings.md`](../../procedures/recommended-settings.md) §1・§2・§3・§8 — 現行レシピ
  - [`procedures/divergence-and-startup.md`](../../procedures/divergence-and-startup.md) — 段階起動
  - [`procedures/inlet-profile.md`](../../procedures/inlet-profile.md) — `inletProfile` CSV の与え方 (流入 BL)
  - [`case/48.flat_plate_cooled_m4/README.md`](../../case/48.flat_plate_cooled_m4/README.md) — 等温壁 × 低 Re SST の起動レシピと実測 (本計画の直接の先行。**引用は量ごとに限定する** — §3.3)
  - [`case/37.pintle_nozzle/README.md`](../../case/37.pintle_nozzle/README.md) — FreeCAD → Salome (SMESH ViscousLayers) → msh4.1 → forge の 3D チェーン
- **related_plans**:
  - [`tooling-nozzle-isothermal-wall-chain.md`](tooling-nozzle-isothermal-wall-chain.md) — 等温壁 × 低 Re SST の検証根拠 (case/48) を与える計画
  - [`turbulence-node-sst-wallfunction.md`](../accepted/turbulence-node-sst-wallfunction.md) / [`turbulence-node-wf-representative-point.md`](../accepted/turbulence-node-wf-representative-point.md) — node 壁関数の代表点。**3D 角線は対象外**のため本計画は壁関数を使わない (§4.5)
- **created**: `2026-09-19`
- **owner**: `sano`

## 1. 目的

超音速流 (既定 M5) が流れる**断熱平板**に開いた**環状の深いキャビティ** (既定: 外径 50 mm・
内部円柱径 45 mm・深さ 50 mm → すきま 2.5 mm) の**内部温度場**を求める。キャビティ壁は等温
(既定 500 K)。完了時に得られる状態:

- (a) キャビティ内部の温度分布 — 深さ方向プロファイル (周方向の代表方位別)、対称面コンター
- (b) キャビティ壁 (外筒 / 内円柱側面 / 床) の熱流束と総入熱、内円柱上面は**別枠**で報告
- (c) 高温ガスの侵入深さ (定義は §4.7)、開口部の流入/流出質量流束
- (d) 上記が収束解・準定常であることのツール VERDICT、格子・領域・流入 BL の感度

**可変前提 (ユーザ指示 2026-09-19)**: 内部円柱は**下流側に偏心**させる場合があり (すきまが周方向に
不均一になる)、**すきま量**と**マッハ数**も変更しうる。したがって形状・条件・メッシュ・BC・評価を
**単一ソースのパラメータから駆動**するチェーンとして作る (§3.1・§4.1)。

## 2. スコープ

- **やる**:
  - パラメトリック形状生成 (FreeCAD): すきま・偏心・深さ・領域・M を設定ファイルで振れる
  - Salome (NETGEN + SMESH ViscousLayers) の非構造メッシュ (tet + prism) と品質ゲート
  - 流入乱流 BL を **2D 前駆計算 + `inletProfile`** で与える (§4.3)
  - forge (node, SLAU, 低 Re SST) の定常 RANS を段階起動で収束させる
  - 評価ツール (温度場・熱流束・侵入深さ・開口流束・収支) と、収束/準定常/感度の定量判定
  - 熱量的不完全性の決着: CPG と semi-perfect (TP) を Stage A で比較し**生産 EOS を確定**する
  - 定常・半割 RANS の妥当性確認として、Stage A メッシュで **URANS (dual-time) の非定常性チェック**
- **やらない** (別計画・将来):
  - 生産格子での本格 URANS / DES (Stage A の非定常性チェックで必要と判明したら別計画へ)
  - 放射伝熱、固体側熱伝導 (CHT)、壁温分布の連立 (等温を与件とする)
  - 実在気体 (解離・電離)。M5・$T_t \le 1300$ K の空気では不要
  - 形状パラメータの系統スイープ (チェーンは振れるようにするが、本計画で回すのは既定形状 + 感度のみ)

## 3. 関連 docs と前提

### 3.1 作動条件 — **パラメータであり、既定値は仮定である**

ユーザ指定は「M5」「平板断熱」「キャビティ等温 500 K」のみで、**高度 (静圧・静温) は決まらない**。
既定値として「高度 20 km 標準大気・M5」を採るが、これは**解析の仮定**であって導出ではない
(codex M7 指摘を採用)。値は `case/49.plate_annular_cavity_m5/conditions.json` を単一ソースとし、
`setup.py` が自由流・前駆計算の目標・BC 値・EOS 係数を導出する。M や高度を変えたらここだけ直す。

物性は forge 実装に合わせる: **Sutherland は `viscMethod: 1` の実装定数 $\mu_0$=1.716e-5,
$T_0$=**273.0** K, $S$=**111.0** K** ([`gasProperties_d.cu:59`](../../solver_density_cuda/cuda_forge/gasProperties_d.cu)。
計画の旧版は 273.15/110.4 と書いていた: codex m11)。`physProp.visc` は `viscMethod: 0` 専用なので
`physProp.visc` / `thermCond` は **`viscMethod: 1` でも parser が必須**とする
(値は Sutherland のとき使われないが、欠けると起動時エラー。旧版は「書かない」と誤記していた。
2026-09-19 実機確認)。自由流値を入れておく。$R = c_p(\gamma-1)/\gamma$ = 287.0 (case/48 と同じ)。

| 量 | CPG (γ=1.4, cp=1004.5) | semi-perfect (TP, 乾燥空気) |
| --- | --- | --- |
| $T_\infty$ | 216.65 K (20 km, 仮定) | 同 |
| $P_\infty$ | 5474.72 Pa (同) | 同 |
| $\rho_\infty$ | 0.0880483 kg/m³ | 同 (低温では TP≈CPG: $c_p$(216.65)=1002.5) |
| $a_\infty$ / $U_\infty$ | 295.042 / 1475.21 m/s | 同 (γ(216.65)=1.4012) |
| $\mu_\infty$ | 1.42178e-5 Pa·s | 同 |
| $Re/m$ | 9.1357e6 /m | 同 |
| $T_t$ | 1299.90 K | **1222.73 K** |
| $T_{aw}$ ($r=Pr^{1/3}$=0.8963) | 1187.55 K | **1126.33 K** |
| $c_p$(1300 K) | 1004.5 (定数) | **1188.31** (= +18.27 % vs 298 K), γ=1.3185 |
| $T_w/T_{aw}$ (500 K 壁) | 0.421 | 0.444 |

> 上表の TP 値は本計画で [`design/forge_design/gas/semiperfect.py`](../../design/forge_design/gas/semiperfect.py) の
> NASA-9 係数を使って再計算した実測値 (乾燥空気 $Y$ = N2 0.75518 / O2 0.23139 / AR 0.012885 / CO2 0.000545)。
> **CPG は高温端を $T_{aw}$ で +61 K (+5.4 %) 過大に出す**。本ケースの成果物が温度そのものなので、
> これは「後で 1 run 確認する感度」ではなく **EOS の選定問題**として Stage A で決着させる (§4.10, codex M8)。

その他の仮定 (結論に効くので明記する):

- 流入 BL は**全乱流**とし、$\delta$ を与件として選ぶ (§4.3)。自然遷移は扱わない。
- 内円柱上面 (`cyl_top`) は等温 500 K (温調体という解釈)。**侵入深さに効かないという根拠は無い**ので
  断熱版の感度 run を必須にする (§6, codex M7)。
- 入口乱流強度 TI 0.5 % ($\mu_t/\mu\approx10$)。前駆計算の入口値であり、キャビティ位置の $k,\omega$ は
  前駆解から受け継ぐ。

### 3.2 ツールチェーンと、その既知の制約 (実機確認済み 2026-09-19)

| 役割 | 実体 | 備考 |
| --- | --- | --- |
| 形状 | FreeCAD 1.1.1 `/home/sano/opt/squashfs-root/usr/bin/freecadcmd` | `QT_QPA_PLATFORM=offscreen`。**例外時に stdout が飛ぶので全 print に `flush=True`** |
| メッシュ | SALOME 9.14.0 `.../SALOME-9.14.0-native-UB24.04-SRC/salome --keep-paths -t` | `LD_LIBRARY_PATH=/home/sano/opt/salome/syslibs/usr/lib/x86_64-linux-gnu` 必須 |
| MED→msh4.1 | case/37 の `med_to_msh41.py` を case/49 用に一般化 (グループ→physID を JSON 化) | 体要素は全て単一 `fluid` に入るので**幾何が複数ソリッドでも可** |
| 取り込み | `solver_density_cuda/build-native/convertGmshToForge` | cwd の `solverConfig.yaml`/`bcondConfig.yaml` を読む。dual 閉性をログに出す |
| 計算 | **AWS g5.xlarge (A10G 23 GB, RAM 15 GB, 4 vCPU, 空き 44 GB) `ubuntu@52.3.224.27`** | ローカルは RAM 11 GB で 3D が落ちる ([aws-p1-instance-state])。**メッシュ生成はローカル** (Salome は AWS 未導入)、h5 を転送して計算 |
| aws CLI | `~/.local/bin/aws` 2.36.49 (2026-09-19 導入) | **認証情報未設定** (`~/.aws` 無し) のため IP 取得・start/stop は不可。IP はユーザから受ける |

**本計画で確認した既存ツールの制約** (codex C1/C2/M6 を実機検証):

1. **`interp_field.py` は 2 次元最近傍**: [`interp_field.py`](../../solver_density_cuda/tools/interp_field.py) の
   `centroids()` が座標を `[:, :2]` に切って cKDTree を作る。深さ 50 mm のキャビティでは z が無視され
   $z=-49$ mm と $z=-1$ mm が同じ供給点を拾う。**→ 本計画は 3D cross-mesh restart を使わない**
   (Stage A と Stage B はそれぞれ段階起動で一様 IC から立てる。§4.6)。
2. **低 Re 壁の壁面ダンプに熱流束が入らない**: `outputHDFflg: 1` の `res_wall_<id>_*.h5` を実測すると
   (case/48 run_0005) `qwall` / `utau` / `ypls` は**全点 0**、`twall_x/y/z` (壁せん断) と `Ps`/`Ts` は有効。
   **→ $q_w$ は場の温度勾配から評価する** (case/48 `tools/cooled_plate_eval.py` と同じ方式)。$\tau_w$ は
   ダンプの `twall_*` を使える。
3. **`check_mesh_quality.py` の 3D 指標は辺長比 AR と面内角 skew のみ**で、体積が潰れた sliver を
   検出しない (codex の人工 tet 例: 厚さ 1e-6 で AR 1.414 / skew 0.250 で PASS)。**→ 体積・sliver 判定を
   case 側ツールで足す** (§4.4)。

### 3.3 先行 case/48 の到達点 — 引用は量ごとに限定する (codex M10 採用)

case/48 の run は `check_convergence.py` で **`NOT CONVERGED (stalled/plateau)`** であり
(run_0004 / 0005 / 0007 / 0008 / 0011 を codex が再実行して確認)、旧版計画が書いた
「VD-II/SU2 と ±3 %/±1 %」という一括引用は**成立しない**。本計画が引用してよいのは:

- 冷却壁 (300 K) の $C_f$/VD-II 0.97–0.99、$2St/C_f$=1.16、SU2 比 $C_f$ 1.000–1.001 / $q_w$ 1.008 /
  $\theta$ 0.998 / $\delta^*$ 1.037–1.047 (run_0011 ↔ su2_B, 準定常 STEADY の量)
- $y_1^+$ 掃引: $q_w$ は $y_1^+\le1$ で −1.3 %、$\le1.8$ で −3.7 %、$\le3.6$ で −7.0 % (run_0007–0009)
- 起動レシピ (層流暖機 → SST soft → mid → 2 次ランプ → 本段) が M4.19 冷却壁で通ること

**「残差プラトーでも準定常な積分量は SU2 と一致する」**という限定付きの実績として扱う。

## 4. 設計方針

### 4.1 形状 — パラメトリック (実装済み・検証済み)

`case/49.plate_annular_cavity_m5/cad/build_geom.py` (FreeCAD)。`geom_config.json` で上書き。

- 座標系: $x$=流れ方向、$z$=平板法線 (流体 $z>0$)、$y$=スパン。外筒軸を原点、平板面 $z=0$。
- パラメータ: `Ro` (キャビティ半径, 既定 25)、`Ri` (内円柱半径, 既定 22.5)、**`x_off`** (内円柱の
  $x$ 偏心, 既定 0)、`depth` (既定 50)、領域 (`x_in/x_plate/x_out/y_max/z_top`)、`r_patch` (平板の
  細分パッチ半径)。すきまは $R_o-R_i$、偏心時は周方向に $R_o-R_i-|x_{off}|$ 〜 $R_o-R_i+|x_{off}|$。
- **$x$ 方向の偏心は $y=0$ 対称を保つ**ので半割 ($y\ge0$) を維持できる。$y$ 偏心は半割不可
  (スクリプトは $x$ のみ許す)。
- 構築: 環状リングを箱内へ `protrude` 突出させて fuse (同一平面 $z=0$ の coincident-face fuse 回避)、
  箱を `x_plate` で 2 分割 (助走 slip / no-slip 平板を別面にする)、`r_patch` 円柱を箱内部に fuse
  (体積不変で $z=0$ 面に seam を残し、開口周りだけ細分できる面 `plate_in` を作る)、最後に $y<0$ を cut。
  **`removeSplitter` は呼ばない** (呼ぶと seam が消えて面を分けられない)。
- **検証 (codex M4 採用、実装済み)**:
  - 単一閉ソリッド・bbox・総表面積が解析値と一致 (内部面の残留を検出)
  - **グループ別面積が解析値と一致** — 総面積一致では誤タグを検出できないため必須
  - 未分類 0 / 全 face がちょうど 1 グループ
  - 面分類は**重心を使わない**: 円筒面は `Surface`(Part.Cylinder) の**半径と軸**、平面は法線と
    支持位置、$z=0$ 面の内外は**頂点の最大半径**で判定する (半円筒面の重心半径は $2R/\pi$ で
    円筒上に無く、曲線境界を持つ平面の重心も境界円の内側に寄るため)
  - 実測: 既定形状 / 偏心 `x_off`=1 (すきま 1.5–3.5 mm) / すきま変更 `Ri`=23.5 + 助走ありの 3 通りで
    **全グループ面積の相対誤差 ≤ 1.2e-8** (円筒面の面積は tessellation 由来でこの桁)

### 4.2 計算領域と境界条件

流入 BL を前駆計算から与える (§4.3) ので、BL 発達用の長い助走を持たない。既定領域 [mm]:
$x\in[-80,120]$, $y\in[0,70]$, $z\in[0,60]$ + キャビティ ($z\in[-50,0]$)。

| physID | グループ | 位置 | BC (本段) | 値 |
| --- | --- | --- | --- | --- |
| 1 | `inlet` | $x=-80$ | `inlet_uniformVelocity` + `ints: {inletProfile: 1}` | ro/Ux/Ps/k/omega を $z$ 分布で (§4.3) |
| 2 | `outlet` | $x=+120$ | `outlet_statPress` | Ps=$P_\infty$ + 逆流 Pt=$P_\infty$/Tt=$T_\infty$ |
| 3 | `top` | $z=+60$ | `slip` | — |
| 4 | `side` | $y=+70$ | `slip` | — |
| 5 | `sym` | $y=0$ | `slip` | — |
| 6 | `plate` | $z=0$, $r>r_{patch}$ | `wall` (断熱 no-slip) | — |
| 7 | `plate_in` | $z=0$, $R_o<r<r_{patch}$ | `wall` (断熱 no-slip) | 開口周りの細分用に分離 |
| 8 | `cav_outer` | $r=R_o$, $-50<z<0$ | `wall_isothermal` | Ts=500 |
| 9 | `cav_floor` | $z=-50$ の環 | `wall_isothermal` | Ts=500 |
| 10 | `cyl_side` | 内円柱側面 | `wall_isothermal` | Ts=500 |
| 11 | `cyl_top` | $z=0$, 内円柱上面 | `wall_isothermal` | Ts=500 (断熱版の感度 run あり) |
| (12) | `runup` | $z=0$, $x<x_{plate}$ | `slip` | 既定では作らない (`x_plate`=`x_in`) |

- **変換は「最終形の壁タグ」で行う** (codex M5 採用): SST の `wall_dist` は変換時の no-slip 壁から
  作られるので、`plate`/`plate_in`/キャビティ 4 面を **`wall`/`wall_isothermal` にしたまま変換**し、
  段階起動 S0 の全面 slip は**実行時に `bcondConfig.yaml` を差し替えて**実現する。
- **`space.pRef` = 5474.72** (= $P_\infty$)。非直交・float32 の自由流保存に必須 (旧版で欠落: codex M5)。
- 出口は超音速だが node は壁列が常に亜音速なので `outflow` でなく $P_s$ 一致の `outlet_statPress`
  ([node-supersonic-exit-outflow], case/48 と同じ)。
- 領域寸法の根拠と**その限界**: M5 のマッハ角 $\mu$=11.54° から、キャビティ起源の**微小**擾乱は
  出口 ($\Delta x$=120) で $z=y=$24 mm に届く → `top` 60 / `side` 70 には余裕 2.5–3 倍。
  **ただしマッハ角は微小擾乱の角度**であり、有限強度の圧縮波や亜音速 BL 経由の上流影響は
  これでは排除できない (codex M7)。→ **領域感度を必須ゲートにする** (§6: `top` 60→90、
  `side` 70→100、`x_out` 120→180 の 3 通りで結論量の変化 ≤ 2 %)。

### 4.3 流入境界層 — 2D 前駆計算 + `inletProfile`

平板 BL を 3D 領域内で発達させると、(a) 平板面積が大きくなり prism 節点数が予算を超え、
(b) すきま 2.5 mm と外部 BL 5 mm で必要な prism 総厚が両立しない (§4.4)。よって:

1. **2D 前駆** (`case/49/precursor/`): case/48 と同じ平面 2D 構造メッシュ・同じ起動レシピで、
   同一自由流・**断熱**平板の M5 乱流 BL を計算する (node, 低 Re SST, $y_1$=6 µm)。
2. $\delta_{99}\approx$ **5 mm** (= 既定すきま幅の 2 倍) になる station $x_0$ で法線分布を抜く。
3. 3D は $x=-80$ mm からこの BL を流す。

- **CSV は case 側ツールが直接書く** (codex M3 採用、実機確認): `gen_inlet_profile.py gen --table` は
  `Ps` を出力列から落とす (`Ps` は `Tt`/`M` 換算経路でしか再生成されない) ため、そのまま渡すと
  **圧力分布が BC の一様値に戻る**。`tools/extract_profile.py` が
  **`z ro Ux Uy Uz Ps k omega` を直接** `inlet_profile_1.csv` に書く。
  - 2D 前駆の**壁法線速度 $V$ を 3D の `Uz`** に写す (発達 BL では $V\ne0$。旧版は欠落)。
    スパン方向は `Uy`=0。
  - 起動ログの `[applyInletProfiles] ... applied:` に **`ro Ux Uy Uz Ps k omega` が全て並ぶこと**を確認し、
    さらに `gen_inlet_profile.py verify` 相当で**入口断面の実際の $U,T,k,\omega$ が前駆解と一致**することを
    確認する (亜音速部の閉包は CPG と TP で違うので「内面で上書きするから成立する」で済ませない)。
- **3D 側の BL 検査は「下流発達を織り込んで」行う** (codex M4 採用):
  - 比較相手は前駆解の**同じ物理位置** — 3D の $x=-30$ mm は前駆の $x_0+50$ mm と比べる
    (旧版は $x_0$ と比べており、物理的な発達分を解像度誤差に混ぜていた。case/48 の断熱 run では
    $x$ 0.30→0.35 m で $\delta^*$ +11.9 % / $\theta$ +13.7 % 発達する)。
  - まず**キャビティを塞いだ同一方式の 3D メッシュ** (`cad/build_geom.py --no-cavity` 相当) で
    入口移送と BL 発達だけを検証する。ここで合わなければキャビティ計算に進まない。
  - 比較量は厚さだけでなく **$U,T,k,\omega$ の法線分布と壁摩擦 $\tau_w$** (壁ダンプの `twall_*` が使える)。
  - **外層・prism→tet 遷移を細分した差**を取って離散誤差を分離する。
  - **「BL 内 19 点」を先験的に合格とも不合格ともしない** — 上の比較で採否を決める (codex M4)。
- $\delta$/すきま幅は**支配パラメータ**なので、既定 2.0 に加えて **$\delta$≈2.5 mm (比 1.0)** の
  感度 run を行う (§5.1)。
- CSV は **run ごとに置く** (restart でも毎回読む)。

### 4.4 メッシュ方針

手順: `cad/build_geom.py` → STEP → `cad/mesh_salome.py` (NETGEN 1D2D3D + `StdMeshers_ViscousLayers`)
→ MED → `cad/med_to_msh41.py` → `forge.msh` → `convertGmshToForge` → `forge.h5`。
面グループは Salome 側でも `KindOfShape` (円筒の軸・半径) と頂点最大半径で判定し (§4.1 と同一規則)、
**グループ別面積を解析値と照合**してから mesh を切る。`geom_used.json` の派生量 (すきま min/max) を
読んでサイズを決めるので、偏心・すきま変更でサイズ設定が自動追従する。

- **prism 境界層 (ViscousLayers)**: Salome は 1 メッシュ 1 設定なので、**すきま最小値が総厚を縛る**:
  総厚 $T_{VL} \le 0.4\,g_{min}$ (両壁から張って中心に 20 % のコアを残す)。既定 $g_{min}$=2.5 mm →
  $T_{VL}\le1.0$ mm。採用: **第一層 10 µm、stretch 1.25、14 層 → 総厚 0.869495 mm、最終層 0.181899 mm**
  (等比和。旧版の 0.853/0.186 は転記誤り: codex m10。**層厚などの派生量は manifest から生成**して転記を無くす)。
  - $y_1$=10 µm の根拠: case/48 (M4.19, $T_w$=300 K, $y_1$=3 µm → $y_1^+$=0.44–0.48) を
    $Re/m$ 1.84 倍・$T_w$ 500 K ($\rho_w/\mu_w$ が 300 K 比 0.41 倍) でスケールすると
    $y_1^+\approx0.32\,(y_1/3\,\mu m)$ → 10 µm で ≈1。**実測で確認**し、$y_1^+>2$ なら細くする。
  - **$y_1^+$ は壁ノードではなく第一内部ノードで評価する** (node の壁ノードは $u$=0・$T$=$T_w$ に
    ピンされ $y_1^+$ の定義に使えない。codex M6)。
  - 偏心で $g_{min}$ が 1.5 mm に減ると $T_{VL}\le0.6$ mm → 層数を 12 に落とす (自動計算)。
- **すきま方向の解像**: 両壁 14 層 (10 µm→0.1819 mm, 計 1.739 mm) + コア tet 0.76 mm → 幅 2.5 mm を **計 ~30 セル**。
  深さ 50 mm を ~50 セル、周方向 (半周 74.6 mm) を ~65 セル。
- **外部 BL の法線解像**: prism は 0.8695 mm ($y^+\approx$87) までしか
  覆わない。残り ($\delta$=5 mm まで) は `plate_in` の面サイズ 0.8 mm の tet が担い、**BL 内 計 ~19 点**
  (prism 14 + tet 5, prism→tet の厚さ比 0.1819→0.8 = 4.4) になる。**入口側の `plate` は 2.5 mm なので
  そこはさらに粗い**。十分かどうかは先験的に決めず:
  - §4.3 の**塞ぎメッシュでの前駆比較 (分布・$\tau_w$ 込み) と細分差**で採否を決める
  - 足りなければ `plate_in` を 0.5 mm に細分 (節点増) か、**箱とすきまを別ソリッドに分割して
    ViscousLayers を 2 系統にする** (SMESH の sub-mesh。§5.1 #10 に退避策として置く)
- **VL 総厚 ≤ 壁面の接線セルサイズ (必須ガード、2026-09-19 実測 4 例)**: 凸角 (開口リップ) まわりの
  接線セルサイズが VL 総厚を下回ると**層が自己交差し、NETGEN の tet 充填が無言で失敗する**
  (`Compute: True` なのに tet 0 / 体積の 2.2 % しか埋まらない。警告は層厚未達のみ)。

  | リップ エッジ | 壁面サイズ | VL 総厚 | 結果 |
  | --- | --- | --- | --- |
  | 0.30 mm | 0.8–1.2 mm | 0.8695 mm | **tet 0 (失敗)** |
  | 0.50 mm | 1.2–2.0 mm | 0.4257 mm | 成功 |
  | 1.00 mm | 0.8–1.2 mm | 0.8695 mm | 成功 |
  | 0.90 mm | **0.5 mm** | 0.8964 mm | **tet 0 (失敗)** ← エッジだけの下限では不十分 |

  → `setup.py` が **壁面サイズ全部 (`size_plate_in`,`size_cav`,`size_floor`,`size_cyl_top`,
  `size_edge_opening`) に下限 = VL 総厚**を課し、引き上げた面を manifest の `size_bumped_by_vl` に残す。
  **帰結 (トレードオフ)**: キャビティ内面を細かくするほど VL 総厚も薄くせざるを得ず、
  **継ぎ目比 = 最終層厚/壁面サイズ ≤ (s−1)/s** で頭打ちになる。第一層 $y_1$ は独立なので $y^+$ は保てる。
- **継ぎ目の滑らかさ (ユーザ要件 2026-09-19)**: prism 最終層と隣の tet のサイズ差を manifest の
  `junction_ratio` (= 最終層厚 / 最小壁面サイズ) で管理し、NETGEN の `SetGrowthRate`
  (manifest `growth_rate`, 既定 0.3 → **0.15**) で壁から離れるときの粗大化を緩める。
  実績: Stage C 0.145 (7 倍段差) → Stage D **0.239** + growth 0.15 (AR 最大 382→128, 平均 6.1)。
- **サイズ設定 (既定)**: 大域 maxh 6 mm / minh 0.02 mm、`plate` 2.5 mm、`plate_in` 0.8 mm、
  開口円エッジ ($r=R_o, R_i$ @ $z$=0) 0.3 mm、`cav_outer`/`cyl_side` 1.2 mm、`cav_floor` 0.8 mm、
  `cyl_top` 1.2 mm。
- **節点数**: Stage A ~2e5 (VL 8 層・粗サイズ)、Stage B 目標 5e5〜8e5。**上限は Stage A の実測
  (GPU メモリ・host メモリ・s/step) から確定する** (codex M6。A10G 23 GB / RAM 15 GB が上限)。
- **品質ゲート** (codex M9 採用: 要素種別ごとに指標・閾値・不合格条件を確定する):
  - `check_mesh_quality.py` (3D モード): skew ≤ 0.9、AR ≤ 1000。**AR 緩和 ≤ 5000 は
    「壁法線 prism 層のセルだけ」に適用**し (tet には適用しない)、緩和対象セル数を台帳に書く。
  - **追加 (case 側 `tools/check_mesh_extra.py`)**:

    | 対象 | 指標 | 閾値 / 不合格条件 |
    | --- | --- | --- |
    | 全要素 | 符号付き体積 | **> 0** (1 個でもあれば FAIL) |
    | tet | 正規化体積 $V/(\frac{1}{6\sqrt2}L_{rms}^3)$ (= 正四面体で 1) | **≥ 0.01** を全数。sliver 判定はこちら (旧版の $V/L_{max}^3$ は正常な薄 prism も小さくなり区別できない: codex M9) |
    | prism | 局所 Jacobian の**符号が 6 頂点で一様**、底面 skew ≤ 0.9、層厚比 (隣接層) ≤ 1.5 | 符号不一致 1 個で FAIL |
    | dual CV | **体積 > 0** 全数、**CV ごとの閉性** $\|\sum \mathbf{S}_f\|/\sum\|\mathbf{S}_f\|$ ≤ 1e-6 | 変換器のログは**全体最大面積で正規化**するので微小 CV の局所閉性を保証しない (codex M9)。CV ごとに測る |
    | 層あり率 | **壁 6 グループ** (`plate`,`plate_in`,`cav_outer`,`cav_floor`,`cyl_side`,**`cyl_top`**) | 100 % |

### 4.4.1 全ヘキサ (構造化) メッシュ — **生産メッシュはこちらに切り替える** (2026-09-19)

ユーザ質問「gmsh でヘキサに切れるか」への回答として作り、**tet+prism より全面的に良い**と判明した。
`cad/build_hex_mesh.py` (gmsh Python API)。

**構成**: z=0 面を四角形でブロッキングして押し出すだけ。

| 領域 | ブロック | 押し出し |
| --- | --- | --- |
| すきま (環状スリット) | 4 セクタ × 半径 2 層 = 8 面 | **下向き** (深さ方向。開口側と床側で別々に細分) |
| 内円柱上面側の外部流 | バタフライ (コア + 3 面) | 上向き |
| 開口まわり | 円→角の O グリッド 4 面 | 上向き |
| 遠方 | 矩形 3 面 | 上向き |

環状面を上下の押し出しが共有するので**開口は内部面**になり自動的に整合する。偏心は
外筒軸からのレイが内円と交わる半径 $r_i(\phi)$ を使うだけでトポロジは不変。

**実測比較 (同一形状・既定条件)**:

| | tet+prism (Stage D) | **全ヘキサ** |
| --- | --- | --- |
| 節点 / 要素 | 741k / tet 1.51M + prism 0.88M | 1.00M / **hex 961k (100 %)** |
| skew 最大 / 平均 | 0.849 / 0.214 | **0.500 / 0.050** |
| AR 最大 | 127.6 | 302.0 (壁法線の構造層。緩和枠内) |
| hex 正規化体積 (最小) | — | **0.744** (ほぼ直方体) |
| メッシュ生成 | 約 8 分 | **2.4 秒** |
| 1 step (node SST) | 39–64 ms (741k) | **24.9 ms (1.00M)** = 節点あたり約 2 倍速 |
| すきま横断 | ~28 セル (VL + tet コア、継ぎ目に段差) | 24 セル (両壁 20 µm から連続) |
| CV 閉性 max | 5.11e-6 (PASS) | 1.21e-5 (PASS) |
| 開口の正味/片道 質量流量 | 3.3〜7.3 % | **0.73 %** |
| CV エネルギー収支 残差 | 22 % | **2.3 %** (ゲート 5 % 内) |

**保存性はどちらも取れている**: CV 閉性はどちらも PASS で、差は**後処理の面積分の精度**に出る
(§4.8 の「収支残差の出どころ」)。ヘキサは評価面が押し出しの節点層と一致するため積分が素直。

**ユーザ確認済み (2026-09-19)**: 「ただヘキサで良いです」→ **生産メッシュはヘキサで確定**。
tet+prism 経路 (FreeCAD → SALOME) は残置するが以後の生産には使わない。偏心スイープ
`sweep_offset.py` も既定をヘキサ経路に差し替えた (FreeCAD/SALOME 不要 = **AWS 上で完結**)。

**tet+prism で抱えていた制約が構造的に消える**:

- レイヤー↔空間メッシュの継ぎ目が無い (半径方向の分布を直接指定)
- 「VL 総厚 ≤ 接線セルサイズ」の制約が無い (層が自己交差しない) → **内面をいくらでも細かくできる**
- **系統的な格子細分列が作れる** (`--scale` で全方向の分割数を一律倍、第一セルを 1/倍)。
  tet+prism では VL 総厚と接線サイズが結合していて独立に振れず、Stage A〜D が
  系統列にならなかった (§6.4 の格子感度が判定不能だった原因)

**注意点**: 床の最小 CV (µm 級) では float32 のメッシュ量で CV 閉性が 1.2e-5 まで上がる
(全 100 万 CV 中 48 個、すべて死水域の床)。床側の第一セルを開口側より粗くして (20 µm / 80 µm)
抑えている。閾値の根拠は `tools/check_mesh_extra.py` の docstring。

### 4.5 ソルバ設定 (現行レシピ §1+§2)

```yaml
mesh: {discretization: "node", nodeWallDirichlet: 1, nodeInletCornerWall: 1, meshFileName: mesh.h5, valueFileName: mesh.h5}
solver: "SLAU"
physProp: {thermalMethod: 0, viscMethod: 1, thermCondMethod: 1, prandtlLam: 0.72, cp: 1004.5, gamma: 1.4}   # CPG。TP は §4.10
space: {convMethod: 1, limiter: 2, pRef: 5474.72}
time: {unsteady: 0, timeIntegration: 11, nStepInner: 4,
       deltaT: {control: 1, cfl: C, cfl_pseudo: C, implicitRelax: 0.7, blockDPLUR: 1, detectNaN: 1}}
turbulence: {model: "sst", scalarDiffusion: 1, dilatationCorrection: 2, katoLaunder: 0,
             wallTreatmentSST: 0, turbulentPrandtl: 0.9, kInit: <前駆>, omegaInit: <前駆>}
output: {level: 1}
```

- **壁関数は使わない** (`wallTreatmentSST: 0` 低 Re 壁解像)。本形状はすきま底・開口リップ・
  内円柱上面の縁に**角線 (2 壁交線) が全周にある**。node 壁関数は角線ノードで代表点が取れず
  $u_\tau$=0 に退避する既知の未対応があり ([node-sst-wallfunction-utau-zero]、
  [`turbulence-node-wf-representative-point.md`](../accepted/turbulence-node-wf-representative-point.md) は
  3D 角線を対象外と明記)、壁解像なら node の $C_f$ は SU2 と一致する
  ([node-wallfunction-pk-convention-deficit])。
- `katoLaunder: 0` (剪断層主体では A/B 無影響)。`dilatationCorrection: 2` は現行既定。
- `sstEnergyIncludesK: 0` (既定)。全温は `VALUE/h0` から `total_quantities.py`。
- **`mesh.bndFirstOrder` は書かない** (使用禁止)。
- `implicitRelax: 0.7` は §1 レシピだが緩和は解をわずかに動かす (SERN 平面で 0.05 %,
  [sern-cfl-no-implicit-relax])。本段収束後に **`implicitRelax` 無し ($cfl\_pseudo\le6$) の確認 run** を取る。
- EOS 圧力床: $P_\infty$=5474.7 Pa で既定 `pMin`=1 Pa から 3 桁余裕。開口リップの局所膨張で床に
  当たったら `cfl_pseudo` を下げる ([implicit-cfl-ceiling-eos-floor])。

### 4.6 段階起動と、段間の移行判定

一様 IC + SST は前縁で NaN になる (case/48 実測) ため壁条件と物理を段階的に足す。
**同一メッシュ内の段間引き継ぎは index コピー** (`res` → `mesh.h5` の VALUE 直接コピー)。
**3D cross-mesh restart は使わない** (§3.2-1)。Stage A と Stage B はそれぞれ S0 から立てる。

| 段 | 変更点 | 次数 | cfl_pseudo | nStepInner | step (初期値) |
| --- | --- | --- | --- | --- | --- |
| S0 | 全壁 `slip` (bcond 差し替え)・層流 (`model: none`)。キャビティ内圧の平衡化 | 1 次 | 0.5 | 10 | 2000 |
| S1 | `plate`/`plate_in` + キャビティを no-slip 断熱に | 1 次 | 0.3 | 10 | 2000 |
| S2 | キャビティ 4 面を `wall_isothermal` Ts=500 に | 1 次 | 0.5 | 10 | 2000 |
| S3 | SST 投入 (`kInit`/`omegaInit`) | 1 次 | 0.3 | 10 | 3000 |
| S4 | SST mid | 1 次 | 1.0 | 10 | 3000 |
| S5 | 2 次ランプ 0.5 → 1 → 2 | 2 次 | 0.5/1/2 | 4 | 各 2000 |
| S6 | 本段 (+`implicitRelax: 0.7`) | 2 次 | 2 | 4 | 20000〜 |

- **IC (codex M5 採用)**: 自由流一様 + 壁近傍 tanh ランプ**だけでは不足**。すきま中央の初期速度が
  ほぼ M5 のままになり「深部は準静止」という想定と矛盾し大きな速度過渡を生む。よって
  **キャビティ内 ($z<0$) は $U$=0・$T$=$T_w$・$P$=$P_\infty$ から始め**、開口近傍 ($|z|<$ 1 mm) で
  外部場へ滑らかに繋ぐ。平板上は前駆 BL 分布を $z$ 方向に写す (入口 CSV と同じもの)。
- **段間の移行判定は step 数固定ではない** (固定値は初回の目安):
  全残差が単調に下降していること、$P,T,|U|$ が物理範囲 ($\rho>0$, $T>0$, $P_s\le P_t$) に収まること、
  **キャビティ内の質量とエネルギーの変化率が段内で減少に転じていること**を満たしてから次段へ。
  満たさない段は step を伸ばすか cfl を下げる。
- `detectNaN: 1` は常時。各段投入直後に `residual_history.csv` 先頭 50 行の非有限を確認。

### 4.7 収束・準定常の判定量

「NaN なし」「残差低下」だけでは足りない。**報告する量そのもの**の定常化を
`check_quasisteady.py --series-csv <csv> --series-cols ...` で判定する。CSV は
`tools/cavity_eval.py --series` が全 `res_*.h5` から作る。

**ゼロ近傍の列をそのまま渡さない** (codex M7 採用、実機確認): `check_quasisteady.py` は
`max(|mean|,1e-30)` で正規化するので、**±1e−12 のゼロ平均系列も 1e−10→1e−14 の減衰系列も
`DRIFTING`** になる。よって各列は**非ゼロの基準尺度で無次元化した検査列**を渡し、生値も併記する。

| 列 (検査に渡す) | 定義 | 基準尺度 / 絶対許容 |
| --- | --- | --- |
| `dT_floor` | 床から 1 mm 上・周方向平均のガス温度 − $T_w$ | 尺度 $T_{aw}-T_w$ (626 K)。絶対許容: 末尾窓の drift ≤ **2 K** |
| `dT_mid` | 深さ中央 ($z=-d/2$)・すきま中央・周方向平均 − $T_w$ | 同上 |
| `dT_mouth` | 開口下 2 mm・すきま中央・周方向平均 − $T_w$ | 同上 |
| `dT_up`, `dT_dn` | **流れ方向**上流側 / 下流側 方位 ($\theta$=180°/0°) の `dT_mid` | 同上。**偏心 (流れ方向非対称) の指標であって、半割が禁じる左右非対称の指標ではない** (codex M1) |
| `dT_left`, `dT_right` | $\theta$=±90° 方位の `dT_mid` | **全周 URANS でのみ意味を持つ** (半割では恒等) |
| `zpen` | $T-T_w > \epsilon_T$ を満たす最深点の深さ。$\epsilon_T$ = 10 / **25** / 50 K を併記 | 尺度 $d$ (50 mm)。絶対許容: drift ≤ **0.5 mm** |
| `mdot_in`, `mdot_out` | 開口面の $\rho u_z$ の負側 / 正側の**絶対**積分 | 尺度は自身の末尾平均 (非ゼロ) |
| `mdot_imbalance` | $\dot m_{net}/\max(\dot m_{in},\dot m_{out})$ | **交換流量で正規化した不釣合い** (生の `mdot_net` は定常で 0 へ行くので検査列にしない: codex M7)。合格は \|値\| ≤ **1e-3** |
| `q_outer`,`q_cylside`,`q_floor` | 各等温壁の面積分入熱 [W] (半割) | 尺度 $\sum q$ (非ゼロ)。絶対許容: 各 drift ≤ 総入熱の 1 % |
| `q_cyltop` | 内円柱上面の入熱 [W] | **CV 外**なので別枠 (§4.8) |

- **合否条件は §6.4 と §8 で同一**にする (codex M7: 旧版は §4.9「定常解との差 10 %」/
  §6.4「振幅 10 %」/ §8「全量 STEADY」がばらばらだった)。定常解の合否と周期統計の合否を分け、
  URANS では**平均差・振幅・窓間 drift の 3 つ**を判定する (§6.4 の表)。
- `T\ge500` K を無条件の棄却基準にはしない。物理妥当性は $\rho>0$, $T>0$, $P_s\le P_t$ と収支で見る。

### 4.7.1 熱伝達率の基準温度 (ユーザ指定 2026-09-19)

熱伝達率は **2 通り**出す。深いキャビティでは外部流の回復温度を駆動温度にすると、実際の駆動温度差
(深部で 37 K) に対して 626 K で割ることになり見かけ上ほぼ 0 に落ちるため、**局所の基準温度**を主とする。

| 記号 | 駆動温度 | 定義 |
| --- | --- | --- |
| **`h_ref` (主)** | **基準温度 $T_{0,ref}$** = **すきま中央面の最近傍点の総温**。**中央面の底は「すきま幅の半分」で切り取る** (床の壁点に対する最近傍点が床から gap/2 上になり、側壁の壁点が中央面まで gap/2 なのと整合する) | $h_{ref}=q''/(T_{0,ref}-T_w)$ |
| `h_aw` (従) | 外部流の回復温度 $T_{aw}$ (一定) | $h_{aw}=q''/(T_{aw}-T_w)$ |

- 総温は **`VALUE/h0` から `tools/total_quantities.py` の `total_state` で逆算**する
  (スクリプトで $T+u^2/2c_p$ を組まない: AGENTS.md「出力と後処理の原則」。CPG/TP 両対応)。
- 実測 (Stage B): `h_ref` は深さによらず **70–100 W/m²K** でほぼ一定 (開口部 230–265, 床 41)、
  一方 `h_aw` は深部で 7.8 W/m²K まで落ちる。設計に使うのは `h_ref`。

### 4.8 検査体積とエネルギー収支

**CV の定義は偏心に対応させる** (codex M5 採用): キャビティ検査体積は
$$\{(x,y,z)\;:\;x^2+y^2<R_o^2,\;(x-x_{off})^2+y^2>R_i^2,\;-d<z<0\}$$
であって $R_i<r<R_o$ **ではない**。既定半径・$x_{off}$=1 mm で同心マスクを使うと
真の半環断面の **12.06 % (22.498 mm²) を取りこぼし**、同面積の固体側を誤って含む
(面積の総和が合っても誤りは見えない)。CV マスク・測線・周方向平均の重みは
**manifest から解決した共通関数**を全ツールが使う (§4.11)。

CV の境界は **`cav_outer` + `cyl_side` + `cav_floor` + 開口面 + 対称面 2 枚**。
**`cyl_top` は CV を囲まない外部流の面**なので収支に入れない。

**物理の壁熱流束と、離散の保存検査を別々に実装する** (codex M6 採用、実機確認):

- **ソルバ出力 $q_w$ (一次情報, 2026-09-19 実装)**: 低 Re 壁 (`wallTreatment==0`) では従来 `qwall_b` /
  `utau_b` に誰も書かず壁面ダンプが全点 0 だったため、**[`viscousFlux_d.cu`](../../solver_density_cuda/cuda_forge/viscousFlux_d.cu)
  の壁カーネルで、残差に入れたのと同じ解像熱流束を同じ符号規約 (壁→流体正) で `qwall_b` に、
  $u_\tau$ を `utau_b` に診断出力する**ようにした (mode 1/2 と `sstEnergyWf != 0` では入力側なので
  触れない = 解はビット不変)。これで `res_<群>_<physID>_<step>.h5` から $q_w$・$\tau_w$・$u_\tau$・$y^+$ が
  そのまま取れる。検算: 断熱壁 `plate` は**厳密に 0**、等温壁は負 (= 壁に熱が入る)。
- **物理 $q_w$ (照合用)**: 場から $q_w=\lambda_w\,\partial T/\partial n$ を片側 2 次差分で組む。
  **壁近傍を除外したマスクで評価してはいけない** (2026-09-19 実測でソルバ値の 7 倍になった)。
  差分間隔は VL 第一層に見合わせる。これは**物理量の推定**であって離散保存量ではない。
- **離散保存**: forge は node の壁 Dirichlet で**壁ノードのエネルギー残差を 0 化**し
  ([`nodeWallDirichlet_d.cu`](../../solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu) の
  `zero_res_roe_bplane_d`)、壁流束も後処理の 3 点差分とは別の勾配を使う
  ([`viscousFlux_d.cu`](../../solver_density_cuda/cuda_forge/viscousFlux_d.cu))。よって
  **勾配積分と離散的なエネルギー授受が有限格子で一致する保証はない**。
  離散収支は**非拘束 CV 群 + 壁隣接界面の流束**で別途検算し、
  **勾配入熱との差が細分で減ること**を確認する (減らなければ評価法の欠陥)。
- **検証**: 平板区間に加え、**角と曲面を持つ伝導問題** (解析解のある同心環状すきま) で
  法線・面積重み・共有ノード (等温と断熱が隣接するリップ) の扱いを確かめてから本計算に使う。
- 開口の全エンタルピー流束は**保存量 `VALUE/h0` と採用 EOS の datum** を使う
  (スクリプトで $T+u^2/2c_p$ を組まない: AGENTS.md)。

**収支残差の出どころは離散スキームでなく後処理の面積分** (2026-09-19 実測, ユーザ指摘
「テトラだと各面での収支が取れていない説」への切り分け):

- tet+prism (stageD, 740k 節点) でも **CV ごとの閉性 $\|\Sigma\mathbf{S}\|/\Sigma|\mathbf{S}|$ は
  max 5.11e-6 / p99.9 2.42e-6 で PASS** (`tools/check_mesh_extra.py`)。双対 CV が閉じているので
  面ごとの収支が崩れているのではない。全ヘキサは max 1.21e-5 で同程度。
- 崩れているのは**後処理で張った開口面の面積分**で、これは**質量で直接測れる**:
  底の閉じたキャビティなので定常なら開口を通る**正味**質量流量は厳密に 0。実測の
  正味/片道は tet 3.3 % (z=-0.5 mm) 〜 7.3 % (z=-1 mm) に対し**ヘキサ 0.73 %**。
- 正味エンタルピー流入は「熱い流入」と「冷えた流出」という**大きな 2 項の差**なので、
  同じ絶対積分誤差が相対では一桁増幅され、収支残差 tet 22 % / **ヘキサ 2.3 %** になる。
  ヘキサでは評価面が押し出しの節点層と一致して面積重みがそのまま使えるのが効いている。
- したがって**収支残差はメッシュ型の保存性の指標ではなく、評価面上の積分精度の指標**として読む。
  報告には必ず**質量の正味/片道**を併記し、これが数 % を超える評価面では収支残差を根拠にしない。
  `tools/cavity_eval.py --flux-depth <mm>` で評価面深さを振って感度を見る。

### 4.9 定常・半割 RANS の適用限界と、非定常確認

- **cavity 分類**: 「深さ/すきま幅 20」は分類指標ではない。cavity flow の分類に使う
  **流れ方向開口長 / 深さ**は対称面の各スリットで $2.5/50=0.05$ で、極端に深い
  (narrow deep) キャビティである。**せん断層が開口を跨ぐのが open cavity、床へ再付着するのが
  closed cavity** (旧版は逆に書いていた)。
- **定常擬似時間の収束は物理的振動の不在を証明しない**。半割の `slip` は法線流速を除去するので
  ([`boundaryCond_d.cu`](../../solver_density_cuda/cuda_forge/boundaryCond_d.cu) の鏡像)、
  **半割 URANS では、半割が禁じたモードを検出できない** (codex M1)。
- よって定常・半割 RANS は**条件付きの初期評価**とし、次を**無条件の必須ゲート**にする:

  **(a) Stage A 全周 (360°) URANS** — 半割 URANS では代用しない。
  - 完全対称 IC を無擾乱で回すだけでは検査にならないので、**小さな左右非対称擾乱**
    (開口近傍の $k$ に振幅 1 % の $\sin\theta$ 分布) を初期に与える。
  - 測るのは左右の `dT_left`/`dT_right`・開口交換流量の左右差、すなわち**反対称成分の成長/減衰**。
  - `dT_up`/`dT_dn` は流れ方向の指標なのでこの検査には使えない。

  **(b) URANS の設定** (codex M2 採用、実機確認):
  ```yaml
  time: {unsteady: 1, dualTime: 1, timeIntegration: 11,
         last: {nStepOuter: N}, nSubIterDualTime: 10,
         deltaT: {control: 0, dt: <固定>, cfl_pseudo: 12, blockDPLUR: 1, detectNaN: 1}}   # implicitRelax は書かない
  ```
  - **`deltaT.control: 0` (固定物理 dt) が必須**。定常設定の `control: 1` を継承すると
    `main.cpp` が例外終了する (`Dual-time implicit requires time.deltaT.control=0`)。
    `dualTime: 1` と `blockDPLUR: 1` も必須。
  - 物理 dt は物理 CFL ≲ 12 (Stage A の最小セル 40 µm → dt ≈ 3e-7 s) を**初期値**とし、
    **dt 半減・`nSub` 倍増・開口せん断層の細分**をそれぞれ 1 本ずつ回して、
    **平均・振幅・卓越周波数の変化が結論量の 2 % 以下**であることを確認する
    (`nSub` 10 の根拠は別ケース限定で、現行手順も nSub 倍増比較を要求している)。
  - 観測窓 2.5 ms (キャビティ音響往復 223 µs の ~11 倍) は**最初の確認時刻**にすぎない。
    音響時間は深部の交換・熱緩和時間を保証しないので、**窓を倍にして深部温度・入熱・交換量の
    統計 (平均・振幅) が安定するまで延長**する。
  - 粗い Stage A で振動が減衰しても**空間・時間離散化による減衰**を排除できないので、
    細分側の比較 (上記) を必ず付ける。
  - 判定: (i) 反対称成分が減衰し、かつ結論量の振幅が §6.4 の許容内 → 定常・半割 RANS を正本にする。
    (ii) リミットサイクル → **平均±振幅**で報告し、定常解との差を明記する。
    (iii) 反対称成分が成長 → 半割仮定は棄却し、生産計算を全周に変える (別計画に上げる)。

### 4.10 生産 EOS — **semi-perfect (TP) を既定とする**

§3.1 の通り CPG は $T_{aw}$ を +61 K (+5.4 %) 過大にする。成果物が温度そのものなので、
**生産 EOS は semi-perfect (TP) を既定**とする (codex M8 の趣旨を採用して旧版の
「Stage A で 3 % 判定して CPG を選ぶ」を撤回した。数値側に格子 5 % / 収支 5 % / drift 2 K の
許容がある中で「3 % 未満だから CPG で良い」は示せない)。

- TP 設定: `physProp: {thermalMethod: 2, species: [MIXDRY], speciesDBFile: species_db.yaml,
  thermoHrefTemp: 298.15}`。**`MIXDRY` は内蔵種でないので DB が必須**
  ([`speciesDB.cpp`](../../solver_density_cuda/input/speciesDB.cpp) は未定義種を拒否する)。
  DB は `design/forge_design/gas/semiperfect.mixture_pseudo_species` で乾燥空気 1 擬似種として生成し、
  `conditions.json` の `dry_air_Y` を単一ソースにする。**前駆 2D・入口 CSV も TP で作る**
  (EOS ごとに整合させる)。
- TP の `cfl_pseudo` は 0.5〜2 から上げる。`implicitRelax: 0.7` を付ければ 6〜8 まで可 (§3 レシピ)。
- **CPG は「パイプライン疎通と段階起動レシピの確立」にだけ使う** (安く速い)。CPG 場から TP へ
  移すときは **`roe` をコピーしない** ($P,T,\mathbf{U},Y$ を保持して TP の EOS で $\rho$・
  内部エネルギー・保存量を組み直す。datum が違う)。
- **CPG/TP の差は Stage A で測って記録する** (採否の判定ではなく、既往結果との比較用の情報):
  `dT_mid`・`dT_mouth`・`zpen`・**深さ方向温度プロファイル**・**壁別入熱**を比較する。
  もし TP が Stage A で安定に回らなければ CPG に落とし、**偏り (+61 K 級) を明記**して報告する。

### 4.11 単一入力 (manifest) — 偏心・すきま・マッハ数の可変前提に対応

ユーザ要件 (偏心・すきま量・マッハ数が変わる) と codex M5 を受け、**最上位入力を 1 つ**にする。

```
case.json  (= conditions.json + geometry を束ねた最上位入力)
   └─ setup.py --resolve --> manifest.json  (解決済みの全派生量)
          ├─ 自由流 (CPG/TP), 入口 BC 値, 前駆の目標 delta
          ├─ 形状 (Ro, Ri, x_off, depth, 領域) と派生: gap_min/max/nom, グループ別面積
          ├─ メッシュ: VL 第一層/stretch/層数/総厚 (= f(gap_min)), 面別サイズ
          └─ 評価: CV マスク式, 測線 (偏心円との交点), 周方向平均の重み, 測点深さ, 基準尺度
```

- CAD・メッシュ・IC・BC・評価の**全段が manifest だけを読む**。層厚などの派生量を
  文書や別スクリプトに転記しない (codex m10 の転記誤り対策)。
- **共通関数**を `tools/geom_common.py` に置く: CV マスク (§4.8 の偏心対応式)、すきま中央線
  (偏心 2 円の中点軌跡)、周方向平均の重み、$\delta/g_{nom}$ の基準幅。
- manifest には**適用範囲と再検証条件**も書く: M・すきま・偏心を変えたら、
  (i) 前駆 2D をやり直す、(ii) VL 層数を $g_{min}$ から引き直す、(iii) 領域感度を取り直す。

## 5. 実装ステップ

1. ✅ `case/49.plate_annular_cavity_m5/` と `cad/build_geom.py` (パラメトリック形状 + グループ面積検証。
   偏心・すきま変更・助走の 3 通りで通過)。
2. ✅ `conditions.json` + `setup.py` (自由流・CPG/TP 両方の $T_t,T_{aw}$ を導出) →
   **`case.json` + `manifest.json` に統合する** (§4.11)。
3. ✅ `precursor/` の 2D メッシュと段階起動 (`gen_mesh.py`/`gen_runs.py`)。TP 版を追加する。
4. `tools/geom_common.py` (manifest から CV マスク・測線・重み) と `tools/extract_profile.py`
   (前駆 → 入口 CSV を**直書き**、$\delta_{99},\delta^*,\theta,Re_\theta$ と時系列)。
5. `tools/cavity_eval.py` (§4.7 の列・§4.8 の物理 $q_w$ と離散収支の**別実装**・時系列 CSV) と
   `tools/check_mesh_extra.py` (§4.4 の表)。**Stage A より前に用意し、解析解のある伝導問題で検算**。
6. `cad/mesh_salome.py` + `cad/med_to_msh41.py` → **塞ぎメッシュ** (キャビティ無し) と Stage A メッシュ。
7. **塞ぎ 3D で入口移送と BL 発達を検証** (§4.3。ここが通らなければ先に進まない)。
8. `gen_runs.py`: S0–S6 (config 生成・bcond 差し替え・index コピー・移行判定・IC)。
9. Stage A: 起動レシピ通し、GPU/host メモリと s/step 実測、$y_1^+$ 実測、CPG↔TP 差の記録、
   **全周 URANS** (dt 半減・nSub 倍増・細分の 3 本込み)、評価ツールの検算。
10. Stage B メッシュ → 品質ゲート → S0 から段階起動 → S6 を準定常まで延長。
11. 感度: 格子 (両方向)・領域・流入 BL $\delta$・`cyl_top` 断熱・`implicitRelax` 無し。
12. 結果まとめ (図・表・README run 一覧)、codex result レビュー。

### 5.1 残作業 (優先順)

**進捗 (2026-09-19)**: ✅ = 完了。CPG で Stage A/B/C/D の 4 格子まで通し、熱流束・熱伝達率・
すきま内部の流れまで出せる状態。

| # | 項目 | 状態 / 内容 |
| --- | --- | --- |
| 1 | manifest 統合 | ✅ `case.json` → `setup.py --resolve` → `manifest.json`。CAD/メッシュ/IC/BC/評価が全部これを読む。`CASE49_MANIFEST` で差し替え可 (塞ぎ形状用) |
| 2 | `tools/geom_common.py` | ✅ 偏心対応 CV マスク (同心式は 12 % 取りこぼすことを MC で再現)・すきま中央線・周方向重み。**方位は上流基準 θ=0** (2026-09-19 ユーザ指定) |
| 3 | 入口 CSV 直書き | ✅ `tools/extract_profile.py`。起動ログで `applied: ro Ux Uy Uz Ps k omega` の 7 列反映を確認 |
| 4 | 評価ツール | ✅ `tools/cavity_eval.py` (温度・侵入深さ・開口流束・**ソルバ出力 q_w**・**h_ref**)、`tools/plot_cavity_fields.py` (流れ・熱流束・熱伝達率・底面)、`tools/compare_runs.py` (格子比較)。⏳ `check_mesh_extra.py` (体積/sliver/dual) は未着手 |
| 5 | 塞ぎ 3D の BL 検証 | ⏳ 形状・メッシュ (410k 節点, PASS) と変換まで完了。**run と前駆比較が未実施** |
| 6 | Salome メッシュ | ✅ Stage A 69k / B 330k / C 589k / D 741k。品質は全て PASS。VL/サイズの整合ガード実装済 |
| 7 | 段階起動 | ✅ S0–S6 が 4 格子とも NaN 0 で完走。6.2 / 21 / 27 / 64 ms/step |
| 8 | TP 化 | ✅ 乾燥空気 1 擬似種 DB 生成・`--gas TP`・IC の datum 組み直し・R_tp での自由流。⏳ **TP の本 run 未実施** |
| 9 | **全周 URANS** | ⏳ 未着手 (plan §4.9 の必須ゲート) |
| 10 | 予算確定 | ✅ 実測: 589k 節点で GPU 2.27 GB / 12 GB。**100 万節点級まではローカルで回せる** |
| 11 | 感度 run 群 | ⚠ **格子感度は判定不能だった**: Stage A–D は y₁・層数・面サイズ・成長率を同時に変えており**系統列でない** (tet+prism では VL 総厚と接線サイズが結合し独立に振れない)。C↔D で総入熱が 20.6 % 違う。→ **全ヘキサ (§4.4.1) の `--scale` で系統列を作り直す**。`h_ref` は C↔D で 0.16 % 一致しており格子に頑健。領域・流入 BL δ・`cyl_top` 断熱・`implicitRelax` 無しは未着手 |
| 12 | BL 解像が不足なら | (退避策。まだ不要) |
| 13 | 3D 最近傍 restart | (必要になっていない。各 Stage は S0 から立てている) |
| 14 | **偏心・すきま変更 run** | ⏳ `sweep_offset.py` を用意 (CAD→メッシュ→計算→評価を偏心量ごとに通す。VL 層数は最小すきまから自動)。格子収束の確認後に投入 |

## 6. 検証

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| plan | `2026-09-19` | [`notes/reviews/2026-09-19-case-plate-annular-cavity-m5-plan.md`](../../notes/reviews/2026-09-19-case-plate-annular-cavity-m5-plan.md) | **NO-GO**, C2/M8/m1 | **全件採用** (根拠は自分で再検証済み): C1 (interp_field 2D) → §3.2-1 + 3D restart を使わない設計 + #12 / C2 (CV 誤り・qwall 0) → §4.8 + #2 / M3 (cavity 分類と定常 RANS の正当化) → §4.9 + #7,#11 / M4 (重心分類の不成立) → §4.1 (実装済み・検証済み) / M5 (変換時の壁タグ・pRef・IC) → §4.2,§4.6 / M6 (sliver・$y_1^+$・予算) → §4.4 + #2,#8 / M7 (仮定の明示・領域感度) → §3.1,§4.2 + #9 / M8 (TP cp +18.27 %) → §3.1,§4.10 + #6 / M9 (準定常の列定義) → §4.7 / M10 (case/48 の引用限定・中止条件) → §3.3,§6.4 / m11 (参照・Sutherland 定数) → §3.1 |
| plan-2 | `2026-09-19` | [`notes/reviews/2026-09-19-case-plate-annular-cavity-m5-plan-2.md`](../../notes/reviews/2026-09-19-case-plate-annular-cavity-m5-plan-2.md) | **GO-with-changes**, C0/M9/m1 | **全件採用** (根拠は自分でコード確認済み): M1 (半割 URANS では半割仮定を検証できない — `slip` が法線流速を除去) → §4.9(a) 全周 URANS を無条件ゲート化 + §5.1 #9 / M2 (`deltaT.control:0` 必須・nSub 倍増・窓延長) → §4.9(b) / M3 (`gen --table` が `Ps` を落とす) → §4.3 + #3 / M4 (BL 比較に下流発達を織り込む・塞ぎメッシュ検証・19 点は先験判断しない) → §4.3 + #5 / M5 (偏心 CV マスク 12 % 誤り・manifest) → §4.8, §4.11 + #1,#2 / M6 (物理 $q_w$ と離散収支の別実装) → §4.8 + #4 / M7 (ゼロ近傍の準定常判定・合否条件の不整合) → §4.7, §6.4 / M8 (EOS 3 % 判定は不可) → §4.10 で **TP を生産既定に変更** / M9 (品質ゲートの定量化) → §4.4 の表 / m10 (参照・VL 数値 0.869495/0.181899) → メタ・§4.4 |

### 6.2 メッシュ

**§4.4 の品質ゲート表**をそのまま合否条件にする (`check_mesh_quality.py` + `tools/check_mesh_extra.py`)。
層あり率は**壁 6 グループ** (`cyl_top` を含む) で 100 %。節点数は Stage A 実測から定めた上限以内。

### 6.3 計算

- **NaN**: 各段投入直後と最終で `residual_history.csv` 全列と `res_*.h5` の `VALUE/*` の非有限を確認。
- **収束**: `check_convergence.py` の VERDICT を貼る。**全列** (`rms_ro`,`rms_roUx/y/z`,`rms_roe`,
  `rms_roK`,`rms_roOmega`)。`rms_ro` 単独で判断しない。**PASS と準定常は別ゲート**として維持する。
- **準定常**: §4.7 の CSV 列すべてが `STEADY` (+ $\Delta T$ 列は絶対 K の drift 許容も満たす)。
  `DRIFTING`/`TRANSIENT-UNSETTLED` を「定常」と書かない。`OSCILLATING` は平均±振幅で報告。

### 6.4 物理の妥当性ゲートと、中止・診断条件

参照解 (実験・文献) は本形状に存在しない。よって**内部整合ゲート**で妥当性を主張する:

| ゲート | 基準 |
| --- | --- |
| キャビティ CV 収支 (§4.8) | $\sum_{3壁}q_w + \Phi_{開口}$ が総入熱の 5 % 以内で閉じる |
| 断熱平板の回復温度 | $T_w$ が $T_{aw}$ (採用 EOS の値) の ±3 %。**評価区間はキャビティ・入口・出口の影響を避けた $x\in[-60,-40]$ mm に限定** (codex M7) |
| 流入 BL の移送・発達 (§4.3) | **塞ぎ 3D** の $x=-30$ mm が **前駆の $x_0+50$ mm** と: $\delta_{99},\delta^*,\theta$ 5 % 以内、$U,T,k,\omega$ の法線分布と $\tau_w$ が 5 % 以内、かつ外層・prism→tet 遷移の細分で差が縮む |
| $y_1^+$ | 第一内部ノードで平板・キャビティ壁とも ≤ 2 (case/48 の $q_w$ 許容 −3.7 % 以内) |
| 格子感度 | $y_1$・`plate_in`・すきま接線を**細化側も含めて**振り、`dT_mid`/`zpen`/総入熱の変化 ≤ 5 % |
| 領域感度 | top/side/x_out を広げて結論量の変化 ≤ 2 % |
| 非定常性 (§4.9) | **全周 URANS** で (i) 反対称成分が減衰、(ii) 結論量 (`dT_mid`,`dT_mouth`,`zpen`,総入熱) の**振幅が平均の 10 % 未満**、(iii) dt 半減 / nSub 倍増 / 開口細分での**平均・振幅・卓越周波数の変化が 2 % 以下**、(iv) 観測窓を倍にして統計が変わらない。この 4 つが揃って初めて定常・半割 RANS を正本にする |
| 準定常 (§4.7) | 検査列 (尺度化済み) が全て `STEADY` **かつ** 絶対許容 ($\Delta T$ drift ≤ 2 K, `zpen` ≤ 0.5 mm, 壁入熱 drift ≤ 総入熱の 1 %, \|`mdot_imbalance`\| ≤ 1e-3) を満たす。**§4.7・§6.3・§8 はこの 1 つの条件を指す** |

**中止・診断条件** (codex M10): S6 を 20000 step 延長しても (a) 残差が下がらず、かつ
(b) `dT_*`/`zpen` の drift が 10k step あたり 5 % を超え続ける場合は、step を足すのをやめて
診断に切り替える — 出力間隔を 5 step にして最初の非有限/床当たりの位置を特定
([detectnan-reports-symptom-not-cause])、EOS 床・開口リップの膨張・角線ノードを確認、
`cfl_pseudo` と `nStepInner` を振る。**「反復を増やす」以外の手を打つ**。

## 7. 影響範囲

- **新規のみ**: `case/49.plate_annular_cavity_m5/` (cad / precursor / tools / gen_runs.py / README)。
- **ソルバコード・共有ツールの変更は無い**。`interp_field.py` の 2D 制約は本計画では回避し
  (§3.2-1)、共有ツールには手を入れない (必要になれば別 plan)。`check_mesh_quality.py` に
  sliver 判定が無い件も case 側ツールで補い、昇格は別 plan とする。
- `plans/README.md` の active 一覧に追記 (済)。

## 8. 完了条件

- [ ] 形状・メッシュが §6.2 の全ゲート PASS (偏心・すきま変更でもスクリプトが通ることを確認済)
- [ ] 段階起動が全段 NaN 0 で完走し、S6 が `check_convergence.py` の VERDICT を満たす
- [ ] §6.4 の「準定常」ゲート (尺度化列が `STEADY` + 絶対許容) を満たす
- [ ] §6.4 の妥当性ゲート (収支・回復温度・BL 維持・$y_1^+$・格子/領域/BL 感度・非定常性) を満たす
- [ ] 生産 EOS = semi-perfect (TP) で回し切る (§4.10)。CPG に落とした場合は偏りを明記
- [ ] **全周 URANS の 4 条件** (§6.4) をクリアし、定常・半割 RANS を正本にできると示す
- [ ] 温度場の図・深さプロファイル・等温壁 3 面 + `cyl_top` の熱流束表を README に記載
- [ ] codex レビューを §6.1 に記録し Critical/Major の採否を §5.1 に反映 (plan 段は 2 巡)
- [ ] `status: done` + §9 変更ログ + `plans/accepted/` へ移動 + `plans/README.md` 同期

## 9. 変更ログ

- `2026-09-19` — 初稿 (20 km M5・境界 11 グループ・VL 13 層・段階起動 S0–S6)。
- `2026-09-19` — **codex plan レビュー #1 (NO-GO, C2/M8/m1) を全件採用して全面改稿**。主な変更:
  ① 3D cross-mesh restart を廃止 (`interp_field.py` が 2 次元最近傍であることを実機確認)、
  ② キャビティ検査体積を 3 壁 + 開口に正し `cyl_top` を別枠に、$q_w$ は場の勾配から
  (壁ダンプの `qwall`/`utau`/`ypls` が低 Re で全点 0 であることを実機確認)、
  ③ cavity 分類の記述を訂正し定常・半割 RANS を条件付き初期評価に格下げ、URANS 非定常性チェックを
  必須ゲート化、④ 面分類を曲面種別・軸・半径 + 頂点最大半径に変更しグループ別面積照合を追加
  (実装・3 形状で検証済)、⑤ 変換は最終壁タグ・S0 は bcond 差し替え・`pRef` 追加・キャビティ静止 IC、
  ⑥ sliver/体積ゲートと $y_1^+$ の第一内部ノード評価、⑦ 仮定 (高度・全乱流・`cyl_top`) の明示と
  領域・BL 感度の必須化、⑧ TP の $c_p$ 差を再計算 (+18.27 %, $T_{aw}$ 1126 K vs CPG 1188 K) し
  生産 EOS の決定を Stage A のゲートに、⑨ 準定常の列を $T-T_w$・流入/流出別・侵入深さ定義に、
  ⑩ case/48 の引用を量ごとに限定し中止・診断条件を追加、⑪ 参照と Sutherland 定数を実装に合わせた。
- `2026-09-19` — **ユーザ追加要件**: 内円柱の下流偏心 (すきま不均一)・すきま量変更・**マッハ数変更**が
  ありうるため、形状/条件/メッシュ/BC/評価を単一ソースのパラメータ駆動にする方針を §1・§3.1・§4.1 に反映。
  形状スクリプトは偏心・すきま変更・助走有無の 3 通りで検証済み。
- `2026-09-19` — **codex plan レビュー #2 (GO-with-changes, C0/M9/m1) を全件採用**。主な変更:
  ① 半割 URANS では半割仮定を検証できない (`slip` が法線流速を除去) ため **Stage A 全周 360° URANS +
  左右非対称擾乱を無条件の必須ゲート**に、② URANS は `unsteady:1`+`dualTime:1`+**`deltaT.control:0`**
  (定常の `control:1` を継承すると例外終了することをコードで確認) とし dt 半減・nSub 倍増・細分で
  2 % 以内・観測窓は統計が安定するまで倍増、③ `gen_inlet_profile.py gen --table` は **`Ps` を落とす**ので
  入口 CSV は case 側で直書き (`Uz` = 前駆の壁法線速度も追加)、④ BL 比較は**下流発達を織り込み**
  前駆の $x_0+50$ mm と比べる + **キャビティを塞いだ 3D** で移送を先に検証、19 点は先験判断しない、
  ⑤ 偏心時の CV マスクは同心式では 12 % 誤るので**最上位 manifest と共通関数**に集約 (§4.11 新設)、
  ⑥ 物理 $q_w$ (勾配) と**離散保存検査を別実装**し細分で差が縮むことを確認 (壁ノードのエネルギー残差が
  0 化されることをコードで確認)、⑦ ゼロ近傍で `check_quasisteady` が破綻するので**尺度化列 + 絶対許容**に、
  合否条件を §6.4 の 1 つへ統一、⑧ **生産 EOS を semi-perfect (TP) に確定** (数値許容 5 % の中で
  「差 3 % 未満だから CPG」は示せない)、⑨ メッシュ品質ゲートを要素種別ごとの定量表に、
  ⑩ 参照と VL 数値 (総厚 0.869495 mm / 最終層 0.181899 mm) を修正。
- `2026-09-19` — **設計変更 (自発)**: 流入 BL を 3D 領域内で発達させる案を廃し、**2D 前駆計算 +
  `inletProfile`** に変更して領域を $x\in[-250,200]$ → $[-80,120]$ に縮小 (§4.3)。理由: すきま 2.5 mm が
  prism 総厚を ≤1 mm に縛る一方で外部 BL は 5 mm あり、Salome の単一 ViscousLayers では両立しない。
  縮小で平板面積が 3 分の 1 になり、かつ「流入 BL が既知」なので 3D の平板解像度を前駆解との
  $\delta^*,\theta$ 比較で定量検査できる。
