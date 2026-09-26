# 軸対称の固体伝導 (`fem2d`) — ノズル壁の共役熱伝達

## メタ

- **area**: `boundary`
- **status**: `draft`
- **related_docs**:
  - [`methods/boundary.md`](../../methods/boundary.md) 「共役熱伝達 (CHT)」— 軸対称は現状拒否 (2026-09-27) と、本計画の設計 (「軸対称の `fem2d` (計画中)」)
  - [`methods/architecture/overview.md`](../../methods/architecture/overview.md) — 軸対称の $r$ 重み幾何 (`axisymMethod: 0`)
- **related_plans**:
  - [`../accepted/boundary-conjugate-heat-transfer.md`](../accepted/boundary-conjugate-heat-transfer.md) — **親** (平面 2D の `fem2d`、界面契約・G-if・G-cons)。§5.1 #95 を本計画に切り出した
  - [`tooling-nozzle-isothermal-wall-chain.md`](tooling-nozzle-isothermal-wall-chain.md) — ノズル等温壁チェーン (壁温を仮定している側)
- **created**: `2026-09-27`
- **owner**: `sano`

## 1. 目的

`mode: fem2d` を**軸対称 (`mesh.isAxisymmetric: 1`, `axisymMethod: 0`)** で使えるようにし、
**ノズル壁の共役熱伝達** (壁厚 + 外表面に熱伝達率 $h$ と $T_c$ を貼る = 固体の `hole` Robin) を解けるようにする。
いま `case/40` 等は壁温を**仮定**して生産値を出しているので (壁温 1400 ± 15 K を仮定)、それを**解**にするのが目的。
完了時: 円筒殻の解析解 (対数抵抗) に一致する軸対称 CHT がソルバ内で回り、平面 2D の結果は変わらない。

## 2. スコープ

- **やる**:
  - 固体の弱形式に $r$ 重みを入れる (要素剛性・Robin 辺・界面の集中量)。C++ (`conjugate/solidFem2d.cpp`) と参照オラクル (`tools/solid_fem2d.py`) の両方
  - 軸対称の起動時検査 (`axisymMethod: 0` のみ受理、固体節点 $r\ge0$、界面の集中量 $A_i^r>0$)
  - `iface_q_eff` (熱流束の診断・`local1d` の入力) の面積を軸対称で $r$ 重みにする — **現状は平面面積で割っており $r$ 倍ずれる** (§3)
  - 検証: 固体単体の円筒殻解析解、ソルバ内連成の同心円環 (case/52 の軸対称版)、平面 2D の回帰
- **やらない**:
  - `axisymMethod: 1` (SU2 流 planar + $1/y$ ソース) — 流体の界面荷重の単位が別物になるので**拒否のまま**。必要になったら別 plan
  - `local1d` の軸対称 (薄肉の対数抵抗化) — **拒否のまま**。ノズル壁は `fem2d` で扱う
  - ノズル生産ケース (`case/40` 等) への適用 — 検証後に別 plan (壁厚・外表面 $h$ の出典決めが要る)
  - 3D (`fem2d` の 3D、親 plan §5.1 #96)、放射、dual-time

## 3. 関連 docs と前提 (2026-09-27 に実測・コードで確認)

- **軸対称は黙って間違っていた → 拒否を入れた (本計画 §5.1 #1、済)**。軸対称メッシュも $z\equiv0$ なので
  `initSolidFem2d` の平面ガード (`conjugateWall.cpp`) を素通りする。`axisymMethod: 0` では
  `variables.cpp`:530 が面ベクトル・面積・体積を $r$ 倍するので、流体の残差 `ifaceRraw` と積分済み荷重 `iface_Qf_eff` は
  **[W/rad]** ($2\pi$ を掛けない)。一方
  - `fem2d` の固体は平面 [W/m] の集中辺長で組む → **荷重が節点ごとに $r$ 倍ずれる**
  - `iface_q_eff` = `iface_Qf_eff` / `msh.planes[ip].surfArea` の分母は**ホストの平面面積** (`variables.cpp`:506 で $r$ を掛けるのはデバイス配列だけ) → **$q$ が $r$ 倍ずれる**。`local1d` はこれを使う
  → `initConjugateWalls` で**全モードの軸対称を拒否**した。ローカル 1 step で拒否を確認 (case/58 入力に `isAxisymmetric: 1`、rc=1)、平面はガードを通過。
- 親 plan の記述「軸対称の $r$ 重みは host・double で行う」(`methods/boundary.md`) は**未実装だった** (訂正済み)。
- 流体の $r$ は `y` 座標 (軸は $x$ 軸、$r=y\ge0$)。`axisRFloor` 以下の面は $r$ を床に上げる。
- 流体は半割面の**面重心の $r$** を掛ける。直線面で $r$ が線形なら**これは厳密な積分**で、固体の FEM 集中量 $\int N_i r\,ds$ との差は**積分領域と重み関数の違い**による (「一次近似」ではない。2026-09-27 codex diagnose)。例 $L=1, r_a=1, r_b=2$: 流体の半辺積分 (0.625, 0.875)、固体の集中量 (2/3, 5/6) — **総和は一致するが、一様な物理熱流束でも $Q_f/A_i^r$ は両端で 0.9375$q$ / 1.05$q$**。**保存性は保てても局所整合性は別に検証が要る** (§6 V-ax2b)。これを理由に荷重を面積で掛け直すことはしない。

## 4. 設計方針

### 4.1 単位 — ラジアンあたりに揃える

流体の $r$ 重み幾何は $2\pi$ を掛けない**ラジアンあたり**なので、固体も [W/rad]・[m²/rad] で組む。
界面荷重 `iface_Qf_eff` は**そのまま**渡す (親 plan の「積分済み荷重を面積で割って掛け直さない」を維持)。

### 4.2 固体の弱形式 ($r=y$)

$$\int_\Omega k\,\nabla T\cdot\nabla v\,r\,dA \;+\; \int_{\Gamma_R} h\,(T-T_c)\,v\,r\,ds \;=\; \sum_i Q_{f,i}\,v_i$$

線形三角形・要素内 $k_e$ = 節点 $k_s(T)$ の平均 (現行規約) として

- **剛性**: 勾配が要素内一定なので $K^e_{ij}=k_e\,\bar r_e\,A_e\,\nabla N_i\cdot\nabla N_j$ ($\bar r_e$ = 要素重心の $r$)。**厳密**。
- **Robin 辺** ($r$ は辺上で線形): $M^e=\frac{hL}{12}\begin{pmatrix}3r_a+r_b & r_a+r_b\\ r_a+r_b & r_a+3r_b\end{pmatrix}$、
  荷重 $b^e=\frac{hT_cL}{6}\,(2r_a+r_b,\;r_a+2r_b)$。**厳密**。
- **界面の集中量**: $A_i^r=\sum_{e\ni i}\frac{L_e}{6}(2r_i+r_j)$ [m²/rad] (= $\int N_i\,r\,ds$)。
  使い道は親 plan と同じく**熱流束への換算と $D_f=g_f A_i^r$ だけ** (荷重には使わない)。
- **最終積分に $r$ を掛けるだけでは足りない** (親 plan §4.3 の codex M6)。要素剛性と Robin 辺の両方に入れる。
- 平面 ($r\equiv1$ と置いたもの) と**同じコード経路**にし、`isAxisymmetric == 0` ではビット同一にする
  (重みを掛ける分岐を入れず、平面では重み配列を作らない)。

### 4.3 界面熱流束の診断 (`iface_q_eff`) — **`axisymMethod: 0` に限る**

`r` 重みは `axisymMethod: 0` だけ (`variables.cpp`:530)。`1` の幾何は平面のまま (`:593`) で、界面残差の意味が別物になる。

- `axisymMethod: 0`: 分母を**流体面の $r$ 重み面積** $\sum_f |S_f|\,\max(r_f, r_{\rm floor})$ (`variables.cpp` と同じ床) にする。
  $Q$ 側 (`iface_Qf_eff`) は変えない。
- `axisymMethod: 1`: `iface_q_eff` / `iface_q_eff_raw` / `iface_Qf_eff` を**有効な値として出さない** (NaN を書き、ログと h5 属性に「未検証」と残す)。
  **診断単独の軸対称利用** (CHT 以外の `interfaceDiag: 1`) にも効く — CHT の拒否ガードはここまで届かないので、この処置が要る。
- **平面では変えない**。
- 熱流束の比較は**どちらの面積か明記**する (G-if ①・$D_f$ は固体側 $A_i^r$、壁熱流束の報告は流体側)。

### 4.4 起動時検査

受理: `isAxisymmetric: 1` かつ `axisymMethod: 0` かつ `mode: fem2d`。拒否:

- `axisymMethod: 1`、`mode: local1d`
- 固体節点に $r<0$
- **界面節点の半径を直接検査**し $r_i\le r_{\rm floor}$ (軸上) を拒否 ($A_i^r>0$ は判定にならない — 軸上端点でも $A_i^r=L/6>0$)
- **連結成分ごとに Robin 拘束** $\sum_{\Gamma_R\cap\text{成分}} h\,L\,\bar r>0$ を要求 (剛性単独は半正定値。軸上だけに Robin 辺がある成分は零空間)

固体**内部**の軸上節点 ($r=0$) は許す (要素重心は $r>0$)。

### 4.4b 固体の出力熱量

`q_hole` (Robin の持ち去り) は現在 `Lh(T−Tc)/2` の**平面・集中**積分 (`solidFem2d.cpp`:301)。軸対称では**組立てと同じ consistent 行列** $M^e(T-T_c)$ で出す
(組立てだけ直すと出力収支が合わない)。平面は現行のまま (節点値は違うが総和は consistent と一致し、既存出力を変えない)。

### 4.5 参照オラクル

`tools/solid_fem2d.py` の `Fem2DOperator` に `axisym=False` 引数を足し、同じ式で組む。
`test_solid_fem2d_cpp.py` に軸対称の組立て・求解の突き合わせを足す (§6 V-ax1)。

## 5. 実装ステップ

1. `conjugateWall.cpp` — 軸対称の拒否 (**済**、2026-09-27)
2. `conjugate/solidFem2d.{hpp,cpp}` — 重み ($r$) を持たせ、剛性・Robin・集中量に適用。`solid_fem2d_tool` に軸対称フラグ
3. `tools/solid_fem2d.py` / `test_solid_fem2d_cpp.py` — オラクルと突き合わせ
4. `conjugateWall.cpp` — 拒否を §4.4 の受理条件に置き換え、`fillInterfaceDiagnostics` の分母を §4.3 に
5. `case/59.conjugate_annulus/` — 同心円環の検証ケース (メッシュ生成・固体生成・評価器)
6. `methods/boundary.md` を「計画中」から現行仕様へ

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | **軸対称の拒否** (済 2026-09-27) | `initConjugateWalls` で `isAxisymmetric == 1` を全モード拒否。ローカル 1 step で rc=1 を確認、平面はガード通過。**V-ax0〜3 がすべて合格するまで解除しない** | O |
| 2 | **§4 / §6 を上位に諮る** | 判断: 2026-09-27、codex (diagnose) に諮った: [`notes/reviews/2026-09-27-cht-axisym-design-diagnose.md`](../../notes/reviews/2026-09-27-cht-axisym-design-diagnose.md) — Critical 0 / Major 5 / Minor 1、**全件採用・却下 0** → §3・§4.3・§4.4・§4.4b・§6・§8 を改訂。次は codex plan レビュー | F |
| 3 | 固体の $r$ 重み (C++ + オラクル) + 出力熱量 | §5 ステップ 2–3 と §4.4b。合格: §6 V-ax1 (a)〜(e) | O |
| 4 | 受理条件・`iface_q_eff` の分母 (method 0 限定)・method 1 の無効化 | §5 ステップ 4、§4.3・§4.4。合格: §6 V-ax0 | O |
| 5 | 同心円環 (`case/59`) と判定 | §5 ステップ 5。合格: §6 V-ax2 | O |
| 6 | 同軸円板 (半径が変わる界面) | 合格: §6 V-ax2b。**FAIL なら止めて諮る** | O→F |
| 7 | 既存の軸対称 `interfaceDiag` 報告値の洗い出し | `grep` で `isAxisymmetric: 1` かつ `interfaceDiag: 1` の run と、`iface_q_eff` を読むツール・報告を列挙。誤っていた値の撤回要否を決める (§8-1) | O→F |
| 8 | 平面の回帰 | 合格: §6 V-ax3 | O |

## 6. 検証 (事前登録。結果を見る前に書いた — 2026-09-27、codex diagnose を反映して改訂)

合否はすべて **全域 FP64 ビルド**・`FORGE_CUDA_BLOCKSIZE=128` で取る (親 plan #94)。run は AWS。
**軸対称 CHT の拒否は V-ax0〜V-ax3 がすべて合格するまで解除しない** (円筒殻が合った段階で解除しない)。

### V-ax0 起動時検査 (負例、各 1 step)

`case/59` の受理入力を 1 か所ずつ変える: `axisymMethod: 1` / `mode: local1d` / 固体節点 1 つを $r<0$ / **界面の軸上端点のみ** / **界面の軸上辺** /
Robin を持たない連結成分。すべて rc≠0 と固有メッセージ。受理入力は起動する。`axisymMethod: 1` + `interfaceDiag: 1` (CHT なし) で
`iface_q_eff` が NaN・属性「未検証」になることも確認する。

### V-ax1 固体単体 (流体なし、`solid_fem2d_tool` と Python)

- (a) C++ ↔ Python: $K$・$b$ の相対差 ≤1e-12、同一荷重での解 ≤1e-9 K
- (b) **独立な辺積分**: 半径が変わる辺 ($r_a=1,r_b=2$ など 3 通り、斜め辺を含む) で Robin 行列・荷重・$A_i^r$・`q_hole` を
  **閉形式の厳密値** (本 plan §4.2 の式を評価器に独立に書く) と照合、相対 ≤1e-13。**C++↔Python の一致だけでは共通の誤りを排除できない**ため
- (c) **厚肉円筒殻**: 内面 $r_1$ に一様熱流束 $q_1$ (界面荷重)、外面 $r_2=2r_1$ に Robin。$T(r)=T_c+q_1r_1\left[\frac{1}{hr_2}+\frac{\ln(r_2/r)}{k}\right]$。
  半径方向 4/8/16 層で**収束率 ≥1.8**、16 層で内面温度誤差 ≤0.1 % of 殻の温度降下。**平面のまま組むと FAIL** することも確認 (検出力)
- (d) **半径が変わる界面・Robin 辺の厳密試験 — 軸対称の円板**: 固体 $x\in[0,t]$、$r\in[r_1,r_2]$ ($r_2=4r_1$)、界面 $x=0$ に一様熱流束、
  Robin $x=t$。解は $x$ の 1 次関数で $r$ に依らない (線形要素で厳密) → **全節点温度の誤差 ≤1e-9 K**、`q_hole` の節点値 = consistent 積分の厳密値 (≤1e-12)
- (e) 収支: 界面入熱 = Robin 持ち去り (相対 ≤1e-12)、(c)(d) とも

### V-ax2 ソルバ内連成 — 同心円環 (`case/59.conjugate_annulus/run_*_annulus`)

case/52 の軸対称版。**幾何・物性 (登録値)**: 静止ガス $r_a$=5 mm ≤ $r$ ≤ $r_b$=10 mm、$k_f$=0.0241 W/mK 一定 (`viscMethod: 0`・`thermCondMethod: 0`)、
$p_0$=1013.25 Pa (case/52 と同じく緩和を速める)、内壁 $r_a$ 等温 $T_a$=350 K (非連成)、外壁 $r_b$ = 共役壁、固体殻 $r_b$ ≤ $r$ ≤ $r_c$=20 mm、
$k_s$=0.027 W/mK、外面 Robin $h$=1e8、$T_c$=300 K。軸方向長さ $L_x$=2 mm、両端 slip。**軸を含まない**。流体格子 半径方向 32 × 軸方向 4 (一様)、固体 16 層。
**解析解** (単位軸長・ラジアンあたり): $R'_f=\ln(r_b/r_a)/k_f$=28.761、$R'_s=\ln(r_c/r_b)/k_s$=25.672、$R'_h=1/(hr_c)$=5e-7 [m·K/W·rad]、
$Q'=(T_a-T_c)/\sum R'$=**0.91855 W/(m·rad)**、$T_{w,*}=T_c+Q'(R'_s+R'_h)$=**323.581 K**、固体の温度降下 **23.581 K** (0.5 % = 0.118 K)、
$q_*=Q'/r_b$=91.86 W/m²。**節点荷重の総和との比較は $Q'L_x$ = 1.8371e-3 W/rad**。平面近似 ($t/k_s$ を $r_b$ で) との差は **+44 %** (検出力)。

合格 (すべて):

- $|T_w-T_{w,*}|\le$ **0.5 % of 23.581 K** (全界面節点の max)。事前予測 ≤0.2 %
- **G-cons**: 流体の内壁放熱・固体入熱 $\sum Q_f$・Robin 持ち去り $\sum$`q_hole` の相互差 ≤0.5 % of $Q'L_x$
- **連成保存**: $\max_i|Q_{{\rm sol},i}-Q_{f,i}|/A_i^r\le$ 0.1 % of $q_*$ (0.092 W/m²)
- **G-if** (`check_cht_interface.py`、全域、`conjugate.gate` に登録): `eps_abs_Wm2` **0.46** (= 0.5 % of $q_*$) / `eps_rel` **1e-3** /
  `dT_K` **2e-3** (≈1e-4 × 降下) / `tol_solid` **1e-9 W/rad** (節点荷重 ≈9e-5 W/rad の約 1e-5) / `n_consec` **80**
- **準定常**: 壁温は **$T_w-T_c$** を系列として渡す (平均 ≈ 降下なので `--drift/--osc 0.001` = 降下の 0.1 %)、熱流束は $q$ そのもの。
  `check_quasisteady.py --series-csv`、`--tail 0.5`、帯平均でなく**全界面節点**。**絶対温度を渡さない** (codex の合成試験: 絶対温度では降下の 1.97 % のドリフトが STEADY になった)
- **流体残差**: `check_convergence.py` が `NOT CONVERGED` のときは、次を**すべて**満たす場合に限り「静止場の丸め床」として報告項目にする —
  (i) 全残差列が上昇していない、(ii) 末尾 20 % で流体・固体の温度場の最大変化 ≤1e-3 × 降下、(iii) `rms_*` が丸め床の推定 (`tools/rounding_floor.py`) の 10 倍以内。
  1 つでも外れたら**未収束**として不合格
- **感度**: 固体 8/16 層、流体 半径方向 16/32 で $T_w$ 差 ≤0.1 % of 降下

### V-ax2b ソルバ内連成 — 同軸円板 (半径が変わる界面の局所整合)

静止ガス層 $x\in[0,H]$ ($H$=5 mm)、$r\in[r_1,r_2]$=[5, 20] mm (軸を含まない、$r$ 端は slip)、$x=0$ 等温 350 K、$x=H$ = 共役壁、
固体 $x\in[H,H+t]$ ($t$=5 mm、$k_s$=0.027)、$x=H+t$ Robin ($h$=1e8、$T_c$=300 K)。解は $x$ の 1 次関数で**$r$ に依らない**:
$q=(350-300)/(H/k_f+t/k_s)$、$T_w=T_c+q\,t/k_s$。**界面・Robin 辺の半径が 4 倍変わる**ので、§3 の流体半辺積分と固体集中量の違いが
局所の壁温の歪みとして出るならここで出る。合格: **全界面節点で $|T_w-T_{w,*}|\le$ 0.5 % of 降下**、かつ V-ax2 と同じ G-cons・連成保存・G-if・準定常。
**FAIL したら荷重の再配分を入れずに止めて諮る** (荷重を面積で掛け直すと親 plan の Mark II 後縁 +29.9 % の歪みを再生産する)。

### V-ax3 平面の回帰

- **固体単体**: 平面の組立て・求解は新旧コードで**ビット一致** (`solid_fem2d_tool`、C3X・Mark II・case/58 の固体)。平面の演算順序は変えない
- **GPU 連成**: 基準バイナリ = `1a389fef` の FP64 ビルド、比較量 = 最終の界面壁温 `Ts` と `conjugate_history.csv` の `Tw_mean`、ノルム = 全界面節点の max、
  入力 = `case/58` の `run_0012_v6p_df5_i50_nlog` と同一 (10000 step)。**先に同一バイナリ 2 回のノイズ床 $\epsilon_0$ を測り**、
  合格 = 新旧差 ≤ max(2$\epsilon_0$, 1e-9 K)。この式は結果を見る前に固定した

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |

## 7. 影響範囲

- `solver_density_cuda/conjugate/solidFem2d.{hpp,cpp}`、`conjugate/solid_fem2d_tool.cpp`、`conjugateWall.cpp`
- `solver_density_cuda/tools/solid_fem2d.py`、`test_solid_fem2d_cpp.py`、`check_cht_interface.py` (必要なら)
- 新規 `case/59.conjugate_annulus/`
- 平面 2D の既存ケース (case/48・52・53・54・58) は §4.2 の設計で不変 (V-ax3 で確認)
- `methods/boundary.md`、`methods/index.md`

## 8. 未確定事項

1. ~~既存の軸対称 `interfaceDiag` 報告値の誤り可能性を本 plan で扱うか~~ **決着 (2026-09-27, codex diagnose)**: **本 plan で扱う** (§5.1 #7)。
   `interfaceDiag` 全項目が誤りとは限らないので、`iface_q_eff`・`iface_q_eff_raw` とその利用先 (軸対称 run の報告・ツール) を洗い出して誤りの有無を決める。

## 9. 変更ログ

- `2026-09-27` — codex (diagnose) の Major 5 / Minor 1 を全件採用して §3・§4・§6・§8 を改訂 (軸上判定を半径の直接検査に、診断の分母を method 0 限定、出力熱量の consistent 化、半径が変わる辺の独立試験と同軸円板、準定常を温度降下基準に、G-if の登録値、平面回帰の許容を事前固定)。
- `2026-09-27` — 初稿 (親 plan §5.1 #95 から切り出し、ユーザ指示「起票していいよ」)。軸対称が**黙って間違っていた**ことをコードで確定し、全モードで拒否 (§5.1 #1)。§4/§6 は上位に諮る前の草稿。
