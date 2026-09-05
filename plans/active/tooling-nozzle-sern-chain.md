# ノズル設計 SERN チェーン (⑤): 燃焼器出口 starting line + 2D 平面最大推力逆設計 (key point dv) + 多作動点 RANS 評価 + MOO

## メタ

- **area**: `tooling / optimization`
- **status**: `in_progress`  <!-- 2026-09-04 起票 (branch feature/sern-design)。S0–S6 完了 (同日、S6 は Euler)。node+SST は解決 (真因 = stage 間 interp 移植)。残 = §5.1 の残作業表 (最優先 = 作動点のサイクル値化、外部流ブロック、3D SST の後縁 3 重点) -->
- **related_docs**:
  - [`methods/design/overview.md`](../../methods/design/overview.md) 「SERN チェーン」節 (現在仕様。本計画と同時に起草)
  - [`design/CAPABILITIES.md`](../../design/CAPABILITIES.md) (問題タイプ `sern_2d` を 📋 で登録)
- **related_plans**:
  - 親: [`tooling-nozzle-design-tool.md`](tooling-nozzle-design-tool.md) §4.6 ⑤ (本計画で差し替え。旧「壁圧 Bézier dv + 局所帰還 + 3D FFD」は撤回)
  - 資産元: [`../accepted/tooling-nozzle-axismach-chain.md`](../accepted/tooling-nozzle-axismach-chain.md) (MOC カーネル・semi-perfect ガス・CFD-in-the-loop の運用)、[`../accepted/tooling-nozzle-moo-loop.md`](../accepted/tooling-nozzle-moo-loop.md) (`opt/` の DOE→Kriging→NSGA-II→EHVI)、[`../accepted/tooling-nozzle-axismach-viscous-deltastar.md`](../accepted/tooling-nozzle-axismach-viscous-deltastar.md) (NS 場からの δ\* 抽出)
  - 出典調査: [`notes/investigations/sern-design-method-survey.md`](../../notes/investigations/sern-design-method-survey.md) (2026-09-04。方針の根拠は全てここ)
- **created**: `2026-09-04`
- **owner**: `sano`

## 1. 目的

⑤ SERN (single expansion ramp nozzle) の設計チェーンを、調査 (上記ノート §4) の推奨どおり
**「燃焼器出口 starting line → 2D 平面最大推力理論の key point を目標にした逆 MOC → forge 2D RANS を
作動点セット (設計 NPR + オフデザイン) で評価 → 既存 `opt/` の MOO」** として実装する。
完成時には、問題定義 YAML (`type: sern_2d`) から ランプ/カウル形状・多作動点の $C_T, C_L, C_M$・
パレート (設計点推力 / オフデザイン推力 / ピッチモーメント / 長さ) が一気通貫で出て、
NASA TM X-71972 の基準形でカウル長・カウル角の傾向が再現される状態を得る。
壁圧規定は主線から外し、剥離制約の判定量と二段膨張オプション (S8) に限定する。

## 2. スコープ

- **やる**
  - 問題定義 `sern_2d` (入口 = 燃焼器出口の超音速一様状態を既定、音速スロート給気は接続点として選択可)
  - 2D 平面・非対称 2 壁 (ランプ + カウル) MOC (前進) と、平面最大推力理論 (Guderley–Hantsch/Rao の平面版) の **key point** $(M_c,\theta_c,\dot m_c/\dot m)$ を目標にしたランプ壁の逆生成
  - カウル後縁以降の自由境界 (等圧せん断層) と外部流圧の整合 (設計点は無波、オフデザインは CFD に委ねる)
  - 3 ブロック構造メッシュ (内部 / カウル下外部流 / 下流プリューム) の決定的生成 (msh4.1 直書き) と品質ゲート
  - forge 2D 平面 RANS (SST, node) を作動点セットで回す runner・力積分メトリクス ($C_T,C_L,C_M$, 剥離位置)・ゲート
  - δ\* 一発補正 (forge NS 場から抽出 → 法線オフセット → 再評価)
  - 多作動点束ねの MOO (既存 `opt/` の目的関数ベクトル化)
  - 2D パレート数点の 3D RANS 確認 (側壁・隅 R・有限スパン) と、3D 設計 (流線追跡 / FFD) の要否判定
- **やらない**
  - 壁圧 $p_w(x)$ を dv にする逆設計・局所 $p\to\theta$ 帰還 (撤回。理由は調査ノート §3)
  - NS 帰還ループ (壁を反復更新する帰還エンジン)。δ\* は一発補正のみ
  - 3D FFD in-loop (S7 で乖離が大きい場合に別 plan で再考)
  - 化学非平衡 MOC (Cain §5 の弱結合法)。凍結 semi-perfect で開始、非平衡は [`chemistry-finite-rate-h2.md`](chemistry-finite-rate-h2.md) 以降
  - 非一様入口 (接触不連続) 用の回転流 MOC (S5 オプション。初版は質量平均一様状態)
  - TBCC over/under のモード遷移・可変幾何

## 3. 前提と既存資産の照合

| 部品 | 既存資産 | ギャップ・判断 |
| --- | --- | --- |
| 平面 MOC 単位過程 | `geometry/moc_kernel.py` (`delta=0` で平面、semi-perfect ガス対応、放射源流/平面単純波則で検証済み) | **再利用**。非対称 2 壁 (上壁 = ランプ、下壁 = カウル、両角部の Prandtl–Meyer 扇、カウル後縁の等圧自由境界) の march は新規 (`moc_sern.py`) |
| 最大推力理論 | なし (③ベルは TOP 幾何 dv で理論を使っていない) | **新規** `rao_planar.py`。平面版は本文未入手のため Rao 1958 / Guderley–Hantsch 1955 / Cain 式 4.1–4.3 から自前導出し §6 の解析検算で担保する |
| 逆 MOC (目標 → 壁流線) | `geometry/moc_inverse.py` (軸目標 → C⁻ 後退 → 壁流線)。①用で軸対称・軸目標 | 構造は流用するが目標が「最終特性線上の状態」に変わる。ランプ壁は kernel 終端 C⁻ と目標 C⁻ の間を前進 MOC で埋めた場の壁流線 (質量流量法) として得る |
| starting line | `feedback/cfd_anchor.py` (CFD 場から線状態抽出)、`transonic.py` (Sauer/Hall) | 既定は一様燃焼器出口 (解析)。音速スロート給気モードでは既存 Hall + kernel MOC で $M\approx1.3$–1.5 の C⁻ まで内部ノズルを作り接続 |
| ガス | `gas/` (CPG / semi-perfect NASA-9, MOC・forge とも対応) | 再利用。凍結組成 |
| メッシュ | `meshing/mesh2d.py` (軸対称 1 ブロック、msh4.1 直書き)、case/23 の plume 6 ブロック `.geo` | **新規** `mesh_sern.py`: 3 ブロック TFI (内部 / カウル下外部流 / 下流)。mesh2d の station 配置と壁クラスタリング関数を流用 |
| 評価 runner | `evaluate/runner.py` (③ベル: 2 段起動・VERDICT・warm seed)、`runner_axismach.py` (TP ガス配管) | **新規** `runner_sern.py`: 作動点ごとに run を作る。段階起動は [`procedures/divergence-and-startup.md`](../../procedures/divergence-and-startup.md) 準拠。外部流入口 + 超音速ノズル入口の 2 入口 |
| メトリクス | `metrics/extract.py::thrust_metrics` (出口面積分)、`metrics/deltastar.py::deltastar_from_run` (NS 場から質量収支 δ\*) | **新規** `metrics/sern_forces.py`: ランプ・カウル内外面の $p,\tau_w$ 積分から $C_T,C_L,C_M$ (基準点指定)、剥離位置 ($\tau_w$ 符号)。出口面積分は検算用。δ\* 抽出は再利用 |
| MOO | `opt/` (DOE / Kriging / NSGA-II / EHVI / driver / polish) | 再利用。1 評価 = N run (作動点) の束ねと目的ベクトル化を driver に追加 |
| 3D | 3D median-dual (`discretization-median-dual-3d.md`)、case/37 Netgen ブリッジ | S7 確認 run のみ。設計側の 3D 生成は押し出し + 側壁 + 隅 R (Gmsh) |

## 4. 設計方針

### 4.1 座標・幾何 (無次元: 燃焼器出口高さ $H=1$)

平面 2D、$x$ = 流れ方向、$y$ = 上向き。入口面 $x=0$、$0\le y\le1$。**ランプ = 上壁** ($y=1$ から角部で
$\theta_{r0}$ だけ上向きに折れて膨張)、**カウル = 下壁** ($y=0$ から $\theta_{c0}$、下向き正)。カウル後縁
$x=L_{\rm cowl}$ 以降は外部流とのせん断層 (等圧自由境界 $p=p_{\rm ext}$)。ランプ後縁 $x=L_{\rm ramp}$。
幾何包絡は spec: $L_{\rm ramp}^{\max}$、後胴線 (ランプが越えてはいけない $y_{\max}(x)$)、$L_{\rm cowl}^{\max}$。
モーメント基準点 $(x_{\rm ref},y_{\rm ref})$ は spec (機体 CG 相当)。符号は頭上げ正。

### 4.2 問題定義 (`type: sern_2d`)

| 区分 | 内容 |
| --- | --- |
| spec | 入口: `inflow.mode: supersonic` ($M_{\rm in}$, $p_{\rm in}$, $T_t$, 組成 — 既定) または `sonic_throat` ($P_t,T_t$, スロート諸元 → 内部対称ノズルで $M_{\rm hand}$ まで)。作動点リスト `operating_points[]` = {名前, $M_\infty$, $p_\infty$, $T_\infty$, 入口状態の上書き, 重み $w_k$}。幾何包絡、モーメント基準点、$\delta^*_{\rm in}$ (入口境界層排除厚、既定 0) |
| derived | 設計点の $p_e/p_a$、出口高さ $H_e$ (逆設計の結果)、$C_{T,\rm ideal}$ (入口状態から $p_a$ まで等エントロピー膨張の推力 = 正規化基準) |
| dv ($d=6$) | **key point** $M_c$、$\theta_c$、質量流量比 $\dot m_c/\dot m$ / ランプ初期角 $\theta_{r0}$ / カウル初期角 $\theta_{c0}$ / カウル長 $L_{\rm cowl}$。任意固定可 |
| 目的 | $f_1=-C_T^{\rm design}$、$f_2=-\sum_k w_k C_T^{(k)}$ (オフデザイン束ね)、$f_3=C_M$ (または目標値との差)、$f_4=L_{\rm ramp}$ |
| 制約 | 幾何包絡 (逆設計結果が超えたら INFEASIBLE)、最低 NPR 作動点で剥離位置 $x_{\rm sep}/L_{\rm ramp}\ge$ 許容値、$C_L$ 符号 (任意) |

### 4.3 平面最大推力理論と key point 逆設計 (S1 で導出・検証)

制御面を「ランプ後縁 $e$ から出る最終 C⁻」に取り、質量流量一定・長さ (= $e$ の位置) 固定で推力
$\int(p+\rho u^2)\,dy$ を最大化する Lagrange 問題。Cain 式 4.1–4.3 の平面版 ($y$ 依存項が落ちる) から、
制御面上で $M,\theta$ を結ぶ乗数関係と、縁 $e$ での
$\tfrac12\rho_e w_e^2\sin2\theta_e=(p_e-p_a)\cot\mu_e$ が得られる。平面流では C⁻ 上の適合関係
($\theta+\nu=$ const) と併せて **制御面上の状態が一様** になることを S1 で確認する (成立すれば
制御面は直線 Mach 線)。kernel (入口一様流 + ランプ角部扇 + カウル角部扇 + カウル壁) の中で乗数関係を
満たす点の軌跡を求め、$c$ を「$c$–$e$ 間の質量流量 = 指定比」で決める、というのが順設計。

**逆設計 (NUAA 型)**: 順設計で反復して求める $c$ の状態 $(M_c,\theta_c,\dot m_c/\dot m)$ を **dv として
先に与え**、(i) kernel を前進 MOC で作り、(ii) $c$ を kernel 内の指定質量流量比の C⁺ 上で $M=M_c$ の点として
探し、(iii) $c$ から目標状態の C⁻ を張り、(iv) kernel 終端 C⁻ と目標 C⁻ の間を前進 MOC で埋め、
(v) 質量流量法で壁流線 = ランプ壁を抽出する。$\theta_c$ と $M_c$ から縁の関係式で設計 $p_e/p_a$ が従属に
決まる (= 設計 NPR は dv の関数)。$c$ が kernel 内に存在しない dv は INFEASIBLE として optimizer に返す。

**文献の対応 (2026-09-05 更新)**。Cain 2010 (RTO-EN-AVT-185 Lecture 12) の**式 4.1–4.3 は軸対称版**である
(式 4.2 に半径 $y$ が入る)。平面版では $\lambda_3$ の $y$ 重みが落ち、それが「制御面上の状態が一様」の根拠になる。
Cain 自身も、流線追跡で弧角が半径依存になると Rao の解析解は成立しないと書いており、$y$ 重みの扱いが要点である。

**平面版・カウル切り詰め・外部流の一次資料が見つかった (2026-09-05)**。当初「本文未入手のため自前導出」と書いていたが、
NASA が本チェーンとほぼ同一の問題を扱っている:

| 文献 | 何が書いてあるか |
| --- | --- |
| **Shyne 1988, NASA TM-100955** (88 p, `papers/nozzle_design/`) | Rao 法を **2 次元** に修正し、**カウルを切り詰めた scarf ノズル**を**外部流条件込み**で最適化。出口 M 6.0 / 外部流 M 5.0 の例。本文・導出・結果が揃った本体 |
| **Shyne & Keith 1990, NASA TM-103175 = AIAA-90-2222** (14 p, 同) | 同じ内容の会議論文版。変分の第一変分から乗数条件 (式 11–15) までが短くまとまっており**読むならこちらが先** |
| Cain 2010 §4.2 (同) | 軸対称の乗数条件 (式 4.1–4.3) の最も読みやすい記述 |
| Rao 1958, *Jet Propulsion* 28(6) 377–382 | 本体 (有償: `10.2514/8.7324`)。上 2 件が引用元として使えるので必須ではなくなった |
| Nickerson 1982, "The Rao Optimum Nozzle Program", SEA Report 6/82/800.1 | Rao 法の実装。Shyne が使った 2 次元修正の出所 |
| Zucrow & Hoffman, *Gas Dynamics* Vol. II (1977) | 乗数条件の導出が式レベルで追える古典 (Shyne の参考文献 6) |

Shyne の定式化と本実装の対応 (**要照合**): 制御面 CE 上で $I=\int(f_1+\lambda_2 f_2+\lambda_3 f_3)\,dl$ を最大化し、
$f_1 = [(p-p_a) + \rho v^2 \sin(\phi-\theta)\cos\theta/\sin\phi]$, $f_2 = \rho v \sin(\phi-\theta)/\sin\phi$,
$f_3 = \cot\phi$ (長さ拘束)。第一変分 = 0 から DE 上の 3 条件 (式 11–13) と縁 E の条件 (式 14)、
$f_3$ が $M,\phi$ に依らないことから式 15 が出る。**本実装 (`rao_planar.theta_lip`) は Cain 式 4.3 の平面形しか持たない**ので、
Shyne の式 11–15 と突き合わせて (a) 平面縮約が一致するか、(b) 式 15 が制御面の一様性に対応するかを確認する (§5.1-10)。
さらに Shyne は**カウル切り詰め位置に最適値がある** (それを超えると推力が落ちる) と結論しており、
本チェーンの `L_cowl` 上限拡大 (§8-8) と直接比較できる。

NUAA 系 (key point パラメータ化の出所) は引き続き要旨のみ: **Lv 2023 *PAS* 143** (総説) → Lv 2017 *AST* 66 →
Yu 2020 *Acta Astro.* 166 / *AST* 105。平面 MLN の既知解 (§6 検証 (ii) の $\theta_{\max}=\nu(M_e)/2$) は
Argrow & Emanuel 1988, *J. Fluids Eng.* 110 283–288。

### 4.4 カウルと外部流

カウルは $\theta_{c0}$ の直線 (S1) → 必要なら短い放物線 (S6 で dv 追加可)。カウル後縁からの波は、
設計点では後縁圧 $=p_{\rm ext}$ となるよう $p_{\rm ext}$ を derived にする (無波) か、指定 $p_{\rm ext}$ に
対する等圧自由境界 (膨張扇 / 斜め衝撃は MOC の外 → INFEASIBLE 警告) として扱う。オフデザインの
波・剥離・外部流干渉は全て forge が解く。**設計に外部流を入れない** のは NASA 1974 と同じ割り切り。

### 4.5 starting line とガス

既定は一様燃焼器出口。非一様入口 (`inflow.profile:` CSV) は S5 オプションで、初版は質量平均した
一様状態を使い、非一様の影響は forge の入口分布 BC ([`boundary-inlet-profile.md`](../accepted/boundary-inlet-profile.md))
で評価側だけに入れる。ガスは semi-perfect 凍結 (`gas.model: semiperfect`, `thermo_href_temp: 298.15` 必須)。

### 4.6 δ\* 一発補正

非粘性設計壁で forge RANS (設計点) を 1 回回し、`metrics/deltastar.py::deltastar_from_run` で
ランプ・カウル両壁の $\delta^*(x)$ を質量収支定義で抽出 → 法線オフセット → 再評価。帰還は回さない。
入口境界層は spec の $\delta^*_{\rm in}$ を入口分布 BC に反映する。効果は $\Delta C_T$ として台帳に残す。

### 4.7 評価 (forge)

- メッシュ: 2 バンド (noz: ランプ–カウル間+プルーム / bot: カウル下の外部流、`mesh_sern.py`) に、**ランプ側外部流 (top + wake、§4.11)** を足した 4 ブロック。
  TFI、壁クラスタリング ($y^+\approx30$–80 + `wallTreatmentSST=1`、AR ≤ 1000)、`check_mesh_quality.py` VERDICT。
  node config で変換 (RANS: no-slip 壁 bcond 必須)。
- 境界: ノズル入口 = 超音速 Dirichlet (全量指定)、外部入口 = 自由流、出口 = `outlet_statPress` + 逆流 Pt/Tt、
  遠方 = 自由流。`space.pRef` = 外部静圧。
- 段階起動 (**3 段、2026-09-05 確定**): **層流暖機** (`turbulence: none` + 粘性あり、1 次、cfl 0.2、2000) →
  **SST soft** (1 次、cfl 0.2、2000) → **SST 本段** (2 次、cfl 0.5、6000)。`run_staged(warm_lam_steps=, warm_lam_cfl=, soft_cfl=)`、
  YAML は `opt.warm_lam_steps / warm_lam_cfl / soft_cfl` と `evaluate.cfl_main`。作動点間は warm restart (同一メッシュ = index コピー)。
- ゲート: `check_convergence.py` PASS、`check_quasisteady.py --quantity` (力係数・剥離位置)。低 NPR の
  RSS/FSS は `OSCILLATING` 前提 → 平均±振幅で報告し、振幅が閾値超なら `SUSPECT`。
- メトリクス: $C_T=F_x/(p_{\rm in}A_{\rm in}\gamma M_{\rm in}^2)$ 系の無次元 (実装時に $C_{T,\rm ideal}$ 正規化を選ぶ)、
  $C_L$、$C_M$ (基準点)、剥離位置。NASA 流の帳簿: ランプ + カウル内面 + カウル外面 (外部流側) を含め、
  制御体積を明示する。

### 4.8 MOO

`opt/driver.py` を「1 評価 = 作動点数分の run」に拡張し、目的ベクトル $(f_1..f_4)$ と制約を返す。
EHVI は 2 目的閉形式のみなので、3 目的以上は NSGA-II 側の hypervolume 近似 (既存 `moo.py` の扱いに従う)。
DOE は $10d=60$ 点、infill 30–50 点。1 評価の実測時間 (2D RANS × 3 作動点) を S6 冒頭で取る。

### 4.9 3D 確認 (S7)

2D パレートから 2–3 点を押し出し + 側壁 + 隅 R (spec) で 3D メッシュ化 (Gmsh)、3D RANS (node) で
$C_T, C_L, C_M$ を 2D と比較。推力差 < 2% なら 2D 設計を正、揚力/モーメントの差は 3D 補正表として持つ。
乖離が大きければ流線追跡 (親場に横方向膨張が要る) か FFD の plan を別途起票。

### 4.10 作動点の定義 (サイクル値・2026-09-05 決定)

作動点は **`external` (飛行条件) だけを振ってはいけない**。飛行 $M_\infty$ が変わればインレット圧縮と燃焼加熱が
変わり、ノズル入口 (= 燃焼器出口) 状態 $(M_{\rm in}, p_{\rm in}, T_{\rm in}, \gamma)$ が従属して変わる。
初版 (run_0010 / 0017 / 0019) は `spec.inflow` を全作動点で固定し NPR を外部圧側で作っていたため非物理だった
(特に「飛行 $M_\infty 1.5$ で燃焼器出口 $M 2.5$」の低 NPR 点)。

**アンカー = NASA TM X-71972 TABLE 1** (p.31、`papers/nozzle_design/`。S4 で傾向照合に使った当の文書)。
$\bar q_\infty = 71850\ \mathrm{N/m^2}$ (1500 psf)、$\alpha = 2°$ の定動圧上昇経路上で、station 1 (インレット前) と
**station 3 (燃焼器出口 = ノズル入口)** の $p, T, V$・当量比 $\phi$・作動モードが与えられている。$\phi = 0$ 行は燃料遮断。

| $M_\infty$ | $\phi$ | mode | $p_1$ [Pa] | $T_1$ [K] | $V_1$ [m/s] | $p_3$ [Pa] | $T_3$ [K] | $V_3$ [m/s] |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 4 | 1.0 | ramjet | 10366 | 244 | 1153 | 157717 | 2343 | 1093 |
| 4 | 0 | off | 〃 | 〃 | 〃 | 17888 | 338 | 1072 |
| 6 | 1.0 | scramjet | 6248 | 282 | 1753 | 101027 | 2328 | 1621 |
| 6 | 0 | off | 〃 | 〃 | 〃 | 16327 | 448 | 1655 |
| 10 | 1.5 | scramjet | 3840 | 353 | 2982 | 57935 | 2222 | 2837 |
| 10 | 0 | off | 〃 | 〃 | 〃 | 14172 | 772 | 2831 |

表は $M_3$ と $\gamma$ を与えないので **CEA2 で埋める** (`.venv-cea/nasa_cea/FCEA2`、`tp` 問題、H$_2$-air 平衡、
$\phi = 0$ 行は空気)。$p_\infty = \bar q_\infty / (0.7 M_\infty^2)$:

| $M_\infty$ | $\phi$ | $p_\infty$ [Pa] | **NPR** $= p_3/p_\infty$ | $\gamma$ (CEA) | $R$ [J/kgK] | $a$ [m/s] | **$M_3 = V_3/a$** |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 4 | 1.0 | 6415 | **24.6** | 1.187 | 340.1 | 972.5 | **1.12** (ramjet) |
| 4 | 0 | 6415 | **2.8** | 1.399 | 287.1 | 368.4 | **2.91** |
| 6 | 1.0 | 2851 | **35.4** | 1.183 | 340.3 | 968.1 | **1.67** |
| 6 | 0 | 2851 | **5.7** | 1.392 | 287.1 | 423.0 | **3.91** |
| 10 | 1.5 | 1026 | **56.4** | 1.226 | 391.0 | 1032.2 | **2.75** |
| 10 | 0 | 1026 | **13.8** | 1.357 | 287.1 | 548.3 | **5.16** |

$\gamma$ は CEA の平衡 GAMMAs。凍結 $\gamma$ にすると $M_3$ は 1〜3 % 動く (semi-perfect TP テーブル化とあわせて確定する)。

ここから決まること:

1. **初版 smoke 値はずれている**: $M_\infty 6$ 巡航で $M_{\rm in}$ 2.5 → 実際 **1.67**、$p_{\rm in}$ 20 kPa → **101 kPa**、
   NPR 20 → **35.4**、$\gamma$ 1.4 → **1.18**。
2. **低 NPR 点の正体は「燃料遮断」であって低速飛行ではない**。TM 本文 p.21 も fuel shutdown を明示的に扱う
   (engine module drag + 上下壁圧の同時低下)。power-off 行は $M_3 = 2.9/3.9/5.2$ と全て超音速なので
   `inflow.mode: supersonic` のまま書ける = **音速スロート実装 (§4.5) を待たずに剥離評価の作動点が張れる**。
3. **音速スロートが要るのは $M_\infty 4$ powered だけ** ($M_3 = 1.12$)。TM 本文「M4 では熱閉塞で純スクラムジェットは
   化学量論比にできずラムジェット性能を使った」と整合。§8-1 の判断根拠はこれで確定。
4. **$\gamma$ が作動点ごとに違う** (powered 1.18 / power-off 1.39) → **`operating_points[].gas` の上書きが必要** (未実装)。
   逆設計 kernel の $\gamma$ は**設計点の値で固定**する (形状は設計点で 1 つに決まるため)。

**生産構成 (既定)**: 設計点 = $M_\infty 6$ powered。作動点 = M6 powered ($w$ 0.5) / M10 powered (0.3) / M4 power-off (0.2)。
`external` は自由流のまま。NPR 35.4 / 56.4 / 2.8 で設計点をまたぐ (過膨張〜不足膨張)。

**サイクル計算の位置づけ**: 表は 3 飛行点しかないので、$\bar q_\infty$ や $\phi$ を変える・$M_\infty$ を刻むには
1D サイクル (定動圧経路 → インレット全圧回復 → Rayleigh 加熱 → station 3) が要る。ただし
**検定条件 = TM の 6 行を再現できること** (全圧回復率・燃焼効率をチューニングパラメータにする)。
アンカー無しのサイクル自作は相関の当否を確かめられないので、順序を逆にしないこと。

### 4.11 ランプ側の外部流ブロック (2026-09-05 決定)

**問題** (ユーザ指摘): 現行 2 バンドメッシュはランプ後縁の先 (`top_out`) を流れに平行な線で閉じ、そこを `outlet_statPress`
(2D 既定) か slip にしている。どちらも**ランプ側に外気が存在しない**。過膨張 (m4_off, NPR 2.8) でランプ上の $p$ が
$p_\infty$ を下回っても、後縁から境界層を通って上流へ伝わる外圧が無いので剥離 (RSS/FSS) が**起きようがない**。
run_0032 の「剥離なし・滑らかな再圧縮」は領域の産物で物理ではない。カウル側は下バンド (`inlet_ext`) で外気を持つので非対称。

**決定**: 機体を「ランプの上に載る有限厚の板」としてモデル化し、その上に自由流バンドを足す **4 ブロック構造** にする
(`mesh.ext_top: 1`)。x station は既存 2 バンドと共有。

| ブロック | 範囲 | 境界 |
| --- | --- | --- |
| bot (既存) | 遠方 → カウル外面/せん断層 | `inlet_ext` / `bottom` / `cowl_out` / `outlet` |
| noz (既存) | カウル内面/せん断層 → ランプ/プルーム上線 | `inlet_nozzle` / `cowl_in` / `ramp` / `outlet`。**プルーム上線 ($x>L_{\rm ramp}$) は内部線になる** |
| wake (新) | プルーム上線 → それに平行な高さ $h_{\rm base}$ の線、$x \ge L_{\rm ramp}$ | 左端 = **base** (鉛直、slip) / `outlet` |
| top (新) | 機体上面線 $y_{\rm veh}$ ($x \le L_{\rm ramp}$) / wake 上線 ($x > L_{\rm ramp}$) → +`top_depth` | `inlet_ext` (自由流、下バンドと同じタグ) / **vehicle** 上面 (slip) / `top_out` (outlet か slip) / `outlet` |

- $y_{\rm veh} = \max(\text{ランプ } y) + $ `vehicle_clearance`、$h_{\rm base} = y_{\rm veh} - y_e$。設計ランプは $\theta_e<0$ で
  後縁が頂点より下がるので base は常に有限 ($h_{\rm base} \le 10^{-9}$ なら wake ブロックを省き top の $j=0$ をプルーム線に直結)。
- **機体上面は水平・slip**: 上面の流れは自由流 $(M_\infty, p_\infty)$ をそのまま後縁に運ぶ。後縁の外圧を $p_\infty$ にする
  **最も中立なモデル** (機体上面形状は未知なので圧縮も膨張もさせない)。後縁 (base 角) で外気は $\theta_e$ 方向へ角膨張して
  プルーム境界に沿う。
- **力の帳簿**: `vehicle` (上面 + base) は機体の力なので $C_T, C_L, C_M$ に**入れない** (NASA TM と同じく nozzle force =
  ramp + cowl 内外面)。base 圧は診断として壁出力する。
- ノード番号・セル順序は既存 2 バンドの**後ろに追加** (`sern_deltastar._structured_upper` 互換)。
- 既定パラメータ: `top_depth` 2H, `nj_ext_top` 41, `nj_wake` 9, `vehicle_clearance` 0.02H, `first_top_frac` 0.02H。
  `SernMeshParams.ext_top` の既定は False (既存 YAML・テスト不変)、cycle3op 生産構成で ON。
- 期待される効き: m4_off で後縁からの再圧縮が上流へ入り、SST で `sep_frac_ramp > 0`。Euler は剥離しないが後縁圧が
  $p_\infty$ に張り付く (run_0032 との差で外圧の到達を確認する)。3D (`mesh_sern3d`) は未対応 (§5.1)。
- **base → テーパに変更 (2026-09-05, run_0034 の結果)**: 鉛直 base 版は node Euler で soft 段 step 5 に発散し、NaN は base 角
  (x 11.07–11.9, y 2.68–2.85) に限局した (step 0 の低密度ノード分布は run_0032 と一致し、余分な 30 点が全て base 下流)。
  原因は **node の 90° 二重 slip 角** (ランプ∩base、上面∩base): カウル TE の 2 壁は ~180° で問題ないが、直交 2 壁を 1 ノードが
  持つ形は node の slip 射影と整合しない。→ `vehicle_taper` (既定 2H): 機体厚 $t_v(x) = (y_{\rm veh} - y_{\rm ramp})\times$
  係数 (x ≤ L_ramp − taper で 1、TE で 0) で後縁を**厚さ 0 に絞り TE を共有** = カウルと同じ構造。base/wake ブロックは無くなり
  3 ブロック。機体上面の流れは後端で緩い boat-tail (θ ≈ −1.7°) 膨張を受け、TE の外圧は $p_\infty$ よりわずかに低い
  (base 圧 $< p_\infty$ の実機と同じ向き)。鉛直 base 版 (`vehicle_taper: 0`) はコードに残す (cell 用)。
- **カウル TE の SST ω 発散 (run_0035)** は ext_top と無関係 (NaN は x 1.20, y −0.11 = カウル TE)。板厚 0 の node 双子ノード
  問題 → `cowl_thickness: 0.002` で解決 (run_0037 完走)。生産 YAML の既定に。
- **結果 (2026-09-05, run_0038/0039)**: 外気込みの node SST でも m4_off (NPR 2.8) は**剥離しない** (`sep_frac_ramp` 0、τ_w>0 全点、
  L_cowl 1.2 と dv 上限 2.5 の両方)。理由は **cowl TE 圧が 1.35–1.59 p∞ で不足膨張**なこと: NPR 2.8 ではカウル側衝撃が存在せず、
  ランプの過膨張 (最小 0.57–0.79 p∞) はランプ自身の転向によるもので、後流のプルーム境界からの圧縮で穏やかに回復する。
  ランプ側外気は後縁 0.1H の圧だけを変える (1.32 → 0.93 p∞)。**RSS/FSS が出るのは cowl TE 圧 < p∞、すなわち
  NPR ≲ 1/(p_cowlTE/p_in) ≈ 2.1 以下**で、TM X-71972 の 6 点にはその作動点が無い (§8-6)。

## 5. 実装ステップ

| Step | 内容 | 主要ファイル | 規模 |
| --- | --- | --- | --- |
| **S0** | `sern_2d` 問題定義 (spec/derived/dv/作動点スキーマ、過拘束検査)・幾何コンテナ (ランプ/カウル折れ線、包絡検査)・CAPABILITIES 更新 | `probdef.py`, `geometry/sern_geometry.py`, `design/CAPABILITIES.md` | 小 |
| **S1** | 平面最大推力理論の導出と key point 逆設計: 非対称 2 壁前進 MOC (角部扇・カウル後縁自由境界)、乗数関係・縁条件、$c$ 探索、目標 C⁻ 張り、壁流線抽出。解析検算 (§6) | `geometry/moc_sern.py`, `geometry/rao_planar.py`, `design/tests/run_sern_moc_tests.py` | 中 |
| **S2** | 3 ブロック TFI メッシュ生成 (msh4.1 直書き)・品質ゲート・node 変換の配管 | `meshing/mesh_sern.py` | 中 |
| **S3** | 評価 runner (作動点ごとの run 生成、段階起動、warm restart)・力積分メトリクス・ゲート | `evaluate/runner_sern.py`, `metrics/sern_forces.py`, `evaluate/health.py` | 中 |
| **S4** | 検証: (a) 設計点 Euler で MOC の $C_T$ と forge の一致、(b) NASA TM X-71972 基準形 (ランプ 20°/18.5H、カウル 6°/3.12H) の M10/M4 傾向再現、(c) 外部流ブロックの妥当性 (case/23 プリューム資産と突合) | `case/46.sern_design/` (新設、README run 一覧) | 中 |
| **S5** | δ\* 一発補正の配管 (抽出 → オフセット → 再評価) と $\Delta C_T$ 台帳。非一様入口の評価側 BC | `feedback/deltastar_sern.py` (薄いラッパ), `runner_sern.py` | 小 |
| **S6** | MOO: 多作動点束ね・目的ベクトル・制約・DOE → パレート。1 評価時間の実測と予算更新 | `opt/driver.py`, `opt/moo.py`, `design/problem_sern_*.yaml` | 中 |
| **S7** | 3D 確認 run (側壁・隅 R)・2D との差の表・3D 設計要否の判定 | `case/46.sern_design/mesh3d/`, README | 中 |
| S8 (任意) | 二段膨張オプション: 基部の壁圧プラトー規定 (④延長部の壁圧帰還と共通機構) + 延長部は最大推力 | 別 plan に切り出す | 中 |

S0→S1 は CFD 不要で先行できる。S2–S3 は S1 と並行可 (固定形状で配管)。

### 5.1 残作業 (優先順・2026-09-05 時点)

**本節が残作業の正本**。作業ログ側 ([`notes/sessions/2026-09-05-handover-sern-design.md`](../../notes/sessions/2026-09-05-handover-sern-design.md))
は使い方・罠の記録に限り、優先順はここを見る。

| # | 項目 | 内容 |
| --- | --- | --- |
| 1 | ~~作動点のサイクル値化~~ **実装・検証済 (2026-09-05, §10)**。残 = **MOO の再取得** (第 1 回 run_0040 は起動不安定で破棄、確定レシピ §4.12 で run_0054 を投入済) | `problem_moo_sst_node_cycle3op.yaml` (ext_top・板厚 2e-3・cm_min −7.0)。前提の §4.11 / 1b / §8-6 は全て決着。run_0019 の結論は非物理な作動点が駆動したので破棄 |
| 1b | ~~`opt.cm_min` の再設定~~ **決着 (2026-09-05): −7.0** | dv 箱 324 点の MOC スイープ (`case/46/sweep_moc_dv_cycle3op.{py,log}`, 成立 246): 設計点 $C_M$ min −13.8 / p10 −10.7 / med −6.9 / p90 −2.5 / max +1.5、$C_T$ 0.81–0.965。重み付き $C_M$ は設計点の ≈0.7 倍 (m10_on −6.1、m4_off +1.0 で緩む) なので −7.0 で最悪 ~15 % を落とす。旧 −2.5 は候補の 9 割を落とす |
| 2 | ~~ランプ後縁側の外部流ブロック~~ **実装・検証済 (2026-09-05, §4.11)**。残 = 剥離が出る作動点の追加可否 (§8-6) | `mesh.ext_top` + `vehicle_taper` で領域は閉じた。m4_off では形状によらず剥離しない (cowl TE 不足膨張)。剥離制約を効かせるなら NPR ≲ 2 の燃料遮断点 (M∞ ≲ 3) を外挿で足す |
| 3 | **3D SST の後縁 3 重点** | カウル後縁 ∩ 側壁後縁で $\omega \to \infty$ (run_0028)。候補 = 後縁に板厚 (2 本の剥離線を分離) / `FORGE_FREEZE_TURB=1` で soft 段を凍結乱流 / 局所 $\omega$ 上限。Euler は同格子で成立 (run_0027) |
| 4 | **側壁長 $L_{\rm sw}$ を dv/仕様に** | 3D 効果 (推力 −4.5 %、揚力反転) の主因 |
| 5 | 亜音速外部流のチャンバー構成 | 静止・地上試験用。case/23 方式。超音速外部流なら不要 |
| 6 | カウル輪郭の設計 | 今は直線固定 (Lv & Xu 2021 はカウルにも最適化理論) |
| 7 | 音速給気のスロート接続 (`sonic_throat`) | §4.10-3 より**対象は $M_\infty 4$ powered のみ**。Hall + kernel MOC で $M$ 1.3–1.5 まで対称内部ノズル |
| 8 | 非一様入口 (回転流 MOC)・凍結 $\gamma$ → 有限速度化学 | 入口分布 BC は評価側に既にある |
| 9 | CFD 領域の縮小 | `mesh.bot_depth` 3H → 超音速外部流なら 1H 程度 |
| 10 | **Shyne の式 11–15 と突き合わせ、厳密 Rao 最適の位置を確かめる** | まず NASA TM-100955 / TM-103175 (§4.3) の平面版乗数条件と `rao_planar` を照合する。次に `rao_planar.lip_residual` ($\theta_c-\theta_{\rm lip}$、0 が Rao 最適) は今は診断値で、dv には課していない。残差 0 になる $(M_c,f)$ を箱の中で求め、**MOO のパレートがその点を含むか**を見る。理論の最適点は推力側の端点として出てくるはずで、出てこなければ平面縮約か実装のどこかが違うという判定になる (§4.3 の依拠箇所の検算) |

## 6. 検証

- **単体 (S1)**: (i) 平面単純波則 ($\theta+\nu$ 保存) 1e-10、(ii) 対称極限 (カウル = ランプの鏡像) で
  平面最小長ノズルの既知解 ($\theta_{\max}=\nu(M_e)/2$、Argrow–Emanuel) を再現、(iii) $p_e/p_a=1$ で
  制御面が一様平行流、(iv) 制御面上の推力積分 = 出口面の推力積分 (保存則、<0.1%)、(v) 逆設計で
  張った目標 C⁻ を前進 MOC で再計算した状態が目標と $M$ 0.1% / $\theta$ 0.05° 以内、(vi) 格子収束。
- **CFD (S4)**: 設計点 Euler (`visc: 0.0`, `thermCond: 0.0`) で $C_T$ が MOC 値と 1% 以内、RANS で
  全作動点 `check_convergence.py` PASS (低 NPR は `check_quasisteady.py` の統計)。NASA 基準形で
  「カウル短縮 → $C_T$ 低下・頭下げ $C_M$」「カウル角 6° 近傍で $C_M$ 最良」の傾向再現。
- **メッシュ**: 各 run で `check_mesh_quality.py` VERDICT PASS を README に記録。
- **δ\* (S5)**: 補正前後で $C_T$ の変化が符号・桁で妥当 (壁摩擦込みで −1〜−3% 程度) であること。
- **3D (S7)**: 2D vs 3D の $C_T$ 差 < 2% を「2D 設計で足りる」判定基準とする。

## 7. 影響範囲

- 新規: `design/forge_design/geometry/{sern_geometry,moc_sern,rao_planar}.py`, `meshing/mesh_sern.py`,
  `evaluate/runner_sern.py`, `metrics/sern_forces.py`, `design/tests/run_sern_moc_tests.py`, `case/46.sern_design/`
- 変更: `probdef.py` (`KNOWN_TYPES` に `sern_2d`), `opt/driver.py` / `opt/moo.py` (多作動点・多目的), `design/CAPABILITIES.md`
- forge 本体: 変更なし (2D 平面 RANS・超音速入口・外部流は既存機能)。不足が出たら別 plan
- 文書: `methods/design/overview.md` 「SERN チェーン」節を実装に同期、親 plan §4.6 ⑤ / §5 Phase 4

### 4.12 SST の起動安定性 (2026-09-05 確定)

サイクル値作動点では**カウル後縁のせん断層**で `roOmega` が発散する。排気と外部流の密度・温度比が大きい
(m10_on: ρ 0.067 vs 0.012、T 2222 vs 228 K) ためで、3 作動点すべてで同じ場所 ($x \approx L_{\rm cowl}$)。効かなかった対策と効いた対策:

| 対策 | 結果 |
| --- | --- |
| カウル板厚 | **不可**。作動点で要求が逆転する (m4_off は 2e-3 で通り 5e-3 で落ち、m10_on は 5e-3 で通り 2e-3 で落ちる)。MOO では成立しない |
| soft 段 CFL 0.5 → 0.2 | 入口面∩カウル外面の角 (M∞10 で顕在化、run_0041) には有効。後縁には無効 (step 1719 に遅れるだけ、run_0042) |
| **層流暖機段**を soft の前に入れる | **有効**。平均場が立ってから乱流方程式を入れる。3 作動点とも通り、力係数は暖機なしの成功例と一致 (run_0048) |
| 本段 CFL 1.0 → 0.5 | **必要**。1.0 は限界的で、同一メッシュ・同一 config でも run により発散する (run_0048 成功 / run_0049 失敗)。0.5 で m10_on 2 回とも同一値 (run_0050/0051) |

確定レシピ = 板厚 5e-3 + 層流暖機 (cfl 0.2, 2000) + SST soft (cfl 0.2, 2000) + SST 本段 (cfl 0.5, 6000)。
検証値 (dv 既定): m6_on C_T 0.9089 (摩擦込み 0.9018) / m10_on 0.9246 (0.9175) / m4_off 0.9759 (0.9701)、全て STEADY。

### 4.13 評価の受理方針と再試行梯子 (2026-09-05 決定)

dv 箱を動かすと**作動点 × 形状の組み合わせごとに固い場所が変わる** (ランプ膨張角部 → §4.12 の丸めで解決、
カウル後縁 → m4_off × 短カウルで残存)。全点に効く単一設定は見つからなかったので、**梯子**にする。

1. **標準レシピ** (§4.12 + `mesh.ramp_fillet: 0.1`): 層流暖機 2000 @cfl 0.2 → SST soft 2000 @0.2 → SST mid (2 次) 2000 @0.2 → 本段 6000 @0.5
2. 落ちたら **緩レシピ**で 1 回だけ再試行: 各段 step 倍・cfl 半分 (4000 @0.1 ×3 → 本段 8000 @0.25)。別 run ディレクトリ (`*_retry`) に作る
3. それでも `rc != 0` の場合でも、**力係数の履歴が 4 スナップショット以上あり 3 成分とも `STEADY`** なら採用し、
   台帳に `degraded: true` を記録する。カウル後縁の破綻は**力積分が収束しきった後に起きる**ため
   (run_0067: C_T 0.95766→0.95741, trend 3.6e-5 の 6 点で `STEADY` の後に発散)、積分値は使える。
4. 上を満たさなければ FAIL (最適化側には INFEASIBLE として返る)

**`degraded` の扱い**: パレート上の点が `degraded` なら、その点は**採用前に単独で回し直して確認する**。
台帳と `pareto.json` に必ず残すこと (「収束した」とは呼ばない)。残差は元々このケースでプラトーなので、
収束判定は力係数の `steadiness` が正本 (§4.7)。

## 8. 未確定事項 (ユーザ確認)

1. ~~**入口の既定**~~ **決着 (2026-09-05, §4.10-3)**: 超音速給気 (燃焼器出口 $M_{\rm in}$ 指定) を既定で続ける。TM X-71972 の作動点では音速スロートが要るのは $M_\infty 4$ powered ($M_3 = 1.12$、ラムジェット) の 1 点だけで、残り 5 点は $M_3 \ge 1.67$ の超音速。`sonic_throat` は §5.1-7 に後回しでよい。
2. ~~**作動点セット**~~ **決着 (2026-09-05, §4.10)**: NASA TM X-71972 TABLE 1 (定動圧 1500 psf 経路) をアンカーにサイクル値で定義する。既定 = 設計点 $M_\infty 6$ powered、作動点 M6 powered ($w$ 0.5) / M10 powered (0.3) / M4 power-off (0.2)。低 NPR は「低速飛行」ではなく「燃料遮断」で作る。
3. **モーメント基準点と目標**: $C_M$ を最小化するのか目標値に合わせるのか。
4. **隅 R・側壁**: S7 で spec として評価するだけ。設計変数にするかは S7 の結果で判断。
5. **外部流を設計段階に入れるか**: 入れない (NASA 流) を既定。
6. ~~**剥離が出る作動点をセットに入れるか**~~ **決着 (2026-09-05, ユーザ判断): (b) 入れない**。TM X-71972 の 3 作動点で
   MOO を回し、剥離制約は持たず `sep_frac_ramp` を台帳記録に留める。根拠: 剥離しないのは**作動点の性質** (最低 NPR 2.8 でも
   cowl TE が 1.35–1.59 p∞ の不足膨張) であって形状の選び方ではなく (L_cowl 1.2 と dv 上限 2.5 の両方で `sep_frac` 0、
   §4.11 結果)、制約として設計を弁別する情報量が無い。NPR ≲ 2 の点を入れるには TM アンカー外の外挿が要り、
   そこまでして得るものが無いと判断した。降下・アボート側の要件が案件で立ったときに再検討する。

7. **外部流のガスが排気ガスのまま** (2026-09-05 発見、未決): `gas_states` は単一 CPG の $(\gamma, R)$ を排気にも外部流にも使う。
   §4.10 で作動点ごとに排気ガス (燃焼生成物) を入れた結果、**外部流 (空気) が燃焼生成物の $\gamma, R$ で計算される**。
   $M_\infty, p_\infty, T_\infty$ は指定どおりだが $\rho_\infty, u_\infty$ がずれる: m10_on で $\gamma R$ = 479.5 (空気 401.8) →
   $u_\infty$ +9.3 %、$\rho_\infty$ −27 %。m6_on は $\gamma R$ = 402.5 で偶然ほぼ一致 (誤差 0.1 %)、m4_off は空気そのもので厳密。
   影響: 外部流の斜め衝撃角・Prandtl–Meyer が変わり、カウル TE のせん断層が実際より強くなる (m10_on の数値不安定の一因の疑い)。
   選択肢: (a) 現状維持 (ノズル性能は排気ガスで決まるので設計目的には妥当) + 制限として明記、(b) 外部流を空気に合わせる
   (排気の膨張が狂うので不可)、(c) `gas.model: semiperfect` の多成分化 (§4.5) — 別 plan 規模。**推奨 = (a)** だが、
   m10_on の安定性が (c) を要求するなら再考する。

8. **`L_cowl` の上限 2.5 が狭すぎる** (2026-09-05 発見、要判断): 現行 dv は 0.6–2.5 H だが、**NASA TM X-71972 の基準形は
   カウル長 3.12 H** で箱の外にある。C_T は上限で飽和しておらず単調増加で、最適化が上限に張り付く:

   | 出典 | カウル長 [H] | $C_T$ |
   | --- | --- | --- |
   | S4 run_0003 (NASA 平板, M∞10) | 2.0 | 0.9605 |
   | S4 run_0004 (**NASA 基準形**) | 3.12 | 0.9704 |
   | S4 run_0005 | 4.5 | 0.9720 (ほぼ飽和) |
   | MOC dv スイープ (中央値) | 0.60 / 1.55 / 2.50 | 0.8749 / 0.9060 / 0.9161 |

   スイープの $C_T$ 上位 8 点は**全て `L_cowl` = 2.5 (上限)**、キャンペーン run_0069 の高 $C_T$ 点も
   doe_004 (2.30) / doe_014 (2.06)。→ **上限を 4.5 H まで広げるべき** (S4 で飽和が確認できている値)。
   副次的に、数値的に固い短カウル域 (`L_cowl` ≲ 0.7 で m4_off / m10_on が落ちる、§4.13) から離れる利点もある。
   **決着 (2026-09-05, ユーザ判断): 打ち切って上限 4.5 H で回し直す**。run_0069 は 22 評価 / PASS 14 で打ち切り
   (最良 doe_004 C_T_w 0.9434 @ Lc 2.30)。広い箱の MOC スイープでは C_M が Lc とともに緩む (中央値 −7.25 @0.6 →
   −2.02 @4.5) ので `cm_min` −7.0 は据え置き。`x_max_kernel` は 14 → 20。

## 9. 完了条件

- [ ] `methods/design/overview.md` の SERN 節を実装と同期
- [ ] S0–S7 完了、§6 の判定を満たす (S8 は別 plan)
- [ ] `design/CAPABILITIES.md` の `sern_2d` を ✅ に更新 (検証ケース = case/46)
- [ ] `status: done` にして §10 に変更ログ、`plans/accepted/` へ移動、`plans/README.md` 同期

## 10. 変更ログ

- `2026-09-04` — 初稿。調査ノート [`sern-design-method-survey.md`](../../notes/investigations/sern-design-method-survey.md) の推奨 (§4) を計画化。親 plan §4.6 ⑤ の壁圧 Bézier dv・局所帰還・3D FFD in-loop を撤回し本計画に差し替え。branch `feature/sern-design`。
- `2026-09-04` — **S0–S1 実装** (`geometry/{sern_geometry,rao_planar,moc_sern}.py`, `probdef.KNOWN_TYPES` に `sern_2d`,
  `tests/run_sern_moc_tests.py` 全 PASS)。平面 MOC は逆 (格子) 法 + 角部扇の解析閉包 + **TE 扇と自由境界反射の
  C⁺ レイ束** (格子だけでは扇が最初の station で 1 セルに潰れ伝播しない・楔域が古い値のまま残る、の 2 件を実測して
  導入)。閉包判定の罠: 「K⁻ が入口値のまま」は TE 扇の下流でも真になるので TE 先頭レイより上流に限定した。
  §6 解析検算: 一様流機械精度 / ランプ扇 M,θ 誤差 0 (解析閉包) / カウル反射後 M(ν_in+2θ_r) 一致 /
  **対称 MLN 極限 (M_in 1.5→M_e 3)**: θ_c=0.0000°, 出口高さ = 等エントロピー面積比 (0.01%), c–e 流量 = f (0.01%),
  K⁻ 一定 1e-12, 上半分総推力 = 出口運動量流束 (0.03%) / カウル付き SERN (M_in 2.5, 15°/5°, L_cowl 1, p_ext 0.05):
  格子倍で C_T 変化 0.001%・L_ramp 0.04%。dv 掃引 33 点 (f 0.35–0.6, M_c 3.2–4.4): 等長候補間で C_T 最大点が
  縁条件残差最小側 — 粗いので最適性の確証は別途 (fine sweep か等長拘束の直接最適化)。
  **速度**: kernel march 14 s (nj 301, dx 2e-3, x 9H)、逆設計 1 本 <1 s。図: scratchpad `sern/` (artifact 化)。
  **未対応**: semi-perfect ガス (圧力比写像)、非一様入口、p_ext > p_TE の衝撃。
- `2026-09-04` — **S2–S3 実装 + S4(a) 合格**。`meshing/mesh_sern.py` (2 バンド構造 quad、カウルはスリット =
  中間線ノードを上下 2 重、**TE 点は共有** — 重複させると TE 直後に幅 0 の隙間 = 未タグ境界辺が出る、
  `tests/run_sern_mesh_tests.py` で検出・修正)、`evaluate/runner_sern.py` (逆設計→メッシュ→品質ゲート→領域別一様 IC→
  段階起動 soft 1 次 cfl0.5 3000 step → 本段 2 次 cfl4。staging の場移植は `interp_field.py`)、
  `metrics/sern_forces.py` (壁面出力 `res_<name>_<id>_<step>.h5` の CONNE 面ごとに (p−p_a)·n を積分。cell = 面値 / node = 節点平均)。
  **S4(a)**: case/46 run_0002 (cell Euler slip, 56k セル, 品質 PASS) で C_T 0.9660 vs MOC 0.9666 (−0.05 %)、
  C_L 0.1427 vs 0.1392、C_M −0.939 vs −0.946。壁圧分布はランプ・カウル内面とも MOC に重なる
  (x/H 6–7 の後縁扇反射位置に小差)。残差は 1.4 桁プラトー (`NOT CONVERGED`) だが力係数は STEADY (1e-5)。
  `geometry.mode: straight` (平板ランプを切るだけ、NASA 照合用) を追加。S4(b) は run_0003–0007 (カウル長 2/3.12/4.5H、
  カウル角 3/6/12°, M∞10, γ1.3) を投入。
- `2026-09-04` — **S4(b) 合格** (case/46 run_0003–0007, `geometry.mode: straight`, cell Euler, M∞10/γ1.3): 内面の力は
  forge と MOC が 5 形状すべてで C_T +0.001 以内・C_M 0.05 以内で一致。カウル長 2.0/3.12/4.5H で C_T 0.9605/0.9704/0.9720
  (短縮で大きく減、延長の利得は小 = NASA)、前方基準 (−20H) の C_M −1.31/−0.68/−0.36 (短いカウルで大きな頭下げ = NASA)。
  カウル角 3/6/12° で C_T 0.9786/0.9704/0.9314 (単調減 = NASA)。**カウル外面の衝撃圧は MOC に無い寄与** (12° で C_T −0.035) で、
  評価器で必ず取る。C_M の傾向は基準点に依存するため、問題定義の `moment_ref` は機体 CG を必ず与える (§8-3 の回答)。
  詳細は case/46 README。
- `2026-09-04` — **RANS 化・S5・S6 (Euler) 完了** (case/46 run_0008–0013, run_0010_moo_euler_2op):
  (i) node + SST は本段でカウル角部直下流の内面側から ω が発散 (cfl 4/1 とも、run_0008/0009) → **未解決**。評価器は
  cell + SST (壁関数, cfl 2) で成立 (run_0011: C_T(p) 0.9685 / 摩擦 −0.0085 / C_L 0.154 / C_M −0.980、力 STEADY)。
  twall の符号は「流体に働く traction」(viscousFlux_d.cu L527/L679) で確認し `sern_forces` に反映。
  (ii) S5: 単独抽出の δ* は膨張扇を欠損と誤認 (0.24 H) → **Euler 基準の質量流束欠損** `deltastar_sern_vs_euler` を採用
  (ランプ 0.009→0.11 H)。オフセット壁 RANS (run_0013) で壁圧が MOC に全域一致、C_L 0.146 / C_M −0.930 (Euler 0.143 / −0.939)、
  C_T は不変。(iii) S6: `opt/driver_sern.py` で 2 作動点 (cruise w0.6 / accel w0.4) Euler MOO、18 評価 17 PASS、HV 1.179、
  パレート 8 点 (L 3.9–12.5 H で C_T,w 0.960–0.977)。1 評価 ≈ 135 s。RANS 化した MOO は評価器の切替 (`evaluate.model: sst`,
  `discretization: cell`) だけで回るが、1 評価 ≈ 3× のコスト。作動点 accel (p_ext/p_in 0.15) は Euler では剥離が出ないため
  RANS で再評価が要る。
- `2026-09-04` — **node + SST 発散の真因を解決** (case/46 run_0014–0016): 真因は solver ではなく、段階起動の stage 間移植
  `interp_field.py` (座標最近傍) が**スリットの座標一致双子壁ノードを同じ元ノードに写す**こと (合成場で 134/134 station 誤写像を確認)。
  排気側壁ノードが外部流の圧力を持ち 2 次で爆発していた。`runner_sern.restart_by_index` (index コピー) に替えて板厚 0 のまま完走
  (run_0016: C_T(p) 0.9691 / 摩擦 −0.0065 / C_L 0.155 / C_M −0.981、cell と一致)。板厚オプション `mesh.cowl_thickness` は任意
  (入口側で 0 に絞ると 4 境界ノードができて発散するので入口から一定にする)。twall の符号規約は node (壁に働く力) と cell (流体に働く力)
  で逆 → `sern_forces` で離散化ごとに切替。**残観察**: node の rms_roOmega が本段で 3e18 一定 (壁ノード ω ピン留め残差の混入、場は健全) と
  境界角の単ノード圧力外れ (入口角・TE)。**教訓**: 座標一致ノードを持つメッシュ (スリット/薄板) では最近傍補間の restart を使わない。
- `2026-09-04` — **S6 RANS 版 (node SST) 完了** (case/46 run_0017): 2 作動点 (cruise / accel p_ext 0.2) + C_M ≥ −2.5 制約、
  目的は摩擦込み C_T。14 評価 10 PASS、HV 0.953、パレート 5 点 (L 5.5–9.7 H, C_T,w 0.960–0.967)。1 点 76 s。
  剥離指標 (`sep_frac_ramp`: 壁に働く接線力が逆向きの長さ割合) を `sern_forces` に追加したが、この作動点範囲では全点 0 —
  剥離を評価するには p_ext/p_in ≳ 0.5 の作動点が要る (未実施)。RANS でも前線の形は Euler と同じ。
- `2026-09-05` — **設計点固定バグ修正**: 作動点ごとに kernel (自由境界圧) を再設計していた (run_0010/0017 は作動点間で形状不一致)。
  `design_from_problem(p, design_external)` で spec.external に固定。**低 NPR 作動点** (M∞1.5, p_ext/p_in 0.6) を単点評価 (run_0018:
  C_T 0.92, C_L −0.32, C_M +1.8 頭上げ、剥離なし・滑らかな再圧縮) → 第 3 作動点 (w 0.2) にした **3 作動点 node SST MOO** (run_0019,
  10 PASS): 低 NPR が支配的で長いランプは C_T 0.83–0.90 に落ち、前線は最短ランプ (L 3.8H) に退化、剥離割合 sep_frac が最大 0.10 で発火。
  チャンバー: 超音速外部流 (M∞>1) では不要、亜音速/静止 (地上試験) には未実装 (case/23 方式が要る)。
  **S7 3D**: `meshing/mesh_sern3d.py` (2 バンド押し出し hex、カウル/側壁スリット、TE/側壁後縁共有、ランプ線も内外 2 重)、
  `evaluate/runner_sern3d.py` (index IC、quad 力積分 [CONNE = 5 整数/面])。**外側空間なし基準 (run_0023) は 2D と 1 % 以内で一致**。
  外側空間あり (run_0020) は横端トポロジで soft 段 NaN が未解決 (3 バグ修正済み: 入口タグ / 2 種入口の共有ノード / IC)。
- `2026-09-05` — **kernel の打ち切り** (ユーザ指摘: 解析領域が必要以上に広い): `PlanarMOC.march(stop_at=(f, M_c))` で f 流線を
  同時積分し key point 到達の 0.3H 先で march を止める (x_max は上限)。設計はビット同一、領域 9H → x_c+0.3。CFD 側の
  `bot_depth` 3H も超音速外部流なら 1H 程度で足りる (次キャンペーンで縮める)。
  **3D 外側空間あり**: Euler は 1 次 soft 段を完走、2 次本段でノズル幅外のランプ角部 (M∞6 が 15° 凸角で p/10 に膨張) から発散。
  加速作動点 (M∞3.5) で Euler/SST を再試行中 (run_0024/0025)。
- `2026-09-05` — **S7 3D: Euler で外側空間あり構成が成立** (run_0027, 加速点, 52.5 万セル, top_out = slip):
  幅内ランプ T −0.005/L −0.019 (2D +0.014/+0.085)、幅外ランプ T −0.025/L −0.104 → C_T 0.932 (2D 0.976), C_L −0.20 (2D 0.00),
  C_M +1.07。側壁がカウル後縁で終わると排気が横に逃げてランプ圧が 1/3 に落ちる = 3D 効果は推力 −4.5 %・揚力反転で大きい。
  3D の落とし穴 (修正済み): 入口面の幅外は inlet_ext / ランプ∩側壁の共有ノードは 2 種入口 → 2 重化 / index IC のカウル外面 /
  暖機コピーは状態量のみ (wall_dist を潰さない) / top_out 静圧出口は流れ平行で不安定 → slip。M∞6 では幅外ランプ角部の真空膨張で
  2 次が落ちる (加速点 M∞3.5 で回避)。**SST 3D は後縁 3 重点で ω → inf (run_0028) が未解決**。
  壁面出力 CONNE (3D) = [5, n0..n3]。kernel は key point で打ち切り。
- `2026-09-05` — **作動点をサイクル値で定義し直す方針を決定** (§4.10 新設、§8-1 / §8-2 決着、§5.1 に残作業の正本を新設)。
  初版の作動点は `spec.inflow` を全点で固定し NPR を外部圧側で作っていて非物理だった。アンカーに
  **NASA TM X-71972 TABLE 1** (定動圧 1500 psf 経路の station 1 / station 3 状態、$\phi$、作動モード) を採り、
  CEA2 (`tp`, H$_2$-air 平衡) で $\gamma, R, a$ を補って $M_3$ と NPR を確定した。結果: $M_\infty 6$ 巡航は
  $M_{\rm in} = 1.67$ / $p_{\rm in} = 101$ kPa / NPR 35.4 / $\gamma = 1.18$ (初版は 2.5 / 20 kPa / 20 / 1.4)。
  **低 NPR 点の正体は燃料遮断** (power-off, $M_3 = 2.9$–5.2 で全て超音速) であり低速飛行ではない。
  実装は未着手 (`operating_points[].gas` 上書きが必要)。既存 MOO (run_0010/0017/0019) の結論は作動点が
  非物理なので破棄・再取得とする。
- `2026-09-05` — **作動点サイクル値化を実装・検証** (§4.10)。`runner_sern.design_snapshot` で設計点 (入口・外部流・ガス) を
  作動点適用前に控え、`design_from_problem(p, design=)` が**それで逆設計を固定**する (従来は `p.spec["inflow"]["M_in"]` を
  読んでいたので、作動点が inflow を上書きした瞬間に形状が作動点依存になっていた)。`select_operating_point` が
  `operating_points[].gas` (gamma/cp) を上書き。`collect` の `cfd_vs_moc` は設計点 run 限定 (`on_design_point`)。
  `runner_sern3d` も同形。`driver_sern` は無変更で通る。導出スクリプト `case/46/cea/tmx_operating_points.py`。
  検証: 単体テスト ALL PASS / run_0030 (prepare ×3) で 3 作動点の**輪郭 sha 一致** (L_ramp 11.36 H) / run_0031 (m6_on, node Euler)
  **C_T 0.9059 = MOC 0.9073 (−0.15 %)**, C_L 0.3135 (0.3105), C_M −8.23 (−8.17), 力係数 STEADY, `check_quasisteady` ALL STEADY,
  残差は 1.0–1.7 桁プラトーで NOT CONVERGED (run_0002 と同性格) / run_0032 (m4_off) C_T 0.974, C_M +0.96 / run_0033 (m10_on)
  C_T 0.923, C_M −6.10。**dv 範囲を張り直した**: NPR 35.4・γ 1.183 では完全膨張 M_e 3.59 で key point は M_c ≲ 3.4 まで、
  M_c 3.2 で既に L_ramp 16–20 H (旧 3.2–4.4 は不成立) → M_c 2.4–3.4, L_cowl 上限 2.5, ni_noz 160。
  **`opt.cm_min: -2.5` は使えない** (§5.1-1b)。
- `2026-09-05` — **ランプ側外部流ブロック (§4.11) を実装・Euler で検証**。`mesh_sern.py` に `ext_top` (top バンド + `vehicle` 壁、
  `_add_ext_top`)、`vehicle_taper` (後端を厚さ 0 に絞り TE 共有; 鉛直 base + wake 版は node で step 5 発散 run_0034)、
  `_geom_start`。`runner_sern`: `vehicle` bcond (slip, 壁出力)、`paste_region_ic` が `y_top` で燃焼器領域を閉じる、YAML キー
  `mesh.ext_top/top_depth/nj_ext_top/nj_wake/vehicle_clearance/first_top_frac/vehicle_taper`。mesh テスト 27 項目 ALL PASS
  (既存 9 + ext_top 12 + taper 6; 順序互換・閉境界・辺数・CCW)。**run_0036** (node Euler m4_off, ext_top テーパ): MESH PASS
  (AR 443, skew 0.464, 60k cells)、C_T 0.9736 / C_M +0.973 (run_0032 と同値)、STEADY、残差 0.8 桁 still converging。
  **後縁圧 p/p∞ 1.317 (run_0032, 領域の産物) → 0.932 (機体上面 boat-tail 膨張後の外気)** で外気の到達を確認。
  ただし m4_off (NPR 2.8) ではランプ圧の主因は**カウル側**プルーム境界 (0.5–7H で 0.84 p∞ の軽い過膨張 → 7H 以降カウル側からの
  圧縮で 1.3 p∞ へ) で、ランプ側外気は後縁 0.1H だけを変える。剥離が出るかは cowl TE 衝撃の強さ (= L_cowl 依存) の問題で、
  L_cowl 1.2 の smoke 形状は m4_off で剥離しない可能性が高い (SST run_0038 で確認中)。
  **node SST のカウル TE ω 発散** (run_0035, 板厚 0) は ext_top と無関係で、`cowl_thickness: 0.002` で解決 (run_0037 完走、
  sep_frac 0) → 生産 YAML の既定に。
- `2026-09-05` — **`opt.cm_min` を −7.0 に決定** (§5.1-1b)。dv 箱の MOC スイープで設計点 $C_M$ の分布を取り (CFD と 0.7 % 一致するので代理可)、重み付き比 ≈0.7 を掛けて最悪 ~15 % を落とす値にした。$C_T$ 上位は M_c 3.0 / θ_r0 22° / L_cowl 2.5 (L_ramp 12–16 H, C_T 0.957–0.965) に集中 —
  θ_r0 と L_cowl が上限に張り付くので、MOO 前に範囲の妥当性 (θ_r0 上限 22°、L_cowl 2.5) を見直す余地あり。
- `2026-09-05` — **S6 MOO 再取得を投入** (case/46 `run_0040_moo_cycle3op`)。§8-6 をユーザ判断で (b) に決着させ、
  剥離制約は持たず `sep_frac_ramp` は台帳記録のみとした。構成: `problem_moo_sst_node_cycle3op.yaml` (設計点 M∞6 powered、
  3 作動点 m6_on 0.5 / m10_on 0.3 / m4_off 0.2、node SST、ext_top テーパ、`cowl_thickness` 2e-3、`cm_min` −7.0)、
  dv 5 変数、DOE 50 + infill 10×2 = 70 評価。実測 1 評価 ≈ 5 分 (SST 1 run ≈ 80 s × 3 作動点) → 約 6 時間、ディスク ≈ 6 GB。
  HV 基準点 = (−0.75, 20.5) (旧既定 −0.90 は新しい C_T 帯の半分を切り捨てる)。
- `2026-09-05` — **SST 起動レシピを確定** (§4.12) し **MOO 第 2 回を投入** (case/46 `run_0054_moo_cycle3op`)。第 1 回 (run_0040) は
  m10_on が soft 段 step 3 で `roOmega` 発散して全滅 → 切り分け (run_0041–0053): 発散は**カウル後縁のせん断層**が本体で、
  板厚は作動点で要求が逆転して使えず、CFL は場所を移すだけ。**層流暖機段**を `run_staged(warm_lam_steps=)` に追加し、
  本段 CFL も 1.0 → 0.5 に下げた (1.0 は同一メッシュ・同一 config で結果が割れる限界状態)。3 作動点の検証値は暖機なしの
  成功例と一致。DOE 40 + infill 8×2 = 56 評価、1 評価 ≈ 8 分。
- `2026-09-05` — **`L_cowl` 上限を 2.5 → 4.5 H に拡大し、キャンペーンを AWS g5 (A10G) へ移して再投入** (§8-8 決着)。
  旧箱は NASA 基準形 (3.12 H) を含まず C_T が上限で張り付いていた。run_0069 は 22 評価 (PASS 14) で打ち切り。
  併せて driver の容量節約が `.xmf` を消さず「中身の無い xmf」を残していたのを修正。
- `2026-09-05` — **平面版 Rao の一次資料を発見** (§4.3)。当初「Rao 本文未入手のため自前導出」としていたが、
  **Shyne 1988 NASA TM-100955** と **Shyne & Keith 1990 NASA TM-103175 (AIAA-90-2222)** が Rao 法の 2 次元修正 +
  **カウル切り詰め (scarf) + 外部流**を扱っており、本チェーンとほぼ同一の問題設定だった。両方 NTRS から取得して
  `papers/nozzle_design/` に配置。乗数条件は式 11–15 の形で書かれており、`rao_planar` との照合を §5.1-10 に積んだ。
  Shyne は**カウル切り詰め位置に最適値がある**と結論しており §8-8 (L_cowl 上限拡大) の直接の比較対象になる。
