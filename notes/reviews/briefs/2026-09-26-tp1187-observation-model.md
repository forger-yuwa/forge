# 諮問: TP-1187 (case/56) の観測モデルと D1/D4 判定の確定前

AGENTS.md「モデル分担とエスカレーション」の条件 1 (§6 検証計画を変える) / 3 (事前登録値との比較が FAIL) /
4 (§5.1 に無い修正へ進む前) / 7 (result 段の解釈を確定する前) に当たるため諮る。

- plan: `plans/active/case-hypersonic-gap-heating-validation.md` (§5.1 の項目 #60 が 3D 項)
- 台帳: `case/56.gap_tp1187/acceptance.json` (ゲート T4-3D)、`case/56.gap_tp1187/conditions.json`
- 先行レビュー: `notes/reviews/2026-09-26-case-hypersonic-gap-heating-validation-result.md`
  (NO-GO, Critical 0 / Major 9 / Minor 1。**Major 9 件はすべて採用済み**)
- 原報 PDF: `papers/gap_heating/Avery_1978_NASA-TP-1187_gap_heating_laminar_turbulent.pdf`
  (PyMuPDF で読める。本文は PDF page index 3-13、Fig 5 は index 27)

## 読んでよいもの (巨大ファイルは読まない: `res_*.h5`, `*.log`, `residual_history.csv`, `plans/README.md`)

- `case/56.gap_tp1187/acceptance.json`
- `case/56.gap_tp1187/conditions.json`
- `case/56.gap_tp1187/tools/skin_forward3d.py` (観測モデル、約 200 行)
- `case/56.gap_tp1187/tools/gap3d_eval.py`
- `case/56.gap_tp1187/tools/extract_inlet_table.py`
- 原報 PDF の page index 4,5,6,7,9,10,11,12 と 27 (fitz の `get_text()` / `get_pixmap()`)

---

## 1. 観測事実

### 1-a. 幾何と試験条件 (原報)
θ=90° タイル配列、W=0.18 cm、L=15.24 と 30.48 cm、タイルのエッジ半径 r=0.25 cm、
薄肉 304SS 0.08 cm、タイル深さ d=6.35 cm。90° 配列はパネルホルダ位置 II (前縁から 188 cm)。
**位置 II では層流データが取れず (遷移的)、90° 配列は乱流データのみ**。
参照対は run 8 と run 14 (反復対): Pt_c 17.2/17.3 MPa、q_dyn 59.01/60.18 kPa、
Re' 4.26/4.33e6 /m、α 0.2°、鋭利前縁、**δ*(位置 II) = 1.62 cm**、q_FP = 63.44 kW/m²。

### 1-b. 実測 (Table III(c)/(f)、q/q_FP、run8 / run14)
L=15.24: TC87(深さ3.81cm) 0.00/0.01, TC88(2.54) 0.05/0.10, TC90(0.76) 2.44/2.77,
TC91(0.51) 3.21/3.59, TC92(0.25) 4.06/4.40, TC93 4.26/4.60, TC94 3.86/4.04, TC95 2.75/2.57
L=30.48 (番号 +25): 112 0.01/0.02, 114(1.52) 0.39/0.37, 115(0.76) 1.39/1.49,
116(0.51) 2.19/2.35, 117(0.25) 2.93/3.15, 118 3.17/3.17, 119 2.88/2.93, 120 2.24/2.33
L30/L15 比の実測区間: TC92 0.716-0.722, TC91 0.655-0.682, TC90 0.538-0.570,
TC93 0.689-0.744, TC94 0.725-0.746

### 1-c. forge の結果 (実測分母 63.44 で割った値、現行の観測モデル適用後)
TC93 1.538 / TC94 1.200 / TC92 0.779 / TC91 0.277 / TC90 0.098 / TC88 0.00 / TC87 0.00。
深さ 0.25→0.76 cm の減衰が実測 1.7 倍に対し forge 7.9 倍。**D1 は不合格**。
codex が現行観測モデルで再計算した L30/L15 比: TC92 0.6342 / TC91 0.1483 / TC90 0.01969
(生 qwall 比は 0.320 / 0.0236 / 0.0122)。

### 1-d. 壁前の全温 (`tools/tt_audit.py`、`total_quantities.total_state` 経由)
前向き壁に達するガスは冷たい。T0 − T_wall = 191 / 96 / 43 / 0.6 / 0 K
(深さ 0.25 / 0.51 / 0.76 / 1.52 / 2.54 cm)。

### 1-e. 本セッションで新たに測った 2 件 (どちらも D1 の原因候補を減らす)
- **入口 BL は合っている**: `case/56.gap_tp1187/run_0002_fp_t8_long` (2D 平板、res_150000) の
  圧縮性 δ* は x=1.8797 m (= 位置 II の 188 cm) で **1.6625 cm**、原報の 1.62 cm に対し **+2.6 %**。
  3 枚のスナップショットで 1.5176-1.5177 cm (入口断面 x=1.7075 m) と定常。
  **ただし原報 p.8 の δ* は ref.6 由来の計算値**で直接測定ではない。
- **分母のずれは forge の欠陥ではない**: forge の 2D 平板 q_FP = 77.34 vs 実測 63.44 kW/m² = +21.9 %。
  原報 p.8 が「turbulent boundary-layer theory to be approximately 30 percent higher than the
  interpolated experimental data」と自分で書いており、forge の +21.9 % はその内側 (実測寄り)。
  forge 自身の分母で正規化すると TC93 は 1.261 とさらに下がるので、**D1 不合格は分母のアーチファクトではない**。

---

## 2. 期待値と出典 (すべて原報の逐語)

- p.11: 「on the 90° walls, the heating is 4.3 times greater [than flat plate]」→ 実測ピーク 4.26-4.60 と整合。
- p.11: 「At all test conditions the turbulent heating has reached 0.1 q_FP or less at 60 percent of the
  tile depth」「the constant heating only occurs over the upper 5 percent of the gap depth」
  (θ=90°, W=0.18 はこの「その他」側)。
- p.11 (Fig 17): 「For W = 0.10 cm the maximum heating occurred at the tangency point of the edge radius
  and the top surface. As W was increased, the maximum heating rate moved around the edge radius to
  **approximately mid arc for W = 0.18 cm**」
- p.7 (データ整理): 「the time of data analysis ... was chosen as the time when the model reached the
  tunnel center line」「The model traverses 210 cm to the tunnel center line in a little over 2 s」
  「The temperature rise rate dT/dt was calculated by **averaging the model temperatures over 0.25-s
  intervals and then determining the temperature rise rate between each interval**」
  「the tile temperature response shows that **the maximum temperature rise rate occurs before the model
  reaches the tunnel center line**」「calculations made to account for conduction and radiation effects
  indicated these heat losses were **less than 2 percent** at the time data were taken」
- p.10 (層流 θ=60/75 の等値線の文脈): 「The uniformity of the heating rates in the transverse direction
  for W ≈ 0.18 cm indicates that the gap flow is **primarily two-dimensional** in nature and that the
  longitudinal gap flow has a **negligible effect** on these heating rates for gap widths less than 0.18 cm」
- p.11 (乱流 θ=60/75/90 の等値線): 「The isometric heating rate contours indicate that the flow is
  **three-dimensional**. ... The region affected by the three-dimensional flow extends at least
  **three gap widths** in the transverse and longitudinal directions」
- p.12: 「the effect of longitudinal gap length is **not as obvious**」「for θ = 90° the impingement heating
  rate for the longer length gap **does not increase** to levels higher than the corresponding levels for
  the shorter length gap. This inconsistency with data at θ = 60° and 75° is **not fully understood**
  because of limited data; however, the inconsistency **may be attributed to the thicker boundary layer
  which is the only difference in the data**」「Unpublished results ... indicate that a **periodic heating
  pattern exists along the length of longitudinal gaps** (gap lengths were 15, 48, and 84 cm ...).
  The fluctuations in the impingement heating rate indicated in figure 20 may be due to the same
  unexplained phenomenon」
- p.5: 「two instrumented stainless-steel thin-wall (0.08 cm) tiles」
- Fig 5 キャプション (index 27): 「Instrumentation for 90° thin wall tiles. L = 15.24 cm; dimensions are
  in cm. (For thermocouple designation where L = 30.48 cm, add 25 to each thermocouple number.)」
- p.4 記号表: 「s — surface distance from tile radius mid arc, cm」
- p.5-6: 熱電対の平均縦中心線からのばらつきは **±0.03 cm**、模型設置時に「thermocouple center line
  with the upstream longitudinal gap center line」を合わせている。

---

## 3. 再現条件

forge node 離散化 / SLAU / blockDPLUR / timeIntegration 11 / unsteady 0 (局所擬似時間) /
SST `wallTreatmentSST: 0` / `qAccumulatorFP64: 1` / `FORGE_CUDA_BLOCKSIZE=128`
(既定 512 は node SLAU のレジスタ上限 481 を超え起動不能。**この値は結果に効く**)。
commit `c52d4eb8` (branch `feature/gap-heating-precision`)。
入口は 2D 平板 `run_0002_fp_t8_long` の x=1.7076 m 断面から作った `inletProfile` CSV。

---

## 4. 実施済みの操作と結果

| run | 内容 | VERDICT |
|---|---|---|
| `case/56.gap_tp1187/run_0061_prod3d_cfl4` (AWS) | 生産 4,038,516 セル、梯子 74.5k + cfl_pseudo 4/implicitRelax 0.7 で +150k | quasisteady **ALL STEADY** (drift 0.0 %)、convergence **NOT CONVERGED (stalled/plateau)** |
| `case/56.gap_tp1187/run_0052_ab3d_open_long` | L=15.24 粗 1,426,740 セル | ALL STEADY |
| `case/56.gap_tp1187/run_0064_ab3d_L30_long` | L=30.48 粗 1,560,000 セル | ALL STEADY、mesh SOFT-PASS |
| `case/56.gap_tp1187/run_0065_L15_v2` | 梯子連結を修正 (M1) + y1/z1 分離した L=15.24、60k step | res_nan 0、梯子 7 段すべて `restart_field.py` が `VERDICT: OK`、壁解像 **FAIL** |
| `case/56.gap_tp1187/run_0066_L30_v2` | 同条件の L=30.48 | **実行中** (梯子 3 段目) |

### 潰した候補 (**結論ではなく潰した証拠**)
- **格子**: 粗 (1.43M) → 細 (4.04M) で TC92 1.22 倍 / TC91 1.43 倍 / TC90 1.63 倍上がる (深いほど改善)。
  減衰比 10.5 → 7.9。2 格子の Richardson 外挿 (仮定次数 p=1/2) で TC90 は 0.14-0.19 止まりで実測 2.44 に届かない。
  **ただし先行レビュー M6 の指摘どおり 2 点では次数も漸近域も決まらず (p=0.04 なら外挿 2.82)、
  これは「上限の証明」にならないと認めた**。第 3 格子は未実施。
- **入口 BL**: 上記 1-e のとおり δ* が +2.6 % で一致。格子非依存 (y=1 mm で 1106 K)。
- **腕 A 幾何 (上流開放)**: 定常点を持たない (124k step cfl 1.5 / 120k cfl 0.5 で ±145 % 振れ、
  12 枚中 6 枚に交差部の超音速節点)。生産は腕 B。
- **梯子未連結 (先行レビュー M1)**: 実際に起きていた (`run_0063` の `res_0` は `roUy` 全点ゼロ)。
  修正して `run_0065_L15_v2` で連結を確認済み。**ただし修正後の D1 の値はまだ出していない**。
- **壁解像 (先行レビュー M5)**: `check_wall_resolution.py --target 1 --over-frac 0` は FAIL。
  y1 と z 方向の第一刻みを `--z1` で分離して平板 >1 を 98.1 %→0.8 %、gap 12.0 %→7.3 % に改善したが
  まだ FAIL。**ツールの `over` は面積割合ではなく点数割合** (ツール内部で
  `100*count(yp>target)/good.sum()`)。自作の節点数え上げでは gap 438/49,494 節点 (0.9 %)・**最大 3.3**、
  超過はすべて「円弧/口」領域、上位は x=−149.90 mm, y=0, z≈0.9-1.0 mm = **縦すきま口の鋭いエッジ**
  (z=W/2)。ツールの最大 243.7 との差は、鋭いエッジで壁法線が一意でないため「法線方向の第一内点」が
  対岸の節点を拾うためと見ている (ツール自身が「幾何的特異点では最大値は格子収束しないので判定に
  使わない」と明記)。**場所については両者一致**。
- **縦すきま口のエッジ半径 (先行レビュー M7)**: 生成器 `gen_mesh_gap3d.py` が省略している
  (r=0.25 cm > W=0.18 cm なので微小でない)。横すきま側の半径は入っている。**未検証**。

### D4 の参照比について自分で一度誤った判断をしかけた (記録として)
p.12 の「δ* の違いが唯一の差」を読んで「D4 の参照比は δ* 交絡で無効」と書きかけたが、**誤り**だった。
p.5「two instrumented thin-wall tiles」+ Fig 5「L=30.48 は番号 +25」から、
**1 枚のパネルに L=15.24 と L=30.48 の計測タイルが両方載り、同一 tunnel run (8/14) で同時に測られている**。
p.12 の当該文は **90° vs 60/75° の比較**についての記述。したがって L30/L15 比は同一 δ*・同一 q_FP・
同一 run の統制された比較である。

---

## 5. 仮説 (どれも確定させていない)

(a) 観測モデルの時間演算が違う。原報は解析時刻 ≈2 s・0.25 s 区間平均の差分なのに、
    実装 (`skin_forward3d.py:115`) は t=0.25 s まで加熱した瞬時値。先行レビューの感度試算では
    0.25→2 s で TC90 の L30/L15 比が 0.0197 → 約 0.537 に動き、**実測区間 0.538-0.570 に入る**。
    304SS の横方向拡散長 √(α t) は 0.25 s で 0.97 mm、2 s で 2.7 mm で、熱電対間隔 2.5 mm と同程度。
(b) 縦すきま口のエッジ半径の省略が取り込み流量を減らしている (M7)。
(c) 定常 RANS の枠外 (原報が言う「縦すきま長に沿う周期的加熱パターン」が本質)。
(d) 壁解像不足が深部の q を過小にしている (M5)。

---

## 6. 諮りたいこと (推奨は 1 つに絞ってほしい。両論併記は避けてほしい)

1. **観測モデルの時間演算をどう定義すべきか。**
   原報の手順は「傾斜露出 (2 s で中心線到達、最大 dT/dt は中心線到達前) → 中心線到達時が解析時刻 →
   0.25 s 区間平均の区間間差分」。露出履歴が未確定 (Fig 9 の挿入履歴しかない) なので、先行レビューは
   「加熱履歴・解析時刻・平均化窓を別々に定義し、未確定な履歴は感度幅として残す」と提案している。
   **forge の定常 q 場を入力として、この手順を再現する最小限の定義は何が妥当か。**
   特に (i) 露出をステップとみなして t=2 s まで加熱する と (ii) ランプ (0→2 s で線形に q を立ち上げる)
   の差を感度幅として出す の、どちらを主にすべきか。

2. **仮説 (a) の一致を信用してよいか。**
   TC90 の比が時間を 2 s にすると実測区間にぴったり入るのは、**都合のよい丸めに見える**。
   懸念は「時間を自由パラメタにすれば 3 点のどれか 1 つは必ず区間に入る」こと。
   これを都合のよい一致でなく物理として主張できる条件 (**どの独立な量が同時に合わなければならないか**)
   を挙げてほしい。

3. **D1 と D4 のどちらを一次のゲートにすべきか。**
   D4 (L30/L15 比) は分母 q_FP が約分され、同一パネル上の 2 枚の計測タイルを同一 tunnel run で測った比
   なので試験条件の不確かさが消える強みがある。一方で原報自身が θ=90° の L 依存を「未解明」
   「周期パターンの可能性」と書いている。D1 は絶対値なので forge の分母のずれ (+21.9 %) を被るが、
   上記 1-e のとおりそのずれは原報自身が書く理論-実測差 (+30 %) の内側である。

4. **壁解像ゲートの扱い。**
   超過が縦すきま口の鋭いエッジ (幾何的特異点、かつ M7 で省略を指摘された箇所) に局在することが
   分かった状況で、(i) 例外として文書化して先へ進む、(ii) エッジ半径を入れた生成器 Version B を作って
   M5 と M7 を同時に解消する、のどちらを取るべきか。Version B は断面が z 依存で変形するため
   生成器の作り直しに相当する (数日規模)。

5. **原報 p.10 の「W ≤ 0.18 cm では縦すきま流の影響は無視でき、流れは主に 2 次元的」は、
   乱流 θ=90° にも及ぶと読むべきか。**
   p.11 の乱流等値線は「流れは 3 次元的」「3 次元の影響域は横・縦に少なくともすきま幅 3 個分」と
   逆のことを言っているように見える。もし乱流 θ=90° W=0.18 でも実質 2 次元なら、
   **2D 計算で 4.26 q_FP を出せるかどうかがはるかに安い決定実験**になる
   (2D run は既にある: `case/56.gap_tp1187/run_0013*`, `run_0014*`)。

各論点について、根拠は原報の逐語または `ファイル:行` / run の数値で示してほしい。
