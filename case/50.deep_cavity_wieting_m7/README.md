# 50. Wieting 深キャビティ (NASA TN D-5908) — 層流・冷壁の一次検証

計画: [`plans/active/case-hypersonic-gap-heating-validation.md`](../../plans/active/case-hypersonic-gap-heating-validation.md) (T0/T1)。
文献基盤: [`notes/investigations/hypersonic-gap-cavity-heating-survey.md`](../../notes/investigations/hypersonic-gap-cavity-heating-survey.md)。

**問い**: 平板基準 $q/q_{fp}$ に対して、深いキャビティ内部の熱伝達はどうなるか。
case/49 (M5 環状すきま $D/W$=20) の外部基準を作るための、**2 次元・層流・冷壁**の一次検証。

## 条件と幾何 (正本は JSON)

- [`conditions.json`](conditions.json) — Table IV 全 12 系列 (300 dpi 描画から目視。OCR は行がずれるので使わない)
- [`geometry.json`](geometry.json) — 図 1 / Table I–III。**後壁は前縁から 6.24 in = 158.5 mm** (旧稿は圧力孔 3.5 in と誤認)、
  前縁公称 $r$=0.76 mm、深さ 20.32 mm、幅 1.270/4.293/7.772/10.643 mm
- 試験気体は**メタン–空気燃焼生成物**。[`tools/gas_model.py`](tools/gas_model.py) が Tt から当量比を逆算し
  (Tt 1639 K → φ 0.589, 燃料 3.30 %)、NASA-9 + Chapman–Enskog/Wilke で cp・γ・μ・λ・Pr を作る。
  **空気で自己検算 PASS** (μ ±2 %, Pr 0.69 vs 0.707)。
- 壁 294 K 等温 (薄肉過渡法の計測時壁温)

## チェーン

```
gen_mesh.py --case T0|T1 [--wd 0.063]      # 平面 2D 構造 (node)。T0/T1 は同一ブロック構成
tools/conditions.py                        # ゲート A (状態の再現) / ゲート B (q_fp の照合)
tools/make_case.py --run ... --mesh ...    # TP 5 種 + 層流 + 段階起動 (soft→mid→2次ランプ→本段)
tools/plate_eval.py RUN [--series]         # 壁 q_w(x) を場から抽出し解析値・文献値と照合
```

- forge 設定: node / SLAU / `thermalMethod: 2` (NASA-9, 5 種) / **`viscMethod: 2`** (Chapman–Enskog + Wilke/Mason–Saxena
  = Python 側と同じモデル族) / 層流 (`turbulence.model: none`) / 陰解法 `timeIntegration 11` + `blockDPLUR`
- **出口前に slip バッファ** (模型長 197 mm → 出口 230 mm)。無いと 2 次段で出口境界層 (x=200, y≈1 mm) から発散する

## 結果 (2026-09-19 時点)

- **ゲート A**: 12 系列すべてで $(T_\infty,p_\infty,\rho_\infty,U_\infty)$ が一意に解けた。
  再構成 $p_t$ / 実測 $p_t$ = 0.919–0.981 (差はノズル全圧損失 + 物性差。**合わせに行かない**)
- **ゲート B**: 解析 $q_{fp}$ は Table IV より **+15〜20 %** → **FAIL** (判定ライン ±5 %)。切り分け:
  気体モデル (燃焼生成物 vs 空気 + Pr 0.75) で約 7 点、残り約 10 点は未説明。
  **有力仮説 = 前縁鈍化 ($r$=0.76 mm) のエントロピー層** (解析式は鋭前縁、CFD なら再現できる) → 未検証
- **T0 (run_0003)**: forge $q_w$ = 34.56 kW/m² vs 解析 34.59 → **forge/解析 = 0.999**。
  → **抽出器と物性の整合は取れている**。+18 % は forge 側ではなく「我々の物性・式 ↔ 文献の報告値」の差

## T1 (深キャビティ w/d = 0.063) の所見

### 三者比較 (run_0006, 80000 step; 分母は文献 Table IV の $q_{fp}$)

| $x/d$ | 実験 (Fig 6a) | 理論 (Burggraf, B4) | forge | 実験/理論 | forge/理論 |
| --- | --- | --- | --- | --- | --- |
| 0.010 | 0.55 | 0.469 | 0.829 | 1.17 | 1.77 |
| 0.030 | 0.42 | 0.187 | 0.221 | 2.25 | 1.18 |
| 0.050 | 0.34 | 0.113 | 0.066 | 3.02 | 0.59 |
| 0.100 | 0.17 | 0.053 | 0.0052 | 3.24 | 0.098 |
| 0.200 | 0.07 | 0.023 | 0.0001 | 3.05 | 0.004 |
| ≥0.40 | 0.02 (計測限界帯) | 0.004–0.010 | ~0 | 2–5 | ~0 |

開口面積平均: forge **0.580** (文献分母) / 0.492 (forge 分母) vs 実験 **1.07**。

- **forge は上端 0.6 すきま幅までは理論と一致** (x/d 0.03 で 1.18 倍) し、そこから急落する。
- **実験は全域で理論の約 3 倍**。W70 はこれを「粘性コアだから」と説明するが、
  **forge は粘性計算なのに理論より下**に出る。
- 理論の適用条件 (B1) $N^*_{Re,L,cr}=240(L/w)^{4/3}(1+d/w)$ を実装すると **2.52e6** となり、
  試験値 1.4–1.8e5 より大きい → **粘性コア = 理論は適用外** (W70 本文と同じ結論)。
  **副産物: (B1) が幾何の独立確認になる** — $L$=157.2 mm で W70 の表値 (2.50) に一致、
  $L$=88.9 mm (圧力孔位置) なら 1.18e6 で合わない。

### 理論解の作り方と壁温依存 (Appendix B)

$$q_s/q_{fp} = 0.6\,Q(X) = \frac{0.21}{\sqrt{1+d/w}}\left[\zeta\!\left(\tfrac12, \tfrac{s}{2(w+d)}\right)
 - \zeta\!\left(\tfrac12, \tfrac{s+w}{2(w+d)}\right)\right]$$

- 形状だけの分布 $Q(X)$ (一般化 Riemann ゼータ) に、**無限深キャビティの平均値 = Chapman の 0.60 $q_{fp}$** を掛けて絶対値にする。
- **$Pr=1$ を仮定**しているので、その仮定の下では**壁温は比から完全に落ちる**
  (コアは全エンタルピー、平板の回復エンタルピーも $Pr$=1 なら全エンタルピー)。
- 実際は $Pr\approx0.69$–0.75 なので落ちず、比は $(h_t-h_w)/(h_{aw}-h_w)$ 分だけ壁温に依存する。
  本条件 ($T_t$ 1639 / $T_{aw}$ 1499 / $T_w$ 294 K) で **+12 %** 程度。冷壁なので効きは小さい。
- **壁が熱いほど効きは増大**し、$T_w\to T_{aw}$ で発散する。**case/49 は平板が断熱なので $q_{fp}=0$、
  この「$q$ 比」正規化自体が定義できない** → case/49 側は $h$ 比 (TN D-8233 形式) を使う。
  実装は [`tools/burggraf.py`](tools/burggraf.py)。

### 圧力比 (G3) — 分母問題の影響を受けない検証点

| 量 | forge (run_0006) | 実験 (W70 Fig 4) |
| --- | --- | --- |
| $p_c/p_m$ (床+1.27 mm ÷ 前縁 88.9 mm の面圧) | **0.979** | 1.0 ± 0.1 (Re・スパン長に依らず) |
| $p_m/p_\infty$ | 1.050 | — (前縁の粘性干渉で僅かに高いのは妥当) |

**G3 PASS**。キャビティ内部がほぼ模型面圧に等しいという開口キャビティの基本挙動は再現できている。

### 収束状態 (この数字の扱い)

- `check_quasisteady` = **TRANSIENT-UNSETTLED** (drift 4.2 %/tail)。$Q_c$ は 47.0 → 21.6 W/m と単調減で、
  増分減衰から**漸近値 ≈ 21.5 W/m** (平均比 0.578)。**定常解は「内部が壁温と等温 = 深部 q→0」に収束する**
  (静止気体を 294 K の壁が囲むので当然)。
- 一方で**残差はリミットサイクル** (rms_roe 0.06–2.66、主要周期 ~9700 pseudo-step)。
  → 定常解が無い可能性があるため **URANS (dual-time) を run_0007 で実施中**。

### 実験との差の候補

1. **計測側**: [`tools/skin_conduction_estimate.py`](tools/skin_conduction_estimate.py)。
   感温面は 0.305 mm の 304SS 薄板で、挿入 0.8 s の間にスキン内を熱が拡散する長さは
   $\sqrt{\alpha_s t}$ = **1.74 mm = 1.4 すきま幅 = $x/d$ 0.086** — **実験が理論を上回り始める深さ帯と一致**。
   横方向伝導が作る見かけ熱流束は ΔT 5–22 K で $0.08$–$0.7\,q_{fp}$ と、実測の超過と同じ桁。
   W70 本文も「表面伝導・放射の補正はしていない」と明記。
2. **非定常**: **否定された** (run_0007)。URANS (dual-time BDF2, dt 5e-8 s, 物理 195 µs = せん断層 115 周期) で
   $Q_c$ = 21.54 W/m — 定常 21.60 W/m と **0.3 % 差**。せん断層振動による深部輸送の増強は無い。
3. **未検証**: 低マッハ前処理、すきま方向の格子細分。

## 計算 run 一覧

| `run_*` | 目的・主要設定差分 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `run_0001_T0_A1` | T0 初回 (バッファ無しメッシュ `t0_wd0.063_y4um`, 本段 cfl 4) | **発散**: 2 次段 step 444 で出口境界層 (x=200 mm, y≈0.96 mm) が非有限。`res_nan_*` あり | 破棄予定 (記録として保持) |
| `run_0002_T0_A1` | T0 本番 (`t0_v2` = slip バッファ付き, ランプ 0.5/1/2 → 本段 cfl 2 + relax 0.7, 8000 step) | 完走。最終場は `mesh.h5` (段間は `interp_field` で引き継ぎ、res は削除) | active (run_0003 の IC) |
| `run_0003_T0_A1_main` | run_0002 の本段を継続 (8000 step, res 1000 step 毎) = **分母診断 (G2)** | `q_w`(x=157.9 mm) = **34.56 kW/m²**、解析比 **0.999**、文献比 1.180。`check_quasisteady --series-csv q_series.csv` = **STEADY** (drift 0.5 %/tail)。`wall_q.csv` / `q_series.csv` | active (**T0 基準**) |
| `run_0004_T1_wd0063_A1` | T1 = 深キャビティ $w/d$=0.063 ($d/w$=16), 同一条件・同一ブロック構成 (`t1_wd0063_v2`), 段階起動 + 本段 10000 step | 完走 (発散なし)。最終場は `mesh.h5` | active (run_0005 の IC) |
| `run_0005_T1_wd0063_main` | 本段継続 10000 step (res 1000 step 毎) = 最初の評価 | `check_quasisteady` **DRIFTING** ($Q_c$ 47.0→31.7 W/m, −12.8 %/tail)。暫定比較は上表。`cavity_eval.json` / `cavity_rear.csv` / `compare_rear.png` | active (**未収束**) |
| `run_0006_T1_wd0063_long` | run_0005 から 80000 step (res 5000 step 毎) — 内部の緩和を追い込む | $Q_c$ 21.6 W/m (漸近 ≈21.5)、平均比 **0.580** (文献分母)。`check_quasisteady` TRANSIENT-UNSETTLED (4.2 %/tail)。残差はリミットサイクル。三者比較表は上。`compare_rear.png` | active (**T1 主結果**) |
| `run_0007_T1_urans` | URANS (dual-time BDF2, dt 5e-8 s = 物理 CFL 4.4, 20 subiter, 4000 step = 195 µs) — 定常のリミットサイクルが物理的非定常かを見る | $Q_c$ = **21.54 W/m** (定常 21.60 と **0.3 % 差**)。**非定常による増強は無い** | active (**非定常の否定**) |
| `run_0008_T1_wd0211` (+`_long`) / `run_0009_T1_wd0383` (+`_long`) | **w/d 掃引** (0.211 / 0.383)。W70 で理論と実験が一致する広い側で forge を検証する | 計算中 | 投入中 |

## 既知の注意

- 段間引き継ぎで `res_*` を消すので、**評価用の時系列は「本段を継続した run」で取る** (run_0003 の形)
- `check_convergence.py` は継続 run では低下桁数が小さく出る (本段の収束場から始まるため)。
  **量の定常性は `check_quasisteady.py --series-csv` で見る** (plan §4.7)
