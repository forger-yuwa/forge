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

## 計算 run 一覧

| `run_*` | 目的・主要設定差分 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `run_0001_T0_A1` | T0 初回 (バッファ無しメッシュ `t0_wd0.063_y4um`, 本段 cfl 4) | **発散**: 2 次段 step 444 で出口境界層 (x=200 mm, y≈0.96 mm) が非有限。`res_nan_*` あり | 破棄予定 (記録として保持) |
| `run_0002_T0_A1` | T0 本番 (`t0_v2` = slip バッファ付き, ランプ 0.5/1/2 → 本段 cfl 2 + relax 0.7, 8000 step) | 完走。最終場は `mesh.h5` (段間は `interp_field` で引き継ぎ、res は削除) | active (run_0003 の IC) |
| `run_0003_T0_A1_main` | run_0002 の本段を継続 (8000 step, res 1000 step 毎) = **分母診断 (G2)** | `q_w`(x=157.9 mm) = **34.56 kW/m²**、解析比 **0.999**、文献比 1.180。`check_quasisteady --series-csv q_series.csv` = **STEADY** (drift 0.5 %/tail)。`wall_q.csv` / `q_series.csv` | active (**T0 基準**) |
| `run_0004_T1_wd0063_A1` | T1 = 深キャビティ $w/d$=0.063 ($d/w$=16), 同一条件・同一ブロック構成 (`t1_wd0063_v2`), 10000 step | 計算中 | 投入中 |

## 既知の注意

- 段間引き継ぎで `res_*` を消すので、**評価用の時系列は「本段を継続した run」で取る** (run_0003 の形)
- `check_convergence.py` は継続 run では低下桁数が小さく出る (本段の収束場から始まるため)。
  **量の定常性は `check_quasisteady.py --series-csv` で見る** (plan §4.7)
