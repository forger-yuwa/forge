# case/59 — Cary TN D-5863 冷却平板 (M6, 空気) = **T4-0a (判別)**

計画: [`plans/active/case-hypersonic-gap-heating-validation.md`](../../plans/active/case-hypersonic-gap-heating-validation.md) §4.12 / §6 G14 / §5.1 #61
精読メモ: [`notes/investigations/tp1187-test-system-reading.md`](../../notes/investigations/tp1187-test-system-reading.md)

## 目的

TP-1187 (case/56) の乱流分母 q_FP の出所 TN D-7275 は、実測が Eckert 参照温度法より 10–30 % 低い。
forge は理論側に乗る。**このずれが冷壁に特有か (モデル側) / 水準だけか** を、壁温比だけを振った実験で判別する。
合格条件は「実測に一致」ではなく、完全乱流域の R = St_forge / St_exp の**壁温比依存**を測ること (plan §6 G14)。

## 実験 (Cary 1970, NASA TN D-5863)

- Langley 20-inch hypersonic tunnel、**空気**、M∞ 6.02 ± 0.02、Tt 533 K、T∞ 65 K
- 鋭い前縁 (厚み 0.0038 cm) の平板、AISI 405、薄板 0.076 cm、液体窒素で内部冷却、両側に端板、α = 0° (Me 6.0)
- Pt 3.55 MPa (Re 2.6–3.1 × 10⁵/cm) と 1.83 MPa (Re 1.4 × 10⁵/cm)
- **Tw/Tt 0.2–0.7 の 9 系列** (Table II、転記済み `conditions.json`)。自然遷移で、加熱ピーク (遷移終端) は x = 14–35 cm
- St_∞ = q / (ρ∞ u∞ cp (Taw − Tw))、乱流の回復係数 0.89

## 事前登録

`acceptance.json` (T4-0a-0: 比較点 17 点、原点 A/B、|D|>10 % / <10 % / またぐ の読み方)。結果を見る前に commit `f1e56cac`。

## 計算 run 一覧

| run | 目的・主要設定差分 | 主要結果・成果物 | 状態 |
|---|---|---|---|
| `run_0001_tw02_re027` | **T4-0a-0** (比較原点だけの A/B)。系列 Re0.27_Tw0.2 (Tw 106.6 K、Re 2.7×10⁵/cm → p∞ 2231 Pa、逆算 Pt 3.59 MPa は原報 3.55 と 1 %)。完全乱流 SST (Tu 0.5 %、μt/μ 10)、メッシュ `fp_y1_1.5um_nx900` (y1 1.5 µm、平板 0.6 m、AR 957・skew 0 で `VERDICT: PASS`)、段階起動 lam→soft→mid→2 次ランプ→本段 20k (cfl 2、relax 0.7、restart_field.py で引き継ぎ全段 `VERDICT: OK`)。AWS `~/forge56` (ソルバは 04531b4 と同一ソース) | 完走・NaN 0。区間 ramp0→main で `NOT CONVERGED (still converging)` (全列 2.0–2.6 dec 低下中) → run_0002 で延長 | active (延長元) |
| `run_0002_tw02_re027_ext` | run_0001 の res_20000 から `restart_field.py` (`VERDICT: OK`) で本段と同一設定 +40,000 step (出力 2000 毎) | `NOT CONVERGED (stalled/plateau)` (延長区間で 0.1–0.4 dec、roK のみ 2.8)。比較点の St 34 列 `check_quasisteady` **ALL STEADY** (drift 0.0 %、`_st_series.csv`)。壁解像 `--target 1 --over-frac 0` **FAIL** (y1+ 平均 0.541、>1 が面積 0.4 %、最大 3.94 は前縁)。**T4-0a-0**: D = R_B/R_A − 1 = +0.091〜+0.258 (平均 +0.146、15/17 点で >10 %) → 事前登録の読みは**判定不能** (下流 2 点が 10 % 未満)。R_A 平均 1.128 [0.974, 1.338]、R_B 1.288 [1.125, 1.472] | active |
| `run_0003_prec_f32` | **精度 A/B の A** (`acceptance.json` T4-0a-P)。run_0002 の res_40000 から restart_field (md5 は B と同一)、float32 + qAccumulatorFP64、+20k | `NOT CONVERGED (stalled/plateau)` (0.0–0.1 dec)。St 系列 ALL STEADY | ref (精度 A/B) |
| `run_0004_prec_f64` | **精度 A/B の B**。同じ初期場・同じ設定で**全域 FP64 ビルド** (`~/forge56-double`、`flowFormat.hpp` の typedef のみ変更、git 4687c3cd) | **`PASS (converged)`** (全列 5.0–6.2 dec)。壁 q_w の A との差 最大 4.8e-5 (相対)、R_A/R_B/D は小数 3 桁まで同じ → **プラトーは丸め床で報告量に効かない**。T4-0a-0 は判定不能のまま (壁解像 FAIL: 前縁 0.72 mm) | active (収束した基準場) |
