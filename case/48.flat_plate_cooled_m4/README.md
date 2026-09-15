# 48. 超音速冷却壁平板 (等温壁 × 低 Re SST の検証)

計画: [`plans/active/tooling-nozzle-isothermal-wall-chain.md`](../../plans/active/tooling-nozzle-isothermal-wall-chain.md) §4.3。
ノズル設計チェーンを断熱壁から等温壁 ($T_w$ 指定) に切り替える前に、**冷却された超音速乱流境界層**の摩擦・熱流束・積分厚さを
forge (node, 低 Re SST `wallTreatmentSST: 0`, 壁関数なし) が理論/経験式と SU2 に対して再現するかを確かめる。

## 条件 (case/44 va3 風洞の試験部壁を模す)

- 空気 CPG (γ 1.4, R 287, Sutherland μ₀ 1.716e-5 / T₀ 273 / S 111, Pr 0.72, Pr_t 0.9)
- 自由流: M 4.19, P 5037.4 Pa, T 283 K (ρ 0.0620 kg/m³, U 1412.9 m/s, Re/m 4.96e6), $T_{aw}$ ≈ 1174 K (r = Pr^{1/3})
- 入口乱流: k 75 m²/s², ω 26000 1/s (TI 0.5 %, μt/μ 10.1)
- 領域 x∈[−0.1, 1.0] m, y∈[0, 0.2] m。x<0 は slip 助走、x≥0 が平板 (`wall` / `wall_isothermal` Ts)。上面 slip、入口 `inlet_uniformVelocity` (超音速 Dirichlet)、出口 `outlet_statPress` (Ps = P∞)
- メッシュ: `gen_mesh.py` (平面 2D quad, node 変換)。streamwise 1000 セル (前縁側 0.22 mm → 2.7 mm)、壁法線 第一セル y₁ = 3 / 6 / 12 / 24 µm (等比 1.09)。品質: 3/6 µm SOFT-PASS (外れ値 <0.1 % = 上流 slip 助走の粗いセル)、12/24 µm PASS
- 壁温: A 断熱 / B 300 K ($T_w/T_{aw}$ 0.26) / C 700 K (0.60)

## 起動レシピ (実測, 2026-09-12)

一様 IC + SST は前縁 (x 0〜3 mm, y<60 µm) で step 114 に NaN。**層流暖機** (`model: none`, 1 次, cfl 0.2, 2000) → SST soft (1 次, cfl 0.3, 2000)
→ mid (1 次, cfl 1.0, 2000) → 2 次ランプ cfl 0.5 / 1 / 2 (各 2000) → 本段 cfl 2 (+ `implicitRelax 0.7`)。本段 cfl 4 は前縁 x≈1〜3 mm で P 床 → NaN。
IC は壁近傍を tanh ランプ (δ₀ 300 µm) で 0 に落とす。壁温切替 (warm start) は soft/mid + ランプ (`--stages soft`)。`gen_runs.py` が全段を実行する。

## 評価

`tools/cooled_plate_eval.py RUN [--ref RUN_A] [--su2 DIR] [--series] [--closure] [--out prefix]`:
壁量 ($C_f$ 2 次片側差分, $q_w=\lambda_w\,dT/dn$, $T_w$, y₁⁺)・積分厚さ (ρu 0.995 縁, δ\*, θ, H)・van Driest II (K–S 基準)・$2St/C_f$・
CONTUR 平面積分 (`deltastar_integral.flat_plate_integral`)・温度–速度関係 (Walz / Duan–Martín)・van Driest 変換分布・エネルギー閉合・時系列 drift。

## 計算 run 一覧

| `run_*` | 目的・主要設定差分 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `run_0001_A_ad_y3` | A 断熱, y₁ 3 µm, 段階起動 (層流暖機→soft→mid→ランプ 0.5/1/2) + 本段 cfl 2, 12000 step | 起動レシピの確立。cfl 4 は NaN (記録 `_failed_logs/`)。T_w 未発達 (881 K @x=0.95) | 中継 (run_0004 の IC) |
| `run_0004_A_ad_y3_cont` | A 断熱の本計算: run_0001 から +48000 step (cfl 2) | `check_convergence` NOT CONVERGED (平坦床 rms_ro 1.3e-7 / roe 0.16, k/ω falling)。T_w 1157–1160 K (回復係数 0.883), $C_f$ 1.55/1.40/1.31e-3 @x 0.3/0.6/0.9 (VD-II 比 0.87–0.92, CONTUR 比 0.90–0.96), θ CONTUR 比 −4.5 %, δ\* +7 %, H 9.3 (CONTUR 8.3)。series: C_f/δ\*/θ STEADY, T_w drift 1–1.7 % (末尾でまだ +0.5 %/12k) | active (**断熱基準**) |
| `run_0005_B_tw300_y3` | B $T_w$ 300 K, y₁ 3 µm (**y₁⁺ 0.44–0.48**, 断熱 0.08 の ×5.7), run_0004 から壁温切替 (soft/mid/ランプ) + 48000 step | 残差床は A と同水準。**$C_f$/VD-II 0.97–0.99, $2St/C_f$ 1.16 (全ステーション), δ\* CONTUR 比 +3.5〜6 %, θ ±0 %, H 4.67 (CONTUR 4.47), エネルギー閉合 1.017, VD 対数則差 −0.3, T(u) は Walz 2.6 % / Duan–Martín 4.7 %**。series ALL STEADY。B/A 比: $C_f$ 1.50 (VD-II 比 1.35–1.39 → +8 %), δ\* 0.69–0.75 (CONTUR 0.77–0.80), θ 1.48 (CONTUR 1.43)。図 `run_0005_B_tw300_y3/cooled_plate_eval.png` | active (**冷却壁 主対象**) |
| `run_0006_C_tw700_y3` | C $T_w$ 700 K ($T_w/T_{aw}$ 0.60), 同 (y₁⁺ 0.15) | $C_f$/VD-II 0.91–0.95, $2St/C_f$ 1.15, δ\* CONTUR 比 +8〜16 %, θ ±3 %, H 6.9–7.3 (CONTUR 6.2), 閉合 1.018, T(u) Walz 1.8 % / DM 3.7 %。C/A 比: $C_f$ 1.20 (VD-II 1.15–1.16), δ\* 0.84–0.89 (CONTUR 0.885–0.893), θ 1.19 (1.18)。series ALL STEADY。図 `run_0006_C_tw700_y3/cooled_plate_eval.png` | active (中間温度) |
| `run_0007_B_tw300_y6` / `run_0008_B_tw300_y12` / `run_0009_B_tw300_y24` | B の y₁⁺ 掃引 (y₁ 6 / 12 / 24 µm = **y₁⁺ 0.9 / 1.8 / 3.6**), run_0005 から cross-mesh (`interp_field`) + soft/mid/ランプ + 36000 step | run_0005 (y₁⁺ 0.46) 比 [x 0.3–0.9 平均]: **y₁⁺ 0.9: $C_f$ −1.0 %, $q_w$ −1.3 %, δ\* −1.0 %, θ −1.0 % / y₁⁺ 1.8: −2.6 %, −3.7 %, −2.1 %, −2.0 % / y₁⁺ 3.6: −3.1 %, −7.0 % ($2St/C_f$ 1.12), −2.2 %, −2.2 %**。→ δ\*/θ (設計チェーンの量) は y₁⁺ ≤ 3.6 でも 2 % 内、$q_w$ (熱負荷) は y₁⁺ ≤ 1 で 1.5 % 内・≤ 2 で 4 %。残差床・series は全て run_0005 と同水準 (STEADY) | active (**y₁⁺ 許容の根拠**) |
| `run_0010_Aplain_ad_y3` / `run_0011_Bplain_tw300_y3` | 素 SST (`dilatationCorrection 0, katoLaunder 0`; SU2 比較の基準対), run_0004 / run_0005 から +36000 step | 生産 SST 比: A-plain $C_f$ +1.6 %, δ\* +1.0 %, θ +1.4 % / B-plain $C_f$ +1.0 %, δ\* +0.7 %, θ +0.8 % (平板では dilatation 補正の効きは 1〜2 %)。**B-plain vs SU2-B (積分量)**: CD 2.286e-3 vs 2.268e-3 (**+0.8 %**), HF 103.86 vs 104.68 kW/m (**−0.8 %**)。A-plain CD 1.547e-3 vs SU2-A 1.565e-3 (−1.2 %, SU2 まだ上昇中)。B-plain の $C_f$/VD-II 0.98–1.00, $2St/C_f$ 1.16。series: B STEADY | active (**SU2 基準対**) |
| `su2_A_ad` / `su2_B_tw300` | SU2 v8.5 RANS-SST (V2003m), 同一 `.su2` (y₁ 3 µm), 断熱 / `MARKER_ISOTHERMAL 300`, 超音速入口/出口, CFL 適応 ≤30, ~1 it/s | **中間 (it ≈3000)**: B の積分量は静定 (rms[ρ] −11): SU2 CD (=∫Cf dx) 2.268e-3 / HF (=∫q_w dx) 104.68 kW/m vs forge run_0005 (生産 SST) 2.264e-3 (**−0.2 %**) / 102.8 kW/m (**−1.8 %**)。A: SU2 CD 1.565e-3 (まだ +0.06 %/300 it で上昇中) vs forge run_0004 1.523e-3 (−2.7 %)。**ステーション比較 (it 5000 の restart, B-plain run_0011 ↔ SU2-B, 同じ評価関数)**: $C_f$ 1.000–1.001, $q_w$ 1.008, θ 0.998, **δ\* 1.037–1.047** (forge H 4.66 vs SU2 4.48; θ と Cf が一致して δ\* だけ +4 % = 壁近傍の密度分布差。T(u) の Walz 形からの差 forge 2.6 % / SU2 2.2 %)。**A-plain run_0010 ↔ SU2-A** (it 5000, SU2 は $T_w$ 1145 K でまだ上昇中): $C_f$ 0.997–1.000, θ 0.98–1.02, δ\* 0.96–1.04 | active (計算中) |
| `run_0012_B_tw300_y3_nx1500` | 流れ方向格子感度: nx 1000 → 1500 (`fp_y1_3um_nx1500`, SOFT-PASS), B 300 K, run_0005 から cross-mesh + soft/mid/ランプ + 36000 step | run_0005 比 $C_f$ +0.2 %, δ\* ±0.1 %, θ ±0.1 % → 1〜3 % 比較の前提を満たす | active (格子感度) |
