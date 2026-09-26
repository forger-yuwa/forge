# 諮問ブリーフ: Cary 冷却平板の比較原点 A/B (T4-0a-0) が判定不能・収束ゲート未達。次に何をするか

- plan: `plans/active/case-hypersonic-gap-heating-validation.md` §4.12 (T4-0a-0)、§6 G14、§5.1 #61 の末尾
- 事前登録: `case/59.flatplate_cary_m6/acceptance.json` (結果を見る前に commit f1e56cac)
- 台帳: `case/59.flatplate_cary_m6/conditions.json` (Cary TN D-5863 Table II 転記)、`case/59.flatplate_cary_m6/README.md`
- 生成器: `case/59.flatplate_cary_m6/gen_runs.py`、`gen_mesh.py`、比較: `case/59.flatplate_cary_m6/tools/cary_compare.py`
- 前回の諮問: `notes/reviews/2026-09-26-t4-0-calibration-plates-diagnose.md`
- run はリモート (AWS) にあり、このリポジトリには無い。数値は下に全部書く。読まないこと: plan の §4.12・§6・#61 以外。

## 1. 観測事実

- `run_0001_tw02_re027`: Cary 系列 Re0.27_Tw0.2 (Tw 106.6 K、Re 2.7×10⁵/cm、p∞ 2231 Pa、T∞ 64.6 K、U∞ 970 m/s)。2D node、低 Re SST (Tu 0.5 %、μt/μ 10、Prₜ 0.9、dilatation 2、Kato-Launder)、熱量的完全空気、Sutherland、Pr 0.72。
  メッシュ y1 1.5 µm・平板 0.6 m・900 セル (`check_mesh_quality` PASS、AR 957)。段階起動 lam→soft→mid→2 次 0.5/1/2→本段 cfl 2 relax 0.7 × 20k、段の引き継ぎは `restart_field.py` 全段 OK。
  区間 ramp0→main: `NOT CONVERGED (still converging)`、全列 2.0–2.6 dec 低下中。
- `run_0002_tw02_re027_ext`: run_0001 の res_20000 から同一設定 +40k。延長区間単独: `NOT CONVERGED (stalled/plateau)`、rms_ro 1.44e-7→9.0e-8 (0.2 dec)、roe 0.2、roOmega 0.4、roK 2.8 dec。
  比較点 St の 34 系列 (20 枚): `check_quasisteady` ALL STEADY、drift 0.0 %。
  壁解像 `check_wall_resolution --target 1 --over-frac 0`: FAIL、y1+ 平均 0.541・p99 0.717・最大 3.94 (前縁 index 0)、>1 が面積 0.4 %。
- CFD の St は前縁で最大 (完全乱流)。x=1 cm 1.29e-3、5 cm 1.04e-3。
- A/B (事前登録どおり): A = x_CFD = x_exp、B = x_CFD = x_exp − 22.22 cm。

| x_exp cm | St_exp | R_A | R_B | D |
|---|---|---|---|---|
| 27.94 | 8.23e-4 | 0.982 | 1.236 | +0.258 |
| 29.85 | 7.84e-4 | 1.022 | 1.244 | +0.217 |
| 34.61 | 7.71e-4 | 1.017 | 1.181 | +0.161 |
| 38.41 | 6.74e-4 | 1.146 | 1.297 | +0.132 |
| 42.21 | 6.41e-4 | 1.189 | 1.324 | +0.114 |
| 46.04 | 5.62e-4 | 1.338 | 1.472 | +0.100 |
| 47.00 | 5.81e-4 | 1.290 | 1.417 | +0.099 |
| 49.20 | 6.00e-4 | 1.241 | 1.355 | +0.091 |

  (17 点全体: R_A 平均 1.128 [0.974, 1.338]、R_B 1.288 [1.125, 1.472]、D 平均 +0.146 [+0.091, +0.258]、15/17 点で D>10 %)
- 実測の St は x 27.9→49.2 cm で 8.2e-4→6.0e-4 (−27 %)。CFD は同区間 A で −8 %、B (x−22.22) で −20 %。

## 2. 期待値と出典

- 事前登録の読み: 全点 |D|>10 % → 原点感度は小さいを棄却 / 全点 <10 % → 原点の影響は小さい / またぐ → 判定不能。前提ゲートは収束 PASS・St 系列 STEADY・壁解像 PASS。

## 3–4. 実施済み

上のとおり。原因帰属はしていない。

## 5. 仮説 (未確認)

- 実測の下流の急な減少は遷移末端の加熱ピーク (overshoot) の後の緩和で、完全乱流の CFD はこれを持たない → 比較点を下流に寄せるほど R が上がる。
- 残差プラトーは float32 の丸め床か有界振動で、St は定常 (drift 0.0 %) — このリポジトリの平板 (case/48・case/56) でも同型のプラトーが出ている。
- 壁解像の超過は前縁の特異点付近だけ。

## 6. 聞きたいこと

1. **この A/B の結果から何を言ってよいか** (判定不能のまま、何を記録するか)。「D が 9–26 % で比較座標の選び方が結果を大きく動かす」と書いてよいか。
2. 収束ゲート (プラトー) と壁解像ゲート (前縁 0.4 %) を**どう満たすか / 満たせない場合にどう扱うか**。丸め床の確認 (倍精度の対照や `rounding_floor.py`) を先にやるべきか。
3. **次の一手を 1 つ**: 比較座標の問題を避けるには、遷移を再現する (forge には γ–Reθt 遷移モデル LM2009 が実装済み) か、比較点を「加熱ピークから十分下流」に限定するか、実測から仮想原点を当てはめるか、別の方法か。
