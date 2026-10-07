# 諮問ブリーフ: Euler の全温の異常の plan の結果の解釈 (result 段の前)

- 日付: 2026-10-07
- plan: `plans/active/verification-case45-euler-total-enthalpy.md` (§9 の末尾「§8 の整理」が確定したい解釈)
- エスカレーション条件: 7 (result 段の解釈を確定する前)

## 1. 観測事実 (すべて plan §9 に run パスと出力つきで記録済み)

- E1: 等エントロピーの IC (A) は全節点 \|T₀ − 1600\| ≤ 0.011 K。soft 段の後 (B、run_0153 の nozzle.h5) で 1245〜1915 K。本段の res_0 (C) と B の差は roe・roY* の最大 1 ULP、他はビット一致 (`_band_ab/euler_t0_e1_BC_diff.json`)。
- E2: 同じ壁・同じ IC・同じ step で配点だけを変えた A (G1、壁際細分化) と B (全域 0.005)。B は全領域 ≤ 0.103 K・時間の幅の前提成立、A は数百 K の超過が窓の中でも動き前提不成立 → 登録上「判別不能」。
- E3 (`mesh_euler` の実装)・E4 段 1 (判別不能)・E4V (合格、較正値 6.8825e-6)。
- E5: NS の 3 条件の最終場を固定し Euler 参照の窓 13 枚を差し替えた δ_E/δ_C の幅 1.6e-5〜4.1e-5 (≤ 5e-4 で PASS)。
- E6: run_0062 (初期線の元、1100 × 65・全域 0.005) の res_6000 は T₀ − 1600 が −0.058〜+0.073 K。初期線の量の 5 枚の安定性は既存記録で 1e-4 以下。
- NS (上流の多項式化の plan §9): 新しい較正で NS の出口コア M 5.9985 (旧 5.99887)。ゲートはユーザ決定で ± 0.05 % に。

## 2. 確定したい解釈 (plan §9「§8 の整理」の ①〜③)

- ① 発生する段は「G1 格子の soft 段の計算の中」。原因の候補は、薄い高 AR セルでの離散化か float32 の幾何の精度、または G1 での整定しない収束。どちらかは未分離。
- ② 初期線 run_0062 は影響なし (E6)。出口較正の変化 −3.70e-4 は格子の違い全体の結果で、全温の異常だけの寄与とは扱わない。単調壁の E′ は取り消さない (ユーザ決定)。MOC の V5/V5b は V5d で置き換え。
- ③ 出口較正の採用 (E4V) と δ_E・C2 の参照の採用 (E5) は別のゲート。run_0164・run_0174 は両方を満たす。

## 3. 諮りたいこと

1. ①〜③ の書き方に、証拠を超えた主張・見落とし・誤りはないか (特に ① の「発生する段」の範囲と、② の run_0062 の扱い)。
2. この状態で `status: done` にして accepted へ移してよいか (完了の区分は「原因の候補と未確定の範囲での完了」。§8 はそれを認めている)。足りないものがあれば挙げてほしい。

## 読んでよいファイル

- `plans/active/verification-case45-euler-total-enthalpy.md`、`case/45.isobutane_m6_d155/_band_ab/euler_t0_e1_BC_diff.json`・`euler_ref_de_sensitivity.json`・`euler_t0_e2_eval.json`・`e4_recal_eval.json`、`case/45.isobutane_m6_d155/README.md`、`plans/active/tooling-nozzle-upstream-poly-and-throat-sizing.md` §9
