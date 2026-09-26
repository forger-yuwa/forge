# 諮問ブリーフ: Cary 冷却平板で「比較座標」をどう固定するか (T4-0a の次の一手)

- plan: `plans/active/case-hypersonic-gap-heating-validation.md` §4.12 (T4-0a、T4-0a-0)、§6 G14、§5.1 #61 の末尾 (2026-09-26 の精度 A/B まで)
- 事前登録と結果: `case/59.flatplate_cary_m6/acceptance.json` (T4-0a-0 の outcome、T4-0a-P の outcome)、`case/59.flatplate_cary_m6/README.md`
- 実験データ: `case/59.flatplate_cary_m6/conditions.json` (Cary TN D-5863 Table II、9 系列 × 31 点)、Table IV (遷移開始 x_b・終端 x_tr) は `papers/gap_heating/Cary_1970_NASA-TN-D-5863_turbulent_heat_transfer_cooling_M6.pdf` p.28 相当
- 遷移モデル: `procedures/recommended-settings.md` の LM2009 の節、`plans/*/` の遷移モデル plan (γ–Reθt、T3A で SU2 と一致、翼では入口粘性比次第)
- 前回の諮問: `notes/reviews/2026-09-26-cary-origin-ab-result-diagnose.md`
- run は AWS。数値は下に書く。読まないこと: plan の §4.12・§6・#61 以外。

## 1. 観測事実 (前回からの差分)

- 精度 A/B: 同じ初期場から float32 (`run_0003_prec_f32`) は plateau、全域 FP64 (`run_0004_prec_f64`) は `PASS (converged)` 5.0–6.2 dec。
  壁 q_w の差は相対 5e-5 以下、R_A/R_B/D は小数 3 桁まで同じ。**残差プラトーは丸め床で報告量に効かない**。
- T4-0a-0 (比較原点 A/B) は事前登録どおり判定不能のまま: D = St_CFD(x−22.22)/St_CFD(x) − 1 = +0.091〜+0.258 (15/17 点で >10 %)。
  R_A (x_CFD = x_exp) 平均 1.128 [0.974, 1.338]、R_B (x_CFD = x_exp − 22.22 cm) 1.288 [1.125, 1.472]。
- 実測 St は 27.9→49.2 cm で −27 %、CFD は A で −8 %、B で −20 %。実測の加熱ピーク (遷移終端) はこの系列で 22.22 cm、Table IV の x_tr は 18.8–21.6 cm (Tw/Tt 0.19–0.24)。
- 壁解像: y1+>1 は前縁 x ≤ 0.72 mm の 4 列のみ (最大 3.94 @ x=0、長さ重み 0.16 %)、比較点 28–49 cm は 0.48–0.50。
  **並行して y1 0.35 µm のメッシュで FP64 の再計算を始める** (G14 の前提ゲートを満たすため。前縁除外は事後導入しない)。
- CFD は前縁から完全乱流 (SST、Tu 0.5 %、μt/μ 10)、St は前縁で最大。

## 2. 期待値と出典

- Cary 原報 p.13: 原点の選び方だけで予測が 7–18 % 変わる。Cary は仮想原点 = 加熱ピーク、Re_v > 10⁶ で比較。
- 事前登録: 比較座標・比較点・共通 Re_v 範囲・重みは G14 本体の前に `acceptance.json` に固定する。結果に合わせて選ばない。

## 5. 候補 (どれも未実施)

- (a) **遷移を再現する**: LM2009 (γ–Reθt) を ON にして、実測の遷移位置 (Table IV の x_b, x_tr) と加熱ピークを再現できるかを先に確かめ、合えば同じ x で比較する。
  リスク: Tw/Tt 0.2 の極冷壁・M6 での LM2009 の実績なし、dilatation 補正との組合せ未検証、入口乱流量の依存。
- (b) **乱流を実測の遷移終端から始める**: 前縁から x_tr まで層流 (SST の生成を切る / 乱流量ゼロ) にして、そこから完全乱流。x_CFD = x_exp で比較。forge で区間ごとに乱流を切れるかは未確認。
- (c) **遷移の影響が消える下流だけで比べる**: 事前に「x_tr から δ の何倍」などの物理的基準で比較点を決める。ただし Cary の下流端 49.2 cm までで遷移の影響が消えている保証はない。
- (d) **運動量厚さ (Re_θ) で揃える**: 位置でなく局所 Re_θ を比較座標にする (Hopkins らの推奨)。実測に θ は無いので、実測側の Re_θ をどう得るかが問題。

## 6. 聞きたいこと

1. **比較座標の固定方法を 1 つ選んでほしい** (上の a–d か別案)。事前登録に書く形 (比較点・座標・許容の読み方) まで。
2. その方法で G14 (壁温比 0.2–0.7 の 9 系列の R) に進めるか、それとも先に 1 系列で方法そのものの妥当性を確かめる段が要るか。
3. 以後の Cary run を FP64 で回すべきか (float32 との報告量の差は 5e-5)。
