# 引き継ぎ ①: θ_r の反復の不確かさと、設計チェーンの許容差 (2026-10-10)

共通ルールは [`2026-10-10-handoff-common.md`](2026-10-10-handoff-common.md)。正本は `plans/active/time_integration-line-implicit-speed.md` の §5.1 #22 (元のセッションが反映する)。

## 問い

case/45 の NS の定常解は、反復の方式 (ライン / point) によって止まる値が違う。θ_r(40/70/94 r_t) で約 0.1 %、Q_w で約 0.01〜0.07 % で、どちらも残差は下がりきっていない (速度 plan §6.15・§6.20・§6.22)。
**設計チェーン (ノズル壁の設計・δ_r の補正・判定) が、θ_r・δ_r・Q_w にどれだけの許容差を要求しているか**を調べ、反復の不確かさ ±0.1 % がそれより小さいかを判定する。

## 手がかり

- 上位の指摘: `notes/reviews/2026-10-10-line-implicit-speed-results-diagnose.md` の「解釈への指摘 (a)」(W5 で 7e-6 r_t が合否を分けた件)。
- メモ: `design-venv-opt-interpreter`。system python3 の numpy 1.26 で δ_r が 7e-6 r_t 動き、W5 が誤って FAIL した教訓。
- 設計チェーン: `design/forge_design/` (feedback・metrics・evaluate)、`case/45.isobutane_m6_d155/` の問題 YAML (`problem_d155_ns_prod*.yaml`)。関連 plan は `plans/active/tooling-nozzle-isothermal-wall-chain.md`、`plans/active/verification-m6-axis-wave-mesh-su2.md`。
- θ_r の抽出は `case/45.isobutane_m6_d155/cold_xcheck.py` (帯は `common_yb`)。

## やること

1. 設計チェーンの中で、θ_r・δ_r・Q_w (と、それらから作る量) に掛かる許容差・判定の閾値を洗い出す (ファイル:行と値、相対か絶対か)。
2. 反復の不確かさ (θ_r で ±0.1 %、Q_w で ±0.07 %) を、それぞれの閾値と同じ単位に直して比べる。
3. 「ラインで止まった値を設計に使ってよいか」の結論を、根拠つきで書く。不確かさが閾値に近い・超える箇所があれば、それも挙げる。

- GPU は不要。コードも変えない。
- 結果は `notes/investigations/2026-10-xx-theta-iteration-tolerance.md` に書き、元のセッションに知らせる。
