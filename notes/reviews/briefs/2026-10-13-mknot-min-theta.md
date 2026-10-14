# 諮問ブリーフ: M_knot の自動選択を最大壁角の最小化で (plan tooling-nozzle-mknot-auto §4・§6 の書き直し)

AGENTS.md 条件 1 (plan §4 の方針を変える)。判断役のスイッチが codex。plan 段のレビュー (`notes/reviews/2026-10-13-tooling-nozzle-mknot-auto-plan.md`、GO-with-changes M5/m2) の M1 で
「M6 の 2.6・2.7 は余裕最大でなく最大壁角が最小の点」と訂正され、ユーザは目的を「最大壁角の最小化 (M6 と同じ)」に決めた (2026-10-13)。

## 観測事実 (plan §3)

- M5 デモ (case/70、ns2、r_t 固定、n_axis_inv 2400): θ_max は 1.98: 15.2994°・2.00: 15.2942°・2.05: 15.2886°・2.08: 15.2900°・2.10: 15.2934°・2.25: 15.363°・3.0: 16.47°。
  ≤ 1.97 は壁 QA (スロートの曲率の不整合 ≥ 15 %) で不成立。θ_max の最小 2.05〜2.08 は成り立たなくなる境目 (約 1.975〜1.98) から 0.07〜0.1。
- 余裕最大の規則での解像度の A/B (2.40〜2.60、N 2400/4800) は選択値が変わらなかった (山の頂上が平ら)。θ_max の最小化での解像度の振る舞いは未測定。

## 書き直した方針 (plan §4)

- 目的: θ_max を最小。制約: 設計チェーン・網の反転/交差/重なり 0・最小 sin 角 ≥ 0.02・余裕 ≥ 1°・軸則の M″ の符号反転 ≤ 1。
- **境目からの距離 (新規)**: 選ぶ点は格子の上で不成立の最も近い点から 0.05 以上。θ_max の最小がそれより近ければ境目から 0.05 の点を選ぶ。
- 同点 (最小との差 1e-3° 以内) は境目から最も遠い点。全局所極小を詰める。最終の (r_t, M_K, 実効 geometry) の組で解像度のゲート (N と 2N で選択値の変化 ≤ 0.01)。

## 問い

**この書き直しのうち、実装後に「選んだ M_knot が境目の近くで不安定になる、または M6 の規則を再現しない」形で失敗する可能性が最も高い点を 1 つ挙げ、防ぐ変更を 1 つ推奨してほしい。**
候補 (縛られなくてよい): (a) ΔM_safe 0.05 の根拠 (境目の位置自体が解像度・r_t で動く量を測っていない)、(b) 境目を「格子の上の不成立の点」で定義すると格子の刻みに依る、
(c) 同点の規則 (境目から最も遠い) が θ_max の平らな谷で M6 の 2.6/2.7 からずれる、(d) θ_max は `wall_qa` の壁点列の最大で、M″ の制約・余裕と同じ解像度の問題を持つ。

## 禁止事項 (厳守)

- ファイルを変更しない (read-only サンドボックスで動いている)。
- **`*.log`, `residual_history*.csv`, `res_*.h5`, `*.vtu`, `*.msh`, `plans/README.md` を読まない** (巨大)。
- 下に列挙した `sed -n 'A,Bp' <file>` 以外のファイル読みをしない。grep は可。
- 両論併記で逃げず、**推奨は 1 つに絞る**。根拠は `ファイル:行` か提示した数値で示す。

## 読んでよいもの

- `sed -n '1,170p' plans/active/tooling-nozzle-mknot-auto.md`
- `cat case/70.m5_air_wt_demo/_mknot_sweep/low_end_N2400.json case/70.m5_air_wt_demo/_mknot_sweep/low_end2_N2400.json`
- `sed -n '85,125p' case/42.isobutane_wt/sweep_axislaw_A_MK.py`
- `grep -n "robust_min_theta\|robust_max_margin" case/42.isobutane_wt/sweep_axislaw_A_MK_summary.csv | head -20`
- `grep -n "kappa0_R_ratio\|0.15\|スロート曲率不整合" design/forge_design/geometry/*.py | head -10`
