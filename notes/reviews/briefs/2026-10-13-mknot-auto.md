# 諮問ブリーフ: M_knot を長さごとに自動で選ぶルール (plan tooling-nozzle-mknot-auto §4・§6)

AGENTS.md 条件 1 (新規 plan の §4・§6)。判断役のスイッチが codex。ユーザ指示 (2026-10-13)「m6 のやり方と同じように、knot の最適位置を定めるようなやり方をルールにしてくれるかな」。

## 観測事実 (plan §3)

- M6 (knot plan): 長さごとに壁での余裕 min(μ_w − θ_w) が M_K に対して単峰。最短ロバストの探索の制約は hard・単峰・余裕 ≥ 1°・網の最小 sin 角 ≥ 0.02。CFD では M_K 2〜3 で軸 M のずれ 0.040 % (鈍感)。
- M5 デモ (case/70、長さ 1.2 m を from_length で固定): 余裕 2.0: 4.19°・2.25: 4.32°・2.5: 4.37°・2.75: 4.30°・3.0: 4.14°・3.5: 3.64°; 1.6・1.8 は壁 QA (スロートの曲率の不整合) で不成立。
  曲率・最大壁角は M_K が小さいほど小さい。L_c は全点 21.49 で同じ。
- 今の要求カードは M_knot = 0.45 × M_design (M6 の最短ロバスト点からの外挿)。

## 方針 (plan §4)

- 目的: 余裕を最大。制約: 設計チェーンが通る・網の反転 0・交差 0・重なり 0・最小 sin 角 ≥ 0.02・余裕 ≥ 1°。成立点なしは例外。
- 探索: (M_A + 0.1, M_d − 0.1) を 0.1 刻み → 最大の格子点の ±0.1 を 0.01 刻み。極大 2 つ以上は警告して大域の最大。差 1e-3° 以内は小さい M_K。
- 物理スロート径のカード: 最初の r_t で選んで固定点の間は固定、収束後に詰めをやり直し、0.01 を超えて動けば固定点をやり直す (最大 2 回)。
- カードの既定を自動にする。数値の override は今と同じ。記録を YAML と card_log に。

## 問い

**このルールのうち、実装後に「選んだ M_knot が本当は最適でない、または設計ごとに不安定に飛ぶ」形で失敗する可能性が最も高い点を 1 つ挙げ、防ぐ変更を 1 つ推奨してほしい。**
候補 (縛られなくてよい): (a) 余裕の山が平ら (M5 で 2.0〜2.75 が 4.19〜4.37°) で、離散化 (逆 MOC の解像度) の雑音で最大の位置が飛ぶ、(b) from_length の L_c の解き (許容 0.02 r_t) が M_K ごとに別の解に落ち、余裕に段差が入る、
(c) 網の検査の seam の除外の設定を M6 用からそのまま使う、(d) 余裕が「作れるか」の指標で流れの質と結びつかない、(e) r_t の固定点との入れ子。

## 禁止事項 (厳守)

- ファイルを変更しない (read-only サンドボックスで動いている)。
- **`*.log`, `residual_history*.csv`, `res_*.h5`, `*.vtu`, `*.msh`, `plans/README.md` を読まない** (巨大)。
- 下に列挙した `sed -n 'A,Bp' <file>` 以外のファイル読みをしない。grep は可。
- 両論併記で逃げず、**推奨は 1 つに絞る**。根拠は `ファイル:行` か提示した数値で示す。

## 読んでよいもの

- `sed -n '1,140p' plans/active/tooling-nozzle-mknot-auto.md`
- `sed -n '1,80p' case/70.m5_air_wt_demo/mknot_sweep.py`
- `cat case/70.m5_air_wt_demo/_mknot_sweep/mknot_sweep_L1200.json`
- `sed -n '50,120p' case/42.isobutane_wt/sweep_axislaw_A_MK.py`
- `sed -n '53,95p' case/42.isobutane_wt/compare_axislaw_D.py`
- `sed -n '440,470p;530,560p' design/forge_design/evaluate/runner_axismach.py`
- `sed -n '451,465p' design/forge_design/geometry/moc_diagnostics.py`
