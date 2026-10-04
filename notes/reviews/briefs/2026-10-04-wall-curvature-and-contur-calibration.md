# 諮問ブリーフ: 壁の 2 階微分の波打ちの直し方と、CONTUR を CFD 場で較正する案

plan: `plans/active/verification-m6-axis-wave-mesh-su2.md` §5.1 #8・§9 (2026-10-04 末尾の各項)。報告ページ https://claude.ai/artifact/6bo7GNmcBAGhG4YnA3mgdU 。
読んでよいファイル: `design/forge_design/geometry/wall_axismach.py` (PhysicalNozzleWall L233–、特に L271–340)、`design/forge_design/evaluate/runner_axismach.py` (prepare_ns L671–800)、
`design/forge_design/feedback/deltastar_integral.py` (integral_bl, delta_r_function)、`design/forge_design/metrics/deltastar.py` (smooth_delta_quintic)、
`design/forge_design/feedback/deltastar_loop.py` (extract_and_merge)、`design/forge_design/meshing/mesh2d.py`、`case/45.isobutane_m6_d155/viz_band_ab.py`。h5 は読まない。

## ユーザの関心 (最優先)

「二回微分がガタガタしているのが一番気になる。どうにかしないといけない」。加えてユーザ提案「CFD で得られた場を使って CONTUR の算出式 (N など) を修正する方がきれいな壁ができるかもしれない」。

## 観測事実

- 物理壁 (`PhysicalNozzleWall`, offset=radial): 設計壁を 3000 点で評価 → δ_r(x) を **np.interp (直線補間)** で足す (δ_r の表は 1250 点、間隔 0.0904 r_t = 前回 NS の格子断面) →
  真のスロート探索 → スロート下流は**全点を通る補間 5 次 B-spline** (クランプ; ノット 2662、平均間隔 0.036 r_t)。上流は 5 次 Hermite。
- 2 階微分は数学的に連続 (B1 の壁でノット左右の r″・r‴ が一致を確認)。しかし r″ − r_inv″ は周期 ≈0.100 r_t (δ 表の間隔 0.090 と整合) で振動し、x=20 で p-p 5.8e-4 [1/r_t]。
  高周波成分の最大: x=0.5〜6 で約 2e-3、下流で 4〜11e-4 [1/r_t] (旧最終壁・B0・B1 の 3 枚とも同程度; 帯の選び方と無関係)。設計壁 |r″| は x=3〜90 で最大 1.4e-2。
- CFD 格子 (1250 断面、x 間隔 0.0904 = δ 表と同じ) の壁節点から取った離散曲率 − 設計壁の離散曲率: 高周波最大 1.6〜3.0e-4 (各区間)、なめらかな成分は 1.5〜2.4e-4 (x<6 で 1.1e-3)。
- 積分法 (CONTUR) の初期壁 run_0022 も、積分法の出力 1500 点を `delta_r_function` (np.interp) で載せていたので同じ鋸歯がある (前ページの図で最大 6e-3)。
- 帯修正の CFD A/B: 方式 E の壁 (B1, run_0046) で軸の局所の山谷は消えた (b70 0.36→0.013 %pt) が、P(0.1) 0.079 % が残り判定保留 (x≈48 の小山 +0.035 % と x 56→80 の −0.04 % の傾き)。
  B1 の場の C− で壁 x=13 → 軸 x=47.6。B1 は固定点 (δ_r/δ_in 1.0017)。新/旧 Euler 参照の差 ≤0.01 %。
- 抽出 (方式 E)/CONTUR の比: x≥5 で旧バイナリ 1.003 (x=5) → 1.000 (20) → 0.994 (60) → 0.983 (94)、新バイナリ 1.080 → 1.052 → 1.038 → 1.028、x<5 は 1.04〜1.29 (δ 0.001〜0.02 r_t)。なめらかで単調、局所の凸凹なし。

## 仮説

- H1: 2 階微分の波打ちは「δ_r の表の直線補間 → 折れ目を通す補間 5 次スプライン」による。δ_r をなめらかな関数 (P-spline の係数そのもの、または表を高次で補間) として壁に渡せば消える。
  ただし補間 5 次スプライン (全 3000 点を通す) 自体が、入力の微小ノイズを増幅する可能性もある。
- H2: B1 に残る P(0.1) 0.08 % の一部は、この壁曲率の波打ち (格子間隔と同じ波長) が流れに与える影響かもしれない (未確認)。
- H3: CONTUR を CFD 場で較正 (N の表・C_f の少数係数) すれば、壁は常微分方程式の解なので構造的になめらかになる。ただし較正パラメータの選び方次第で x 方向の形が合わない可能性。

## 聞きたいこと

1. 2 階微分の波打ちの直し方を 1 つ推薦してほしい (壁の作り方: δ_r を関数で渡す / 平滑化スプラインで壁を作る / 3000 点の評価をやめてノットを粗くする / その他)。
   既存の A14 (制約付き LSQ B-spline 壁、不採用; `plans/accepted/tooling-nozzle-axismach-physical-throat.md` に記録あり) との関係も。オフラインで合否を決める基準 (例: 壁 r″ の高周波最大、格子の離散曲率) を数値で。
2. その修正で B1 の P(0.1) の残りが減るかを確かめる価値はあるか (H2)。あるなら最小の A/B を 1 つ。
3. CONTUR 較正案 (H3) の方式: どの係数を動かすか (N(Re_δ) の表のずらし・C_f の倍率・その組合せ)、合わせる x 範囲、基準 (抽出 δ か軸 M か)。1 と組み合わせる順番。
4. 次の一手を 1 つ。
