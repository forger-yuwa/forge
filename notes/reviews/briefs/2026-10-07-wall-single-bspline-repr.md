# 諮問: 物理壁を全域 1 本の x の 5 次 B-spline で「表現だけ」作り直し、STEP で渡す方針 (案 A) の §4 と §6

日付 2026-10-07。諮問先 codex (diagnose)。エスカレーション条件 1 (plan §4・§6 を書き直した)。
plan: `plans/active/tooling-nozzle-wall-single-bspline.md` (全文を読むこと。§9 に案 B を取り下げた経緯)。作業ツリー `/home/sano/work/forge-integ-1005`。
別の実装担当が同じツリーで `design/forge_design/geometry/moc_*.py`・`evaluate/runner_axismach.py` を編集中 (MOC の plan)。読むのは自由だが、未 commit の差分が混ざっている。HEAD (934d3086 以降) の版を基準に読むこと。

## 経緯 (前回の諮問・レビューとの関係)

- 前回の諮問 (`notes/reviews/2026-10-07-wall-single-bspline-diagnose.md`) と plan 段レビュー (`notes/reviews/2026-10-07-tooling-nozzle-wall-single-bspline-plan.md`) は、**案 B** (δ_r の平滑化そのものを x の B-spline に替える) に対するもの。plan 段レビューの独立試算で、事前に決めた条件を満たす λ が 242 候補中 0 だったので、案 B を取り下げた。
- ユーザ決定「A でいいかな」「STEP でいいよ、それで進めて」: 設計の中身 (MOC・当てはめ・δ_r の平滑化・ランプ) は変えず、今の物理壁を 1 本の B-spline に作り直す。CAD には STEP で渡す。

## 観測事実

- **今の物理壁** (r_t 単位、run_0147、`PhysicalNozzleWall` の解析経路 `r = design.r + s·δ_r`、`wall_axismach.py:497-500`):
  - 直管 6.485 ([−12.5, −12))、H(x) ([−12, −11))、H + s·δ_r ([−11, −6))、H + δ_r ([−6, 0))、S(x) + δ_r(x) ([0, 95.245])。
  - H は 5 次の多項式 (`UpstreamThroatPoly`)、s は 5 次 smoothstep、S は 5 次 B-spline (内部ノット 253)。
- **δ_r は x の 5 次 B-spline** (初稿で「x の区分多項式ではない」と書いたのは誤りで、訂正済み):
  - 第 1 段 `smooth_delta_quintic` (ξ = √(x + 16.5) の P-spline) の値を 1500 点 (間隔 0.0719 r_t) で取り、第 2 段 `delta_r_from_table` (`runner_axismach.py:909`) が `make_interp_spline(x, d, k=5)` で通す。内部ノットは表の点 x[3:−3] (1494 個)。
  - 壁が使うのは第 2 段 (`integral_delta_r` の joint 壁の分岐)。
- **事前試算** (`case/45.isobutane_m6_d155/wall_single_bspline_probe.py` → `_band_ab/wall_single_bspline_probe.json`、CFD 0 step):
  - 再構成した物理壁と run_0147 の `wall_physical.csv` の差: 1.5e-15 m。
  - ノット: 両端重複度 6、継ぎ目 (−12・−11・−6・0) 重複度 3、[−11, 0) に δ_r の内部ノット、(0, x_e) に S と δ_r の内部ノットの和集合。係数 1747、異なるノット 1735、最小のノット間隔 7.6e-5 r_t。
  - 係数は全ノット区間の Gauss 8 点での最小二乗 (scipy の `lstsq`、密行列)。所要約 30 s。
  - 元の壁との差の最大 (各ノット区間の内部 39 点):

    | 区間 | 半径 | 1 階 | 2 階 | 3 階 |
    | --- | --- | --- | --- | --- |
    | 直管 | 2.4e-14 | 1.7e-13 | 2.6e-12 | 5.5e-11 |
    | H | 2.4e-14 | 8.0e-14 | 9.8e-13 | 6.5e-12 |
    | ランプ (当てはめ) | 7.5e-14 | 4.6e-12 | 3.7e-10 | 7.5e-8 |
    | H + δ_r | 2.0e-14 | 1.4e-12 | 2.7e-10 | 2.1e-8 |
    | S + δ_r | 1.1e-13 | 3.5e-12 | 1.2e-9 | 4.1e-7 |

  - 継ぎ目の跳び: 値・1〜2 階微分 ≤ 1.3e-12。3 階微分は 0.065 (−12)・0.003 (−11)・0.0035 (−6)・0.185 (0) で、元の壁の跳びと同じ設計上のもの。
- **メッシュ**: `meshing/mesh2d.py:147` が `wall.r(xs)` で壁節点を置く。座標は float32 (`geom_float`)。r ≈ 0.775 m で 1 ulp ≈ 6e-8 m。スロート付近の壁節点間隔は 0.0155 r_t (run_0117 の `nozzle.h5`、run_0147 も同じメッシュ設定)。
- **STEP**: FreeCAD 1.1.1 (OpenCascade) の `Part.BSplineCurve.buildFromPolesMultsKnots` → `exportStep` で 5 次の `B_SPLINE_CURVE_WITH_KNOTS` が書けることを確かめた (小さな例で)。

## 方針 (plan §4)

- §4.1: 上の事前試算の方法で 1 本にする。許容誤差 (半径 ≤ 1.3e-7 r_t、1 階 ≤ 1e-7、2 階 ≤ 1e-5) を超えたら止める。
- §4.2: 新しい壁クラスは `PhysicalNozzleWall` の下流向け属性をすべて持つ。報告は保存した係数から評価する。
- §4.3: STEP は平面の B-spline 曲線 (x(u) = u をグレビル点の制御点で厳密に表す)、mm。CAD の形と CFD の形の関係 (CFD は節点を直線でつないだ多角形で、節点間で最大約 1.2 µm 内側) を添え書きに書く。
- §4.4: `geometry.physical_wall_repr` (`legacy` | `single_bspline`、既定 `legacy` でビット同一)。
- §6: W0 (既定のビット同一)、W1 (精度)、W2 (スロート量)、W3 (メッシュ座標の差 ≤ 2 ulp なら CFD のやり直し不要)、W4 (下流の道具)、W5 (STEP の読み直しと独立評価器)。

## 問い

1. 厳密化の方法 (ノットの和集合 + Gauss 点の最小二乗) は妥当か。近接ノット (最小間隔 7.6e-5) や係数 1747 個で、数値上・CAD 上の問題は出ないか。最小二乗でなくノット挿入で厳密に組むべき理由はあるか。
2. 事前試算の結果は許容誤差より 5 桁以上小さい。許容誤差は事前登録として据え置いたが、据え置きでよいか (結果を見て締めるのは「結果を見てから合格条件を作る」に当たるか)。
3. W3 (メッシュ座標の差 ≤ 2 float32 ulp なら NS・凝縮をやり直さない) は、CFD を回さない根拠として十分か。壁の表現が CFD の入力に入る経路は、メッシュの節点座標以外にあるか (`wall_dist`・初期値 `ic.py`・δ* の抽出・境界条件 など)。
4. §4.2 の壁クラスの置き換えで、見落としやすい下流の依存は何か。スロート量を 1 本の B-spline から求め直すことの影響。
5. STEP: 制御点約 1750 個・最小ノット間隔 0.0059 mm の曲線を CAD に渡すことに実務上の問題はあるか。W5 の検査 (OCC で読み直し、独立の de Boor 評価器、位置 ≤ 1e-6 mm 等) で足りるか。グレビル点の x 制御点で x(u) = u が浮動小数点でどこまで再現されるか。
6. §4.3 に書いた「CAD の形と CFD の形の関係」(弦の内側へのずれ (節点間隔)²·曲率/8 ≈ 1.2 µm、float32 の丸め約 0.06 µm、どちらも加工公差より 1 桁以上小さい) は正しいか。
7. MOC の軸処理の plan (`discretization-moc-axis-limit-and-corrector.md`) が生産に入ると設計壁 S が µm 級で変わる。本 plan は「その時点の生産壁を 1 本にする」道具なので独立だと考えているが、順序や交絡の問題はあるか。

## 読んでよいもの

- 上記 plan、`plans/active/discretization-moc-axis-limit-and-corrector.md`、`plans/accepted/tooling-nozzle-throat-monotone-r2.md`
- 前回の諮問・レビュー: `notes/reviews/2026-10-07-wall-single-bspline-diagnose.md`、`notes/reviews/2026-10-07-tooling-nozzle-wall-single-bspline-plan.md`
- `design/forge_design/geometry/wall_axismach.py` (`PhysicalNozzleWall`・`JointFitCFDWall`)、`design/forge_design/evaluate/runner_axismach.py` (`delta_r_from_table`・`integral_delta_r`・`prepare_ns`)、`design/forge_design/evaluate/ic.py`、`design/forge_design/meshing/mesh2d.py`、`design/forge_design/report/nozzle_report.py`
- `case/45.isobutane_m6_d155/wall_single_bspline_probe.py` と `_band_ab/wall_single_bspline_probe.json`
- run_0147 の保存物は別の作業ツリー `/home/sano/work/forge/case/45.isobutane_m6_d155/run_0147_ns_mono_final/` (`delta_r_initial.csv`・`wall_physical.csv`・`prepare_info.json`)

編集は禁止。
