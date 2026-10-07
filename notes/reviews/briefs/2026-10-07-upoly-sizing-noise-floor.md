# 諮問: 寸法の逆算 (U2・U2b) が CFD 前の δ_r の雑音の床で判定できない — 許容差・δ_r の経路・`solve_rt` の扱い

日付 2026-10-07。諮問先 codex (diagnose)。エスカレーション条件 3 (事前登録の比較が FAIL) と 4 (plan に無い修正へ進む前)。
plan: `plans/active/tooling-nozzle-upstream-poly-and-throat-sizing.md` (§4.2・§6 U2・U2b・§9 の 2026-10-07 の記録)。作業ツリー `/home/sano/work/forge-integ-1005` (HEAD 2c09d681)。

## 観測事実

- 実装 commit a382f37d。`design/forge_design/feedback/deltastar_loop.py` の `solve_rt_throat`・`solve_rt` (CFD 前の経路を `build_physical_wall` と共通の δ_r 経路に揃えた)、`SIZING_TOL_M = 1e-9`、`SizingNotConverged`。検証の出力 `case/45.isobutane_m6_d155/_band_ab/upoly/{U2.json, U2b.json, u2_solve.json}`。
- CFD 前の δ_r は `integral_bl` (CONTUR の積分法、RK45 rtol 1e-6) → 5 次 P-spline の平滑化 → `delta_r_from_table`。r_t を 1e-12 m 動かすと δ_r(x_e) が約 1e-5 r_t、δ_r(0) が約 1e-7 r_t 跳ぶ (r_t に対して滑らかでない)。
- 残差の雑音の床 (トレンドを除いた散らばり): スロート半径で 0.7〜1.4e-8 m、出口半径で 1.1〜2.6e-6 m。事前登録の合格条件は \|差\| ≤ 1e-9 m。
- U2 (スロート径から): B (共通経路) は 3 目標とも往復の差 0 / −3.6e-10 / −4.1e-12 m で形式上 PASS だが、床の中で 6〜7 回目に偶然入っただけ。A (`integral_bl` 直接) は生の δ_r で poly のゲートを通れず (\|Q″ − H″\| 4.9e-2 > 5e-3)、想定の「往復の差で外れる」形ではなく失敗。NS 後の経路 (δ_E の全分布) は 3 回で ≤ 4.4e-10 m に収束して PASS。
- U2b (出口径から、CFD 前): 3 目標とも 30 回で収まらず `SizingNotConverged` (最後の残差 7.2e-7 / −5.5e-7 / 9.5e-7 m)。**今のコードでは CFD 前の `solve_rt` は常に例外**になり、手順書 (`procedures/nozzle-design-workflow.md` ④) の CFD 前の寸法の決め方が止まる。
- 実寸の目安: case/45 の出口半径 0.775 m、スロート半径 0.0767 m。加工公差は数十 µm と見込んでいるが、指定値は未確認。メッシュの座標は float32 (出口で 1 ulp 約 0.06 µm)。
- 生産の寸法は NS 後の経路 (C2 方式、抽出した δ_E を Re^−0.2 で換算) で決めている。CFD 前の寸法は初回の見積もり。

## 問い

1. 許容差 1e-9 m は何を守るためのものか (数値解法の閉じ具合)。CFD 前の経路で、許容差を雑音の床に合わせる (例: スロート 1e-7 m・出口 1e-5 m) ことは「結果を見て合格条件を緩める」に当たるか。当たるなら、どう登録し直すのが筋か (雑音の床を測ってから許容差を決める手順を、新しい判定として登録する、など)。
2. δ_r を r_t に対して滑らかにする手 (RK45 の rtol を締める・固定刻みの積分にする・r_t 依存を解析的に分離する (Re のスケーリング) など) のうち、どれが妥当か。生産の δ_r (NS の前の初期壁) を変えることになるので、その影響と検証の条件は何か。
3. CFD 前の `solve_rt` が常に例外になる今の状態をどう扱うか (すぐに直す・旧経路に戻す・許容差を変える・警告にする)。
4. U2 の判別 A/B (A が「壁を作れない」で外れた) の結論はどう読むか。共通経路の採用は支持されるか。

## 読んでよいもの

- 上記 plan、`notes/reviews/2026-10-07-upstream-poly-throat-sizing-diagnose.md`、`notes/reviews/2026-10-07-tooling-nozzle-upstream-poly-and-throat-sizing-plan.md`
- `design/forge_design/feedback/{deltastar_loop.py, deltastar_integral.py}`、`design/forge_design/evaluate/runner_axismach.py` (`integral_delta_r`・`build_physical_wall`)、`design/forge_design/metrics/deltastar.py` (`smooth_delta_quintic`)
- `case/45.isobutane_m6_d155/{upoly_verify.py, upoly_solve_rt_record.py}`、`_band_ab/upoly/` の JSON
- `procedures/nozzle-design-workflow.md`

編集は禁止。
