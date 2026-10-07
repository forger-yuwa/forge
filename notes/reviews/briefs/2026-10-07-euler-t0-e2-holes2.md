# 諮問 (1 回で採否): E2 の登録の穴 2 つ — 近零の主指標の相対 STEADY と、本段 18000 step では still converging

日付 2026-10-07。諮問先 codex (diagnose)。エスカレーション条件 1。結果はまだ誰も見ていない (E2 は未投入。スクリプトは乾式確認済み)。ユーザから「時間がかかりすぎ」との声があるので、下の案の採否を 1 回で決めてほしい。
plan: `plans/active/verification-case45-euler-total-enthalpy.md` §6 E2 (commit f158bac2 の版)。前回の諮問 `notes/reviews/2026-10-07-euler-t0-e2-convergence-precondition-diagnose.md`。HEAD 548778bc 以降。

## 観測事実 (実装担当、`case/45.isobutane_m6_d155/euler_t0_e2_eval.py` 未 commit・乾式確認済み)

1. **相対 STEADY が近零の主指標で成り立たない**: `check_quasisteady.py:280–286` は drift と fluct を |平均| で割る。E2 の主指標は T₀ − 1600 (領域ごとの最大・最小・99 % 点) で、正常な領域では復元の床 (IC で max\|dev\| 0.0107 K、壁域の 99 % 点 0.0045 K) の近くに乗る。本物のツールで: 平均 0.07 K・幅 0.018 K → DRIFTING、平均 0.01 K・幅 0.0036 K → DRIFTING、平均 0.5 K・幅 0.016 K → STEADY。どれも幅 ≤ 0.1 K は満たす。P2 は 27 列 × 2 窓 × 両腕すべての STEADY を要求するので、B が仮説どおり正常でも不成立になる見込み。
2. **本段 18000 step では A が P1 を満たさない見込み**: run_0153 (腕 A と同じ壁・格子・IC の手順・段) の本段区間は `NOT CONVERGED (still converging)` (全列 falling、2.6〜2.7 dec)。P1 は「PASS か、不合格の理由が停滞だけ」。延長の run_0160 (通算 54000 step) は stalled/plateau になった (V5b)。Euler G1 は 2 本並列で数百 step/s。

## 案

- Q1 (近零): (a) `check_quasisteady` には T₀ − 1600 でなく T₀ そのもの (≈ 1600 K) を渡し、整定の判別は時間方向の幅 ≤ 0.1 K に任せる (STEADY は記録)。(b) \|平均\| が小さい列 (例 < 1 K) は幅の条件だけで判定。(c) 別案。
- Q2 (長さ): 本段を 54000 step に延ばし、評価窓を本段 42000〜54000 の 13 枚と末尾 5 枚にする (V5b で ISEN が plateau に達した長さ)。soft 3000 は同じ。固定の予算で、収束の予測ではない。延長の再延長はしない。(別案: P1 に still converging を許す)

## 問い

1. Q1・Q2 の採否と、plan に書く文言。
2. 実装担当の読み (結果を見る前、plan §9 に記録する): 評価点は列ごとに A と B の節点の η の和集合・区分線形の補間 / 「異常」はスロートの 3 領域の主指標の max(\|最大\|, \|最小\|) > 100 K、13 枚すべてで / 「差 ≤ 1 K」は最大・最小・99 % 点に当て割合には当てない / 境界は ≤、割合は狭義の > / 固定の条件の不成立は判別不能 / soft 段の後の時点の欠けは記録の欠け、評価窓の欠けは前提不成立。これでよいか。

読んでよいもの: 上記 plan と前回の諮問、`case/45.isobutane_m6_d155/{euler_t0_e2_eval.py, euler_t0_e2.py, run_euler_t0_e2.sh, test_euler_t0_e2_eval.py}` (未 commit)、`solver_density_cuda/tools/check_quasisteady.py`、`case/45.isobutane_m6_d155/README.md` (run_0153・run_0160 の行)。編集は禁止。
