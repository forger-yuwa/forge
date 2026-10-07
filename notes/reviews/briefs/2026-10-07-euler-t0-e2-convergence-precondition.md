# 諮問 (短い): E2 の収束の前提を「PASS」にするか「PASS か既知の plateau」にするか

日付 2026-10-07。諮問先 codex (diagnose)。エスカレーション条件 5 (諮問の Major の一部を採らないかの判断)。結果はまだ誰も見ていない (E2 は未実行)。
plan: `plans/active/verification-case45-euler-total-enthalpy.md` §6 E2 (改訂版、「判定の前提」の項)。前回の諮問 `notes/reviews/2026-10-07-euler-t0-e2-fix-diagnose.md` は「本段の `check_convergence` が両側 PASS」を前提にした。

## 観測事実

- case/45 の Euler の本段の `check_convergence` は、全温が正常な run も含めて、手元で確認できたものはすべて NOT CONVERGED (stalled/plateau): run_0062 (1100 × 65、全温 +0.073 K、最終 step 5999、rms_ro 0.3 dec)、run_0114 (G1、+241 K)、run_0143〜0145・0150〜0160 (G1、本段・延長とも 0.1〜0.9 dec)。PASS の Euler の run は見つかっていない。
- MOC の plan の V5 の前提は「pass か既知の plateau」だった。

## 問い

1. E2 の収束の前提を「PASS」のままにすると、両側とも plateau で判別不能になる見込みが高い。plan の改訂版のとおり「PASS か既知の plateau (RISING・DIVERGED なし)」に緩め、代わりに主指標の時系列の STEADY (全 13 枚と末尾 5 枚) で整定を見ることは妥当か。
2. 妥当でないなら、E2 の前に何をすべきか (収束させる設定の探索など)。

読んでよいもの: 上記 plan と前回の諮問、`solver_density_cuda/tools/{check_convergence.py, check_quasisteady.py}`。編集は禁止。
