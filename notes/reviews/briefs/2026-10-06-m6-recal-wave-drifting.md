# 諮問: 最終 NS の波 η0.1 が値は閾値内だが check_quasisteady で DRIFTING (延長後も) — ③ 合格としてよいか、④ へ進めるか

日付 2026-10-06。codex (diagnose)。エスカレーション条件 3 (事前登録「③④ の未達は延長 1 回、なお未達なら諮問」)。
plan: `plans/active/tooling-nozzle-cfd-pinned-initial-line.md` §5.1 #11〜#11f、§9 末尾 (NS 連鎖の結果)。

## 観測事実

- ③ `case/45.isobutane_m6_d155/run_0116_ns_recal_final` (細分格子、Md_moc_offset +3.770e-4、k_f 1.054129、r_t 76.6539 mm、段階起動 → cfl 1・60000 step) と延長 `run_0117_ns_recal_final_ext` (restart_field ビット一致、cfl 1・60000)。
- 時系列 (`exitM_sampling_ab.py`、各 run 12 枚 5000〜60000、Euler 参照 run_0114):
  - run_0116: 出口コア M 5.998860 (幅 3.0e-5) STEADY; δ_E/δ_C 0.99985 STEADY; オーバーシュート η0.1 0.00752 % ± 0.00054 OSCILLATING; 波 η0.1 0.00660 % (幅 0.00049) DRIFTING (drift 5.2 %/tail・fluct 7.4 %)。
  - run_0117: 出口コア M **5.998887** (幅 9.2e-6) STEADY; δ_E/δ_C 0.99983 STEADY; オーバーシュート η0.1 **0.00794 %** (幅 0.00025) STEADY; 波 η0.1 **0.00651 %** (幅 0.00053) **DRIFTING** (drift 5.4 %/tail・fluct 8.2 %、"extremum at tail-end")。
- 登録ゲート: 出口コア M 6.000 ± 0.02 % (5.9988〜6.0012) → 合格 (余裕 0.00009); 波 η0.1 ≤ 0.01 %・オーバーシュート η0.1 ≤ +0.035 % (末尾 5 枚 STEADY); 出口半径 0.775 ± 0.1 mm (0.7749995) 合格; δ_E/δ_C 1 ± 0.5 % 合格; 壁解像 PASS (3.6 %)。check_convergence 両 run NOT CONVERGED (plateau)。
- 波 η0.1 の定義: `nozzle_report.metrics` — 試験窓 [x_E+2, x_F−1] の 100(M/6−1) から P-spline (knot 10) を引いた残差の最大絶対値。0 近傍ではない (0.0065 %) が、最大値なので位置が跳ぶと値が段で変わる。
- 過去の同指標: run_0109 (前の壁) 0.0066 % STEADY、run_0094 (粗格子) 0.0036 %。

## 問い

Q1. 値が閾値の 2/3 で幅 0.0005 %pt の DRIFTING を、登録ゲート「≤ 0.01 % (STEADY)」の未達とみなすべきか。みなすなら次の一手 (さらなる延長は登録外)。
Q2. ③ を合格扱いにして ④ (凝縮 ON、cfl 1・18000 step) へ進めてよいか。進める場合の記録の仕方 (DRIFTING を残す)。
Q3. 最大値型の指標 (波) の準定常判定として妥当な代替 (例: 最大値の位置の固定・上位 N 点平均) を将来の登録に入れるべきか。

## 読んでよいもの

上記 plan、`case/45.isobutane_m6_d155/{README.md,exitM_sampling_ab.py,run_recal_chain.sh}`、`design/forge_design/report/nozzle_report.py` (metrics)、`solver_density_cuda/tools/check_quasisteady.py` (classify)。run は AWS のみ。編集禁止。
