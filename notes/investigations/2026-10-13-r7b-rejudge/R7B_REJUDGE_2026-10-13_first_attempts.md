# R7b ④-0 既存 3 本の再判定 (rejudge-r7b-v1) — 機械的な表示 (解釈は呼び出し側)

- 生成 2026-10-13T03:06:03+0000、git 4829694052e1bedb53db25228b619701434535ca、現行の受理検査の版 `sern-acceptance-v1`
- **結果 B: 欠落・判定不能・不合格の run がある (版の付け替えだけで救済しない)**
  - `run_1088_tewake_m6_L10`: 現行の受理検査で FAIL (JUDGMENT_UNSETTLED): 判定区間: 窓・定常性が未達: C_L: a 1.37e-05・d 9.08e-05 (許容 5e-05)・STEADY; 判定区間: 窓・定常性が未達: C_M: a 2.90e-04・d 1.71e-03 (許容 5e-04)・STEADY / 判定区間 (#2・#3) が PASS でない (UNSETTLED)
  - `run_1097_tewake_lsw10_m6_cont20k`: 現行の受理検査で FAIL (JUDGMENT_UNSETTLED): 判定区間: 窓・定常性が未達: C_T_with_shear: a 2.64e-06・d 6.85e-05 (許容 5e-05)・STEADY; 判定区間: 窓・定常性が未達: C_M: a 5.47e-04・d 1.97e-03 (許容 5e-04)・STEADY; 判定区間: 最低温度が未定常: ALL の最低温度 TRANSIENT-UNSETTLED / 判定区間 (#2・#3) が PASS でない (UNSETTLED)
  - `run_1085_tewake_g4_L10_cont20k`: 現行の受理検査で FAIL (JUDGMENT_UNSETTLED): 判定区間: 窓・定常性が未達: C_L: a 2.43e-06・d 5.87e-05 (許容 5e-05)・STEADY; 判定区間: 窓・定常性が未達: C_M: a 6.35e-05・d 1.72e-03 (許容 5e-04)・STEADY / 判定区間 (#2・#3) が PASS でない (UNSETTLED)

## 対象の表

判定区間の列は #1 床事象 / #2 温度 / #3 力の窓。力は判定区間の末尾 W step の「平均 (a / d)」(a = max|q − 平均|、d = 前窓との差)。

| run | L_sw | 作動点 | 格子 | 状態 | 現行 verdict | fail_class | 判定区間 #1/#2/#3 | C_T | C_T_with_shear | C_L | C_M | 欠けた入力・判定不能の理由 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `run_1088_tewake_m6_L10` | 0.8 | m6_on | g3 | JUDGED | FAIL | JUDGMENT_UNSETTLED | PASS/PASS/UNSETTLED | 0.8891113 (4.3e-06 / 7.3e-06) | 0.8704129 (4.2e-06 / 1.2e-06) | 0.06209784 (1.4e-05 / 9.1e-05) | -1.626624 (2.9e-04 / 1.7e-03) | — |
| `run_1097_tewake_lsw10_m6_cont20k` | 1.0 | m6_on | g3 | JUDGED | FAIL | JUDGMENT_UNSETTLED | PASS/UNSETTLED/UNSETTLED | 0.8926328 (2.9e-06 / 2.6e-05) | 0.8725508 (2.6e-06 / 6.8e-05) | 0.07149329 (1.7e-05 / 4.6e-05) | -1.849267 (5.5e-04 / 2.0e-03) | — |
| `run_1085_tewake_g4_L10_cont20k` | 0.8 | m10_on | g4 | JUDGED | FAIL | JUDGMENT_UNSETTLED | PASS/PASS/UNSETTLED | 0.9162843 (3.6e-07 / 6.2e-06) | 0.8965902 (6.4e-07 / 4.0e-05) | 0.06728739 (2.4e-06 / 5.9e-05) | -1.776773 (6.4e-05 / 1.7e-03) | — |

## 参考: g3 の設計差 (L_sw 1.0 − 0.8、現行の verdict が PASS の行だけ、q = C_T_with_shear)

| 作動点 | run 0.8 | run 1.0 | 量 | q(0.8) | q(1.0) | Δ | 時間幅 U = u₀.₈ + u₁.₀ |
|---|---|---|---|---|---|---|---|
- 未算出: m10_on: g3 の L_sw 0.8 (対応表に無い: 未実施) / 1.0 (対応表に無い: 未実施) の両方が現行で PASS でない
- 未算出: m6_on: g3 の L_sw 0.8 (run_1088_tewake_m6_L10: FAIL) / 1.0 (run_1097_tewake_lsw10_m6_cont20k: FAIL) の両方が現行で PASS でない
- 未算出: m4_off: g3 の L_sw 0.8 (対応表に無い: 未実施) / 1.0 (対応表に無い: 未実施) の両方が現行で PASS でない
- ΔJ: 未算出 (3 作動点の対がそろわない)
- m10_on の格子感度 G₁₀ = |Δ₁₀,g4 − Δ₁₀,g3| + U₁₀,g4 + U₁₀,g3: **未確定** — Δ₁₀,g4 に要る g4・L_sw 1.0・m10_on の判定区間が無い (g4 は L_sw 0.8 だけ)

## 使った問題 YAML (te_wake の YAML を基本に、evaluate の判定区間のキーだけ生産 YAML の値)

| run | te_wake の YAML (sha256・最終コミット) | 生産 YAML (sha256・最終コミット) | 生産 YAML から取った値 (te_wake の YAML の値) | 他に違うキーの数 |
|---|---|---|---|---|
| `run_1088_tewake_m6_L10` | problem_3d_prod_3op_wallres_lswx08_tewake_g3_L10.yaml (d281fd40292d・7c88aa35) | problem_3d_prod_3op_wallres_lswx08.yaml (df1969f2bf1e・fa816325) | judgment_steps 20000 ((無い))・judgment_extend_once 1 ((無い))・floor_events 1 (1)・floor_events_window_steps 10000 (10000)・outStepInterval 500 (500) | 0 |
| `run_1097_tewake_lsw10_m6_cont20k` | problem_3d_prod_3op_wallres_lswx10_tewake_g3_L10.yaml (7a259984523f・7f056898) | problem_3d_prod_3op_wallres_lswx10.yaml (9d8134c794b6・fa816325) | judgment_steps 20000 ((無い))・judgment_extend_once 1 ((無い))・floor_events 1 (1)・floor_events_window_steps 10000 (10000)・outStepInterval 500 (500) | 0 |
| `run_1085_tewake_g4_L10_cont20k` | problem_3d_prod_3op_wallres_lswx08_tewake_g4_L10.yaml (c6f01b51a87f・7c88aa35) | problem_3d_prod_3op_wallres_lswx08.yaml (df1969f2bf1e・fa816325) | judgment_steps 20000 ((無い))・judgment_extend_once 1 ((無い))・floor_events 1 (1)・floor_events_window_steps 10000 (10000)・outStepInterval 500 (500) | 8 |

## 各 run の記録

- `run_1088_tewake_m6_L10`: `/home/ubuntu/forge-r8/case/46.sern_design/run_1088_tewake_m6_L10/REJUDGE_2026-10-13` (REJUDGE_INFO.json・AB_GATES.md・REJUDGE_SUMMARY.txt・metrics.json・CONFIG_DIFF.txt) 終了 2026-10-13T02:49:45+0000・git 4829694052e1bedb53db25228b619701434535ca
- `run_1097_tewake_lsw10_m6_cont20k`: `/home/ubuntu/forge-r8/case/46.sern_design/run_1097_tewake_lsw10_m6_cont20k/REJUDGE_2026-10-13` (REJUDGE_INFO.json・AB_GATES.md・REJUDGE_SUMMARY.txt・metrics.json・CONFIG_DIFF.txt) 終了 2026-10-13T02:56:29+0000・git 4829694052e1bedb53db25228b619701434535ca
- `run_1085_tewake_g4_L10_cont20k`: `/home/ubuntu/forge-r8/case/46.sern_design/run_1085_tewake_g4_L10_cont20k/REJUDGE_2026-10-13` (REJUDGE_INFO.json・AB_GATES.md・REJUDGE_SUMMARY.txt・metrics.json・CONFIG_DIFF.txt) 終了 2026-10-13T03:06:03+0000・git 4829694052e1bedb53db25228b619701434535ca
