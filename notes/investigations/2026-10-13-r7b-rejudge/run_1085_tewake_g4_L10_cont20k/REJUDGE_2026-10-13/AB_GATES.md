# A/B: 受理検査だけを旧版 (保存された metrics.json) / 現行版 (sern-acceptance-v1) で変える — run_1085_tewake_g4_L10_cont20k

- 旧: `/home/ubuntu/forge-r8/case/46.sern_design/run_1085_tewake_g4_L10_cont20k/metrics.json` (sha256 6da44d7eb0387364602514a52ad70366ac5917429de99a2afe5cc58b28d9ce4e, mtime 2026-10-10T19:13:10+0000) — run の時点の collect の出力
- 現行: `/home/ubuntu/forge-r8/case/46.sern_design/run_1085_tewake_g4_L10_cont20k/REJUDGE_2026-10-13/metrics.json` — `runner_sern3d.collect(problem_rejudge.yaml, run, out_dir=...)`

| 項目 | 旧 | 現行 | 同じ |
|---|---|---|---|
| verdict | PASS | FAIL | ≠ |
| fail_class | None | JUDGMENT_UNSETTLED | ≠ |
| acceptance_version | (なし) | sern-acceptance-v1 | ≠ |
| rc | 0 | 0 | = |
| 保存場 field | OK | OK | = |
| 残差 residual | OK NOT CONVERGED (stalled/plateau) | OK NOT CONVERGED (stalled/plateau) | = |
| 残差の入力 residual_evidence | input_ok None・欠けた列 None・時系列 None 行 (step None..None、最終場 None) | input_ok True・欠けた列 []・時系列 20000 行 (step 0..19999、最終場 res_20000.h5) | ≠ |
| 4 量の頭打ち steadiness | OK C_T_with_shear STEADY・C_T STEADY・C_L STEADY・C_M STEADY | OK C_T_with_shear STEADY・C_T STEADY・C_L STEADY・C_M STEADY | = |
| 最終場の床 floors | OK | OK | = |
| #1 床事象 floor_events | OK PASS 区間 [10000, 20000] 更新 10000 必須 True | OK PASS 区間 [10000, 20000] 更新 10000 必須 True | = |
| 力の部品 force_parts | (なし) | OK 不完全 0 | ≠ |
| 判定区間 judgment (#2・#3) | (なし) | UNSETTLED JUDGMENT_UNSETTLED | ≠ |
| #2 温度 temperature | (なし) | PASS 窓 [10000, 20000] 標本 20 | ≠ |
| #3 力の窓 force_window | (なし) | UNSETTLED 窓 [10000, 20000] 標本 20/20 | ≠ |
| 残差の桁 residual_scale | OK | OK | = |
| 理由 reasons | (なし) | 判定区間: 窓・定常性が未達: C_L: a 2.43e-06・d 5.87e-05 (許容 5e-05)・STEADY; 判定区間: 窓・定常性が未達: C_M: a 6.35e-05・d 1.72e-03 (許容 5e-04)・STEADY | ≠ |

判定区間 #3 の窓 (現行、末尾 W step の平均 m・a = max|q − m|・前窓との差 d):

| 量 | 平均 | a | d | 許容 | 定常 | 合否 |
|---|---|---|---|---|---|---|
| C_T | 0.9162843 | 3.62e-07 | 6.18e-06 | 5e-05 | STEADY | OK |
| C_T_with_shear | 0.8965902 | 6.42e-07 | 3.96e-05 | 5e-05 | STEADY | OK |
| C_L | 0.06728739 | 2.43e-06 | 5.87e-05 | 5e-05 | STEADY | NG |
| C_M | -1.776773 | 6.35e-05 | 1.72e-03 | 5e-04 | STEADY | NG |
