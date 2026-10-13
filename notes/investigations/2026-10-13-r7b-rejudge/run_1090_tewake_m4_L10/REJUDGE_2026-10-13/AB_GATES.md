# A/B: 受理検査だけを旧版 (保存された metrics.json) / 現行版 (sern-acceptance-v1) で変える — run_1090_tewake_m4_L10

- 旧: `/home/ubuntu/forge-r8/case/46.sern_design/run_1090_tewake_m4_L10/metrics.json` (sha256 9c737c981c70c391d5061b5ecf9b5ecf9c6320619874d2528d5498ffe2759b5b, mtime 2026-10-11T03:09:25+0000) — run の時点の collect の出力
- 現行: `/home/ubuntu/forge-r8/case/46.sern_design/run_1090_tewake_m4_L10/REJUDGE_2026-10-13/metrics.json` — `runner_sern3d.collect(problem_rejudge.yaml, run, out_dir=...)`

| 項目 | 旧 | 現行 | 同じ |
|---|---|---|---|
| verdict | PASS | PASS | = |
| fail_class | None | None | = |
| acceptance_version | (なし) | sern-acceptance-v1 | ≠ |
| rc | 0 | 0 | = |
| 保存場 field | OK | OK | = |
| 残差 residual | OK NOT CONVERGED (stalled/plateau) | OK NOT CONVERGED (stalled/plateau) | = |
| 残差の入力 residual_evidence | input_ok None・欠けた列 None・時系列 None 行 (step None..None、最終場 None) | input_ok True・欠けた列 []・時系列 20000 行 (step 0..19999、最終場 res_20000.h5) | ≠ |
| 4 量の頭打ち steadiness | OK C_T_with_shear STEADY・C_T STEADY・C_L STEADY・C_M STEADY | OK C_T_with_shear STEADY・C_T STEADY・C_L STEADY・C_M STEADY | = |
| 最終場の床 floors | OK | OK | = |
| #1 床事象 floor_events | OK PASS 区間 [10000, 20000] 更新 10000 必須 True | OK PASS 区間 [10000, 20000] 更新 10000 必須 True | = |
| 力の部品 force_parts | (なし) | OK 不完全 0 | ≠ |
| 判定区間 judgment (#2・#3) | (なし) | PASS | ≠ |
| #2 温度 temperature | (なし) | PASS 窓 [10000, 20000] 標本 20 | ≠ |
| #3 力の窓 force_window | (なし) | PASS 窓 [10000, 20000] 標本 20/20 | ≠ |
| 残差の桁 residual_scale | OK | OK | = |
| 理由 reasons | (なし) | (なし) | = |

判定区間 #3 の窓 (現行、末尾 W step の平均 m・a = max|q − m|・前窓との差 d):

| 量 | 平均 | a | d | 許容 | 定常 | 合否 |
|---|---|---|---|---|---|---|
| C_T | 0.9791903 | 2.34e-07 | 1.21e-06 | 5e-05 | STEADY | OK |
| C_T_with_shear | 0.9676994 | 2.53e-07 | 1.86e-06 | 5e-05 | STEADY | OK |
| C_L | -0.02486542 | 1.90e-06 | 1.33e-05 | 5e-05 | STEADY | OK |
| C_M | 0.5040016 | 3.88e-05 | 3.93e-04 | 5e-04 | STEADY | OK |
