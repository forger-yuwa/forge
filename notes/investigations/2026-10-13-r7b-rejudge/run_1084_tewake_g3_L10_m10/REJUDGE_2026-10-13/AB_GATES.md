# A/B: 受理検査だけを旧版 (保存された metrics.json) / 現行版 (sern-acceptance-v1) で変える — run_1084_tewake_g3_L10_m10

- 旧: `/home/ubuntu/forge-r8/case/46.sern_design/run_1084_tewake_g3_L10_m10/metrics.json` (sha256 8a8e02e491a64bcc58d6b38c45af4adda9137dfa2d5c61adc57826e3992c3dff, mtime 2026-10-10T17:31:27+0000) — run の時点の collect の出力
- 現行: `/home/ubuntu/forge-r8/case/46.sern_design/run_1084_tewake_g3_L10_m10/REJUDGE_2026-10-13/metrics.json` — `runner_sern3d.collect(problem_rejudge.yaml, run, out_dir=...)`

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
| C_T | 0.9163372 | 1.61e-06 | 1.28e-06 | 5e-05 | STEADY | OK |
| C_T_with_shear | 0.8969471 | 1.42e-06 | 1.39e-06 | 5e-05 | STEADY | OK |
| C_L | 0.06748933 | 6.31e-06 | 1.96e-05 | 5e-05 | STEADY | OK |
| C_M | -1.783215 | 1.40e-04 | 4.27e-04 | 5e-04 | STEADY | OK |
