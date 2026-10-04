# main 統合 (2026-10-05): feature/gap-heating-precision

統合ブランチ `integrate/main-2026-10-05` = origin/main (77318d0e、2026-10-04 統合 PR #3 の後) ← feature/gap-heating-precision (d7ec3e14)。
feature 側の main から分かれた後の 50 commit は、M6 ノズル設計 (plan `verification-m6-axis-wave-mesh-su2`: 抽出の帯選択 E、壁入力の 5 次補間、
CONTUR 較正つまみ、ノズル設計の標準出力と pptx 報告 `design/forge_design/report/`、手順書 `procedures/nozzle-design-outputs.md`) と、
gap heating (case/56・59・60 の記録、`interp_field.py` の 3D 最近傍修正) の文書・Python ツール。**CUDA ソルバ (`solver_density_cuda/cuda_forge/` ほか C++) の変更は無い**。

## 衝突の解消 (2 件、両側の追加)
- `design/forge_design/evaluate/runner_axismach.py` `prepare_ns` の `prepare_info`: main 側の種ごとの輸送物性の記録 (`info["transport"]`) と、feature 側の
  実際の nj・`axis_gap_frac` の記録の両方を残す。
- `plans/README.md` active 一覧: 両側の追加行を両方残す。

## 事前登録 (確認の前に commit)
ソルバの変更が無いので CUDA の回帰は行わない (main の回帰結果がそのまま有効)。確認は Python と設計チェーンの再現:
1. Python 単体試験: `solver_density_cuda/tests/unit/test_*.py` が main (77318d0e) と同じ判定 (main で PASS のものが統合版でも PASS)。
   `design/tests/run_deltastar_tests.py` ALL PASS、`solver_density_cuda/tools/test_interp_field_3d.py` PASS、`solver_density_cuda/tools/test_gate_bad_input.py` PASS。
2. `check_plans.py` の FAIL 件数が main から増えない。
3. 設計チェーンの再現: 統合版で最終設計 run_0051 の壁を作り直す (`prep_contur_cal.py`、`problem_d155_ns_c2final.yaml`、k_f 1.0257341596、IC なし) と、
   `nozzle.msh` が run_0051 と md5 一致、`wall_physical.csv` の最大差 ≤ 1e-12 m。
4. 報告ツール: 統合版の `nozzle_report` を run_0051 に回し、`report.json` の評価量 (波・オーバーシュート・出口コア M・ṁ 比) が統合前と一致 (相対 1e-9)、pptx が生成される。
FAIL のときは main へ入れない (原因を切り分ける)。main へは PR で入れる (ユーザ指示 2026-10-05「main 取り込むとともにこの成果を main に入れて」)。
