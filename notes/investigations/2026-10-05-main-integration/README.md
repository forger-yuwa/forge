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

## 結果 (2026-10-05)
1. **PASS**: `solver_density_cuda/tests/unit/test_*.py` 21 本が main (77318d0e) と統合版で同じ判定 (両方とも全件 rc 0)。
   `design/tests/run_deltastar_tests.py` ALL PASS、`test_interp_field_3d.py` ALL PASS、`test_gate_bad_input.py` PASS。
2. **PASS**: `check_plans.py` の FAIL は main 28 件 (21/49 OK)、統合版 28 件 (22/50 OK; 1 件増えたのは feature 側の plan `verification-m6-axis-wave-mesh-su2` で OK 側)。
3. **FAIL**: 統合版で最終設計の壁を作り直すと run_0051 と一致しない。
   - まず main では semiperfect TP の NS/SST に**化学種ごとの輸送物性 (`gas.transport`、viscMethod 2) が必須**になっており、feature 側で作った
     `problem_d155_ns_c2final{,_cond}.yaml` (Sutherland の時代の rt77p02 から複製) は prepare で止まる。main が既存の問題ファイルに入れたのと同じ
     `transport` 行を 2 ファイルに追加した (この統合ブランチの変更)。壁の形は積分法 (粘性は内部の Sutherland) で決まり transport には依存しない。
   - そのうえで壁の差は最大 **5.0e-7 m (0.5 µm、相対 6.5e-7)**、`nozzle.msh` の md5 不一致。差の出どころは**設計壁 (逆 MOC)**: 軸 M の目標 3e-9、
     設計壁 2.3e-7 m (出口側)。main 側の気体熱物性の変更 (`gas/semiperfect.py`・`gas/composition.py`: NASA-9 区間可変化と組成の扱い) で物性が 1e-8 級に動いたもの。
     feature 側の変更や衝突解消の誤りではない。形状としては無視できる大きさだが、事前登録の基準 (≤1e-12 m) は満たさない。
   - **より重要な含意**: main の NS は粘性モデルが変わっている (Sutherland → 種ごとの CEA/IAPWS)。最終設計の検証 run (run_0051/0052) は Sutherland で
     回したもので、main で同じ設計を回すと境界層 (出口 δ) が変わりうる → 出口較正の k_f と r_t を main の物性で取り直す必要がある (未実施)。
4. **PASS**: 統合版の `nozzle_report` を run_0051 に回した評価量 29 項目が統合前と完全一致 (相対差 0)。統合用の作業ツリーには `.venv-pptx` が無いので pptx は作られない (図と json のみ; 想定どおり)。

→ 事前登録どおり 3 が FAIL なので、main への取り込みはユーザ判断を待つ (原因は特定済み: main 側の熱物性・輸送物性の変更)。

**ユーザ判断 (2026-10-05)**: 3 の 0.5 µm の差は許容として記録し、PR で main に取り込む。取り込み後、main の物性 (種ごとの輸送物性) で最終設計
(NS + 凝縮 ON) を回し直し、出口較正の k_f と r_t を確認する。
