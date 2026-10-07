# main 統合 (2026-10-07): feature/nozzle-wall-fit-and-pipeline

統合ブランチ `integrate/main-2026-10-07` = origin/main (77318d0e、2026-10-04 統合 PR #3) ← feature/nozzle-wall-fit-and-pipeline。
feature 側は 2026-10-05 の統合ブランチ `integrate/main-2026-10-05` (feature/gap-heating-precision の取り込み、main へは未マージ) を含み、
その後に case/45 の M6 ノズルの設計チェーンの改良を積んだ。main はこのブランチの祖先なので、統合は衝突なしの merge commit で、
統合ブランチのファイルの中身は feature の先端と同一になる。

含むもの (accepted の plan):
- `plans/accepted/tooling-nozzle-throat-monotone-r2.md` (単調壁)
- `plans/accepted/tooling-nozzle-wall-single-bspline.md` (全域 1 本の B-spline 壁と STEP)
- `plans/accepted/tooling-nozzle-upstream-poly-and-throat-sizing.md` (上流の 5 次多項式・スロート径からの寸法の逆算)
- `plans/accepted/discretization-moc-axis-limit-and-corrector.md` (MOC の軸上の解析極限・収束する予測修正。**axis-Mach の問題の MOC の既定を analytic + converge に変更**)
- `plans/accepted/verification-case45-euler-total-enthalpy.md` (Euler 専用の格子 `mesh_euler`)
- 2026-10-05 の統合ブランチの内容 (`notes/investigations/2026-10-05-main-integration`)

**CUDA ソルバ (`solver_density_cuda/cuda_forge/` ほか C++) の変更は無い** (`git diff --stat origin/main..HEAD -- solver_density_cuda` は tools の Python だけ)。

## 事前登録 (確認の前に commit)
1. Python の単体試験: `solver_density_cuda/tests/unit/test_*.py` と `solver_density_cuda/tools/test_*.py` が main (77318d0e) と同じ判定 (main で PASS のものが統合版でも PASS)。
2. design の試験 (`design/.venv-opt` の Python): `design/tests/run_*.py` の FAIL が main と同じか少ない。新しく FAIL になるものが無い (既知の 3 本 `run_sern_gates_tests.py`・`run_sern_moc_tests.py`・`run_species_attrs_ic_tests.py` は main でも FAIL のはず — main で確かめる)。
3. `check_plans.py` の FAIL 件数が main から増えない。
4. 設計チェーンの再現: 生産の問題 `case/45.isobutane_m6_d155/problem_d155_ns_prod.yaml` の準備が生産の run (run_0167_ns_n012_N2) の入力を再現する (`prod_confirm.py`、AWS、feature の c600b507 で PASS 済み)。統合ブランチの中身が feature の先端と同一 (`git diff` が空) であることを確かめて、この結果を引き継ぐ。
5. MOC の既定の切り替えの波及 (`case/45.isobutane_m6_d155/moc_default_sweep.py`、CFD 0 step): MOC のキーを書いていない axis-Mach の問題 (case/41〜45) で新しい既定のゲートが不合格になるものがあれば、ユーザに諮ってから統合する (統合の前に問題 YAML に旧方式を明示するか、そのままにするか)。

FAIL のときは main へ入れない (原因を切り分ける)。main へは PR (merge commit) で入れる。gh は未認証なので PR はブラウザで作る。
