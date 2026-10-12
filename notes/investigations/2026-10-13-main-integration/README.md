# main 統合 (2026-10-13): feature/nozzle-wall-fit-and-pipeline

ユーザ決定 (2026-10-13): 「main には統合してもらってよい」。統合ブランチ `integrate/main-2026-10-13` = origin/main (90a09685、PR #6) ← feature/nozzle-wall-fit-and-pipeline。
試しの merge は衝突なしで、統合した木のファイルの中身は feature の先端と同一 (PR #6 の調整のルールの中身は feature に入っている)。main へは PR (merge commit) で入れる。gh は未認証なので PR はブラウザで作る (またはユーザが `gh auth login`)。

## 含む変更 (前回の統合 PR #5 の後)
- **CUDA ソルバの変更 34 commit** (`git log origin/main..HEAD -- solver_density_cuda/cuda_forge solver_density_cuda/input solver_density_cuda/main.cpp ...`)。前回の統合 (2026-10-07) はソルバの変更が無かったが、今回はある。
- design: CONTUR v2 (等温壁の熱の閉包)、`mesh2d` の第一層の厚さの表・密度・壁法線の近壁層のオプション、`msh_digits` ほか (8 ファイル)。

### 既定の振る舞いが変わるもの (検証の記録)
| 変更 | 既定 | 影響 | 検証の記録 |
| --- | --- | --- | --- |
| 軸対称の面ベクトルを双対の区間ごとの r で重み付け (`mesh.axisSegmentRWeight`、0b4dff4e) | 1 (修正あり) | 軸対称の計算の値が変わる | plans/active/axisymmetric-freestream-hoop-gauge.md §4.10〜§4.12 (閉性の修正の検証)、§4.13〜§4.16 |
| 変換器が幾何を double で書く (cb6a0d8b)・float のビルドが double の幾何を読む (段 1〜4、1da5593d〜e93361b8) | 常に | float のビルドの値が変わる。変換後の HDF5 は以前の版とビット同一にならない | plans/active/architecture-float-state-double-geometry.md §6 (V0〜V4) |
| double のビルドで格子の座標を double で読む (c096d66c) | double のビルドだけ | FP64 のビルドの値が変わる (冷却壁の第一層) | plans/active/tooling-nozzle-isothermal-wall-chain.md §5.1 #13〜#16 |
| ライン Thomas の配列の並び LAYOUT2 を既定に (fd2d4c3c) | 既定 | 速さだけ (補正と因子のビット比較で同一) | plans/active/time_integration-line-implicit-speed.md §6.7 |
| 前処理の block カーネルの node の行の処理 (e224ae11) と `speciesPrecondDt` | 既定 1 | node で `lowMachPrecond ≥ 2` の計算だけが変わる (境界の凍結の修正) | float plan §5.1 #26・§6.57・§6.59 |
| 方向別 dt の上限・熱の Jacobian・`lineViscCoupling 2/3`・各種の診断 | 0 / off | 変わらない (opt-in) | 各 plan |

## 事前登録 (確認の前に commit)
1. Python の単体試験: `solver_density_cuda/tests/unit/test_*.py` と `solver_density_cuda/tools/test_*.py` が main (90a09685) と同じ判定 (main で PASS のものが統合版でも PASS)。新しい試験は PASS。
2. design の試験 (`design/.venv-opt` の Python): `design/tests/run_*.py` の FAIL が main と同じか少ない。新しく FAIL になるものが無い。
3. `check_plans.py` の FAIL 件数が main から増えない。
4. 設計チェーンの準備の再現 (CFD 0 step、変換の前で比べる): 生産の問題 `case/45.isobutane_m6_d155/problem_d155_ns_prod.yaml` の準備が出す `.msh`・`wall_design.csv`・`delta_r_initial.csv`・`bcondConfig.yaml` が、前回の統合の版 (c600b507) とバイト同一。変換後の HDF5 は変換器の変更 (cb6a0d8b) で設計上ビット同一にならないので比べない。
5. CUDA のビルド: float と FP64 の両方が通る。CUDA の単体試験 (`tests/unit/test_precond_node_rows.cu`・`test_quantize_state.cu`・`test_qacc_axis.cu`) が PASS。統合の木のソルバは feature の先端と同一なので、先端 (e224ae11 以降) でのビルドと試験の記録を使ってよい。
6. 移植の回帰 R1〜R3 (float plan §6.59): R1・R3 合格、R2 結果 A。
7. 上の「既定の振る舞いが変わるもの」の表を PR の説明に載せる。

FAIL のときは main へ入れない (原因を切り分ける)。

## 実行の記録
- **4: 実行不能 (入力が消えていた、2026-10-13)**。生産の問題の準備は Euler の参照 run `run_0062_euler_wallfit_fit_r1_ext6k` の `nozzle.h5` (初期線を作るのに使う) を読むが、A のこの run の `nozzle.h5` は 2026-10-10 の後片付けで消えていた (`case/45/_disk_cleanup_2026-10-10.tsv` の 11 行目、3 ファイル・27 MB)。継続の元の `run_0056` にも `nozzle.h5` は無い。A で前回の版 (c600b507) と先端 (b155b39c) の両方で準備を試し、両方ともこの読み込みで止まった (`~/integ-p4-{old,tip}.log`)。作り直すには Euler の連鎖を回し直す必要があり、今回は行わない。**残る危険**: 前回の統合の後の design の変更 (CONTUR v2、`mesh2d` の第一層の表・密度・壁法線の近壁層、`msh_digits`) が、新しいキーを書かない生産の問題の準備を変えないことは、確認 2 の design の試験 (`run_contur_v2_tests.py`・`run_mesh_*` ほか) で確かめた範囲に限る。基準をゆるめたのではなく、確かめられなかったこととして残す。
