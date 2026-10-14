# 諮問ブリーフ: G0c (縮流部を粗くした格子) と B0 (生産の陰解法) を設計チェーンに入れる設計と、チェーンの開始の場からの検証 (2026-10-13)

- 依頼者: 主セッション (エスカレーション条件 1: plan の §4 設計方針・§6 検証計画を新規に書く)
- plan: `plans/active/tooling-nozzle-isothermal-wall-chain.md` §4.8・§5.1 #29〜#31・§6 V-g0c-chain (起案、諮問の前)。関連 `plans/active/architecture-float-state-double-geometry.md` §5.1 #25・#28、§6.58・§6.61 の結果
- 前の諮問: `notes/reviews/2026-10-13-speed-s-result-diagnose.md` (G0c・FP64・上限 50 を登録の条件に限って採用。次は「設計チェーンが実際に作る NS の開始の場から 1 本」、一様場から本段へ直接は不可、開始の場は生成の手順と入力の sha で定義)
- ユーザ (2026-10-13): 「よい」(チェーンが G0c を作れるようにして、事前登録してから 1 本確かめる段取りを承認)。目的は case/45 の NS を速くして、壁温を変える設計へ進むこと。

## 観測事実 (2026-10-13 の調査、file:line は作業ツリー `/home/sano/work/forge-integ-1005`)

- G0c の作り方: `case/45.isobutane_m6_d155/make_g0c.py` が forge-25 の `core_grid_mesh.py` (commit 96197821) を `git show` で読み、`use_grefc_stations(True)` で `mesh2d._x_stations` を差し替える。区間 `GREFC_SEG=(−8, −7.5, −2.5, −2)`、`ni_ref` 4719 固定、密度は `x_density_table` の線形補間だけ (スロートのガウスは無視、この YAML では throat_refine 1.0)、20,001 点で積分、d_n = d·(1 − w(1 − max(1, d/2)/d))、n_new = round(589 × 0.5729) = 337、ni = 4719 − 589 + 337 = 4467。区間の外の station は元の配列のまま (`make_g0c.py:38-39` で `array_equal`)。`.msh` は 17 桁、sha256 285f0fa4ecc0ed1c。
- `mesh2d._x_stations` (`design/forge_design/meshing/mesh2d.py:103-122`) は密度を 8,001 点で取り全域で正規化して ni 等分の逆写像。区間の外を変えないオプションは無い。半径方向・壁法線の層は列ごとに x だけで決まる。
- チェーンの NS: `prepare_ns` (`design/forge_design/evaluate/runner_axismach.py:1252`) は準一次元の場を貼り、`ic_from` があれば `interp_field.py` で移す (1421-1454)。段階起動 `run_staged_ns` (1539) full = S1_soft (1 次・cfl 0.5・nStepInner 10・3,000 step) → S2_mid (1 次・cfl 1・3,000 step) → 本段 (2 次・`cfl_main`) (1585-1613)。本段の設定は `runner_wt._config_sst_node` (247-275)。**B0 のキーはチェーンに無い**: case/45 の `cold_cfl.py prep` (122-157) が足し、`implicitRelax 0.7` は `prepare_ns(implicit_relax=0.7)` (1413)。B0 の定義は `plans/active/time_integration-line-implicit-speed.md:450`。
- δ* ループの各 pass は `ic_from` = 前の pass、段階起動 full (`design/forge_design/feedback/deltastar_loop.py:363-381`)。
- 生産の場の系譜: `run_0179_ns_n012_N2_ext` res_20000 → (`cold_pair.py prep`、格子をまたいで移し、段階起動 full、本段 100,000 step、cfl 1、relax 0.7、FP64) `run_0182` → `run_0183` (+100,000) → (`hp7.sh`、B0、格子を変換し直し) `run_0483` (到達 145,000) → Q32 → G0c の起点。`run_0179` res_20000 と `_band_ab/prod_confirm/prep` は AWS A にある。
- 停止規則 (float plan §6.11、`case/45/m9_watch.py`) はチェーンに無い。
- forge-25 の core-grid plan は「メッシャのコードは変えない、生産の格子の置き換えはしない」(試験用の格子として扱う)。

## 設計の案 (plan §4.8)

1. メッシャに `mesh.x_respace: [{x0, x1, blend, factor, ni_ref}]`。元の station を ni_ref で作り、区間の中だけ上の式で置き直す。結果の ni は YAML の `ni` と照合 (違えば止める)。密度は `_x_stations` と同じ式から取る。合格: case/45 の tw300 の YAML に (−8, −2, 0.5, 2, 4719)・ni 4467 を書いた .msh が既存の G0c とバイト一致、`x_respace` の無い YAML は今と同じ出力。
2. `evaluate.main_line: {line: dir, cap: 50, itj: 5, cfl: 4, relax: 0.7}` で本段だけに B0 のキーを書く (S1・S2 は今のまま)。合格: 本段の solverConfig が `run_0717_c50_d1` と step 数・出力の間隔・ファイル名以外でキーごとに一致。キーの無い問題は今と同一。
3. G0c の YAML は `ar_max: 10000` と書き、run の README に AR 緩和を明記 (一般の基準の変更ではない)。
4. 停止規則はチェーンに入れず、検証の台本で本段の区間だけに `m9_watch.py` を当てる。

## 検証の案 (§6 V-g0c-chain)

- チェーンの `prepare_ns` を上の YAML で、`ic_from = run_0179_ns_n012_N2_ext` (res_20000)、`delta_r_csv = _band_ab/prod_confirm/prep/delta_r_initial.csv`・offset radial で呼ぶ → 段階起動 full → 本段 (G0c・FP64 7e95e5bb・B0) を最大 20 万 step・§6.11 の停止規則。
- 停止の出力を T2a・T2b と M1〜M7 (float plan §6.58 と同じ式・許容) で比べる。対照の出力は `run_0804_sp_g0c` (再投入しない)。総時間と本段の時間を分ける。A = 到達かつ全比較が許容内 → 適用の拡大を支持。B = 未到達・発散・許容超過 → 拡大を退ける。欠損は判定不能。

## 問い

1. メッシャのオプション (一般化した `x_respace`) と「チェーンの出力に変換を掛ける」のどちらがよいか。一般化のパラメータ (x0・x1・blend・factor・ni_ref) と、密度を `_x_stations` と同じ式から取る (forge-25 版は表だけ) ことの扱い。バイト一致の合格条件は妥当か。
2. 本段だけ B0、S1・S2 は今のまま、でよいか。B0 の値を YAML に持たせる形 (`main_line`) でよいか。
3. 開始の場として `run_0179` res_20000 (生産の系譜と同じ、別の格子・早い段の NS) + 段階起動 full でよいか。準一次元の場から (ic_from 無し) も要るか。予算 20 万 step は本段だけに数えるか、段階起動を含めるか。
4. 判定 (M1〜M7 対 T2、停止規則を本段の区間だけに当てる) で足りないもの。段階起動の各段のゲート・本段の区間の `check_convergence --segment` の扱い。
5. この検証が A でも言えないこと (壁温を変えた問題、δ* ループの pass の間、float のビルド)。

## 読んでよいファイル

- `plans/active/tooling-nozzle-isothermal-wall-chain.md` §4.8・§5.1 #29〜#31・§6
- `plans/active/architecture-float-state-double-geometry.md` §6.11・§6.58・§6.61 と結果、§5.1 #25・#28
- `design/forge_design/meshing/mesh2d.py`、`design/forge_design/evaluate/runner_axismach.py` (prepare_ns・run_staged_ns)、`design/forge_design/evaluate/runner_wt.py` (`_config_sst_node`)
- `case/45.isobutane_m6_d155/make_g0c.py`・`cold_cfl.py`・`cold_pair.py`・`hp7.sh`・`m9_watch.py`・`problem_d155_ns_prod_coldmesh_tw300.yaml`
- forge-25 の `core_grid_mesh.py`: `git show 96197821:case/45.isobutane_m6_d155/core_grid_mesh.py`
