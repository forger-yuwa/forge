# 設計チェーンの起動の経路の棚卸し (2026-10-14)

plan [tooling-generator-defaults-and-launch-gate](../../plans/active/tooling-generator-defaults-and-launch-gate.md) §4.5・§5.1 #9 の棚卸し。
`design/forge_design` から forge (`run_case.sh`) を起動する全経路について、段ごとに「段の名前・最終の solverConfig の作り方・段間の引き継ぎ・
履歴の退避・`StageManifest` の記録」を表にした。読んだのは `feature/generator-defaults-launch-gate` の 570790e8 (main PR #9 を取り込み済み)。
「変更後」の列は同日の #9 の実装 (この表と同じ commit) の後の状態。

## 0. 要約

- 起動の経路は 15 (入口の関数・CLI の単位。§1 の表の行)。forge を起動するのは全経路とも `run_forge` → `run_case.sh` (runner.py・runner_sern.py の 2 つの実装)
  か、停止規則の見張り (`stop_rule.watch_main_stage` が `run_case.sh` を子プロセスで持つ) だけ。
- **`StageManifest` を書いていなかった経路**: runner_wt の段階起動 (WT モード F)、walldriven、SERN 2D・3D (driver_sern を含む) — **#9 で 3 つとも書くようにした**。
  1 段だけの経路 (runner.py の CPG・TP の単段、opt/driver の CPG、polish の CPG、euler_loop) は今も書かない (段が 1 つで区間の判定が要らない。
  ただし driver の CPG の soft は別ディレクトリの別の run で、本段とは区間がつながらない)。
- **同一メッシュの段間に `interp_field.py` を使っていた所**: runner_wt.run_staged、runner_walldriven.run_staged — **#9 で段終了ゲート + `restart_field.py` に替えた**。
  **未対応 (範囲外として残した)**: euler_loop のチャンクの継続 (`cont_<j>`、同じ格子の nozzle.h5 を写して `interp_field.py`)、
  opt/driver の CPG の soft → 本段と opt/polish の CPG (同じ幾何から作り直した別ディレクトリの格子へ `interp_field.py`。格子の生成が決定的なら同一メッシュ)。
- **段終了ゲート (必要保存量・有限・ρ>0) の無い段間**: SERN の `restart_by_index` (化学種の照合はするが、非有限・負の密度は見ない)。範囲外として残した (§3)。
- 段の config を**正規表現・文字列置換**で作っている経路: runner_wt・walldriven・SERN・opt/driver の CPG (§3.1 のとおり)。axis-Mach・runner.py の TP は YAML の位置
  (`_cfg_set`) か関数の引数。書き換えの方式の統一は plan §5.1 #11。

## 1. 経路の表

略号: R = `runner.py`、W = `runner_wt.py`、A = `runner_axismach.py`、D = `runner_walldriven.py`、S = `runner_sern.py`、S3 = `runner_sern3d.py`。
「config」の列は最終の solverConfig (forge が読む文字列) の作り方。`prepare` の品質検査用の cell 変換 (一時コピー、起動しない) は除く。

| # | 経路 (入口) | 機種 | 段 (tag) | 最終の config の作り方 | 段間の引き継ぎ | 履歴の退避 | `StageManifest` (前 → 後) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | `R.main` (CPG、`--stage` なし) = `R.prepare` → `R.run_forge` | ③ベル CPG (NS SST・node) | main だけ | `R._solver_config` (cfl 4・2 次・`nStepInner 5`・`katoLaunder 1`) | node は `--warm-from` 必須: 別の run の res を `interp_field.py` (cross-mesh) | なし (`residual_history.csv` のみ) | なし → なし |
| 2 | `R.main --stage S` (TP) = `R.prepare(stage)` → `R.run_forge` | ③ベル TP (NS・Euler) | soft / mid / main のどれか 1 つ | `R._solver_config_tp(p, n, oi, stage)` (`TP_STAGES`) | `--warm-from` は `interp_field.py` (cross-mesh)、冷間は soft だけ | なし | なし → なし |
| 3 | `R.main --staged` = `R.run_staged_tp` (`opt/driver` の TP 評価も) | ③ベル TP | `S1_soft` → `S2_mid` → `S3_main` (`evaluate.stages`) | 段ごとに `R._solver_config_tp(p, n, oi, stage)` | 段終了ゲート (`A.stage_gate`) + `restart_field.py` (`A._restart_same_mesh`) | `residual_history_<tag>.csv` | あり → あり |
| 4 | `R.main --continue-from` = `R.continue_tp` (`opt/polish` の TP も) | ③ベル TP | 親の段 + `S<k>_main_p` | `R._solver_config_tp(p, n, oi, "main")` | 親の最終 res を `restart_field.py` | `residual_history_<tag>.csv` (親の履歴も写す) | あり (親の記録に足す) → あり |
| 5 | `opt/driver` CPG の評価 (`_evaluate`、`<tag>_stgA` → `<tag>`) | ③ベル CPG | soft (別ディレクトリ `<tag>_stgA`) → main (`<tag>`) | soft = `R.prepare(nsteps=3000)` の config を文字列置換 (cfl 4→0.5・`convMethod 1→0`・`outStepInterval`)、main = `R.prepare` の config | soft: `--warm-from` 相当の `interp_field.py` (最近傍 PASS 評価の res)。main: 同じ幾何から作り直した格子へ soft の res を `interp_field.py` | 2 つの run に分かれる (段の履歴の退避なし) | なし → なし |
| 6 | `opt/polish` CPG (`polish_run`、`<tag>_p<k>`) | ③ベル CPG | main の継続 (チャンク) | `R.prepare` の config | 同じ幾何から作り直した格子へ前のチャンクの res を `interp_field.py` (docstring も「名前に反して同一メッシュの継続ではない」) | run ごと | なし → なし |
| 7 | `W.main` = `W.prepare(euler=True)` → `W.run_staged` | ① 風洞 モード F (Euler、**cell**) | soft → main | main = `W._config_euler` (cell・cfl 4・`nStepInner 5`)、soft = 文字列置換 (cfl 0.5・`convMethod 0`・`nStepInner 10`・3000 step) | **前: `interp_field.py` (同一メッシュ) → 後: 段終了ゲート + `restart_field.py`** | **前: なし → 後: `residual_history_<tag>.csv`** | **なし → あり** |
| 8 | `W.prepare(euler=False, ic_from)` → `W.run_staged` (NS v1、design の中に呼び出し元なし。case のスクリプトから) | ① 風洞 NS v1 (node SST) | soft → main | main = `W._config_sst_node`、soft = 行 7 と同じ置換 | `ic_from` は `interp_field.py` (cross-mesh) + ω 床。段間は行 7 と同じ | 行 7 と同じ | **なし → あり** |
| 9 | `feedback/euler_loop.run_loop` (`pass_<k>`、`cont_<j>`) | ① 風洞 B10 の帰還 (Euler、node 既定 / cell) | パスごとに main 1 段 + 継続チャンク | `W._config_euler_node` (node) か `W._config_euler` (cell)、CHUNK step | パス ≥ 2: 前パスの res を `interp_field.py` (壁が動くので cross-mesh で正しい)。**継続チャンク `cont_<j>`: 同じ格子 (nozzle.h5 を写す) へ `interp_field.py` (同一メッシュ)** | ディレクトリごと | なし → なし |
| 10 | `A.main` (Euler) = `A.prepare` → `A.run_staged` | ① 風洞 axis-Mach Euler (node、CPG / TP) | soft → [mid] → main (`--stages full/soft/none`、mid は `evaluate.mid_stage` か TP) | main = `W._config_euler_node` → `A._apply_gas_to_config` (TP の physProp) → `--cfl` / `--implicit-relax` の置換。soft = `A._first_order(A._cfg_set(cfl 0.5))`、mid = `_cfg_set(cfl 1.0)`、段の長さ 3000 (YAML の位置で書き換え) | 段終了ゲート + `restart_field.py` | `residual_history_<tag>.csv` | あり → あり |
| 11 | `A.prepare_ns` → `A.run_staged_ns` (`feedback/deltastar_loop` の 3 か所、case/45 のスクリプト群) | ① 風洞 axis-Mach NS (node SST、TP) | `full`: `S1_soft` → `S2_mid` → main / `ramp`: `R<i>_cfl<c>` → main / `none`: main | main = `W._config_sst_node` → `_apply_gas_to_config(viscous=True)` → `implicit_relax` の挿入 (+ `evaluate.main_line` の B0 は `A._apply_main_line`、前段は `PRESTAGE_CONFIG`)。前段は `A.ns_stage_configs` (`_ns_prestage`: 1 次・cfl 0.5/1.0・`nStepInner 5→10`) | 段終了ゲート + `restart_field.py`。`ic_from` は `interp_field.py` (cross-mesh) + ω 床 | `residual_history_<tag>.csv`、`MAIN_STAGE_STARTED`、停止規則の記録 | あり → あり |
| 12 | `D.main` = `D.prepare` → `D.run_staged` | ① 風洞 wall-driven (Euler、node) | soft → main | main = `W._config_euler_node` (cfl 4)、soft = 文字列置換 (cfl 0.5・`convMethod 0`・3000 step。`nStepInner` は 5 のまま) | **前: `interp_field.py` (同一メッシュ) → 後: 段終了ゲート + `restart_field.py`** | **前: なし → 後: `residual_history_<tag>.csv`** | **なし → あり** |
| 13 | `S.main` = `S.prepare` → `S.run_staged` (`opt/driver_sern` の評価も) | ⑤ SERN 2D (Euler / SST、node / cell、CPG / frozen_tp) | [`warm_adapt`] → [`warm_cfl<c>` …] → `soft_cfl<c>` … → [`mid_cfl<c>` …] → `main` (`opt:` の梯子。`--stages main/none`) | main = `S._solver_config` (`solverConfig_main.yaml`)、前段は正規表現・文字列置換 (cfl・`convMethod`・層流暖機の `turbulence`・段の長さ) | `restart_by_index` (同一メッシュの index コピー・化学種の照合。段終了ゲートなし)、warm start は `warm_from_run` (作動点間の相似リマップ) | `residual_<tag>.csv` (前: 暖機・soft・mid だけ → 後: 適応段・本段も) | **なし → あり** |
| 14 | `S3.main` = `S3.prepare` → `S.run_staged` | ⑤ SERN 3D | 行 13 と同じ (`--soft-steps`、`opt:` の梯子。`warm_lam_ramp`・`soft_ramp` は渡さない) | main = `S._solver_config` (2D と同じ関数)、前段は行 13 と同じ | 行 13 と同じ。別 run の場の流用は `S3.warm_from_same_mesh` (index コピー、試験からだけ呼ぶ) | 行 13 と同じ | **なし → あり** |
| 15 | `opt/driver_sern._eval_op` (標準 → 緩レシピの梯子) | ⑤ SERN 2D | 行 13 と同じ (warm start は `operating_points[].warm_from`) | 行 13 と同じ (緩レシピは step 倍・cfl 半分) | 行 13 と同じ | 行 13 と同じ | **なし → あり** |

起動しない経路 (参考): `export/deliverables.py` は品質検査用に cell 変換するだけ。`intake/` は問題 YAML を作るだけ。`solver_density_cuda/tools/rerun_conditions.py`
は参照 run から新しい run を作る (`restart_field.py` + 場の検査) が起動はしない。case 配下のスクリプト (case/45 の `tw370_run.py`・`h370_run.py`・`cold_pair.py`
など) は `A.prepare_ns`・`A.run_staged_ns` を呼ぶので §1 の表の行 11 と同じ経路 (`cold_cfl.py prep` のように自前で config を書き換えるものは棚卸しの対象外)。

### 1.1 経路とプロファイル (2026-10-14、plan §5.1 #6 の後)

各経路の段を書き写したプロファイル (`solver_density_cuda/tools/recipes/recommended_settings.yaml`、`recipe_settings.py list`)。
`RECIPE.json` (宣言) はまだどの生成器も書かない (共通の組み立て関数で書かせるのは plan §5.1 #11)。手で作る run は
`recipe_settings.py declare` で書ける。

| 行 | 経路 | プロファイル (段) | 書き写しに使った代表の問題 |
| --- | --- | --- | --- |
| 1・5・6 | ③ベル CPG (`runner.main`・opt/driver・polish) | `bell_cpg_ns_sst` (soft → main) | case/40 `problem_bell_top.yaml` (opt/driver の CPG 評価) |
| 2・3・4 | ③ベル TP (`--stage`・`--staged`・`--continue-from`) | `bell_tp_ns_sst` / `bell_tp_euler` (soft → mid → main) | case/68 `problem_bell_tp_m_ns_l100_corr.yaml` / `problem_bell_tp_m_euler_l100.yaml` |
| 7 | ① 風洞 モード F (`runner_wt.main`) | `wt_modef_euler_cell` (soft → main) | case/41 `problem_m4.yaml` |
| 8 | ① 風洞 NS v1 (`runner_wt.prepare(euler=False)`) | なし (design の中に呼び出し元が無い) | — |
| 9 | ① 風洞 B10 の帰還 (`euler_loop`) | `wt_feedback_euler_node` (main だけ) | case/41 `problem_m4_b10.yaml` |
| 10 | ① axis-Mach Euler (`runner_axismach.main`) | `axismach_euler_tp` (soft → mid → main) | case/45 `problem_d155_euler_pin_G1_recal_mono.yaml` |
| 11 | ① axis-Mach NS (`prepare_ns` → `run_staged_ns`) | `axismach_ns_sst_tp` / `axismach_ns_sst_tp_b0` (S1_soft → S2_mid → main) | case/45 `problem_d155_ns_prod_coldmesh_tw300.yaml` / `..._tw370_cgref_h_b0.yaml` (main_line の B0)。`stages=ramp` は書き写していない |
| 12 | ① wall-driven (`runner_walldriven.main`) | `walldriven_euler_cpg` (soft → main) | case/41 `problem_m4_walldriven.yaml` |
| 13・15 | ⑤ SERN 2D (`runner_sern.main`・driver_sern) | `sern2d_finite_base_frozen_tp_sst` (warm_lam → soft → mid → main) | case/46 `problem_moo_frozen_tp_cycle3op.yaml` (既定の作動点)。warm start の適応段は書き写していない |
| 14 | ⑤ SERN 3D (`runner_sern3d.main`) | `sern3d_frozen_tp_sst` (warm_lam → soft → mid → main) | case/46 `problem_3d_prod_m6on_g2.yaml` |

## 2. 段の作り方の詳細 (書き写しの対象)

| 生成器 | 本段の値 (このコミット) | 前段の作り方 |
| --- | --- | --- |
| `R._solver_config` (CPG) | node・SLAU・`timeIntegration 11`・`blockDPLUR 1`・cfl 4・`nStepInner 5`・`convMethod 1`・`limiter 2`・SST (`katoLaunder 1`、`turbulentPrandtl`・`wallTreatmentSST` を書かない)・`nodeInletCornerWall` を書かない・`implicitRelax` を書かない | opt/driver の soft だけ (cfl 0.5・1 次、`nStepInner 5` のまま) |
| `R._solver_config_tp` (TP) | node・cfl `evaluate.cfl_main` (既定 4)・`implicitRelax 0.7`・`nStepInner 4`・`nodeInletCornerWall 1` (NS)・SST `turbulentPrandtl 0.9`・`output.level 1` | `TP_STAGES` (soft・mid: `convMethod 0`・`limiter 0`・cfl 0.5 / 1.0・`nStepInner 10`、2000 step) — §1・§1.2 どおりの唯一の生成器 |
| `W._config_euler` (cell) | **cell**・cfl 4・`nStepInner 5` | 文字列置換 (行 7) |
| `W._config_euler_node` | node・cfl 4 (axis-Mach は `--cfl`・`evaluate.cfl_main`)・`nStepInner 5`・`implicitRelax` を書かない | axis-Mach: `_first_order`・`_cfg_set` (行 10)、walldriven: 文字列置換 (行 12) |
| `W._config_sst_node` | node・cfl (axis-Mach NS は `evaluate.cfl_main` 既定 1.0)・`nStepInner 5`・`turbulentPrandtl 0.9`・`wallTreatmentSST 0`・`nodeInletCornerWall` を書かない | `A.ns_stage_configs` (行 11) / 文字列置換 (行 7・8) |
| `S._solver_config` | 離散化は `mesh.discretization` (**既定 cell**)・cfl `evaluate.cfl_main` (既定 4)・`nStepInner 5`・`limiterScaled`・`venkatK`・`pRef` = 外気圧・SST は `turbulentPrandtl` を書かない・`nodeInletCornerWall` は `mesh.node_inlet_corner_wall` のときだけ・`implicitRelax` は `evaluate.implicit_relax` のときだけ | 正規表現・文字列置換 (行 13) |

## 3. 範囲外として残したもの (#9 の後)

1. **euler_loop の継続チャンク** (§1 の表の行 9): 同じ格子への `interp_field.py`。AGENTS.md の規則 (同一メッシュは `restart_field.py`) に当たる。
   初期場が変わる修正なので、置き換えるなら #9 と同じ restart の検査を付ける。どの項目で扱うかは親が決める。
2. **opt/driver・opt/polish の CPG** (§1 の表の行 5・6): 同じ幾何から作り直した格子への `interp_field.py`。格子の生成が決定的なら同一メッシュで、
   `restart_field.py` が座標の一致を確かめられる。CPG のベルを今も使うか (生産は TP の行 3) と合わせて判断が要る。
3. **SERN の `restart_by_index`** (§1 の表の行 13〜15): 段終了ゲート (必要保存量・有限・ρ>0) が無い。化学種の照合と座標一致ノードの扱い
   (2026-09-04 の `interp_field.py` の誤写像の対策) は持っている。ゲートを足すかは plan §5.1 #11 で `restart_field.py` に寄せるときに合わせて決める。
4. **段の config の書き換えの方式** (正規表現・文字列置換): YAML の位置の方式 (`yaml_strict.replace_scalars`) に揃える (plan §5.1 #11)。#9 では config の文字列を変えていない
   (変更前とのバイト同一を `design/tests/run_launch_path_tests.py` で固定)。

**#19 の後 (2026-10-15)**: 1〜3 は段終了ゲート + `restart_field.py` (`runner_wt._staged_handoff`) に揃えた。1 は `euler_loop.start_cont_chunk`
(パスの格子の写しへ)、2 は `runner.continue_cpg` (継続元の変換済みの格子を写して使い、作り直さない。opt/driver の CPG の soft → 本段と
opt/polish の CPG のチャンク継続)、3 は `runner_sern.restart_by_index` (次の段の config を渡し、層流暖機 → SST で足す roK・roOmega の
初期値を前の段のゲートと分けて検査する: `runner_axismach.next_stage_init_problems`)。負例は `design/tests/run_launch_path_tests.py` (g)〜(i)。
残る同一格子の index コピーは `runner_sern3d.warm_from_same_mesh` (試験からだけ呼ぶ、段終了ゲートなし)。

## 4. 確かめたこと

- `design/tests/run_launch_path_tests.py` (forge の代わりに実際の HDF5 を書く偽の forge): runner_wt (Euler cell・NS node SST)・walldriven・SERN
  (層流暖機・soft・mid の梯子、CFL の梯子、warm start、`stages=main`・`none`) の段の config が変更前 (570790e8) とバイト同一、`stage_manifest.json` の段の並び、
  段の履歴、`restart_field.py` による初期場のビット一致と `wall_dist` の保持、restart の負例 (保存量 roe が無い・roUx に Inf・ρ<0・NS で roOmega が無い・
  節点数の違う res) で本段を起動しないこと。
