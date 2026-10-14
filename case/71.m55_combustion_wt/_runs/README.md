# case/71 の run の記録の写し (AWS A から回収、h5 は無い)

2026-10-15、plan [`tooling-nozzle-card-local-mesh`](../../../plans/active/tooling-nozzle-card-local-mesh.md) §5.1 #3b。
codex 諮問 [`2026-10-15-m55-ns-unsettled-diagnose.md`](../../../notes/reviews/2026-10-15-m55-ns-unsettled-diagnose.md) の「証拠不足」(全残差履歴・
stage_manifest・起動時の設定・停止規則の記録・終了状態・実効の設定・バイナリの精度・判定ツールの VERDICT・壁解像) に応えるための写し。
本体は AWS A の `~/forge-cardmesh/case/71.m55_combustion_wt/run_0002_ns`・`run_0001_euler` (保存場 `res_*.h5`・格子 `nozzle.h5` はそちらだけ)。

**実行のバイナリ** (両 run 共通、`RUN_PROVENANCE.txt` の forge_bin): FP64 (`flowFormat.hpp` の `flow_float`・`geom_float` = double) の
`/home/ubuntu/forge-integ-dflow2-fp64/solver_density_cuda/build/forge`、sha256 `3a9e347bfac51d003c4bfc9d384f8740b85f3058df9aa486aff4806a089231f8`、
mtime 2026-10-14T13:18:15Z、ツリーの HEAD 2fd36b0c (ソルバの源の最後の commit 75cb99e7 はビルドより前、ビルドより新しいソースのファイル 0 件)、
Release・sm_86、`FORGE_CUDA_BLOCKSIZE=128`。生成器 (run の作成・停止規則の見張り) は `~/forge-cardmesh` の b0d528bf (git_dirty 0 件)。

## 中身

| パス | 中身 |
| --- | --- |
| `SHA256SUMS_collect.txt` | 回収した元のファイル (A) の sha256 と、gzip の中身の sha256 (一致を確かめた) |
| `run_0002_ns/` | NS pass 0 の記録: `residual_history.csv.gz` (本段 0〜183,998 step の全行、内反復を含む)・`residual_history_S1_soft.csv.gz`・`residual_history_S2_mid.csv.gz`、`stage_manifest.json` (S1・S2 だけ: 本段の後の記録の前にドライバが OSError で落ちた)、`solverConfig.yaml` (本段、nStepOuter 200,000)、`bcondConfig.yaml`、`prepare_info.json`、`MAIN_STAGE_STARTED` (本段・S1・S2 の設定の sha256)、`RUN_PROVENANCE.txt`、`forge_launches.jsonl` (3 回の起動の cfg_fnv)、`CONVERGENCE_VERDICT.txt`、`stop_rule.json` (見張りの全記録; status RUNNING、行は 1,000〜182,000)、`stop_rule_band/` (帯 A の抽出の出力)、`run_0002_ns.driver.log` (ドライバ全体)、`forge_run.head300.log`・`forge_run.tail200.log`、`restart_field.log`・`run_case_stdout.log`・`MESH_QUALITY.txt`・化学種の記録・`delta_r_initial.*`、`_absent.txt` (`STOP_RULE_STATE.json` は無い・`stop_rule.json.tmp` は 0 byte) |
| `run_0002_ns/judge/` | 判定ツールの出力の原文 (A で 2026-10-15 に実行; 実行の台本 `run_0002_ns/judge_0002.sh`、最終の場の検査 `final_field_check.py`、実効の設定の照合 `config_regen_check.py`)。`stop_rule_series.csv.gz` は `check_quasisteady --series-csv` に渡した停止規則の系列 |
| `run_0002_ns/theta_band_ab.json` | θ_r の測定帯の A/B (台本 `run_0002_ns/theta_band_ab.py`)。帯 B の抽出の出力は `run_0002_ns/stop_rule_band_B_183000/` (`delta_r_edge.csv.gz` の y_b 列が帯、`delta_r_extract_summary.json`) |
| `run_0001_euler/` | Euler の参照の記録: 残差履歴 (soft・mid・本段・`--segment` の連結) の gz、`stage_manifest.json`・`solverConfig.yaml`・`bcondConfig.yaml`・`prepare_info.json`・`RUN_PROVENANCE.txt`・`forge_launches.jsonl`・`CONVERGENCE_VERDICT.txt`・`metrics.json`・`euler_driver.log`・`achieved_vs_target.csv.gz`・forge の log の先頭と末尾 |

## 判定の結果 (2026-10-15、AWS A、forge 0 step; 原文は `run_0002_ns/judge/`)

| 項目 | 結果 (ツールの原文の要点) |
| --- | --- |
| 区間 | 記録の `stage_manifest.json` には本段が無い (ドライバが本段の後の `_record(cfg_main, "main")` の前に OSError で落ちた) ので、`check_convergence --segment` をそのまま回すと **S1_soft → S2_mid (1 次) を判定してしまう** (`check_convergence_segment_asis.txt`: NOT CONVERGED、last step 5999)。本段を `StageManifest.add` で足した写し (`stage_manifest_with_main.json`; 本段の cfg_fnv `d903391a098b4612` は forge_launches の 3 回目と一致) では区間 1 = main だけ (convMethod 1 で S2 と別) → 判定区間は本段 0〜183,998 |
| `check_convergence --segment` (本段) | **NOT CONVERGED (still converging — run more steps)**: rms_ro 1.8 桁・roUx 2.4・roUy 2.4・roe 2.2・roK 2.9・roOmega 2.7・roY0/Y1 1.8 桁、全列 falling (`check_convergence_segment_with_main.txt`; 本段だけの通常判定 `check_convergence_main_only.txt` も同じ数字)。参考 (ツールの外、`residual_history.csv.gz` の outer_end の行から): 各列の最小は 157,795〜172,251 step (rms_ro 6.66e−7 @ 163,988)、直近 2 万 step (164k〜184k) の log10 の平均はその前の 2 万 step (144k〜164k) より 0.07〜0.44 桁高い (ro +0.15・roUx +0.32・roUy +0.44・roe +0.22・roK +0.07・roOmega +0.09) |
| `check_quasisteady --series-csv` (停止規則の系列 1,000〜183,000、tail 既定 0.4) | **DRIFTING**: r DRIFTING (drift 107.6 %/tail)、θ_r 0.4・0.7 TRANSIENT-UNSETTLED (drift 3.7・3.5 %/tail)、θ_r 0.95 STEADY (2.1 %/tail; 既定の許容 5 %)、ṁ_in・ṁ_out・出口コア M STEADY (`check_quasisteady_series.txt`)。系列の 183,000 の行は stop_rule.json に無い (保存の失敗) ので腕 A で計算し直した行 (同じ関数で 100k〜182k が記録と差 0.0 で一致) |
| 最終の場 res_183000 (`final_field_check.json`) | 非有限 0 (全 23 量)。ρ・P・T > 0 (最小 0.0193 kg/m³・2,693.7 Pa・202.8 K)、k ≥ 0 (k = 0 の 4,800 点は全て壁の行 j = 91)。静圧 ≤ 全圧 (`total_quantities.total_state`、TP-NASA9 frozen): **298 点で P > P0、最大 P/P0 = 1 + 5.8e−8** — 全て壁の行 (U = 0、T = 350 K) で、P0 = P になるはずの点の逆算の差の水準。ほか: P0 > Pt (1 + 1e−6) が 50,441 点 (壁から 27 点以上離れた主流〜軸、x −1〜8 m; P0/Pt の 90 % 点 1.0000055・99 % 点 1.00045・最大 1.00236)、T0 > Tt + 1 K が入口の列 (i = 0、壁から 19〜27 点) の 9 点 (T0 最大 1,253.4 K) |
| 壁解像 `check_wall_resolution --step 183000 --groups wall --weight area` | **PASS (壁解像)**: no-slip 壁は bcond の `wall` (physID 3、wall_isothermal 350 K) だけ。y1+ 平均 0.485・p99 0.821・最大 1.004 (x −1 m の入口の端)、y1+ > 1 の面積 0.0 %、評価できた面積 100 %。壁の出力 `res_wall_3_183000.h5` は全場と同じ step にある。**停止の step ではない** (停止規則は未到達) ので plan §6 V4 の判定ではない |
| 実効の設定 (`config_regen_check.json`) | **一致 (生成器どおり)**: 問題 YAML (`problem_m55_final.yaml`、sha256 f1ae7ee1…、手元の commit と同じ) から生成器で作り直した本段の solverConfig (nStepOuter を上限 200,000 に書き換えた後) が run の `solverConfig.yaml` と文字列で一致 (sha256 9aefd8aa… = MAIN_STAGE_STARTED)、S1・S2 の設定の sha256 が MAIN_STAGE_STARTED・prepare_info と一致、forge_launches の 3 回の cfg_fnv が S1・S2・本段と一致、bcondConfig・species_meta が一致 (bcond_fnv caea65b4…)。BC: 入口 inlet_Pressure Pt 3.5 MPa・Tt 1,250 K・Y0 0.9142/Y1 0.0858・k 1・ω 18,000、出口 outlet_statPress Ps 2,770.74 Pa (= spec.p_ambient)、壁 wall_isothermal Ts 350 K、軸 axis。停止規則は既定の要求 (問題 YAML に `evaluate.stop_rule` 無し → `explicit` false)、euler_ref はディレクトリだけの指定 (固定のスナップショットでない; 帯 A の抽出は当時の最新 res_12000 を読んだ = `stop_rule_band/delta_r_extract_summary.json`)。forge の log の先頭 (`forge_run.head300.log`) に `thermoFloat 1 (hybrid float Newton + double polish)`・輸送物性の表 (float) が出ている (既定のまま) |
| θ_r の測定帯の A/B (`run_0002_ns/theta_band_ab.json`) | 下の節 |
| Euler の参照 `run_0001_euler` の `--segment` | NOT CONVERGED (stalled/plateau)、区間 mid → main、2.1〜2.5 桁 (`check_convergence_segment_euler.txt`) |

## θ_r の測定帯の A/B (forge 0 step)

同じ保存場 (res_100000〜180000 の 1 万 step ごと + 181,000〜183,000、12 枚)・同じ Euler の参照 (`run_0001_euler/res_12000.h5`、
sha256 93cb88b0…)・同じ 3 断面 (x = 32.03・56.05・76.10 r_t) で、積分帯の外縁 y_b だけを変えた。A = `stop_rule.json` の固定帯 (本段の
最初の出力 res_1000 から抽出)、B = res_183000 から同じ関数 (`stop_rule.fix_band`、同じ抽出の条件) で抽出した帯。腕 A の再計算は
見張りの記録 (stop_rule.json の行) と 100k〜182k で差 0.0。D_θ は主の窓 = 100k〜183k (12 点、前半・後半の平均の step の差 52,286) を
2 万 step あたりに換算。

| 断面 | D_θ A [%/2 万 step] | D_θ B | 閾値 | 判定 (事前に固定) |
| --- | --- | --- | --- | --- |
| x/x_w 0.4 | 0.969 | 2.332 | 0.05 | 両腕とも超過 |
| 0.7 | 0.912 | 1.960 | 0.05 | 両腕とも超過 |
| 0.95 | 0.538 | 1.748 | 0.05 | 両腕とも超過 |

**判定: 3 断面とも両腕で D_θ > 0.05 % → 「帯の選び方だけで未到達を説明できる」を棄却。** 副の見方 (判定に使わない) も同じ側:
1 万 step の等間隔 9 点で A 0.55〜1.00・B 1.81〜2.44、最小二乗の傾きで A 0.53〜0.95・B 1.72〜2.30、腕 A の 1,000 step ごとの見張りの系列
(100k〜182k、83 点) で 0.56〜1.02 (%/2 万 step)。θ_r は両腕とも 100k〜183k で単調に増えている (例 0.4: A 0.039617 → 0.041133、B 0.046519 → 0.051045 r_t)。

帯そのもの: 後期の帯 B は A より 1.40〜5.52 倍 (使用域 [0, 80.07] r_t の中央値 2.82) 広く、後期の縁 y* (= y_b,B / 1.25) は使用域の全列
(100 %) で帯 A の外にある (3 断面の y_b: A 0.492・0.893・0.970、B 0.753・1.419・1.800 r_t)。つまり帯 A の θ_r は後期の場では境界層の
一部しか積分していない。これは測定の妥当性の記録で、上の判定 (帯を広げても未到達は消えない) とは別。

## 容量 (AWS A、2026-10-15 08:13 JST; 原文 `run_0002_ns/judge/disk.txt`)

| 項目 | 値 |
| --- | --- |
| 空き容量 (`/`) | 16,231,657,472 B (16 GB、使用 93 %) |
| 全場の出力 1 枚 `res_183000.h5` (441,600 節点、FP64、23 量) | 100,605,412 B |
| 壁の出力 1 枚 `res_wall_3_183000.h5` (4,800 点) | 1,170,416 B |
| 出口の出力 1 枚 `res_outlet_2_183000.h5` | 28,896 B |
| 格子 + 初期場 `nozzle.h5` | 187,637,372 B |
| `residual_history.csv` (本段 183,999 step、内反復の行を含む 551,997 行 + 見出し) | 109,468,246 B (gzip 16,035,021 B) |
| `forge_run.log` (本段) | 19,755,555 B |
| run dir の合計 | run_0002_ns 2,799,718,212 B (全場 22 枚・壁/出口の出力 183 枚ずつ)、run_0001_euler 1,640,406,895 B (Euler の全場 1 枚 106,077,312 B) |

1,000 step ごとに全場を書くと 1 万 step で約 1.0 GB、壁の出力は 1 万 step で約 12 MB。
