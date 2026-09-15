# codex レビュー: condensation-float-speedup (result)

- **plan**: [`plans/active/condensation-float-speedup.md`](../../plans/active/condensation-float-speedup.md)
- **stage**: `result` (diff base `3b77cc4b`)
- **date**: 2026-09-13
- **commit**: `a1cd6afa` (feature/perf-3d-speedup)
- **codex**: effort `high`, 5.4 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M1/m2
- **focus**: 凝縮 float 化の result レビュー 3 回目。2 回目 (NO-GO M4/m2) の採用: 範囲外セルの dry 判定を旧 double 飽和圧で行う修正 (反例を device 試験に追加)、周期保存と cell h0 絶対ゲートを『未達 (既存挙動)』として condensation-followups F-cf1/F-cf2 に転記、最終バイナリ SHA 8c7096927b03b5d6 での dry/起動/発達場の再計測、各ラベルの VERDICT ファイル保存、文書同期。A10G は未計測で §8 に明記。この状態を『ローカル検証完了・A10G 保留』として accepted に移してよいか、それとも in_progress のまま A10G を待つべきかを判定してください。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**

最終バイナリでのローカル高速化と、前回の範囲外 dry 判定修正は裏付けられました。  
ただし A10G の性能目標は依然として完了条件であり、「ローカル検証完了・A10G 保留」での `accepted` 移動は支持しません。  
保存ゲートの既存不具合と高速化の回帰は区別できていますが、判定表示・文書同期には残件があります。

1. **Major — 本来の性能目標が未検証で、完了条件を満たしていない。**

   [plan §8:230](/home/sano/work/forge-perf/plans/active/condensation-float-speedup.md:230) は、A10G の発達場 `≤50`、起動区間 `≤45 ms/step` を完了条件として残しています。§5.1 #13/#15 に未計測と明記したことは適切ですが、条件の達成にはなりません。

   実バイナリの SHA は `8c7096927b03b5d68cf836ffc707ee3d3f1bcbb7e37dd0e10ed4a6343a718190`。同 SHA のログと基準ログを照合した結果は次のとおりです。単位は `ms/step`、各3反復です。

   | 証拠ディレクトリ（リポジトリ相対） | 基準 | 最終版 |
   |---|---:|---:|
   | `case/16.nozzle_wys/run_0455_perf_big_norcm_bench/` | 62.24–62.73 | 61.76–62.54 |
   | `case/16.nozzle_wys/run_0459_perf_cond3d_local_bench/` | 194.79–195.97 | 85.92–88.06 |
   | `case/16.nozzle_wys/run_0463_perf_cond3d_local_dev_bench/` | 215.20–220.53 | 110.45–111.91 |

   対象は各ディレクトリの `mF_161240_r{1,2,3}_n100`／`fF_161240_r{1,2,3}_n100`、成果物は `bench_*.log`。新版ログの反転失敗警告は0です。**ローカル改善は実在しますが、A10G の絶対時間を証明しません。**

   **対案:** `in_progress` を維持し、A10G で同一最終バイナリによる dry・起動・発達場の交互計測を完了する。目標を満たさなければ、追加最適化を残作業に戻してください。

2. **Minor — `h0` の最終 VERDICT が、絶対ゲートと非劣化ゲートを識別できない。**

   [cond_axis_h0.py:97](/home/sano/work/forge-perf/solver_density_cuda/tools/cond_axis_h0.py:97) は `--ref` 指定時に絶対上限超過を失敗条件から外し、[112行](/home/sano/work/forge-perf/solver_density_cuda/tools/cond_axis_h0.py:112) で単一の `VERDICT: PASS` を返します。再実行でも次を確認しました。

   | run／ラベル | 新版／基準の偏差 [J/kg] | 絶対上限300 | 現ツール |
   |---|---:|---|---|
   | `case/34.arthur_n2_nozzle/run_0100_merge_regress_cell_air/condf`／`merged` | 428.542／434.125 | 超過 | `PASS` |
   | `case/34.arthur_n2_nozzle/run_0103_condf_regress_n2_cell/condf`／`merged` | 482.370／481.952 | 超過 | `PASS` |

   両比較とも末尾は `STEADY`。途中に `EXCEED` は表示されるため数値を隠してはいませんが、VERDICTだけを読むと絶対保存基準も通ったように見えます。

   **対案:** `ABSOLUTE_VERDICT: FAIL` と `NONDEGRADATION_VERDICT: PASS` を明示し、保存ファイルも再生成する。既存誤差の修正は、転記済みの [followups F-cf2:54](/home/sano/work/forge-perf/plans/active/condensation-followups.md:54) で扱う整理を支持します。

3. **Minor — 「文書同期済み」に対し、前回指摘された不整合が残っている。**

   根拠は以下です。

   - [plan:153](/home/sano/work/forge-perf/plans/active/condensation-float-speedup.md:153) は蒸発 `λ≤1e-5`。実試験は [test_cond_float.cpp:168](/home/sano/work/forge-perf/solver_density_cuda/tests/unit/test_cond_float.cpp:168) の**絶対誤差 `≤5e-4`**です。
   - [methods/condensation.md:908](/home/sano/work/forge-perf/methods/condensation.md:908) は蒸発端 `|Δsj|dt≤1e-6`。実試験は [device試験:148](/home/sano/work/forge-perf/solver_density_cuda/tests/unit/test_cond_float_device.cu:148) の `1e-4` です。
   - [solver-settings.md:174](/home/sano/work/forge-perf/procedures/solver-settings.md:174) は平衡形の EOS を「常に double」と説明します。しかし [dependentVariables_d.cu:140](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/dependentVariables_d.cu:140) が別扱いするのは `condEquilibrium==2`。`==1` の湿潤 TP セルはハイブリッド反転へ進みます。
   - [plan変更ログ:242](/home/sano/work/forge-perf/plans/active/condensation-float-speedup.md:242) に「場はノイズ床内」が残っています。`run_0460` の `rog_0` 差は `1.31e-5`、ノイズ床は `5.71e-6`。合格理由は絶対基準です。

   **対案:** §4.3・現在仕様・試験の閾値と分岐を統一し、変更ログにも訂正を明示する。完了するまで §5.1 #14 を未完了に戻してください。

実装については、範囲外の旧 double 飽和圧による判定と委譲を [condensationSourceKernels_d.cuh:353](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/condensationSourceKernels_d.cuh:353) で確認しました。double source のセル本体と double clamp 本体は変更前と同一です。`case/16.nozzle_wys/run_0458_frozen1step_node2d_cond/{merged,new_cf0}/res_1.h5` の T・P・凝縮診断6量もビット一致しました。

`case/16.nozzle_wys/run_0460_condf_phys_node2d_cond/{merged,condf}/` は、再実行で双方 **`NOT CONVERGED (stalled/plateau)`**、`pmax/machmax` は **`OVERALL: ALL STEADY`**、場比較は **32/32 PASS**でした。これは未収束状態の回帰合格です。[case README の run 一覧](/home/sano/work/forge-perf/case/16.nozzle_wys/README.md:75) の訂正後の説明を支持します。

周期箱 `case/09.Taylor-Green/run_0053_cond_periodic_conservation/` の20 step後のモーメント収支も、新版約 `1.111e-3`、基準約 `1.143e-3` を再現しました。保存ゲートは未達ですが、[F-cf1:53](/home/sano/work/forge-perf/plans/active/condensation-followups.md:53) への転記は確認でき、前回の「未転記」は解消しています。

**推奨は、`in_progress` を維持して A10G 実測を待つことです。** 待機中に判定表示と文書を修正し、性能ゲートの結果を得てから完了判定してください。

ファイルは変更しておらず、指摘は **plan 未反映**です。device 単体試験の再実行は CUDA driver/runtime 不整合で実行不能でした。今回、GPU単体試験の PASS を独立再現したとは扱っていません。

指摘数: Critical 0 / Major 1 / Minor 2
