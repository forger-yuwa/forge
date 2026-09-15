# codex レビュー: condensation-source-limiter-steady (result)

- **plan**: [`plans/active/condensation-source-limiter-steady.md`](../../plans/active/condensation-source-limiter-steady.md)
- **stage**: `result` (diff base `f1e40f41`)
- **date**: 2026-09-16
- **commit**: `0dac3473` (feature/sern-design)
- **codex**: effort `high`, 7.3 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M3/m2
- **focus**: 1 回目 result レビュー (NO-GO M1–M7/m1) の採否が §6.1 の記載どおりに実装・検証へ反映されているか (蒸発の一様 ṙ 形, condEquilibrium 0 限定, volumePartial, 全クランプ診断, 1 セル固定点/密度変化/多分散の単体, 凝縮固有量 series STEADY, node/cell × cfl, eq=1 回帰, Wys series, Arthur 反復基準)。§9 の 2026-09-16 変更ログの主張を run (case/44 run_0137–0146, case/34 lim1b*, case/16 run_0470 lim1b/lim0b/t*, case/09 run_0060/0061) と単体テスト tests/unit/test_cond_limiter_steady.cu で突き合わせ、accepted にしてよいか

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

# 判定: **NO-GO**

主要な実装修正と最終場の差分値は確認できました。  
ただし、受入条件の未達、継続計算の入力不備、必須試験の欠落が残っています。  
前回指摘の「全反映済み・検証完了」という結論は成立しません。

## 指摘一覧

### 1. **Major — 凝縮固有量が `STEADY` でも、残差の受入条件を満たしていない**

**根拠:** 対象 run に `check_convergence.py` を再実行し、保存済み判定と同じ **`NOT CONVERGED (stalled/plateau)`** を確認しました。

[plan §6](/home/sano/work/forge/plans/active/condensation-source-limiter-steady.md:133) の代替ゲートは、**全モーメント残差の3桁低下・固有量の `STEADY`・補正量ゼロの全条件**を要求しています。

| case/44 の run | `rog / roQ2 / roQ1 / roQ0` の低下桁数 |
|---|---|
| [run_0137：node CFL 2](/home/sano/work/forge/case/44.vitiated_air_wt/run_0137_va3_M4.19_Lc8_noneq_inletTt_lim1b_cfl2/CONVERGENCE_VERDICT.txt) | **3.4 / 3.1 / 2.8 / 2.7** |
| [run_0141：cell CFL 2](/home/sano/work/forge/case/44.vitiated_air_wt/run_0141_va3_M4.19_Lc8_noneq_inletTt_cell_lim1_cfl2/CONVERGENCE_VERDICT.txt) | **2.2 / 1.6 / 1.4 / 1.1** |
| [run_0146：cell CFL 0.5](/home/sano/work/forge/case/44.vitiated_air_wt/run_0146_va3_M4.19_Lc8_noneq_inletTt_cell_lim1_cfl05_cont/CONVERGENCE_VERDICT.txt) | **2.3 / 2.0 / 1.7 / 1.5** |

これらの `cond_series.csv` は、指定の `--drift 0.002 --osc 0.002` で **`ALL STEADY`** を再現しました。しかし、それで残差条件は代替できません。特に cell の `rms_roe` 最終値は **26.8 / 6.27** で、plan が例示する node の床 **0.42** とも異なります。

**対案:** §5.1 の回帰を未完了に戻し、残差床の原因を調べて受入ゲートを満たしてください。現状の比較は「最終保存場の差が許容範囲内」と報告し、「定常解の一致を検証完了」としないこと。

### 2. **Major — `run_0146` は凝縮場を引き継いでおらず、「72000 step 継続」は誤り**

**根拠:**  
`case/44.vitiated_air_wt/run_0146_va3_M4.19_Lc8_noneq_inletTt_cell_lim1_cfl05_cont/` の入力・初期出力を照合しました。

- `nozzle.h5` には `rog_0` 最大 **0.001173936**、`roQ0_0` 最大 **7.821119e13** が存在。
- 一方、`roY1`／`Y1` が欠落。
- `res_0.h5` では **`Y1=0`、4モーメントすべて0**。
- 初期 `condClampCorr_0` 最大 **0.016189324**。引き継いだ液相が初期クランプで削除されています。

欠落した化学種は [variables.cpp:700](/home/sano/work/forge/solver_density_cuda/variables.cpp:700) で第1種以外を0に初期化され、[condensationTransport_d.cu:84](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationTransport_d.cu:84) の `rog≤roY_w` で液相が消えます。

さらに [interp_field.py:90](/home/sano/work/forge/solver_density_cuda/tools/interp_field.py:90) は、転送先に存在しない凝縮モーメントを新設しますが、**存在しない `roY*` は新設しません**。今回の入力不備を生む経路がコードにもあります。

**対案:** 化学種を含む全保存量を引き継ぎ、起動直後に入力と `res_0.h5` の保持を確認して再検証してください。現在の run は「流れを引き継ぎ、化学種・凝縮場を再発達させた48000 step」と訂正すべきです。

### 3. **Major — 前回 M7 と全クランプ診断の必須試験が、まだ完了していない**

**根拠:**

- [plan:142](/home/sano/work/forge/plans/active/condensation-source-limiter-steady.md:142) が要求する **node/cell × `condFloat 0/1` ごとの CFL 比較**に対し、cell の `run_0141/0142/0146` はすべて `condFloat` 省略＝1。**cell × `condFloat 0` がありません**。
- 追加された固定点試験は [test_cond_limiter_steady.cu:265](/home/sano/work/forge/solver_density_cuda/tests/unit/test_cond_limiter_steady.cu:265) で source と update limiter を呼びますが、**後段の `cond_realizability_clamp` を通していません**。
- 同試験の合否は状態の比較と `theta` が中心で、[同:279](/home/sano/work/forge/solver_density_cuda/tests/unit/test_cond_limiter_steady.cu:279) に**未加工残差・全硬クランプ補正量の判定がありません**。
- 更新試験では `corQ` を確保して渡すだけで、[同:118](/home/sano/work/forge/solver_density_cuda/tests/unit/test_cond_limiter_steady.cu:118) で読み戻して検証していません。

**対案:** 完了扱いを解除し、cell の両精度比較と、上限超過・各 Q の負値・消滅を含む**実際の更新全経路**の試験を追加してください。固定点試験には、残差と補正量の許容値を明示した判定を入れるべきです。

### 4. **Minor — 旧経路の補正診断が「このステップの量」になっていない**

**根拠:** 診断のリセットは新 limiter の [condensationUpdateLimiter_d.cuh:80](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationUpdateLimiter_d.cuh:80) にあります。旧経路は [condensationTransport_d.cu:398](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationTransport_d.cu:398) を通り、実現可能性クランプの補正を累積し続けます。

実測でも `case/16.nozzle_wys/run_0470_limiter_regress_wys/lim0b/` の `condClampCorr_0` 総和は、

- step 4000：**0.025118401**
- step 24000／48000：**0.025119125**

となり、起動時の補正が最終出力まで残ります。新経路との診断比較を誤らせます。

**対案:** ステップ開始時に共通でリセットし、旧更新の floor も記録するか、旧経路ではこの診断を未対応として明示してください。

### 5. **Minor — 現在仕様・台帳・残作業表に前回の誤記が残っている**

**根拠:**

- [case/44 README:598](/home/sano/work/forge/case/44.vitiated_air_wt/README.md:598) は、依然として `run_0131` を「正本」と記載。
- [methods/condensation.md:833](/home/sano/work/forge/methods/condensation.md:833) も、未収束の `run_0131` との「一致」を現行検証の根拠に使用。
- [followups F-cf7](/home/sano/work/forge/plans/active/condensation-followups.md:59) は「未実装・未着手」のまま。
- [対象 plan §5.1](/home/sano/work/forge/plans/active/condensation-source-limiter-steady.md:125) では、登録済みの F-cf8 が未完了表示で、実際に未達の検証が完了扱い。
- Wysłouzil の記録で使う [compare_condfix.py:79](/home/sano/work/forge/case/16.nozzle_wys/compare_condfix.py:79) の判定基準は、指定された `check_quasisteady.py` の0.2%基準とは別です。

**対案:** 現行の根拠 run と未達条件を正本へ同期し、正式ツールの判定を保存してください。`methods/index.md` は既存ファイルの更新なので、変更不要です。

## 確認できた修正・実測

前回指摘がすべて未対応という評価ではありません。

- **蒸発の一様 ṙ 形、`condEquilibrium==0` 限定、部分体積によるソース積分**は実装を確認。case/09 `run_0060/0061` の seam 比も、旧 **2／4／8** → 新 **0.99999976～1.0000001** を再現しました。
- case/44 の最終場差は記載どおり。`0139−0137` は g 相対 L1 **0.0018964**、最大温度差 **0.21033 K**。`0146−0141` は **0.00067980／0.11574 K**。これは上記の収束条件未達とは別の確認です。
- Arthur `lim1b*` の **28/28 PASS** と、N2 cell の反復間 `ro` 差 **2.42e-4～4.63e-4** を再現。ただし、この PASS は場の回帰比較で、収束 VERDICT は **NOT CONVERGED** です。
- Wysłouzil `lim1b/lim0b` は、報告量を抽出して正式な `check_quasisteady.py --series-csv` に渡し直しても、0.2%基準で両方 **ALL STEADY**。g の正規化最大差 **1.1973e-5** も再現しました。
- eq=1 のモード差は反復差以内。単独速度計測の約 **+6%** もログと整合します。今回確認した全保存場・全残差列に NaN/Inf はありませんでした。

## 推奨

**`in_progress` を維持し、入力引き継ぎを修正したうえで、残差ゲートと不足試験を満たしてから再レビューしてください。** 現時点で `accepted` へ移すべきではありません。

レビュー範囲は `f1e40f41...0dac3473`。CUDA 単体試験の再実行は driver/runtime 不整合により実行不能で、計算結果の不合格とは区別しました。ファイルは変更しておらず、本指摘は **plan 未反映**です。

指摘数: Critical 0 / Major 3 / Minor 2
