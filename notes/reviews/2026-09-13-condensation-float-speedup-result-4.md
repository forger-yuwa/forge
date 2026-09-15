# codex レビュー: condensation-float-speedup (result)

- **plan**: [`plans/active/condensation-float-speedup.md`](../../plans/active/condensation-float-speedup.md)
- **stage**: `result` (diff base `3b77cc4b`)
- **date**: 2026-09-13
- **commit**: `01b5034c` (feature/perf-3d-speedup)
- **codex**: effort `high`, 6.3 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M2/m1
- **focus**: 凝縮 float 化の result レビュー 4 回目。3 回目 (NO-GO M1/m2) の採用: A10G 実測 (dry 33.8 不変 / 起動 89.4→44.0 [目標 ≤45 達成] / 発達場 98.0→52.3 [目標 ≤50 に 4.6 % 未達, 残りは followups F-cf3])、h0 二重 VERDICT、閾値の文書同期。証拠は case/16.nozzle_wys/_aws_perf_evidence/cond/ (bench ログ・nsys CSV・BINARY_SHA256)、3D 300 step の場比較と分布統計は §5.1 #13。発達場の目標未達 (4.6 %) を残作業に戻した上で accepted に移してよいか判定してください。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**

`F-cf3` への転記と `h0` 二重判定の修正は確認できました。  
しかし、回収済みベンチログは最終バイナリのものではなく、44.0／52.3 ms/step を裏付けません。  
さらに最終変更には、計画対象外の二相音速の変更と、その誤差を見逃す試験上の穴があります。

1. **Major — 最終版の性能・3D回帰を、別バイナリの証拠で完了扱いしている。**

   [plan §5.1 #13](/home/sano/work/forge-perf/plans/active/condensation-float-speedup.md:193) は `forge_condf2` の100 step×3交互計測としています。しかし、`case/16.nozzle_wys/_aws_perf_evidence/cond/` の18本を再集計すると、次の結果です。

   | 対象テンプレート | 基準 `forge_merged` | 回収された新版 `forge_condf` |
   |---|---:|---:|
   | `run_0400_perf_baseline` | 33.76–33.77 | 33.68–33.71 |
   | `run_0420_perf_cond_baseline` | 89.39–89.40 | **45.61–45.63** |
   | `run_0421_perf_regress_cond3d` | 98.01–98.02 | **55.12–55.13** |

   単位は `ms/step`。新版ログのSHAはすべて `eb9d571e995ada5e…` です。[起動ログ:1](/home/sano/work/forge-perf/case/16.nozzle_wys/_aws_perf_evidence/cond/run_0420_perf_cond_baseline_bench__f_072807_r1_n100.log:1)、[発達場ログ:239](/home/sano/work/forge-perf/case/16.nozzle_wys/_aws_perf_evidence/cond/run_0421_perf_regress_cond3d_bench__f_072807_r1_n100.log:239) で直接確認できます。

   一方、[BINARY_SHA256.txt:3](/home/sano/work/forge-perf/case/16.nozzle_wys/_aws_perf_evidence/cond/BINARY_SHA256.txt:3) の `forge_condf2` は `96c037c11bae5c48…`。最終版名の nsys CSV はカーネル時間合計 **51.8386 ms/step** を裏付けますが、これは壁時計の交互計測を代替しません。

   また、指定フォルダには300 stepの場、`CMP_VERDICT`、分布統計の集計結果がなく、ローカルに `run_0421_perf_regress_cond3d/` もありません。したがって「30/30 PASS」「分布統計一致」は今回独立検証できません。

   **対案:** #13/#15と§8の完了扱いを戻し、最終版の交互計測ログ、基準・新版の300 step比較資料、入力HDF5を含むSHAを回収してください。**現状では「残り4.6%」という移管対象の数値自体が確定していません。**

2. **Major — 二相音速への表微分導入は、Newton高速化だけではない。接続点で既定挙動が変わる。**

   [dependentVariables_d.cu:219](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/dependentVariables_d.cu:219) は、音速に使う従来の
   `(L(T+0.1)−L(T−0.1))/0.2`
   を表の解析微分に置換しています。これは反転の収束速度だけでなく、`gamma`・`sonic`、さらに陰解法に渡す係数を変えます。

   特に [condensationProperties_d.cuh:252](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/condensationProperties_d.cuh:252) は液相エンタルピーを373.15 Kでクランプします。コードの係数・表構築・float32評価をPythonで再評価した合成状態では、以下を得ました。

   - `T=373.15 K`：従来の `dL/dT ≈ −218.609`、表微分 `≈ +1889.944 J/(kg·K)`。
   - `g=0.02`、`cp_allvap=1020`、`Rmix=290 J/(kg·K)`：二相音速は **379.906→383.023 m/s、約0.82%差**。
   - 二相比熱の差は約4.1%。コードコメントの「相対誤差 ≤1e-8」は成立しません。

   これはGPU runの観測値ではなく、実装式の独立再評価です。[単体試験:50](/home/sano/work/forge-perf/solver_density_cuda/tests/unit/test_cond_float.cpp:50) は373.15 K近傍を除外し、参照微分の刻みも本番旧式と異なります。既存試験のPASSではこの変更を検証できません。[現在仕様:890](/home/sano/work/forge-perf/methods/condensation.md:890) の「二相音速はdoubleのまま」とも不整合です。

   **対案:** 本planでは**音速用微分を従来のdouble差分に戻す**ことを推奨します。表微分はNewtonの傾きに限定し、修正後の最終版で性能と回帰を再確認してください。

3. **Minor — 残作業表と正式な検証条件が、最新の完了主張に追随していない。**

   [§5.1 #15:195](/home/sano/work/forge-perf/plans/active/condensation-float-speedup.md:195) は旧SHAを「最終バイナリ」と呼び、依然「A10Gは未計測」。[§6:202](/home/sano/work/forge-perf/plans/active/condensation-float-speedup.md:202) は短期回帰をAWS実施としていますが、#5ではローカル予備・AWS正式取得未了のままです。[§8:231](/home/sano/work/forge-perf/plans/active/condensation-float-speedup.md:231) の完了チェックと両立しません。

   `F-cf1/F-cf2/F-cf3` は [followups §5.1](/home/sano/work/forge-perf/plans/active/condensation-followups.md:53) に存在します。ただし、転記は最終版の未検証項目を消す根拠にはなりません。

   **対案:** 各検証にコミット・バイナリSHA・実施環境・成果物を対応付け、ローカル結果を正式採用する項目と未完了項目を明示してください。#13–#15、§6、§8、`plans/README.md` を同じ状態に揃える必要があります。

前回修正のうち、以下は裏付けられました。

- double sourceのセル本体とdouble clamp本体は変更前と同一。`case/16.nozzle_wys/run_0458_frozen1step_node2d_cond/{merged,new_cf0}/res_1.h5` のT・P・凝縮診断6量もビット一致しました。
- `case/16.nozzle_wys/run_0460_condf_phys_node2d_cond/{merged,condf}/` を再判定すると、双方 **`NOT CONVERGED (stalled/plateau)`**、`pmax/machmax` は **`OVERALL: ALL STEADY`**、場比較は **32/32 PASS**。未収束状態の回帰合格という記述は適切です。
- `case/34.arthur_n2_nozzle/run_0100_merge_regress_cell_air/condf/` と `run_0103_condf_regress_n2_cell/condf/` のh0偏差は **428.542／482.370 J/kg**。再実行でも双方 **`ABSOLUTE_VERDICT: FAIL`／`NONDEGRADATION_VERDICT: PASS`** となり、二重表示の修正を確認しました。

**推奨は、`in_progress` を維持することです。** 音速の変更を戻し、その最終版の性能・回帰証拠を揃えてから再レビューしてください。現状で `accepted` への移動は支持しません。

ファイルは変更しておらず、指摘は **plan未反映**です。既存host試験バイナリはALL PASSでしたが、device試験はCUDA driver/runtime不整合で再実行できませんでした。

指摘数: Critical 0 / Major 2 / Minor 1
