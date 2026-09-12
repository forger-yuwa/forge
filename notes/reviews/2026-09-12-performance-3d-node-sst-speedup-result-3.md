# codex レビュー: performance-3d-node-sst-speedup (result)

- **plan**: [`plans/active/performance-3d-node-sst-speedup.md`](../../plans/active/performance-3d-node-sst-speedup.md)
- **stage**: `result` (diff base `0512823d`)
- **date**: 2026-09-12
- **commit**: `f55918af` (feature/perf-3d-speedup)
- **codex**: effort `high`, 5.5 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M2/m3

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**NO-GO**

高速化は生ログで確認でき、提示された回帰も現行比較ツールでは PASS を再現しました。  
ただし、長時間比較の証拠と未達事項の管理が不足しており、`accepted` 移行の根拠はまだ揃っていません。確認した差分から Critical 相当の数値実装誤りは見つかっていません。

1. **Major — 長時間比較の必須証拠が不足したまま、残作業 #11 を完了扱いにしている**

   [plan:113](/home/sano/work/forge-perf/plans/active/performance-3d-node-sst-speedup.md:113) は両者の準定常確認を必須とし、116 行では「未収束なら比較を確定しない」としています。一方、[残作業 #11:160](/home/sano/work/forge-perf/plans/active/performance-3d-node-sst-speedup.md:160) は「済」です。

   再確認した結果は次のとおりです。

   - `case/16.nozzle_wys/_aws_perf_evidence/run_0416_final_sweep4/`：`check_convergence.py` は **NOT CONVERGED**。`rms_roUy=4.09e-10`、`rms_roK=7.09e-8` が plateau。残差履歴に NaN/Inf はありません。
   - 同ディレクトリの `check_quasisteady.py --quantity pmax,machmax`：**NO mesh / NOT ALL STEADY**。保存された準定常 VERDICT も見当たりません。これは非定常性の証明ではなく、証拠不足です。
   - 基準 `run_0234` の残差履歴・VERDICT は今回のツリーにありません。`_aws_ic/run_0234_res_12000.h5` 一枚では代替できません。

   [壁圧判定ファイル](/home/sano/work/forge-perf/case/16.nozzle_wys/_aws_perf_evidence/run_0416_final_sweep4/WALL_PP0_SERIES_VERDICT.txt:1) の輪郭壁の末尾差 **2.2118e-4** は CSV から再現しました。ただし、これは独自の隣接二枚判定であり、要求されたツールによる両者の定常性確認とは別です。

   **対案:** #11 を未完了へ戻し、基準の残差履歴、両者の準定常判定出力、再計算可能な入力・時系列を回収する。目的と完了条件を「有限時間継続の非劣化」に統一し、収束解保証は未達として明記してください。

2. **Major — 前回指摘された node 軸拘束の回帰が残作業から抜けている**

   [timeIntegration_d.cu:1001](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1001) は、今回 `cached` に応じて軸の拘束行処理を変えています。

   [plan:180](/home/sano/work/forge-perf/plans/active/performance-3d-node-sst-speedup.md:180) の軸対称回帰は **cell** です。追加された case/39 は node 周期・line-implicit であり、node の軸拘束分岐を通りません。これは前回レビュー M4 の対案にも含まれていましたが、§5.1 に未完了項目として残っていません。

   なお、`case/39.periodic_hills/run_0100_perf_regress_node_lineimplicit_dualtime/` の比較は再実行でも **PASS 18/18** でした。この対応自体は確認できています。

   **対案:** node 軸対称の回帰を §5.1 に戻す。最終版で軸の `roUy/dq_roUy` 拘束、全保存量残差、近軸の場を基準版と比較し、残置する対角キャッシュの off/on も確認してください。

3. **Minor — 温度反転試験のドリフト集計が、依然として異なる量・単位を混在させている**

   [test_thermo_float.cpp:72](/home/sano/work/forge-perf/solver_density_cuda/tools/test_thermo_float.cpp:72) は純 float 反転の **絶対温度差 `[K]`** を `maxdrift` に入れています。同じ変数へ 93 行でハイブリッド再格納の **相対温度差** を入れ、108 行で相対量 `maxdriftD` と比較しています。

   float 保存量の往復処理は追加されていますが、表示される `drift/T` と合否の集計は正しく分離されていません。

   **対案:** 純 float の絶対ドリフトと、ハイブリッド・double の相対ドリフトを別変数にする。後二者を同じ単位で比較し、修正版の試験結果で文書の数値を更新してください。

4. **Minor — 回帰ツールが場の欠落・形状不一致を黙って除外し、PASS にできる**

   [perf_regress.py:91](/home/sano/work/forge-perf/solver_density_cuda/tools/perf_regress.py:91) は比較先にない場を、93 行は形状の違う場を `continue` で除外します。さらに、比較失敗時も非ゼロ終了コードを返しません。

   ファイルを作らないメモリ上の再現確認で、基準八場から比較先の `ro` を欠落させ、`P` の形状を変えると、残った **6/6 場が合格**しました。

   **対案:** 必須場の欠落・形状不一致・非有限値を明示的な FAIL にし、一件でも失敗したら非ゼロ終了する。今回確認した PASS を直ちに否定する問題ではありませんが、「全場」の検証基盤として修正が必要です。

5. **Minor — 現在仕様と報告値に、修正済みとされた不整合が残っている**

   - [plan:181](/home/sano/work/forge-perf/plans/active/performance-3d-node-sst-speedup.md:181) と [case README:73](/home/sano/work/forge-perf/case/16.nozzle_wys/README.md:73) の凝縮回帰「全場 ≤1.6e-5」は誤りです。`run_0456_perf_regress_node2d_cond/base_r1` 対 `final` は、`h0=2.09e-5`、`vis_turb=1.08e-3`。現行閾値では **PASS 32/32** ですが、報告値は直っていません。
   - [performance.md:66](/home/sano/work/forge-perf/methods/architecture/performance.md:66) と plan §4.3 は絶対基準を記述する一方、実装と §6 は「絶対基準 **or** 2×ノイズ床」です。
   - [recommended-settings.md:35](/home/sano/work/forge-perf/procedures/recommended-settings.md:35) に「同一収束」が残っています。同ファイル47行と `solver-settings.md` は研磨を一段と説明していますが、実装は最大三段です。
   - [methods/index.md:15](/home/sano/work/forge-perf/methods/index.md:15) は表の途中へ箇条書きを挿入し、表を分断しています。

   **対案:** 数値・採否基準・研磨段数を最終仕様へ統一し、「同一収束」を未収束の有限時間比較という記述へ直す。索引は表の行として追加してください。`check_plans.py` の **PASS** は確認しましたが、これらの内容矛盾は検出しません。

推奨は、**`in_progress` を維持し、長時間比較の証拠回収と node 軸回帰を最優先で完了させ、判定ツール・文書を修正して再レビューすること**です。`run_0416` の **32.62 ms/step**、RCM の **33.87/33.88 → 32.74/32.74 ms/step** は生ログで確認できており、高速化そのものを撤回する理由はありません。

ファイルは変更していません。本レビューは **plan 未反映**です。

指摘数: Critical 0 / Major 2 / Minor 3
