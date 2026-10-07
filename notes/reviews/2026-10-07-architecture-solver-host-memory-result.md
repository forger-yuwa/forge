# codex レビュー: architecture-solver-host-memory (result)

- **plan**: [`plans/active/architecture-solver-host-memory.md`](../../plans/active/architecture-solver-host-memory.md)
- **stage**: `result` (diff base `main`)
- **date**: 2026-10-07
- **commit**: `f5cff054` (feature/sern-design)
- **codex**: effort `high`, 4.6 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M2/m2

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**GO-with-changes**

R1–R3 は設計方針に沿っており、確認した差分に数値カーネルの変更や新たな符号・単位・境界処理の誤りはありません。  
SERN g3 のホストメモリ約49%削減は保存報告と整合します。  
ただし、比較器の欠落検出と検証証跡・完了記述を補うまでは `accepted` への移管を勧めません。

`main...HEAD` の指定差分を取得し、本変更は `9c9f623c..93e55957` として切り分けました。`hostCell` の検査順序、ゴーストを含む確保長、出力名の共通化、面配列の退避 fallback は設計と一致しています。

検証については、保存 TSV の **3147＋255量、FAIL 0**、固定幅の SHA256、追加診断の単体試験 **8件 PASS** を確認しました。ただし AWS の `res_*.h5`・残差 CSV・個別 VERDICT はこの作業ツリーに無く、場の再集計はできていません。[run 索引]( /home/sano/work/forge-sern-design/case/66.hostmem_regression/README.md:165)は全 forge run を `NOT CONVERGED` と明記しており、定常解の一致を偽っている記述は確認しませんでした。

1. **Major — 比較対象が欠落していても、追加診断 B が失敗になりません。**

   **根拠:** [`compare_runs.py:562`](/home/sano/work/forge-sern-design/case/66.hostmem_regression/compare_runs.py:562) は「存在する最大ステップ」を最終出力として採用し、予定した N step への到達を検査しません。さらに、存在しないディレクトリを base/new 各3本として `compare_abs` に渡す実行確認で、**「量 0、FAIL 0」**になりました。[集計処理:766](/home/sano/work/forge-sern-design/case/66.hostmem_regression/compare_runs.py:766)もこれを失敗にしません。`runs_for` は終了コードによらず `done` を選ぶため、全反復が同じ地点で失敗・早期終了する場合も検出できません。

   **対案:** 比較前に、構成ごとの予定反復数、実行成功、必要ファイル、最終ステップ・物理時刻、必須列を検査し、欠落・比較量ゼロを明示的な FAIL にしてください。欠落、初期出力のみ、予定 N 未到達の負例を追加し、既存成果物に完全性検査を掛け直すことが移管条件です。今回の3402量が空だったという指摘ではありません。

2. **Major — メモリの「傾き」と負例試験を、提出された証跡から再検証できません。**

   **根拠:** [plan §5.1 #4](/home/sano/work/forge-sern-design/plans/active/architecture-solver-host-memory.md:130)には、傾き **2737→2227→1371 B/節点**と `P`・`roN` の負例成功が書かれていますが、対応する run パス・工程別ログ・停止ログへの参照がありません。掲載された `memlog_table.txt` は変更前の測定です。

   一方、[g3 の工程別報告](/home/sano/work/forge-sern-design/case/66.hostmem_regression/results/2026-10-07_base9c9f623c_vs_new93e55957/sern_g3_memlog.txt:18)は **5157→2630 MiB**を支持しますが、同一格子の反復なので傾きの検証にはなりません。[README:242](/home/sano/work/forge-sern-design/case/66.hostmem_regression/README.md:242)の「GPU 単点で同値 → GPU 傾きの条件を満たす」も論拠が不足しています。

   **対案:** s070/s085 の各段階について、入力・バイナリ、節点数、工程別 RSS/HWM/GPU、回帰式を追跡できる測定表を残してください。負例も入力、変更内容、終了コード、名前付き停止ログを結び付けてください。実施済み試験の記録補完を優先し、原本が失われている項目だけ再検証するのが妥当です。

3. **Minor — checkpoint の「決定的な一致」と、時間発展後の比較が混同されています。**

   **根拠:** [§6.4:244](/home/sano/work/forge-sern-design/plans/active/architecture-solver-host-memory.md:244)は checkpoint を決定的な比較にまとめ、分割／連続も「全判定を通った」範囲に含めています。しかし `c09ckpt100` の最終 `CHECKPOINT/roXiP` は **D_abs＝2.980232e−6**であり、ビット一致ではありません。[比較表:248](/home/sano/work/forge-sern-design/case/66.hostmem_regression/results/2026-10-07_base9c9f623c_vs_new93e55957_abs/all_quantities.tsv:248)

   また、[`split_vs_cont.py:8`](/home/sano/work/forge-sern-design/case/66.hostmem_regression/split_vs_cont.py:8)は明示的に「判定はしない」としています。

   **対案:** 「共通 checkpoint の復元直後はビット一致」「時間発展後は A/B の定義した比較」「分割／連続は差を記録した参考比較」と分けてください。これは checkpoint 実装の破損を示すものではなく、確認済み範囲の記述修正です。

4. **Minor — 現在仕様と残作業表が、移管できる状態まで同期されていません。**

   **根拠:** [`methods/architecture/overview.md:192`](/home/sano/work/forge-sern-design/methods/architecture/overview.md:192)は現在形でホスト全確保を説明し、その後の括弧内で新実装を説明しています。旧メモリ値と「常駐＋約150 MB」も、旧測定であることが不明瞭です。

   [plan §5.1 #3–#4](/home/sano/work/forge-sern-design/plans/active/architecture-solver-host-memory.md:129)には「途中」「未実行」「本格回帰が残り」が残っています。§8 が要求する実測後の規模見積もりも、[起点 plan の B4b-2](/home/sano/work/forge-sern-design/plans/active/tooling-sern-mesh-blocking.md:848)へ未反映です。

   **対案:** 現在仕様を主文にし、旧値は変更前測定として分離してください。残作業表を「完了／範囲限定で免除／後継 plan へ移管」に整理し、NaN ダンプ保留も明示的に処理してください。B4b-2 には傾き・切片と共有環境の実効予算を用いた見積もりを返し、1000万節点は実行保証と区別してください。移管時には `plans/README.md` と参照リンクも同期が必要です。

**推奨は、§6.4 の範囲限定を維持し、上記1→2→3→4の順で修正してから移管することです。** 凝縮 dual-time の固定幅試験は、記録どおり base/new とも4/6本が超過した「判定不能」のまま残してください。今回の受入れを数値同等性の立証に言い換えるべきではありません。

ファイル変更なし。指摘・推奨は plan 未反映です。

指摘数: Critical 0 / Major 2 / Minor 2
