# codex レビュー: species-passive-scalar-unification (plan)

- **plan**: [`plans/active/species-passive-scalar-unification.md`](../../plans/active/species-passive-scalar-unification.md)
- **stage**: `plan`
- **date**: 2026-09-16
- **commit**: `22b7f0ba` (feature/sern-design)
- **codex**: effort `high`, 6.4 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M6/m1
- **focus**: 12 回目 (plan-11 NO-GO M1–M6 の反映; ゲートの完備性のみ, 実装済 commit 済)。変更: passive_gate_common.py に集約 — (1) solver と同じ既定値 (bdfOrder 2, nSub 20 等) で正規化した実効設定、整数の型/範囲検査; config_diff で許可差分 (dt/nStepOuter/nSub/outStepInterval/monitorInterval/valueFileName) 以外の solverConfig 全項目と bcondConfig・メッシュ・IC の md5 を照合; (2) residual_history_check: config 由来の必須列 (流れ・SST・化学種・受動種) の存在、全行の全数値の有限性、期待 step 0..nStepOuter−1 の outer_begin/outer_end/inner_iter、inner 番号 1..nSub−1 (solver は最終 sub-iter を outer_end に書く)、ゼロ始点は全 inner 行ゼロのときだけ受理、全列の低下 ≥2 桁; (3) check_passive_budget: 全 [passive] 行の全数値トークンを読込時点で検査し一度でも非有限/解析不能なら失敗を保持、CSV と res_<nStepOuter>.h5 を必須、mode auto = fct / conservative (全 bcond periodic かつ凝縮なしのみ) / unsteady (保存は判定不能と明示) / steady; (4) 終了場は res_<nStepOuter>.h5 に固定; (5) 実効刻み dt_eff = float32(dt) と名目時刻 nStepOuter×dt_eff を正本に checkpoint (double 積算) と照合; (6) end-to-end 試験 tests/unit/test_passive_gates.py (合成 run で main()→exit code: 正常系 [2 次系列 PASS, dt 非依存 exact PASS], 反例 [1 次系列, 最終場欠落, 古い場のみ, CSV 欠落, nan トークン, 列欠落, 0→1→0, inner 不足, convMethod/bcond 差, KEEP で --expect-fct, 特異モーメント, ξ>1])。問い: 残る偽 PASS/誤 FAIL、§6-2/§6-6 の合否をこれらの VERDICT に委ねてよいか。設計の再指摘は不要。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**  
対象は `22b7f0ba`。ゲートに偽 PASS と誤 FAIL が残り、§6-2／§6-6 の最終合否はまだ委ねられません。  
数値設計の再指摘はせず、コード・既存 run・ファイルを作らない模擬入力で確認した問題に限定します。

1. **Major — sub-iter の初回残差を取り違えている**

   **根拠:** solver は反復番号 **0** の BDF 込み残差を `outer_begin`／`inner_begin` に記録します（[main.cpp:1754](/home/sano/work/forge/solver_density_cuda/main.cpp:1754)、[同:290](/home/sano/work/forge/solver_density_cuda/main.cpp:290)）。しかし、[passive_gate_common.py:212](/home/sano/work/forge/solver_density_cuda/tools/passive_gate_common.py:212) は `inner_iter` だけを保持し、番号 **1** を低下判定の始点にしています。

   実書式に合わせた反例で次を確認しました。

   | 反復番号 0→1→2→3 の残差 | 現行判定 |
   |---|---|
   | `0 → 1 → 0.01 → 0.001` | **3 桁低下として PASS**。実際のゼロ始点を見逃す |
   | `1 → 0.001 → 0.0003 → 0.0001` | **1 桁として FAIL**。実際には 4 桁低下 |

   前者は次数ツールの `main()` でも **`VERDICT: PASS`、exit 0** でした。なお、現行 solver の `outer_end` は最後の `inner_iter` と同じ残差を参照します（[main.cpp:330](/home/sano/work/forge/solver_density_cuda/main.cpp:330)）。

   **対案:** `inner_begin(0)` から `inner_iter(nSub−1)` までを評価してください。`outer_end` との一致、行の重複・欠落・範囲外番号も検査し、ゼロ始点判定に番号 0 を含める必要があります。

2. **Major — 「解析不能・過去の異常を保持する」契約が未完成**

   **根拠:** [check_passive_budget.py:54](/home/sano/work/forge/solver_density_cuda/tools/check_passive_budget.py:54) の正規表現は、数値として認識できる部分だけを抽出します。`BROKEN` のような数値欄の異常文字列は抽出されません。さらに、解析結果は最新行で上書きされ、`per-step` 等を含む行は無条件に読み飛ばせます（[同:59](/home/sano/work/forge/solver_density_cuda/tools/check_passive_budget.py:59)、[同:81](/home/sano/work/forge/solver_density_cuda/tools/check_passive_budget.py:81)）。

   以下を再現しました。

   - 最終行の `per-step ... lo BROKEN`：**`main()` が PASS、exit 0**。
   - 過去の `nonfinite 1` の後に正常な最終記録：**`main()` が PASS、exit 0**。
   - 過去の累積絶対補正が `BROKEN`、最終行は正常：`parse_lines → evaluate` が **True**。

   **対案:** 各書式の数値欄を全文解析し、その場で変換・有限性・非有限フラグを検査してください。解析失敗と異常フラグは最新値への上書きと独立して保持し、部分文字列による無条件受理を廃止してください。

3. **Major — メッシュ・IC・BC が欠落しても同一扱いになり、外部入力の差も見逃す**

   **根拠:** [_md5:65](/home/sano/work/forge/solver_density_cuda/tools/passive_gate_common.py:65) はファイル欠落時に `None` を返します。[config_diff:99](/home/sano/work/forge/solver_density_cuda/tools/passive_gate_common.py:99) は戻り値の等値比較だけなので、両側欠落は一致になります。

   また、内容比較は `bcondConfig.yaml`・メッシュ・IC に限定されています。solver が別途読む `inlet_profile_<physID>.csv`（[boundaryCond.cpp:249](/home/sano/work/forge/solver_density_cuda/boundaryCond.cpp:249)）や `speciesDBFile` の内容は対象外です。

   模擬系列で、次の両方が次数ツールの **`VERDICT: PASS`、exit 0** になりました。

   - 全 run からメッシュ・IC・`bcondConfig.yaml` を欠落させる。
   - 設定上の参照名を同じにして、run ごとに入口プロファイル・物性 DB の内容を変える。

   **対案:** 必須入力の存在を比較前に検査してください。そのうえで、有効な BC が参照する入口プロファイルと外部物性 DB も内容照合に含めてください。

4. **Major — 化学種の時間 1 次誤差を、既定の次数ゲートが検出しない**

   **根拠:** [analyze_moment_order.py:43](/home/sano/work/forge/case/44.vitiated_air_wt/analyze_moment_order.py:43) の既定評価量には化学種がありません。検査は指定された `fields` だけで、config 由来の必須評価量との照合もありません（[同:86](/home/sano/work/forge/case/44.vitiated_air_wt/analyze_moment_order.py:86)）。残差列の存在確認は、化学種の時間次数確認を代替しません。

   非一様な `roY0/roY1` に刻み比例の誤差を与え、他の評価量を刻み非依存とした系列を通しました。

   - **既定評価量:** `VERDICT: PASS`、exit 0。
   - **`roY0/roY1` を追加:** 両種の観測次数 **1.000**、`VERDICT: FAIL`、exit 1。

   **対案:** §6-6 用の正式ゲートでは、config から全化学種・トレーサ・全凝縮種の必要評価量を生成してください。`--fields` による部分評価を正式な全体 PASS と区別し、`exact` は時間次数を実証した量として数えないでください。

5. **Major — 保存が「判定不能」でも、収支ゲート全体が PASS を返す**

   **根拠:** [check_passive_budget.py:147](/home/sano/work/forge/solver_density_cuda/tools/check_passive_budget.py:147) は `unsteady` モードで保存判定を省略しますが、最終結果は通常の **PASS／exit 0** です（[同:228](/home/sano/work/forge/solver_density_cuda/tools/check_passive_budget.py:228)）。

   非 FCT・非周期の模擬 run で、初期総量 **1**、最終総量 **2**、境界流束・ソースの収支記録なしでも、`main()` は **`VERDICT: PASS`、exit 0** でした。総量増加自体が誤りとは限りませんが、保存を検証した証拠にはなりません。

   **対案:** 保存判定不能は `INDETERMINATE` 等の独立した結果と非成功終了コードにしてください。floor／lim／場の部分合格は別表示とし、§6-2 の合格には収支検証済みを必須にしてください。

6. **Major — checkpoint restart の開始時刻を無視し、正当な継続計算を拒否する**

   **根拠:** [passive_gate_common.py:60](/home/sano/work/forge/solver_density_cuda/tools/passive_gate_common.py:60) は終了時刻を常に `nStepOuter × dt_eff` とします。しかし solver は履歴復元時に checkpoint の時刻を引き継ぎます（[main.cpp:1035](/home/sano/work/forge/solver_density_cuda/main.cpp:1035)）。

   実 run `case/09.Taylor-Green/run_0095_passiveC_m3_restart100/` では、ログ上の開始時刻は **0.2000000477**、追加計算は **100 step × 約 0.002**。終了 checkpoint は **0.400000453** ですが、現行場ゲートは **0.2000000095** を期待して FAIL でした。旧成果物の丸め差とは別に、開始時刻約 **0.2** を丸ごと落としています。

   **対案:** 実際に復元された開始時刻を正本とし、`t_start + nStepOuter × dt_eff` を照合してください。履歴復元成功と BDF1 再開を区別し、現在の double 積算による正常 restart を試験に追加してください。

7. **Minor — 実効設定を作っても、比較には未正規化 YAML を使っている**

   **根拠:** [passive_gate_common.py:92](/home/sano/work/forge/solver_density_cuda/tools/passive_gate_common.py:92) は `ca['raw']`／`cb['raw']` を比較します。`bdfOrder: 2` の明示と省略は同じ実効設定ですが、実行すると **`solverConfig time.bdfOrder: 2 vs None`** が差分として返ります。

   **対案:** solver と同じ既定値を反映した設定を比較してください。元 YAML は差分説明用に保持する構成が適切です。

**推奨は、数値設計を維持し、上記の順でゲート修正を完了してから正式検証を判定することです。** 特に、部分検査の PASS と §6 全体の合格を結び付ける判定を明示してください。

既存の `test_check_passive_budget.py` は **ALL PASS** でした。今回の反例はファイル入出力をメモリに置き換え、判定関数と `main()` を実行して確認しています。

実成果物も再検査しました。`case/44.vitiated_air_wt/run_0321_passiveG_order3_bdf2_dt1.6e-5_nsub40/`、`run_0322_passiveG_order3_bdf2_dt8e-6_nsub40/`、`run_0323_passiveG_order3_bdf2_dt4e-6_nsub40/`、`run_0324_passiveG_order3_bdf2_dt8e-6_nsub80/` の終了場を比較した現行次数ゲートは **`VERDICT: FAIL`**。`Q0/Q1/Q2` の次数は **1.374／1.626／1.652** でした。これはゲート修正版 solver の合格実績ではありません。run 索引は [case README](/home/sano/work/forge/case/44.vitiated_air_wt/README.md) です。

ファイル変更・新規 CFD 実行はしていません。提案は **plan 未反映**です。

指摘数: Critical 0 / Major 6 / Minor 1
