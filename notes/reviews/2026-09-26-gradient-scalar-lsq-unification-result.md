# codex レビュー: gradient-scalar-lsq-unification (result)

- **plan**: [`plans/active/gradient-scalar-lsq-unification.md`](../../plans/active/gradient-scalar-lsq-unification.md)
- **stage**: `result` (diff base `36d8ba03`)
- **date**: 2026-09-26
- **commit**: `1cfd5623` (feature/sern-design)
- **codex**: effort `high`, 5.3 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M5/m2
- **focus**: Phase 1 (opt-in) の result 1 回目 (plan §5.1 #5r)。S2 は codex diagnose で『条件付き・未合格』と判定済み: 起点床の収束ゲートが設定差のある 4 case で gg・lsq とも不成立 (case/39 は両方 PASS、#5a で chi 単独説明は棄却)、FCT smoke は旧上限超過。これらを保存したまま、(1) 実装 diff (7ebf6d7d・f99f236d) に欠陥が無いか、(2) 物理ゲートと lsq−gg の差・S3 の主張が run の記録で裏付けられているか、(3) Phase 2 (既定化 #6) に進む前に何が不足か、を判定してほしい。結果の一次記録は 4 case README の run_095x_sglsq_* 行と case/09.Taylor-Green/_g0_lsq_seam/ の txt。
- **extra**: `case/16.nozzle_wys/README.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**
LSQ の中核実装に明確な欠陥は見つかりませんでしたが、検証基盤に誤合格・監視漏れを再現しました。
S2 は未合格のままであり、Phase 2 の既定化と `accepted` への移行は承認できません。

`36d8ba03...1cfd5623` の指定 diff を確認しました。係数共有、差分形、境界 incidence の skip、最大4変数のチャンク、周期二重 gather の回避は設計と整合しています。既定 `gg` と cell の GG 固定も維持されています。S0/S1 の保存記録は訂正後の判定を支持します。一方、S2 の対象 run はローカルに存在せず、AWS の残差 CSV・HDF5・VERDICT 原本からの再計算はできていません。以下では、再現した不具合と記録上の未達を区別します。

1. **Major — `check_floor_ratio.py` が「判定不能」や壊れた参照床を PASS にする。**

   [check_floor_ratio.py:39](/home/sano/work/forge-sern-design/solver_density_cuda/tools/check_floor_ratio.py:39) は `CC.analyze()` の総合判定を捨て、入力エラーのうち `(入力)` だけを引き継ぎます。列ごとの「末尾窓の代表値が0」は無視されます。また、[同ファイル:93](/home/sano/work/forge-sern-design/solver_density_cuda/tools/check_floor_ratio.py:93) は参照系列を検証せず床を計算します。

   ファイルを書かないメモリ上の CSV 試験で、必須5列をそろえて次を再現しました。

   | 入力 | 通常判定・参照状態 | 床比判定 |
   |---|---|---|
   | 前半80点が `1e-6`、末尾20点が0 | `check_convergence`: 判定不能 | **PASS** |
   | 参照が全点 Inf、対象が有限の定数 | 参照床 Inf | **PASS** |

   **対案:** 低下桁数・plateau の免除と入力不正を分離し、列単位の判定不能を必ず不合格へ伝播する。参照系列の NaN/Inf・必須列・床の正値性も検査する。この負の試験を追加し、既存の床比判定を再確認する。

2. **Major — S2 の収束ゲート未達を解消する根拠が、追加診断後も不足している。**

   [plan §5.1 #5a–#5c](/home/sano/work/forge-sern-design/plans/active/gradient-scalar-lsq-unification.md:103) の記録では、case/48 の `run_0956_sglsq_s2_gg_chi0` も旧床比 `rms_roUy=9.08` で、chi 単独説明は棄却されています。0476 系は継続後も `--from-floor` 不成立、0482 系も通常収束ゲート未達です。[case/40 の run 索引](/home/sano/work/forge-sern-design/case/40.nozzle_design_tool/README.md:112) では `rms_roOmega` が旧床の約56倍です。

   物理量の差が許容内という記録と、収束ゲート未達という記録は両立します。**両双子が同様に外れることは、共通の回帰不具合を排除しません。** 現行の「条件付き・未合格」という整理は維持すべきです。

   **対案:** 起点床を差し替えず、実装直前版 `36d8ba03` と現行 `gg` を同一入力・実効設定で比較して共通変化を切り分ける。判定区間付きの残差 CSV、床比出力、準定常系列・VERDICT、provenance をレビュー可能な証拠一式として保存する。新基準を採る場合も独立した妥当性確認と事前規則を必要とする。

3. **Major — FCT の保存性と差の妥当性は、#5d の測定だけでは確認できない。**

   [plan:106](/home/sano/work/forge-sern-design/plans/active/gradient-scalar-lsq-unification.md:106) は、流入流出との収支照合が未実施、`nSub=15` の十分性も未確認と明記しています。[case/16 の run 索引](/home/sano/work/forge-sern-design/case/16.nozzle_wys/README.md:393) にある `run_0954/0955_sglsq_s2_fct_*` の差は、旧上限に対して密度約11倍、`roY` 約4.4倍、`roUy` 約230倍です。

   floor 補正が小さいこと、有界性、gg/lsq の総量が近いことは、各 run の保存収支が閉じる証明にはなりません。

   **対案:** `check_passive_budget.py` の FCT 判定で境界流束・ソース・ピン交換を含む収支を確認する。さらに同一物理時刻で sub-iteration 数への感度を調べ、反復不足と作用素変更による差を分離する。旧上限超過は保存し、0476 の定常差を代替上限にしない。

4. **Major — S3 の競合監視が実際の測定バイナリを数えない。**

   [s3_perf.py:22](/home/sano/work/forge-sern-design/case/09.Taylor-Green/_g0_lsq_seam/s3_perf.py:22) が起動するのは `forge_f99f236d` ですが、監視は `pgrep -xc forge`、走行中の汚染判定は件数 `>1` です。このプロセス名を持つ子プロセスを使い、前者が **0件**を返すことを再現しました。測定プロセス自身が数えられないため、別の `forge` が1本加わっても検出条件を満たしません。

   保存結果の比 `1.0097`／`1.0165` と `VERDICT S3: PASS` は計算上整合しています。ただし「走行中監視により競合を排除した」という裏付けが成立しません。根拠は [case/48 の記録](/home/sano/work/forge-sern-design/case/09.Taylor-Green/_g0_lsq_seam/S3_case48.txt:10) と [case/39 の記録](/home/sano/work/forge-sern-design/case/09.Taylor-Green/_g0_lsq_seam/S3_case39.txt:10) です。

   **対案:** GPU 使用 PID と測定対象 PID を照合する監視へ修正し、競合検出の試験を通す。独占実行を保証できる記録がなければ S3 を取り直す。今回の測定が実際に汚染されたとまでは断定しません。

5. **Major — provenance の不足は、既定化前の明示 gg→lsq 切替にも影響する。**

   [stage_manifest.py:225](/home/sano/work/forge-sern-design/solver_density_cuda/tools/stage_manifest.py:225) の hard key は `scalarGradient` を含みません。メモリ上で、他の設定が同一の明示 `gg`／`lsq` 二段を作ると、実際に `segments()` が **1区間**へ連結しました。異なる作用素の残差を一つの収束区間として扱えます。

   [plan §4.4](/home/sano/work/forge-sern-design/plans/active/gradient-scalar-lsq-unification.md:65) の既知課題ですが、指定された `plans/active/tooling-stage-manifest-launch-binding.md` は存在せず、S4 も未完了です。

   **対案:** #2g を起票・完了し、実効 `scalarGradient`、起動順、バイナリ識別を区間判定へ反映する。S4 に同一バイナリでの明示 gg→lsq も加える。修正完了までは両者を同一区間として判定しない。

6. **Minor — 現在仕様の文書に、実装済み opt-in と矛盾する記述が残る。**

   [methods/discretization.md:624](/home/sano/work/forge-sern-design/methods/discretization.md:624) は「LSQ は NS だけ」、[methods/gradient.md:176](/home/sano/work/forge-sern-design/methods/gradient.md:176) は適用条件外のスカラーを一律「片側 GG」と記載しています。実装済みの LSQ 選択時には後者は片側 LSQ です。`procedures/solver-settings.md` にも `scalarGradient` の説明がありません。

   **対案:** 既定 gg／opt-in lsq／cell 固定を明確に分け、「実装中」を「実装済み・検証未完了」へ更新する。#5e の chi 既定値不整合も担当 plan と同期する。`methods/index.md` の既存リンクに問題はありません。

7. **Minor — 追加 run と未解決事項の恒久索引が追いついていない。**

   [case/48 README:55](/home/sano/work/forge-sern-design/case/48.flat_plate_cooled_m4/README.md:55) に chi=0 の `run_0956_sglsq_s2_gg_chi0` がなく、[case/16 README:393](/home/sano/work/forge-sern-design/case/16.nozzle_wys/README.md:393) は0476継続 run の名前だけで、#5b の最終結果を反映していません。#5a は診断完了ですが、床移動の原因は未解決です。

   **対案:** 各 run の結果・成果物・状態を索引へ追記し、床移動の残る切り分けを §5.1 の独立した未完了項目にする。

**推奨は、既定 `gg`・`in_progress` を維持し、検証基盤を修正したうえで不足する証拠をそろえ、Phase 1 の result を再審査することです。** 今回の結果で Phase 2 へ進めず、`accepted` へも移さないでください。

読み取り専用の依頼に従い、ファイルは変更していません。**plan 未反映**です。

指摘数: Critical 0 / Major 5 / Minor 2
