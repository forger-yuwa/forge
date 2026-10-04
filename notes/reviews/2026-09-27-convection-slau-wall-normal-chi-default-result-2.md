# codex レビュー: convection-slau-wall-normal-chi-default (result)

- **plan**: [`plans/active/convection-slau-wall-normal-chi-default.md`](../../plans/active/convection-slau-wall-normal-chi-default.md)
- **stage**: `result` (diff base `22976398`)
- **date**: 2026-09-27
- **commit**: `8556373a` (feature/sern-design)
- **codex**: effort `high`, 4.4 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M2/m1
- **focus**: result 2 回目。1 回目 NO-GO の 7 件は §5.1 #11 で対応、B1-a 等温壁の壁圧は出口角 1 点をユーザ決定で除外 (#12)。一次記録は case/46.sern_design/chi_evidence/ (JUDGEMENTS.txt・INDEX.md) と各 case の CHI_DEFAULT_*.txt。chi 関係の変更 (三値解決・起動記録・stage_manifest の ambiguous・runner_sern/runner_sern3d/driver_sern の来歴と選別・docs) に絞って、accepted にしてよいか判断してほしい。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**  
三値解決・来歴の選別・残差の受入判定は概ね裏付けられました。  
ただし、前回修正で SERN の形状生成が壊れており、出口角除外後の合格値も保存証拠から再計算できません。

1. **Major — 前回修正で `_dv` が削除され、2D・3D の設計チェーンが停止します。**

   commit `34aa2335` は `_last_launch_chi` の修正と同時に `_dv` の定義を削除しています。しかし [runner_sern.py:210](/home/sano/work/forge-sern-design/design/forge_design/evaluate/runner_sern.py:210) には呼び出しが残っています。3D も [runner_sern3d.py:114](/home/sano/work/forge-sern-design/design/forge_design/evaluate/runner_sern3d.py:114) から同じ関数を使用します。

   実際に `problem_moo_frozen_tp_cycle3op.yaml` を `load_problem` で読み、`design_from_problem(p)` を呼ぶと、**`NameError: name '_dv' is not defined`** を再現しました。保存済み B1 の成功は、修正後 HEAD の正常動作を保証していません。

   **対案:** `_dv` を復元し、実際の問題 YAML を使う形状生成の回帰試験を追加してください。2D・3D の準備経路と来歴収集まで確認し、plan §5.1 #4 の完了根拠を更新する必要があります。現在の選別試験23件は通りますが、この経路を検査していません。

2. **Major — 出口角除外後の `L∞ = 0.031%` を再現できる証拠・判定処理がありません。**

   出口角を除外するユーザ決定自体は異議の対象にしません。問題は、その決定を適用した実測判定が保存されていないことです。

   [chidef_b1a.py:120](/home/sano/work/forge-sern-design/case/40.nozzle_design_tool/chidef_b1a.py:120) は依然として全壁ノードの最大差を取り、除外処理がありません。[CHI_DEFAULT_B1a_isoT_ext.txt:17](/home/sano/work/forge-sern-design/case/40.nozzle_design_tool/CHI_DEFAULT_B1a_isoT_ext.txt:17) の記録も **`1.3583% → FAIL`**、総合は **`VERDICT B1-a 量: FAIL`** のままです。

   一次記録の [chidef_b1a_wall.csv:1](/home/sano/work/forge-sern-design/case/46.sern_design/chi_evidence/40/0959_chidef_b1a_isoT_omit_ext/chidef_b1a_wall.csv:1) は壁圧3点の時系列だけです。全221点の座標・圧力・除外マスクがなく、ローカルに対象 HDF5 もないため、「超過は出口角1点だけ」「残る220点は最大0.031%」を独立に検算できません。

   **対案:** `run_0958/0959_chidef_b1a_isoT_*_ext` の全壁点について、ノードID・座標・両 run の圧力・境界所属を保存してください。壁∩出口で除外点を特定する処理を追加し、**除外前後の最大差・除外点数・新しい VERDICT** を併記すべきです。元の FAIL 記録は履歴として残してください。

3. **Minor — 前回指摘の文書同期が一部未完了です。**

   [plan:117](/home/sano/work/forge-sern-design/plans/active/convection-slau-wall-normal-chi-default.md:117) は完了記述の末尾に「m4_off・m10_on は走行中」「(iii) 未」を残しています。[plan:179](/home/sano/work/forge-sern-design/plans/active/convection-slau-wall-normal-chi-default.md:179) も、決着した #12 を「判断待ち」としています。

   **対案:** 残作業表とレビュー対応欄を同期し、上記1・2を未完了項目として登録してください。各 case の run 索引と `plans/README.md` の更新は確認できました。

確認できた範囲では、chi の流束式・符号・単位・ゼロ割閾値への変更はありません。manifest 試験 **17件**、設計 DB 選別試験 **23件**は PASS。同一 `cfg_fnv` の chi `[0,1]` は `ambiguous` となり、`check_convergence --segment` が判定不能・終了コード1になることも再現しました。試験用ファイルはメモリに置き換えています。

保存残差21本を `check_convergence` の判定関数で再評価すると、通常判定は **NOT CONVERGED、非有限値・RISINGなし**。全8組の床比判定は **PASS**、最大末尾平均比は **1.141214倍**でした。保存 CSV の準定常判定も再現でき、SERN の `m4_off` は一部 **OSCILLATING**、`m10_on` は **STEADY**、5列の差区間は双方とも帯内です。これは plan のプラトー許容による受入を支持しますが、通常基準での収束を意味しません。

**推奨は `in_progress` を維持し、1 → 2 → 3 の順に対応してから再判定することです。** ファイルは変更しておらず、本指摘は plan 未反映です。

指摘数: Critical 0 / Major 2 / Minor 1
