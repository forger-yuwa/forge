# codex レビュー: convection-slau-wall-normal-chi-default (result)

- **plan**: [`plans/active/convection-slau-wall-normal-chi-default.md`](../../plans/active/convection-slau-wall-normal-chi-default.md)
- **stage**: `result` (diff base `22976398`)
- **date**: 2026-09-27
- **commit**: `76eec0d1` (feature/sern-design)
- **codex**: effort `high`, 4.2 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M5/m2
- **focus**: node の slauWallNormalChi 既定 auto (=1) 化の result。検証は §5.1 #6–#9 と §6 (B0・B1 通常域・B1 追加域)、一次記録は case/46・16・39・40 の CHI_DEFAULT_*.txt。diff には同時期の別 plan (スカラー勾配 LSQ) の変更も含まれるので、chi に関係する部分 (solverConfig の三値解決・起動記録・stage_manifest・runner_sern・docs) に絞って判断してほしい。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**  
三値解決とカーネルへの受け渡しは設計どおりですが、検証証拠の不足、軸対称の未実施項目、来歴管理の欠陥が残っています。  
「B0・B1 全件完了、除外域なし」として `accepted` に移す根拠は不足しています。

指定 diff を取得し、chi 関連部分を確認しました。新たな流束演算・符号・単位・絶対閾値の変更はありません。既存の manifest 試験 **15/15**、診断試験 **22/22** は再実行して PASS。ただし、以下の反例もファイル変更なしで再現しました。

1. **Major — B1 の実測を独立に検証できる成果物が不足しています。**

   対象の12分岐 run に `check_convergence.py` を実行した結果、この環境ではすべて `NO residual_history.csv`、総合は `OVERALL: CHECK FAILURES ABOVE` でした。`CONVERGENCE_VERDICT.txt`、設定原本、起動記録、派生量 CSV も対象 run として参照できません。

   [CHI_DEFAULT_B1ii_all.txt:13](/home/sano/work/forge-sern-design/case/46.sern_design/CHI_DEFAULT_B1ii_all.txt:13) にはプラトー許容 PASS と通常判定 `NOT CONVERGED` が併記されていますが、全列の末尾平均比や比較元、実行引数がありません。「省略側 ≤ 2× 明示0側」を再計算できません。接続模型についても [README.md:347](/home/sano/work/forge-sern-design/case/46.sern_design/README.md:347) は場ファイル削除済みとしています。

   **対案:** 設定・起動記録・全残差 CSV・派生量 CSV・省略していない判定出力・入力とバイナリのハッシュを、run ごとに保存してください。遠隔に残っているなら所在と取得手順を明記し、レビュー可能にする必要があります。記載数値が誤りとは断定しませんが、現状は裏付け未確認です。

2. **Major — B1-a は登録済みの等温壁熱流束検証を満たしていません。**

   [plan:152](/home/sano/work/forge-sern-design/plans/active/convection-slau-wall-normal-chi-default.md:152) は等温壁 `q_w(x)`、相対 L2 ≤ 1%、3点の準定常確認を要求しています。一方、[CHI_DEFAULT_B1a.txt:21](/home/sano/work/forge-sern-design/case/40.nozzle_design_tool/CHI_DEFAULT_B1a.txt:21) は「断熱壁のため対象外」として PASS にしています。

   `run_0954/0955_chidef_b1a_*` の壁圧差 **0.3765%** や推力・流量の結果は、この欠測を補いません。抽出器で事前に対象外と書いていても、plan の受入条件は変更されていません。

   **対案:** 等温壁の軸対称 A/B を追加し、登録した熱流束基準まで検証してください。それまでは #8 を未完了に戻し、「軸対称を含め全域 PASS」を撤回すべきです。

3. **Major — manifest は同じ設定の異なる起動を区別できません。既存試験もその反例を検査していません。**

   [stage_manifest.py:152](/home/sano/work/forge-sern-design/solver_density_cuda/tools/stage_manifest.py:152) は `cfg_fnv` ごとの最後の値だけを保持し、[同:288](/home/sano/work/forge-sern-design/solver_density_cuda/tools/stage_manifest.py:288) で全対応段を上書きします。

   実際に同じ YAML の起動記録を chi `[0, 1]` とすると、結果は **実効値 `[1, 1]`、1区間**でした。[試験:40](/home/sano/work/forge-sern-design/solver_density_cuda/tools/test_stage_manifest_wall_normal_chi.py:40) は「同一 YAML」と書きながら別ハッシュ `"aa"` / `"bb"` を与えているため、この欠陥を逃しています。

   **対案:** 起動結合の本実装が[別 plan で延期済み](/home/sano/work/forge-sern-design/plans/active/tooling-stage-manifest-launch-binding.md:6)であることを、本 plan §4.3・§5.1 にも反映してください。当面は chi の既定変更をまたぐ継続を禁止し、対応が曖昧な起動記録では区間判定を停止する必要があります。「来歴対応完了」のままにはできません。

4. **Major — 設計 DB は chi の異なる評価・実効値不明の評価を混ぜられます。**

   [driver_sern.py:70](/home/sano/work/forge-sern-design/design/forge_design/opt/driver_sern.py:70) は日付と `scalarGradient` だけを確認します。`status=PASS`、現行日付、`scalarGradient=lsq` の行で、chi を **`None` / `0` / `1`** としたところ、すべて `_learnable=True` でした。明示0を許す現行 runner では、日付だけで chi 方針を区別できません。

   さらに [runner_sern.py:61](/home/sano/work/forge-sern-design/design/forge_design/evaluate/runner_sern.py:61) は最後の起動が壊れている・chi を欠く場合に過去の値へ遡ります。反例では過去の `0` を返し、不正値 `2` も受理しました。

   **対案:** chi の取得を既存の `_last_launch_value(..., allowed=("0","1"))` と同じ厳格な処理に統一してください。学習・Pareto 選別では、全作動点の chi が campaign の期待値と一致することを要求し、不明値を除外してください。

5. **Major — 推奨設定の正本に相反する運用指示が残っています。**

   [recommended-settings.md:70](/home/sano/work/forge-sern-design/procedures/recommended-settings.md:70) は auto=1 を推奨しますが、同じ節の [87行](/home/sano/work/forge-sern-design/procedures/recommended-settings.md:87) は診断不能な新構成で1を禁止し、[92行](/home/sano/work/forge-sern-design/procedures/recommended-settings.md:92) は **「現行運用は0」「1にしてよい根拠ではない」** としています。

   **対案:** 旧方針を履歴として明確に分離し、現行本文を auto=1 と「0へ変更する判断材料」に統一してください。衝撃足の未評価という制限は残すべきです。

6. **Minor — 3D runner の来歴出力が実装されていません。**

   [runner_sern3d.py:329](/home/sano/work/forge-sern-design/design/forge_design/evaluate/runner_sern3d.py:329) は chi 実効値・`flag_policy` を追加せず `metrics.json` を書きます。[plan #4](/home/sano/work/forge-sern-design/plans/active/convection-slau-wall-normal-chi-default.md:112) の `runner_sern*.py` 対応完了とは一致しません。また B1(i) は `cad/run_junction_model.py` 経由であり、3D runner 自体の実行検証を代替しません。

   **対案:** 2D・3D 共通の来歴取得処理を使い、3D runner の変換から計算・収集までの確認を #4 に残してください。

7. **Minor — 残作業表・結果節・run 索引が完了状態と一致していません。**

   [plan:115](/home/sano/work/forge-sern-design/plans/active/convection-slau-wall-normal-chi-default.md:115) は完了記述の後に「走行中」「未」が残り、[§6.2](/home/sano/work/forge-sern-design/plans/active/convection-slau-wall-normal-chi-default.md:177) は「未実施」です。[case/46 README:348](/home/sano/work/forge-sern-design/case/46.sern_design/README.md:348) も基準場が走行中のまま。case/16・39・40 には対象 chi run の索引行がありません。

   **対案:** 上記未解決事項を §5.1 に登録し、§6.2 と各 case の run 一覧を同期してください。`plans/README.md` の「draft・planレビュー前」も更新が必要です。

なお、`NOT CONVERGED` や `OSCILLATING` の存在だけを不合格理由にはしていません。plan はプラトーと差区間による受入を認めており、記載された m4_off の振動列も差区間は帯内です。問題は、その受入条件を再検証できる証拠と、未実施項目の扱いです。

**推奨は、`in_progress` を維持し、証拠の保存・再判定、B1-a の熱流束検証、来歴と DB 選別の修正、文書同期を終えてから result レビューを再実施することです。** ファイルは変更しておらず、本指摘は plan 未反映です。

指摘数: Critical 0 / Major 5 / Minor 2
