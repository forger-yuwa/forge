# codex レビュー: convection-slau-wall-normal-chi (result)

- **plan**: [`plans/accepted/convection-slau-wall-normal-chi.md`](../../plans/accepted/convection-slau-wall-normal-chi.md)
- **stage**: `result` (diff base `39526328`)
- **date**: 2026-09-25
- **commit**: `1795e3d0` (feature/sern-design)
- **codex**: effort `high`, 5.3 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M0/m3
- **focus**: 4 回目の result 段。(1) 対応表の『限定』行 (R2-M2/R2-M3/R3-M2/R3-M4 = 固定点・格子感度) は受入条件から外して §5.1 #11 (既定化) の前提へ移した (§2 に線引きを明記)。この移動が opt-in・既定 0 の受入として妥当かを判定してほしい (再測定を求めるなら、どの量をどの長さで)。(2) §6.2『受入の範囲』の各文が対応表の『閉』行のみを根拠にしているか。(3) #14 の判定列を差系列・自 run 偏差から run の量 (絶対 x_foot・窓内壁圧) に変えた理由 (operand の型) が結果依存でないか。(4) V1-c (全接続面収支の凍結近傍求根) と §4.2 (単面近似、合否に使わない) の区別が読めるか。(5) V3 表『収束の悪化』行: SERN 2D は登録量 (本段の低下桁数、flag1 ≥ flag0−0.2 桁) で PASS、分岐 2 件 (case/16・case/48) は代替量という扱いが受入根拠として成立するか。前 3 回の型は『測れた範囲を超えて書いた』。
- **extra**: `notes/reviews/2026-09-25-convection-slau-wall-normal-chi-result.md`, `case/46.sern_design/_r3_m1m3/V1B_b_VERDICT.txt`, `case/46.sern_design/_r3_m1m3/run0437_budget.txt`, `case/46.sern_design/_r3_m1m3/PERRUN_QS_flag0.txt`, `case/46.sern_design/_r3_m1m3/PERRUN_QS_flag1.txt`, `case/46.sern_design/_r3_m1m3/PROVENANCE.txt`, `case/46.sern_design/_r3_m1m3/sern2d_conv/CONV_DETERIORATION_VERDICT.txt`, `case/46.sern_design/_r3_m1m3/sern2d_conv/CONV_flag0.txt`, `case/46.sern_design/_r3_m1m3/sern2d_conv/CONV_flag1.txt`, `solver_density_cuda/tools/test_stage_manifest_wall_normal_chi.py`, `case/16.nozzle_wys/v3_case16/B_10b_RULE.txt`, `case/16.nozzle_wys/v3_case16/C_CONVERGENCE_RATIO.txt`, `case/16.nozzle_wys/v3_case16/QUASISTEADY_VERDICT.txt`, `case/16.nozzle_wys/_v2/V2_VERDICT.txt`, `case/16.nozzle_wys/_v2/SIDE_new0_vs_new1.txt`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
既定0の opt-in として、検証したケース・設定・有限区間に限定した受入は支持します。  
固定点・格子収束を #11 へ移す判断も妥当ですが、「害が無い」という一般化と、以下の記録不整合は移動前に直してください。

HEAD `1795e3d0` で指定 diff を取得しました。質量流束だけの変更、再構成後の面速度、壁 mask、境界半割面の除外、config の拒否条件は設計と整合しています。追加式に明白な符号・単位・ゼロ割の誤りは見つかりませんでした。

再実行では、流束単体試験 `ALL PASS`、V5 保存ダンプの P0–P5 `PASS`、`stage_manifest` 7/7 `PASS`、case/16 V2 `VERDICT: PASS` を確認しました。SERN の #14 は両側 `ALL STEADY`。SERN 本段と case/16 分岐の残差判定は、記載どおり両側 **`NOT CONVERGED`** です。

1. **Minor — #14 の分布全体の根拠を、既存の `selfL2_pct` から補ってください。**

   **根拠:** [抽出コード:55](/home/sano/work/forge-sern-design/case/46.sern_design/v3sern_foot_perrun.py:55) は窓内34点の分布偏差を計算していますが、判定対象は窓平均と3点です。この4量の `STEADY` だけでは、点間で相殺する分布変動を排除できません。

   ただし、保存 CSV の末尾29 dump を再計算すると、末尾平均分布からの最大相対 L2 は flag0 **0.031285 %**、flag1 **0.017736 %**。三角不等式による任意の2 dump 間の上限は **0.062570 %／0.035472 %**です。同じ正規化で線形回帰ドリフトの上限も **0.090619 %／0.051373 %**となり、0.2 %以内です。

   **対案:** この分布全体の上限を #14 と対応表へ追記してください。偏差系列を自身の平均で割る判定から、物理量の系列へ変更した理由は妥当です。**追加の solver run は不要**で、既存 CSV から根拠を補えます。

2. **Minor — case/48 の「収束悪化」根拠が、末尾水準と終端瞬時値を混同しています。**

   **根拠:** [受入結論:1225](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi.md:1225) の `rms_ro` **1.83e−7 対1.94e−7** は終端値です。`case/48.flat_plate_cooled_m4/_v3/v3_flag{0,1}` の末尾40 %平均では、flag1/flag0 は `rms_ro` **1.0108**、非ゼロ全列で **0.9760–1.0114**でした。悪化が小さいという結論は維持できますが、瞬時値から「低い」と説明すべきではありません。

   SERN はツールと同じ末尾中央値で再計算して、低下桁数差の最悪値が **−0.04370桁**。登録条件の **−0.2桁以上**を全8列で満たしました。case/16 の末尾平均比 **1.0003–1.0024**も再現しています。

   **対案:** case/48 も全列の末尾平均比で記録し、事後評価と明記する。case/16 の種は `NOT CONVERGED` なので、「収束場からの分岐」は「残差がプラトーの場からの分岐」に直してください。これらは有限区間の回帰根拠として成立しますが、漸近的な収束率の保証にはなりません。

3. **Minor — 撤回済みの結論が本文・台帳に残り、対応表の「閉」と矛盾しています。**

   **根拠:** 以下が残っています。

   - [plan:1107](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi.md:1107)：固定点「済・PASS、0.01 %一致」
   - [plan:763](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi.md:763)：「5桁目まで一致」
   - [plan:1257](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi.md:1257)：固定点を V3 で確認する旧方針
   - [case/16 README:377](/home/sano/work/forge-sern-design/case/16.nozzle_wys/README.md:377)：「反復ノイズ以下」
   - [case/46 README:342](/home/sano/work/forge-sern-design/case/46.sern_design/README.md:342)：`run_0439` が「実行中」

   **対案:** 現役の要約・台帳を最新結論へ同期し、過去の記述には撤回と後継節を付けてください。§6.2 は「受入根拠」と「未解決・適用限界」を分ければ、「限定」行の観測を合格根拠に混ぜずに残せます。

重点への判断として、**固定点・格子収束の再測定を今回の受入条件には要求しません**。ただし受入文は「試験した条件・区間・判定量で回帰許容内」としてください。既定0であることだけでは、flag1 使用時の無害性は証明できません。

V1-c の全接続面・凍結近傍求根と、§4.2 の単面近似の区別は現在の記述で読めます。V1-b の保存記録も、起点の正味流入と最大 **0.0603 %/dump、`VERDICT: PASS`**を示しています。ただし、AWS 側の元 HDF5・V6/V7・SERN の力係数時系列は今回独立再計算できておらず、保存記録と plan の記載による確認です。

**推奨は、上記1→2→3を反映し、既定0と #11 の未解決事項を維持したまま `accepted` へ移すことです。** ファイルは変更しておらず、レビュー提案は plan 未反映です。

指摘数: Critical 0 / Major 0 / Minor 3
