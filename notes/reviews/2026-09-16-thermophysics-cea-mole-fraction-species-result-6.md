# codex レビュー: thermophysics-cea-mole-fraction-species (result)

- **plan**: [`plans/active/thermophysics-cea-mole-fraction-species.md`](../../plans/active/thermophysics-cea-mole-fraction-species.md)
- **stage**: `result` (diff base `4b357400~1`)
- **date**: 2026-09-16
- **commit**: `417a3d2f` (feature/sern-design)
- **codex**: effort `high`, 3.6 min, rc=0
- **判定**: **GO**, 指摘 C0/M0/m0
- **focus**: 6 回目 (result-5 の M1/m1 反映後; §6.1 の採否表と §5.1 #30)。前回までの再指摘は不要。実物で確認: M1 interp_field.py の Xi→roXi 転送 (tools/interp_field.py; tests/unit/test_interp_field_tracer.py 3 件) と転送配列の完備検査, m1 plans/README と case/44 README の記述。accepted に移してよいか。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

## 判定: **GO**

result-5 の M1/m1 は修正済みです。転送処理・失敗分岐・run の再判定で、今回の受入れを妨げる追加指摘はありません。  
**改定済みの node/Euler 検証範囲で `accepted` へ移してよい**と判断します。

## 指摘一覧

新規指摘なし。修正済みの過去指摘は再掲しません。

### 確認した根拠

- **`Xi → roXi` 転送:** [interp_field.py:141](/home/sano/work/forge/solver_density_cuda/tools/interp_field.py:141) の補完処理を実行し、宛先の `roXi/ρ` が **0.8 → 0.3** になることを確認しました。[追加試験](/home/sano/work/forge/solver_density_cuda/tests/unit/test_interp_field_tracer.py:24)の3項目は **ALL PASS**。`roXi` が存在する場合は、その保存量が優先されます。
- **転送配列の完備検査:** [interp_field.py:152](/home/sano/work/forge/solver_density_cuda/tools/interp_field.py:152) は宛先を書込み用に開く前に検査します。トレーサ欠落、および入力検査を通っても転送配列に `roY1` が不足する条件を追加確認し、**拒否・宛先不変**を確認しました。
- **文書:** [plans/README.md:27](/home/sano/work/forge/plans/README.md:27) の古い残作業記述と、[case/44 README:645](/home/sano/work/forge/case/44.vitiated_air_wt/README.md:645) の固定点に関する断定は修正済みです。`methods/index.md` の参照先も存在します。plan の構造検査は **`VERDICT: PASS`**。

### run の再検証

| 対象・run 索引 | 確認結果 |
|---|---|
| `case/16.nozzle_wys/run_0471_tp5_node_euler_va3/`、`run_0472_tp2_node_euler_va3_lumped/` — [run 一覧](/home/sano/work/forge/case/16.nozzle_wys/README.md:373) | 両方 **`PASS (converged)`**。`machmax,pmax` は **`OVERALL: ALL STEADY`**。保存済み `MESH_QUALITY.txt` は **`VERDICT: PASS`** |
| case/44 の `run_0195`・`0196`・`0197`・`0199` — [対象一覧](/home/sano/work/forge/case/44.vitiated_air_wt/README.md:635) | 全て **`NOT CONVERGED`**。`cond_series.csv` は0.2%閾値で **`OVERALL: ALL STEADY`** |
| case/46 の `run_0100`–`0103_species_regress_euler_*` — [run 一覧](/home/sano/work/forge/case/46.sern_design/README.md:144) | 全て **`NOT CONVERGED`**。`force_history.csv` の `C_T,C_L,C_M` は0.1%閾値で **`OVERALL: ALL STEADY`** |

再判定は保存済み `CONVERGENCE_VERDICT.txt` と整合しました。上記10 run の全保存スナップショット、計166枚の `VALUE/*` に **NaN/Inf はありません**。

小型ケースの最終場比較も再計算し、相対差は `ro` **1.02e-6**、`P` **1.73e-6**、`T` **1.04e-6**、`Y_H2O` **3.85e-6**。plan の記載を裏付けます。case/44・46 の評価は、改定ゲートによる報告量の比較に限定されます。

## 推奨

**今回のレビューを §6.1 に記録し、通常の完了処理として `status: done`・`accepted` への移動・索引同期を行ってください。** トレーサ拡散と粘性二流体試験は、明記済みの §5.1 #8 `F-sp1` として残してください。

対象は `HEAD=417a3d2f`。指定 diff を取得して確認しました。試験の一時入出力はメモリへ置換し、新規ビルド・CFD投入は行っていません。ファイル変更なし、**plan 未反映**です。

指摘数: Critical 0 / Major 0 / Minor 0
