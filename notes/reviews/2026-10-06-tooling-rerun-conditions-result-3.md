# codex レビュー: tooling-rerun-conditions (result)

- **plan**: [`plans/active/tooling-rerun-conditions.md`](../../plans/active/tooling-rerun-conditions.md)
- **stage**: `result` (diff base `main`)
- **date**: 2026-10-06
- **commit**: `7936b118` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 8.3 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M3/m2

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**NO-GO**

主要な run 数値と準定常判定は再現できましたが、実装には入力解釈と段終了検査の欠落が残っています。  
`accepted` への移動は保留し、以下の Major 3 件を先に修正してください。ファイルは変更していません。

1. **Major — 重複 YAML キーを受理し、ソルバと異なる保存量集合をスケールする**

   根拠: [rerun_conditions.py:365](/home/sano/work/forge-integ-1005/solver_density_cuda/tools/rerun_conditions.py:365)、[同:600](/home/sano/work/forge-integ-1005/solver_density_cuda/tools/rerun_conditions.py:600)。`yaml.safe_load` は重複キーの後勝ちですが、ソルバの `yaml-cpp` は先勝ちです。この差は既存の [migrate_solver_config.py:17](/home/sano/work/forge-integ-1005/solver_density_cuda/tools/migrate_solver_config.py:17) にも明記されています。

   `run_0094` の設定をメモリ上で `turbulence: {model: "sst", model: "none", …}` に変更すると、`build_plan` は Pt 0.8 倍のスケールを受理し、必要保存量から `roK`・`roOmega` を除外しました。ソルバは `sst` と解釈するため、生成場では ρ だけが 0.8 倍になり、k・ω は意図せず **1.25 倍**になります。「T・U・Y・k・ω を保つ」という契約を破ります。

   **対案:** 全階層の重複キーを作成前に拒否する共通 YAML ローダーを導入し、`rerun_conditions` と runner の設定検査に適用してください。既存の重複検出方式を利用できます。

2. **Major — NS/Euler 分類が粘性・熱伝導・種拡散を見ていない**

   根拠: [rerun_conditions.py:678](/home/sano/work/forge-integ-1005/solver_density_cuda/tools/rerun_conditions.py:678)。分類は壁種別と `sst*` だけです。しかし、全壁 `slip`・乱流なしでも、`viscMethod: 2` なら内部の粘性・種拡散は残ります。種拡散の実行条件は [speciesTransport_d.cu:960](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/speciesTransport_d.cu:960) です。

   読み取り内容だけを差し替えた再現では、`viscMethod: 2` と `transport` を残した全壁 `slip` の入力が **Euler と判定され、`stages: none` を推奨**しました。生成 config の本段 CFL 5・6000 step も警告だけで通ります。これは腕 E の非粘性検証を適用できる入力ではありません。

   **対案:** v1 の Euler 対応を、乱流・粘性・熱伝導・種拡散が無効であることを確認できる設定に限定してください。全壁 `slip` でも輸送が有効な入力は、未対応として作成前に拒否するのが適切です。

3. **Major — Euler の段階起動は Inf・負密度を段終了ゲートで止めない**

   根拠: [runner_axismach.py:732](/home/sano/work/forge-integ-1005/design/forge_design/evaluate/runner_axismach.py:732)。Euler の `_stage` は終了コードと保存 step だけを確認し、そのまま restart します。NS 側の [同:1159](/home/sano/work/forge-integ-1005/design/forge_design/evaluate/runner_axismach.py:1159) にある `stage_gate` 呼び出しがありません。

   メモリ内 HDF5 の `ro=-1`／`ro=Inf` は実際の `stage_gate` で不合格になります。一方、forge と書込みをモックした Euler `full` の制御試験では、**gate 呼出し 0 回、solver 呼出し 2 回、戻り値 0**でした。`restart_field` のビット一致検査は Inf・負密度を排除しません。[restart_field.py:118](/home/sano/work/forge-integ-1005/solver_density_cuda/tools/restart_field.py:118)

   **対案:** Euler にも restart 前の `stage_gate` を適用し、Inf・負密度で次段を呼ばない回帰試験を追加してください。新たに `full` を推奨する実行経路として、段別履歴と manifest も残すべきです。

4. **Minor — Pt＋壁温変更で、受理結果と「禁止」という説明が矛盾する**

   根拠: スケール拒否条件には Tw が含まれませんが、推奨文は Tw 変更を禁止条件と説明します。[rerun_conditions.py:568](/home/sano/work/forge-integ-1005/solver_density_cuda/tools/rerun_conditions.py:568)、[同:689](/home/sano/work/forge-integ-1005/solver_density_cuda/tools/rerun_conditions.py:689)、[overview.md:999](/home/sano/work/forge-integ-1005/methods/design/overview.md:999)

   `--Pt 4400000 --Ps 1789.6 --Tw 350 --scale-ic pt` はメモリ上の等温壁 fixture で受理され、同時に「scale-ic pt は禁止条件に当たる」と記録されました。

   **対案:** plan §4.4 に合わせ、Pt＋Tw の IC 変換は許容しつつ、説明を「複合条件の起動・整定は未検証で、推奨対象外」に統一してください。

5. **Minor — run 台帳に「残り約 0.1 %」という未確定の整定余量が残る**

   根拠: [case README:97](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/README.md:97)。plan 変更ログでは簡易予想値との差へ訂正されていますが、台帳は依然として「残り約 0.1 %」です。

   `run_0138_rerun_fullpath_blk4/quantities_series.csv` の流量末尾平均は **17136.9285**、直前窓差は **+0.821875**、再判定は **DRIFTING**。17150 との差は真の漸近値までの距離ではありません。

   **対案:** 台帳も「簡易予想値との差約 0.1 %、整定余量は未確定」に訂正し、§5.1 #11 を参照させてください。

検証結果については、主ツリーの `case/45.isobutane_m6_d155/` にある17 run に `check_convergence.py` と、plan 所定閾値の `check_quasisteady.py` を再実行しました。主な結果は次のとおりです。

| run | 再確認した数値 | VERDICT |
|---|---|---|
| `run_0119_rerun_ctrl` | δ_E 0.725282185、出口 M 5.999261065 | 残差 `NOT CONVERGED`、対象3量 `STEADY` |
| `run_0120_rerun_euler_pt08` | 出口 M 5.999997716、流量 13344.9650 | 残差 `NOT CONVERGED`、M・流量 `STEADY` |
| `run_0132_rerun_pt08_scale_cfl1_ext` | δ_E 0.749448414、出口 M 5.992124197 | 残差 `NOT CONVERGED`、対象3量 `STEADY` |
| `run_0134_rerun_euler_tt1500` | 出口 M 6.021834240、直前窓差 +0.020794759 | 残差 `NOT CONVERGED`、M・流量 `DRIFTING` |
| `run_0139_rerun_pt08_noscale_cfl1_ext` | 出口 M 5.993976711、最終場の Ux<0 は779節点 | 残差 `NOT CONVERGED`、M・流量 `DRIFTING` |

数値は各 run の `quantities_series.csv`、逆流節点数は `res_6000.h5` を確認しています。A3 と B3 の結果は、限定された Pt 変更条件でスケール IC を推奨する根拠になります。**同一定常解への到達や解の一意性を示した結果ではありません。** plan の最新の限定表現は妥当です。

回帰確認では `run_mesh_params_tests.py` は **FAIL 0**、plan lint は **PASS**。ファイル作成を伴う `test_rerun_conditions.py` 全体は再実行せず、上記の検査をメモリ上で実施しました。§5.1 #9・#11には物理受入れと整定長の未解決事項が残っており、この持越し自体は明示されています。

**推奨は、移動を保留し、指摘1→2→3の順に実装修正と回帰試験を行い、その後に文書を同期して再レビューすることです。** 既存の Pt スケール検証結果を撤回する必要はありません。本レビューによる変更提案は、read-only 指定に従い **plan 未反映**です。

指摘数: Critical 0 / Major 3 / Minor 2
