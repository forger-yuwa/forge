# codex レビュー: tooling-nozzle-deltastar-band-root-zone (result)

- **plan**: [`plans/active/tooling-nozzle-deltastar-band-root-zone.md`](../../plans/active/tooling-nozzle-deltastar-band-root-zone.md)
- **stage**: `result` (diff base `2dda16c5`)
- **date**: 2026-10-13
- **commit**: `8acac69f` (feature/nozzle-deliverables)
- **codex**: effort `high`, 8.0 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M3/m2
- **extra**: `case/70.m5_air_wt_demo/README.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**

主要修正の効果は確認できましたが、参照流束の検査に端点の漏れがあります。
M5 の V5 と成果物の寸法検証も未完了で、現在の完了条件のまま `accepted` へ移すことには反対します。

1. **Major — `nonpos_ref` が帯下端 y_b を検査していない**

   根拠: [deltastar.py:276](/home/sano/work/forge-deliv/design/forge_design/metrics/deltastar.py:276)。検査対象は `qref[y <= y_b]` の格子点だけです。y_b が格子点間にある場合、その端点で参照流束が負でも通過します。D < 0 では累積流量の検査にも入りません。

   読み取り専用の合成試験で、既定の `refine=4`、r_w=1、y_b=0.10003、q_E=1、帯の参照 q_ref=−0.001−100(y−y_b) とすると、次を再現しました。

   ```text
   q_ref(y_b) = -0.0010000000000005926
   core_rms   = 0
   band_rms   = 1.34e-15
   delta_r    = -0.12195253943731954
   reason     = ['negative_deficit']
   ```

   `negative_deficit` は soft なので、[deltastar.py:589](/home/sano/work/forge-deliv/design/forge_design/metrics/deltastar.py:589) 以降で平滑化に使われます。実 run での発生は未確認ですが、「探索域の参照が非正なら hard」という仕様を満たしていません。

   **対案:** 符号分岐の前に y_b での参照値も評価し、閉区間全体の有限性・正値性を検査する。格子点間でゼロを横切る負例を追加し、`nonpos_ref` → NaN → held → 寸法停止まで確認してください。

2. **Major — V5 は未達。保存された `STEADY` は専用受入れ判定になっていない**

   根拠: [qs_series.py:43](/home/sano/work/forge-deliv/case/70.m5_air_wt_demo/qs_series.py:43) は末尾5枚だけを、既定の drift 5%・osc 10% で判定します。[plan:170](/home/sano/work/forge-deliv/plans/active/tooling-nozzle-deltastar-band-root-zone.md:170) の「変化・推定残差とも ≤2.5e−4」を実装していません。

   以下はすべて `case/70.m5_air_wt_demo/` 配下です。最終 HDF5 から出口 M を再計算し、CSV と一致することを確認しました。

   | run | 最終出口 M | 末尾40%の M の最大−最小 |
   |---|---:|---:|
   | `run_0014_ns_ns2` | 5.002198385 | 1.87045e−3 |
   | `run_0014b_ns_ns2_cont60k` | 5.001334078 | 3.26299e−4 |
   | `run_0014c_ns_ns2_cont120k` | 4.999110293 | 6.73212e−4 |
   | `run_0014d_ns_ns2_cfl5` | 4.997856084 | 1.11755e−3 |

   監査用に末尾40%、相対許容5e−5（M≈5で絶対2.5e−4相当）として `check_quasisteady` で再判定すると、**4 run とも `DRIFTING`** でした。

   残差については `run_0014` の `PASS (converged)` を再現しました。継続 run の b/c は、規定の `--from-floor run_0014_ns_ns2` なら **`PASS (within 1.5x of reference floor)`**、d は **`NOT CONVERGED`** です。したがって b/c の問題は、継続区間の低下桁数ではなく、目的量が動き続けることです。

   現 plan の「V5 判定保留」は実測と整合しています。しかし [完了条件:194](/home/sano/work/forge-deliv/plans/active/tooling-nozzle-deltastar-band-root-zone.md:194) は V1〜V5 合格のままです。−0.19%→+0.04% を定常解の改善実績として採用できません。

   **対案:** M5 は未達として保持し、実装修正の完了範囲から分離する。V5 用には絶対許容・判定窓・残りの変化の推定を明示した判定結果を保存してください。

3. **Major — 「寸法ゲート PASS」の根拠が、要求された成果物ゲートではない**

   根拠: [problem_m5_air_ns2.yaml.card_log.json:699](/home/sano/work/forge-deliv/case/70.m5_air_wt_demo/problem_m5_air_ns2.yaml.card_log.json:699) は、記録を「設計壁から測った寸法」「最終の合格は成果物（STEP・メッシュ）から測る」と明記しています。

   拡大部長 **1.199522039 m、誤差 −0.477961 mm** と設計寸法の PASS は確認できました。しかし `run_0014*` の4ディレクトリには、いずれも `deliverables/` と `DIMENSIONS.json` がありません。[plan:173](/home/sano/work/forge-deliv/plans/active/tooling-nozzle-deltastar-band-root-zone.md:173) の成果物寸法ゲートを満たす証拠はありません。

   **対案:** 記述を「カード生成時の設計寸法 PASS」に限定する。成果物まで完了範囲に含めるなら、対象 run の STEP・メッシュを生成して出荷寸法ゲートを通し、run 一覧に根拠を追加してください。

4. **Minor — 運用手順が新経路と再現性の変更に追随していない**

   根拠: [nozzle-design-workflow.md:88](/home/sano/work/forge-deliv/procedures/nozzle-design-workflow.md:88) は従来の `solve_rt_throat --prev` の案内だけで、実寸を維持する `wt_card --ns-correction` → `--delta-r-csv` の経路がありません。また [同:213](/home/sano/work/forge-deliv/procedures/nozzle-design-workflow.md:213) の「同じ壁を再現したいだけなら YAML を変えない」は今回の変更後には不十分です。

   保存 A/B では、同じ YAML でも出口が c2b で **+0.845760 mm**、case/45 生産で **+5.353145 mm** 動きます。これは意図した修正で、plan と `methods` には記録されています。

   **対案:** 新 CLI の手順、壁記録の必須条件、旧 run の再現にはコード版・保存壁の確認が必要なことを手順書へ反映する。`methods/index.md` の既存リンクは有効でした。

5. **Minor — 比較条件と残作業表に古い記述が残っている**

   根拠: [plan:171](/home/sano/work/forge-deliv/plans/active/tooling-nozzle-deltastar-band-root-zone.md:171) は「同じ IC の A/B ではない」と訂正していますが、174行には「同じ IC・段取り・step 数」が残っています。実際の IC は [run_0010/prepare_info.json:636](/home/sano/work/forge-deliv/case/70.m5_air_wt_demo/run_0010_ns_c2b_final/prepare_info.json:636) が `run_0009`、[run_0014/prepare_info.json:606](/home/sano/work/forge-deliv/case/70.m5_air_wt_demo/run_0014_ns_ns2/prepare_info.json:606) が `run_0010b` です。

   §5.1 の #10・#11 に未解決事項が残っている点は適切です。一方、#8/#12 の result 段が重複し、#12 の「諮問未実施」は §6.1 の記録と矛盾します。

   **対案:** 比較は異なる初期場からの計算であると統一し、残作業を「実装の修正」「検証未達」「後続課題」に整理してください。

再確認できた肯定的な結果は、V1 の解析解との差 **1.29e−7**、既存 `run_deltastar_tests.py` と壁単調性試験の **ALL PASS**、case/45 の4対の最終 CSV のバイト同一です。単調性 A/B の33壁・旧判定不合格14壁という集計も一致しました。対象7 run の最終 HDF5 に NaN/Inf はありません。ソルバ本体の diff はなく、境界半割面・周期 seam・軸の離散式を変更した差分はありません。

**推奨は、現 plan を `paused` のまま保持し、端点検査を修正したうえで、検証済みの実装修正を独立した完了範囲として再レビューすることです。** M5 の定常性能・成果物検証・帯選択の課題は未完として引き継いでください。今回はファイルを変更しておらず、この推奨は plan 未反映です。

指摘数: Critical 0 / Major 3 / Minor 2
