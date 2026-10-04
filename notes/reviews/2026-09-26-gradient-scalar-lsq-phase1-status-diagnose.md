# codex 諮問 (diagnose): gradient-scalar-lsq-phase1-status

- **brief**: [`notes/reviews/briefs/2026-09-26-gradient-scalar-lsq-phase1-status.md`](../../notes/reviews/briefs/2026-09-26-gradient-scalar-lsq-phase1-status.md)
- **plan**: [`plans/active/gradient-scalar-lsq-unification.md`](../../plans/active/gradient-scalar-lsq-unification.md)
- **date**: 2026-09-26
- **commit**: `3da913c8` (feature/sern-design)
- **codex**: effort `high`, 2.2 min, rc=0
- **結論**: **Phase 1 を `in_progress`・既定 `gg` のまま保ち、result 再審査の前に FCT smoke の実装直前版／現行版 gg の200物理 step A/Bを行う。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表です。参照する plan は [gradient-scalar-lsq-unification.md](/home/sano/work/forge-sern-design/plans/active/gradient-scalar-lsq-unification.md) です。

| 対象 | 採否・重大度 | 根拠と対案 |
|---|---|---|
| Phase 1 合格を求める result 再審査 | **却下・Major** | S2 の必須条件が残る。0476 継続の `rms_roY1` は起点床の約64倍、case/40 の `rms_roOmega` は約56倍で、`NOT BACK ON FLOOR`（[FLOOR_RECHECK.txt:83](/home/sano/work/forge-sern-design/case/09.Taylor-Green/_g0_lsq_seam/FLOOR_RECHECK.txt:83)、[同:169](/home/sano/work/forge-sern-design/case/09.Taylor-Green/_g0_lsq_seam/FLOOR_RECHECK.txt:169)）。未合格の進捗審査として提出することは可能だが、今それを次の一手には推奨しない。 |
| FCT smoke を Phase 1 完了条件から外す | **却下・Major** | 200 step の `E/D` は `ro` 1.975、`roUy` 2.032、`roY0` 0.229 で全て規定の0.1未満を満たさず、`VERDICT #5h(b): B`（[M3_fct.txt:18](/home/sano/work/forge-sern-design/case/09.Taylor-Green/_g0_lsq_seam/M3_fct.txt:18)）。さらに収支 FAIL が記録されている。**試験を実施したことと、検証対象を合格させたことを分け、#5h は「試験完了・ゲート未合格」とする。** |
| 床移動の真因調査を別 plan に分ける | **採用。ただし S2 の免除は却下・Major** | 旧版比較は case/48 の **gg 対 gg**。他ケースにも lsq にも一般化できない。実際、case/48 延長段の起点床比は `rms_roe` が gg 1.876、lsq 2.662（[FLOOR_RECHECK.txt:28](/home/sano/work/forge-sern-design/case/09.Taylor-Green/_g0_lsq_seam/FLOOR_RECHECK.txt:28)）。§5.1 #5l に調査先をリンクし、#5・§6 の未合格条件を残す。§4.4 には「調査の移管は Phase 1 GO を意味しない」と明記する。 |
| `in_progress` を維持して #2g を起票する | **採用** | #2g は既定化の独立した前提（[plan:67](/home/sano/work/forge-sern-design/plans/active/gradient-scalar-lsq-unification.md:67)）。起票・準備は可能。ただし Phase 2 の既定変更への移行とは区別する。今回の診断として優先するのは下記 A/B。 |
| 「物理量を変えず、性能劣化 ≤2%」という主張 | **却下・Major** | 差が許容内であることは不変の証明ではない。性能の事前上限も5%。記録上は `VERDICT S3: PASS`、中央値比0.9974／1.0184（[S3_case48.txt:10](/home/sano/work/forge-sern-design/case/09.Taylor-Green/_g0_lsq_seam/S3_case48.txt:10)、[S3_case39.txt:10](/home/sano/work/forge-sern-design/case/09.Taylor-Green/_g0_lsq_seam/S3_case39.txt:10)）。**「測定した2ケースで S3 を通過したが、S2 は収束・FCT ゲート未達で Phase 1 は未合格」**を提出文とする。 |

結論: **Phase 1 を `in_progress`・既定 `gg` のまま保ち、result 再審査の前に FCT smoke の実装直前版／現行版 gg の200物理 step A/Bを行う。**

第 1 仮説: FCT の収支 FAIL は、本 plan の実装以前から再現する共通の計算手順または診断上の問題である。 **確度: 中。ただし実装以前からの再現は未確認。**  
　根拠: [case/16.nozzle_wys/README.md:393](/home/sano/work/forge-sern-design/case/16.nozzle_wys/README.md:393) は gg／lsq × nSub 15／30 の4本とも `remainder 8.45e-2`、収支判定 FAIL と記録している。これは共通要因を疑う根拠にはなるが、定常場からの切替が原因だとは示さない。  
　反証条件: 入力・実効設定・診断定義を揃えた旧版 gg が保存ゲートを通過し、現行 gg だけが同じゲートで失敗する。

第 2 仮説: 本 plan の commit 群が gg にも共通する経路を変え、FCT の FAIL を生んだ。**確度: 低・未確認。** case/48 の旧版比較では FCT 経路を検査していないため、除外できない。

判別 A/B: **変更因子はバイナリだけ**（`36d8ba03`／`f99f236d`）。`run_0954_sglsq_s2_fct_gg` の開始入力を新規 run に複製し、両方 gg、nSub 15、dt・IC・BC・GPU・BLOCKSIZE を固定して各200物理 step。旧版に存在しない `scalarGradient` キーは両方で省略し、現行の実効 gg を確認する。同じ checker・同じ閾値で、収支閉合、独立総量照合、remainder と有界性を比較する。  
　→ **A: 旧版も同じ remainder 項で FAIL**なら、「今回の commit 群が初めてこの FAIL を生んだ」という仮説を棄却する。ただし lsq の精度ゲートは未合格のまま。  
　→ **B: 旧版 PASS・現行 FAIL**なら、既存問題として切り離す判断を棄却し、commit 群の共通経路へ戻る。  
　診断列の意味・実効設定が揃わない場合は判定不能とし、数値の近さで代用しない。

やらない方がよいこと: 未合格項目を Phase 2 に移して Phase 1 を GO にすること、現行 gg の床への差し替え、nSub を増やして合格が出るまで試すこと、case/48 の結果で全ケースの床移動を除外すること。

呼び出し側の前提への異議: **「gg でも同値」は「スカラー勾配とは独立」の証明ではない。** 共通の診断誤差や、gg／lsq 双方に作用する変更も残る。また、`E_q/D_q ≈ 2` は反復感度が大きい証拠であり、lsq−gg 差の原因を定量的に説明した割合ではない。

不足情報: 指定範囲の [M2_case48.txt:34](/home/sano/work/forge-sern-design/case/09.Taylor-Green/_g0_lsq_seam/M2_case48.txt:34) は条件(1)(2)だけを判定し、STEADY・物理ゲート・NaN/RISING は別出力で確認するよう明記している。その原判定出力、および FCT 4本の `check_passive_budget` 原出力・remainder の定義と閾値が必要。今回の閲覧範囲ではコード上の真因まで確認していない。

**plan 未反映**。変更禁止のため、上記の反映は呼び出し側で行う。
