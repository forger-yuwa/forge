# codex レビュー: time_integration-line-viscous-jacobian (result)

- **plan**: [`plans/active/time_integration-line-viscous-jacobian.md`](../../plans/active/time_integration-line-viscous-jacobian.md)
- **stage**: `result` (diff base `4d394a71`)
- **date**: 2026-10-10
- **commit**: `3de780a4` (feature/faceh-audit-viscjac-close)
- **codex**: effort `high`, 3.3 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M3/m1
- **focus**: この plan を「本線不採用・診断として残置」で accepted に移してよいかを判定する。対象は §6.17 (閉じるときの処置: 試験ごとの処置・診断コードの残置の範囲・再開の条件・U0 と FP32 の扱い) と、§5.1 #4・#5g・§6.3 の訂正。methods/time_integration/implementation.md の lineViscCoupling 2・3 の節と procedures/solver-settings.md の lineViscCoupling の記述は作業ツリーの未 commit の変更なので、ファイルを直接読むこと。諮問 notes/reviews/2026-10-09-line-viscous-faceh-diagnose.md が閉じる前に書くよう求めた項目 (不合格・判定不能・未実施の処置、診断コードの残置範囲、再開条件、U0・FP32 の適用範囲と移管先、古い要約の訂正) を満たしているか。U0 の代わりに挙げた速度 plan §6.0 の B1 の短期の一致 (run_0276〜0281、AB_abB.json は AWS 上にありリポジトリには無い。値は §6.17 に転記) を代わりの証拠と呼んでよいか。総時間の引用 (速度 plan §6.14・§6.15・§6.20) が原文と合っているか。コードの差分は 5ab83056・6d49738e・4fdc2c0e・846727be・5197e00e を git show で見る (4d394a71..HEAD には速度 plan の変更も混ざる)。コードは変えない方針。
- **extra**: `notes/reviews/2026-10-09-line-viscous-faceh-diagnose.md`, `plans/active/time_integration-line-implicit-speed.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**

「本線不採用・診断として残置」で閉じる方針は妥当です。
ただし、検証の合格範囲と古い要約に不整合が残っています。以下を優先順に訂正してから `accepted` へ移してください。

1. **Major — U-J の PASS は、登録した列別精度の合格を示していません。**  
   [事前登録](/home/sano/work/forge-faceh/plans/active/time_integration-line-viscous-jacobian.md:115)は「列ごとの相対誤差 ≤ 1e-6」ですが、[試験実装](/home/sano/work/forge-faceh/solver_density_cuda/tools/test_line_visc_jacobian.cpp:69)は全列の誤差を行列全体の最大値 `sD/sK` で割っています。小さい列の相対誤差は、この判定では保証されません。  
   **対案:** コードを変えずに閉じるなら、§6.0・§6.17 を「実装した全行列正規化の試験は PASS、登録した列別判定は未確認」と訂正し、列別照合を再開時の課題に残してください。失敗と断定する必要もありません。

2. **Major — U2 の「収束した」を裏付ける VERDICT が提示されていません。**  
   [§6.0](/home/sano/work/forge-faceh/plans/active/time_integration-line-viscous-jacobian.md:146)、[§6.17](/home/sano/work/forge-faceh/plans/active/time_integration-line-viscous-jacobian.md:436)、[case/52 の索引](/home/sano/work/forge-faceh/case/52.conjugate_slab/README.md:49)は、`run_0011_u2_lvc2_cfl50` の温度誤差を根拠に「収束」としています。確認できる記録は、950 step で 0.05 K、5000 step で 1.4e−6 K という解析解誤差です。対象 run の全残差に対する判定は確認できませんでした。  
   **対案:** 「温度の解析解誤差が上記水準まで低下。250 step の基準は不合格、残差収束は未確認」に統一してください。「収束」を残す場合は、判定区間付きの `check_convergence` VERDICT が必要です。

3. **Major — 本文で限定した因果結論が、恒久索引と残作業表では断定のままです。**  
   [case/45 README:173](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/README.md:173)は E1 を「熱伝導の近傍 K を除くと破綻を避けられる」と記載していますが、[§6.9](/home/sano/work/forge-faceh/plans/active/time_integration-line-viscous-jacobian.md:285)の正式な扱いは「比較無効・探索的な観測」です。  
   また、[§5.1 #5i](/home/sano/work/forge-faceh/plans/active/time_integration-line-viscous-jacobian.md:103)と README:174 の「密度の列は必要でない」は、[§6.11 の訂正](/home/sano/work/forge-faceh/plans/active/time_integration-line-viscous-jacobian.md:313)にある「同じ破綻機構かは未確認」という限定を落としています。  
   **対案:** E1 は比較無効を明記し、密度列の試験は「除去しても 566 step で非有限。同一機構への必要性は未確定」に統一してください。今回の #5g・§6.3 の訂正自体は適切です。

4. **Minor — 総時間の引用には、比較条件と未収束の留保も残してください。**  
   [§6.17](/home/sano/work/forge-faceh/plans/active/time_integration-line-viscous-jacobian.md:422)の数値は原文と整合します。再計算でも L5 は 1.17875 h、L0 は 1.15967 h、差は **68.7 秒**です。  
   ただし、[速度 plan §6.14](/home/sano/work/forge-faceh/plans/active/time_integration-line-implicit-speed.md:391)には L5 の両段が `NOT CONVERGED`、到達窓の θ_r が `TRANSIENT-UNSETTLED` とあります。また、[§6.20 の B0](/home/sano/work/forge-faceh/plans/active/time_integration-line-implicit-speed.md:555)は **2 出力連続**の到達条件です。  
   **対案:** 「専有単価で換算した到達時間の推定であり、収束解の同等性比較ではない」「B0 は到達条件が異なる」を引用先にも添えてください。

指定の 5 commit、現在のカーネル、未 commit の文書を確認した範囲では、D/K の符号、壁拘束、マスク、値 0/1 との分岐に新たな明白な実装誤りは見つかりませんでした。周期 seam・軸を含む実行時の健全性まで確認した意味ではありません。

**B1 は、値 0・キー 5・方向別・上限 50 に限った補助証拠として引用できます。** [比較スクリプト](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/checkB.sh:23)と commit の包含関係は説明どおりです。ただし `AB_abB.json` は未閲覧なので、転記された合格結果の独立検証はできていません。U0 全体や現在の HEAD の回帰合格には拡張せず、U0 未実施の扱いを維持してください。

面エンタルピーについては、別作業ツリーに残る [集計 JSON](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_jprobe/fh_judge.json:1)で有効条件と相対差比 **0.39934** を確認できました。ただし元配列からの再計算はしていません。対象 run の残差・VERDICT が本作業ツリーにないため、今回 `check_convergence`／`check_quasisteady` の再実行はできていません。

**推奨は、上記の記録を訂正して「本線不採用・診断残置」の判断文書として閉じることです。** U0・FP32・E2 の未実施、方向微分の判別不能、残差精度監査の移管は維持します。移動時には `status`・変更ログ・result 採否・`plans/README.md` と相互リンクを同期してください。ファイルは変更していません。plan 未反映です。

指摘数: Critical 0 / Major 3 / Minor 1
