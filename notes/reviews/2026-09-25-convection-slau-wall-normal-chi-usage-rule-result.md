# codex レビュー: convection-slau-wall-normal-chi-usage-rule (result)

- **plan**: [`plans/accepted/convection-slau-wall-normal-chi-usage-rule.md`](../../plans/accepted/convection-slau-wall-normal-chi-usage-rule.md)
- **stage**: `result` (diff base `093bcdca`)
- **date**: 2026-09-25
- **commit**: `910d1514` (feature/sern-design)
- **codex**: effort `high`, 5.6 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M4/m2
- **focus**: 結論が測った範囲 (m6_on・生産格子・末尾窓の帯判定 / 1 dump・1 構成の診断可能性 / 衝撃足は両作動点で判定保留) を超えていないか、判定保留のまま閉じる扱いと t=2 mm vs cowl_thickness 0.5 mm の食い違いの記録が適切か、procedures/recommended-settings.md §1.0a の規則文と注記が §6.2 の結果と整合するか。
- **extra**: `case/46.sern_design/_rule_usage/BAND_VERDICT.txt`, `case/46.sern_design/_rule_usage/APPLICABILITY_run0437_res0.txt`, `case/46.sern_design/_rule_usage/QS_flag0.txt`, `case/46.sern_design/_rule_usage/QS_flag1.txt`, `case/46.sern_design/_rule_usage/sf_diag_m6.txt`, `case/46.sern_design/_rule_usage/sf_diag_m4.txt`, `case/46.sern_design/_rule_usage/m6_on_VERDICT.json`, `case/46.sern_design/_rule_usage/_rule_usage_sf/m4_off_VERDICT.json`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
m6_on の5列の帯判定は元 CSV から再現でき、衝撃足を「判定保留」とする結論も維持できます。  
ただし、診断ツールの誤合格、証拠を超えた適用条件の断定、未完了の再判定を直してから移すべきです。

1. **Major — 診断ツールが NaN・照合面ゼロを合格にする。**  
   根拠: [diag_applicability.py:139](/home/sano/work/forge-sern-design/case/46.sern_design/cad/diag_applicability.py:139) 以降は `r > 1` のみで失敗を判定し、非有限値や照合件数を検査していません。ファイルを書かない入力模擬試験で、**流束 NaN と接続面0件の両方が exit 0、`VERDICT: 診断可能`、最大誤差0.000** になりました。  
   また、142行で境界半割面を除外しています。既存の [run0437_budget.txt:4](/home/sano/work/forge-sern-design/case/46.sern_design/_r3_m1m3/run0437_budget.txt:4) では対象各 CV は6面ですが、新しい照合は各5面です。「全接続面」の確認にはなっていません。  
   **対案:** 状態・面積・許容差・誤差の有限性、正の密度・音速、対象 CV の妥当性、照合面数を必須検査にする。境界半割面は境界種別と質量流束ゼロを確認して明示的に扱い、NaN・空集合・欠落面の否定試験を追加する。

2. **Major — §4.1.1 の設定検査を実装し切れていない。**  
   根拠: [diag_applicability.py:67](/home/sano/work/forge-sern-design/case/46.sern_design/cad/diag_applicability.py:67) はログが無くても継続し、95行の `thermalMethod` は無条件で OK。周期境界の検査もありません。模擬試験では、**YAML の `thermalMethod: 2` とログの `0` の不一致、ログ欠落、旧配置 `physProp.isAxisymmetric: 1` がすべて `bad=[]`** でした。旧配置は実際に [solverConfig.cpp:881](/home/sano/work/forge-sern-design/solver_density_cuda/input/solverConfig.cpp:881) が受理します。  
   **対案:** 実ソルバと同じ優先順位で実効設定を解決し、周期・軸対称を拒否する。ログ欠落と不一致は診断不能にする。エコーを出さないキーは、その制約と代替の確認根拠を記録し、「両方から確認済み」と扱わない。

3. **Major — 適用規則の根拠が測定範囲を超えている。**  
   根拠: [plan:221](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi-usage-rule.md:221) は1 dump の検査から「診断3条件が成立」としています。しかし [§4.1:55](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi-usage-rule.md:55) の条件(i)(ii)には3 dump以上の密度減少・正味流出の持続が必要です。今回の表が裏付けるのは、その瞬間の密度比、面流束照合、条件(iii)です。  
   また、[recommended-settings.md:85](/home/sano/work/forge-sern-design/procedures/recommended-settings.md:85) の「2D では排出の3条件が成立しない」は、力係数の帯判定からは導けません。  
   **対案:** 結論を「1 dumpで診断可能性と条件(iii)を確認」に限定する。2Dの規則は「現行運用は0、1を適用する根拠は未確認」とし、帯判定には **m6_on・生産格子・step 22000–36000** の限定を付ける。過去の時系列を根拠にするなら、その flag・設定・3 dump の数値を直接参照する。

4. **Major — 膨張角窓の再判定が未完了なのに、残作業 #2 が完了になっている。**  
   根拠: [plan:99](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi-usage-rule.md:99) と [#2:155](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi-usage-rule.md:155) は、`C_T` **と膨張角窓**の再計算を要求しています。結果と前 plan の今回の訂正は `C_T` のみで、[前 plan:980](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi.md:980) には撤回対象の振幅加算による **0.857 %、判定不変** が残っています。  
   **対案:** 既存 dump で膨張角窓の全 dump 対の分布比較を再計算し、現行の受入記述と [case README:376](/home/sano/work/forge-sern-design/case/46.sern_design/README.md:376) を同期する。それまでは #2 の該当部分を未完了に戻す。旧定義の「カウル板厚2 mm」も、**実厚0.5 mmとは異なる登録済み評価長さ**と注記する。今回の登録値を事後変更しなかった判断自体は妥当です。

5. **Minor — 衝撃足が検出された場合の分岐に、登録手順との相違がある。**  
   根拠: [v3sern_shockfoot.py:130](/home/sano/work/forge-sern-design/case/46.sern_design/v3sern_shockfoot.py:130) で末尾40 %に絞った後、[同:102](/home/sano/work/forge-sern-design/case/46.sern_design/v3sern_shockfoot.py:102) で再び `--tail 0.4` を適用します。72 dumpなら準定常判定は29点ではなく12点です。人工系列で、29点では `DRIFTING`、再切出しした12点では `STEADY` を再現しました。さらに [同:70](/home/sano/work/forge-sern-design/case/46.sern_design/v3sern_shockfoot.py:70) は非有限点を除外するため、50点中1点だけ有限でも唇由来率1.0になります。  
   **対案:** 窓選択を一度に統一し、50点の有効性を確認する。148行の未登録の「半数以上検出」条件も登録手順に合わせる。今回の両側0候補という結論は、これらの比較分岐を通らないため影響しません。

6. **Minor — `solver` の省略時挙動について plan と実装が食い違う。**  
   根拠: [plan:159](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi-usage-rule.md:159) は「省略 ≡ 既定 SLAU」を試験済みとしていますが、[solverConfig.cpp:253](/home/sano/work/forge-sern-design/solver_density_cuda/input/solverConfig.cpp:253) では必須キーです。[stage_manifest.py:114](/home/sano/work/forge-sern-design/solver_density_cuda/tools/stage_manifest.py:114) の説明は正しく、実行した10試験も `PASS` でした。省略と SLAU の同一性は試験していません。  
   **対案:** 実装を維持し、planを「必須キー、SLAU／ROE／SLAU2の区間分離と大文字小文字の正規化を確認」に訂正する。

検証上、5列の平均差・差区間は提示値を再現しました。`check_quasisteady.py` は flag0 の `C_T` 系を `STEADY`、`C_L` 系・`C_M` を `OSCILLATING`、flag1を全列 `STEADY` と判定しました。保存された `run_0450–0452` の収束判定はすべて **`NOT CONVERGED (stalled/plateau)`** です。ローカルで再実行できた m6_on の残差判定は延長前の11999 stepまでで、延長 run の残差と衝撃足・面流束の生データは独立再計算できていません。

**推奨は、測定範囲を限定した受入を維持し、上記を直してから `accepted` に移すことです。** 優先順は **#1・#2の診断ゲート修正 → #3の断定訂正 → #4の未完了再判定 → #5・#6の同期**。衝撃足は判定保留のまま閉じて構いませんが、§5.1に保留状態と§7の再開条件を残してください。3D格子収束・固定点の `R5o-chi` への委譲は記録されています。指定 diff に CUDA カーネルや既定値の変更はなく、ファイルも変更していません。

指摘数: Critical 0 / Major 4 / Minor 2
