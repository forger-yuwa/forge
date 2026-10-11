# codex レビュー: time_integration-line-viscous-jacobian-faceh (result)

- **plan**: [`plans/active/time_integration-line-viscous-jacobian-faceh.md`](../../plans/active/time_integration-line-viscous-jacobian-faceh.md)
- **stage**: `result` (diff base `3de780a4`)
- **date**: 2026-10-11
- **commit**: `274359b7` (feature/faceh-audit-viscjac-close)
- **codex**: effort `high`, 6.9 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M1/m1
- **focus**: plan を閉じて accepted/ へ移す前の result 段。§6.2・§6.4・§6.6・§6.8・§6.10・§6.13 の判定と変更ログの主張が、判定の記録 (case/45.isobutane_m6_d155/_band_ab/cold_pair/lvc*_judge.json・*_pregate.json) と判定器のコードで裏付けられているか。§5.1 #8〜#10 は後継 plan time_integration-line-viscous-jacobian-dt へ移した。既定はこの plan では決めない (ユーザ決定 2026-10-11)。実装の diff は timeIntegration_d.cu の診断用の切替 (FORGE_LVC_TERMS など) と試験・判定の道具
- **extra**: `case/45.isobutane_m6_d155/_band_ab/cold_pair/lvcgeom_judge.json`, `case/45.isobutane_m6_d155/_band_ab/cold_pair/lvc75_judge.json`, `plans/active/time_integration-line-viscous-jacobian-dt.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**

主要な結論は、確認できた実データ・判定記録と整合しています。通常ビルドの数値処理に今回の差分による変更はありません。accepted への移動前に、判定器の検査漏れと文書の旧記述を修正してください。

1. **Major — 監査判定器が、非有限の sweep 配列を検査済みとして通す**

   [lvcaudit_judge.py:137](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcaudit_judge.py:137) は `rhs_s*`・`dqnew_s*`・`dqold_s*` の存在を確認しますが、156 行の有限性検査には含めていません。さらに [同:366](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcaudit_judge.py:366) の T 判定は、同じビット列の NaN も一致として通します。

   保存済み入力をメモリ上で差し替え、監査・対照・再実行の `rhs_s1` に同一の NaN を入れると、**S は不合格 0、T の `rhs_s1` は `ok: true`** のままでした。入力不備を `INVALID` にする登録条件を満たしていません。

   実際の `run_0550`〜`0552` の対象配列には非有限値がなく、今回の FAIL・T 保留という結論は変わりません。ただし、[plan §6.8:297](/home/sano/work/forge-faceh/plans/active/time_integration-line-viscous-jacobian-faceh.md:297) の「次から検査する」は未実装で、後継の残作業表にも明示されていません。

   **対案:** 3 本すべての sweep 配列に期待形状 `(N, 5)` と有限性の検査を追加し、不備は比較前に `INVALID` にする。NaN・空配列の反例と保存済み正常入力で確認し、対応を残作業表・変更ログへ記録する。

2. **Minor — 現在仕様と索引が、今回の実装を反映していない**

   [implementation.md:543](/home/sano/work/forge-faceh/methods/time_integration/implementation.md:543) は、単体試験を「行列全体の最大値で割って判定」「列ごとの相対誤差は保証しない」と説明しています。現在の実装は列別判定を行い、`--ref mp100 --eta 1e-8` による参照精度の比較も備えています。

   また、[plans/README.md:33](/home/sano/work/forge-faceh/plans/README.md:33) は初期の面エンタルピー A/B だけを説明し、「コードは変えない」としたままです。

   **対案:** 現在仕様には列別判定・多倍長参照・保証範囲を記載する。索引には監査出力、幾何精度 A/B、マスク 7/5 比較までの結果を反映し、「通常ビルドの数値処理は変更しない」と正確に書く。

照合結果は次のとおりです。run はすべて `case/45.isobutane_m6_d155/` 以下で、恒久索引は [README.md:191](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/README.md:191) です。

| 対象 | 確認結果 |
|---|---|
| 実装・回帰 | `timeIntegration_d.cu` の追加された監査用 `#if` 16 ブロックを除去すると、`3de780a4` の内容と完全一致。既定値変更なし。 |
| §6.2 | JSON のゲート 45 項目合格、4 本 DIVERGED、ログの検知 step 123／122／123／123 は本文と一致。 |
| §6.4・§6.6 | 旧 FAIL の記録を保持。多倍長 CSV は 3560 列すべて合格、最大誤差 2.408e−9、旧不一致 13 列もすべて合格。本文と一致。 |
| §6.8 | 保存済み host 配列を使った Python 側の再照合でも、C の不一致 125、D の判別不能 3、T の外れ 9。係数差の最大 10.785 % も一致。 |
| §6.10 | `run_0560`〜`0563` の CSV 最初の非有限 step は 292／356／198／234。破綻前の保存場も確認し、本文と一致。 |
| §6.13 | `run_0570`・`0572` は CSV step 347／279 で非有限。`run_0571`・`0573` は各 2000 step 有限、各 21 保存時点の検査に異常なし。事前ゲートの証拠 53 ファイルもハッシュ一致。 |
| 残作業 | #8〜#10 の後継への移管を確認。C・D・T と実残差との整合は未解決として保持。上記の sweep 検査漏れだけ処置が不足。 |

マスク 5 の判定を保存済み区間 CSV（main、step 0〜1999）から再実行した結果は、両方とも **`NOT CONVERGED (still converging — run more steps)`** でした。`rms_roK` の末尾 500 step の傾きも +0.179／+0.228 桁で、plan の「有限だが未収束」という限定は妥当です。

確認範囲の制限として、`run_0540`〜`0543` の実体は本環境に見つからず、§6.2 は JSON・README の照合までです。C++ 単体試験はコードと保存出力・CSV を照合し、再ビルドはしていません。ファイルは変更していません。

**推奨:** 上記 1 → 2 の順で修正し、既存証拠による再確認を記録してから、診断の区切りとして accepted へ移す。追加の GPU 計算や既定値の決定は、このレビューへの対応には不要です。

指摘数: Critical 0 / Major 1 / Minor 1
