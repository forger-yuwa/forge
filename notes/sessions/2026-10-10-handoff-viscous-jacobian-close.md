# 引き継ぎ ④: 粘性ヤコビアン plan を閉じる (2026-10-10)

共通ルールは [`2026-10-10-handoff-common.md`](2026-10-10-handoff-common.md)。対象は `plans/active/time_integration-line-viscous-jacobian.md` (このセッションが担当してよい)。

## 状況

- 上位の諮問 `notes/reviews/2026-10-09-line-viscous-faceh-diagnose.md` の推奨は、§6.15 を限定つきの支持として記録し、値 2・3 の探索を終えて、既定 (値 0) を維持すること。
  - plan を閉じる前に次を書くよう求めている: 各試験の不合格・判定不能・未実施の処置、診断コードの残置の範囲、再開の条件、U0 と FP32 の扱い (適用範囲の制限と移管先)、§5.1 の古い要約 (「実装経路は正しい」「1 ulp」など) の訂正。
- その後の総時間の結果 (速度 plan §6.14・§6.15・§6.20):
  - 粘性入りのライン (値 3・マスク 5、`FORGE_LVC_TERMS=5`、上限なし) は、水準まで 11.5 万 step・1.18 h。
  - 値 0 + 上限 50 は 12.0 万 step・1.16 h で、判別できない。
  - 本番は値 0 のまま (環境変数に頼る経路を本番に置かない)。この結果を本 plan の結論に引用する。
- 走っている run は無い。文書の作業だけ。

## やること

1. 上の諮問の採否を §5.1・§6 に反映し、未実施・不合格の項目の処置を書く (消さない)。
2. `python3 solver_density_cuda/tools/codex_review.py plans/active/time_integration-line-viscous-jacobian.md --stage result` を回し、指摘の採否を §6.1 に書く。
3. 採否が済んだら、status を done (結論は「本線不採用・診断として残置」) にして `plans/accepted/` へ移す。`plans/README.md` の両方の表を同期する。
4. `methods/time_integration/implementation.md` の `lineViscCoupling` の記述 (値 2・3) と整合させる。

- コードは変えない (診断コードは残置)。
