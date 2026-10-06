# codex レビュー: tooling-nozzle-throat-monotone-r2 (result)

- **plan**: [`plans/active/tooling-nozzle-throat-monotone-r2.md`](../../plans/active/tooling-nozzle-throat-monotone-r2.md)
- **stage**: `result` (diff base `c821b71e`)
- **date**: 2026-10-07
- **commit**: `9bc4d67f` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 5.5 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M1/m3
- **focus**: 2 回目の result 段レビュー。1 回目 (notes/reviews/2026-10-07-tooling-nozzle-throat-monotone-r2-result.md) の M1・M2・m3・m4 への対応 (commit b405c19f・9bc4d67f、§6.1・§9 2026-10-07) が足りているか、accepted に移してよいかを判断してほしい
- **extra**: `case/45.isobutane_m6_d155/README.md`, `notes/reviews/2026-10-07-tooling-nozzle-throat-monotone-r2-result.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**

前回 M1 の実入口準備 A/B は対応を支持します。単調拘束の実装と形状ゲートも確認できました。
ただし、M2 の壁取り違え検査に抜けが残っています。以下を修正してから `accepted` に移してください。

1. **Major — 同じ非単調壁を A/B 両腕に使っても「採用」になる。**

   [判定器:75](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/throat_mono_practical_eval.py:75) は他腕との照合の `status == "ok"` だけを確認します。しかし、[壁照合処理:124](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/throat_mono_ab.py:124) は、壁が同一でも `ok` を返します。

   メモリ上の負例で、両腕の spline・壁座標を同じ非単調壁にし、B の設定記録だけ `[0, 1.5]` にすると、**実座標差 0、識別可能節点 0、前提不成立 0 件で「単調壁を候補形状として採用」**を再現しました。今回の実 run が取り違えられていたという指摘ではありません。

   **対案:** この A/B では識別可能節点が存在し、それらが自腕の壁に一致することを必須にする。加えて、B の保存 spline 自体に S1 を適用するか、S1 合格済み spline と照合してください。同一壁・壁取り違えの負例を最終判定まで通すテストが必要です。

2. **Minor — 手順書が、撤回済みの比較結果を再び主張している。**

   [nozzle-design-workflow.md:102](/home/sano/work/forge-integ-1005/procedures/nozzle-design-workflow.md:102) の「Euler・NS・凝縮で旧壁との差は許容幅内」は、plan §9 の訂正と不整合です。N は新壁の絶対ゲート、K は凝縮量の準定常判定であり、新旧差の許容幅を検証したものではありません。

   **対案:** 「Euler は E′ 合格、dry NS は N 合格、凝縮 NS は K 合格。新旧差は参考値」に統一してください。

3. **Minor — IC 差の向きの訂正が JSON の説明文に届いていない。**

   [throat_mono_practical_eval.py:144](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/throat_mono_practical_eval.py:144) の `limits` は依然「α−β」です。更新後の保存 JSON にも残っています。一方、数値欄は `ic_beta_minus_alpha` で、オーバーシュート η0.1 は **β−α = +0.0014334544 %pt** です。

   **対案:** 説明文も「β−α＝番号写像−最近傍」に直し、JSON を再生成してください。

4. **Minor — 将来課題の移管先が未確定のまま「移管」と記載されている。**

   [plan §5.1 #9:162](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-throat-monotone-r2.md:162) の比較条件の訂正は妥当です。ただし、移管先は「#10 の別 plan」とだけ書かれ、#10 も「plan を分ける」という未着手の記述です。前回 m4 の移管先明示は未完了です。

   **対案:** 受け皿となる plan と残作業項目を明示するか、「未移管・本 plan に将来課題として保持」と正確に記載してください。

確認できた裏付けは次のとおりです。

- 指定された `c821b71e...HEAD` の diff を確認。ソルバ本体の変更はなく、単調拘束なしの係数・ノットのビット同一と、形状ゲート **`VERDICT: PASS`** を再現しました。
- 実入口準備 A/B のコードと[保存結果](/home/sano/work/forge/case/45.isobutane_m6_d155/_band_ab/throat_mono_entry_prep_ab.json)は整合しています。B1 は config・幾何、B2 は幾何の一致を記録し、B2 の反復数・出力間隔・緩和設定の違いも明示しています。
- `case/45.isobutane_m6_d155/run_0149_ns_mono_final_ext/` の `quantities_series.csv`：窓 60000〜80000、4 量とも **`OVERALL: ALL STEADY`** を再現。
- `case/45.isobutane_m6_d155/run_0148_ns_mono_final_cond/` の `cond_series.csv`：窓 14000〜18000、4 量とも **`OVERALL: ALL STEADY`** を再現。
- 保存された残差判定は **`NOT CONVERGED (stalled/plateau)`**、メッシュ品質は **`VERDICT: PASS (AR<=5000, skew<=0.90)`**。収束したとは判断しません。[run 索引](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/README.md:21)も確認しました。
- Euler の更新前後の集計数値は不変で、最大 `(D＋2SE)/Δq = 0.417795`。ただし、Euler 原系列と新規 NS の残差 CSV・保存場はローカルにないため、原系列からの再集計・NaN 再検査は未実施です。
- 単調拘束・メッシュ設定・物理壁・判定関数のテストは合格。排除厚さテストも、ローカル NumPy の API 差をメモリ内で補って合格しました。

**推奨は、1→2→3→4 の順に修正し、限定した生産採用として閉じることです。** 判定器を実データで再評価し、§6.1・§5.1 と移動時の索引を同期してください。今回の指摘だけを理由とする CFD 再計算は不要です。ファイル変更・forge 起動は行っておらず、本レビューは **plan 未反映**です。

指摘数: Critical 0 / Major 1 / Minor 3
