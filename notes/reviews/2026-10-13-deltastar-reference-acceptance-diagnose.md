# codex 諮問 (diagnose): deltastar-reference-acceptance

- **brief**: [`notes/reviews/briefs/2026-10-13-deltastar-reference-acceptance.md`](../../notes/reviews/briefs/2026-10-13-deltastar-reference-acceptance.md)
- **plan**: [`plans/active/tooling-nozzle-deltastar-band-root-zone.md`](../../plans/active/tooling-nozzle-deltastar-band-root-zone.md)
- **date**: 2026-10-13
- **commit**: `c13be1fa` (feature/nozzle-deliverables)
- **codex**: effort `high`, 1.9 min, rc=0
- **結論**: ?

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表：

| 指摘 | 採否 | 根拠と処置 |
|---|---|---|
| M1 — **Major**：hard／held の伝播 | 採用 | [`deltastar_loop.py:85`](/home/sano/work/forge-deliv/design/forge_design/feedback/deltastar_loop.py:85) は重みをゼロにしても平滑化値を全点へ返し、同:335 の寸法入力は有限性しか検査しない。使用域の hard／held で寸法解きを停止し、抽出→平滑化→CSV→寸法の伝播試験を入れる。 |
| M2 — **Major**：負の参照流束 | 採用 | [`deltastar.py:270`](/home/sano/work/forge-deliv/design/forge_design/metrics/deltastar.py:270) は負の D を直接 `_negative_delta_r` に渡し、同:126 は非正の参照を `1e-30` に置換する。符号分岐前に有限性・正値性・探索域を検査する。正常な負の δ_r は維持する。 |
| M3 — **Major**：端点の片側許容 | 採用 | レビュー記録の `run_0009_ns_c2_final/` は出口超過 `2.70614e-6 r_t`。出典照合を維持して両側 `|Δ| ≤ ε_f32` とし、負座標の入口も試験する。 |
| M4 — **Major**：入力 δ_E の準定常 | 採用・現状未達 | [`deltaE_series.csv`](/home/sano/work/forge-deliv/case/70.m5_air_wt_demo/_deltaE_series_run_0010/deltaE_series.csv:2) の x=48 は40k→60kで `2.1960e-4 r_t` 単調減少。x=3 の報告 VERDICT も `DRIFTING`。末尾2枚の平均では受入れない。 |
| M5 — **Major**：出口 M の許容 | 採用 | [`check_quasisteady.py:543`](/home/sano/work/forge-deliv/solver_density_cuda/tools/check_quasisteady.py:543) の既定許容は目標精度より大きい。末尾変化・推定残変化とも `≤2.5e-4`、単調なら漸近値も合否範囲内を要求する。推定不能は合格にしない。 |
| m6 — **Minor**：同じ IC という記述 | 採用 | [`plan:102`](/home/sano/work/forge-deliv/plans/active/tooling-nozzle-deltastar-band-root-zone.md:102) の新 NS は `run_0010` を初期値にする。「各設計の受入れ済み定常解どうし」に修正する。 |
| m7 — **Minor**：実装手順・文書 | 採用 | [`plan:111`](/home/sano/work/forge-deliv/plans/active/tooling-nozzle-deltastar-band-root-zone.md:111) は最終 NS まで含まない。`methods` の先行更新に加え、カードから寸法照合までの試験・最終 NS・result レビューを手順と影響範囲に含める。 |

**結論: 最終 NS は保留し、既存の NS 5枚を使って「帯位置を毎回再選択／50k の帯位置に固定」だけを変える抽出 A/B を先に行う。**

**第 1 仮説:** 50k→55k の局所的な δ_E の段差は、`edge` の選択分岐の切替が主成分である。 **確度: 中**

根拠: [`deltastar.py:327`](/home/sano/work/forge-deliv/design/forge_design/metrics/deltastar.py:327) は正値だけを選別し、同:336 は最初に閾値を満たす交点を選ぶ。交点位置は補間するため、単純な「離散帯への丸め」ではないが、有効点・交点の枝は切り替わり得る。さらに同:366 の平滑化で影響が周囲へ広がる。`run_0010_ns_c2b_final/` の x=3 は50k→55kで `+1.7137e-4 r_t`、ブリーフの x=5 は約 `+5.1e-4 r_t`。**実際に枝が変わった記録は未確認。**

反証条件: 帯位置を固定しても、x=5 の50k→55kの差が `1e-4 r_t` を超えて残るなら、「帯の切替でこの段差を許容内まで説明できる」という仮説を棄却する。

**第 2 仮説:** 帯の切替とは別に、NS の抽出量に緩慢な変化が残っている。確度: 中。x=48 の5枚は一貫して減少しており、局所段差だけを消しても入力全体の受入れは成立しない。

**第 3 仮説:** Euler の局所場の変化が帯選択を通じて増幅されている。確度: 中。「float32 の床の雑音」という原因までは未確認。

**判別 A/B:** forge **0 step**。Euler は `case/70.m5_air_wt_demo/run_0012_euler_c2b/` に固定し、`run_0010_ns_c2b_final/` の40k・45k・50k・55k・60kを再抽出する。

- **A:** 各時刻で通常どおり `edge` を再選択。
- **B:** 50k の最終帯位置 `y_b(x)` を全時刻で固定。流束のフィット・積分・δ_E の平滑化は各時刻で再計算する。
- 見る量は、x=5 の50k→55kの差、使用域全体の最大差、`y_star`・`y_b`・有効梯子点・選択交点・hard／held の変化。
- **Aで段差が再現し、Bでは x=5 の差が `≤1e-4 r_t`** → 固定帯でも同程度の段差が残る仮説を棄却し、帯選択の寄与を支持する。
- **Bでも `>1e-4 r_t`** → 帯選択だけで説明する仮説を棄却する。これだけで流れの過渡とは断定しない。

Bは診断用であり、その補正表をそのまま設計入力に採用しない。

**やらない方がよいこと:** 使用域への限定だけで Euler を合格にすること、55k・60kだけを選んで平均すること、原因を確かめる前に CFL と緩和を同時に変えて回し直すこと。

受入れの推奨は次の一組に絞る。

- **Euler：現時点では受け入れない。** [`README.md:26`](/home/sano/work/forge-deliv/case/70.m5_air_wt_demo/README.md:26) の記録は A が `NOT CONVERGED`、B が `NOT CONVERGED (stalled/plateau)`。既登録の「収束 PASS＋対象量 STEADY＋延長感度の閾値」を維持する。同一実効設定の A本段→B継続を連結した区間でも収束判定を確認し、未合格の A を `--from-floor` の収束済み参照には使わない。感度閾値は `|ΔM|≤2.5e-4`、流量相対差 `≤1e-4`、全分布 `max|Δδ_E|≤1e-4 r_t` のままとする。
- **入力 δ_E：受入れを通った、事前固定の5枚・20k step窓の点ごとの平均を使う。** 受入れ済み Euler と同じ抽出手順を固定し、5k間隔で評価する。`check_quasisteady --series-csv` は窓全体を対象にし、VERDICT を残す。加えて使用域全体で `maxₓ(maxₜδ_E−minₜδ_E)≤1e-4 r_t`、hard／held なしを要求する。単調な量は推定残変化も同じ許容内、推定不能なら保留。**現在の40k〜60k窓は不合格であり、平均で不合格を消さない。**

**呼び出し側の前提への異議:**

- **Major —「縮流部を除けば Euler を受け入れられる」は数値に反する。** 許可された A/B の JSON を直接比較すると、x=10 の δ_E は `0.0317892442510 → 0.0319364854074`、差は **`1.47241156e-4 r_t`**。使用域内でも閾値超過。対案は上記の受入れ保留。
- **Major —「x<0 の δ_r を壁に直接使わない」から「設計へ影響しない」は導けない。** [`wall_axismach.py:1721`](/home/sano/work/forge-deliv/design/forge_design/geometry/wall_axismach.py:1721) は壁の区分式の記録であり、抽出・平滑化の依存範囲を保証しない。`edge` は x方向の平滑化を行う。対案は、領域限定を将来提案する場合も、下流の補正とスロート接続への伝播を確認してから事後変更として登録すること。
- **Major —「plateau＝float32 の床」は未立証。** B単独の低下0.1桁と密度残差 `9.6e-8` は停滞の観測であり、精度限界の証明ではない。対案は同一設定区間の判定と局所抽出の感度を分けて記録すること。

**不足情報:** A本段→B継続を連結した判定区間・VERDICT、段差前後の帯選択記録、使用域全点の hard／held。今回、VERDICT は README・レビュー・ブリーフの記録を照合したもので、禁止された残差・HDF5は読まず、判定ツールも再実行していない。

ファイル変更なし。受入れ条件と A/B の提案は **plan 未反映**。呼び出し側で対象 plan §5.1 #5・#6、§6へ事前登録する。
