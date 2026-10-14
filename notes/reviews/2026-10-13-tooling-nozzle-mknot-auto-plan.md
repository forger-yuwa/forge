# codex レビュー: tooling-nozzle-mknot-auto (plan)

- **plan**: [`plans/active/tooling-nozzle-mknot-auto.md`](../../plans/active/tooling-nozzle-mknot-auto.md)
- **stage**: `plan`
- **date**: 2026-10-13
- **commit**: `19e5abc7` (feature/nozzle-cooled-wall-settling)
- **codex**: effort `high`, 5.6 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M5/m2

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**

自動選択の目的と CFD を含めないスコープは妥当ですが、M6 の根拠に誤読があり、探索・最終確定・検証条件にも修正が必要です。  
保存済み M5 データを再集計し、両解像度で選択値 `2.44` になることは確認しました。ただし、これだけではカード生成全体の合格を裏付けません。

1. **Major — M6 の「余裕最大」と「最大壁角最小」を取り違えている。**

   §1・§3の「`L_c=45` で `2.6`、`50` で `2.7` が余裕最大」は誤りです。保存 CSV を同じ `n_axis=600` で集計すると次になります。

   | `L_c` | 最大壁角を最小にする `M_K` | その余裕 | 余裕を最大にする `M_K` | その余裕 |
   |---|---:|---:|---:|---:|
   | 45 | 2.6 | 1.844668° | 3.0 | 2.018833° |
   | 50 | 2.7 | 2.535646° | 3.1 | 2.767106° |

   根拠: [sweep_axislaw_A_MK_summary.csv:91](/home/sano/work/forge-deliv/case/42.isobutane_wt/sweep_axislaw_A_MK_summary.csv:91)、[同:95](/home/sano/work/forge-deliv/case/42.isobutane_wt/sweep_axislaw_A_MK_summary.csv:95)、[同:118](/home/sano/work/forge-deliv/case/42.isobutane_wt/sweep_axislaw_A_MK_summary.csv:118)、[同:122](/home/sano/work/forge-deliv/case/42.isobutane_wt/sweep_axislaw_A_MK_summary.csv:122)。計画が引用する高解像度の `14.045° / 1.841°` も、CSV では明確に `robust_min_theta` です（[同:177](/home/sano/work/forge-deliv/case/42.isobutane_wt/sweep_axislaw_A_MK_summary.csv:177)）。

   **対案:** 余裕最大という今回の目的は維持し、「従来の代表点を再現する規則」とは区別すること。§1・§3を訂正し、M6 の `L_c=45/50` を目的関数の取り違えを検出する回帰試験に追加してください。

2. **Major — 「M6 と同じ制約」になっていない。軸則の単峰条件が欠けている。**

   旧探索は `hard_gate_pass` に加えて `axis_gates["mpp_quality_ok"]` を必須にしています（[sweep_axislaw_A_MK.py:91](/home/sano/work/forge-deliv/case/42.isobutane_wt/sweep_axislaw_A_MK.py:91)）。これは **軸方向の `M″` の符号反転数 ≤ 1** であり、`M_K` に対する目的関数の山の数とは別です。

   `KnotQuinticAxisLaw.gates()` はこの品質値を返しますが、`violations` に追加するのは `M′<0` だけです（[axis_law.py:263](/home/sano/work/forge-deliv/design/forge_design/geometry/axis_law.py:263)）。したがって、`design_chain` が正常終了しても単峰条件の合格にはなりません。§4.1・§4.2の列挙にも、この条件がありません。

   **対案:** `mpp_quality_ok == True` を候補の必須条件に追加し、`mpp_sign_changes` を記録すること。V2 に「他条件は全部合格、軸則の単峰条件だけ不合格」の棄却試験を追加してください。

3. **Major — 最良の粗格子点だけを細分化しても、多峰時の最大は選べない。**

   [plan:62](/home/sano/work/forge-deliv/plans/active/tooling-nozzle-mknot-auto.md:62) は粗格子の勝者の前後だけを細分化します。一方、次行と V2(b) は多峰時にも大域最大を選ぶとしています。この保証は成立しません。

   合成関数  
   `m(x) = 1 + max(3 exp(−((x−2)/0.2)²), 3.5 exp(−((x−3.04)/0.03)²))`  
   を全点制約合格として実行確認しました。

   - `0.1` 格子では山を二つ検出し、勝者は `x=2.0, m=4.0`。
   - 計画どおり勝者の ±0.1 を細分化しても `4.0`。
   - 全域の `0.01` 格子では `x=3.04, m=4.5`。

   これは合成試験であり、実際の MOC にこの形が発生したという主張ではありません。しかし、計画自身の V2(b) に対する反例です。

   **対案:** 検出した全局所極大の周辺を細分化して比較すること。粗探索・細探索・解像度比較で共通の格子原点を定義し、同点規則は「評価済み最大値との差 ≤ `1e-3°` の集合から最小 `M_K`」と明文化してください。保証する範囲も「評価済み候補中の最大」と限定すべきです。

4. **Major — 最終 YAML に対する解像度ゲートと `r_t` 固定点の成立条件が閉じていない。**

   §4.2で解像度を確認した後、§4.3では `r_t` を変更して再選択します。しかし、**最終 `r_t` でも N/2N の確認をやり直すか**、選択値がちょうど `0.01` 動いた場合にどちらを出力するかが未定義です。

   さらに既存カード生成器は、長さ合わせの停滞時に `n_axis_inv` を倍増し、再停滞時には元へ戻して `L_total_tol` を緩和します（[wt_card.py:1324](/home/sano/work/forge-deliv/design/forge_design/intake/wt_card.py:1324)）。場が短ければ `x_end_margin` も変更します（[同:1357](/home/sano/work/forge-deliv/design/forge_design/intake/wt_card.py:1357)）。初回の選択・解像度確認を、そのまま最終出力の根拠にはできません。

   **対案:** 最終採用する `(r_t, M_K, 実効 geometry)` を一組として確定し、その組で寸法残差・全制約・N/2N ゲートを確認すること。設定変更時は関連する評価を無効化して再評価してください。V4 には YAML の値の一致だけでなく、物理スロート径・長さの公差、最終設定でのゲート成立を加える必要があります。

5. **Major — V1 の比較元が別条件で、「同じ設計の再計算」になっていない。**

   V1 は `card_m5_air_demo_c2b.yaml` から §3 の数値を相対 `1e-9` で再現するとしています。しかし保存結果は異なります。

   | 保存結果、`M_K=2.25` | `r_t` [m] | `L_c` | 最大壁角 |
   |---|---:|---:|---:|
   | `problem_m5_air_c2b` | 0.0249735612732 | 21.6992444848 | 15.2323803065° |
   | `problem_m5_air_ns2` | 0.0250732576586 | 21.4902330499 | 15.3632487897° |
   | §3 の掃引 JSON | 記録なし | 21.4902330499 | 15.3632487897° |

   根拠: [c2b の記録:479](/home/sano/work/forge-deliv/case/70.m5_air_wt_demo/problem_m5_air_c2b.yaml.card_log.json:479)、[ns2 の記録:512](/home/sano/work/forge-deliv/case/70.m5_air_wt_demo/problem_m5_air_ns2.yaml.card_log.json:512)、[掃引 JSON:29](/home/sano/work/forge-deliv/case/70.m5_air_wt_demo/_mknot_sweep/mknot_sweep_L1200.json:29)。掃引値は `ns2` の記録に対応しています。さらに、自動選択後に `r_t` を解き直すカード試験と、`r_t` 固定の掃引は別の検証です。

   **対案:** 固定入力の再現試験とカード全体の統合試験を分けること。前者は入力 YAML・`r_t`・実効設定・環境・ハッシュを固定して未丸め JSON と比較し、後者は選択条件と寸法公差で判定してください。既定の `axismach_explore_v1` も統合試験に含めるべきです。

6. **Minor — V3 の「時刻以外はバイト同一」は比較対象の定義が不足している。**

   生成 YAML の `_card.r_t_history` には実行時間 `time_s` が入り、`card_log` にも設計時間・総時間が入ります（[wt_card.py:1442](/home/sano/work/forge-deliv/design/forge_design/intake/wt_card.py:1442)、[同:1513](/home/sano/work/forge-deliv/design/forge_design/intake/wt_card.py:1513)、[同:1563](/home/sano/work/forge-deliv/design/forge_design/intake/wt_card.py:1563)）。生成日時だけを除外しても再現しません。

   **対案:** 日時・所要時間・パスなどの非決定的フィールドを明示的に正規化した比較に変更し、数値 override 時に探索関数が呼ばれないことを別途検査してください。

7. **Minor — 必須の解像度ゲートが完了条件から漏れている。**

   [plan:123](/home/sano/work/forge-deliv/plans/active/tooling-nozzle-mknot-auto.md:123) は「V1〜V4 が合格」で完了できます。既存 JSON の再集計では V5 の `2.44 / 2.44` を確認できましたが、これは今後実装するゲートの動作試験ではありません。

   **対案:** 完了条件を V1〜V5 に修正し、選択移動 `>0.01`、高解像度側だけ制約違反、探索端に最大がある場合の棄却試験も必須にしてください。

**推奨は、余裕最大の方針を維持し、上記を修正してから実装することです。** 実装前の優先順は **①根拠・単峰制約の訂正 → ②探索保証の修正 → ③最終設定の確定手順 → ④検証入力と合格条件の修正** とします。

既存の knot・要求カードの accepted plan は部品を提供していますが、この自動選択は未実装です。campaign plan の [#15](/home/sano/work/forge-deliv/plans/active/tooling-design-problem-campaign-recipe.md:189) は Euler/C2 を含む探索なので、今回の設計計算だけの選択器は独立した小機能として実装できます。CUDA、node/cell、block-DPLUR、周期境界の処理変更はなく、それらの CFD 回帰を今回の必須条件にする必要はありません。

今回の確認は保存データの再集計・コード照合・合成反例の実行です。CFD の収束・定常性は判定しておらず、§3の「CFD は鈍感」を承認根拠には採用していません。ファイル変更なし、**plan 未反映**です。

指摘数: Critical 0 / Major 5 / Minor 2
