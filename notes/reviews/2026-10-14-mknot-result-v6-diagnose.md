# codex 諮問 (diagnose): mknot-result-v6

- **brief**: [`notes/reviews/briefs/2026-10-14-mknot-result-v6.md`](../../notes/reviews/briefs/2026-10-14-mknot-result-v6.md)
- **plan**: [`plans/active/tooling-nozzle-mknot-auto.md`](../../plans/active/tooling-nozzle-mknot-auto.md)
- **date**: 2026-10-14
- **commit**: `88675245` (feature/nozzle-cooled-wall-settling)
- **codex**: effort `high`, 2.7 min, rc=0
- **結論**: **V6 の限定的な合格を採用し、同点幅 1e-3° を維持したまま「制約→最小→同点」の順序を確定して、#3 のカード組み込みへ進む。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（Critical なし）

| 対象 | 採否・重大度 | 根拠と推奨 |
|---|---|---|
| V6 の合格と #3 への移行 | **採用** | 事前登録は同点規則込みの選択差を判定する（[plan:124](/home/sano/work/forge-deliv/plans/active/tooling-nozzle-mknot-auto.md:124)）。提示値の 2.07／2.07、距離 0.098 は合格条件を満たす。**この固定入力・N/2N での安定性**として採用し、下記を組み込み仕様に反映する。 |
| 穴1：端数の丸め | **却下・Minor** | 原点は `M_A + coarse`（[mknot_select.py:386](/home/sano/work/forge-deliv/design/forge_design/evaluate/mknot_select.py:386)）。YAML は評価済みの値を保持し、表示だけ丸める。丸めた別の点を未評価で出力しない。 |
| 穴2：処理順序 | **採用・Major** | 現状は距離制約前の同点集合で早期 return する（[mknot_select.py:326](/home/sano/work/forge-deliv/design/forge_design/evaluate/mknot_select.py:326)、[同:332](/home/sano/work/forge-deliv/design/forge_design/evaluate/mknot_select.py:332)）。**全制約・安全距離で候補を絞る → 最小 → 同点規則**に統一する。「必ず境目＋0.05を選ぶ」という文言は削除する。 |
| 穴4：2N の境目が保守側へ移動 | **採用・Major** | 両解像度の保守的境目を使うのが既定方針（[plan:63](/home/sano/work/forge-deliv/plans/active/tooling-nozzle-mknot-auto.md:63)）。共通の安全領域で各解像度の選択を独立にやり直し、選択差・制約・距離を再判定する。 |
| 穴5：入口の MOC ゲート | **採用・Major** | 選択時は記録だけ（[mknot_select.py:172](/home/sano/work/forge-deliv/design/forge_design/evaluate/mknot_select.py:172)）だが、`prepare` は必須判定する（[runner_axismach.py:334](/home/sano/work/forge-deliv/design/forge_design/evaluate/runner_axismach.py:334)）。同じ関数の判定を選択制約へ組み込む。`fixed2` は `pass=None` を許す既存仕様も保持する。 |

結論: **V6 の限定的な合格を採用し、同点幅 1e-3° を維持したまま「制約→最小→同点」の順序を確定して、#3 のカード組み込みへ進む。**

第 1 仮説: **［Major］現行の処理順序は、安全距離で除外される点の θ_max に最終選択が左右されるため、意図した制約付き最小化と一致しない。V6 ではこの欠陥が表面化していない。** 確度: 高

  根拠: `apply_rule` の早期 return は、距離制約を適用した候補集合で最小値を取り直さない。許可されたコードから `boundary_distance`・`_pick`・`apply_rule` を抽出し、成立判定などを合成入力用に差し替えた読み取り専用の実行で、以下を確認した。

  下側境目 1.98、ΔM_safe=0.05、tie=0.001° とする。

  | M_K | θ_max A [°] | θ_max B [°] | 安全距離 |
  |---|---:|---:|---|
  | 2.00 | 15.0000 | 15.0020 | 不合格 |
  | 2.04 | 15.0008 | 同左 | 合格 |
  | 2.05 | 15.0015 | 同左 | 合格 |

  **実行結果：現行コードは A→2.04、B→2.05。制約を先に適用すると両方→2.05。** 安全候補の目的値を一切変えていないのに選択が変わる。これは M5 の実測値ではなく、処理順序の反例である。

  反証条件: この合成入力で、安全距離外の点の θ_max を変更しても選択が変わらず、両方 2.05 になること。現行コードでは満たさなかった。

第 2・第 3 仮説: 追加しない。提示された V6 は、前回懸念した「この固定入力で境目への接近が選択を不安定化する」を再提示する根拠にはならない。

判別 A/B: **上表の安全距離外の点の θ_max だけを変える合成 A/B を推奨し、本レビューで実施済み。** CFD・MOC の追加計算は不要。選択が不変なら第1仮説を棄却、変われば支持する。実行結果は後者だった。修正後は両方 2.05 を要求する。

やらない方がよいこと: 素の最小の移動 0.02 を理由に、事前登録済みの V6 を不合格へ読み替えること。また、2.08 を同点集合へ入れるために、結果を見て tie を拡大すること。

呼び出し側の前提への異議: **［Minor］「θ_max の相対差が小さいから同じ設計」は未立証。** θ_max は壁形状全体・曲率・流れの質を代表しない。「目的関数について設定した許容差内」と表現すべきである。一方、2N の 2.07 は最小との差が提示値で 0.00079°となり、既定の同点条件内にある。段の存在だけでは、同点規則の失敗とも判定できない。

不足情報: V1a・V6 の元 JSON は指定に従い未読で、数値と PASS は提示情報に基づく。全試験の再実行もしていない。最終の `r_t`・実効 `geometry` での V3/V4 は今後必要であり、今回の合格からカード全体の検証完了までは主張できない。
