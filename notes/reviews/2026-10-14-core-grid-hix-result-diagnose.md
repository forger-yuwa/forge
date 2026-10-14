# codex 諮問 (diagnose): core-grid-hix-result

- **brief**: [`notes/reviews/briefs/2026-10-14-core-grid-hix-result.md`](../../notes/reviews/briefs/2026-10-14-core-grid-hix-result.md)
- **plan**: [`plans/active/tooling-nozzle-core-grid.md`](../../plans/active/tooling-nozzle-core-grid.md)
- **date**: 2026-10-14
- **commit**: `4308df75` (feature/nozzle-core-grid)
- **codex**: effort `high`, 2.0 min, rc=0
- **結論**: **入口マッハ数の下限は今回設けず、現形状を保つ既存の `x_respace` 案の設計チェーンへの組み込み・受入確認を優先する。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論: **入口マッハ数の下限は今回設けず、現形状を保つ既存の `x_respace` 案の設計チェーンへの組み込み・受入確認を優先する。**

採否表:

| 重大度 | 判断対象 | 採否・根拠・対案 |
|---|---|---|
| Major | 「入口マッハ数を上げると B0 の停滞が消える」 | **却下**。H でも残差は横ばいになり、継続区間の判定は `NOT CONVERGED`。形状・格子・境界層の履歴も同時に変わっている。**「今回の上流形状変更で停滞水準が大幅に低下した」**とする。[plan:1418](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:1418)、[plan:1369](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:1369) |
| Major | §4.101 の結果 A と登録した解釈 | **採用**。判定 JSON の資格4項目はすべて真。縮流部の比と二量の `STEADY` は登録した読みを支持する。ただし排除できるのは「改善が θ_r・Q_w の過渡中だけだった」であり、全場の定常化ではない。[plan:1415](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:1415)、[plan:1417](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:1417) |
| Major | θ_r +6.0 %・Q_w +7.7 % を形状変更の確定した影響量とする | **却下**。C は未収束なので、差には形状の効果と解き残しが混在する。同設定の再現差より大きくても、その分離にはならない。**指定窓の観測差として残し、設計効果・誤差上限・許容差には使わない。**[plan:1419](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:1419) |
| Major | 入口 M ≈ 0.0565 を設計上または計算上の下限にする | **却下**。これは改善した一条件であり、必要条件でも境界値でもない。形状を変えない Grefc でも大きな改善がある。まず格子による対処を選ぶ。[plan:759](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:759) |
| Major | `x_respace` の有効性から生産採用まで済んだと扱う | **要再検証**。対処の優先順位としては採用するが、Grefc の既存試験は物理量の参照値への到達判定を満たしていない。G0c の別セッションの報告も、今回の許可範囲では判定資料を直接確認できていない。**既存の受入証拠を照合し、不足分だけ検証する。**[plan:760](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:760)、同 plan §5.1 #13(dd) |

ユーザに言える範囲は次の4文。

> 今回の入口半径半減を含む形状変更では、H の継続区間でも縮流部の残差が C の固定窓に対し、ρで約1/50、ρEで約1/140〜1/170に低下しました。H の θ_r・Q_w は `STEADY` となり、この改善が二量の過渡中だけの現象だったという説明は退けられます。一方、H の残差にも横ばいが残り、マッハ数だけの効果は分離できていないため、「停滞が消えた」「入口マッハ数に下限が必要」とは言えません。θ_r +6.0 %・Q_w +7.7 % は未収束の C との窓平均の差であり、確定した設計効果ではありません。

根拠となる区間と判定:

- H: `case/45.isobutane_m6_d155/run_0916_hi_h_cont/`、本段通算60k〜100k。保存済み `check_convergence` 出力は **`NOT CONVERGED (stalled/plateau)`**。本段＋継続は plan 記載で **`PASS`**。両方を併記し、都合のよい区間だけを選ばない。
- 同区間の `check_quasisteady`: θ_r **`STEADY`**、単調増加・推定漸近値 **0.121686**（最終比 +0.041 %）。Q_w **`STEADY`**、単調減少・推定漸近値 **9.96603×10⁶**（−0.015 %）。
- C: `case/45.isobutane_m6_d155/run_0912_hi_c_main/`、比較窓40k〜60k。本段の判定は **`NOT CONVERGED`**。
- 判定資料: [hix_judge.json](/home/sano/work/forge-coregrid/case/45.isobutane_m6_d155/_band_ab/core_grid/hix_judge.json)。今回は指示どおり生の計算出力を読まず、保存された判定を確認した。

第 1 仮説: **現形状のまま縮流部の x の間隔を変更することで、入口マッハ数の下限を導入せずに、実用上必要な残差低減と物理量の精度を得られる。** 確度: **中**

  根拠: `case/45.isobutane_m6_d155/run_0625_cg_grefc/` は、対照 `run_0619_cg_gref_b0_ctl/` に対し、縮流部のρの残差比0.027〜0.028、ρEの比0.016。形状変更は改善の必須条件ではない。ただし両 run の `check_convergence` は `NOT CONVERGED`、物理量の到達判定も未達なので、実用精度までの主張は仮説に留める。[plan:755](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:755)

  反証条件: 実装・測定の資格を満たした `x_respace` の試験で、固定予算内の残差改善条件を満たさない、または準定常になった物理量が適格な同形状参照から許容差を超えて外れる。

第 2 仮説: H の改善には局所マッハ数の上昇も寄与するが、その大きさは格子・境界層の変化から分離できていない。確度: **低〜中、寄与は未確認**。

第 3 仮説: 設けない。既存の候補を増やす証拠がない。

判別 A/B: **追加 run が受入確認に必要な場合だけ、現形状・B0 を固定して `mesh.x_respace` の有無を比較する。** 半径の追加系列は回さない。これは低マッハ機構の解明ではなく、選んだ対処の受入試験とする。

- **構成・長さ:** 同じ物理壁・BC・バイナリ・数値設定・共通ドナーを使い、本段各60,000 step、1,000 step ごとに出力。既存の対照がこの条件を満たせば流用する。W1=(20k,40k]、W2=(40k,60k] を固定する。
- **主指標:** §4.36 と同じ縮流部の体積重み付き残差 q。結果Aは両窓のρ・ρEで「有／無」≤1/3、有効側の W2/W1≤1.1、他成分・領域外・`throatwall` の監視比≤1.5。結果Bは両窓のρ・ρEで比≥2/3、かつ有効側の W2/W1≥0.8。中間は中間のまま残す。これらは**診断用の提案閾値**である。
- **資格:** 意図した局所粗化、物理壁・範囲外座標・第一層条件の保持、入力と実効設定、移送、メッシュ品質、出力数、有限性・正値を確認する。**残差の低下継続を資格違反にしない。**
- **物理量は別ゲート:** `check_quasisteady` と全保存量の `check_convergence` を別に記録する。参照との精度判定は双方の定常性を確認してから、§4.36の既存予算（不確かさ込み相対0.1 %、軸・η=0.1の M は0.005 %pt）を使う。未定常なら残差の結果Aを取り消さず、**生産採用は未判定**とする。

→ **Aかつ物理量の受入条件成立なら、現条件で形状変更・入口 M 下限が必要という説を退け、格子対処を採る。Bなら、この `x_respace` による固定予算内の大幅改善という仮説を退ける。Bでもマッハ数下限の必要性が証明されたことにはならない。**

既存の G0c 検証がこれを満たしていれば、新しい A/B は重複させない。

やらない方がよいこと: M=0.0565 の下限化、半径0.75倍などの系列追加、H の単純延長による原因確定、入口切断を同等の問題とみなすこと。`inlet_Pressure` は速度分布を直接入力できず、切断後の状態再現には別の確認が必要。[inlet-profile.md:16](/home/sano/work/forge-coregrid/procedures/inlet-profile.md:16)

呼び出し側の前提への異議: **「顕著な差が出たので下限の決め方へ進める」という推論は受け入れない。** 差の存在は確認できたが、マッハ数の因果効果と下限の必要性は確認できていない。また「形を変えない格子対処」も離散化誤差を変えるため、精度確認は必要である。

不足情報: G0c／T2 の根拠 run、実効設定、判定区間、収束・準定常の VERDICT、0.1 % 比較の対象量と不確かさ。この不足は下限導入を見送る判断を妨げないが、生産採用済みと認定するには必要である。

**plan 未反映**。読み取り専用の諮問として出力した。採用する場合の反映先は `plans/active/tooling-nozzle-core-grid.md` §4.101 の後続節と §5.1 #13。
