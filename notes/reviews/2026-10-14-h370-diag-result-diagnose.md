# codex 諮問 (diagnose): h370-diag-result

- **brief**: [`notes/reviews/briefs/2026-10-14-h370-diag-result.md`](../../notes/reviews/briefs/2026-10-14-h370-diag-result.md)
- **plan**: [`plans/active/tooling-nozzle-isothermal-wall-chain.md`](../../plans/active/tooling-nozzle-isothermal-wall-chain.md)
- **date**: 2026-10-14
- **commit**: `f6e6d7df` (HEAD)
- **codex**: effort `high`, 3.4 min, rc=0
- **結論**: **Eulerの資格FAILと受け渡し保留を維持し、既存17出力に対して帯位置を共通化した参照交換A/Bを、追加CFD 0 stepで実施してください。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（Critical 0 / Major 3）

数値とVERDICTはplanの記録を引用しています。対象run・結果JSON・NPZがこの作業ツリーにないため、実測の独立再判定はできていません。

| 論点 | 採否・重大度 | 根拠と対案 |
|---|---|---|
| 感度は結果B | **条件付き採用** | 記録のD = 0.0288 mm、全17出力で0.01 mm超なら、登録どおり `B_D` です。[判定実装](/home/sano/work/forge-review-h370/case/45.isobutane_m6_d155/h370_judge.py:1374)も最大差で分類しており、合成入力で確認しました。ただし実結果の確定にはJSON・NPZの照合が必要です。 |
| 原因はスロート下流のEuler場の時間的振動 | **要再検証・Major 1** | 今回交換したのは、異なるCFLの**固定参照各1枚**です。測ったのは参照間差であり、各参照の時間変動ではありません。さらに帯位置も参照ごとに再計算します。[抽出実装](/home/sano/work/forge-review-h370/case/45.isobutane_m6_d155/h370_judge.py:1292)。対案は、帯位置を共通にした参照交換で寄与を切り分けることです。 |
| 試験部とスロート下流を分けて記述する | **採用。ただし資格の分割は却下・Major 2** | 「今回の固定参照交換では、40〜94 r_tの最大差は0.0011 mm、0.46〜6.67 r_tで0.01 mm超」と書けます。[結果記録](/home/sano/work/forge-review-h370/plans/active/tooling-nozzle-isothermal-wall-chain.md:488)。これは試験部の絶対精度やEulerの適格性を保証しません。全使用域の条件を維持し、区間別予算・補正停止・平均参照への変更は、誤差の裏付けを伴う別登録にしてください。 |
| CFL 0.125の継続、局所リミッタ調査、時間平均参照 | **次の一手としては採らない** | CFL依存は残っていますが、既存の収束VERDICTは両腕 `NOT CONVERGED` です。[記録](/home/sano/work/forge-review-h370/plans/active/tooling-nozzle-isothermal-wall-chain.md:469)。まず既存場で抽出の寄与を分けます。平均化で変動が小さくなっても、共通の偏りや残差条件の未達は解消したことになりません。 |
| 壁解像を別課題として扱う | **採用。ただし受け渡し阻害条件として残す・Major 3** | 壁では超過面積2.9 % > 許容2 %、出口群は判定不能で、記録された全体判定は **`INDETERMINATE`** です。[結果記録](/home/sano/work/forge-review-h370/plans/active/tooling-nozzle-isothermal-wall-chain.md:487)。壁側の未達と出口側の測定不足を分けて記録してください。Grefの第一層や入口縮小による速度上昇が原因とは未確定です。局所のy₁・接線壁応力・ρ_w・μ_wと境界分類を確認してから、必要な格子修正を決めます。 |

結論: **Eulerの資格FAILと受け渡し保留を維持し、既存17出力に対して帯位置を共通化した参照交換A/Bを、追加CFD 0 stepで実施してください。**

第 1 仮説: **今回の0.01 mm超過には、参照交換に伴う帯位置の移動が必要であり、帯を共通化すれば参照場の差だけでは予算を超えない。** 確度: **低・未確認**。  
根拠: [deltastar.py:483](/home/sano/work/forge-review-h370/design/forge_design/metrics/deltastar.py:483) はδ_rの帯位置依存から最初の閾値交差を選び、平滑化して帯を決めます。[同:350](/home/sano/work/forge-review-h370/design/forge_design/metrics/deltastar.py:350)と[同:382](/home/sano/work/forge-review-h370/design/forge_design/metrics/deltastar.py:382)では、帯位置が参照フィットと欠損積分域の両方に効きます。報告された帯位置差1.34 mmは候補を支持しますが、その最大位置x = 11.2はδ_E差の最大位置x = 3.93と異なり、因果の証明にはなりません。  
反証条件: 共通帯で、有効な全17出力・全使用域の最大参照間差が **0.01 mmを超えること**。

第 2 仮説: **帯を共通化してもEuler参照の局所的なρUx分布差だけで予算を超える。** 確度: **低・未確認**。その差が時間振動か、CFL間の持続的な差かは現データの要約から判別できません。

判別 A/B: **変更する入力はEuler参照だけ。帯は両腕で同じものを使います。**

- NS入力：`case/45.isobutane_m6_d155/run_0951_h370_ns_b0/` の **112,500〜152,500、2,500 step間隔の全17出力**。比較列は既存の2,938列から縮めません。
- 腕A：`case/45.isobutane_m6_d155/run_0952_h370_euler_cfl2/res_40000.h5`
- 腕B：`case/45.isobutane_m6_d155/run_0953_h370_euler_cfl05/res_40000.h5`
- 共通帯：既存NPZの**腕Aのy_b(k,x)**を、結果を見る前に固定して両腕に適用します。各NS時刻には対応する帯を使い、参照交換時の再選択だけを止めます。[`y_b_fixed`の経路](/home/sano/work/forge-review-h370/design/forge_design/metrics/deltastar.py:339)を使えます。
- 一定xへの補間、帯幅、平滑化、緩和、検査条件は維持します。腕Aが元の抽出を再現することを確認し、生のδ_rと、登録と同じ未緩和δ_Eを記録します。
- 主判定量は **D_fixed = maxₖ,ₓ |δ_E,B − δ_E,A|**。帯の移動が寄与した位置も、既存の差分との比較で残します。

→ **結果A：全入力有効、D_fixed ≤ 0.01 mm**なら、「参照場の差だけで今回の予算超過を説明できる」を退け、帯の再選択による増幅を支持します。  
→ **結果B：D_fixed > 0.01 mm**なら、第1仮説を退けます。帯の選択だけを直しても解決しません。  
→ 固定帯で欠測・非有限・抽出不成立が出た場合は**判定不能**です。

これは原因を分ける診断です。結果Aでも、本番の抽出を固定帯へ変更したり、Eulerを適格としたりはしません。

やらない方がよいこと: **時間平均による資格の代用、スロート下流だけの事後的な予算緩和、問題区間のδ*補正の無根拠な凍結、残差床の低下率からδ_E誤差を比例推定すること。** 適格参照を得た後も、登録どおり保存した全履歴から停止規則・保持・2窓条件を再評価する必要があります。

呼び出し側の前提への異議: **δ_in = NS壁 − Euler壁は幾何量です。** 同じ格子・壁を使う今回のCFL変更でδ_in自体は動きません。[定義](/home/sano/work/forge-review-h370/design/forge_design/metrics/deltastar.py:338)。動き得るのはEulerの流れ場、それを使う参照フィット、選択される帯です。また、記録上のNS収束 `PASS`（本段）と準定常 `ALL STEADY`（保持窓）は、全使用域の参照感度や壁解像の合格を意味しません。

不足情報: 対象runの保存場・残差CSV・実効設定・各VERDICT、`h370_diag_sens.json`・NPZ、壁解像の境界別詳細がありません。コードはHEAD `f6e6d7df` を確認し、指定の3ファイルは実行版 `bc72c383` から差分なしでした。

ファイル変更・forge起動なし。**plan未反映**。呼び出し側の反映先は `plans/active/tooling-nozzle-isothermal-wall-chain.md` §6の結果解釈・新規診断登録、§5.1の残作業表です。
