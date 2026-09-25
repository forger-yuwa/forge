# codex レビュー: convection-slau-wall-normal-chi (result)

- **plan**: [`plans/accepted/convection-slau-wall-normal-chi.md`](../../plans/accepted/convection-slau-wall-normal-chi.md)
- **stage**: `result` (diff base `39526328`)
- **date**: 2026-09-23
- **commit**: `254ceae3` (feature/sern-design)
- **codex**: effort `high`, 6.7 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M5/m2
- **focus**: 前回 (2026-09-23) の result 段 NO-GO (C0/M5/m2) への対応を検証してほしい。M1: 診断ツールに chi_n を入れ単体試験でカーネル照合し、run_0439 (ALL STEADY) の最後 3 dump で V1-b/c を再判定 (観測/予測 1.000)。M2: 欠けていた outflow+flag0 (run_0440) を実施し 2x2 を揃えた (排出の軌跡が 5 桁一致・同じ step 815 で NaN)。M3: V0 単体・衝撃足 (事前固定ノルムで L2 0.842 %)・cfl 0.4 固定点 (累積 CFL を揃えて 0.01 %)・格子感度 2 水準を実施。case/16・周期軸対称・格子感度 3 水準目は未実施と明記。M4: 4.40 桁に訂正しノイズ床の断定を削除。M5: 剥離縁限定を撤回し dp 比例に訂正 (単体試験で実証)。特に見てほしいのは (1) 撤回と再判定の記述が履歴として適切か (取り消し線で残している) (2) 格子感度 2 水準での結論「定性的に格子によらない・細かいほど早い」が 2 点から言い過ぎでないか (3) 未実施 3 項目を残したまま in_progress で止める判断の妥当性 (done にはしない) (4) C_L 0.433 % と衝撃足 0.842 % を既定化の条件に回した扱い。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**  
`chi_n` の実装は設計と整合し、診断修正と追加実験には前進があります。ただし、受入条件を満たさない比較に PASS が残っています。  
**既定 0・`in_progress` で留める判断は妥当です。現時点で `accepted` へ移すことは支持しません。**

指定 diff を HEAD `5b7bc6d0` まで確認しました。質量流束だけの変更、再構成・補正後の面速度、壁 mask、境界半割面の除外、config の拒否条件は整合しています。追加部分に明白な符号・単位・ゼロ割ガードの誤りは見つかりませんでした。周期・軸対称は未検証です。

独立に再現できた結果は次のとおりです。

- `test_diag_wall_cv_budget.py`: **`ALL PASS`**。
- case/48 `_v2/v2_{old,new}_{a,b,c}`: **`VERDICT: PASS`**。ノイズ床比 0.86–1.36 を再現。
- case/48 `_v3/v3_flag{0,1}`: 両側 **`NOT CONVERGED (stalled/plateau)`**。一方、各11 dump・5地点の `Cf`・`qw` は、既存評価関数と `check_quasisteady.py` の判定関数で指定閾値の **`STEADY`**。全 dump の `VALUE/*` に非有限値なし。

case/46 の `run_0439`・`run_0440`・`run_0442–0444` と `_v3sern` の元 CSV・VERDICT はローカルで確認できませんでした。以下では、**コードで確認したことと、plan 記載値から判断できることを区別**しています。

1. **Major — SERN の準定常ゲートが未完了で、受入条件が「既定化条件」へ移っています。**

   **根拠:** [plan:248](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi.md:248) は `C_T` と衝撃足に厳しい準定常条件を要求しています。しかし [plan:521](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi.md:521) の **0.842 % は最終壁 dump 1枚の比較**です。力係数についても、各側3本の反復はありますが、指定閾値での時系列 VERDICT がありません。反復間の再現性は、共通して残る時間方向の偏りを除去しません。

   また、[§5.1:201](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi.md:201) は case/16 を「既定化の必須条件」としていますが、[V3:241](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi.md:241) は明確に**現在の受入ゲート**です。残作業は未実施3項目だけではありません。

   **対案:** `C_T`・衝撃足壁圧・衝撃足位置の両側時系列を指定閾値で判定し、それまでは「終端差は許容内、受入判定は保留」とする。case/16、格子感度、周期・軸対称に加え、この準定常確認と累積 CFL あたりの収束比較を §5.1 に戻してください。

   **`C_L` 0.433 % は、事前の `C_T` ゲート違反ではありません。** 独立した精度評価を既定化の条件にするのは妥当です。ただし、flag 1 を設計最適化へ使う際にも必要です。衝撃足 0.842 % は既に現在の受入対象なので、その準定常確認まで将来へ送れません。

2. **Major — 累積 CFL を揃えただけでは、固定点の一致を証明できません。**

   **根拠:** [plan:531](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi.md:531) は `run_0439` と `run_0442` の比較を PASS としていますが、`run_0442` の収束・準定常 VERDICT が提示されていません。`run_0439` も記載上は `NOT CONVERGED` で、密度系列に許した準定常閾値は既定の drift 5 %／osc 10 % です。これは **0.01 % 規模の固定点比較を支える精度ではありません**。

   同じ累積 CFL で近い場になることは、同じ過渡を追っている場合にも起こります。また表の密度 L2 差 **0.011 %** は、本文の「0.01 %以内」でもありません。

   **対案:** 現結果を「同じ累積 CFL での終端場差」と記録し、固定点 PASS を保留する。両 CFL で比較精度より十分小さい時間変化を確認し、全保存量の収束判定と、独立に延長した末尾区間の比較で固定点を判定してください。

3. **Major — 2水準の起動試験から、格子非依存性と発現速度の一般則を導いています。**

   **根拠:** [plan:554](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi.md:554) で確認したのは、110万／245万節点における「flag 0 は NaN、flag 1 は所定 step まで完走」です。粗い側の flag 1 は6000 stepまでで、定常解・壁圧・収束率の格子比較はありません。これは**試した2格子で起動改善が再現した証拠**であり、格子感度の懸念が解消した証拠ではありません。

   [plan:573](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi.md:573) の体積比1.98と NaN step 比3.72についても、局所 `Δt` の実測がありません。更新には `Σṁ × Δt/V` が関わり、格子変更で面流束・時間刻み・初期補間誤差も変わります。体積だけでは説明できません。

   **対案:** 結論を「試した2格子では同じ成否。採用格子の方が少ない反復で発散」に限定する。3水準目と定常量・収束率の比較を残し、発現速度の機序を述べるなら、対応する局所 CV の面収支・`Δt/V`・初期状態を測定してください。

4. **Major — 診断の `chi_n` 漏れは直りましたが、V0・V1-b/c の全面合格を支える記録が不足しています。**

   **根拠:** [単体試験:23](/home/sano/work/forge-sern-design/case/46.sern_design/cad/test_diag_wall_cv_budget.py:23) の `kernel_mdot` は、CUDA カーネルを実行せず式を Python で書き直したものです。代数確認として有用ですが、float32、実際の壁 mask、非対象面のビット不変性は試していません。それでも [plan:589](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi.md:589) は「非対象面」まで試験済みにしています。

   [診断:104](/home/sano/work/forge-sern-design/case/46.sern_design/cad/diag_wall_cv_budget.py:104) の `roMin`・`pMin` は床到達数の判定用で、§6で要求した前処理の再現ではありません。[同:197](/home/sano/work/forge-sern-design/case/46.sern_design/cad/diag_wall_cv_budget.py:197) の平衡求根も、通常の1回の実行では最後の dump だけを処理します。

   [結果表:418](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi.md:418) には起点直後の正味流入、`|Σṁ|Δt/V/ρ` の数値、3 dump 分の平衡予測が揃っていません。153797 の1点も欠落しています。最大残差 `3.94e−13` を起点 `1.49e−9` と比べた縮小は **3.58桁**で、表全体を「5桁縮小」とは言えません。

   **対案:** 「Python による式の確認」と「実カーネルの確認」を分け、実効設定・前処理を含む面流束照合を追加する。前処理が不活性なら、その証拠を示す。起点直後の収支と、3 CV × 3 dump の予測・観測・正規化収支・実行コマンドを保存して判定してください。現段階では、観測／予測1.000を独立に追認できません。

5. **Major — M5 の撤回が本文に反映しきれておらず、誤った機序説明が現役の結論として残っています。**

   **根拠:** [plan:286](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi.md:286) は依然として「剥離縁だけを面数で確認」、同295行以降は付着境界層では `χ_n ≈ χ`、[同:507](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi.md:507) は「剥離縁のみ」「元々 `χ=0`」と述べています。いずれも壁 mask と面数比からは導けません。3D の面数比0.53 %を2Dの力係数変化の説明に転用する点も未修正です。

   一方、[V1-c の旧判定:401](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi.md:401) を取り消し線で残し、理由と再判定を続けた形式は適切です。M5にも同じ処置が必要です。

   **対案:** 上記旧説明を明示的に撤回し、現在の説明を「全壁隣接内部面が対象。凍結した状態での流束差は `Δχ × Δp` に比例」に統一する。影響位置を主張するなら、各ケース自身の `Δχ`・圧力差・流束差の分布で示してください。

6. **Minor — 2×2 の追加は有効ですが、「5桁一致」「独立と確定」は強すぎます。**

   **根拠:** [plan:370](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi.md:370) の step 800 は `2.6535e−5` 対 `2.6542e−5`、相対差 **0.0264 %**です。5有効桁の一致ではありません。両者が step 815 で発散する記録は「出口変更だけでは初期の壁排出を防げない」を強く支持しますが、長時間の相互作用まで否定しません。

   また、撤回するとした step 800 対1000の比較が [plan:385](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi.md:385) に残っています。両側とも新 BC の試験から「出口既定変更の回帰も実証」とする [同:513](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi.md:513) も、変更前後の比較にはなっていません。

   **対案:** 「本条件の初期排出に対する BC の影響は小さく、両変更が必要だった」に限定する。旧比較を撤回し、出口既定変更の回帰は flag を固定した別判定として扱ってください。

7. **Minor — 前回指摘した台帳・数値・検証範囲の不整合が残っています。**

   **根拠:**

   - [plans/README.md:62](/home/sano/work/forge-sern-design/plans/README.md:62) は「本命は未確定・plan レビューが先」のまま。
   - [case/46 README:340](/home/sano/work/forge-sern-design/case/46.sern_design/README.md:340) は撤回済みの旧 V1-c PASS、同342行は `run_0439` 実行中のまま。
   - [plan:196](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi.md:196) に4.8桁が残っています。記載値から計算すると **4.399桁**です。
   - case/48 の保存済み `cooled_plate_eval.json` から再計算した最大差は **`Cf` 0.001804 %／`qw` 0.002388 %**。plan の0.000 %／0.011 %とは異なります。許容内という判断は変わりません。
   - case/48 の起点 `run_0025` 自身も再判定は **`NOT CONVERGED`**。「収束済み場・収束床」という説明は未修正です。
   - [methods:377](/home/sano/work/forge-sern-design/methods/convection/theory.md:377) に非周期・非軸対称という検証範囲がありません。`methods/index.md` の既存リンクは有効です。

   **対案:** 現在の結論・履歴・未完了を同期し、数値には評価定義と成果物の実パスを付ける。`outflow` への既定変更は別 plan と procedures に明記されていますが、差分全体に「既定0だから既存挙動への影響なし」は適用しないでください。

**推奨は、既定0の実験機能として `active / in_progress` に留めることです。** 優先順は、V0/V1の証拠を揃える → SERN の準定常と CFL 固定点を再判定する → 未実施の回帰・格子試験を完了する → 文書を同期して再レビュー、です。ファイル変更・plan 反映はしていません。

指摘数: Critical 0 / Major 5 / Minor 2
