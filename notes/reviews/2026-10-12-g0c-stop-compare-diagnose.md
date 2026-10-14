# codex 諮問 (diagnose): g0c-stop-compare

- **brief**: [`notes/reviews/briefs/2026-10-12-g0c-stop-compare.md`](../../notes/reviews/briefs/2026-10-12-g0c-stop-compare.md)
- **plan**: [`plans/active/architecture-float-state-double-geometry.md`](../../plans/active/architecture-float-state-double-geometry.md)
- **date**: 2026-10-12
- **commit**: `894e5a43` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 3.8 min, rc=0
- **結論**: **現案の採用判定は保留し、G0 と G0c を同じ FP64・`lineDtDirectionalCap: 1` で継続する、格子だけの A/B に組み直す。** G0 の参照更新もこの比較に含める。
- **extra**: `case/45.isobutane_m6_d155/cold_xcheck.py`, `case/45.isobutane_m6_d155/t2_judge.py`, `case/45.isobutane_m6_d155/_band_ab/cold_pair/t2_judge.json`, `plans/active/axisymmetric-freestream-hoop-gauge.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論: **現案の採用判定は保留し、G0 と G0c を同じ FP64・`lineDtDirectionalCap: 1` で継続する、格子だけの A/B に組み直す。** G0 の参照更新もこの比較に含める。

採否表（Critical 0 / Major 5 / Minor 1）:

| 重大度 | 対象・採否 | 根拠と対案 |
|---|---|---|
| Major | **未収束の G0 を固定参照として採用判定する案：却下** | `case/45.isobutane_m6_d155/run_0752_t2_d1`・`run_0754_t2_d2` は、保存された判定でともに `OVERALL: CHECK FAILURES ABOVE`。再開後 122,500〜142,500 step の θ_r 全断面・Q_w は `DRIFTING`。θ_r の変化は +0.038〜0.045 %/2 万 step。再実行差の小ささは残る反復誤差を制限しない。[t2_judge.json:21](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/cold_pair/t2_judge.json:21)。**参照も継続し、両格子の収束・準定常を比較の前提にする。** |
| Major | **G0c・上限 50 対 G0・上限 1 を格子の試験とする案：却下** | 上限は平均流と SST の両方の局所時間刻みに効く。[solver-settings.md:75](/home/sano/work/forge-integ-1005/procedures/solver-settings.md:75)。未収束の比較では格子変更と反復経路の変更が混ざる。**まず両腕を上限 1 に固定する。** 上限 50 を含む設定一式の性能比較は、この格子試験とは別の問い。 |
| Major | **粗くした区間だけ q_w の許容を 3 % にする案：却下** | 根拠が「格子間隔が 2 倍」に留まる。[ブリーフ:24](/home/sano/work/forge-integ-1005/notes/reviews/briefs/2026-10-12-g0c-stop-compare.md:24)。解像度変更の倍率は、許してよい熱負荷誤差を決めない。**現案を修正するなら全区間 1 % を暫定基準とし、3 % への緩和には設計側の許容根拠を求める。** θ_r・δ_loc・Q_w の 0.1 % も今回の運用上の許容であり、絶対精度の保証ではない。 |
| Major | **`reduce_fields` の出力をそのまま物理的な熱負荷・壁応力と扱う案：要再検証** | `CP=1360`・`PR=0.72` と Sutherland 則が固定され、q_w はその λ で計算される。また τ_w は絶対値で、逆向きのせん断を区別できない。[cold_xcheck.py:31](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/cold_xcheck.py:31)、[同:80](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/cold_xcheck.py:80)。**実効の物性設定との整合を確認し、不一致なら共通の物性評価で再抽出する。** 比較座標・補間法・ゼロ近傍の絶対許容も事前に固定する。 |
| Major | **θ_r・Q_w の準定常だけで全指標を採用する案：却下** | δ_loc・局所 q_w・出口分布の定常性は別途必要。抽出器は入口出口の質量・エンタルピー流量も返すが、案では収支が採否から抜けている。[cold_xcheck.py:94](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/cold_xcheck.py:94)。**これらの時系列、壁応力・局所 y₁⁺、出口分布の局所差を加える。** 軸を除くコア評価は採用できるが、既知の軸の飛びは別指標として残す。[hoop plan:251](/home/sano/work/forge-integ-1005/plans/active/axisymmetric-freestream-hoop-gauge.md:251)。 |
| Minor | **数値収集は並行、速度は占有で測る案：採用** | 既存の採否と整合する。[float plan:1459](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:1459)。ただし速度比較には両腕の同じ到達条件・起点の扱いが必要。共有 GPU の経過時間や異なる停止条件の step 数で加速率を出さない。 |

第 1 仮説: **現案の 0.1 % 判定は残る過渡に敏感で、参照を継続すると採否が変わり得る。** 確度: 中  
根拠: 上記 T2 の θ_r は判定許容の約 4 割に相当する変化を 2 万 step で残している。さらに停止規則は θ_r のドリフトと流量欠損で決まり、Q_w・δ_loc・局所熱流束の準定常は条件になっていない。[t2_judge.py:21](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/t2_judge.py:21)。これは残りの誤差が 0.04 % という意味ではなく、その上限が分からないという問題。  
反証条件: 両腕の収束・準定常を確認した比較でも、参照更新による変化と時間窓による差が事前の反復誤差予算内に収まり、採否が変わらないこと。

第 2 仮説: **粗格子化による局所的な熱流束・せん断の変化が、下流の積分量では隠れる。** 確度: 低・未確認。G0c の定常な壁分布比較が提示されていないため、現時点では除外できない。

判別 A/B: **変更する条件は G0 / G0c の格子だけ。**

- 両腕とも同じ FP64 バイナリ・実効設定・境界条件・上限 1。`run_0754_t2_d2/res_142500.h5` を共通の診断起点とし、G0 は同一格子 restart、G0c は cross-mesh 移植を使う。これは**未収束場からの診断**であり、収束済み参照とは呼ばない。移植後・初期化後の保存量と物性を照合する。
- **各 60,000 step 固定、2,500 step ごとに保存**。最初の 20,000 step を移植・再開の応答として分け、20,000〜40,000 と 40,000〜60,000 の二窓を登録する。既存の停止値を通過しても止めない。
- θ_r・δ_loc・Q_w、壁の q_w と符号付きせん断、出口 M 分布、質量・エネルギー収支を見る。帯の外縁・出口断面・平均の重みを両格子で共通に固定する。出口 M はコア平均だけでなく分布の最大差も見る。
- **比較前提**は両腕の全保存量の `check_convergence` 合格と、採否に使う各量の `check_quasisteady` 合格。0.1 % の比較量には、各窓の drift/osc を 0.01 % 以下とする暫定予算を提案する。ただし、窓内ドリフトを残りの反復誤差の上限と同一視しない。単調なら漸近値の推定も併記する。
- **結果 A**：前提を満たし、両窓とも全指標が許容内 → この条件で「粗格子化により許容超過が起きる」を支持しない。G0c・上限 50 の速度・安定性まで合格にはしない。
- **結果 B**：前提を満たし、同じ指標が両窓で許容超過 → この条件での格子変更を不採用。
- 収束・準定常が未達、窓で判定が逆転、必要データが欠ける場合は**判別不能**。60,000 step は診断予算であり、収束を保証する長さではない。

やらない方がよいこと: G0c の旧 Q32 起点から同じ 20 万 step 試験を二重に回し、動いている T2 停止値との近さで採用すること。`DRIFTING` を限界事項として添えるだけでは、採用判定の欠けた前提は補えない。既存 T2 の結果 B や、その停止規則も書き換えない。

呼び出し側の前提への異議: **「7 桁の再現性」「登録停止規則への到達」「収束解の格子非依存性」は別である。** §6.11 自身も停止規則は定常解の同等性を保証しないと明記している。[float plan:492](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:492)。Gref/Grefc の近い数値も、両者の判定区間と VERDICT が揃うまでは補助観測に留める。AR 緩和の判断範囲は、この G0c・この条件に限定する。

不足情報: 対象四 run の実ディレクトリ・残差 CSV・HDF5 はこのワークスペースに存在しなかった。今回は保存された判定 JSON とコードを確認し、FP64 再実行差を再計算したが、判定ツールの再実行や場の独立検査はできていない。実効の輸送物性、局所 q_w の許容根拠、出口評価域の厳密な定義も必要。**ファイル変更なし・plan 未反映。反映先は呼び出し側の float plan §5.1 #25 と試験 G の事前登録節。**
