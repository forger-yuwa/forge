# codex 諮問 (diagnose): t2-result

- **brief**: [`notes/reviews/briefs/2026-10-12-t2-result.md`](../../notes/reviews/briefs/2026-10-12-t2-result.md)
- **plan**: [`plans/active/architecture-float-state-double-geometry.md`](../../plans/active/architecture-float-state-double-geometry.md)
- **date**: 2026-10-12
- **commit**: `de9a3c15` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 5.9 min, rc=0
- **結論**: **修正版 N3 を採り、T2 を結果 B として閉じ、double 廃止を保留したうえで、生産の停止値の設定依存性をユーザへ報告する。**
- **extra**: `case/45.isobutane_m6_d155/t2_judge.py`, `case/45.isobutane_m6_d155/_band_ab/cold_pair/t2_judge.json`, `case/45.isobutane_m6_d155/_band_ab/cold_pair/rcap_judge.json`, `plans/active/time_integration-line-implicit-speed.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（Critical なし）

| 重大度 | 対象 | 採否・根拠・対案 |
|---|---|---|
| Major | §6.56 の結果 B | **採用。** JSON から再計算すると、ρu 残差の中央値は float／FP64 が **2.715〜2.955 倍**。θ_r だけなら片方は許容内だが、**両方の float が登録要件を外している**。「上限 1 だけでは今回の運用要件を満たさなかった」と記す。「float は使えない」への一般化は却下。[判定値](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/cold_pair/t2_judge.json:172) |
| Major | Q1・N1：float だけ延長する | **提示された判別方法は却下。** FP64 自身も動いているため、float が過去の FP64 到達値を通過しても、同じ定常値に向かう証明にならない。差÷末尾ドリフトは時間尺度の換算にすぎない。遅れと持続的な偏りの識別は未了として残す。[結果記録](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:1430) |
| Major | Q2：許容境界なので判定不能か | **却下。停止時比較の判定は B、定常精度の判断は未確定。** 停止ドリフト ≤0.05 % は、残りの反復誤差 ≤0.05 % を保証しない。ただし、これは登録済み判定を覆す理由ではない。将来の許容変更は別の運用判断として記録し、T2 の B を保持する。θ_r の許容だけを緩めても残差フラグは残る。[登録基準](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:490) |
| Major | Q3：丸めの床なので運用上問題ない | **要再検証。** 比べているのは異なる場の残差であり、中央値の約 3 倍だけでは、残差評価の丸め・更新による場の偏り・残存過渡を分離できない。「末尾窓の残差高値」と記す。欠損とドリフトの合格を、独立した残差要件の代用にしない。[中央値と判定](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/cold_pair/t2_judge.json:104) |
| Major | Q4・N2：生産比較へ進む | **懸念の共有は採用、N2 の即時実施は却下。** −1.219 % は「上限 1 の FP64 停止値と旧起点の差」であり、生産設定の真の誤差ではない。N2 は精度と上限を同時に変えるため、候補一式の運用比較には使えても原因診断には使えない。既存の Gref の上限比較も踏まえ、生産側の検証方針と統合する。[起点比較](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:1435)、[既存調査](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:184) |
| Major | N3：「精度差が 20〜40 分の 1」とまとめる | **N3 を修正採用。** 旧比較は共通起点からの変化量・20k〜30k の窓、新比較は各停止時の値・FP64 停止値を分母とする。したがって、20〜40 倍を同条件での改善率にしない。「短期の大差を抑える効果は確認できたが、停止時の登録要件は未達」とまとめる。[旧指標の定義](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:1187) |
| Minor | 準定常判定の記録 | **補正を採用。** 実際の閾値は θ_r が `1e-4`、Q_w が `2e-4`。また JSON は最後の総合判定行だけを保存しており、各量の `DRIFTING` 行を保持していない。量別判定・区間を保存する。[判定コード](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/t2_judge.py:114) |

結論: **修正版 N3 を採り、T2 を結果 B として閉じ、double 廃止を保留したうえで、生産の停止値の設定依存性をユーザへ報告する。**

第 1 仮説: ρu の約 3 倍の高値には、同じ状態でも現れる残差評価経路の精度差が大きく寄与している。**確度: 低（未確認）**。  
根拠: `case/45.isobutane_m6_d155/run_0753_t2_f1`・`run_0755_t2_f2` の中央値は 2.20e−4・2.04e−4、FP64 の `run_0752_t2_d1`・`run_0754_t2_d2` は 7.43e−5・7.50e−5。ただし状態を固定した比較ではない。  
反証条件: 同じ状態で残差を評価すると精度間の差が小さく、FP64 でも float の高い残差が残る。

第 2 仮説: 高い残差は、精度依存の更新によって形成された場の反復誤差を主に反映している。未確認。  
第 3 仮説: θ_r の差は有限期間の進み方の違いで、最終的には小さくなる。未確認。**短い延長だけで漸近値の一致・不一致を確定することはできない。**

判別 A/B: **今回は追加投入せず、原因調査を再開する場合は、同じ float 停止場からの「更新前の残差評価」の精度 A/B を一つだけ行う。**

- 変えるのはビルド精度だけ。上限 1・格子・設定を固定し、各 1 step、再実行付き。残差組立時の全保存量が同じことを前提とし、初期化で変われば判定不能。
- 観測量は同一領域・同一定義の `res_roUx` の RMS と場の差。事前登録案として、再実行幅を考慮しても float／FP64 が **2 倍以上なら A**、**1.25 倍以下なら B**、中間は判定不能。
- **A なら**残差評価の精度依存性を支持する。**B なら**「評価精度だけが約 3 倍の主因」を退け、場の違い・更新履歴を調べる。どちらも θ_r の定常偏差や生産採用の判定には転用しない。

やらない方がよいこと: float だけを旧 FP64 値に近づくまで延長すること、T2 を通すために閾値を変更すること、上限 1 または 50 の停止値を真値と置くこと。共有 GPU の経過時間を採用判断用の速度差にすることも避ける。

呼び出し側の前提への異議: **「止まった値」「丸めの床」「生産の精度の問題」は、観測より強い表現である。** 保存済みの正式判定は全 4 本で次のとおり。

```text
check_convergence: OVERALL: CHECK FAILURES ABOVE
check_quasisteady（θ_r・Q_w）: OVERALL: NOT ALL STEADY (see above)
```

準定常の判定区間は再開後、FP64 が **122,500〜142,500 step**、float が **120,000〜140,000 step**。報告はこの区間と停止時の観測に限る。run と成果物の索引は [case README の T2 行](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/README.md:224)。

不足情報: 対象 run の HDF5・CSV・見張り系列・量別 VERDICT はこの環境にない。今回はコード、計画、保存済み JSON を照合して比較値を再計算した。全期間の有限性や初期化後の同一性は独立に再検証できていない。**ファイル変更なし・forge 起動なし・plan 未反映**。採否と今後の判断は、呼び出し側で float plan §6.56・§5.1 に反映する。
