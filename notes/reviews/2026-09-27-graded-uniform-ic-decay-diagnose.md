# codex 諮問 (diagnose): graded-uniform-ic-decay

- **brief**: [`notes/reviews/briefs/2026-09-27-graded-uniform-ic-decay.md`](../../notes/reviews/briefs/2026-09-27-graded-uniform-ic-decay.md)
- **plan**: [`plans/active/axisymmetric-graded-grid-static-gas.md`](../../plans/active/axisymmetric-graded-grid-static-gas.md)
- **date**: 2026-09-27
- **commit**: `2ebc46a3` (feature/cht-axisym-fem2d)
- **codex**: effort `high`, 3.4 min, rc=0
- **結論**: **A を同一数値設定・非連成のまま累積100000 stepまで延長することを推奨するが、これは静止保持への到達試験であり、ソルバ欠陥の有無を二択で確定する試験にはしない。**
- **extra**: `case/62.conjugate_disk/README.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

| 重大度 | 判断 | 根拠と対案 |
|---|---|---|
| **Major** | 「居座る」の訂正を**採用。今行う** | `run_0017_hold_A_uniform` は残差 **`NOT CONVERGED (still converging)`**、準定常は **全6系列 `DRIFTING`**。定常的な異常状態は観測していない。記述を「20000 step では大きな起動過渡が残り、減衰中。最終到達状態は未確認」に替える。 |
| **Major** | 「延長で合格すればソルバの欠陥ではない」を**却下** | 初期場変更で振幅が小さくなることと、不自然な増幅・遅い減衰の原因が正常であることは別。[流体側 plan §5.1 #2](/home/sano/work/forge-cht/plans/active/axisymmetric-graded-grid-static-gas.md:77) の結果も後者を証明しない。合格時の結論は「この条件では反復延長で静止保持基準に到達可能」に限定する。 |
| **Major** | 「非一様では減衰が約500倍遅い」を**却下** | [plan §3](/home/sano/work/forge-cht/plans/active/axisymmetric-graded-grid-static-gas.md:35) の `0.231 / 4.3e-4 ≈ 537` は同一反復数での**振幅比**。減衰率比ではない。一様側の時系列から同じ量の対数減衰率を求めるまで、倍率を原因説明に使わない。 |
| **Minor** | B の連結残差 PASS を**採用、記録は補完** | `case/62.conjugate_disk/run_0018_hold_B_ext40k/CONVERGENCE_VERDICT.txt` は継続区間単独の `NOT CONVERGED`。一方、旧区間と連結して正式判定関数をメモリ上で再実行すると **`PASS`、8.0–8.3桁低下**を再現した。累積30000–40000の準定常も **72系列 `ALL STEADY`**を再現。連結判定の成果物を別名で残し、区間を明示する。 |

結論: **A を同一数値設定・非連成のまま累積100000 stepまで延長することを推奨するが、これは静止保持への到達試験であり、ソルバ欠陥の有無を二択で確定する試験にはしない。**

第 1 仮説: **A の20000 step時点の大振幅は、起動で励起された反復過渡が主体で、最終的には静止保持の許容範囲へ減衰する。** 確度: **中**  
  根拠: `case/62.conjugate_disk/run_0017_hold_A_uniform/static_hold_series.csv` を再集計すると、10000→20000で `max|U|` は **1.479→0.23094 m/s**、`qerr_cj` は **14079.8→3125.3 %**。全保存量の有効残差列は `falling`、準定常の VERDICT は **`DRIFTING`**。ただし、`max|U|` の対数回帰による e-fold は10000–15000で **4441 step**、15000–20000で **6892 step**、18000–20000で **7894 step**と長くなっている。一定率の指数減衰は実証されていない。  
  反証条件: 非連成・同一設定で、残差 `PASS` と対象全系列 `STEADY` を満たしながら、静止保持の閾値を超える状態に落ち着くこと。**100000 stepでまだ減衰中、だけでは反証にならない。**

第 2 仮説: **非一様格子と局所擬似時間・陰解法更新の組合せが、遅い反復モードを作っている。** 確度: **低、未確認**。局所刻みが `cfl_pseudo` で決まることは設定リファレンスにあるが、両格子の実効刻み分布と減衰モードの対応は未提示。

第 3 仮説: **境界処理または離散化に起因する、閾値を超える別の定常状態への接近。** 確度: **低、未確認**。現在の減衰だけでは排除できないが、既知の slip 欠陥と同一とは認めない。

判別 A/B: **既存の初期場 A/B の A 側だけを延長する。新しい設定変更の A/B は重ねない。**

- A は `run_0017_hold_A_uniform` の最終場から、新規 run に **+80000 step**。B は静止保持を確認した `run_0017_hold_B_conduct`＋`run_0018_hold_B_ext40k` を参照とする。変更する実験因子は反復長のみ。
- 保存間隔は100 step。既存6系列に加え、**固定節点IDごとの両壁の符号付き熱流束**を保存する。
- 最終判定窓を**累積90000–100000**に事前固定する。その窓を含む累積80000–100000の系列に、既存の `--tail 0.5 --drift 0.001 --osc 0.001` を適用。閾値は速度、相対圧力差、両壁それぞれの全節点熱流束誤差の既登録値を維持する。残差はAの旧区間と継続区間を連結して通常判定する。
- **静止保持の全条件を満たすなら**、この構成で「異常状態に居座り続ける」を棄却する。**残差 PASS・全系列 STEADY なのに閾値超過なら**、「過渡が消えれば許容状態へ入る」を棄却する。**それ以外は判別未了**で止める。NaN/Infなら早期停止する。
- `warmup` は**実際に使われる絶対 step の終了値より大きく**する。現値100000のまま累積100000まで延ばしても非連成とは限らない（[設定:44](/home/sano/work/forge-cht/case/62.conjugate_disk/run_0017_hold_A_uniform/solverConfig.yaml:44)、[plan §5.1 #3](/home/sano/work/forge-cht/plans/active/axisymmetric-graded-grid-static-gas.md:78)）。非連成条件の維持を確認する。

**100000は試験上限であり、収束予測ではない。** `qerr_cj` の末尾 e-fold約9000 stepを楽観的に固定しても、100000時点は約 **0.43 %**で、ようやく振幅閾値付近になるだけ。厳しい準定常基準まで通る保証はない。この証拠から「安い有限長の試験で必ず二択が決着する」とは設計できない。

やらない方がよいこと: 延長と同時にCFL・次数・slipを変更すること、未達を直ちに欠陥と断定すること、合格するまで判定窓を動かすこと。今回の結果だけで純伝導IC限定のCHT保証範囲を広げることもしない。

呼び出し側の前提への異議: **訂正は延長結果を待たずに行う。** [流体側 plan §5.1 #4](/home/sano/work/forge-cht/plans/active/axisymmetric-graded-grid-static-gas.md:79)、[methods/boundary.md:500](/home/sano/work/forge-cht/methods/boundary.md:500)、[発注元 plan §6:217](/home/sano/work/forge-cht/plans/accepted/boundary-cht-axisymmetric-fem2d.md:217)、case README・報告ページの「居座る」を、上記の観測記述へ訂正する。過去の失敗、再登録の履歴、純伝導IC限定の保証範囲は維持する。また「全量が単調」は全期間には成立せず、CSVはstep 100から始まり、`qerr_hot` は10000–15000にも増加区間が10回ある。**末尾で減衰している**という限定なら正しい。

不足情報: 一様格子側の同じ量の時系列、両格子の実効局所擬似時間刻み、系列抽出処理と固定節点IDの対応、報告ページの所在。これらがないため、減衰速度差の機構は未確定。**ファイルは変更しておらず、plan未反映。**
