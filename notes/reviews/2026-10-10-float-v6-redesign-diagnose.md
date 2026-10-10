# codex 諮問 (diagnose): float-v6-redesign

- **brief**: [`notes/reviews/briefs/2026-10-10-float-v6-redesign.md`](../../notes/reviews/briefs/2026-10-10-float-v6-redesign.md)
- **plan**: [`plans/active/architecture-float-state-double-geometry.md`](../../plans/active/architecture-float-state-double-geometry.md)
- **date**: 2026-10-10
- **commit**: `5002e423` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 3.8 min, rc=0
- **結論**: ?
- **extra**: `case/48.flat_plate_cooled_m4/README.md`, `plans/active/axisymmetric-freestream-hoop-gauge.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（Critical 0／Major 6／Minor 1）。**推奨は、判定器と実効設定を直して case/48 の精度A/Bを先行すること。V6の回帰責務を全面的に撤去する案は却下する。**

| 重大度 | 対象 | 採否・根拠・対案 |
|---|---|---|
| Major | V6をfloat対FP64中心にする | **目的は採用、従来検証の全面代替は却下。** 段①〜④で通らない化学種・受動種・凝縮の経路をV6で検証する約束が残る（[plan:304](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:304)）。新floatと新FP64の一致だけでは、両者共通の退行を検出できない。精度比較を主軸にしつつ、Sodの厳密解、Taylor–Greenの保存則・時刻歴、凝縮の成分収支、未検証経路の回帰を残す。非定常ケースを「同じ所に止まるか」で評価しない。 |
| Major | case/48でcase/45の原因を二分する | **case/48先行は採用、提示された因果解釈は却下。** 平面化以外にもTP、積分法、圧力勾配、格子・壁解像などが変わる（[brief:32](/home/sano/work/forge-integ-1005/notes/reviews/briefs/2026-10-10-float-v6-redesign.md:32)）。差が出れば「軸対称がなくても精度差が現れる」まで。差が出なくても、冷却壁SSTと他条件との相互作用は除外できない。「SST・低Re壁が原因」とは確定できない。 |
| Major | `run_0048`を生産設定として複製 | **要再検証。** コピー元は`run_0954_sglsq_s2_f1probe_gg`であり（[fg2.sh:32](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/fg2.sh:32)、[fg8.sh:37](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/fg8.sh:37)）、現行node既定は`lsq`（[設定:506](/home/sano/work/forge-integ-1005/procedures/solver-settings.md:506)）。実YAMLが手元にないため`gg`指定の残存は未確認。生産採用を問う今回は両腕で`scalarGradient: lsq`を明示し、`slauWallNormalChi`、`limiterScaled`、`implicitSolvePrecision`なども固定・記録する。定常の実効CFLは`cfl_pseudo: 2`。旧起点からの変更は新しい判定区間として扱う。 |
| Major | `cooled_plate_eval.py --json --closure --series`をそのまま判定に使う | **却下。** `series()`は系列を保存せず真偽値だけを返し、許容幅は厚さ1 %・摩擦／熱流束1.5 %（[評価器:240](/home/sano/work/forge-integ-1005/case/48.flat_plate_cooled_m4/tools/cooled_plate_eval.py:240)）。今回の許容差より15〜20倍緩い。CD・HFには`--integrals`が必要で、それも最終場だけ（[評価器:277](/home/sano/work/forge-integ-1005/case/48.flat_plate_cooled_m4/tools/cooled_plate_eval.py:277)）。各保存時刻の全判定量をCSVに出して正式ツールへ渡す処理が必要。HFは端部を除く積分、閉合の壁熱量は別の積分点列なので定義も固定する。閉合の全エンタルピーは現在の自前再構成（[評価器:194](/home/sano/work/forge-integ-1005/case/48.flat_plate_cooled_m4/tools/cooled_plate_eval.py:194)）から`VALUE/h0`へ改める。 |
| Major | 「許容差以内、かつ再実行差の10倍以下なら一致」 | **却下。** [brief:41](/home/sano/work/forge-integ-1005/notes/reviews/briefs/2026-10-10-float-v6-redesign.md:41)の式では、再実行差がゼロに近いほど合格が難しく、雑音が大きいほど通りやすい。再実行差は許容差を広げる根拠ではなく、比較の不確かさとして加える。2本ずつは初回診断には使えるが、信頼区間や再現変動の上限を保証する標本数ではない。 |
| Major | `fl1_judge.py`を信頼して§6.20を解釈する | **判定器の修正・再評価を採用。** 全期間の残差検査ではなく末尾のみ（[判定器:48](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/fl1_judge.py:48)）。総合判定は正式ツール実行前に確定する（[判定器:81](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/fl1_judge.py:81)）。読み取りだけの合成入力検証で、`Q_w=NaN`でも`reach()`が`REACH, 25000`を返し、`abs(dv)>0.1`も偽になることを確認した。全期間の必須量・時刻・残差の有限性、出力欠損、ツール実行失敗を先に検査し、不完全なら判定不能にする。既定規則で許した`NOT CONVERGED`とツール失敗は区別する。 |
| Minor | §6.19のQ_wの不合格理由 | **訂正を採用。** 掲載値ではAの移動量だけでなく、効果／再実行差の条件も落ちる。`|Δ_B−Δ_A|=1.340 %`に対し`10|Δ_A2−Δ_A|=2.250 %`（[plan:614](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:614)、[判定式:36](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/v6ab_judge.py:36)）。「判別不能」は維持し、Aの静止条件だけを緩めれば支持になるとは考えない。 |

**結論:** 呼び出し側でV6の目的・実効設定・評価器・判定規則を事前登録し直し、§6.20の完了を待たずにcase/48のfloat対FP64比較を行う。

**第1仮説:** 現行設定のcase/48にも、時間変動と再実行差では説明できない、許容差を超える精度依存差が現れる。**確度: 低（case/48では未確認）**。  
根拠: case/45では精度切替に対するθ_rの応答が記録されているが、総合判定は「判別不能」であり、case/48への一般化はまだ仮説にすぎない（[plan:604](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:604)）。  
反証条件: 下記の有効な比較で、全対象量が時間変動・再実行差込みで許容内に入る。この場合、**このcase/48構成について**仮説を棄却する。

**第2仮説:** case/45の大きな差には、収縮部・軸対称・TP・積分法・格子条件との相互作用が必要である。確度: 低。case/48で明瞭な差が出れば、「これらがなければ精度差は出ない」という必要条件説を棄却できる。ただしcase/45と同じ機構とは限らない。

**判別A/B:** 同一ソース・同一double幾何・同一初期保存量・同一実効設定で、**ビルド精度だけ**を変える。

- **A＝float×2、B＝FP64×2、各48,000 step固定。** 初回診断としてこの本数と長さを採用する。2,000 stepごとに出力し、共通の24,000〜48,000 step、13点を判定窓とする。48,000 stepで静定しなければ判別不能とし、結果を見て窓を選び直さない。
- 初期場は`run_0025_B_tw300_y3_fx05/res_48000.h5`を7保存量の共通Q32として移す。FP64側はその値を正確に広げる。同じ入力HDF5を使い、節点対応・境界所属・保存量一致・格子品質を確認する。後処理にも共通double座標と同じステーション節点列を使い、境界層端`je`の飛びを記録する。
- 主判定量はx＝0.3／0.6／0.9のθ・δ*・Cf・q_w、およびCD・HF。相対許容差τはθ・δ*に`0.0005`、Cf・q_w・CD・HFに`0.001`。これは**精度間比較の暫定許容差**であり、物理的正確さの保証ではない。
- 各runの窓を`check_quasisteady.py --series-csv … --tail 1.0`で判定する。厚さは`--drift 0.0001 --osc 0.0001`、その他は`0.0002`。全対象の`STEADY`を要求し、単調なら漸近値も併記する。
- 各精度の2本・全窓点の平均をμ₃₂、μ₆₄とし、D＝|μ₃₂−μ₆₄|/|μ₆₄|。各精度内で平均からの最大偏差をb₃₂、b₆₄とし、E＝(b₃₂＋b₆₄)/|μ₆₄|とする。**全量でD＋E≤τなら「登録窓の対象量は許容内」**。いずれかでD−E>τなら「許容外の精度差を検出」。その間は判別不能。D−E≥0.005なら「0.5 %以上の明瞭な差」と追加表示する。これは観測した4本の変動幅による判定で、統計的信頼区間とは呼ばない。
- 全残差の`check_convergence`も併記する。`NOT CONVERGED`なら、その事実を残し、**「同じ定常解に収束」「FP64を全面廃止できる」とは結論しない**。

→ **許容内なら**、このcase/48で許容外の差が残る第1仮説を棄却する。  
→ **変動幅を超える許容外差なら**、軸対称・TP・ライン陰解法・収縮部が精度差の必須条件という説を棄却する。  
→ **非定常・変動過大・欠損なら**、原因について結論を出さない。

**やらない方がよいこと:** 旧floatとの回帰をすべて削除する、case/48の結果だけでSSTまたは軸対称を犯人と決める、再実行差に比例して許容差を緩める、`qAccumulatorFP64`や残差のdouble化を真因確定前に実装すること。

**呼び出し側の前提への異議:** `run_0025`は記録上`NOT CONVERGED (stalled/plateau)`であり、「floatで収束した場」は不正確。また、既存の`SERIES VERDICT: STEADY`は今回より緩い別判定である。`run_0028/0029`はCHT構成、`run_0030/0031`は層流・別壁解像、S2は`gg/lsq`の設定差なので、それらの数値を今回の精度比較の雑音幅として流用しない。

**不足情報:** 対象run、実YAML、HDF5、判定JSONはこの作業ツリーに存在せず、実測値とVERDICTの再集計はできなかった。数値は提供記録として扱い、コード上の欠陥と合成入力での再現を分けて判断した。ファイル変更・forge起動は行っていない。**plan未反映**。呼び出し側で`architecture-float-state-double-geometry.md`の§5.1 #10・§6 V6へ反映する。
