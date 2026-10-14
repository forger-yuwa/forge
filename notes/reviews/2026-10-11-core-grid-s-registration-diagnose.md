# codex 諮問 (diagnose): core-grid-s-registration

- **brief**: [`notes/reviews/briefs/2026-10-11-core-grid-s-registration.md`](../../notes/reviews/briefs/2026-10-11-core-grid-s-registration.md)
- **plan**: [`plans/active/tooling-nozzle-core-grid.md`](../../plans/active/tooling-nozzle-core-grid.md)
- **date**: 2026-10-11
- **commit**: `7060402a` (feature/nozzle-core-grid)
- **codex**: effort `high`, 3.1 min, rc=0
- **結論**: **基準窓を固定したまま Gref・Grefx に各 6 万 step の point 継続を一度だけ行い、後半 4 万の資格と、s に point 側の変動幅も加えた予算で判定する案を推奨する。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（Critical なし）

| 対象 | 重大度・採否 | 根拠と対案 |
| --- | --- | --- |
| (a) 起点・基準窓の固定 | **採用** | `run_0608`・`run_0610` の `res_60000` から、同一格子用の `restart_field.py` で継続する。基準の平均・U・分母は変更しない。[plan:551](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:551) の扱いと整合する。 |
| (a) 4 万 step と point 側の資格 | **Major・修正採用** | 提案では point 側の d・U・準定常判定の窓が未定義。既存判定は「最初の 2 万を除外し、後半 4 万で資格、最後の 2 万で平均」である。[判定コード:157](/home/sano/work/forge-coregrid/case/45.isobutane_m6_d155/cg_gref_judge.py:157)。**今回は初めから各 6 万 step を上限に登録**し、この窓構成を維持する。4 万 step が本質的に不足すると断定したものではない。 |
| (b) 平均差 s だけを加える合否式 | **Major・修正採用** | s を窓平均の差に変更すると、point 側の変動幅が最終予算に入らない。既存コードの予算は基準側の U と p だけである。[判定コード:223](/home/sano/work/forge-coregrid/case/45.isobutane_m6_d155/cg_gref_judge.py:223)。**s は平均差として報告し、最終式には両格子の U_point も足す**。基準側の U を s 内に再加算しない。 |
| 判定対象を 8 量に限定する実装 | **Major・要修正** | [plan:582](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:582) は [40,50] を除外したが、コードは同窓を生成し、CSV の全量を資格・合否に使う。[判定コード:23](/home/sano/work/forge-coregrid/case/45.isobutane_m6_d155/cg_gref_judge.py:23)、[181](/home/sano/work/forge-coregrid/case/45.isobutane_m6_d155/cg_gref_judge.py:181)。**必須の 8 量を明示し、欠測なら判定不能、[40,50] は記録専用**にする。今回の結果 A が誤りだという指摘ではない。 |
| (c) 資格不足時の延長 | **採用：上限で止める** | 登録した 6 万 step で資格がなければ判定不能。自動延長しない。[plan:553](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:553) と同じ扱いにする。 |
| (b) 「生産の候補」の意味 | **Minor・修正採用** | この試験で通すのは精度・反復頑健性の条件。[plan:174](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:174) の総時間条件は別に残る。28 対 41.8 ms/step だけで生産置換まで確定しない。 |

結論: **基準窓を固定したまま Gref・Grefx に各 6 万 step の point 継続を一度だけ行い、後半 4 万の資格と、s に point 側の変動幅も加えた予算で判定する案を推奨する。**

第 1 仮説: Grefx は point 継続後にも資格を満たし、反復感度と追加の変動幅を含めて許容差内に残る。確度: **低**。  
根拠: `case/45.isobutane_m6_d155/run_0608_cg_gref_cont_k1/` と `case/45.isobutane_m6_d155/run_0610_cg_grefx_ext_k1/` の比較では、θ_r [65,75] の残予算は 0.684 %、Q_w は 0.788 %ある。[plan:565](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:565)、[570](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:570)。ただし現構成の s は未測定。旧 `run_0262` の θ_r(70) −0.081 % は、格子・バイナリ・指標が異なり、今回の上限には使えない。  
反証条件: 登録した point 窓でどれかの必須量が資格を欠く、または下記の最終予算が τ を超える。その場合、「この予算・窓で候補性を確認できる」を退ける。真の格子誤差が τ を超えるという意味ではない。

第 2・第 3 仮説: 追加しない。方式固有の効果と自然な継続変化の分離は、この試験の目的に含めない。

判別 A/B: **既存の B0 窓を A、新しい point 継続を B とする比較を、両格子で一度だけ行う。変更因子は反復プリセットである。**

- 起点はそれぞれ `run_0608`・`run_0610` の `res_60000`。新しい run に保存量をビット一致で移す。
- point は `time.deltaT` 配下の `lineImplicit: 0`、`lineDtDirectional: 0`、`implicitThermalJacobian: 0`。`cfl_pseudo: 4`・`implicitRelax: 0.7`、格子・境界条件・空間離散化・`mesh.axisSegmentRWeight: 1`・バイナリは固定する。この比較から個々のスイッチの寄与は主張しない。
- 各 60,000 step、2500 ごとに出力。point の μ は **40k〜60k の台形時間平均**、d は **20k〜40k と 40k〜60k の平均差**、U_point は **20k〜60k の μ からの最大偏差**とする。
- 基準の μ は両 run の **40k〜60k**、U は既存の **20k〜60k**から固定する。相対量はすべて固定した D = |μ_Gref,base| で % に換算する。軸・η 0.1 の M は %pt の差を使う。
- point の資格は、必須 8 量について `check_quasisteady.py` の **20k〜60k、`--tail 1.0` の VERDICT: STEADY**、d ≤ τ/8、U_point ≤ τ/4、漸近値と最終値の差 ≤ τ/8。drift・osc の換算は現コードと同様に行う。必要な 17 出力・全必須量・有限値が揃わなければ判定不能。残差の判定も point の区間で取得し、既存 run の PASS を引き継がない。

量ごとに、s_G = |μ_G,point − μ_G,base| として、最終条件を次にする。

```text
B_total = |E| + U_ref,base + U_x,base + p_ref + p_x
              + s_ref + s_x + U_ref,point + U_x,point
B_total ≤ τ
```

U_point を資格だけで制限する案は採らない。例えば θ_r [65,75] で s が両方ゼロでも、両 U_point が 0.37 %なら各々 τ/4 = 0.375 %以内だが、合計予算は 0.816 + 0.37 + 0.37 = **1.556 %**になる。これは実測ではなく、提案式で取りこぼす場合の例である。

→ **結果 A**：両 point が資格を満たし、全 8 量で B_total ≤ τなら、登録条件内の運用上の頑健性を支持し、精度面の候補としてユーザに諮る。  
→ **結果 B**：資格を満たすがどれかで超過なら、今回の予算では候補化を見送る。資格不足なら別扱いの**判定不能**として止める。いずれも自動延長しない。

やらない方がよいこと: 基準窓・分母の取り直し、符号付き s 同士の相殺、合格するまでの延長、point の残差 PASS だけでの採用、G0・旧バイナリの s の流用。

呼び出し側の前提への異議: 「資格あり」は「収束済み」と同義ではない。記録上、Gref は **NOT CONVERGED**、Grefx は `--from-floor` で **PASS**であり、両者の量が STEADY という別の条件で比較している。[plan:516](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:516)、[558](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:558)。また、s は自然な継続変化も含む運用差として扱い、真の反復誤差の上限や方式単独の因果効果とは呼ばない。

不足情報: point の実効設定・初期保存量の一致・現バイナリでの費用は投入時に確認が必要。今回、生データの読み取りと判定再実行は禁止範囲のため行っておらず、数値・VERDICT は plan の記録を根拠とした。**ファイル変更なし、plan 未反映**。採用時の反映先は `plans/active/tooling-nozzle-core-grid.md` の s 登録節・§5.1 #13・§6。
