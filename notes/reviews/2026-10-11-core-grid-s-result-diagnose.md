# codex 諮問 (diagnose): core-grid-s-result

- **brief**: [`notes/reviews/briefs/2026-10-11-core-grid-s-result.md`](../../notes/reviews/briefs/2026-10-11-core-grid-s-result.md)
- **plan**: [`plans/active/tooling-nozzle-core-grid.md`](../../plans/active/tooling-nozzle-core-grid.md)
- **date**: 2026-10-11
- **commit**: `e00e39fc` (feature/nozzle-core-grid)
- **codex**: effort `high`, 3.2 min, rc=0
- **結論**: **候補 1 を選び、旧判定を「判定不能」のまま残して、Gref の point だけを追加 6 万 step 継続する新しい固定長試験を一度登録する。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（根拠は許可された plan の記録。生データの再検証はしていない）

| 重大度 | 主張・候補 | 採否と理由・対案 |
| --- | --- | --- |
| Major | 「Gref は反復方式で 1 % 動く」 | **要再検証**。観測されたのは「point 継続中の θ_r の窓平均が B0 基準から約 1 % 動いた」。`check_quasisteady` は **DRIFTING / TRANSIENT-UNSETTLED** で、定常解の差でも方式単独の効果でもない。この限定した記述を採用する。根拠: [plan:603](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:603)、[plan:596](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:596)。 |
| Major | 「Grefx は方式に依らず、Gref より信頼できる」 | **一般化は却下**。「試した B0→point の継続では必須 8 量が STEADY、s_x ≤ 0.0015 %」は採用。反復に対する再現性の証拠は Gref より強いが、空間精度の優劣は示さない。根拠: [plan:603](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:603)、[plan:408](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:408)。 |
| Major | Grefx を今、生産候補として諮る | **却下**。参照側の資格が欠け、登録した予算が評価できない。θ_r [65,75] は s・point の変動を加える前で 0.816 %。残りは最大 **0.684 %** なので、s_ref ≈ 1.10 % のまま定常化しても通らない。資格の回復と予算内への復帰は別条件である。根拠: [plan:566](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:566)、[plan:594](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:594)、[plan:604](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:604)。 |

結論: **候補 1 を選び、旧判定を「判定不能」のまま残して、Gref の point だけを追加 6 万 step 継続する新しい固定長試験を一度登録する。**

第 1 仮説: Gref の point に残る変動は減衰途中であり、point 通算 12 万 step までに、登録した量の資格を回復できる。**確度: 低**。追加 6 万で足りるという予測は未確認であり、これは期限を設けて反証する仮説である。  
　根拠: `case/45.isobutane_m6_d155/run_0611_cg_gref_point/` の θ_r [65,75] は −1.12 % まで動いた後、−1.07 % へ戻っている。一方、Q_w はまだ下がり、rms_roOmega は 9.1e1 → 2.6e3 に増えている。減衰過渡とは整合するが、持続振動と区別できていない。[plan:602](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:602)、[plan:605](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:605)。  
　反証条件: 通算 80k〜120k の固定窓で、必須 8 量のどれかが下記の資格を満たさない。その場合、「追加 6 万で適格化する」を棄却する。過渡一般を否定したとは扱わない。

第 2 仮説: 現行 point 設定では、Gref に予算を超える持続振動または長時間のドリフトが残る。確度: 低、未確認。ω 残差の増加だけでは原因を乱流方程式に特定できない。  
第 3 仮説: なし。許可された資料から、実装の特定箇所を原因として挙げる根拠はない。

判別 A/B: **変えるのは継続長だけ（point 通算 60k → 120k）**。方式の因果切り分けではなく、参照側の資格が回復するかを調べる試験とする。

- `run_0611_cg_gref_point/res_60000.h5` から、同一メッシュへ `restart_field.py` で再開する。新しい run に追加 60,000 step、出力間隔 2500。格子・バイナリ・BC・point の設定は固定する。Grefx は再実行しない。
- 判定窓は **point 通算 80k〜120k の 17 出力**。μ は 100k〜120k の平均、d は 80k〜100k との平均差、U は 80k〜120k の μ からの最大偏差。B0 の基準窓・分母・Grefx の既存評価は固定する。
- 資格は §4.26 を維持する。必須 8 量が `check_quasisteady` **STEADY**、d ≤ τ/8、U ≤ τ/4、漸近値と最終値の差 ≤ τ/8、欠測・非有限なし。θ_r・δ_loc・Q_w ではそれぞれ **0.1875 %、0.375 %、0.1875 %**。M は登録済みの単位と換算を維持する。
- `check_convergence` は同一設定の point 通算 0〜120k を対象に取り直す。現時点の point 区間は **両格子とも NOT CONVERGED** であり、Grefx の B0 の PASS を引き継がない。[plan:602](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:602)。
- **結果 A：資格を回復** → 新しい s_ref・U_ref,point で §4.26 の B_total を再計算する。全 8 量が τ 以下なら精度・反復の頑健性の候補化へ進む。超過なら現行予算での候補化を見送る。総時間の条件は別途必要。
- **結果 B：資格を回復しない** → 第 1 仮説を棄却し、Gref をこの構成で参照として使う判断を保留する。延長は重ねない。非有限・発散は途中で停止し、正常でも追加 60k で必ず止める。

やらない方がよいこと: Grefx の値へ「近づいた」ことを合格条件にすること、Gref の約 1 % を格子差から差し引くこと、戻り始めた θ_r から漸近値を決めること。新しい細格子も今は追加しない。参照の反復状態が未確定のまま比較対象を増やしても、今回の判定不能は解消しない。

呼び出し側の前提への異議: **「参照側が怪しいから候補側を採れる」は成立しない。** Gref の問題は Grefx の不合格を意味しないが、適格な物差しを不要にもしていない。また、Q_w の s_ref = 0.013 % は、正負をまたぐ過渡の窓平均で小さくなった可能性があり、感度が小さい証拠にはならない。対案は上記の資格確認を先行すること。重大度 **Major**、根拠は [plan:604](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:604)・[plan:605](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:605)。

不足情報: 同じ起点から B0 を続けた対照がないため、方式変更と自然な継続の寄与は分離できない。今回の提案もその因果関係は確定しない。指定の閲覧制限に従い、数値・VERDICT は plan 記録の引用である。ファイル変更・forge 実行は行っていない。**plan 未反映**。採用する場合は呼び出し側が §4.28 と §5.1 #13 に新試験として記録し、§4.27 の判定を保存する。
