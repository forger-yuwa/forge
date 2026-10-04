# codex 諮問 (diagnose): twophase-g3a-result

- **brief**: [`notes/reviews/briefs/2026-10-04-twophase-g3a-result.md`](../../notes/reviews/briefs/2026-10-04-twophase-g3a-result.md)
- **plan**: [`plans/active/condensation-two-phase-default.md`](../../plans/active/condensation-two-phase-default.md)
- **date**: 2026-10-04
- **commit**: `4fbbecfa` (feature/species-transport)
- **codex**: effort `high`, 5.3 min, rc=0
- **結論**: **既定化は保留し、保存済み状態で「蒸発の Q1 ソースに半径 0 の粒子数を含める／含めない」だけを変える 0 step A/B を行う。**
- **extra**: `solver_density_cuda/cuda_forge/condensationSource_d.cuh`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論: **既定化は保留し、保存済み状態で「蒸発の Q1 ソースに半径 0 の粒子数を含める／含めない」だけを変える 0 step A/B を行う。**

採否表:

| 重大度 | 対象 | 採否・根拠・対案 |
|---|---|---|
| **Major** | (a) 領域外で Q0 を下げ、射影を置き換える | **却下**。破れている Q1·Q3 ≥ Q2² は Q0 に依存しない。Q0 をどう変えても、この違反は修復できない。現行射影も Q3 を固定して Q1・Q2 を変更している。[`condensationRealizability_d.cuh:263`](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/condensationRealizability_d.cuh:263)。対案は、**領域外に出る前の蒸発作用素で消滅を扱う**こと。 |
| **Major** | S_Q0 < 0 だけを追加する試算・修正 | **そのままでは却下**。H_up の瞬間微分には S_Q0 が現れない。必要なのは、消滅する粒子を取り除く操作と、その後の正半径分布から作るソースの整合。現行 S_Q1 は q0 全量を使う。[`condensationSource_d.cuh:348`](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/condensationSource_d.cuh:348)。対案は下記 A/B。 |
| **Major** | Q2・Q1・Q0 の FAIL をすべて「射影との釣り合い」で説明 | **Q1・Q2 は支持、Q0 は要再検証**。`case/16.nozzle_wys/run_0604_g3a_op_561_g2_onprod_r1/` の Q0 は abs(ΣR)=3.5814e13、B=7.8455e13、閾値=5.4903e13。残差単独は閾値内で、**B だけでも閾値を超える**。また更新診断の Q0 射影量は 0。[`G3A_pair0.txt:22`](/home/sano/work/forge-species/notes/investigations/2026-10-04-twophase-g3/evidence/G3A_pair0.txt:22)、[`G3B_RESULT.txt:116`](/home/sano/work/forge-species/notes/investigations/2026-10-04-twophase-g3/evidence/G3B_RESULT.txt:116)。Q0 は「規定上 FAIL、射影原因とは未確定」と分けて記録する。 |
| **Major** | 最終場 3 本の ΔF を事前登録の時間平均と同等に扱う | **診断用には採用、正式な代替としては却下**。事前登録は本番出口流束系列そのものの準定常判定を要求している。[`condensation-two-phase-default.md:105`](/home/sano/work/forge-species/plans/active/condensation-two-phase-default.md:105)。反復間の散らばりでは、全反復に共通する時間ドリフトを検出できない。現結果は「最終場平均を用いた変更後評価」と明記し、正式評価では保存済み時系列から本番面流束を再評価する。 |
| **Major** | 対処 (b) の方向と作業範囲 | **蒸発の r=0 境界処理を整合させる方向を採用。ただし実装方式は未確定**。これは二相拡散の既定値変更ではなく、ON/OFF 共通の蒸発モデル変更。[`condensation-source-limiter-steady.md:57`](/home/sano/work/forge-species/plans/accepted/condensation-source-limiter-steady.md:57)。**別 plan** に分離し、本 plan の前提条件にする。OFF の変化もモデル修正として検証し、旧結果との一致を合格条件にしない。 |

第 1 仮説: **実現可能領域の上側境界にあるゼロ半径の重みを、現行 S_Q1=q0ṙ が蒸発する粒子数として数えているため、蒸発ソースが境界の外を向く。** 確度: **高〔局所機構〕／中〔観測全体の主因〕**

  根拠: 保存量を q₀,q₁,q₂、q₃=ρg/(4πρ_l/3) と置く。温度・物性を固定し、q1e・q2e の上限が非作動なら、現行コードは

  S₁=q₀ṙ、S₂=2q₁ṙ、S₃=3q₂ṙ。

  D=q₁q₃−q₂² とすると、**Ḋ=ṙ(q₀q₃−q₁q₂)**。境界 D=0 では、w₁=q₁²/q₂、w₀=q₀−w₁ より、

  **Ḋ=ṙq₃w₀、Ḣ_up=ṙw₀/q₁。**

  したがって ṙ<0、w₀>0 なら、**無限小のソース段階ですでに外向き**になる。有限刻みだけの問題ではない。式の根拠は [`condensationSource_d.cuh:323`](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/condensationSource_d.cuh:323)。

  保存された `run_0604` の集計も整合する。全成分正の 4831 節点中、蒸発分岐は 4810。Ḣ_up のソース寄与は中央値 −3.1615e5 /s、負の割合 99.3 %、二相拡散は中央値 +9.1309e3 /s。[`BPUSH_run_0604…txt:1`](/home/sano/work/forge-species/notes/investigations/2026-10-04-twophase-g3/evidence/BPUSH_run_0604_g3a_op_561_g2_onprod_r1.txt:1)。

  **ただし「実際の液滴がゼロ半径に到達した履歴」を確認したわけではない。** 現在のモーメントがそのような分布を表すことと、そこへ至った物理・数値過程は分ける。

  反証条件: 正確な ρ_l、実際のソース値、上限処理の作動を照合した同一節点で、下記 B にしても外向きソース寄与が A の半分以上残るなら、この機構を主因とする説を棄却する。

第 2・第 3 仮説: **現時点では追加しない。** 対角不一致は `case/16.nozzle_wys/run_0592_pj_update_A/` 対 `run_0593_pj_update_B/` で射影比・外向き割合とも B/A=1.0000、`VERDICT: REJECT`。主因候補へ戻す根拠はない。[`PJ_VERDICT.txt:317`](/home/sano/work/forge-species/notes/investigations/2026-10-04-twophase-g3/evidence/PJ_VERDICT.txt:317)。

判別 A/B: **新規 run・更新なし。既存 `tp_operator.h5` に対するソースの反実仮想評価を 1 組だけ行う。**

- 対象は `case/16.nozzle_wys/run_0604_g3a_op_561_g2_onprod_r1/` と同じ ON の残り 2 本。`run_0592` の射影集合を固定し、全成分正・蒸発分岐・上側境界近傍を抽出する。境界判定は ρ_l=1000 近似をやめ、射影と同じ物性で **abs(D)/(q₁q₃+q₂²) ≤ 1e−6**。q1e/q2e の上限作動点は別集計する。
- **A:** 記録された現行ソース。
- **B:** **S_Q1 だけ**を `S_Q1ᴮ = S_Q1ᴬ·q₁²/(q₀q₂)` に置換する。これは同じ ṙ のまま q₀ を正半径側の重み w₁ に替える試算。格納状態・S_Q0・S_Q2・S_g・輸送項は固定する。**本番修正案ではない。**
- 指標は同一集合の Ḣ_up と、外向き寄与 `N=ΣVq₁·max(−Ḣ_up,0)`。ソース単独と全残差の両方を出す。対象外節点の件数と N も残す。
- **事前登録:** 対象集合が元の外向きソース寄与 N の 90 % 以上を覆い、3 本とも N_B/N_A ≤0.1 → ゼロ半径の重みの寄与が主因という仮説を支持し、(b) の別 plan へ進む。N_B/N_A ≥0.5 → 主因説を棄却。中間・対象の不足・記録不整合は判別不能。全残差に外向きが残れば、ソース機構の支持と問題全体の解消を分ける。

やらない方がよいこと: **Q0 だけによる領域外修復、S_Q0 だけを足して H_up の符号変化を期待する試験、射影撤去、機構を確認しない延長、旧 λ スケールへの即時復帰。** 旧形は定常残差に Δτ を入れるため、既存の設計判断を逆戻りさせる。文献でも消滅流束は分布の境界情報を要する閉鎖問題として扱われている。単なる負の数ソース追加とは区別すべきである。[Pollack et al., 2019](https://iris.polito.it/handle/11583/2775692)。

呼び出し側の前提への異議:

- **H_up の段分解は、物性固定の方向微分として採用する。** ρ は確かに消える。ただし実現可能性の完全な余裕は `ln(q₁q₃/q₂²)` であり、温度変化まで含む時間微分には **−d lnρ_l(T)/dt** が必要。現在の [`boundary_push.py:38`](/home/sano/work/forge-species/notes/investigations/2026-10-04-twophase-g3/boundary_push.py:38) はこれを評価していない。**Minor:** 「物性固定の残差方向」と明記し、張り付き判定には正確な物性を使う。
- **G3-a の記録は FAIL のまま維持する。** 保存された 3 組はすべて `(1) assembly closure: ON PASS, OFF PASS`、変更後の ΔF 評価では `(2) ON all: w UNDETERMINED, g PASS, Q2 FAIL, Q1 FAIL, Q0 FAIL`。Q1・Q2 は誤差上界を引いても閾値超過なので、実質的な未達は明瞭。一方、**Q0 の FAIL は同じ因果説明にまとめない**。
- **総水分 w の判定不能は採用するが、「ON/OFF で同じ」は却下する。** 4.33e−9 の差に対し反復間誤差が 1.87e−8 なので「この評価では差を検出できない」が正しい。総水分の保存性を確認する絶対基準は、差 ΔF が検出不能なときにも機能するものとして別途事前登録が必要。
- **G3-a だけでは「射影が釣り合いを取る」とは証明できない。** Q1・Q2 については G3-b の更新量と射影量の打ち消しが補強する。作用素残差と、DPLUR・局所刻みを経た更新補正を、そのまま同じ単位の収支として足さない。

不足情報: 対象の `run_0561–0569`・`run_0592–0593`・`run_0604–0609` の実体はローカルに存在せず、今回は保存された `evidence/*.txt`・CSV とコードを照合した。**HDF5 からの再集計、G2 の収束・準定常 VERDICT 原票の独立確認は未実施**。次の A/B には AWS 側の節点状態・正確な物性・ソース記録が必要。本番出口流束時系列の準定常 VERDICT も不足している。

ファイル変更・forge 起動は行っていない。**plan 未反映**。呼び出し側で §5.1 #4g3n と §6 に本判断を反映する。
