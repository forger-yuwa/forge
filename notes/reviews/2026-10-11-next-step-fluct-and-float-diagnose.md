# codex 諮問 (diagnose): next-step-fluct-and-float

- **brief**: [`notes/reviews/briefs/2026-10-11-next-step-fluct-and-float.md`](../../notes/reviews/briefs/2026-10-11-next-step-fluct-and-float.md)
- **plan**: [`plans/active/architecture-float-state-double-geometry.md`](../../plans/active/architecture-float-state-double-geometry.md)
- **date**: 2026-10-11
- **commit**: `461963f4` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 3.1 min, rc=0
- **結論**: **既存結果の全期間の有限性と判定ツールの正常実行を確認したうえで、`FORGE_FREEZE_TURB=1` の FP64・float・float 再実行を比較し、「SST の更新を止めても精度差が残るか」を次の一手とする。**
- **extra**: `plans/active/axisymmetric-freestream-hoop-gauge.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論: **既存結果の全期間の有限性と判定ツールの正常実行を確認したうえで、`FORGE_FREEZE_TURB=1` の FP64・float・float 再実行を比較し、「SST の更新を止めても精度差が残るか」を次の一手とする。**

採否表（Critical なし）:

| 重大度 | 対象・採否 | 根拠と対案 |
|---|---|---|
| **Major** | 「凍結すると k・ω・μt が固定される」：**却下** | フラグが飛ばすのは SST 更新ブロック（[main.cpp:2380](/home/sano/work/forge-integ-1005/solver_density_cuda/main.cpp:2380)）。原始量は毎回 `k=roK/ρ`、`ω=roOmega/ρ` と再計算される（[dependentVariables_d.cu:307](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/dependentVariables_d.cu:307)）。壁の ω・ρω の上書きも残り、μt は密度・速度勾配・物性にも依存する（[ransBoundary_d.cu:63](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/ransBoundary_d.cu:63)、[turbulent_viscosity_d.cu:198](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/turbulent_viscosity_d.cu:198)）。**「SST の輸送更新を停止する試験」と定義し直す。** |
| **Major** | 「差が消えれば発生源は k・ω の式」：**却下**。凍結試験自体は**採用** | 更新停止は、SST 内で発生する誤差と、平均流の誤差を SST が増幅する経路の両方を切る。差が消えても両者は分離できない。結論は「この精度差の発達には SST の更新・フィードバックが必要だった」に限る。 |
| **Major** | 現在の結果判定を無条件に確定：**要再検証** | `tr_an.py` は主判定を先に確定し、後で実行する収束ツールの失敗・発散を判定へ戻さない（[tr_an.py:59](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/tr_an.py:59)、[同:104](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/tr_an.py:104)）。`snapshot()` の有限性検査も ρ・速度・T・k に限定される（[cold_series.py:29](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/cold_series.py:29)）。**全残差列・必要保存量の有限性、ツールの実行成否を先に確認して再判定する。未収束と発散・判定不能は区別する。** |
| **Major** | 「G1・G1x の揺れは格子そのものの性質」：**要再検証** | 記録は異なる格子・履歴・判定期間の比較であり、格子だけの効果を取り出していない（[hoop plan:230](/home/sano/work/forge-integ-1005/plans/active/axisymmetric-freestream-hoop-gauge.md:230)）。**「その実行条件で G0 より大きかった」に留める。**変換の主因説を退ける範囲も、試験済みの G0・キー0・登録窓に限定する。 |
| **Minor** | キー0で静まるという観察：**採用**。逆切替案：**設計の骨格は採用、今回の投入は後回し** | 登録判定「支持しない」と、後半の J が小さい事実は両立する（[hoop plan:240](/home/sano/work/forge-integ-1005/plans/active/axisymmetric-freestream-hoop-gauge.md:240)）。ただし逆切替後の **J の比だけ**では、新しい過渡や対照側の偶然小さい分母を拾う。実施時は傾きを除いた RMS・絶対振幅・固定した複数窓も事前条件に含め、元の判定を書き換えない。 |

第 1 仮説: **float による θ_r の偏りの発達には、SST の状態更新を含むフィードバックが大きく関与している。誤差の発生源が SST 内かは未確定。**　確度: **中**

  根拠: [float plan §6.27](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:722) に記録された表を計算し直すと、θ_r の六条件はすべて成立する。ただし W1 の効果／再実行差は約 **10.7〜11.3倍**で、10倍基準に対する余裕は小さい。収縮部の k・ω 系に差が残る記録は仮説に整合するが、先に離れた変数が誤差の発生源とは限らない。

  反証条件: SST 更新停止下でも、両 float 腕が凍結 FP64 腕から D の方向へ **0.5|D| 以上**離れ、再実行差を十分上回る状態が登録した両窓・三断面で続くこと。この場合、「観測規模の偏りには SST 更新が必要」という仮説を棄却する。

第 2 仮説: 平均流・物性・乱流粘性の再評価など、**凍結後にも残る経路**の精度差だけで偏りが発達する。確度: 中。発生箇所は未確認。

判別 A/B:

- **共通条件**: §6.26 と同じ `_tr/q32_init.h5`・格子・キー1・実効設定を使い、三腕とも `FORGE_FREEZE_TURB=1`。既存試験から追加する介入はこのフラグだけ。化学種の凍結は混ぜない。
- **腕**: A＝FP64、B＝float、B′＝float 再実行。各 **30,000 step 固定、500 step ごと**。窓は従来と同じ W1＝20,000〜25,000、W2＝25,000〜30,000。
- **見る量**: θ_r 三断面を主判定、Q_w を独立した副判定とする。各腕の起点からの移動 Δを計算し、**δB＝ΔB−ΔA、δB′＝ΔB′−ΔA** を比較する。凍結 FP64 自身の移動を差し引くため、従来の `|ΔA|≤0.1|D|` は要求しない。ただし ΔA 自体と場の変化は必ず記録する。
- **差が残る場合**: 両窓・三断面すべてで、δB・δB′が D の方向へ 0.5|D| 以上、かつ平均効果が `10|δB−δB′|` を超える → **第1仮説の必要性を棄却し、第2仮説を支持**。
- **差が抑えられる場合**: 両窓・三断面すべてで `max(|δB|,|δB′|)+10|δB−δB′|≤0.1|D|` → **SST 更新停止下でも従来規模の差が発達する説を退ける**。SST 内の演算が発生源だとは断定しない。
- **中間・窓間不一致・欠損・非有限**: 判別不能。結果を見て期間を延長しない。上の閾値は診断用の暫定値であり、信頼区間ではない。

凍結が効いた証拠として、境界ピン対象を除く内部の `roK`・`roOmega` の変化を確認する。k・ω・μt の変化は別に記録する。**μt 一定を合格条件にはしない。**

やらない方がよいこと: **SST 状態・commit・ソース・拡散を、現時点で順番に double 化し始めること。** `S_lostfrac` は符号を持たない更新誤差の集計で、累積する偏りを示さない。次の実装箇所は、同じ入力での「残差」「要求更新」「実反映更新」のどこに説明力のある差があるかで決める。拡散についても、旧 V3 の局所診断を長期軌跡全体の除外証拠に拡張しない。

呼び出し側の前提への異議: **`ALL STEADY` を「十分に静まった」の根拠には採らない。** 今回の呼出しは既定閾値で、ドリフト5%、変動幅10%である（[check_quasisteady.py:543](/home/sano/work/forge-integ-1005/solver_density_cuda/tools/check_quasisteady.py:543)）。plan 記載の VERDICT は `check_convergence: CHECK FAILURES`、`check_quasisteady: ALL STEADY（既定閾値）`。これは小さな反復変動や定常解の同等性を保証しない。

不足情報: 対象の `case/45.isobutane_m6_d155/run_0485_kab_k1/`〜`run_0490_cv_old_k0/`、判定 JSON、残差・保存場・正式 VERDICT はローカルに無く、**実測の独立再集計は未実施**。所在は [case README の run 一覧](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/README.md:208) に記録されている。特に `cv_an.py` 自体は準定常ツールを呼ばないため、その判定記録も必要。ファイル変更・forge 起動は行っていない。**plan 未反映**であり、反映は依頼どおり呼び出し側が行う。
