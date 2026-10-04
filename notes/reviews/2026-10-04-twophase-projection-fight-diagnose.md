# codex 諮問 (diagnose): twophase-projection-fight

- **brief**: [`notes/reviews/briefs/2026-10-04-twophase-projection-fight.md`](../../notes/reviews/briefs/2026-10-04-twophase-projection-fight.md)
- **plan**: [`plans/active/condensation-two-phase-default.md`](../../plans/active/condensation-two-phase-default.md)
- **date**: 2026-10-04
- **commit**: `02a3220f` (feature/species-transport)
- **codex**: effort `high`, 4.5 min, rc=0
- **結論**: **既定化を保留し、同じ ON 保存場から、液・Q の前処理対角だけを共通化する診断 A/B を各 1 更新行う。**
- **extra**: `plans/active/condensation-two-phase-transport.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論: **既定化を保留し、同じ ON 保存場から、液・Q の前処理対角だけを共通化する診断 A/B を各 1 更新行う。**

| 重大度 | 判断対象 | 採否・根拠・対案 |
|---|---|---|
| **Major** | 射影が更新を打ち消していても、報告量が定常なら既定化する | **却下**。`case/16.nozzle_wys/run_0561_g2_onprod_r1/` の保存済み判定は **`NOT CONVERGED`**、7 報告量は **`OVERALL: ALL STEADY`**。射影込みの写像には固定点があり得るが、射影前の残差がゼロとは限らない。既存契約は離散方程式を満たすことなので、既定化は保留する。根拠：[収束判定](/home/sano/work/forge-species/notes/investigations/2026-10-04-twophase-g2/evidence/run_0561_g2_onprod_r1/CONVERGENCE_VERDICT.txt:2)、[準定常判定](/home/sano/work/forge-species/notes/investigations/2026-10-04-twophase-g2/evidence/run_0561_g2_onprod_r1/QUASISTEADY_7q.txt:11)、[更新の契約](/home/sano/work/forge-species/methods/condensation.md:1138)。 |
| **Major** | 「共通 θ の追加」「g=0 なら Q も消す」を主対策にする | **却下**。θ は既に共通。さらに射影は **ρg>0 かつ ρQ0>0** の節点でしか作動しないので、g=0 の観測は大きな射影補正の直接原因を説明しない。まず正の液が残る節点の更新方向を調べる。根拠：[共通 θ と commit](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/twoPhaseDiffusion_d.cuh:190)、[射影の作動条件](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/condensationRealizability_d.cuh:257)。 |
| **Major** | 再正規化 H3 を主因とする | **却下。ただし収支への影響は別途残す**。液・全 Q に同じ正の係数を掛ける処理は、温度固定・丸めなしでは Hankel 不等式の符号を変えない。ブリーフでも違反の発生は再正規化前。撤去する根拠はない。根拠：[共通係数の適用](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/speciesTransport_d.cu:2181)。 |
| **Major** | G3-b の帳尻が閉じたので数値補正も許容する | **却下**。帳簿の閉鎖と補正の許容性を別判定にする。ON の正式補正ゲートは射影により **FAIL**。G2 の既存 PASS は維持するが、G3 全体は未達として扱う。根拠：[G2 の正式ゲート記録](/home/sano/work/forge-species/notes/investigations/2026-10-04-twophase-g2/evidence/G2_VERDICT_v2.txt:1)、[G3 の設計](/home/sano/work/forge-species/plans/active/condensation-two-phase-default.md:93)。 |

第 1 仮説: **H1 のうち、蒸発域で液と Q の前処理が異なることが更新方向を歪め、実現可能領域の外へ押している。** 確度: **中**。

根拠:

- 蒸発分岐では全ソース対角をゼロ初期化した後、`sj_g` だけを設定し、Q の対角はゼロのまま戻る。[condensationSourceKernels_d.cuh:321](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/condensationSourceKernels_d.cuh:321)、[同:395](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/condensationSourceKernels_d.cuh:395)。
- DPLUR はその異なる分母を実際に使う。一方、共通 θ は増分の大きさを変えるだけで、実現可能性を保証しない。[twoPhaseDiffusion_d.cuh:148](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/twoPhaseDiffusion_d.cuh:148)、[condensationTransport_d.cu:887](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/condensationTransport_d.cu:887)。
- `case/16.nozzle_wys/run_0580_g3b_update_on/tp_update.h5` の保存済み集計では、**符号相殺を避けた絶対量でも**、射影量／制限後増分は全域で Q1 約 **72 %**、Q2 約 **69 %**。既存の局所領域 `x=35–70 mm、壁距離<0.4 mm` では約 **91 %、97 %**。単なる全域の符号付き和の偶然ではない。[G3B_RESULT.txt:106](/home/sano/work/forge-species/notes/investigations/2026-10-04-twophase-g3/evidence/G3B_RESULT.txt:106)、[同:231](/home/sano/work/forge-species/notes/investigations/2026-10-04-twophase-g3/evidence/G3B_RESULT.txt:231)。

反証条件: 対角を共通化しても、正の液を持つ主要補正節点で、**増分の大きさで正規化した領域外向き成分と射影比がほぼ残る**なら、「対角差が主因」は棄却する。

第 2 仮説: **対角差を除いても、蒸発ソースと輸送を合わせた残差、またはその有限増分化が実現可能性を破る。** 確度: 中。共通 θ だけでは防げない反例がある。記載された一様 ṙ のソースを、温度・液密度固定の無次元単分散状態 `(q0,q1,q2,q3)=(1,1,1,1)` に前進更新すると、`(1,1−h,1−2h,1−3h)` となり、両 Hankel 行列式は **−h²**。これは実 run の原因確定ではないが、「共通係数なら保証される」は数学的に成立しない。[methods/condensation.md:821](/home/sano/work/forge-species/methods/condensation.md:821)。

第 3 仮説: **g の床・アンダーフローによるモーメント塵は併発しているが、大きな射影収支の主因ではない。** 確度: 低。g=0 の節点は射影を通らず、後段の消滅処理で扱われるため、両集合の混同を避ける。

判別 A/B:

- 共通入力は `case/16.nozzle_wys/run_0561_g2_onprod_r1/res_48000.h5`。ON のまま、同一の更新開始状態・残差・物性・面流束・ピン条件から**各 1 更新**。
- **A:** 現行の成分別対角。**B:** 液・Q2・Q1・Q0 の対角だけを、節点ごとの `max(Dg,DQ2,DQ1,DQ0)` に共通化する。蒸気の対角、残差、sweep 数、緩和、θ、射影は維持する。これは既存 YAML の切替ではなく、診断用変更である。
- 両者で候補・floor 後・射影前後を記録する。Q3 は **本番と同じ `ρ_l(T)`** で算出し、g>0 の射影対象と g=0 の塵を別集計する。
- 指標は、既存の全域・局所領域ごとの成分別 `ΣV|C_projection| / ΣV|δq_limited|` と、Hankel 条件に対する更新方向。**射影絶対量の減少だけでは支持としない**。
- **A で再現し、B で Q1・Q2 の射影比がともに A の 0.1 以下となり、外向きの更新方向も解消** → 第 1 仮説を支持。  
  **B でも正規化した外向き成分・射影比が同程度に残る** → 対角差の主因説を棄却し、第 2 仮説を優先する。中間結果は判別不能とする。

やらない方がよいこと: 射影の撤去、許容の緩和、g=0 の節点だけを直して完了扱い、機構を見ない長時間延長。**最小修正の方向は液・Q の更新を一体として実現可能にすること**であり、共通対角化そのものを完成修正とはしない。射影を避けるため θ=0 に張り付く修正も不可。

呼び出し側の前提への異議:

- 「毎更新ほぼ打ち消す」は、現時点では **1 更新の直接計測＋末尾窓の補正ゲート FAIL** まで。各更新で同じ相殺率とは未確認。
- 「レシピ非依存」ではなく、既存記録どおり「レシピ感度が ON−OFF 差の 10 % 以下」。
- **G3 への反映案:** G3-b の帳簿閉鎖は維持し、射影を含む数値補正の許容判定を追加する。既存 κ の全域ゲートに加え、事前固定した局所領域で、各更新・各成分の `ΣVΣ操作|C| / ΣVq_start` を評価する。丸め上界・物理的消滅・境界供給は別計上する。G3-a の残差／出口流束差の判定も必要で、擬似時間の補正量を単純に Δτ で割って物理流束の代用にしない。**局所ゲートは今回の提案であり、過去の G2 判定を遡って変更しない。**

不足情報: 生の `tp_update.h5` はローカルに無く、節点別の対角・増分・温度・射影作動集合を独立再集計できなかった。したがって **更新経路の問題という診断は支持できるが、対角差が真因との確定はできない**。ファイル変更・forge 実行はしていない。**plan 未反映**。呼び出し側で `condensation-two-phase-default.md` §5.1 #4g3・#4g3b／§6 に採否を反映すること。
