# codex 諮問 (diagnose): linedir-divergence-2

- **brief**: [`notes/reviews/briefs/2026-10-09-linedir-divergence-2.md`](../../notes/reviews/briefs/2026-10-09-linedir-divergence-2.md)
- **plan**: [`plans/active/tooling-nozzle-isothermal-wall-chain.md`](../../plans/active/tooling-nozzle-isothermal-wall-chain.md)
- **date**: 2026-10-09
- **commit**: `36b95969` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 8.1 min, rc=0
- **結論**: **次の実行は、`run_0203` の条件から `FORGE_FREEZE_TURB=1` だけを変える最大200 stepのA/Bに絞る。事前に既存帳簿の状態位相を揃え、平均流だけで成長するのか、SSTの更新が必要なのかを判別する。**
- **extra**: `notes/reviews/briefs/2026-10-08-linedir-divergence.md`, `notes/reviews/2026-10-08-linedir-divergence-diagnose.md`, `plans/accepted/time_integration-line-implicit-viscous-v2.md`, `solver_density_cuda/cuda_forge/setDT_d.cu`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論: **次の実行は、`run_0203` の条件から `FORGE_FREEZE_TURB=1` だけを変える最大200 stepのA/Bに絞る。事前に既存帳簿の状態位相を揃え、平均流だけで成長するのか、SSTの更新が必要なのかを判別する。**

第 1 仮説: **方向別Δτで時間対角が弱くなり、近似LHSが表していない壁近傍の平均流の結合が露出している。特に、等温壁の密度更新に伴うエネルギー・圧力応答が有力。** 確度: **中**。

  根拠:

- **[Major] 壁の「温度固定」と、陰解法の「エネルギー補正ゼロ」は同じ拘束ではない。**  
  [timeIntegration_d.cu:1053](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1053) は壁の `δroe=0` を課し、密度行は残す。一方、[nodeWallDirichlet_d.cu:83](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu:83) は更新後の密度から `roe=ρe(Tw,Y)`、`P=ρR(Y)Tw` を作る。固定組成・静止壁なら、実際の拘束の微分は **δroe=e_w δρ、δP=R Tw δρ**。現行ライン解はこのエネルギー応答を隣接点へ同時に伝えていない。温度ピンは外側更新の直後にも行われる。[main.cpp:2516](/home/sano/work/forge-integ-1005/solver_density_cuda/main.cpp:2516)  
  **対案:** 壁密度・壁圧力と第一内点の法線運動量の位相を調べ、この結合が成長モードに参加しているか確認する。不整合の存在は確認できるが、今回の主因とはまだ断定しない。

- **[Major] `convMethod: 0` でも、SLAUの残差とFVSのLHSは一致しない。**  
  ラインの結合は [timeIntegration_d.cu:890](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:890) のFVSから抽出される。SLAUの質量流束は [convectiveFlux_slau_d.inc.cuh:584](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:584) の別の式である。  
  **対案:** `run_0207` が退けるのは「二次再構成が発散に必要」という説までとし、**流束・壁拘束を含むLHS/RHSの不整合は候補に残す**。

- 方向別Δτは [setDT_d.cu:161](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/setDT_d.cu:161) で内部ライン面の制約を除き、[timeIntegration_d.cu:827](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:827) の V/Δτ を弱める。提供記録の「精度変更・sweep増加・一次化でも成長が残る」は、この機構と整合する。**Thomasが厳密に解くのは組み立てた近似行列であって、実際の残差と境界処理の全結合ではない。**

  反証条件: SST更新を止め、内部の保存乱流量が凍結されたことを確認した状態で平均流の成長が消えるなら、「平均流側の不整合だけで今回の成長を説明できる」という第1仮説は退け、SSTとの結合を主候補に上げる。成長が残っても、等温壁の拘束が主因と確定するわけではない。

第 2 仮説: **粘性・熱伝導の近似LHSが不足している。** 確度: 中〜低。既定では粘性はスカラー対角のみで、`lineViscCoupling=1` でも完全な粘性・熱伝導Jacobianにはならない。またLHSは `cfg.visc`、実際の粘性流束は局所 `vis_lam` を使う。TPでは両者の値の照合が必要。[timeIntegration_d.cu:937](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:937)、[1496](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1496)、[viscousFlux_d.cu:164](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/viscousFlux_d.cu:164)

第 3 仮説: **SSTの分離更新と平均流のフィードバックが成長に必要。** 確度: 低〜中、未検証。SSTはライン行列の外で更新されるため、`nStepInner` の増加ではこの結合を解き直せない。[main.cpp:2374](/home/sano/work/forge-integ-1005/solver_density_cuda/main.cpp:2374)

判別 A/B: **A＝元の条件、B＝`FORGE_FREEZE_TURB=1` のみ追加。最大200 step。**

- `run_0203` と同じ開始場・リミッタ基準・CFL・精度・次数・sweep数・緩和を使う。実装上、この環境変数は定常経路でSST更新を実際に止める。[main.cpp:2348](/home/sano/work/forge-integ-1005/solver_density_cuda/main.cpp:2348)、[2376](/home/sano/work/forge-integ-1005/solver_density_cuda/main.cpp:2376)
- **事前予測は「Bでも成長が残る」**。見る量は全保存量残差、壁ρ・P、第一内点から30層目までのρ・P・法線速度。初期の固定したρ・P・音速で正規化し、単一点ではなく対象帯の振幅と成長率を比較する。
- **成長が残る場合:** SSTの時間更新が発散に必要という第3仮説を外す。
- **Bだけで成長が消える場合:** 第1仮説の「平均流だけで十分」を退け、SST更新との結合を支持する。200 stepまで有限で、末尾の32 step窓3本で振幅が増大しないことを確認する。発散が遅れるだけなら保留。
- 凍結の実効性は、壁・入口ピンを除く `roK/roOmega` の不変性で確認する。**これは `k/ω/μt` の完全固定ではない**。密度・速度勾配の変化や壁境界の更新は残る。

要求された候補の優先順位は次のとおり。**今回回すのは1位だけ**とする。

| 優先 | 単因子 | 事前予測と判別範囲 |
|---|---|---|
| 1 | `FORGE_FREEZE_TURB=1` | 成長が残る予測。SST更新が必要かを分ける |
| 2 | `time.deltaT.lineViscCoupling: 1` | 単独では止まらない予測。止まれば粘性LHSの変更に感度がある。ただし対角も2α→αへ変わるので、非対角追加だけの効果とは呼ばない |
| 3 | `time.deltaT.implicitRelax: 0.7 → 0.3` | 成長の鈍化・周期の伸長を予測し、停止は未確定。安定化しても、壁・対流・SSTのどれが主因かは分からない |

CFDを回さずに行う帳簿解析は、このA/Bの比較指標を作る作業として次に限定する。

1. **状態の位相を揃える。** `after_eos_bc` を基本にし、速度は同じ保存量からも再計算する。`call` はstep番号ではなく残差組立ての呼び出し番号なので対応を確認する。
2. 同じcallで、対流＝`res_after_conv − res_before_conv`、粘性＝`res_after_viscous − res_after_sources`、射影＝`res_final − res_after_viscous` を取る。絶対値だけでなく符号と相殺を見る。ソース差分は複数処理の合計なので、直ちにhoop単独へ帰属させない。[main.cpp:1993](/home/sano/work/forge-integ-1005/solver_density_cuda/main.cpp:1993)
3. 壁—第一内点面の質量流束を、圧力差項 **ṁ_p=−ss·χ_mass·P_del/(2c_hat)** と残部に分ける。`lowMachPrecond=0` の今回ならこの式でよい。面の向きは `ic0/ic1` に合わせる。`slauWallNormalChi=0` の直接効果は、保存済みの `chi_mass` を `chi` に置き換えて評価できる。ただし、その後の状態変化までは予測できない。
4. 残差の最大項を「成長の原因」と即断しない。陰解法の実際の補正は残差を行列で解いた結果であり、単純な `R/V` 更新ではない。

やらない方がよいこと:

- **等温→断熱を最初の診断にしない。** 方程式の境界条件と熱場の過渡を変えるため、壁Jacobianだけの検査にならない。
- `slauWallNormalChi: 0` を本命視しない。低Machでは両χとも1に近い可能性があり、まず面帳簿で変更強度を測れる。
- 壁境界面のΔτ制約解除や、第一層だけの場当たり的な除外をしない。
- コードによる安定化策を選ぶなら、診断後に **Δτ=min(Δτ_directional, R·Δτ_base)** という全ライン内点共通の倍率上限を、第一層だけの例外処理より優先する。これは安定性と加速余地を測る試験であり、真因の修正とは区別する。壁の時間制約は維持する。
- 200 stepの有限性を高速化成功としない。step単価がpointの約2倍なら、同じ収支・収束ゲートまでのstep数を半分未満にする必要がある。今回の短い診断だけでは採算は決まらない。

呼び出し側の前提への異議:

- **[Major]「Uyが先行した」は再解析が必要。** `entry` はEOS再評価前で、保存量と原始量の位相が混在し得る。補正カーネルが更新するのは保存量であり、内部点の `U/P/T` はそこで更新しない。[main.cpp:1879](/home/sano/work/forge-integ-1005/solver_density_cuda/main.cpp:1879)、[update_d.cu:301](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/update_d.cu:301)  
  **対案:** `after_eos_bc` で比較し直す。なお壁温ピンはstep末尾にもあるため、`entry→after_eos_bc` の差だけではピンの全作用を測れない。

- **[Major]「軸対称hoopのJacobianがLHSに無い」は却下。** 圧力・hoopの近似Jacobianが組み込まれている。[timeIntegration_d.cu:973](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:973)  
  **対案:** 問うなら「実際のTP残差とどこまで整合しているか」であり、「欠落」ではない。

- **[Major] 面帳簿はFP64ではない。** 面値の格納先は `float*`、CSVも9桁である。[convectiveFlux_slau_d.inc.cuh:648](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:648)、[convectiveFlux_d.cu:725](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/convection/convectiveFlux_d.cu:725)  
  **対案:** 圧力差は丸め済みの `Pf_R−Pf_L` から作り直さず、差を取ってから保存した `P_del` を使う。帳簿から倍精度級の閉合を要求しない。

- 精度主因・単なるsweep不足・二次再構成必須説を優先候補から外す判断は、提供記録の範囲で採用する。**Bの「隣接Δτの比そのものが原因」と、Cの「旧格子と同じモード」は未立証**。壁近傍と下流全断面では発生場所も異なるため、共通するのは現時点では「時間対角を弱めると不安定になる」まで。

不足情報: 対象の `case/45.isobutane_m6_d155/run_0183*・run_0191*・run_0201〜0207*` と帳簿は、この作業ツリーには存在しない。したがって実効設定、初期成長率、残差・HDF5の非有限検査、各VERDICTは直接再検証できていない。特に `lineViscCoupling`、局所Δτ、局所物性、実行バイナリとの対応が必要。コード上の指摘と、提供資料に記録された観測を区別した診断である。**ファイル変更・forge起動なし。plan 未反映（依頼どおり呼び出し側が反映）。**
