# codex 諮問 (diagnose): twophase-drop-nonneg-theta

- **brief**: [`notes/reviews/briefs/2026-10-02-twophase-drop-nonneg-theta.md`](../../notes/reviews/briefs/2026-10-02-twophase-drop-nonneg-theta.md)
- **plan**: [`plans/active/condensation-two-phase-transport.md`](../../plans/active/condensation-two-phase-transport.md)
- **date**: 2026-10-02
- **commit**: `856dbc43` (feature/species-transport)
- **codex**: effort `high`, 6.7 min, rc=0
- **結論**: ?
- **extra**: `notes/investigations/2026-10-02-twophase-diffusion-kernel-design.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（対象 HEAD: `856dbc43`。run の数値はブリーフ・plan 記載値であり、今回再測定していない）

| 論点 | 採否 | 理由・対案 |
|---|---|---|
| Q1：定常擬似時間の試行値で負を許す | **採用。ただし診断用 opt-in に限定** | 契約は元の離散方程式の固定点なので、中間値の射影は検討できる。ただし「射影後に動かない」は残差ゼロを意味しない。独立残差の条件を維持する。[plan:88](/home/sano/work/forge-species/plans/active/condensation-two-phase-transport.md:88) |
| Q1：非負 θ だけを外し、既存クランプへ任せる | **却下 — Major、最重要の穴** | 現行 commit は蒸気が負になると **液を削らず総水分を増やす**。OFF 側ではこの処理を変更し、総水分の floor・再正規化を別計上した上で、固定した総水分に対して液を射影する。[twoPhaseDiffusion_d.cuh:213](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/twoPhaseDiffusion_d.cuh:213) |
| Q2：負値が物理評価へ渡らない保証 | **要再検証 — Major** | EOS 自体は負の g を排除しない。局所的な `max` の有無だけでなく、更新から次の物理評価までの呼出順を確認する必要がある。 |
| Q3：既存 `[cond-corr]` に κ を適用 | **そのままでは却下 — Major** | 更新ごと・成分ごとの補正を測れていない。下記の定義へ変更する。κ は再正規化の丸め尺度であり、一般のクランプ誤差の保証にはならない。 |

**結論:** OFF 側の commit と補正計測を下記の契約に揃え、既定値を維持したまま、DPLUR 同士で非負制限 ON/OFF の 2000 step A/B を一組行う。

**第 1 仮説:** 共通の非負 θ による局所停止を解除すれば、元の方程式を補正で代替せずに、二相系の反復を改善できる。**確度: 中**

- **根拠:** `ρg=0、δρg<0` では [非負制限](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/twoPhaseDiffusion_d.cuh:195) が θ をゼロにし、蒸気・液・Q の全増分を止める。ブリーフによれば `case/16.nozzle_wys/run_0516_twophase_dplur/` の末尾200更新では液の非負制限による停止が20977件あり、閾値制限による停止はない。
- **反証条件:** OFF 側で停止を解除しても、後述の独立残差の改善基準を満たさない、または補正が許容を超える。この場合、「制限解除だけで実用的に改善できる」という仮説を棄却する。

**第 2 仮説:** 停止は全域未収束の主因ではなく、解除しても凝縮域下流の Q 残差が残る。確度: 中。`run_0511_twophase_relax05_diag3/` では直接支配説が200/200更新で棄却されたという記録がある。ただし点対角・緩和0.5の結果を、DPLUR の停止セルへそのまま移せない。[plan:130](/home/sano/work/forge-species/plans/active/condensation-two-phase-transport.md:130)

**判別 A/B:**

A は `condensation.condTwoPhaseNonnegLimit: 1`、B は `0`。両側とも二相拡散 ON、DPLUR、同じビルド・初期保存量 `run_0482 res_48000`・実効設定で2000 step。変える設定はこのキーだけ、評価窓は更新1800–1999とする。

試験前に、次を実装契約・事前基準として固定する。

1. **OFF 側の補正経路を明確にする。**  
   `dg_max`・`dT_max` は残す。`vround` の `wnew=gnew` を一般の蒸気負値に適用しない。総水分の floor・再正規化後に `0≤ρg≤ρY_w` を保証し、Q の floor・射影も含め、物理評価へ渡る状態を非負にする。既存 ON 側の数値動作は維持する。

2. **補正は各更新・各成分で測る。**  
   q ∈ {ρY_w, ρv, ρg, ρQ2, ρQ1, ρQ0} に対して、

   `C_q,n = Σ_i Σ_a |Δq_i,n^(a)| V_i / Σ_i q_i,n^開始 V_i`

   とする。a は commit 内補正・floor・上限クランプ・モーメント射影などの補正段階。各段階の**実際の格納値の差を double で評価**し、相殺させず、重複計上しない。分母には負を含み得る試行値を使わない。分母ゼロなら補正も厳密ゼロを要求し、非有限は不合格とする。

   B の数値補正は `max_n C_q,n ≤ κ = 4.7683716e−7` を要求する。これは今回の事前許容であり、保存の証明ではない。再正規化は既存 `[renorm-gate]` で別判定する。液滴消滅も理由別に残し、負値の修復を「物理的消滅」に紛れ込ませない。

3. **改善は θ の件数ではなく独立残差で判定する。**  
   末尾200更新の Q2・Q1・Q0 の独立残差最大値が、それぞれ A の末尾最大値と自身の初期値の両方の0.1倍以下、または既存の独立残差許容以下になること。総水分・蒸気・液の残差は、許容を超える範囲で A より悪化しないこと。各更新の補正後状態は有限・非負であること。

   再正規化は、A が κ 以下の項目は B も κ 以下、A が超過する項目は B で悪化しないことを診断条件とする。**正式受入では全項目 κ 以下を要求する。**

→ **B が全診断条件を満たせば**第1仮説を支持する。**満たさなければ**「この制限解除だけで十分」を棄却する。2000 step の結果で正式収束・定常解の一致は主張しない。正式受入には、元の独立残差監査と `check_convergence` の PASS、残した閾値制限・θ_src の解除も必要である。

**やらない方がよいこと:** θ=0 件数の消失や、クランプ後の非負だけを成功条件にすること。補正が小さいという理由で元の独立残差許容を緩めること。今回の結果を見る前に既定を OFF にすること。

**呼び出し側の前提への異議:**

- **「既存クランプで受ける」は現行 commit の説明になっていない。**  
  θ_thr=1、`ρY_w=0.01、ρg=0.0099、δρv=−0.0002、δρg=0` の float32 代数検算では、試行総水分は約0.0098だが、現行 `vround` は約0.0099へ増やす。補正は初期総水分の約1％で、丸め補正ではない。forge は実行していない。

- **既存監視は提案ゲートを直接評価できない。**  
  `Qcut` は次元の異なる Q2・Q1・Q0 の補正を一つに加算している。[twoPhaseDiffusion_d.cuh:207](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/twoPhaseDiffusion_d.cuh:207) また、更新後は再正規化、`passiveBounds` の順であり、後段の実現可能性クランプだけを見ると先行補正を落とす。[condensationTransport_d.cu:1193](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/condensationTransport_d.cu:1193) `[cond-corr]` の区間累積を総量で割る方式も、更新ごとの窓内最大値とは異なる。[同:642](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/condensationTransport_d.cu:642)

- **液クランプが崩す保存の種類は区別すべき。**  
  実現可能性クランプ自体は総水分・全エネルギーを変更せず、液量・モーメントを変更する。したがって問題は液・モーメント方程式への人工補正と EOS 応答である。一方、上記 `vround` は総水分そのものを変更する。

**不足情報:**

負値の消費経路について、許可された読取範囲で確認できたのは以下までである。

- EOS の `cond_twophase_resid` は g を直接使用し、非負化しない。[condensationEOS_d.cuh:351](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/condensationEOS_d.cuh:351)
- `cond_vapor_state_f` は `Yw−g` を非負化するが、負の g 自体は防がない。負の g が来れば蒸気量を過大評価し得る。[condensationSourceF_d.cuh:64](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/condensationSourceF_d.cuh:64)
- 実現可能性クランプも、総水分が非負であることを前提にしている。負の `gmax` では下限処理後に再び負を代入し得る。[condensationRealizability_d.cuh:219](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/condensationRealizability_d.cuh:219)
- 二相拡散は蒸気・液・Q の入力を直接使う。[twoPhaseDiffusion_d.cuh:72](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/twoPhaseDiffusion_d.cuh:72)

主セッションでは、更新後から EOS・核生成／成長／蒸発・気相組成による物性評価・次の移流拡散まで、補正済み保存量と primitive／ghost が渡る呼出順を確認する必要がある。**「全消費箇所で安全」は今回未確認。**

ファイル変更・forge 実行は行っていない。**plan 未反映**。反映先は `plans/active/condensation-two-phase-transport.md` §4.2、§5.1、§6。
