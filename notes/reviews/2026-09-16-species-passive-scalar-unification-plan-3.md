# codex レビュー: species-passive-scalar-unification (plan)

- **plan**: [`plans/active/species-passive-scalar-unification.md`](../../plans/active/species-passive-scalar-unification.md)
- **stage**: `plan`
- **date**: 2026-09-16
- **commit**: `15e2ec6d` (feature/sern-design)
- **codex**: effort `high`, 7.2 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M4/m2
- **focus**: 3 回目 (result レビュー 2 回目 NO-GO の M1–M3 への設計対応)。新設した §4.6 (モーメント非負化を成分ごとに; 共通 θ 停止の撤回), §4.7 (dual-time 限定の Zalesak 型 FCT: 低次予測 ρφ_L* = ρφ^k + Δt/(aV) R_low, 局所極値 + [0,1]/≥0, P±/Q±/R±, 面ごとの α_f を両 CV で共有, 周期 node の gather 順, 定常・RK は対象外とする理由), §4.8 (周期 node のリミッタ 2 段化 [極値 gather → ψ 評価 → min gather], 流れ側 limiter_r1_d も同じ) の設計と §5.1 #20–#24 / §6 の試験計画が妥当か。特に: (a) dual-time の sub-iter 反復形 FCT で収束時に有界性と保存が両立するか、a (BDF 係数) の扱い、ソース込み低次予測で核生成の新極値を許す扱い; (b) 定常に FCT を掛けない判断 (Δτ 依存の固定点) の是非と、定常の起動緩和 θ_b を非保存のまま残す扱い; (c) 周期 group の極値 gather に float の順序保存 int 化 atomicMax/Min を使う点; (d) 前制限 (Zalesak prelimiting) の採否; (e) 見落としている整合点 (ピン行、node 境界半割面、φ_N·δρ 項との順序、passiveImplicitCoupling 0/1 両方での適用)。既存部分 (§4.1–4.5) の再指摘は不要。実装はまだ (設計のみ)。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**

§4.6 の共通停止撤回と §4.8 の周期群統合は妥当です。しかし §4.7 は、低次予測の有界性と BDF2 の履歴項を処理できておらず、保存・有界性の保証が成立しません。  
追加設計の実装前に、FCT の定式化と境界・更新処理の契約を修正してください。

`plans/README.md` と関連 `accepted/` を確認しました。F-sp1／F-cf8 の未完部分を扱う目的は適切で、今回の FCT・周期リミッタ修正を解決済みとする別計画は見当たりません。以下は追加設計への指摘です。

1. **Major — 低次予測が有界とは限らず、BDF2 では低次方程式を厳密に解いても負値が出る。**

   **根拠:** [plan:174](/home/sano/work/forge/plans/active/species-passive-scalar-unification.md:174) の `Δt/(aV)` は、[main.cpp:1713](/home/sano/work/forge/solver_density_cuda/main.cpp:1713) の BDF 係数と整合しています。ただし、BDF 残差を代入すると
   \[
   q_L^*=\frac{bq^n-cq^{n-1}}a+
         \frac{\Delta t}{aV}R_{\mathrm{space,low}}(q^k)
   \]
   です。これは低次陰解法を解いた値ではありません。BDF1 の初回反復でも、周期一次風上、`ρ=V=1`、物理 CFL=2、`q=[1,0,0,0]` なら予測は **`[-1,2,0,0]`** になります。`Q<0` をゼロにしても、逸脱した予測を領域内へ戻す補正は保証されません。

   さらに、独立した double の反例計算で、64セル周期一次風上、先頭8セルだけ `q=1`、残りゼロ、BDF1 起動、物理 CFL=4 とすると、**BDF2 第2 step の厳密解で `min(q)=-0.0181718904`** でした。行列の M 行列性だけでは、負の履歴係数を持つ右辺の非負性まで保証しません。

   ソースでも同じ問題があります。`dq/dt=-λq`、`λΔt=4`、`q⁰=1` なら BDF1 で `q¹=0.2`、次の BDF2 は **`q²=-0.0181818`**。空間反拡散 `A_f=0` なので、今回の面制限では直せません。「ソース込み予測を局所極値へ含める」は核生成の正の新極値を許しますが、負の予測を救済しません。

   **対案:** 有界な backward-Euler 低次解を基準として実際に求め、**BDF2 との差も制限対象に含める保存的 FCT**へ定式化を改めてください。ソースの非負性条件も明示する必要があります。現在の恒等式だけを有界性の証明として採用することには反対です。

2. **Major — FCT・再ピン・密度補正・成分別非負化の処理順が未整合。**

   **根拠:** [plan:173](/home/sano/work/forge/plans/active/species-passive-scalar-unification.md:173) はピン除去後に FCT 補正を加え、ピン行では `R±=1` としています。そのまま両 CV に加算すると、**ピン残差が復活**します。現行のピン除去と残差記録は [main.cpp:1738](/home/sano/work/forge/solver_density_cuda/main.cpp:1738) にあり、追加位置の指定が必要です。

   node 境界半割面では、[passiveKernels_d.cuh:25](/home/sano/work/forge/solver_density_cuda/cuda_forge/passiveKernels_d.cuh:25) のとおり `Pface` が未書込みです。全輸送面について一律に `F_H=ṁPface` を作ると、境界で stale 値を読みます。

   また、現行はモーメント更新後に `φ_N·δρ` を加えます（[condensationTransport_d.cu:322](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationTransport_d.cu:322)）。§4.6 の規則をその位置へ入れると、例えば `N=1, ρ_pre=1, ρ_new=0.5, θ_ud=-2` では、成分クリップでゼロにした後、密度補正で **−0.5** になります。「確定値0」という契約と一致しません。

   **対案:** §4.7 に次を処理契約として追加してください。

   - node 境界半割面は現行境界流束を維持し、当面 `A_f=0`。周期の架空半割面も補正対象外。
   - 補正の独立バッファを周期 gather → 残差へ加算 → **再ピン → 補正後残差を記録**。ピンへの交換量は境界収支として計上。
   - 成分非負化は密度補正済みの基点 `φ_Nρ_new` と整合させる。
   - 同じ補正済み残差を `passiveImplicitCoupling 0/1` の両方へ渡す。

3. **Major — §4.8 が実際の流れ側リミッタ経路を指定しておらず、非負専用 atomic も流用できない。**

   **根拠:** 流れの `ro,Ux,Uy,Uz,P` は、[limiter_d.cu:284](/home/sano/work/forge/solver_density_cuda/cuda_forge/limiter_d.cu:284) で **`limiter_r1_fused5_d`** を呼びます。`limiter_r1_d` の修正だけでは、狙っている `limiter_ro` の周期不整合が残ります。

   また [plan:183](/home/sano/work/forge/plans/active/species-passive-scalar-unification.md:183) の「非負 float」の順序保存は速度には使えません。単純な signed int ビット解釈では、`−2` と `−1` の大小関係が逆になります。FCT の低次予測も、指摘1のとおり負になり得ます。

   実測の出発点は確認できました。`case/09.Taylor-Green/run_0089_passiveC_m2_rerun0068/res_500.h5` のノード4532／2982は、`Xi`・`dXidx` が同値でも `limiter_Xi=0.985149／1.0` でした。

   **対案:** `fused5` を含む実呼出し経路を #22 に明記し、周期 node では全対象を「極値統合 → ψ評価 → min統合」に通してください。整数キーは符号付き有限 float 全域で順序を保存する変換を使い、±0・NaN/Inf の扱いも定義すること。試験には負速度、辺・角の多重周期群を含めてください。

4. **Major — #20–#22 の検証では、新しい FCT と成分別クランプの主要な失敗を検出しきれない。**

   **根拠:** [plan:188](/home/sano/work/forge/plans/active/species-passive-scalar-unification.md:188) の追加試験は主にステップ移流と平行移動です。一方、既存の時間次数実績には KEEP・一次再構成の系列があり、**FCT が作動する SLAU＋SFR 2 の時間次数を代替しません**。

   §4.6 は補正量・作動数を記録しますが、#20 が要求する「停止セルの残差寄与」とは別です。成分クリップでも `N_k=0, d_k<0` の当該成分は停止します。共通停止の撤回だけで、残差床と次数未達の解消は保証されません。

   再計測した `run_0257.../res_200.h5` では、共通 θ=0 の663ノード中661ノードに **`Q1=0` かつ `Q2>0, g>0`** がありました。非負性だけでなく、モーメント間の実現可能性も検査すべきです。現行の [condensationRealizability_d.cuh:29](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationRealizability_d.cuh:29) は、一般のモーメント不等式を保証する投影ではありません。

   **対案:** §6 の必須ゲートへ以下を追加してください。

   - **FCT有効**の滑らかな移流・核生成・蒸発で、BDF1/BDF2、`passiveImplicitCoupling 0/1`、3刻み＋sub-iter倍増。
   - `g,Q0` に加えて **`Q1,Q2`** の誤差、制限作動集合の未制限残差、`Q1²≤Q0Q2`・`Q2²≤Q1Q3`。
   - 密度が変化する圧縮流、非周期入口、周期辺・角での収支。
   - 「残差2桁低下」だけでなく、最終 BDF 残差の積分から見積もる保存誤差が `1e-6` 以下。

5. **Minor — float32 の更新に、double 集計だけで保存誤差 `1e-12` は要求できない。**

   **根拠:** [plan:188](/home/sano/work/forge/plans/active/species-passive-scalar-unification.md:188) に対し、[flowFormat.hpp:6](/home/sano/work/forge/solver_density_cuda/flowFormat.hpp:6) の保存量は float です。共有面流束を等量逆符号で加えても、各 CV への書戻しで別々に丸められます。単純な2セル計算でも、double 集計した総量差は **`7.45e-9`** でした。

   **対案:** `1e-12` は double 参照実装の代数的保存試験に限定し、本番 CUDA/float32 は `1e-6` と丸め誤差の測定で判定してください。非周期試験は総量一定ではなく、境界流束・ソースを含む収支で評価します。

6. **Minor — 前制限は採用してよいが、`P±` 集計より前に確定するべき。**

   **根拠:** [plan:178](/home/sano/work/forge/plans/active/species-passive-scalar-unification.md:178) は `P±` 計算後に前制限します。破棄予定の反拡散も分母へ数えるため、残す流束を必要以上に絞ります。

   記載された符号条件は、`A_f>0` を `i→j` の補正流束と読む限り妥当です。ただし「受け側／出し側」は **`ṁ` ではなく `A_f` の符号**で選ぶ必要があります。

   **対案:** 前制限 → `P±` 集計 → `R±` → 共有 `α_f` の順に固定してください。前制限は指摘1の低次有界性問題を解決するものではありません。

定常へ今回の物理時間 FCT を適用しない判断は支持します。ただし、定常の `θ_b` は [passiveKernels_d.cuh:124](/home/sano/work/forge/solver_density_cuda/cuda_forge/passiveKernels_d.cuh:124) でゼロになれるため、「固定点では無作用」は保証ではなく検証条件です。**補正が止まったことと、未制限残差が消えたことを分けて確認**してください。RK を対象外にするスコープも妥当ですが、非SSP RKでは FCT が不可能という一般論にはできません。RK4 の積分流束を制限する先行例があります。[Chaplin–Colella](https://arxiv.org/abs/1506.02999)

既存の非保存問題も再現しました。`case/09.Taylor-Green/run_0098_passiveC_m2_step_seam_s3_final/` の `res_0.h5 → res_500.h5` は、周期代表点・合併体積で **83.71696370 → 83.26320518、損失0.542015%**。両 run の `check_convergence.py` 通常判定は **`NOT CONVERGED`** でした。これは非定常の収支・停止状態の確認であり、定常解の主張ではありません。run 索引は [case/09 README](/home/sano/work/forge/case/09.Taylor-Green/README.md) と [case/44 README](/home/sano/work/forge/case/44.vitiated_air_wt/README.md) です。

**推奨は、§4.7 を「有界な低次解＋時間差分を含む保存的補正」に再設計し、上記1→2→3→4の順で計画を直してから実装することです。** §4.6・§4.8 の方向性は維持できますが、現行の FCT 式の実装着手は認められません。ファイルは変更しておらず、レビュー提案は **plan 未反映**です。

指摘数: Critical 0 / Major 4 / Minor 2
