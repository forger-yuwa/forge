# codex レビュー: species-passive-scalar-unification (plan)

- **plan**: [`plans/active/species-passive-scalar-unification.md`](../../plans/active/species-passive-scalar-unification.md)
- **stage**: `plan`
- **date**: 2026-09-16
- **commit**: `15e2ec6d` (feature/sern-design)
- **codex**: effort `high`, 7.4 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M6/m1
- **focus**: 4 回目 (plan-3 NO-GO の M1/M2/M4/m1/m2 への再設計)。§4.7 を「物理 step 末尾の保存的 FCT 補正」に書き換えた: sub-iter が収束した完全陰的 HO 解 q_H はそのまま、同じ BDF 履歴・質量流束・凍結ソース・拡散で移流だけ 1 次風上にした低次陰解 q_L を Jacobi sweep で実際に解き、2 つの離散式の差が厳密に面流束 A'_f = ṁ(P_face − φ_L,up) + [J(q_H) − J(q_L)] の和になることを使って、Zalesak (前制限 → P± → 局所極値 {φ^n, φ_L} ∩ 物理限界 → R± → 面共有 α_f) で q^{n+1} = q_H − Δt/(aV) Σ (1−α_f) A'_f とする (α=1 の面では q_H 不変 = 時間 2 次の完全陰的スキームは不変、保存、q_L の解誤差は限界にしか効かない)。BDF2 の負係数・負ソースで q_L が僅かに負になり得る点は「保証の穴」として floor 収支で定量化。処理契約 (内部面のみ、補正後の再ピンと境界収支、floor、実現可能性、coupling 0/1 共通)。§6-2/§6-6 のゲートも plan-3 M4/m1 を反映 (S3 での次数試験、Q1/Q2、実現可能性不等式、float 1e-6)。§4.8 は実装済 (§5.1 #22, fused5 含む, CAS 比較 atomicMax/Min)。問い: (a) この post-step 補正で有界性 (q_L が有界な範囲) と保存、時間 2 次の両立の論理に穴が無いか; (b) 最終 sub-iter の状態 (P_face, ṁ, S) と収束状態の差 (sub-iter 誤差) が α=1 の恒等式に入る扱い; (c) q_L の Jacobi 反復の有界性の主張 (M 行列, 凸結合) と圧縮性 (Σ流入≠Σ流出) での限界; (d) 局所極値に φ^n と φ_L を使う選択 (q_H を入れない理由: q_H が制限対象); (e) モーメントで上限を局所極値にする扱いと核生成の新極値 (ソース込み q_L が自身の限界に入る); (f) 実装上の落とし穴 (周期 gather 順, ピン行, 拡散差分の面計算)。§4.1–4.6, §4.8 の再指摘は不要。実装はこれから (設計のみ)。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**  
目的と統一方針は妥当ですが、§4.7 の「有界な低次解に基づく補正」という証明が、前制限・反復誤差・境界面の扱いで成立していません。  
`α=1` で `q_H` を保持することと、内部面補正の総量保存は成立します。それだけでは有界性と時間2次の両立を保証できません。

`plans/README.md`、`accepted/`、設定・仕様・検証手順を確認しました。F-sp1／F-cf8 の受け皿として重複はなく、周期箱と凝縮ノズルの選択も適切です。cell 検証の除外はユーザー指定として扱います。以下は今回の §4.7 と、その受入条件に限定した指摘です。

1. **Major — 前制限で `A'_f` を上書きすると、高次解と低次解の差の恒等式が壊れる。**

   **根拠:** [plan:185](/home/sano/work/forge/plans/active/species-passive-scalar-unification.md:185) は前制限で `A'_f=0` とし、[plan:189](/home/sano/work/forge/plans/active/species-passive-scalar-unification.md:189) はその流束を使って `q_H` から補正を引きます。しかし、`q_H−q_L` を表すのは**前制限前の流束差**です。消した成分が `q_H` に残り、限界計算からだけ抜けます。

   読み取り専用のPythonで、FCT段の3CV周期・BE・CFL=1の代数試験を実行しました。密度・体積・時間刻みを1とし、履歴 `[1,0,0]`、固定面値 `[0,0.5,0]` とすると、

   - 厳密低次解: `[4/7, 2/7, 1/7]`
   - 高次解: `[1, −1/2, 1/2]`
   - 記載どおり前制限・補正した結果: `[5/7, −3/14, 1/2]`

   **低次解は有界で、線形解誤差もないのに負値が残ります。** これはforge本体のrunではなく、補正式の単体反例です。

   **対案:** `Araw` を保存し、前制限後・Zalesak制限後の流束を `Alim` として、補正を
   \[
   q_C=q_H-M^{-1}B(A_{\rm raw}-A_{\rm lim}),\qquad
   M_{ii}=aV_i/\Delta t
   \]
   と定義してください。`B` は面流束を両CVへ逆符号で加える作用素です。前制限で棄却した面では `Alim=0` でも、引き戻す `Araw` は残します。

2. **Major — `q_L` の解誤差は「限界にしか効かない」のではなく、補正後状態の限界逸脱に直接残る。**

   **根拠:** [plan:178](/home/sano/work/forge/plans/active/species-passive-scalar-unification.md:178)、[plan:183](/home/sano/work/forge/plans/active/species-passive-scalar-unification.md:183)、[plan:189](/home/sano/work/forge/plans/active/species-passive-scalar-unification.md:189)。

   同じ履歴・ソースを `h` に含め、離散残差を
   \[
   r_X=h+BF_X-Mq_X
   \]
   と定義すると、正しい関係は
   \[
   q_C=q_L+M^{-1}BA_{\rm lim}+M^{-1}(r_L-r_H)
   \]
   です。Zalesakが制限する最初の2項に、**制限されない残差差**が加わります。限界ちょうどのCVでは、小さな誤差でも逸脱します。

   実装も [main.cpp:1730](/home/sano/work/forge/solver_density_cuda/main.cpp:1730) で固定回数反復し、残差評価後に状態を更新して終了します。終了状態での再評価・収束受入判定はありません。また、[speciesTransport_d.cu:1423](/home/sano/work/forge/solver_density_cuda/cuda_forge/speciesTransport_d.cu:1423) の増分制限が残るため、更新停止を未制限HO方程式の収束と見なせません。

   **対案:** 最終更新後の同一状態で `massflux`・`Pface`・拡散係数・ソース・BDF残差を再評価し、固定してください。実際の基点
   \[
   q_B=q_H-M^{-1}BA_{\rm raw}
   \]
   の許容範囲を確認し、その基点から制限を組む契約に変更することを推奨します。基点が不適格なら追加反復し、解消しなければそのstepを受け入れない設計が必要です。

   `α=1` で補正がゼロになること自体は反復誤差に依存しません。ただし、その状態の時間精度にはHO反復誤差が残ります。

3. **Major — 非周期境界の `A'_f=0` は、記載された低次陰解と両立しない。**

   **根拠:** [plan:184](/home/sano/work/forge/plans/active/species-passive-scalar-unification.md:184) は境界半割面を除外します。しかし [passiveKernels_d.cuh:29](/home/sano/work/forge/solver_density_cuda/cuda_forge/passiveKernels_d.cuh:29) の境界流束は `ṁ q_i/ρ_i` です。HO・LOで同じ式でも、`q_H≠q_L` なら流束は異なります。特に非ピンの出口で差が残ります。

   2CVのBE代数試験でも、入口流束1、内部HO流束0.75、出口は自セル値とすると、`q_H=[0.25,0.375]`、`q_L=[0.5,0.25]`。内部面だけを完全に低次へ戻しても `[0.5,0.125]` となり、低次解を回復しません。

   拡散の符号も明文化が必要です。既存コードの `J` は [passiveKernels_d.cuh:177](/home/sano/work/forge/solver_density_cuda/cuda_forge/passiveKernels_d.cuh:177) でCV0へ **加算**されます。この変数を再利用するなら、CV0→CV1向きの全流束は `ṁφ−J` です。

   **対案:** 境界流束差も離散式の差に含め、境界交換として収支に計上してください。ピン行の反力・再ピン補正も別に記録します。内部面・境界面・拡散で共通の符号規約を定義し、`α=0` で低次解を回復する単体試験を必須にしてください。

4. **Major — M行列だけでは、圧縮性流れのprimitive上限とJacobi各反復の有界性を保証できない。**

   **根拠:** [plan:178](/home/sano/work/forge/plans/active/species-passive-scalar-unification.md:178) の「上限は流入側の凸結合」「BEでは厳密に有界」は条件不足です。

   低次行列を `L`、右辺を `f` とすると、トレーサの `0≤q_L≤ρ^{n+1}` を保証する十分条件は、M行列性に加えて
   \[
   0\le f\le L\rho^{n+1}
   \]
   です。BE・ソースなしなら、同じ質量流束による離散連続式が満たされることで上限を導けます。`Σ流入≠Σ流出` 自体が問題なのではなく、**密度の時間変化との整合が必要**です。Jacobiの初期値も許容範囲内でなければなりません。

   BDF2では負値だけでなく上限側も破れます。無流束・密度1・履歴 `q^n=0.9, q^{n−1}=0.5` なら、右辺は正でも `q_L=1.033333…` です。

   また、周期一様移流のJacobi誤差減衰率は `C/(1+C)`。CFL=100では100反復後も初期誤差の約37%が残ります。「相対変化≤1e-6」は残差や限界誤差の保証ではありません。

   **対案:** 非負の面係数、連続式との整合、初期値、上下限双方の成立条件を明記し、残差による受入判定を置いてください。周期群では近傍和だけでなく面対角も合算し、合併体積のBDF項を一度だけ加える契約まで低次solveに記載する必要があります。

5. **Major — 成分別FCTと既存クランプでは、モーメントの実現可能性を保証できない。**

   **根拠:** [plan:186](/home/sano/work/forge/plans/active/species-passive-scalar-unification.md:186) の成分別上下限は、モーメント間の不等式を含みません。[condensationRealizability_d.cuh:28](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationRealizability_d.cuh:28) の既存クランプも、非負化・液相上限・液滴消滅を処理するだけで、`Q1²≤Q0 Q2`、`Q2²≤Q1 Q3` を保証しません。

   無次元モーメントで、実現可能な単分散状態 `[1,1,1,1]` と `[1,2,4,8]` の成分別範囲内にある `[1,2,1,8]` は、`Q1²=4>Q0 Q2=1` です。各成分の上下限を守ることでは防げません。

   実測でも `case/44.vitiated_air_wt/run_0257_passiveD_order_bdf2_dt8e-6_nsub40/res_200.h5` に第1不等式の相対違反が73ノードありました。ただし `Q1≤5.45e-31` の微小量であり、有意な物理誤差の証拠とは扱いません。同runの `check_convergence.py` 再実行は **`NOT CONVERGED`** です。

   **対案:** モーメントベクトルの許容領域と `ρg≤ρY_w` を面補正の制約に含めてください。単に4本で同じ係数を使うだけでも一般の多面更新の保証にはならず、更新後ベクトルの許容性を保証する必要があります。後段クランプによる各保存量の変更も、現在の合否式から抜けている `condClampCorr` を含めて収支判定してください。

6. **Major — `g` のpost-step補正後に二相EOSを更新する契約が抜けている。**

   **根拠:** [plan:193](/home/sano/work/forge/plans/active/species-passive-scalar-unification.md:193) の後処理は「クランプ→primitive→mirror」です。しかし [condensationTransport_d.cu:138](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationTransport_d.cu:138) のprimitive処理は、二相の `T/P` を更新しません。温度反転は [dependentVariables_d.cu:156](/home/sano/work/forge/solver_density_cuda/cuda_forge/dependentVariables_d.cu:156) にあり、`g` に依存します。

   `g` を補正した後、旧 `T/P` を残すと、出力・消滅判定・次のソース評価で状態が不整合になります。また、補正後のソースと流れの残差は、補正前の完全陰的解の残差とは別です。

   **対案:** 確定保存量から二相EOS・物性・必要な境界状態を更新し、その後に出力とcheckpointを確定してください。補正後のソース差と連成残差を測定し、**実際に `α<1` が生じる**凝縮試験で時間次数を判定することが必要です。

7. **Minor — 局所限界の逸脱量を `floorCorr` で代用できない。**

   **根拠:** [plan:196](/home/sano/work/forge/plans/active/species-passive-scalar-unification.md:196) は限界逸脱量を「floorが処理した量」と同一視しています。しかし [passiveKernels_d.cuh:58](/home/sano/work/forge/solver_density_cuda/cuda_forge/passiveKernels_d.cuh:58) のfloorは物理限界のみです。局所上限0.4に対する値0.5は、トレーサの物理範囲内なので検出されません。

   **対案:** floor前の局所上下限逸脱、物理限界逸脱、補正後逸脱を分けて計測してください。モーメント不等式は適切な無次元化と絶対・相対許容値を併用し、最終場だけでなく各物理stepで確認します。

**推奨は、post-step FCTの方向を維持し、上記の離散式契約を修正してから実装することです。** 優先順は、①生流束差と境界を含む恒等式、②残差を含めた有界な基点の認証、③モーメント許容領域とEOS更新、④その失敗を検出する単体・次数試験です。

`{φ^n,φ_L}` を局所限界に使い、制限対象の `q_H` を除く選択は妥当です。ソース込み `φ_L` を含めれば、低次解に現れた核生成の新極値も許容できます。ただし、滑らかな極値で高次精度が維持される証明にはなりません。引用された [Chaplin–Colella](https://arxiv.org/html/1506.02999v1) も、滑らかな極値の限界を別途設計しています。現計画の次数ゲートは維持し、未達を許容しないことを推奨します。

ファイル変更なし。以上の修正提案は **plan未反映** です。

指摘数: Critical 0 / Major 6 / Minor 1
