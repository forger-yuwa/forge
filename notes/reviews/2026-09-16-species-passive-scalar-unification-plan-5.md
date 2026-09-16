# codex レビュー: species-passive-scalar-unification (plan)

- **plan**: [`plans/active/species-passive-scalar-unification.md`](../../plans/active/species-passive-scalar-unification.md)
- **stage**: `plan`
- **date**: 2026-09-16
- **commit**: `fbefd706` (feature/sern-design)
- **codex**: effort `high`, 6.6 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M5/m1
- **focus**: 5 回目 (plan-4 NO-GO M1–M6/m1 の反映)。§4.7 v3 のみ: (1) 生流束 A^raw を保持し q_C = q_H − M^{-1}B(A^raw − A^lim) (前制限で棄却した面も引き戻す); (2) 終了状態で assembleResidual を再評価して massflux/P_face/ソースを固定し r_H を診断、基点 q_B = q_H − M^{-1}B A^raw の物理限界逸脱を計測しクリップ、局所限界は {φ^n, φ_B} なので q_B は常に局所限界内; (3) 境界半割面の差 ṁ(φ_H − φ_L) を恒等式に含め α=1 固定 (補正には入らない)、全流束 ṁφ − J の規約; (4) 低次陰解の有界性条件 (f≥0, f≤Lρ^{n+1} [BE は離散連続式], 初期値クリップ, BDF2 の破れは記録) と線形残差による受入; (5) モーメント実現可能性: 不等式違反時に (Q0, g) 保存の単分散射影 + 違反数 log + condClampCorr を合否式へ; (6) g 補正後の依存変数再評価。問い: 残る論理の穴 (特に (2) の基点クリップの扱いと (5) の射影の選択 [Q1,Q2 を単分散に戻す = Q1,Q2 を変える] が妥当か、より良い最小射影があるか)、単体試験 (i)–(iv) と §6-2/§6-6 のゲートで失敗を検出できるか。§4.1–4.6, §4.8 の再指摘は不要。実装は plan-4 の提案どおり並行して進めている (レビューで変更が要れば反映する)。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**  
`A^raw` の保持と終了状態の再評価は妥当です。しかし、基点クリップと境界固定寄与の扱いでは有界性が成立せず、厳密解でも反例が出ます。単分散射影も、微小な違反に対して過大な変更を加えます。

対象は指定どおり §4.7 v3 と対応する検証ゲートです。`plans/README.md`・accepted の関連計画を確認し、今回の課題は F-sp1/F-cf8 の未解決部分と整合しています。以下の反例はホスト上で計算し、境界の反例は有理数演算で残差ゼロまで確認しました。作業中コードの v3 未反映部分を、そのまま新たな指摘として数えてはいません。

1. **Major — 限界計算だけの基点クリップでは、有界性と保存性を両立できない**

   **根拠:** [plan:187](/home/sano/work/forge/plans/active/species-passive-scalar-unification.md:187) は限界計算で `q_B` をクリップしますが、[補正式:193](/home/sano/work/forge/plans/active/species-passive-scalar-unification.md:193) の基点は未クリップの `q_B` です。

   流束ゼロで `q_B=-0.01` なら、限界計算をどうクリップしても `q_C=-0.01` のままです。実際の基点もゼロに置き換えれば、今度は保存量を追加します。「自身を局所極値に含める」ことは、物理限界との交差後に自身が許容内である保証にはなりません。

   BDF2 の例 `q^n=0.9, q^{n-1}=0.5` でも、無流束・無ソースの厳密解は `31/30` です。反復を収束させても解消しません。また、`r_H` は記録のみ、低次線形残差の未達も警告のみなので、現在の記述には実効的な受入拒否条件がありません。

   **対案:** 実際に加算する基点をクリップせず、その許容性を物理 step の受入条件にしてください。不成立なら履歴・時刻を進めず、原因が反復誤差なら反復を継続し、BDF 履歴・ソース由来なら step を拒否します。残差は成分別に評価し、局所誤差 `M⁻¹(r_L−r_H)` と物理限界逸脱にも閾値を置くべきです。`f=0` 用の絶対許容値も必要です。

   **検証:** 無流束の上下限逸脱、BDF2 の不適合履歴、意図的な HO/LO 未収束を単体へ追加し、**floor 前の失敗と履歴未更新**を検査してください。FCT の有界性が許容な基点に依存する点は、[Kuzmin–Möller–Turek の解析 §6](https://diamhomes.ewi.tudelft.nl/~mmoller/files/publications/TecRep221.pdf)とも整合します。

2. **Major — 境界面を `α=1` に固定するなら、その寄与を制限予算の基点へ含める必要がある**

   **根拠:** [plan:185–196 の起点](/home/sano/work/forge/plans/active/species-passive-scalar-unification.md:185)。全輸送面で `q_B` を作る一方、`P±` は内部面だけです。境界寄与が無制限に戻るため、実際の補正は

   $$
   q_C=
   \underbrace{q_B+M^{-1}B_{\partial}A^{raw}_{\partial}}_{q_*}
   +M^{-1}B_I A_I^{lim}
   $$

   となります。内部面の予算は `q_B` ではなく **`q_*`** に対して作る必要があります。node 境界で自セル値を使う実装は [passiveKernels_d.cuh:25](/home/sano/work/forge/solver_density_cuda/cuda_forge/passiveKernels_d.cuh:25) で確認しました。

   **厳密な反例:** 3 CV、`ρ=1`、BE、`M=0.5 I`、内部質量流束 `+1`、左境界 `−1`、右境界 `+1`、ソース・拡散なし。左境界はピンされていない逆流境界とします。

   | 量 | CV 0 | CV 1 | CV 2 |
   |---|---:|---:|---:|
   | `q^n` | 0.5 | 0.3 | 0.5 |
   | `q_H` | 0.9 | 0.7 | 0.5 |
   | `q_L=q_B` | 0.5 | 13/30 | 41/90 |
   | v3 による `q_C` | **1.3** | 7/18 | 37/90 |

   HO 内部面値を `(0.7, 0.5)` とすると **`r_H=r_L=0`**。`A_I^raw=(1/5,1/15)`、前制限と Zalesak 後は `(0,1/45)` となり、総量 `2.1` を保存したまま上限を破ります。反復誤差や基点クリップが原因ではありません。

   **対案:** 境界固定方針を維持するなら、`q_*` で許容性・局所極値・`Q±` を統一してください。`q_*` が不適合なら指摘1の受入拒否へ接続します。

   **検証:** この反例を非周期単体として追加してください。周期箱だけでは検出できません。単体(iii)も、内部面だけ `α=0` なら期待値は `q_*` であり、一般には `q_L` ではありません。

3. **Major — 単分散射影は実現可能な状態を作るが、微小違反に対して不連続かつ過大**

   **根拠:** [plan:200](/home/sano/work/forge/plans/active/species-passive-scalar-unification.md:200)。

   `Q0,Q3>0` を固定し、

   $$
   r_{30}=(Q_3/Q_0)^{1/3},\qquad
   x=\frac{Q_1}{Q_0r_{30}},\quad
   y=\frac{Q_2}{Q_0r_{30}^2}
   $$

   とすると、不等式は `x²≤y, y²≤x`、単分散状態は `(1,1)` です。

   `(x,y)=(0.5,√0.5+10⁻⁷)` では違反は約 `1.41×10⁻⁷` ですが、提案は `x` を **0.5**、`y` を約 **0.293** 変更します。一方、`y` を `2×10⁻⁷` 下げるだけで両不等式を厳密に満たせます。

   `Q1` は平均半径・成長速度、`Q2` は液相生成ソースに使われています。[condensationSourceF_d.cuh:131](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationSourceF_d.cuh:131)  
   したがって `g,Q0` を保存しても、次 step の物理応答を大幅に変えます。

   **対案:** 無次元 `(x,y)` 上で、重みを明記した**許容誤差付きの最小補正**を採用してください。正常な状態は保持し、特異境界は次の指摘の整合条件で扱います。単分散への一律置換は採用しません。補正予算を超える場合は step を拒否します。

   **検証:** 違反を `10⁻³→10⁻⁷` と縮小したとき補正量も縮小すること、多分散状態が保持されることを追加してください。時間次数だけでは、変更された物理モデルに対して高次収束してしまう可能性があります。

4. **Major — 2本の不等式だけでは「モーメント実現可能性」を判定できない**

   **根拠:** [plan:200](/home/sano/work/forge/plans/active/species-passive-scalar-unification.md:200) と [最終場の判定:295](/home/sano/work/forge/plans/active/species-passive-scalar-unification.md:295)。

   `(Q0,Q1,Q2,Q3)=(1,1,1,8)` は両不等式を通ります。しかし、

   $$
   \int(r-1)^2\,d\mu=Q_2-2Q_1+Q_0=0
   $$

   なので分布は `r=1` に集中し、必ず `Q3=1`。`Q3=8` とは両立しません。**違反数ゼロでも非実現可能です。**

   **対案:** 非負性に加え、2つの Hankel 行列の半正定値性と、特異時の整合条件

   $$
   (Q_2,Q_3)^T\in
   \operatorname{Range}
   \begin{pmatrix}Q_0&Q_1\\Q_1&Q_2\end{pmatrix}
   $$

   を検査してください。これは [Curto–Fialkow, Theorem 5.1](https://www.math.uh.edu/~hjm/v017n4/0603CURTO.pdf) の4モーメントの場合です。float32 では無次元化した判定と許容誤差を定義します。

   `Q0=0,g>0` は両者を保存する射影自体が不可能です。`S>1` だから見逃す契約にはできません。また、`Q3` の換算に使う `ρ_l(T)` は、**g 補正後の EOS 更新と同じ温度**で確定させてください。

   **検証:** 上の特異反例、`Q0=0,g>0`、単分散・多分散・ゼロ状態を含めます。現在の §6-6 の「2不等式違反ゼロ」では、この穴を検出できません。

5. **Major — `condClampCorr` の体積積分だけでは、射影によるモーメント収支を監査できない**

   **根拠:** [plan:202](/home/sano/work/forge/plans/active/species-passive-scalar-unification.md:202) は射影を `condClampCorrQ` に記録し、合否には `condClampCorr` を追加しています。しかし既存実装では、

   - `condClampCorr`：`|Δ(ρg)|/ρ`
   - `condClampCorrQ`：各モーメントの最大相対補正

   です。[condensationRealizability_d.cuh:38](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationRealizability_d.cuh:38)  
   どちらも `Q1/Q2` の符号付き・絶対保存量補正の積算ではありません。

   実測でも、[run_0257 の res_200.h5](/home/sano/work/forge/case/44.vitiated_air_wt/run_0257_passiveD_order_bdf2_dt8e-6_nsub40/res_200.h5) は全 `VALUE` が有限ですが、`Q0>0,g>0` の6275点中、`Q1²≤Q0Q2` の相対違反が42点、`Q1=0,Q2>0` が661点あります。`condClampCorr_0` の最大は **0**、`condClampCorrQ_0` は **1** でした。これは瞬時場の監査であり、定常値の主張ではありません。再実行した `check_convergence.py` の VERDICT は **`NOT CONVERGED`** です。run 索引は [case/44 README](/home/sano/work/forge/case/44.vitiated_air_wt/README.md) です。

   **対案:** `g,Q0,Q1,Q2` ごとに、射影前後の `∫Δq dV` と `∫|Δq|dV` を周期 root のみで、各 step・全期間について積算してください。異なる単位の成分は別々に正規化し、ゼロ総量の扱いも定義します。

   **検証:** §6-2/§6-6 は全 step の実状態と、この成分別収支で判定してください。BDF 収支に使う境界流束・ソースが「補正時に凍結した値」なのかも固定し、後処理後の状態との収支を検査する必要があります。

6. **Minor — 単体(i)と double 参照の合否条件が、アルゴリズムの性質より強い**

   **根拠:** [plan:210](/home/sano/work/forge/plans/active/species-passive-scalar-unification.md:210)。滑らかであるだけでは、極値付近で前制限・FCT が無作用とは限りません。

   また、現行の [double 参照:71](/home/sano/work/forge/solver_density_cuda/tests/unit/test_passive_fct.cu:71) は float の `q_H` を入力しながら、[判定:220](/home/sano/work/forge/solver_density_cuda/tests/unit/test_passive_fct.cu:220) では物理履歴からの保存を `10⁻¹²` で要求しています。保存的補正が保存するのは入力 `q_H` の総量なので、HO の反復誤差まで除去することはできません。

   **対案:** 無作用試験は「前制限なし・全 `α=1` を満たす入力」でビット一致を検査してください。double の `10⁻¹²` は補正前後の代数的保存に適用し、物理 step 全体の収支は HO/LO 残差込みで別判定します。

**推奨は、post-step FCT を維持しつつ、実際の基点に対する受入判定を必須にする修正です。** 優先順は、①境界固定寄与込みの基点と step 拒否、②実現可能性判定と無次元最小補正、③成分別補正収支、④上記反例を含む単体、⑤§6-2/§6-6 の再検証です。現行のゲートだけでは十分ではありません。

ファイル変更・新規 CFD 計算は行っていません。レビュー提案は **plan 未反映**です。

指摘数: Critical 0 / Major 5 / Minor 1
