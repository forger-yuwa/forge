# codex レビュー: species-passive-scalar-unification (plan)

- **plan**: [`plans/active/species-passive-scalar-unification.md`](../../plans/active/species-passive-scalar-unification.md)
- **stage**: `plan`
- **date**: 2026-09-16
- **commit**: `fbefd706` (feature/sern-design)
- **codex**: effort `high`, 5.3 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M5/m1
- **focus**: 6 回目 (plan-5 NO-GO M1–M5/m1 の反映)。§4.7 v4 のみ: 核心は BDF2 を流束形の BE に書き換える点 — 前 step の増分率 I^n = (V/Δt)(q^n − q^{n−1}) を面流束 G^n と局所残り H^n に厳密分解 (確定状態から作るので恒等式は常に厳密) し、BDF2 = BE + F^eff = (1/a)F + (c/a)G^n、S^eff = (1/a)S + (c/a)H^n とする。低次 BE 陰解 q_L は有界 (BE, M 行列)、履歴の負係数は G^n として反拡散流束 A^raw = F^eff_H − F_L(q_L) に入り Zalesak で制限される。境界半割面も Zalesak に入れ外部側 R=1。基点 q_B = q_H − M^{-1}B A^raw (全面) の物理限界逸脱は反復誤差と負の S^eff (蒸発・局所 H) だけ → 記録し事後の合否で拒否 (物理 step の再試行は枠組みに無い)。実現可能性: 許容領域 x≤1, x²≤y≤√x (ε 内側余裕で特異整合を回避) への最小補正、成分別収支。問い: (a) 流束形 BDF2 の恒等式 (b=a+c) と G/H の分解・確定状態からの再構成・restart (G/H 欠落時は全部 H) に穴が無いか; (b) 低次 BE と流れ側 BDF2 の質量流束の不整合 (上限 ρ^{n+1} の O(Δt²) の差) の扱い; (c) 局所 H^n の制限不能な寄与 (前 step のクランプ・ソース) の扱いと記録で十分か; (d) 実現可能性の許容領域と ε 射影の妥当性 (H1, H2 ⪰ 0 + 特異整合); (e) 単体 (i)–(vii) と §6-2/§6-6 のゲート。§4.1–4.6, §4.8 の再指摘は不要。実装は並行して進めている。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**  
固定刻みでの BDF2 の書き換えと `G/H` 分解は成立します。しかし、低次 BE の有界性と実現可能性射影には反例があります。  
§4.7 v4 の設計修正が必要です。並行実装中の未完成部分そのものは指摘に数えていません。

1. **Major — 流れ側 BDF2 と低次 BE の不整合は、反復誤差ゼロでも基点を密度上限の外へ出す**

   **根拠:** [plan:185](/home/sano/work/forge/plans/active/species-passive-scalar-unification.md:185) は不整合を認識していますが、[plan:192](/home/sano/work/forge/plans/active/species-passive-scalar-unification.md:192) は基点逸脱の原因を「反復誤差と負の `S_eff`」に限定しています。これは誤りです。

   有理数演算で、2 CV、`V=Δt=1`、内部質量流束 `0→1` に `1`、ソース・拡散なし、密度履歴 `(1,1)` を確認しました。BDF2 の密度解と一様トレーサ `φ=1` は次のとおりです。

   | 量 | CV 0 | CV 1 |
   |---|---:|---:|
   | `ρ_H=q_H` | `1/3` | `5/3` |
   | `q_L=q_B` | `1/4` | `7/4` |
   | `q_B/ρ_H` | `3/4` | **`21/20=1.05`** |

   `G^n=H^n=0`、HO/LO 残差ともゼロです。上限逸脱 `1/12` は反復を増やしても消えません。最終補正がこの例を救う場合でも、**有界な基点を前提とする証明と基点逸脱ゲートは成立しません**。`O(Δt²)` は符号条件や `1e-6` 以下を保証しません。

   **対案:** 密度の BE 形連続式と整合する低次面質量流束を構成し、`Lρ^{n+1}−f≥0` を成立条件にしてください。この圧縮流反例と、一様 `φ=1` の保持を単体・§6-2 に追加すべきです。

2. **Major — node の逆流境界では、低次作用素が M 行列である前提も破れる**

   **根拠:** [plan:186](/home/sano/work/forge/plans/active/species-passive-scalar-unification.md:186) の境界半割面は、符号付き `ṁ/ρ_own` を対角へ加えます。[実装:108](/home/sano/work/forge/solver_density_cuda/cuda_forge/passiveFct_d.cuh:108) も同じです。

   ピンされていない逆流境界を持つ1 CVで、`V=Δt=1`、`ρ^n=ρ^{n−1}=1/5`、流入質量流束 `2` とすると、残差ゼロの BDF2 密度は `23/15`。ところが低次 BE は、

   \[
   L=1-\frac{2}{23/15}=-\frac7{23},\qquad
   f=\frac15>0,\qquad q_L=-\frac{23}{35}.
   \]

   **非負 RHS でも負解になります。** 境界の外部側 `R=1` は反拡散制限の修正であり、この低次作用素の問題を解消しません。M 行列性は陽性保証の前提です。[Kuzmin–Möller–Turek §4.1](https://diamhomes.ewi.tudelft.nl/~mmoller/files/publications/TecRep221.pdf)

   **対案:** 低次の逆流境界を、許容な外部状態による流入 RHS と非負の流出対角に分け、HO との差を境界反拡散へ含めてください。低次対角を小さい正数へ置換して解いたことにする処理は不可です。単体(v)には、上限だけでなく低次作用素の成立確認も必要です。

3. **Major — ε 射影は区間が空になり、`x=0` の非実現可能状態も残す**

   **根拠:** [plan:203](/home/sano/work/forge/plans/active/species-passive-scalar-unification.md:203) の式をそのまま評価しました。

   `x=1, ε=10⁻⁶` では `x'=0.999999` に対して、

   - 下限：`x'²(1+ε)=0.999998999999`
   - 上限：`√x'(1−ε)=0.999998500000375`

   となり、**下限が上限を約 `5×10⁻⁷` 上回ります**。内側区間が非空となるには、少なくとも
   \[
   x'\le\left(\frac{1-\varepsilon}{1+\varepsilon}\right)^{2/3}
   \]
   が必要です。

   また `Q0,Q3>0, Q1=0` では `x'=0` のまま、射影は `Q2'=0` にします。しかし非負半径分布で `Q1=0` なら `Q3=0` が必要なので、依然として非実現可能です。特異整合条件はεを付けるだけでは回避できません。[Curto–Fialkow, Theorem 5.1](https://www.math.uh.edu/~hjm/v017n4/0603CURTO.pdf)

   これは実入力にも存在します。`case/44.vitiated_air_wt/run_0257_passiveD_order_bdf2_dt8e-6_nsub40/res_200.h5` を再読込すると、`g,Q0>0` の6275点中、**`Q1=0,Q2>0` が661点**でした。全 `VALUE` は有限です。再実行した `check_convergence.py` は **`NOT CONVERGED`**。これは瞬時場の監査であり、定常性の主張ではありません。[case/44 run 索引]( /home/sano/work/forge/case/44.vitiated_air_wt/README.md)

   **対案:** 正常な単分散・多分散状態を保持し、特異整合を明示判定する射影にしてください。補正時には `x,y` の両方を動かし、距離・重み・許容誤差を定義します。現在の逐次クランプは「最小補正」の定義にもなっていません。単体(vii)に `x=0`、単分散、単分散近傍を追加してください。

4. **Major — 実現可能性射影と EOS 更新の順序が、同じ温度を使う契約を満たさない**

   **根拠:** [plan:199](/home/sano/work/forge/plans/active/species-passive-scalar-unification.md:199) は「実現可能性 → EOS 再評価」の順でありながら、射影の `ρ_l(T)` は再評価後の温度を使うとしています。現行の [main.cpp:1833](/home/sano/work/forge/solver_density_cuda/main.cpp:1833) も射影が EOS 更新に先行しています。

   射影時に `Q2²=Q1Q3` でも、後続 EOS 更新で `ρ_l` が2%増えれば `Q3` は `1/1.02` 倍になり、最終状態は
   \[
   \frac{Q_2^2}{Q_1Q_3}=1.02
   \]
   です。ε余裕では吸収できません。

   **対案:** `g` の上限・消滅処理を確定 → EOS 更新 → その温度で `Q1,Q2` を射影 → 最終実現可能性検査、の順に分離してください。再度 `g` を変える処理があるなら EOS も再更新します。`G/H` の構築は、この最終検査後です。

5. **Major — `G/H` 欠落 restart の「全部 H」は、制限不能な影響を1 stepに限定できない**

   **根拠:** [plan:179](/home/sano/work/forge/plans/active/species-passive-scalar-unification.md:179) の代入は恒等式を保ちますが、制限可能性を保ちません。

   2 CV、`ρ=V=Δt=1`、履歴を `(0.5,0.5)→(0.9,0.1)`、現在の流束・ソースをゼロとした厳密計算では、

   - 正しい履歴面流束 `G^n=−0.4, H^n=0`：基点 `(0.9,0.1)`。履歴反拡散を制限できる。
   - 欠落時 `G^n=0,H^n=(0.4,−0.4)`：基点 **`(31/30,−1/30)`**。面制限では修復できず、floor の絶対補正は総量比 **`1/15≈6.67%`**。

   さらに無ソース・無補正でも `H^{n+1}=(c/a)H^n=H^n/3` と残ります。「その1 stepだけ制限不能」ではありません。

   **対案:** FCT 有効時の `G/H` を checkpoint の必須履歴に含め、欠落時は既存契約と同じく**全系を BDF1 へ揃えて再開**してください。§6-6 に「状態履歴だけ存在し、FCT 履歴だけ欠ける restart」を追加します。ログだけで継続する利点は、この失敗を許容する費用に見合いません。

6. **Minor — BDF 残差の係数と単体(iii)の期待値が不整合**

   **根拠:** [plan:184](/home/sano/work/forge/plans/active/species-passive-scalar-unification.md:184) と [診断実装:64](/home/sano/work/forge/solver_density_cuda/cuda_forge/passiveFct_d.cuh:64) の `r_H` は元の BDF 残差です。この定義なら正しい差分恒等式は、

   \[
   M(q_H-q_L)=BA^{raw}+r_L-\frac{r_H}{a}
   \]

   です。[plan:191](/home/sano/work/forge/plans/active/species-passive-scalar-unification.md:191) は `1/a` が欠けています。また単体(iii)の `q_B−q_L` は低次線形残差だけでは決まりません。

   **対案:** `r_H/a` を BE 換算残差として定義し直すか、式に係数を明記してください。意図的な HO 未収束入力で恒等式を検査し、残差判定はスケールが大きく異なる各成分について個別に行います。

**推奨は、流束形 BDF2 と post-step FCT を維持し、実際の基点の許容性を保証する設計へ修正することです。** 優先順は、①密度整合・逆流境界を含む低次作用素、②実現可能性射影と EOS 順序、③FCT 履歴の一括復元、④上記反例の単体化、⑤§6-2/§6-6 の統合検証です。

目的は既存 accepted plan の残件に対応しており、重複実装ではありません。境界面を Zalesak に含める修正と成分別収支も妥当です。ただし、局所 `H` を含む基点逸脱に対する**事後拒否は結果の採用条件にはできても、有界性保証の代わりにはなりません**。再試行機構は不要ですが、逸脱・補正予算超過を必須の `FAIL` 判定に結び付ける必要があります。

ファイル変更・新規 CFD run は行っていません。提案は **plan 未反映**です。

指摘数: Critical 0 / Major 5 / Minor 1
