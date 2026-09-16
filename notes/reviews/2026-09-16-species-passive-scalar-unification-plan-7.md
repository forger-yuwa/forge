# codex レビュー: species-passive-scalar-unification (plan)

- **plan**: [`plans/active/species-passive-scalar-unification.md`](../../plans/active/species-passive-scalar-unification.md)
- **stage**: `plan`
- **date**: 2026-09-16
- **commit**: `fbefd706` (feature/sern-design)
- **codex**: effort `high`, 6.4 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M6/m1
- **focus**: 7 回目 (plan-6 NO-GO M1–M5/m1 の反映; §4.7 v5)。実装済みで単体試験は通っている (tests/unit/test_passive_fct.cu (i)–(vi), test_cond_limiter_steady (l))。変更点: (1) 低次 BE 作用素の質量流束を流れの BE 形連続式と整合する ṁ^eff = (1/a)ṁ^{n+1} + (c/a)ṁ^eff,n にし φ≡1 を厳密解にする (plan-6 M1 の反例を解消); (2) node 逆流境界は許容な外部状態 φ^n_own の定数 RHS + 流出のみ対角 (M 行列を保つ; 差は境界反拡散で外部側 R=1); (3) 実現可能性: 相対許容 ε で判定し境界へ、退化 (x≤δ または y≤δ²) は (Q0,g) 保存の単分散再初期化、他は許容領域 A={x≤1, x²≤y≤√x} へのユークリッド最近点 (x,y 両方を動かす); 特異境界の整合条件は float では判定できず ε 内の境界点は近傍の実現可能状態と区別しないと明記; (4) 後処理は g クランプ → EOS 更新 → その T で Q1,Q2 射影 → G/H/ṁ^eff; (5) G/H/ṁ^eff を FCT 有効時の必須 checkpoint 履歴にし欠落は全系 BDF1; (6) r_H/a。問い: (a) ṁ^eff の再帰 (流れの sub-iter 残差の局所残り H_ρ を無視) で上限 q_L≤ρ^{n+1} の主張はどこまで正しいか、H_ρ の扱い; (b) 逆流境界の外部状態に φ^n_own を使う妥当性 (代替: 境界値 bvar); (c) 退化しきい δ=1e-3 と最近点射影 (重みなし) の妥当性; (d) 残る論理の穴と §6-2/§6-6 ゲート (check_passive_budget.py の FAIL 判定を含む) の十分性。§4.1–4.6, §4.8 の再指摘は不要。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**  
有効質量流束・逆流境界 RHS・必須 checkpoint 履歴の方向は妥当です。  
ただし、収支ゲートの偽 PASS、実現可能性射影の欠陥、EOS 後処理の不整合が残り、現行 v5 を実装確定の根拠にはできません。

目的は正しく同定されています。[plans/README.md:26](/home/sano/work/forge/plans/README.md:26) と accepted plan を確認しました。F-sp1／F-cf8 の未解決部分を引き受けており、重複実装ではありません。node 限定も今回の指定として妥当です。以下は §4.7 とその検証に限定します。

1. **Major — `check_passive_budget.py` が、拒否すべき入力を PASS にする**

   **根拠:** [check_passive_budget.py:52](/home/sano/work/forge/solver_density_cuda/tools/check_passive_budget.py:52) は補正項を個別に比較し、計画の「合計 ≤1e-6」を検査していません。有限性の検査もなく、FCT／クランプ記録の欠落は無視します。

   実際のスクリプトに、ファイルを作らずメモリ上のログを与えて確認しました。

   | 入力 | 実測結果 |
   |---|---|
   | `floor_rel=NaN`, `lim_rel=NaN` | **PASS、exit 0** |
   | `floor_rel=6e-7`, `lim_rel=6e-7`、合計 `1.2e-6` | **PASS、exit 0** |
   | floor 行だけで FCT 記録なし | **PASS、exit 0** |

   また [speciesTransport_d.cu:1450](/home/sano/work/forge/solver_density_cuda/cuda_forge/speciesTransport_d.cu:1450) は総量がゼロなら補正比を **0** にします。消滅した成分に補正量が残る場合も、これでは予算超過を検出できません。

   **対案:** 必須成分・必須診断・最終 step の完備性と有限性を検査し、成分ごとに補正の**絶対量の合計**を判定してください。分母ゼロで補正が非ゼロなら FAIL。ログに印字された丸め済み相対値ではなく、十分な桁数の生積算値から計算し、上記を失敗系試験に追加すべきです。

2. **Major — 非周期の「境界・ソース込み保存 ≤1e-6」を、現在の境界診断では検証できない**

   **根拠:** [passiveFct_d.cuh:349](/home/sano/work/forge/solver_density_cuda/cuda_forge/passiveFct_d.cuh:349) は実現流束を
   \[
   G=F_L+\alpha A^{pre}
   \]
   としますが、356 行の境界診断は **`abs(α Apre)` のみ**です。符号、低次流束、HO から引き戻した量を保持していません。

   境界で `α=0, Araw≠0` なら、FCT は総量を変えますが `bndExchange` は **0** です。さらに [check_passive_budget.py:62](/home/sano/work/forge/solver_density_cuda/tools/check_passive_budget.py:62) は境界交換を表示するだけで、保存残差を計算していません。

   **対案:** 実際の `G` の符号付き境界積分、局所ソース、ピン交換、後処理補正、初期・最終総量を独立に記録し、閉じた収支を検査してください。`H := 増分−divG` をそのまま「ソース」として収支を閉じると、反復誤差まで恒等的に吸収してしまいます。**物理ソースと残差・補正由来の局所項を分ける**必要があります。

3. **Major — 退化判定が `(x,y)=(0,0)` を見逃し、`δ=1e-3` は微小違反を巨大補正に変える**

   **根拠:** [condensationRealizability_d.cuh:18](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationRealizability_d.cuh:18) は、不等式違反がなければ退化判定の前に return します。そのため `Q0,Q3>0, Q1=Q2=0` は**変更なし**になります。非負半径分布では `Q1=0` なら `Q3=0` であり、この状態は実現不能です。特異ケースには半正定値性に加えて整合条件が必要です。[Curto–Fialkow, Theorem 5.1](https://hjm.math.uzh.ch/v017n4/0603CURTO.pdf)

   同じ分岐を評価すると、以下も確認できます。

   - `x=1e-3, y=1e-6(1−2e-6)`：下側境界からの絶対違反は約 **`2e-12`**ですが、出力は **`(1,1)`**。
   - `x=1.001e-3` の同程度の違反：通常の最近点射影。

   閾値の両側で補正が不連続です。「小さい正の `x`」は実現不能の根拠にもなりません。例えば `x=5e-4,y=1e-4` は領域の内部です。単体 [test_cond_limiter_steady.cu:375](/home/sano/work/forge/solver_density_cuda/tests/unit/test_cond_limiter_steady.cu:375) は `(0,0)` と閾値両側を含みません。

   **対案:** 厳密ゼロ・アンダーフロー由来の特異不整合を通常の不等式検査より先に扱い、小さい正の状態には最近点射影を使ってください。単分散再初期化は明示した修復方針として限定し、収支超過なら拒否するべきです。無次元座標での重みなしユークリッド距離自体は合理的な選択ですが、**`δ=1e-3` による強制再初期化は採用を勧めません**。

4. **Major — EOS 後の処理が「Q1/Q2 だけの射影」になっていない**

   **根拠:** [main.cpp:1855](/home/sano/work/forge/solver_density_cuda/main.cpp:1855) は旧温度で `condensationPrimitive_d_wrapper` を実行し、EOS 更新後の 1861 行で同じ関数を再実行します。

   しかしその関数は [condensationTransport_d.cu:169](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationTransport_d.cu:169) から、**`g` の上限・消滅処理を含む**クランプを呼びます。[condensationRealizability_d.cuh:114](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationRealizability_d.cuh:114) 以降では更新後の `T,P` によって消滅判定が変わり、125 行で `g` を再びゼロにできます。その後に EOS 再更新はありません。

   したがって「2 回目では `g` は変わらない」というコメントは保証されません。旧温度で先に射影することも、最終温度で一度だけ射影する操作とは異なります。

   **対案:** **`g` の制約・消滅確定 → EOS → Q1/Q2 専用射影 → 最終検査 → G/H** に関数を分離してください。EOS 後に `g` が変わる場合は、EOS の再整合が必要です。温度変化で消滅条件を跨ぐ試験を追加すべきです。

5. **Major — 線形解・HO 残差の合成ノルムが、トレーサと小さいモーメントの未収束を隠す**

   **根拠:** [passiveFct_d.cuh:77](/home/sano/work/forge/solver_density_cuda/cuda_forge/passiveFct_d.cuh:77)、155 行、188 行は、異なる単位・桁の全受動種を同じ二乗和へ加えます。[speciesTransport_d.cu:1708](/home/sano/work/forge/solver_density_cuda/cuda_forge/speciesTransport_d.cu:1708) はその単一ノルムで Jacobi を終了します。

   例えば `rhs_Q0=1e15, rhs_Xi=1, residual_Xi=1` では、トレーサの相対残差が **1** でも合成相対残差は **`1e-15`** です。

   実際に
   `case/44.vitiated_air_wt/run_0257_passiveD_order_bdf2_dt8e-6_nsub40/res_200.h5`
   を再読込すると、最大値は `roQ0_0=5.53e13`、`rog_0=1.18e-3`。桁差は実データにもあります。この run の再実行した `check_convergence.py` の判定は **`NOT CONVERGED`** であり、ここでは瞬時場のスケール確認にだけ使用しています。[run 索引](/home/sano/work/forge/case/44.vitiated_air_wt/README.md)

   **対案:** 成分ごとに残差を正規化して**全成分が条件を満たすこと**を要求してください。低次解の未達と HO 残差の最大値は全物理 step を通じて保持し、最終 step の表示だけで済ませないこと。§6-6 の全量判定と同じ粒度が必要です。

6. **Major — `ṁeff` は改善だが、記載された上限成立条件はまだ不十分**

   **根拠:** [plan:192](/home/sano/work/forge/plans/active/species-passive-scalar-unification.md:192) の「`φⁿ≤1` かつ `Hρ≥0` で上限成立」は、受動種側の有効ソース・履歴項を落としています。

   一定 `Δt`、非ピン行で、密度残差を
   \[
   r_\rho=B\dot m-M(a\rho^{n+1}-b\rho^n+c\rho^{n-1})
   \]
   とし、
   \[
   E_\rho^{n+1}=M(\rho^{n+1}-\rho^n)-B\dot m^{eff,n+1}
   \]
   と置くと、
   \[
   E_\rho^{n+1}=\frac caE_\rho^n-\frac1a r_\rho^{n+1}.
   \]
   **現在 step の残差がゼロでも、過去の残りは消えません。**

   上限の十分条件は境界 RHS も含めた **`f ≤ Lρ`** です。`Hρ≥0` だけでは足りません。現行単体 [test_passive_fct.cu:296](/home/sano/work/forge/solver_density_cuda/tests/unit/test_passive_fct.cu:296) 自身が、密度一定・無流束で、正の局所履歴 `H=0.4` により **`qL=31/30>ρ=1`** となる例を持っています。

   **対案:** `Eρ` と `Lρ−f` を実際の作用素から評価し、上限・下限の成立条件を明記してください。正の局所履歴による上限逸脱も対象です。事後拒否を選ぶことは可能ですが、**有界性保証とは区別し、M1 の修正済みゲートへ接続する必要があります**。検証には、変動密度での `φ≡1`、符号の異なる密度残差、過去残差だけが残る複数 step を追加してください。

7. **Minor — v5 の式が実装と同期していない**

   **根拠:** [plan:195](/home/sano/work/forge/plans/active/species-passive-scalar-unification.md:195) の展開式は、低次移流を `ṁeff` に変更した後も旧式のままです。HO と LO の upwind 側も、両質量流束の符号が異なれば一致しません。196–197 行には依然として `rH/a` の係数がありません。

   **対案:** 展開式を廃して、
   \[
   A^{raw}=\frac1aF_H+\frac caG^n-F_L(\dot m^{eff},q_L),\qquad
   M(q_H-q_L)=BA^{raw}+r_L-\frac{r_H}{a}
   \]
   を正本にしてください。実装 [passiveFct_d.cuh:224](/home/sano/work/forge/solver_density_cuda/cuda_forge/passiveFct_d.cuh:224) の直接差分は、この方向になっています。

**推奨は、流束形 BDF2＋post-step FCT を維持し、上記を修正してから §6-2／§6-6 の統合検証へ進むことです。** 優先順は、①偽 PASS と保存収支、②退化射影と EOS 処理分離、③成分別残差と密度整合条件、④追加反例・FCT 有効の時間次数・履歴欠落 restart です。#21 にこれらを明示し、単体 PASS だけで完了扱いにしないことを勧めます。

逆流境界の `φⁿ_own` は、現行の Neumann 型境界に対する低次補助状態として妥当です。指定流入値がある境界ではその `bvar` とピンを使い、境界反拡散が作動した場合の交換を M2 の符号付き収支で検証してください。FCT 履歴欠落時の全系 BDF1 化は、コードで確認できました。

ファイル変更・新規 CFD run は行っていません。CUDA は `cuInit=100` でデバイスを認識できず、GPU 単体試験の ALL PASS は再確認できていません。提案は **plan 未反映**です。

指摘数: Critical 0 / Major 6 / Minor 1
