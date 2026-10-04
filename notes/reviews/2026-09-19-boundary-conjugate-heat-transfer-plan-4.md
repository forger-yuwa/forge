# codex レビュー: boundary-conjugate-heat-transfer (plan)

- **plan**: [`plans/active/boundary-conjugate-heat-transfer.md`](../../plans/active/boundary-conjugate-heat-transfer.md)
- **stage**: `plan`
- **date**: 2026-09-19
- **commit**: `f8fca224` (feature/sern-design)
- **codex**: effort `high`, 6.4 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M5/m1
- **focus**: 3 巡目。2 巡目 (notes/reviews/2026-09-19-boundary-conjugate-heat-transfer-plan-3.md, NO-GO, C0/M6/m1) を全件採用して改訂した。まず 2 巡目の #1-#7 が実際に解消されたかを 1 件ずつ検証してほしい: #1 依存診断 (tooling-energy-balance-diagnostics) の dual-time/周期の解除マイルストーンと過渡の拘束反力式 C=Dt(VE)-Rraw、#2 D_f を初期推定値へ格下げ + 未緩和残差による受理判定と退避 (交番モード含む) が反例に耐えるか、#3 fem2d の (K_s+E^T D_f E)u=b_s+E^T[Q_f+D_f E u] と荷重転送・共有角 1 回計上、#4 陰解法フックを advanceImplicitSteady/implicitNonlinearUpdate 経路にしたこと、#5 G-cons/G-if の数値仕様 (sum|Q| と絶対床での規格化、局所ノルム、連続反復数、準定常閾値を比較許容の 1/5) と Phase 1 前倒し、#6 V5 の 3 段分解 (実測壁温→流体 / 固体単独 / CHT) と不確かさ項目、#7 索引・親 plan の同期。残る穴があれば挙げ、実装に進んでよいかを判定すること。判定が GO でないなら、あと何を確定させれば GO になるかを具体的に列挙してほしい。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
固定点保存式・界面熱量の符号・陰解法フックは妥当です。  
ただし、受理判定には再試行から抜けられない反例があり、診断・共有角・合格基準にも未確定事項があります。以下を実装前に仕様へ反映してください。

ファイルは変更していません。コード照合、CPU 上の数値検算、既存2 run の収束再判定を実施しました。新規 CFD 計算は実行していません。

前回 #1–#7 の再評価は次のとおりです。「解消」は計画上の評価です。

| 前回 | 再評価 |
|---|---|
| #1 診断依存・過渡反力 | **部分解消**。`C=Dt(VE)−Rraw` と解除順序は正しい。ただし提供側 plan に解除作業がなく、採取契約も未確定。今回 #2。 |
| #2 `D_f`・受理判定 | **未解消**。上界主張の撤回は適切。しかし最大ノルムによる受理と一律倍増では停止する反例がある。今回 #1。 |
| #3 `fem2d`・荷重転送 | **数式・単位・試験追加は解消**。`Eᵀ` による積分荷重転送は妥当。共有角の所有規則には穴が残る。今回 #4。 |
| #4 陰解法フック | **解消**。[`main.cpp:1753`](/home/sano/work/forge/solver_density_cuda/main.cpp:1753) は実際に `implicitNonlinearUpdate` を呼ぶ。同一状態の採取と次回組立て前の更新も適切。 |
| #5 ゲート | **部分解消**。絶対値和による規格化と Phase 1 前倒しは適切。振動許容と具体的な判定仕様が不足。今回 #3。 |
| #6 V5 | **分解方針・不確かさ項目は解消**。ただし (a)(b) の入力・比較対象・合格条件が未定義。今回 #5。 |
| #7 索引・親 plan | **指定箇所は解消**。[索引:28](/home/sano/work/forge/plans/README.md:28) と[親 plan:165](/home/sano/work/forge/plans/active/tooling-nozzle-isothermal-wall-chain.md:165) は修正済み。周辺の旧記述は残る。今回 #6。 |

目的は妥当です。`plans/accepted/` に完成済み CHT の重複実装は見当たらず、親 plan からの移管も確認しました。外部連成から始め、node の既存等温壁経路を使う方針は支持します。検証を node 主体とすることも [`procedures/verification/README.md`](/home/sano/work/forge/procedures/verification/README.md) と整合しています。

**指摘一覧**

1. **Major — 残差最大ノルムの単調減少を要求すると、安定な伝導系でも全候補を棄却します。**

   **根拠:** [plan:103](/home/sano/work/forge/plans/active/boundary-conjugate-heat-transfer.md:103) の受理判定を、線形問題で検算しました。以下は forge の実測ではなく、計画の更新式に対する CPU 検算です。

   \[
   A_s=0.1I,\quad
   H=-Q_f'=\begin{pmatrix}2&-1\\-1&2\end{pmatrix},\quad
   D_f=s\,\mathrm{diag}(1,10),\quad r=(1,1)^T.
   \]

   `A_s` と `H` はともに対称正定値です。`H` は各節点から固定温度へのコンダクタンス1、節点間コンダクタンス1の伝導網に相当します。

   試行後の残差は
   \[
   r_{\rm trial}=[I-(A_s+H)(A_s+D_f)^{-1}]r
   \]
   となります。

   | 倍率 `s` | 試行後の最大ノルム（試行前は1） |
   |---:|---:|
   | 1 | 1.701170 |
   | 2 | 1.371713 |
   | 4 | 1.191533 |
   | 1024 | 1.000771 |

   第2成分は
   \[
   1+\frac1{0.1+s}-\frac{2.1}{0.1+10s}>1\qquad(s\ge1)
   \]
   なので、**何度倍増しても受理できません**。一方、`s=1` の固定点反復のスペクトル半径は **0.960450 < 1** です。発散する反復を止める問題ではなく、収束可能な反復を受理規則が止めています。交番モードだけの試験では検出できません。

   **対案:** 局所最大ノルムは最終合格ゲートに残し、反復の受理には、降下方向を確認する残差メリット関数と line search を定義してください。伝導 SPD 問題では、比較中に重みを固定した
   \(\Phi=r^T(A_s+D_f)^{-1}r\)
   が候補です。本反例では `s=1` で **1.008101 → 0.883107** と減少します。一般の流体応答には SPD 保証がないため、降下しない場合の方向変更、再試行上限、量子化による停滞の扱いも必要です。**停滞を合格にしないこと**を明記し、V1b に非一様 `D_f` と本反例を追加してください。

2. **Major — 診断の解除マイルストーンが、提供側の作業と採取仕様に接続していません。**

   **根拠:** [CHT plan:141](/home/sano/work/forge/plans/active/boundary-conjugate-heat-transfer.md:141) は解除順序を記載しましたが、[診断 plan:32](/home/sano/work/forge/plans/active/tooling-energy-balance-diagnostics.md:32) は引き続き dual-time・周期を拒否し、[残作業表:148](/home/sano/work/forge/plans/active/tooling-energy-balance-diagnostics.md:148) にも対応する解除項目がありません。

   また、`Rraw` を空間残差とするか BDF 込みとするかは、属性だけでなく計算式を左右します。実装は [`main.cpp:1814`](/home/sano/work/forge/solver_density_cuda/main.cpp:1814) で空間残差を組み、その後 BDF 項を加えます。壁エネルギー残差はそれ以前に射影されるため、後から `res_roe` を読むだけでは復元できません。

   **対案:** 提供側にも次を登録してください。

   - `Rraw` は壁残差射影前の空間残差と固定する。
   - 同一評価状態について、実際の BDF 係数・履歴・体積から `Dt(VE)` を計算し、`C=Dt(VE)−Rraw` とする。
   - 非零の蓄積を持つ試験で、BDF1/BDF2 の符号と単位を検証する。
   - 周期は root 単位の一回集計を採用し、部分 CV・壁反力との対応を出力する。

   初版診断の範囲を維持したまま、後続マイルストーンとして追加して構いません。また、[V6:366](/home/sano/work/forge/plans/active/boundary-conjugate-heat-transfer.md:366) の「壁関数併用」は、§4.3 の対象外宣言に合わせて**解除後の追加試験**へ分離してください。

3. **Major — 準定常判定はまだ比較精度を保証せず、G-if も再現可能な合格式になっていません。**

   **根拠:** [plan:374](/home/sano/work/forge/plans/active/boundary-conjugate-heat-transfer.md:374) が厳しくするのは `drift` だけです。実装は [`check_quasisteady.py:293`](/home/sano/work/forge/solver_density_cuda/tools/check_quasisteady.py:293) で drift と振動幅を別々に判定します。

   同ツールの判定関数に、末尾が `[99,101,101,99]` の系列を与えた結果は次のとおりでした。

   ```text
   drift=0.0004, osc=0.10:
   STEADY — tail mean=100, drift=0.0%, fluct=2.0%

   drift=0.0004, osc=0.0004:
   OSCILLATING
   ```

   **0.2%比較のため drift を0.04%にしても、±1%の振動が合格します。**

   また、[G-if:368](/home/sano/work/forge/plans/active/boundary-conjugate-heat-transfer.md:368) は「事前登録する」とだけあり、相対残差の分母、絶対許容との論理関係、連続回数が未定義です。積分荷重 `[W]` の絶対閾値だけでは、細分化で節点荷重が小さくなる効果も混入します。

   **対案:** `drift` **と** `osc` の両方を比較許容の1/5以下にし、比較する局所温度・局所熱量にも適用してください。G-if は局所面積に対応した尺度、絶対・相対条件の組合せ、必要連続回数、欠損・非有限・量子化停滞時の不合格を数式化してください。少なくとも V1 の具体的な値を実装仕様の基準として固定する必要があります。

4. **Major — 共有角で、非連成の等温壁による上書きを防げません。**

   **根拠:** [plan:187](/home/sano/work/forge/plans/active/boundary-conjugate-heat-transfer.md:187) が拒否するのは「異なる `conjugateGroup` の共有ノード」です。しかし共有相手が通常の `wall_isothermal` なら、この条件では捕捉できません。

   コードは [`nodeWallDirichlet_d.cu:154`](/home/sano/work/forge/solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu:154) で各 bcond の温度ピンを順番に適用します。角ノードは複数 bcond に重複します（[`gmshReader.hpp:2292`](/home/sano/work/forge/solver_density_cuda/mesh/gmshReader.hpp:2292)）。例えば連成側800 K、非連成側300 Kなら、最後に適用した側が勝ちます。

   **対案:** 初版は、共有 CV に接する**全温度拘束**を走査し、非連成等温壁との共有を起動時に拒否してください。同一グループ内は global CV ID ごとに固体 DOF を一意に対応させ、
   \[
   Q_j=\sum_{\text{当該CVの壁面}}F^E-C_j
   \]
   を一度だけ転送すると明記します。正常な共有角保存試験に加え、連成壁と非連成等温壁の競合試験が必要です。

5. **Major — V5 の (a)(b) は、まだ「通った」を判断できません。**

   **根拠:** [V5:365](/home/sano/work/forge/plans/active/boundary-conjugate-heat-transfer.md:365) は (a)(b) の合格を前提にしますが、(a) の `h` の定義・比較点・許容、(b) の外周境界条件・比較量がありません。冷却孔の Robin 条件だけでは、固体単独試験の外周条件が閉じません。

   **対案:** 同じ翼・同じ条件で、以下に固定することを勧めます。

   - **(a)** 実測 `Tw` を課し、原典と同じ温度基準で定義した `h`／熱流束と壁圧を比較する。
   - **(b)** 実測 `Tw` を外周 Dirichlet、公開冷却条件を孔 Robin として解き、外周反力熱流束を比較する。これは**原典のデータ処理の再現検査**として位置づけ、入力した壁温への一致を成果にしない。
   - **(c)** 壁温を未知に戻して CHT を解く。

   各段の比較点・ノルム・許容帯の算定方法と、不確かさの合成方法を仕様化してください。「別々に積む」だけでは、単純和か統計的合成か、相関をどう扱うかが決まりません。

   今回は [NASA CR-168015 の PDF](https://ntrs.nasa.gov/api/citations/19830020105/downloads/19830020105.pdf) の直接取得が403で、表原画像の再確認はできていません。原典数値を新たに検証済みとはしていません。この指摘は plan 自体の試験定義不足に基づきます。

6. **Minor — 同期済みの記述と、まだ残る旧方針を整理してください。**

   **根拠:** [対象 plan:297](/home/sano/work/forge/plans/active/boundary-conjugate-heat-transfer.md:297) は索引・移管表を未修正としていますが、その指定箇所は修正済みです。一方、[親 plan:138](/home/sano/work/forge/plans/active/tooling-nozzle-isothermal-wall-chain.md:138) は「弱ループが収束しない＝軸方向伝導が支配的」という旧判断を残しています。反復不収束からモデル不足は判定できません。

   **対案:** 指定箇所の同期を完了扱いにし、親 plan の判断基準を本計画 §4.8 のモデル感度・厚さ方向近似の評価へ置き換えてください。残作業表も、追加した重要項目 #20–#26 を実際の依存順に並べ直すべきです。

既存 run の確認結果も記します。各 run の本段 `residual_history.csv` 全体を現行 `check_convergence.py` で再判定し、保存済み `CONVERGENCE_VERDICT.txt` と結論は同じでした。

| run | VERDICT |
|---|---|
| `case/48.flat_plate_cooled_m4/run_0011_Bplain_tw300_y3/` | `NOT CONVERGED (stalled/plateau)` |
| `case/44.vitiated_air_wt/run_0114_va_ns_fine_iso300_samewall_cont/` | `NOT CONVERGED (stalled/plateau)`。`rms_roUy`・`rms_roOmega` は `RISING` |

**推奨は、A→B の方針を維持し、上記 #1→#5 の順に契約を確定してから Phase 0/1 の実装へ進むことです。** GO に必要なのは、①停止しない受理・退避仕様、②提供側と一致する診断解除契約、③再現可能なゲート、④共有角の唯一の所有者、⑤V5各段の入力と合格条件です。Phase 2 は引き続き V1–V3 合格後としてください。

レビュー提案は **plan 未反映**です。

指摘数: Critical 0 / Major 5 / Minor 1
