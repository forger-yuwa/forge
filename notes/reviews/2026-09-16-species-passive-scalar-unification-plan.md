# codex レビュー: species-passive-scalar-unification (plan)

- **plan**: [`plans/active/species-passive-scalar-unification.md`](../../plans/active/species-passive-scalar-unification.md)
- **stage**: `plan`
- **date**: 2026-09-16
- **commit**: `e7ae9fda` (feature/sern-design)
- **codex**: effort `high`, 5.9 min, rc=0
- **判定**: **NO-GO**, 指摘 C1/M7/m1
- **focus**: 設計方針 §4 (受動種を化学種カーネルに乗せる §4.1、S3 [2 次 face 移流] の本番化と LHS/RHS defect-correction 不整合の解消策 §4.2 [implicitRelax・deferred-correction 係数・scalarCflMax]、凝縮モーメントの移流置換 §4.3、dual-time の化学種修正移植と受動種 BDF §4.4) と検証計画 §6 の妥当性・見落とし。特に: (1) deferred-correction 分割で固定点が変わらない主張と擬似時間安定性の根拠、(2) 受動種の面値クリップ (非正規化) の保存性、(3) 化学種の再正規化・トレーサのクランプが同一 run 内の差 1.8e-4 の原因という仮説 (§4.0) の妥当性、(4) dual-time 移植で chem ブランチ固有部分を除外する境界、(5) node 周期・ピン・境界の受動種対応の漏れ、(6) 検証ゲートの実効性 (case/28 の S3 発散再現、case/44 cfl 6 + relax 0.7)。cell は対象外 (ユーザ指示)。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

# 判定: **NO-GO**

受動スカラの輸送共通化と dual-time 修正は必要です。しかし、§4.2 の固定点不変の主張は数式上誤りで、S3 発散の原因説明も実際の更新方式と一致しません。周期処理・有界性・検証ゲートを修正してから実装すべきです。

対象 HEAD: `1784924a`。コード、移植元 `e296f0d0`、保存済み run を確認し、収束・準定常判定を再実行しました。ファイル変更・新規 CFD 計算はしていません。

## 指摘一覧

### 1. **Critical — deferred-correction の記載どおりに実装すると固定点が変わる**

**根拠:** [plan:91](/home/sano/work/forge/plans/active/species-passive-scalar-unification.md:91)。

記載された残差は

\[
R_\omega=R_1+\omega_{dc}(R_2-R_1)
\]

です。この更新が停止する条件は **\(R_\omega=0\)** であり、一般に \(R_2=0\) ではありません。例えば \(\omega_{dc}=0.5\) では、一次・二次残差を半分ずつ混ぜた別の離散方程式を解きます。「定常で差分項が全量残る」ための履歴変数・反復式は計画にありません。

また、一次 LHS・二次 RHS という事実だけでは発散は証明できません。必要なのは実際の反復作用素の安定性です。残差を全量保持しても、緩和係数 0.7 が安定性を保証するわけではありません。

**対案:** **RHS は常に完全な \(R_2\) とし、緩和は増分と擬似時間刻みに限定する**方針へ変更してください。`speciesDeferredCorrection: 0.5` は撤回する。凍結流れの小型問題で反復の減衰を確認し、収束した \(R_2\) の解が緩和・擬似 CFL に依存しないことを検証します。

### 2. **Major — case/28 の発散例は「未緩和の point-implicit」ではない**

**根拠:** 実在する発散 run は

`case/28.cutler_coaxial_jet/run_0052_rycl_B_s3_cfl4/`

です。[solverConfig.yaml:33](/home/sano/work/forge/case/28.cutler_coaxial_jet/run_0052_rycl_B_s3_cfl4/solverConfig.yaml:33) に、すでに以下があります。

- `implicitRelax: 0.7`
- `speciesImplicitCoupling: 1`
- `speciesFaceReconstruction: 2`

この更新は近傍補正を持つ scalar-DPLUR で、緩和も適用済みです。[speciesTransport_d.cu:339](/home/sano/work/forge/solver_density_cuda/cuda_forge/speciesTransport_d.cu:339)、[同:829](/home/sano/work/forge/solver_density_cuda/cuda_forge/speciesTransport_d.cu:829)。

再確認では **最初の NaN/Inf は step 433**、`check_convergence.py` は **`DIVERGED (NaN/Inf)`**。さらに、この config は `discretization` 未指定であり、読込時の既定は `cell` です。[solverConfig.cpp:185](/home/sano/work/forge/solver_density_cuda/input/solverConfig.cpp:185)。

**対案:** node 上で、同一 IC・BC に対する S2/S3 と coupling 0/1 の比較を**実装前の原因確認**へ移してください。過去の cell 結果は参考証拠として扱い、node の発散再現と原因を同一視しない。今回 cell の追加検証は要求しません。

### 3. **Major — 面クリップの保存性を誤解しており、有界性の保証もない**

**根拠:** [plan:102](/home/sano/work/forge/plans/active/species-passive-scalar-unification.md:102) の「面クリップで局所非保存」は誤りです。[speciesTransport_d.cu:565](/home/sano/work/forge/solver_density_cuda/cuda_forge/speciesTransport_d.cu:565) は、同一面で作った一つの流束を隣接 CV に逆符号で加えます。**共有面値をクリップしても、内部面の保存性は丸めを除いて維持されます。** 正規化の有無とは別問題です。

一方、面値の非負性は更新後の非負性を保証しません。現行 MUSCL は `phiC + limiter·grad·dx`。[convectiveFlux_common_d.cuh:80](/home/sano/work/forge/solver_density_cuda/cuda_forge/convection/convectiveFlux_common_d.cuh:80)。

例えば一様密度・正方向流れで、ステップ直前の値が `…0,0,1,1…`、`limiter_ro=1` なら、最後のゼロ点から出る面値は中央勾配により `0.25` になり得ます。流入面値はゼロなので、その点の増分は負です。面を `[0,1]` にクリップしても防げません。現行更新の floor は、この負値を切って保存量を増やします。[scalarTransport_d.cu:287](/home/sano/work/forge/solver_density_cuda/cuda_forge/scalarTransport_d.cu:287)。

また `Xi` は輸送に読む原始量なので、「表示だけクリップ」ではありません。

**対案:** 保存性の説明を訂正し、**面再構成の制限と保存量更新の補正を別々に診断**してください。保存的な流束制限を設計する必要があります。「新規リミッタなし」で有界性まで保証する条件は撤回すべきです。試験では `Xi` 表示値ではなく生の `roXi/ro`、総トレーサ量、floor 補正量を判定し、モーメントは非負性に加えて相互の実現可能性も確認してください。

### 4. **Major — node 周期は名前ベース登録だけでは成立しない**

**根拠:** [periodicNode_d.cu:84](/home/sano/work/forge/solver_density_cuda/cuda_forge/periodicNode_d.cu:84) の追加 gather は `res_*` だけです。`transport_diag_*` は集約されません。勾配 gather も固定リストで、化学種・トレーサ・凝縮モーメントの勾配は含まれていません。[同:152](/home/sano/work/forge/solver_density_cuda/cuda_forge/periodicNode_d.cu:152)。

再利用予定の `species_gradient_d` は全 plane を走査し、周期半割面を除外する引数もありません。[speciesTransport_d.cu:585](/home/sano/work/forge/solver_density_cuda/cuda_forge/speciesTransport_d.cu:585)。既存の周期勾配 gather は、その除外を前提としています。

さらに保存量 mirror は平均流更新時に呼ばれる経路があり、その**後に更新する受動量**の同期まで保証しません。[main.cpp:1447](/home/sano/work/forge/solver_density_cuda/main.cpp:1447)。

**対案:** §4.1 に周期処理の順序を明記してください。

- 周期半割面を除いた勾配積算と gather。
- 輸送対角の集約、ソース Jacobian の体積と整合した集約。
- 受動量更新後の状態 mirror。
- 空間残差 gather 後に、合併体積を使って BDF 項を一度だけ追加。

node 周期移流・拡散と、一様な凝縮ソースの試験を必須にしてください。継ぎ目・辺・角で保存量と更新速度を確認します。

### 5. **Major — 1.8e-4 の原因は未確定で、共通化だけでは一致を保証できない**

**実測:** `case/46.sern_design/run_0104_species_regress_euler_lumped_tracer/` の保存量から再計算しました。

| 出力 step | max \(|roXi/ro-roY0/ro|\) | 平均 |
|---:|---:|---:|
| 0 | 0 | 0 |
| 500 | 2.33948e-4 | 4.01995e-7 |
| 6000 | 1.80900e-4 | 4.61209e-7 |

step 6000 の最大差点では `Xi=0.9998191`、`Y0=1.0`。差の存在は再現できますが、保存済み出力から補正前の再正規化・クランプ寄与は復元できません。収束判定も **`NOT CONVERGED (stalled/plateau)`** です。

**根拠:** 化学種は更新後に非負化・再正規化されます。[speciesTransport_d.cu:142](/home/sano/work/forge/solver_density_cuda/cuda_forge/speciesTransport_d.cu:142)。提案はこの処理を残し、受動種からは除外します。したがって輸送残差を共通化しても更新写像は異なります。coupling 1/2 なら更新方式も異なります。

**対案:** 「残る差の原因である」を**検証する仮説**に戻す。凍結した共通流れ・同一 IC/BC で、残差、更新前後、再正規化前後、クランプ前後を段階比較してください。一致ゲートは同じ更新方式・同じ拡散係数の制御試験に限定し、別々の流れを解く `full`/`lumped` 間に同じ機械精度一致を要求しないでください。

### 6. **Major — dual-time 移植の境界条件・更新分岐・時間精度ゲートが不足**

**根拠:** 移植元 `e296f0d0` を確認すると、以下が残ります。

- `species_add_unsteady_d` は全内部 CV に BDF を追加し、ピンを考慮しない。同 commit の `speciesTransport_d.cu:817`。
- 移植元の dual-time は coupling 2 以外を `speciesTimeIntegration_d_wrapper` に流し、**coupling 1 の DPLUR 分岐を持たない**。同 commit の `main.cpp:1391`。
- 現行の入口ピン残差除去は `assembleResidual` 内、BDF 追加はその後です。[main.cpp:1249](/home/sano/work/forge/solver_density_cuda/main.cpp:1249)、[同:1570](/home/sano/work/forge/solver_density_cuda/main.cpp:1570)。単純移植では BDF がピン残差を再導入します。
- `check_convergence.py` は物理 step の `outer_end` を読むため、**各物理 step 内のサブ反復収束を直接判定するものではありません**。[check_convergence.py:32](/home/sano/work/forge/solver_density_cuda/tools/check_convergence.py:32)。

CMC、反応 block LU、反応用 `dt_local_sp` を除外する判断は妥当ですが、それだけで移植範囲は確定しません。

**対案:** coupling 0/1/2 の処理順、BDF 追加後のピン拘束、更新後の入口値再適用、履歴初期化・restart を明記してください。物理時間項は `scalarCflMax` で変更しない。

時間精度試験は滑らかな解を使い、同一最終時刻で \(\Delta t,\Delta t/2,\Delta t/4\) の三水準比較にする。全化学種・受動量の BDF 込み残差をサブ反復ごとに確認し、反復数倍増で結果が変わらないことも必須です。凝縮は純移流だけでなく、ソースが作動する試験を追加してください。

### 7. **Major — S3 本番化のゲートを通らずに完了できる**

**根拠:** [plan:173](/home/sano/work/forge/plans/active/species-passive-scalar-unification.md:173) の完了条件は「検証 1–4」のみで、核心の **5: S3 安定性、6: 原因切り分け**が抜けています。残作業表にもステップ 2b がありません。

また `passiveScalarScheme: 1` だけでは、既定 `speciesFaceReconstruction: 0` のまま一次輸送です。§6.2 の A/B には S3 を明示していないため、二次化を試さず通る余地があります。

今回の再判定は次のとおりです。

| run | `check_convergence.py` | `check_quasisteady.py` |
|---|---|---|
| `case/44.vitiated_air_wt/run_0170_va3_M4.19_Lc8_noneq_inletTt_lim1e_cfl2/` | **NOT CONVERGED** | `cond_series.csv`: **ALL STEADY** |
| `case/44.vitiated_air_wt/run_0183_va3_M4.19_Lc8_noneq_inletTt_lim1e_cfl6_relax07_from0170/` | **NOT CONVERGED** | `cond_series.csv`: **ALL STEADY** |

準定常判定は指定六量、`--drift 0.002 --osc 0.002`。この結果は量の安定を示しますが、**離散方程式の固定点に到達した証拠ではありません**。run 索引は [case/44 README:708](/home/sano/work/forge/case/44.vitiated_air_wt/README.md:708)。

**対案:** 検証 1–6 を完了条件に入れ、全比較の SFR・coupling・緩和・実効 scalar CFL を固定してください。固定点不変は小型の `PASS` ケースで検証し、case/44・46 は準定常回帰として別ゲートにする。「完走」「S2 と同程度の残差床」だけでは本番化を認めないでください。

### 8. **Major — F-sp1 を閉じるための拡散検証がない**

**根拠:** [accepted plan:231](/home/sano/work/forge/plans/accepted/thermophysics-cea-mole-fraction-species.md:231) は、**粘性二流体・差動拡散試験を F-sp1 に延期**しています。今回の計画は F-sp1 を閉じますが、主要なトレーサ比較は Euler で、追加する `passive_diffusion_d` が実行されません。

化学種拡散は定数 Sc と混合平均拡散を区別し、さらに \(\Sigma J=0\) 補正を行います。[speciesTransport_d.cu:248](/home/sano/work/forge/solver_density_cuda/cuda_forge/speciesTransport_d.cu:248)。固定 Sc のトレーサと差動拡散する化学種は、一般には同じ輸送方程式ではありません。

**対案:** node の解析解付き拡散試験を追加し、層流 `Sc`、乱流 `Sc_t`、無流束壁、入口、周期を検証してください。`nSpecies==1` でも動作させ、モーメントには拡散が入らないことも確認する。化学種との一致は等拡散係数条件に限定し、差動拡散では非一致を正しい挙動として評価します。

### 9. **Minor — 前提訂正と設定仕様が本文全体へ反映されていない**

**根拠:**

- §1・§3 は依然「化学種は二次」「`limiter_Y` を使用」と説明し、§4.0・§4.1 と矛盾。[plan:20](/home/sano/work/forge/plans/active/species-passive-scalar-unification.md:20)。
- `implicitRelax` のコード既定は **1.0**。0.7 は推奨レシピです。[solverConfig.cpp:329](/home/sano/work/forge/solver_density_cuda/input/solverConfig.cpp:329)。
- 移植元の履歴カウンタはシフト後に増えるため、`g_speciesLevelShifts < 2` は通常**最初の一 step**を BDF1 にします。「最初の二 step」と異なります。
- 共通スカラ更新へ無条件に緩和を追加すると、`passiveScalarScheme: 0` の旧挙動維持と両立しません。化学種も別途影響を受けます。

**対案:** §1・§3・設定表・旧経路の契約を一括修正し、既定値、推奨値、適用対象を明確にしてください。

## 推奨

**完全な二次残差を保持し、共通化 → BDF 整備 → S3 本番化の順に段階実装する案に絞ります。**

実装前の修正優先順は次のとおりです。

1. §4.2 の残差ブレンドを撤回し、node で発散原因と更新方式を確認する。
2. 保存性・有界性・周期・ピン・BDF の処理順を設計へ追記する。
3. 一次経路で輸送共通化と差の原因を検証する。
4. 拡散・時間精度・S3 の必須ゲートを追加し、全ゲート通過後に既定変更する。

目的自体は未解決の F-sp1/F-cf8 と整合し、重複実装ではありません。ただし S3 は既存 accepted plan が本番採用を退けた機能です。今回の計画では、その判断を覆す証拠がまだ不足しています。

**plan 未反映** — read-only の依頼に従い、上記はレビュー提案として提示しています。

指摘数: Critical 1 / Major 7 / Minor 1
