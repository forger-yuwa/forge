# codex レビュー: condensation-kantrowitz-gamma-twophase-sonic (plan)

- **plan**: [`plans/active/condensation-kantrowitz-gamma-twophase-sonic.md`](../../plans/active/condensation-kantrowitz-gamma-twophase-sonic.md)
- **stage**: `plan`
- **date**: 2026-09-10
- **commit**: `5db69d93` (feature/sern-design)
- **codex**: effort `high`, 6.5 min, rc=0
- **判定**: **NO-GO**, 指摘 C1/M5/m1
- **focus**: §4.1 Kantrowitz の γ を凝縮種 (蒸気) の γ_v にする根拠と、§4.2 二相 frozen 音速 c²=γ_2φ R_eff T の導出 (c_p,2φ = c_p^全蒸気 − g L'(T), c_v,2φ = c_p,2φ − R_eff, R_eff = R_mix − g R_w) が一温度二相 EOS と整合しているか。per-cell gamma を γ_2φ に置換して block-DPLUR の κ=γ−1 と χ=c²−κh を組むことの妥当性。検証計画 (5 run, 単体 FD テスト) の十分性
- **extra**: `methods/condensation.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**  
蒸気 `γ_v` の採用と、固定組成・固定 `g` に対する二相音速の導出は妥当です。  
しかし、N2 で `c²<0` になる既存物性、Jacobian の固定変数の混同、検証の不足があり、全凝縮ケースへの既定適用は認められません。

1. **Critical — N2 では、提案式が負の音速二乗を生成する。**

   根拠は [`condensationProperties_d.cuh:67`](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationProperties_d.cuh:67) の `n2_latent`。コードの多項式と計画指定の中心差分を倍精度で評価すると、次になります。**run の実測ではなく、物性式の数値評価です。**

   | 入力・評価量 | 値 |
   |---|---:|
   | `T`, `g` | `45.2 K`, `0.2` |
   | `cp_allvap`, `R` | `1038.8`, `296.8 J/(kg·K)` |
   | `L′` | `+4159.75 J/(kg·K)` |
   | `cp_2φ` | `206.85 J/(kg·K)` |
   | `cv_2φ` | **`−30.59 J/(kg·K)`** |
   | `c²` | **`−72572.3 m²/s²`** |

   この温度は実装が掲げる物性有効域内で、`g=0.2` も plan が想定する値です。現行の温度反転は [`condensationEOS_d.cuh:81`](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationEOS_d.cuh:81) で Newton 分母を正値化しますが、これは EOS の熱力学的安定性を保証しません。そこから生の微分を音速に使えば `sqrt(c²)` が NaN になります。

   **対案:** N2 への既定適用を今回から外す。N2 は潜熱・内部エネルギーの整合した修正と、有効温度・湿り度範囲全体の `cv_2φ>0`, `c²>0` 検証を先行させる。音速だけを床で隠す修正は不可です。

2. **Major — `γ_2φ` の代入で「分離解法の Jacobian が厳密に整合する」という説明は成立しない。CPG の Jacobian も変わる。**

   固定 `g,Y` なら、

   \[
   de/dT=c_p^{全蒸気}-R_{\rm eff}-gL',\qquad
   \kappa=R_{\rm eff}/c_{v,2\phi}
   \]

   であり、§4.2 の音速式は正しいです。一方、実装は NS 保存量を更新した後で `rog` を別に更新します。根拠は [`main.cpp:1315`](/home/sano/work/forge/solver_density_cuda/main.cpp:1315)、[`main.cpp:1363`](/home/sano/work/forge/solver_density_cuda/main.cpp:1363)。**固定 `g` と固定 `ρg` は異なります。**

   固定組成について、

   \[
   \xi_g=\left.\frac{\partial p}{\partial(\rho g)}\right|_{\rho,\rho e}
        =-R_wT+\kappa(L-R_wT),
   \]
   \[
   \chi_{\rho g}=\underbrace{c_f^2-\kappa h}_{\chi_g}-g\xi_g .
   \]

   提案が与えるのは固定 `g` の近似ブロックであり、`ρg` 列を含む保存系の厳密な Jacobian ではありません。`condEquilibrium:2` は [`dependentVariables_d.cu:119`](/home/sano/work/forge/solver_density_cuda/cuda_forge/dependentVariables_d.cu:119) で `g` 自体を再決定するため、さらに別の微分です。

   また、CPG 分岐は [`block_dplur_jacobian_d.cuh:58`](/home/sano/work/forge/solver_density_cuda/cuda_forge/block_dplur_jacobian_d.cuh:58) で `c/(γ−1)` を使います。**`gamma` が一定でも `sonic` を変えれば Jacobian は変化します。** plan §4.2 の「変更なし」は誤りです。

   **対案:** 今回の LHS を「固定 `g,Y` の frozen 近似」と明記し、厳密整合の主張を削除する。実際の `accumulate_split_jacobian_cf` を使い、固定 `g` の流束 FD、固定 `ρg` との差、float32 での行列作用を検証する。CPG・平衡拘束形への展開は別途検証する。

3. **Major — `dependentVariables` 一箇所では、境界まで二相音速に統一できない。**

   [`boundaryCond_d.cu:612`](/home/sano/work/forge/solver_density_cuda/cuda_forge/boundaryCond_d.cu:612) の TP 静圧出口はセル `gamma` を特性構成に使いますが、[`同:660`](/home/sano/work/forge/solver_density_cuda/cuda_forge/boundaryCond_d.cu:660) では ghost の温度・エネルギー・音速を全蒸気 EOS から再構築します。提案後は内部側が二相、ghost 側が全蒸気になります。

   node 等温壁も [`nodeWallDirichlet_d.cu:78`](/home/sano/work/forge/solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu:78) で `P/roe/sonic` を全蒸気状態に上書きします。これらは既存の制約ですが、全凝縮ケースへの既定変更を保証できない直接の理由です。5 run の `outflow`・断熱壁では検出できません。

   **対案:** 境界による状態再構築を影響範囲表に追加し、二相状態を保持する境界試験を設ける。今回未検証の境界構成へ新モデルを自動適用しない。

4. **Major — 現行ツールを呼ぶだけでは、凝縮の収束と報告量の定常性を判定できない。**

   [`check_convergence.py:27`](/home/sano/work/forge/solver_density_cuda/tools/check_convergence.py:27) は検査列を NS 5 本と SST 2 本に固定し、**`rms_rog_*`, `rms_roQ*` を検査しません。**

   実際に再実行した根拠 run は  
   `case/16.nozzle_wys/run_0230_user_node_sst_cond_outflow/` です。

   - `check_convergence.py`: **`NOT CONVERGED (still converging — run more steps)`**
   - `rms_roe`: 最終 `6.87e-3`、ピークから **1.8 桁**低下
   - `rms_roOmega`: 最終 `3.54e2`、**2.5 桁**低下
   - `check_quasisteady.py --quantity pmax,machmax`: **両方 `STEADY`**

   `pmax` は全域最大値を取るだけです（[`check_quasisteady.py:96`](/home/sano/work/forge/solver_density_cuda/tools/check_quasisteady.py:96)）。この run の最終最大圧位置は `x=−37.57 mm` であり、凝縮 onset の監視になっていません。さらに既定の drift 許容値は **5%**。計画が評価したい壁圧差 `≲1%` より大きいです。

   **対案:** 凝縮・化学種の保存量残差を検査対象に追加する。onset、壁圧偏差、出口 `g/M/c`、`h0` 誤差そのものを時系列判定し、期待差より十分小さい許容幅を設定する。比較は残差判定と対象量の定常判定を満たすまで延長する。12000 step の同時刻比較は過渡比較として扱う。

5. **Major — 5 run は変更要因の分離には有効だが、既定変更の検証として不足する。**

   [`plan:136`](/home/sano/work/forge/plans/active/condensation-kantrowitz-gamma-twophase-sonic.md:136) の全 run が TP carrier・非平衡・node・同一境界です。cell、凝縮 OFF、pure TP、CPG、平衡拘束形を通りません。共有コードに対する node/cell 両方の検証は [`procedures/verification/README.md:39`](/home/sano/work/forge/procedures/verification/README.md:39) の必須事項です。

   単体試験も `c²` のスカラー FD だけでは、`gam_array` 配線、CPG 分岐、境界上書き、陰解法固有系を確認できません。`L′` の差分幅依存性、物性接続点、湿り度上限、float32 格納後の `gamma−1` 精度も未検証です。

   **対案:** 5 run を維持し、今回対応する構成について wet cell 回帰、独立した凝縮 OFF 回帰、実カーネル出力の照合を追加する。単体試験は温度・組成・湿り度を掃引し、FD 刻み依存性と物性接続点を確認する。未検証構成の既定変更は完了条件から外す。

   なお、既存 `run_0230` のメッシュ品質記録は **`VERDICT: PASS`**、最大 AR `724.2`、skewness `0.203`。IC の index コピー元も記録されており、この前提選択は妥当です。

6. **Major — 回帰条件と物理効果の合否条件が誤っている。**

   [`plan:148`](/home/sano/work/forge/plans/active/condensation-kantrowitz-gamma-twophase-sonic.md:148) の「新モデルでも dry セルの場が一致」は保証できません。**同一入力状態に対する dry EOS の一致**と、凝縮帯から圧力・粘性・陰解法を介して影響を受ける**計算場の一致**は別です。

   また、相対差 `≤1e-6` はビット同等ではありません。既存 [`case/34 README:43`](/home/sano/work/forge/case/34.arthur_n2_nozzle/README.md:43) にも同一設定の run 間差が記録されています。

   [`plan:152`](/home/sano/work/forge/plans/active/condensation-kantrowitz-gamma-twophase-sonic.md:152) の「onset ≲1 mm、壁圧 ≲1%という見積りと整合」を合否条件にするのも不適切です。これは未検証の仮説で、正しい実装が予想外の応答を示す可能性があります。

   **対案:** dry の厳密一致は同一状態を与える単体試験で判定する。run 回帰は同一バイナリの反復差を測り、ゼロ近傍も扱える絶対・相対許容値を定める。物理効果の見積りは観測項目に移し、合否は EOS 微分・保存性・有界性・収束性で決める。

7. **Minor — `γ_v` の選択は支持するが、carrier 中の物理誤差の説明を限定すべき。**

   純蒸気形の係数は

   \[
   \frac{2(\gamma_v-1)}{\gamma_v+1}
   =\frac{R_v}{c_{v,v}+R_v/2}
   \]

   と書け、蒸気の熱容量を使う説明は整合します。既存 NASA 係数を評価した値も、200 K で `γ=1.33212`、240 K で `1.33122`、300 K で `1.32888`。定数 `1.331` はこの温度範囲では妥当です（[`thermo_d.cu:104`](/home/sano/work/forge/solver_density_cuda/cuda_forge/thermo_d.cu:104)）。

   ただし、これは**純蒸気形近似内の修正**です。carrier による熱除去を含めたモデルの誤差が「15%」と確定するわけではありません。carrier 衝突を加えると非等温抑制が弱まることは一次研究にも示されています。[Wedekind et al., 式9–10](https://arxiv.org/pdf/0804.1516)

   **対案:** §1・§4.1 を「純蒸気形近似の係数修正」と明記する。上下限は **抑制量 `θ` の上下限**と明記し、核生成率では逆に `J_pure ≤ J_carrier ≤ J_iso` となることを示す。この条件付き関係から onset の実験誤差改善までは保証しない。

**推奨は、今回を「TP carrier・非平衡 H2O の限定導入」に縮小して再レビューすることです。** 優先順は、① N2・平衡拘束形・未対応境界への自動適用を外す、② frozen 近似 Jacobian の定義を修正する、③ 凝縮残差・対象量の判定を整える、④ node/cell と実カーネル試験を追加する、です。

課題自体は未解決で、既存 accepted plan と重複する完成済み機能ではありません。ただし、既存の loose coupling を今回だけで厳密化したとは扱えません。ファイルは変更しておらず、以上は **plan 未反映**です。

指摘数: Critical 1 / Major 5 / Minor 1
