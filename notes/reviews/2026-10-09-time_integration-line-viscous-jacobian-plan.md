# codex レビュー: time_integration-line-viscous-jacobian (plan)

- **plan**: [`plans/active/time_integration-line-viscous-jacobian.md`](../../plans/active/time_integration-line-viscous-jacobian.md)
- **stage**: `plan`
- **date**: 2026-10-09
- **commit**: `dc3bf1cc` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 7.0 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M5/m2

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**GO-with-changes**  
粘性・熱伝導でライン内の結合を補う方針は妥当で、薄層モデル内の D・K の符号も整合しています。  
ただし、壁拘束との不整合、実装に必要な入力経路、検証・採用基準を修正してから着手すべきです。

目的は既存計画と完全には重複しません。[既存 v2](/home/sano/work/forge-integ-1005/plans/accepted/time_integration-line-implicit-viscous-v2.md:21) はスカラー結合であり、今回の速度・温度・粘性仕事のブロック結合は未実装です。既存の `storeLU` → block-Thomas に載せる構造も妥当です。node 限定は現行運用と整合し、cell の追加計算を要求する必要はありません。

確認の限界として、参照 run の本体はローカルにありません。`check_convergence.py` を `run_0217`・`run_0223`・`run_0224` に実行した結果はすべて `NO residual_history.csv` でした。以下の時系列の再判定は、保存済みの抽出 JSON に対するものです。ファイルは変更していません。

1. **Major — 等温壁の ΔT=0 という前提が、採用予定のキー 5 では成立しない。**

   **根拠:** [計画:52](/home/sano/work/forge-integ-1005/plans/active/time_integration-line-viscous-jacobian.md:52) は壁側の熱伝導 K を ΔT_w=0 として消します。一方、[実装:1091](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1091) は、ビット 2 が無ければ Δ(ρE)_w=0 とし、Δρ_w は残します。キー 5 にはビット 2 がありません。

   静止壁では ΔT_w=−e_wΔρ_w/(ρ_wc_v)。CPG・300 K・Δρ/ρ=1% の例なら、線形補正は **−3 K** です。更新後の温度ピンで戻しても、線形系内部の拘束は一致しません。壁側 K の消去と壁行の単位行化を、そのまま一貫した Dirichlet 線形化とは呼べません。

   **対案:** mode 2 では強制等温壁に **Δ(ρE)_w−e_wΔρ_w=0** を課し、その拘束下で壁側の粘性 K を消してください。ビット 2 は面ごとの熱伝導処理から分離し、「ライン外の面にだけ効く」という記述の対象外にします。壁密度を動かす単体試験で ΔT_w=0 を確認することが先です。

2. **Major — mode 2 単独で動かすための入力経路と対応範囲が不足している。**

   **根拠:** [実装:1562](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1562) は `implicitThermalJacobian & 1` が偽なら `thermCond`・`cp`・`fx` に `nullptr` を渡します。計画どおり分岐だけを追加すると、`lineViscCoupling: 2` 単独では必要な配列が届きません。また、現カーネルの層流粘性入力は `cfg.visc` ですが、残差は [面平均の `vis_lam`](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/viscousFlux_d.cu:164) を使います。

   さらに、[計画:66](/home/sano/work/forge-integ-1005/plans/active/time_integration-line-viscous-jacobian.md:66) の設定検査は `heatCorrSU2` を制限しません。しかし [残差:274](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/viscousFlux_d.cu:274) の非ゼロ経路では、熱伝導の係数が計画の δ/dcc と異なります。壁関数も [残差:224](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/viscousFlux_d.cu:224) で traction を変更します。

   **対案:** §4.3 に wrapper・引数の変更まで記載し、mode 2 単独でも物性配列を渡してください。初版は導出済みの **`heatCorrSU2: 0`・低 Re・強制壁処理**に限定し、未対応の組合せは明示的に拒否します。mode 2＋thermal キー 0 の起動試験を必須にしてください。

3. **Major — U2 とビルド確認だけでは、今回追加する Jacobian を検証できない。**

   **根拠:** [U2](/home/sano/work/forge-integ-1005/plans/active/time_integration-line-viscous-jacobian.md:87) は静止した純伝導です。せん断結合、高速流の温度交差微分、粘性仕事の誤実装を十分に検出できません。FP32 は [§5:72](/home/sano/work/forge-integ-1005/plans/active/time_integration-line-viscous-jacobian.md:72) でコンパイルしか要求していません。Thomas が内部 double でも、D・K は [入力段階で `flow_float` に保存](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1824) されます。

   今回、計画の薄層流束を独立に実装して 50 組の状態で中心差分と照合したところ、D=−∂R/∂Q_i、K=∂R/∂Q_j の列ごとの相対誤差は最大約 **1.7×10⁻⁷**でした。これは式の確認であり、CUDA 実装の保証にはなりません。

   **対案:** U2 の前に、実装する共通関数を対象とした有限差分試験を置いてください。CPG/TP、高速接線流、両面方向、壁拘束を含め、組み上げた短いラインを密行列解と照合します。続いてせん断試験と FP32 実行を追加し、線形解の後退誤差・factor 失敗数・正値性を確認します。U0 も、維持を約束している **mode 1** を含めてください。

4. **Major — 安定性・収束・派生量の定常性を判定する手順が未確定。**

   **根拠:** [V-n1:90](/home/sano/work/forge-integ-1005/plans/active/time_integration-line-viscous-jacobian.md:90) の「NaN なし・開始の 10 倍以内」は短期の有界性検査です。末尾で成長中のケースや、物性床に張り付いたケースも通り得ます。§6 には判定ツールの実行区間・閾値がありません。

   保存済み JSON の末尾 2 万 step を `check_quasisteady.py` の系列入力へ渡し、`--tail 1 --drift 0.0005`、既定の最低 4 点で再判定しました。

   | 対象 run（`case/45.isobutane_m6_d155/` 配下） | 区間 | VERDICT |
   |---|---:|---|
   | `run_0217_ns_coldmesh_tw300_cfl4_ext3` | 180000–200000 | **TRANSIENT-UNSETTLED**：3 点しかない |
   | `run_0224_ns_coldmesh_tw300_linedir_tj5_cap50_ext` | 40000–60000 | **DRIFTING**：θ_r・Q_w とも |

   `run_0217` の判定は「未定常を証明」ではなく、**この窓では確認不足**です。参照値 0.11477 を確定した定常解として扱う根拠は揃っていません。run 索引は [case README:141](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/README.md:141) にあります。

   **対案:** 短期安定性、収束判定、採用判定を分けて事前登録してください。全残差の末尾成長、ρ/P/T の正値性、床・クリップ作動、`check_convergence` の区間と VERDICT、θ_r・Q_w の `check_quasisteady` を明示します。欠損は符号付き和だけでなく **絶対値と局所残差ノルム**も使い、相殺を除外してください。

5. **Major — V-n2 は「高速化していない変更」を合格にできる。**

   **根拠:** [V-n2:92](/home/sano/work/forge-integ-1005/plans/active/time_integration-line-viscous-jacobian.md:92) は、同じ step の欠損が同等で、θ_r がわずかに参照値へ近づき、ms/step が 10% 増えても合格です。到達 step 数がほぼ変わらなければ、壁時計では遅くなります。

   また、比較対象の `run_0223` は旧バイナリです。同時進行の [速度計画:73](/home/sano/work/forge-integ-1005/plans/active/time_integration-line-implicit-speed.md:73) が同じカーネルを変更するため、過去の ms/step との比較では今回の費用を分離できません。

   **対案:** 同一バイナリ・同一 IC・専有 GPU・native 環境で mode 0/2 を再比較し、主指標を **事前登録した同じ状態水準までの壁時計**にします。固定 step の比較は途中診断、10% は費用上限として扱ってください。まず R=50 で評価し、その後に上限なしを試す順序を推奨します。

6. **Minor — P=I＋⅓nnᵀ は、現残差の「勾配を凍結した厳密微分」ではない。**

   **根拠:** [残差:205](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/viscousFlux_d.cu:205) では直接の速度差項は βΔu で、転置・発散は面勾配から入ります。その勾配を固定すると運動量の直接微分は βI です。[計画:43](/home/sano/work/forge-integ-1005/plans/active/time_integration-line-viscous-jacobian.md:43) の P は、それらの法線寄与を薄層モデルとして近似的に戻しています。

   **対案:** P は維持してよいですが、「現離散残差の厳密 Jacobian」ではなく **薄層近似の前処理行列**と明記してください。単体試験では薄層モデルとの微分一致と、実残差との差を別々に測ります。

7. **Minor — 「run_0222 で SST の寄与を確認」は断定が強い。**

   **根拠:** [関連計画:199](/home/sano/work/forge-integ-1005/plans/active/time_integration-implicit-thermal-jacobian.md:199) 自身が、凍結 run の `roK` に最大 **2.9%** の変化があり、壁・入口ピンの除外は未確認と記録しています。事前に要求した凍結の実効性確認が完了していません。

   **対案:** §4.2 は「SST 更新を止めた試行で増幅が弱まったが、寄与の分離は未完了」に修正してください。凍結対象外を明示した `roK`・`roOmega` と μ_t の監査が終わるまで、原因を確定扱いしないことです。

推奨は、**壁拘束と対応範囲を確定した mode 2 を opt-in で実装し、局所 Jacobian／ライン解の検証 → U2・せん断・FP32 → R=50 の同一条件 A/B → 上限なし**の一本です。実装前の修正優先順は **1 → 2 → 3 → 4 → 5**。方針は撤回不要ですが、現状の §6 では成功も失敗も原因を判別しきれません。以上はレビュー提案であり、read-only 指定に従い plan 未反映です。

指摘数: Critical 0 / Major 5 / Minor 2
