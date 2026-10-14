# codex レビュー: time_integration-implicit-thermal-jacobian (plan)

- **plan**: [`plans/active/time_integration-implicit-thermal-jacobian.md`](../../plans/active/time_integration-implicit-thermal-jacobian.md)
- **stage**: `plan`
- **date**: 2026-10-09
- **commit**: `9f8eaf04` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 4.6 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M4/m1
- **extra**: `notes/reviews/2026-10-09-linedir-divergence-2-diagnose.md`, `solver_density_cuda/cuda_forge/timeIntegration_d.cu`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
温度微分をエネルギー行へ入れる方針は妥当ですが、提案する係数は対象ケースの実際の熱伝導率と整合していません。  
係数・適用範囲・検証条件を修正してから、対角ブロックだけを変える限定試験に進むことを推奨します。

指定ファイルを全文確認し、添付された plan・コードが作業ツリーと一致することも確認しました。ただし、対象の `case/45.isobutane_m6_d155/run_0181*・run_0183*・run_0191*・run_0203〜0211*` はこの作業ツリーにありません。発散ステップや帳簿の数値は提供記録として扱い、残差・HDF5・VERDICTを直接再検証した結果とは区別します。

1. **[Major] 提案する Λᵀ は、対象ケースの熱伝導残差を凍結線形化した係数になっていません。**

   根拠: [plan:59](/home/sano/work/forge-integ-1005/plans/active/time_integration-implicit-thermal-jacobian.md:59) は定数 `physProp.visc` と `prandtlLam` から係数を作ります。一方、[対象問題:33](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/problem_d155_ns_prod_coldmesh_tw300.yaml:33) は種別輸送物性を指定し、[生成器:152](/home/sano/work/forge-integ-1005/design/forge_design/evaluate/runner_axismach.py:152) は `viscMethod: 2` に変更して `thermCondMethod: 1` を除去します。[gasProperties_d.cu:107](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/gasProperties_d.cu:107) の λ は輸送モデルから独立に求まり、μcp/Pr ではありません。

   さらに、実際の内部面係数は [viscousFlux_d.cu:264](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/viscousFlux_d.cu:264) の **面補間した層流熱伝導率＋面cp×面μt/Prt** です。セル値だけの提案式とは、定Prの場合も一般には一致しません。この状態でV1が失敗しても、温度Jacobianの仮説を棄却できません。

   **対案:** 運動量・連続行の従来係数は維持し、熱伝導行には **Λᵀ = (k_eff,f / c_v,i)·δ/dcc、c_v,i = cp_i/γ_i** を使ってください。既存の `thermCond`・`cp`・`vis_turb`・面補間重みを参照し、物性値と勾配補正は凍結する近似と明記します。物性の温度微分や熱伝導非対角まで追加する必要はありません。

2. **[Major] `node` 専用というスコープを、記載された実装条件では保証できません。**

   根拠: [plan:29](/home/sano/work/forge-integ-1005/plans/active/time_integration-implicit-thermal-jacobian.md:29) は `cell` を対象外としていますが、[plan:54](/home/sano/work/forge-integ-1005/plans/active/time_integration-implicit-thermal-jacobian.md:54) の `!(isNode && !has_nbr)` は、`isNode == 0` なら境界面を含めて真になります。[設定検査:63](/home/sano/work/forge-integ-1005/plans/active/time_integration-implicit-thermal-jacobian.md:63) にも離散化の条件がありません。

   **対案:** 有効化時に `mesh.discretization == "node"` を必須とし、熱伝導分岐自体も `isNode && has_nbr` で限定してください。対象外設定の起動拒否を検証項目へ追加します。現在の検証方針どおり、cellの計算回帰を新設する必要はありません。

3. **[Major] 検証がノズルの挙動に偏り、Jacobianの実装誤りや線形解法の停止を見逃します。**

   根拠: [plan:74](/home/sano/work/forge-integ-1005/plans/active/time_integration-implicit-thermal-jacobian.md:74) はFP32をコンパイル確認だけとし、[§6:90](/home/sano/work/forge-integ-1005/plans/active/time_integration-implicit-thermal-jacobian.md:90) に温度微分・行列作用の検証がありません。V1の「有限かつ振動が増えない」は、更新が止まっても成立します。実際、[Thomas因子分解:1816](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1816) は失敗フラグを立て、[solve:1850](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1850) は前回補正を保持します。

   **対案:** V0より前に、次の小規模検証を置いてください。

   - 固定組成・凍結物性で、CPG/TPの温度微分を有限差分と照合。非零速度、異なるエネルギー基準、ST=float/doubleを含める。
   - 実装した対角ブロックの作用と、独立に組んだ参照行列を照合。壁行の上書きとラインへの保存も確認する。
   - FP32で熱伝導を有効にした小さなnodeケースを実行し、全保存量残差・温度・線形solve失敗件数を確認する。

   V0はpointだけでなく、**新バイナリ・キー0で元のdirectional条件も再現**してください。V1にはsolve失敗ゼロと、全残差・物理状態が悪化していない条件を足します。FP64ビルドとST=doubleは別軸なので、検証表にも分けて記載すべきです。

4. **[Major] V2・V3の「欠損が小さい」では、質量収支改善と高速化を誤判定できます。**

   根拠: [plan:95](/home/sano/work/forge-integ-1005/plans/active/time_integration-implicit-thermal-jacobian.md:95) は符号付きの `Σ res_ro×2π` を比較します。負側への大きな行き過ぎも「小さい」を満たし、局所残差の正負相殺も隠れます。また、同じstepでの改善と `ms/step` の記録だけでは、同じ誤差水準に達する時間の優劣が定まりません。

   **対案:** 合否量を **|2πΣres_ro|/ṁ_in** とし、局所残差のノルムも併記してください。高速化は、事前登録した同一の誤差水準までの壁時計時間で判定します。20,000 step内に双方が達しない場合は、共通の経過時間で比較し、「過渡の改善」に結論を限定します。

   §6には `check_convergence.py` と `check_quasisteady.py` の実行・対象区間・VERDICT保存を明記してください。短期試験で `NOT CONVERGED` や `DRIFTING` を許容すること自体は妥当ですが、安定化・過渡加速と、収束・定常解一致を別の判定にする必要があります。

5. **[Minor] 証拠が支持する範囲より、熱伝導への原因帰属が強く書かれています。**

   根拠: [plan:46](/home/sano/work/forge-integ-1005/plans/active/time_integration-implicit-thermal-jacobian.md:46) は粘性段のエネルギー残差増大を熱伝導の増大として扱っています。しかし [viscousFlux_d.cu:336](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/viscousFlux_d.cu:336) は粘性仕事を作り、その後に熱伝導を加えます。段差分だけでは両者を分離できません。断熱化も温度場と壁拘束を同時に変えています。

   **対案:** 「粘性段の増大が熱伝導Jacobian不足と整合する」という仮説表現に留め、V1では熱伝導と粘性仕事を別々に記録してください。また、非対角と壁拘束を残す今回の変更が失敗しても、「熱伝導の陰的扱い全般」を棄却せず、今回の対角近似では不十分だったと判断します。

**推奨は、上記修正後の「実際の面熱伝導率を凍結した、対角のみのopt-in試験」に絞ります。** 実装前の優先順は、①係数の整合、②node限定の保証、③局所検証と失敗検知、④誤差・時間による合否条件の確定です。

固定組成の単相気体に対する ∂e/∂Q の式と加算符号は妥当です。既存の一般EOS planは対流Jacobian、line-v2は保存量へのスカラー粘性結合なので、本件は重複していません。近傍熱結合と壁拘束の同時変更を避ける限定スコープにも合理性があります。ただし「RHS不変」は同じ離散残差の根を維持する根拠であり、その根への到達保証にはなりません。

ファイル変更・forge起動は行っていません。**plan未反映**です。

指摘数: Critical 0 / Major 4 / Minor 1
