# codex 諮問 (diagnose): linedir-divergence

- **brief**: [`notes/reviews/briefs/2026-10-08-linedir-divergence.md`](../../notes/reviews/briefs/2026-10-08-linedir-divergence.md)
- **plan**: [`plans/active/tooling-nozzle-isothermal-wall-chain.md`](../../plans/active/tooling-nozzle-isothermal-wall-chain.md)
- **date**: 2026-10-08
- **commit**: `6c124d46` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 4.6 min, rc=0
- **結論**: **まず実効の `time.deltaT.implicitSolvePrecision` を確認し、0 ならこれだけを 1 に変えた A/B を行う。FP64 ビルドという説明だけでは、陰解法内部の単精度誤差を除外できない。**
- **extra**: `plans/accepted/time_integration-line-implicit-viscous-v2.md`, `solver_density_cuda/cuda_forge/setDT_d.cu`, `methods/time_integration/implementation.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論: **まず実効の `time.deltaT.implicitSolvePrecision` を確認し、0 ならこれだけを 1 に変えた A/B を行う。FP64 ビルドという説明だけでは、陰解法内部の単精度誤差を除外できない。**

第 1 仮説: **単精度で組み立てたライン行列の誤差が、方向別 Δτ によって時間対角の保護を失い、近壁の成長振動として現れる。** 確度: **中。ただし実効 `implicitSolvePrecision=0` が条件で、当該 run では未確認。**

  根拠:

- **[Major] FP64 ビルドでも行列組立ては自動では倍精度にならない。** [timeIntegration_d.cu:1531](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1531) は、設定が 1 の場合だけ `<double>`、それ以外は明示的に `<float>` を起動する。`flow_float` の typedef 変更では、この分岐は変わらない。既定 0 は [implementation.md:235](/home/sano/work/forge-integ-1005/methods/time_integration/implementation.md:235) にも記載されている。
- **座標差を取る前に丸めている。** [timeIntegration_d.cu:926](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:926) は `ST(cc_neighbor) − ST(cc_self)` で粘性対角用の幾何を作る。`ST=float` なら、倍精度メッシュでも極薄層の距離が再び単精度に量子化される。plan §5.1 #16 が問題にした「第一層が数 ulp」という条件を、LHS 内で再導入し得る。
- Thomas の消去自体は double だが、入力は既に組み立てた D・K である。[timeIntegration_d.cu:1794](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1794) の D − Kprev·W は、失われた入力精度を回復しない。方向別 Δτ は [setDT_d.cu:161](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/setDT_d.cu:161) で内部ライン面を除外し、[timeIntegration_d.cu:827](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:827) の V/Δτ を弱める。

  対案: **次数・CFL・緩和を固定したまま、組立て精度だけを切り替える。** これはコード変更不要で、空間離散化も維持できる。

  反証条件: 元の実効設定が既に 1 なら、この仮説は対象 run に当てはまらない。0→1 でも同じ場所で同程度の初期成長率・周期を持つ振動が残れば、「単精度組立てが主因」を棄却する。NaN の到達時刻だけが遅れる場合は支持としない。

第 2 仮説: **近似 LHS と実際の残差・壁処理の不整合が、Δτ 拡大で顕在化する。** 確度: 中。FVS の近傍結合を抽出しているだけで、2 次 SLAU の残差全体を厳密に微分した行列ではない（[timeIntegration_d.cu:896](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:896)）。粘性結合も別スイッチであり、「壁法線方向をすべて厳密に陰化した」という前提は成立しない。

第 3 仮説: **5 sweep では off-line の補正結合を十分処理できない。** 確度: 低〜中。off-line のみ旧補正を読む実装は確認できる（[timeIntegration_d.cu:866](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:866)）が、この格子での sweep 履歴がない。周期約16 stepという外側更新の観測だけでは、線形反復不足と第2仮説を区別できない。

判別 A/B: **変更は `time.deltaT.implicitSolvePrecision: 0 → 1` の一点。**

- 共通条件は `run_0203_ns_coldmesh_tw300_linedir_diag` の開始場・固定リミッタ基準・`cfl_pseudo=4`・`nStepInner=5`・`implicitRelax=0.7`・2次・directional ON。同一バイナリで、新規 run に分岐する。
- 最大 **200 step**、非有限で停止。既存の527節点の帳簿を毎 step 保存する。
- 見る量は全残差列、ρ・P・壁法線速度の増分、非有限の初発位置。増分は両腕共通の初期状態で正規化し、**16 step 窓の振幅と窓間の増幅率**を比較する。小さい局所 Ux だけを分母にしない。
- **A の単精度側で既知の成長が再現され、B の倍精度側では200 stepまで有限、末尾3窓の振幅が増えない**なら、第1仮説を支持する。ただし幾何の量子化と行列演算の丸めの寄与までは分離できない。
- **B でも同じ位置・周期・成長傾向が残る**なら、第1仮説を主因から外す。両腕とも再現しない場合や、成長が遅れるだけの場合は判定保留。

やらない方がよいこと:

- 精度を確認する前に、`convMethod`・`implicitRelax`・CFLを同時に変えること。安定化しても機構が分からなくなる。
- 壁境界面の時間刻み制約を外すこと。今回の切り分けにならず、既往の発散記録もある。
- 直ちに Δτ 上限や第一層だけの例外処理を実装すること。いずれもコード変更を伴う安定化策で、今回の単因子診断の後に判断する。
- 200 stepの有限性を、収束や高速化の証明として扱うこと。

呼び出し側の前提への異議:

- **[Major]「CFLによらない」は未立証。** [case README:131](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/README.md:131) の記録は、4で65 step、1で129 stepの発散。同じ場所という観測はあるが、時間刻みに依存しない証拠ではない。対案は「試した2値では安定化しなかった」と記述し、実効 `cfl_pseudo` も確認すること。
- **[Major] 壁と第一内点の Δτ の差だけでは、拘束の不整合を説明できない。** 壁の運動量・エネルギー行は単位行化される一方、密度行は残る（[timeIntegration_d.cu:1029](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1029)、[1053](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1053)）。対案は、壁密度・圧力と第一内点の補正の位相を見ること。壁全体を固定 Dirichlet と扱わない。
- **[Minor] 全ラインのNaNは発生機構の指紋ではない。** Thomas の前進・後退代入で非有限がライン内に伝播し得る。対案はNaN直前の全断面分布より、初期の有限な成長を優先すること。
- **追究する価値は、この短い診断一回にはある。高速化の採用判断はまだできない。** 約2倍のstep単価なら、同じ収支・収束ゲートへの到達step数を半分未満にする必要がある。旧格子の結果をこの格子の限界に一般化しないことは、[既存 plan:222](/home/sano/work/forge-integ-1005/plans/accepted/time_integration-line-implicit-viscous-v2.md:222) とも整合する。

不足情報: 対象 `run_0183`・`run_0191`・`run_0200〜0203` のディレクトリは、この作業ツリーには存在しなかった。したがって実効設定、帳簿、残差、`CONVERGENCE_VERDICT.txt` は直接検証できていない。必要なのは特に **`implicitSolvePrecision`、`cfl_pseudo`、`lineViscCoupling`、実行バイナリとソースの対応**。上記はコードに基づく診断仮説であり、run 数値は提供資料の記録として扱った。ファイル変更・forge起動は行っていない。**plan 未反映**（依頼どおり、反映は呼び出し側）。
