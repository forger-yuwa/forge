# codex 諮問 (diagnose): core-grid-mechanism

- **brief**: [`notes/reviews/briefs/2026-10-11-core-grid-mechanism.md`](../../notes/reviews/briefs/2026-10-11-core-grid-mechanism.md)
- **plan**: [`plans/active/tooling-nozzle-core-grid.md`](../../plans/active/tooling-nozzle-core-grid.md)
- **date**: 2026-10-11
- **commit**: `96197821` (feature/nozzle-core-grid)
- **codex**: effort `high`, 3.0 min, rc=0
- **結論**: **F1 を採用し、上限 50 の B0 で `FORGE_FREEZE_TURB` だけを変える 2 万 step の比較を行う。ただし、判別するのは「SST の更新が停滞を支えているか」であり、「μt の遅れが原因」とはまだ判定しない。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論: **F1 を採用し、上限 50 の B0 で `FORGE_FREEZE_TURB` だけを変える 2 万 step の比較を行う。ただし、判別するのは「SST の更新が停滞を支えているか」であり、「μt の遅れが原因」とはまだ判定しない。**

第 1 仮説: 平均流と SST の分離反復が、縮流部の温度境界層外縁にある遅いモードを維持している。確度: **中**
  
根拠:
- 平均流の block 更新後、SST は組立済みの残差から別に更新される。[main.cpp:2378](/home/sano/work/forge-coregrid/solver_density_cuda/main.cpp:2378)
- SST の分母は `V/dt_local + V·src_jac + transport_diag` で、方向別 dt を使う一方、更新式には隣接節点の未知増分がない。[update_d.cu:403](/home/sano/work/forge-coregrid/solver_density_cuda/cuda_forge/update_d.cu:403)
- 記録上、上限 1 は残差集中域の dt をほぼ変えず（比の中央値 0.977）、同域の残差を 1/25〜1/160 に下げた。局所の dt だけでなく、壁側との結合を調べる根拠になる。ただし、その結合が SST 経由か平均流経由かは未確定。[tooling-nozzle-core-grid.md:746](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:746)

反証条件: 凍結の実効性を確認したうえで、後述の両窓で縮流部の ρ・ρE 残差が対照の 2/3 以上に残り、同じ場所のモードが持続すること。この場合、「SST 更新を止めれば 2 万 step 内に停滞が大幅に緩和する」という限定した仮説を退ける。

第 2 仮説: 平均流側のライン外結合・熱伝導 Jacobian・有限回の内部反復だけでも、細かい x 格子上の遅いモードが残る。確度: **中、未確認**。凍結後も停滞が残る場合の候補。
第 3 仮説: 追加しない。

判別 A/B: **変更は `FORGE_FREEZE_TURB=0/1` の一点、比較長は各 20,000 step。**

- 起点は `case/45.isobutane_m6_d155/run_0608_cg_gref_cont_k1/res_60000.h5`。同一格子への `restart_field.py` を使い、保存量のビット一致を確認する。B0・上限 50・バイナリ・境界条件は固定する。
- 非凍結の対照は、同じ起点・構成が記録され、領域残差も保存済みの `case/45.isobutane_m6_d155/run_0619_cg_gref_b0_ctl/` の **0〜20k** を使う。凍結側だけ新しい、番号の衝突しない `run_*` に作る。`run_0615` は全域残差の補助照合に使える。
- 判定窓は **10k〜15k と 15k〜20k**。全域の平均流・化学種の全残差列と、縮流部 x/r_t ∈ [−8, −2] の `q = √(ΣV(R/V)² / ΣV)` を比較する。領域指標は対照と同じ保存時刻・同じ集計方法にそろえる。
- 凍結側では、拘束境界を除く `roK`・`roOmega` の変化と、`vis_turb`・`dt_local` の変化も記録する。凍結が効かなければ判定不能。

→ **結果 A:** 両窓で全域の対象残差列、領域の ρ・ρE がすべて対照の **1/3 以下**。SST 更新への強い依存を支持し、「SST 更新を止めても同程度の停滞が残る」を退ける。μt の遅れや x 格子依存の機構まで確定しない。  
→ **結果 B:** 両窓で領域の ρ・ρE がともに **2/3 以上**で、同じ場所のモードが残る。SST の内部更新なしでも停滞が持続すると判断し、第 2 仮説を優先する。  
→ 中間・成分間の不一致・新しい過渡への置換は判定不能。非有限は別扱い。自動延長しない。倍率は診断用であり、収束判定ではない。

やらない方がよいこと: F2、追加の格子変更、CFL・ITJ・内部反復数の同時変更。凍結で ω の更新が止まることを、ω の方程式が解けた証拠にしない。また、凍結側の θ_r・Q_w が point の値に近づくことを合否に使わない。解いている系が変わるためである。

呼び出し側の前提への異議:

- **[Major] 「凍結すると k・ω、したがって μt が実質固定」は未成立。** 止まるのは `ρk`・`ρω` の更新で、原始量はそれらを密度で割って再計算する。[dependentVariables_d.cu:307](/home/sano/work/forge-coregrid/solver_density_cuda/cuda_forge/dependentVariables_d.cu:307) μt も密度・ひずみ・F2 に依存する。[turbulent_viscosity_d.cu:235](/home/sano/work/forge-coregrid/solver_density_cuda/cuda_forge/turbulent_viscosity_d.cu:235) **対案:** 「SST 保存量の内部更新を止める試験」と定義し、μt の変化を実測する。
- **[Major] F1 の成功だけから「分離更新のずれが原因」とは言えない。** 凍結は更新方式の改善ではなく、SST 方程式を解く操作の停止である。失敗しても平均流には乱流粘性の影響が残る。根拠は更新ゲートの範囲。[main.cpp:2380](/home/sano/work/forge-coregrid/solver_density_cuda/main.cpp:2380) **対案:** 上の A/B の限定した読み方を採る。
- **[Major] 「前線に沿った自由度が増えるから遅い」はまだ説明になっていない。** Grefc は列数だけでなく、離散演算子・局所 dt・壁側の結合も変える。§4.37 が示すのは局所粗化への感度である。記載された `check_convergence` は対照・Grefc とも **NOT CONVERGED**、副判定も到達未成立。[tooling-nozzle-core-grid.md:757](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:757) **対案:** 今回は SST 更新への依存を一段切り分け、節点数との比例則や機構の確定を保留する。

不足情報: 実行バイナリと読んだソースの対応、凍結時の境界上書き・μt の変化量、細／粗格子でのモード減衰率と内部線形反復の解き残し。禁止された生データは読まず、数値と VERDICT は依頼文・plan の記録を根拠にした。ファイル変更・forge 起動はしていない。**plan 未反映**。採用時の反映先は呼び出し側の §4.38。
