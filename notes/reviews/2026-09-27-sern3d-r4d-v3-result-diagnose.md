# codex 諮問 (diagnose): sern3d-r4d-v3-result

- **brief**: [`../forge-sern-design/notes/reviews/briefs/2026-09-27-sern3d-r4d-v3-result.md`](../../../forge-sern-design/notes/reviews/briefs/2026-09-27-sern3d-r4d-v3-result.md)
- **plan**: [`../forge-sern-design/plans/active/tooling-nozzle-sern-3d.md`](../../../forge-sern-design/plans/active/tooling-nozzle-sern-3d.md)
- **date**: 2026-09-27
- **commit**: `f255f889` (feature/gap-heating-precision)
- **codex**: effort `high`, 3.4 min, rc=0
- **結論**: **(b)を先行し、g3の同一格子で `side_far` だけをslipから検証済みの特性型遠方境界へ変更するA/Bを1組行う。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表です。コードは依頼元の `/home/sano/work/forge-sern-design/` を確認しました。**対象runはローカルに存在せず、数値とVERDICTは依頼文・planの報告値であり、原票を独立確認した結果ではありません。**

| 重大度 | 対象・採否 | 根拠と対案 |
|---|---|---|
| **Major** | **分岐A「現側方領域は不合格」：採用** | [plan §5.1 R4d](/home/sano/work/forge-sern-design/plans/active/tooling-nozzle-sern-3d.md:1942)の事前規則に対し、報告されたDは `C_L` が許容の **1.48倍**、`C_M` が **4.36倍**。「現BC・g3において、この側方拡大への感度が許容超」と記録する。拡大側の領域が十分であることや、反射が真因であることまでは認定しない。 |
| **Major** | **節点Mₙから境界流入・反射経路を認定：却下** | [slip実装](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/boundaryCond_d.cu:68)は境界用速度を Uᵦ = Uᵢ − (Uᵢ·n)n とする。[境界流束](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_boundary_d.inc.cuh:192)はその速度から質量流束を計算する。内部節点のMₙ≠0は、slipを通る流入の証拠ではない。対案は、内部状態のMₙと実際の境界流束を分け、面積重み付きの統計を出すこと。 |
| **Major** | **次は幅の系列を一括投入：却下、(b)を先行** | [runner:63](/home/sano/work/forge-sern-design/design/forge_design/evaluate/runner_sern3d.py:63)で `side_far` は閉じたslip境界に固定されている。まず同一格子で、この閉塞条件を外す影響を分離する。`top_out`も同時変更すると原因分離ができない。 |
| **Major** | **領域誤差を注記するだけでC_L・C_Mを生産の目的・制約に使用：却下** | [§4.46](/home/sano/work/forge-sern-design/plans/active/tooling-nozzle-sern-3d.md:1845)の格子差と今回のDの単純和でも `C_L ≈ 0.002121`、`C_M ≈ 0.062295` と総許容を超える。しかもDは真の誤差上限ではない。**探索用の参考値に留め、設計の合否・最終順位を確定する用途は保留**する。正式なGは同じ窓の平均・振幅で再計算する。 |
| **Major** | **V4–V6を現BCで続行：保留** | 今後BCを変えるなら、その領域試験は最終仕様の認定に使えない。**側方BCの診断→遠方BCの確定→最終BCでV3–V6と同時拡大確認**の順にする。`top_out`は別途確認し、Mₙ≈0の`bottom`もBC適合性の確認対象に残す。 |

結論: **(b)を先行し、g3の同一格子で `side_far` だけをslipから検証済みの特性型遠方境界へ変更するA/Bを1組行う。**

第1仮説: **現側方位置でのslipによる閉塞が、C_L・C_Mに許容超の偏りを与えている。** 確度: **中**。  
　根拠: 上記コードで質量を通さない境界であることを確認した。また、報告では `case/46.sern_design/run_0986_r4d_g3_base/` と `case/46.sern_design/run_0988_r4d_g3_zapp075_cont20k/` の領域変更で、C_L・C_MのDが許容を超えている。ただし「下流のプルーム反射がランプへ戻る」という具体的経路は未確認。  
　反証条件: 下記の比較が有効に成立し、C_L・C_MともD≤εとなること。「同じ位置で側方BCだけを変えれば許容超に動く」という仮説を棄却する。あらゆる反射の不存在を意味しない。

第2仮説: **旧境界節点が内部節点になる際の双対CV・離散化の変化と、追加した粗い外側セルが、V3の差に寄与した。** 未確認。共通節点の座標・保存量のビット一致だけでは、旧境界付近の離散演算まで同一とは保証できない。

第3仮説: 設定しない。

判別A/B:

- **変更点は `side_far` のBCだけ。** Aはslip、Bは特性型遠方境界。格子は元のg3、遠方面z/H=2.5に固定し、`top_out`・`bottom`・`outlet`は維持する。
- 両者とも `run_0986_r4d_g3_base/` の同じ最終場から、同一格子restartで新規runへ分岐する。**同一の新バイナリ**を使い、Aを対照にする。
- Bは流入特性に外気状態を与え、流出特性は内部から取る。TP熱力学・組成・SSTと整合し、自由流保持と波の出入りを検証した実装を前提とする。開発費は別途必要であり、単なるYAML変更では成立しない。
- 各20000 step、500 step間隔。末尾10000 stepでDを算出し、前10000 stepとの平均差≤0.1εを要求する。未達なら延長。εは既定どおり、`C_T`・`C_T_with_shear`・`C_L`が0.0005、`C_M`が0.005。
- 4量の`check_quasisteady`、床・NaN検査を確認し、`check_convergence`も併記する。残差プラトーだけでは、既存のユーザ決定に反して却下しない。

→ **結果A：C_LまたはC_MでD>εなら、側方BCの影響を支持する。** BCを確定してから領域系列へ進む。ただし、これだけでBの精度や反射経路を証明したことにはならない。  
→ **結果B：両量ともD≤εなら、第1仮説の上記定量形を棄却する。** V3の差をslip閉塞だけで説明せず、第2仮説を候補に残す。  
準定常条件やBC実装の検証を満たさなければ判定不能。

やらない方がよいこと: **slip→outflowを「非反射化」と呼んで生産採用すること。** `_bcond_config`の独立実行でも、`top_out_kind: outflow`は実際にはslipになった。側方と上方を同時変更すること、今回の結果に合わせて許容を緩めること、C_TのV3合格だけで全体を認定することも避ける。

呼び出し側の前提への異議: 「節点速度を射影しない」という説明は、読んだslip実装と整合する。一方、**「流入51.7%」は内部節点の速度符号の割合**と表記すべきで、実流入率・流入面積率ではない。また、面別の最終snapshotは原因探索の補助に留め、「差はほぼランプ」を確定するなら同じ比較窓の面別力時系列が必要。

不足情報: 対象runの実効設定、比較窓の係数時系列、判定区間付きVERDICT原票、境界統計の抽出方法。報告上の判定は **`GATES: PASS`／4量`STEADY`／`check_convergence: NOT CONVERGED`** であり、残差収束は主張しない。ファイル変更・forge実行はしていない。**plan未反映**。採用時の反映先は `plans/active/tooling-nozzle-sern-3d.md` §5.1 R4d・§6・§8。
