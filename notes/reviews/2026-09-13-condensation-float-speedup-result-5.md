# codex レビュー: condensation-float-speedup (result)

- **plan**: [`plans/active/condensation-float-speedup.md`](../../plans/active/condensation-float-speedup.md)
- **stage**: `result` (diff base `3b77cc4b`)
- **date**: 2026-09-13
- **commit**: `861245bf` (feature/perf-3d-speedup)
- **codex**: effort `high`, 6.9 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M1/m2
- **focus**: 凝縮 float 化の result レビュー 5 回目。4 回目 (NO-GO M2/m1) の採用: 二相音速の dL/dT を double 差分に戻した (1551099d)、最終バイナリ forge_condf3 (sha256 f3d03e2d) で A10G を再計測 (dry 33.9 / 起動 44.0 / 発達場 52.3 ms/step, 100 step ×3 交互)、run_0421 の 300 step CMP_VERDICT と分布統計と全 sha256 を case/16.nozzle_wys/_aws_perf_evidence/cond/ に回収、§5.1 #13/#15・§6・§8 をコミット/SHA/環境付きに整理。発達場の目標 ≤50 に対する 4.6 % 未達は followups F-cf3 に移管した上で、accepted に移してよいか判定してください。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**

前回の二相音速修正と、最終SHAによる44.0／52.3 ms/stepの計測は確認できました。  
ただし、N2物性表の上端を跨ぐJacobian摂動に退避漏れがあり、陰的係数を大きく変えます。  
性能の4.6%未達を `F-cf3` へ移すこと自体は妥当ですが、この実装不具合は受け入れ前に修正が必要です。

1. **Major — 範囲内セルの `T+0.1 K` が表外へ出ても、Jacobian評価がdoubleへ退避しない。**

   [condensationSourceKernels_d.cuh:353](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/condensationSourceKernels_d.cuh:353) は現在温度 `Td` だけで退避を判定します。その後、[同:445](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/condensationSourceKernels_d.cuh:445) は `Td+0.1f` を無条件にfloat表で評価します。N2表の上端は125.6 Kなので、摂動先の物性は [condensationTables_d.cuh:49](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/condensationTables_d.cuh:49) で125.6 Kへクランプされます。

   コードの物性式・Hermite係数構築・float32評価をPythonで独立再評価した反例です。**GPU runの実測ではありません。**

   条件：pure N2、`condKantrowitz=1`、`T=125.6 K`、`P=3305011 Pa`、`ρ=88.658295 kg/m³`、`g=Q0=Q1=Q2=0`、`dt=1e-7 s`。

   | 評価量 | 従来double | float表経路 |
   |---|---:|---:|
   | 摂動先の過飽和度 | 約0.99654 | 約1.00100 |
   | 摂動先の質量ソース `Sg` | 0 | 約1.217×10⁵ |
   | `src_jac_g` | 約1.276×10⁵ s⁻¹ | 約3.508×10³ s⁻¹ |

   **Jacobianを約97%過小評価**し、`|Δsj|dt≈1.24e-2` になります。核生成率は約2.72×10²⁷で上限未満、律速係数も1です。上限セルに対する意図的変更では説明できません。

   [device試験:173](/home/sano/work/forge-perf/solver_density_cuda/tests/unit/test_cond_float_device.cu:173) の通常N2状態は110 Kまでで、さらに10 bar超を除外しているため、この境界を検査していません。

   **対案：** 本体だけでなく、**全摂動評価点が物性表の有効範囲内であること**をfloat経路の条件にしてください。範囲を跨ぐセルは既存doubleセル関数へ丸ごと委譲し、125.5–125.7 Kの境界試験を追加することを推奨します。§5.1に未解決項目として追加が必要です。

2. **Minor — 最終版の証拠は揃ったが、要約に旧版の数値・名前が残っている。**

   [plan §5.1 #13:193](/home/sano/work/forge-perf/plans/active/condensation-float-speedup.md:193) と、回収された [CMP_VERDICT:3](/home/sano/work/forge-perf/case/16.nozzle_wys/_aws_perf_evidence/cond/run_0421_CMP_VERDICT_condf3.txt:3) は次のように食い違います。

   | 正規化最大差 | plan | 最終 `condf3` 証拠 |
   |---|---:|---:|
   | `P` | 4.9e-3 | 5.68e-3 |
   | `T` | 1.6e-2 | 1.73e-2 |
   | `rog_0` | 9.8e-2 | 1.10e-1 |

   `PASS 30/30` は維持されますが、「基準ノイズの0.8–0.9倍」は修正が必要です。また、[case README:75](/home/sano/work/forge-perf/case/16.nozzle_wys/README.md:75) は最終版をまだ `forge_condf2` としています。`nsys` CSVも `condf2` 名で、最終版の壁時計計測とは区別すべきです。

   **対案：** planとrun索引を回収済み数値へ同期し、旧版プロファイルには取得版を明記してください。300 stepの分布比較は「当該時点の統計比較」と記し、定常性確認済みの意味で「一致」と表現しないでください。

3. **Minor — 現在仕様と分岐表に、既に訂正された説明が残っている。**

   [plan:114](/home/sano/work/forge-perf/plans/active/condensation-float-speedup.md:114) はdouble反転を「成功時はビット不変」としますが、[dependentVariables_d.cu:167](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/dependentVariables_d.cu:167) はdouble反転後にも研磨します。旧Newtonが返した有限温度も追加補正され得ます。

   また、[methods/condensation.md:890](/home/sano/work/forge-perf/methods/condensation.md:890) とplanの分岐表は平衡形反転を一括してdoubleとしていますが、`condEquilibrium=1` の湿潤TPセルはハイブリッドへ進みます。末尾の補足と本文が矛盾しています。

   **対案：** 本文・分岐表を `procedures/solver-settings.md` の修正済み説明に統一してください。ビット一致の保証範囲はsource／clampの凍結評価に限定します。

今回、以下は裏付けられました。

- `_aws_perf_evidence/cond/` の18本：最終SHA `f3d03e2d…` で、dry **33.90–33.95**、起動 **44.01–44.02**、発達場 **52.30–52.31 ms/step**。
- `case/16.nozzle_wys/run_0460_condf_phys_node2d_cond/{merged,condf}/` の再判定：双方 **`NOT CONVERGED (stalled/plateau)`**、`pmax/machmax` は **`OVERALL: ALL STEADY`**。`h0` は新旧差 **+0.068 J/kg**、**`ABSOLUTE_VERDICT: PASS`／`NONDEGRADATION_VERDICT: PASS`**。
- `run_0458_frozen1step_node2d_cond/{merged,new_cf0}/res_1.h5` のT・P・凝縮診断6量はビット一致。double sourceセル本体とclamp本体も変更前と同一。
- `F-cf1/F-cf2/F-cf3` は [followups §5.1](/home/sano/work/forge-perf/plans/active/condensation-followups.md:53) に存在します。

**推奨は、`in_progress` を維持し、指摘1の退避漏れを修正・検証してから受け入れることです。** 性能目標の残り2.3 msより、既定経路のJacobian互換性を優先してください。

ファイルは変更しておらず、指摘は **plan未反映**です。既存host試験バイナリはALL PASSでしたが、device試験はCUDA driver/runtime不整合で再実行できませんでした。

指摘数: Critical 0 / Major 1 / Minor 2
