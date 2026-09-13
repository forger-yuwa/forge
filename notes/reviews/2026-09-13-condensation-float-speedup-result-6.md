# codex レビュー: condensation-float-speedup (result)

- **plan**: [`plans/active/condensation-float-speedup.md`](../../plans/active/condensation-float-speedup.md)
- **stage**: `result` (diff base `3b77cc4b`)
- **date**: 2026-09-13
- **commit**: `709cb9e0` (feature/perf-3d-speedup)
- **codex**: effort `high`, 3.7 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M0/m2
- **focus**: 凝縮 float 化の result レビュー 6 回目。5 回目 (NO-GO M1/m2) の採用: 表範囲判定に src_jac の摂動点 T+0.1 K を含め跨ぐセルは double 委譲 (境界 125.45–125.7 K の device 試験追加)、最終バイナリ forge_condf4 で A10G を再計測 (dry 33.9 / 起動 44.0 / 発達場 52.3 ms/step; 証拠 _aws_perf_evidence/cond/ を最終 SHA で回収、300 step CMP_VERDICT と統計は condf4)、§4.2-5・methods §9・README の文言を訂正。発達場 ≤50 の 4.6 % 未達は followups F-cf3 に移管済み。accepted に移してよいか判定してください。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**

5 回目の Jacobian 退避漏れは修正され、`condf4` の A10G 計測証拠も確認できました。  
新たな Critical / Major は見つかりませんでした。ただし、「訂正済み」とされた文書に旧記述が残っており、移動前の同期が必要です。

1. **Minor — 最終バイナリと検証範囲の記録が統一されていない。**

   [plan §5.1 #15:195](/home/sano/work/forge-perf/plans/active/condensation-float-speedup.md:195) と [§8:233](/home/sano/work/forge-perf/plans/active/condensation-float-speedup.md:233) は、最終版を依然として `1551099d / forge_condf3` としています。「最終差分は二相音速の復帰のみ」も、今回の範囲判定変更を含みません。

   実際の [BINARY_SHA256.txt:2](/home/sano/work/forge-perf/case/16.nozzle_wys/_aws_perf_evidence/cond/BINARY_SHA256.txt:2) は `forge_condf4 = ef3c6fae`、SHA256 `3604fd4e824e46b6…`。18本の計測ログもこの SHA に対応しています。

   **対案：** #15・§8 を `condf4` に同期し、旧 SHA のローカル検証については、今回の変更が適用される温度範囲と境界試験による補完を明記してください。#5 の「AWS 正式取得未了」も、ローカル正式採用への変更と整合させる必要があります。

2. **Minor — 平衡形と `condFloat: 0` の現在仕様に訂正漏れがある。**

   [methods/condensation.md:890](/home/sano/work/forge-perf/methods/condensation.md:890) は、温度反転について平衡形を一括して「double のまま」としています。しかし [dependentVariables_d.cu:140](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/dependentVariables_d.cu:140) の専用分岐は `condEquilibrium == 2` であり、`1` の湿潤 TP セルはハイブリッドへ進みます。末尾の補足だけが正しい状態です。

   また、[plan:125](/home/sano/work/forge-perf/plans/active/condensation-float-speedup.md:125) の「失敗時のみ挙動差」は、成功時にも double Newton 後の研磨が入る実装と矛盾します。

   **対案：** `methods` 本文・設定表と plan 分岐表を、修正済みの `solver-settings.md` に統一してください。ソースの退避条件も、現在温度だけでなく **`T+0.1 K` を含む**と明記してください。

判断の裏付けとして、以下を確認しました。

- **実装：** [退避条件:355](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/condensationSourceKernels_d.cuh:355) は摂動点を含み、範囲外は旧 double セル関数へ委譲します。旧 source／clamp 本体は、変更前と空白を除いて同一でした。面潜熱は引き続き面カーネル内で評価しています。
- **速度：** `case/16.nozzle_wys/_aws_perf_evidence/cond/` の最終版ログで、dry **33.92–33.94**、起動 **44.01–44.04**、発達場 **52.30–52.34 ms/step**。発達場は基準から約 **46.6%短縮**、50 ms 目標には約 **4.6%未達**です。dry は約 +0.15 ms の小幅増で、厳密な不変ではありません。
- **物理検証：** `case/16.nozzle_wys/run_0460_condf_phys_node2d_cond/{merged,condf}/` を再判定し、双方 **`NOT CONVERGED (stalled/plateau)`**、`pmax/machmax` は **`OVERALL: ALL STEADY`**。比較は **`PASS 32/32`**、h0 は新旧差 **+0.068 J/kg**、**`ABSOLUTE_VERDICT: PASS / NONDEGRADATION_VERDICT: PASS`**。収束解の一致を証明した結果ではありません。[run索引](/home/sano/work/forge-perf/case/16.nozzle_wys/README.md:76)
- **互換性・未達項目：** `run_0458_frozen1step_node2d_cond/{merged,new_cf0}/res_1.h5` の T・P・凝縮診断6量のビット一致を再確認しました。周期 `run_0053` の20 step収支誤差は旧 **1.143×10⁻³**、新 **1.111×10⁻³**で、保存ゲート未達は事実です。cell h0 絶対ゲートと性能残件も含め、[followups §5.1:53](/home/sano/work/forge-perf/plans/active/condensation-followups.md:53) に `F-cf1/F-cf2/F-cf3` が残っています。

**推奨は、上記1→2の順で記録を訂正し、今回のレビューを §6.1 に残してから `accepted` へ移すことです。** 性能残件と既存の保存誤差は、明記された followups で継続する判断を支持します。

既存 host 試験バイナリは `ALL PASS`。device 試験は CUDA driver/runtime 不整合により、この環境では再検証できませんでした。ファイルは変更しておらず、指摘は **plan 未反映**です。

指摘数: Critical 0 / Major 0 / Minor 2
