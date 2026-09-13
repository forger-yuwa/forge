# codex レビュー: performance-3d-node-sst-speedup (result)

- **plan**: [`plans/active/performance-3d-node-sst-speedup.md`](../../plans/active/performance-3d-node-sst-speedup.md)
- **stage**: `result` (diff base `0512823d`)
- **date**: 2026-09-12
- **commit**: `703cb24e` (feature/perf-3d-speedup)
- **codex**: effort `high`, 4.6 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M1/m2
- **focus**: 4 回目の result レビュー。前回 (result-3) の M2 で指摘された node 軸対称 SST TP (case/44 run_0200) の T 3.8e-5 の決着 (§5.1 #13, §4.3 の double 対照, tools/perf_regress.py --truth, 証拠 case/44.vitiated_air_wt/run_0200_perf_regress_node_axisym_sst_tp/bisect_slau_summary.txt と cmp_final_truth.txt) が妥当か、accepted/ に移してよいかを判定してください。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**GO-with-changes**

case/44 の温度差を理由に SLAU を旧演算順へ戻す必要はありません。double 対照への T 最大差は旧版 `3.58e-5` → 最終版 `3.14e-5` で、3構成の `PASS 18/18` を再現しました。  
ただし、比較ツールの誤合格経路と、完了条件・レビュー記録の不整合は移行前に修正が必要です。

1. **Major — `--noise` / `--truth` がデータ不備を合格に変えてしまう**

   根拠は [perf_regress.py:115](/home/sano/work/forge-perf/solver_density_cuda/tools/perf_regress.py:115) と [同:126](/home/sano/work/forge-perf/solver_density_cuda/tools/perf_regress.py:126)。欠落・形状不一致を表す `inf` が、そのまま許容差や比較基準に使われます。ファイルを作らない再現確認で、次の両方が **PASS 8/8、終了コード0** になりました。

   - ノイズ用データから T を欠落させると `noise=inf` になり、新版の T が基準の100倍でも合格。
   - 基準の `ro` が長さ2、新版と double 対照が長さ1の場合、`truth-dist base=inf new=0` により **SHAPE エラーが合格**。

   **対案:** 基準・新版・ノイズ・double 対照の必須場、形状、有限性を先に検査し、不備は無条件で非ゼロ終了する。数値判定へ進めるのは検査済みデータだけにしてください。この2例を退行試験に追加すべきです。

   なお、今回の case/44 の7ラベルは、別途確認して**場の構成・形状が同じ、非有限値0**でした。この不具合によって今回の `PASS` が生じたわけではありません。

2. **Minor — 「有限時間の非劣化」と「収束解保証」が完了条件で混在している**

   [plan:115](/home/sano/work/forge-perf/plans/active/performance-3d-node-sst-speedup.md:115) は「未収束なら比較を確定しない」とする一方、[残作業 #11:162](/home/sano/work/forge-perf/plans/active/performance-3d-node-sst-speedup.md:162) は完了扱いです。

   再実行した `check_convergence.py` の判定は、case/16 の基準・最終版とも **NOT CONVERGED**。case/44 の `base_r1`・`final`・`base_dbl_r1` も **NOT CONVERGED** でした。したがって double 結果は、同一入力からの有限時間計算に対する高精度の**参照結果**です。[performance.md:69](/home/sano/work/forge-perf/methods/architecture/performance.md:69) の「真値」という呼称は保証範囲を広く見せます。

   **対案:** §1・§4.3・§8 を「指定ケース・継続時間における非劣化と速度改善」に統一し、収束解不変は未検証と明記する。double 対照は「高精度参照」と呼んでください。また、[plan:138](/home/sano/work/forge-perf/plans/active/performance-3d-node-sst-speedup.md:138) に残る却下済みのセル前計算案も、実装どおり面状態評価へ修正してください。

3. **Minor — result-3 の採否がレビュー表に残っていない**

   [result-3 の記録](/home/sano/work/forge-perf/notes/reviews/2026-09-12-performance-3d-node-sst-speedup-result-3.md:1) は存在しますが、[plan §6.1:194](/home/sano/work/forge-perf/plans/active/performance-3d-node-sst-speedup.md:194) には result-3 の行がありません。#13 の説明だけでは、前回 M1 の証拠回収や各 Minor の対応状況を追えません。

   **対案:** result-3 の各指摘の採否・対応先を記録し、今回の修正事項を §5.1 に追加する。修正後の再集計結果を記録してから完了扱いにしてください。`check_plans.py` は **PASS** でしたが、この内容不足は検出していません。

**推奨は、SLAU の現実装を維持し、上記1→2→3を修正してから `accepted/` へ移すことです。**

重点の #13 を支持する実測は十分あります。`case/44.vitiated_air_wt/run_0200_perf_regress_node_axisym_sst_tp/` の最大温度差の節点12813では、旧版／最終版／double 対照がそれぞれ **838.542175／838.501587／838.507836 K**。対照との差は **0.034339 → 0.006249 K** に縮小しています。軸上2401節点の `roUy` は、対角キャッシュ有効版を含む全7ラベルで **0** でした。これは運動量流束3行の二分結果と整合します。証拠の索引は [case/44 README](/home/sano/work/forge-perf/case/44.vitiated_air_wt/README.md:522) です。

前回不足していた case/16 の基準残差と両者の準定常判定も回収されています。`_aws_perf_evidence/ref_run_0234_evidence/` と `run_0416_final_sweep4/` の保存判定は、pmax/machmax とも **OVERALL: ALL STEADY**。最終版の **32.62 ms/step** は生ログで確認しました。収束解保証と区別すれば、高速化の採用根拠になります。

ファイルは変更していません。本レビューは **plan 未反映**です。

指摘数: Critical 0 / Major 1 / Minor 2
