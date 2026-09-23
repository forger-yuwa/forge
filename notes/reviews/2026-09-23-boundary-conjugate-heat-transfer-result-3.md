# codex レビュー: boundary-conjugate-heat-transfer (result)

- **plan**: [`plans/active/boundary-conjugate-heat-transfer.md`](../../plans/active/boundary-conjugate-heat-transfer.md)
- **stage**: `result` (diff base `main`)
- **date**: 2026-09-23
- **commit**: `7c1ac6c4` (feature/cht-phase2-fem2d)
- **codex**: effort `high`, 6.9 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M4/m1
- **focus**: 2 巡目 (notes/reviews/2026-09-23-boundary-conjugate-heat-transfer-result-2.md) の Major 5 + Minor 1 への対処を最優先で検証する。(1) 物理残差 r=K(u)u-b-E^T Qf を直接組む形が量子化停滞を防ぐか、(2) ゲート窓を動かさない実装が末尾の不良を捨てられなくなったか、(3) content_sha1 照合・更新位相の保持・最終 step 強制保存、(4) cht_loop/check_cht_balance の積分済み荷重統一、(5) 0 からの発散検知。再判定は run_0154_r2_avg42 (C3X)・run_0155/0156 (再開同値, interval 非倍数で分割)・run_0036_r2_avg42 (Mark II)

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**NO-GO**

C3X の G-if は再実行でも PASS でした。前回の温度精度混用・末尾行の除外・ゼロ始点の発散検知は改善しています。  
ただし、再開時刻、収支ゲートの比較状態、流束選択に欠陥が残っており、`accepted` への移動は見送るべきです。

1. **Major — 最終保存が1ステップ早く、再開位相も一部の設定で失われます。**

   [main.cpp:1808](/home/sano/work/forge-cht/solver_density_cuda/main.cpp:1808) は出力へ `iStep+1` を渡しますが、[conjugateWall.cpp:993](/home/sano/work/forge-cht/solver_density_cuda/conjugateWall.cpp:993) は最終ステップを `nStepOuter-1` と判定しています。

   **実測:** `case/53.c3x_vane_cht/run_0156_r2_legB/` には流体の `res_19975.h5` がある一方、固体は `res_solid_5_19974.h5` までで、チェックポイントの累積 `step` も **39999** です。最終壁ダンプは **15000** までです。前半の `run_0155_r2_legA/` は出力間隔4005が20025を割り切るため、この欠陥を回避していました。

   また、[step の復元](/home/sano/work/forge-cht/solver_density_cuda/conjugateWall.cpp:509)は平均バッファが存在する分岐内だけなので、既定の `flux_avg=1` では位相を復元しません。[warmup 判定](/home/sano/work/forge-cht/solver_density_cuda/conjugateWall.cpp:849)も再開後のローカルステップを使い、暖機を繰り返します。

   **対案:** 完了ステップ数と保存位相を統一し、流体・壁・固体を同じチェックポイントとして保存してください。累積ステップは平均バッファと独立に復元し、`flux_avg=1`、`warmup>0`、出力間隔の非倍数での分割を回帰試験に含める必要があります。

2. **Major — G-cons は現在の固体状態を検査していません。**

   [check_cht_balance.py:146](/home/sano/work/forge-cht/solver_density_cuda/tools/check_cht_balance.py:146) は、最新の流体荷重に対して、**初期入力の `wall_profile_5.csv` から復元した固体**を比較しています。保存された `SOLID/T` は読みません。前回 M4 の指摘が残っています。

   **再現:** `case/53.c3x_vane_cht/run_0154_r2_avg42/` で報告の **0.0958%／PASS** は再現しましたが、その固体放熱量は **43430.64 W/m** でした。実際の保存固体から計算すると **43473.52 W/m**、同時刻の流体荷重は **43472.29 W/m**、不釣合いは約 **0.0028%** です。今回の結論は許容内ですが、現行ゲートは対象の固体状態を検証していません。

   **対案:** ソルバ内 `fem2d` 用の入力経路を設け、対応時刻・固体ハッシュを確認した保存温度から Robin 放熱を計算してください。平均荷重と瞬時荷重のどちらを比較したかも出力に残してください。

3. **Major — 外部ループが `--flux` の指定を無視します。**

   [cht_loop.py:319](/home/sano/work/forge-cht/solver_density_cuda/tools/cht_loop.py:319) は `iface_Qf_eff` が存在すると無条件に採用し、[同:341](/home/sano/work/forge-cht/solver_density_cuda/tools/cht_loop.py:341)で選択済みの流束を上書きします。

   **再現:** 該当コードを合成入力で実行すると、`q_eff=10`、`q_compact=3`、`q_recon=7` に対して、どの `--flux` 指定でも渡す荷重は **10** でした。履歴には指定した流束名が記録されるため、A/B の条件と実際の計算が食い違います。

   **対案:** 積分済み荷重を優先する条件を **`a.flux == "q_eff"`** に限定してください。他の指定では選択した流束を使い、平均処理にも同じ条件を適用してください。

4. **Major — 準定常の保存系列が、依然として別の熱流束を評価しています。**

   [cht_wall_series.py:40](/home/sano/work/forge-cht/solver_density_cuda/tools/cht_wall_series.py:40) は既定が `q_compact` で、選択肢にも `q_eff` がありません。前回指摘後も、対象 run の `wall_series.csv` はこの系列でした。

   **実測:** C3X の最終 `q_total` は同 CSV では **42246.9 W/m**、連成の積分済み荷重は **43472.3 W/m** です。Mark II は、登録許容 `drift=osc=0.002` で、保存形の瞬時荷重が **DRIFTING**、移動平均された連成荷重が **STEADY** となりました。この違いは明示が必要です。

   **対案:** `iface_Qf_eff` の直接集計に対応し、瞬時値と連成に使う移動平均値を別列として再判定してください。`ALL STEADY` には対象系列・平均窓・許容を併記してください。

5. **Minor — 現在仕様と結果記録の同期が未完了です。**

   [methods/boundary.md:462](/home/sano/work/forge-cht/methods/boundary.md:462) の例には必須の `tol_solid` がなく、[同:476](/home/sano/work/forge-cht/methods/boundary.md:476)には旧来の `D_f` 凍結説明が残っています。[plans/README.md:30](/home/sano/work/forge-cht/plans/README.md:30)も `draft`・ソルバ内薄肉シェルの説明です。

   また、現在の保存状態による `run_0154` 対 `run_0156` の差は、平均 **+0.001139 K**、局所 RMS **0.001496 K**、最大 **0.003499 K** で、[C3X 台帳](/home/sano/work/forge-cht/case/53.c3x_vane_cht/README.md:780)の値を再現しませんでした。

   **対案:** 設定例・更新式・再開手順・比較成果物を同期し、V2・V3・V6などの未完了事項と、Phase 2 の承認範囲を明確にしてください。

再判定の主要結果は以下です。準定常許容は壁温 `0.0006`、熱量 `0.002` です。

| run | G-if | 壁温最大 | 保存形の瞬時総熱量 | 平均後の連成総熱量 |
|---|---|---|---|---|
| `case/53.c3x_vane_cht/run_0154_r2_avg42/` | **PASS** | **STEADY** | **STEADY** | **STEADY** |
| `case/54.markii_vane_cht/run_0036_r2_avg42/` | **NOT CONVERGED** | **DRIFTING** | **DRIFTING** | **STEADY** |

両翼の流体残差は `check_convergence.py` で **NOT CONVERGED (stalled/plateau)** です。確認した保存場に NaN/Inf はありませんでした。`test_solid_shell.py` と `test_solid_fem2d.py` はともに **PASS (all)**。ゲートの合成試験では、163更新の待機不足と末尾のバッファ再初期化を正しく不合格にしました。前回の4節点・`Df=1e5` の温度精度混用例も、現在は相対残差約 **1.0** として検出されます。

推奨は、**`active` に留め、上記1〜4を修正・回帰検証した後、同一時刻の成果物で V4b を再判定すること**です。ファイルは変更しておらず、提案は **plan 未反映**です。

指摘数: Critical 0 / Major 4 / Minor 1
