# codex レビュー: time_integration-line-viscous-jacobian-dt (result)

- **plan**: [`plans/active/time_integration-line-viscous-jacobian-dt.md`](../../plans/active/time_integration-line-viscous-jacobian-dt.md)
- **stage**: `result` (diff base `eea98f4b^`)
- **date**: 2026-10-11
- **commit**: `274359b7` (feature/faceh-audit-viscjac-close)
- **codex**: effort `high`, 4.6 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M3/m2
- **focus**: plan を閉じて accepted/ へ移す前の result 段。§6.2 (方向別 dt の A/B)・§6.4 (上限 50 の 3 腕の比較)・§6.6 (全部入りの再実行) の判定と変更ログの主張が、判定の記録 (case/45.isobutane_m6_d155/_band_ab/cold_pair/lvcdt_judge.json・lvcmp_judge.json・lvcmp_rerun_judge.json・*_pregate.json) と判定器のコード (lvcdt_judge.py・lvcmp_judge.py・lvcmp_rerun_judge.py) で裏付けられているか。既定はこの plan では決めない (ユーザ決定 2026-10-11)。§5.1 #3・#4 は未解決のまま残し、再開の条件 (§6.4) を置いて閉じる。ソルバのコードの変更は無い (台本と判定器だけ)
- **extra**: `case/45.isobutane_m6_d155/_band_ab/cold_pair/lvcmp_judge.json`, `case/45.isobutane_m6_d155/_band_ab/cold_pair/lvcmp_rerun_judge.json`, `case/45.isobutane_m6_d155/_band_ab/cold_pair/lvcdt_judge.json`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**

保存済み JSON は、§6.2 の「回避」、§6.4 の「登録幅以上の到達差を検出せず」、§6.6 の「到達時点の差の再現」を支持する。
ただし、本文の定常性・解精度に関する主張と、判定器の証拠検査を修正してから `accepted/` へ移すべき。

指定 diff を取得し、ソルバ・`design/`・`methods/`・`procedures/` に変更がないことを確認した。既定変更もない。プロンプト内の plan・JSON は作業ツリーのファイルと一致する。

検証には制限がある。対象 run の CSV・HDF5・`CONVERGENCE_SEGMENT.txt` はこの作業ツリーに存在せず、生データからの独立再計算はできていない。以下の実測値は保存済み判定記録との照合結果である。

1. **Major — 未定常の到達時点を「止まる状態」「一致」と解釈している。**

   [plan:261](/home/sano/work/forge-faceh/plans/active/time_integration-line-viscous-jacobian-dt.md:261) は L0/L5 を「一致」、[plan:347](/home/sano/work/forge-faceh/plans/active/time_integration-line-viscous-jacobian-dt.md:347) は「止まる状態のずれ」の再現とする。しかし、[`lvcmp_rerun_judge.json`:263](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/_band_ab/cold_pair/lvcmp_rerun_judge.json:263) 以降の記録は次のとおり。

   - L7・L7r の `Q_w`: **DRIFTING**
   - L0 の `theta_r_70/94`、L7r の `theta_r_40`: **TRANSIENT-UNSETTLED**
   - 4 本すべて: **NOT CONVERGED**

   これは §6.5 自身の「止まる状態そのものの再現は言わない」という制限にも反する。

   **対案:** §5.1・§6.4・§6.6・case README を「同じ停止規則の到達時点で観測された差」に統一する。L0/L5 は「到達時点の値の差が 0.03 % 以内」とし、定常解の一致や停止状態の不変性を主張しない。

2. **Major — 残差の大小から解への距離・θ_r の正しさを推論している。**

   [plan:263](/home/sano/work/forge-faceh/plans/active/time_integration-line-viscous-jacobian-dt.md:263) の「解から遠い」「θ_r の正しい値は生産の側が有力」は、この比較では裏付けられない。同じ残差作用素でも、残差ノルムの大小だけでは状態誤差や派生量誤差の順序は保証されない。

   記録で確認できるのは、到達窓の `rms_ro` 中央値が L7/L0 = **3.298 倍**、L7r/L0 = **3.360 倍**という事実まで。構成間の θ_r 差を、そのまま全部入りの反復誤差とも扱えない。

   **対案:** 「到達窓の密度残差中央値が高い」「構成間差が許容目安 0.1 % を超えた」に限定する。解精度の順位は保留する。今回の plan を閉じるために point の追加計算は必要ない。

3. **Major — 判定器が途中の不正な監視値と判定失敗を通す。**

   [`lvcmp_judge.py`:207](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcmp_judge.py:207) は系列の step を検査するが、量の有限性を直接確認するのは到達行だけである。同ファイル [241 行](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcmp_judge.py:241) は準定常判定の失敗を受けても `REACHED` を返す。

   実関数をメモリ上の模擬ファイルに適用し、次を確認した。

   - 途中行の `Q_w = NaN` または `theta_r_40 = Inf` → **REACHED**
   - 準定常判定が `rc=2, overall=None` → **REACHED**
   - 残差 CSV がなくても、判定テキストがあれば収束記録の検査を通過

   `lvcmp_rerun_judge.py` もこの処理を共有する。今回のデータが不正だった証拠ではないが、「ゲート全件合格」だけでは排除できない。

   **対案:** 全監視行の必須量・有限性を検査し、準定常判定の実行失敗は `INVALID` にする。正常な `DRIFTING` 判定は観測として保持する。収束判定を対象 CSV・設定・区間のハッシュと結び付け、保存済み証拠で再判定する。

4. **Minor — `CENSORED` の残差窓が 1 step はみ出す。**

   [`lvcmp_rerun_judge.py`:138](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcmp_rerun_judge.py:138) は `floor_median(..., 150000)` を呼び、step 150000 の `outer_begin` を必須にする。しかし、[ソルバのループ:3482](/home/sano/work/forge-faceh/solver_density_cuda/main.cpp:3482) は 0〜149999 で終了する。

   正常な末尾 CSV を模した入力でも、`20000 / 20001` 行として拒否されることを確認した。今回の `REACHED` 結果には影響しない。

   **対案:** 完走時の残差 step と出力番号の対応を明記し、`CENSORED` では実在する末尾窓を使う。窓の定義変更を plan に記録する。

5. **Minor — 再開条件と影響範囲に古い記述が残る。**

   [plan:282](/home/sano/work/forge-faceh/plans/active/time_integration-line-viscous-jacobian-dt.md:282) の再開条件は「流れ方向の細分化」だが、根拠となる同 [295 行](/home/sano/work/forge-faceh/plans/active/time_integration-line-viscous-jacobian-dt.md:295) は「密度を 1/2」にした粗格子化で、方向が逆である。また、[§7:373](/home/sano/work/forge-faceh/plans/active/time_integration-line-viscous-jacobian-dt.md:373) の「2000 step に限る」は追加された長時間比較を反映していない。

   **対案:** 「粗格子化」に訂正し、§7・`plans/README.md` を短期 A/B と最大 150000 step の到達比較に分けて更新する。

推奨は、**追加 GPU 計算を行わず、上記の記述・判定器を修正し、保存済み証拠を再判定してから閉じる**ことである。優先順は指摘 1・2 → 3・4 → 5。§5.1 #3・#4 を未解決のまま残し、既定を決めない方針は妥当。

`check_plans.py` は **PASS**、Python 6 本の構文確認と shell 4 本の `bash -n` も通過した。ファイル変更は行っていない。

指摘数: Critical 0 / Major 3 / Minor 2
