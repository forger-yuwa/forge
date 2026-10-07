# codex レビュー: tooling-nozzle-upstream-poly-and-throat-sizing (result)

- **plan**: [`plans/active/tooling-nozzle-upstream-poly-and-throat-sizing.md`](../../plans/active/tooling-nozzle-upstream-poly-and-throat-sizing.md)
- **stage**: `result` (diff base `main`)
- **date**: 2026-10-07
- **commit**: `de7a5589` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 6.7 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M2/m1
- **focus**: U0〜U4 の結果 (§9) が事前登録を満たすか。U4 は N0・N1 (+延長・凝縮・N0 の凝縮の延長) とユーザ決定 2 件 (出口 M ± 0.05 %・零に近いオーバーシュートの絶対の幅) の扱い、主張の限定 (N1 − N0 は参考値) が証拠の範囲か。生産採用はユーザ判断
- **extra**: `case/45.isobutane_m6_d155/ns_n012_eval.py`, `case/45.isobutane_m6_d155/ns_n012_cond_ext.py`, `case/45.isobutane_m6_d155/README.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**

壁生成・寸法逆算の実装と、U0・U1・U3・U2d の保存記録は概ね整合しています。  
U4 のユーザ決定２件と「N1 − N0 は参考値」という限定も妥当です。ただし、合否判定の欠陥を直して再評価するまで `accepted` への移動は保留してください。

1. **Major — N0 凝縮延長の判定器が「判定不能」を合格にする**

   根拠: [ns_n012_cond_ext.py:151](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/ns_n012_cond_ext.py:151)。`res_ok` は見出しの存在と `RISING`・`DIVERGED` の不在だけで決まり、`判定不能` を排除していません。

   メモリ上の負例で、`rms_ro` だけを持つ履歴を `check_convergence.analyze` に渡すと、必須４列の欠損を「判定不能」と報告します。しかし延長側の判定式は `residual_no_rising=True` になります。４派生量が `STEADY` なら K 合格にできます。

   また [同ファイル:117](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/ns_n012_cond_ext.py:117) の `judge` は、親子の対応・実効設定・バイナリ・`RUN_RC`・`NAN_SCAN` を確認しません。通常の U4 評価器にある前提検査が延長判定から抜けています。

   **対案:** 判定処理を共通化し、欠損・判定不能・未知の出力は不合格にする。凝縮モーメントを含む必要列と親子の実行条件を確認し、`run_0168_ns_n012_N0_cond/`＋`run_0180_ns_n012_N0_cond_ext/` を再判定してください。既存記録が誤合格だったと断定するものではありませんが、現状の判定器では保証できません。

2. **Major — 「NaN の全走査」が中間 HDF5 を検査していない**

   根拠: [ns_n012.py:855](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/ns_n012.py:855)。残差履歴は全期間を読みますが、HDF5 は `rs[-1]` だけを開きます。[run_ns_n012.sh:201](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/run_ns_n012.sh:201) はこれを「NaN の全走査」としています。

   実関数に模擬入力を渡し、中間 `res_1000.h5` の密度を NaN、最終 `res_2000.h5` を正常にすると、開くのは最終ファイルだけで `CLEAN` を返すことを確認しました。中間場を含む検査を要求するリポジトリ規則を満たしません。

   **対案:** `res_0`・中間・最終の全保存場を検査し、検査済みファイル一覧と最初の異常位置を残す。U4 の対象 run を再走査して保存記録を更新してください。証拠の再検査で済むため、この指摘だけで CFD の再計算は不要です。

3. **Minor — 残作業表と計画索引が現在の到達点を表していない**

   根拠: [plan §5.1:119](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-upstream-poly-and-throat-sizing.md:119) は、訂正済みの「U2 は判定として成り立たない」「扱いを諮る」を現在の残作業として残しています。#6 も U4 実施後の状態に更新されず、[plans/README.md:29](/home/sano/work/forge-integ-1005/plans/README.md:29) は `draft・諮問待ち` のままです。

   **対案:** 過去の FAIL・棄却記録を保持しつつ、現在の状態を「U2d により初期見積もり用途を支持」「U4 は上記再検査後に確定」「生産採用はユーザ判断」に更新する。完了区分を **機能実装の完了** と明記し、P3 の再調査条件への参照も残してください。

検証記録については、次の範囲を確認しました。

| 対象 | 確認結果 |
|---|---|
| U0・U1 | 保存 JSON は明示 `ramp` の互換性、下流のビット同一、形状ゲートの主張と整合 |
| U2〜U2d | U2b の FAIL・U2c の棄却を保持。U2d の往復誤差は出口 −1.83e−6 m、スロート 5.36e−9 m で、限定した許容差内 |
| U3 | 保存壁を独立に読み直し、100,001 点で比較。半径・１階・２階微分の最大差は 1.07e−14・9.97e−14・2.49e−11。登録の許容差内 |
| U4 | 保存 JSON 上、dry の残差は `NOT CONVERGED (stalled/plateau)`、オーバーシュートは `DRIFTING` を保持して絶対幅条件で合格。N0 延長の凝縮４量は `STEADY`、残差は `NOT CONVERGED (stalled/plateau)` |

U4 の数値は [通常評価 JSON](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/ns_n012_eval.json)・[N0 延長評価 JSON](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/ns_n012_cond_ext_run_0180_ns_n012_N0_cond_ext.json)・[run 一覧](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/README.md)との照合です。対象 run の生の CSV・HDF5 はこの環境に無く、全列の再判定・場の再走査は未実施です。保存ログの既存テスト結果と、今回独立に実行した検査も区別しています。

**推奨は、指摘 1 → 2 → 3 の順に修正・再評価し、機能実装完了として `accepted` に移すことです。** 生産採用、流れの同等性、性能改善の確定には広げないでください。ファイルは変更しておらず、本レビューの指摘は plan 未反映です。

指摘数: Critical 0 / Major 2 / Minor 1
