# codex レビュー: discretization-moc-axis-limit-and-corrector (result)

- **plan**: [`plans/active/discretization-moc-axis-limit-and-corrector.md`](../../plans/active/discretization-moc-axis-limit-and-corrector.md)
- **stage**: `result` (diff base `main`)
- **date**: 2026-10-07
- **commit**: `de7a5589` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 6.3 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M1/m2
- **focus**: V0〜V5d・V5′ の結果 (§9) が事前登録を満たすか。V5d は形式上保留 (前提不成立)、NS の V5′ (N2) は dry・凝縮のゲートに合格 (ユーザ決定 2 件を適用)。主張の限定 (N2 − N1 は参考値、性能改善は確定しない) が証拠の範囲か。生産採用は V5d の保留を明示した限定採用としてユーザに諮る
- **extra**: `case/45.isobutane_m6_d155/ns_n012_eval.py`, `case/45.isobutane_m6_d155/moc_v5d_eval.py`, `case/45.isobutane_m6_d155/README.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**

MOC の数式と主要な精度改善は独立試験で再現できました。  
ただし、不合格の MOC を計算準備で拒否する処理と、結果・残作業の記述を直す必要があります。  
V5d は保留のままです。NS の登録ゲート合格を、性能改善や生産採用条件の全面達成には読み替えられません。

1. **Major — MOC の不合格が計算投入経路で強制されない**

   [runner_axismach.py:784](/home/sano/work/forge-integ-1005/design/forge_design/evaluate/runner_axismach.py:784) の `prepare` は、`design_chain` の結果から `moc.gate.pass` を確認せずメッシュ生成へ進みます。`prepare_ns` も同様です。また、[ns_n012.py:571](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/ns_n012.py:571) はゲートを記録しますが、N2 の合格を要求していません。

   メモリ上の負例で `applicable=True・pass=False・n_iter_fail=1` を渡しても、`prepare` がメッシュ生成関数を呼ぶことを確認しました。§4.2 の「反復失敗が1対でもあれば検証・生産は不合格」が徹底されません。

   **対案:** 診断を取得できる設計経路は残し、計算準備・起動の境界で `converge` の合格を必須にしてください。N2 の検査にも追加し、不合格・診断欠損を拒否する試験を置くべきです。今回の N2 が実際に反復失敗していたという指摘ではありません。

2. **Minor — 「悪化の向きの量は無い」は保存数値と矛盾する**

   [plan §9:253](/home/sano/work/forge-integ-1005/plans/active/discretization-moc-axis-limit-and-corrector.md:253) と [case README:112](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/README.md:112) にこの記述があります。しかし [V5d 評価 JSON:426](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/moc_v5d_eval.json:426) の `exit_M_dev` は、腕 B の 4.12e−6 から腕 M の 4.49e−5 へ増えています。差は +4.07e−5 です。

   **対案:** 「参考計算では全量の D＋2SE が許容幅内。ただし出口 Mach 誤差は増加方向。前提不成立のため正式判定は保留」と訂正してください。許容幅内と、悪化方向の差がないことは別です。

3. **Minor — 残作業表と索引が最新の到達点を反映していない**

   [plan §5.1:111](/home/sano/work/forge-integ-1005/plans/active/discretization-moc-axis-limit-and-corrector.md:111) は終了済みの plan レビューを未完了のまま残し、#5 は実施済みの V5d・V5′ を今後の実行手順として記載しています。[plans/README.md:28](/home/sano/work/forge-integ-1005/plans/README.md:28) も V5c 後の待機状態です。

   **対案:** 現在の残作業を「指摘1の修正・検証、本レビューの採否記録、V5d 保留を明示した生産採用判断」に更新してください。`accepted` に移す場合は、§8・§9 に完了区分を **機能実装の完了** と明記し、生産採用と既定値変更は未決定として残してください。

検証の裏付けは次の範囲です。

| 対象 | 確認結果 |
|---|---|
| V0・V2 | 4解像度を独立に再実行。第1段の誤差比・壁誤差・最細区間の次数 **1.9902426** を再現。反復失敗0 |
| V1 | 変更前 `170b0d75` と比較。case/45 の点群609×4・係数259・ノット265、および放射源流の網が完全一致 |
| V3 | 放射源流の再出発試験を再実行して許容差内。case/45 は保存記録との照合 |
| V4・V6 | case/45 を再生成し、第1点の角度差 **0.00544057°**、反復最大34、ゲート合格を再現。形状ゲート・時間比1.48は保存記録と整合 |
| V5・V5b・V5c | 保存 JSON は、それぞれ保留・保留・「出力復元経路だけが原因」を棄却する記録と整合 |
| V5d | **保留（前提不成立）**。腕 M の全温時間幅0.113〜0.126 Kは上限0.1 K超過。`armM_M6` も正式には **判別不能** |
| V5′ | 保存 JSON 上、N1・N2 はユーザ決定2件を適用して登録ゲート合格。ただし下記の元判定を保持する必要あり |

V5′ の対象は `case/45.isobutane_m6_d155/` の dry `run_0166_ns_n012_N1`＋`run_0178_ns_n012_N1_ext`、`run_0167_ns_n012_N2`＋`run_0179_ns_n012_N2_ext`、凝縮 `run_0169_ns_n012_N1_cond`・`run_0170_ns_n012_N2_cond` です。

[保存評価 JSON](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/ns_n012_eval.json) の判定は、残差が **`NOT CONVERGED (stalled/plateau)`**、dry のオーバーシュートが **`DRIFTING`**、凝縮4量が **`STEADY`**。dry は通算80000〜100000、凝縮は14000〜18000の窓です。したがって「N2自身は変更後の登録ゲート合格、N2−N1は参考値、性能改善は未確定」という限定は妥当です。

最新 CFD run の生の CSV・HDF5 はこの環境にないため、全列の独立再判定・場の再走査は未実施です。評価器の SHA256 は保存 JSON と一致しました。恒久索引は [case README の run 一覧](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/README.md:112) です。

**推奨は、指摘1 → 2 → 3の順に修正し、機能実装完了として `accepted` に移すことです。** 生産採用は、V5d の保留を残した限定採用としてユーザに諮ってください。ファイルは変更しておらず、本指摘は plan 未反映です。

指摘数: Critical 0 / Major 1 / Minor 2
