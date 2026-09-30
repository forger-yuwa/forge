# codex レビュー: boundary-cht-conjugate-flat-plate (plan)

- **plan**: [`plans/active/boundary-cht-conjugate-flat-plate.md`](../../plans/active/boundary-cht-conjugate-flat-plate.md)
- **stage**: `plan`
- **date**: 2026-10-01
- **commit**: `b1573788` (feature/cht-conjugate-benchmarks)
- **codex**: effort `high`, 3.5 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M2/m2

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
B1 を最初に行う方針は妥当です。ただし、診断量の採取・判別方法と、C 全体の最終合否条件を実装前に固定してください。  
現時点では slip を真因とする根拠はなく、ソルバ修正への着手は支持しません。

目的の同定と計画の分割は適切です。発注元は **A の限定結果で閉鎖**され、C は明示的に未完了として移されています（[発注元:27](/home/sano/work/forge-cht/plans/accepted/boundary-cht-conjugate-benchmarks.md:27)）。重複実装ではありません。`warmup` による更新停止も累積 step の実装と整合します。今回の適用範囲は node・FP64・平面・定常であり、cell・float32・周期・軸対称への保証は含めない扱いが正しいです。

確認した実測は以下です。対象は `case/65.conjugate_flat_plate/run_0007_c1_n64/`、索引は [case README](/home/sano/work/forge-cht/case/65.conjugate_flat_plate/README.md:11) です。

- `conjugate_history.csv` を再集計：末尾400更新の `res_abs_Wm2` は **13.714～741.761、平均182.779**。計画の記載を確認できました。
- 保存済み流体判定は **`NOT CONVERGED`**。ただし、手元での再実行は **`NO residual_history.csv`** で、全残差からの独立再判定はできません。
- 保存された評価窓内77節点の温度・熱流束系列を正式ツールで再判定：ともに **`OVERALL: ALL STEADY`**。窓外の定常性は示しません。
- メッシュ品質の再実行：**`VERDICT: PASS`**、最大 AR **149.7**、最大 skewness **0.000**。

1. **Major — B1 の観測・判別手順が、現状では再現可能な事前登録になっていません。**

   **根拠:** [対象 plan:47](/home/sano/work/forge-cht/plans/active/boundary-cht-conjugate-flat-plate.md:47) の「細かく」「代表値」「局所変動も減衰」が未定義です。さらに、B は [conjugateWall.cpp:1081](/home/sano/work/forge-cht/solver_density_cuda/conjugateWall.cpp:1081) で戻るため、更新処理内の `conjugate_history.csv`・界面節点ログも出ません（[同:1056](/home/sano/work/forge-cht/solver_density_cuda/conjugateWall.cpp:1056)）。

   `check_quasisteady.py <run>` の既定対象は `shock/asym/machmax/pmax` で、必要な界面熱流束や局所 `P/T/U_n` ではありません（[check_quasisteady.py:381](/home/sano/work/forge-cht/solver_density_cuda/tools/check_quasisteady.py:381)）。絶対圧力を既定の相対尺度で評価すれば、101325 Pa に対する数 Pa の振動を容易に `STEADY` と判定します。

   **対案:** A/B 共通で壁 HDF5 の `iface_q_eff` を採取してください。この出力は更新停止中も独立に計算されます（[output.cpp:299](/home/sano/work/forge-cht/solver_density_cuda/output/output.cpp:299)）。併せて以下を登録します。

   - 連成間隔50に対する採取位相、出力間隔、前後縁の対象節点集合と法線。
   - 残差は `outer_end` の各20000 step区間の算術平均、横ばいは例えば区間平均の `max/min ≤ 1.1`。
   - 局所変動の振幅定義、B/A の減衰比、絶対尺度。対象系列を `--series-csv` と明示的な閾値で判定する。
   - A の再現確認には、低下桁数だけでなく元 run の残差水準・変動位置との比較を加える。

2. **Major — 発注元から、どの格子で最終合否を決めるかという未解決事項を引き継いでいます。**

   **根拠:** [対象 plan:72](/home/sano/work/forge-cht/plans/active/boundary-cht-conjugate-flat-plate.md:72) は発注元 §6 をそのまま参照しますが、そこには「各3水準」「格子間差が減る」とあり、合格を要求する格子が明記されていません（[発注元:150](/home/sano/work/forge-cht/plans/accepted/boundary-cht-conjugate-benchmarks.md:150)）。この曖昧さは発注元自身も [同:140](/home/sano/work/forge-cht/plans/accepted/boundary-cht-conjugate-benchmarks.md:140) で認めています。細格子ほど G-if が悪化した今回、判定対象の事後選択は結論を変えます。

   **対案:** 今回は **C1/C2 × n16/32/64 の6本すべてに、前提ゲートと主判定の合格を要求する**と明記してください。B1 は診断用であり、この6本の代替にはしません。C2 の軸方向伝導効果、不確かさ、固体指標の準定常性も再判定項目として残作業表に列挙してください。

3. **Minor — 同時刻 restart の成立確認を、B1 投入前の作業として明記してください。**

   **根拠:** 指定された `res_600000.h5`・`conjugate_state_5.h5` と元の残差 CSV は手元にありません。また、固体状態の読込みは `st.u` を復元しますが、壁温は別経路の `wallProfile` です（[conjugateWall.cpp:645](/home/sano/work/forge-cht/solver_density_cuda/conjugateWall.cpp:645)、[methods/boundary.md:516](/home/sano/work/forge-cht/methods/boundary.md:516)）。固体状態をコピーするだけでは壁温復元を保証しません。

   **対案:** §5.1 の B1 前に、一次データ回収、バイナリ照合、累積 step・メッシュ・保存量 dtype の確認、`wallProfile: 1` の設定、固体界面温度と壁温分布の照合を追加してください。短い再開確認で **A/B の初期壁温が同じこと、B の壁温が固定されること**を測ってから長時間 run に進めます。

4. **Minor — ④の丸めは、判定値だけでなく固体更新にも入る可能性があります。**

   **根拠:** 非界面の診断残差は `Au−b`（[conjugateWall.cpp:953](/home/sano/work/forge-cht/solver_density_cuda/conjugateWall.cpp:953)）ですが、更新補正の右辺も `−(Au−rhs)` です（[同:975](/home/sano/work/forge-cht/solver_density_cuda/conjugateWall.cpp:975)）。非界面では同じ相殺を含みます。「診断の丸めによる偽 FAIL」だけを調べると、更新への影響を見落とします。ただし、このコードだけから大きな界面振動の原因とは断定できません。

   **対案:** 同一の保存温度を固定し、現行式・温度差形式・高精度参照で、非界面残差に加えて補正量 `Δu` も比較してください。一様310 Kだけでなく、実際の非一様固体場を使います。残差評価式を変更するなら、判定用と更新用の両経路を検証対象にします。

**推奨は、B1 を補完して実行することです。** 投入前の優先順は、③の再開成立確認、①の観測・判別仕様固定、②の完了条件明文化です。④は残差評価の変更に先立って行ってください。B1 の結果が出るまでは `Df_scale`・slip・CFL・閾値を変更しない方針を維持します。

ファイルは変更していません。以上の提案は **plan 未反映**です。

指摘数: Critical 0 / Major 2 / Minor 2
