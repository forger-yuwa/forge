# codex レビュー: time_integration-implicit-thermal-jacobian (plan)

- **plan**: [`plans/active/time_integration-implicit-thermal-jacobian.md`](../../plans/active/time_integration-implicit-thermal-jacobian.md)
- **stage**: `plan`
- **date**: 2026-10-10
- **commit**: `3de780a4` (feature/faceh-audit-viscjac-close)
- **codex**: effort `high`, 3.5 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M5/m1
- **focus**: §6.2 (2026-10-10 に追加した事前登録: point 仕上げでの面エンタルピーの精度の A/B) だけを点検する。§4・§6.0 は既に実装・検証済みで対象外。点検してほしいこと: (1) 主判定の共通の評価器を d (面エンタルピーを double で評価) にした理由が妥当か (諮問の記録と引き継ぎメモ notes/sessions/2026-10-10-handoff-face-enthalpy-audit.md は「A の設定の 1 step だけの run」を示唆していた)。(2) 判定の閾値 (10 %、ノイズ = max(再評価の差, 3 本の軌道の幅) の 3 倍、範囲が分かれる条件、床の判定の傾き −0.01 桁/500 step) と 3 本ずつの反復で足りるか。(3) ゲート (step 0 の表示と同じ状態の評価の一致 1e-6、設定の一致、状態の sha256) の抜け。(4) fh_floor.sh と fh_floor_judge.py が §6.2 の文言どおりか (台本は AWS の ~/forge-wallfit/case/45.isobutane_m6_d155 で動く。cold_pair.py・nozzle 系のモジュールはリポジトリの同じディレクトリにある)。出発の run_0354 の末尾の rms_roe は step ごとに 6.51〜6.54 を行き来する。point の 1 step は約 19 ms、2000 step は約 40 s。
- **extra**: `notes/reviews/2026-10-09-line-viscous-faceh-diagnose.md`, `case/45.isobutane_m6_d155/fh_floor.sh`, `case/45.isobutane_m6_d155/fh_floor_judge.py`, `case/45.isobutane_m6_d155/cold_cfl.py`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
§6.2 の A/B と、d を共通評価器にする方針は妥当です。ただし、現行のゲートと判定分岐では、不正な入力や増加中の残差を受理できます。以下を実行前に修正してください。§4・§6.0 は評価対象外です。

1. **Major — 終了失敗・途中の非有限値を判定から排除できない。**

   根拠: [fh_floor.sh:42](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/fh_floor.sh:42) と同52行は終了コードを記録するだけです。[fh_floor_judge.py:100](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/fh_floor_judge.py:100) は最終step番号と最終場の `ro/P/T` を確認しますが、全期間の残差の有限性・行の欠落・重複・`RUN_RC` を確認しません。評価runもstep 0の存在だけです。

   判定器本体をメモリ上の人工入力で実行すると、**step 100 の `rms_roOmega` がNaNでも `gates_ok=true`、`SUPPORT`** になりました。

   対案: 全runの正常終了、軌道の0〜1999行の一意な連続性、評価runの所定行、全 `rms_*` の有限性・非負性を先に検査する。最終場も全保存量の有限性を確認する。失敗時は計算を続けず `INVALID` を出し、評価runのHDF5削除は検査後に移してください。比・対数のゼロ入力も明示的に扱う必要があります。

2. **Major — 軌道6本が指定した `res_40000.h5` から始まる保証がない。**

   根拠: [fh_floor.sh:59](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/fh_floor.sh:59) はSを `res_40000.h5` に固定しますが、軌道の準備は[同40行](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/fh_floor.sh:40)でrunディレクトリを渡します。[cold_cfl.py:115](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/cold_cfl.py:115) は、その時点の**最後のres**を選びます。判定器のSHA照合は評価runだけで、軌道にはありません。

   また、設定照合は `solverConfig.yaml` のみです。`bcondConfig.yaml` などをrun_0183から複製するため、「run_0354と同じ設定」の保証として不足しています。

   対案: 軌道も固定した `_fh_states/s0` を入力とし、6本すべての `parent_res`・SHAを指定ファイルの実ハッシュと照合する。境界条件・解決済み物性・メッシュも照合対象にしてください。step 0のRMS一致は補助検査であり、状態の同一性の代用にはできません。

3. **Major — 「減衰中ではない」を「頭打ち」と判定している。**

   根拠: [plan:295](/home/sano/work/forge-faceh/plans/active/time_integration-implicit-thermal-jacobian.md:295)、[判定器:163](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/fh_floor_judge.py:163)。人工入力で3本とも傾き **+0.1桁/500 step** にすると、結果は `FLAT` でした。1本だけ増加する場合も、全体を頭打ちに分類できます。

   −0.01桁は500 stepで約**2.28%低下**です。これ未満の継続的な低下も、床の証明にはなりません。Bだけを見ており、Aの悪化による相対差も区別できません。

   対案: 各腕について「減衰」「増加」「振動・反復間不一致」「登録した許容内」を分ける。±0.01桁を使うなら、全反復の傾きと複数の末尾窓の水準が許容内であることを要求し、結論は「この窓で有意なドリフトを検出せず」に限定する。既存の1000・1500・2000 stepの共通評価も整合確認に使い、Sからの変化を併記してください。自動延長しない方針は維持できます。

4. **Major — 「支持しない」側がノイズを無視する。**

   根拠: [判定器:143](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/fh_floor_judge.py:143)。支持側は範囲分離とノイズ倍率を要求しますが、否定側は `r > 0.9` だけです。

   人工入力でAを `[10,10,10]`、Bを `[7,9.5,12]` とすると、**r=0.95、noise=0.526でも `NOT_SUPPORT`** になります。これは10%の効果を識別できないデータです。

   対案: 10%境界をノイズ込みの不確かさ範囲がまたぐ場合は、両側とも `INDETERMINATE` にする。各腕3本は大きな効果を探す探索試験には使えますが、標本範囲の3倍は統計的な「3σ」ではありません。同一状態の再評価差は固定状態での実行ばらつきであり、step間の振動も別に記録してください。10%は実用上の検出対象として事前登録する値で、3本から精度保証が得られる値ではありません。

5. **Major — 共通評価で10%未満なら「表示の定義の違いだけ」とするのは過剰な帰属。**

   根拠: [plan:293](/home/sano/work/forge-faceh/plans/active/time_integration-implicit-thermal-jacobian.md:293)、[判定器:178](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/fh_floor_judge.py:178)。

   共通評価で9%改善し、自方式の表示で11%改善した場合も「表示だけ」となります。9%の状態改善は消えていません。さらに、共通評価は2000 stepの単一状態、自方式は末尾500 stepの中央値で、比較対象も異なります。

   対案: この分岐は「自方式では10%以上、共通評価では登録した10%基準に未達」と記述する。「表示だけ」と言うには、共通評価の状態差がノイズ内であることと、**同じ保存状態に対するf/dの差**で表示差を説明できることを確認してください。

6. **Minor — §6.2で要求した区間付きVERDICTの記録が自動化されていない。**

   根拠: [plan:299](/home/sano/work/forge-faceh/plans/active/time_integration-implicit-thermal-jacobian.md:299) は各腕の `check_convergence --segment` を要求します。一方、[run_case.sh:64](/home/sano/work/forge-faceh/solver_density_cuda/tools/run_case.sh:64) は通常判定であり、今回の台本・判定器には区間判定の追加実行・記録がありません。

   対案: 6本の区間0〜1999についてVERDICTを保存し、集計から参照する処理を追加する。短期試験の `NOT CONVERGED` は許容し、`DIVERGED`・判定不能とは分けてください。

**推奨は、dを主評価器とした各腕3本・2000 stepの探索試験を維持し、上記1〜5を修正してから実行することです。**

dの選択自体を戻す必要はありません。既定経路はFP64ビルドでも入力のT・Yをfloatへ変換し、floatを返します。[thermo_d.cuh:1079](/home/sano/work/forge-faceh/solver_density_cuda/cuda_forge/thermo_d.cuh:1079)。dはその丸め経路を避けるため、両腕の状態を測る共通の物差しとして合理的です。ただし「dで評価した離散残差の改善」と限定し、物理精度や収束床の改善まで一般化しないこと。fとの結果不一致も報告してください。

引き継ぎの「Aの設定」は簡便な方法の提案であり、元の諮問は共通評価器をfに固定していません。d採用を今回の事前登録で明示すれば整合します。自己一致の1e−6も保守的な停止基準として維持できますが、不合格だけで「別の残差を読んだ」と断定せず、再評価ばらつきと切り分けてください。

確認範囲は依頼全文、掲載コードと実ファイルの一致、関連実装、シェル構文検査、人工入力による判定器の検査です。AWSの対象runとバイナリはローカルに存在せず、実測値の独立再計算はしていません。ファイル変更・forge起動はなし。**plan未反映**です。

指摘数: Critical 0 / Major 5 / Minor 1
