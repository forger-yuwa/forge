# codex レビュー: tooling-sern-te-wake-grid (result)

- **plan**: [`plans/active/tooling-sern-te-wake-grid.md`](../../plans/active/tooling-sern-te-wake-grid.md)
- **stage**: `result` (diff base `main`)
- **date**: 2026-10-12
- **commit**: `4257ae4d` (feature/sern-design)
- **codex**: effort `high`, 6.0 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M2/m2
- **focus**: 2026-10-11 の result レビュー (NO-GO、M1〜M3) の指摘への対応 (§5.1 R1〜R3) と、その後の codex diagnose 2 回 (2026-10-12) で追加した R2a〜R2h5・R2c (実 GPU の延長経路 run_1100) が、指摘を塞いでいるか。生産の標準経路で、必須の証拠が欠けた評価が受理・延長・学習に入らないか。移管した項目 (R2h(6)〜(8)・R15・2D plan・R7b) の追跡先と適用範囲の書き方が十分か。status done にして accepted へ移してよいか。
- **extra**: `case/46.sern_design/README.md`, `notes/investigations/2026-10-12-r2c/run_1100_r2c_lsw10_m6_judgment/R2C_SUMMARY.txt`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**  
前回の R1〜R3 と `run_1100` の延長成功は裏付けられます。  
ただし、化学種残差の欠落・残差時系列の切り詰めを受理する経路が残っており、`accepted` への移動は認めません。

1. **Major — TP の化学種残差が欠けても受理・延長できる**

   **根拠:** [`check_convergence.py:171`](/home/sano/work/forge-sern-design/solver_density_cuda/tools/check_convergence.py:171) の必須集合は流体・SST・遷移の列だけです。化学種残差は存在するときだけ検査します。R2h2 の修正は保存場の `roY*` には適用されていますが、残差には及んでいません。

   `run_1100` の実効 config は２化学種を要求します。保存された残差 fixture から `rms_roY0`・`rms_roY1` を削除し、実際の CSV ローダと `residual_health` に通すと、**`ok=True`、`unverifiable=False`** でした。

   他の証拠を保存結果に固定した接続試験では、１回目相当は延長可、２回目相当は **`GATES PASS`**。対照の `rms_roOmega` 欠落は正しく拒否されました。したがって、修正した入力検査の対象集合に穴が残っています。

   **対案:** 実効 config から化学種・トレーサの必須残差列も決め、欠落を `JUDGMENT_UNVERIFIABLE` にしてください。単独列欠落・全化学種列欠落について、受理と延長の両方を拒否する試験を追加すべきです。

2. **Major — 残差が１行だけでも「プラトー」として受理できる**

   **根拠:** [`check_convergence.py:189`](/home/sano/work/forge-sern-design/solver_density_cuda/tools/check_convergence.py:189) は空データを拒否しますが、時系列の必要範囲・欠測・終了位置を検査しません。[`sern_gates.py:179`](/home/sano/work/forge-sern-design/design/forge_design/metrics/sern_gates.py:179) も `last_step` を保存するだけで、判定段との照合に使っていません。

   同じ fixture を **先頭の `step=0` の１行だけ**、または **末尾の `step=19999` の１行だけ**にすると、いずれも次の結果でした。

   - 残差：`ok=True`、`unverifiable=False`、`NOT CONVERGED (stalled/plateau)`
   - 他の証拠を固定した１回目相当：延長可
   - 同じく２回目相当：`GATES PASS`

   最終場や床カウンタが完全でも、欠けた残差区間の NaN・上昇を検査したことにはなりません。

   **対案:** 判定段の設定と残差出力の時点規約から、必要な step・phase の範囲を固定し、欠測・重複・切り詰めを判定不能にしてください。残差は今回 `0…19999`、場は最終 `20000` なので、その対応を明示して照合する必要があります。`require_residual_pass=True` への変更では解決しません。

3. **Minor — run 索引が未達の `run_1097` まで合格と読める**

   **根拠:** [`case/46.sern_design/README.md:408`](/home/sano/work/forge-sern-design/case/46.sern_design/README.md:408) は `run_1097` を含む行に「全 run … STEADY・窓 OK」と記載しています。

   保存系列の再判定では、`run_1097_tewake_lsw10_m6_cont20k/` は **力の窓 `UNSETTLED`、温度も `UNSETTLED`**。延長後の `run_1099_tewake_lsw10_m6_cont40k/` が両方 `PASS` です。

   **対案:** 「1096・1098・1099 が合格、1097 は窓・温度未達で延長」と明記してください。恒久索引でも plan と同じ区別が必要です。

4. **Minor —「既存 YAML・run には影響しない」は適用範囲を過小に書いている**

   **根拠:** [`plan:174`](/home/sano/work/forge-sern-design/plans/active/tooling-sern-te-wake-grid.md:174) は、blend の既定値が０であることから変更全体を無影響としています。しかし R2f〜R2h は、blend の指定にかかわらず、必須量・実効 config・力の部品が欠ける既存評価の判定を変えます。[`evaluate_gates:813`](/home/sano/work/forge-sern-design/design/forge_design/metrics/sern_gates.py:813)

   **対案:** 「既定の格子形状は不変」と「受理検査は既存経路にも適用」を分け、再判定で不合格になり得る入力を明記してください。

確認できた修正・結果は以下です。

| 対象 | 今回の確認 |
|---|---|
| R1・R2b | 学習条件が床記録欠落・`judgment=UNSETTLED` を拒否 |
| R3 | 同一変換器識別の下で `L_sw=0.8/1.0` は同じ方式 ID、旧０格子・band は別 ID |
| 既定格子 | 小型構成で導入前版と座標がバイト一致、接続・境界面も一致 |
| `run_1100` | 保存系列の再判定で１回目 `UNSETTLED`、２回目 `PASS` を再現 |
| 移管項目 | 共有ツール #12〜#14、親 R15・R7b、2D plan に追跡先あり。固定形状への限定も記載あり |

`case/46.sern_design/run_1100_r2c_lsw10_m6_judgment/` の[保存記録](/home/sano/work/forge-sern-design/notes/investigations/2026-10-12-r2c/run_1100_r2c_lsw10_m6_judgment/)では、２回目の力・温度 CSV に `check_quasisteady.py` を再実行し、両方 **`OVERALL: ALL STEADY`**。床記録の単体監査も **`PASS`、区間 `(10000, 20000]`、10000更新、T/ρ/P 各０件**でした。残差の保存判定は **`NOT CONVERGED (stalled/plateau)`**で、`rms_roK` は `falling`。収束や全列プラトーとは認定しません。run 索引は [README:409](/home/sano/work/forge-sern-design/case/46.sern_design/README.md:409) です。

`run_1086 − run_1084` の末尾窓平均差も、C_T −5.26786e−5、C_T_with_shear −3.56800e−4、C_L −1.98803e−4、C_M +6.34811e−3 と plan を再現しました。局所化格子の結果を撤回する根拠はありません。

**推奨は、局所化方式と R2c の結果を維持し、指摘１→２の残差完全性検査を塞いでから再レビューすることです。** この２点は現在の生産経路を通るため、移管だけで完了扱いにはできません。

レビュー対象は `4257ae4d`。欠落試験はメモリ内で実施しました。AWS の実 HDF5・完全な残差系列・GPU 回帰は独立再実行していません。ファイルは変更しておらず、指摘は **plan 未反映**です。

指摘数: Critical 0 / Major 2 / Minor 2
