# codex レビュー: tooling-settle-convergence-and-monitor (result)

- **plan**: [`plans/active/tooling-settle-convergence-and-monitor.md`](../../plans/active/tooling-settle-convergence-and-monitor.md)
- **stage**: `result` (diff base `48cab832`)
- **date**: 2026-10-16
- **commit**: `b283aa48` (feature/card-local-mesh)
- **分類**: `ok`
- **backend**: `codex` — 選んだ理由: 既定 (auto: codex を試し、unavailable のときだけ claude)
- **要求したモデル**: `既定 (~/.codex/config.toml)`
- **実効のモデル**: `gpt-6-astra` (codex のログの見出し)
- **CLI の版**: `codex-cli 0.162.1`
- **依頼側との関係**: 不明 (`--requester-model` 未指定)
- **permission_denials**: 対象外 (codex は read-only サンドボックス)
- **実行**: effort `high`, 5.1 min, rc=0
- **試行のログ**: `notes/reviews/2026-10-16-tooling-settle-convergence-and-monitor-result-5.codex.log` (codex: ok)
- **判定**: **GO-with-changes**, 指摘 C0/M1/m1
- **focus**: 4 回目の result 段 (notes/reviews/2026-10-16-tooling-settle-convergence-and-monitor-result-4.md、NO-GO、M1・試走 A の証拠不足) への対応 (§5.1 #16 の続き、#13 の試走 A の一次記録 case/71.m55_combustion_wt/_runs/qs_trial_a_evidence/ と手元の case/90.qs_trial_a/) を確かめ、§8 の完了の範囲で §6 の各 V の判定が事前登録どおりかを見る

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

判定: GO-with-changes

前回の M1 は修正され、試走 A の主要な一次記録も確認できました。
ただし、旧形式の v3 記録が診断キーの追加だけで成果物側では「未評価」になる回帰が残っています。
以下を直してから `accepted` へ移してください。ファイルは変更していません。

1. **Major — 試走 A の正常な停止記録を、現行の成果物処理が拒否する。**

   [stop_rule.py:1349](/home/sano/work/forge-cardmesh/design/forge_design/evaluate/stop_rule.py:1349) は量の辞書のキー集合に完全一致を要求します。一方、[settle.py:1221](/home/sano/work/forge-cardmesh/solver_density_cuda/tools/settle.py:1221) が追加した `required`・`required_source` は、試走 A の停止時記録にはありません。

   次の両 run で、入力ハッシュ・基準・起動記録は一致し、現行 CLI は終了コード0・`VERDICT: SETTLED` でした。しかし成果物が使う `_settled_convergence` は、次の不一致だけで **`result: 未評価`** を返しました。

   ```
   quantities[0].values: キーが違う
   再計算だけ ['required', 'required_source']、記録だけ []
   ```

   - `case/90.qs_trial_a/run_0003_ns_p0/`：停止42,060
   - `case/90.qs_trial_a/run_0005_ns_p1/`：停止30,088

   [wall_tables.py:730](/home/sano/work/forge-cardmesh/design/forge_design/export/wall_tables.py:730) がこの不一致を成果物へ伝えるため、既存結果から報告を再生成できません。前回の原本保全修正では、この互換性問題は解消しません。

   **対案:** 原本のハッシュ照合と現行の残差被覆検査を維持し、旧形式の `residual_finite` に限った明示的な比較の正規化を追加してください。キーの欠落を一律に無視せず、この2 run を回帰試験へ追加し、正常入力は通過・入力変更時は拒否を確認してください。§5.1にも未解決事項として記録が必要です。

2. **Minor — 試走 A の run 索引が実績と一致していない。**

   [case/90.qs_trial_a/README.md:11](/home/sano/work/forge-cardmesh/case/90.qs_trial_a/README.md:11) と13行目は、p0・p1とも結果欄が空で状態が「予定」です。実際には両方の `REACHED`・`SETTLED`・関門通過を確認できました。

   **対案:** 手元の run 索引に停止step、判定区間、主要成果物、保存状態を記載してください。証拠として保存した原本は保持し、索引からその所在へリンクしてください。

確認できた主要な証拠は次のとおりです。

| 対象 | 確認結果 |
|---|---|
| 前回 M1：V4・V6 | 通常CLIの保存をメモリ上で実行。原本移行後は照合成功。流量CSV変更後・原本喪失時は、成果物の照合と受け渡しの照合がともに拒否 |
| `case/71.m55_combustion_wt/run_0011_ns_v2/` | 現行CLIで **`VERDICT: SETTLED`**。区間36,000〜46,000、停止時観測46,202 |
| `case/71.m55_combustion_wt/run_0014_chain_v2/` | 同じく **`SETTLED`**。区間36,000〜46,000、停止時観測46,171。保存済み関門は `proceed=true`、公開表と停止時抽出表はバイト同一 |
| V3(b)：run_0012・run_0013 | 保存記録はともに **`A`**、46,000〜126,000。planの数値と整合 |
| `case/71.m55_combustion_wt/run_0017_ns_float_free130/` | 既存の `settle_virtual_stop.analyse(write=False)` で **`A` を独立再現**。32,000〜112,000の41点で違反なし。最大比は出口コアM **0.05832**、線 **0.12448**、流れ角 **0.01032** |
| 試走 A：p0・p1 | 現行CLIでともに **`SETTLED`**。判定区間は32,000〜42,000／20,000〜30,000。保存済み関門・現行の受け渡し照合も通過。ただし指摘1あり |

V4・V6・run_0017・試走 A の両 run は、全残差履歴と最終場で非有限0、ρ・P・Tは正でした。[case/71のrun索引](/home/sano/work/forge-cardmesh/case/71.m55_combustion_wt/README.md:35)とも照合しました。

実装差分、既定の `wt_nozzle_v3`、文書の判定方式の分岐、小さいノズル検証の後続計画への移管は整合しています。V1の元のビット一致条件は判定不能のままで、改訂後の限定された検証結果と区別されています。今回、GPU計算と書き込みを伴う試験一式は再実行していません。

**推奨:** 指摘1の互換性修正と実データ回帰試験を先に完了し、次に指摘2の索引を更新してから、§8の限定された完了範囲で `accepted` へ移してください。閾値は暫定のまま維持します。

証拠: 足りる

指摘数: Critical 0 / Major 1 / Minor 1
