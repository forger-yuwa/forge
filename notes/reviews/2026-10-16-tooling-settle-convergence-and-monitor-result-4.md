**未成立 (incomplete)** — incomplete: 証拠不足: §8の完了条件に含まれる `case/90.qs_trial_a/` のp0・p1。手元には[手順書の結果表](/home/sano/work/forge-cardmesh/procedures/design-chain-quickstar …

# codex レビュー: tooling-settle-convergence-and-monitor (result)

- **plan**: [`plans/active/tooling-settle-convergence-and-monitor.md`](../../plans/active/tooling-settle-convergence-and-monitor.md)
- **stage**: `result` (diff base `48cab832`)
- **date**: 2026-10-16
- **commit**: `61fc3d59` (feature/card-local-mesh)
- **分類**: `incomplete`
- **backend**: `codex` — 選んだ理由: 既定 (auto: codex を試し、unavailable のときだけ claude)
- **要求したモデル**: `既定 (~/.codex/config.toml)`
- **実効のモデル**: `gpt-6-astra` (codex のログの見出し)
- **CLI の版**: `codex-cli 0.162.1`
- **依頼側との関係**: 不明 (`--requester-model` 未指定)
- **permission_denials**: 対象外 (codex は read-only サンドボックス)
- **実行**: effort `high`, 6.1 min, rc=0
- **試行のログ**: `notes/reviews/2026-10-16-tooling-settle-convergence-and-monitor-result-4.codex.log` (codex: incomplete)
- **判定**: **NO-GO**, 指摘 C0/M1/m0
- **focus**: 3 回目の result 段 (notes/reviews/2026-10-16-tooling-settle-convergence-and-monitor-result-3.md、NO-GO、M1・m1・m2・証拠不足) への対応 (§5.1 #16 = 61fc3d59、#17・§8 の完了の範囲 = ユーザ決定、README の訂正、#13 の追加の証拠 nozzle.h5・run_0017 の res_130000.h5) を確かめ、§8 の完了の範囲で §6 の各 V の判定が事前登録どおりかを見る

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

判定: NO-GO

追加されたデータから V3(c) の `VERDICT: A` を独立再現でき、前回の証拠不足は解消しました。
ただし、停止記録の修正により、旧形式の run で入力変更後の照合まで通る問題が生じています。
今回の完了条件に含まれる試走 A は、一次記録が手元になく確認できません。ファイルは変更していません。

1. **Major — 旧形式の互換処理が、停止時の証拠との不一致を通常の再判定で解消してしまう。**

   [stop_rule.py:1528](/home/sano/work/forge-cardmesh/design/forge_design/evaluate/stop_rule.py:1528) は、停止時のファイルハッシュが一致しなくても、判定・step・量ごとの合否などの要約が同じなら最新の判定を代用します。[deltastar_loop.py:1265](/home/sano/work/forge-cardmesh/design/forge_design/feedback/deltastar_loop.py:1265) は、それを現在の入力から再計算した結果と比較します。停止時の入力ハッシュ・量の値との照合が失われます。

   `case/71.m55_combustion_wt/run_0011_ns_v2/` で、流量CSVの流入・流出の絶対値を各1 kg/s増やす変更と、CLIの保存を**メモリ上だけで模擬**しました。ファイルサイズ、起動記録、停止記録、停止時の全場は維持しています。

   | 操作 | 結果 |
   |---|---|
   | 入力変更後、CLI実行前 | `SETTLED` だが `match_reeval.matches = false`。入力ハッシュの不一致を検出 |
   | 通常の `check_settled.py` | 終了コード0、`SETTLED` |
   | 保存後の照合 | `matches = true`、`form = overwritten`、`exact = false` |
   | 受け渡しの照合 | `_v2_binding.pass = true`、問題0 |

   これは丸め差への許容ではなく、変更後の入力を変更後の記録と照合しているためです。[plan:157](/home/sano/work/forge-cardmesh/plans/active/tooling-settle-convergence-and-monitor.md:157) に弱点を記載していても、§4.9の停止時記録との照合契約は満たしません。

   **対案:** 旧形式は、停止時のハッシュと一致する原本を通常CLIの上書き前に専用記録へ移行してください。原本が失われた場合は照合不能として拒否します。「入力変更→通常CLI→成果物・受け渡し」を負例に追加し、§5.1 #16を未完了に戻してください。新形式の専用ファイル方式は維持できます。

確認した結果を示します。run の接頭辞は `case/71.m55_combustion_wt/` です。

| 対象 | 確認結果 |
|---|---|
| #16の通常操作 | V4・V6とも、入力を変えなければCLI後の照合は成功。新旧両形式をメモリ上で確認 |
| V4 `run_0011_ns_v2/` | 現行CLIで `VERDICT: SETTLED`。区間36,000〜46,000、停止時観測46,202 |
| V6 `run_0014_chain_v2/` | 同じく `SETTLED`。区間36,000〜46,000、停止時観測46,171。保存済み関門は `proceed=true`、公開表と停止時抽出表はバイト同一 |
| V3(c) `run_0017_ns_float_free130/` | 既存の `settle_virtual_stop.analyse(write=False)` で **A**。t_s＝32,000、検証区間32,000〜112,000。最大比は出口コアM **0.05832**、線 **0.1245**、流れ角 **0.01032** |
| float の `run_0015`・`run_0016` | 再生はともに t_s＝32,000、必要期間の不足で **未成立**。記載と整合 |
| 非有限・出力 | V4・V6・run_0017の全残差履歴と最終場で非有限0、ρ・P・Tは正。V4・V6の `STOPPED_AT.files` の欠落・サイズ不一致0 |
| 前回の文書指摘 | [run索引](/home/sano/work/forge-cardmesh/case/71.m55_combustion_wt/README.md:35) の容量訂正・float結果を確認。V4の監視出力は再集計でも **68,855,382 bytes** |

小さいノズルの検証を後続計画へ移した決定は、[今回の完了範囲](/home/sano/work/forge-cardmesh/plans/active/tooling-settle-convergence-and-monitor.md:226) と[移管先の#11](/home/sano/work/forge-cardmesh/plans/active/tooling-design-orchestrator.md:224)で整合しています。既定の `wt_nozzle_v3` と暫定閾値の扱いも実装・文書で確認しました。

**推奨:** `active` に留め、指摘1の互換処理を修正し、試走 A の一次記録を補完してから再レビューしてください。

証拠不足: §8の完了条件に含まれる `case/90.qs_trial_a/` のp0・p1。手元には[手順書の結果表](/home/sano/work/forge-cardmesh/procedures/design-chain-quickstart.md:755)だけがあり、`REACHED`・`SETTLED`・関門通過・非有限0を裏付ける停止記録、判定基準、VERDICT、残差履歴、判定に使用した観測データと最終場を確認できません。

指摘数: Critical 0 / Major 1 / Minor 0
