**未成立 (incomplete)** — incomplete: 証拠不足: `run_0015_ns_float_a/`・`run_0016_ns_float_b/`・`run_0017_ns_float_free130/` の再生に使用した `nozzle.h5`。保存記録のハッシュと一致する入力が必要です。特に … / `run_0017_ns_float_free130/res_130000.h5`。残差の非有限0は確認できましたが、plan の「最終場の非有限0・ρ・T・P > 0」は独立に確認できていません。

# codex レビュー: tooling-settle-convergence-and-monitor (result)

- **plan**: [`plans/active/tooling-settle-convergence-and-monitor.md`](../../plans/active/tooling-settle-convergence-and-monitor.md)
- **stage**: `result` (diff base `48cab832`)
- **date**: 2026-10-16
- **commit**: `8cf2a284` (feature/card-local-mesh)
- **分類**: `incomplete`
- **backend**: `codex` — 選んだ理由: 既定 (auto: codex を試し、unavailable のときだけ claude)
- **要求したモデル**: `既定 (~/.codex/config.toml)`
- **実効のモデル**: `gpt-6-astra` (codex のログの見出し)
- **CLI の版**: `codex-cli 0.162.1`
- **依頼側との関係**: 不明 (`--requester-model` 未指定)
- **permission_denials**: 対象外 (codex は read-only サンドボックス)
- **実行**: effort `high`, 9.2 min, rc=0
- **試行のログ**: `notes/reviews/2026-10-16-tooling-settle-convergence-and-monitor-result-3.codex.log` (codex: incomplete)
- **判定**: **NO-GO**, 指摘 C0/M1/m2
- **focus**: 2 回目の result 段 (notes/reviews/2026-10-16-tooling-settle-convergence-and-monitor-result-2.md、未成立・NO-GO、M1〜M3・m1・証拠不足) への対応 (§5.1 #9〜#13、e92f6b1e・de5a425f・552d0a31) と、#14 の V3(c) の結果 (run_0017、§6) を確かめ、§6 の各 V の判定が事前登録どおりかを見る。証拠の写しは case/71.m55_combustion_wt/_runs/settle_evidence/・_runs/run_0017_ns_float_free130/ と手元の run_0011・0014・0015・0016・0017 の写し

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

判定: NO-GO

前回の照合拒否・残差欠損検査の修正は確認できましたが、通常の再判定で停止時の証拠を上書きし、成果物・受け渡しの照合を壊す問題が残っています。  
V4・V6 の `SETTLED` は再現しました。V3(c) の保存記録は事前登録と整合しますが、独立再生に必要な入力が不足しています。ファイルは変更していません。

1. **Major — 通常の再判定が停止時の証拠を上書きし、その後の照合を失敗させる。**

   [check_settled.py:185](/home/sano/work/forge-cardmesh/solver_density_cuda/tools/check_settled.py:185) は既定で `SETTLE_VERDICT.json` を上書きします。一方、[stop_rule.py:1405](/home/sano/work/forge-cardmesh/design/forge_design/evaluate/stop_rule.py:1405) と [deltastar_loop.py:1250](/home/sano/work/forge-cardmesh/design/forge_design/feedback/deltastar_loop.py:1250) は、停止時に記録した量のハッシュとの完全一致を要求します。今回追加した許容付き比較に到達する前に拒否されます。

   `run_0011_ns_v2`・`run_0014_chain_v2` で、通常 CLI の書き込みをメモリ上に捕捉して再読込を模擬すると、両方とも次を再現しました。

   - 実行前：`match_reeval.matches = true`
   - CLI：終了コード0、`SETTLED`
   - 保存後：`matches = false`、受け渡しの照合も失敗

   したがって、移動した run を通常のコマンドで再判定すると、成果物側は [wall_tables.py:724](/home/sano/work/forge-cardmesh/design/forge_design/export/wall_tables.py:724) で「未評価」になります。

   **対案:** 停止時の再評価を変更されない専用ファイルに保存し、`STOP_RULE_STATE`・受け渡し・成果物はそれを参照する。通常 CLI の最新判定は別に保存する。「通常 CLI → 成果物 → 受け渡し」の連続操作を、相対パス・移動後の run で回帰試験にしてください。

2. **Minor — 完了条件と残る検証の関係が未整理。**

   [plan:156](/home/sano/work/forge-cardmesh/plans/active/tooling-settle-convergence-and-monitor.md:156) の #15 は小さいノズルの確認が未完了です。一方、[§6:163](/home/sano/work/forge-cardmesh/plans/active/tooling-settle-convergence-and-monitor.md:163) の V3(b) は小さいノズルと case/71 の両方を要求し、[§8:221](/home/sano/work/forge-cardmesh/plans/active/tooling-settle-convergence-and-monitor.md:221) は「V1〜V5 がすべて合格」を完了条件としています。

   **対案:** 今回完了とする検証範囲を §8 に明記し、#15 が完了条件に残るのか、後続計画へ移すのかを確定してください。暫定プリセットの採用と一般的な妥当性の確認は区別する必要があります。

3. **Minor — 恒久的な run 索引に訂正が反映されていない。**

   [case README:35](/home/sano/work/forge-cardmesh/case/71.m55_combustion_wt/README.md:35) は V4 の容量を依然 **76.4 MB** と記載しています。今回、保存された一覧から集計すると、節点・境界の HDF5/XMF と流量 CSV は **68,855,382 bytes** でした。また、同 README:39 の float 2 run は投入時の記載だけです。

   **対案:** 容量の集計対象を明記して訂正し、float 2 run の終了結果・再生結果・証拠へのリンクも索引へ反映してください。

確認結果は以下のとおりです。run の接頭辞はすべて `case/71.m55_combustion_wt/` です。

| 対象 | 確認結果 |
|---|---|
| 前回 M1／#9 | `stage: main → S1` の差し替えを通常 CLI が `INDETERMINATE`・終了コード3で拒否 |
| 前回 M2／#10 | 最終1行・1列・列なし・欠測・重複を拒否。非有限は `DIVERGED` |
| 前回 M3／#11 | 原本を保持した比較は成功。ただし指摘1の保存経路が残る |
| V4 `run_0011_ns_v2/` | `check_settled.py --no-write`：**SETTLED**。区間36,000〜46,000、停止時観測46,202。オンライン13件との不一致0 |
| V6 `run_0014_chain_v2/` | 同じく **SETTLED**。区間36,000〜46,000、停止時観測46,171。オンライン13件との不一致0。公開表と停止時抽出表はバイト同一 |
| V4・V6 の出力 | `STOPPED_AT.files` の欠落・サイズ不一致0。全残差履歴と最終場の非有限0、ρ・P・Tは正 |
| V3(c) `run_0017_ns_float_free130/` | 保存記録は **A**：t_s＝32,000、区間32,000〜112,000の41点、最大比0.05832／0.12448／0.01032。基準を現行生成器で再生成すると保存基準と一致。残差390,000行の非有限0 |

V4・V6 の容量条件は訂正後も満たします。V3(c) の保存結果を反証する数値は見つかっていませんが、手元での再生は `nozzle.h5` の欠落により**未成立**でした。これは物理量の超過による結果 B ではありません。

**推奨:** `active` に留め、指摘1の修正と連続操作の回帰試験、完了範囲・run 索引の整理、以下の証拠補完を済ませてから結果レビューを再実施してください。

証拠不足: `run_0015_ns_float_a/`・`run_0016_ns_float_b/`・`run_0017_ns_float_free130/` の再生に使用した `nozzle.h5`。保存記録のハッシュと一致する入力が必要です。特に run_0017 は、この欠落で目的量3種類が `INDETERMINATE` になります。

証拠不足: `run_0017_ns_float_free130/res_130000.h5`。残差の非有限0は確認できましたが、plan の「最終場の非有限0・ρ・T・P > 0」は独立に確認できていません。

指摘数: Critical 0 / Major 1 / Minor 2
