**未成立 (incomplete)** — incomplete: 証拠不足: §5.1 #8 の `v3_replay`・`v3b_replay` の基準ファイルと判定記録。手元に float の `run_0015_ns_float_a/`・`run_0016_ns_float_b/` もなく、記載された停止候補・ … / V4・V6 の `STOPPED_AT.files` に列挙された停止時の入口・出口・軸の HDF5/XMF と壁XMF。手元の写しでは欠落しています。V4 の「全場4枚・全場以外76.4 MB」を裏付ける全出力の一覧・サイズ記録も必要です …

# codex レビュー: tooling-settle-convergence-and-monitor (result)

- **plan**: [`plans/active/tooling-settle-convergence-and-monitor.md`](../../plans/active/tooling-settle-convergence-and-monitor.md)
- **stage**: `result` (diff base `48cab832`)
- **date**: 2026-10-16
- **commit**: `b64dbfdd` (feature/card-local-mesh)
- **分類**: `incomplete`
- **backend**: `codex` — 選んだ理由: 既定 (auto: codex を試し、unavailable のときだけ claude)
- **要求したモデル**: `既定 (~/.codex/config.toml)`
- **実効のモデル**: `gpt-6-astra` (codex のログの見出し)
- **CLI の版**: `codex-cli 0.162.1`
- **依頼側との関係**: 不明 (`--requester-model` 未指定)
- **permission_denials**: 対象外 (codex は read-only サンドボックス)
- **実行**: effort `high`, 5.0 min, rc=0
- **試行のログ**: `notes/reviews/2026-10-16-tooling-settle-convergence-and-monitor-result-2.codex.log` (codex: incomplete)
- **判定**: **NO-GO**, 指摘 C0/M3/m1
- **focus**: 前回の result 段 (notes/reviews/2026-10-16-tooling-settle-convergence-and-monitor-result.md、未成立・NO-GO) の指摘 C0/M3/m2 への対応 (§5.1 #7c・#8、V4・V6 の一次データの手元の写し) を確かめ、§6 の各 V の判定が事前登録どおりかを見る

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

判定: NO-GO

前回の修正は進んでいますが、照合の回避経路と `wt_nozzle_v3` の入力検査不足が残っています。V4・V6 の停止時 `SETTLED` は再現できました。`accepted` への移行は保留です。ファイルは変更していません。

1. **Major — `invalid` と判定された基準でも CLI が合格を返す。**  
   [check_settled.py:130](/home/sano/work/forge-cardmesh/solver_density_cuda/tools/check_settled.py:130) は、照合の問題を拒否する対象を `main`・`unconfirmed` に限定しています。`run_0011_ns_v2` の基準をメモリ上で `grid.segment.stage: main → S1` と変更すると、選択処理は `mode: invalid`・問題ありになりますが、通常の CLI は **終了コード0・`VERDICT: SETTLED`** を返しました。基準の SHA256 が本段記録と変わったことも検査されません。  
   **対案:** `invalid` または照合の `problems` がある場合は、判定処理に進まず `INDETERMINATE`・終了コード3にする。この1項目差し替えを CLI の回帰試験に追加してください。

2. **Major — `wt_nozzle_v3` は残差監視に必要な行・列が欠けても合格する。**  
   [settle.py:1199](/home/sano/work/forge-cardmesh/solver_density_cuda/tools/settle.py:1199) の `residual_finite` は、存在する数値だけを調べ、最終行の存在だけで被覆を認めています。`run_0011_ns_v2` の入力と生成器が解決した v3 基準を使い、残差を①最終1行だけ、②`rms_ro` 列だけにすると、いずれも全体が **`SETTLED`** になりました。さらに、残差列が一つもない `step,inner,phase` の1行も見張り単体では合格します。  
   **対案:** 整定の閾値を復活させず、必要な残差列と監視区間の行の被覆・重複・単調性を必須検査にする。欠落は `INDETERMINATE`、非有限は `DIVERGED` にしてください。

3. **Major — 判定結果のハッシュがパス表記と演算の微小差に依存し、同じ計算結果を拒否する。**  
   [stop_rule.py:1309](/home/sano/work/forge-cardmesh/design/forge_design/evaluate/stop_rule.py:1309) は、入力パスや浮動小数点の計算値を含む量の辞書全体をハッシュ化し、[wall_tables.py:714](/home/sano/work/forge-cardmesh/design/forge_design/export/wall_tables.py:714) が完全一致を要求します。  
   同一環境・同一データでも、run を相対パス／絶対パスで渡すだけでハッシュ不一致を再現しました。AWS からの写しには残差統計の丸め差もあり、V4・V6 とも再判定は `SETTLED`・終了コード0なのに、成果物側の収束判定が **「未評価」** になります。  
   **対案:** 入力内容・基準・起動のハッシュ照合は維持し、所在地のパスを量の同一性判定から外す。再計算値は明示した数値許容で比較し、各検査の合否は一致を要求してください。

4. **Minor — 現行仕様と未完了の検証が整理されていない。**  
   [plan:89](/home/sano/work/forge-cardmesh/plans/active/tooling-settle-convergence-and-monitor.md:89) の既定は依然 `wt_nozzle_v1`、§4.9 は v2、実装は v3です。[§5.1 #8](/home/sano/work/forge-cardmesh/plans/active/tooling-settle-convergence-and-monitor.md:149) は完了扱いですが、float の8万 step検証は未成立で、「ユーザの判断待ち」と後続決定も同じ欄に残っています。別条件の小さいノズルの確認も独立した残作業になっていません。  
   **対案:** 現行仕様を v3 に統一し、過去仕様は履歴と明示する。実装完了・検証未成立・適用範囲外を分け、残る検証を§5.1に独立項目として残してください。

確認できた証拠は次のとおりです。

| 対象 | 確認結果 |
|---|---|
| 前回の照合指摘 | 起動記録の欠落・FNV不一致は拒否。ただし指摘1の経路が残る |
| 前回の出力指摘 | XMF の flush・close 検査、HDF5 の flush、完了確認への伝播をコードで確認 |
| 前回の流量指摘 | 出口ゼロ・逆流・符号反転は `NOT_SETTLED`、入口ゼロは `INDETERMINATE`。例外なし |
| V3・v2再生 | 保存記録は両 run とも `VERDICT: A`、t_s＝46,000、区間46,000〜126,000。参考値0.08724の訂正も整合 |
| V4・V6 | 両方ともオンライン13件と再生の verdict 不一致0 |

V4 の `case/71.m55_combustion_wt/run_0011_ns_v2/`、V6 の `case/71.m55_combustion_wt/run_0014_chain_v2/` では、`check_settled.py --no-write` がともに **`VERDICT: SETTLED（閾値は暫定）`**。目的量の区間は36,000〜46,000、停止時観測はそれぞれ46,202・46,171です。全残差履歴138,606行・138,513行と最終場の非有限は0。保存されたメッシュ品質判定はともに **`PASS（AR≤5000、skew≤0.90）`** です。V6 の抽出表と公開表はバイト同一で、成果物の実ファイルも確認できました。run 索引は [case README](/home/sano/work/forge-cardmesh/case/71.m55_combustion_wt/README.md:21) です。

**推奨:** `active` に留め、指摘1→2→3の修正と負例検証、残作業整理、以下の証拠補完を終えてから結果レビューをやり直してください。

証拠不足: §5.1 #8 の `v3_replay`・`v3b_replay` の基準ファイルと判定記録。手元に float の `run_0015_ns_float_a/`・`run_0016_ns_float_b/` もなく、記載された停止候補・変化比を照合できません。

証拠不足: V4・V6 の `STOPPED_AT.files` に列挙された停止時の入口・出口・軸の HDF5/XMF と壁XMF。手元の写しでは欠落しています。V4 の「全場4枚・全場以外76.4 MB」を裏付ける全出力の一覧・サイズ記録も必要です。今回再現できた `SETTLED` と、出力一式・容量条件の検証は別です。

指摘数: Critical 0 / Major 3 / Minor 1
