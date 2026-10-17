# codex レビュー: tooling-nozzle-card-local-mesh (result)

- **plan**: [`plans/accepted/tooling-nozzle-card-local-mesh.md`](../../plans/accepted/tooling-nozzle-card-local-mesh.md)
- **stage**: `result` (diff base `48cab832`)
- **date**: 2026-10-17
- **commit**: `af3c6081` (feature/card-local-mesh)
- **分類**: `ok`
- **backend**: `codex` — 選んだ理由: 既定 (auto: codex を試し、unavailable のときだけ claude)
- **要求したモデル**: `既定 (~/.codex/config.toml)`
- **実効のモデル**: `gpt-6-astra` (codex のログの見出し)
- **CLI の版**: `codex-cli 0.162.1`
- **依頼側との関係**: 不明 (`--requester-model` 未指定)
- **permission_denials**: 対象外 (codex は read-only サンドボックス)
- **実行**: effort `high`, 5.1 min, rc=0
- **試行のログ**: `notes/reviews/2026-10-17-tooling-nozzle-card-local-mesh-result-4.codex.log` (codex: ok)
- **判定**: **NO-GO**, 指摘 C0/M1/m0
- **focus**: 3 回目の result 段 (notes/reviews/2026-10-17-tooling-nozzle-card-local-mesh-result-3.md、NO-GO、M1 成果物の生成元) への対応 (§5.1 #6c2、7bc9a664、run_0014 と試走 A の作り直しと accepted) を確かめ、§8 の完了の範囲で accepted に移してよいかを見る

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

判定: NO-GO

前回の「別 run の成果物で accepted になる」経路は修正されています。
ただし、配布する `nozzle.msh` と計算に使った `nozzle.h5` の対応が保証されず、生成元照合に穴が残っています。

1. **Major — 実際に書き出すメッシュが生成元照合の対象から漏れている**

   [deliverables.py:99](/home/sano/work/forge-cardmesh/design/forge_design/export/deliverables.py:99) が記録する格子は `nozzle.h5` ですが、配布メッシュ・Fluent・寸法検査の入力は [同:365](/home/sano/work/forge-cardmesh/design/forge_design/export/deliverables.py:365) の `nozzle.msh` です。

   両者を結ぶ検査はありません。`prepare_info.json` に記録された msh のハッシュとの不一致も、[同:359](/home/sano/work/forge-cardmesh/design/forge_design/export/deliverables.py:359) では `msh_sha256_matches: false` と記録するだけです。[wt_qualify.py:592](/home/sano/work/forge-cardmesh/design/forge_design/feedback/wt_qualify.py:592) はその不一致を検査しません。

   ファイルを変更せず、次を確認しました。

   - `run_0014` の壁を保ち、径方向配点だけ変えたメッシュをメモリ上で作成。内部節点は最大 **1.918 mm** 動きますが、構造の自己検査・規則検査を通ります。msh の出典不一致も拒否されません。
   - 保存済み MANIFEST の `msh_sha256_matches` をメモリ上で `false` にした負例でも、実際の `qualify()` は **`accepted`・生成元 `MATCH`・`missing: []`** を返しました。

   したがって、格子を作り直して h5 を更新し忘れた場合、古い計算結果と新しい配布メッシュを組み合わせて合格にできます。

   **対案：** §5.1 #6c2 を再オープンし、配布元 msh のハッシュと、計算用 h5 を生成した際の入力 msh の記録を照合してください。不一致・照合不能は `incomplete` にし、同じ壁で内部配点だけ異なる負例を追加してください。

確認できた正常系の結果：

| run（リポジトリルート相対） | 判定区間／停止 step | 再実行した `check_settled.py --no-write` |
|---|---|---|
| `case/71.m55_combustion_wt/run_0014_chain_v2/` | 36000–46000／46171 | `VERDICT: SETTLED` |
| `case/90.qs_trial_a/run_0003_ns_p0/` | 32000–42000／42060 | `VERDICT: SETTLED` |
| `case/90.qs_trial_a/run_0005_ns_p1/` | 20000–30000／30088 | `VERDICT: SETTLED` |

3 本とも全残差履歴・停止時の `res_<step>.h5` の `VALUE/*` に非有限値はなく、ρ・P・T は正でした。

`run_0014` 自身の成果物では **`accepted`・`MATCH`** を再現しました。出口コア M **5.502841913**、η=0.03 の最大偏差 **0.062305 %**・片振幅 **0.057285 %**、壁解像 `PASS`。`run_0010` の古い成果物を渡す負例は **`incomplete`・`ABSENT`** になりました。試走 A の作り直し後の `accepted`・`MATCH` は保存記録で確認しています。索引は [case/71 README](/home/sano/work/forge-cardmesh/case/71.m55_combustion_wt/README.md:23) です。

手組み格子の保存記録は品質 `PASS`、AR 最大 **4849.03**、AR 緩和違反 **0**。V8 の一次記録は記載値を裏付け、V9 の保存データから出口位置差最大 **1.9938×10⁻⁶ r_t** を再確認しました。生成器統合などの後続計画への移管も記録されています。

回帰スイートの再実行は read-only 環境でのキャッシュ・一時ファイル作成制約により完了していません。上記の再判定では描画モジュールのみ代替し、評価・照合は実コードを使用しました。

**推奨：msh と計算用 h5 の対応照合を完成させ、負例の拒否と正常系の accepted 維持を確認してから、§8 の範囲で accepted へ移してください。**

証拠: 足りる

指摘数: Critical 0 / Major 1 / Minor 0
