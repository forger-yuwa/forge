# codex レビュー: tooling-nozzle-card-local-mesh (result)

- **plan**: [`plans/accepted/tooling-nozzle-card-local-mesh.md`](../../plans/accepted/tooling-nozzle-card-local-mesh.md)
- **stage**: `result` (diff base `48cab832`)
- **date**: 2026-10-17
- **commit**: `519ef3d1` (feature/card-local-mesh)
- **分類**: `ok`
- **backend**: `codex` — 選んだ理由: 既定 (auto: codex を試し、unavailable のときだけ claude)
- **要求したモデル**: `既定 (~/.codex/config.toml)`
- **実効のモデル**: `gpt-6-astra` (codex のログの見出し)
- **CLI の版**: `codex-cli 0.162.1`
- **依頼側との関係**: 不明 (`--requester-model` 未指定)
- **permission_denials**: 対象外 (codex は read-only サンドボックス)
- **実行**: effort `high`, 4.5 min, rc=0
- **試行のログ**: `notes/reviews/2026-10-17-tooling-nozzle-card-local-mesh-result-5.codex.log` (codex: ok)
- **判定**: **GO-with-changes**, 指摘 C0/M0/m1
- **focus**: 4 回目の result 段 (notes/reviews/2026-10-17-tooling-nozzle-card-local-mesh-result-4.md、NO-GO、msh と h5 の対応) への対応 (§5.1 #6c2 の再オープンの分: msh の照合と h5 の座標・接続の照合、f7d3c8ff) を確かめ、§8 の完了の範囲で accepted に移してよいかを見る。試走 A の p1 の判定し直し (accepted・MATCH) は _runs/qs_trial_a_evidence/run_0005_ns_p1/WT_QUALIFICATION_msh_2026-10-17.json

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

判定: GO-with-changes

前回の Major（配布 msh と計算用 h5 の対応）は解消されています。
正常系の合格維持と不一致の拒否を確認しました。ただし、終端 C⁺ の変更に追随していない回帰試験が 1 件残っています。

1. **Minor — CONTUR の回帰試験が、新しい設計壁を旧 run の参照値と比較して失敗する**

   [run_two_coef_initializer_tests.py:33](/home/sano/work/forge-cardmesh/design/tests/run_two_coef_initializer_tests.py:33) は現在の MOC で設計を作り、[同:75](/home/sano/work/forge-cardmesh/design/tests/run_two_coef_initializer_tests.py:75) 以降で旧 `run_0003_ns_cont` 由来の表と比較しています。

   再実行では、項目 (4) の最大差が **0.00245 mm > 許容 0.002 mm** となり、終了コード 1。他の項目は合格しました。同じ環境で `_design_cplus` だけを `eabbab97` の親コミットの実装に戻した対照では、最大差 **0.000419 mm、ALL PASS** でした。

   終端 C⁺ の意図した変更と、過去の参照値との比較が混在しています。CONTUR の式・導関数の不具合を示す結果ではありませんが、現行 HEAD で試験合格とは扱えません。

   **対案：** 過去 run との比較には、その run に対応する固定の設計形状を使ってください。現行 MOC を通す検証は別に維持し、許容値を単に緩めず、この試験の合格を確認して §5.1・§6 の記録を更新してください。

確認した根拠：

- [wt_qualify.py:600](/home/sano/work/forge-cardmesh/design/forge_design/feedback/wt_qualify.py:600) の msh 照合と、同ファイルの `_h5_binding` により、準備記録・run・配布 msh、および h5 の座標・接続が照合されます。
- `run_0014` の正常系は **`accepted`・`MATCH`**。メモリ上で不一致フラグや h5 座標ハッシュの不一致を注入すると、実際の `qualify()` は **`incomplete`・`MISMATCH`**。準備記録の欠落も `UNCHECKABLE` になりました。
- 再実行した `check_settled.py --no-write` の結果：

| run（リポジトリルート相対） | 判定区間／停止 step | VERDICT |
|---|---|---|
| `case/71.m55_combustion_wt/run_0014_chain_v2/` | 36000–46000／46171 | SETTLED |
| `case/90.qs_trial_a/run_0003_ns_p0/` | 32000–42000／42060 | SETTLED |
| `case/90.qs_trial_a/run_0005_ns_p1/` | 20000–30000／30088 | SETTLED |

3 本とも全残差履歴と停止時の `res_<step>.h5` の `VALUE/*` に非有限値はなく、ρ・P・T は正でした。

`run_0014` の出口コア M **5.502842**、η=0.03 の最大偏差 **0.0623 %**・片振幅 **0.0573 %**、壁解像 `PASS` を再現しました。試走 A p1 の **`accepted`・`MATCH`** は、指定された `WT_QUALIFICATION_msh_2026-10-17.json` でも確認しました。run 索引は [case/71 README](/home/sano/work/forge-cardmesh/case/71.m55_combustion_wt/README.md:23) です。

格子の保存記録は品質 `PASS`、AR 最大 **4849.03**、AR 緩和違反 **0**。V8 の一次記録は記載値を裏付け、V9 の保存データから出口位置差最大 **1.9938×10⁻⁶ r_t** を再確認しました。生成器統合などの後続計画への移管も記録されています。

**推奨：上記の回帰試験の参照条件を修正して合格を確認し、レビュー記録・移動に伴う索引とリンクを同期したうえで、§8 の範囲で accepted へ移してください。**

ファイル変更はありません。試験では Numba のディスクキャッシュを無効化し、設計判定では描画モジュールだけを代替しました。回帰スイート全体の再実行はしていません。

証拠: 足りる

指摘数: Critical 0 / Major 0 / Minor 1
