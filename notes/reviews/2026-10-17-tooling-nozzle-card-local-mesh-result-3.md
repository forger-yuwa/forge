# codex レビュー: tooling-nozzle-card-local-mesh (result)

- **plan**: [`plans/accepted/tooling-nozzle-card-local-mesh.md`](../../plans/accepted/tooling-nozzle-card-local-mesh.md)
- **stage**: `result` (diff base `48cab832`)
- **date**: 2026-10-17
- **commit**: `f8ecc12d` (feature/card-local-mesh)
- **分類**: `ok`
- **backend**: `codex` — 選んだ理由: 既定 (auto: codex を試し、unavailable のときだけ claude)
- **要求したモデル**: `既定 (~/.codex/config.toml)`
- **実効のモデル**: `gpt-6-astra` (codex のログの見出し)
- **CLI の版**: `codex-cli 0.162.1`
- **依頼側との関係**: 不明 (`--requester-model` 未指定)
- **permission_denials**: 対象外 (codex は read-only サンドボックス)
- **実行**: effort `high`, 4.1 min, rc=0
- **試行のログ**: `notes/reviews/2026-10-17-tooling-nozzle-card-local-mesh-result-3.codex.log` (codex: ok)
- **判定**: **NO-GO**, 指摘 C0/M1/m0
- **focus**: 2 回目の result 段 (notes/reviews/2026-10-17-tooling-nozzle-card-local-mesh-result-2.md、NO-GO、M1・M2・m1) への対応 (§5.1 #7a・#7b の再オープンの分、methods の断定の訂正) を確かめ、§8 の完了の範囲で accepted に移してよいかを見る

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

判定: NO-GO

前回の M1・M2 と文書の断定は修正され、負例の拒否も確認できました。
ただし、別 run の成果物を組み合わせて `accepted` にできる経路が残っています。
計算結果そのものの合格を否定する指摘ではなく、設計判定ツールの誤受理です。

1. **Major — 成果物と評価対象の run・停止 step が結び付いていない**

   [wt_qualify.py:595](/home/sano/work/forge-cardmesh/design/forge_design/feedback/wt_qualify.py:595) は成果物内のファイルを `MANIFEST.json` のハッシュと照合しますが、その生成元を現在の評価対象と照合しません。[同:607](/home/sano/work/forge-cardmesh/design/forge_design/feedback/wt_qualify.py:607) は run 名の不一致も警告だけで通します。

   ファイルを変更せず、実データで次を再現しました。

   - 評価対象：`case/71.m55_combustion_wt/run_0014_chain_v2/`、停止 step **46171**
   - カード：`case/71.m55_combustion_wt/card_m55_3coef.yaml`
   - 渡した成果物：`case/71.m55_combustion_wt/run_0010_ns_b0/deliverables/`
   - 結果：**`accepted`、`missing: []`、成果物・寸法とも `PASS`**

   渡した成果物の [MANIFEST.json:912](/home/sano/work/forge-cardmesh/case/71.m55_combustion_wt/run_0010_ns_b0/deliverables/MANIFEST.json:912) は、CFD 表の生成元を **run_0010、step 86000、status `diagnostic`** と記録しています。両 run の `wall_repr.json`・`nozzle.h5` の SHA-256 も異なります。単なる保存場所の移動ではありません。

   **対案：** 成果物の生成元の壁・メッシュ・停止場・設定のハッシュと step を記録し、評価対象との一致を必須にしてください。別ディレクトリへのコピーは内容の一致で認め、不一致・照合不能は `incomplete` にします。今回の組み合わせを負例に追加し、§5.1 #6b に残作業として記録してください。

確認できた修正・検証結果は以下です。

- **前回 M1：解消。** 通常の run_0014 は停止時の記録との照合が通ります。η=0.06 の評価、基準の読み取り値だけの変更、基準内容の変更、停止時判定の欠落は、それぞれ `incomplete` または `INDETERMINATE` になります。
- **前回 M2：解消。** `--groups` 省略時に未出力の壁を境界設定へメモリ上で追加すると、`INDETERMINATE`・終了コード2。存在しない step も拒否し、通常入力は `wall/3@46171`・`PASS` でした。
- **前回 Minor：解消。** [methods/design/overview.md:828](/home/sano/work/forge-cardmesh/methods/design/overview.md:828) は壁補正の滑らかさと対象 run の実測に限定した記述になっています。

`check_settled.py --no-write` の再実行結果：

| run（リポジトリルート相対） | 判定区間／停止観測 | VERDICT |
|---|---|---|
| `case/71.m55_combustion_wt/run_0014_chain_v2/` | 36000–46000／46171 | `SETTLED` |
| `case/90.qs_trial_a/run_0003_ns_p0/` | 32000–42000／42060 | `SETTLED` |
| `case/90.qs_trial_a/run_0005_ns_p1/` | 20000–30000／30088 | `SETTLED` |

これは各基準に登録された量・区間についての判定です。3本とも全残差履歴と停止 step の `res_<step>.h5` の `VALUE/*` に非有限値はなく、ρ・P・T は正でした。

run_0014 を自身の成果物で評価すると **`accepted`** を再現しました。出口コア M **5.502841913**、η=0.03 の最大偏差 **0.062305 %**・片振幅 **0.057285 %**、壁解像も `PASS` です。対応は [case/71 の run 一覧](/home/sano/work/forge-cardmesh/case/71.m55_combustion_wt/README.md:23) で追跡できます。試走 A の最終 `accepted` は保存記録で確認しました。

手組み格子の保存記録は品質 `PASS`、AR最大 **4849.03**、AR緩和違反 **0**、第一層／要求厚さ最大 **0.99999965**。V8 の一次記録も存在し、V9 の保存データからは2400／4800の出口位置差最大 **1.9938×10⁻⁶ r_t** を再確認しました。生成器統合などの後続計画への移管は記録されており、今回へ戻す必要はありません。

回帰試験スイートは read-only 制約による一時ファイル作成失敗で実行できませんでした。上記は書き込み不要の判定とメモリ上の負例による確認です。`wt_qualify` の実行では描画モジュールだけを代替し、評価・照合処理は実コードを使用しました。

**推奨：成果物の生成元照合を修正し、上記負例の拒否と通常の `accepted` 維持を確認してから、§8 の範囲で accepted へ移してください。**

証拠: 足りる

指摘数: Critical 0 / Major 1 / Minor 0
