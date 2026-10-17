# codex レビュー: tooling-nozzle-card-local-mesh (result)

- **plan**: [`plans/accepted/tooling-nozzle-card-local-mesh.md`](../../plans/accepted/tooling-nozzle-card-local-mesh.md)
- **stage**: `result` (diff base `48cab832`)
- **date**: 2026-10-17
- **commit**: `126e00c8` (feature/card-local-mesh)
- **分類**: `ok`
- **backend**: `codex` — 選んだ理由: 既定 (auto: codex を試し、unavailable のときだけ claude)
- **要求したモデル**: `既定 (~/.codex/config.toml)`
- **実効のモデル**: `gpt-6-astra` (codex のログの見出し)
- **CLI の版**: `codex-cli 0.162.1`
- **依頼側との関係**: 不明 (`--requester-model` 未指定)
- **permission_denials**: 対象外 (codex は read-only サンドボックス)
- **実行**: effort `high`, 6.5 min, rc=0
- **試行のログ**: `notes/reviews/2026-10-17-tooling-nozzle-card-local-mesh-result-2.codex.log` (codex: ok)
- **判定**: **NO-GO**, 指摘 C0/M2/m1
- **focus**: 前回の result 段 (notes/reviews/2026-10-17-tooling-nozzle-card-local-mesh-result.md、NO-GO、M1・M2・m1・m2・V8 の証拠不足) への対応 (§5.1 #7a・#7b、V7 の訂正、methods/design/overview.md、_runs/v8_float_fp64/) を確かめ、§8 の完了の範囲で §6 の判定が事前登録どおりかを見る

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

判定: NO-GO

前回の負例への修正、V7 の訂正、V8 の証拠追加は確認できました。ただし、未監視の量や未評価の壁を合格扱いできる経路が残っています。`accepted` への移動は保留を推奨します。

1. **Major — 停止判定と現在の設定の照合がなく、未監視の線を保証済みにできる**

   [wt_qualify.py:225](/home/sano/work/forge-cardmesh/design/forge_design/feedback/wt_qualify.py:225) は `STOP_RULE_STATE` の `REACHED` と `SETTLED` だけで計算成立を認め、[同:407](/home/sano/work/forge-cardmesh/design/forge_design/feedback/wt_qualify.py:407) は現在の `settle_criteria.yaml` と評価線が一致すれば準定常確認を省略します。停止時の基準・設定との同一性を検証していません。

   `case/71.m55_combustion_wt/run_0014_chain_v2/` で再現しました。実際の監視線 η=0.03 に対して η=0.06 を評価すると、全場系列不足で正しく `incomplete` になります。しかし、**基準ファイルの読み取り値だけをメモリ上で η=0.06 に変更**すると、古い停止記録のまま `accepted` になります。ファイルは変更していません。

   **対案:** [deltastar_loop.py:1152](/home/sano/work/forge-cardmesh/design/forge_design/feedback/deltastar_loop.py:1152) の `_v2_binding` など既存の照合処理を再利用し、停止時の判定・基準・入力が結び付いた場合だけ監視済みと認める。不一致・記録欠落は `incomplete` とし、基準変更後の負例を追加してください。

2. **Major — `--groups` 省略時、壁ダンプが一度もない壁を見落とす**

   [check_wall_resolution.py:330](/home/sano/work/forge-cardmesh/solver_density_cuda/tools/check_wall_resolution.py:330) の壁名・`physID` の照合は `--groups` 指定時だけです。省略時は存在するダンプから対象を列挙するため、未出力の壁が検査対象に入りません。[deltastar_loop.py:683](/home/sano/work/forge-cardmesh/design/forge_design/feedback/deltastar_loop.py:683) の受け渡し判定は、この省略経路を使用しています。

   同じ `run_0014_chain_v2/` の境界設定の読み取り値に、ダンプのない `missing_wall` をメモリ上で追加して確認しました。

   - `--step 46171`：`wall/3@46171` だけを評価し、終了コード 0・`PASS`。
   - `--step 46171 --groups wall,missing_wall`：終了コード 2・`INDETERMINATE`。

   **対案:** 省略時も `bcondConfig.yaml` から評価必須の壁集合を作り、全対象の指定 step のダンプを要求してください。「存在する壁の指定 step が欠ける」負例に加え、「対象壁のダンプが一度もない」負例が必要です。

3. **Minor — 現在仕様に、検証で保証していない無波動の断定が残る**

   [methods/design/overview.md:828](/home/sano/work/forge-cardmesh/methods/design/overview.md:828) の「なめらかさを保つので軸の波を作らない」は、壁形状の滑らかさから流れの無波動を保証する表現です。V7 を未収束による判定保留へ訂正した内容とも釣り合いません。

   **対案:** 「壁補正の滑らかさを保つ」と記述し、軸付近の Mach 偏差は対象 run・評価線・判定区間に限定した実測として示してください。

確認できた証拠は次のとおりです。

- `check_settled.py --no-write` の再実行結果：

  | run（リポジトリルート相対） | 判定区間 | VERDICT |
  |---|---:|---|
  | `case/71.m55_combustion_wt/run_0014_chain_v2/` | 36000–46000、停止 step 46171 の確認を含む | `SETTLED` |
  | `case/90.qs_trial_a/run_0003_ns_p0/` | 32000–42000、停止 step 42060 の確認を含む | `SETTLED` |
  | `case/90.qs_trial_a/run_0005_ns_p1/` | 20000–30000、停止 step 30088 の確認を含む | `SETTLED` |

  これは各基準に設定された量・区間についての判定です。各 run の全残差履歴と最終 `res_*.h5` の `VALUE/*` に非有限値はなく、ρ・P・T は正でした。

- `run_0014_chain_v2/` は実際の設定で `accepted`、停止 step の壁解像も `PASS`。前回の「存在しない step・壁名」の負例も拒否されます。上記 Major は、これとは別の残存経路です。
- V7 の `run_0010_ns_b0` の保存判定は `NOT CONVERGED`。plan の判定保留への訂正は妥当です。
- V8 の追加記録は [v8_float_fp64](/home/sano/work/forge-cardmesh/case/71.m55_combustion_wt/_runs/v8_float_fp64) で確認しました。抽出 CSV から再計算した δ_E 最大差は、FP64 同士 `8.78×10⁻⁶ mm`、float 同士 `0.05359 mm`。前回不足していた比較証拠は補われています。これは指定 step の比較であり、定常解の一致認定ではありません。
- V9 の保存データから、網 2400／4800 の x_F 最大差 `1.994×10⁻⁶ r_t` を再確認しました。手動格子の保存品質判定も `PASS（AR≤5000、skew≤0.90）` でした。
- run の対応は [case/71 の run 一覧](/home/sano/work/forge-cardmesh/case/71.m55_combustion_wt/README.md:23) で追跡できます。試走 A の最終 `accepted` は保存記録で確認しましたが、手元の case/90 には成果物一式がなく、総合判定の再現は `incomplete` です。

回帰試験スイートは、一時ディレクトリを作れない read-only 制約で起動できませんでした。上記は保存記録の照合、書き込み不要の判定実行、メモリ上の入力差し替えによる確認です。

**推奨は一つです:** §5.1 の #7a・#7b を再オープンし、上記二つの誤合格経路と文書の断定を修正してから `accepted` へ移してください。生成器への統合など、明示的に別 plan へ移した作業を今回へ戻す必要はありません。

証拠: 足りる

指摘数: Critical 0 / Major 2 / Minor 1
