# codex レビュー: condensation-h2o-latent-supercritical (result)

- **plan**: [`plans/active/condensation-h2o-latent-supercritical.md`](../../plans/active/condensation-h2o-latent-supercritical.md)
- **stage**: `result` (diff base `HEAD~1`)
- **date**: 2026-09-14
- **commit**: `9daedec7` (feature/perf-3d-speedup)
- **codex**: effort `high`, 4.6 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M3/m4

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**

Watson 外挿・646.15 K での凍結・低温側の潜熱維持は確認できました。
ただし、高温の安全性テストに検証漏れがあり、場の回帰も計画に書かれた合格基準を満たしていません。以下を修正してから `accepted` に移してください。

確認対象は `580dd8be → 9daedec7` の指定 diff です。境界半割面・周期 seam・軸の幾何処理に変更はありません。既存ビルドの host テスト6本は再実行で PASS、低温潜熱の32,816点の A/B 出力も同一でした。GPU テスト2本の PASS は本環境では再確認できませんでした（`test_cond_float_device` は CUDA driver/runtime 不整合で失敗）。

1. **Major — 高温の有限性テストが、修正したゼロ除算経路を通っていません。**

   **根拠:** [test_cond_air.cpp:270](/home/sano/work/forge-perf/solver_density_cuda/tests/unit/test_cond_air.cpp:270) は蒸発条件で `cond_growth` を呼んでいます。しかし `growthModel=1` は [condensationSource_d.cuh:138](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/condensationSource_d.cuh:138) で潜熱評価前に `return 0.0`。H2O の `growthModel=0` は潜熱を使わない Hertz–Knudsen 式です。つまり、両モデルとも問題の `1/L²` を検証していません。

   今回追加したガードは別関数の [cond_evap_rate:249](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/condensationSource_d.cuh:249) にあります。「最大35.0 m/s、全点有限」という出力は再現しましたが、これをゼロ除算修正の証拠にはできません。

   **対案:** 指定した高温・半径の組み合わせで `cond_evap_rate` を直接呼び、負かつ有限であることを検査してください。float 経路も対象にし、成長側は過飽和条件で `1/L²` の枝を通してください。計画の存在しない関数名 `cond_growth_rate` も訂正が必要です。

2. **Major — 約束した高温・高液相率の二相温度反転試験が未実装です。**

   **根拠:** [plan:146](/home/sano/work/forge-perf/plans/active/condensation-h2o-latent-supercritical.md:146) は `g=0.1/0.5/0.9`、`T=380/500/640 K` の往復と失敗時の `roe` 保護確認を要求しています。追加された `(i)` 節は `Tsat` 試験で終了しており、この項目がありません。

   既存の [test_cond_float.cpp:197](/home/sano/work/forge-perf/solver_density_cuda/tests/unit/test_cond_float.cpp:197) も `g≤0.198`、高温側は500/900/1500 Kで、640 Kや凍結接続点を検査していません。一方、CPG 呼び出し側の `ok=false` 時に `roe` を更新しない分岐自体は [dependentVariables_d.cu:267](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/dependentVariables_d.cu:267) に存在します。

   **対案:** §6-7を実装し、646.15 Kの両側も追加してください。温度誤差・エネルギー残差・`ok` を検査し、意図的な反転失敗では実際の呼び出し経路で `roe` 保護を確認すること。§5.1の「実装＋単体完了」は、それまで未完了へ戻してください。

3. **Major — 回帰の PASS は、計画の基準とも「ノイズ床以下」という説明とも一致しません。**

   **根拠:** [plan:153](/home/sano/work/forge-perf/plans/active/condensation-h2o-latent-supercritical.md:153)・[plan:158](/home/sano/work/forge-perf/plans/active/condensation-h2o-latent-supercritical.md:158) は `VALUE/*` 全配列のビット同一を要求しています。実測は非同一です。

   以下の各 `res_300.h5` を読み、新側3反復の全ペア最大を配列ごとのノイズ床として再計算しました。

   - node: `case/16.nozzle_wys/run_0464_h2olw_node2d_{base,new,new_r2,new_r3}/`
   - cell: `case/16.nozzle_wys/run_0465_h2olw_cell2d_{base,new,new_r2,new_r3}/`

   | 系統・配列 | base–new 正規化最大差 | 新側反復のノイズ床 |
   |---|---:|---:|
   | node `Uy` | 9.8752e−6 | 1.4264e−5 |
   | node `P` | **1.7313e−6** | 1.1986e−6 |
   | node `roe` | **5.3651e−6** | 4.6070e−6 |
   | cell `Uy` | 1.0334e−3 | 1.3407e−3 |
   | cell `omega` | **1.7935e−4** | 1.4384e−4 |
   | cell `roe` | **4.2920e−5** | 3.8171e−5 |

   `Uy` の記載値は正しいですが、全配列が床以下ではありません。これだけで物性変更による劣化とは断定できない一方、「物性変更が場を動かさない」との断定もできません。

   **対案:** 非決定性を考慮する判定へ正式に変更し、[followups:65](/home/sano/work/forge-perf/plans/active/condensation-followups.md:65) の既存規約である全ペア最大床と `diff_res.py --tolfile --factor 2` に揃えて全配列を判定してください。ゼロ近傍量の絶対許容も明記し、§6・§8.1・READMEを同じ基準に統一してください。

4. **Minor — 検証成果物が不足しています。未収束の説明自体は正しいです。**

   **根拠:** 上記8 runには `CONVERGENCE_VERDICT.txt` と `residual_history.png` がありません。今回ツールを再実行した結果、全8 runとも次の判定でした。

   - `check_convergence.py`: **NOT CONVERGED (stalled/plateau)**
   - `check_quasisteady.py --quantity pmax`: **TRANSIENT-UNSETTLED**

   全 `res_300.h5` の `VALUE/*` に NaN/Inf はなく、記載された温度範囲と373.15 K超のセル数0も確認しました。[case README:53](/home/sano/work/forge-perf/case/16.nozzle_wys/README.md:53) が回帰専用と限定している点は適切です。

   **対案:** VERDICT、残差図、全配列比較結果を各 run に保存し、READMEから参照してください。短時間の回帰対照として扱う限り、未収束そのものを理由に長時間計算へ延ばす必要はありません。

5. **Minor — 凍結導入前の説明と、現在の誤差評価が混在しています。**

   **根拠:** [plan:93](/home/sano/work/forge-perf/plans/active/condensation-h2o-latent-supercritical.md:93) は647.095 K付近の急勾配・反転失敗を説明していますが、現行実装は646.15 K以上で定数です。この位置で `±0.1 K` の差分を取っても `dL/dT=0` になります。

   また、[methods/condensation.md:629](/home/sano/work/forge-perf/methods/condensation.md:629) の最大表誤差7.1e−6に対し、既存テストの再実行値は646.04 Kで **2.374e−5**。許容1e−4には収まりますが、同じ誤差として扱えません。`dL/dT` の7.427e−5も、[test_cond_float.cpp:60](/home/sano/work/forge-perf/solver_density_cuda/tests/unit/test_cond_float.cpp:60) で645.8–646.5 Kを除いた正規化指標です。

   **対案:** 凍結前の数値は旧案の再現例と明記し、現在の誤差には入力温度のfloat丸め・除外域・正規化式を付記してください。傾きの単位 `J/(kg·K²)` も `J/(kg·K)` に訂正してください。

6. **Minor — 「N2/空気の経路は無関係」という影響範囲は誤りです。**

   **根拠:** [plan:179](/home/sano/work/forge-perf/plans/active/condensation-h2o-latent-supercritical.md:179) は変更を `h2o_latent` のみに限定していますが、実際には種共通の [cond_Tsat:366](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/condensationProperties_d.cuh:366) と `cond_evap_rate` も変更しています。`cond_Tsat` の勾配変更にはH2O限定の分岐がありません。

   **対案:** §7に実際の変更ファイル・関数を列挙し、「N2の潜熱相関は不変、共通の飽和温度反転は変更」と記載してください。N2についても対応する飽和温度の往復結果を残してください。

7. **Minor — レビュー状態と後続課題の記録が先走っています。**

   **根拠:** [condensation-followups.md:57](/home/sano/work/forge-perf/plans/active/condensation-followups.md:57) は既に「plan/result レビュー済み・完了」としています。一方、対象planは result レビュー待ちです。[methods/condensation.md:632](/home/sano/work/forge-perf/methods/condensation.md:632) が後続表へ残すとしている高温湿潤状態の課題も、その表にありません。`plans/README.md:26` は依然 `draft` 表記です。

   **対案:** 本レビューの採否と修正作業を§5.1・§6.1へ記録し、高温湿潤状態の扱いを後続表へ明示的に移管してください。完了表記と索引は修正・再検証後に同期してください。`methods/index.md` の既存リンクは有効で、今回の変更による追加修正は不要です。

推奨は、**Watson＋凍結の実装方針を維持し、上記1→2→3を優先して検証を補正した後、4〜7の記録を整えて `accepted` へ移すこと**です。現状の「実装・検証完了」は承認できません。ファイル変更禁止に従い、本レビュー内容は plan 未反映です。

指摘数: Critical 0 / Major 3 / Minor 4
