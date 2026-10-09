# codex レビュー: time_integration-line-implicit-speed (plan)

- **plan**: [`plans/active/time_integration-line-implicit-speed.md`](../../plans/active/time_integration-line-implicit-speed.md)
- **stage**: `plan`
- **date**: 2026-10-09
- **commit**: `84fe2cee` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 3.6 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M5/m1
- **focus**: §6.4 (Thomas を float にする #8 の事前登録) と、その実装 (solver_density_cuda/cuda_forge/timeIntegration_d.cu の lineThomasFactorT_d・lineThomasSolveT_d・line_implicit_t・lineF32Mode・line_f32・launchLineFactor/launchLineSolve の f32 の分岐・比較の経路) と case/45.isobutane_m6_d155/jf32.sh・inv_cmp_judge.py・f32_time_judge.py だけを点検してほしい。§6.2 の点検 (notes/reviews/2026-10-09-time_integration-line-implicit-speed-plan-2.md) で指摘された型の穴が残っていないか、精度の組み合わせの表と実装が一致しているか、判定の基準・分岐に穴がないか

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**GO-with-changes**

§6.4 の精度表は、ISP 0・`flow_float=double` の登録条件では実装と整合しています。
ただし、F32cs の非有限検査と採否ゲートに誤合格の穴があり、現状では実験結果を採否に使えません。
対象の `run_0339〜0346` は手元に存在せず、以下はコード監査と書き込みを抑止した再現試験の結果です。

1. **Major — `jf32.sh` は最初の比較判定で停止する。**

   **根拠:** [`jf32.sh:18`](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/jf32.sh:18) のループは `m` と `r` だけを設定し、:22 で未設定の `$arm` を参照します。`set -u` が有効なので、通常の環境では `arm: unbound variable` で終了します。同じ展開を実行して終了コード **127** を確認しました。環境から `arm` を継承している場合も、両モードを同じ腕として判定してしまいます。

   **対案:** ループ内で `m=1 → arm=F32c`、`m=2 → arm=F32cs` を明示し、対応する `PASS` を設定してください。solver を呼ばない模擬実行で、両腕の判定と後段の選別まで確認すべきです。

2. **Major — F32cs の因子を検査しておらず、「非有限 0 件」を保証できない。**

   **根拠:** F32cs の factor は [`timeIntegration_d.cu:2754`](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:2754) で `line_f32::W/LU` に書きます。しかし比較は [:2820](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:2820) で、直前にゼロ化した別の **double 配列** `msh.line_LU_d/line_W_d` を読みます。「差は意味がない」というログ注記では、実因子の非有限検査まで失われる問題を解消できません。

   また、float 分解の [:2174](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:2174) は `+Inf` のピボットを通します。有限行列の末尾 2×2 ブロックを `[[3e38,3e38],[-3e38,3e38]]`、rhs をゼロとした host float32 再現では、**分解成功・LU に Inf 1 件・解は有限のゼロ**になりました。したがって、dq と η だけでは因子の非有限を必ず検出できません。

   **対案:** 実際に使った float の `LU/W/y` を、ライン上の有効要素に限定して検査してください。両腕の D/K/rhs・因子・補正の非有限を独立に集計し、1 件でも失格にします。分解では非有限も失敗扱いにしてください。

   なお、plan の「NaN も失敗、double 版と同じ」は誤記です。従来版の [:1845](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1845) は `pa < 1e-30` だけで、NaN を拒否しません。既定経路の変更有無も明記する必要があります。

3. **Major — 「判別不能」を記録しながら、終了コード 0 で合格にする。**

   **根拠:** [`inv_cmp_judge.py:32`](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/inv_cmp_judge.py:32) は LU が閾値を超えた solve を候補側の失格判定から除きます。しかし :36 の `ok` に「判別不能がない」という条件がありません。

   実スクリプトへ模擬ログを入力した結果は次のとおりです。

   | 入力 | 実際の判定 |
   |---|---|
   | F32cs η=1e-2、LU η=1e-8 | `pass: true`、終了コード 0 |
   | 最終 solve の LU に非有限 1 本 | `pass: true`、終了コード 0 |

   後者では、非有限累計の印字が η 評価より前にあることも効きます（[`timeIntegration_d.cu:2877`](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:2877)）。最終 η で増えた件数は、次の solve の累計に載りません。

   **対案:** `PASS / FAIL / INDETERMINATE` を分離し、後段へ進めるのは完全な `PASS` だけにしてください。非有限は腕を問わず無条件失格、有限な LU の η 超過は判別不能とし、ログの累計に依存せず各 η 行も集計します。

4. **Major — 全 solve・全ライン・全計測 step の完備性を確認していない。**

   **根拠:** [`inv_cmp_judge.py:18`](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/inv_cmp_judge.py:18) は各種類の行が一つでもあれば進みます。腕間の solve 番号の一致、期待する factor/solve 数、`n_eval` を合格条件に使っていません。

   再現では、**factor 1 回・solve 1 回だけでも合格**、さらに **評価ライン数 0 でも合格**しました。これでは「20 step の全ライン・全 solve」という登録条件を検証できません。

   [`f32_time_judge.py:7`](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/f32_time_judge.py:7) も件数しか確認しません。899 行あるものの、実際は **450 種類の step・449 step 欠落**という重複ログを与えても「速い」と判定しました。

   **対案:** 今回の条件では factor 20 回・solve 100 回、各 solve の両腕について 4719 本の評価を要求してください。番号の欠落・重複・腕間不一致は判定不能にします。性能ログは step 集合が厳密に `101…999` で、各 step が一度だけ現れることを検査してください。

5. **Major — 後処理の失敗後も進み、古い判定結果を利用できる。**

   **根拠:** [`jf32.sh:35`](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/jf32.sh:35) は性能判定の失敗をログに書くだけで、:37 で固定名の JSON を読みます。過去の JSON が残っていれば、今回の判定失敗後でも古い「速い」を使います。

   同様に [:45](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/jf32.sh:45) の時系列生成失敗も継続し、:51 は到達判定の成否を確認せず JSON をコピーします。そこで記録する `$?` は判定器ではなく **`cp` の終了コード**です。収束検査も戻り値を記録するだけで、許容した `NOT CONVERGED` と `DIVERGED` を分岐していません。

   **対案:** 判定器の実行失敗は即停止し、今回の run・バイナリに紐づく結果だけを読み込んでください。到達判定までの必須データが欠けた場合は候補選別を保留します。`NOT CONVERGED` の許容と、非有限・`DIVERGED` の失格を明確に分けてください。

6. **Minor — F32c/F32cs の時間差だけでは、演算とメモリの律速を確定できない。**

   **根拠:** [`plan:195`](/home/sano/work/forge-integ-1005/plans/active/time_integration-line-implicit-speed.md:195) は短縮比から主因を断定しています。しかし F32c は保存幅を保っても、演算型、変換命令、一時変数の幅が変わります。F32cs では保存幅に加えて変換も変わるため、純粋な帯域だけの比較ではありません。

   **対案:** この比較は「各変更の実効的な短縮量」の測定として扱い、律速の解釈は仮説に留めてください。#10 の実施判断は、配列変更そのものの A/B 計測で行うべきです。

**推奨は、既定 double を維持し、上記を修正した後に F32c/F32cs の限定実験を進めることです。** 修正の優先順は **①腕名の設定 → ②実因子の非有限検査 → ③判別不能の分離 → ④ログの完備性 → ⑤後処理失敗の遮断 → ⑥律速解釈の修正**です。

Thomas 本体では、`T=float` による積和・LU・代入、`S=double/float` による保存の使い分けを確認しました。緩和の乗算は `flow_float`、今回なら double です。この段も精度表に追記すると契約が明確になります。最終採用を #9 の品質・総時間評価に残す方針は妥当です。

ファイルは変更していません。修正案は **plan 未反映**です。

指摘数: Critical 0 / Major 5 / Minor 1
