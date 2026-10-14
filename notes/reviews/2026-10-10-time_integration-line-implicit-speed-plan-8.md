# codex レビュー: time_integration-line-implicit-speed (plan)

- **plan**: [`plans/active/time_integration-line-implicit-speed.md`](../../plans/active/time_integration-line-implicit-speed.md)
- **stage**: `plan`
- **date**: 2026-10-10
- **commit**: `a87d3bf5` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `medium`, 1.2 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M3/m1
- **focus**: §6.19 (LAYOUT2 を既定にする前の経路の確認と既定化の作り) と case/45.isobutane_m6_d155/lay3.sh だけを点検してほしい。(1) 確かめる経路の選び方に、既定化で壊れうる経路の抜けがないか (timeIntegration_d.cu の lineLayoutEnabled・line_layout::ensure・lineThomasFactorLP_d/lineThomasSolveLP_d・lineThomasFactor_d_wrapper/lineThomas_d_wrapper を grep で見て)。(2) 合格の規則 (η の判定不能の扱い) が誤合格を生まないか。(3) 既定化の作り (環境変数が無いときは 2、診断スイッチで従来へ、区画 > 1.5 × 被覆で従来へ) の穴。(4) lay3.sh の投入前の不具合。短く
- **extra**: `case/45.isobutane_m6_d155/lay3.sh`, `case/45.isobutane_m6_d155/line_cmp_judge.py`, `solver_density_cuda/mesh/mesh.cpp`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**

§6.19 の数値経路の選定は概ね妥当です。実装前に、台本の失敗伝播、実効モードの決定方法、メモリ不足時の退避を補ってください。ファイル変更・計算投入はしていません。

1. **Major — `lay3.sh` が失敗を終了状態に反映しない。**  
   根拠: `case/45.isobutane_m6_d155/lay3.sh:15` の `judge()` は最後の `echo` により成功を返します。判定器を終了コード 1 に置き換えた模擬実行でも `judge_return=0` でした。実行失敗もログだけで続行し、26 行で成果物を削除、41 行で `lay3.done` を作ります。`prep` は既存 run を拒否するため、再投入時は古いログを判定する危険もあります（`cold_cfl.py:106`）。  
   **対案:** 実行成功と今回生成された判定 JSON を必須にし、全件の合格を集約して終了コードに反映する。η の例外は `fail` が空で、`indeterminate` が指定理由だけの場合に限定する。失敗時の成果物は残す。

2. **Major — 新設する退避分岐の検証と、モードの一貫性の契約が不足。**  
   根拠: `plans/active/time_integration-line-implicit-speed.md:530` の確認には、診断スイッチ併用・区画比による退避がありません。現コードは有効判定が静的キャッシュ（`timeIntegration_d.cu:2884`）、配列確保が別処理（2901 行）、カーネル選択が別のモード参照（3103・3154 行）です。`ensure()` だけで退避すると、確保先と実行経路が食い違います。  
   **対案:** メッシュ情報を含めた実効モードを確保前に一度決め、factor・solve・比較・表示で共有する。診断スイッチ各種の「未指定なら退避／明示 1・2 なら拒否」、区画比の境界・超過、`lineImplicit=1` でライン 0 本を検証表に追加する。

3. **Major — 区画比 ≤1.5 でも、既定化によるメモリ不足は防げない。**  
   根拠: plan:528–529 は旧配列を保持します。`timeIntegration_d.cu:2908` 以降の追加確保はピボット込み **645 bytes/区画**で、等長ラインでも必要です。現在の `cudaMalloc` は失敗すると停止するため、従来動いたケースが既定変更だけで起動不能になり得ます。  
   **対案:** 自動選択時は追加確保に失敗したら部分確保を解放して従来へ退避する仕様を追加し、確保失敗を模擬して確認する。

4. **Minor — η 超過の説明が参照尺度の欠落に偏っている。**  
   根拠: `timeIntegration_d.cu:3036` は `flow_float` に保存された緩和後の補正から `dq/relax` を復元して η を測ります。float ビルドでは、この保存丸めでも 1e−11 を超え得ます。  
   **対案:** plan:523 にこの制約を追記する。厳密なビット一致と失敗ゼロを満たす場合の例外は、**配置変更の同値性確認としては妥当**ですが、線形解精度や収束の合格とは区別してください。

推奨は、**1 → 2 → 3 を計画と台本に反映してから既定化へ進む**ことです。上記はレビュー提案であり、plan 未反映です。

指摘数: Critical 0 / Major 3 / Minor 1
