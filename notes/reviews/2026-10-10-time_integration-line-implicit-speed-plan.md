# codex レビュー: time_integration-line-implicit-speed (plan)

- **plan**: [`plans/active/time_integration-line-implicit-speed.md`](../../plans/active/time_integration-line-implicit-speed.md)
- **stage**: `plan`
- **date**: 2026-10-10
- **commit**: `ce548732` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 4.0 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M2/m4
- **focus**: §6.6 (ncu の結果の読み) と §6.7 (Thomas の配列の並べ替え #10 の事前登録)、その実装 (solver_density_cuda/cuda_forge/timeIntegration_d.cu の lineThomasFactorL_d・lineThomasSolveL_d・lineLIdx・lineLayoutEnabled・line_layout・launchLineFactor/launchLineSolve の layout の分岐・比較の経路) と case/45.isobutane_m6_d155/jlay.sh・line_cmp_judge.py (--bitwise)・time_pairs_judge.py・ncu_wrap_lay.sh だけを点検してほしい。演算が既定の lineThomasFactor_d<false>・lineThomasSolve_d<false> と要素ごとに同じか (ビット一致の前提が成り立つか)、判定の穴、ncu の読みの過大な主張

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**GO-with-changes**  
配列の並べ替えを次に試す方針は妥当です。既定版と演算順序は対応していますが、判定器に偽 PASS の穴があり、`ncu` 表にも異なる sweep の混在があります。  
以下を試験投入前に修正してください。ファイルは変更していません。

確認した範囲では、`lineLIdx` の添字変換、D の読込み、Kprev·W の積和、LU・ピボット、前進・後退代入、緩和、分解失敗時の処理は既定版に対応しています。Kprev のコピーも、既存の `storeLU` による凍結と整合します。根拠は [timeIntegration_d.cu:1993](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1993)、[同:2169](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:2169)、[同:2924](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:2924) です。**ソース上、数値を変える変更は見つかりませんでした。ただし GPU でのビット一致は未確認**です。ローカルには `run_0348` 以降の実測結果がありません。

1. **Major — `line_cmp_judge.py` は必要な証拠が欠けても PASS になる。**

   根拠：[line_cmp_judge.py:20](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/line_cmp_judge.py:20) は factor の番号だけを数え、失敗比較の完全な行や `[lineNonfinite]` の存在・一意性を要求していません。同ファイル:36 で取得する `n_over` も判定に使いません。

   書込みをモック化した実行で、次の偽 PASS を再現しました。

   - `[lineNonfinite]` を全削除 → `PASS`
   - factor 行を番号だけに置換し、失敗比較を削除 → `PASS`
   - η 最大 `1.000e-11`、超過件数 **1 本** → `PASS`
   - η 最大 `nan`、非有限件数 0 → `PASS`

   3 番目は実際にも起こり得ます。CUDA 側は最大値を `%.3e` で出すため、例えば **1.0004e−11 が閾値と同じ表示**になります。[timeIntegration_d.cu:2915](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:2915)

   **対案：** factor・solve の完全な比較行と非有限検査行を番号ごとにちょうど 1 本要求し、欠落・重複は `INDETERMINATE`。数値は明示的に有限性を検査する。今回の閾値では `n_over` を採否に使い、最大値も十分な桁数で出力してください。

2. **Major — §6.6 の block カーネル行は、異なる sweep の測定値を混ぜている。**

   根拠：[plan:228](/home/sano/work/forge-integ-1005/plans/active/time_integration-line-implicit-speed.md:228) と、`case/45.isobutane_m6_d155/_band_ab/cold_pair/ncu/run_0347_{details,raw}.csv` の ID 0・3。

   | 測定 ID | 時間 | DRAM 使用率 | 読み／書き | load sectors/request |
   |---|---:|---:|---:|---:|
   | 0：最初の block | 2.47 ms | 89.54% | 0.614／0.714 GB | 10.04 |
   | 3：後続の block | 0.776 ms | 88.57% | 0.232／0.180 GB | 10.86 |

   plan は ID 3 の時間・使用率・占有率に、ID 0 の通信量・セクタ比・scoreboard を接合しています。したがって、この行は実在する 1 回の実行を表しません。Thomas の factor・solve の主要数値は CSV と対応していました。

   **対案：** block を初回と後続に分け、各行に測定 ID と sweep を記載する。LAYOUT 側も同じ位相どうしで比較してください。

3. **Minor — `--bitwise` はビット比較ではなく、浮動小数点の数値比較。**

   根拠：[timeIntegration_d.cu:2821](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:2821) の不一致判定は `x[i] != ref[i]` です。**+0 と −0 は異なるビット列でも通ります**。また LAYOUT の因子比較は、実際の並べ替えバッファを比較していません。[同:3039](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:3039)

   **対案：** `dq` を元の型幅の整数ビット列として比較する。数値不変の切り分けを確実にするため、LU・W・ピボットも `lineLIdx` で対応付けて比較してください。現在の検査結果を「ビット一致」と呼ぶのは不正確です。

4. **Minor — 非連続アクセスの証拠は強いが、律速原因の確定と実験失敗の解釈は強すぎる。**

   根拠：[plan:230](/home/sano/work/forge-integ-1005/plans/active/time_integration-line-implicit-speed.md:230)、[同:248](/home/sano/work/forge-integ-1005/plans/active/time_integration-line-implicit-speed.md:248)。

   31.28／31.24 sectors/request はアクセス効率の悪さを支持します。ただし scoreboard の 19.1／47.1 は時間百分率ではありません。CSV の発行間隔で割ると long scoreboard は factor 約 **48%**、solve 約 **77%**です。キャッシュや依存連鎖を含む待ちを、変更対象の配列だけに帰属できません。[NVIDIA の定義](https://docs.nvidia.com/nsight-compute/ProfilingGuide/index.html#l1-tex-cache)でも、セクタ比は L1 アクセス効率の指標です。

   **対案：** 「非連続アクセスを有力原因として #10 を検証する」に修正する。速くならなければ「この実装では改善を得られなかった」と判定してください。D・Knext などは未変更で、コピー費用も増えるため、「律速の読みが外れた」とは限りません。

5. **Minor — 性能判定は、測った腕と専有条件を検証していない。**

   根拠：[time_pairs_judge.py:9](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/time_pairs_judge.py:9) は時間行だけを読み、[jlay.sh:29](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/jlay.sh:29) は GPU プロセス数を記録するだけです。

   模擬ログでは、**両腕の実効モードが LU でも**時間差があれば「速い」、B の全時間が **0 ms**でも「速い」になりました。

   **対案：** 各 run の factor・solve の実効モードが期待する LU／LAYOUT であること、比較・診断が無効であること、時間が有限かつ正であることを確認する。GPU 専有条件を満たさない組は採否から除外してください。

6. **Minor — `ncu` が失敗しても正常終了扱いになる。**

   根拠：[jlay.sh:38](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/jlay.sh:38) は終了コードを記録するだけで、その後に出力削除と `jlay.done` 作成へ進みます。`ncu_wrap_lay.sh` のカーネル選択は LAYOUT の名前に対応していますが、取得成功を保証する検査がありません。

   **対案：** `ncu` の非ゼロ終了、レポート欠落、対象カーネルの内訳不一致を失敗として伝播する。factor 1 回・solve 5 回・block 5 回と実効 LAYOUT を確認してから、完了記録と削除へ進めてください。

**推奨は、#10 を opt-in のまま継続することです。** 優先順は、①比較器の偽 PASS とビット比較、②`ncu` 表・因果解釈、③性能・プロファイル実行の検査、です。その修正後に登録済み A/B を実施し、既定化前には §5.1 #5 の可変長・部分被覆・凍結・周期・拘束行の確認を完了させてください。長期到達試験の省略は、数値不変を確認できた適用範囲に限定すべきです。レビュー提案は **plan 未反映**です。

指摘数: Critical 0 / Major 2 / Minor 4
