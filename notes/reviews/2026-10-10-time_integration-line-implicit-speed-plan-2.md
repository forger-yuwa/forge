# codex レビュー: time_integration-line-implicit-speed (plan)

- **plan**: [`plans/active/time_integration-line-implicit-speed.md`](../../plans/active/time_integration-line-implicit-speed.md)
- **stage**: `plan`
- **date**: 2026-10-10
- **commit**: `df1fd6e3` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 4.4 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M7/m1
- **focus**: §6.9 (#9 仕上げまでの総時間の事前登録) と case/45.isobutane_m6_d155/m9.sh・m9_watch.py・m9_unit.py・m9_judge.py だけを点検してほしい。P と L0 は既存の系列の事後の集計であることを明記している。判定の規則・総時間の組み立て・E2・見張りの止め方と消し方・分岐に穴がないか、事後の集計の扱いが妥当か

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
P／L0 を事後集計として使うこと自体は妥当です。ただし、現規則では L0 が E2 不合格となり、主比較が成立しません。  
単価計測・品質判定・停止／再実行にも修正が必要です。ファイルは変更せず、保存 JSON の再計算とコードの読み取り・メモリ上の模擬実行で確認しました。

1. **Major — E2 の参照状態が不適切で、L0 は既に比較対象から落ちる**

   `m9_judge.py` の関数を保存系列に適用すると、以下を再現しました。

   | 項目 | 再計算結果 |
   |---|---:|
   | P の水準到達 | 570,000 step |
   | L0 のライン水準到達 | 120,000 step |
   | L0 の実際の切替 | 135,000 step |
   | L0 の point 段の水準到達 | 35,000 step |
   | L0 の通算到達 | 170,000 step |

   ところが、L0 到達時の P 到達時に対する E2 は、θ_r(40/70/94) が **+0.1821／+0.1590／+0.1479%**、Q_w が **−0.2180%**。全量が上限0.1%を外れます。現在の実装は正しく「比べない」を返すため、**L5 の結果によらず `L5_vs_L0` の速度判定ができません**。

   さらに P 自身も57万→64万 stepで θ_r(40) が+0.1110%、Q_w が−0.1290%変化します。E2 超過を「別の解」と断定する根拠はありません。参照側の過渡でも超えます。根拠: [m9_judge.py:28](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/m9_judge.py:28)、[E2 と比較分岐:49](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/m9_judge.py:49)。

   **対案:** 現 E2 の結果は不合格として残す。L5 投入前に、定常性を確認できた point 参照窓と、E1・E2・品質条件を同時に満たす終了規則を新しい版として登録する。閾値を L0 に合わせて緩めず、E2 超過の表現は「同等性未確認」とする。

2. **Major — 単価に最終出力の時間が混入する**

   `m9_unit.py` はログの101〜999を採用します。しかし forge のログは0始まりで、出力には `iStep+1` を渡した後に `monitor.report(iStep)` を呼びます。したがって **ログの999には `res_1000.h5` と境界出力の費用が入ります**。lineK の commit `00da938b` でも同じでした。

   根拠: [m9_unit.py:16](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/m9_unit.py:16)、[main.cpp:2535](/home/sano/work/forge-integ-1005/solver_density_cuda/main.cpp:2535)、[main.cpp:3458](/home/sano/work/forge-integ-1005/solver_density_cuda/main.cpp:3458)。

   この単価を数十万 stepへ掛けるため、出力1回の費用を過大に積算します。

   **対案:** 実際の出力に対応するログ行を除外する。今回なら999を外す。出力費用は別に測り、登録した出力回数に応じて加算する。

3. **Major — 「総時間」の組み立てと不確かさが登録内容に一致しない**

   起動時間の実装は「起動から1 step目」ではなく、**ラッパー全実行時間−ログの全 step時間**です。初期出力・終了処理なども含みます。また、腕別に起動時間を測っているのに、全腕へ P の起動単価だけを適用しています。根拠: [m9_unit.py:21](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/m9_unit.py:21)、[m9_judge.py:45](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/m9_judge.py:45)。

   起動時間の有限性・非負性も未検査です。メモリ上の模擬実行では、`wall.txt=nan` を与えても正常終了し、3腕すべての `median_startup_s` に NaN を出しました。

   さらに `uncert_s` は出力間隔だけです。単価3本のばらつき、開始場から終盤までの費用変化、旧バイナリの step 数を転用する不確かさは含みません。§6.0 は「再実行差の3倍以内」であり、§6.9 の前提にある長期軌道のビット一致を証明していません。

   **対案:** 指標を「共通環境の単価で換算した推定総時間」と定義し、計算・出力・腕別起動の内訳を固定する。全入力の有限性・非負性とログ全区間の完全性を検査し、単価のばらつきも比較区間へ反映する。旧系列の転用は明示的なモデル仮定として残す。

4. **Major — 品質ゲートと非有限チェックが最終判定につながっていない**

   `m9.sh` は `check_convergence` の終了コードを記録するだけで、DIVERGED／判定不能でも先へ進みます。`check_quasisteady` は呼ばれず、`m9_judge.py` も両 VERDICT を読みません。根拠: [m9.sh:43](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/m9.sh:43)、[m9.sh:54](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/m9.sh:54)。

   見張りの非有限判定は `cold_series.snapshot` に依存しますが、対象は `ro/Ux/Uy/T/k` だけで、ω・保存エネルギー・全残差列を覆いません。`reach()` は系列の `nonfinite` を無視します。模擬系列では、θ_r が `1→2→1→2→1` と振動し、途中に `nonfinite=10` があっても20,000 stepで到達扱いでした。根拠: [cold_series.py:25](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/cold_series.py:25)、[m9_judge.py:20](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/m9_judge.py:20)。

   P 到達直前2万 stepの保存系列も正式ツールへ渡しましたが、3標本しかなく、結果は **`TRANSIENT-UNSETTLED`、`OVERALL: NOT ALL STEADY`** でした。これは発散の意味ではなく、指定窓の判定資料が足りないという意味です。

   **対案:** 全腕について品質 VERDICT と対象区間を判定入力にする。DIVERGED・非有限・証拠欠落は比較不可。NOT CONVERGED を許すなら、その比較は「登録した状態水準への到達」に限定する。保存系列を CSV 入力として準定常判定する処理と閾値を明記する。

5. **Major — 停止・異常終了・打切りを区別できない**

   見張りは forge が見つからないだけで `ENDED` にします。正常に上限まで進んだか、起動失敗か、途中クラッシュかを確認しません。シェル側のバックグラウンド処理も、最後の `echo` で forge の失敗を覆い隠し、判定器の失敗後にも `m9.done` を作れます。失敗までの時間も判定 JSON に残りません。

   根拠: [m9_watch.py:81](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/m9_watch.py:81)、[m9.sh:37](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/m9.sh:37)、[m9.sh:55](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/m9.sh:55)。

   到達時も `SIGTERM` の送信後20秒待つだけで、終了を確認せず後続出力を削除します。その後の restart は「最大 stepのファイル」を選ぶため、削除後に出力が残れば登録した到達場と異なる場を使う余地があります。根拠: [m9_watch.py:65](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/m9_watch.py:65)、[cold_cfl.py:111](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/cold_cfl.py:111)。

   **対案:** 起動したプロセスを識別して終了を待ち、`REACHED / CENSORED / DIVERGED / EXEC_ERROR` を分離する。停止確認後に削除し、restart は `reach_step` とハッシュで指定する。監視失敗時にも子プロセスを回収する。

6. **Major — 読取失敗と再実行で、停止不能・古い証拠の再利用が起きる**

   読めない HDF5 は例外を握って再試行するため、forge 終了後の破損ファイルでも永久に待ちます。また、既存系列を読み込んだ step は処理済みとして飛ばす一方、状態は `running` に初期化します。系列と状態ファイルは別々の非原子的な保存なので、中断位置によって到達を再判定できません。

   根拠: [m9_watch.py:17](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/m9_watch.py:17)、[保存処理:44](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/m9_watch.py:44)、[例外処理:55](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/m9_watch.py:55)。

   通常の間引きでは系列の保存前に HDF5 を削除するため、その間の停止で一次データと未保存の集計行を両方失います。

   **対案:** 読取再試行に上限を設け、終了後の破損は `DATA_ERROR` にする。系列・状態を整合する形で原子的に確定してから削除する。再開時には既存系列を検証・再判定し、run の識別情報が合わなければ拒否する。

7. **Major — 実効設定とキャッシュの検証が不足する**

   `m9.sh` は `FORGE_LINE_NOOP`・`FORGE_LINE_DEBUG_POINT`・`FORGE_LINE_MONO`・`FORGE_FREEZE_TURB` を解除しません。これらは数値処理や費用を変える実在のスイッチです。単価判定は実効モードを確認せず、既存の `m9_unit.json` は存在するだけで再利用します。

   根拠: [m9.sh:7](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/m9.sh:7)、[キャッシュ分岐:16](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/m9.sh:16)、[m9_unit.py:15](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/m9_unit.py:15)。

   **対案:** 許可した環境変数だけで起動し、ログから実効モード・マスク5・設定・バイナリを検証する。キャッシュには設定・入力場・計測器のハッシュを結び付ける。`nvidia-smi` 失敗も「0プロセス」と区別する。

8. **Minor — 初回判定窓と step表記を揃える必要がある**

   見張りは step 0を除外するため、5000間隔では最初の2万 step窓が成立するのは25,000 stepです。既存 L0 の系列は step 0を含み、20,000 stepから評価できます。根拠: [m9_watch.py:53](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/m9_watch.py:53)。

   また §6.9 の「切り戻しの pointで17万 step」は、保存系列とコードでは「point 3.5万、通算17万」です。

   **対案:** step 0の扱いを統一し、各段の局所 stepと通算 stepを併記する。

**推奨は、上記を反映した §6.9 の改訂版を固定してから L5 を投入することです。** 優先順は、①E2参照と終了条件、②単価・総時間、③品質ゲート、④停止・削除・再開、⑤実効設定の検証です。

P／L0 の事後集計という表示、L0 に実際の13.5万 stepを計上する処理、上限の異なる「構成同士」を比較する整理は維持してよいです。ただし、P の57万 stepは保存系列で確認できる最初の到達であり、それ以前の欠けた系列まで含めた初到達の証明ではありません。手元には L5・単価計測の実測成果物がなく、速度の採否は未判定です。修正提案は read-only 指示に従い **plan 未反映**です。

指摘数: Critical 0 / Major 7 / Minor 1
