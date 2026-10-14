# codex レビュー: time_integration-line-implicit-speed (plan)

- **plan**: [`plans/active/time_integration-line-implicit-speed.md`](../../plans/active/time_integration-line-implicit-speed.md)
- **stage**: `plan`
- **date**: 2026-10-10
- **commit**: `63b3ff39` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 3.1 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M4/m1
- **focus**: §6.11 (#9 の事前登録 v2) と case/45.isobutane_m6_d155/m9.sh・m9_watch.py・m9_unit.py・m9_judge.py・m9_ref.json だけを点検してほしい。前回の点検 (notes/reviews/2026-10-10-time_integration-line-implicit-speed-plan-2.md) の指摘が反映されているか、M2 の却下の根拠 (実測で出力は表示の step 1000 に入る) が妥当か、新たな穴がないか

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
参照 R の変更と前回 M2 の却下は妥当です。ただし、現スクリプトは正常な専有環境で起動できず、point 段の監視も例外で停止します。  
品質証拠の欠落を通す判定と再開処理も修正してから L5 を投入してください。

1. **Major — GPU が空いていると、単価計測を開始できない**

   根拠: [m9.sh:16](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/m9.sh:16)、[起動前検査:35](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/m9.sh:35)。

   `gpu_others()` の戻り値が末尾の `grep -c .` に依存しています。プロセスがゼロなら表示は `0` ですが、終了コードは **1** です。元の関数を使った模擬実行で、正常・空の `nvidia-smi` 出力に対して **`g=0, rc=1`** を再現しました。したがって、必要条件である専有状態で「nvidia-smi が失敗」と終了します。逆に他プロセスが存在すると起動し、9 本を回した後で `m9_unit.py` が拒否します。

   **対案:** コマンドの成功・失敗と件数を分離し、取得成功なら件数ゼロでも正常終了させる。`g != 0` は投入前に拒否する。空・他プロセスあり・取得失敗の3条件を確認すること。

2. **Major — point 段は最初の出力で `KeyError` になり、forge が回り続ける**

   根拠: [m9_watch.py:92](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/m9_watch.py:92)、[m9_ref.json:1](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/m9_ref.json:1)、[m9.sh:61](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/m9.sh:61)。

   E2 の内包表記が `REF` の全キーを走査しますが、参照 JSON には説明文字列の `source` があります。実ファイルの式と参照値をメモリ上で実行すると **`KeyError: 'source'`** になりました。これは読取例外の `try` の外なので、状態保存も `stop()` も実行されません。シェルは監視の非ゼロ終了を受けても `wait $job` に進むため、到達時停止・間引きが働かなくなります。

   また、モード不一致経路も SIGTERM 後に期限なしで待ち、登録した SIGKILL への移行を使っていません（[m9.sh:54](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/m9.sh:54)）。

   **対案:** E2 は数値4項目を明示し、型・有限性を起動前に検査する。監視例外・モード不一致にも共通の停止処理を適用し、終了確認と失敗状態の保存まで保証する。

3. **Major — 品質記録や E2 の再検証に失敗しても「速い」と判定できる**

   根拠: [m9_judge.py:46](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/m9_judge.py:46)、[L5 集計:74](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/m9_judge.py:74)、[比較条件:79](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/m9_judge.py:79)。

   - 収束記録の欠落・ツール異常は「記録なし」となり、比較を通ります。
   - 保存状態が `REACHED` なら、`reach(rows)` が未到達・非有限を返しても `total_s` を設定します。
   - 準定常ツールの実行失敗も比較条件に入りません。

   実際の `cmp()` に「収束記録なし・E2なし・準定常処理エラー」の模擬入力を与えると、**「L5 が速い」**を返しました。前回 M4 は未解消です。

   **対案:** 時間比較の前に、到達の再計算結果と `reach_step` の一致、必要な収束 VERDICT、準定常判定の正常実行を必須にする。証拠欠落・処理異常・非有限は比較不可とする。正常に得られた `NOT CONVERGED`／`NOT ALL STEADY` は、登録どおり「水準と E2 への到達時間」に限って扱い、収束解の同等性とは区別する。

4. **Major — 原子的保存だけでは再開できず、保存済みの到達を失う**

   根拠: [m9_watch.py:78](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/m9_watch.py:78)、[到達処理:99](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/m9_watch.py:99)、[m9.sh:28](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/m9.sh:28)。

   到達時に保存する最初の JSON は、到達行を含んでも状態が `running` のままです。停止処理中や最終状態の保存前に中断すると、再開時には保存済み行を処理済みとして飛ばし、到達を再判定しません。到達条件を満たす25,000 stepまで保存した状態を模擬すると、再開後は **`EXEC_ERROR, last_step=25000`** になりました。

   シェル全体の再実行も、最初に既存の `m9_unit.json` を削除し、既存の単価計測ディレクトリを検出して終了します。前回 M6 の「再開時の再判定」は実装されていません。

   **対案:** 起動前に保存系列から終了条件を再評価し、停止・削除・最終状態保存を繰り返しても同じ結果になる処理にする。新規実行と再開を分け、再開ではバイナリ・設定・入力場・参照値の一致を検査して有効な成果物を再利用する。

5. **Minor — `--tail 20000` の単位がツール仕様と違う**

   根拠: [m9_judge.py:43](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/m9_judge.py:43)、[check_quasisteady.py:542](/home/sano/work/forge-integ-1005/solver_density_cuda/tools/check_quasisteady.py:542)。

   `--tail` は step 数ではなく、標本数に対する割合です。今回は CSV を先に2万 step窓へ切り出しているため、巨大な割合を指定して偶然その全体を評価できています。

   **対案:** 窓の切出しは維持し、`--tail 1.0` に修正する。§6.11 の記載も揃える。

前回指摘のうち、次は確認できました。

- **M1 の参照変更は妥当。** `m9_ref.json` の4数値は保存系列の `run_0263_ns_coldmesh_tw300_cfl4_ext4`、40,000 stepと完全一致しました。同 run の20,000〜40,000 stepを正式ツールへ渡した結果は **`OVERALL: ALL STEADY`** でした。
- **M2 の却下は正しい。** [main.cpp:593](/home/sano/work/forge-integ-1005/solver_density_cuda/main.cpp:593) は表示に `iStep + 1` を使い、出力処理後に `report()` を呼びます。旧 commit `00da938b` でも同じです。前回レビューは内部添字と表示番号を取り違えています。ただし、記載の **34.65／274.75 ms** は手元に生ログがなく、直接再検算していません。
- **到達値を再現。** P は通算610,000 step、L0 は `run_0262_ns_coldmesh_tw300_cutback_point` の35,000 step。E2 は全4量で0.1%以内です。ただし正式ツールでは、P の590,000〜610,000 step窓は標本不足、L0 の15,000〜35,000 step窓は θ_r 3量が `TRANSIENT-UNSETTLED` で、**両方とも `OVERALL: NOT ALL STEADY`** でした。

**推奨は、v2 の評価規則を維持し、①起動・監視の確定的な不具合、②判定の証拠検査、③中断からの復旧、④`--tail` の順に直してから L5 を投入することです。** P／L0 の事後集計と推定総時間という位置付けは維持してよいです。

ファイルは変更していません。修正提案は **plan 未反映**です。

指摘数: Critical 0 / Major 4 / Minor 1
