# codex レビュー: time_integration-line-implicit-speed (plan)

- **plan**: [`plans/active/time_integration-line-implicit-speed.md`](../../plans/active/time_integration-line-implicit-speed.md)
- **stage**: `plan`
- **date**: 2026-10-10
- **commit**: `689c694f` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `medium`, 1.9 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M3/m0
- **focus**: §6.18 (S2・M64・M64+S2 の総時間、2 出力連続、損益分岐で打ち切り) と case/45.isobutane_m6_d155/tt.sh・tt_judge.py・m9_watch.py の --consec 追加だけを、投入前の不具合 (起動できない・誤合格・forge が監視なしで回る・予算の計算違い・判定の誤り・ディスク・既定の挙動の変化) に絞って短く点検してほしい。m9_watch.py の既定 (--consec 1) の挙動が従来と同じかも確かめてほしい
- **extra**: `case/45.isobutane_m6_d155/tt.sh`, `case/45.isobutane_m6_d155/tt_judge.py`, `case/45.isobutane_m6_d155/m9_watch.py`, `case/45.isobutane_m6_d155/scr_check.py`, `notes/reviews/2026-10-10-line-implicit-speed-results-diagnose.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**GO-with-changes**

§6.18 の予算式は記載どおりです。`--consec 1` は旧版との合成系列2,000件で一致し、B0 の2回連続到達も125000 stepを再現しました。ただし投入前に次の3点を修正してください。

1. **Major — ディスク監視が空きゼロまで発火しない。**  
   根拠: [tt.sh:76](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/tt.sh:76)。`df -BG` は切り上げなので、空き500 MBでも `1G` となり、`-lt 1` は偽です。また監視開始は最初の腕の起動から480秒後です。  
   **対案:** `-B1` でバイト比較し、最初の腕の起動時から監視する。容量取得失敗も停止条件にする。

2. **Major — 比較不可の結果が性能判定を通る。**  
   根拠: [tt_judge.py:70](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/tt_judge.py:70)。`CENSORED` を先に処理するため、収束判定が `DIVERGED`／「記録なし」でも「遅い」になります。条件分岐の実行で再現しました。また、収束ツールの「判定不能」は総括が `NOT CONVERGED` になる場合があり、現状では比較対象に入ります。  
   **対案:** 発散・記録欠落・入力の判定不能を先に除外し、その後で到達／予算打ち切りを分類する。

3. **Major — 終端状態では `--consec 2` の整合性検査を飛ばす。**  
   根拠: [m9_watch.py:75](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/m9_watch.py:75)。既存状態が `REACHED` なら、`phase/budget/consec` の検証前に正常終了します。旧来の1回到達を2回連続の結果として再利用でき、[tt_judge.py:69](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/tt_judge.py:69) も区別しません。  
   **対案:** 引数・保存設定の検証を終端判定より前に移し、最終判定器でも `phase=line`・`consec=2`・登録予算との一致を必須にする。

推奨は、**上記3点を修正してから、登録済みの3腕を投入すること**です。ファイルは変更していません。構文確認は通過しましたが、AWS側の `m9_unit.json` と未計測の単価がないため、実際の予算値・起動は未検証です。レビュー提案は **plan 未反映**です。

指摘数: Critical 0 / Major 3 / Minor 0
