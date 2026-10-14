# codex レビュー: time_integration-line-implicit-speed (plan)

- **plan**: [`plans/active/time_integration-line-implicit-speed.md`](../../plans/active/time_integration-line-implicit-speed.md)
- **stage**: `plan`
- **date**: 2026-10-10
- **commit**: `dfa779fb` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `medium`, 1.2 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M3/m0
- **focus**: §6.13 (implicitRelax 0.85・1.0 の 4 腕の事前登録) と case/45.isobutane_m6_d155/m9r.sh・m9r_judge.py・cold_cfl.py の --relax・m9_watch.py の残す出力の変更だけを、投入前の不具合 (起動できない・誤合格・forge が監視なしで回る・ディスク) に絞って短く点検してほしい

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**

§6.13 の比較方針は妥当ですが、投入前に以下の3点を修正すべきです。
構文確認は合格。判定器の誤判定はメモリ上の模擬入力で再現しました。ファイル変更・forge 起動はしていません。

1. **Major — 判定記録が欠けても「速い」と判定する。**  
   根拠: [m9r_judge.py:60](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/m9r_judge.py:60)。新規腕の `DIVERGED` だけを拒否し、`記録なし` と比較相手の収束判定を検査していません。実際の `cmp()` に模擬入力を渡すと、**両段が「記録なし」でも「速い」**、**比較相手が `DIVERGED` でも「速い」**となりました。比較相手の準定常判定が欠けても通ります。  
   **対案:** 既存 `m9_judge.py` の `evidence()` と同じ前提検査を両腕に適用し、必要な VERDICT の欠落・判別不能・発散は比較不可にする。

2. **Major — 監視開始前の中断で forge が取り残される。**  
   根拠: [m9r.sh:28](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/m9r.sh:28)。forge をバックグラウンド起動し、180秒待ってから監視を開始します。`phase` を実行するシェルに終了時の停止処理がなく、この間にそのシェルへ `SIGTERM` が届くと、子の forge は監視なしで残り得ます。監視の非ゼロ終了への対処だけでは、この経路を覆えません。  
   **対案:** 各腕の起動前に `EXIT/TERM/INT/HUP` の停止処理を登録し、子プロセスの終了まで確認する。監視は起動直後から開始し、180秒後のモード確認も監視下で行う。

3. **Major — 「最新だけ残す」では空き7.6 GB以内を保証できない。**  
   根拠: [m9_watch.py:121](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/m9_watch.py:121) の削除対象は、解析済みの **step > 0** だけです。`res_0`、各段の `nozzle.h5`、未処理・書き込み中の出力は残ります。[cold_cfl.py:188](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/cold_cfl.py:188) は段ごとに入力を複製するため、4腕の仕上げ段追加でも容量が増えます。空き容量の事前検査・監視はありません。**満杯になるとは断定できませんが、登録した容量対策だけでは不足です。**  
   **対案:** 実ファイルサイズから全8段の固定保持量と同時出力のピークを見積もり、起動前の容量ゲートと実行中の下限停止を追加する。既に動いている監視プロセスにも新しい保持規則が適用されているか確認する。

**推奨:** 上記を **1 → 2 → 3** の順に修正・模擬確認してから、登録済み4腕を投入する。`--relax` の置換・差分検査には、指定した `0.85/1.0` を阻む不具合は見つかりませんでした。実 run がローカルに無いため、実容量と起動成功は未検証です。レビューのみのため **plan 未反映**です。

指摘数: Critical 0 / Major 3 / Minor 0
