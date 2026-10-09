# codex レビュー: time_integration-line-implicit-speed (plan)

- **plan**: [`plans/active/time_integration-line-implicit-speed.md`](../../plans/active/time_integration-line-implicit-speed.md)
- **stage**: `plan`
- **date**: 2026-10-10
- **commit**: `fabf2d4c` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 2.3 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M3/m1
- **focus**: §6.11 の末尾の「追加の腕 L5L」と case/45.isobutane_m6_d155/m9_branch.sh・m9_watch.py (--phase line_e2)・m9_judge.py (L5L の集計と比較) だけを点検してほしい。L5 のライン段は走行中 (まだ水準に未到達) で、到達の前に入れた追加。穴があれば到達の前に直せる範囲で

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
L5L を追加する意義はあります。ただし、見張りの起動不能と、L5L の到達時間を過大評価する判定があります。  
対象 run 本体はローカルに無く、以下はコードとメモリ内の再現試験による指摘です。実測の収束・優劣は判定していません。

1. **Major — `--budget` がコメントに取り込まれ、見張りが起動しない。**

   根拠: [m9_watch.py:19](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/m9_watch.py:19)。`line_e2` の説明コメント以降に `ap.add_argument("--budget", ...)` があり、引数登録されません。実際の引数定義を抽出して実行すると、`line`・`point`・`line_e2` の全てが `unrecognized arguments: --budget 200000`、終了コード `2` になりました。

   [m9_branch.sh:35](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/m9_branch.sh:35) はその失敗後も `wait $job` に進むため、forge が無監視で走り続けます。走行中の旧見張りとは別に、**次に起動される L5 の point 段も同じ問題を受けます**。

   **対案:** `--budget` 登録を独立した行へ戻し、3 phase の引数解析を確認する。併せて見張りの非ゼロ終了時は forge を停止して `EXEC_ERROR` を残す。参照読込など、現在 `try` の外にある初期化の失敗も同じ停止処理に含める。

2. **Major — 同じライン構成の継続なのに履歴を捨て、L5L に最低 25000 step を余計に要求する。**

   根拠: [m9_watch.py:61](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/m9_watch.py:61) は現在の run 内で 20000 step の幅を要求し、[m9_judge.py:93](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/m9_judge.py:93) も `run_0363` だけで再判定します。step 0 を除外するため、5000 間隔では最早到達が **25000 step** です。

   合成入力で、分岐元から全量が参照値と等しく、欠損 `0.01`・ドリフト `0` の系列でも、継続後 20000 step までは未到達、25000 step で初めて到達となりました。分岐時点で既に E を満たしていても、追加計算を必須にします。これは「ラインを続けた場合の最初の E 到達」を測れていません。

   **対案:** `line_e2` は、同一の実効設定を確認したうえで `run_0353` の末尾履歴を絶対 step で引き継ぐ。分岐点でも E2 を評価し、既に満たせば追加 step を `0` とする。見張り・再集計・準定常窓・起動回数をこの定義に統一する。point への切替は設定が変わるため、従来どおり切替後の窓で判定する。

3. **Major — point 側の状態ファイルが無いと、独立した L5L の比較まで落ちる。**

   根拠: [m9_judge.py:75](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/m9_judge.py:75)。共通ライン段が `REACHED` なら、L5L の処理より先に `run_0354/m9_watch.json` を無条件で読みます。欠落を模擬すると `FileNotFoundError` となり、L5L の集計に到達しません。point 側の起動・見張り失敗によって、正常終了した L5L と L0・P の比較まで失われます。

   **対案:** 各腕を独立して読み、欠落・破損はその腕の「比較不可」として記録する。共通ライン段の証拠が有効なら、point 側の失敗にかかわらず L5L の集計と L0・P との比較を完了させる。

4. **Minor — `NOT CONVERGED` の場合に必要な限定文が抜ける。**

   根拠: [m9_judge.py:115](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/m9_judge.py:115)。限定文は準定常判定だけを見ています。全段 `NOT CONVERGED`、準定常 `ALL STEADY` の合成入力では、限定なしの「L5L が速い」が返りました。§6.11 の品質条件と不整合です。

   **対案:** 比較対象の段に `NOT CONVERGED` が一つでもあれば、「水準と E2 への到達時間の比較に限る」を必ず付ける。

推奨は、**L5L の追加を維持し、上記 1 → 2 → 3 → 4 を分岐前に修正してから実行すること**です。特に 2 の到達定義は結果を見る前に §6.11 へ反映してください。ファイルは変更しておらず、本レビューの対案は plan 未反映です。

指摘数: Critical 0 / Major 3 / Minor 1
