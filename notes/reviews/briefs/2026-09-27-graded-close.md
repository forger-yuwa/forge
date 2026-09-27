# 諮問ブリーフ: 流体側 plan (非一様格子の静止ガス) を閉じてよいか (2026-09-27)

AGENTS.md 条件 **7** (解釈の確定) と **1** (plan の結論・方針)。

## 読んでよいファイル
- plan `plans/active/axisymmetric-graded-grid-static-gas.md` (全文、§5.1 #2〜#4 に結果)
- `case/62.conjugate_disk/run_0020_hold_A_ext100k/`: `static_hold_series.csv`・`QUASISTEADY_hold.txt`・`CONVERGENCE_VERDICT_concat.txt`
- `case/62.conjugate_disk/README.md` の run 一覧
- 発注元 `plans/accepted/boundary-cht-axisymmetric-fem2d.md` §6 V-ax2b の再登録文 (非一様は純伝導 IC 限定)

## 結果 (登録どおり)
- A (一様 IC、非一様格子、非連成) を累積 100000 まで延長: 窓 90000–100000 で閾値 4 項目 PASS (max|U| 3.64e-4 m/s、熱流束誤差 ≤2.8e-3 %、相対圧力差 2.8e-7)、準定常 72 系列 ALL STEADY、残差 PASS (連結 13–14 桁)。**累積 40000 で既に静止状態**に達し以後不変。
- B (純伝導 IC) は 40000 で同じく全条件 PASS (max|U| 5.18e-4)。
- → 「異常状態に居座り続ける」を棄却、「この条件では反復延長で静止保持基準に到達可能」。
- 連成あり (run_0014、一様 IC、warmup 20000) の失敗は、step 20000 で連成を始めた時点で過渡がまだ大きかった (max|U| 0.23 m/s、壁熱流束の市松 ±4000 W/m²) ことと整合する。

## 私の提案 (未実施)
1. 流体側 plan を「**起動過渡の減衰が遅い (機構は未確認) が、有限の反復で静止保持に届く。対処は起動手順**」で閉じる (done → accepted)。
   起動手順の注意書き: 加熱壁のある静止ガスの CHT は、伝導平衡 IC を使うか、非連成の warmup で流体が静止保持基準に届いてから連成する (本ケースで 40000 step)。methods/boundary.md に反映。
2. 発注元 CHT plan の保証範囲 (非一様は純伝導 IC 限定) は**今は広げない**。広げるなら、一様 IC・warmup 60000 の CHT run を 1 本、V-ax2b と同じ基準で事前登録して回す (別作業)。
3. 減衰が非一様格子で遅い機構 (局所擬似時間刻みの分布など) は追わず、未確認として記録する。

## 諮りたいこと (推奨を 1 つに)
1. この結論で plan を閉じてよいか。書き方の修正があれば。
2. 閉じる前に要るもの (result レビュー 1 回で足りるか、plan レビューを飛ばしている点の扱い — 本 plan は §4/§6 を diagnose で決めたが plan 段の codex レビューは回していない)。
3. 提案 2 (保証範囲を広げる追加 run) をやる価値があるか。
