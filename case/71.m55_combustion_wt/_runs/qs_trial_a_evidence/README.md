# 手順書の試走 A (case/90.qs_trial_a) の一次記録の写し (2026-10-16、M6ノズルの設計)

plan `tooling-settle-convergence-and-monitor` §8 の完了の範囲に入れた「試走 A」の証拠 (result 段レビュー 4 回目の証拠不足への対応、§5.1 #13)。
原本は AWS B の `~/qs_trial/forge/case/90.qs_trial_a/` (手順書 `procedures/design-chain-quickstart.md` の付録、HOME を `/home/ubuntu/qs_trial` にした素の clone)。
コードは p0 が a36a87ff、p1 の §2 から eabbab97、成果物のやり直しで 77aa8abf、公差 0.1 % の判定で 6c7bccc3 (B の clone で `git pull`)。

| ファイル | 中身 |
| --- | --- |
| `run_0003_ns_p0/`・`run_0005_ns_p1/` | 停止規則・判定・関門の記録 (`STOP_RULE_STATE.json`・`SETTLE_VERDICT.json`・`settle_criteria.yaml`・`CONVERGENCE_VERDICT.txt`・`HANDOFF_GATE.json`・`STOPPED_AT`・`MAIN_STAGE_STARTED`・`forge_launches.jsonl`・`boundary_massflux.csv` ほか)、較正・判定の出力 (`contur_calibrate/`・`wt_qualify*/`) |
| `listing_<run>.tsv` | run の全ファイルの一覧 (大きさ・時刻) |
| `nonfinite_check.txt` | 停止の step の場と残差の履歴の非有限の検査 (手元の写しで計算) |
| `card_p0.yaml`・`card_p1.yaml`・`README.md`・`_log/` | 試走のカード・case の README・各段のログ |

手元の大きい写し (git に入れない): `case/90.qs_trial_a/` (停止の step の場 `res_42060.h5`・`res_30088.h5`、`nozzle.h5`、残差の履歴、節点・境界の出力)。
case の番号 90 は試走の名前で、リポジトリの case としては使っていない。
