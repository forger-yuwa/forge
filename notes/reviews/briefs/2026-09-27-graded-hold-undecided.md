# 諮問ブリーフ: 静止保持 A/B が「判別未了」— 次の一手 (2026-09-27)

AGENTS.md 条件 **3** (事前登録した比較が判定不能) と **4**。

## 読んでよいファイル
- plan `plans/active/axisymmetric-graded-grid-static-gas.md` の §5.1 #3 の行 (登録と結果)
- `case/62.conjugate_disk/run_0017_hold_B_conduct/static_hold_series.csv` (100 step ごと、全 200 行)、`CONVERGENCE_VERDICT.txt`、`solverConfig.yaml`
- 評価器 `case/62.conjugate_disk/eval_static_hold.py`

## 事実
- B (純伝導 IC): 窓 10000–20000 で閾値 4 項目 PASS (max|U| 5.18e-4 m/s、熱流束誤差 hot 4.98e-3 % / cj 5.65e-3 %、相対圧力差 2.88e-7)、残差 PASS (6.2–7.6 桁)。
  準定常 (`--tail 0.5 --drift 0.001 --osc 0.001`): Umax・Uymax・dPrel・checker・qerr_hot は STEADY、**qerr_cj だけ DRIFTING (drift 0.6 %/tail、fluct 0.8 %)**。
- qerr_cj の末尾: 18900 5.65243e-3 → 20000 5.65259e-3 %、100 step あたりの増分 1.8e-8 → 1.0e-8 (単調増加・増分減衰)。値そのものは許容 0.5 % の 1/90。
- A (一様 IC) は異常を再現 (max|U| 1.48 m/s、全系列 DRIFTING、残差 NOT CONVERGED)。
- 登録の決定則: 「閾値内でも DRIFTING → 判別未了で止める」。

## 私の提案 (未実施)
B を最終状態 (step 20000) から同一設定で +20000 step 継続し (流体 `restart_field.py --keep-src-dtype`、連成しないので固体・壁状態は不要)、
**同じ閾値・同じ準定常基準**を、継続区間の末尾半分 (累積 30000–40000) に当てる。A は判定済み (異常再現) なので延長しない。

## 諮りたいこと (推奨を 1 つに)
1. この延長で再判定するのは妥当か (事後の緩和に当たらないか)。代案があれば 1 つ。
2. 値が許容の 1/90 の量に相対 drift 0.1 % を課すのは適切か — 今回は変えないが、今後の登録で「値の絶対 drift」を併用すべきか。
