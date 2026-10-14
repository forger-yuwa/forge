# 諮問ブリーフ: V-h370 の判定の基準の具体化 (Euler 参照の資格・抽出の参照の A/B の式) の確認 (2026-10-14、投入前)

- 依頼者: 主セッション `M6ノズルの設計` (エスカレーション条件 1・3: §6 の事前登録の基準を、結果を見る前に確定する)
- plan: `plans/active/tooling-nozzle-isothermal-wall-chain.md` §6「V-h370 の事前登録」、§5.1 #36e・#36f
- 前の諮問: `notes/reviews/2026-10-14-h370-registration-diagnose.md` (登録の形を決めた)
- 実装: commit `b5b53a1b` (`case/45.isobutane_m6_d155/h370_run.py`・`h370_judge.py`・`test_h370_launch.py`、手元の試験 83 項目 ALL PASS)。forge はまだ回していない。

## 観測事実 (実装者が登録の文言を具体化した点)

### F1 Euler 参照 (`run_0950_h370_euler_ref`) の資格
- 登録の文言: 「全残差の `check_convergence.py` PASS、NaN/Inf なし、出口 M・流量・参照場の時系列の `check_quasisteady.py` の VERDICT。資格を満たした時点の res と sha を固定。参照の時間の変動が抽出の差の判定を左右するなら資格未達」。
- 具体化 (`h370_judge.py --euler-qualify`):
  - RUN_RC 0、投入前の検査の合格、本段の最後の出力が res_60000 (Euler の経路 run_staged には早く止める仕組みが無いので常に 60,000 step 回る)。
  - 全段の残差・全出力に NaN/Inf なし、ρ・P・T > 0。`check_convergence.py --segment --drop 3` PASS。
  - `check_quasisteady.py` (CSV 系列モード、**既定の閾値** tail 0.4・drift 5 %・osc 10 %) で、出口コア M・流量 (stop_rule.recorded_quantities と同じ定義) と、**「参照場」= 軸と η 0.1 の線の M (x 40・70・94 r_t)** がすべて STEADY。OSCILLATING も未達。
  - `nozzle_report.metrics` の η 0.1 の評価量 (wave・overshoot・slope・range) は平均が 0 に近く相対の閾値が意味を持たないので記録だけ。
  - 固定するのは本段の最後の出力 (res_60000) とその sha。σ_ref を測るため末尾の窓の始まりと中ほど (37000・49000) の出力の sha も記録。
- 懸念: 既定の閾値 (相対 5 %・10 %) は 0.01 mm の寸法の予算に比べて緩い。参照の時間の変動の効き方は F2 の σ_ref で別に測る設計。

### F2 抽出の参照の A/B の式
- 登録の文言: 「結果 A: 共通の使用域の全断面・判定窓で max|δ_E,B − δ_E,A| ≤ 0.01 mm かつ抽出の合否が同じ。結果 B: 0.01 mm 超か合否が変わる。時間の変動や参照の未整定で判別できなければ判定不能」。
- 具体化 (`h370_judge.py`):
  - D_k = 2 窓の各出力 k での、両腕が抽出した共通の列の max|δ_E,B − δ_E,A| (外挿しない)。腕 B の値は stop_rule.json の見張りの記録 (sha が合わなければ抽出し直す)。
  - σ_ref = 停止の出力を、Euler の 37000・49000 の出力で抽出し直したときの max|Δδ_E| (参照の時間の変動の効き)。
  - **結果 A**: max_k D_k + σ_ref ≤ 0.01 mm かつ全出力で抽出の合否が同じ。
  - **結果 B**: どれかの出力で抽出の合否が違う、または min_k D_k − σ_ref > 0.01 mm。
  - それ以外 (窓の中で 0.01 mm をまたぐ、σ_ref でまたぐ、σ_ref を測れない、記録の欠け) は判定不能。
  - 抽出の合否 = 「例外なし・共通の列に held なし・値が有限」。
  - 腕 A の次の壁は、共通の列だけ腕 A の値に替えた表で作る (記録だけ)。

## 期待値と出典
- 0.01 mm は検証済みの精度でなく暫定の予算 (登録)。

## 問い
1. F1 の具体化でよいか。特に (a) 「参照場」を軸と η 0.1 の線の M (3 断面) とすること、(b) `check_quasisteady.py` の既定の閾値のままでよいか (σ_ref で効き方を別に測るので資格は緩くてよい、と読んでよいか)、(c) 固定を res_60000 にすること。
2. F2 の式でよいか (σ_ref の定義と、A・B・判定不能の境界)。
3. 結果を見る前に直すべき点があれば、最小で。

## 読んでよいファイル
- `case/45.isobutane_m6_d155/h370_judge.py`・`h370_run.py`・`test_h370_launch.py`
- `plans/active/tooling-nozzle-isothermal-wall-chain.md` §6 V-h370・§4.9 の codex diagnose 3
- `solver_density_cuda/tools/check_quasisteady.py`、`design/forge_design/evaluate/stop_rule.py` (`recorded_quantities`)、`design/forge_design/report/nozzle_report.py`
