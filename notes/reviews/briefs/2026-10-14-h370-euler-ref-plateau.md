# 諮問ブリーフ: V-h370 の Euler 参照 run_0950 が予備の資格の残差の条件で FAIL (残差は 2 桁で有界な横ばい、量は STEADY) (2026-10-14)

- 依頼者: 主セッション `M6ノズルの設計` (エスカレーション条件 3: §6 に事前に書いた条件の FAIL、条件 4: 登録に無い対処へ進む前)
- plan: `plans/active/tooling-nozzle-isothermal-wall-chain.md` §6「V-h370 の事前登録」・§5.1 #36f・#36i
- 前の諮問: `notes/reviews/2026-10-14-h370-registration-diagnose.md` (Euler の資格: 全残差の `check_convergence.py` PASS ほか)、`2026-10-14-h370-judge-criteria-diagnose.md`、`2026-10-14-h370-prelaunch-review-diagnose.md` (σ_ref は 17 × 全 20 出力、使用域の σ_ref,use ≤ 0.01 mm を独立の必要条件、判定不能は更新不可; #36i で実装中)

## 観測事実 (AWS B、`case/45.isobutane_m6_d155/run_0950_h370_euler_ref/`、台本 fdf3b728、FP64 7e95e5bb)

- 段: soft (1 次、3,000 step、19 s) → mid (2 次、cfl 1、3,000) → main (2 次、cfl 2、60,000、339 s)。RUN_RC 0、全段の残差・60 出力に NaN/Inf なし。投入前の検査は合格 (品質 PASS、変換後の δ_in 最大 7.96e−4 mm・全列被覆)。
- `check_convergence.py --segment --drop 3` (判定区間 = mid → main、stage_manifest の hard キーが同じ区間): **NOT CONVERGED (stalled/plateau)**。全列の低下は 2.0〜2.2 桁 (rms_ro 2.47e−4 → 2.27e−6、roUx 4.78e−1 → 3.51e−3、roUy 8.29e−2 → 1.03e−3、roe 4.74e2 → 3.68、roY0・roY1 も 2.1 桁)。
- 本段の残差の推移 (outer_end の行、log10 の線形回帰の傾き):

  | 列 | [0, 5000) | [5000, 20000) | [20000, 40000) | [40000, 60000) |
  |---|---|---|---|---|
  | rms_ro (min / max / 傾き dec/1e4) | 1.27e−6 / 4.00e−6 / −0.251 | 1.26e−6 / 2.73e−6 / +0.003 | 1.33e−6 / 2.72e−6 / −0.003 | 1.26e−6 / 2.81e−6 / −0.002 |
  | rms_roUx | 1.68e−3 / 4.91e−3 / −0.151 | 1.61e−3 / 4.33e−3 / +0.003 | 1.93e−3 / 4.30e−3 / −0.003 | 1.75e−3 / 4.37e−3 / −0.002 |
  | rms_roe | 2.08 / 6.77 / −0.280 | 2.06 / 4.47 / +0.003 | 2.17 / 4.43 / −0.003 | 2.06 / 4.58 / −0.002 |

  → 本段の 5,000 step 以降は約 2 倍の幅の有界な振動で、傾きは 0 (±0.003 dec/1e4)。
- 予備の資格の量 (`check_quasisteady.py --series-csv`、既定の閾値、60 出力): 出口コア M・流量・軸の M (40・70・94)・η 0.1 の M (40・70・94) の 8 列すべて **STEADY** (drift 0.0 %/tail、fluct 0.0 %、表示の桁で)。
- 判定器は残差の条件で止まり、固定 (pin) は書いていない。NS (run_0951) は投入していない。
- **旧い Euler 参照も同じ**: `run_0174_euler_v5d_M_r1` (生産の δ* チェーンの参照に使ってきた) は、plan `discretization-moc-axis-limit-and-corrector` の V5d の記録で「全 run RUN_RC 0、本段 NOT CONVERGED stalled/plateau」。残差のファイルは掃除で消えていて今は再判定できない。
- 前の諮問 (registration) で「Euler: 全残差の `check_convergence.py` PASS」を資格に入れたのは呼び出し側のブリーフの案に codex が同意した形で、Euler の残差が構造的に 3 桁下がらないことは考慮していなかった (呼び出し側の見落とし)。

## 期待値と出典
- 登録 (§6 V-h370): Euler の資格 = 全残差の `check_convergence.py` PASS、NaN/Inf なし、出口 M・流量・参照場の時系列の `check_quasisteady.py`、参照の時間の変動が抽出の差の判定を左右するなら資格未達。強化 (prelaunch-review M3): 使用域の σ_ref,use ≤ 0.01 mm (NS の 17 出力 × Euler の末尾の窓の全 20 出力) を独立の必要条件。
- AGENTS.md: 残差のプラトーは収束ではない。量の定常化を併せて確かめる。リミットサイクルの可能性。

## 仮説 (呼び出し側の案)
- H1: 2 次・リミッタつきの Euler のこの格子では、残差の床は有界な振動 (リミッタの切り替えなど) で、場の量は定常。設計の入力として必要なのは「抽出への参照の時間の変動の効き」(σ_ref,use) で、それは登録の強化で直接測る。
- 案 A: 登録を改める (結果を見た後なので、その旨を明記): Euler の参照の資格の残差の条件を「本段の後半の残差が有界 (傾き 0 の有意性・振幅の上限) + 全量 STEADY」に置き換え、最終の資格は σ_ref,use ≤ 0.01 mm (17 × 20、既に必要条件) で決める。
- 案 B: 登録どおり FAIL とし、Euler の設定を変える (例: リミッタの凍結・cfl・1 次) か、別の判定の道具で残差の床を測る (`rounding_floor.py`、倍精度の対照) まで NS を回さない。
- 案 C: その他。

## 問い
1. 案 A・B・C のどれか。A なら、置き換える残差の条件の具体 (窓・傾きの有意性・振幅の上限) を、今回の観測に合わせて作る形にならないよう、どう決めるべきか。
2. 案 A のとき、NS (run_0951) を今の固定の出力 (res_60000) で投入してよいか (最終の資格は NS の後の σ_ref,use が決める)。
3. 旧い参照の腕 A (run_0174) も同じ残差の状態なので、A/B の解釈に影響するか。

## 読んでよいファイル
- `plans/active/tooling-nozzle-isothermal-wall-chain.md` §6 V-h370・§5.1 #36e〜#36i
- `case/45.isobutane_m6_d155/h370_judge.py` (`euler_prequalify`)、`solver_density_cuda/tools/check_convergence.py`・`check_quasisteady.py`
- `plans/accepted/discretization-moc-axis-limit-and-corrector.md` (V5d の記録)
- `notes/reviews/2026-10-14-h370-registration-diagnose.md`・`2026-10-14-h370-prelaunch-review-diagnose.md`
