# 諮問ブリーフ: M5.5 の NS (case/71 run_0002_ns) が停止規則に届かずディスク満杯で 184k で停止 — 次の一手 (2026-10-15)

- 依頼者: 主セッション `M6ノズルの設計` (エスカレーション条件 2・4: 未収束の対処、plan に無い再実行の前)
- plan: `plans/active/tooling-nozzle-card-local-mesh.md` §5.1 #3a・#3b、関連 `plans/active/tooling-nozzle-cooled-wall-ns-settling.md` (case/70 の冷却壁の入口の直管の質量の過渡で停止規則が未到達だった件; この作業ツリーに無ければ main の版)
- 入力: `case/71.m55_combustion_wt/problem_m55_final.yaml` (M5.5・燃焼ガス・Pt 3.5 MPa・Tt 1250 K・等温壁 350 K・r_t 0.0999 m・入口の直管 0.5 m + 縮流部 0.5 m・拡大部 8 m、手組みの格子 ni 4800 × nj 92 = 441,600 点、`axis_cap_frac` 0.03)、Euler 参照 `run_0001_euler` (res_12000、出口の軸の M 5.5002、残差 2.1〜2.5 桁で停滞 = `check_convergence --segment` NOT CONVERGED)。バイナリは FP64 (2fd36b0c)、生成器は feature/card-local-mesh b0d528bf。

## 観測事実

- `deltastar_loop --init-integral --stages full` で S1 soft → S2 mid → 本段 (停止規則: 窓 2 万 step・質量の収支の平均 ≤ 0.05 %・振幅 ≤ 0.5 %・θ_r (x 0.4/0.7/0.95 L) の変化 ≤ 0.05 %・振幅 ≤ 0.5 %、到達の後 4 万 step 保持、上限 20 万)。本段は約 14 ms/step。
- 本段の 183,998 step で AWS A のディスクが満杯になり (全場の出力 1 枚 100 MB × 1,000 step ごと、主セッションの見積もりの誤り)、forge が Write failed、ドライバも OSError で終了。停止規則の状態は RUNNING のまま (到達なし)。途中の全場の出力は 1 万 step ごとと 180,000〜183,000 だけ残した。
- 停止規則の系列 (`stop_rule.json` の rows、抜粋):

  | step | 質量の収支 r | ṁ_in | ṁ_out | θ_r 0.4 | θ_r 0.7 | θ_r 0.95 | 出口コア M |
  | --- | --- | --- | --- | --- | --- | --- | --- |
  | 50,000 | 0.01672 | 122.844 | 120.790 | 0.03893 | 0.06177 | 0.06452 | 5.52532 |
  | 100,000 | 0.01046 | 122.193 | 120.914 | 0.03962 | 0.06288 | 0.06518 | 5.52396 |
  | 150,000 | 0.00551 | 121.740 | 121.070 | 0.04056 | 0.06435 | 0.06601 | 5.52316 |
  | 180,000 | 0.00340 | 121.586 | 121.172 | 0.04109 | 0.06510 | 0.06650 | 5.52212 |
  | 182,000 | 0.00331 | 121.578 | 121.175 | 0.04112 | 0.06514 | 0.06652 | 5.52214 |

  r は 3 万 step で約 0.62 倍 (時定数 約 6 万 step) で単調に減り、ṁ_in が下がり ṁ_out が上がって近づいている。θ_r は単調に増え続け、最後の窓の変化 D は 0.4〜0.75 % (条件 0.05 %)。
- 残差 (本段): ρ 1.3 桁・ρE 1.7 桁・運動量 1.9〜2.0 桁・k・ω 1.7〜2.4 桁の低下で、ρ・ρE は直近 2 万 step で横ばい (+0.1 桁)。
- NaN・Inf は全出力で 0 (停止規則の監査)。

## 期待値と出典

- 停止規則の条件 (上) と、plan cooled-wall-ns-settling §6 S-stop の「冷却壁の case/70 は入口の直管の質量の過渡で未到達」。
- 出口の M は δ* の補正の前なので 5.5 からずれてよい (次の pass で直す)。

## 呼び出し側の仮説・案

1. 遅い過渡は入口の直管の質量の再分配 (冷却壁で密度が上がる境界層の発達) で、時定数 約 6 万 step。183,000 の場から同じ格子・同じ設定で続け (`restart_field.py`、出力の間隔を 5,000 step に)、r ≤ 0.0005 まで約 12 万 step (約 30 分) 伸ばせば停止規則に届く見込み。
2. 続きの run の停止規則は窓の数え直しになる (到達・保持の判定は続きの run の中で)。判定区間は続きの run の本段だけ (同じ設定なので本段 + 続きを連結してもよいか)。

## 問い

1. 案 1 (続けて伸ばす) でよいか。r の減り方 (指数に近い) から到達の見込みと上限 (例 15 万 step) をどう置くか。伸ばす前に切り分けるべき原因 (入口の BC・初期場・冷却壁) はあるか。
2. 続きの run の判定 (停止規則・収束の区間・`--from-floor`) の置き方。
3. 出口コア M 5.522 (+0.4 %) と θ_r の値は、最初の pass として妥当か (次の δ* の補正に渡してよい状態か) — 判断に要るものは何か。

## 読んでよいファイル

- `case/71.m55_combustion_wt/` の README・問題 YAML・`_mesh/`
- `design/forge_design/evaluate/stop_rule.py`、`design/forge_design/evaluate/runner_axismach.py` (`run_staged_ns`)、`design/forge_design/feedback/deltastar_loop.py`
- `plans/active/tooling-nozzle-card-local-mesh.md`、`plans/active/tooling-nozzle-cooled-wall-ns-settling.md` (あれば)
- run の記録は AWS A にあり手元に無い (上の表が抜粋)。
