# 諮問ブリーフ: ③ベル (M) の δ* 補正 — 理論式を NS で較正する設計 (#7b)

- 日付: 2026-10-12
- 依頼元: 主セッション (Opus)。エスカレーション条件 1 (§4・§6 の設計)。このセッションの諮問は diagnostician (Fable)
- 作業ツリー: `/home/sano/work/forge-thruster` (HEAD `17d8442a`)。NS・Euler の最終場は手元に写した (`case/68.bell_two_mode/run_0018_tp_m_euler_l100_pin/`・`run_0019_tp_m_ns_l100/`、git 追跡外)
- 対象: plan `plans/active/tooling-nozzle-bell-two-mode.md` §4.4 の「δ* 補正の扱い (ユーザ決定 2026-10-12)」と「δ* と物理寸法」・§5.1 #7・§6 V3。前回の記録 `notes/reviews/2026-10-12-bell-m-deltastar-route-diagnostician.md`

## 1. 前提 (ユーザ決定 2026-10-12、覆さない)

- δ* 補正は (M) の**設計者の選択肢** (補正する / しない)。
- 補正の方法は **δ* の理論式 (積分法など) を NS の計算結果で較正する形**。前回のあなたの推奨 (C を生産・A を棚上げ) のうち、(A) の棚上げは不採用。補正しない (C) も選択肢として残る。
- 補正で狙いに近づくのは壁圧・壁沿いのマッハ数・出口の分布 (1 次) であることを説明したうえでの決定 (④・海面作動では壁圧が効く)。

## 2. 観測事実

**(M) の NS 4 本 (補正なし、AWS B、全ゲート PASS — 壁解像 y1+ > 1 は 0 %、推力系 ALL STEADY、残差 STALLED 併記)**: (D) の MOO と同じ格子族・recipe。

| L_k | run | η_exit (NS) | MOC η_cs (非粘性) | (D) の Pareto (同じ L で内挿、記録のみ) |
| --- | --- | --- | --- | --- |
| L70 (15.86) | run_0020 | 0.970936 | 0.977298 | ≈ 0.96973 |
| L80 (18.13) | run_0021 | 0.977196 | 0.984210 | ≈ 0.97556 |
| L90 (20.39) | run_0022 | 0.980506 | 0.988133 | ≈ 0.97961 |
| L100 (22.66) | run_0019 | 0.982155 | 0.990345 | Pareto の最長は 21.69 で 0.98083 |

- NS の ṁ/ṁ_1D は 4 本とも 0.99075 (スロートの境界層で Euler 0.99263 より −0.19 %)。

**δ*(x) の抽出** (`metrics/deltastar.py` の `deltastar_from_core_matched_euler`、band 法、NS = run_0019、Euler 基準 = run_0018 [同じピン版 L100 の壁の Euler]、全 1201 断面で ok、band_rms 最大 4.5e-3 [線 5e-3]、core_rms 最大 0.129): 
x 0.5 で 0.00178 r_t (r_w 1.134)、1.0 で 0.00317、2.0 で 0.00620、4.0 で 0.01374、6.0 で 0.02314、8.0 で 0.03260、12.0 で 0.05683、16.0 で 0.08222、20.0 で 0.11175、22.5 で 0.13318 (r_w 7.057) → **出口で δ*_e/r_e ≈ 1.9 %**。前回の見積もりで置いた 0.5 % の約 4 倍。抽出には③の run の prepare_info に `x_E` が無くて流量比の段で止まる欠け (手でその段を外して回した) — 道具の小さな直しが要る。
- 抽出値の配列は `/tmp/claude-1000/-home-sano-work-forge/9e48a22d-63d9-46f8-a767-e8fdf4578030/scratchpad/dstar_l100.npz` (x・d・rw)。

**理論式の手持ち**: `feedback/deltastar.py` (乱流平板の相関 + Eckert の参照温度、断熱、CPG の γ・cp を引数に取る `dstar_flatplate(s, M, Pt, Tt, γ, cp)`・`deltastar_offset`)、`feedback/deltastar_integral.py` (CONTUR 型の積分法、CPG。plan §4.4 は TP 化を求めている)。①風洞の教訓: 積分法・経験式はスロートの δ* を NS の実効値の 3〜12 倍に出した。

## 3. 問い (推奨を 1 つずつ)

1. **前回の「補正しない損は 2 次で ≈ 1e-4」の見積もりは、δ*_e/r_e 1.9 % でどう変わるか**。補正の効きの見込みと、それを確かめる測定 (例: L100 で NS の δ* をそのまま使った補正壁を NS で評価する一発測定) の事前登録。
2. **理論式と較正の形**: どの理論式を基にするか (平板の相関 / 積分法 — TP 化の要否)、較正の形 (定数倍 / x の関数 / 無次元量の関数)、較正に使う NS の本数と形 (例: L100 で較正し L70 で検証)、合格の線 (検証形状で δ* の差 / 補正後の η や壁圧の差)。
3. **非粘性壁の ε と出口半径の合わせ方、スロートの扱い**: 前回の推奨 (ε_inv = (√ε − δ*_e/r_t)²、スロートは補正しない) を、補正が設計者の選択肢になった前提で見直すか。生成器を ε_inv で回し直す費用 (1 本 約 25 分) を含めて。
4. **V3 の比較の組み方**: (M) は「補正あり」「補正なし」の 2 本立てで (D) と比べるか、どちらかを主にするか。設計時間の表にどう入れるか (補正ありは NS の較正が一度きりの開発費か、設計ごとの費用か)。
5. 次の一手 1 つとその判定。

## 4. 読んでよいもの

- `plans/active/tooling-nozzle-bell-two-mode.md` (§4.2・§4.4・§4.5・§5.1 #7・#8・§6 V3)
- `notes/reviews/2026-10-12-bell-m-deltastar-route-diagnostician.md` (前回)
- `design/forge_design/metrics/deltastar.py` (先頭 30 行、`deltastar_from_core_matched_euler` 394〜、`massflow_ratio` 111〜)
- `design/forge_design/feedback/deltastar.py`、`design/forge_design/feedback/deltastar_integral.py` (先頭 60 行・200〜260 行)
- `case/68.bell_two_mode/rao_m/rao_pin_lk.csv`、`case/68.bell_two_mode/README.md` の run 一覧
- python で抽出値 (上の npz) や理論式を呼んで確かめてよい。ファイルは書かない。禁止: `*.log`, `res_*.h5` の丸読み, `*.pdf`, `plans/README.md`
