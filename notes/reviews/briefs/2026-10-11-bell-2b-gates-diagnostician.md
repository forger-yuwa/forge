# 諮問ブリーフ: ③ベル共通評価器 (#2b) のゲートが基準ケースの全流れ解を通さない件

- 日付: 2026-10-11
- 依頼元: 主セッション (Opus)。エスカレーション条件 3 (事前に書いた合格条件との比較が判定不能・不合格) と 4 (plan に無い修正の前)。ユーザ指示によりこのセッションの諮問は diagnostician (Fable)
- 作業ツリー: `/home/sano/work/forge-thruster` (HEAD `b630ca99`)
- 対象 plan: `plans/active/tooling-nozzle-bell-two-mode.md` §4.1 (推力・ゲート (i)〜(v))、§5.1 #2b

## 1. 観測事実 (主セッションで再確認したもの: 試験 ALL PASS、run_0003 の `thrust_history.csv` 末尾、4 run の `CONVERGENCE_VERDICT.txt`)

- 実装: `metrics/extract.py` の `thrust_cv` (全境界の Π_x、node の半辺の r 重み、外向き法線を 2 通りで照合、壁の traction は node で −twall)、`evaluate/runner.py` の `run_staged_tp`・`continue_tp` (同一メッシュで restart_field + StageManifest)、`opt/driver.py` の `tp_gates` (i)〜(v)、`opt/polish.py` の `polish_run_tp`。試験 `design/tests/run_bell_cv_gate_tests.py` ALL PASS。
- 基準ケース (LOX/LH2・10 MPa・O/F 6・凍結 TP、TOP 形状 θn 32°・θe 9°・L 18 r_t、ε 50、断熱、真空、推奨 TP レシピ):

| run | 内容 | CV 残差 | η (出口面) | η (壁) | ṁ/ṁ_1D |
| --- | --- | --- | --- | --- | --- |
| `case/68.bell_two_mode/run_0002_tp_staged_cv` | soft 2000 → mid 2000 → 本段 12000、561×121 | −6.6e-4 | 0.97978 | — | 0.99341 |
| `run_0003_tp_cv_continue` | run_0002 の本段を同一メッシュで +24000 (本段計 36000) | **−5.2e-7** | 0.979220 | 0.979220 | 0.99255 |
| `run_0004_tp_staged_wallfine` | 1201×121 (AR 979.5)、スロートの第一セル 0.97 µm、本段 36000 | −7.3e-5 | 0.97921 | — | 0.99264 |

- NaN/Inf は 3 本とも 0。`check_convergence --segment`: run_0003 (区間 S3_main → S4_main_p、本段 36000) は **NOT CONVERGED (stalled)** — rms_ro 2.5・rms_roUy 2.2・rms_roe 2.5・rms_roOmega 1.6 桁で横ばい。run_0004 (S3_main) も stalled (rms_roUy 2.7・rms_roOmega 1.7 桁、他は 3 桁以上)。
- `check_quasisteady --series-csv thrust_history.csv`: 壁会計・出口面会計・η・ṁ は **STEADY** (単調、漸近値までの残り 0.03 % 以下)。**cv_residual は DRIFTING** (末尾窓の相対ドリフト run_0003 225 %、run_0004 96 %)。
- 局所 y1+ (`check_wall_resolution.py --groups wall`): run_0003 は最大 5.68 (スロート全点 > 1)、run_0004 は最大 2.55 (入口直管の始端)、スロート最大 1.60 で点の 68.7 % が > 1。実装者の探索: 現行の格子 (半径方向の線が x 一定で拡大部の近壁セルにスキューがあり AR 5000 の例外に当たらない) では ni ≤ 1201 で y1+ ≤ 1 と AR ≤ 1000 が両立しない。ni 1601 で y1+ ≤ 0.91 なら AR ≈ 1400。
- CPG の case/40 (出口 `outlet_statPress`、Ps 2e4 Pa) の 16 run で CV 残差 +3.4e-3〜+5.5e-3。出口の壁際の節点 P は 6〜8e4 Pa (実装者は境界流束と節点状態の不一致を疑うが未確認)。配線回帰 (旧 `thrust_metrics` と新 `thrust_cv` 出口面会計の η 差) は最大 1.5e-5。
- 収縮部の P・T の超過 (2a の引き継ぎ (e)): 収束後は P > Pc が 0 ノード、T > Tc が最大 +0.35〜0.51 K (壁から 50 µm 以内、x/r_t −3.5〜−1.8)。

## 2. 前例 (リポジトリ内)

- ⑤ SERN では「残差プラトーは受理の阻害要因としない。力は 7 桁で定常・派生量は ALL STEADY・場の動きは再現性幅の 3〜62 倍、`NOT CONVERGED` は消さず併記」とユーザが決定 (`plans/active/tooling-nozzle-sern-3d.md` の §4.39 付近、R5g の行)。
- AGENTS.md「収束確認」: 残差プラトーは収束ではなく、積分量の定常化と場の発達を併せて確認する。プラトーが停滞か有界振動かは測る (メモリ residual-plateau-bounded-oscillation: 保存量の同期・傾きの有意性・振幅・CFL を伸ばす)。

## 3. 問い (各 1 つに絞って推奨)

- **(B) ゲート (i)**: 推奨 TP レシピで本段が 2.2〜2.7 桁 (ω 1.6〜1.7) で横ばい、推力量は STEADY。(a) 横ばいの原因を切り分ける (何を 1 つ測るか)、(b) ゲート (i) を SERN と同型の「プラトー併記 + 推力系 STEADY + 場の動きの上限」に変える、のどちらを先にするか。(b) なら場の動きの上限の測り方。
- **(A) cv_residual の定常判定**: 0 を横切る系列に相対ドリフトは使えない。絶対値の帯 ($F_{1D}$ 単位) にするなら値と窓。
- **(C) スロート y1+ ≤ 1 と AR**: 壁法線に沿った近壁格子にして AR 5000 の例外を使う / ni ≥ 2000 程度 / 目標を緩める (緩めるなら格子感度の裏付けが要る — AGENTS.md「壁解像確認」)、から 1 つ。
- **(D) η の定義**: 出口面会計と壁会計のどちらを主にするか。
- **(E) ゲート (ii) の閾値** (今は drift 5 %・osc 10 % の既定): η を 1e-3 の精度で比べる用途に対して。
- **(F) `space.pRef`**: 推奨レシピ §5 は軸対称で動作静圧を書くよう求めるが、ノズルは 1e7→1.4e4 Pa と幅がある。書くなら何の値か、書かないならその根拠。
- **(G) CPG の CV 残差 3〜5e-3**: TP 側の受理に関係するか、別件として記録だけでよいか。

## 4. 読んでよいもの

- `plans/active/tooling-nozzle-bell-two-mode.md`
- `design/forge_design/metrics/extract.py` (`thrust_cv`:205 付近)、`design/forge_design/opt/driver.py` (97〜212 行のゲート)、`design/forge_design/evaluate/runner.py` (`run_staged_tp`:585 付近、TP の config)
- `case/68.bell_two_mode/README.md`、各 run の `CONVERGENCE_VERDICT.txt`・`CONVERGENCE_SEGMENT_VERDICT.txt`・`QUASISTEADY_THRUST_VERDICT.txt`・`WALL_RESOLUTION_VERDICT.txt`・`GATES.json`・`thrust_history.csv`・`solverConfig.yaml`
- `procedures/recommended-settings.md` (§2・§5)、`procedures/solver-settings.md` (pRef の節)
- `residual_history.csv` は統計量だけ python で取るなら可 (全量を出力しない)
- 禁止: `*.log`, `res_*.h5` の全量ダンプ, `*.pdf`, `plans/README.md`。ファイルは変更しない。forge は起動しない
