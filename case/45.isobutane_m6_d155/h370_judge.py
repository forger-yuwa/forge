"""plan tooling-nozzle-isothermal-wall-chain §5.1 #36f・#36g・#36i・§6「V-h370 の事前登録」(§4.9 codex diagnose 3、判定の基準の諮問
`notes/reviews/2026-10-14-h370-judge-criteria-diagnose.md`、投入前のレビューの諮問 `notes/reviews/2026-10-14-h370-prelaunch-review-diagnose.md`)
の判定。B の case dir で回す (h370_run.py が呼ぶ)。壁は更新しない (固定壁の pass)。

--euler-prequalify: (i) Euler 参照 run_0950_h370_euler_ref の**予備の資格** (prequalify; §6 V-h370 の「Euler」、§5.1 #36g (1)):
  - RUN_RC 0、投入前の検査 (H370_PRECHECK.json) の合格、本段が上限 60,000 step まで回った (最後の出力 res_60000.h5)。
  - NaN/Inf なし: 全段の残差の履歴 (residual_history_{soft,mid,main}.csv) の rms_* と、本段の全出力の VALUE の浮動小数点の量。ρ・P・T > 0。
  - 全残差: check_convergence.py --segment --drop 3 (stage_manifest の最後の区間 = 方程式・BC・離散化が同じ本段の区間、全列) の PASS。
  - 出口 M・流量・参照場の時系列 (本段の全出力) の check_quasisteady.py の VERDICT (CSV 系列モード、閾値は既定のまま): 量は
    出口コア M・無次元の流量 (stop_rule.recorded_quantities = nozzle_report.metrics・qs_series.py と同じ定義) と、参照場として軸 (η = 0) と
    η = 0.1 の線の M を x = 40・70・94 r_t で (nozzle_report.eta_line、停止規則の断面と同じ位置)。すべて STEADY であること
    (OSCILLATING も参照の時間の変動なので未達)。nozzle_report.metrics の η = 0.1 の線の評価量 (wave・overshoot・slope・range、%) は
    平均が 0 に近く相対の閾値が意味を持たないので記録だけ。閾値を細かくした (1e-4) 判定も記録だけ。
  これは 6 点の M・出口 M・流量の既定の相対閾値の準定常で、抽出に使う参照場全体の資格には足りない (codex の Major)。**最終の資格**は
  NS の判定の σ_ref・σ_ref,use (下の 5) で決め、σ_ref を測れない・判定不能・参照の時間の変動が A/B の判定を左右する (参照の未整定)・
  σ_ref,use > 0.01 mm・A/B が判定不能なら `update_conditions_met` を偽にする。
  満たせば本段の最後の出力とその sha256 を _band_ab/half_inlet/h370_euler_ref_pin.json に固定する (既にあって違えば止める; 消さない)。
  σ_ref の代わりの参照 (§5.1 #36i (M1)) は Euler の末尾の窓 = 本段の最後の 20,000 step [40,000, 60,000] の**全**出力 (出力 1,000 ごと、
  固定の 60,000 を含め 21 出力のうち代わりは 40,000〜59,000 の 20 出力; 出力の間隔と窓から列挙し、在る出力だけから作らない) とその sha。
  1 つでも欠ければ予備の資格を判定不能にする。旧い代表 (5,000 ごとの 40,000・45,000・50,000・55,000) と末尾の窓の全出力の場の差
  (`tail_field_dev`) は記録だけで、資格の判定には使わない。スクリーニングの出力 (check_quasisteady の末尾の窓の始まりと中ほど、停止の出力
  だけに当てる) も記録する。
  → _band_ab/half_inlet/h370_euler_prequalify.json。終了コード 0 予備の資格あり / 1 未達 / 2 判定不能。

(既定): (ii) NS run_0951_h370_ns_b0 の §6 V-h370 の判定 → _band_ab/half_inlet/h370_judge.json:
  1. 記録: 投入前の検査 (H370_PRECHECK.json)・本段の直前の照合 (H370_MAIN_START.json)・基盤の停止 (PRECHECK_FAILED・WATCH_ABORT・
     RESOURCE_STOP・DISK_STOP・WRITE_FAILED)・RUN_PROVENANCE.txt のバイナリと環境変数・段の残差の履歴の非有限。
  2. 停止の状態 (STOP_RULE_STATE.json): 状態・停止の出力 s・到達と中断と保持 (stop_rule.json の hold)・要求の断面と実際の列・停止の出力の sha。
  3. 受け渡しの関門 (deltastar_loop.handoff_gate、HANDOFF_GATE.json): 停止の状態 REACHED・check_convergence.py --segment --drop 3 の PASS・
     保存則・局所 y₁⁺ (check_wall_resolution.py --step s)・2 窓 [s − 40000, s − 20000]・[s − 20000, s] の δ_E・δ_target (全使用域で
     最大 − 最小 ≤ 0.01 mm、窓の平均の差 ≤ 0.01 mm; 入るべき 17 出力を出力の間隔から列挙)・Euler の参照と ω・抽出の全域の被覆・
     使用域の held と欠測なし・次の壁の幾何。上書き (allow_not_reached・allow_unconverged) は使わない。合格なら関門が抽出の表を run の直下に
     写す (次の pass は自動でつながない)。
  4. check_quasisteady.py の VERDICT (記録; 既定の相対閾値で寸法の予算の合格を代用しない): 2 窓の r・θ_r (断面ごと)・出口コア M・流量・
     δ_E・δ_target (x = 40・70・94 r_t) の時系列。
  5. 抽出の参照の A/B と σ_ref (同じ NS の保存系列に対して Euler の参照だけを変える; §5.1 #36g (2)〜(5)):
     - 判定窓の完全性: 停止の時刻 s・出力の間隔 (stop_rule.json の settings.out_step_interval、2,500)・2 窓から入るべき出力
       (stop_rule.required_window_steps、17 出力) を列挙し、stop_rule.json の行と res の有無と照合する。行が無くても res があれば
       抽出し直して補い、res が無ければ判定不能。
     - 比べる列 F を**抽出の前に**座標だけで決める: NS の使用域の列 (stop_rule.json の delta_series.x) のうち、腕 A (旧い
       run_0174_euler_v5d_M_r1、最後の出力を sha で固定) と腕 B (固定した run_0950) の Euler の補間できる範囲 (共通の関数
       deltastar.euler_cover_range、受け渡しの抽出と同じ端の許容つき = 同一の端点の丸めを吸収する許容 (端の外側の微小区間への
       定数の延長、端の列の間隔 × 1e−4)) の両方に入る列。F の中の非有限・欠測・held は比較を
       合格にしない (後から有限の列だけを選ばない)。旧い参照の範囲の外 (出口の区間) は F に入らないので比べない (外挿しない)。
     - 出力 k ごと: 各腕の合否 pf (抽出の例外なし・F の中に held・非有限・欠測なし)、両腕とも pf のとき
       D_k = F の max|δ_E,B − δ_E,A| [mm]。抽出の条件は受け渡しと同じ (edge・一定 x の断面への補間・端は欠測・ω 0.5・平滑化の λ・
       昇格の上限; stop_rule.json の extract)。
     - σ_ref (腕 B の参照の時間の変動; §5.1 #36i (M1)・(M3)): 必要な全出力 k (17) について、固定した出力の代わりに Euler の末尾の窓の
       **全**出力 (40,000〜59,000 の 20 出力 = `sigma_ref_required`、登録の定数から列挙) で抽出し直し (17 × 20 = 340 抽出)、
       σ_ref = max_k max_e max_F |δ_E(k, e) − δ_E(k, 固定)| [mm] (比べる列 F)、
       σ_ref,use = max_k max_e max_U |δ_E(k, e) − δ_E(k, 固定)| [mm] (使用域の全体 U = 固定の参照で抽出できる使用域の全列 = 受け渡しの
       関門の coverage_check と同じ使用域 [0, x_e]; F ⊂ U)。最大を与えた (NS の出力・Euler の出力・x) を両方とも記録する。
       抽出の前に、予備の資格の記録の代わりの参照が 20 出力ちょうどで、各出力が在り sha が記録と同じことを照合する (欠け・不一致は
       σ_ref を測れない → 最終の資格なし; 抽出は回さない)。代わりの参照で F・使用域の held・欠測・非有限の列 (列ごとの印) や合否が
       変われば σ_ref は判定不能、抽出が失敗すれば測れない。停止の出力 × スクリーニングの 2 出力の差は「スクリーニング」として記録だけ
       (σ_ref の上限とは扱わない)。
     - 分類 (classify_ab): 必要な出力の欠け・抽出の失敗・σ_ref を測れない・判定不能 → 判定不能。どれかの出力で pf_A ≠ pf_B → 代わりの
       参照の全部でもその違いが保たれれば 結果 B (保たれなければ判定不能)。両腕とも pf でない出力 (F の中の欠測) → 判定不能。
       それ以外: max_k D_k + σ_ref ≤ 0.01 mm → 結果 A、min_k D_k − σ_ref > 0.01 mm → 結果 B、σ_ref を足し引きして初めて 0.01 mm を
       またぐ → 判定不能 (参照の未整定)、D_k だけで窓の中で 0.01 mm をまたぐ → 判定不能 (時間の変動)。
     - 参照の最終の資格 (reference_final; §5.1 #36i (M2)・(M3)): 次の**すべて**のとき資格あり — σ_ref・σ_ref,use を測れて判定できる、
       **σ_ref,use ≤ 0.01 mm** (A/B の分類とは独立の必要条件: 旧い参照との差が明確なことと、新しい参照が設計の入力として落ち着いている
       ことは別)、A/B が結果 A か結果 B (B_D・B_pf) に分類された (σ_ref が分類を左右しない)。σ_ref・σ_ref,use を測れない・判定不能・
       σ_ref,use > 0.01 mm・参照の未整定 → 資格なし。A/B が判定不能 (窓の中で 0.01 mm をまたぐ時間の変動 U_time・出力の欠け・抽出の
       失敗・F の中の欠測・合否の違いを確かめられない) → 資格を確かめられない。どちらも更新に進めない。
     停止の出力 s では両腕の全体の抽出 (run の H370_AB/{A,B}) で帯の位置 (y_b・y*)・感度・held・δ_target の差・次の壁 (腕 B の表、腕 A は
     共通の列だけ腕 A の δ_target に替えた表 = 出口の区間を延ばさない) の真のスロートを記録する。
     抽出の回数と所要の目安 (V-h370、§5.1 #36i): 腕 A 17 + σ_ref 17 × 20 = 340 + スクリーニング 2 (+ 記録の無い出力の腕 B の抽出し直し)
     = 359 回を 1 つの作業の一覧にしてプロセス並列 (`run_snaps`、既定の上限 3; --jobs 1 の直列も子プロセス 1 つで回す: BLAS のスレッド数で
     丸めが変わる (実寸で δ_E 最大 6.9e−18 r_t) ので、子の BLAS を 1 スレッドにそろえて直列と並列をビット同一にする) で回し、停止の出力の全体の抽出
     2 回はこのプロセスで直列。1 抽出は
     実寸の合成の場で 116〜123 s (直列なら 361 回 ≈ 11.6〜12.3 h)。手元 (12 論理 CPU、他の負荷あり) の 3 本同時は 1 本 157 s (1 本だけの
     とき 123 s) で 359 × 157 / 3 + 2 × 123 ≈ 5.3 h。B (g5.xlarge の 4 vCPU = 2 物理コア × 2) では並列の伸びがこれより小さいことがある。
     1 抽出のプロセスの最大メモリは 181〜187 MiB (SNAP_PEAK_MIB 256 で見積もる)。抽出の本体 (帯の梯子 = band_local_deficit) は Euler の参照に
     依る (deltastar.py の euler_q_at を梯子の各段が呼ぶ) ので、代わりの参照の間で使い回せる NS 側の処理は読み込みと一定 x への写像だけで
     1 抽出の約 0.05 % (0.06 s) — 抽出器を変えずに分けられる分は無視できるので使い回しはしない。
  6. 設計を終える条件の量 (記録だけ; 最初の pass の判定には使わない): 使用域の未緩和の max|δ_E − δ_in|。
  判定: 「更新に進める条件」= 関門の proceed かつ 参照の最終の資格 (と 1 の基盤の不成立が無いこと)、抽出の参照の A/B は
  結果 A / 結果 B / 判定不能。
  **診断専用の run** (印 H370_DIAG_ONLY.json か投入前の検査の phase ns-diag; §6 V-h370-diag) はこの判定をせず (判定不能、rc 2)、受け渡しの
  関門を呼ばない。

--diag-sens (§5.1 #36l・§6「V-h370-diag の事前登録」、codex diagnose `notes/reviews/2026-10-14-h370-diag-ns-diagnose.md`、CFD 0 step):
  診断の NS (h370_run.py --phase ns-diag の run_0951) で Euler の参照の入れ替えへの δ_E の感度を測る → _band_ab/half_inlet/h370_diag_sens.json
  (配列は h370_diag_sens_arrays.npz)。受け渡しの関門・正式の判定・壁の更新は呼ばない。update_conditions_met は常に偽。
  - 入力: 停止の時刻 s (REACHED は停止の出力、CENSORED は最後の出力 = 「打切りの時点の診断の窓」) の [s − 40,000, s] の 2,500 ごとの全 17 出力
    (登録の窓・間隔から列挙し、在ることと終了の時点の sha の一覧 H370_DIAG_SHA256.json との一致を確かめる; 欠ければ判定不能、縮めない)。
    Euler の参照は _band_ab/half_inlet/h370_diag_euler_ref.json に固定した腕 A = run_0952_h370_euler_cfl2/res_40000.h5・腕 B =
    run_0953_h370_euler_cfl05/res_40000.h5 (sha を照合)。補助 = 診断の NS の参照 run_0950 の res_60000 (記録だけ)。抽出の条件は停止規則の記録の
    extract (= #36e で固定した受け渡しの条件: edge・remap_constant_x・端は欠測・ω 0.5・平滑化の λ 1・ノット 2・昇格の上限 ×1000) で、#36e と違えば
    判定不能。抽出器・見張りの版 (forge_design の該当ファイルの sha) が投入の時点と違えば判定不能。
  - 比べる列 F は抽出の前に座標と腕 A・B の Euler の補間できる範囲から決める (fixed_columns)。帯の選び方 (edge) は固定、選ばれる帯の位置は
    参照ごとに抽出器が計算し直す (diag_snapshot が y_b・y* を記録)。
  - 主の量 D = max over (17 出力 × F) |δ_E,B − δ_E,A| [mm] (未緩和の δ_E)。分類 (classify_diag; classify_ab の時刻ごとの U_time・σ_ref は
    使わない): 入力の不足・抽出の例外 → 判定不能、入れ替えで合否 (F の中の held・欠測・非有限) が変わる → 結果 B、D > 0.01 mm → 結果 B
    (最大で決める)、両腕とも合否 否の出力が残る → 判定不能、それ以外 (全入力が有効・D ≤ 0.01 mm・合否が同じ) → 結果 A。
  - 抽出は 17 × 3 = 51 回 (run_snaps、既定の並列 3)。終了コード 0 結果 A・結果 B / 2 判定不能。

--diag-sens-fixed-band (§5.1 #36n・§6「V-h370-diag の帯を共通にした A/B の事前登録」、codex diagnose
  `notes/reviews/2026-10-14-h370-diag-result-diagnose.md`、CFD 0 step): 帯の位置を共通にした参照の入れ替えの A/B → _band_ab/half_inlet/
  h370_diag_sens_fixedband.json (配列は h370_diag_sens_fixedband.npz)。原因を分ける診断。受け渡しの関門・正式の判定・壁の更新は呼ばない。
  - 入力の照合は --diag-sens と同じ (diag_inputs) に、--diag-sens の記録 (h370_diag_sens.json・h370_diag_sens_arrays.npz) との照合 (npz の sha・
    17 出力・使用域の 2,938 列・比べる列 F・参照の sha・抽出の条件・版) を足す。
  - 共通の帯: npz の腕 A の y_b(k, x) を両腕に当てる。固定は deltastar.edge_band_positions の戻りの y_b の差し替え (fixed_band; 抽出器はその
    y_b を各列の band_local_deficit(y_b_fixed=...) に渡す。抽出器の本体は変えない)。npz は使用域 [0, x_w] の列だけなので、使用域の外の列
    (縮流部・入口の直管) は腕 A が同じ出力で選ぶ帯を使い、腕 B には腕 A が使った全列の帯を渡す (腕 A → 腕 B の 2 段、17 + 17 = 34 抽出)。
  - 腕 A の再現: 固定の帯の δ_E が npz の A_delta_E と一致 (同じ有限の列・差の最大 ≤ 1e−12 r_t; 最大差を記録)。不成立なら判定不能。
  - D_fixed = max over (17 出力 × F) |δ_E,B − δ_E,A| [mm] (未緩和の δ_E)。分類 (classify_fixed、登録どおり): 欠測・非有限・抽出の不成立 (F の中の
    held を含む)・腕 A の再現の不成立 → 判定不能、D_fixed > 0.01 mm → 結果 B、それ以外 → 結果 A。区間ごと (0〜1・1〜3・3〜5・5〜10・10〜40・
    40〜94・94〜出口) の D_fixed と --diag-sens の D の最大を並べ、帯の移動 (腕 B が選び直した帯 − 共通の帯) と生の δ_r も残す。
usage: python3 h370_judge.py --euler-prequalify     python3 h370_judge.py [--no-ab] [--jobs N]     python3 h370_judge.py --diag-sens [--jobs N]
       python3 h370_judge.py --diag-sens-fixed-band [--jobs N]
       (--no-ab は A/B と σ_ref を回さない = 参照の最終の資格を確かめないので、更新に進める条件は満たさない)
       (--jobs は抽出の並列数の上限 1〜4、既定 3。実際の並列数は min(上限, 空きメモリ / 1 抽出の最大メモリ SNAP_PEAK_MIB))
       (--case-dir・--design・--tools は試験用; 既定は B の ~/forge-g0cchain、無ければこの作業ツリー)
"""
import argparse
import contextlib
import csv
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
HOOKS = {}            # 試験用の差し替え (handoff_gate・delta_snapshot・extract_and_merge・next_wall_check)。本番では空
PLAN = "plans/active/tooling-nozzle-isothermal-wall-chain.md §5.1 #36f・#36g・§6 V-h370 (§4.9 codex diagnose 3、判定の基準の諮問)"
EU_RUN, NS_RUN = "run_0950_h370_euler_ref", "run_0951_h370_ns_b0"
NS_PROBLEM = "problem_d155_ns_prod_coldmesh_tw370_cgref_h_b0.yaml"
OLD_EULER = "run_0174_euler_v5d_M_r1"
PIN = "_band_ab/half_inlet/h370_euler_ref_pin.json"
QUAL = "_band_ab/half_inlet/h370_euler_prequalify.json"   # 予備の資格の記録 (§5.1 #36g (1); 旧名 h370_euler_qualify.json)
OUT = "_band_ab/half_inlet/h370_judge.json"
EU_STEPS = 60000
EU_OUT = 1000                                 # Euler の本段の出力の間隔 (h370_run.EU_OUT・problem_d155_euler_h_ref.yaml と同じ)
EU_TAIL = 20000                               # σ_ref の Euler の末尾の窓 (本段の最後の 20,000 step、出力 1,000 ごと 21 出力)
SIGMA_REF_REP_EVERY = 5000                    # 記録だけ: 旧い代表の間隔 (40,000・45,000・50,000・55,000)。資格の判定には使わない (§5.1 #36i M1)
JOBS_CAP = 3                                  # 抽出の並列数の上限の既定 (B = g5.xlarge の 4 vCPU を他のセッションの計算と共有する)
JOBS_MAX = 4                                  # 並列数の上限の最大 (g5.xlarge の vCPU)
SNAP_PEAK_MIB = 256                           # 1 抽出のプロセスの最大メモリの見積もり [MiB] (実寸の合成の場 NS 4719 × 160・Euler 2000 × 97
                                              # で 1 抽出を回したプロセスの ru_maxrss 181 MiB、3 抽出を続けた子 187 MiB を切り上げ; 並列数を空きメモリで抑える)
BIN_SHA = "7e95e5bbd8f50c9e"
XS = (40.0, 70.0, 94.0)                       # 参照場・記録の量の位置 [r_t] (停止規則の要求の断面 = x_fracs × x_w)
QUAL_COLS = ["exit_core_M", "mdot_rt"] + [f"M_axis_{x:g}" for x in XS] + [f"M_eta0.1_{x:g}" for x in XS]
REC_COLS = ["wave_eta0.1", "overshoot_eta0.1", "slope_eta0.1", "range_eta0.1"]
QS_TIGHT = ["--drift", "0.0001", "--osc", "0.0001"]       # 記録だけ (既定の tail 0.4)
VERDICT_RE = re.compile(r"^\s+(\S+)\s*:\s(.*?)\s*(STEADY|DRIFTING|OSCILLATING|TRANSIENT-UNSETTLED|NONFINITE)\s*$")
TOL_MM = 0.01                                 # 抽出の参照の A/B (§6 V-h370) と 2 窓の予算
AB_DIR = "H370_AB"
MARKERS = ("PRECHECK_FAILED", "WATCH_ABORT", "RESOURCE_STOP", "DISK_STOP", "WRITE_FAILED")
# --- 診断の感度の判定 (--diag-sens; §5.1 #36l・§6 V-h370-diag)。名前は h370_run.py の DIAG_* と同じ (試験 L が照合) ------------------
PLAN_DIAG = ("plans/active/tooling-nozzle-isothermal-wall-chain.md §5.1 #36l・§6 V-h370-diag (codex diagnose "
             "notes/reviews/2026-10-14-h370-diag-ns-diagnose.md)")
DIAG_REF = "_band_ab/half_inlet/h370_diag_euler_ref.json"   # 診断専用の Euler の参照と感度の腕 A・B の固定 (h370_run.py --phase ns-diag が書く)
DIAG_MARK = "H370_DIAG_ONLY.json"                 # 診断専用の run の印
DIAG_END = "H370_DIAG_END.json"                   # 診断の NS の終了の記録
DIAG_SHA = "H370_DIAG_SHA256.json"                # 終了の時点の run の全ファイルの sha256 の一覧
DIAG_OUT = "_band_ab/half_inlet/h370_diag_sens.json"
DIAG_NPZ = "_band_ab/half_inlet/h370_diag_sens_arrays.npz"   # 出力 × 腕 × 使用域の列の δ_E・δ_target・held・帯の位置 (JSON に入れない配列)
DIAG_WINDOW, DIAG_OUT_INT, DIAG_N_REQ = 20000, 2500, 17     # 登録の窓 (2 窓 × 20,000)・出力の間隔・判定の窓の出力の数
DIAG_ARMS = ("A", "B")                            # 主の判定の腕 (A = run_0952 cfl 2、B = run_0953 cfl 0.5)。補助 aux = run_0950 (記録だけ)
VERSION_FILES = ("evaluate/stop_rule.py", "metrics/deltastar.py", "feedback/deltastar_loop.py", "evaluate/runner_axismach.py",
                 "evaluate/runner.py", "metrics/delta_in_check.py")
# --- 帯を共通にした A/B (--diag-sens-fixed-band; §5.1 #36n・§6「V-h370-diag の帯を共通にした A/B の事前登録」) ----------------------------
PLAN_FIXED = ("plans/active/tooling-nozzle-isothermal-wall-chain.md §5.1 #36n・§6「V-h370-diag の帯を共通にした A/B の事前登録」(codex diagnose "
              "notes/reviews/2026-10-14-h370-diag-result-diagnose.md)")
FIXED_OUT = "_band_ab/half_inlet/h370_diag_sens_fixedband.json"
FIXED_NPZ = "_band_ab/half_inlet/h370_diag_sens_fixedband.npz"
REPRO_TOL_RT = 1e-12                              # 腕 A の再現の許容 [r_t] (浮動小数の丸め程度: 7.7e−11 mm = 0.01 mm の 1e−8。同じ機械・子の
                                                  # BLAS 1 スレッドならビット同一のはず; 最大差は記録する)
FIXED_SECTIONS = ((0.0, 1.0), (1.0, 3.0), (3.0, 5.0), (5.0, 10.0), (10.0, 40.0), (40.0, 94.0), (94.0, None))   # 区間 [lo, hi) [r_t]、最後は出口まで


def now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256_file(p) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def jload(p):
    try:
        return json.loads(Path(p).read_text())
    except (OSError, ValueError):
        return None


def jwrite(p: Path, d) -> None:
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = Path(str(p) + ".tmp")
    tmp.write_text(json.dumps(d, indent=1, ensure_ascii=False, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o)))
    os.replace(tmp, p)


def res_steps(d: Path) -> list:
    return sorted((int(f.stem.split("_")[1]) for f in Path(d).glob("res_[0-9]*.h5") if re.fullmatch(r"res_\d+\.h5", f.name)))


def residual_nonfinite(p: Path):
    """残差の履歴 CSV の rms_* 列の非有限の数と最初の step (tw370_judge.py と同じ)。無ければ None。"""
    if not p.is_file():
        return None
    n_bad, first, rows, last = 0, None, 0, None
    with open(p, newline="") as f:
        r = csv.reader(f)
        head = next(r, [])
        cols = [i for i, c in enumerate(head) if c.strip().startswith("rms_")]
        for row in r:
            if not row:
                continue
            rows += 1
            last = row[0]
            for i in cols:
                try:
                    v = float(row[i])
                except (ValueError, IndexError):
                    v = float("nan")
                if not math.isfinite(v):
                    n_bad += 1
                    if first is None:
                        first = row[0]
    return {"file": p.name, "rows": rows, "rms_cols": len(cols), "nonfinite": n_bad, "first_nonfinite_step": first, "last_step": last}


def run_qs(tools: Path, series_csv: Path, cols: list, extra=()) -> dict:
    """check_quasisteady.py の CSV 系列モード。列ごとの VERDICT を出力の行から読む (行の無い列は missing)。"""
    cmd = [sys.executable, str(tools / "check_quasisteady.py"), "--series-csv", str(series_csv), "--series-cols", ",".join(cols), *extra]
    r = subprocess.run(cmd, capture_output=True, text=True)
    v = {}
    for ln in r.stdout.splitlines():
        m = VERDICT_RE.match(ln)
        if m and m.group(1) in cols:
            v[m.group(1)] = {"verdict": m.group(3), "detail": m.group(2).strip()}
    return {"cmd": " ".join(cmd[1:]), "rc": r.returncode, "verdicts": v, "missing": [c for c in cols if c not in v],
            "overall": next((ln.strip() for ln in r.stdout.splitlines() if ln.startswith("OVERALL")), None),
            "stderr": r.stderr.strip()[-500:]}


def write_series(path: Path, rows: list, cols: list) -> None:
    with open(path, "w") as fh:
        fh.write("step," + ",".join(cols) + "\n")
        for r in rows:
            fh.write(f"{r['step']}," + ",".join(f"{float(r.get(c, float('nan'))):.12g}" for c in cols) + "\n")


# ===================================================================================================== (i) Euler の資格
def euler_series(run: Path, steps: list) -> tuple:
    """本段の出力ごとの量 (資格の量 QUAL_COLS・記録の量 REC_COLS) と、全 VALUE の非有限・ρ・P・T の最小。戻り (行, 場の問題)。"""
    import h5py
    from forge_design.evaluate import stop_rule as SR
    from forge_design.report import nozzle_report as NR
    info = json.loads((run / "prepare_info.json").read_text())
    ni, nj, S, x_E = int(info["mesh"]["ni"]), int(info["mesh"]["nj"]), float(info["scale_m"]), float(info["x_E"])
    with h5py.File(run / "nozzle.h5", "r") as h:
        nc = h["/MESH/COORD"][:].reshape(-1, 3)
    Xf, Rf = (nc[:, 0] / S).reshape(ni, nj), (nc[:, 1] / S).reshape(ni, nj)       # stop_rule.RunGeometry の Xf・Rf と同じ式
    rows, probs = [], []
    for st in steps:
        name = f"res_{st}.h5"
        row = {"step": st}
        with h5py.File(run / name, "r") as h:
            V = {k: np.asarray(d[:]) for k, d in h["VALUE"].items()}
        nf = {k: int(np.count_nonzero(~np.isfinite(a))) for k, a in V.items() if a.dtype.kind == "f"}
        row["nonfinite_all"] = int(sum(nf.values()))
        for k in ("ro", "P", "T"):
            if k not in V:
                probs.append(f"{name}: VALUE/{k} が無い")
                continue
            row[f"{k}_min"] = float(np.nanmin(V[k])) if V[k].size else float("nan")
        if row["nonfinite_all"]:
            probs.append(f"{name}: 非有限 {row['nonfinite_all']} 個 ({ {k: n for k, n in nf.items() if n} })")
            rows.append(row)
            continue
        if any(row.get(f"{k}_min", 1.0) <= 0 for k in ("ro", "P", "T")):
            probs.append(f"{name}: ρ・P・T の最小 {[row.get(f'{k}_min') for k in ('ro', 'P', 'T')]} (> 0 であること)")
        sh = (ni, nj)
        row.update(SR.recorded_quantities(Xf, Rf, {k: V[k].reshape(sh) for k in ("ro", "Ux", "Uy", "sonic")}, x_E))
        F = NR.load_field(run, name)
        for x in XS:
            row[f"M_axis_{x:g}"] = float(NR.eta_line(F, "M", 0.0, np.array([x]))[0])
            row[f"M_eta0.1_{x:g}"] = float(NR.eta_line(F, "M", 0.1, np.array([x]))[0])
        try:
            m = NR.metrics(run, F)["eta0.1"]
            row.update({"wave_eta0.1": m["wave_pct"], "overshoot_eta0.1": m["overshoot_pct"], "slope_eta0.1": m["slope_pct"],
                        "range_eta0.1": m["range_pct"]})
        except Exception as e:  # noqa: BLE001 — 記録だけの量
            row["metrics_error"] = repr(e)
        rows.append(row)
    return rows, probs


def sigma_ref_required(last=None, tail=None, every=None) -> list:
    """σ_ref の代わりの参照に要る Euler の出力の step (§5.1 #36i (M1)): 末尾の窓 [last − tail, last] の出力の間隔 every ごとの全出力から
    固定の last を除いたもの (V-h370: 40,000〜59,000 の 20 出力)。登録の定数 (既定 EU_STEPS・EU_TAIL・EU_OUT、呼ぶ時点の値) から列挙する
    (在る出力・予備の資格の記録からは作らない)。"""
    last, tail, every = (EU_STEPS if last is None else last), (EU_TAIL if tail is None else tail), (EU_OUT if every is None else every)
    lo = last - tail
    return [lo + k * every for k in range(tail // every) if lo + k * every < last]


def sigma_ref_selection(steps: list, last=None, tail=None, every=None) -> dict:
    """σ_ref の代わりの参照の選び方 (§5.1 #36i (M1); 予備の資格の時点 = NS の判定の前に固定する)。末尾の窓 [last − tail, last] の
    **全**出力 (出力の間隔 every ごと、固定の last を除く; V-h370 は 40,000〜59,000 の 20 出力) を代わりの参照にする。
    旧い代表 (SIGMA_REF_REP_EVERY ごと) は記録だけ。
    戻り {"tail_window", "tail_steps" (窓の在る出力), "required" (要る出力の列挙), "alternatives" (要る出力のうち在るもの), "missing"
    (要る出力のうち無いもの), "representatives_record" (旧い代表、記録だけ), "rule", "why"}。"""
    last, tail, every = (EU_STEPS if last is None else last), (EU_TAIL if tail is None else tail), (EU_OUT if every is None else every)
    lo = last - tail
    tail_steps = [s_ for s_ in steps if lo <= s_ <= last]
    want = sigma_ref_required(last, tail, every)
    reps = [lo + k * SIGMA_REF_REP_EVERY for k in range(tail // SIGMA_REF_REP_EVERY) if lo + k * SIGMA_REF_REP_EVERY < last]
    return {"tail_window": [lo, last], "tail_steps": tail_steps, "required": want, "alternatives": [s_ for s_ in want if s_ in tail_steps],
            "missing": [s_ for s_ in want if s_ not in tail_steps], "representatives_record": reps,
            "rule": f"末尾の窓 [{lo}, {last}] の全出力 ({every} step ごと、固定の {last} を除く {len(want)} 出力) を代わりの参照にする",
            "why": (f"代表 ({SIGMA_REF_REP_EVERY} step ごと) に間引くと代表の間の出力だけの参照の変動を見ない (codex diagnose "
                    "2026-10-14-h370-prelaunch-review の反例: NS 82,500 × Euler 42,000 だけ 0.030 mm 動く場を代表 4 出力は σ_ref 0 と見た)。"
                    f"NS の 17 出力 × {len(want)} = {17 * len(want)} 抽出は NS の判定でプロセス並列にして回す。旧い代表と場の差 tail_field_dev は"
                    "記録だけ (資格の判定に使わない)")}


def tail_field_dev(run: Path, steps: list, pin_step: int, reps: list) -> dict:
    """記録だけ (資格の判定には使わない、§5.1 #36i (M1)): 末尾の窓の各出力の ρU_x の固定の出力との差 max|ρU_x − ρU_x,pin| / max|ρU_x,pin|。
    旧い代表 (reps) の中の最大と代表の外の最大を並べる。閾値・抽出の誤差の上界を持たないので σ_ref の保証にはならない (σ_ref は全出力で
    抽出し直して測る)。"""
    import h5py

    def q(st):
        with h5py.File(run / f"res_{st}.h5", "r") as h:
            return np.asarray(h["/VALUE/ro"][:], dtype=float) * np.asarray(h["/VALUE/Ux"][:], dtype=float)
    q0 = q(pin_step)
    den = float(np.max(np.abs(q0))) or 1.0
    per = [{"step": st, "max_rel": float(np.max(np.abs(q(st) - q0)) / den)} for st in steps if st != pin_step]
    rep = [p_["max_rel"] for p_ in per if p_["step"] in reps]
    oth = [p_ for p_ in per if p_["step"] not in reps]
    o = max(oth, key=lambda p_: p_["max_rel"]) if oth else None
    return {"per_output": per, "max_rel_representatives": max(rep) if rep else None,
            "max_rel_others": o["max_rel"] if o else None, "step_max_others": o["step"] if o else None,
            "note": "記録だけ (合否に使わない; 閾値も抽出の誤差の上界も無い)。σ_ref は末尾の窓の全出力で抽出し直して測る (§5.1 #36i M1)"}


def euler_prequalify(case: Path, tools: Path) -> int:
    """(i) Euler の参照の予備の資格と固定 (上の docstring)。終了コード 0 予備の資格あり / 1 未達 / 2 判定不能。"""
    run = case / EU_RUN
    R = {"plan": PLAN, "created": now(), "run": EU_RUN, "tools": str(tools), "undecidable": [], "failures": [], "prequalified": False,
         "kind": "予備の資格 (prequalify)。最終の資格は NS の判定の σ_ref (h370_judge.py の既定の経路の reference_final) が決め、"
                 "update_conditions_met に繋ぐ (§5.1 #36g (1))",
         "definitions": {"qualifying_cols": QUAL_COLS, "record_cols": REC_COLS, "positions_rt": list(XS),
                         "quasisteady": "check_quasisteady.py --series-csv (既定の --tail 0.4・--drift 0.05・--osc 0.10・--min-snaps 4)、全列 STEADY",
                         "exit_core_M・mdot_rt": "forge_design.evaluate.stop_rule.recorded_quantities (nozzle_report.metrics・qs_series.py と同じ定義)",
                         "M_axis・M_eta0.1": "forge_design.report.nozzle_report.eta_line(F, 'M', η, x) (η = 0 は軸、0.1 は metrics の判定の線)",
                         "record_cols": "nozzle_report.metrics の eta0.1 (wave_pct・overshoot_pct・slope_pct・range_pct)、記録だけ",
                         "pin": "本段の最後の出力 (上限 60,000 step) とその sha256",
                         "sigma_ref_snapshots": (f"σ_ref の代わりの参照: 末尾の窓 (最後の {EU_TAIL} step) の全出力 ({EU_OUT} step ごと、固定を除く "
                                                 f"{len(sigma_ref_required())} 出力; §5.1 #36i M1)"),
                         "sigma_ref_representatives_record": f"旧い代表 ({SIGMA_REF_REP_EVERY} step ごと) と tail_field_dev は記録だけ (資格の判定に使わない)",
                         "screening_snapshots": "スクリーニング (check_quasisteady の末尾の窓の始まりと中ほど; NS の停止の出力だけに当てる記録、σ_ref の上限とは扱わない)"}}
    und, fails = R["undecidable"], R["failures"]
    qual_p, pin_p = case / QUAL, case / PIN

    def finish(rc):
        R["rc"] = rc
        jwrite(qual_p, R)
        print(f"== Euler の予備の資格 ({EU_RUN}): {'予備の資格あり → ' + PIN if rc == 0 else ('未達' if rc == 1 else '判定不能')}"
              + (f" — {'; '.join(und + fails)[:1500]}" if und or fails else ""))
        print(f"→ {qual_p}")
        return rc
    if not run.is_dir() or not (run / "prepare_info.json").is_file() or not (run / "nozzle.h5").is_file():
        und.append(f"{EU_RUN} が無いか prepare_info.json・nozzle.h5 が無い")
        return finish(2)
    rc_txt = (run / "RUN_RC").read_text().strip() if (run / "RUN_RC").is_file() else None
    R["run_rc"] = rc_txt
    if rc_txt != "0":
        fails.append(f"RUN_RC {rc_txt!r} (0 であること)")
    pc = jload(run / "H370_PRECHECK.json") or {}
    R["precheck_ok"] = pc.get("ok")
    if not pc.get("ok"):
        fails.append("投入前の検査 (H370_PRECHECK.json) が無いか不成立")
    for m in MARKERS:
        if (run / m).is_file():
            fails.append(f"基盤の停止・不成立の印 {m}: {(run / m).read_text()[:200]}")
    steps = [s for s in res_steps(run) if s > 0]
    R["outputs"] = {"n": len(steps), "first": steps[0] if steps else None, "last": steps[-1] if steps else None}
    if not steps:
        und.append("本段の出力 res_<N>.h5 (N > 0) が無い")
        return finish(2)
    if steps[-1] != EU_STEPS:
        fails.append(f"本段の最後の出力 res_{steps[-1]}.h5 が上限 {EU_STEPS} step でない")
    R["residual_nonfinite"] = {t: residual_nonfinite(run / f"residual_history_{t}.csv") for t in ("soft", "mid", "main")}
    for t, v in R["residual_nonfinite"].items():
        if v is None:
            fails.append(f"residual_history_{t}.csv が無い")
        elif v["nonfinite"]:
            fails.append(f"residual_history_{t}.csv の rms_* に非有限 {v['nonfinite']} 個 (最初 step {v['first_nonfinite_step']})")
    # 全残差 (同じ設定の区間)
    cmd = [sys.executable, str(tools / "check_convergence.py"), "--segment", "--drop", "3", str(run)]
    cv = subprocess.run(cmd, capture_output=True, text=True)
    (run / "H370_EULER_CONVERGENCE.txt").write_text(cv.stdout + cv.stderr)
    R["convergence"] = {"cmd": " ".join(cmd[1:]), "rc": cv.returncode,
                        "verdict": [ln.strip() for ln in cv.stdout.splitlines() if ln.startswith("===") or ln.startswith("OVERALL")]}
    if cv.returncode != 0:
        fails.append(f"check_convergence.py --segment --drop 3 が PASS でない (rc {cv.returncode}: {R['convergence']['verdict'][-1:] })")
    # 時系列
    try:
        rows, probs = euler_series(run, steps)
    except Exception as e:  # noqa: BLE001
        und.append(f"出力の量を作れない ({type(e).__name__}: {e})")
        return finish(2)
    fails += probs
    sc = run / "h370_euler_series.csv"
    write_series(sc, [r for r in rows if not r.get("nonfinite_all")], QUAL_COLS + REC_COLS)
    q = run_qs(tools, sc, QUAL_COLS)
    R["quasisteady"] = q
    R["quasisteady_record"] = {"record_cols (既定の閾値)": run_qs(tools, sc, REC_COLS), "tight (1e-4、記録だけ)": run_qs(tools, sc, QUAL_COLS + REC_COLS, QS_TIGHT)}
    if q["missing"]:
        und.append(f"check_quasisteady.py の出力に列 {q['missing']} の VERDICT が無い ({q['stderr'][-200:]})")
    bad = {c: v["verdict"] for c, v in q["verdicts"].items() if v["verdict"] != "STEADY"}
    if bad:
        fails.append(f"check_quasisteady.py の VERDICT が STEADY でない: {bad}")
    last = rows[-1]
    R["last_values"] = {c: last.get(c) for c in QUAL_COLS + REC_COLS}
    # σ_ref の代わりの参照 (末尾の窓の全出力; NS の判定の前のここで固定する) と末尾の窓の全出力の sha (§5.1 #36i M1)
    sel = sigma_ref_selection(steps)
    R["sigma_ref_selection"] = {k: sel[k] for k in ("tail_window", "required", "alternatives", "missing", "representatives_record", "rule", "why")}
    R["tail_outputs"] = [{"res": f"res_{s}.h5", "step": s, "sha256": sha256_file(run / f"res_{s}.h5")} for s in sel["tail_steps"]]
    R["sigma_ref_snapshots"] = [o for o in R["tail_outputs"] if o["step"] in sel["alternatives"]]
    R["sigma_ref_representatives_record"] = [o for o in R["tail_outputs"] if o["step"] in sel["representatives_record"]]
    if sel["missing"]:
        # 本段が上限まで回っていれば出力の欠け (消された) で判定不能、上限より前に終わった run は既に未達 (その帰結)
        (und if steps[-1] == EU_STEPS else fails).append(f"σ_ref の代わりの参照 (末尾の窓の全出力) のうち出力 {sel['missing']} が無い "
                                                         "(σ_ref を測れない)")
    elif steps[-1] == EU_STEPS:
        try:
            R["tail_field_dev"] = tail_field_dev(run, sel["tail_steps"], EU_STEPS, sel["representatives_record"])
        except Exception as e:  # noqa: BLE001 — 記録だけ
            R["tail_field_dev"] = {"error": f"{type(e).__name__}: {e}"[:300]}
    # スクリーニング (check_quasisteady の末尾の窓 = 既定 tail 0.4 の始まりと中ほど; 停止の出力だけに当てる記録)
    n = len(steps)
    k = max(3, int(math.ceil(0.4 * n)))
    tail = steps[-k:]
    alt = sorted({tail[0], tail[len(tail) // 2]} - {steps[-1]})
    R["screening_snapshots"] = [{"res": f"res_{s}.h5", "step": s, "sha256": sha256_file(run / f"res_{s}.h5")} for s in alt]
    R["quasisteady_tail_window"] = [tail[0], tail[-1]]
    if und:
        return finish(2)
    if fails:
        return finish(1)
    pin = {"plan": PLAN, "created": now(), "run": EU_RUN, "run_abs": str(run.resolve()), "res": f"res_{steps[-1]}.h5", "step": steps[-1],
           "sha256": sha256_file(run / f"res_{steps[-1]}.h5"), "prequalify": QUAL,
           "note": "NS の euler_ref は {run: <case dir>/run_0950_h370_euler_ref, res, sha256} (h370_run.py read_pin)。予備の資格で固定し、"
                   "最終の資格は NS の判定の σ_ref"}
    R["pin"] = {k: pin[k] for k in ("run", "res", "step", "sha256")}
    old = jload(pin_p)
    if old is not None:
        if (old.get("run"), old.get("res"), old.get("sha256")) != (pin["run"], pin["res"], pin["sha256"]):
            und.append(f"{PIN} が既にあり、今回の固定 ({pin['res']} {pin['sha256'][:16]}) と違う ({old.get('res')} {str(old.get('sha256'))[:16]}) — 消さない")
            return finish(2)
        R["pin_note"] = "既にある同じ固定をそのまま使う"
    else:
        pin_p.parent.mkdir(parents=True, exist_ok=True)
        fd = os.open(str(pin_p), os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        os.write(fd, json.dumps(pin, indent=1, ensure_ascii=False).encode())
        os.close(fd)
    R["prequalified"] = True
    return finish(0)


# ===================================================================================================== (ii) NS
@contextlib.contextmanager
def snapshot_view(run: Path, step: int):
    """run の prepare_info.json・nozzle.h5・res_<step>.h5 だけを symlink した一時ディレクトリ (抽出器が最新の出力を読むので 1 出力に絞る)。"""
    with tempfile.TemporaryDirectory(prefix="h370_view_") as td:
        v = Path(td) / run.name
        v.mkdir()
        for name in ("prepare_info.json", "nozzle.h5", f"res_{step}.h5"):
            (v / name).symlink_to((run / name).resolve())
        yield v


def _arr(d, k):
    return np.asarray(d.get(k, []), dtype=float)


# A/B の分類の理由の符号 (classify_ab の code)。参照の最終の資格 (reference_final) は SIGMA_FAIL なら資格なし、REF_OK (結果 A・結果 B に
# 分類できたもの) なら σ_ref,use ≤ 0.01 mm を併せて資格あり、それ以外 (判定不能: U_time を含む) は確かめられない (更新に進めない)。
# U_time (窓の中の時間の変動で D が 0.01 mm をまたぐ) は §5.1 #36i (M2) で REF_OK から外した (判定不能は更新不可に統一)
SIGMA_FAIL = ("U_sigma_unmeasurable", "U_sigma_undecidable", "U_sigma_straddle")
REF_OK = ("A", "B_D", "B_pf")


def fixed_columns(x_use, eu_runs: dict, edge_tol: bool) -> dict:
    """比べる列の集合 F を抽出の前に座標だけで決める (§5.1 #36g (3)): NS の使用域の列 x_use のうち、各腕の Euler の補間できる範囲
    (`deltastar.euler_cover_range`; edge_tol は受け渡しの抽出の remap_edge == "missing" と同じ) の全部に入る列。範囲は各 run の格子
    (prepare_info.json・nozzle.h5 の列の x) だけで決まり、場・抽出の結果を見ない。戻り {"mask" (bool の配列), "n", "x_range", "cover",
    "excluded" (腕ごとに範囲の外で外した列の数と x の範囲)}。"""
    from forge_design.metrics.delta_in_check import wall_columns_from_run
    from forge_design.metrics.deltastar import euler_cover_range, euler_covers
    xu = np.asarray(x_use, dtype=float)
    F = np.ones(len(xu), dtype=bool)
    out = {"cover": {}, "excluded": {}, "edge_tol": bool(edge_tol)}
    for arm, rd in eu_runs.items():
        rng = euler_cover_range(wall_columns_from_run(rd)["x"], edge_tol=edge_tol)
        c = euler_covers(xu, rng)
        out["cover"][arm] = rng
        out["excluded"][arm] = {"n": int((~c).sum()), "x_range": [float(xu[~c].min()), float(xu[~c].max())] if (~c).any() else None}
        F &= c
    out.update(mask=F, n=int(F.sum()), n_use=int(len(xu)), x_range=[float(xu[F].min()), float(xu[F].max())] if F.any() else None)
    return out


def arm_check(d, F) -> dict:
    """1 腕 1 出力の抽出 d (`stop_rule.delta_snapshot` の戻り) の、列の集合 F の中の合否 pf: 抽出の例外なし・F の中に held なし・δ_E・δ_target が
    有限 (抽出しなかった列 = 欠測は NaN なのでここで落ちる)。例外・形の違いは {"error"} (合否ではなく判定不能の理由)。"""
    if not isinstance(d, dict) or "error" in d or "delta_E" not in d:
        return {"error": (d or {}).get("error", "δ_E が無い") if isinstance(d, dict) else "記録が無い"}
    e, t = _arr(d, "delta_E"), _arr(d, "delta_target")
    h = np.asarray(d.get("held", []), dtype=bool)
    if not (len(e) == len(t) == len(h) == len(F)):
        return {"error": f"列の数 {len(e)}・{len(t)}・{len(h)} が比べる列の表 {len(F)} と違う"}
    miss = np.zeros(len(F), dtype=bool)
    miss[[j for j in (d.get("missing") or []) if 0 <= int(j) < len(F)]] = True
    nonfin = F & ~(np.isfinite(e) & np.isfinite(t))
    held = F & h
    return {"pf": bool(not nonfin.any() and not held.any()), "n_nonfinite": int(nonfin.sum()), "n_held": int(held.sum()),
            "n_missing": int((F & miss).sum()), "i_bad": [int(j) for j in np.flatnonzero(nonfin | held)[:8]]}


def ab_pair(step: int, dA: dict, dB: dict, mm: float, src_b: str, x_use, F) -> dict:
    """出力 1 つの腕 A・B の比較 (§5.1 #36g (3)): 比べる列は先に決めた F (抽出の後に有限の列だけを選ばない)。各腕の合否 pf_A・pf_B は F の中で
    見る。両腕とも pf のときだけ D_mm = F の max|δ_E,B − δ_E,A| [mm] (Dt_mm は δ_target)。どちらかの抽出の例外は error_<腕>。"""
    r = {"step": step, "B_source": src_b, "n_fixed": int(np.sum(F))}
    xu = np.asarray(x_use, dtype=float)
    for arm, d in (("A", dA), ("B", dB)):
        c = arm_check(d, F)
        if "error" in c:
            r[f"error_{arm}"] = str(c["error"])[:300]
        else:
            r[f"pf_{arm}"] = c["pf"]
            r[f"check_{arm}"] = {k: c[k] for k in ("n_nonfinite", "n_held", "n_missing")}
            if c["i_bad"]:
                r[f"x_bad_{arm}"] = [float(xu[j]) for j in c["i_bad"]]
    if any(k.startswith("error") for k in r) or not np.any(F):
        if not np.any(F):
            r["error_fixed"] = "比べる列 F が空"
        return r
    if r["pf_A"] and r["pf_B"]:
        eA, eB, tA, tB = _arr(dA, "delta_E"), _arr(dB, "delta_E"), _arr(dA, "delta_target"), _arr(dB, "delta_target")
        d = np.abs(eB - eA)[F] * mm
        i = int(np.argmax(d))
        r["D_mm"], r["i_at_D"] = float(d[i]), int(np.flatnonzero(F)[i])
        r["Dt_mm"] = float(np.max(np.abs(tB - tA)[F]) * mm)
    else:
        r["D_mm"] = None
    return r


def _sigma_norm(sigma) -> dict:
    """classify_ab の sigma を {"status", "mm", "per_step", "why"} にそろえる (float = 測れた・per_step なし、None = 測れない)。"""
    if isinstance(sigma, dict):
        return {"status": sigma.get("status"), "mm": sigma.get("mm"), "per_step": sigma.get("per_step"), "why": sigma.get("why")}
    if sigma is None:
        return {"status": "unmeasurable", "mm": None, "per_step": None, "why": "σ_ref が無い"}
    return {"status": "ok", "mm": float(sigma), "per_step": None, "why": None}


def classify_ab(per: list, sigma, required, tol: float = TOL_MM) -> dict:
    """A/B の結果 (§6 V-h370・§5.1 #36g (4)(5)): 戻り {"result" (結果 A / 結果 B / 判定不能), "why", "code"}。
    per: 出力ごとの `ab_pair` の戻り。required: 判定窓に入るべき出力の step (`stop_rule.required_window_steps`; 無ければ判定不能)。
    sigma: `sigma_ref_measure` の戻り (status ok / unmeasurable / undecidable、mm、per_step) か、試験用に float (測れた・per_step なし)・None。
    順に: 必要な出力の欠け → 判定不能、抽出の失敗 → 判定不能、σ_ref を測れない・判定不能 → 判定不能、どれかの出力で pf_A ≠ pf_B →
    代わりの参照の全部 (per_step の alts の pf_F) でも違いが保たれれば 結果 B・確かめられなければ判定不能、両腕とも pf でない出力
    (F の中の欠測・held・非有限) → 判定不能、それ以外は D: max D + σ ≤ tol → 結果 A、min D − σ > tol → 結果 B、σ を足し引きして初めて
    tol をまたぐ → 判定不能 (参照の未整定)、D だけで窓の中で tol をまたぐ → 判定不能 (時間の変動)。"""
    sg = _sigma_norm(sigma)

    def res(result, code, why):
        return {"result": result, "code": code, "why": why}
    if not required:
        return res("判定不能", "U_required", "判定窓に入るべき出力の集合が無い (出力の間隔・窓から列挙できない)")
    if not per:
        return res("判定不能", "U_empty", "比べる出力が無い")
    have = {r["step"] for r in per if "error_missing" not in r}
    lack = [k for k in required if k not in have]
    if lack:
        return res("判定不能", "U_incomplete", f"判定窓に入るべき {len(required)} 出力のうち比べていない出力 (記録も res も無く再抽出で"
                                              f"補えない・渡されていない) {len(lack)}: {lack[:8]}")
    per = [r for r in per if r["step"] in set(required)]
    err = [r["step"] for r in per if any(k.startswith("error") for k in r)]
    if err:
        return res("判定不能", "U_error", f"抽出・比較の失敗 (出力 {err[:8]})")
    if sg["status"] == "unmeasurable" or (sg["status"] == "ok" and sg["mm"] is None):
        return res("判定不能", "U_sigma_unmeasurable", f"参照の時間の変動 σ_ref を測れない ({sg.get('why')})")
    if sg["status"] != "ok":
        return res("判定不能", "U_sigma_undecidable", f"参照の時間の変動 σ_ref が判定不能 — 代わりの参照で held・欠測・合否が変わる ({sg.get('why')})")
    sig = float(sg["mm"])
    pfd = [r for r in per if r["pf_A"] != r["pf_B"]]
    if pfd:
        ps = sg.get("per_step") or {}
        unconf = []
        for r in pfd:
            alts = (ps.get(r["step"]) or ps.get(str(r["step"])) or {}).get("alts") or []
            if not alts or any(a.get("pf_F") is None or a["pf_F"] == r["pf_A"] for a in alts):
                unconf.append(r["step"])
        detail = "; ".join(f"{r['step']}: A {'合' if r['pf_A'] else '否'}・B {'合' if r['pf_B'] else '否'}" for r in pfd[:6])
        if unconf:
            return res("判定不能", "U_pf_unconfirmed", f"抽出の合否が変わる出力 ({detail}) の違いを、代わりの参照の時刻で確かめられない・保たれない "
                                                     f"(出力 {unconf[:6]})")
        return res("結果 B", "B_pf", f"抽出の合否が変わる出力 {len(pfd)} ({detail}); 代わりの参照の全部の時刻でも違いが保たれる "
                                    "— 旧い参照の流用を退ける")
    inc = [r["step"] for r in per if r.get("D_mm") is None]
    if inc:
        return res("判定不能", "U_fixed_cols", f"両腕とも比べる列 F の中に欠測・held・非有限がある出力 {inc[:8]} — 比較を合格にしない")
    D = [r["D_mm"] for r in per]
    dmax, dmin = max(D), min(D)
    why = f"max_k D_k {dmax:.4g} mm・min_k D_k {dmin:.4g} mm・σ_ref {sig:.4g} mm (許容 {tol} mm、出力 {len(per)})"
    if dmax + sig <= tol:
        return res("結果 A", "A", f"比べる列 F の全断面・判定窓の全出力で max|δ_E,B − δ_E,A| + σ_ref ≤ {tol} mm かつ合否が同じ ({why}) — "
                                "「旧い参照が共通区間で予算を超えるずれを生む」を退ける (旧い参照の出口の不足は残る)")
    if dmin - sig > tol:
        return res("結果 B", "B_D", f"全出力で max|δ_E,B − δ_E,A| − σ_ref > {tol} mm ({why}) — 旧い参照の流用を退ける")
    if dmax <= tol or dmin > tol:
        return res("判定不能", "U_sigma_straddle", f"参照の未整定: σ_ref を足し引きして初めて {tol} mm をまたぐ ({why})")
    return res("判定不能", "U_time", f"窓の中の時間の変動で D_k が {tol} mm をまたぐ ({why})")


def reference_final(ab, do_ab: bool, tol: float = TOL_MM) -> dict:
    """参照の最終の資格 (§5.1 #36g (1)・#36i (M2)・(M3)): 次のすべてのとき資格あり — σ_ref・σ_ref,use を測れて判定できる (status ok)、
    **σ_ref,use ≤ tol** (使用域の全体の参照の変動; A/B の分類とは独立の必要条件)、A/B が結果 A・結果 B に分類された (REF_OK: σ_ref が
    分類を左右しない)。σ を測れない・判定不能・σ_ref,use > tol・参照の未整定 (SIGMA_FAIL) → 資格なし (False)。A/B が判定不能
    (U_time を含む) → 確かめられない (None)。False・None はどちらも更新に進めない。
    戻り {"qualified", "why", "sigma_ref_mm", "sigma_ref_at", "sigma_ref_use_mm", "sigma_ref_use_at", "sigma_status", "ab_code", "tol_mm"}。"""
    if not do_ab:
        return {"qualified": None, "why": "A/B と σ_ref を回していない (--no-ab) — 参照の最終の資格を確かめていない"}
    if not ab or "classify" not in ab:
        return {"qualified": None, "why": (ab or {}).get("why") or "A/B を回せない — 参照の最終の資格を確かめていない"}
    c, sg = ab["classify"], ab.get("sigma_ref") or {}
    use = sg.get("mm_use")
    out = {"sigma_ref_mm": sg.get("mm"), "sigma_ref_at": sg.get("at_F"), "sigma_ref_use_mm": use, "sigma_ref_use_at": sg.get("at_use"),
           "sigma_status": sg.get("status"), "ab_code": c["code"], "tol_mm": tol}
    if sg.get("status") != "ok" or use is None or sg.get("mm") is None:
        out.update(qualified=False, why=f"資格なし: σ_ref・σ_ref,use を測れない・判定不能 (status {sg.get('status')}: {sg.get('why')})"[:600])
    elif use > tol:
        at = sg.get("at_use") or {}
        out.update(qualified=False, why=(f"資格なし: 使用域の全体の参照の変動 σ_ref,use {use:.4g} mm > {tol} mm (NS の出力 {at.get('step')}・Euler の出力 "
                                         f"{at.get('euler_res')}・x {at.get('x')}) — 新しい参照が設計の入力として落ち着いていない "
                                         f"(A/B の分類 {c['result']} に関係なく、§5.1 #36i M3)"))
    elif c["code"] in SIGMA_FAIL:
        out.update(qualified=False, why=f"資格なし: {c['why']}")
    elif c["code"] in REF_OK:
        out.update(qualified=True, why=(f"資格あり: σ_ref,use {use:.4g} mm ≤ {tol} mm、σ_ref {sg.get('mm'):.4g} mm は A/B の分類 "
                                        f"({c['result']}) を左右しない"))
    else:
        out.update(qualified=None, why=f"確かめられない: A/B が判定不能 ({c['code']}: {c['why']}) — 判定不能は更新に進めない (§5.1 #36i M2)")
    return out


def sigma_ref_alternatives(spec_b: dict, qual) -> tuple:
    """σ_ref の代わりの参照の照合 (§5.1 #36i (M1); 抽出の前): 固定の出力が本段の上限の出力 (res_<EU_STEPS>.h5) で、予備の資格の記録の
    sigma_ref_snapshots が末尾の窓の全出力 (`sigma_ref_required`、20 出力) ちょうどで、各出力が腕 B の run に在り sha が記録と同じ。
    戻り (alts (step の順), check {"ok", "required", "recorded", "missing_in_record", "extra_in_record", "missing_files", "sha_mismatch", "why"})。
    ok でなければ σ_ref を測れない (抽出は回さない)。"""
    req = sigma_ref_required()
    raw = [a for a in ((qual or {}).get("sigma_ref_snapshots") or []) if isinstance(a, dict)]
    rec = {}
    for a in raw:
        try:
            rec[int(a.get("step"))] = a
        except (TypeError, ValueError):
            pass
    chk = {"required": req, "recorded": sorted(rec), "missing_in_record": [k for k in req if k not in rec],
           "extra_in_record": sorted(k for k in rec if k not in req), "missing_files": [], "sha_mismatch": [], "ok": False}
    why = []
    if spec_b.get("res") != f"res_{EU_STEPS}.h5":
        why.append(f"固定の出力 {spec_b.get('res')} が本段の上限の出力 res_{EU_STEPS}.h5 でない")
    if chk["missing_in_record"] or chk["extra_in_record"]:
        why.append(f"予備の資格の記録の代わりの参照が末尾の窓の全 {len(req)} 出力でない (無い {chk['missing_in_record'][:8]}・余分 "
                   f"{chk['extra_in_record'][:8]})")
    rd = Path(spec_b["run"])
    alts = []
    for k in req:
        a = rec.get(k)
        if a is None:
            continue
        if a.get("res") != f"res_{k}.h5":
            chk["missing_files"].append(k)
            continue
        f = rd / a["res"]
        if not f.is_file():
            chk["missing_files"].append(k)
        elif sha256_file(f) != a.get("sha256"):
            chk["sha_mismatch"].append(k)
        else:
            alts.append({"res": a["res"], "step": k, "sha256": a["sha256"]})
    if chk["missing_files"]:
        why.append(f"代わりの参照の出力が腕 B の run に無い {chk['missing_files'][:8]}")
    if chk["sha_mismatch"]:
        why.append(f"代わりの参照の出力の sha が予備の資格の記録と違う {chk['sha_mismatch'][:8]}")
    chk["ok"] = not why
    if why:
        chk["why"] = "; ".join(why)
    return alts, chk


def mem_available_mib():
    """/proc/meminfo の MemAvailable [MiB] (読めなければ None)。"""
    try:
        for ln in Path("/proc/meminfo").read_text().splitlines():
            if ln.startswith("MemAvailable:"):
                return int(ln.split()[1]) / 1024.0
    except (OSError, ValueError, IndexError):
        pass
    return None


def snap_jobs(requested=None) -> dict:
    """抽出の並列数 (§5.1 #36i 計算量): min(上限 (既定 JOBS_CAP = 3、--jobs で 1〜JOBS_MAX), 空きメモリ / SNAP_PEAK_MIB)。空きメモリが
    読めなければ 1。戻り {"jobs", "cap", "mem_available_mib", "per_extraction_peak_mib", "by_memory", "cpu_count"}。"""
    cap = JOBS_CAP if requested is None else int(requested)
    if not 1 <= cap <= JOBS_MAX:
        raise ValueError(f"--jobs は 1〜{JOBS_MAX} ({requested!r})")
    avail = mem_available_mib()
    by_mem = int(avail // SNAP_PEAK_MIB) if avail is not None else 1
    return {"jobs": max(1, min(cap, by_mem)), "cap": cap, "mem_available_mib": None if avail is None else round(avail),
            "per_extraction_peak_mib": SNAP_PEAK_MIB, "by_memory": by_mem, "cpu_count": os.cpu_count()}


def _call_snap(snap, task) -> dict:
    """抽出 1 回 (task = (NS の run, step, Euler の参照の指定, 抽出の条件, 使用域の列 x[, 足す引数 …]); 足す引数は帯を共通にした A/B
    (§5.1 #36n、`diag_snapshot_fixed`) の固定の帯だけ)。戻り {"ok": 抽出} か {"exc": 例外の文}。
    直列でも並列でも同じこの関数を通す (例外の扱いを揃える)。"""
    run, k, spec, ext, x_use, *more = task
    try:
        return {"ok": snap(run, k, spec, ext, x_use, *more)}
    except Exception as e:  # noqa: BLE001 — sha の不一致・出力が無いを含む
        return {"exc": f"{type(e).__name__}: {e}"}


def _snap_worker(task) -> dict:
    """プロセス並列の子 (spawn) の抽出: 本物の stop_rule.delta_snapshot (sys.path は親から引き継ぐ)。"""
    from forge_design.evaluate import stop_rule as SR
    return _call_snap(SR.delta_snapshot, task)


def _diag_snap_worker(task) -> dict:
    """プロセス並列の子 (spawn) の抽出: 診断の感度の判定の `diag_snapshot` (§6 V-h370-diag)。"""
    return _call_snap(diag_snapshot, task)


def _diag_fixed_snap_worker(task) -> dict:
    """プロセス並列の子 (spawn) の抽出: 帯を共通にした A/B の `diag_snapshot_fixed` (§5.1 #36n)。帯の固定はこの子の中だけで効く。"""
    return _call_snap(diag_snapshot_fixed, task)


def _process_worker(snap):
    """snap が本物の抽出なら子プロセスで回す関数 (stop_rule.delta_snapshot → _snap_worker、diag_snapshot → _diag_snap_worker、
    diag_snapshot_fixed → _diag_fixed_snap_worker)、それ以外 (試験の差し替え) は None (このプロセスで直列)。"""
    from forge_design.evaluate import stop_rule as SR
    if snap is SR.delta_snapshot:
        return _snap_worker
    if snap is diag_snapshot:
        return _diag_snap_worker
    if snap is diag_snapshot_fixed:
        return _diag_fixed_snap_worker
    return None


@contextlib.contextmanager
def _single_thread_env():
    """子プロセスの BLAS・OpenMP を 1 スレッドにする (並列の子どうしでコアを取り合わない)。spawn の子は起動の時点の環境変数を読む。"""
    keys = ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")
    old = {k: os.environ.get(k) for k in keys}
    os.environ.update({k: "1" for k in keys})
    try:
        yield
    finally:
        for k, v in old.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v


def run_snaps(tasks: list, snap, jobs: int = 1, progress=None) -> tuple:
    """抽出の作業をまとめて回す (§5.1 #36i 計算量)。snap が本物の stop_rule.delta_snapshot なら **jobs によらず**子プロセス (spawn、
    max_workers = jobs、子は BLAS・OpenMP 1 スレッド) で回す: jobs 1 (直列) も子 1 つで回す。BLAS のスレッド数で行列積・解の丸めが変わり
    (実寸の合成の場で親プロセス (BLAS 既定のスレッド数) と子 (1 スレッド) の抽出が 2,938 列のうち 771 列で違った: δ_E 最大 6.9e−18 r_t、
    δ_target 最大 2.8e−17 r_t、held は同じ; 2026-10-14 の手元の実測)、
    直列と並列の値をビット同一にするには計算する環境をそろえる必要があるため。試験の差し替え (snap が別の関数) は子に渡せないので
    このプロセスで直列。結果は作業の順に並べる (抽出は互いに独立なので、並列でも記録の順・値は直列と同じ)。
    戻り (結果 [{"ok": 抽出} | {"exc": 文}], 記録 {"n", "mode", "jobs", "seconds", "children_maxrss_mib"})。progress(i, n, 経過秒) は進み。
    診断の感度の判定の `diag_snapshot` (§6 V-h370-diag) も本物の抽出なので同じく子プロセスで回す (`_process_worker`)。"""
    import resource
    import time
    n = len(tasks)
    info = {"n": n, "mode": "in-process (serial)", "jobs": 1}
    t0 = time.time()
    out = [None] * n
    worker = _process_worker(snap)
    if n and worker is not None:
        import multiprocessing as mp
        from concurrent.futures import ProcessPoolExecutor, as_completed
        w = max(1, min(int(jobs), n))
        info.update(mode="process (spawn)", jobs=w, blas_threads=1)
        with _single_thread_env(), ProcessPoolExecutor(max_workers=w, mp_context=mp.get_context("spawn")) as ex:
            futs = {ex.submit(worker, t): i for i, t in enumerate(tasks)}
            for m, f in enumerate(as_completed(futs), 1):
                try:
                    out[futs[f]] = f.result()
                except Exception as e:  # noqa: BLE001 — 子の異常終了 (BrokenProcessPool など) は抽出の失敗として残す
                    out[futs[f]] = {"exc": f"{type(e).__name__}: {e}"}
                if progress:
                    progress(m, n, time.time() - t0)
    else:
        if worker is None:
            info["note"] = "抽出が差し替え (試験) なので子に渡さずこのプロセスで直列"
        for i, t in enumerate(tasks):
            out[i] = _call_snap(snap, t)
            if progress:
                progress(i + 1, n, time.time() - t0)
    info["seconds"] = round(time.time() - t0, 1)
    info["children_maxrss_mib"] = round(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024.0, 1)
    return out, info


def _masks(d, n: int):
    """抽出 d の使用域の列ごとの印: held と、δ_E・δ_target が有限 (欠測の列は NaN なので偽)。"""
    e, t = _arr(d, "delta_E"), _arr(d, "delta_target")
    h = np.asarray(d.get("held", []), dtype=bool)
    if not (len(e) == len(t) == len(h) == n):
        return None, None
    return h, np.isfinite(e) & np.isfinite(t)


def sigma_ref_eval(steps: list, alts: list, d0_by: dict, dalt: dict, x_use, F, mm: float) -> dict:
    """σ_ref・σ_ref,use (§5.1 #36g (2)・#36i (M1)・(M3)) を抽出の結果から測る。d0_by[k] = 固定の参照の抽出 (腕 B)、dalt[(k, res)] =
    代わりの参照の抽出の結果 ({"ok": 抽出} か {"exc": 文})。
    σ_ref = max_k max_e max_F |Δδ_E| [mm] (比べる列 F)、σ_ref,use = max_k max_e max_U |Δδ_E| [mm] (U = 固定の参照で δ_E・δ_target が
    有限な使用域の列 = 使用域の全体; F ⊂ U)。最大を与えた (NS の出力・Euler の出力・x) を at_F・at_use に残す。
    status: ok / unmeasurable (代わりの参照が無い・抽出の例外 = sha の不一致を含む・固定の抽出が無い) / undecidable (代わりの参照で
    F・使用域の合否、または使用域の列ごとの held・有限の印が変わる)。ok でなければ mm・mm_use は None (観測値は *_observed)。"""
    xu = np.asarray(x_use, dtype=float)
    out = {"status": "ok", "mm": None, "mm_use": None, "at_F": None, "at_use": None, "alts": [a.get("res") for a in alts], "per_step": {},
           "why": None, "seconds": 0.0}
    if not alts:
        return {**out, "status": "unmeasurable", "why": "代わりの参照 (予備の資格の sigma_ref_snapshots) が無い"}
    FU = np.ones(len(F), dtype=bool)
    worst, worst_u, unmeas, undec = 0.0, 0.0, [], []
    for k in steps:
        d0 = d0_by.get(k)
        c0F, c0U = arm_check(d0, F), arm_check(d0, FU)
        if "error" in c0F:
            unmeas.append(f"{k}: 固定の参照の抽出が無い ({c0F['error'][:80]})")
            continue
        e0 = _arr(d0, "delta_E")
        h0, f0 = _masks(d0, len(F))
        rec = {"pf_F_pin": c0F["pf"], "pf_use_pin": c0U["pf"], "n_use_extractable": int(f0.sum()), "alts": []}
        for a in alts:
            r = dalt.get((k, a.get("res")))
            if r is None or "exc" in r:
                msg = "抽出していない" if r is None else r["exc"]
                unmeas.append(f"{k} × {a.get('res')}: {msg}"[:200])
                rec["alts"].append({"res": a.get("res"), "error": msg[:200]})
                continue
            d = r["ok"]
            out["seconds"] += float(d.get("seconds") or 0.0)
            cF, cU = arm_check(d, F), arm_check(d, FU)
            if "error" in cF:
                unmeas.append(f"{k} × {a.get('res')}: {cF['error'][:80]}")
                rec["alts"].append({"res": a.get("res"), "error": cF["error"][:200]})
                continue
            e1 = _arr(d, "delta_E")
            h1, f1 = _masks(d, len(F))
            both = f0 & f1
            dd = np.abs(e1 - e0) * mm
            ar = {"res": a.get("res"), "pf_F": cF["pf"], "pf_use": cU["pf"], "max_mm_F": None, "x_at_F": None, "max_mm_use": None,
                  "x_at_use": None, "change": {q: [c0U[q], cU[q]] for q in ("n_held", "n_missing", "n_nonfinite") if c0U[q] != cU[q]}}
            diff = (h0 != h1) | (f0 != f1)
            if diff.any():
                ar["change"]["x_mask_diff"] = [float(v) for v in xu[diff][:8]]
            for key, m in (("F", F & both), ("use", both)):
                if m.any():
                    j = int(np.flatnonzero(m)[int(np.argmax(dd[m]))])
                    ar[f"max_mm_{key}"], ar[f"x_at_{key}"] = float(dd[j]), float(xu[j])
            if cF["pf"] != c0F["pf"] or cU["pf"] != c0U["pf"] or ar["change"]:
                undec.append(f"{k} × {a.get('res')}: 合否 F {c0F['pf']}→{cF['pf']}・使用域 {c0U['pf']}→{cU['pf']} {ar['change']}")
            if ar["max_mm_F"] is not None and (out["at_F"] is None or ar["max_mm_F"] > worst):
                worst = ar["max_mm_F"]
                out["at_F"] = {"step": k, "euler_res": a.get("res"), "x": ar["x_at_F"], "mm": worst}
            if ar["max_mm_use"] is not None and (out["at_use"] is None or ar["max_mm_use"] > worst_u):
                worst_u = ar["max_mm_use"]
                out["at_use"] = {"step": k, "euler_res": a.get("res"), "x": ar["x_at_use"], "mm": worst_u}
            rec["alts"].append(ar)
        out["per_step"][k] = rec
    out["mm"], out["mm_use"] = (worst, worst_u)
    if unmeas:
        out.update(status="unmeasurable", why="; ".join(unmeas[:4]) + (" …" if len(unmeas) > 4 else ""))
    elif undec:
        out.update(status="undecidable", why="; ".join(undec[:4]) + (" …" if len(undec) > 4 else ""))
    if out["status"] != "ok":
        out["mm_observed"], out["mm_use_observed"], out["mm"], out["mm_use"] = out["mm"], out["mm_use"], None, None
    return out


def sigma_ref_measure(run: Path, steps: list, spec_b: dict, alts: list, ext: dict, x_use, F, dB: dict, mm: float, snap, jobs: int = 1) -> dict:
    """σ_ref を測る (抽出して `sigma_ref_eval`): 判定窓の全出力 k × alts (代わりの参照 {res, step, sha256}) を、腕 B の固定の Euler の出力の
    代わりにして抽出し直す (`run_snaps`、jobs > 1 ならプロセス並列)。dB[k] は固定の参照の抽出。"""
    tasks, keys = [], []
    for k in steps:
        for a in alts:
            tasks.append((run, k, {"run": spec_b["run"], "res": a["res"], "sha256": a["sha256"]}, ext, x_use))
            keys.append((k, a["res"]))
    res, info = run_snaps(tasks, snap, jobs)
    out = sigma_ref_eval(steps, alts, dB, dict(zip(keys, res)), x_use, F, mm)
    out["parallel"] = info
    return out


def ab_full(run: Path, s: int, specs: dict, ext: dict, problem: Path, x_w: float, mm: float) -> dict:
    """停止の出力 s の両腕の全体の抽出 (run/H370_AB/<腕>) と、帯の位置・感度・held・δ_target・次の壁の比較 (記録)。"""
    from forge_design.feedback import deltastar_loop as DL
    ext_fn = HOOKS.get("extract_and_merge", DL.extract_and_merge)
    nw_fn = HOOKS.get("next_wall_check", DL.next_wall_check)
    out, tabs = {}, {}
    for arm, spec in specs.items():
        od = run / AB_DIR / arm
        if od.exists():
            shutil.rmtree(od)
        od.mkdir(parents=True)
        try:
            with snapshot_view(run, s) as v:
                summ = ext_fn(v, spec, out_dir=od, **ext)
        except Exception as e:  # noqa: BLE001
            out[arm] = {"error": f"{type(e).__name__}: {e}"[:600]}
            continue
        ne = summ.get("not_extracted") or {}
        out[arm] = {"n_stations": summ.get("n_stations"), "n_held": summ.get("n_held"), "lam": (summ.get("smooth") or {}).get("lam"),
                    "not_extracted": {"n": ne.get("n"), "x_range": [min(ne["x"]), max(ne["x"])] if ne.get("x") else None,
                                      "reasons": sorted(set(ne.get("reason") or []))}}
        tabs[arm] = {"next": np.loadtxt(od / "delta_r_next.csv", delimiter=",", skiprows=1, ndmin=2),
                     "equiv": np.loadtxt(od / "delta_r_equiv.csv", delimiter=",", skiprows=1, ndmin=2),
                     "edge": (np.loadtxt(od / "delta_r_edge.csv", delimiter=",", skiprows=1, ndmin=2) if (od / "delta_r_edge.csv").is_file() else None),
                     "head": (od / "delta_r_next.csv").read_text().splitlines()[0]}
    if set(tabs) != {"A", "B"}:
        out["compare"] = {"error": "どちらかの腕の抽出が失敗"}
        return out
    A, B = tabs["A"], tabs["B"]
    xa, xb = A["next"][:, 0], B["next"][:, 0]
    use = lambda x: (x >= 0.0) & (x <= x_w + 1e-9)  # noqa: E731
    ca = np.isin(xa, xb) & use(xa)
    cb = np.isin(xb, xa) & use(xb)
    cmp = {"n_common_use": int(ca.sum()), "x_common": [float(xa[ca].min()), float(xa[ca].max())] if ca.any() else None,
           "B_use_not_in_A": {"n": int((use(xb) & ~np.isin(xb, xa)).sum()),
                              "x_range": ([float(xb[use(xb) & ~np.isin(xb, xa)].min()), float(xb[use(xb) & ~np.isin(xb, xa)].max())]
                                          if (use(xb) & ~np.isin(xb, xa)).any() else None)}}
    if ca.any():
        ia, ib = np.flatnonzero(ca), np.flatnonzero(cb)
        for name, col in (("delta_E", 3), ("delta_target", 1)):
            d = np.abs(A["next"][ia, col] - B["next"][ib, col]) * mm
            cmp[f"{name}_max_diff_mm"] = float(d.max())
            cmp[f"{name}_x_at_max"] = float(xa[ia][int(np.argmax(d))])
        cmp["held_common"] = {"A": int(A["next"][ia, 4].sum()), "B": int(B["next"][ib, 4].sum())}
        ea, eb = A["equiv"], B["equiv"]
        ja, jb = np.flatnonzero(np.isin(ea[:, 0], xa[ia])), np.flatnonzero(np.isin(eb[:, 0], xb[ib]))
        if len(ja) == len(jb) and len(ja):
            sa = np.max(np.abs(ea[ja][:, 7:9] - ea[ja][:, [1]]), axis=1) / np.maximum(np.abs(ea[ja][:, 1]), 1e-30)
            sb = np.max(np.abs(eb[jb][:, 7:9] - eb[jb][:, [1]]), axis=1) / np.maximum(np.abs(eb[jb][:, 1]), 1e-30)
            cmp["band_sensitivity_rel"] = {"A_max": float(np.nanmax(sa)), "B_max": float(np.nanmax(sb)), "max_abs_diff": float(np.nanmax(np.abs(sa - sb)))}
        if A["edge"] is not None and B["edge"] is not None:
            ga, gb = A["edge"], B["edge"]
            ka, kb = np.flatnonzero(np.isin(ga[:, 0], xa[ia])), np.flatnonzero(np.isin(gb[:, 0], xb[ib]))
            if len(ka) == len(kb) and len(ka):
                for name, col in (("y_b", 5), ("y_star", 2)):
                    d = np.abs(ga[ka, col] - gb[kb, col])
                    fin = np.isfinite(d)
                    cmp[f"band_{name}_max_diff_mm"] = float(np.max(d[fin]) * mm) if fin.any() else None
                    cmp[f"band_{name}_max_rel"] = (float(np.max(d[fin] / np.maximum(np.abs(gb[kb, col][fin]), 1e-30))) if fin.any() else None)
                    cmp[f"band_{name}_n_nonfinite"] = int((~fin).sum())
        # 未緩和の固定点の残差 (設計を終える条件の量、記録だけ): 腕 B の使用域の max|δ_E − δ_in|
        ub = use(xb)
        r = np.abs(B["next"][ub, 3] - B["next"][ub, 2]) * mm
        cmp["B_fixed_point_residual_mm"] = {"max": float(r.max()), "x_at_max": float(xb[ub][int(np.argmax(r))]), "note": "0.05 mm は設計を終える条件 (ここでは記録だけ)"}
    out["compare"] = cmp
    # 次の壁: 腕 B の表 / 腕 A の共通の列の δ_target だけを替えた表 (出口の区間は腕 B のまま = 延ばさない)
    nw = {}
    try:
        nw["B"] = nw_fn(problem, run / AB_DIR / "B" / "delta_r_next.csv")
        hyb = B["next"].copy()
        if ca.any():
            hyb[np.flatnonzero(cb), 1] = A["next"][np.flatnonzero(ca), 1]
        hd = run / AB_DIR / "A_on_common"
        hd.mkdir(parents=True, exist_ok=True)
        np.savetxt(hd / "delta_r_next.csv", hyb, delimiter=",", comments="", header=B["head"] + " [H370: 共通の列の delta_target を腕 A に替えた表]",
                   fmt="%.18e")
        nw["A_on_common"] = nw_fn(problem, hd / "delta_r_next.csv")
        ta, tb = (nw["A_on_common"] or {}).get("throat"), (nw["B"] or {}).get("throat")
        if ta and tb:
            nw["throat_diff"] = {"x_rt": ta["x"] - tb["x"], "r_mm": (ta["r"] - tb["r"]) * mm}
    except Exception as e:  # noqa: BLE001
        nw["error"] = f"{type(e).__name__}: {e}"[:600]
    out["next_wall"] = nw
    return out


def ab_compare(case: Path, run: Path, s: int, rec: dict, spec_b: dict, qual: dict, problem: Path, x_w: float, jobs: int = 1) -> dict:
    """抽出の参照の A/B と σ_ref・σ_ref,use (上の docstring の 5、§5.1 #36g (2)〜(5)・#36i)。腕 B の抽出し直し・腕 A・σ_ref・スクリーニングの
    抽出は 1 つの作業の一覧にまとめて `run_snaps` で回す (jobs > 1 ならプロセス並列; 結果は作業の順に並べるので記録は直列と同じ)。"""
    from forge_design.evaluate import stop_rule as SR
    from forge_design.feedback import deltastar_loop as DL
    snap = HOOKS.get("delta_snapshot", SR.delta_snapshot)
    st, ds = rec.get("settings") or {}, rec.get("delta_series") or {}
    W = int(st.get("window", 0) or 0)
    x_use = np.asarray(ds.get("x") or [], dtype=float)
    S = ds.get("scale_m")
    out = {"tol_mm": TOL_MM, "arm_B": spec_b, "windows": [[s - 2 * W, s - W], [s - W, s]], "result": "判定不能", "why": None}
    ext = {k: v for k, v in (ds.get("extract") or DL.handoff_extract_settings()).items() if k != "source"}
    out["extract"] = ext
    old = case / OLD_EULER
    ors = res_steps(old) if old.is_dir() else []
    if not ors or not len(x_use) or S is None or not W:
        out["why"] = f"腕 A の {OLD_EULER} の出力が無いか、stop_rule.json に窓・δ の列が無い"
        return out
    mm = float(S) * 1e3
    # --- (4) 判定窓の完全性: 入るべき出力を出力の間隔と 2 窓から列挙する (記録にある行だけから作らない)
    oi = st.get("out_step_interval")
    try:
        required = SR.required_window_steps(int(s), W, oi)
    except (TypeError, ValueError) as e:
        out["why"] = f"判定窓に入るべき出力を列挙できない (settings.out_step_interval {oi!r}: {e})"
        return out
    rows = {r["step"]: r for r in rec.get("rows") or [] if isinstance(r.get("step"), int)}
    comp = {"required": required, "n_required": len(required), "out_step_interval": int(oi),
            "rows_missing": [k for k in required if k not in rows], "res_missing": [k for k in required if not (run / f"res_{k}.h5").is_file()]}
    out["completeness"] = comp
    spec_a = {"run": str(old.resolve()), "res": f"res_{ors[-1]}.h5", "sha256": sha256_file(old / f"res_{ors[-1]}.h5")}
    out["arm_A"] = spec_a
    # --- (3) 比べる列 F を抽出の前に座標と各腕の Euler の補間できる範囲から決める
    try:
        fx = fixed_columns(x_use, {"A": old, "B": Path(spec_b["run"])}, edge_tol=(ext.get("remap_edge") == "missing"))
    except Exception as e:  # noqa: BLE001
        out["why"] = f"比べる列を決められない ({type(e).__name__}: {e})"[:400]
        return out
    F = fx.pop("mask")
    out["columns"] = {**fx, "note": "抽出の前に座標だけで決めた比べる列 (NS の使用域の列 ∩ 両腕の Euler の補間できる範囲)。この中の欠測・held・"
                                    "非有限は比較を合格にしない"}
    if not F.any():
        out["why"] = "比べる列 F が空 (両腕の Euler の範囲が重ならない)"
        return out
    # --- 抽出の作業の一覧 (腕 B の抽出し直し → 腕 A → σ_ref → スクリーニングの順; 直列のときもこの順に回す)
    tasks, keys, plan = [], [], {}
    for k in required:
        f = run / f"res_{k}.h5"
        if not f.is_file():
            plan[k] = None
            continue
        sha_k = sha256_file(f)
        dB = (rows.get(k) or {}).get("delta")
        if k not in rows:
            src, dB = "抽出し直し (stop_rule.json に行が無い)", None
        elif isinstance(dB, dict) and "delta_E" in dB and len(dB["delta_E"]) == len(x_use) and dB.get("res_sha256") == sha_k:
            src = "記録 (stop_rule.json の delta)"
        else:
            src, dB = "抽出し直し (記録が無いか出力の sha が違う)", None
        plan[k] = {"src": src, "dB": dB}
        if dB is None:
            tasks.append((run, k, spec_b, ext, x_use))
            keys.append(("B", k, None))
        tasks.append((run, k, spec_a, ext, x_use))
        keys.append(("A", k, None))
    present = [k for k in required if plan.get(k) is not None]
    # σ_ref の代わりの参照の照合 (§5.1 #36i M1: 末尾の窓の全 20 出力ちょうど・在る・sha が記録と同じ)。合わなければ σ_ref の抽出は回さない
    alts, achk = sigma_ref_alternatives(spec_b, qual)
    if achk["ok"]:
        for k in present:
            for a in alts:
                tasks.append((run, k, {"run": spec_b["run"], "res": a["res"], "sha256": a["sha256"]}, ext, x_use))
                keys.append(("S", k, a["res"]))
    scr_snaps = (qual or {}).get("screening_snapshots") or []
    for sn in scr_snaps:
        tasks.append((run, s, {"run": spec_b["run"], "res": sn.get("res"), "sha256": sn.get("sha256")}, ext, x_use))
        keys.append(("scr", s, sn.get("res")))
    prog_p = run / "h370_judge_progress.log"

    def progress(m, n, el):
        if m == n or m % 5 == 0:
            with open(prog_p, "a") as fh:
                fh.write(f"{now()} 抽出 {m}/{n} (経過 {el:.0f} s、残りの見込み {el / m * (n - m):.0f} s、並列 {jobs})\n")
    results, pinfo = run_snaps(tasks, snap, jobs, progress=progress)
    got = dict(zip(keys, results))

    def as_d(r):
        return r["ok"] if "ok" in r else {"error": r["exc"][:300]}
    per, dB_by = [], {}
    for k in required:
        pk = plan.get(k)
        if pk is None:
            per.append({"step": k, "error_missing": f"出力 res_{k}.h5 が無い (再抽出で補えない)"})
            continue
        dB = pk["dB"] if pk["dB"] is not None else as_d(got[("B", k, None)])
        dA = as_d(got[("A", k, None)])
        dB_by[k] = dB
        per.append(ab_pair(k, dA, dB, mm, pk["src"], x_use, F))
    comp["reextracted_B"] = [r["step"] for r in per if str(r.get("B_source", "")).startswith("抽出し直し")]
    out["per_output"] = per
    # --- (2) σ_ref・σ_ref,use: 判定窓の全出力 × 末尾の窓の全 20 出力
    if achk["ok"]:
        sig = sigma_ref_eval([k for k in required if k in dB_by], alts, dB_by, {(k, r_): got[(t_, k, r_)] for (t_, k, r_) in keys if t_ == "S"},
                             x_use, F, mm)
    else:
        sig = {"status": "unmeasurable", "mm": None, "mm_use": None, "at_F": None, "at_use": None, "alts": [a["res"] for a in alts],
               "per_step": {}, "seconds": 0.0, "why": f"σ_ref の代わりの参照が揃わない (抽出していない): {achk.get('why')}"}
    sig["alternatives_check"] = achk
    sig["selection"] = (qual or {}).get("sigma_ref_selection")
    sig["note"] = ("判定窓の全出力を、固定した Euler の出力の代わりに末尾の窓の全出力 (予備の資格の時点で固定した 20 出力と sha) で抽出し直した "
                   "δ_E の差 (mm は比べる列 F、mm_use は使用域の全体 U = 固定の参照で抽出できる使用域の全列; at_F・at_use は最大を与えた "
                   "NS の出力・Euler の出力・x)")
    out["sigma_ref"] = sig
    # スクリーニング (停止の出力 × check_quasisteady の末尾の窓の始まり・中ほど): 記録だけ、σ_ref の上限とは扱わない
    scr = []
    d_s = dB_by.get(s)
    for sn in scr_snaps:
        r = got[("scr", s, sn.get("res"))]
        if "exc" in r:
            scr.append({"res": sn.get("res"), "error": r["exc"][:300]})
            continue
        d = r["ok"]
        e0, e1 = _arr(d_s or {}, "delta_E"), _arr(d, "delta_E")
        m = F & np.isfinite(e0) & np.isfinite(e1) if len(e0) == len(e1) == len(F) else np.zeros(0, bool)
        scr.append({"res": sn["res"], "max_mm_F": float(np.max(np.abs(e1 - e0)[m]) * mm) if m.any() else None})
    out["screening"] = {"per_snapshot": scr, "label": "スクリーニング (停止の出力 × 2 出力だけ。σ_ref の上限とは扱わない、記録だけ)"}
    ok = [r for r in per if r.get("D_mm") is not None]
    if ok:
        i_max = max(ok, key=lambda r: r["D_mm"])
        out["D_mm"] = {"max": i_max["D_mm"], "min": min(r["D_mm"] for r in ok), "n_outputs": len(per), "n_with_D": len(ok),
                       "x_at_max": float(x_use[i_max["i_at_D"]]), "step_at_max": i_max["step"], "x_fixed": out["columns"]["x_range"]}
    out["classify"] = classify_ab(per, sig, required)
    out["result"], out["why"] = out["classify"]["result"], out["classify"]["why"]
    n_a = sum(1 for r in per if "error_missing" not in r)
    out["extractions"] = {"arm_A": n_a, "arm_B_reextracted": len(comp["reextracted_B"]),
                          "sigma_ref": sum(1 for (t_, _k, _r) in keys if t_ == "S"), "screening": len(scr),
                          "sigma_ref_seconds": round(sig.get("seconds", 0.0), 1), "parallel": pinfo,
                          "note": "停止の出力の全体の抽出 (at_stop、両腕) は別に 2 回 (直列)"}
    try:
        out["at_stop"] = ab_full(run, s, {"A": spec_a, "B": spec_b}, ext, problem, x_w, mm)
    except Exception as e:  # noqa: BLE001 — 記録だけ
        out["at_stop"] = {"error": f"{type(e).__name__}: {e}"[:600]}
    return out


def gate_summary(g: dict) -> dict:
    c = g.get("checks") or {}
    pick = lambda d, ks: {k: d.get(k) for k in ks} if isinstance(d, dict) else d  # noqa: E731
    return {"proceed": g.get("proceed"), "blocked_by": g.get("blocked_by"), "stop_rule": g.get("stop_rule"),
            "convergence": pick(g.get("convergence") or {}, ("kind", "rc", "verdict", "cmd")),
            "conservation": pick(c.get("conservation"), ("pass", "r_mean", "r_amp", "n", "why")),
            "wall_resolution": pick(c.get("wall_resolution"), ("pass", "kind", "verdict", "cmd")),
            "two_windows": pick(c.get("two_windows"), ("pass", "windows", "n", "worst", "why")),
            "euler_ref": c.get("euler_ref"), "omega": c.get("omega"),
            "coverage": pick(c.get("coverage"), ("pass", "n_use", "n_held_use", "n_nonfinite_use", "not_extracted", "why")),
            "next_wall": pick(c.get("next_wall"), ("pass", "throat", "upstream_gate_pass", "validate", "why")),
            "made": (g.get("extract") or {}).get("made"), "override": g.get("override")}


def ns_judge(case: Path, tools: Path, do_ab: bool, out_p: Path, jobs=None) -> int:
    from forge_design.evaluate import stop_rule as SR
    from forge_design.feedback import deltastar_loop as DL
    run = case / NS_RUN
    R = {"plan": PLAN, "created": now(), "run": NS_RUN, "case_dir": str(case), "tools": str(tools), "undecidable": [], "records": {},
         "note": "壁は更新しない (固定壁の pass)。0.05 mm の固定点の条件は設計を終える用で、ここでは記録だけ"}
    und = R["undecidable"]

    def finish(rc):
        jwrite(out_p, R)
        print(f"== VERDICT V-h370: {R.get('verdict')}")
        print(f"→ {out_p}")
        return rc
    if not run.is_dir():
        und.append(f"{NS_RUN} が無い")
        R["verdict"] = "判定不能: " + "; ".join(und)
        return finish(2)
    if (run / DIAG_MARK).exists() or (jload(run / "H370_PRECHECK.json") or {}).get("phase") == "ns-diag":
        # 診断専用の run (§6 V-h370-diag): 正式の判定をしない・受け渡しの関門を呼ばない (感度の判定は --diag-sens)
        und.append(f"診断専用の run ({DIAG_MARK}、§6 V-h370-diag) — 正式の判定をせず受け渡しの関門を呼ばない (感度の判定は --diag-sens)")
        R["update_conditions_met"] = False
        R["verdict"] = "判定不能: " + "; ".join(und)
        return finish(2)
    pin, qual = jload(case / PIN), jload(case / QUAL)
    spec = {"run": str(case / pin["run"]), "res": pin["res"], "sha256": pin["sha256"]} if pin else None
    R["euler_ref"] = spec
    if not pin or not (qual or {}).get("prequalified"):
        und.append("Euler の参照の固定・予備の資格の記録が無い")
    # --- 1. 記録 ---
    pc, ms = jload(run / "H370_PRECHECK.json") or {}, jload(run / "H370_MAIN_START.json") or {}
    R["records"]["precheck_ok"], R["records"]["main_start_ok"] = pc.get("ok"), ms.get("ok")
    if not pc.get("ok"):
        und.append("投入前の検査 (H370_PRECHECK.json) が無いか不成立")
    if not ms.get("ok"):
        und.append("本段の直前の照合 (H370_MAIN_START.json) が無いか不成立")
    for m in MARKERS:
        if (run / m).is_file():
            und.append(f"基盤の停止・不成立の印 {m}: {(run / m).read_text()[:200]}")
    R["records"]["run_rc"] = (run / "RUN_RC").read_text().strip() if (run / "RUN_RC").is_file() else None
    prov = (run / "RUN_PROVENANCE.txt").read_text() if (run / "RUN_PROVENANCE.txt").is_file() else ""
    fsha = next((ln.split(":", 1)[1].strip() for ln in prov.splitlines() if ln.startswith("forge_sha256")), None)
    envs = [ln.split(":", 1)[1].strip() for ln in prov.splitlines() if ln.startswith("env")]
    R["records"]["provenance"] = {"forge_sha256": fsha, "env": envs}
    if not (fsha or "").startswith(BIN_SHA) or "FORGE_CUDA_BLOCKSIZE=128" not in envs:
        und.append(f"RUN_PROVENANCE.txt のバイナリ ({str(fsha)[:16]}) か FORGE_CUDA_BLOCKSIZE が登録と違う")
    R["records"]["residual_nonfinite"] = {t: residual_nonfinite(run / f"residual_history_{t}.csv") for t in ("S1_soft", "S2_mid", "main")}
    # --- 2. 停止の状態 ---
    state = SR.read_state(run) or {}
    rec = jload(run / SR.RECORD) or {}
    status = state.get("status") or "MISSING"
    s = state.get("stop_step") if state.get("stop_step") is not None else rec.get("stop_step")
    rows = [r for r in rec.get("rows") or [] if isinstance(r.get("step"), int)]
    last = rows[-1] if rows else {}
    R["stop"] = {"status": status, "stop_step": s, "last_output": state.get("last_output"), "forge_rc": state.get("forge_rc"),
                 "hold": rec.get("hold"), "note": state.get("note"), "sections_requested": rec.get("sections"),
                 "sections_chosen_x": {k: v for k, v in last.items() if k.startswith("x_theta_r_")}, "recorded_at_stop": rec.get("recorded_at_stop")}
    if status == "REACHED" and s is not None:
        f = run / f"res_{s}.h5"
        R["stop"]["stop_sha256_match"] = f.is_file() and sha256_file(f) == state.get("stop_sha256")
        if not R["stop"]["stop_sha256_match"]:
            und.append(f"停止の出力 res_{s}.h5 が無いか sha が STOP_RULE_STATE.json と違う")
        # 壁の出力の名前は投入前の検査の記録 (wall_output_check) から、無ければ case/45 の bcond の壁 (wall・physID 3)
        wname = (((pc.get("items") or {}).get("N3_wall_output") or {}).get("wall_dump_name") or "res_wall_3_<step>.h5").replace("<step>", str(s))
        R["stop"]["wall_dump_at_stop"] = sorted(p.name for p in run.glob(f"res_*_*_{s}.h5"))
        R["stop"]["wall_dump_required"] = wname
        if not (run / wname).is_file():
            und.append(f"停止の出力 s = {s} の壁の出力 {wname} が無い (局所 y₁⁺ を判定できない)")
    # --- 3. 受け渡しの関門 ---
    gate = None
    if spec is not None:
        try:
            gate = HOOKS.get("handoff_gate", DL.handoff_gate)(run, problem=case / NS_PROBLEM, euler_ref=spec, omega=DL.HANDOFF_OMEGA,
                                                             allow_not_reached=False, allow_unconverged=False, write=True)
            R["gate"] = gate_summary(gate)
        except Exception as e:  # noqa: BLE001
            und.append(f"受け渡しの関門が例外 ({type(e).__name__}: {e})"[:400])
    # --- 4. check_quasisteady (記録) ---
    W = int((rec.get("settings") or {}).get("window", 0) or 0)
    if s is not None and W and rows:
        from forge_design.evaluate.stop_rule import theta_key
        xs_th = [theta_key(f) for f in (rec.get("settings") or {}).get("x_fracs") or []]
        x_use = np.asarray((rec.get("delta_series") or {}).get("x") or [], dtype=float)
        ser = []
        for r in rows:
            if not (s - 2 * W <= r["step"] <= s) or r.get("nonfinite_all"):
                continue
            o = {"step": r["step"], "r": r.get("r", float("nan")), **{k: r.get(k, float("nan")) for k in xs_th},
                 "exit_core_M": r.get("exit_core_M", float("nan")), "mdot_rt": r.get("mdot_rt", float("nan"))}
            d = r.get("delta")
            for x in XS:
                ok = isinstance(d, dict) and "delta_E" in d and len(d["delta_E"]) == len(x_use) and len(x_use)
                o[f"dE_{x:g}"] = float(np.interp(x, x_use, d["delta_E"])) if ok else float("nan")
                o[f"dT_{x:g}"] = float(np.interp(x, x_use, d["delta_target"])) if ok else float("nan")
            ser.append(o)
        cols = ["r"] + xs_th + ["exit_core_M", "mdot_rt"] + [f"d{q}_{x:g}" for q in ("E", "T") for x in XS]
        if ser:
            sc = run / "h370_qs_series.csv"
            write_series(sc, ser, cols)
            R["quasisteady"] = {"window": [s - 2 * W, s], "n": len(ser), "default": run_qs(tools, sc, cols),
                                "tight (1e-4)": run_qs(tools, sc, cols, QS_TIGHT),
                                "note": "記録だけ (§6 V-h370: 既定の相対閾値で寸法の予算の合格を代用しない; 判定は関門の 2 窓)"}
    # --- 5. 抽出の参照の A/B ---
    if do_ab and spec is not None and s is not None and status in ("REACHED", "CENSORED") and rows:
        info = jload(run / "prepare_info.json") or {}
        x_w = float(((info.get("physical_wall") or {}).get("domain") or [0.0, float("inf")])[1])
        try:
            pj = snap_jobs(jobs)
            R["parallel_plan"] = {**pj, "rule": f"並列数 = min(上限 {pj['cap']}, 空きメモリ {pj['mem_available_mib']} MiB / 1 抽出の最大 "
                                                 f"{SNAP_PEAK_MIB} MiB = {pj['by_memory']}) (§5.1 #36i)"}
            R["ab"] = ab_compare(case, run, int(s), rec, spec, qual, case / NS_PROBLEM, x_w, jobs=pj["jobs"])
            if status != "REACHED":
                R["ab"]["note"] = f"停止の状態 {status} の窓で行った診断の抽出 (更新の許可とは別)"
        except Exception as e:  # noqa: BLE001
            R["ab"] = {"result": "判定不能", "why": f"A/B が例外 ({type(e).__name__}: {e})"[:400]}
    else:
        R["ab"] = {"result": "判定不能" if do_ab else "未実施 (--no-ab)",
                   "why": f"停止の状態 {status}・停止の出力 {s} では A/B の窓を作れない" if do_ab else None}
    # --- 参照の最終の資格 (§5.1 #36g (1)・#36i: σ_ref・σ_ref,use を測れない・判定不能・σ_ref,use > 0.01 mm・参照の未整定・A/B が判定不能
    #     なら更新に進めない) ---
    rf = reference_final(R["ab"], do_ab)
    R["reference_final"] = rf
    # --- 判定 ---
    proceed = bool(gate and gate.get("proceed"))
    met = proceed and rf.get("qualified") is True
    blk = ([] if proceed else ((gate or {}).get("blocked_by") or ["関門の記録が無い"])) + \
          ([] if rf.get("qualified") is True else [f"Euler の参照の最終の資格 (σ_ref・σ_ref,use): {rf.get('why')}"])
    upd = ("更新に進める条件をすべて満たす (関門と参照の最終の資格; 壁の更新は自動でつながない)" if met
           else "更新に進める条件を満たさない: " + "; ".join(blk)[:900])
    R["update_conditions_met"] = bool(met and not und)
    abr = R["ab"].get("result")
    if und:
        R["verdict"] = ("判定不能: " + "; ".join(und)[:900] + f" (記録: 関門 {('進める' if proceed else '止める') if gate else 'なし'}、"
                        f"参照の最終の資格 {rf.get('qualified')}、A/B {abr})")
        return finish(2)
    R["verdict"] = f"{upd} / 抽出の参照の A/B: {abr}" + (f" ({R['ab'].get('why')})" if R["ab"].get("why") else "")
    return finish(0)


# ===================================================================================================== 診断の感度の判定 (§6 V-h370-diag)
def module_versions() -> dict:
    """抽出器・見張り・準備の版: 読み込んだ forge_design の VERSION_FILES の sha256 (h370_run.module_versions と同じ)。"""
    import forge_design
    root = Path(forge_design.__file__).resolve().parent
    return {n: (sha256_file(root / n) if (root / n).is_file() else None) for n in VERSION_FILES}


def diag_snapshot(run_dir, step: int, euler_ref, extract: dict, x_use) -> dict:
    """診断の感度の判定の抽出 1 回 (§6 V-h370-diag): `stop_rule.delta_snapshot` と同じ抽出 (prepare_info.json・nozzle.h5・その res だけを
    symlink した一時ディレクトリ、条件 `extract` = 受け渡しの固定の条件 (edge・一定 x の断面・端は欠測・ω・平滑化の λ・昇格の上限)、
    使用域の列 x_use への x の一致の対応) に、選ばれた帯の位置 (delta_r_edge.csv の y_b・縁 y*) を足して返す。δ_E・δ_target・held・missing
    は delta_snapshot と同じ値 (試験 L で照合)。抽出の表は一時ディレクトリに書いて捨てる (受け渡しの表として残さない)。帯の位置は参照
    ごとに抽出器が計算し直す (帯の選び方 edge は固定)。"""
    import time
    from forge_design.feedback.deltastar_loop import extract_and_merge, read_delta_r_next
    run_dir = Path(run_dir).resolve()
    res = run_dir / f"res_{step}.h5"
    t0 = time.time()
    kw = {k: v for k, v in extract.items() if k != "source"}
    with tempfile.TemporaryDirectory(prefix="h370_diag_") as td:
        src = Path(td) / run_dir.name
        src.mkdir()
        for name in ("prepare_info.json", "nozzle.h5"):
            (src / name).symlink_to(run_dir / name)
        (src / res.name).symlink_to(res)
        summ = extract_and_merge(src, euler_ref, out_dir=Path(td) / "out", **kw)
        nx = read_delta_r_next(Path(td) / "out" / "delta_r_next.csv")
        ep = Path(td) / "out" / "delta_r_edge.csv"
        edge = np.loadtxt(ep, delimiter=",", skiprows=1, ndmin=2) if ep.is_file() else None
    x_use = np.asarray(x_use, dtype=float)
    xe = np.atleast_1d(np.asarray(nx["x_rt"], dtype=float))
    if len(xe) == 0:
        raise ValueError(f"res_{step}.h5 の抽出の表が空")
    k = np.clip(np.searchsorted(xe, x_use), 0, len(xe) - 1)
    hit = xe[k] == x_use
    dE = np.where(hit, np.atleast_1d(nx["delta_E"])[k], np.nan)
    dT = np.where(hit, np.atleast_1d(nx["delta_target"])[k], np.nan)
    held = np.where(hit, np.atleast_1d(nx["held"])[k] != 0, False)
    out = {"delta_E": dE.tolist(), "delta_target": dT.tolist(), "held": [int(v) for v in held],
           "missing": [int(j) for j in np.where(~hit)[0]], "n_missing": int((~hit).sum()), "n_held": int(held.sum()),
           "lam": summ["smooth"]["lam"], "res_sha256": sha256_file(res), "seconds": round(time.time() - t0, 2),
           "not_extracted_n": (summ.get("not_extracted") or {}).get("n")}
    if edge is not None and len(edge):
        o = np.argsort(edge[:, 0], kind="stable")
        xg = edge[o, 0]
        j = np.clip(np.searchsorted(xg, x_use), 0, len(xg) - 1)
        hb = xg[j] == x_use
        out["y_b"] = np.where(hb, edge[o, 5][j], np.nan).tolist()
        out["y_star"] = np.where(hb, edge[o, 2][j], np.nan).tolist()
    else:
        out["y_b"] = out["y_star"] = [float("nan")] * len(x_use)
        out["band_note"] = "delta_r_edge.csv が無い (帯の位置を記録できない)"
    return out


def diag_pair(step: int, dA: dict, dB: dict, mm: float, x_use, F) -> dict:
    """出力 1 つの腕 A・B の比較 (§6 V-h370-diag): 合否 pf と D_mm (両腕とも pf のとき F の max|δ_E,B − δ_E,A| [mm]) は `ab_pair` と同じ。
    足すもの: D_partial_mm (どちらかが pf でないとき、両腕とも有限・held でない F の列だけの最大; D の下限なので 0.01 mm を超えれば結果 B の
    根拠になる)、δ_target の差、帯の位置 (y_b・y*) の差、held・欠測の列の数。"""
    r = ab_pair(step, dA, dB, mm, "抽出 (診断)", x_use, F)
    if any(k.startswith("error") for k in r):
        return r
    xu = np.asarray(x_use, dtype=float)
    eA, eB, tA, tB = _arr(dA, "delta_E"), _arr(dB, "delta_E"), _arr(dA, "delta_target"), _arr(dB, "delta_target")
    hA, hB = np.asarray(dA.get("held", []), dtype=bool), np.asarray(dB.get("held", []), dtype=bool)
    ok = F & np.isfinite(eA) & np.isfinite(eB) & ~hA & ~hB
    r["n_valid_both"] = int(ok.sum())
    if ok.any():
        d = np.abs(eB - eA) * mm
        j = int(np.flatnonzero(ok)[int(np.argmax(d[ok]))])
        r["D_valid_mm"], r["x_at_D_valid"] = float(d[j]), float(xu[j])
        okt = ok & np.isfinite(tA) & np.isfinite(tB)
        if okt.any():
            dt = np.abs(tB - tA) * mm
            jt = int(np.flatnonzero(okt)[int(np.argmax(dt[okt]))])
            r["Dt_valid_mm"], r["x_at_Dt_valid"] = float(dt[jt]), float(xu[jt])
        if r.get("D_mm") is None:
            r["D_partial_mm"] = r["D_valid_mm"]
    if r.get("D_mm") is not None:
        r["x_at_D"] = float(xu[r["i_at_D"]])
    for name in ("y_b", "y_star"):
        a, b = _arr(dA, name), _arr(dB, name)
        if len(a) == len(b) == len(F):
            m = F & np.isfinite(a) & np.isfinite(b)
            r[f"band_{name}_max_diff_mm"] = float(np.max(np.abs(b - a)[m]) * mm) if m.any() else None
            r[f"band_{name}_n_nonfinite_F"] = {"A": int((F & ~np.isfinite(a)).sum()), "B": int((F & ~np.isfinite(b)).sum())}
    return r


def classify_diag(per: list, required, tol: float = TOL_MM) -> dict:
    """診断の感度の分類 (§6 V-h370-diag; `classify_ab` (時刻ごとの U_time・σ_ref) は流用しない): 戻り {"result", "code", "why", "D_max"}。
    per: 出力ごとの `diag_pair` の戻り (無い出力は {"step", "error_missing"})。required: 判定の窓に入るべき 17 出力。順に:
    入るべき出力の集合が無い・欠ける → 判定不能 (入力の不足)、抽出の例外 → 判定不能、どれかの出力で入れ替えにより合否 (pf: 例外なし・
    比べる列 F の中に held・欠測・非有限なし) が変わる → 結果 B、D = 全出力 × F の max|δ_E,B − δ_E,A| が tol を超える (どちらかが pf でない
    出力は両腕とも有効な列だけの最大 = D の下限で見る) → 結果 B、両腕とも pf でない出力がある (D を全列で出せない) → 判定不能、
    それ以外 (全出力・全列で両腕とも有効、D ≤ tol、合否が同じ) → 結果 A。**最大で決める**: 窓の中で D が tol をまたぐ時刻があっても
    max D > tol なら結果 B (時刻ごとに判定不能にしない)。"""
    def res(result, code, why, dmax=None):
        return {"result": result, "code": code, "why": why, "D_max": dmax}
    if not required:
        return res("判定不能", "U_required", "判定の窓に入るべき出力の集合が無い")
    have = {r["step"] for r in per if "error_missing" not in r}
    lack = [k for k in required if k not in have]
    if lack:
        return res("判定不能", "U_incomplete", f"判定の窓に入るべき {len(required)} 出力のうち比べられない出力 {len(lack)}: {lack[:8]} (縮めない)")
    per = [r for r in per if r["step"] in set(required)]
    err = [r["step"] for r in per if any(k.startswith("error") for k in r)]
    if err:
        return res("判定不能", "U_error", f"抽出の例外・比較の失敗 (出力 {err[:8]})")
    vals = [(r["D_mm"], r["step"], r.get("x_at_D")) for r in per if r.get("D_mm") is not None] + \
           [(r["D_partial_mm"], r["step"], r.get("x_at_D_valid")) for r in per if r.get("D_mm") is None and r.get("D_partial_mm") is not None]
    top = max(vals, key=lambda v: v[0]) if vals else None
    dmax = {"mm": top[0], "step": top[1], "x": top[2]} if top else None
    pfd = [r for r in per if r["pf_A"] != r["pf_B"]]
    over = [v for v in vals if v[0] > tol]
    if pfd or over:
        why = []
        if pfd:
            why.append("入れ替えで抽出の合否が変わる出力 " + "; ".join(f"{r['step']}: A {'合' if r['pf_A'] else '否'}・B {'合' if r['pf_B'] else '否'}"
                                                           for r in pfd[:6]))
        if over:
            why.append(f"D = max|δ_E,B − δ_E,A| {top[0]:.4g} mm > {tol} mm (出力 {top[1]}・x {top[2]}; tol を超える出力 {len(over)} / {len(per)})")
        return res("結果 B", "B_pf" if pfd else "B_D", "; ".join(why) + " — 「今回の参照の入れ替えが寸法の予算を超える」を退けない", dmax)
    both = [r["step"] for r in per if r.get("D_mm") is None]
    if both:
        return res("判定不能", "U_both_invalid", f"両腕とも比べる列 F の中に held・欠測・非有限がある出力 {both[:8]} (D を全列で出せず、両腕の有効な列の"
                                                  f"最大 {top[0] if top else None} mm ≤ {tol} mm)", dmax)
    return res("結果 A", "A", f"全 {len(per)} 出力・比べる列 F の全列で両腕とも有効、D = max|δ_E,B − δ_E,A| {top[0]:.4g} mm ≤ {tol} mm (出力 {top[1]}・"
                              f"x {top[2]})、合否が同じ — 「今回の参照の入れ替えが寸法の予算を超える」を退ける", dmax)


def diag_inputs(case: Path, R: dict, und: list):
    """診断の感度の判定の入力の照合 (§6 V-h370-diag の 1〜6; `diag_sens` (#36l) と帯を共通にした A/B `diag_sens_fixed_band` (#36n) で共通):
    診断の run であること・記録・終了の状態と停止の時刻 s・停止規則の記録の窓・抽出の条件・使用域の列・抽出器と見張りの版・判定の窓の 17 出力の
    完全性と sha の一覧・参照 A・B・aux の固定と sha。記録 (records・stop・window_kind・required・extract・versions・completeness・refs) を R に
    書き、不成立は und に足す。run が無ければ None、それ以外は {"run", "rec", "required", "ext", "x_use", "mm", "specs", "rdirs"}
    (und が空でなければ呼び出し側が判定不能で終える)。"""
    from forge_design.evaluate import stop_rule as SR
    from forge_design.feedback import deltastar_loop as DL
    run = case / NS_RUN
    if not run.is_dir():
        und.append(f"{NS_RUN} が無い")
        return None
    # --- 1. 診断の run であること・記録 ---
    mark, pc = jload(run / DIAG_MARK), jload(run / "H370_PRECHECK.json") or {}
    ms, end = jload(run / "H370_MAIN_START.json") or {}, jload(run / DIAG_END)
    R["records"] = {"diag_mark": bool(mark), "precheck_phase": pc.get("phase"), "precheck_ok": pc.get("ok"), "main_start_ok": ms.get("ok"),
                    "diag_end": {k: (end or {}).get(k) for k in ("end_state", "stop_step", "window", "main_outputs", "formal_records", "run_rc")}}
    if not mark or pc.get("phase") != "ns-diag":
        und.append(f"診断の印 {DIAG_MARK} が無いか投入前の検査が診断の経路 (ns-diag) のものでない — 診断の経路で回した run でない")
    if not pc.get("ok"):
        und.append("投入前の検査 (H370_PRECHECK.json) が無いか不成立")
    if not ms.get("ok"):
        und.append("本段の直前の照合 (H370_MAIN_START.json) が無いか不成立")
    for m in MARKERS:
        if (run / m).is_file():
            und.append(f"基盤の停止・不成立の印 {m}: {(run / m).read_text()[:200]}")
    prov = (run / "RUN_PROVENANCE.txt").read_text() if (run / "RUN_PROVENANCE.txt").is_file() else ""
    fsha = next((ln.split(":", 1)[1].strip() for ln in prov.splitlines() if ln.startswith("forge_sha256")), None)
    envs = [ln.split(":", 1)[1].strip() for ln in prov.splitlines() if ln.startswith("env")]
    R["records"]["provenance"] = {"forge_sha256": fsha, "env": envs}
    if not (fsha or "").startswith(BIN_SHA) or "FORGE_CUDA_BLOCKSIZE=128" not in envs:
        und.append(f"RUN_PROVENANCE.txt のバイナリ ({str(fsha)[:16]}) か FORGE_CUDA_BLOCKSIZE が登録と違う")
    # --- 2. 終了の状態と停止の時刻 s ---
    state = SR.read_state(run) or {}
    status = state.get("status") or "MISSING"
    s = state.get("stop_step")
    R["stop"] = {"status": status, "stop_step": s, "last_output": state.get("last_output"), "hold": state.get("hold"), "note": state.get("note")}
    if status not in ("REACHED", "CENSORED") or not isinstance(s, int):
        und.append(f"終了の状態 {status} (停止の時刻 {s!r}) — 停止の時刻が決まらない (REACHED か CENSORED であること)")
    elif end and (end.get("end_state"), end.get("stop_step")) != (status, s):
        und.append(f"終了の記録 {DIAG_END} の状態・停止の時刻 ({end.get('end_state')}, {end.get('stop_step')}) が STOP_RULE_STATE.json と違う")
    R["window_kind"] = ("停止の時刻の窓 [s − 40,000, s]" if status == "REACHED" else
                        "打切りの時点の診断の窓 (上限の打切り CENSORED: 最後の出力で終わる末尾の 17 出力; §6 V-h370-diag で先に登録)"
                        if status == "CENSORED" else None)
    # --- 3. 停止規則の記録の窓・抽出の条件・使用域の列 ---
    rec = jload(run / SR.RECORD) or {}
    st, ds = rec.get("settings") or {}, rec.get("delta_series") or {}
    if (st.get("window"), st.get("out_step_interval")) != (DIAG_WINDOW, DIAG_OUT_INT):
        und.append(f"停止規則の窓・出力の間隔 ({st.get('window')}, {st.get('out_step_interval')}) が登録 ({DIAG_WINDOW}, {DIAG_OUT_INT}) と違う")
    required = SR.required_window_steps(int(s), DIAG_WINDOW, DIAG_OUT_INT) if isinstance(s, int) and s >= 0 else []
    R["required"] = required
    if required and len(required) != DIAG_N_REQ:
        und.append(f"判定の窓の出力が {len(required)} ({DIAG_N_REQ} であること)")
    ext = {k: v for k, v in (ds.get("extract") or {}).items() if k != "source"}
    R["extract"] = ext
    if ext != DL.handoff_extract_settings():
        und.append(f"停止規則の記録の抽出の条件 {ext} が #36e で固定した受け渡しの条件 {DL.handoff_extract_settings()} と違う")
    try:
        uc = SR.use_columns(run)
        x_use = np.asarray(uc["x"], dtype=float)
        mm = float(uc["scale_m"]) * 1e3
        if [float(v) for v in ds.get("x") or []] != x_use.tolist():
            und.append("停止規則の記録の使用域の列 (delta_series.x) が run の格子の使用域の列と違う")
    except Exception as e:  # noqa: BLE001
        und.append(f"使用域の列を読めない ({type(e).__name__}: {e})"[:300])
        x_use, mm = np.zeros(0), None
    # --- 4. 抽出器・見張りの版 (投入の時点と同じ) ---
    vnow, vrun = module_versions(), pc.get("versions")
    R["versions"] = {"now": vnow, "at_launch": vrun, "same": vnow == vrun}
    if vnow != vrun:
        und.append("抽出器・見張りの版 (forge_design の該当ファイルの sha) が投入の時点 (H370_PRECHECK.json の versions) と違う")
    # --- 5. 判定の窓の 17 出力の完全性と sha の一覧 ---
    shal = jload(run / DIAG_SHA) or {}
    files = shal.get("files") or {}
    comp = {"required": required, "missing_res": [], "sha_mismatch": [], "not_in_sha_list": [], "sha_list": DIAG_SHA if shal else None}
    for k in required:
        f = run / f"res_{k}.h5"
        if not f.is_file():
            comp["missing_res"].append(k)
        elif f.name not in files:
            comp["not_in_sha_list"].append(k)
        elif sha256_file(f) != files[f.name]:
            comp["sha_mismatch"].append(k)
    R["completeness"] = comp
    if not shal:
        und.append(f"run の sha の一覧 {DIAG_SHA} が無い (診断の NS の終了の記録が書かれていない)")
    if comp["missing_res"]:
        und.append(f"判定の窓に入るべき出力 {comp['missing_res'][:8]} が無い (入力の不足; 縮めない)")
    if comp["not_in_sha_list"] or comp["sha_mismatch"]:
        und.append(f"判定の窓の出力の sha が終了の時点の一覧と合わない (一覧に無い {comp['not_in_sha_list'][:8]}・違う {comp['sha_mismatch'][:8]})")
    # --- 6. 参照の固定 (DIAG_REF) と sha の照合 ---
    dref = jload(case / DIAG_REF) or {}
    specs, rdirs = {}, {}
    for arm, r_ in (("A", (dref.get("sens_refs") or {}).get("A")), ("B", (dref.get("sens_refs") or {}).get("B")), ("aux", dref.get("euler_ref"))):
        if not isinstance(r_, dict) or not r_.get("run") or not r_.get("res") or not r_.get("sha256"):
            und.append(f"{DIAG_REF} に参照 {arm} の固定 {{run, res, sha256}} が無い")
            continue
        g = case / r_["run"] / r_["res"]
        if not g.is_file() or sha256_file(g) != r_["sha256"]:
            und.append(f"参照 {arm} の {r_['run']}/{r_['res']} が無いか sha256 が {DIAG_REF} の {r_['sha256'][:16]}… と違う")
            continue
        specs[arm] = {"run": str((case / r_["run"]).resolve()), "res": r_["res"], "sha256": r_["sha256"]}
        rdirs[arm] = case / r_["run"]
    R["refs"] = {k: {**v, "run": Path(v["run"]).name} for k, v in specs.items()}
    info = jload(run / "prepare_info.json") or {}
    eref = info.get("euler_ref") or {}
    if "aux" in specs and (Path(str(eref.get("run", ""))).name, eref.get("res"), eref.get("sha256")) != \
            (Path(specs["aux"]["run"]).name, specs["aux"]["res"], specs["aux"]["sha256"]):
        und.append(f"診断の NS の prepare_info の euler_ref {eref} が {DIAG_REF} の固定 (run_0950) と違う")
    return {"run": run, "rec": rec, "required": required, "ext": ext, "x_use": x_use, "mm": mm, "specs": specs, "rdirs": rdirs}


def diag_sens(case: Path, tools: Path, out_p: Path, jobs=None) -> int:
    """診断の感度の判定 (§6 V-h370-diag、CFD 0 step): 診断の NS (run_0951、印 H370_DIAG_ONLY.json) の停止の時刻 s (CENSORED は最後の出力) の
    [s − 40,000, s] の 2,500 ごとの全 17 出力 (列挙して完全性・sha の一覧と照合; 欠ければ判定不能、縮めない) に、Euler の参照の
    {run, res, sha256} だけを替えて (腕 A = run_0952 の res_40000・腕 B = run_0953 の res_40000、どちらも DIAG_REF の sha を照合) 抽出する
    (`diag_snapshot`、条件は停止規則の記録の extract = 受け渡しの固定の条件、比べる列 F は抽出の前に座標と両腕の Euler の補間できる範囲から)。
    D = max over (17 出力 × F) |δ_E,B − δ_E,A| [mm] (未緩和の δ_E) で `classify_diag`。δ_target・帯の位置・held・欠測・合否も記録する。補助:
    run_0950 (res_60000、診断の NS の参照) で抽出した値との差 (記録だけ; 主の判定に使わない)。受け渡しの関門・正式の判定は呼ばない。
    update_conditions_met は常に偽。終了コード 0 結果 A・結果 B / 2 判定不能。"""
    R = {"plan": PLAN_DIAG, "created": now(), "run": NS_RUN, "case_dir": str(case),
         "kind": "診断 (§6 V-h370-diag の感度の判定)。資格・受け渡し・壁の更新には使わない",
         "scope": ("この固定の系列に対する感度の結論で、NS の収束・Euler の適格性・異なる参照で選ばれる停止の時刻の同等性を意味しない。"
                   "結果 A でも残差の条件は外さない (§6 V-h370-diag のやらないこと)"),
         "update_conditions_met": False, "result": "判定不能", "code": None, "why": None, "undecidable": [], "tol_mm": TOL_MM}
    und = R["undecidable"]
    arrays = {}

    def finish(rc):
        R["update_conditions_met"] = False
        if arrays:
            npz = case / DIAG_NPZ
            npz.parent.mkdir(parents=True, exist_ok=True)
            np.savez_compressed(npz, **arrays)
            R["arrays"] = {"file": DIAG_NPZ, "sha256": sha256_file(npz), "keys": sorted(arrays)}
        if und:
            R["result"], R["code"] = "判定不能", R.get("code") or "U_input"
        if R["result"] == "判定不能" and not R.get("why"):
            R["why"] = "; ".join(und)[:900]
        jwrite(out_p, R)
        print(f"== VERDICT V-h370-diag (感度の診断、更新・受け渡しには使わない): {R['result']} ({R.get('code')}) — {str(R.get('why'))[:600]}")
        print(f"→ {out_p}")
        return 0 if R["result"] in ("結果 A", "結果 B") else 2
    inp = diag_inputs(case, R, und)                  # 入力の照合 1〜6 (帯を共通にした A/B (#36n) と共通)
    if inp is None or und:
        return finish(2)
    run, rec, required, ext, x_use, mm, specs, rdirs = (inp[k] for k in ("run", "rec", "required", "ext", "x_use", "mm", "specs", "rdirs"))
    # --- 7. 比べる列 F (抽出の前に座標と両腕の Euler の補間できる範囲から) ---
    try:
        fx = fixed_columns(x_use, {a: rdirs[a] for a in DIAG_ARMS}, edge_tol=(ext.get("remap_edge") == "missing"))
        fx_aux = fixed_columns(x_use, {"aux": rdirs["aux"]}, edge_tol=(ext.get("remap_edge") == "missing"))
    except Exception as e:  # noqa: BLE001
        und.append(f"比べる列を決められない ({type(e).__name__}: {e})"[:400])
        return finish(2)
    F = fx.pop("mask")
    F_aux = fx_aux.pop("mask")
    R["columns"] = {**fx, "aux_cover (記録)": {"n": int((F_aux & F).sum()), "excluded": fx_aux["excluded"]},
                    "note": "抽出の前に座標だけで決めた比べる列 (NS の使用域の列 ∩ 腕 A・B の Euler の補間できる範囲)。この中の held・欠測・非有限は合格にしない"}
    if not F.any():
        und.append("比べる列 F が空")
        return finish(2)
    # --- 8. 抽出 (17 出力 × 腕 A・B + 補助 17) ---
    snap = HOOKS.get("diag_snapshot", diag_snapshot)
    tasks, keys = [], []
    for k in required:
        for arm in DIAG_ARMS + ("aux",):
            tasks.append((run, k, specs[arm], ext, x_use))
            keys.append((arm, k))
    try:
        pj = snap_jobs(jobs)
    except ValueError as e:
        und.append(str(e))
        return finish(2)
    R["parallel_plan"] = pj
    prog_p = run / "h370_diag_sens_progress.log"

    def progress(m, n, el):
        if m == n or m % 3 == 0:
            with open(prog_p, "a") as fh:
                fh.write(f"{now()} 抽出 {m}/{n} (経過 {el:.0f} s、残りの見込み {el / m * (n - m):.0f} s、並列 {pj['jobs']})\n")
    results, pinfo = run_snaps(tasks, snap, pj["jobs"], progress=progress)
    got = {kk: (r_["ok"] if "ok" in r_ else {"error": r_["exc"][:300]}) for kk, r_ in zip(keys, results)}
    R["extractions"] = {"n": len(tasks), "parallel": pinfo}
    # --- 9. 出力ごとの比較・分類 ---
    per = [diag_pair(k, got[("A", k)], got[("B", k)], mm, x_use, F) for k in required]
    R["per_output"] = per
    c = classify_diag(per, required)
    R["result"], R["code"], R["why"], R["D_max"] = c["result"], c["code"], c["why"], c["D_max"]
    Ds = [r["D_mm"] for r in per if r.get("D_mm") is not None]
    R["D_mm"] = {"max": max(Ds) if Ds else None, "min": min(Ds) if Ds else None, "n_with_D": len(Ds), "n_outputs": len(per),
                 "n_over_tol": sum(1 for d in Ds if d > TOL_MM), "per_output": {r["step"]: r.get("D_mm") for r in per},
                 "note": "D の判定は全出力の最大で決める (窓の中で 0.01 mm をまたいでも時刻ごとに判定不能にしない)"}
    # --- 10. 補助: run_0950 (診断の NS の参照) で抽出した値との差・停止規則の記録の δ との照合 (記録だけ) ---
    rows = {r_["step"]: r_ for r_ in rec.get("rows") or [] if isinstance(r_, dict) and isinstance(r_.get("step"), int)}
    aux = []
    for k in required:
        d0 = got[("aux", k)]
        a_ = {"step": k}
        if "error" in d0:
            a_["error"] = d0["error"]
        else:
            e0 = _arr(d0, "delta_E")
            for arm in DIAG_ARMS:
                d_ = got[(arm, k)]
                if "error" in d_:
                    continue
                e1 = _arr(d_, "delta_E")
                m = F & np.isfinite(e0) & np.isfinite(e1)
                a_[f"max_abs_vs_{arm}_mm"] = float(np.max(np.abs(e1 - e0)[m]) * mm) if m.any() else None
            rd = (rows.get(k) or {}).get("delta")
            if isinstance(rd, dict) and len(rd.get("delta_E") or []) == len(e0):
                er = _arr(rd, "delta_E")
                m = np.isfinite(er) & np.isfinite(e0)
                a_["recorded_row_max_abs_diff_mm"] = float(np.max(np.abs(er - e0)[m]) * mm) if m.any() else None
                a_["recorded_row_res_sha_same"] = rd.get("res_sha256") == d0.get("res_sha256")
        aux.append(a_)
    R["aux_run_0950"] = {"per_output": aux, "note": "補助の記録だけ (主の判定に使わない)。recorded_row は停止規則の見張りが同じ参照で記録した δ との差"}
    # --- 11. 配列 (出力 × 腕 × 使用域の列) ---
    arrays.update(steps=np.asarray(required, dtype=np.int64), x_use=x_use, F=F)
    for arm in DIAG_ARMS + ("aux",):
        for q in ("delta_E", "delta_target", "held", "y_b", "y_star"):
            arrays[f"{arm}_{q}"] = np.array([(_arr(got[(arm, k)], q) if (len(_arr(got[(arm, k)], q)) == len(x_use))
                                              else np.full(len(x_use), np.nan)) for k in required], dtype=float)
    return finish(0)


# ===================================================================================================== 帯を共通にした A/B (§5.1 #36n)
@contextlib.contextmanager
def fixed_band(band: dict):
    """抽出器の帯の位置 y_b を固定する (§5.1 #36n; 抽出器の本体は変えない): band_select "edge" の帯の位置を決める
    `deltastar.edge_band_positions` の戻りの y_b だけを `band` の値に差し替える。抽出器はその y_b を各列の
    `band_local_deficit(y_b_fixed=...)` に渡す (deltastar.py の edge の経路そのまま)。梯子・縁 y* の推定・感度の欄 (y_b × 0.8・1.25)・検査・
    平滑化・緩和は抽出器がそのまま行う (case/70 deltaE_series.fix_band と同じ差し替え)。
    band = {"x": 列の x [r_t], "y_b": 帯の位置 [r_t], "fill_own": 真偽}: 抽出器が帯を決める列 (Euler の範囲の中の全列) の x を表の x と
    **完全一致**で引く (同じ NS の格子の同じ式の x なので丸めの許容は置かない)。表に無い列は fill_own なら抽出器が選んだ帯のまま、そうでなければ
    ValueError。表から引いた値が非有限なら ValueError (抽出の不成立)。
    yield する dict に入るもの: calls (差し替えが呼ばれた回数)、x・y_b (固定した帯)・y_b_natural (抽出器が選び直した帯 = 固定しなければ使った帯)、
    n_table・n_own・x_own_range (表から引いた列・自分の帯の列)、raw (抽出器の戻りの x・生の δ_r・hard_ok; `deltastar_loop` が呼ぶ抽出器を
    包んで値を写すだけ)。"""
    import forge_design.metrics.deltastar as DSm
    from forge_design.feedback import deltastar_loop as DLm
    tx = np.asarray(band["x"], dtype=float)
    ty = np.asarray(band["y_b"], dtype=float)
    if tx.ndim != 1 or tx.shape != ty.shape:
        raise ValueError(f"固定する帯の表の形が違う (x {tx.shape}・y_b {ty.shape})")
    if len(np.unique(tx)) != len(tx):
        raise ValueError("固定する帯の表の x に重複がある")
    o = np.argsort(tx, kind="stable")
    tx, ty = tx[o], ty[o]
    fill_own = bool(band.get("fill_own"))
    cap = {"calls": 0}
    orig_pos, orig_ext = DSm.edge_band_positions, DLm.deltastar_from_core_matched_euler

    def pos(x, *a, **k):
        out = dict(orig_pos(x, *a, **k))
        xx = np.asarray(x, dtype=float)
        nat = np.asarray(out["y_b"], dtype=float)
        if len(tx):
            j = np.clip(np.searchsorted(tx, xx), 0, len(tx) - 1)
            hit = tx[j] == xx
            yb = np.where(hit, ty[j], np.nan)
        else:
            hit, yb = np.zeros(len(xx), dtype=bool), np.full(len(xx), np.nan)
        if not hit.all() and not fill_own:
            raise ValueError(f"固定する帯の表に無い列 {int((~hit).sum())} (x {xx[~hit][:4].tolist()} …) — 帯を共通にできない")
        bad = hit & ~np.isfinite(yb)
        if bad.any():
            raise ValueError(f"固定する帯の値が非有限の列 {int(bad.sum())} (x {xx[bad][:4].tolist()} …) — 抽出の不成立")
        yb = np.where(hit, yb, nat)
        cap.update(calls=cap["calls"] + 1, x=xx.copy(), y_b=yb.copy(), y_b_natural=nat.copy(), n_table=int(hit.sum()), n_own=int((~hit).sum()),
                   x_own_range=[float(xx[~hit].min()), float(xx[~hit].max())] if (~hit).any() else None)
        out["y_b"] = yb
        return out

    def ext(*a, **k):
        d = orig_ext(*a, **k)
        cap["raw"] = {"x": np.asarray(d["x"], dtype=float).copy(), "delta_r_raw": np.asarray(d["delta_r_raw"], dtype=float).copy(),
                      "hard_ok": np.asarray(d["hard_ok"], dtype=bool).copy()}
        return d

    DSm.edge_band_positions, DLm.deltastar_from_core_matched_euler = pos, ext
    try:
        yield cap
    finally:
        DSm.edge_band_positions, DLm.deltastar_from_core_matched_euler = orig_pos, orig_ext


def diag_snapshot_fixed(run_dir, step: int, euler_ref, extract: dict, x_use, band: dict) -> dict:
    """帯を共通にした A/B の抽出 1 回 (§5.1 #36n): #36l の `diag_snapshot` (stop_rule.delta_snapshot と同じ抽出に帯の位置を足したもの) を
    `fixed_band(band)` の中で呼ぶ (帯の位置 y_b だけを固定; 抽出の条件・一定 x の断面・使用域の列への x の一致の対応は #36l と同じ)。
    戻り: diag_snapshot の量 (y_b は固定して使った帯 = 抽出の表 delta_r_edge.csv の値、y_star は縁の推定) に足すもの: y_b_natural (抽出器が
    選び直した帯)・delta_r_raw (生の δ_r) を使用域の列で、band_all (帯を決めた全列の x・固定した y_b・選び直した y_b; 腕 B に渡す)、
    band_fill (表から引いた列・自分の帯の列の数と範囲、抽出の表の帯の位置が固定した値と同じか = 差し替えが表まで届いたか)。"""
    with fixed_band(band) as cap:
        d = diag_snapshot(run_dir, step, euler_ref, extract, x_use)
    if cap["calls"] != 1 or "raw" not in cap:
        raise RuntimeError(f"帯の位置の固定が 1 回だけ効いていない (edge_band_positions の呼び出し {cap['calls']} 回; band_select が edge でない?)")
    xu = np.asarray(x_use, dtype=float)

    def on_use(xs, v):
        o = np.argsort(xs, kind="stable")
        xs, v = xs[o], v[o]
        j = np.clip(np.searchsorted(xs, xu), 0, len(xs) - 1)
        return np.where(xs[j] == xu, v[j], np.nan)
    yb_fix = on_use(cap["x"], cap["y_b"])
    d["y_b_natural"] = on_use(cap["x"], cap["y_b_natural"]).tolist()
    d["delta_r_raw"] = on_use(cap["raw"]["x"], cap["raw"]["delta_r_raw"]).tolist()
    d["band_all"] = {"x": cap["x"].tolist(), "y_b": cap["y_b"].tolist(), "y_b_natural": cap["y_b_natural"].tolist()}
    d["band_fill"] = {"n_table": cap["n_table"], "n_own": cap["n_own"], "x_own_range": cap["x_own_range"],
                      "edge_table_matches_fixed": bool(np.array_equal(np.asarray(d["y_b"], dtype=float), yb_fix, equal_nan=True))}
    return d


def classify_fixed(per: list, required, repro_ok: bool, tol: float = TOL_MM) -> dict:
    """帯を共通にした A/B の分類 (§6「V-h370-diag の帯を共通にした A/B の事前登録」のとおり; #36l の classify_diag の「合否の変化 → 結果 B」・
    「D の下限 → 結果 B」は登録に無いので使わない): 戻り {"result", "code", "why", "D_max"}。順に: 入るべき出力の集合が無い・欠ける →
    判定不能、腕 A の再現の不成立 → 判定不能、抽出の例外 (帯の固定の不成立を含む) → 判定不能、どちらかの腕の比べる列 F の中に held・欠測・
    非有限 (D_mm なし) → 判定不能、D_fixed = max over (出力 × F) |δ_E,B − δ_E,A| > tol → 結果 B、それ以外 (全入力が有効・D_fixed ≤ tol) → 結果 A。"""
    def res(result, code, why, dmax=None):
        return {"result": result, "code": code, "why": why, "D_max": dmax}
    if not required:
        return res("判定不能", "U_required", "判定の窓に入るべき出力の集合が無い")
    have = {r["step"] for r in per if "error_missing" not in r}
    lack = [k for k in required if k not in have]
    if lack:
        return res("判定不能", "U_incomplete", f"判定の窓に入るべき {len(required)} 出力のうち比べられない出力 {len(lack)}: {lack[:8]} (縮めない)")
    if not repro_ok:
        return res("判定不能", "U_repro", "腕 A の再現の不成立 (固定の帯で抽出した腕 A の δ_E が #36l の記録 A_delta_E と一致しない)")
    per = [r for r in per if r["step"] in set(required)]
    err = [r["step"] for r in per if any(k.startswith("error") for k in r)]
    if err:
        return res("判定不能", "U_error", f"抽出の例外・帯の固定の不成立・比較の失敗 (出力 {err[:8]})")
    part = [(r["D_partial_mm"], r["step"], r.get("x_at_D_valid")) for r in per if r.get("D_mm") is None and r.get("D_partial_mm") is not None]
    bad = [r["step"] for r in per if r.get("D_mm") is None]
    if bad:
        top_p = max(part, key=lambda v: v[0]) if part else None
        return res("判定不能", "U_invalid", f"比べる列 F の中に held・欠測・非有限がある出力 {bad[:8]} (固定の帯で欠測・非有限・抽出の不成立; "
                                            f"両腕とも有効な列だけの最大 {top_p[0] if top_p else None} mm は記録だけ)",
                   {"mm": top_p[0], "step": top_p[1], "x": top_p[2], "partial": True} if top_p else None)
    top = max(((r["D_mm"], r["step"], r.get("x_at_D")) for r in per), key=lambda v: v[0])
    dmax = {"mm": top[0], "step": top[1], "x": top[2]}
    n_over = sum(1 for r in per if r["D_mm"] > tol)
    if top[0] > tol:
        return res("結果 B", "B_D", f"D_fixed = max|δ_E,B − δ_E,A| {top[0]:.4g} mm > {tol} mm (出力 {top[1]}・x {top[2]}; tol を超える出力 "
                                    f"{n_over} / {len(per)}) — 帯の選び方だけを直しても解決しない", dmax)
    return res("結果 A", "A", f"全 {len(per)} 出力・比べる列 F の全列で両腕とも有効、D_fixed = max|δ_E,B − δ_E,A| {top[0]:.4g} mm ≤ {tol} mm "
                              f"(出力 {top[1]}・x {top[2]}) — 「参照の場の差だけで今回の予算の超過を説明できる」を退ける", dmax)


def section_max(x, Dkx, F, steps, sections=FIXED_SECTIONS) -> list:
    """区間ごとの最大: Dkx (出力 × 使用域の列) の比べる列 F の中・区間 [lo, hi) (hi None は出口まで) の最大とその出力・x。
    非有限は除く (区間に有限の値が無ければ max_mm None)。"""
    x, Dkx, F = np.asarray(x, dtype=float), np.asarray(Dkx, dtype=float), np.asarray(F, dtype=bool)
    out = []
    for lo, hi in sections:
        m = F & (x >= lo) & ((x < hi) if hi is not None else True)
        r = {"x_range": [lo, hi if hi is not None else (float(x[m].max()) if m.any() else None)], "n_cols": int(m.sum()), "max_mm": None}
        if m.any():
            sub = np.where(np.isfinite(Dkx[:, m]), Dkx[:, m], -np.inf)
            if np.isfinite(sub).any():
                i, j = np.unravel_index(int(np.argmax(sub)), sub.shape)
                r.update(max_mm=float(sub[i, j]), step=int(steps[i]), x=float(x[m][j]))
        out.append(r)
    return out


def diag_sens_fixed_band(case: Path, tools: Path, out_p: Path, jobs=None) -> int:
    """帯の位置を共通にした Euler 参照の入れ替えの A/B (§5.1 #36n・§6「V-h370-diag の帯を共通にした A/B の事前登録」、CFD 0 step):
    #36l の記録 (h370_diag_sens.json・h370_diag_sens_arrays.npz) の**腕 A の帯 y_b(k, x)** を両腕に当て (各 NS 出力 k にその出力の帯)、腕 A
    (run_0952 res_40000)・腕 B (run_0953 res_40000) で抽出し直す (`diag_snapshot_fixed`、run_snaps の子プロセス・並列)。
    1. 入力の照合は #36l と同じ (`diag_inputs`) に、#36l の記録との照合を足す: npz の sha、17 出力・使用域の列・比べる列 F (座標と両腕の
       Euler の範囲から決め直して npz の F と同じ、列の数も記録と同じ)・参照 A・B の sha・抽出の条件・版が #36l と同じ。
    2. 腕 A (17 抽出): 帯 = npz の A_y_b(k) (使用域の列)。使用域の外の列 (縮流部・入口の直管; #36l の記録は使用域だけなので npz に無い) は
       腕 A が同じ出力で選ぶ帯 (fill_own; #36l と同じ抽出器の選び方)。**腕 A の再現**: 固定の帯の δ_E が npz の A_delta_E と同じ有限の列で、
       差の最大 ≤ REPRO_TOL_RT [r_t] (最大差を記録)。不成立・腕 A の抽出の例外なら判定不能で、腕 B は回さない。
    3. 腕 B (17 抽出): 帯 = 腕 A が同じ出力で使った帯 (帯を決めた全列; 腕 B が帯を決める列が表に無ければ抽出の不成立)。
    4. D_fixed = max over (17 出力 × F) |δ_E,B − δ_E,A| [mm] (未緩和の δ_E)。分類は `classify_fixed` (登録どおり)。区間ごと (FIXED_SECTIONS) の
       D_fixed と #36l の D (帯を選び直したとき、npz の B_delta_E − A_delta_E) の最大を並べ、帯の移動 (腕 B が選び直した帯 − 共通の帯) も記録する。
       生の δ_r・未緩和の δ_E・δ_target・held・帯の位置は配列 (FIXED_NPZ)。
    受け渡しの関門・正式の判定・壁の更新は呼ばない。update_conditions_met は常に偽。終了コード 0 結果 A・結果 B / 2 判定不能。"""
    R = {"plan": PLAN_FIXED, "created": now(), "run": NS_RUN, "case_dir": str(case),
         "kind": "診断 (§6 V-h370-diag の帯を共通にした A/B; 原因を分ける)。資格・受け渡し・壁の更新には使わない",
         "scope": ("結果 A でも本番の抽出を固定の帯に変えない・Euler を適格としない (§6 の登録)。この固定の系列に対する診断で、NS の収束・"
                   "Euler の適格性を意味しない"),
         "update_conditions_met": False, "result": "判定不能", "code": None, "why": None, "undecidable": [], "tol_mm": TOL_MM,
         "repro_tol_rt": REPRO_TOL_RT, "source_36l": {"json": DIAG_OUT, "npz": DIAG_NPZ}}
    und = R["undecidable"]
    arrays = {}

    def finish(rc):
        R["update_conditions_met"] = False
        if arrays:
            npz = case / FIXED_NPZ
            npz.parent.mkdir(parents=True, exist_ok=True)
            np.savez_compressed(npz, **arrays)
            R["arrays"] = {"file": FIXED_NPZ, "sha256": sha256_file(npz), "keys": sorted(arrays)}
        if und:
            R["result"], R["code"] = "判定不能", R.get("code") or "U_input"
        if R["result"] == "判定不能" and not R.get("why"):
            R["why"] = "; ".join(und)[:900]
        jwrite(out_p, R)
        print(f"== VERDICT V-h370-diag 帯を共通にした A/B (診断、更新・受け渡しには使わない): {R['result']} ({R.get('code')}) — "
              f"{str(R.get('why'))[:600]}")
        print(f"→ {out_p}")
        return 0 if R["result"] in ("結果 A", "結果 B") else 2
    inp = diag_inputs(case, R, und)                  # 入力の照合 1〜6 (#36l と同じ)
    if inp is None or und:
        return finish(2)
    run, required, ext, x_use, mm, specs, rdirs = (inp[k] for k in ("run", "required", "ext", "x_use", "mm", "specs", "rdirs"))
    R["mm_per_rt"] = mm
    # --- 7. #36l の記録との照合 ---
    src = R["source_36l"]
    r36, npz_p = jload(case / DIAG_OUT), case / DIAG_NPZ
    if not r36 or not npz_p.is_file():
        und.append(f"#36l の記録 {DIAG_OUT} か配列 {DIAG_NPZ} が無い")
        return finish(2)
    sha_npz = sha256_file(npz_p)
    src.update(created=r36.get("created"), result=r36.get("result"), code=r36.get("code"), D_max=r36.get("D_max"), npz_sha256=sha_npz)
    if (r36.get("arrays") or {}).get("sha256") != sha_npz:
        und.append(f"#36l の配列 {DIAG_NPZ} の sha256 が記録 ({str((r36.get('arrays') or {}).get('sha256'))[:16]}…) と違う")
    with np.load(npz_p) as zf:
        z = {k: zf[k] for k in zf.files}
    keys = ("steps", "x_use", "F", "A_y_b", "A_delta_E", "A_delta_target", "A_held", "B_delta_E", "B_y_b")
    lack = [k for k in keys if k not in z]
    if lack:
        und.append(f"#36l の配列に {lack} が無い")
        return finish(2)
    n_req, n_x = len(required), len(x_use)
    if [int(v) for v in z["steps"]] != list(required) or r36.get("required") != list(required):
        und.append(f"#36l の 17 出力 (配列 {[int(v) for v in z['steps']][:3]}…・記録 {str(r36.get('required'))[:40]}…) が今の判定の窓 {required[:3]}… と違う")
    if not np.array_equal(z["x_use"], x_use):
        und.append("#36l の配列の使用域の列 x_use が run の使用域の列と違う")
    if any(z[k].shape != (n_req, n_x) for k in keys[3:]):
        und.append(f"#36l の配列の形 {[z[k].shape for k in keys[3:]]} が (出力 {n_req}, 列 {n_x}) と違う")
    for arm in DIAG_ARMS:
        a36, now_ = (r36.get("refs") or {}).get(arm) or {}, R["refs"].get(arm) or {}
        if (a36.get("run"), a36.get("res"), a36.get("sha256")) != (now_.get("run"), now_.get("res"), now_.get("sha256")):
            und.append(f"参照 {arm} の固定 {now_} が #36l の記録 {a36} と違う")
    if r36.get("extract") != ext:
        und.append(f"抽出の条件 {ext} が #36l の記録 {r36.get('extract')} と違う")
    if (r36.get("versions") or {}).get("now") != R["versions"]["now"]:
        und.append("抽出器・見張りの版が #36l の記録の版と違う")
    if und:
        return finish(2)
    # --- 8. 比べる列 F (座標と両腕の Euler の補間できる範囲から決め直し、#36l の F と照合) ---
    try:
        fx = fixed_columns(x_use, {a: rdirs[a] for a in DIAG_ARMS}, edge_tol=(ext.get("remap_edge") == "missing"))
    except Exception as e:  # noqa: BLE001
        und.append(f"比べる列を決められない ({type(e).__name__}: {e})"[:400])
        return finish(2)
    F = fx.pop("mask")
    F36 = np.asarray(z["F"], dtype=bool)
    same_F = bool(np.array_equal(F, F36)) and (r36.get("columns") or {}).get("n") == int(F.sum())
    R["columns"] = {**fx, "same_as_36l": same_F, "n_36l": (r36.get("columns") or {}).get("n"),
                    "note": "比べる列は #36l と同じ (座標と両腕の Euler の補間できる範囲から決め直して照合)。この中の held・欠測・非有限は判定不能"}
    if not same_F:
        und.append(f"比べる列 F ({int(F.sum())} 列) が #36l の F ({int(F36.sum())} 列・記録 {(r36.get('columns') or {}).get('n')}) と違う")
    if not F.any():
        und.append("比べる列 F が空")
    yA = np.asarray(z["A_y_b"], dtype=float)
    nf = int((~np.isfinite(yA[:, F])).sum())
    if nf:
        und.append(f"共通の帯 (#36l の腕 A の y_b) が比べる列の中で非有限 {nf} 点 (欠測)")
    if und:
        return finish(2)
    # --- 9. 腕 A (17 抽出): 帯 = npz の A_y_b(k)、使用域の外は腕 A が選ぶ帯 ---
    snap = HOOKS.get("diag_snapshot_fixed", diag_snapshot_fixed)
    try:
        pj = snap_jobs(jobs)
    except ValueError as e:
        und.append(str(e))
        return finish(2)
    R["parallel_plan"] = pj
    prog_p = out_p.with_name(out_p.stem + "_progress.log")

    def progress(tag):
        def f(m, n, el):
            if m == n or m % 3 == 0:
                with open(prog_p, "a") as fh:
                    fh.write(f"{now()} 腕 {tag} の抽出 {m}/{n} (経過 {el:.0f} s、残りの見込み {el / m * (n - m):.0f} s、並列 {pj['jobs']})\n")
        return f
    tA = [(run, k, specs["A"], ext, x_use, {"x": x_use, "y_b": yA[i], "fill_own": True, "source": "#36l の npz の A_y_b (使用域の列)"})
          for i, k in enumerate(required)]
    resA, infoA = run_snaps(tA, snap, pj["jobs"], progress=progress("A"))
    gotA = {k: (r_["ok"] if "ok" in r_ else {"error": r_["exc"][:300]}) for k, r_ in zip(required, resA)}
    R["extractions"] = {"A": infoA}
    # 腕 A の再現 (固定の帯の δ_E が #36l の A_delta_E と一致する)
    rep = []
    for i, k in enumerate(required):
        d = gotA[k]
        if "error" in d:
            rep.append({"step": k, "ok": False, "error": d["error"]})
            continue
        e, e0 = _arr(d, "delta_E"), np.asarray(z["A_delta_E"][i], dtype=float)
        if len(e) != len(e0):
            rep.append({"step": k, "ok": False, "error": f"列の数 {len(e)} が #36l の {len(e0)} と違う"})
            continue
        fe, f0 = np.isfinite(e), np.isfinite(e0)
        m = fe & f0
        dmax = float(np.max(np.abs(e - e0)[m])) if m.any() else None
        t, t0 = _arr(d, "delta_target"), np.asarray(z["A_delta_target"][i], dtype=float)
        mt = np.isfinite(t) & np.isfinite(t0) if len(t) == len(t0) else np.zeros(0, dtype=bool)
        h, h0 = np.asarray(d.get("held", []), dtype=float), np.asarray(z["A_held"][i], dtype=float)
        ybn, yb = _arr(d, "y_b_natural"), _arr(d, "y_b")
        mb = np.isfinite(ybn) & np.isfinite(yA[i]) if len(ybn) == n_x else np.zeros(0, dtype=bool)
        bf = d.get("band_fill") or {}
        rep.append({"step": k, "ok": bool(np.array_equal(fe, f0) and dmax is not None and dmax <= REPRO_TOL_RT and bf.get("edge_table_matches_fixed")),
                    "finite_pattern_same": bool(np.array_equal(fe, f0)), "max_abs_delta_E_rt": dmax,
                    "max_abs_delta_E_mm": None if dmax is None else dmax * mm, "bitwise_delta_E": bool(np.array_equal(e, e0, equal_nan=True)),
                    "max_abs_delta_target_rt": float(np.max(np.abs(t - t0)[mt])) if mt.any() else None,
                    "held_same": bool(len(h) == len(h0) and np.array_equal(h, h0)),
                    "band_used_equals_npz": bool(len(yb) == n_x and np.array_equal(yb, yA[i], equal_nan=True)),
                    "band_natural_vs_npz_max_abs_rt": float(np.max(np.abs(ybn - yA[i])[mb])) if mb.any() else None,
                    "band_natural_bitwise_npz": bool(len(ybn) == n_x and np.array_equal(ybn, yA[i], equal_nan=True)), "band_fill": bf})
    dm = [r_["max_abs_delta_E_rt"] for r_ in rep if r_.get("max_abs_delta_E_rt") is not None]
    repro_ok = len(rep) == n_req and all(r_["ok"] for r_ in rep)
    R["repro_A"] = {"ok": repro_ok, "tol_rt": REPRO_TOL_RT, "max_abs_delta_E_rt": max(dm) if dm else None,
                    "max_abs_delta_E_mm": max(dm) * mm if dm else None, "bitwise_all": all(r_.get("bitwise_delta_E") for r_ in rep),
                    "band_natural_bitwise_all": all(r_.get("band_natural_bitwise_npz") for r_ in rep), "per_output": rep,
                    "note": ("腕 A を固定の帯 (#36l の npz の A_y_b、使用域の外は腕 A が選ぶ帯) で抽出した δ_E と #36l の A_delta_E の照合。"
                             "band_natural_* は腕 A が選び直した帯と npz の帯の照合 (記録)")}

    def stack(got, q):
        return np.array([(_arr(got[k], q) if len(_arr(got[k], q)) == n_x else np.full(n_x, np.nan)) for k in required], dtype=float)
    arrays.update(steps=np.asarray(required, dtype=np.int64), x_use=x_use, F=F)
    qs = ("delta_E", "delta_target", "held", "y_b", "y_b_natural", "y_star", "delta_r_raw")
    for q in qs:
        arrays[f"A_{q}"] = stack(gotA, q)
    if not repro_ok:
        errA = [r_["step"] for r_ in rep if "error" in r_]
        R["code"] = "U_error" if errA else "U_repro"
        und.append(f"腕 A の抽出の例外 (出力 {errA[:8]})" if errA else
                   f"腕 A の再現の不成立: 固定の帯の δ_E が #36l の A_delta_E と一致しない (最大差 {R['repro_A']['max_abs_delta_E_rt']} r_t、"
                   f"許容 {REPRO_TOL_RT}; 不成立の出力 {[r_['step'] for r_ in rep if not r_['ok']][:8]}、うち有限の列が違う出力 "
                   f"{[r_['step'] for r_ in rep if not r_.get('finite_pattern_same', True)][:8]}) — 腕 B は回さない")
        return finish(2)
    # --- 10. 腕 B (17 抽出): 帯 = 腕 A が同じ出力で使った帯 (帯を決めた全列) ---
    tB = [(run, k, specs["B"], ext, x_use, {"x": np.asarray(gotA[k]["band_all"]["x"], dtype=float), "y_b": np.asarray(gotA[k]["band_all"]["y_b"], dtype=float),
                                             "fill_own": False, "source": "腕 A が同じ出力で使った帯 (この判定の腕 A の抽出の band_all)"})
          for k in required]
    resB, infoB = run_snaps(tB, snap, pj["jobs"], progress=progress("B"))
    gotB = {k: (r_["ok"] if "ok" in r_ else {"error": r_["exc"][:300]}) for k, r_ in zip(required, resB)}
    R["extractions"].update(B=infoB, n=len(tA) + len(tB))
    for q in qs:
        arrays[f"B_{q}"] = stack(gotB, q)
    # --- 11. 出力ごとの比較・分類 ---
    per = []
    D36 = np.abs(np.asarray(z["B_delta_E"], dtype=float) - np.asarray(z["A_delta_E"], dtype=float)) * mm
    for i, k in enumerate(required):
        dA, dB = gotA[k], gotB[k]
        r_ = diag_pair(k, dA, dB, mm, x_use, F)
        if not any(kk.startswith("error") for kk in r_):
            for arm, d in (("A", dA), ("B", dB)):
                if not (d.get("band_fill") or {}).get("edge_table_matches_fixed"):
                    r_[f"error_band_{arm}"] = "抽出の表の帯の位置が固定した値と違う (帯の固定が効いていない)"
            ya, yb_ = _arr(dA, "y_b"), _arr(dB, "y_b")
            r_["band_common"] = bool(len(ya) == len(yb_) == n_x and np.array_equal(ya, yb_, equal_nan=True))
            if not r_["band_common"]:
                r_["error_band_common"] = "両腕の帯の位置 (使用域の列) が同じでない"
            ybn = _arr(dB, "y_b_natural")
            mv = F & np.isfinite(ybn) & np.isfinite(ya) if len(ybn) == len(ya) == n_x else np.zeros(n_x, dtype=bool)
            if mv.any():
                dmv = np.abs(ybn - ya) * mm
                j = int(np.flatnonzero(mv)[int(np.argmax(dmv[mv]))])
                r_["band_move_B_max_mm"], r_["x_at_band_move_B"] = float(dmv[j]), float(x_use[j])
            r_["band_fill_B"] = dB.get("band_fill")
        m36 = F & np.isfinite(D36[i])
        if m36.any():
            j = int(np.flatnonzero(m36)[int(np.argmax(D36[i][m36]))])
            r_["D_36l_mm"], r_["x_at_D_36l"] = float(D36[i][j]), float(x_use[j])
        per.append(r_)
    R["per_output"] = per
    c = classify_fixed(per, required, repro_ok)
    R["result"], R["code"], R["why"], R["D_max"] = c["result"], c["code"], c["why"], c["D_max"]
    Ds = [r_["D_mm"] for r_ in per if r_.get("D_mm") is not None]
    R["D_fixed_mm"] = {"max": max(Ds) if Ds else None, "min": min(Ds) if Ds else None, "n_with_D": len(Ds), "n_outputs": len(per),
                       "n_over_tol": sum(1 for d in Ds if d > TOL_MM), "per_output": {r_["step"]: r_.get("D_mm") for r_ in per},
                       "D_36l_per_output": {r_["step"]: r_.get("D_36l_mm") for r_ in per},
                       "note": "D_fixed は帯を共通にしたとき、D_36l は #36l (帯を参照ごとに選び直したとき) の同じ出力の最大"}
    # --- 12. 区間ごとの最大 (D_fixed と #36l の D を並べる)・帯の移動・列ごとの最大 ---
    Dfix = np.abs(arrays["B_delta_E"] - arrays["A_delta_E"]) * mm
    move = np.abs(arrays["B_y_b_natural"] - arrays["A_y_b"]) * mm
    s_fix, s_36, s_mv = (section_max(x_use, a_, F, required) for a_ in (Dfix, D36, move))
    R["sections"] = [{"x_range": a_["x_range"], "n_cols": a_["n_cols"], "D_fixed": {k: a_.get(k) for k in ("max_mm", "step", "x")},
                      "D_36l": {k: b_.get(k) for k in ("max_mm", "step", "x")}, "band_move_B": {k: m_.get(k) for k in ("max_mm", "step", "x")}}
                     for a_, b_, m_ in zip(s_fix, s_36, s_mv)]
    R["sections_note"] = ("区間 [lo, hi) [r_t] (最後は出口まで) の比べる列 F の中の 17 出力の最大 [mm]。D_fixed = 帯を共通にした |δ_E,B − δ_E,A|、"
                          "D_36l = #36l (帯を選び直した) の同じ量、band_move_B = 腕 B が選び直した帯 − 共通の帯 (腕 A の帯) の |差| [mm]")

    def colmax(a_):
        a_ = np.where(np.isfinite(a_), a_, -np.inf).max(axis=0)
        return np.where(np.isfinite(a_), a_, np.nan)
    arrays.update(D_fixed_kx_mm=Dfix, D_36l_kx_mm=D36, band_move_B_kx_mm=move, D_fixed_x_mm=colmax(Dfix), D_36l_x_mm=colmax(D36),
                  band_move_B_x_mm=colmax(move))
    o36, ofx = F & (arrays["D_36l_x_mm"] > TOL_MM), F & (arrays["D_fixed_x_mm"] > TOL_MM)
    R["over_tol_columns"] = {k: {"n": int(m_.sum()), "x_range": [float(x_use[m_].min()), float(x_use[m_].max())] if m_.any() else None}
                             for k, m_ in (("D_36l", o36), ("D_fixed", ofx), ("D_36l_only", o36 & ~ofx), ("both", o36 & ofx))}
    xs0 = [np.asarray(gotA[k]["band_all"]["x"], dtype=float) for k in required]
    if all(np.array_equal(xs0[0], v) for v in xs0):
        arrays["band_all_x"] = xs0[0]
        arrays["A_band_all_y_b"] = np.array([gotA[k]["band_all"]["y_b"] for k in required], dtype=float)
        arrays["A_band_all_y_b_natural"] = np.array([gotA[k]["band_all"]["y_b_natural"] for k in required], dtype=float)
        if all("band_all" in gotB[k] and np.array_equal(np.asarray(gotB[k]["band_all"]["x"], dtype=float), xs0[0]) for k in required):
            arrays["B_band_all_y_b"] = np.array([gotB[k]["band_all"]["y_b"] for k in required], dtype=float)
            arrays["B_band_all_y_b_natural"] = np.array([gotB[k]["band_all"]["y_b_natural"] for k in required], dtype=float)
    else:
        R["band_all_note"] = "帯を決めた全列の x が出力で違う (全列の帯の配列は残さない)"
    return finish(0)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--euler-prequalify", action="store_true", help="(i) Euler の参照の予備の資格の判定と固定")
    ap.add_argument("--diag-sens", action="store_true",
                    help="§6 V-h370-diag: 診断の NS (run_0951、ns-diag) の Euler 参照の入れ替えの感度 (CFD 0 step; 関門・正式の判定は呼ばない)")
    ap.add_argument("--diag-sens-fixed-band", action="store_true",
                    help="§5.1 #36n: --diag-sens の記録の腕 A の帯 y_b(k, x) を両腕に当てた参照の入れ替えの A/B (CFD 0 step; 関門・正式の判定は呼ばない)")
    ap.add_argument("--no-ab", action="store_true",
                    help="抽出の参照の A/B と σ_ref を回さない (A/B は未実施、参照の最終の資格を確かめないので更新に進める条件は満たさない)")
    ap.add_argument("--case-dir", default=str(HERE))
    ap.add_argument("--design", default=None, help="forge_design の design/ (既定: ~/forge-g0cchain の、無ければこの作業ツリーの)")
    ap.add_argument("--tools", default=None, help="forge の tools (既定: ~/forge-g0cchain の、無ければこの作業ツリーの)")
    ap.add_argument("--out", default=None)
    ap.add_argument("--jobs", type=int, default=None,
                    help=f"抽出の並列数の上限 (1〜{JOBS_MAX}、既定 {JOBS_CAP}; 実際は min(上限, 空きメモリ / {SNAP_PEAK_MIB} MiB))")
    a = ap.parse_args(argv)
    case = Path(a.case_dir).resolve()
    design = Path(a.design) if a.design else Path.home() / "forge-g0cchain/design"
    if not design.is_dir():
        design = HERE.parents[1] / "design"
    tools = Path(a.tools) if a.tools else Path.home() / "forge-g0cchain/solver_density_cuda/tools"
    if not tools.is_dir():
        tools = HERE.parents[1] / "solver_density_cuda/tools"
    if str(design) not in sys.path:
        sys.path.insert(0, str(design))
    if a.euler_prequalify:
        return euler_prequalify(case, tools)
    if a.diag_sens:
        return diag_sens(case, tools, Path(a.out) if a.out else case / DIAG_OUT, jobs=a.jobs)
    if a.diag_sens_fixed_band:
        return diag_sens_fixed_band(case, tools, Path(a.out) if a.out else case / FIXED_OUT, jobs=a.jobs)
    return ns_judge(case, tools, not a.no_ab, Path(a.out) if a.out else case / OUT, jobs=a.jobs)


if __name__ == "__main__":
    sys.exit(main())
