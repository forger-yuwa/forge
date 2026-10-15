r"""NS の本段の停止規則 (plan tooling-nozzle-cooled-wall-ns-settling §4.6・plan tooling-nozzle-isothermal-wall-chain §4.9・§5.1 #36a・#36e、2026-10-14)。

**由来と統合** (2026-10-14、main への統合 1 回目で 2 つの版を 1 つにした):
forge-f5 の feature/nozzle-cooled-wall-settling (66364fff に θ_r の断面の基準を設計壁の出口に直した 60de8125 を重ねた版、plan
tooling-nozzle-cooled-wall-ns-settling §4.6、ユーザ決定・諮問 stop-rule を全件採用) を、M6 の feature/nozzle-wall-fit-and-pipeline が
写して変えた版 (plan tooling-nozzle-isothermal-wall-chain §5.1 #36a・#36e) を土台にし、使うかどうかの既定を forge-f5 の版に戻した。
  - **既定で使う** (forge-f5、ユーザ決定): 問題 YAML に `evaluate.stop_rule` が無ければ使う。使えない (euler_ref が無い・見つからない・
    出力の間隔・上限が出力の間隔の倍数でない・設計壁の出口が無い) ときは forge-f5 の版の動き = UNAVAILABLE を記録し、本段を `fallback_steps`
    (`--steps` > `evaluate.ns_steps` > 60,000) の長さで回す。
  - **`evaluate.stop_rule: true` を明示** (または関数の引数 stop_rule=True) したときだけ、使えない設定を `prepare_ns` が run dir を作る前に・
    `run_staged_ns` が forge を起動する前に例外にする (M6 の版の動き; 記録 `explicit`)。`evaluate.stop_rule: false`・CLI `--no-stop-rule` は使わない。
  - 到達の後の継続 `evaluate.stop_rule_hold_steps` (M6 の版で足した) の**既定は `HOLD_STEPS_DEFAULT` = 40,000 step** (ユーザ決定、forge-f5 経由;
    値はこの 1 か所の定数)。最初に条件を満たした出力 s₀ の後も回し、s₀ + H までの全ての出力で条件を保ったら REACHED (止める出力は s₀ + H
    以後の最初の出力)。途中で外れたら到達を取り消して `hold.breaks` に記録し、次の到達から数え直す。上限 (`cap`) は到達 + 保持を含む本段の
    step 数で、保持の途中で上限に届いたら CENSORED。0 を書けば最初の到達で止める (forge-f5 の版の動き)。
  - M6 の版で足したものは残す: 終了の状態の記録 `STOP_RULE_STATE.json` (`write_state`・`read_state`; `stop_rule.json` は見張りの途中も
    書き換わる詳しい記録、こちらは本段が終わった後に 1 度だけ書く状態の印)、問題 YAML の `evaluate.stop_rule_thresholds` (`thresholds` と
    同じキー; 例 case/45 の試験部の断面 x = 40・70・94 r_t は `x_fracs` に x / x_w を書く。既定は forge-f5 の版のまま)、上限の倍数の検査、
    (§5.1 #36e) 帯の位置の抽出を受け渡しの抽出と同じ条件 (`deltastar_loop.handoff_extract_settings`: 一定 x の断面への補間
    `remap_constant_x=True`・端の補間範囲の不足は欠測 `remap_edge="missing"`) で行うこと、Euler の参照の固定 {run, res, sha256}
    (`deltastar.resolve_euler_ref`、sha が合わなければ使わない)、(§5.1 #36e (7)) 見張りの中で出力ごとに δ_E・δ_target (抽出器の量、使用域
    [0, x_w] の全列) を受け渡しと同じ条件で抽出して `stop_rule.json` の行の `delta` に残すこと (**判定には使わない記録の量**。受け渡しの
    関門の 2 窓の判定がこの記録を読み、計算の後に抽出し直さない)。抽出は別スレッド 1 本で行い (`DeltaRecorder`)、停止の判定を待たせない。
    追いつかないときは最新の出力から 2 窓 (2W) より古い出力を飛ばす (`skipped`; 停止の時刻 s ≥ 最新の出力なので、[s − 2W, s] に入らない
    出力だけ)。要求の断面の位置 (x_fracs × x_w) は `sections` に残す。
  - θ_r の断面の基準は設計壁の出口 `physical_wall.domain[1]` (60de8125 の修正; 軸 Mach 則の終点 x_E ではない)。判定・量・見張り・forge の
    止め方は forge-f5 の版のまま。
  - (§5.1 #36g) 帯の生成と δ の記録は抽出器 (remap_edge="missing") を通るので、Euler の参照の補間できる範囲は共通の
    `deltastar.euler_cover_range` (端の外側だけに幾何の許容 = 端の列の間隔 × 1e−4、許容の内側で端を越える列は端の列の値) になる。
    2 窓に入るべき出力は `required_window_steps` (出力の間隔から列挙; 記録にある行だけから作らない) で決める。

**新しい規則であり、旧 `m9_watch.py` (§6.11、`case/45.isobutane_m6_d155/m9_watch.py`) と同等とは言わない**: m9_watch は窓の両端の差で傾向と振動を
混ぜて判定し、断面・量・連続の到達の条件もここと違う (codex diagnose 2026-10-14、`notes/reviews/2026-10-14-half-inlet-370k-design-diagnose.md`)。
ここは M6 の考え方を移したうえで、**平均の変化と振幅を分けて**判定する別の規則として登録する。m9_watch の到達の step・結果をこの規則の結果として読み替えない。
REACHED は「この規則の量が停止の条件を満たした」の意味に限り、収束・設計の量の定常・壁の更新に進む条件 (§4.9) の合格ではない。

判定量 (出力 res_<N>.h5 ごと、N > 0):
  - 質量の収支 $r = (\dot m_{in} - \dot m_{out})/\dot m_{in}$: 入口 (列 i = 0)・出口 (列 i = ni − 1) の列の折れ線を通る流束
    $2\pi\int(\rho u_x\,dr - \rho u_r\,dx)\,r$ (M6 の `cold_xcheck.reduce_fields` と同じ式)。
  - 境界層の運動量厚さ $\theta_r$ を $x/x_w$ = 0.4・0.7・0.95 の断面で ($x_w$ は設計壁 (物理壁) の出口 `prepare_info.physical_wall.domain[1]`、r_t 単位;
    `prepare_info.x_E` は軸 Mach 則の終点で別物 — 2026-10-14 まで誤って x_E を使っていた、M6ノズルの設計・forge-25 の指摘)。各行 j を x 方向に補間した
    一定 x の断面で、帯の外縁 $y_b$ の $(\rho_e, u_e)$ を縁にした面積で等価な厚さ
    $\theta_r = r_w - \sqrt{r_w^2 - 2\int_{r_b}^{r_w}\rho u(u_e-u)\,r\,dr/(\rho_e u_e^2)}$ (`reduce_fields` と同じ式、断面は最も近い列)。
  - $y_b(x)$ は**本段の最初の出力**を固定 Euler 参照 (`prepare_info.euler_ref`) で `extract_and_merge(.., band_select="edge")` と同じ抽出にかけた
    帯の位置 (縁 y* の 1.25 倍、さらに掛けない; 諮問 stop-rule (a)) を本段の間は固定する (`fix_band`)。
  - 記録だけ (判定に使わない; 諮問 (d)): 出口コア M (最後の列の η = r/r_w ∈ [0.05, 0.7] の M の平均) と無次元の流量
    (2π∫ρU_x r dr の x ∈ (−2.5, x_E) の中央値) — `report/nozzle_report.py` の `metrics` と `case/70.m5_air_wt_demo/qs_series.py` と同じ定義。
    非有限の数 (VALUE の全ての浮動小数点の量)。

判定 (`judge_rows`、閾値・窓・断面は引数、既定は §4.6):
  - 窓 = 直近 W = 2 万 step の出力 (step ∈ [s − W, s]、出力の間隔に依らない)。本段の履歴が実際に W を覆い (最初の出力 ≤ s − W)、窓に 9 点以上。
  - r: |窓の平均| ≤ 0.05 %・振幅 (最大 − 最小)/2 ≤ 0.5 %。
  - θ_r (各断面): $D_\theta$ = |後半の平均 − 前半の平均| / |窓の平均| × W / (後半の平均の step − 前半の平均の step) ≤ 0.05 % (2 万 step あたりに換算;
    諮問の Major)、振幅 (最大 − 最小)/2 / |窓の平均| ≤ 0.5 %。前半・後半は窓の中点 s − W/2 で分け、**中点ちょうどの出力は後半に入れる**
    (点数が奇数で等間隔のとき。諮問の計算と同じ割り方。換算するので線形の系列では割り方に依らない)。換算しない値も `D_raw` として記録する。
  - 満たした最初の出力で REACHED (連続は要らない; 諮問 (b))。ただし hold_steps (既定 `HOLD_STEPS_DEFAULT`) > 0 のときは上の「到達の後の継続」。
    非有限で DIVERGED。REACHED は「指定の量の停止の条件を満たした」の意味に限り、収束・設計の量の定常の VERDICT とは別。

状態: REACHED / CENSORED (上限まで正常に回って届かない) / DIVERGED (非有限、または forge の detectNaN) / EXEC_ERROR (forge が上限より前に異常終了・
見張りの例外) / DATA_ERROR (出力が読めない・判定量が有限でない) / UNAVAILABLE (停止規則を使えない: euler_ref が無い・抽出できない・出力の間隔が
窓に 9 点を入れられない — 本段は ns_steps の長さで回す)。記録は `<run>/stop_rule.json` (一時ファイル → rename で原子的に書く)。
明示 (`evaluate.stop_rule: true`) の run は準備・起動の時点で使えないと分かったら例外で止めるので、UNAVAILABLE は見張りの途中で帯の位置を
決められないとき (本段は prepare の nStepOuter = fallback_steps の出力で止める) だけ出る。既定 (キー無し) の run は準備・起動の時点で使えなくても
UNAVAILABLE を記録して fallback_steps の長さで回す。本段が終わったら `run_staged_ns` が終了の状態を `<run>/STOP_RULE_STATE.json` に
1 度だけ書く (`write_state`)。`run_staged_ns` の終了コード (REACHED・CENSORED・一部の UNAVAILABLE で 0) で合否を決めず、この状態を直接見る
(`feedback.deltastar_loop.handoff_gate`)。

見張り (`watch_main_stage`): `run_case.sh` を子プロセス (新しいプロセスグループ) で持ったまま、書き終わった res を判定し、REACHED・DIVERGED で forge に
SIGTERM → 待って SIGKILL (止まったことを /proc で確かめる)、止めた出力より後の res を消す。forge の終了コードは `run_case.sh` の
`[run_case] forge exit=N` の行から読む (`FORGE_STRICT_CONVERGENCE` の終了コードと取り違えない)。forge だけを止めるので `run_case.sh` の後処理
(収束判定・残差の png) はそのまま走る; 走らなかったら同じ後処理を呼ぶ。
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import signal
import subprocess
import sys
import tempfile
import time
from datetime import datetime
from pathlib import Path

import h5py
import numpy as np

# §4.6 の既定 (閾値は M6 の実績を移しただけ; 風洞の形状での妥当性は §6 S-stop (c) で見る)
WINDOW = 20000
MIN_POINTS = 9
R_MEAN_MAX = 5e-4          # |r の窓の平均| ≤ 0.05 %
R_AMP_MAX = 5e-3           # r の振幅 ≤ 0.5 %
THETA_DRIFT_MAX = 5e-4     # D_θ ≤ 0.05 % / 2 万 step
THETA_AMP_MAX = 5e-3       # θ_r の振幅 / |平均| ≤ 0.5 %
X_FRACS = (0.4, 0.7, 0.95)
CAP_DEFAULT = 200_000      # 停止規則を使うときの本段の上限 (諮問 (c)、M6 と同じ)。到達 + 保持を含む本段の step 数
HOLD_STEPS_DEFAULT = 40_000   # 到達の後の継続の既定 (2026-10-14 ユーザ決定、forge-f5 経由)。evaluate.stop_rule_hold_steps で上書き (0 = 最初の到達で止める)
STATES = ("REACHED", "CENSORED", "DIVERGED", "EXEC_ERROR", "DATA_ERROR", "UNAVAILABLE")
STATE = "STOP_RULE_STATE.json"      # 本段の終了の状態 (run_staged_ns が本段の後に 1 度だけ書く; M6 の版で足した)
PLAN = ("plans/active/tooling-nozzle-cooled-wall-ns-settling.md §4.6・plans/active/tooling-nozzle-isothermal-wall-chain.md §4.9・§5.1 #36a "
        "(2026-10-14 統合: 既定で使う・明示のときだけ使えない設定を例外・保持の既定 40,000)")
RECORD = "stop_rule.json"
BAND_DIR = "stop_rule_band"
_RES_RE = re.compile(r"res_(\d+)\.h5")
_OUT_RE = re.compile(r"res_(?:(?!nan_).+_)?(\d+)\.(?:h5|xmf)")   # res_<N>.h5/.xmf と境界の res_<名前>_<id>_<N>.h5/.xmf (res_nan_* は除く)
_FORGE_EXIT_RE = re.compile(r"\[run_case\] forge exit=(-?\d+)")


class Unavailable(Exception):
    """停止規則を使えない (例外として外へ出さず、UNAVAILABLE として記録する)。"""


def thresholds(**kw) -> dict:
    """判定の閾値・窓・断面 (既定は §4.6)。未知のキーは拒否する。値の検査 (M6 の版で足した、問題 YAML の `evaluate.stop_rule_thresholds` も
    ここを通る): window・min_points は正の整数 (min_points ≥ 2)、4 つの閾値は有限の正の数、x_fracs は (0, 1] の有限の数の空でない
    重複の無いリスト (設計壁の出口 x_w に対する比; 断面は x_fracs × x_w に最も近い列)。"""
    d = dict(window=WINDOW, min_points=MIN_POINTS, r_mean_max=R_MEAN_MAX, r_amp_max=R_AMP_MAX,
             theta_drift_max=THETA_DRIFT_MAX, theta_amp_max=THETA_AMP_MAX, x_fracs=list(X_FRACS))
    bad = set(kw) - set(d)
    if bad:
        raise ValueError(f"停止規則の未知の設定 {sorted(bad)}")
    d.update(kw)
    for k in ("window", "min_points"):
        v = d[k]
        if isinstance(v, bool) or not isinstance(v, int) or v <= 0 or (k == "min_points" and v < 2):
            raise ValueError(f"停止規則の {k} は正の整数{' (≥ 2)' if k == 'min_points' else ''} ({v!r})")
    for k in ("r_mean_max", "r_amp_max", "theta_drift_max", "theta_amp_max"):
        v = d[k]
        if isinstance(v, bool) or not isinstance(v, (int, float)) or not np.isfinite(v) or v <= 0:
            raise ValueError(f"停止規則の {k} は有限の正の数 ({v!r})")
    xf = d["x_fracs"]
    if (not isinstance(xf, (list, tuple)) or not xf
            or any(isinstance(v, bool) or not isinstance(v, (int, float)) or not np.isfinite(v) or not 0.0 < v <= 1.0 for v in xf)
            or len({float(v) for v in xf}) != len(xf) or len({theta_key(v) for v in xf}) != len(xf)):
        raise ValueError(f"停止規則の x_fracs は (0, 1] の数の空でない重複の無いリスト ({xf!r})")
    d["x_fracs"] = [float(v) for v in xf]
    return d


def theta_key(frac: float) -> str:
    return f"theta_r_{float(frac):g}"


def required_window_steps(step: int, window: int, out_int: int, n_windows: int = 2) -> list:
    """停止の出力 s で終わる n_windows 窓 [s − n W, s] (閉区間) に入るべき出力の step を、出力の間隔から列挙する (昇順)。
    plan tooling-nozzle-isothermal-wall-chain §5.1 #36g (codex diagnose 2026-10-14 の Major「判定窓の完全性」): 記録 (stop_rule.json の行) に
    ある出力だけから窓を作ると、行そのものの欠けを検出できない。forge は outStepInterval の倍数で res を書き、s もその倍数なので、
    必要な出力は s − k·Δ (k = 0..⌊n W/Δ⌋)。V-h370 (W 20,000・Δ 2,500・2 窓) では 17 出力。引数が整数でなければ ValueError。"""
    vals = {"step": step, "window": window, "out_int": out_int, "n_windows": n_windows}
    for k, v in vals.items():
        if isinstance(v, bool) or not isinstance(v, (int, np.integer)) or (v < 0 if k == "step" else v <= 0):
            raise ValueError(f"required_window_steps: {k} は {'0 以上' if k == 'step' else '正'}の整数 ({v!r})")
    s, n, d = int(step), int(n_windows) * int(window), int(out_int)
    return sorted(s - k * d for k in range(n // d + 1))


# ------------------------------------------------------------------ 設定 (prepare_ns が prepare_info.json の stop_rule に書く)
EVALUATE_KEYS = ("stop_rule", "stop_rule_cap", "stop_rule_thresholds", "stop_rule_hold_steps")


def check_evaluate_keys(evaluate: dict) -> tuple:
    """問題 YAML の `evaluate.stop_rule` (true / false)・`evaluate.stop_rule_cap` (正の整数)・`evaluate.stop_rule_thresholds`
    (`thresholds` のキーの mapping)・`evaluate.stop_rule_hold_steps` (0 以上の整数) を検査して (stop_rule, cap, 閾値の上書き, hold) を返す
    (無ければ None・None・{}・`HOLD_STEPS_DEFAULT`; hold は null も既定)。不正値は ValueError (prepare_ns は run dir を作る前に呼ぶ)。"""
    req = cap = None
    if "stop_rule" in evaluate:
        req = evaluate["stop_rule"]
        if not isinstance(req, bool):
            raise ValueError(f"evaluate.stop_rule は true / false ({req!r})")
    if evaluate.get("stop_rule_cap") is not None:
        cap = evaluate["stop_rule_cap"]
        if isinstance(cap, bool) or not isinstance(cap, int) or cap <= 0:
            raise ValueError(f"evaluate.stop_rule_cap は正の整数 ({cap!r})")
    over = evaluate.get("stop_rule_thresholds") or {}
    if not isinstance(over, dict):
        raise ValueError(f"evaluate.stop_rule_thresholds は mapping ({over!r})")
    thresholds(**over)                          # 未知のキー・不正値はここで ValueError
    hold = evaluate.get("stop_rule_hold_steps")
    if hold is None:
        hold = HOLD_STEPS_DEFAULT
    if isinstance(hold, bool) or not isinstance(hold, int) or hold < 0:
        raise ValueError(f"evaluate.stop_rule_hold_steps は 0 以上の整数 ({hold!r})")
    return req, cap, dict(over), int(hold)


def is_explicit(sr: dict | None) -> bool:
    """prepare_info の stop_rule が**明示の要求** (`evaluate.stop_rule: true` か引数 stop_rule=True) か。明示なら使えない設定は例外、
    既定 (キー無し) なら UNAVAILABLE を記録して fallback_steps で回す。統合 (2026-10-14) より前の記録 (`explicit` が無い) は
    requested_source で決める ("evaluate.stop_rule"・"引数" は明示、"既定 (使う)" は既定)。"""
    if not isinstance(sr, dict) or not sr.get("requested"):
        return False
    if "explicit" in sr:
        return bool(sr["explicit"])
    return sr.get("requested_source") in ("evaluate.stop_rule", "引数")


def resolve_settings(evaluate: dict, *, nsteps_arg, n_main: int, n_main_source: str, out_int: int, euler_ref,
                     requested: bool | None = None) -> dict:
    """停止規則の設定を決める。requested: None なら `evaluate.stop_rule` (**無ければ使う** = 既定)、True / False は関数の引数・
    CLI の --no-stop-rule。`explicit` は明示の要求 (引数 True か `evaluate.stop_rule: true`) — 明示なら呼ぶ側 (`prepare_ns`・`run_staged_ns`)
    は使えない設定を例外にし、既定なら `use` = False と理由を記録して本段を fallback_steps で回す。
    上限: 引数 nsteps_arg (`--steps`) があればそれ、無ければ `evaluate.stop_rule_cap` (既定 200,000)。上限は到達 + 保持を含む本段の step 数。
    使わない・使えないときの本段の長さは n_main (prepare_ns の `--steps` > `evaluate.ns_steps` > 60,000) — `fallback_steps` に記録する。
    使えない条件 (`use` = False と理由): euler_ref が無い、出力の間隔が窓に MIN_POINTS 点を入れられない (間隔 > W/(MIN_POINTS − 1))、
    上限が出力の間隔の倍数でない (forge は outStepInterval の倍数でしか res を書かないので、上限の出力が無く CENSORED を見分けられない)。"""
    from ..feedback.deltastar_loop import handoff_extract_settings
    from ..metrics.deltastar import resolve_euler_ref
    yaml_req, yaml_cap, th_over, hold = check_evaluate_keys(evaluate)
    th = thresholds(**th_over)
    if requested is None:
        if yaml_req is not None:
            requested, src, explicit = yaml_req, "evaluate.stop_rule", bool(yaml_req)
        else:
            requested, src, explicit = True, "既定 (使う)", False
    else:
        requested, src = bool(requested), "引数"
        explicit = requested
    # euler_ref: run dir (str) か固定のスナップショット {run, res, sha256} (§5.1 #36e (6); ここでは形だけ検査、sha は prepare_ns・起動前に照合)
    eref = resolve_euler_ref(euler_ref, verify=False)["record"] if euler_ref else None
    rec = {"requested": requested, "requested_source": src, "explicit": explicit,
           "fallback_steps": int(n_main), "fallback_steps_source": n_main_source,
           "out_step_interval": int(out_int), "euler_ref": eref, **th,
           "thresholds_source": ({"evaluate.stop_rule_thresholds": sorted(th_over)} if th_over else "既定"),
           "hold_steps": hold,
           "hold_steps_source": ("evaluate.stop_rule_hold_steps" if evaluate.get("stop_rule_hold_steps") is not None
                                 else f"既定 {HOLD_STEPS_DEFAULT}"),
           "extract": {**handoff_extract_settings(), "source": "feedback.deltastar_loop.handoff_extract_settings (§5.1 #36e)"}}
    if not requested:
        rec.update(use=False, reason="使わない指定", cap=None, cap_source=None)
        return rec
    if nsteps_arg:
        cap, cap_src = int(nsteps_arg), "引数 (--steps)"
    elif yaml_cap is not None:
        cap, cap_src = yaml_cap, "evaluate.stop_rule_cap"
    else:
        cap, cap_src = CAP_DEFAULT, f"既定 {CAP_DEFAULT}"
    why = []
    if not euler_ref:
        why.append("euler_ref が無い (y_b を決められない)")
    if int(out_int) * (th["min_points"] - 1) > th["window"]:
        why.append(f"出力の間隔 {int(out_int)} step では窓 {th['window']} step に {th['min_points']} 点が入らない "
                   f"(evaluate.outStepInterval ≤ {th['window'] // (th['min_points'] - 1)})")
    if int(cap) % int(out_int):
        why.append(f"上限 {int(cap)} step が出力の間隔 {int(out_int)} step の倍数でない (上限の出力が書かれず CENSORED を見分けられない)")
    rec.update(use=not why, reason=("; ".join(why) or None), cap=cap, cap_source=cap_src)
    if cap < th["window"] + int(out_int) + hold:
        rec["note"] = (f"上限 {cap} では最初の出力から窓 {th['window']} step を覆って保持 {hold} step を終えられない (必ず CENSORED)")
    return rec


# ------------------------------------------------------------------ 量
def _read_coords(run_dir: Path):
    with h5py.File(run_dir / "nozzle.h5", "r") as h:
        return np.array(h["/MESH/COORD"][:]).reshape(-1, 3)


def read_values(path) -> dict:
    """res の VALUE の全データセット (保存の型のまま)。"""
    with h5py.File(path, "r") as h:
        return {k: np.array(d[:]) for k, d in h["VALUE"].items()}


def nonfinite_count(V: dict) -> int:
    return int(sum(np.count_nonzero(~np.isfinite(a)) for a in V.values() if a.dtype.kind == "f"))


def mass_flows(X, R, RO, UX, UY, S: float) -> tuple:
    """入口・出口の列の折れ線を通る質量流量 [kg/s] (`cold_xcheck.reduce_fields` の式; X, R は r_t 単位の (ni, nj))。"""
    out = []
    for i in (0, X.shape[0] - 1):
        xx = X[i] * S; rr = R[i] * S
        a = RO[i] * UX[i] * rr; b = RO[i] * UY[i] * rr
        out.append(2 * np.pi * float(np.sum(0.5 * (a[1:] + a[:-1]) * np.diff(rr) - 0.5 * (b[1:] + b[:-1]) * np.diff(xx))))
    return out[0], out[1]


def theta_r_sections(X, R, RO, UX, yb_x, yb, xs_list) -> list:
    """一定 x の断面 (各行 j を x 方向に補間) で、帯の外縁 y_b の (ρ_e, u_e) を縁にした運動量厚さ θ_r (`cold_xcheck.reduce_fields` の式、
    断面は xs に最も近い列)。戻り [(x_断面, θ_r)]。y_b が軸を越える断面は NaN。"""
    xt = X[:, 0]
    nj = X.shape[1]
    out = []
    for xs in xs_list:
        i = int(np.argmin(np.abs(xt - xs)))
        x0 = float(xt[i])
        rr = np.array([np.interp(x0, X[:, j], R[:, j]) for j in range(nj)])
        ro = np.array([np.interp(x0, X[:, j], RO[:, j]) for j in range(nj)])
        ux = np.array([np.interp(x0, X[:, j], UX[:, j]) for j in range(nj)])
        rw = rr[-1]; rb = rw - float(np.interp(x0, yb_x, yb))
        if not np.isfinite(rb) or rb <= rr[0]:
            out.append((x0, float("nan")))
            continue
        rf = np.linspace(rb, rw, 4001)
        rho = np.interp(rf, rr, ro); u = np.interp(rf, rr, ux)
        re_, ue_ = rho[0], u[0]
        qq = np.trapezoid(rho * u * (ue_ - u) * rf, rf)
        out.append((x0, float(rw - np.sqrt(max(rw ** 2 - 2.0 * qq / (re_ * ue_ ** 2), 0.0)))))
    return out


def recorded_quantities(Xf, Rf, V: dict, x_E: float) -> dict:
    """記録だけの量: 出口コア M と無次元の流量 (`nozzle_report.metrics`・`qs_series.py` と同じ定義; 座標は nozzle_report.load_field と同じ
    float32 の COORD / S、VALUE は float64)。"""
    ux = V["Ux"].astype(float); uy = V["Uy"].astype(float)
    M = np.hypot(ux, uy) / V["sonic"].astype(float)
    eta = Rf[-1] / Rf[-1, -1]; core = (eta >= 0.05) & (eta <= 0.7)
    q = V["ro"].astype(float) * ux
    m = 2 * np.pi * np.trapezoid(q * Rf, Rf, axis=1)
    xs = Xf[:, 0]; sel = (xs > -2.5) & (xs < x_E)
    return {"exit_core_M": float(M[-1][core].mean()), "mdot_rt": float(np.median(m[sel]))}


class RunGeometry:
    """run の格子 (prepare_info・nozzle.h5) を一度だけ読む。"""

    def __init__(self, run_dir):
        self.run_dir = Path(run_dir)
        info = json.loads((self.run_dir / "prepare_info.json").read_text())
        self.info = info
        self.ni, self.nj = int(info["mesh"]["ni"]), int(info["mesh"]["nj"])
        self.S = float(info["scale_m"]); self.x_E = float(info["x_E"])
        pw = info.get("physical_wall") or {}
        dom = pw.get("domain") if isinstance(pw, dict) else None
        if not (isinstance(dom, (list, tuple)) and len(dom) == 2):
            raise ValueError(f"{self.run_dir}: prepare_info.physical_wall.domain が無い (θ_r の断面の基準 = 設計壁の出口を決められない)")
        self.x_wall_exit = float(dom[1])     # θ_r の断面の基準 (設計壁の出口)。x_E (軸 Mach 則の終点) ではない
        nc = _read_coords(self.run_dir)
        if nc.shape[0] != self.ni * self.nj:
            raise ValueError(f"{self.run_dir}: 節点数 {nc.shape[0]} が ni × nj = {self.ni * self.nj} でない (構造格子の node run でない)")
        # θ_r・r は M6 と同じ float64 の座標、記録の量は nozzle_report.load_field と同じ float32 の COORD / S
        self.X = (nc[:, 0].astype(float) / self.S).reshape(self.ni, self.nj); self.R = (nc[:, 1].astype(float) / self.S).reshape(self.ni, self.nj)
        self.Xf = (nc[:, 0] / self.S).reshape(self.ni, self.nj); self.Rf = (nc[:, 1] / self.S).reshape(self.ni, self.nj)


def snapshot(geo: RunGeometry, step: int, band: dict | None, x_fracs=X_FRACS, values: dict | None = None) -> dict:
    """出力 res_<step>.h5 の 1 行。band が None (停止規則を使えない) なら θ_r は出さない。非有限があれば量は出さない (DIVERGED の行)。"""
    V = values if values is not None else read_values(geo.run_dir / f"res_{step}.h5")
    for k in ("ro", "Ux", "Uy", "sonic"):
        if k not in V or V[k].shape[0] != geo.ni * geo.nj:
            raise ValueError(f"res_{step}.h5 の VALUE/{k} が無い・長さが違う")
    row = {"step": int(step), "nonfinite_all": nonfinite_count(V)}
    if row["nonfinite_all"]:
        return row
    sh = (geo.ni, geo.nj)
    RO = V["ro"].astype(float).reshape(sh); UX = V["Ux"].astype(float).reshape(sh); UY = V["Uy"].astype(float).reshape(sh)
    mi, mo = mass_flows(geo.X, geo.R, RO, UX, UY, geo.S)
    row.update(mdot_in=mi, mdot_out=mo, r=(mi - mo) / mi if mi != 0 else float("nan"))
    if band is not None:
        xs = [f * geo.x_wall_exit for f in x_fracs]
        for f, (x0, th) in zip(x_fracs, theta_r_sections(geo.X, geo.R, RO, UX, np.asarray(band["x"]), np.asarray(band["y_b"]), xs)):
            row[theta_key(f)] = th; row[f"x_{theta_key(f)}"] = x0
    row.update(recorded_quantities(geo.Xf, geo.Rf, {k: V[k].reshape(sh) for k in ("ro", "Ux", "Uy", "sonic")}, geo.x_E))
    return row


# ------------------------------------------------------------------ 帯の位置 y_b (本段の最初の出力 × Euler 参照)
def resolve_path(p) -> Path | None:
    """prepare_info の euler_ref (相対なら cwd、無ければリポジトリ直下からの相対) を解決する。無ければ None。固定のスナップショットの指定
    {run, res, sha256} (§5.1 #36e) は run のディレクトリを解決する (sha の照合は `check_euler_ref`)。"""
    if isinstance(p, dict):
        p = p.get("run")
    if not p:
        return None
    cand = [Path(p)]
    if not Path(p).is_absolute():
        cand.append(Path(__file__).resolve().parents[3] / p)
    for c in cand:
        if c.is_dir():
            return c.resolve()
    return None


def euler_ref_spec(euler_ref):
    """prepare_info の euler_ref を抽出器に渡す指定にする: run dir (str、`resolve_path` で解決) か固定のスナップショット
    {run (解決した絶対パス), res, sha256}。見つからなければ None。"""
    eu = resolve_path(euler_ref)
    if eu is None:
        return None
    if isinstance(euler_ref, dict):
        return {"run": str(eu), "res": euler_ref.get("res"), "sha256": euler_ref.get("sha256")}
    return eu


def check_euler_ref(euler_ref) -> dict:
    """euler_ref を解決して固定のスナップショットなら sha256 を照合する (§5.1 #36e (6))。見つからない・sha が違えば ValueError。
    戻り: `deltastar.resolve_euler_ref` の戻り。"""
    from ..metrics.deltastar import resolve_euler_ref
    spec = euler_ref_spec(euler_ref)
    if spec is None:
        raise ValueError(f"euler_ref が無い・見つからない ({euler_ref!r})")
    return resolve_euler_ref(spec, verify=True)


def band_extract_kw(sr: dict) -> dict:
    """帯の位置の抽出 (`fix_band`) に渡す条件: 受け渡しの抽出と同じ一定 x の断面への補間と端の扱い (prepare_info の stop_rule.extract、
    無ければ今の既定 `handoff_extract_settings`)。帯の位置は抽出器の縁の推定 (δ_r の平滑化・緩和の前) で決まるので、ω・平滑化・昇格の上限は
    渡さない (昇格の上限で帯の位置が使えなくならないように)。"""
    from ..feedback.deltastar_loop import handoff_extract_settings
    ex = sr.get("extract") or handoff_extract_settings()
    return {k: ex[k] for k in ("remap_constant_x", "remap_edge") if k in ex}


def fix_band(run_dir, step: int, euler_ref, out_dir, **extract_kw) -> dict:
    """本段の最初の出力 res_<step>.h5 を Euler 参照で `extract_and_merge(.., band_select="edge")` と同じ抽出にかけ、帯の位置 y_b(x) [r_t] を返す
    (`delta_r_edge.csv` の y_b = 縁 y* の edge_c (1.25) 倍; さらに掛けない)。抽出は run の他の出力を読まないよう、prepare_info.json・
    nozzle.h5・その 1 つの res だけを一時ディレクトリに symlink して行う (`_load_structured` は最新の res を読むため)。
    使えない (euler_ref が無い・抽出が失敗・y_b が有限でない) ときは Unavailable。euler_ref が固定のスナップショット {run, res, sha256} なら
    その res を読み、sha が合わなければ Unavailable (§5.1 #36e (6))。extract_kw は `band_extract_kw` (一定 x の断面への補間・端の扱い)。"""
    run_dir = Path(run_dir).resolve(); od = Path(out_dir)
    eu = euler_ref_spec(euler_ref)
    if eu is None:
        raise Unavailable(f"euler_ref が無い・見つからない ({euler_ref!r})")
    res = run_dir / f"res_{step}.h5"
    t0 = time.time()
    try:
        from ..feedback.deltastar_loop import extract_and_merge
        with tempfile.TemporaryDirectory(prefix="stop_rule_band_") as td:
            src = Path(td) / run_dir.name
            src.mkdir()
            for name in ("prepare_info.json", "nozzle.h5"):
                (src / name).symlink_to(run_dir / name)
            (src / res.name).symlink_to(res)
            extract_and_merge(src, eu, band_select="edge", out_dir=od, **extract_kw)
        edge = np.loadtxt(od / "delta_r_edge.csv", delimiter=",", skiprows=1, ndmin=2)
    except Exception as e:  # noqa: BLE001 — 抽出の失敗は「使えない」として記録する (例外として外へ出さない)
        raise Unavailable(f"抽出が失敗: {type(e).__name__}: {e}") from e
    x, yb = edge[:, 0], edge[:, 5]
    if len(x) < 2 or not np.all(np.isfinite(yb)) or not np.all(yb > 0):
        raise Unavailable(f"帯の位置 y_b が有限・正でない (点 {len(x)}、非有限 {int(np.count_nonzero(~np.isfinite(yb)))})")
    return {"x": x.tolist(), "y_b": yb.tolist(), "step": int(step), "res_sha256": _sha256(res),
            "euler_ref": (eu if isinstance(eu, dict) else str(eu)),
            "extract": {"function": "feedback.deltastar_loop.extract_and_merge", "band_select": "edge", **{k: repr(v) for k, v in extract_kw.items()}},
            "out_dir": str(od), "seconds": round(time.time() - t0, 1)}


# ------------------------------------------------------------------ 出力ごとの δ_E・δ_target の記録 (§5.1 #36e (7); 判定には使わない)
DELTA_WAIT_S = 7200.0      # 本段の後に残りの抽出を待つ上限 [s] (抽出 1 本 ~2 分 × 追いつかない分; 飛ばす規則で 2W/出力の間隔 + 1 本まで)


def use_columns(run_dir) -> dict:
    """δ の記録の列: NS の列の x (抽出器の `_load_structured` と同じ式 COORD / S) のうち使用域 [0, x_w] (x_w = 設計壁の出口
    physical_wall.domain[1]) の列。戻り {"x" (r_t、列の昇順), "i" (列の番号), "x_w", "scale_m"}。"""
    run_dir = Path(run_dir)
    info = json.loads((run_dir / "prepare_info.json").read_text())
    ni, nj = int(info["mesh"]["ni"]), int(info["mesh"]["nj"])
    S = float(info["scale_m"])
    xw = float((info.get("physical_wall") or {})["domain"][1])
    with h5py.File(run_dir / "nozzle.h5", "r") as h:
        nc = h["/MESH/COORD"][:].reshape(-1, 3)
    x = (nc[:, 0].reshape(ni, nj) / S)[:, 0]
    sel = np.where((x >= 0.0) & (x <= xw + 1e-9))[0]
    return {"x": x[sel], "i": sel, "x_w": xw, "scale_m": S}


def delta_snapshot(run_dir, step: int, euler_ref, extract: dict, x_use) -> dict:
    """出力 res_<step>.h5 1 つを受け渡しの抽出と同じ条件 (`extract` = `handoff_extract_settings`: edge・一定 x の断面への補間・端は欠測・
    ω・平滑化の λ・昇格の上限、固定した Euler の参照) で抽出し、使用域の列 x_use の δ_E (未緩和の平滑化抽出)・δ_target (緩和・再平滑化の
    後の次の壁の入力)・held を返す (r_t 単位、抽出しなかった列は NaN と `missing`)。抽出は run の他の出力を読まないよう一時ディレクトリに
    prepare_info.json・nozzle.h5・その res だけを symlink して行い、出力の表は捨てる (行の値だけを stop_rule.json に残す)。"""
    from ..feedback.deltastar_loop import extract_and_merge, read_delta_r_next
    run_dir = Path(run_dir).resolve()
    res = run_dir / f"res_{step}.h5"
    t0 = time.time()
    kw = {k: v for k, v in extract.items() if k != "source"}
    with tempfile.TemporaryDirectory(prefix="stop_rule_delta_") as td:
        src = Path(td) / run_dir.name
        src.mkdir()
        for name in ("prepare_info.json", "nozzle.h5"):
            (src / name).symlink_to(run_dir / name)
        (src / res.name).symlink_to(res)
        summ = extract_and_merge(src, euler_ref, out_dir=Path(td) / "out", **kw)
        nx = read_delta_r_next(Path(td) / "out" / "delta_r_next.csv")
    x_use = np.asarray(x_use, dtype=float)
    xe = np.atleast_1d(np.asarray(nx["x_rt"], dtype=float))
    if len(xe) == 0:
        raise ValueError(f"res_{step}.h5 の抽出の表が空")
    # 列の対応は x の一致で取る (抽出器の x と use_columns の x は同じ式 COORD / S。表の x は '%.18e' で往復しても同じ値)
    k = np.clip(np.searchsorted(xe, x_use), 0, len(xe) - 1)
    hit = xe[k] == x_use
    dE = np.where(hit, np.atleast_1d(nx["delta_E"])[k], np.nan)
    dT = np.where(hit, np.atleast_1d(nx["delta_target"])[k], np.nan)
    held = np.where(hit, np.atleast_1d(nx["held"])[k] != 0, False)
    return {"delta_E": dE.tolist(), "delta_target": dT.tolist(), "held": [int(v) for v in held],
            "missing": [int(j) for j in np.where(~hit)[0]], "n_missing": int((~hit).sum()), "n_held": int(held.sum()),
            "lam": summ["smooth"]["lam"], "res_sha256": _sha256(res), "seconds": round(time.time() - t0, 2)}


class DeltaRecorder:
    """出力ごとの δ_E・δ_target の抽出 (`delta_snapshot`) を別スレッド 1 本で回す (停止の判定を抽出で待たせない)。submit した順に処理し、
    最新の submit から span (= 2 窓 2W) より古い出力は抽出せず `skipped` にする (停止の時刻 s は最新の出力以後なので [s − 2W, s] に入らない)。
    例外は行の `error` に残す (記録の量なので見張りは止めない)。"""

    def __init__(self, run_dir, euler_ref, extract: dict, x_use, span: int) -> None:
        import queue
        import threading
        self.run_dir, self.euler_ref, self.extract, self.x_use, self.span = Path(run_dir), euler_ref, dict(extract), list(x_use), int(span)
        self._q = queue.Queue()
        self._lock = threading.Lock()
        self._done: dict = {}
        self._latest = None
        self._stop = False
        self._t = threading.Thread(target=self._work, name="stop_rule_delta", daemon=True)
        self._t.start()

    def submit(self, step: int) -> None:
        with self._lock:
            self._latest = int(step) if self._latest is None else max(self._latest, int(step))
        self._q.put(int(step))

    def _work(self) -> None:
        while True:
            step = self._q.get()
            try:
                if step is None:
                    return
                with self._lock:
                    latest, stop = self._latest, self._stop
                if stop:
                    out = {"error": "見張りが終わったので抽出しない"}
                elif latest is not None and step < latest - self.span:
                    out = {"skipped": f"最新の出力 {latest} から 2 窓 {self.span} step より古い (受け渡しの 2 窓に入らない)"}
                else:
                    try:
                        out = delta_snapshot(self.run_dir, step, self.euler_ref, self.extract, self.x_use)
                    except Exception as e:  # noqa: BLE001 — 記録の量の抽出の失敗は行に残す (見張りは止めない)
                        out = {"error": f"{type(e).__name__}: {e}"}
                with self._lock:
                    self._done[step] = out
            finally:
                self._q.task_done()

    def collect(self) -> dict:
        """終わった抽出 {step: 結果} を取り出す (取り出した分は消す)。"""
        with self._lock:
            out, self._done = self._done, {}
        return out

    def finish(self, wait_s: float, abandon: bool = False) -> dict:
        """残りの抽出を wait_s まで待って止める (abandon なら待たずに未着手を捨てる)。戻り: 取り出していない結果。"""
        if abandon:
            with self._lock:
                self._stop = True
        self._q.put(None)
        self._t.join(timeout=max(float(wait_s), 0.0))
        out = self.collect()
        out["_unfinished"] = bool(self._t.is_alive())
        return out


def _attach_deltas(rows: list, got: dict) -> int:
    """`DeltaRecorder` の結果を行の `delta` に付ける。付けた数を返す。"""
    by = {r["step"]: r for r in rows}
    n = 0
    for st, d in got.items():
        if isinstance(st, int) and st in by:
            by[st]["delta"] = d
            n += 1
    return n


# ------------------------------------------------------------------ 判定
def judge_rows(rows: list, k: int | None = None, *, window=WINDOW, min_points=MIN_POINTS, r_mean_max=R_MEAN_MAX, r_amp_max=R_AMP_MAX,
               theta_drift_max=THETA_DRIFT_MAX, theta_amp_max=THETA_AMP_MAX, x_fracs=X_FRACS, relative_keys=()) -> dict:
    """rows[k] (step 昇順、本段の出力、step > 0) の時点の判定と内訳。ok なら停止の条件を満たす。
    relative_keys: θ_r の断面のほかに、θ_r と同じ式・閾値 (D ≤ theta_drift_max・振幅/|平均| ≤ theta_amp_max) で判定する行のキー
    (結果は out["relative"]、ok に含める)。停止規則 (見張り・first_stop) は渡さない (空 = 判定は今のまま)。計算の後の判定
    (`judge_run.py` の criterion stop_rule) が出口コア M をこの関数で判定するために足した (定義を二重に持たない;
    plan tooling-convergence-and-wall-resolution-gates §5.1 #12b (2)・#14a (2))。"""
    k = len(rows) - 1 if k is None else int(k)
    s_end = rows[k]["step"]
    win = [r for r in rows[:k + 1] if r["step"] >= s_end - window]
    out = {"step": s_end, "n_window": len(win), "history_span": s_end - rows[0]["step"], "covers": s_end - rows[0]["step"] >= window,
           "ok": False}
    if not out["covers"]:
        out["why"] = f"本段の履歴 {out['history_span']} step が窓 {window} step を覆わない"
        return out
    if len(win) < min_points:
        out["why"] = f"窓の点数 {len(win)} < {min_points}"
        return out
    s = np.array([w["step"] for w in win], float)
    second = s >= s_end - window / 2.0           # 中点ちょうどの出力は後半
    first = ~second
    if not first.any() or not second.any():
        out["why"] = "窓の前半か後半に出力が無い"
        return out
    rv = np.array([w.get("r", np.nan) for w in win], float)
    out["r_mean"] = float(rv.mean()); out["r_amp"] = float((rv.max() - rv.min()) / 2.0)
    out["r_ok"] = bool(np.isfinite(rv).all() and abs(out["r_mean"]) <= r_mean_max and out["r_amp"] <= r_amp_max)
    ds = float(s[second].mean() - s[first].mean())

    def _relative(key):
        # θ_r の式 (前半・後半の平均の差を 2 万 step あたりに換算した D と、振幅 / |平均|)
        v = np.array([w.get(key, np.nan) for w in win], float)
        mean = float(v.mean())
        d_raw = float(abs(v[second].mean() - v[first].mean()) / abs(mean)) if mean else float("nan")
        amp = float((v.max() - v.min()) / 2.0 / abs(mean)) if mean else float("nan")
        D = d_raw * window / ds
        return {"mean": mean, "D": D, "D_raw": d_raw, "amp": amp,
                "ok": bool(np.isfinite(v).all() and D <= theta_drift_max and amp <= theta_amp_max)}
    th = {theta_key(f): _relative(theta_key(f)) for f in x_fracs}
    out["theta"] = th; out["half_step_gap"] = ds
    rel = {key: _relative(key) for key in relative_keys}
    if relative_keys:
        out["relative"] = rel
    out["ok"] = bool(out["r_ok"] and all(t["ok"] for t in th.values()) and all(t["ok"] for t in rel.values()))
    if not out["ok"]:
        bad = ((["r"] if not out["r_ok"] else []) + [k_ for k_, t in th.items() if not t["ok"]]
               + [k_ for k_, t in rel.items() if not t["ok"]])
        out["why"] = "条件を満たさない: " + ", ".join(bad)
    return out


def row_data_problem(row: dict, x_fracs=X_FRACS) -> str | None:
    """場が有限なのに判定量が有限でない (y_b が軸を越える・流量 0 など) → 理由。問題なければ None。"""
    keys = ["r"] + [theta_key(f) for f in x_fracs]
    bad = [k for k in keys if not np.isfinite(row.get(k, np.nan))]
    return f"判定量が有限でない: {', '.join(bad)}" if bad else None


def _hold_step(hs: dict, ok: bool, step: int, hold: int, why=None) -> bool:
    """到達の後の継続 (`hold_steps`; M6 の版で足した、既定 `HOLD_STEPS_DEFAULT`) の状態 hs = {start, reached, breaks} を出力 1 つぶん
    進め、REACHED にするかを返す。
    条件を満たした最初の出力 s₀ から数え、step − s₀ ≥ hold になった出力で True。途中で条件が外れたら到達を取り消して breaks に記録し、
    次に満たした出力から数え直す。hold 0 は条件を満たした出力で True (forge-f5 の版の「最初の到達で止める」と同じ)。"""
    if not ok:
        if hs.get("start") is not None:
            hs["breaks"].append({"start": hs["start"], "broken_at": int(step), "why": why})
            hs["start"] = None
        return False
    if hs.get("start") is None:
        hs["start"] = int(step)
        hs["reached"].append(int(step))
    return int(step) - hs["start"] >= int(hold)


def first_stop(rows: list, hold_steps: int = 0, **th) -> dict | None:
    """保存済みの系列の最初の終わり {status: REACHED / DIVERGED / DATA_ERROR, k, step, judge}。無ければ None。
    hold_steps > 0 は到達の後 hold_steps の間すべての出力で条件を保った出力を REACHED にし、`hold` に経過を付ける。**関数の既定 0 は
    最初の到達** (系列の再生用); 見張り (`watch_main_stage`) と同じ結果にするには prepare_info の stop_rule.hold_steps (既定
    `HOLD_STEPS_DEFAULT`) を渡す。"""
    fr = th.get("x_fracs", X_FRACS)
    hs = {"start": None, "reached": [], "breaks": []}
    for k, row in enumerate(rows):
        if row.get("nonfinite_all"):
            return {"status": "DIVERGED", "k": k, "step": row["step"], "judge": None}
        p = row_data_problem(row, fr)
        if p:
            return {"status": "DATA_ERROR", "k": k, "step": row["step"], "judge": None, "note": p}
        j = judge_rows(rows, k, **th)
        if _hold_step(hs, j["ok"], row["step"], hold_steps, j.get("why")):
            out = {"status": "REACHED", "k": k, "step": row["step"], "judge": j}
            if hold_steps:
                out["hold"] = {"hold_steps": int(hold_steps), **hs}
            return out
    return None


# ------------------------------------------------------------------ 記録・後処理
def _sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def save_record(run_dir, rec: dict) -> Path:
    """`<run>/stop_rule.json` を原子的に書く (一時ファイル → rename)。"""
    p = Path(run_dir) / RECORD
    tmp = p.with_name(p.name + ".tmp")
    tmp.write_text(json.dumps(rec, indent=1, ensure_ascii=False, default=_json_default))
    os.replace(tmp, p)
    return p


def write_state(run_dir, state: dict) -> Path:
    """本段の終了の状態 `<run>/STOP_RULE_STATE.json` を原子的に書く (一時ファイル → rename; M6 の版で足した)。state は少なくとも status を持つ。
    `run_staged_ns` が本段の後に 1 度だけ書く。見張りの途中の記録 (stop_rule.json、status は RUNNING から変わる) と違い、このファイルが
    あれば本段は終わっている (無いのに prepare_info が停止規則を要求していれば、本段が終わっていないか記録が失われた)。"""
    if state.get("status") not in STATES + ("NOT_USED",):
        raise ValueError(f"停止の状態 {state.get('status')!r} は {STATES + ('NOT_USED',)} のどれか")
    p = Path(run_dir) / STATE
    tmp = p.with_name(p.name + ".tmp")
    tmp.write_text(json.dumps({"version": 1, "plan": PLAN, **state}, indent=1, ensure_ascii=False, default=_json_default))
    os.replace(tmp, p)
    return p


def read_state(run_dir) -> dict | None:
    """`<run>/STOP_RULE_STATE.json` を読む (無ければ None)。"""
    p = Path(run_dir) / STATE
    return json.loads(p.read_text()) if p.exists() else None


def _json_default(o):
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, Path):
        return str(o)
    raise TypeError(type(o))


def _now() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def res_steps(run_dir) -> list:
    return sorted(int(m.group(1)) for f in Path(run_dir).iterdir() if (m := _RES_RE.fullmatch(f.name)))


def delete_after(run_dir, step: int) -> list:
    """止めた出力 step より後の出力 (res_<N>.h5/.xmf と境界の res_<名前>_<id>_<N>.h5/.xmf、N > step) を消す。消した名前を返す。"""
    gone = []
    for f in sorted(Path(run_dir).iterdir()):
        m = _OUT_RE.fullmatch(f.name)
        if m and int(m.group(1)) > step:
            f.unlink(); gone.append(f.name)
    return gone


def postprocess(run_dir) -> dict:
    """`run_case.sh` と同じ後処理 (check_convergence.py → CONVERGENCE_VERDICT.txt、plot_implicit_residuals.py → residual_history.png)。
    計算の後の判定 (judge_run.py) はここでは回さない: 本段の総合は run_staged_ns が停止の状態の記録の後に `judge.finalize` で確定する。"""
    from .runner import FORGE_TOOLS, _ENV
    run_dir = Path(run_dir)
    with open(run_dir / "CONVERGENCE_VERDICT.txt", "w") as fh:
        r = subprocess.run(["python3", str(FORGE_TOOLS / "check_convergence.py"), "."], cwd=run_dir, env=_ENV, stdout=fh, stderr=subprocess.STDOUT)
    with open(run_dir / "CONVERGENCE_VERDICT.txt", "a") as fh:
        fh.write(f"CONVERGENCE_EXIT: {r.returncode}\n")
    subprocess.run(["python3", str(FORGE_TOOLS / "plot_implicit_residuals.py"), "--input", "residual_history.csv", "--output", "residual_history.png"],
                   cwd=run_dir, env=_ENV, capture_output=True)
    return {"by": "watcher", "convergence_exit": r.returncode}


# ------------------------------------------------------------------ forge のプロセス
def _proc_stat(pid: int):
    try:
        s = Path(f"/proc/{pid}/stat").read_text()
    except OSError:
        return None
    rest = s[s.rfind(")") + 2:].split()
    return rest[0], int(rest[1])          # 状態, ppid


def _alive(pid: int) -> bool:
    st = _proc_stat(pid)
    return st is not None and st[0] not in ("Z", "X")


def _children(pid: int) -> list:
    out = []
    for d in Path("/proc").iterdir():
        if d.name.isdigit():
            st = _proc_stat(int(d.name))
            if st is not None and st[1] == pid:
                out.append(int(d.name))
    return out


def forge_exe() -> Path:
    """`run_case.sh` が起動する forge (FORGE_BIN、無ければ <FORGE_ROOT>/solver_density_cuda/build/forge)。"""
    from .runner import FORGE_ROOT, _ENV
    return Path(_ENV.get("FORGE_BIN") or (FORGE_ROOT / "solver_density_cuda" / "build" / "forge"))


def _is_forge(pid: int, exe: Path) -> bool:
    try:
        argv = Path(f"/proc/{pid}/cmdline").read_bytes().split(b"\0")
        cwd = os.readlink(f"/proc/{pid}/cwd")
    except OSError:
        return False
    want = os.path.realpath(exe)
    for a in argv[:2]:
        if a and os.path.realpath(os.path.join(cwd, a.decode(errors="replace"))) == want:
            return True
    try:
        return os.path.realpath(os.readlink(f"/proc/{pid}/exe")) == want
    except OSError:
        return False


def _forge_exit(run_dir) -> int | None:
    p = Path(run_dir) / "run_case_stdout.log"
    if not p.exists():
        return None
    m = _FORGE_EXIT_RE.findall(p.read_text(errors="replace"))
    return int(m[-1]) if m else None


def stop_forge(bash_pid: int, run_dir, term_wait_s: float = 180.0, kill_wait_s: float = 60.0, find_wait_s: float = 10.0) -> dict:
    """`run_case.sh` (bash_pid) の子の forge だけに SIGTERM → term_wait_s 待って残れば SIGKILL → kill_wait_s 待ち、止まったことを確かめる。
    forge だけを止めるので run_case.sh の後処理はそのまま走る。forge が見つからず、まだ終わった記録も無いときは (find_wait_s 待った後)
    プロセスグループごと止める (後処理は呼ぶ側が補う)。"""
    exe = forge_exe()
    t0 = time.time()
    pids = [p for p in _children(bash_pid) if _is_forge(p, exe)]
    while not pids and _forge_exit(run_dir) is None and _alive(bash_pid) and time.time() - t0 < find_wait_s:
        time.sleep(0.2)
        pids = [p for p in _children(bash_pid) if _is_forge(p, exe)]
    rec = {"forge_pids": pids, "signals": [], "group_killed": False}
    if not pids:
        if _forge_exit(run_dir) is None and _alive(bash_pid):
            for sig in (signal.SIGTERM, signal.SIGKILL):
                try:
                    os.killpg(bash_pid, sig)
                except OSError:
                    pass
                rec["signals"].append(f"killpg {sig.name}")
                t1 = time.time()
                while _alive(bash_pid) and time.time() - t1 < kill_wait_s:
                    time.sleep(0.2)
                if not _alive(bash_pid):
                    break
            rec["group_killed"] = True
        rec["stopped"] = not _alive(bash_pid) or _forge_exit(run_dir) is not None
        return rec
    for sig, wait in ((signal.SIGTERM, term_wait_s), (signal.SIGKILL, kill_wait_s)):
        live = [p for p in pids if _alive(p)]
        if not live:
            break
        for p in live:
            try:
                os.kill(p, sig)
            except OSError:
                pass
        rec["signals"].append(sig.name)
        t1 = time.time()
        while any(_alive(p) for p in pids) and time.time() - t1 < wait:
            time.sleep(0.2)
    rec["stopped"] = not any(_alive(p) for p in pids)
    return rec


# ------------------------------------------------------------------ 見張り (run_staged_ns の本段)
def watch_main_stage(run_dir, sr: dict, *, poll_s: float = 10.0, settle_s: float = 60.0, term_wait_s: float = 180.0,
                     kill_wait_s: float = 60.0, post_wait_s: float = 1800.0, read_tries: int = 3, delta_wait_s: float = DELTA_WAIT_S) -> dict:
    """本段の forge を `run_case.sh` の子プロセスで起動し、停止規則で見張る。sr は prepare_info.json の stop_rule (use = True、cap が
    solverConfig の nStepOuter に書いてあること)。戻り: {status, stop_step, forge_rc, rc}。rc は呼ぶ側に返す終了コード
    (REACHED・CENSORED・UNAVAILABLE で ns_steps まで回った = 0、それ以外は非 0)。記録は `<run>/stop_rule.json`。

    書き終わった出力: それより後の出力がある、forge が終わった、または最後に書かれてから settle_s 以上。読めない出力は read_tries 回まで
    次の見回りで読み直す (それでも読めなければ DATA_ERROR)。"""
    from .runner import FORGE_TOOLS, _ENV
    run_dir = Path(run_dir).resolve()
    th = {k: sr[k] for k in ("window", "min_points", "r_mean_max", "r_amp_max", "theta_drift_max", "theta_amp_max", "x_fracs")}
    cap = int(sr["cap"]); fallback = int(sr["fallback_steps"])
    # 到達の後の継続 (0 = 最初の到達で止める)。prepare_info の stop_rule.hold_steps (resolve_settings の既定 HOLD_STEPS_DEFAULT)。
    # 統合 (2026-10-14) より前に準備した run の記録にキーが無ければ 0 (準備した時の動きのまま)
    hold = int(sr.get("hold_steps") or 0)
    pre = [n for n in res_steps(run_dir) if n > 0]
    if pre:
        raise RuntimeError(f"{run_dir}: 本段の起動前に出力 res_{pre[0]}.h5 ほか {len(pre)} 個がある — 停止規則の系列に混ざるので起動しない")
    geo = RunGeometry(run_dir)
    rec = {"version": 1, "plan": PLAN, "settings": sr, "status": "RUNNING",
           "started": _now(), "band": None, "unavailable": None, "rows": [], "judge_last": None, "tries": {},
           "timing": {"poll_s": poll_s, "settle_s": settle_s, "term_wait_s": term_wait_s, "kill_wait_s": kill_wait_s},
           "hold": ({"hold_steps": hold, "start": None, "reached": [], "breaks": []} if hold else None),
           # 要求の断面の位置 (x_fracs × x_w)。実際に選ばれた列の x は各行の x_theta_r_<比> (§6 V-h370)
           "sections": [{"key": theta_key(f), "x_frac": float(f), "x_requested": float(f) * geo.x_wall_exit} for f in th["x_fracs"]]}
    # 出力ごとの δ_E・δ_target の記録 (§5.1 #36e (7); 判定には使わない)。列は使用域 [0, x_w] の NS の列、条件は受け渡しの抽出と同じ
    from ..feedback.deltastar_loop import handoff_extract_settings
    ex = sr.get("extract") or {**handoff_extract_settings(), "source": "既定 (prepare_info の stop_rule に extract が無い)"}
    uc = use_columns(run_dir)
    rec["delta_series"] = {"x": uc["x"].tolist(), "use_domain": [0.0, uc["x_w"]], "scale_m": uc["scale_m"], "extract": ex,
                           "euler_ref": sr.get("euler_ref"), "span_steps": 2 * int(th["window"]),
                           "note": "判定には使わない記録の量。行の delta に δ_E (未緩和の平滑化抽出)・δ_target (緩和・再平滑化の後)・held [r_t]"}
    recorder = DeltaRecorder(run_dir, euler_ref_spec(sr.get("euler_ref")), ex, uc["x"], 2 * int(th["window"]))
    save_record(run_dir, rec)
    t_launch = time.time()
    # FORGE_JUDGE_DEFER=1: run_case.sh の judge_run.py は本段の結果だけを記録し、総合の判定は run_staged_ns が停止の状態の記録の後に
    # 確定する (plan tooling-convergence-and-wall-resolution-gates §4.6-2・§5.1 #12)
    from .judge import defer_env
    with open(run_dir / "run_case_stdout.log", "w") as log:
        proc = subprocess.Popen([str(FORGE_TOOLS / "run_case.sh"), str(run_dir)], env=defer_env(_ENV), stdout=log, stderr=subprocess.STDOUT,
                                start_new_session=True)
    rec["run_case_pid"] = proc.pid
    stop = None                      # {"status", "step", ...}: 見張りが止めると決めた出力
    rows = rec["rows"]
    seen = set()

    def complete(n, steps_now, done):
        if done or any(m > n for m in steps_now) or _forge_exit(run_dir) is not None:
            return True
        try:
            return time.time() - (run_dir / f"res_{n}.h5").stat().st_mtime >= settle_s
        except OSError:
            return False

    try:
        while True:
            if _attach_deltas(rows, recorder.collect()):
                save_record(run_dir, rec)
            done = proc.poll() is not None
            steps_now = res_steps(run_dir)
            for n in [m for m in steps_now if m > 0 and m not in seen]:
                if not complete(n, steps_now, done):
                    break
                try:
                    V = read_values(run_dir / f"res_{n}.h5")
                    if rec["band"] is None and rec["unavailable"] is None:
                        try:
                            rec["band"] = fix_band(run_dir, n, sr.get("euler_ref"), run_dir / BAND_DIR, **band_extract_kw(sr))
                        except Unavailable as e:
                            rec["unavailable"] = {"reason": str(e), "step": n,
                                                  "action": f"停止規則を使わず fallback_steps {fallback} step で止める"}
                        save_record(run_dir, rec)
                    row = snapshot(geo, n, rec["band"], th["x_fracs"], values=V)
                except Exception as e:  # noqa: BLE001
                    rec["tries"][str(n)] = rec["tries"].get(str(n), 0) + 1
                    rec["last_read_error"] = f"res_{n}.h5: {type(e).__name__}: {e}"
                    if rec["tries"][str(n)] >= read_tries:
                        stop = {"status": "DATA_ERROR", "step": n, "note": f"出力が読めない ({rec['last_read_error']})",
                                "read_error": True, "forge_done": done or _forge_exit(run_dir) is not None}
                    save_record(run_dir, rec)
                    break
                seen.add(n)
                if row["nonfinite_all"]:
                    stop = {"status": "DIVERGED", "step": n, "note": f"非有限 {row['nonfinite_all']} 個"}
                elif rec["unavailable"] is None:
                    p = row_data_problem(row, th["x_fracs"])
                    if p:
                        stop = {"status": "DATA_ERROR", "step": n, "note": p}
                rows.append(row)
                if not row["nonfinite_all"]:
                    recorder.submit(n)
                if stop is None and rec["unavailable"] is None:
                    j = judge_rows(rows, **th)
                    row["judge"] = j
                    rec["judge_last"] = j
                    if hold:
                        if _hold_step(rec["hold"], j["ok"], n, hold, j.get("why")):
                            stop = {"status": "REACHED", "step": n,
                                    "note": f"最初の到達 {rec['hold']['start']} の後 {hold} step のすべての出力で条件を保った"}
                    elif j["ok"]:
                        stop = {"status": "REACHED", "step": n}
                if stop is None and rec["unavailable"] is not None and n >= fallback:
                    stop = {"status": "UNAVAILABLE", "step": n, "note": f"停止規則を使えないので fallback_steps {fallback} の出力で止めた"}
                save_record(run_dir, rec)
                if stop is not None:
                    break
            if stop is not None:
                break
            if done and all(m in seen for m in res_steps(run_dir) if m > 0):
                break
            time.sleep(poll_s)
    except BaseException as e:
        # 見張りの例外でも forge を止めて状態を残す
        _attach_deltas(rows, recorder.finish(0.0, abandon=True))
        rec["stop_forge"] = stop_forge(proc.pid, run_dir, term_wait_s, kill_wait_s) if proc.poll() is None else None
        rec.update(status="EXEC_ERROR", note=f"見張りの例外: {type(e).__name__}: {e}", finished=_now())
        save_record(run_dir, rec)
        if not isinstance(e, Exception):
            raise
        return {"status": "EXEC_ERROR", "stop_step": None, "forge_rc": _forge_exit(run_dir), "rc": 1}

    t_stop = time.time()
    if stop is not None and proc.poll() is None:
        rec["stop_forge"] = stop_forge(proc.pid, run_dir, term_wait_s, kill_wait_s)
        if not rec["stop_forge"]["stopped"]:
            stop = {"status": "EXEC_ERROR", "step": stop["step"], "note": f"forge を止められない (止める理由: {stop['status']})"}
    try:
        proc.wait(timeout=post_wait_s)
    except subprocess.TimeoutExpired:
        for sig in (signal.SIGTERM, signal.SIGKILL):
            try:
                os.killpg(proc.pid, sig)
            except OSError:
                pass
            try:
                proc.wait(timeout=30)
                break
            except subprocess.TimeoutExpired:
                continue
        rec["run_case_killed"] = True
    rec["run_case_rc"] = proc.returncode
    forge_rc = _forge_exit(run_dir)
    rec["forge_rc"] = forge_rc
    # run_case.sh の後処理が止めた後にも走ったか (CONVERGENCE_VERDICT.txt が止めた後に書かれたか)。走っていなければ同じ後処理を呼ぶ
    vp = run_dir / "CONVERGENCE_VERDICT.txt"
    ran = (vp.exists() and vp.stat().st_mtime >= (t_stop if stop is not None else t_launch) - 1.0
           and "CONVERGENCE_EXIT" in vp.read_text(errors="replace"))
    rec["postprocess"] = {"by": "run_case.sh"} if ran else postprocess(run_dir)

    last = max([n for n in res_steps(run_dir) if n > 0], default=0)
    if stop is not None and stop.get("read_error") and stop.get("forge_done") and forge_rc not in (None, 0):
        # forge が異常終了した後に読めない出力 (書きかけ) が残った: 異常終了を DATA_ERROR と取り違えない
        status, step = "EXEC_ERROR", last
        rec["note"] = f"forge が上限 {cap} より前に異常終了 (forge rc {forge_rc}); 書きかけの出力が読めない ({stop['note']})"
    elif stop is not None:
        status, step = stop["status"], stop["step"]
        if stop.get("note"):
            rec["note"] = stop["note"]
    elif forge_rc == 0 and last >= cap:
        status, step = ("UNAVAILABLE" if rec["unavailable"] is not None else "CENSORED"), last
        if status == "CENSORED" and hold and rec["hold"]["start"] is not None:
            rec["note"] = f"到達 {rec['hold']['start']} の後の {hold} step の保持が上限 {cap} までに終わらない (保持の途中で上限)"
    elif forge_rc not in (None, 0) and sorted(run_dir.glob("res_nan_*.h5")):
        status, step = "DIVERGED", None
        nan_files = sorted(run_dir.glob("res_nan_*.h5"))
        rec["note"] = "forge の detectNaN が止めた (" + ", ".join(f.name for f in nan_files) + ")"
        st_nan = [int(m.group(1)) for f in nan_files if (m := re.fullmatch(r"res_nan_(\d+)\.h5", f.name))]
        rec["fail_step"] = min(st_nan) if st_nan else None
    else:
        status, step = "EXEC_ERROR", last
        rec["note"] = (f"forge が上限 {cap} より前に終わった (forge rc {forge_rc}、最後の出力 {last})"
                       + (f"; 未読の出力あり ({rec.get('last_read_error')})" if rec.get("last_read_error") else ""))
    if status in ("REACHED", "DIVERGED", "UNAVAILABLE") and step is not None and stop is not None:
        rec["deleted_after_stop"] = delete_after(run_dir, step)
    if step is not None and (run_dir / f"res_{step}.h5").exists() and status in ("REACHED", "DIVERGED", "UNAVAILABLE", "CENSORED"):
        rec["stop_sha256"] = _sha256(run_dir / f"res_{step}.h5")
    # 残りの δ の抽出を待って行に付ける (止めた出力より後の出力は行に無い。抽出の時間は delta_series に)
    left = recorder.finish(delta_wait_s)
    rec["delta_series"]["unfinished_at_end"] = left.pop("_unfinished")
    _attach_deltas(rows, left)
    secs = [r["delta"]["seconds"] for r in rows if isinstance(r.get("delta"), dict) and "seconds" in r["delta"]]
    rec["delta_series"]["seconds"] = ({"n": len(secs), "mean": float(np.mean(secs)), "max": float(np.max(secs))} if secs else {"n": 0})
    rec.update(status=status, stop_step=step, last_output=last, finished=_now())
    if status == "REACHED":
        rec["recorded_at_stop"] = window_summary(rows, step, th["window"])
    save_record(run_dir, rec)
    ok = status in ("REACHED", "CENSORED") or (status == "UNAVAILABLE" and (stop is not None or forge_rc == 0))
    return {"status": status, "stop_step": step, "forge_rc": forge_rc, "rc": 0 if ok else (forge_rc if forge_rc not in (None, 0, 143, 137) else 1)}


def window_summary(rows: list, step: int, window: int = WINDOW) -> dict:
    """記録だけの量 (出口コア M・流量) の、止めた出力までの窓の平均 ± 振幅 (判定には使わない; §4.6 の記録)。"""
    win = [r for r in rows if step - window <= r["step"] <= step and "exit_core_M" in r]
    out = {}
    for k in ("exit_core_M", "mdot_rt"):
        v = np.array([r[k] for r in win], float)
        if len(v):
            out[k] = {"mean": float(v.mean()), "amp": float((v.max() - v.min()) / 2.0), "n": int(len(v))}
    return out
