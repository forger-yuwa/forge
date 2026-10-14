#!/usr/bin/env python3
"""NS の本段の停止規則 (`evaluate/stop_rule.py`) と受け渡しの制御 (`feedback/deltastar_loop.handoff_gate`) の試験
(plan tooling-nozzle-isothermal-wall-chain §4.9・§5.1 #36a)。

停止規則は forge-f5 の feature/nozzle-cooled-wall-settling 60de8125 から取り込んだ (同 `design/tests/run_stop_rule_tests.py` の (a)(b) と
断面の基準の試験を、この作業ツリー向けに直して写し、閾値の YAML・到達の後の継続・(c)(d) を足した)。旧 `m9_watch.py` (§6.11) と同等とは
言わない新しい規則として試す。

  (a) 合成の系列 (`judge_rows`・`first_stop`、forge 0 step) と設定 (`resolve_settings`・`runner_axismach._stop_rule_settings`): 取り込み元と同じ
      判定の数値 (REACHED の出力・平均の変化の換算・振幅・r・非有限・窓の点数・DATA_ERROR)。設定は**こちらの既定 = 使わない**
      (evaluate.stop_rule が無ければ prepare_info にも書かない)、上限が出力の間隔の倍数でない・euler_ref が無い・stop_rule_cap だけの指定は
      run dir を作る前に例外。量の式: 合成の構造格子で θ_r・入口と出口の流量が M6 の `cold_xcheck.reduce_fields` と一致、出口コア M・流量が
      `nozzle_report.load_field` からの定義と一致。断面の基準は設計壁の出口 physical_wall.domain[1] (x_E ではない; 無ければ例外)、
      `evaluate.stop_rule_thresholds.x_fracs` で case/45 の試験部の断面 x = 40・70・94 r_t を指定できる。閾値の YAML の不正値は例外。
      到達の後の継続 (`hold_steps`): 保てば到達 + H で REACHED、外れたら取り消して数え直す、保ち切れなければ届かない。
  (b) 偽の forge (`FORGE_BIN` を res を書く台本に差し替え、本物の `run_case.sh` 経由) で `run_staged_ns` の本段の見張り: REACHED (止めて後の出力を
      消す・SIGTERM・SIGTERM を無視する forge は SIGKILL)・CENSORED (rc 0 でも状態は CENSORED)・EXEC_ERROR・DIVERGED (非有限・detectNaN)・
      DATA_ERROR・UNAVAILABLE (見張りの途中で帯の位置が決まらない; rc 0)。**どの場合も状態を STOP_RULE_STATE.json と戻り値の `stop_rule` で
      直接返す** (rc 0 の CENSORED・UNAVAILABLE を REACHED と取り違えない)。stages="full" では S1・S2 に当てず本段だけ。到達の後の継続:
      保ったら到達 + H で止める、途中で外れて上限まで届かなければ CENSORED (外れた記録つき)。要求したのに
      euler_ref が見つからない run は forge を起動する前に例外。stop_rule=False は NOT_USED を記録して今までの経路。
  (c) 受け渡しの制御 (`handoff_gate` と `run_from_table`・`run_pass`・`run_pass0_integral`・CLI): 状態が REACHED で収束の判定が PASS のときだけ
      抽出 (band_select="edge" を明示) と fixed_point.json に進む。REACHED でも収束の判定が不合格なら止める (`allow_unconverged` で進む)。
      CENSORED・UNAVAILABLE・NOT_USED は止める (`allow_not_reached` で進む)。DIVERGED・EXEC_ERROR・DATA_ERROR・記録なし (MISSING)・
      収束の判定の DIVERGED・判定不能は上書きでも進まない。`run_pass` は前の pass の run の関門が不合格なら抽出・準備の前に止める
      (新しい run dir を作らない)。判定の結果は HANDOFF_GATE.json・HANDOFF_CONVERGENCE.txt に残る。CLI は止めたら終了コード 3。
  (d) 指定の無い問題の回帰: prepare_info に stop_rule の無い run の `run_staged_ns` は今までの経路 (run_forge、本段の nStepOuter は準備のまま、
      stop_rule.json・STOP_RULE_STATE.json を書かない、戻りは int と等しい)。

実行 (作業ツリーの直下で): /home/sano/work/forge/design/.venv-opt/bin/python design/tests/run_stop_rule_tests.py [--quick]
所要: (a) < 5 s、(b)〜(d) 約 1〜2 分 (偽の forge 約 25 本)。--quick は (b)〜(d) を省く。出力は一時ディレクトリ (全部合格なら消す)。
"""
from __future__ import annotations

import importlib.util
import json
import math
import os
import shutil
import stat
import sys
import tempfile
import types
from pathlib import Path

import h5py
import numpy as np

DESIGN = Path(__file__).resolve().parents[1]
ROOT = DESIGN.parent
QUICK = "--quick" in sys.argv[1:]
TMP = Path(tempfile.mkdtemp(prefix="stop_rule_tests_"))
FAKE_FORGE = TMP / "fake_forge.py"
os.environ["FORGE_BIN"] = str(FAKE_FORGE)            # runner._ENV は import 時の環境を写すので import の前に
os.environ.setdefault("FORGE_NO_DESIGN_REPORT", "1")
sys.path.insert(0, str(DESIGN))

from forge_design.evaluate import runner_axismach as RA  # noqa: E402
from forge_design.evaluate import stop_rule as SR  # noqa: E402
from forge_design.feedback import deltastar_loop as DL  # noqa: E402

M6_XCHECK = ROOT / "case" / "45.isobutane_m6_d155" / "cold_xcheck.py"     # M6 の式 (読むだけ)
# 偽の forge の run の solverConfig の雛形 (取り込み元は case/70 の run の写しを読んだ; この作業ツリーには無いので文字列で持つ)
TEMPLATE_CFG = '''mesh: {discretization: "node", isAxisymmetric: 1, axisCentroidShift: 1, nodeWallDirichlet: 1, meshFileName: "nozzle.h5", valueFileName: "nozzle.h5"}
gpu: 1
solver: "SLAU"
physProp: {thermalMethod: 0, viscMethod: 1, visc: 1.8e-5, thermCond: 0.0257, prandtlLam: 0.72, cp: 1004.5, gamma: 1.4}
time:
  unsteady: 0
  dualTime: 0
  last: {nStepOuter: 10000}
  deltaT: {control: 1, dt: 1e-8, cfl: 4.0, cfl_pseudo: 4.0, dt_min: 1e-9, dt_max: 0.001, blockDPLUR: 1, implicitRelax: 0.7, lowMachPrecond: 0, detectNaN: 1}
  outStepStart: 0
  outStepInterval: 2500
  timeIntegration: 11
  nStepInner: 5
space: {convMethod: 1, limiter: 2}
turbulence: {model: "sst", scalarDiffusion: 1, dilatationCorrection: 2, katoLaunder: 1, wallTreatmentSST: 0, turbulentPrandtl: 0.9}
initial: "uniform_p101325_u10"
output: {level: 1}
'''
TEMPLATE_BC = "wall: {physID: 3, type: wall_isothermal, Tw: 370.0}\n"
FAIL = 0


def check(name, cond, detail=""):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name + (f"  [{detail}]" if detail else ""), flush=True)
    if not cond:
        FAIL += 1


def skip(name, why):
    print(f"skip {name}  [{why}]", flush=True)


def raises(fn, exc=Exception, must=""):
    """fn() が exc を投げ、文言に must を含むか。戻り (ok, 文言)。"""
    try:
        fn()
    except exc as e:  # noqa: BLE001
        return (must in str(e)), str(e)
    return False, "例外にならない"


# ====================================================================================== (a) 合成の系列
def rows_of(steps, theta=None, r=None, nonfinite=None):
    """系列の行。theta(step) → 3 断面共通の θ_r、r(step) → r。"""
    out = []
    for s in steps:
        row = {"step": int(s), "nonfinite_all": int(nonfinite(s)) if nonfinite else 0, "r": float(r(s)) if r else 0.0}
        for f in SR.X_FRACS:
            row[SR.theta_key(f)] = float(theta(s)) if theta else 1.0e-3
        out.append(row)
    return out


def stop_step(rows, **th):
    d = SR.first_stop(rows, **th)
    return None if d is None else (d["status"], d["step"])


def test_synthetic():
    print("== (a) 合成の系列 (取り込み元と同じ数値)", flush=True)
    s1k = list(range(1000, 60001, 1000))
    s25 = list(range(2500, 60001, 2500))
    check("一定 (間隔 1000) → 履歴が 2 万 step を覆う最初の出力 21,000 で REACHED", stop_step(rows_of(s1k)) == ("REACHED", 21000),
          str(stop_step(rows_of(s1k))))
    check("一定 (間隔 2500、窓に 9 点) → 22,500 で REACHED", stop_step(rows_of(s25)) == ("REACHED", 22500), str(stop_step(rows_of(s25))))
    j = SR.judge_rows(rows_of(s1k), 19)
    check("20,000 (履歴 19,000 step) は窓を覆わない", not j["ok"] and not j["covers"], j.get("why", ""))

    lin = lambda rate: (lambda s: 1e-3 * (1.0 + rate * (s - 1000) / 20000.0))  # noqa: E731
    check("θ_r の平均の線形の変化 0.1 %/2 万 step → 届かない", stop_step(rows_of(s1k, theta=lin(1e-3))) is None)
    j = SR.judge_rows(rows_of(s1k, theta=lin(1e-3)), 20)
    t = j["theta"][SR.theta_key(0.4)]
    check("  その D_θ ≈ 0.1 % (換算あり)", abs(t["D"] - 1e-3 / 1.0005) < 1e-9, f"D {t['D']:.6e}")

    for n_pts, steps in ((9, list(range(2500, 22501, 2500))), (21, list(range(1000, 21001, 1000)))):
        th8 = lambda s, s0=steps[0]: 1e-3 * (1.0 + 8e-4 * (s - s0) / 20000.0)  # noqa: E731
        rows = rows_of(steps, theta=th8)
        j = SR.judge_rows(rows)
        t = j["theta"][SR.theta_key(0.7)]
        want_raw = {9: 0.044982e-2, 21: 0.041983e-2}[n_pts]
        check(f"0.08 %/2 万 step の単調な系列 ({n_pts} 点、窓 {j['n_window']} 点) → 届かない (D_θ = {t['D'] * 100:.6f} %)",
              stop_step(rows) is None and j["n_window"] == n_pts and abs(t["D"] - 8e-4 / 1.0004) < 1e-10)
        check(f"  換算しない前後半の差 D_raw = {t['D_raw'] * 100:.6f} % ≤ 0.05 % (換算しなければ止めてしまう)",
              abs(t["D_raw"] - want_raw) < 5e-9 and t["D_raw"] <= SR.THETA_DRIFT_MAX and t["amp"] <= SR.THETA_AMP_MAX,
              f"振幅 {t['amp'] * 100:.6f} %")
        check(f"  閾値は引数で変えられる (theta_drift_max 1e-3 なら {n_pts} 点の系列で REACHED)",
              stop_step(rows, theta_drift_max=1e-3) == ("REACHED", steps[-1]))

    sinus = lambda a: (lambda s: 1e-3 * (1.0 + a * math.sin(2 * math.pi * ((s // 1000) - 1) / 10.0)))  # noqa: E731
    check("振幅 1 %・平均の変化 0 (周期 1 万 step の正弦) → 届かない", stop_step(rows_of(s1k, theta=sinus(1e-2))) is None)
    check("振幅 0.1 %・平均の変化 0 → 21,000 で REACHED", stop_step(rows_of(s1k, theta=sinus(1e-3))) == ("REACHED", 21000))
    check("r の平均 0.06 % → 届かない", stop_step(rows_of(s1k, r=lambda s: 6e-4)) is None)
    check("r の平均 0.04 % → REACHED", stop_step(rows_of(s1k, r=lambda s: 4e-4)) == ("REACHED", 21000))
    check("r の振幅 0.6 % (交互 ±0.6 %) → 届かない", stop_step(rows_of(s1k, r=lambda s: 6e-3 * (-1) ** (s // 1000))) is None)
    check("r の振幅 0.4 % (交互 ±0.4 %) → REACHED", stop_step(rows_of(s1k, r=lambda s: 4e-3 * (-1) ** (s // 1000))) == ("REACHED", 21000))
    check("非有限 (12,000 で) → DIVERGED", stop_step(rows_of(s1k, nonfinite=lambda s: s == 12000)) == ("DIVERGED", 12000))
    check("非有限 (30,000 で、21,000 で届く系列) → 先に REACHED", stop_step(rows_of(s1k, nonfinite=lambda s: s == 30000)) == ("REACHED", 21000))
    s8 = [s for s in range(2500, 22501, 2500) if s != 12500]
    j = SR.judge_rows(rows_of(s8))
    check("窓の点数 8 (間隔 2500 で 1 つ欠け) → 届かない", stop_step(rows_of(s8)) is None and j["n_window"] == 8, j.get("why", ""))
    s3k = list(range(3000, 90001, 3000))
    check("間隔 3000 (窓に 7 点) → 届かない", stop_step(rows_of(s3k)) is None, SR.judge_rows(rows_of(s3k)).get("why", ""))
    rows = rows_of(s1k)
    rows[5][SR.theta_key(0.95)] = float("nan")
    check("場が有限で判定量 (θ_r) が NaN → DATA_ERROR", stop_step(rows) == ("DATA_ERROR", 6000))

    print("== (a) 設定 (resolve_settings・_stop_rule_settings; こちらの既定 = 使わない)", flush=True)
    kw = dict(nsteps_arg=None, n_main=200000, n_main_source="evaluate.nStepOuter", out_int=2500, euler_ref="case/x/run_euler")
    d = SR.resolve_settings({}, **kw)
    check("キー無し → 使わない (取り込み元は使う)", not d["requested"] and not d["use"] and d["requested_source"] == "既定 (使わない)",
          json.dumps({k: d[k] for k in ("requested", "use", "requested_source")}, ensure_ascii=False))
    d = SR.resolve_settings({"stop_rule": True}, **kw)
    check("evaluate.stop_rule: true → 使う・上限 200,000・閾値は取り込み元のまま", d["use"] and d["requested"] and d["cap"] == 200000
          and d["window"] == 20000 and d["min_points"] == 9 and d["r_mean_max"] == 5e-4 and d["theta_drift_max"] == 5e-4
          and d["x_fracs"] == [0.4, 0.7, 0.95] and d["fallback_steps"] == 200000, json.dumps(d, ensure_ascii=False)[:200])
    check("evaluate.stop_rule_cap 150000 → 上限 150000", SR.resolve_settings({"stop_rule": True, "stop_rule_cap": 150000}, **kw)["cap"] == 150000)
    check("--steps 90000 → 上限 90000 (引数)", SR.resolve_settings({"stop_rule": True}, **{**kw, "nsteps_arg": 90000})["cap"] == 90000)
    check("引数 requested=False → 使わない (YAML が true でも)", not SR.resolve_settings({"stop_rule": True}, **kw, requested=False)["use"])
    d = SR.resolve_settings({"stop_rule": True, "stop_rule_cap": 151000}, **kw)
    check("上限 151,000 が出力の間隔 2,500 の倍数でない → 使えない (こちらで足した)", not d["use"] and "倍数" in d["reason"], d["reason"])
    d = SR.resolve_settings({"stop_rule": True}, **{**kw, "out_int": 5000})
    check("出力の間隔 5000 → 使えない (窓に 9 点が入らない)", not d["use"] and "出力の間隔" in d["reason"], d["reason"])
    for bad, msg in (({"stop_rule": "yes"}, "true / false"), ({"stop_rule_cap": 0}, "正の整数"), ({"stop_rule_cap": 1.5e5}, "正の整数")):
        ok, e = raises(lambda bad=bad: SR.resolve_settings(bad, **kw), ValueError, msg)
        check(f"YAML の不正値 {bad} → ValueError", ok, e)

    print("== (a) 閾値の YAML (evaluate.stop_rule_thresholds)・到達の後の継続 (evaluate.stop_rule_hold_steps)", flush=True)
    xw45 = 95.2974                 # case/45 の修正後の設計壁の出口 [r_t] (§5.1 #34、終点 95.22758 → 95.29740)
    f45 = [40.0 / xw45, 70.0 / xw45, 94.0 / xw45]
    d = SR.resolve_settings({"stop_rule": True, "stop_rule_thresholds": {"x_fracs": f45, "r_mean_max": 1e-3}}, **kw)
    check("stop_rule_thresholds: x_fracs (case/45 の 40・70・94 r_t) と r_mean_max を上書き、他は既定、出所を記録",
          d["use"] and d["x_fracs"] == f45 and d["r_mean_max"] == 1e-3 and d["r_amp_max"] == SR.R_AMP_MAX
          and d["thresholds_source"] == {"evaluate.stop_rule_thresholds": ["r_mean_max", "x_fracs"]}, str(d["x_fracs"]))
    check("キー無しは閾値の出所が既定・hold_steps 0", SR.resolve_settings({"stop_rule": True}, **kw)["thresholds_source"] == "既定 (取り込み元のまま)"
          and SR.resolve_settings({"stop_rule": True}, **kw)["hold_steps"] == 0)
    d = SR.resolve_settings({"stop_rule": True, "stop_rule_thresholds": {"window": 40000}}, **kw)
    check("window 40000 にすると出力の間隔 2500 でも窓に 9 点以上 (使える)", d["use"] and d["window"] == 40000)
    for bad, msg in (({"stop_rule_thresholds": {"x_frac": [0.5]}}, "未知"), ({"stop_rule_thresholds": {"x_fracs": [0.5, 1.2]}}, "x_fracs"),
                     ({"stop_rule_thresholds": {"x_fracs": []}}, "x_fracs"), ({"stop_rule_thresholds": {"x_fracs": [0.5, 0.5]}}, "x_fracs"),
                     ({"stop_rule_thresholds": {"r_mean_max": -1e-3}}, "r_mean_max"), ({"stop_rule_thresholds": {"min_points": 1}}, "min_points"),
                     ({"stop_rule_thresholds": [0.4]}, "mapping"), ({"stop_rule_hold_steps": -1}, "hold_steps"),
                     ({"stop_rule_hold_steps": 1.5}, "hold_steps")):
        ok, e = raises(lambda bad=bad: SR.resolve_settings({"stop_rule": True, **bad}, **kw), ValueError, msg)
        check(f"不正値 {bad} → ValueError", ok, e[:120])
    # 到達の後の継続 (合成の系列; 間隔 1000、一定なら 21,000 で最初に届く)
    st = SR.first_stop(rows_of(s1k), hold_steps=10000)
    check("hold 10,000: 一定の系列は最初の到達 21,000 の後 10,000 step 保って 31,000 で REACHED",
          st["status"] == "REACHED" and st["step"] == 31000 and st["hold"]["reached"] == [21000] and not st["hold"]["breaks"], str(st and st.get("hold")))
    check("hold 0 は取り込み元と同じ (21,000、hold の記録を付けない)", SR.first_stop(rows_of(s1k), hold_steps=0)["step"] == 21000
          and "hold" not in SR.first_stop(rows_of(s1k), hold_steps=0))
    spike = rows_of(s1k, r=lambda s: 2e-2 if s == 25000 else 0.0)        # 25,000 だけ r = 2 % (平均・振幅とも外れる)
    st = SR.first_stop(spike, hold_steps=10000)
    br = st["hold"]["breaks"] if st else None
    check("hold 10,000: 21,000 で届き 25,000 で外れる → 取り消して記録、窓から外れた後に再び届いた出力から数え直して REACHED",
          st is not None and br and br[0]["start"] == 21000 and br[0]["broken_at"] == 25000 and st["step"] == st["hold"]["start"] + 10000
          and st["hold"]["start"] > 25000, str(st and st["hold"]))
    check("hold が系列より長い → 届かない", SR.first_stop(rows_of(s1k[:35]), hold_steps=20000) is None)

    P = lambda ev: types.SimpleNamespace(evaluate=ev)  # noqa: E731
    eu = TMP / "euler_dummy"
    eu.mkdir(exist_ok=True)
    check("_stop_rule_settings: キー無し・引数 None → None (prepare_info に stop_rule を書かない = 今と同じ)",
          RA._stop_rule_settings(P({"nStepOuter": 60000}), None, str(eu), None) is None)
    check("_stop_rule_settings: evaluate.stop_rule: false → 要求なしの記録 (使わない)",
          (lambda d: d is not None and not d["requested"] and not d["use"])(RA._stop_rule_settings(P({"stop_rule": False}), None, str(eu), None)))
    d = RA._stop_rule_settings(P({"stop_rule": True, "nStepOuter": 200000, "outStepInterval": 2500}), None, str(eu), None)
    check("_stop_rule_settings: stop_rule true・nStepOuter 200000・出力 2500 → 使う (上限 200,000・fallback 200,000)",
          d["use"] and d["cap"] == 200000 and d["fallback_steps"] == 200000 and d["out_step_interval"] == 2500, str(d.get("reason")))
    ok_ev = {"stop_rule": True, "nStepOuter": 200000, "outStepInterval": 2500}
    for name, ev, er, msg in (("euler_ref が無い", ok_ev, None, "euler_ref が無い"),
                              ("euler_ref が見つからない", ok_ev, str(TMP / "no_such_euler"), "見つからない"),
                              ("上限が出力の間隔の倍数でない", {**ok_ev, "stop_rule_cap": 151000}, str(eu), "倍数"),
                              ("出力の間隔が窓に 9 点を入れられない (nStepOuter 48000 を割り切らない 2500 → 16000 に直る)",
                               {"stop_rule": True, "outStepInterval": 2500, "stop_rule_cap": 192000}, str(eu), "9 点が入らない"),
                              ("stop_rule_cap だけ (stop_rule が無い)", {"stop_rule_cap": 150000}, str(eu), "stop_rule_cap"),
                              ("stop_rule_thresholds だけ (stop_rule が無い)", {"stop_rule_thresholds": {"x_fracs": [0.5]}}, str(eu),
                               "stop_rule_thresholds")):
        ok, e = raises(lambda ev=ev, er=er: RA._stop_rule_settings(P(ev), None, er, None), ValueError, msg)
        check(f"_stop_rule_settings: {name} → ValueError (prepare_ns は run dir を作る前に呼ぶ)", ok and e.count(";") == 0, e[:140])
    # prepare_ns の本段の step 数・出力の間隔は関数に出す前の式と同じ (回帰)
    def old_steps(ev, nsteps):
        n = nsteps or int(ev.get("nStepOuter", 48000))
        out_int = int(ev.get("outStepInterval", max(n // 6, 1)))
        if n % out_int:
            out_int = next(n // k for k in (3, 2, 4, 6, 1) if n % k == 0)
        return n, out_int
    cases = [({}, None), ({"nStepOuter": 60000}, None), ({"nStepOuter": 200000, "outStepInterval": 2500}, None),
             ({"nStepOuter": 12001}, None), ({"nStepOuter": 50000, "outStepInterval": 3000}, None), ({}, 9000), ({"outStepInterval": 7}, 100)]
    check("_ns_steps (prepare_ns の本段の step 数・出力の間隔) は関数に出す前の式と同じ",
          all(RA._ns_steps(P(ev), ns)[:2] == old_steps(ev, ns) for ev, ns in cases))


# ====================================================================================== 合成の構造格子 (量の式・偽の forge)
NI, NJ, SCALE, X_E = 41, 13, 0.025, 20.0


def synth_geometry(x0=-5.0, x1=30.0):
    x = np.linspace(x0, x1, NI)
    rw = 1.0 + 0.02 * (x - x0)
    X = np.repeat(x[:, None], NJ, axis=1)
    R = rw[:, None] * np.linspace(0.0, 1.0, NJ)[None, :] ** 0.7      # 壁へ寄せる
    return X, R


def synth_base(x0=-5.0, x1=30.0):
    """流量が各列で 1 (同じ折れ線の式で) になる合成の場 (float64; 偽の forge が float32 で書く前の基準)。"""
    X, R = synth_geometry(x0, x1)
    rw = R[:, -1:]
    y = rw - R; delta = 0.15 * rw
    u = 1700.0 * np.minimum(1.0, np.maximum(y / delta, 0.0) ** (1.0 / 7.0))
    I = np.array([2 * np.pi * np.sum(0.5 * ((u[i] * R[i] * SCALE)[1:] + (u[i] * R[i] * SCALE)[:-1]) * np.diff(R[i] * SCALE)) for i in range(NI)])
    ro = np.repeat((1.0 / I)[:, None], NJ, axis=1)
    T = np.full_like(ro, 300.0)
    return {"ro": ro, "Ux": u, "Uy": np.zeros_like(u), "T": T, "P": ro * 287.0 * T, "sonic": np.sqrt(1.4 * 287.0 * T), "k": np.full_like(ro, 1.0)}


def write_mesh(run_dir: Path, x0=-5.0, x1=30.0):
    X, R = synth_geometry(x0, x1)
    coord = np.c_[X.ravel() * SCALE, R.ravel() * SCALE, np.zeros(NI * NJ)].astype(np.float32)
    with h5py.File(run_dir / "nozzle.h5", "w") as h:
        h["MESH/COORD"] = coord.ravel()


def test_quantities():
    print("== (a) 量の式 (合成の構造格子)", flush=True)
    rd = TMP / "q_fixture"
    rd.mkdir()
    write_mesh(rd)
    (rd / "prepare_info.json").write_text(json.dumps({"mesh": {"ni": NI, "nj": NJ}, "scale_m": SCALE, "x_E": X_E, "Md": 5.0,
                                                       "physical_wall": {"domain": [-5.0, X_E]}}))   # 合成では設計壁の出口 = X_E (断面は従来と同じ)
    B = synth_base()
    B["ro"][-1] *= 1.003                                    # 出口の列の流量を 0.3 % 増やす → r = −0.003
    with h5py.File(rd / "res_1000.h5", "w") as h:
        for k, v in B.items():
            h[f"VALUE/{k}"] = v.ravel().astype(np.float32)
    geo = SR.RunGeometry(rd)
    band = {"x": [-5.0, 30.0], "y_b": [0.25, 0.25]}
    row = SR.snapshot(geo, 1000, band)
    check("合成の場: r = (ṁ_in − ṁ_out)/ṁ_in ≈ −0.003 (出口の列の ρ を 1.003 倍)", abs(row["r"] + 0.003) < 2e-6, f"r {row['r']:.7f}")
    check("θ_r の断面は x/x_E = 0.4・0.7・0.95 に最も近い列",
          all(abs(row[f"x_{SR.theta_key(f)}"] - f * X_E) <= 35.0 / (NI - 1) / 2 + 1e-12 for f in SR.X_FRACS))
    from forge_design.report.nozzle_report import load_field
    F = load_field(rd, "res_1000.h5")
    eta = F["R"][-1] / F["R"][-1, -1]; core = (eta >= 0.05) & (eta <= 0.7)
    m = 2 * np.pi * np.trapezoid(F["V"]["ro"] * F["V"]["Ux"] * F["R"], F["R"], axis=1)
    xs = F["X"][:, 0]
    want = {"exit_core_M": float(F["V"]["M"][-1][core].mean()), "mdot_rt": float(np.median(m[(xs > -2.5) & (xs < X_E)]))}
    check("出口コア M・流量が nozzle_report.load_field からの定義とビット同一",
          row["exit_core_M"] == want["exit_core_M"] and row["mdot_rt"] == want["mdot_rt"], f"{row['exit_core_M']!r} {want['exit_core_M']!r}")
    if not M6_XCHECK.exists():
        skip("θ_r・流量が M6 の cold_xcheck.reduce_fields と一致", f"{M6_XCHECK} が無い")
        return
    spec = importlib.util.spec_from_file_location("m6_cold_xcheck", M6_XCHECK)
    XC = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(XC)
    with h5py.File(rd / "nozzle.h5") as h:
        xy = np.array(h["MESH/COORD"], dtype=float).reshape(-1, 3)[:, :2]
    with h5py.File(rd / "res_1000.h5") as h:
        ro, ux, uy, T, k = (np.array(h["VALUE/" + n], dtype=float) for n in ("ro", "Ux", "Uy", "T", "k"))
    o = XC.reduce_fields(xy, ro, ux, uy, T, k, False, NI, NJ, SCALE, np.array(band["x"]), np.array(band["y_b"]))
    th_m6 = [float(o["theta_r"][int(np.argmin(np.abs(o["x"] - f * X_E)))]) for f in SR.X_FRACS]
    th_new = [row[SR.theta_key(f)] for f in SR.X_FRACS]
    check("θ_r が M6 の reduce_fields (同じ y_b、fac 1.0) とビット同一 (式が同じことの確認で、判定の規則の同等ではない)", th_new == th_m6,
          f"{th_new} vs {th_m6}")
    check("入口・出口の流量が M6 の reduce_fields とビット同一",
          row["mdot_in"] == float(o["mdot_in"]) and row["mdot_out"] == float(o["mdot_out"]))


def test_station_reference():
    """θ_r の断面の基準は設計壁 (物理壁) の出口 prepare_info.physical_wall.domain[1] (取り込み元 60de8125 の修正)。x_E (軸 Mach 則の終点) では
    ない。断面の比は evaluate.stop_rule_thresholds.x_fracs で変えられ、case/45 の試験部の 40・70・94 r_t を指定できる。"""
    print("== (a) 断面の基準 (設計壁の出口) と x_fracs の指定", flush=True)
    x0, xw = -12.5, 95.2974                                  # case/45 の H の入口と修正後の設計壁の出口 [r_t] (間隔 約 2.7 r_t の合成の格子)
    rd = TMP / "station_fixture"
    rd.mkdir()
    write_mesh(rd, x0, xw)
    (rd / "prepare_info.json").write_text(json.dumps({"mesh": {"ni": NI, "nj": NJ}, "scale_m": SCALE, "x_E": 39.81, "Md": 6.0,
                                                       "physical_wall": {"domain": [x0, xw]}}))
    with h5py.File(rd / "res_1000.h5", "w") as h:
        for k, v in synth_base(x0, xw).items():
            h[f"VALUE/{k}"] = v.ravel().astype(np.float32)
    geo = SR.RunGeometry(rd)
    band = {"x": [x0, xw], "y_b": [0.25, 0.25]}
    dx = (xw - x0) / (NI - 1)
    row = SR.snapshot(geo, 1000, band)
    check("断面の基準 x_w = physical_wall.domain[1] (95.2974)、x_E (39.81) ではない",
          geo.x_wall_exit == xw and geo.x_E == 39.81, f"x_w {geo.x_wall_exit}")
    got = [row[f"x_{SR.theta_key(f)}"] for f in SR.X_FRACS]
    check("既定の断面は x/x_w = 0.4・0.7・0.95 に最も近い列 (38.1・66.7・90.5 r_t 付近)",
          all(abs(g - f * xw) <= dx / 2 + 1e-9 for g, f in zip(got, SR.X_FRACS)), str([round(g, 3) for g in got]))
    f45 = [40.0 / xw, 70.0 / xw, 94.0 / xw]
    th = SR.thresholds(x_fracs=f45)
    row = SR.snapshot(geo, 1000, band, th["x_fracs"])
    got = [row[f"x_{SR.theta_key(f)}"] for f in th["x_fracs"]]
    check("x_fracs = 40/x_w・70/x_w・94/x_w → 断面は 40・70・94 r_t に最も近い列 (M6 の §6.11 の断面; 規則は同等でない)",
          all(abs(g - x) <= dx / 2 + 1e-9 for g, x in zip(got, (40.0, 70.0, 94.0))), str([round(g, 3) for g in got]))
    (rd / "prepare_info.json").write_text(json.dumps({"mesh": {"ni": NI, "nj": NJ}, "scale_m": SCALE, "x_E": 39.81}))
    ok, e = raises(lambda: SR.RunGeometry(rd), ValueError, "physical_wall.domain")
    check("physical_wall.domain の無い run は RunGeometry が例外 (断面の基準を決めない)", ok, e[:100])


# ====================================================================================== (b) 偽の forge と本段の見張り
FAKE_SRC = r'''#!{py}
"""偽の forge (試験用): solverConfig の nStepOuter・outStepInterval に従って res_<N>.h5 を書く。fake_scenario.json の指定で
出口の列の ρ を交互に揺らす・途中で異常終了・非有限を書く・detectNaN のように止まる・SIGTERM を無視する・残差の減り方 (res_tau)。"""
import json, re, signal, sys, time
import h5py, numpy as np
cfg = open("solverConfig.yaml").read()
n_outer = int(re.search(r"nStepOuter:\s*(\d+)", cfg).group(1))
interval = int(re.search(r"outStepInterval:\s*(\d+)", cfg).group(1))
sc = json.load(open("fake_scenario.json"))
if sc.get("ignore_term"):
    signal.signal(signal.SIGTERM, signal.SIG_IGN)
base = dict(np.load("fake_base.npz"))
nj = int(sc["nj"])
tau = float(sc.get("res_tau", 5000.0))
cols = ["rms_ro", "rms_roUx", "rms_roUy", "rms_roUz", "rms_roe", "rms_roK", "rms_roOmega"]
hist = open("residual_history.csv", "w")
hist.write("step,inner,phase," + ",".join(cols) + "\n")
def write(name, step, k):
    V = {{key: v.copy() for key, v in base.items()}}
    a = float(sc.get("osc_out", 0.0))
    if a and step >= int(sc.get("osc_from", 0)):
        V["ro"][-nj:] *= 1.0 + a * (-1) ** k
    if sc.get("nan_at") == step:
        V["T"][7] = np.nan
    with h5py.File(name, "w") as h:
        for key, v in V.items():
            h["VALUE/" + key] = v.astype(np.float32)
write("res_0.h5", 0, 0)
step, k = 0, 0
print("fake forge start", n_outer, interval, flush=True)
while step < n_outer:
    nxt = min(step + interval, n_outer)
    time.sleep(float(sc.get("dt", 0.2)))
    for s in range(step + 1, nxt + 1, max(interval // 10, 1)):
        hist.write(f"{{s}},-1,outer_end," + ",".join(f"{{1e-3 * np.exp(-s / tau):.6e}}" for _ in cols) + "\n")
    hist.flush()
    step, k = nxt, k + 1
    if sc.get("exit_at") == step:
        print("fake forge: 異常終了", flush=True); sys.exit(int(sc.get("exit_code", 3)))
    if sc.get("detect_nan_at") == step:
        write(f"res_nan_{{step}}.h5", step, k); print("detectNaN", flush=True); sys.exit(1)
    if sc.get("corrupt_at") == step:
        open(f"res_{{step}}.h5", "wb").write(b"not a hdf5 file")       # 書きかけ (読めない出力)
        if sc.get("exit_after_corrupt"):
            print("fake forge: 書きかけで異常終了", flush=True); sys.exit(3)
        continue
    write(f"res_{{step}}.h5", step, k)
print("fake forge end", flush=True)
sys.exit(0)
'''

EULER_DUMMY = TMP / "euler_dummy"


def make_run(name: str, scenario: dict, *, interval=2500, fallback=10000, cap=40000, euler_ref=True, requested=True, with_sr=True,
             hold=0):
    """偽の forge 用の run (prepare_ns の出力の形: prepare_info.json の stop_rule・solverConfig の nStepOuter = fallback)。
    with_sr=False は指定の無い問題の準備 (prepare_info に stop_rule が無い)。戻り (run_dir, prepare_info の中身)。"""
    rd = TMP / name
    rd.mkdir()
    write_mesh(rd)
    np.savez(rd / "fake_base.npz", **{k: v.ravel() for k, v in synth_base().items()})
    (rd / "fake_scenario.json").write_text(json.dumps({"nj": NJ, **scenario}))
    cfg = RA._cfg_set(TEMPLATE_CFG, {"nStepOuter": fallback, "outStepInterval": interval})
    (rd / "solverConfig.yaml").write_text(cfg)
    (rd / "bcondConfig.yaml").write_text(TEMPLATE_BC)
    (rd / "MESH_QUALITY.txt").write_text("VERDICT: PASS (試験の偽物)\n")
    EULER_DUMMY.mkdir(exist_ok=True)
    info = {"mesh": {"ni": NI, "nj": NJ}, "scale_m": SCALE, "x_E": X_E, "physical_wall": {"domain": [-5.0, X_E]},
            "wall_thermal": {"mode": "isothermal", "Tw": 370.0},
            "euler_ref": str(EULER_DUMMY) if euler_ref else None, "nStepOuter": fallback, "cfl_main": 4.0,
            "throat_physical": {}, "dstar_source": "試験", "initializer": None}
    if with_sr:
        sr = SR.resolve_settings({"stop_rule": requested, "stop_rule_hold_steps": hold}, nsteps_arg=None, n_main=fallback,
                                 n_main_source="試験", out_int=interval, euler_ref=str(EULER_DUMMY), requested=None)
        if sr.get("cap") is not None:
            sr["cap"] = cap
        if not euler_ref:
            sr["euler_ref"] = str(TMP / "no_such_euler")         # 準備の後に euler_ref が消えた run
        info["stop_rule"] = sr
    (rd / "prepare_info.json").write_text(json.dumps(info))
    return rd, info


WATCH = {"poll_s": 0.1, "settle_s": 1.0, "term_wait_s": 3.0, "kill_wait_s": 5.0, "post_wait_s": 120.0}


def fake_band(run_dir, step, euler_ref, out_dir, **kw):
    return {"x": [-5.0, 30.0], "y_b": [0.25, 0.25], "step": int(step), "res_sha256": SR._sha256(Path(run_dir) / f"res_{step}.h5"),
            "euler_ref": str(euler_ref), "extract": {"function": "試験の偽物"}, "out_dir": str(out_dir), "seconds": 0.0}


def outputs(rd):
    return SR.res_steps(rd)


def main_nstep(rd):
    man = json.loads((rd / "stage_manifest.json").read_text())
    st = [s for s in man["stages"] if s["tag"] == "main"]
    return int(st[-1]["soft"]["nStepOuter"]) if st else None


def run_case(rd, stages="none", **kw):
    rc = RA.run_staged_ns(rd, stages=stages, stop_rule_watch=WATCH, **kw)
    rec = json.loads((rd / SR.RECORD).read_text()) if (rd / SR.RECORD).exists() else None
    return rc, rec, SR.read_state(rd)


def state_ok(rc, st, status):
    """戻り値の付加情報と STOP_RULE_STATE.json が同じ状態を持つか。"""
    return (st is not None and st["status"] == status and getattr(rc, "stop_rule", None) is not None
            and rc.stop_rule["status"] == status and rc.stop_rule == {k: v for k, v in st.items() if k not in ("version", "plan")})


def test_watcher():
    print("== (b) 偽の forge で本段の見張り・終了の状態", flush=True)
    real_fix = SR.fix_band
    SR.fix_band = fake_band
    try:
        rd, _ = make_run("reach", {"dt": 0.25})
        rc, rec, st = run_case(rd)
        outs = outputs(rd)
        check("REACHED: 22,500 (間隔 2500 で窓が覆う最初の出力) で止める・rc 0",
              rc == 0 and rec["status"] == "REACHED" and rec["stop_step"] == 22500, f"rc {rc}、{rec['status']} {rec.get('stop_step')}")
        check("REACHED: 状態は STOP_RULE_STATE.json と戻り値の stop_rule で同じ (status・stop_step・forge rc・上限)",
              state_ok(rc, st, "REACHED") and st["stop_step"] == 22500 and st["forge_rc"] == 143 and st["cap"] == 40000
              and st["record"] == SR.RECORD and st["stop_sha256"] == SR._sha256(rd / "res_22500.h5"), json.dumps(st, ensure_ascii=False)[:200])
        check("REACHED: 止めた出力より後の res を消した (残る最後の出力 = 22,500)", max(outs) == 22500, f"消した {rec.get('deleted_after_stop')}")
        check("REACHED: forge だけに SIGTERM (止まったことを確かめた)・forge の rc 143",
              rec["stop_forge"]["stopped"] and rec["stop_forge"]["signals"] == ["SIGTERM"] and rec["forge_rc"] == 143)
        vp = rd / "CONVERGENCE_VERDICT.txt"
        check("REACHED: run_case.sh の後処理 (収束判定) が止めた後に走った", rec["postprocess"]["by"] == "run_case.sh" and vp.exists()
              and "CONVERGENCE_EXIT" in vp.read_text(), str(rec.get("postprocess")))
        check("REACHED: solverConfig と stage_manifest.json の本段の nStepOuter = 上限 40,000",
              RA._cfg_value((rd / "solverConfig.yaml").read_text(), "nStepOuter") == 40000 and main_nstep(rd) == 40000, str(main_nstep(rd)))
        ms = json.loads((rd / RA.MAIN_STAGE_STARTED).read_text())
        check("REACHED: 通知 MAIN_STAGE_STARTED の本段の nStepOuter = 上限 (見張りの前に書いた)",
              ms["main_config"]["nStepOuter"] == 40000 and ms["time_unix"] <= Path(rd / "run_case_stdout.log").stat().st_mtime)

        rd, _ = make_run("reach_ignore_term", {"dt": 0.25, "ignore_term": True}, cap=200000)
        rc, rec, st = run_case(rd)
        check("SIGTERM を無視する forge → 待って SIGKILL で止めた (REACHED、forge の rc 137)",
              rc == 0 and state_ok(rc, st, "REACHED") and rec["stop_forge"]["signals"] == ["SIGTERM", "SIGKILL"]
              and rec["forge_rc"] == 137 and max(outputs(rd)) == 22500, json.dumps(rec.get("stop_forge"), ensure_ascii=False))

        rd, _ = make_run("censored", {"dt": 0.05, "osc_out": 0.01}, cap=30000)
        rc, rec, st = run_case(rd)
        check("CENSORED: 上限 30,000 まで正常に回って届かない — **rc は 0 だが状態は CENSORED** (rc で REACHED と取り違えない)",
              rc == 0 and state_ok(rc, st, "CENSORED") and st["stop_step"] == 30000 and st["forge_rc"] == 0 and max(outputs(rd)) == 30000,
              f"{st and st['status']} rc {rc}")

        rd, _ = make_run("exec_error", {"dt": 0.05, "exit_at": 15000, "exit_code": 3})
        rc, rec, st = run_case(rd)
        check("EXEC_ERROR: 上限より前の異常終了 (forge rc 3 を取り違えない)",
              rc == 3 and state_ok(rc, st, "EXEC_ERROR") and st["forge_rc"] == 3 and st["last_output"] == 12500, f"{st}")

        rd, _ = make_run("diverged", {"dt": 0.25, "nan_at": 12500})
        rc, rec, st = run_case(rd)
        check("DIVERGED: 非有限の出力 12,500 で止めて後の出力を消す (rc ≠ 0)",
              rc != 0 and state_ok(rc, st, "DIVERGED") and st["stop_step"] == 12500 and max(outputs(rd)) == 12500, f"{st}")

        rd, _ = make_run("detect_nan", {"dt": 0.05, "detect_nan_at": 10000})
        rc, rec, st = run_case(rd)
        check("DIVERGED: forge の detectNaN (res_nan_10000.h5 + rc 1) は EXEC_ERROR でなく DIVERGED",
              rc == 1 and state_ok(rc, st, "DIVERGED") and "detectNaN" in (st.get("note") or ""), f"{st}")

        rd, _ = make_run("data_error", {"dt": 0.25, "corrupt_at": 7500})
        rc, rec, st = run_case(rd)
        check("DATA_ERROR: 読めない出力 7,500 を 3 回読み直して forge を止める (rc ≠ 0)",
              rc != 0 and state_ok(rc, st, "DATA_ERROR") and st["stop_step"] == 7500 and rec["tries"].get("7500") == 3, f"{st}")
        rd, _ = make_run("crash_partial", {"dt": 0.05, "corrupt_at": 7500, "exit_after_corrupt": True})
        rc, rec, st = run_case(rd)
        check("書きかけの出力を残して forge が rc 3 で終わる → EXEC_ERROR (DATA_ERROR と取り違えない)",
              rc == 3 and state_ok(rc, st, "EXEC_ERROR") and "書きかけ" in (st.get("note") or ""), f"{st}")

        def no_band(*a, **k):
            raise SR.Unavailable("試験: 抽出が失敗")
        SR.fix_band = no_band
        rd, _ = make_run("unavailable_mid", {"dt": 0.25})
        rc, rec, st = run_case(rd)
        check("見張りの途中で帯の位置が決まらない → UNAVAILABLE、fallback_steps 10,000 の出力で止める — **rc は 0 だが状態は UNAVAILABLE**",
              rc == 0 and state_ok(rc, st, "UNAVAILABLE") and st["stop_step"] == 10000 and max(outputs(rd)) == 10000
              and "抽出が失敗" in st["unavailable"]["reason"], f"{st}")
        SR.fix_band = fake_band

        rd, _ = make_run("no_euler_ref", {"dt": 0.02}, euler_ref=False)
        ok, msg = raises(lambda: run_case(rd), RuntimeError, "euler_ref")
        check("要求したのに euler_ref が見つからない → forge を起動する前に RuntimeError (取り込み元は本段を ns_steps で回す)",
              ok and not (rd / "run_case_stdout.log").exists() and SR.read_state(rd) is None, msg[:120])

        rd, _ = make_run("arg_false", {"dt": 0.02})
        rc, rec, st = run_case(rd, stop_rule=False)
        check("stop_rule=False → 使わない (stop_rule.json なし・本段は prepare の長さ 10,000)・状態 NOT_USED を記録",
              rc == 0 and rec is None and max(outputs(rd)) == 10000 and main_nstep(rd) == 10000 and state_ok(rc, st, "NOT_USED"), f"{st}")

        rd, _ = make_run("yaml_false", {"dt": 0.02}, requested=False)
        rc, rec, st = run_case(rd)
        check("prepare で使わない指定 (evaluate.stop_rule: false) → 今までの経路 (記録・状態なし・本段 10,000)",
              rc == 0 and rec is None and st is None and max(outputs(rd)) == 10000 and getattr(rc, "stop_rule", "x") is None)

        rd, _ = make_run("stale_output", {"dt": 0.02})
        (rd / "res_5000.h5").write_bytes(b"")
        ok, msg = raises(lambda: RA.run_staged_ns(rd, stages="none", stop_rule_watch=WATCH), RuntimeError, "起動しない")
        check("本段の起動前に res_<N> (N > 0) がある → forge を起動しない (状態は書かない = 記録なし)",
              ok and not (rd / "run_case_stdout.log").exists() and SR.read_state(rd) is None, msg[:120])

        # --- 到達の後の継続 (hold_steps) ---
        rd, _ = make_run("hold_ok", {"dt": 0.2}, hold=10000)
        rc, rec, st = run_case(rd)
        check("hold 10,000: 22,500 で届いた後 10,000 step 保って 32,500 で REACHED (後の出力を消す)",
              state_ok(rc, st, "REACHED") and st["stop_step"] == 32500 and st["hold"]["reached"] == [22500] and not st["hold"]["breaks"]
              and max(outputs(rd)) == 32500 and st["hold_steps"] == 10000, f"{st and (st['status'], st['stop_step'], st.get('hold'))}")
        rd, _ = make_run("hold_break", {"dt": 0.05, "osc_out": 0.01, "osc_from": 25000}, hold=10000, cap=40000)
        rc, rec, st = run_case(rd)
        br = (st or {}).get("hold", {}).get("breaks") or []
        check("hold 10,000: 22,500 で届き 25,000 から出口の ρ が ±1 % 振れて外れる → 取り消しを記録、上限 40,000 で CENSORED (rc 0)",
              rc == 0 and state_ok(rc, st, "CENSORED") and br and br[0]["start"] == 22500 and br[0]["broken_at"] == 25000
              and st["stop_step"] == 40000, f"{st and (st['status'], st.get('hold'))}")
        rd, _ = make_run("hold_cut", {"dt": 0.05}, hold=20000, cap=30000)
        rc, rec, st = run_case(rd)
        check("hold 20,000 が上限 30,000 までに終わらない → CENSORED (保持の途中と記録)",
              state_ok(rc, st, "CENSORED") and "保持の途中" in (st.get("note") or ""), f"{st and (st['status'], st.get('note'))}")

        # --- stages="full": S1・S2 には当てず本段だけ (S1・S2 は 3,000 step の通常の起動) ---
        orig = (RA.stage_gate, RA._restart_same_mesh)
        RA.stage_gate, RA._restart_same_mesh = (lambda res, cfg: []), (lambda res, mesh: None)
        try:
            rd, _ = make_run("full_stages", {"dt": 0.1})
            rc, rec, st = run_case(rd, stages="full")
        finally:
            RA.stage_gate, RA._restart_same_mesh = orig
        man = json.loads((rd / "stage_manifest.json").read_text())
        steps = [(s["tag"], int(s["soft"]["nStepOuter"])) for s in man["stages"]]
        check("stages=full: S1_soft・S2_mid は 3,000 step のまま、本段だけ上限 40,000 で見張る (REACHED 22,500)",
              steps == [("S1_soft", 3000), ("S2_mid", 3000), ("main", 40000)] and state_ok(rc, st, "REACHED") and st["stop_step"] == 22500,
              f"{steps} {st and st['status']}")
        check("stages=full: 停止規則の系列は本段の出力だけ (最初の行 2,500・前の段の res_3000 を含まない)",
              rec["rows"][0]["step"] == 2500 and rec["band"]["step"] == 2500)
    finally:
        SR.fix_band = real_fix


# ====================================================================================== (c) 受け渡しの制御
CALLS = []


def fake_extract(prev_run, euler_run, omega=0.5, **kw):
    """抽出の偽物: 呼ばれた引数を記録し、delta_r_next.csv の形だけ書く。"""
    CALLS.append({"run": str(prev_run), "euler": str(euler_run), "omega": omega, **kw})
    od = Path(kw.get("out_dir") or prev_run)
    np.savetxt(od / "delta_r_next.csv", np.c_[[0.0, 1.0], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0, 0]], delimiter=",",
               header=",".join(DL.DELTA_R_NEXT_COLUMNS), comments="")
    return {"omega": omega, "n_held": 0, "massflow": {}}


def test_handoff():
    print("== (c) 受け渡しの制御 (handoff_gate・run_from_table・run_pass・run_pass0_integral・CLI)", flush=True)
    real = {"fix": SR.fix_band, "prep": RA.prepare_ns, "staged": RA.run_staged_ns, "collect": RA.collect,
            "extract": DL.extract_and_merge, "card": DL.check_table_against_card}
    scen = {}

    def fake_prepare(problem, run_dir, **kw):
        s = scen["next"]
        rd, info = make_run(Path(run_dir).name, s["scenario"], **s.get("kw", {}))
        if rd != Path(run_dir):
            raise AssertionError(f"{rd} != {run_dir}")
        s["prepare_kw"] = kw
        return info

    def staged(run_dir, **kw):
        return real["staged"](run_dir, stop_rule_watch=WATCH, **kw)

    SR.fix_band = fake_band
    RA.prepare_ns, RA.run_staged_ns, RA.collect = fake_prepare, staged, (lambda problem, run_dir: {"res_file": "試験"})
    DL.extract_and_merge, DL.check_table_against_card = fake_extract, (lambda problem, table: {"table": str(table), "checked": False})
    eu = str(EULER_DUMMY)
    try:
        def from_table(name, scenario, kw=None, **gate_kw):
            scen["next"] = {"scenario": scenario, "kw": kw or {}}
            CALLS.clear()
            rd = TMP / name
            m = DL.run_from_table("試験.yaml", eu, TMP / "table.csv", rd, stages="none", **gate_kw)
            gate = json.loads((rd / DL.HANDOFF_GATE).read_text())
            return rd, m, gate, list(CALLS)

        fast = {"dt": 0.1, "res_tau": 2000.0}         # 残差が 22,500 step で 4.9 桁下がる (収束の判定 PASS)
        slow = {"dt": 0.1, "res_tau": 50000.0}        # 0.2 桁しか下がらない (NOT CONVERGED)
        rd, m, gate, calls = from_table("h_ok", fast)
        check("REACHED + 収束の判定 PASS → 抽出へ進む (band_select='edge' を明示、fixed_point.json あり)",
              gate["proceed"] and gate["stop_rule"]["status"] == "REACHED" and gate["convergence"]["kind"] == "PASS"
              and len(calls) == 1 and calls[0]["band_select"] == "edge" and calls[0]["run"] == str(rd) and (rd / "fixed_point.json").exists()
              and m["handoff_gate"]["proceed"], f"{gate['blocked_by']} {calls}")
        check("  判定の記録: HANDOFF_GATE.json (状態・収束の判定の行・上書きなし)・HANDOFF_CONVERGENCE.txt (check_convergence --segment の出力)",
              gate["override"] == {"allow_not_reached": False, "allow_unconverged": False} and "--segment" in gate["convergence"]["cmd"]
              and (rd / DL.HANDOFF_CONVERGENCE).exists() and "OVERALL" in (rd / DL.HANDOFF_CONVERGENCE).read_text()
              and gate["extract"]["band_select"] == "edge")
        ok_run = rd

        rd, m, gate, calls = from_table("h_unconv", slow)
        check("REACHED でも収束の判定が不合格 (NOT CONVERGED) → 抽出しない・fixed_point.json・delta_r_next.csv を作らない",
              not gate["proceed"] and gate["stop_rule"]["status"] == "REACHED" and gate["convergence"]["kind"] == "NOT_CONVERGED"
              and not calls and not (rd / "fixed_point.json").exists() and not (rd / "delta_r_next.csv").exists()
              and not m["handoff_gate"]["proceed"], f"{gate['blocked_by']}")
        blocked_run = rd
        rd, m, gate, calls = from_table("h_unconv_allow", slow, allow_unconverged=True)
        check("  allow_unconverged=True (明示の上書き) なら抽出へ進み、上書きを記録",
              gate["proceed"] and gate["override"]["allow_unconverged"] and len(calls) == 1 and calls[0]["band_select"] == "edge")

        rd, m, gate, calls = from_table("h_censored", {**fast, "osc_out": 0.01}, kw={"cap": 30000})
        check("CENSORED (rc 0) + 収束 PASS → 止める", not gate["proceed"] and gate["stop_rule"]["status"] == "CENSORED" and not calls,
              f"{gate['stop_rule']['status']} {gate['convergence']['kind']}")
        rd, m, gate, calls = from_table("h_censored_allow", {**fast, "osc_out": 0.01}, kw={"cap": 30000}, allow_not_reached=True)
        check("  allow_not_reached=True なら CENSORED から進む (上書きを記録)", gate["proceed"] and gate["override"]["allow_not_reached"] and len(calls) == 1)

        SR.fix_band = lambda *a, **k: (_ for _ in ()).throw(SR.Unavailable("試験: 抽出が失敗"))
        rd, m, gate, calls = from_table("h_unavail", fast)
        SR.fix_band = fake_band
        check("UNAVAILABLE (rc 0) → 止める", not gate["proceed"] and gate["stop_rule"]["status"] == "UNAVAILABLE" and not calls)

        rd, m, gate, calls = from_table("h_diverged", {**fast, "nan_at": 12500}, allow_not_reached=True, allow_unconverged=True)
        check("DIVERGED は上書き (allow_not_reached・allow_unconverged) でも止める", not gate["proceed"] and gate["stop_rule"]["status"] == "DIVERGED"
              and not calls, f"{gate['blocked_by']}")
        rd, m, gate, calls = from_table("h_exec", {**fast, "exit_at": 15000}, allow_not_reached=True, allow_unconverged=True)
        check("EXEC_ERROR は上書きでも止める", not gate["proceed"] and gate["stop_rule"]["status"] == "EXEC_ERROR" and not calls)

        rd, m, gate, calls = from_table("h_notused", {**fast, "dt": 0.02, "res_tau": 1000.0}, kw={"with_sr": False})
        check("停止規則を使わない run (指定の無い問題、NOT_USED) → 既定は止める (従来は抽出へ進んだ)",
              not gate["proceed"] and gate["stop_rule"]["status"] == "NOT_USED" and not calls and SR.read_state(rd) is None)
        rd, m, gate, calls = from_table("h_notused_allow", {**fast, "dt": 0.02, "res_tau": 1000.0}, kw={"with_sr": False}, allow_not_reached=True)
        check("  allow_not_reached=True なら NOT_USED から進む (収束の判定 PASS のとき)", gate["proceed"] and len(calls) == 1
              and calls[0]["band_select"] == "edge")

        # MISSING: 停止規則を要求した run に状態の記録が無い (本段が終わっていない・記録が失われた)
        miss = TMP / "h_missing"
        shutil.copytree(ok_run, miss)
        (miss / SR.STATE).unlink()
        g = DL.handoff_gate(miss, allow_not_reached=True, allow_unconverged=True)
        check("記録なし (MISSING) は上書きでも止める", not g["proceed"] and g["stop_rule"]["status"] == "MISSING", f"{g['blocked_by']}")
        # 収束の判定不能 (stage_manifest.json が無い) は上書きでも止める
        und = TMP / "h_undetermined"
        shutil.copytree(ok_run, und)
        (und / "stage_manifest.json").unlink()
        g = DL.handoff_gate(und, allow_not_reached=True, allow_unconverged=True)
        check("収束の判定不能 (--segment に stage_manifest.json が無い) は上書きでも止める", not g["proceed"]
              and g["convergence"]["kind"] == "UNDETERMINED", f"{g['convergence']['verdict']}")

        # --- run_pass: 前の pass の run の関門 (抽出・次の壁の作成の前) ---
        CALLS.clear()
        nxt = TMP / "p_from_blocked"
        ok, msg = raises(lambda: DL.run_pass("試験.yaml", eu, blocked_run, nxt, stages="none"), DL.HandoffBlocked, "前の pass")
        check("run_pass: 前の pass の run が関門で不合格 → 抽出・準備の前に止める (抽出を呼ばない・新しい run dir を作らない)",
              ok and not CALLS and not nxt.exists(), msg[:150])
        scen["next"] = {"scenario": fast}
        CALLS.clear()
        nxt = TMP / "p_from_ok"
        info = DL.run_pass("試験.yaml", eu, ok_run, nxt, stages="none", prepare_only=True)
        pi = json.loads((nxt / "prepare_info.json").read_text())
        check("run_pass: 前の pass の run が合格 → 前の run から band_select='edge' で抽出し、その delta_r_next.csv で準備する",
              len(CALLS) == 1 and CALLS[0]["run"] == str(ok_run) and CALLS[0]["band_select"] == "edge"
              and Path(scen["next"]["prepare_kw"]["delta_r_csv"]) == ok_run / "delta_r_next.csv", f"{CALLS}")
        check("run_pass: 前の run の関門の結果を新しい run の prepare_info.json (handoff_from_prev) に残す (前の run には書かない)",
              pi.get("handoff_from_prev", {}).get("proceed") is True and pi["handoff_from_prev"]["run"] == str(ok_run)
              and info.get("handoff_from_prev") is not None)
        CALLS.clear()
        scen["next"] = {"scenario": slow}
        nxt = TMP / "p_run_unconv"
        m = DL.run_pass("試験.yaml", eu, ok_run, nxt, stages="none")
        check("run_pass: 新しい run の本段の後、関門が不合格なら固定点の抽出 (fixed_point.json) へ進まない (前の run からの抽出 1 回だけ)",
              len(CALLS) == 1 and CALLS[0]["run"] == str(ok_run) and not (nxt / "fixed_point.json").exists()
              and not m["handoff_gate"]["proceed"], f"{CALLS}")

        # --- run_pass0_integral: 同じ後段 ---
        CALLS.clear()
        scen["next"] = {"scenario": fast}
        rd = TMP / "p0_ok"
        m = DL.run_pass0_integral("試験.yaml", eu, rd, None, initializer={"model": "contur"}, stages="none")
        check("run_pass0_integral: REACHED + PASS → band_select='edge' で抽出", m["handoff_gate"]["proceed"] and len(CALLS) == 1
              and CALLS[0]["band_select"] == "edge" and (rd / "fixed_point.json").exists())
        CALLS.clear()
        scen["next"] = {"scenario": {**fast, "osc_out": 0.01}, "kw": {"cap": 30000}}
        rd = TMP / "p0_censored"
        m = DL.run_pass0_integral("試験.yaml", eu, rd, None, initializer={"model": "contur"}, stages="none")
        check("run_pass0_integral: CENSORED → 抽出しない", not m["handoff_gate"]["proceed"] and not CALLS)

        # --- CLI: 止めたら終了コード 3、前の run が不合格なら 3 ---
        (TMP / "table.csv").write_text("x_rt,delta_r\n0,0\n")
        scen["next"] = {"scenario": slow}
        rc = DL.main(["--problem", "試験.yaml", "--euler-ref", eu, "--delta-r-csv", str(TMP / "table.csv"), "--run-dir", str(TMP / "cli_unconv"),
                      "--stages", "none"])
        check("CLI --delta-r-csv: 関門が不合格 → 終了コード 3", rc == 3, f"rc {rc}")
        scen["next"] = {"scenario": slow}
        rc = DL.main(["--problem", "試験.yaml", "--euler-ref", eu, "--delta-r-csv", str(TMP / "table.csv"), "--run-dir", str(TMP / "cli_allow"),
                      "--stages", "none", "--allow-unconverged"])
        check("CLI --allow-unconverged → 進む (終了コード 0)", rc == 0, f"rc {rc}")
        rc = DL.main(["--problem", "試験.yaml", "--euler-ref", eu, "--prev", str(blocked_run), "--run-dir", str(TMP / "cli_prev_blocked"),
                      "--stages", "none"])
        check("CLI (run_pass): 前の pass の run が不合格 → 終了コード 3・run dir を作らない", rc == 3 and not (TMP / "cli_prev_blocked").exists())
    finally:
        SR.fix_band = real["fix"]
        RA.prepare_ns, RA.run_staged_ns, RA.collect = real["prep"], real["staged"], real["collect"]
        DL.extract_and_merge, DL.check_table_against_card = real["extract"], real["card"]


# ====================================================================================== (d) 指定の無い問題の回帰
def test_regression():
    print("== (d) 指定の無い問題の回帰 (prepare_info に stop_rule が無い run)", flush=True)
    rd, _ = make_run("reg_no_sr", {"dt": 0.02}, with_sr=False)
    cfg0 = (rd / "solverConfig.yaml").read_text()
    seen = []
    real_forge = RA.run_forge

    def spy(run_dir):
        seen.append((run_dir / "solverConfig.yaml").read_text())
        return real_forge(run_dir)

    RA.run_forge = spy
    try:
        rc = RA.run_staged_ns(rd, stages="none")
    finally:
        RA.run_forge = real_forge
    check("本段は今までどおり run_forge (run_case.sh) で準備の長さ 10,000 まで、solverConfig は準備のまま",
          len(seen) == 1 and seen[0] == cfg0 and max(outputs(rd)) == 10000 and main_nstep(rd) == 10000)
    check("stop_rule.json・STOP_RULE_STATE.json を書かない、戻りは int と等しい (stop_rule None)",
          not (rd / SR.RECORD).exists() and SR.read_state(rd) is None and rc == 0 and isinstance(rc, int) and rc.stop_rule is None
          and json.dumps({"rc": rc}) == '{"rc": 0}' and f"{rc}" == "0")


def main() -> int:
    print(f"出力: {TMP}", flush=True)
    test_synthetic()
    test_quantities()
    test_station_reference()
    if QUICK:
        skip("(b)〜(d) 偽の forge", "--quick")
    else:
        FAKE_FORGE.write_text(FAKE_SRC.format(py=sys.executable))
        FAKE_FORGE.chmod(FAKE_FORGE.stat().st_mode | stat.S_IXUSR)
        test_watcher()
        test_handoff()
        test_regression()
    print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAIL'}")
    if FAIL == 0:
        shutil.rmtree(TMP, ignore_errors=True)
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
