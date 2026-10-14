#!/usr/bin/env python3
"""NS の本段の停止規則 (`evaluate/stop_rule.py`、plan tooling-nozzle-cooled-wall-ns-settling §4.6・§5.1 #1c) の試験 — §6 S-stop (a)(b)。

  (a) 合成の系列 (`judge_rows`・`first_stop`、forge 0 step): 一定 → 窓が 2 万 step を覆い 9 点以上になった最初の出力で REACHED
      (間隔 1000 は 21,000、2500 は 22,500)、θ_r の平均の線形の変化 0.1 %/2 万 step → 届かない、**θ_r が 2 万 step に 0.08 % 単調に増える系列 (r = 0)
      → 9 点・21 点とも届かない (D_θ ≈ 0.07997 %)、換算しない前後半の差 (D_raw) は 0.044982 %・0.041983 % で閾値 0.05 % を下回る = 換算しなければ
      止めてしまう** (諮問 stop-rule の判別 A/B; codex の計算と同じ数値)、振幅 1 %・平均の変化 0 → 届かない、振幅 0.1 %・平均の変化 0 → REACHED、
      r の平均 0.06 % → 届かない (0.04 % なら REACHED)、r の振幅 0.6 % → 届かない、非有限 → DIVERGED、窓の点数 8 → 届かない、
      間隔 3000 (窓に 7 点) → 届かない、場が有限で判定量が NaN → DATA_ERROR、閾値は引数で変えられる。設定 (`resolve_settings`) の既定・上限の出所・
      使えない条件・YAML の不正値。量の式: 合成の構造格子で θ_r・入口と出口の流量が M6 の `cold_xcheck.reduce_fields` (forge-integ-1005、
      読むだけ; 無ければ skip) と一致、出口コア M・流量が `nozzle_report.load_field` からの定義 (`qs_series.py` と同じ式) と一致。
  (b) 偽の forge (`FORGE_BIN` を res を書く台本に差し替え、本物の `run_case.sh` 経由) で `run_staged_ns(stages="none")` の本段の見張り:
      REACHED で止めて後の出力を消す (止めた出力の sha256、forge の終了コード 143、run_case.sh の後処理が走る)、SIGTERM を無視する forge は
      SIGKILL で止める、上限で CENSORED (rc 0)、途中の異常終了で EXEC_ERROR (forge の rc 3 を取り違えない)、非有限で DIVERGED、forge の
      detectNaN (res_nan_* + rc 1) で DIVERGED、読めない出力で DATA_ERROR (書きかけを残した異常終了は EXEC_ERROR)、抽出が使えないとき fallback_steps の出力で止める (UNAVAILABLE)、euler_ref が無い run と
      stop_rule=False は今までの経路 (本段は prepare の長さ)。どの run も solverConfig と stage_manifest.json の本段の nStepOuter が実際の上限。
      帯の位置 (`fix_band`) は偽物に差し替える (本物の抽出は case/70 の replay `stop_rule_replay.py` で通す)。

実行 (forge-deliv の直下で): /home/sano/work/forge/design/.venv-opt/bin/python design/tests/run_stop_rule_tests.py [--quick]
所要: (a) < 5 s、(b) 約 30 s (偽の forge 13 本)。--quick は (b) を省く。出力は一時ディレクトリ (全部合格なら消す)。
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

C70_RUN = ROOT / "case" / "70.m5_air_wt_demo" / "run_0010b_ns_c2b_cont60k"     # solverConfig・bcondConfig の雛形 (読むだけ)
M6_XCHECK = Path("/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/cold_xcheck.py")   # M6 の式 (読むだけ)
FAIL = 0


def check(name, cond, detail=""):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name + (f"  [{detail}]" if detail else ""), flush=True)
    if not cond:
        FAIL += 1


def skip(name, why):
    print(f"skip {name}  [{why}]", flush=True)


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
    print("== (a) 合成の系列", flush=True)
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

    # 0.08 %/2 万 step の単調な系列 (r = 0): 9 点と 21 点 — 諮問の判別 A/B
    for n_pts, steps in ((9, list(range(2500, 22501, 2500))), (21, list(range(1000, 21001, 1000)))):
        th8 = lambda s, s0=steps[0]: 1e-3 * (1.0 + 8e-4 * (s - s0) / 20000.0)  # noqa: E731
        rows = rows_of(steps, theta=th8)
        j = SR.judge_rows(rows)
        t = j["theta"][SR.theta_key(0.7)]
        want_raw = {9: 0.044982e-2, 21: 0.041983e-2}[n_pts]
        check(f"0.08 %/2 万 step の単調な系列 ({n_pts} 点、窓 {j['n_window']} 点) → 届かない (D_θ = {t['D'] * 100:.6f} %)",
              stop_step(rows) is None and j["n_window"] == n_pts and abs(t["D"] - 8e-4 / 1.0004) < 1e-10)
        check(f"  換算しない前後半の差 D_raw = {t['D_raw'] * 100:.6f} % ≤ 0.05 % (換算しなければ止めてしまう; codex の {want_raw * 100:.6f} %)",
              abs(t["D_raw"] - want_raw) < 5e-9 and t["D_raw"] <= SR.THETA_DRIFT_MAX and t["amp"] <= SR.THETA_AMP_MAX,
              f"振幅 {t['amp'] * 100:.6f} %")
        check(f"  閾値は引数で変えられる (theta_drift_max 1e-3 なら {n_pts} 点の系列で REACHED)",
              stop_step(rows, theta_drift_max=1e-3) == ("REACHED", steps[-1]))

    sinus = lambda a: (lambda s: 1e-3 * (1.0 + a * math.sin(2 * math.pi * ((s // 1000) - 1) / 10.0)))  # noqa: E731
    check("振幅 1 %・平均の変化 0 (周期 1 万 step の正弦) → 届かない", stop_step(rows_of(s1k, theta=sinus(1e-2))) is None,
          f"振幅 {SR.judge_rows(rows_of(s1k, theta=sinus(1e-2)), 20)['theta'][SR.theta_key(0.4)]['amp'] * 100:.4f} %")
    j = SR.judge_rows(rows_of(s1k, theta=sinus(1e-3)), 20)
    check("振幅 0.1 %・平均の変化 0 → 21,000 で REACHED", stop_step(rows_of(s1k, theta=sinus(1e-3))) == ("REACHED", 21000),
          f"D {j['theta'][SR.theta_key(0.4)]['D']:.3e}、振幅 {j['theta'][SR.theta_key(0.4)]['amp']:.3e}")
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

    print("== (a) 設定 (resolve_settings)", flush=True)
    kw = dict(nsteps_arg=None, n_main=60000, n_main_source="既定 60000", out_int=1000, euler_ref="case/x/run_euler")
    d = SR.resolve_settings({}, **kw)
    check("既定: 使う・上限 200,000・閾値は §4.6", d["use"] and d["requested"] and d["cap"] == 200000 and d["cap_source"].startswith("既定")
          and d["window"] == 20000 and d["min_points"] == 9 and d["r_mean_max"] == 5e-4 and d["theta_drift_max"] == 5e-4
          and d["x_fracs"] == [0.4, 0.7, 0.95] and d["fallback_steps"] == 60000, json.dumps(d, ensure_ascii=False)[:200])
    check("--steps 90000 → 上限 90000 (引数)", SR.resolve_settings({}, **{**kw, "nsteps_arg": 90000, "n_main": 90000})["cap"] == 90000)
    check("evaluate.stop_rule_cap 150000 → 上限 150000", SR.resolve_settings({"stop_rule_cap": 150000}, **kw)["cap"] == 150000)
    check("evaluate.stop_rule: false → 使わない", not SR.resolve_settings({"stop_rule": False}, **kw)["use"])
    check("引数 requested=False → 使わない (YAML が true でも)", not SR.resolve_settings({"stop_rule": True}, **kw, requested=False)["use"])
    d = SR.resolve_settings({}, **{**kw, "euler_ref": None})
    check("euler_ref が無い → 使えない (理由を記録)", d["requested"] and not d["use"] and "euler_ref" in d["reason"], d["reason"])
    d = SR.resolve_settings({}, **{**kw, "out_int": 5000})
    check("出力の間隔 5000 → 使えない (窓に 9 点が入らない)", not d["use"] and "出力の間隔" in d["reason"], d["reason"])
    check("出力の間隔 2500 → 使える", SR.resolve_settings({}, **{**kw, "out_int": 2500})["use"])
    for bad, msg in (({"stop_rule": "yes"}, "true / false"), ({"stop_rule_cap": 0}, "正の整数"), ({"stop_rule_cap": 1.5e5}, "正の整数")):
        try:
            SR.resolve_settings(bad, **kw); ok, e = False, "例外にならない"
        except ValueError as ex:
            ok, e = msg in str(ex), str(ex)
        check(f"YAML の不正値 {bad} → ValueError", ok, e)


# ====================================================================================== 合成の構造格子 (量の式・偽の forge)
NI, NJ, SCALE, X_E = 41, 13, 0.025, 20.0


def synth_geometry():
    x = np.linspace(-5.0, 30.0, NI)
    rw = 1.0 + 0.02 * (x + 5.0)
    X = np.repeat(x[:, None], NJ, axis=1)
    R = rw[:, None] * np.linspace(0.0, 1.0, NJ)[None, :] ** 0.7      # 壁へ寄せる
    return X, R


def synth_base():
    """流量が各列で 1 (同じ折れ線の式で) になる合成の場 (float64; 偽の forge が float32 で書く前の基準)。"""
    X, R = synth_geometry()
    rw = R[:, -1:]
    y = rw - R; delta = 0.15 * rw
    u = 1700.0 * np.minimum(1.0, np.maximum(y / delta, 0.0) ** (1.0 / 7.0))
    I = np.array([2 * np.pi * np.sum(0.5 * ((u[i] * R[i] * SCALE)[1:] + (u[i] * R[i] * SCALE)[:-1]) * np.diff(R[i] * SCALE)) for i in range(NI)])
    ro = np.repeat((1.0 / I)[:, None], NJ, axis=1)
    T = np.full_like(ro, 300.0)
    return {"ro": ro, "Ux": u, "Uy": np.zeros_like(u), "T": T, "P": ro * 287.0 * T, "sonic": np.sqrt(1.4 * 287.0 * T), "k": np.full_like(ro, 1.0)}


def write_mesh(run_dir: Path):
    X, R = synth_geometry()
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
          all(abs(row[f"x_{SR.theta_key(f)}"] - f * X_E) <= 35.0 / (NI - 1) / 2 + 1e-12 for f in SR.X_FRACS),
          str([round(row[f"x_{SR.theta_key(f)}"], 4) for f in SR.X_FRACS]))
    # 記録の量: nozzle_report.load_field と qs_series.py の式
    from forge_design.report.nozzle_report import load_field
    F = load_field(rd, "res_1000.h5")
    eta = F["R"][-1] / F["R"][-1, -1]; core = (eta >= 0.05) & (eta <= 0.7)
    m = 2 * np.pi * np.trapezoid(F["V"]["ro"] * F["V"]["Ux"] * F["R"], F["R"], axis=1)
    xs = F["X"][:, 0]
    want = {"exit_core_M": float(F["V"]["M"][-1][core].mean()), "mdot_rt": float(np.median(m[(xs > -2.5) & (xs < X_E)]))}
    check("出口コア M・流量が nozzle_report.load_field からの定義 (qs_series.py) とビット同一",
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
    check("θ_r が M6 の reduce_fields (同じ y_b、fac 1.0) とビット同一", th_new == th_m6, f"{th_new} vs {th_m6}")
    check("入口・出口の流量が M6 の reduce_fields とビット同一",
          row["mdot_in"] == float(o["mdot_in"]) and row["mdot_out"] == float(o["mdot_out"]), f"{row['mdot_in']!r} {float(o['mdot_in'])!r}")


# ====================================================================================== (b) 偽の forge と本段の見張り
FAKE_SRC = r'''#!{py}
"""偽の forge (試験用): solverConfig の nStepOuter・outStepInterval に従って res_<N>.h5 を書く。fake_scenario.json の指定で
出口の列の ρ を交互に揺らす・途中で異常終了・非有限を書く・detectNaN のように止まる・SIGTERM を無視する。"""
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
cols = ["rms_ro", "rms_roUx", "rms_roUy", "rms_roUz", "rms_roe", "rms_roK", "rms_roOmega"]
hist = open("residual_history.csv", "w")
hist.write("step,inner,phase," + ",".join(cols) + "\n")
def write(name, step, k):
    V = {{key: v.copy() for key, v in base.items()}}
    a = float(sc.get("osc_out", 0.0))
    if a:
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
        hist.write(f"{{s}},-1,outer_end," + ",".join(f"{{1e-3 * np.exp(-s / 5000.0):.6e}}" for _ in cols) + "\n")
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


def make_run(name: str, scenario: dict, *, interval=2500, fallback=10000, cap=40000, euler_ref="EULER", use=True, requested=True):
    """偽の forge 用の run (prepare_ns の出力の形: prepare_info.json の stop_rule・solverConfig の nStepOuter = fallback)。"""
    rd = TMP / name
    rd.mkdir()
    write_mesh(rd)
    np.savez(rd / "fake_base.npz", **{k: v.ravel() for k, v in synth_base().items()})
    (rd / "fake_scenario.json").write_text(json.dumps({"nj": NJ, **scenario}))
    cfg = (C70_RUN / "solverConfig.yaml").read_text()
    cfg = RA._cfg_set(cfg, {"nStepOuter": fallback, "outStepInterval": interval})
    (rd / "solverConfig.yaml").write_text(cfg)
    shutil.copy2(C70_RUN / "bcondConfig.yaml", rd / "bcondConfig.yaml")
    eu = TMP / "euler_dummy"
    eu.mkdir(exist_ok=True)
    sr = SR.resolve_settings({}, nsteps_arg=None, n_main=fallback, n_main_source="試験", out_int=interval,
                             euler_ref=(str(eu) if euler_ref else None), requested=requested)
    sr["cap"] = cap if sr.get("cap") is not None else None
    (rd / "prepare_info.json").write_text(json.dumps({"mesh": {"ni": NI, "nj": NJ}, "scale_m": SCALE, "x_E": X_E,
                                                       "physical_wall": {"domain": [-5.0, X_E]},
                                                       "wall_thermal": {"mode": "isothermal", "Tw": 350.0}, "stop_rule": sr,
                                                       "euler_ref": sr["euler_ref"]}))
    return rd


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


def run_case(rd, **kw):
    rc = RA.run_staged_ns(rd, stages="none", stop_rule_watch=WATCH, **kw)
    rec = json.loads((rd / SR.RECORD).read_text()) if (rd / SR.RECORD).exists() else None
    return rc, rec


def test_watcher():
    print("== (b) 偽の forge で本段の見張り", flush=True)
    FAKE_FORGE.write_text(FAKE_SRC.format(py=sys.executable))
    FAKE_FORGE.chmod(FAKE_FORGE.stat().st_mode | stat.S_IXUSR)
    real_fix = SR.fix_band
    SR.fix_band = fake_band
    try:
        # --- REACHED: 一定の場 → 22,500 で止めて後の出力を消す ---
        rd = make_run("reach", {"dt": 0.25})
        rc, rec = run_case(rd)
        outs = outputs(rd)
        check("REACHED: 22,500 (間隔 2500 で窓が覆う最初の出力) で止める・rc 0",
              rc == 0 and rec["status"] == "REACHED" and rec["stop_step"] == 22500, f"rc {rc}、{rec['status']} {rec.get('stop_step')}")
        check("REACHED: 止めた出力より後の res を消した (残る最後の出力 = 22,500)、境界の res_* も",
              max(outs) == 22500 and not [f.name for f in rd.iterdir() if SR._OUT_RE.fullmatch(f.name) and int(SR._OUT_RE.fullmatch(f.name).group(1)) > 22500],
              f"消した {rec.get('deleted_after_stop')}、残り {outs[-3:]}")
        check("REACHED: 止めた出力の sha256 を記録", rec.get("stop_sha256") == SR._sha256(rd / "res_22500.h5"))
        check("REACHED: forge だけに SIGTERM (止まったことを確かめた)・forge の rc 143",
              rec["stop_forge"]["stopped"] and rec["stop_forge"]["signals"] == ["SIGTERM"] and rec["forge_rc"] == 143,
              json.dumps(rec.get("stop_forge"), ensure_ascii=False))
        vp = rd / "CONVERGENCE_VERDICT.txt"
        check("REACHED: run_case.sh の後処理 (収束判定) が止めた後に走った", rec["postprocess"]["by"] == "run_case.sh" and vp.exists()
              and "CONVERGENCE_EXIT" in vp.read_text(), str(rec.get("postprocess")))
        check("REACHED: solverConfig と stage_manifest.json の本段の nStepOuter = 上限 40,000",
              RA._cfg_value((rd / "solverConfig.yaml").read_text(), "nStepOuter") == 40000 and main_nstep(rd) == 40000, str(main_nstep(rd)))
        check("REACHED: 系列・判定の内訳・y_b の出所・記録だけの量を残す",
              len(rec["rows"]) >= 9 and rec["rows"][-1]["judge"]["ok"] and rec["band"]["step"] == 2500
              and "exit_core_M" in rec["recorded_at_stop"] and (rd / "stop_rule.json.tmp").exists() is False)

        # --- SIGTERM を無視する forge → SIGKILL ---
        rd = make_run("reach_ignore_term", {"dt": 0.25, "ignore_term": True}, cap=200000)   # 上限まで回ると term_wait より長い
        rc, rec = run_case(rd)
        check("SIGTERM を無視する forge → 待って SIGKILL で止めた (REACHED、forge の rc 137)",
              rc == 0 and rec["status"] == "REACHED" and rec["stop_forge"]["signals"] == ["SIGTERM", "SIGKILL"] and rec["stop_forge"]["stopped"]
              and rec["forge_rc"] == 137 and max(outputs(rd)) == 22500, json.dumps(rec.get("stop_forge"), ensure_ascii=False))

        # --- CENSORED: 出口の列の ρ を ±1 % 交互 → r の振幅 1 % で届かない、上限 30,000 で正常終了 ---
        rd = make_run("censored", {"dt": 0.05, "osc_out": 0.01}, cap=30000)
        rc, rec = run_case(rd)
        check("CENSORED: 上限 30,000 まで正常に回って届かない (rc 0、forge rc 0)",
              rc == 0 and rec["status"] == "CENSORED" and rec["stop_step"] == 30000 and rec["forge_rc"] == 0 and max(outputs(rd)) == 30000,
              f"{rec['status']} {rec.get('stop_step')} rc {rc} forge {rec.get('forge_rc')}; 最後の判定 {rec['judge_last'].get('why')}")
        check("CENSORED: 本段の nStepOuter (solverConfig・manifest) = 上限 30,000", main_nstep(rd) == 30000)

        # --- EXEC_ERROR: 15,000 で異常終了 (rc 3) ---
        rd = make_run("exec_error", {"dt": 0.05, "exit_at": 15000, "exit_code": 3})
        rc, rec = run_case(rd)
        check("EXEC_ERROR: 上限より前の異常終了 (forge rc 3 を取り違えない)",
              rc == 3 and rec["status"] == "EXEC_ERROR" and rec["forge_rc"] == 3 and rec["last_output"] == 12500,
              f"{rec['status']} rc {rc} forge {rec.get('forge_rc')} {rec.get('note')}")

        # --- DIVERGED: 12,500 に NaN ---
        rd = make_run("diverged", {"dt": 0.25, "nan_at": 12500})
        rc, rec = run_case(rd)
        check("DIVERGED: 非有限の出力 12,500 で止めて後の出力を消す (rc ≠ 0)",
              rc != 0 and rec["status"] == "DIVERGED" and rec["stop_step"] == 12500 and max(outputs(rd)) == 12500
              and rec.get("stop_sha256") == SR._sha256(rd / "res_12500.h5"), f"{rec['status']} {rec.get('stop_step')} rc {rc} {rec.get('note')}")

        # --- DIVERGED (forge の detectNaN: res_nan_ + rc 1) ---
        rd = make_run("detect_nan", {"dt": 0.05, "detect_nan_at": 10000})
        rc, rec = run_case(rd)
        check("DIVERGED: forge の detectNaN (res_nan_10000.h5 + rc 1) は EXEC_ERROR でなく DIVERGED",
              rc == 1 and rec["status"] == "DIVERGED" and "detectNaN" in rec.get("note", ""), f"{rec['status']} rc {rc} {rec.get('note')}")

        # --- DATA_ERROR: 読めない出力 (forge は続く) → 3 回読み直して止める ---
        rd = make_run("data_error", {"dt": 0.25, "corrupt_at": 7500})
        rc, rec = run_case(rd)
        check("DATA_ERROR: 読めない出力 7,500 を 3 回読み直して forge を止める (rc ≠ 0)",
              rc != 0 and rec["status"] == "DATA_ERROR" and rec["stop_step"] == 7500 and rec["tries"].get("7500") == 3
              and rec["stop_forge"]["stopped"], f"{rec['status']} {rec.get('stop_step')} rc {rc} {rec.get('note')}")
        # --- 書きかけの出力を残して異常終了 → DATA_ERROR でなく EXEC_ERROR ---
        rd = make_run("crash_partial", {"dt": 0.05, "corrupt_at": 7500, "exit_after_corrupt": True})
        rc, rec = run_case(rd)
        check("書きかけの出力を残して forge が rc 3 で終わる → EXEC_ERROR (DATA_ERROR と取り違えない)",
              rc == 3 and rec["status"] == "EXEC_ERROR" and rec["forge_rc"] == 3 and "書きかけ" in rec.get("note", ""),
              f"{rec['status']} rc {rc} {rec.get('note')}")

        # --- 抽出が使えない → fallback_steps (10,000) の出力で止める ---
        def no_band(*a, **k):
            raise SR.Unavailable("試験: 抽出が失敗")
        SR.fix_band = no_band
        rd = make_run("unavailable_mid", {"dt": 0.25})
        rc, rec = run_case(rd)
        check("抽出が使えない → UNAVAILABLE、fallback_steps 10,000 の出力で止めて後を消す (rc 0)",
              rc == 0 and rec["status"] == "UNAVAILABLE" and rec["stop_step"] == 10000 and max(outputs(rd)) == 10000
              and "抽出が失敗" in rec["unavailable"]["reason"], f"{rec['status']} {rec.get('stop_step')} rc {rc}")
        SR.fix_band = fake_band

        # --- euler_ref が無い run: 今までの経路 (本段は prepare の長さ 10,000)、UNAVAILABLE を記録 ---
        rd = make_run("no_euler_ref", {"dt": 0.02}, euler_ref=None)
        rc, rec = run_case(rd)
        check("euler_ref が無い → 停止規則を使わず本段は prepare の長さ 10,000 (UNAVAILABLE と理由を記録、manifest も 10,000)",
              rc == 0 and rec["status"] == "UNAVAILABLE" and "euler_ref" in rec["reason"] and max(outputs(rd)) == 10000
              and main_nstep(rd) == 10000 and rec.get("forge_rc") == 0, f"{rec['status']} {rec.get('reason')}")

        # --- stop_rule=False (関数の引数): 今までの経路、記録なし ---
        rd = make_run("arg_false", {"dt": 0.02})
        rc, rec = run_case(rd, stop_rule=False)
        check("stop_rule=False → 使わない (stop_rule.json なし・本段は prepare の長さ 10,000)",
              rc == 0 and rec is None and max(outputs(rd)) == 10000 and main_nstep(rd) == 10000)

        # --- prepare で使わない指定 (evaluate.stop_rule: false) → 記録なし ---
        rd = make_run("yaml_false", {"dt": 0.02}, requested=False)
        rc, rec = run_case(rd)
        check("prepare で使わない指定 → 今までの経路 (記録なし・本段 10,000)", rc == 0 and rec is None and max(outputs(rd)) == 10000)

        # --- 起動前に出力が残っている run は起動しない ---
        rd = make_run("stale_output", {"dt": 0.02})
        (rd / "res_5000.h5").write_bytes(b"")
        try:
            RA.run_staged_ns(rd, stages="none", stop_rule_watch=WATCH); ok, msg = False, "例外にならない"
        except RuntimeError as e:
            ok, msg = "起動しない" in str(e) and not (rd / "run_case_stdout.log").exists(), str(e)
        check("本段の起動前に res_<N> (N > 0) がある → forge を起動しない", ok, msg[:120])
    finally:
        SR.fix_band = real_fix


def test_station_reference():
    """θ_r の断面の基準は設計壁 (物理壁) の出口 (prepare_info.physical_wall.domain[1])。prepare_info.x_E (軸 Mach 則の終点) ではない
    (2026-10-14 M6ノズルの設計・forge-25 の指摘: case/45 で x_E 39.81 r_t を基準にして断面が 15.9・27.9・37.8 r_t になっていた)。"""
    import json as _json
    d = ROOT / "case" / "70.m5_air_wt_demo" / "run_0010b_ns_c2b_cont60k"
    if not (d / "nozzle.h5").exists():
        skip("断面の基準", f"{d} の nozzle.h5 が無い")
        return
    info = _json.loads((d / "prepare_info.json").read_text())
    geo = SR.RunGeometry(d)
    xw = float(info["physical_wall"]["domain"][1])
    check("θ_r の断面の基準は設計壁の出口 (physical_wall.domain[1]) で、x_E (軸 Mach 則の終点) ではない",
          geo.x_wall_exit == xw and abs(xw - float(info["x_E"])) > 1.0, f"x_w {xw:.4f}、x_E {float(info['x_E']):.4f} r_t")
    xs = [f * geo.x_wall_exit for f in SR.X_FRACS]
    check("断面は x/x_w = 0.4・0.7・0.95 (case/70 c2b で 19.2・33.6・45.6 r_t 付近)", all(abs(a - b) < 1e-12 for a, b in zip(xs, [0.4 * xw, 0.7 * xw, 0.95 * xw])),
          [round(x, 3) for x in xs])


def main() -> int:
    print(f"出力: {TMP}", flush=True)
    test_synthetic()
    test_quantities()
    test_station_reference()
    if QUICK:
        skip("(b) 偽の forge", "--quick")
    elif not C70_RUN.exists():
        skip("(b) 偽の forge", f"{C70_RUN} が無い (solverConfig の雛形)")
    else:
        test_watcher()
    print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAIL'}")
    if FAIL == 0:
        shutil.rmtree(TMP, ignore_errors=True)
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
