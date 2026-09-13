#!/usr/bin/env python3
"""⑤ SERN 評価ゲート (plan §5.1 R1) の単体テスト (design/.venv-opt/bin/python で実行)。

山場: (1) `steadiness` が NaN 系列を STEADY にしない (codex C1 の再現バグ)、正式ツール classify と同じ判定、
(2) 保存場・残差・力係数の各ゲートが独立に効く、(3) driver が rc != 0 の評価を採用しない (旧 §4.13-3 撤回)、
degraded が台帳と pareto 要約に残る、物理的 INFEASIBLE と数値失敗を混同しない、(4) 3D 帳簿がノズル力と機体力を分ける。
"""
import csv
import json
import math
import subprocess
import sys
import tempfile
import types
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from forge_design.metrics import sern_gates as G  # noqa: E402
from forge_design.metrics.sern_forces import steadiness, write_force_history_csv  # noqa: E402
from forge_design.opt import driver_sern as D  # noqa: E402
import check_quasisteady as CQ  # noqa: E402  (sern_forces が tools を sys.path に足す)

FAIL = 0


def check(name, cond):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name)
    if not cond:
        FAIL += 1


# --- (1) steadiness -------------------------------------------------------------------
nan = float("nan")
check("steadiness [1,1,1,NaN] は STEADY でない (NONFINITE)", steadiness([1, 1, 1, nan])["verdict"] == "NONFINITE")
check("steadiness NaN 個数を残す", steadiness([1, nan, 1, nan, 1])["n_nonfinite"] == 2)
check("steadiness 平坦 10 点 = STEADY", steadiness([0.95] * 10)["verdict"] == "STEADY")
check("steadiness 3 点 = TRANSIENT-UNSETTLED", steadiness([1, 1, 1])["verdict"] == "TRANSIENT-UNSETTLED")
drift = [1.0 + 0.05 * i for i in range(10)]
check("steadiness 単調 5 %/点 = DRIFTING", steadiness(drift)["verdict"] == "DRIFTING")
osc = [1.0 + 0.04 * (-1) ** i for i in range(12)]
check("steadiness ±4 % 振動 = OSCILLATING (drift 2 %/osc 5 %)", steadiness(osc)["verdict"] == "OSCILLATING")
rng = np.random.default_rng(0)
for k in range(5):
    v = 1.0 + 0.01 * np.cumsum(rng.standard_normal(12))
    a = steadiness(v)["verdict"]; b = CQ.classify(np.arange(12.0), v, 0.4, 0.02, 0.05, 4)[0]
    check(f"steadiness == check_quasisteady.classify (乱数系列 {k}: {a})", a == b)
check("mean/amp を返す", abs(steadiness([0.95] * 10)["mean"] - 0.95) < 1e-12)

# --- (1b) force_history.csv → 正式ツール --series-csv ----------------------------------------
tmp = Path(tempfile.mkdtemp(prefix="sern_gates_"))
hist = [{"step": 500 * (i + 1), "C_T": 0.95, "C_T_with_shear": 0.94 if i < 5 else nan, "C_L": 0.1, "C_M": -1.0} for i in range(6)]
csvp = write_force_history_csv(tmp / "force_history.csv", hist)
rows = list(csv.DictReader(open(csvp)))
check("force_history.csv: step + 存在する列だけ", rows[0].keys() >= {"step", "C_T", "C_T_with_shear", "C_L", "C_M"} and "C_T_friction" not in rows[0])
tool = ROOT.parent / "solver_density_cuda" / "tools" / "check_quasisteady.py"
r = subprocess.run([sys.executable, str(tool), "--series-csv", str(csvp), "--series-cols", "C_T,C_T_with_shear", "--drift", "0.02", "--osc", "0.05"],
                   capture_output=True, text=True)
check("check_quasisteady --series-csv: NaN 列は NONFINITE で exit 1", r.returncode == 1 and "NONFINITE" in r.stdout and "C_T_with_shear" in r.stdout)
r = subprocess.run([sys.executable, str(tool), "--series-csv", str(csvp), "--series-cols", "C_T,C_L,C_M"], capture_output=True, text=True)
check("check_quasisteady --series-csv: 平坦列は STEADY で exit 0", r.returncode == 0 and "-> STEADY" in r.stdout)
r = subprocess.run([sys.executable, str(tool), "--series-csv", str(csvp), "--series-cols", "nope"], capture_output=True, text=True)
check("check_quasisteady --series-csv: 未知列は ERROR (黙って STEADY にしない)", r.returncode != 0 and "missing column" in r.stdout)


# --- (2) 保存場・残差ゲート -----------------------------------------------------------------
def make_run(d: Path, roe_nan=False, p_neg=False, nan_dump=False, resid="plateau", rc=0):
    d.mkdir(parents=True, exist_ok=True)
    n = 50
    with h5py.File(d / "res_6000.h5", "w") as f:
        for k, val in (("ro", 1.0), ("roUx", 2.0), ("roUy", 0.0), ("roUz", 0.0), ("roe", 3.0e5), ("P", 1.0e5), ("T", 300.0)):
            a = np.full(n, val)
            if k == "roe" and roe_nan:
                a[7] = nan
            if k == "P" and p_neg:
                a[3] = -1.0
            f.create_dataset(f"VALUE/{k}", data=a)
    if nan_dump:
        (d / "res_nan_12.h5").write_bytes(b"")
    steps = 6000
    with open(d / "residual_history.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["step", "inner", "phase", "rms_ro", "rms_roUx", "rms_roe"])
        for s in range(0, steps + 1, 100):
            t = s / steps
            if resid == "converged":
                v = 10 ** (-1 - 4 * t)
            elif resid == "plateau":
                v = 10 ** (-1 - 1.5 * min(t * 3, 1))
            elif resid == "rising":
                v = 10 ** (-1 - 2 * min(t * 3, 1)) * (1 + 40 * max(t - 0.8, 0))
            elif resid == "nan":
                v = nan if s > steps * 0.9 else 10 ** (-1 - t)
            w.writerow([s, -1, "outer_end", v, v * 10, v * 100])
    (d / "run_case_stdout.log").write_text(f"[run_case] forge exit={rc}\n")


run_ok = tmp / "ok"; make_run(run_ok)
fh = G.field_health(run_ok)
check("field_health 健全 → ok", fh["ok"] and fh["file"] == "res_6000.h5" and fh["min"]["P"] == 1.0e5)
run_nan = tmp / "roe_nan"; make_run(run_nan, roe_nan=True)
check("field_health roe に NaN → 不合格 (個数付き)", not G.field_health(run_nan)["ok"] and G.field_health(run_nan)["nonfinite"] == {"roe": 1})
run_pneg = tmp / "p_neg"; make_run(run_pneg, p_neg=True)
check("field_health P ≤ 0 → 不合格", not G.field_health(run_pneg)["ok"] and G.field_health(run_pneg)["nonpositive"] == {"P": 1})
run_dump = tmp / "dump"; make_run(run_dump, nan_dump=True)
check("field_health res_nan_*.h5 あり → 不合格", not G.field_health(run_dump)["ok"])
check("field_health res 無し → 不合格", not G.field_health(tmp / "nothing_here_")["ok"])

rh = G.residual_health(run_ok)
check("residual_health プラトー → ok (NaN/rising 無し) だが converged でない", rh["ok"] and not rh["converged"] and "stalled" in rh["verdict"])
run_conv = tmp / "conv"; make_run(run_conv, resid="converged")
check("residual_health 4 桁低下 → converged", G.residual_health(run_conv)["converged"])
run_rise = tmp / "rise"; make_run(run_rise, resid="rising")
rr = G.residual_health(run_rise)
check("residual_health 末尾上昇 → rising 列を挙げて不合格", not rr["ok"] and rr["rising"] == ["rms_ro", "rms_roUx", "rms_roe"])
run_rnan = tmp / "rnan"; make_run(run_rnan, resid="nan")
check("residual_health NaN → nan/不合格", not G.residual_health(run_rnan)["ok"] and G.residual_health(run_rnan)["nan"])
check("residual_health csv 無し → MISSING", G.residual_health(tmp / "nothing_here_")["verdict"] == "MISSING")
check("forge_rc_from_log", G.forge_rc_from_log(run_ok) == 0 and G.forge_rc_from_log(tmp / "nothing_here_") is None)

steady_hist = [{"step": 500 * (i + 1), "C_T": 0.95, "C_T_with_shear": 0.94, "C_L": 0.1, "C_M": -1.0} for i in range(8)]
g = G.evaluate_gates(run_ok, steady_hist, 0)
check("evaluate_gates 全合格 → PASS, objective = C_T_with_shear", g["verdict"] == "PASS" and g["fail_class"] is None and g["objective"] == "C_T_with_shear")
check("evaluate_gates Euler (摩擦なし) の objective = C_T", G.evaluate_gates(run_ok, [{k: v for k, v in h.items() if k != "C_T_with_shear"} for h in steady_hist], 0)["objective"] == "C_T")
check("evaluate_gates rc=1 → FAIL/DIVERGED (力係数が STEADY でも)", G.evaluate_gates(run_ok, steady_hist, 1)["fail_class"] == "DIVERGED")
check("evaluate_gates rc 不明 → DIVERGED", G.evaluate_gates(run_ok, steady_hist, None)["fail_class"] == "DIVERGED")
check("evaluate_gates 場 NaN → DIVERGED", G.evaluate_gates(run_nan, steady_hist, 0)["fail_class"] == "DIVERGED")
check("evaluate_gates 残差 rising → RESIDUAL_RISING", G.evaluate_gates(run_rise, steady_hist, 0)["fail_class"] == "RESIDUAL_RISING")
check("evaluate_gates 残差 NaN → DIVERGED", G.evaluate_gates(run_rnan, steady_hist, 0)["fail_class"] == "DIVERGED")
check("evaluate_gates プラトー + require_residual_pass → NOT_CONVERGED", G.evaluate_gates(run_ok, steady_hist, 0, require_residual_pass=True)["fail_class"] == "NOT_CONVERGED")
check("evaluate_gates 履歴なし → NO_FORCES", G.evaluate_gates(run_ok, [], 0)["fail_class"] == "NO_FORCES")
h_nan = [dict(h) for h in steady_hist]; h_nan[-1]["C_T_with_shear"] = nan
gg = G.evaluate_gates(run_ok, h_nan, 0)
check("evaluate_gates 目的量の末尾 NaN → UNSTEADY (NONFINITE)", gg["fail_class"] == "UNSTEADY" and gg["steadiness"]["series"]["C_T_with_shear"]["verdict"] == "NONFINITE")
h_cl = [dict(h) for h in steady_hist]
for i, h in enumerate(h_cl):
    h["C_L"] = 0.1 + 0.02 * i
check("evaluate_gates C_L ドリフト → UNSTEADY (目的量以外も独立に見る)", G.evaluate_gates(run_ok, h_cl, 0)["unsteady" if False else "fail_class"] == "UNSTEADY")
check("evaluate_gates 理由が全ゲート分積まれる (rc + 場 + 残差)", len(G.evaluate_gates(run_rnan, [], 1)["reasons"]) >= 3)


# --- (3) driver: 採用条件・分類・degraded ----------------------------------------------------------
def fake_runner(scenario: dict):
    """R.prepare / run_staged / collect の差し替え。scenario[op] = list of (rc, gate_verdict, fail_class) for ladder rungs。"""
    m = types.SimpleNamespace()
    calls = {"run": []}

    def prepare(prob, rd, op=None, **kw):
        rd = Path(rd); rd.mkdir(parents=True, exist_ok=False)
        info = {"design": {"L_ramp": scenario.get("L_ramp", 10.0), "warnings": scenario.get("warnings", [])}}
        (rd / "prepare_info.json").write_text(json.dumps(info)); return info

    def run_staged(rd, stages, **kw):
        rd = Path(rd); op = rd.name.split("_", 2)[-1].replace("_retry", ""); rung = 1 if rd.name.endswith("_retry") else 0
        calls["run"].append(rd.name)
        rc, _, _ = scenario["ops"][op][min(rung, len(scenario["ops"][op]) - 1)]
        (rd / "res_6000.h5").write_bytes(b""); (rd / "res_5000.h5").write_bytes(b"")
        return rc

    def collect(prob, rd, out_dir=None, rc=None, require_residual_pass=False):
        rd = Path(rd); op = rd.name.split("_", 2)[-1].replace("_retry", ""); rung = 1 if rd.name.endswith("_retry") else 0
        rc_, verdict, fc = scenario["ops"][op][min(rung, len(scenario["ops"][op]) - 1)]
        st = {k: {"verdict": "STEADY", "detail": ""} for k in ("C_T_with_shear", "C_T", "C_L", "C_M")}
        return {"forge_rc": rc, "C_T": 0.95, "C_T_with_shear": 0.94, "C_L": 0.1, "C_M": scenario.get("C_M", {}).get(op, -1.0), "step": 6000,
                "gates": {"verdict": verdict, "fail_class": fc, "reasons": [fc or ""], "objective": "C_T_with_shear",
                          "residual": {"verdict": "NOT CONVERGED (stalled/plateau)"}, "steadiness": {"series": st}}}

    m.prepare, m.run_staged, m.collect, m.calls = prepare, run_staged, collect, calls
    return m


BASE = {"type": "sern_2d", "dv": {k: {"min": 0.0, "max": 1.0} for k in D.DV_ORDER},
        "spec": {"inflow": {"M_in": 2.5, "p_in": 1e5}, "external": {"p_inf": 5e3}, "H_m": 0.1,
                 "operating_points": [{"name": "m6_on", "weight": 0.6, "external": {}}, {"name": "m4_off", "weight": 0.4, "external": {}}]},
        "gas": {"gamma": 1.4, "cp": 1004.5}, "geometry": {"L_ramp_max": 12.0}, "opt": {"cm_min": -7.0}}


def campaign(name, opt_extra=None):
    d = tmp / name; d.mkdir()
    raw = json.loads(json.dumps(BASE)); raw["opt"].update(opt_extra or {})
    import yaml
    (d / "problem.yaml").write_text(yaml.safe_dump(raw))
    return D.SernCampaign(d / "problem.yaml", d / "camp")


X = [0.5] * 5
# (a) 標準レシピ rc=1 (力 STEADY) → 採用しない → 再試行 rc=1 → FAIL/DIVERGED
c = campaign("a"); D.R = fake_runner({"ops": {"m6_on": [(1, "FAIL", "DIVERGED"), (1, "FAIL", "DIVERGED")], "m4_off": [(0, "PASS", None)]}})
row = c.evaluate(X, "doe_000")
check("driver (a) rc≠0 は力 STEADY でも採用しない → FAIL/DIVERGED", row["status"] == "FAIL" and row["fail_class"] == "DIVERGED")
check("driver (a) 標準 → 再試行の 2 段を回した", D.R.calls["run"] == ["doe_000_m6_on", "doe_000_m6_on_retry"])
check("driver (a) 失敗した作動点の要約が台帳に残る (gate/fail_class)", row["ops"]["m6_on"]["gate"] == "FAIL" and row["ops"]["m6_on"]["gate_fail_class"] == "DIVERGED")
# (b) 標準 rc=1 → 再試行 PASS → PASS + degraded、pareto 要約にも残る
c = campaign("b"); D.R = fake_runner({"ops": {"m6_on": [(1, "FAIL", "DIVERGED"), (0, "PASS", None)], "m4_off": [(0, "PASS", None)]}})
row = c.evaluate(X, "doe_000")
check("driver (b) 再試行で通れば PASS + degraded", row["status"] == "PASS" and row["degraded"] and row["degraded_ops"] == ["m6_on"])
s = c.summary()
check("driver (b) pareto 要約に degraded/tag/ops ゲート要約が残る", s["pareto"][0]["degraded"] and s["pareto"][0]["tag"] == "doe_000" and s["n_degraded"] == 1
      and s["pareto"][0]["ops"]["m6_on"]["gate"] == "PASS")
check("driver (b) C_T_w は目的量 (摩擦込み) の加重平均", abs(row["C_T_w"] - 0.94) < 1e-12)
# (c) ゲート UNSTEADY (rc=0) → FAIL/UNSTEADY
c = campaign("c"); D.R = fake_runner({"ops": {"m6_on": [(0, "FAIL", "UNSTEADY"), (0, "FAIL", "UNSTEADY")], "m4_off": [(0, "PASS", None)]}})
row = c.evaluate(X, "doe_000")
check("driver (c) rc=0 でもゲート不合格なら FAIL/UNSTEADY", row["status"] == "FAIL" and row["fail_class"] == "UNSTEADY")
# (d) 設計不成立 → INFEASIBLE/DESIGN (CFD を回さない)
c = campaign("d"); D.R = fake_runner({"warnings": ["key point outside kernel"], "ops": {"m6_on": [(0, "PASS", None)], "m4_off": [(0, "PASS", None)]}})
row = c.evaluate(X, "doe_000")
check("driver (d) 設計 warning → INFEASIBLE/DESIGN、forge を回さない", row["status"] == "INFEASIBLE" and row["fail_class"] == "DESIGN" and D.R.calls["run"] == [])
# (e) 最終輪郭で L_ramp_max 超過 → INFEASIBLE/L_RAMP_MAX
c = campaign("e"); D.R = fake_runner({"L_ramp": 12.14, "ops": {"m6_on": [(0, "PASS", None)], "m4_off": [(0, "PASS", None)]}})
row = c.evaluate(X, "doe_000")
check("driver (e) 最終輪郭 L_ramp > L_ramp_max → INFEASIBLE/L_RAMP_MAX", row["status"] == "INFEASIBLE" and row["fail_class"] == "L_RAMP_MAX")
c = campaign("e2", {"l_ramp_tol": 0.02}); D.R = fake_runner({"L_ramp": 12.14, "ops": {"m6_on": [(0, "PASS", None)], "m4_off": [(0, "PASS", None)]}})
check("driver (e2) opt.l_ramp_tol 2 % で 12.14 は通る", c.evaluate(X, "doe_000")["status"] == "PASS")
# (f) 作動点別 C_M 窓 → INFEASIBLE/CM_WINDOW (加重平均は窓内でも)
c = campaign("f", {"cm_window": {"m4_off": [-1.0, 1.5]}}); D.R = fake_runner({"C_M": {"m6_on": -6.0, "m4_off": 2.2}, "ops": {"m6_on": [(0, "PASS", None)], "m4_off": [(0, "PASS", None)]}})
row = c.evaluate(X, "doe_000")
check("driver (f) 作動点別 C_M 窓 → INFEASIBLE/CM_WINDOW (加重平均 −2.7 は cm_min −7 内)", row["status"] == "INFEASIBLE" and row["fail_class"] == "CM_WINDOW" and "m4_off" in row["note"])
# (g) 想定外例外は INFEASIBLE にしない
c = campaign("g"); D.R = fake_runner({"ops": {"m6_on": [(0, "PASS", None)], "m4_off": [(0, "PASS", None)]}})
def boom(*a, **k):
    raise OSError("disk")
D.R.collect = boom
row = c.evaluate(X, "doe_000")
check("driver (g) I/O 例外 → FAIL/ERROR (INFEASIBLE と混同しない)", row["status"] == "FAIL" and row["fail_class"] == "ERROR" and "OSError" in row["note"])
# (h) summary の status 内訳と HV (PASS のみ)
c = campaign("h"); D.R = fake_runner({"ops": {"m6_on": [(0, "PASS", None)], "m4_off": [(0, "PASS", None)]}})
c.evaluate(X, "doe_000"); D.R = fake_runner({"ops": {"m6_on": [(1, "FAIL", "DIVERGED"), (1, "FAIL", "DIVERGED")], "m4_off": [(0, "PASS", None)]}}); c.evaluate([0.4] * 5, "doe_001")
s = c.summary()
check("driver (h) summary: PASS 1 / FAIL/DIVERGED 1, HV は PASS のみ", s["n_pass"] == 1 and s["status_counts"] == {"PASS": 1, "FAIL/DIVERGED": 1} and s["hv"] > 0)
check("driver (h) _XF は PASS だけ学習に使う", len(c._XF()[0]) == 1)

# --- (4) 3D 帳簿: ノズル力と機体力の分離 ---------------------------------------------------------------
from forge_design.evaluate import runner_sern3d as R3  # noqa: E402
d3 = tmp / "r3d"; d3.mkdir()
P = R3.PHYS_SERN3D


def quad_file(path, quads, ps):
    """quads: list of 4 点 (x,y,z) の面。CONNE は [5, n0..n3] (3D 壁出力の形)。"""
    xyz = np.array([q for f in quads for q in f], float); conn = []
    for k in range(len(quads)):
        conn += [5, 4 * k, 4 * k + 1, 4 * k + 2, 4 * k + 3]
    with h5py.File(path, "w") as f:
        f.create_dataset("MESH/COORD", data=xyz.ravel()); f.create_dataset("MESH/CONNE", data=np.array(conn, np.int64))
        f.create_dataset("VALUE/Ps", data=np.array(ps, float))


# ランプ: z ∈ [0,1] 幅内 2 面 (p=3), 幅外 1 < z ≤ 2 = 機体 (p=2), 2 < z ≤ 3 = 機体幅の外 (p=2)。各面 x∈[0,1], y=1 (法線 +y)
def sq(z0, z1):
    return [(0, 1, z0), (1, 1, z0), (1, 1, z1), (0, 1, z1)]
quad_file(d3 / f"res_ramp_{P['ramp']}_100.h5", [sq(0, 0.5), sq(0.5, 1), sq(1, 2), sq(2, 3)], [3, 3, 2, 2])
quad_file(d3 / f"res_cowl_in_{P['cowl_in']}_100.h5", [[(0, 0, 0), (1, 0, 0), (1, 0, 1), (0, 0, 1)]], [4])
f3 = R3.forces3d(d3, 100, p_a=1.0, F_ideal_per_m=10.0, half_W=1.0, H=1.0, half_W_vehicle=2.0)
# ノズル: ramp 幅内 (3−1)·(0.5+0.5)=2 (+y), cowl_in (4−1)·1=3 (−y) → Fy = −1; 機体: (2−1)·1 = 1 (+y); beyond: 面積 1 は除外
check("forces3d C_L はノズル面だけ (幅内 ramp + cowl_in)", abs(f3["C_L"] - (-1.0 / 10.0)) < 1e-12)
check("forces3d C_L_vehicle は W/2 < z ≤ W_vehicle/2 の機体下面だけ", abs(f3["C_L_vehicle"] - 0.1) < 1e-12 and abs(f3["area_beyond_vehicle"] - 1.0) < 1e-12)
check("forces3d 旧定義 (total_with_vehicle) は両者の和", abs(f3["C_L_total_with_vehicle"] - 0.0) < 1e-12)
f3b = R3.forces3d(d3, 100, p_a=1.0, F_ideal_per_m=10.0, half_W=1.0, H=1.0, half_W_vehicle=None)
check("forces3d W_vehicle 省略 → 幅外を全部機体に (旧挙動; Z_ext 依存)", abs(f3b["C_L_vehicle"] - 0.2) < 1e-12 and f3b["area_beyond_vehicle"] == 0.0)

print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAILED'}")
sys.exit(1 if FAIL else 0)
