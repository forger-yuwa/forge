#!/usr/bin/env python3
"""CFD ピン抽出器 (`feedback/cfd_initial_line.py`) の単体試験 = plan tooling-nozzle-cfd-pinned-initial-line §6 V0 (事前登録)。

case/45 の Euler 格子に Hall 場を載せた**合成 Hall 場**から、Hall の線とアンカーを再現できるかを見る:
  線 max|ΔM| ≤ 2e-4・max|Δθ| ≤ 0.004°・|Δx₀| ≤ 2e-4、アンカー |ΔM| ≤ 2e-5・|ΔM′| ≤ 1e-4・|ΔM″| ≤ 1e-3、
  m* の再現 ≤ 1e-4 (相対)、ν↔M 往復 ≤ 2e-4。
加えて入力の拒否 (亜音速の場・アンカー位置の取り違え) と design_chain のキー検査。
(試作 case/45.isobutane_m6_d155/test_cfd_initial_line.py の移植。)

usage: design/.venv-opt/bin/python design/tests/run_cfd_initial_line_tests.py [格子の run]
  格子の run の既定: $FORGE_CFDPIN_GRID_RUN か /home/sano/work/forge/case/45.isobutane_m6_d155/run_0062_euler_wallfit_fit_r1_ext6k
  (検証 run は git 管理外。無ければ exit 2 で止める — 合格扱いにしない)
"""
import json
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "design"))

from forge_design.evaluate.runner_axismach import _gam_or_gas, design_chain, load_problem  # noqa: E402
from forge_design.feedback.cfd_initial_line import CFDPinnedThroat  # noqa: E402
from forge_design.geometry.moc_inverse import _flux_along, _Pt  # noqa: E402
from forge_design.geometry.moc_kernel import pm_mach, pm_nu  # noqa: E402
from forge_design.geometry.transonic import HallThroat  # noqa: E402

FAIL = 0


def check(name, cond):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name)
    if not cond:
        FAIL += 1


PROB = ROOT / "case/45.isobutane_m6_d155/problem_d155_euler_c2final_n2400.yaml"
run = Path(sys.argv[1] if len(sys.argv) > 1 else os.environ.get(
    "FORGE_CFDPIN_GRID_RUN", "/home/sano/work/forge/case/45.isobutane_m6_d155/run_0062_euler_wallfit_fit_r1_ext6k"))
if not (run / "nozzle.h5").exists():
    print(f"格子の run が無い: {run} (引数か FORGE_CFDPIN_GRID_RUN で与える)")
    sys.exit(2)

import h5py  # noqa: E402

info = json.loads((run / "prepare_info.json").read_text())
S = float(info["scale_m"])
ni = int(info["mesh"]["ni"])
with h5py.File(run / "nozzle.h5") as f:
    nc = f["/MESH/COORD"][:].reshape(-1, 3)
X = (nc[:, 0] / S).reshape(ni, -1)
R = (nc[:, 1] / S).reshape(ni, -1)

# --- V0: 合成 Hall 場 ------------------------------------------------------------
g = 1.2735422978978068
ht = HallThroat(R=2.0, gamma=g)
m = (X[:, 0] > -2.2) & (X[:, 0] < 2.7)
Mh = np.ones(X.shape)
Th = np.zeros(X.shape)
Mh[m] = ht.mach(X[m], R[m])
Th[m] = ht.theta(X[m], R[m])
cfd = CFDPinnedThroat(2.0, g, X, R, Mh, Th)
xc, rc, Mc, tc = cfd.throat_characteristic(n=41)
xh, rh, Mhl, thl = ht.throat_characteristic(n=41)
aC = cfd.axis_anchor(cfd.x0_cfd)
aH = ht.axis_anchor(cfd.x0_cfd)
p = load_problem(PROB)
gas = _gam_or_gas(p)
gg = gas if hasattr(gas, "nu") else g


def flux(x, r, M, t):
    return _flux_along([_Pt(float(a), float(b), float(d), float(pm_nu(float(c), gg)), gg)
                        for a, b, c, d in zip(x, r, M, t)], gg)[-1]


mC, mH = flux(xc, rc, Mc, tc), flux(xh, rh, Mhl, thl)
Ms = np.linspace(1.0005, 1.3, 200)
rt = np.array([float(pm_mach(float(pm_nu(float(v), gg)), gg)) for v in Ms])
out = dict(line_dM=float(np.abs(Mc - Mhl).max()), line_dtheta_deg=float(np.degrees(np.abs(tc - thl)).max()),
           dx0=float(xc[0] - xh[0]), anchor_dM=aC[0] - aH[0], anchor_dMp=aC[1] - aH[1], anchor_dMpp=aC[2] - aH[2],
           mstar_rel=float(mC / mH - 1), nu_roundtrip=float(np.abs(rt - Ms).max()))
gates = dict(line_dM=2e-4, line_dtheta_deg=0.004, dx0=2e-4, anchor_dM=2e-5, anchor_dMp=1e-4, anchor_dMpp=1e-3,
             mstar_rel=1e-4, nu_roundtrip=2e-4)
for k, v in gates.items():
    check(f"V0 {k}: {out[k]:.3e} (≤ {v:g})", abs(out[k]) <= v)

# --- 入力の拒否 --------------------------------------------------------------------
try:
    CFDPinnedThroat(2.0, g, X, R, np.full(X.shape, 0.9), np.zeros(X.shape))
    check("一様 M=0.9 (亜音速) の場を拒否", False)
except ValueError as e:
    check(f"一様 M=0.9 (亜音速) の場を拒否 ({e})", "亜音速" in str(e))
try:
    cfd.axis_anchor(cfd.x0_cfd + 0.01)
    check("軸着地 x0 以外でのアンカー要求を拒否", False)
except ValueError:
    check("軸着地 x0 以外でのアンカー要求を拒否", True)

# --- design_chain のキー -----------------------------------------------------------
for key, val in (("initial_line", "spline"), ("wall_repr", "smooth")):
    pp = load_problem(PROB)
    pp.geometry[key] = val
    try:
        design_chain(pp)
        check(f"geometry.{key}: {val} を拒否", False)
    except ValueError:
        check(f"geometry.{key}: {val} を拒否", True)
pp = load_problem(PROB)
pp.geometry["initial_line"] = "cfd"
try:
    design_chain(pp)
    check("initial_line: cfd で initial_line_run 無しを拒否", False)
except ValueError:
    check("initial_line: cfd で initial_line_run 無しを拒否", True)
pp = load_problem(PROB)
pp.geometry.update({"initial_line": "cfd", "initial_line_res": "res_6000.h5", "Md_moc_offset": -4.16e-4,
                    "initial_line_run": os.path.relpath(run, PROB.parent)})      # problem ファイルからの相対パス
if (run / "res_6000.h5").exists():
    d = design_chain(pp)
    il = d["initial_line"]
    check(f"相対パスの initial_line_run を解決 ({il['run']})", Path(il["run"]).resolve() == run.resolve())
    check(f"報告の Md は spec のまま ({d['Md']}), MOC は較正値 ({d['Md_moc']})",
          d["Md"] == float(pp.spec["M_design"]) and d["Md_moc"] == float(pp.spec["M_design"]) - 4.16e-4)
    check(f"x0 = CFD 線の軸着地 ({il['x0']:.6f})", abs(il["x0"] - d["x_A"]) == 0.0 and il["source"] == "cfd")
else:
    print("skip: res_6000.h5 が無いので design_chain (cfd) の検査を省略")

print(json.dumps(out, indent=1))
print("FAIL 件数:", FAIL)
sys.exit(1 if FAIL else 0)
