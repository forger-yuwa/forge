#!/usr/bin/env python3
"""U1: エネルギー行の熱伝導 Jacobian (`implicitThermalJacobian`) の純伝導の試験
(plan time_integration-implicit-thermal-jacobian §6 U1)。U2・U3 (plan time_integration-line-viscous-jacobian §6) も同じ板で回す:
U2 = ライン陰解法 + 方向別の擬似 dt + `lineViscCoupling` 0/2 の純伝導、U3 = 上壁を z 方向に U で動かす Couette (両壁 300 K)。
U3 の解析解 (定数 μ・k、圧力一様): u_z = U y/H、T = T_w + (μU²/2k)(y/H)(1 − y/H)。

case/52 の流体の板 (一様 20×16 quad、H 0.01 m、W 0.02 m、空気 CPG、定数 k 0.0241、静止、低圧 1013.25 Pa) の
下壁 300 K・上壁 350 K の等温壁で、初期の一様 300 K から回す。解析解は直線 T(y) = 300 + 50·y/H (q = kΔT/H = 120.5 W/m²)。
キー (0 / 1) と cfl_pseudo を変え、温度場の解析解からの最大のずれ (L∞) が 0.05 K (ΔT の 0.1 %) を切る step 数を比べる。

usage:
  python3 u1_thermjac.py mesh                                    (手元: gmsh で mesh/u1_slab.msh を作る)
  python3 u1_thermjac.py prep <run> --key K --cfl C [--steps 20000] [--out 500] --conv <convertGmshToForge>
                               [--line 1] [--dir 1] [--lvc 0|1|2] [--couette U]   (U2・U3)
                                                                 (AWS: 変換・静止の初期場・設定; forge は run_case.sh で回す)
  python3 u1_thermjac.py eval <run> [<run> ...]                  (L∞ の時系列と 0.05 K に入る step → <run>/U1_EVAL.json)
"""
import argparse
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
H, W, NX, NY = 0.01, 0.02, 20, 16
T_BOT, T_TOP, P0 = 300.0, 350.0, 1013.25
CP, GAM, K = 1004.5, 1.4, 0.0241
MU = 1.81e-5
TOL_K = 0.05

GEO = """// U1: 一様構造 quad (physID 1 左 / 2 右 / 3 下 300 K / 4 上 350 K)
W = {W}; H = {H}; nx = {nx}; ny = {ny};
Point(1) = {{0, 0, 0, 1}};  Point(2) = {{W, 0, 0, 1}};
Point(3) = {{W, H, 0, 1}};  Point(4) = {{0, H, 0, 1}};
Line(1) = {{1, 2}}; Line(2) = {{2, 3}}; Line(3) = {{3, 4}}; Line(4) = {{4, 1}};
Line Loop(1) = {{1, 2, 3, 4}}; Plane Surface(1) = {{1}};
Transfinite Line{{1, 3}} = nx + 1;
Transfinite Line{{2, 4}} = ny + 1;
Transfinite Surface{{1}}; Recombine Surface{{1}};
Physical Curve("side_left",  1) = {{4}};
Physical Curve("side_right", 2) = {{2}};
Physical Curve("wall_bot",   3) = {{1}};
Physical Curve("wall_top",   4) = {{3}};
Physical Surface("fluid",    5) = {{1}};
Mesh.MshFileVersion = 4.1;
"""

BCOND = """side_left:  {physID: 1, kind: slip, outputHDFflg: 0, ints: , floats: }
side_right: {physID: 2, kind: slip, outputHDFflg: 0, ints: , floats: }
wall_bot:   {physID: 3, kind: wall_isothermal, outputHDFflg: 1, ints: , floats: {Ux: 0.0, Uy: 0.0, Uz: 0.0, Ts: 300.0}}
wall_top:   {physID: 4, kind: wall_isothermal, outputHDFflg: 1, ints: , floats: {Ux: 0.0, Uy: 0.0, Uz: {uz_top}, Ts: {t_top}}}
"""

SOLVER = """# U1 (plan time_integration-implicit-thermal-jacobian §6): 純伝導、下 300 K / 上 350 K。case/52 template と同じ数値設定。
mesh: {{discretization: "node", nodeWallDirichlet: 1, meshFileName: "mesh.h5", valueFileName: "mesh.h5"}}
gpu: 1
detectNaN: 1
solver: "SLAU"
physProp:
  isAxisymmetric: 0
  thermalMethod: 0
  viscMethod: 0
  visc: 1.81e-5
  thermCond: 0.0241
  cp: 1004.5
  gamma: 1.4
time:
  unsteady: 0
  dualTime: 0
  last: {{nStepOuter: {steps}, time: 1.00}}
  deltaT: {{control: 1, dt: 0.00000001, cfl: 1.0, cfl_pseudo: {cfl}, blockDPLUR: 1, implicitRelax: 0.5, dt_min: 0.00000000001, dt_max: 1.0{key}}}
  outStepStart: 0
  outStepInterval: {out}
  timeIntegration: 11
  nStepInner: 5
space: {{convMethod: 2, limiter: 2}}
turbulence: {{model: "none"}}
output: {{level: 1}}
initial: "uniform_p101325_u10"
"""


def cmd_mesh():
    m = HERE / "mesh"; m.mkdir(exist_ok=True)
    (m / "u1_slab.geo").write_text(GEO.format(W=W, H=H, nx=NX, ny=NY))
    subprocess.run(["gmsh", "-2", "u1_slab.geo", "-o", "u1_slab.msh", "-format", "msh41", "-v", "1"], cwd=m, check=True,
                   stdout=subprocess.DEVNULL)
    print(m / "u1_slab.msh")


def cmd_prep(run: Path, key: int, cfl: float, steps: int, out: int, conv: str, line: int = 0, ddir: int = 0, lvc: int = 0,
             couette: float = 0.0):
    import h5py
    if run.exists():
        sys.exit(f"{run} が既にある — 止める")
    run.mkdir(parents=True)
    shutil.copy(HERE / "mesh" / "u1_slab.msh", run / "u1_slab.msh")
    keytxt = f", implicitThermalJacobian: {key}" if key else ""
    if line:
        keytxt += ", lineImplicit: 1" + (", lineDtDirectional: 1" if ddir else "") + (f", lineViscCoupling: {lvc}" if lvc else "")
    (run / "solverConfig.yaml").write_text(SOLVER.format(steps=steps, cfl=cfl, out=out, key=keytxt))
    (run / "bcondConfig.yaml").write_text(BCOND.replace("{uz_top}", repr(float(couette)))
                                          .replace("{t_top}", repr(T_BOT if couette else T_TOP)))
    (run / "probe.yaml").write_text("outStepInterval: 100\noutStepStart: 0\npoints:\nsurfaces:\n")   # 無いと forge が止まる
    r = subprocess.run([conv, "u1_slab.msh", "mesh.h5"], cwd=run, capture_output=True, text=True)
    (run / "convert.log").write_text(r.stdout + r.stderr)
    if not (run / "mesh.h5").exists():
        sys.exit(f"変換に失敗 (rc {r.returncode})")
    R = CP * (GAM - 1.0) / GAM
    ro = P0 / (R * T_BOT)
    with h5py.File(run / "mesh.h5", "a") as f:
        n = f["CELLS/centCoords"].shape[0] // 3
        dt = f["CELLS/centCoords"].dtype                       # 変換器の精度に合わせる (FP64 は float64)
        g = f.require_group("VALUE")
        for name, val in (("ro", ro), ("roUx", 0.0), ("roUy", 0.0), ("roUz", 0.0), ("roe", ro * CP / GAM * T_BOT)):
            if name in g:
                del g[name]
            g.create_dataset(name, data=np.full(n, val, dtype=dt))
    json.dump({"key": key, "cfl": cfl, "steps": steps, "out": out, "converter": conv, "ic_dtype": str(dt),
               "line": line, "dir": ddir, "lvc": lvc, "couette": couette},
              open(run / "U1_PREP.json", "w"), indent=1)
    print(f"[u1 prep] {run.name}: key {key}, cfl {cfl}, line {line} dir {ddir} lvc {lvc} couette {couette}, {steps} step, 出力 {out} ごと、IC {dt}")


def cmd_eval(runs):
    import h5py
    for run in runs:
        with h5py.File(run / "mesh.h5", "r") as f:
            y = np.asarray(f["MESH/COORD"][:], float).reshape(-1, 3)[:, 1]
        prep = json.load(open(run / "U1_PREP.json"))
        U = float(prep.get("couette", 0.0))
        if U:
            texact = T_BOT + (MU * U * U / (2 * K)) * (y / H) * (1 - y / H)
            uexact = U * y / H
        else:
            texact = T_BOT + (T_TOP - T_BOT) * y / H
        rows = []
        for p in sorted(run.glob("res_[0-9]*.h5"), key=lambda q: int(re.findall(r"res_(\d+)", q.name)[0])):
            st = int(re.findall(r"res_(\d+)", p.name)[0])
            with h5py.File(p, "r") as f:
                T = np.asarray(f["VALUE/T"][:], float)
                uz = np.asarray(f["VALUE/Uz"][:], float) if U else None
            eu = float(np.max(np.abs(uz - uexact))) / U if U else 0.0
            rows.append((st, float(np.max(np.abs(T - texact))), bool(np.all(np.isfinite(T))), eu))
        if U:   # U3: u の L∞ < 0.1 % U かつ T の L∞ < 1 % (T − T_w の最大)
            tol_t = 0.01 * MU * U * U / (8 * K)
            hit = next((st for st, e, fin, eu in rows if fin and e < tol_t and eu < 1e-3), None)
        else:
            hit = next((st for st, e, fin, eu in rows if fin and e < TOL_K), None)
        nan = any(not r[2] for r in rows) or any(run.glob("res_nan_*.h5"))
        out = {"run": run.name, **prep, "linf_series": rows, "first_step_below_tol": hit, "tol_K": TOL_K, "nonfinite": nan}
        (run / "U1_EVAL.json").write_text(json.dumps(out, indent=1))
        print(f"{run.name}: key {prep['key']} cfl {prep['cfl']} → L∞ < {TOL_K} K に入る step: {hit}  非有限: {nan}  "
              f"L∞ の推移 " + " ".join(f"{r[0]}:{r[1]:.3g}" + (f"/u{r[3]:.2e}" if U else "") for r in rows[::max(1, len(rows)//8)]))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest="cmd", required=True)
    sp.add_parser("mesh")
    p = sp.add_parser("prep"); p.add_argument("run"); p.add_argument("--key", type=int, required=True); p.add_argument("--cfl", type=float, required=True)
    p.add_argument("--steps", type=int, default=20000); p.add_argument("--out", type=int, default=500); p.add_argument("--conv", required=True)
    p.add_argument("--line", type=int, default=0); p.add_argument("--dir", type=int, default=0); p.add_argument("--lvc", type=int, default=0)
    p.add_argument("--couette", type=float, default=0.0)
    p = sp.add_parser("eval"); p.add_argument("runs", nargs="+")
    a = ap.parse_args()
    if a.cmd == "mesh":
        cmd_mesh()
    elif a.cmd == "prep":
        cmd_prep(HERE / a.run, a.key, a.cfl, a.steps, a.out, a.conv, a.line, a.dir, a.lvc, a.couette)
    else:
        cmd_eval([HERE / r for r in a.runs])
