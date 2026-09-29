#!/usr/bin/env python3
"""V2c 超音速の斜め衝撃波 (plan boundary-node-farfield-characteristic §6 V2c) の run を作る (AWS で実行)。
3D 薄板 (z 3 層、5 mm、slip)。下面 = 平板 (x ≤ 0.2 m) + 10° ランプ (0.2–0.7 m) + 水平 (0.7 m–、slip、x ≥ 0.2 m が評価線)。
ランプを出口まで伸ばすと流路が 0.3 → 0.124 m に絞られ、slip 上面 (A) で閉塞して発散した (2026-09-29) ので 0.7 m で止める、M 2.5 (CPG 空気 300 K、101325 Pa)。
ランプ角の衝撃波 (β ≈ 31.9°) が上面に当たる。上面:
  --case A: H 0.3 m、slip (反射する対照)      --case B: H 0.3 m、farfield
  --case C: H 0.8 m、slip (上面での反射点が x > 1.2 m になる高さ。共通領域 0 ≤ j ≤ 60 の格子は A/B と同一)
格子: x 0–1.2 m を 240 分割 (Δx 5 mm)。y は 0 ≤ j ≤ 60 で下面から H_A = 0.3 m (水平) までを等分、C はその上に 5 mm で 100 行足す。
入口 xmin = farfield (超音速流入 → 自由流)、出口 xmax = outflow。定常 block-DPLUR。
  python3 setup_v2c.py RUN --case A|B|C [--steps 6000] [--cfl 1.0]
"""
import argparse, math, os, subprocess, sys
import h5py, numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from make_box_msh import write_hex_msh  # noqa: E402

XC, XE, TH, LX, HA, NX, JA, DZ, NZ = 0.2, 0.7, math.radians(10.0), 1.2, 0.3, 240, 60, 0.005, 3


def yb(x):
    return (min(max(x, XC), XE) - XC) * math.tan(TH)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run"); ap.add_argument("--case", choices=("A", "B", "C"), required=True)
    ap.add_argument("--steps", type=int, default=6000); ap.add_argument("--cfl", type=float, default=1.0)
    ap.add_argument("--bin", default=os.path.expanduser("~/forge-pgrad-new/solver_density_cuda/build-ff"))
    a = ap.parse_args()
    run = a.run
    os.makedirs(run)
    JC = 100 if a.case == "C" else 0
    ny = JA + JC
    xs = np.linspace(0.0, LX, NX + 1)

    def coord(i, j, k):
        x = xs[i]
        y = yb(x) + (j / JA) * (HA - yb(x)) if j <= JA else HA + (j - JA) * DZ
        return x, y, k * DZ
    write_hex_msh(os.path.join(run, "ramp.msh"), NX, ny, NZ, coord, (LX, HA + JC * DZ, NZ * DZ))
    P, T, gam, cp, M = 101325.0, 300.0, 1.4, 1004.5, 2.5
    R = cp - cp / gam; ro = P / (R * T); c = math.sqrt(gam * R * T); U = M * c
    open(os.path.join(run, "solverConfig.yaml"), "w").write(f"""mesh: {{discretization: "node", nodeWallDirichlet: 1, meshFileName: "ramp.h5", valueFileName: "ramp.h5"}}
gpu: 1
solver: "SLAU"
physProp: {{thermalMethod: 0, viscMethod: 0, visc: 0.0, thermCond: 0.0, cp: {cp}, gamma: {gam}}}
time:
  unsteady: 0
  dualTime: 0
  last: {{nStepOuter: {a.steps}}}
  deltaT: {{control: 1, dt: 1e-8, cfl: {a.cfl}, cfl_pseudo: {a.cfl}, dt_min: 1e-12, dt_max: 1.0, blockDPLUR: 1, lowMachPrecond: 0, detectNaN: 1}}
  outStepStart: 0
  outStepInterval: {max(a.steps // 6, 1)}
  timeIntegration: 11
  nStepInner: 5
space: {{convMethod: 1, limiter: 2, pRef: {P}, limiterRefLength: 1.0, limiterRoRef: {ro!r}, limiterPRef: {P!r}, limiterARef: {c!r}}}
turbulence: {{model: "none"}}
initial: "uniform_p101325_u10"
""")
    fl = f"ro: {ro!r}, Ux: {U!r}, Uy: 0.0, Uz: 0.0, Ps: {P!r}"
    top = "farfield" if a.case == "B" else "slip"
    kinds = {1: ("xmin", "farfield"), 2: ("xmax", "outflow"), 3: ("ymin", "slip"), 4: ("ymax", top), 5: ("zmin", "slip"), 6: ("zmax", "slip")}
    open(os.path.join(run, "bcondConfig.yaml"), "w").write("".join(
        f"{nm}: {{physID: {p}, kind: {k}, outputHDFflg: 0, ints: , floats: {{{fl if k == 'farfield' else ''}}}}}\n" for p, (nm, k) in kinds.items()))
    open(os.path.join(run, "probe.yaml"), "w").write("outStepInterval: 100000\noutStepStart: 0\npoints:\nsurfaces:\n")
    subprocess.run([os.path.join(a.bin, "convertGmshToForge"), "ramp.msh", "ramp.h5"], cwd=run,
                   stdout=open(os.path.join(run, "convert.log"), "w"), stderr=subprocess.STDOUT)
    if not os.path.exists(os.path.join(run, "ramp.h5")) or "writeInputH5: wrote /VIZMESH" not in open(os.path.join(run, "convert.log")).read():
        raise SystemExit(f"変換に失敗: {run}/convert.log")
    e = P / ((gam - 1.0) * ro)
    with h5py.File(os.path.join(run, "ramp.h5"), "r+") as f:
        V = f["VALUE"]; n = V["ro"].shape[0]
        vals = {"ro": ro, "roUx": ro * U, "roUy": 0.0, "roUz": 0.0, "roe": ro * (e + 0.5 * U * U)}
        for kk, v in vals.items():
            V[kk][...] = np.full(n, v, dtype=V[kk].dtype)
    open(os.path.join(run, "IC_FROM.txt"), "w").write(f"V2c case {a.case}: 一様 M {M} (ρ {ro:.6g}、U {U:.6g})、上面 {top}、H {HA + JC * DZ} m\n")
    print(f"prepared {run}: case {a.case}, nodes ≈ {(NX + 1) * (ny + 1) * (NZ + 1)}")


if __name__ == "__main__":
    main()
