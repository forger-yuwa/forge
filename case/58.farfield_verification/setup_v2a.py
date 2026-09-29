#!/usr/bin/env python3
"""V2a 音響反射 (plan boundary-node-farfield-characteristic §6 V2a) の run を作る (AWS で実行)。
3D 薄板チャネル: x ∈ [0, L] (Δx 5 mm)、y・z 各 3 層 (15 mm、側壁 slip)。一様流 (CPG 空気、300 K、101325 Pa、M) に
右向きの単純音波パルス (p' = δ g(x)、u' = p'/(ρc)、ρ' = p'/c²、δ = 1e-3 P∞、g = ガウス 半値全幅 40 mm、中心 x0 0.4 m) を置く。
左端 (x = 0) は farfield、右端は --end (farfield / slip)。プローブ: x = 0.8 m (評価点)、0.4 m。
dual-time (timeIntegration 11 + dualTime 1、物理 dt 固定)。
リミッタの基準値 (L_ref・ρ・P・a) は自由流で固定する: 既定の自動決定は領域の対角長と初期場の平均を使うので、短・長の領域で
離散化が変わり、その差 (入射パルス通過時に約 1 %) が反射と区別できなくなる (2026-09-29 run_0020/0021 で確認)。
  python3 setup_v2a.py RUN --len 1|3 --mach 0.3 --end farfield|slip [--solver SLAU|SLAU2] [--dt 5e-6] [--nsub 20]
         [--tend 3.0e-3] [--cflp 12] [--bin DIR]
"""
import argparse, math, os, subprocess, sys
import h5py, numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run"); ap.add_argument("--len", type=float, required=True); ap.add_argument("--mach", type=float, required=True)
    ap.add_argument("--end", choices=("farfield", "slip"), required=True); ap.add_argument("--solver", default="SLAU")
    ap.add_argument("--dt", type=float, default=5e-6); ap.add_argument("--nsub", type=int, default=20)
    ap.add_argument("--tend", type=float, default=3.0e-3); ap.add_argument("--cflp", type=float, default=12.0)
    ap.add_argument("--bin", default=os.path.expanduser("~/forge-pgrad-new/solver_density_cuda/build-ff"))
    a = ap.parse_args()
    run = a.run
    os.makedirs(run)
    dx, W, nyz = 0.005, 0.015, 3
    nx = int(round(a.len / dx))
    subprocess.run([sys.executable, os.path.join(HERE, "make_box_msh.py"), os.path.join(run, "chan.msh"), str(nx), str(nyz), str(nyz),
                    repr(a.len), repr(W), repr(W)], check=True)
    P, T, gam, cp = 101325.0, 300.0, 1.4, 1004.5
    R = cp - cp / gam; ro = P / (R * T); c = math.sqrt(gam * R * T); U = a.mach * c
    steps = int(math.ceil(a.tend / a.dt))
    open(os.path.join(run, "solverConfig.yaml"), "w").write(f"""mesh: {{discretization: "node", nodeWallDirichlet: 1, meshFileName: "chan.h5", valueFileName: "chan.h5"}}
gpu: 1
solver: "{a.solver}"
physProp: {{thermalMethod: 0, viscMethod: 0, visc: 0.0, thermCond: 0.0, cp: {cp}, gamma: {gam}}}
time:
  unsteady: 1
  dualTime: 1
  nSubIterDualTime: {a.nsub}
  last: {{nStepOuter: {steps}}}
  deltaT: {{control: 0, dt: {a.dt!r}, cfl: 1.0, cfl_pseudo: {a.cflp!r}, dt_min: 1e-12, dt_max: 1.0, blockDPLUR: 1, lowMachPrecond: 0, detectNaN: 1}}
  outStepStart: 0
  outStepInterval: {steps}
  timeIntegration: 11
  nStepInner: 5
space: {{convMethod: 1, limiter: 2, pRef: {P}, limiterRefLength: 1.0, limiterRoRef: {ro!r}, limiterPRef: {P!r}, limiterARef: {c!r}}}
turbulence: {{model: "none"}}
initial: "uniform_p101325_u10"
""")
    fl = f"ro: {ro!r}, Ux: {U!r}, Uy: 0.0, Uz: 0.0, Ps: {P!r}"
    kinds = {1: ("xmin", "farfield"), 2: ("xmax", a.end), 3: ("ymin", "slip"), 4: ("ymax", "slip"), 5: ("zmin", "slip"), 6: ("zmax", "slip")}
    lines = []
    for p, (nm, k) in kinds.items():
        lines.append(f"{nm}: {{physID: {p}, kind: {k}, outputHDFflg: 0, ints: , floats: {{{fl if k == 'farfield' else ''}}}}}\n")
    open(os.path.join(run, "bcondConfig.yaml"), "w").write("".join(lines))
    yp = W / 3.0
    open(os.path.join(run, "probe.yaml"), "w").write(
        f"outStepInterval: 1\noutStepStart: 0\npoints:\n  eval: {{x: 0.8, y: {yp}, z: {yp}}}\n  mid: {{x: 0.4, y: {yp}, z: {yp}}}\nsurfaces:\n")
    subprocess.run([os.path.join(a.bin, "convertGmshToForge"), "chan.msh", "chan.h5"], cwd=run,
                   stdout=open(os.path.join(run, "convert.log"), "w"), stderr=subprocess.STDOUT)
    if not os.path.exists(os.path.join(run, "chan.h5")) or "writeInputH5: wrote /VIZMESH" not in open(os.path.join(run, "convert.log")).read():
        raise SystemExit(f"変換に失敗: {run}/convert.log")
    delta, fwhm, x0 = 1e-3 * P, 0.04, 0.4
    sig = fwhm / (2.0 * math.sqrt(2.0 * math.log(2.0)))
    with h5py.File(os.path.join(run, "chan.h5"), "r+") as f:
        V = f["VALUE"]; n = V["ro"].shape[0]
        # 節点座標 (node では値の数 = 節点数、MESH/COORD)
        xyz = np.array(f["MESH/COORD"]).reshape(-1, 3)
        if len(xyz) != n:
            raise SystemExit(f"MESH/COORD の数 {len(xyz)} が値の数 {n} と違う (node 変換でない?)")
        x = xyz[:, 0].astype(float)
        pp = delta * np.exp(-0.5 * ((x - x0) / sig) ** 2)
        rr = ro + pp / c ** 2; uu = U + pp / (ro * c); PP = P + pp
        e = PP / ((gam - 1.0) * rr)
        vals = {"ro": rr, "roUx": rr * uu, "roUy": 0.0 * rr, "roUz": 0.0 * rr, "roe": rr * (e + 0.5 * uu * uu)}
        for kk, v in vals.items():
            V[kk][...] = v.astype(V[kk].dtype)
    open(os.path.join(run, "IC_FROM.txt"), "w").write(
        f"V2a 音響: L {a.len} m、M {a.mach}、右端 {a.end}、{a.solver}、dt {a.dt}、nSub {a.nsub}、cfl_pseudo {a.cflp}、"
        f"パルス δ {delta} Pa・FWHM {fwhm} m・中心 {x0} m (右向き単純波)、背景 ρ {ro:.6g} c {c:.6g} U {U:.6g}\n")
    print(f"prepared {run}: nx {nx}, steps {steps}, acoustic CFL {(c + U) * a.dt / dx:.3f}")


if __name__ == "__main__":
    main()
