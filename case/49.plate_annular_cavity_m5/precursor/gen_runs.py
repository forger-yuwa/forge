#!/usr/bin/env python3
"""case/49 前駆 2D 平板 (断熱, M 可変) の forge run を生成して段階起動する。

usage:
  python3 gen_runs.py --run run_0001_precursor_m5 --mesh fp_y1_6um [--main-steps 20000] [--cfl 2]

段階 (case/48 の実測レシピ): 層流暖機 (1 次, cfl 0.2) → SST soft (1 次, cfl 0.3) → mid (1 次, cfl 1.0)
→ 2 次ランプ (cfl 0.5/1/2) → 本段 (2 次, cfl, implicitRelax)。段間は同一メッシュなので
`interp_field.py` (2D 平面なので (x,y) 最近傍 = index 同値) で mesh.h5 の VALUE を更新する。

自由流は ../conditions.json (setup.py) が単一ソース。壁は**断熱** (3D 側の平板と同じ)。
"""
import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
CASE = HERE.parent
ROOT = CASE.parents[1]
TOOLS = ROOT / "solver_density_cuda" / "tools"
ENV = dict(os.environ,
           LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu/hdf5/serial:" + os.environ.get("LD_LIBRARY_PATH", ""))

sys.path.insert(0, str(CASE))
from setup import load as load_conditions  # noqa: E402

D = load_conditions()
IC_DELTA0 = 3.0e-4          # IC の壁近傍速度ランプ幅 [m] (前縁の巨大せん断を避ける)

CFG = """mesh: {{discretization: "node", nodeWallDirichlet: 1, nodeInletCornerWall: 1, meshFileName: "mesh.h5", valueFileName: "mesh.h5"}}
gpu: 1
solver: "SLAU"
physProp: {{thermalMethod: 0, viscMethod: 1, visc: {mu:.6g}, thermCond: {kc:.6g}, thermCondMethod: 1, prandtlLam: {pr}, cp: {cp}, gamma: {gam}}}
time:
  unsteady: 0
  dualTime: 0
  last: {{nStepOuter: {nsteps}}}
  deltaT: {{control: 1, dt: 1e-8, cfl: {cfl}, cfl_pseudo: {cfl}, implicitRelax: {relax}, blockDPLUR: 1, lowMachPrecond: 0, dt_min: 1e-10, dt_max: 1.0, detectNaN: 1}}
  outStepStart: 0
  outStepInterval: {outint}
  timeIntegration: 11
  nStepInner: {ninner}
space: {{convMethod: {conv}, limiter: {lim}, pRef: {pref:.6g}}}
turbulence: {{model: "{model}", scalarDiffusion: 1, dilatationCorrection: 2, katoLaunder: 0, wallTreatmentSST: 0, turbulentPrandtl: {prt}, kInit: {k:.6g}, omegaInit: {om:.6g}}}
output: {{level: 1}}
initial: "uniform_p101325_u10"
"""
BC = """inlet:  {{physID: 1, kind: inlet_uniformVelocity, outputHDFflg: 0, ints: , floats: {{ro: {ro:.6g}, Ux: {u:.6g}, Uy: 0.0, Uz: 0.0, Ps: {p:.6g}, k: {k:.6g}, omega: {om:.6g}}}}}
outlet: {{physID: 2, kind: outlet_statPress, outputHDFflg: 0, ints: , floats: {{Ps: {p:.6g}, Pt: {p:.6g}, Tt: {t:.6g}}}}}
top:    {{physID: 3, kind: slip, outputHDFflg: 0, ints: , floats: }}
wall:   {{physID: 4, kind: wall, outputHDFflg: 1, ints: , floats: {{Ux: 0.0, Uy: 0.0, Uz: 0.0}}}}
sym:    {{physID: 5, kind: slip, outputHDFflg: 0, ints: , floats: }}
"""


def cfg(nsteps, cfl, relax, conv, lim, ninner, outint, lam=False):
    return CFG.format(mu=D["mu_inf"], kc=D["mu_inf"] * D["cp"] / D["prandtl_lam"],
                      pr=D["prandtl_lam"], cp=D["cp"], gam=D["gamma"], nsteps=nsteps, cfl=cfl,
                      relax=relax, outint=outint, ninner=ninner, conv=conv, lim=lim,
                      pref=D["P_inf"], model="none" if lam else "sst", prt=D["prandtl_turb"],
                      k=D["k_inf"], om=D["omega_inf"])


def patch_ic(h5):
    """一様自由流 + 平板上 (x>=0) は壁距離 tanh で速度を落とし、壁ノードは u=0。"""
    with h5py.File(h5, "r+") as f:
        n = f["/VALUE/ro"].shape[0]
        c = f["/MESH/COORD"][:].reshape(-1, 3)
        wd = f["/VALUE/wall_dist"][:]
        ro = np.full(n, D["ro_inf"]); u = np.full(n, D["U_inf"])
        onplate = c[:, 0] >= -1e-9
        u[onplate] = D["U_inf"] * np.tanh(np.maximum(wd[onplate], 0.0) / IC_DELTA0)
        u[(wd <= 0.0) & onplate] = 0.0
        roe = D["P_inf"] / (D["gamma"] - 1.0) + 0.5 * ro * u ** 2
        f["/VALUE/ro"][:] = ro.astype(np.float32)
        f["/VALUE/roUx"][:] = (ro * u).astype(np.float32)
        f["/VALUE/roUy"][:] = np.zeros(n, np.float32)
        f["/VALUE/roUz"][:] = np.zeros(n, np.float32)
        f["/VALUE/roe"][:] = roe.astype(np.float32)
        for name, v in (("roK", ro * D["k_inf"]), ("roOmega", ro * D["omega_inf"])):
            if "/VALUE/" + name in f:
                f["/VALUE/" + name][:] = v.astype(np.float32)
            else:
                f.create_dataset("/VALUE/" + name, data=v.astype(np.float32))
        print("IC: n=%d wall nodes=%d" % (n, int(((wd <= 0.0) & onplate).sum())))


def run_forge(rd):
    r = subprocess.run([str(TOOLS / "run_case.sh"), str(rd)], env=ENV, capture_output=True, text=True)
    (rd / "run_case_stdout.log").write_text(r.stdout + r.stderr)
    return r.returncode


def stage(rd, text, nsteps, tag):
    (rd / "solverConfig.yaml").write_text(text)
    rc = run_forge(rd)
    res = sorted(rd.glob("res_[0-9]*.h5"), key=lambda f: int(f.stem.split("_")[1]))
    if rc != 0 or not res or int(res[-1].stem.split("_")[1]) < nsteps:
        raise SystemExit("stage %s failed rc=%d res=%s" % (tag, rc, [r.name for r in res][-2:]))
    subprocess.run([sys.executable, str(TOOLS / "interp_field.py"), str(res[-1]), str(rd / "mesh.h5")],
                   env=ENV, check=True, capture_output=True, text=True)
    with h5py.File(rd / "mesh.h5", "r+") as f:      # 層流段の res には roK/roOmega が無い
        ro = f["/VALUE/ro"][:].astype(float)
        for name, val in (("roK", D["k_inf"]), ("roOmega", D["omega_inf"])):
            ds = "/VALUE/" + name
            if ds not in f or not np.all(np.isfinite(f[ds][:])) or float(np.abs(f[ds][:]).max()) == 0.0:
                if ds in f:
                    del f[ds]
                f.create_dataset(ds, data=(ro * val).astype(np.float32))
    for f in rd.glob("res_*"):
        f.unlink()
    for nm in ("residual_history.csv", "residual_history.png", "CONVERGENCE_VERDICT.txt", "forge_run.log"):
        if (rd / nm).exists():
            shutil.move(str(rd / nm), str(rd / ("%s_%s%s" % (Path(nm).stem, tag, Path(nm).suffix))))
    print("  stage %s done" % tag)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    ap.add_argument("--mesh", required=True)
    ap.add_argument("--main-steps", type=int, default=20000)
    ap.add_argument("--cfl", type=float, default=2.0)
    ap.add_argument("--relax", type=float, default=0.7)
    ap.add_argument("--out-int", type=int, default=2000)
    ap.add_argument("--ramp", default="0.5,1,2")
    ap.add_argument("--ramp-steps", type=int, default=2000)
    ap.add_argument("--lam-cfl", type=float, default=0.2, help="層流暖機の CFL")
    ap.add_argument("--lam-steps", type=int, default=2000)
    ap.add_argument("--sst-cfl", default="0.3,1.0",
                    help="SST 1 次段の CFL 列 (カンマ区切り)。M9 など厳しい条件では刻む")
    ap.add_argument("--sst-steps", type=int, default=2000)
    ap.add_argument("--dry", action="store_true")
    a = ap.parse_args()

    rd = HERE / a.run
    if rd.exists():
        raise SystemExit("%s exists" % rd)
    rd.mkdir()
    shutil.copy(HERE / "mesh" / (a.mesh + ".h5"), rd / "mesh.h5")
    (rd / "bcondConfig.yaml").write_text(BC.format(ro=D["ro_inf"], u=D["U_inf"], p=D["P_inf"],
                                                   t=D["T_inf"], k=D["k_inf"], om=D["omega_inf"]))
    (rd / "GEN_ARGS").write_text(" ".join(sys.argv[1:]) + "\n")
    (rd / "probe.yaml").write_text("outStepInterval: 100\noutStepStart: 0\npoints:\nsurfaces:\n")
    main_cfg = cfg(a.main_steps, a.cfl, a.relax, 1, 2, 4, a.out_int)
    (rd / "solverConfig.yaml").write_text(main_cfg)
    patch_ic(rd / "mesh.h5")
    if a.dry:
        return
    # **SST 段の CFL は指定できるようにする**。M5 の 0.3 → 1.0 は M9 では `mid` で落ちた
    # (2026-09-19)。マッハ数を上げるときは `--sst-cfl` を細かく刻む。
    stage(rd, cfg(a.lam_steps, a.lam_cfl, a.relax, 0, 0, 10, a.lam_steps, lam=True),
          a.lam_steps, "lam")
    for i, cv in enumerate([float(x) for x in a.sst_cfl.split(",") if x]):
        stage(rd, cfg(a.sst_steps, cv, a.relax, 0, 0, 10, a.sst_steps), a.sst_steps,
              "sst%d_cfl%g" % (i, cv))
    for i, cv in enumerate([float(v) for v in a.ramp.split(",") if v]):
        stage(rd, cfg(a.ramp_steps, cv, a.relax, 1, 2, 4, a.ramp_steps), a.ramp_steps,
              "ramp%d_cfl%g" % (i, cv))
    (rd / "solverConfig.yaml").write_text(main_cfg)
    rc = run_forge(rd)
    print("main rc", rc)
    print((rd / "CONVERGENCE_VERDICT.txt").read_text()[-800:])


if __name__ == "__main__":
    main()
