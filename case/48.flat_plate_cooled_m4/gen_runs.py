#!/usr/bin/env python3
"""case/48 の forge run を生成して段階起動する。

usage: python3 gen_runs.py --run run_0001_A_ad_y3 --mesh fp_y1_3um --wall adiabatic|300|700 [--plain] [--main-steps 12000] [--cfl 4] [--relax 0.7] [--dry]
  --wall adiabatic → kind: wall / 数値 → wall_isothermal Ts=数値
  --plain          → dilatationCorrection 0, katoLaunder 0 (SU2 比較用の素 SST)
段階: soft (1次, cfl 0.5, nStepInner 10, 2000) → mid (1次, cfl 1.0, 2000) → main (2次, cfl, implicitRelax, main-steps)。
段間は interp_field.py (同一メッシュなので index 同値) で mesh h5 の VALUE を更新。
"""
import argparse, os, re, shutil, subprocess, sys
from pathlib import Path
import h5py, numpy as np

ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT / "solver_density_cuda" / "tools"
HERE = Path(__file__).resolve().parent
ENV = dict(os.environ, LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu/hdf5/serial:" + os.environ.get("LD_LIBRARY_PATH", ""))

# 自由流 (case/44 va3 試験部相当, 空気 CPG)
GAM, CP = 1.4, 1004.5; R = CP * (GAM - 1) / GAM
M_INF, P_INF, T_INF = 4.19, 5037.4, 283.0
RO_INF = P_INF / (R * T_INF); A_INF = (GAM * R * T_INF) ** 0.5; U_INF = M_INF * A_INF
K_INF, OM_INF = 75.0, 26000.0     # TI 0.5 %, mu_t/mu ~ 10
IC_DELTA0 = 3.0e-4               # IC の壁近傍速度ランプ幅 [m]

CFG = """mesh: {{meshFormat: "hdf5", discretization: "node", nodeWallDirichlet: 1, nodeInletCornerWall: 1, meshFileName: "mesh.h5", valueFileName: "mesh.h5"}}
gpu: 1
solver: "SLAU"
physProp: {{isCompressible: 1, thermalMethod: 0, viscMethod: 1, ro: 1.2, visc: 1.8e-5, thermCond: 0.0257, thermCondMethod: 1, prandtlLam: 0.72, cp: {cp}, gamma: {gam}}}
time:
  unsteady: 0
  dualTime: 0
  last: {{control: 0, nStepOuter: {nsteps}}}
  deltaT: {{control: 1, dt: 1e-8, cfl: {cfl}, cfl_pseudo: {cfl}, implicitRelax: {relax}, blockDPLUR: 1, lowMachPrecond: 0, dt_min: 1e-10, dt_max: 1.0, detectNaN: 1}}
  outStepStart: 0
  outStepInterval: {outint}
  timeIntegration: 11
  nStepInner: {ninner}
space: {{convMethod: {conv}, limiter: {lim}}}
turbulence: {{model: "sst", scalarDiffusion: 1, dilatationCorrection: {dil}, katoLaunder: {kl}, wallTreatmentSST: 0, turbulentPrandtl: 0.9, kInf: {k}, omegaInf: {om}}}
output: {{level: 1}}
initial: "uniform_p101325_u10"
"""
BC = """inlet:  {{physID: 1, kind: inlet_uniformVelocity, outputHDFflg: 0, ints: , floats: {{ro: {ro:.6g}, Ux: {u:.6g}, Uy: 0.0, Uz: 0.0, Ps: {p:.6g}, k: {k}, omega: {om}}}}}
outlet: {{physID: 2, kind: outlet_statPress, outputHDFflg: 0, ints: , floats: {{Ps: {p:.6g}, Pt: {p:.6g}, Tt: {t:.6g}}}}}
top:    {{physID: 3, kind: slip, outputHDFflg: 0, ints: , floats: }}
wall:   {wallline}
sym:    {{physID: 5, kind: slip, outputHDFflg: 0, ints: , floats: }}
"""


def cfg(nsteps, cfl, relax, conv, lim, ninner, outint, plain, lam=False):
    t = CFG.replace('model: "sst"', 'model: "none"') if lam else CFG
    return t.format(cp=CP, gam=GAM, nsteps=nsteps, cfl=cfl, relax=relax, conv=conv, lim=lim, ninner=ninner,
                      outint=outint, dil=0 if plain else 2, kl=0 if plain else 1, k=K_INF, om=OM_INF)


def wall_line(wall):
    if wall == "adiabatic":
        return "{physID: 4, kind: wall, outputHDFflg: 1, ints: , floats: {Ux: 0.0, Uy: 0.0, Uz: 0.0}}"
    return f"{{physID: 4, kind: wall_isothermal, outputHDFflg: 1, ints: , floats: {{Ux: 0.0, Uy: 0.0, Uz: 0.0, Ts: {float(wall)}}}}}"


def patch_ic(h5, Tw):
    """一様自由流 + 壁ノード (wall_dist=0, x>=0) は u=0・T=Tw (等温) で IC を書く。"""
    with h5py.File(h5, "r+") as f:
        n = f["/VALUE/ro"].shape[0]
        c = f["/MESH/COORD"][:].reshape(-1, 3)
        wd = f["/VALUE/wall_dist"][:]
        ro = np.full(n, RO_INF); u = np.full(n, U_INF); T = np.full(n, T_INF)
        wall = (wd <= 0.0) & (c[:, 0] >= -1e-9)
        # 平板上 (x>=0) は壁距離 tanh(y/δ0) で速度を落とす (前縁の衝撃的立ち上がりと第一ノードの巨大せん断を避ける)
        onplate = c[:, 0] >= -1e-9
        u[onplate] = U_INF * np.tanh(np.maximum(wd[onplate], 0.0) / IC_DELTA0)
        u[wall] = 0.0
        if Tw is not None:
            T[wall] = Tw; ro[wall] = P_INF / (R * Tw)
        roe = P_INF / (GAM - 1) + 0.5 * ro * u ** 2
        f["/VALUE/ro"][:] = ro.astype(np.float32)
        f["/VALUE/roUx"][:] = (ro * u).astype(np.float32)
        f["/VALUE/roUy"][:] = np.zeros(n, np.float32); f["/VALUE/roUz"][:] = np.zeros(n, np.float32)
        f["/VALUE/roe"][:] = roe.astype(np.float32)
        for name, v in (("roK", ro * K_INF), ("roOmega", ro * OM_INF)):
            if f"/VALUE/{name}" in f: f[f"/VALUE/{name}"][:] = v.astype(np.float32)
            else: f.create_dataset(f"/VALUE/{name}", data=v.astype(np.float32))
        print(f"IC: n={n} wall nodes={int(wall.sum())} Tw={Tw}")


def set_wall_T(h5, Tw):
    """warm start 場の壁ノード (x>=0) を T=Tw・u=0 に置き換える (P 不変)。"""
    with h5py.File(h5, "r+") as f:
        c = f["/MESH/COORD"][:].reshape(-1, 3); wd = f["/VALUE/wall_dist"][:]
        ro = f["/VALUE/ro"][:].astype(float); rou = f["/VALUE/roUx"][:].astype(float); rov = f["/VALUE/roUy"][:].astype(float); roe = f["/VALUE/roe"][:].astype(float)
        wall = (wd <= 0.0) & (c[:, 0] >= -1e-9)
        P = (GAM - 1) * (roe[wall] - 0.5 * (rou[wall] ** 2 + rov[wall] ** 2) / ro[wall])
        ro[wall] = P / (R * Tw); rou[wall] = 0.0; rov[wall] = 0.0; roe[wall] = P / (GAM - 1)
        f["/VALUE/ro"][:] = ro.astype(np.float32); f["/VALUE/roUx"][:] = rou.astype(np.float32); f["/VALUE/roUy"][:] = rov.astype(np.float32); f["/VALUE/roe"][:] = roe.astype(np.float32)
        print(f"set_wall_T: {int(wall.sum())} wall nodes -> {Tw} K")


def run_forge(rd):
    r = subprocess.run([str(TOOLS / "run_case.sh"), str(rd)], env=ENV, capture_output=True, text=True)
    (rd / "run_case_stdout.log").write_text(r.stdout + r.stderr)
    return r.returncode


def stage(rd, text, nsteps, tag=""):
    (rd / "solverConfig.yaml").write_text(text)
    rc = run_forge(rd)
    res = sorted(rd.glob("res_[0-9]*.h5"), key=lambda f: int(f.stem.split("_")[1]))
    if rc != 0 or not res or int(res[-1].stem.split("_")[1]) < nsteps:
        raise SystemExit(f"stage failed rc={rc} res={[r.name for r in res][-2:]}")
    subprocess.run([sys.executable, str(TOOLS / "interp_field.py"), str(res[-1]), str(rd / "mesh.h5")], env=ENV, check=True,
                   capture_output=True, text=True)
    # 層流段の res には roK/roOmega が無い → 自由流値を入れ直す (SST 段の IC)
    with h5py.File(rd / "mesh.h5", "r+") as f:
        ro = f["/VALUE/ro"][:].astype(float)
        for name, val in (("roK", K_INF), ("roOmega", OM_INF)):
            if f"/VALUE/{name}" not in f or not np.all(np.isfinite(f[f"/VALUE/{name}"][:])) or float(np.abs(f[f"/VALUE/{name}"][:]).max()) == 0.0:
                if f"/VALUE/{name}" in f: del f[f"/VALUE/{name}"]
                f.create_dataset(f"/VALUE/{name}", data=(ro * val).astype(np.float32))
    for f in rd.glob("res_*"): f.unlink()
    for f in ("residual_history.csv", "residual_history.png", "CONVERGENCE_VERDICT.txt", "forge_run.log"):
        if (rd / f).exists(): shutil.move(str(rd / f), str(rd / f"{Path(f).stem}_{tag}{Path(f).suffix}"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True); ap.add_argument("--mesh", required=True); ap.add_argument("--wall", required=True)
    ap.add_argument("--plain", action="store_true"); ap.add_argument("--main-steps", type=int, default=12000)
    ap.add_argument("--cfl", type=float, default=4.0); ap.add_argument("--relax", type=float, default=0.7)
    ap.add_argument("--out-int", type=int, default=2000); ap.add_argument("--dry", action="store_true")
    ap.add_argument("--stages", default="full", choices=["full", "soft", "none"])
    ap.add_argument("--ramp", default="0.5,1,2", help="本段前の 2 次 cfl ランプ (各 --ramp-steps)。空文字で無し")
    ap.add_argument("--ramp-steps", type=int, default=2000)
    ap.add_argument("--limiter", type=int, default=2)
    ap.add_argument("--ic-from", default=None, help="warm start 元 run (mesh.h5 を interp_field で作る; 同一メッシュ)")
    ap.add_argument("--ic-mesh", default=None, help="場入り mesh h5 をそのまま mesh.h5 に使う (段階起動済み場の再利用)")
    a = ap.parse_args()
    rd = HERE / a.run
    if rd.exists(): raise SystemExit(f"{rd} exists")
    rd.mkdir()
    shutil.copy(HERE / "mesh" / f"{a.mesh}.h5", rd / "mesh.h5")
    Tw = None if a.wall == "adiabatic" else float(a.wall)
    if a.ic_mesh:
        shutil.copy(a.ic_mesh, rd / "mesh.h5"); (rd / "CONTINUED_FROM").write_text(str(a.ic_mesh) + "\n")
    elif a.ic_from:
        src = sorted((HERE / a.ic_from).glob("res_[0-9]*.h5"), key=lambda f: int(f.stem.split("_")[1]))[-1]
        subprocess.run([sys.executable, str(TOOLS / "interp_field.py"), str(src), str(rd / "mesh.h5")], env=ENV, check=True,
                       capture_output=True, text=True)
        (rd / "CONTINUED_FROM").write_text(str(src) + "\n")
        if Tw is not None: set_wall_T(rd / "mesh.h5", Tw)
    else:
        patch_ic(rd / "mesh.h5", Tw)
    (rd / "bcondConfig.yaml").write_text(BC.format(ro=RO_INF, u=U_INF, p=P_INF, t=T_INF, k=K_INF, om=OM_INF, wallline=wall_line(a.wall)))
    main_cfg = cfg(a.main_steps, a.cfl, a.relax, 1, a.limiter, 5, a.out_int, a.plain)
    (rd / "solverConfig.yaml").write_text(main_cfg)
    (rd / "GEN_ARGS").write_text(" ".join(sys.argv[1:]) + "\n")
    (rd / "probe.yaml").write_text("outStepInterval: 100\noutStepStart: 0\npoints:\nsurfaces:\n")
    if a.dry: return
    if a.stages == "full":
        # 層流暖機 (SST なし, 1 次, cfl 0.2) → SST soft (1 次, cfl 0.3) → mid (1 次, cfl 1.0) → 本段
        stage(rd, cfg(2000, 0.2, a.relax, 0, 0, 10, 2000, a.plain, lam=True), 2000, tag="lam")
        stage(rd, cfg(2000, 0.3, a.relax, 0, 0, 10, 2000, a.plain), 2000, tag="soft")
        stage(rd, cfg(2000, 1.0, a.relax, 0, 0, 10, 2000, a.plain), 2000, tag="mid")
    if a.stages == "soft":
        # 壁温切替 (同一メッシュ warm start) 用: SST 1 次 soft → mid だけ (plan §4.4 (a))
        stage(rd, cfg(2000, 0.5, a.relax, 0, 0, 10, 2000, a.plain), 2000, tag="soft")
        stage(rd, cfg(2000, 1.0, a.relax, 0, 0, 10, 2000, a.plain), 2000, tag="mid")
    if a.ramp:
        for i, cv in enumerate([float(v) for v in a.ramp.split(",") if v]):
            stage(rd, cfg(a.ramp_steps, cv, a.relax, 1, a.limiter, 5, a.ramp_steps, a.plain), a.ramp_steps, tag=f"ramp{i}_cfl{cv:g}")
    (rd / "solverConfig.yaml").write_text(main_cfg)
    rc = run_forge(rd)
    print("main rc", rc); print((rd / "CONVERGENCE_VERDICT.txt").read_text()[-600:])


if __name__ == "__main__":
    main()
