#!/usr/bin/env python3
"""case/55 の run 生成 — 乾燥空気 TP (4 種) + 低 Re SST + 段階起動 (case/48 の M4.19 レシピ準拠)。

物性は tools/gas_model.py と**同じ経路**にする (plan §4.3):
  熱力学  : forge `thermalMethod: 2` (NASA-9, speciesDBFile) ← 同じ DB
  輸送    : forge `viscMethod: 2` (Chapman–Enskog + Wilke/Mason–Saxena) ← 同じモデル族
  λ       : viscMethod 2 が Mason–Saxena で直接出す (thermCondMethod は使わない)
usage:
  python3 tools/make_case.py --run run_0001_T0_A1 --mesh t0_wd0.063_y4um --series A1 [--main-steps 8000]
"""
import argparse, json, os, shutil, subprocess, sys
from pathlib import Path
import numpy as np
import h5py, yaml

HERE = Path(__file__).resolve().parent
CASE = HERE.parent
ROOT = CASE.parents[1]
TOOLS = ROOT / "solver_density_cuda" / "tools"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "design"))
from gas_air import dry_air                                                 # noqa: E402
from conditions import solve_state                                          # noqa: E402
from forge_design.gas.composition import ResolvedSpeciesDB                  # noqa: E402

ENV = dict(os.environ, LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu/hdf5/serial:" + os.environ.get("LD_LIBRARY_PATH", ""))
SPECIES = ["N2", "O2", "AR", "CO2"]
IC_DELTA0 = 2.0e-4          # IC の壁近傍速度ランプ幅 [m]

CFG = """mesh: {{discretization: "node", nodeWallDirichlet: 1, nodeInletCornerWall: 1, meshFileName: "mesh.h5", valueFileName: "mesh.h5"}}
gpu: 1
solver: "SLAU"
physProp:
  thermalMethod: 2
  viscMethod: 2
  visc: {mu_inf:.6e}
  thermCond: {lam_inf:.6e}
  cp: {cp_ref:.2f}
  gamma: {gam_ref:.5f}
  species: [{species}]
  speciesDBFile: species_db.yaml
  thermoHrefTemp: 298.15
time:
  unsteady: 0
  dualTime: 0
  last: {{nStepOuter: {nsteps}}}
  deltaT: {{control: 1, dt: 1.0e-9, cfl: {cfl}, cfl_pseudo: {cfl}, implicitRelax: {relax}, blockDPLUR: 1,
    dt_min: 1.0e-12, dt_max: 1.0, detectNaN: 1, monitorInterval: 100}}
  outStepStart: 0
  outStepInterval: {out_int}
  timeIntegration: 11
  nStepInner: {inner}
space: {{convMethod: {conv}, limiter: {lim}}}
turbulence: {{model: "{turb}", scalarDiffusion: 1, dilatationCorrection: 2, katoLaunder: 1,
             wallTreatmentSST: 0, turbulentPrandtl: 0.9, kInit: {kinf:.4f}, omegaInit: {ominf:.2f}}}
initial: "uniform_p101325_u10"
output: {{level: 1, extraFields: [thermCond, vis_lam]}}
"""

BC_PLATE = {
    "adiabatic": '{{physID: 4, kind: wall, outputHDFflg: 1, ints: , floats: {{Ux: 0.0, Uy: 0.0, Uz: 0.0}}}}',
    "isothermal": '{{physID: 4, kind: wall_isothermal, outputHDFflg: 1, ints: , floats: {{Ux: 0.0, Uy: 0.0, Uz: 0.0, Ts: {tw:.2f}}}}}',
}

BC = """inlet:  {{physID: 1, kind: inlet_uniformVelocity, outputHDFflg: 0, ints: , floats: {{ro: {ro:.8e}, Ux: {u:.6f}, Uy: 0.0, Uz: 0.0, Ps: {p:.6f}, k: {kinf:.4f}, omega: {ominf:.2f}{ys}}}}}
outlet: {{physID: 2, kind: outlet_statPress, outputHDFflg: 0, ints: , floats: {{Ps: {p:.6f}, Pt: {p:.6f}, Tt: {t:.4f}}}}}
top:    {{physID: 3, kind: slip, outputHDFflg: 0, ints: , floats: }}
plate:  {plate_bc}
runup:  {{physID: 5, kind: slip, outputHDFflg: 0, ints: , floats: }}
cavity: {{physID: 6, kind: wall_isothermal, outputHDFflg: 1, ints: , floats: {{Ux: 0.0, Uy: 0.0, Uz: 0.0, Ts: {tw:.2f}}}}}
"""


def write_species_db(run_dir):
    db = ResolvedSpeciesDB.builtin()
    out = {s: db[s].to_db_dict() for s in SPECIES}
    (run_dir / "species_db.yaml").write_text(yaml.safe_dump(out, sort_keys=False), encoding="utf-8")


def patch_ic(h5, st, gas, Tw, kinf, ominf, plate_thermal="adiabatic"):
    """一様自由流 + 壁ノード (wall_dist=0, x>=0) は u=0・T=Tw。roY も書く。"""
    with h5py.File(h5, "r+") as f:
        n = f["/VALUE/ro"].shape[0]
        c = f["/MESH/COORD"][:].reshape(-1, 3)
        wd = f["/VALUE/wall_dist"][:]
        ro = np.full(n, st["rho"]); u = np.full(n, st["U"]); T = np.full(n, st["T"])
        onplate = c[:, 0] >= -1e-9
        u[onplate] = st["U"] * np.tanh(np.maximum(wd[onplate], 0.0) / IC_DELTA0)
        wall = (wd <= 0.0) & onplate
        u[wall] = 0.0
        # 平板は断熱、キャビティ壁だけ Tw (y<0 の壁ノード)。IC では区別せず自由流温度のままにし、
        # キャビティ内の壁ノードだけ Tw に落とす (等温 BC が step 0 から効く)
        incav = wall & (c[:, 1] < -1e-9) if plate_thermal == "adiabatic" else wall
        T[incav] = Tw
        ro[incav] = st["p"] / (gas.R * Tw)
        e = np.array([gas.h(t) for t in np.unique(T)])
        emap = {t: ev - gas.R * t for t, ev in zip(np.unique(T), e)}
        e_node = np.array([emap[t] for t in T])
        roe = ro * (e_node + 0.5 * u ** 2)
        f["/VALUE/ro"][:] = ro.astype(np.float32)
        f["/VALUE/roUx"][:] = (ro * u).astype(np.float32)
        f["/VALUE/roUy"][:] = np.zeros(n, np.float32)
        f["/VALUE/roUz"][:] = np.zeros(n, np.float32)
        f["/VALUE/roe"][:] = roe.astype(np.float32)
        for name, val in (("roK", kinf), ("roOmega", ominf)):
            v = (ro * val).astype(np.float32)
            if f"/VALUE/{name}" in f:
                f[f"/VALUE/{name}"][:] = v
            else:
                f.create_dataset(f"/VALUE/{name}", data=v)
        for i, s in enumerate(SPECIES):
            name = f"/VALUE/roY{i}"
            v = (ro * gas.Y[s]).astype(np.float32)
            if name in f:
                f[name][:] = v
            else:
                f.create_dataset(name, data=v)
        print(f"IC: n={n}, 壁ノード {int(wall.sum())}, T∞={st['T']:.2f} K, U∞={st['U']:.1f} m/s, Tw={Tw} K")


def cfg_text(**kw):
    return CFG.format(species=", ".join(SPECIES), **kw)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    ap.add_argument("--mesh", required=True)
    ap.add_argument("--series", default="A1")
    ap.add_argument("--main-steps", type=int, default=8000)
    ap.add_argument("--soft-steps", type=int, default=1500)
    ap.add_argument("--mid-steps", type=int, default=1500)
    ap.add_argument("--soft-cfl", type=float, default=0.5)
    ap.add_argument("--mid-cfl", type=float, default=1.0)
    ap.add_argument("--ic-delta", type=float, default=2.0e-4, help="IC の壁近傍ランプ幅 [m]")
    ap.add_argument("--cfl", type=float, default=2.0)
    ap.add_argument("--ramp", default="0.5,1,2", help="2 次移行の cfl ランプ (空文字で無し)")
    ap.add_argument("--ramp-steps", type=int, default=1000)
    ap.add_argument("--relax", type=float, default=0.7)
    ap.add_argument("--out-int", type=int, default=2000)
    ap.add_argument("--plate-thermal", choices=["adiabatic", "isothermal"], default="adiabatic",
                    help="平板の熱条件。分母 (h_fp) を作る smooth run は isothermal で回す")
    ap.add_argument("--dry", action="store_true")
    a = ap.parse_args()

    cond = json.loads((CASE / "conditions.json").read_text(encoding="utf-8"))
    s = next(x for x in cond["series"] if x["id"] == a.series)
    Tw = cond["wall"]["Tw_cavity"]
    gas = dry_air(s["Tt_tp"], cond["gas"]["Y"])
    T, U, rho, p, mu = s["T_inf"], s["U_inf"], s["ro_inf"], s["P_inf"], s["mu_inf"]
    st = dict(T=T, U=U, rho=rho, p=p, mu=mu)

    rd = CASE / a.run
    rd.mkdir(exist_ok=True)
    for f in rd.glob("res_*"):
        f.unlink()
    for f in ("forge_run.log", "residual_history.csv"):
        if (rd / f).exists():
            (rd / f).unlink()
    shutil.copy(CASE / "mesh" / f"{a.mesh}.h5", rd / "mesh.h5")
    write_species_db(rd)
    (rd / "probe.yaml").write_text("outStepInterval: 100\noutStepStart: 0\npoints:\nsurfaces:\n")
    global IC_DELTA0
    IC_DELTA0 = a.ic_delta
    patch_ic(rd / "mesh.h5", st, gas, Tw, cond["turbulence"]["k_inf"], cond["turbulence"]["omega_inf"],
             plate_thermal=a.plate_thermal)
    ys = "".join(f", Y{i}: {gas.Y[sp]:.8f}" for i, sp in enumerate(SPECIES))
    plate_bc = BC_PLATE[a.plate_thermal].format(tw=Tw)
    (rd / "bcondConfig.yaml").write_text(BC.format(ro=rho, u=U, p=p, t=T, tw=Tw, ys=ys, plate_bc=plate_bc,
                                               kinf=cond["turbulence"]["k_inf"],
                                               ominf=cond["turbulence"]["omega_inf"]))
    (rd / "case_setup.json").write_text(json.dumps(
        dict(series=s, T_inf=T, U_inf=U, rho_inf=rho, p_inf=p, mu_inf=mu, Tw=Tw, phi=gas.phi,
             Y=gas.Y, R=gas.R, mesh=a.mesh, Taw_tp=s["Taw_tp"], plate_thermal=a.plate_thermal),
        indent=2, ensure_ascii=False), encoding="utf-8")

    common = dict(mu_inf=mu, lam_inf=gas.lam(T), cp_ref=gas.cp(T), gam_ref=gas.gamma(T), relax=1.0,
                  kinf=cond["turbulence"]["k_inf"], ominf=cond["turbulence"]["omega_inf"], turb="sst")
    stages = [
        ("lam",  dict(conv=0, lim=0, cfl=0.2, inner=10, nsteps=a.soft_steps,
                      out_int=a.soft_steps, **{**common, "turb": "none"})),
        ("soft", dict(conv=0, lim=0, cfl=a.soft_cfl, inner=10, nsteps=a.soft_steps,
                      out_int=a.soft_steps, **common)),
        ("mid",  dict(conv=0, lim=0, cfl=a.mid_cfl, inner=10, nsteps=a.mid_steps,
                      out_int=a.mid_steps, **common)),
    ]
    for i, c in enumerate([float(x) for x in a.ramp.split(",") if x.strip()]):
        stages.append((f"ramp{i}_cfl{c:g}", dict(conv=1, lim=2, cfl=c, inner=6, nsteps=a.ramp_steps,
                                                 out_int=a.ramp_steps, **common)))
    stages.append(("main", dict(conv=1, lim=2, cfl=a.cfl, inner=4, nsteps=a.main_steps,
                                out_int=a.out_int, **{**common, "relax": a.relax})))
    (rd / "stages.json").write_text(json.dumps([{k: v for k, v in st_.items()} | {"tag": tag}
                                                for tag, st_ in stages], indent=2), encoding="utf-8")
    if a.dry:
        (rd / "solverConfig.yaml").write_text(cfg_text(**stages[0][1]))
        print("dry: 段構成のみ書き出し")
        return
    for tag, kw in stages:
        (rd / "solverConfig.yaml").write_text(cfg_text(**kw))
        print(f"--- stage {tag}: cfl={kw['cfl']} conv={kw['conv']} → step {kw['nsteps']}")
        r = subprocess.run([str(TOOLS / "run_case.sh"), str(rd)], env=ENV, capture_output=True, text=True)
        (rd / f"run_stdout_{tag}.log").write_text(r.stdout + r.stderr)
        res = sorted(rd.glob("res_[0-9]*.h5"), key=lambda f: int(f.stem.split("_")[1]))
        last = int(res[-1].stem.split("_")[1]) if res else -1
        if r.returncode != 0 or last < kw["nsteps"]:
            print(r.stdout[-3000:])
            raise SystemExit(f"stage {tag} failed rc={r.returncode} last={last}")
        # 次段の IC は最終 res を同一メッシュへ index コピー
        subprocess.run([sys.executable, str(TOOLS / "interp_field.py"), str(res[-1]), str(rd / "mesh.h5")],
                       env=ENV, check=True, capture_output=True, text=True)
        for f in rd.glob("res_*"):
            f.unlink()
        for f in ("residual_history.csv", "residual_history.png", "CONVERGENCE_VERDICT.txt", "forge_run.log"):
            if (rd / f).exists():
                shutil.move(str(rd / f), str(rd / f"{Path(f).stem}_{tag}{Path(f).suffix}"))
    print("done:", rd)


if __name__ == "__main__":
    main()
