#!/usr/bin/env python3
"""case/56 の run 生成 — TP-1187 較正パネル相当の 2D 平板 (分母 q_FP を forge で作る)。

熱力学/輸送は case/50 と同じ (`thermalMethod: 2` NASA-9, `viscMethod: 2` Chapman-Enskog +
Wilke/Mason-Saxena)。燃焼ガスの組成・当量比も case/50 の `CombustionProducts` から取る。
自由流は `derived.json` (Gate A: M=7 固定、q から試験部全圧を解いたもの)。

段階起動は case/55 と同じ思想: 層流暖機 → soft → mid → 2 次ランプ → 本段。
"""
import argparse, json, subprocess, sys
from pathlib import Path
import numpy as np, h5py

HERE = Path(__file__).resolve().parent
CASE = HERE.parent
ROOT = CASE.parents[1]
C50 = ROOT / "case" / "50.deep_cavity_wieting_m7" / "tools"
TOOLS = ROOT / "solver_density_cuda" / "tools"
sys.path.insert(0, str(HERE))
# 同名 (make_case) なので明示パスで読む (循環 import を避ける)
import importlib.util as _ilu                                    # noqa: E402
_sp = _ilu.spec_from_file_location("c50_make_case", C50 / "make_case.py")
m50 = _ilu.module_from_spec(_sp); _sp.loader.exec_module(m50)
from gas_htst import htst                                        # noqa: E402

SPECIES = m50.SPECIES
IC_DELTA0 = 2.0e-3

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
output: {{level: 1, extraFields: [thermCond, vis_lam, vis_turb]}}
"""

BC = """inlet:  {{physID: 1, kind: inlet_uniformVelocity, outputHDFflg: 0, ints: , floats: {{ro: {ro:.8f}, Ux: {U:.4f}, Uy: 0.0, Uz: 0.0, Ps: {p:.6f}, k: {kinf:.4f}, omega: {ominf:.2f}, {ymf}}}}}
outlet: {{physID: 2, kind: outlet_statPress, outputHDFflg: 0, ints: , floats: {{Ps: {p:.6f}, Pt: {p:.6f}, Tt: {t:.4f}}}}}
top:    {{physID: 3, kind: slip, outputHDFflg: 0, ints: , floats: }}
plate:  {{physID: 4, kind: wall_isothermal, outputHDFflg: 1, ints: , floats: {{Ux: 0.0, Uy: 0.0, Uz: 0.0, Ts: {tw:.2f}}}}}
slip:   {{physID: 5, kind: slip, outputHDFflg: 0, ints: , floats: }}
"""


def seed_turb(h5, k_inf, om_inf):
    """roK / roOmega が無い or 全ゼロなら自由流値で埋める。

    case/50 の patch_ic は層流ケース用で roK/roOmega を書かない。さらに**層流暖機段の出力には
    そもそも roK/roOmega が存在しない**ので、そのまま SST 段に渡すと k=0 から始まり step 3 で
    発散する (2026-09-20 に踏んだ)。段ごとの引き継ぎの後に必ず呼ぶ。
    """
    with h5py.File(h5, "r+") as f:
        ro = f["/VALUE/ro"][:].astype(float)
        for name, val in (("roK", k_inf), ("roOmega", om_inf)):
            key = f"/VALUE/{name}"
            v = (ro * val).astype(np.float32)
            if key not in f:
                f.create_dataset(key, data=v); print(f"  seed {name} (新規)")
            elif float(np.max(np.abs(f[key][:]))) == 0.0:
                f[key][:] = v; print(f"  seed {name} (全ゼロだった)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    ap.add_argument("--mesh", default="fp_t8")
    ap.add_argument("--series", type=int, default=8, help="TP-1187 の run 番号")
    ap.add_argument("--tw", type=float, default=300.0)
    ap.add_argument("--laminar", action="store_true", help="乱流モデルを使わない (層流 run 用)")
    ap.add_argument("--soft-steps", type=int, default=3000)
    ap.add_argument("--mid-steps", type=int, default=3000)
    ap.add_argument("--ramp", default="0.3,0.6,1.0")
    ap.add_argument("--ramp-steps", type=int, default=2500)
    ap.add_argument("--soft-cfl", type=float, default=0.2)
    ap.add_argument("--mid-cfl", type=float, default=0.5)
    ap.add_argument("--cfl", type=float, default=1.5)
    ap.add_argument("--main-steps", type=int, default=40000)
    ap.add_argument("--out-int", type=int, default=5000)
    ap.add_argument("--dry", action="store_true")
    a = ap.parse_args()

    der = {o["run"]: o for o in json.loads((CASE / "derived.json").read_text(encoding="utf-8"))["series"]}
    o = der[a.series]
    gas = htst(o["Tt_c"])
    T, U, rho, p = o["T_inf"], o["U_inf"], o["ro_inf"], o["p_inf"]
    mu, lam = gas.mu(T), gas.lam(T)
    # 乱流入口: Tu=0.5 %, mu_t/mu = 10
    k_inf = 1.5 * (0.005 * U) ** 2
    om_inf = rho * k_inf / (10.0 * mu)

    rd = CASE / a.run; rd.mkdir(exist_ok=True)
    for f in list(rd.glob("res_*")) + list(rd.glob("*.log")) + list(rd.glob("residual_history*")):
        f.unlink()
    import shutil
    shutil.copy(CASE / "mesh" / f"{a.mesh}.h5", rd / "mesh.h5")
    m50.write_species_db(rd)
    (rd / "probe.yaml").write_text(
        "outStepInterval: 100\noutStepStart: 0\npoints:\nsurfaces:\n", encoding="utf-8")
    ymf = ", ".join(f"Y{i}: {gas.Y[sp]:.8f}" for i, sp in enumerate(SPECIES))
    (rd / "bcondConfig.yaml").write_text(BC.format(
        ro=rho, U=U, p=p, t=T, tw=a.tw, kinf=k_inf, ominf=om_inf, ymf=ymf), encoding="utf-8")
    (rd / "case_setup.json").write_text(json.dumps(
        dict(tp1187_run=a.series, M=o["M"], T_inf=T, U_inf=U, rho_inf=rho, p_inf=p,
             mu_inf=mu, T_aw=o["T_aw"], Tt_c=o["Tt_c"], Tw=a.tw, mesh=a.mesh,
             phi=gas.phi, Y=gas.Y, R=gas.R, laminar=bool(a.laminar),
             x_eval_m={"I": 1.17, "II": 1.88}, x_trip_m=0.13),
        indent=2, ensure_ascii=False), encoding="utf-8")

    turb = "none" if a.laminar else "sst"
    common = dict(mu_inf=mu, lam_inf=gas.lam(T), cp_ref=gas.cp(T), gam_ref=gas.gamma(T),
                  relax=1.0, kinf=k_inf, ominf=om_inf, turb=turb)
    stages = [("lam",  dict(conv=0, lim=0, cfl=0.2, inner=10, nsteps=a.soft_steps,
                            out_int=a.soft_steps, **{**common, "turb": "none"})),
              ("soft", dict(conv=0, lim=0, cfl=a.soft_cfl, inner=10, nsteps=a.soft_steps,
                            out_int=a.soft_steps, **common)),
              ("mid",  dict(conv=0, lim=0, cfl=a.mid_cfl, inner=10, nsteps=a.mid_steps,
                            out_int=a.mid_steps, **common))]
    for i, c in enumerate([float(x) for x in a.ramp.split(",") if x.strip()]):
        stages.append((f"ramp{i}_cfl{c:g}", dict(conv=1, lim=2, cfl=c, inner=6,
                                                 nsteps=a.ramp_steps, out_int=a.ramp_steps, **common)))
    stages.append(("main", dict(conv=1, lim=2, cfl=a.cfl, inner=4, nsteps=a.main_steps,
                                out_int=a.out_int, **common)))
    (rd / "stages.json").write_text(json.dumps(
        [{**st, "tag": tag} for tag, st in stages], indent=2), encoding="utf-8")

    st0 = dict(rho=rho, U=U, T=T, p=p)
    m50.IC_DELTA0 = IC_DELTA0
    m50.patch_ic(rd / "mesh.h5", st0, gas, a.tw)
    seed_turb(rd / "mesh.h5", k_inf, om_inf)
    if a.dry:
        (rd / "solverConfig.yaml").write_text(CFG.format(species=", ".join(SPECIES), **stages[0][1]))
        print("dry: 段構成のみ"); return
    for tag, kw in stages:
        (rd / "solverConfig.yaml").write_text(CFG.format(species=", ".join(SPECIES), **kw))
        print(f"--- stage {tag}: cfl={kw['cfl']} conv={kw['conv']} turb={kw['turb']} → {kw['nsteps']} step")
        subprocess.run([str(TOOLS / "run_case.sh"), str(rd)], check=False)
        for suf in ("forge_run.log", "residual_history.csv", "CONVERGENCE_VERDICT.txt"):
            src = rd / suf
            if src.exists():
                src.rename(rd / f"{Path(suf).stem}_{tag}{Path(suf).suffix}")
        res = sorted(rd.glob("res_[0-9]*.h5"), key=lambda f: int(f.stem.split("_")[1]))
        if res:
            subprocess.run([sys.executable, str(TOOLS / "interp_field.py"), str(res[-1]),
                            str(rd / "mesh.h5")], check=True)
            seed_turb(rd / "mesh.h5", k_inf, om_inf)


if __name__ == "__main__":
    main()
