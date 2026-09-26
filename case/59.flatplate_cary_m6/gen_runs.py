#!/usr/bin/env python3
"""case/59 (Cary TN D-5863 冷却平板, M6 空気) の forge run を生成して段階起動する。

    python3 gen_runs.py --run run_0001_tw02_re027 --series Re0.27_Tw0.2 [--mesh fp_y1_1.5um_nx900] [--dry]

- 自由流: M 6.02、Tt 533 K (γ 1.4 の等エントロピー → T∞ 64.6 K)。**p∞ は系列の単位 Re (Table II) に
  forge 自身の Sutherland μ (gasProperties_d.cu: T0 273 K, μ0 1.716e-5, S 111 K) で合わせる**。
- 壁: 等温 Tw = (Tw/Tt)·533 K。熱量的完全の空気、Pr 0.72 一定、Prₜ 0.9、SST 補正は case/48 と同じで全系列固定。
- 段: lam (1 次, cfl 0.2) → soft (SST 1 次, 0.3) → mid (1 次, 1.0) → 2 次ランプ 0.5/1/2 → main (2 次, cfl 2, relax 0.7)。
  段の引き継ぎは **restart_field.py** (同一メッシュ。case/48 の interp_field.py は AGENTS で禁止された使い方)。
  段ごとの実効設定を `stage_manifest.json` に書く (check_convergence.py --segment 用)。
- `physProp.visc` は剛性見積りにしか入らないので、冷たい自由流ではなく境界層の代表温度 (Taw) の μ を入れる。
"""
import argparse, json, os, shutil, subprocess, sys
from pathlib import Path
import h5py, numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TOOLS = ROOT / "solver_density_cuda" / "tools"
sys.path.insert(0, str(TOOLS))
from stage_manifest import StageManifest                        # noqa: E402

ENV = dict(os.environ, FORGE_CUDA_BLOCKSIZE="128",
           LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu/hdf5/serial:" + os.environ.get("LD_LIBRARY_PATH", ""))
GAM, CP = 1.4, 1004.5
R = CP * (GAM - 1) / GAM
M_INF, TT = 6.02, 533.0
T_INF = TT / (1 + 0.5 * (GAM - 1) * M_INF ** 2)
A_INF = (GAM * R * T_INF) ** 0.5
U_INF = M_INF * A_INF
R_TURB = 0.89                                   # Cary の St 還元に使われた乱流回復係数
IC_DELTA0 = 1.0e-3


def mu_suth(T):
    return 1.716e-5 * (T / 273.0) ** 1.5 * (273.0 + 111.0) / (T + 111.0)


CFG = """mesh: {{discretization: "node", nodeWallDirichlet: 1, nodeInletCornerWall: 1, meshFileName: "mesh.h5", valueFileName: "mesh.h5"}}
gpu: 1
solver: "SLAU"
physProp: {{thermalMethod: 0, viscMethod: 1, visc: {visc:.4e}, thermCond: 0.0257, thermCondMethod: 1, prandtlLam: 0.72, cp: {cp}, gamma: {gam}}}
time:
  unsteady: 0
  dualTime: 0
  last: {{nStepOuter: {nsteps}}}
  deltaT: {{control: 1, dt: 1e-8, cfl: {cfl}, cfl_pseudo: {cfl}, implicitRelax: {relax}, blockDPLUR: 1, lowMachPrecond: 0, dt_min: 1e-11, dt_max: 1.0, detectNaN: 1, qAccumulatorFP64: 1}}
  outStepStart: 0
  outStepInterval: {outint}
  timeIntegration: 11
  nStepInner: {ninner}
space: {{convMethod: {conv}, limiter: {lim}, limiterRoRef: {ro:.8g}, limiterPRef: {p:.8g}, limiterARef: {a:.8g}}}
turbulence: {{model: "{model}", scalarDiffusion: 1, dilatationCorrection: 2, katoLaunder: 1, wallTreatmentSST: 0, turbulentPrandtl: 0.9, kInit: {k:.6g}, omegaInit: {om:.6g}}}
output: {{level: 1}}
initial: "uniform_p101325_u10"
"""
BC = """inlet:  {{physID: 1, kind: inlet_uniformVelocity, outputHDFflg: 0, ints: , floats: {{ro: {ro:.8g}, Ux: {u:.8g}, Uy: 0.0, Uz: 0.0, Ps: {p:.8g}, k: {k:.6g}, omega: {om:.6g}}}}}
outlet: {{physID: 2, kind: outlet_statPress, outputHDFflg: 0, ints: , floats: {{Ps: {p:.8g}, Pt: {p:.8g}, Tt: {t:.8g}}}}}
top:    {{physID: 3, kind: slip, outputHDFflg: 0, ints: , floats: }}
wall:   {{physID: 4, kind: wall_isothermal, outputHDFflg: 1, ints: , floats: {{Ux: 0.0, Uy: 0.0, Uz: 0.0, Ts: {tw:.6g}}}}}
sym:    {{physID: 5, kind: slip, outputHDFflg: 0, ints: , floats: }}
"""


def freestream(re_per_m):
    mu = mu_suth(T_INF)
    ro = re_per_m * mu / U_INF
    p = ro * R * T_INF
    k = 1.5 * (0.005 * U_INF) ** 2                # Tu 0.5 %
    om = ro * k / (10.0 * mu)                     # μt/μ = 10
    return dict(ro=ro, p=p, mu=mu, k=k, om=om)


def patch_ic(h5, fs, Tw):
    with h5py.File(h5, "r+") as f:
        n = f["/VALUE/ro"].shape[0]
        c = f["/MESH/COORD"][:].reshape(-1, 3)
        wd = f["/VALUE/wall_dist"][:]
        onplate = c[:, 0] >= -1e-9
        u = np.full(n, U_INF); T = np.full(n, T_INF)
        s = np.tanh(np.maximum(wd, 0.0) / IC_DELTA0)
        u[onplate] = U_INF * s[onplate]
        T[onplate] = Tw + (T_INF - Tw) * s[onplate]   # 壁で Tw、外で T∞ (圧力は一様)
        ro = fs["p"] / (R * T)
        roe = fs["p"] / (GAM - 1) + 0.5 * ro * u ** 2
        f["/VALUE/ro"][:] = ro.astype(np.float32)
        f["/VALUE/roUx"][:] = (ro * u).astype(np.float32)
        f["/VALUE/roUy"][:] = 0.0
        f["/VALUE/roUz"][:] = 0.0
        f["/VALUE/roe"][:] = roe.astype(np.float32)
        seed_turb(f, fs, force=True)


def seed_turb(f, fs, force=False):
    ro = f["/VALUE/ro"][:].astype(float)
    for name, val in (("roK", fs["k"]), ("roOmega", fs["om"])):
        key = f"/VALUE/{name}"
        v = (ro * val).astype(np.float32)
        if key not in f:
            f.create_dataset(key, data=v)
        elif force or float(np.max(np.abs(f[key][:]))) == 0.0 or not np.all(np.isfinite(f[key][:])):
            f[key][:] = v


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    ap.add_argument("--series", required=True, help="conditions.json の series id (例 Re0.27_Tw0.2)")
    ap.add_argument("--mesh", default="fp_y1_1.5um_nx900")
    ap.add_argument("--main-steps", type=int, default=20000)
    ap.add_argument("--out-int", type=int, default=1000)
    ap.add_argument("--cfl", type=float, default=2.0)
    ap.add_argument("--relax", type=float, default=0.7)
    ap.add_argument("--init-from", default=None,
                    help="収束場の res_*.h5。別メッシュなら interp_field.py で移して段階起動を省き、本段だけを回す")
    ap.add_argument("--forge-tools", default=None,
                    help="run_case.sh のあるディレクトリ (別ビルド、例: 全域 FP64 の worktree)。既定はこのリポジトリ")
    ap.add_argument("--dry", action="store_true")
    a = ap.parse_args()

    cond = json.loads((HERE / "conditions.json").read_text())
    ser = next(s for s in cond["series"] if s["id"] == a.series)
    re_m = ser["Re_per_cm"] * 100.0
    Tw = ser["Tw_over_Tt"] * TT
    fs = freestream(re_m)
    Taw = T_INF * (1 + R_TURB * 0.5 * (GAM - 1) * M_INF ** 2)
    visc = mu_suth(Taw)

    rd = HERE / a.run
    if rd.exists():
        raise SystemExit(f"{rd} は既にある。消さない — 別の run 名にすること")
    rd.mkdir()
    shutil.copy(HERE / "mesh" / f"{a.mesh}.h5", rd / "mesh.h5")
    common = dict(visc=visc, cp=CP, gam=GAM, relax=a.relax, ro=fs["ro"], p=fs["p"], a=A_INF,
                  k=fs["k"], om=fs["om"])
    stages = [("lam", dict(nsteps=2000, cfl=0.2, conv=0, lim=0, ninner=10, outint=2000, model="none")),
              ("soft", dict(nsteps=2000, cfl=0.3, conv=0, lim=0, ninner=10, outint=2000, model="sst")),
              ("mid", dict(nsteps=2000, cfl=1.0, conv=0, lim=0, ninner=10, outint=2000, model="sst"))]
    for i, c in enumerate((0.5, 1.0, 2.0)):
        stages.append((f"ramp{i}_cfl{c:g}", dict(nsteps=2000, cfl=c, conv=1, lim=2, ninner=5, outint=2000, model="sst")))
    stages.append(("main", dict(nsteps=a.main_steps, cfl=a.cfl, conv=1, lim=2, ninner=5, outint=a.out_int, model="sst")))
    bc = BC.format(ro=fs["ro"], u=U_INF, p=fs["p"], t=T_INF, k=fs["k"], om=fs["om"], tw=Tw)
    (rd / "bcondConfig.yaml").write_text(bc)
    (rd / "probe.yaml").write_text("outStepInterval: 100\noutStepStart: 0\npoints:\nsurfaces:\n")
    sm = StageManifest(rd)
    for tag, st in stages:
        sm.add(tag, CFG.format(**common, **st), bc,
               history="residual_history.csv" if tag == "main" else f"residual_history_{tag}.csv")
    sm.write()
    setup = dict(series=a.series, Re_per_m=re_m, Tw_over_Tt=ser["Tw_over_Tt"], Tw=Tw, M=M_INF, Tt=TT,
                 T_inf=T_INF, U_inf=U_INF, **{k + "_inf": v for k, v in fs.items()}, Taw_r089=Taw,
                 visc_stiffness=visc, mesh=a.mesh, cuda_blocksize=128, gen_args=sys.argv[1:])
    (rd / "case_setup.json").write_text(json.dumps(setup, indent=2, ensure_ascii=False))
    (rd / "solverConfig.yaml").write_text(CFG.format(**common, **stages[0][1]))
    if a.init_from:
        stages = stages[-1:]                            # 本段だけ
        sm = StageManifest(rd)
        sm.add("main", CFG.format(**common, **stages[0][1]), bc, history="residual_history.csv",
               restart_from=a.init_from)
        sm.write()
        (rd / "solverConfig.yaml").write_text(CFG.format(**common, **stages[0][1]))
        subprocess.run([sys.executable, str(TOOLS / "interp_field.py"), a.init_from, str(rd / "mesh.h5")],
                       env=ENV, check=True)
        (rd / "CONTINUED_FROM").write_text(f"{a.init_from} (interp_field.py、別メッシュ)\n")
    else:
        patch_ic(rd / "mesh.h5", fs, Tw)
    print(f"{a.run}: Re {re_m:.3g}/m  p∞ {fs['p']:.1f} Pa  ρ∞ {fs['ro']:.4g}  U∞ {U_INF:.1f}  T∞ {T_INF:.2f} K  Tw {Tw:.1f} K  Taw(r0.89) {Taw:.1f} K")
    if a.dry:
        return
    for tag, st in stages:
        (rd / "solverConfig.yaml").write_text(CFG.format(**common, **st))
        print(f"--- {tag}: {st['nsteps']} step, cfl {st['cfl']}", flush=True)
        run_tools = Path(a.forge_tools) if a.forge_tools else TOOLS
        rc = subprocess.run([str(run_tools / "run_case.sh"), str(rd)], env=ENV).returncode
        cur = rd / f"res_{st['nsteps']}.h5"
        if rc != 0 or not cur.exists() or list(rd.glob("res_nan_*.h5")):
            raise SystemExit(f"stage {tag}: rc={rc} res={cur.exists()} → 中断 (古い場を引き継がない)")
        if tag == "main":
            break
        for suf in ("forge_run.log", "residual_history.csv", "CONVERGENCE_VERDICT.txt", "residual_history.png"):
            if (rd / suf).exists():
                (rd / suf).rename(rd / f"{Path(suf).stem}_{tag}{Path(suf).suffix}")
        subprocess.run([sys.executable, str(TOOLS / "restart_field.py"), str(cur), str(rd / "mesh.h5")],
                       env=ENV, check=True)
        with h5py.File(rd / "mesh.h5", "r+") as f:
            seed_turb(f, fs)
        for p in list(rd.glob("res_*")):
            p.rename(rd / f"_{tag}_{p.name}")          # 段の出力は消さずに退避 (main の step 番号と混ざらないように)
    print("main 終了")


if __name__ == "__main__":
    main()
