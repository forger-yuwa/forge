#!/usr/bin/env python3
"""case/49 のメッシュ取り込みと段階起動 (plan §4.2 / §4.5 / §4.6)。

usage:
  python3 gen_runs.py convert --msh cad/forge.msh --out mesh/stageA.h5
  python3 gen_runs.py run --run run_0001_stageA_cpg --mesh mesh/stageA.h5 \
         --inlet-csv precursor/run_0002_.../inlet_profile_1.csv [--gas CPG] [--main-steps 20000]

設定はすべて manifest.json (setup.py --resolve) から。physID は manifest の `phys_id`。

重要 (plan §4.2 / §4.6, codex plan-1 M5):
  - **変換は最終形の壁タグで行う** (SST の wall_dist は no-slip 壁から作られる)。
    段階起動 S0 の「全面 slip」は実行時に bcondConfig.yaml を差し替えて実現する。
  - `space.pRef` = P_inf (非直交 float32 の自由流保存)。
  - IC はキャビティ内 (z<0) を **静止・壁温** から始める (一様 M5 のままだと大きな速度過渡)。
  - 同一メッシュの段間引き継ぎは index コピー (**3D で interp_field.py は使わない**)。
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TOOLS = ROOT / "solver_density_cuda" / "tools"
BUILD = ROOT / "solver_density_cuda" / "build"
# FORGE_BIN は build-native に固定する (別セッションが build/ を使うため。
# 解像壁の qwall/utau 診断出力を入れたバイナリ)。
ENV = dict(os.environ,
           LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu/hdf5/serial:" + os.environ.get("LD_LIBRARY_PATH", ""),
           FORGE_BIN=os.environ.get("FORGE_BIN", str(ROOT / "solver_density_cuda" / "build-native" / "forge")))
sys.path.insert(0, str(HERE))
from setup import load as load_conditions  # noqa: E402

# manifest は環境変数で差し替えられる (塞ぎ形状 manifest_plug.json などの併存用)
MAN = json.loads((HERE / os.environ.get("CASE49_MANIFEST", "manifest.json")).read_text())
PID, G = MAN["phys_id"], MAN["geometry"]
D = load_conditions()

# 群は manifest に**実在するものだけ**にする (塞ぎ形状 plug_cavity ではキャビティ壁が無い)
_iso = ["cav_outer", "cav_floor", "cyl_side"] + (
    ["cyl_top"] if D["cyl_top_thermal"] == "isothermal" else [])
_ad = ["plate", "plate_in"] + ([] if D["cyl_top_thermal"] == "isothermal" else ["cyl_top"])
WALLS_ISO = [g for g in _iso if g in PID]
WALLS_AD = [g for g in _ad if g in PID]
SLIPS = [g for g in (["top", "side", "sym"] + (["runup"] if G["has_runup"] else [])) if g in PID]


# ---------------------------------------------------------------- config
def bcond(stage, inlet_profile=False):
    """stage: 'convert'|'slip'|'adiabatic'|'isothermal'
    convert = 最終形 (wall_dist 用)。slip = S0 用に全壁 slip。"""
    L = []
    prof = "{inletProfile: 1}" if inlet_profile else ""
    L.append("inlet:     {physID: %d, kind: inlet_uniformVelocity, outputHDFflg: 0, ints: %s, "
             "floats: {ro: %.8g, Ux: %.8g, Uy: 0.0, Uz: 0.0, Ps: %.8g, k: %.8g, omega: %.8g}}"
             % (PID["inlet"], prof, D["ro_inf"], D["U_inf"], D["P_inf"], D["k_inf"], D["omega_inf"]))
    L.append("outlet:    {physID: %d, kind: outlet_statPress, outputHDFflg: 0, ints: , "
             "floats: {Ps: %.8g, Pt: %.8g, Tt: %.8g}}"
             % (PID["outlet"], D["P_inf"], D["P_inf"], D["T_inf"]))
    for s in SLIPS:
        L.append("%-10s {physID: %d, kind: slip, outputHDFflg: 0, ints: , floats: }" % (s + ":", PID[s]))
    for w in WALLS_AD + WALLS_ISO:
        if stage == "slip":
            L.append("%-10s {physID: %d, kind: slip, outputHDFflg: 0, ints: , floats: }" % (w + ":", PID[w]))
        elif stage == "adiabatic" or (stage == "convert" and w in WALLS_AD) or \
                (stage in ("convert", "isothermal") and w in WALLS_AD):
            L.append("%-10s {physID: %d, kind: wall, outputHDFflg: 1, ints: , "
                     "floats: {Ux: 0.0, Uy: 0.0, Uz: 0.0}}" % (w + ":", PID[w]))
        else:                                  # convert / isothermal の等温壁
            L.append("%-10s {physID: %d, kind: wall_isothermal, outputHDFflg: 1, ints: , "
                     "floats: {Ux: 0.0, Uy: 0.0, Uz: 0.0, Ts: %.8g}}" % (w + ":", PID[w], D["wall_T"]))
    return "\n".join(L) + "\n"


def solver_cfg(nsteps, cfl, *, conv=1, lim=2, ninner=4, relax=0.7, outint=2000,
               model="sst", mesh_file="mesh.h5", gas=None):
    gas = (gas or D["gas"]).upper()
    kc = D["mu_inf"] * D["cp"] / D["prandtl_lam"]
    if gas == "TP":
        # semi-perfect 乾燥空気 1 擬似種 (plan §4.10)。`thermoHrefTemp` は必須
        # (絶対基準 h だと chi_eos が桁違いになる)。cp/gamma も parser が要求するので残す。
        phys = ('physProp: {thermalMethod: 2, species: ["MIXDRY"], speciesDBFile: "species_db.yaml", '
                'thermoHrefTemp: 298.15, viscMethod: 1, visc: %.8g, thermCond: %.8g, thermCondMethod: 1, '
                'prandtlLam: %s, cp: %s, gamma: %s}'
                % (D["mu_inf"], kc, D["prandtl_lam"], D["cp"], D["gamma"]))
    else:
        phys = ('physProp: {thermalMethod: 0, viscMethod: 1, visc: %.8g, thermCond: %.8g, '
                'thermCondMethod: 1, prandtlLam: %s, cp: %s, gamma: %s}'
                % (D["mu_inf"], kc, D["prandtl_lam"], D["cp"], D["gamma"]))
    turb = ('turbulence: {model: "none"}' if model == "none" else
            'turbulence: {model: "sst", scalarDiffusion: 1, dilatationCorrection: 2, katoLaunder: 0, '
            'wallTreatmentSST: 0, turbulentPrandtl: %.4g, kInit: %.8g, omegaInit: %.8g}'
            % (D["prandtl_turb"], D["k_inf"], D["omega_inf"]))
    return f"""mesh: {{discretization: "node", nodeWallDirichlet: 1, nodeInletCornerWall: 1, meshFileName: "{mesh_file}", valueFileName: "{mesh_file}"}}
gpu: 1
solver: "SLAU"
{phys}
time:
  unsteady: 0
  dualTime: 0
  last: {{nStepOuter: {nsteps}}}
  deltaT: {{control: 1, dt: 1e-9, cfl: {cfl}, cfl_pseudo: {cfl}, implicitRelax: {relax}, blockDPLUR: 1, lowMachPrecond: 0, dt_min: 1e-11, dt_max: 1.0, detectNaN: 1}}
  outStepStart: 0
  outStepInterval: {outint}
  timeIntegration: 11
  nStepInner: {ninner}
space: {{convMethod: {conv}, limiter: {lim}, pRef: {D['P_inf']:.8g}}}
{turb}
output: {{level: 1}}
initial: "uniform_p101325_u10"
"""


# ---------------------------------------------------------------- 取り込み
def cmd_convert(a):
    conv = HERE / "mesh" / "_conv"
    conv.mkdir(parents=True, exist_ok=True)
    # 変換は幾何と wall_dist だけなので EOS に依らない -> CPG で通す (species DB 不要)
    (conv / "solverConfig.yaml").write_text(solver_cfg(10, 0.5, conv=0, lim=0, ninner=5,
                                                       outint=10, mesh_file="m.h5", gas="CPG"))
    (conv / "bcondConfig.yaml").write_text(bcond("convert"))
    r = subprocess.run([str(BUILD / "convertGmshToForge"), str(Path(a.msh).resolve()), "m.h5"],
                       cwd=conv, env=ENV, capture_output=True, text=True)
    (conv / "convert.log").write_text(r.stdout + r.stderr)
    if not (conv / "m.h5").exists():
        print(r.stdout[-3000:], r.stderr[-2000:])
        raise SystemExit("convert failed")
    out = Path(a.out if os.path.isabs(a.out) else HERE / a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    shutil.move(str(conv / "m.h5"), str(out))
    for x in conv.glob("m.xmf"):
        x.unlink()
    for line in (conv / "convert.log").read_text().splitlines():
        if any(k in line for k in ("closure", "dual", "nBconds", "Number of", "bname")):
            print("   ", line)
    q = subprocess.run([sys.executable, str(TOOLS / "check_mesh_quality.py"), str(out)],
                       capture_output=True, text=True)
    (out.parent / ("quality_%s.txt" % out.stem)).write_text(q.stdout + q.stderr)
    print("\n".join(l for l in q.stdout.splitlines() if "VERDICT" in l or "max" in l.lower())[:2000])
    print("wrote", out)


# ---------------------------------------------------------------- IC
def _tp_energy(T, ro, u2):
    """TP (semi-perfect 1 擬似種) の roe。datum は thermoHrefTemp=298.15
    (h(298.15)=0 の基準)。e = h(T) - h(Tref) - R T。"""
    sys.path.insert(0, str(ROOT / "design"))
    from forge_design.gas.semiperfect import GasSemiPerfect
    g = GasSemiPerfect(D["dry_air_Y"], Tt=2000.0)
    href = float(g.h_mass(298.15))
    h = np.array([float(g.h_mass(t)) for t in np.atleast_1d(T)]) - href
    e = h - D["R_tp"] * np.atleast_1d(T)
    return ro * (e + 0.5 * u2)


def patch_ic(h5, inlet_csv, gas="CPG"):
    """外部流は入口 BL 分布を z で写す。キャビティ内 (z<0) は静止・壁温・P_inf。"""
    prof = np.genfromtxt(inlet_csv, names=True)
    zc = np.asarray(prof["z"], float)
    o = np.argsort(zc)
    zc = zc[o]
    get = lambda n: np.asarray(prof[n], float)[o]      # noqa: E731
    with h5py.File(h5, "r+") as f:
        c = np.array(f["MESH/COORD"]).reshape(-1, 3)
        n = f["/VALUE/ro"].shape[0]
        z = c[:, 2]
        ro = np.interp(z, zc, get("ro"))
        ux = np.interp(z, zc, get("Ux"))
        uz = np.interp(z, zc, get("Uz"))
        ps = np.interp(z, zc, get("Ps"))
        kk = np.interp(z, zc, get("k"))
        om = np.interp(z, zc, get("omega"))
        uy = np.zeros(n)
        incav = z < 0.0                                # キャビティ内は静止・壁温
        ro[incav] = D["P_inf"] / (D["R"] * D["wall_T"])
        ux[incav] = 0.0; uz[incav] = 0.0
        ps[incav] = D["P_inf"]
        kk[incav] = D["k_inf"] * 1e-3; om[incav] = D["omega_inf"]
        wd = np.array(f["/VALUE/wall_dist"]) if "/VALUE/wall_dist" in f else np.ones(n)
        wall = wd <= 0.0
        ux[wall] = 0.0; uy[wall] = 0.0; uz[wall] = 0.0
        u2 = ux ** 2 + uy ** 2 + uz ** 2
        if gas.upper() == "TP":
            # T を (ps, ro) から作り直して TP の e で roe を組む (CPG の ps/(γ-1) は datum が違う)
            Tfield = ps / np.maximum(ro * D["R_tp"], 1e-30)
            # 一意な温度だけ評価して展開 (NASA-9 の評価は node 数ぶん回すと遅い)
            Tq = np.round(Tfield, 2)
            uniq, inv = np.unique(Tq, return_inverse=True)
            sys.path.insert(0, str(ROOT / "design"))
            from forge_design.gas.semiperfect import GasSemiPerfect
            g = GasSemiPerfect(D["dry_air_Y"], Tt=2000.0)
            href = float(g.h_mass(298.15))
            e_u = np.array([float(g.h_mass(t)) - href - D["R_tp"] * t for t in uniq])
            roe = ro * (e_u[inv] + 0.5 * u2)
            print("  IC(TP): T %.1f..%.1f K, 一意温度 %d 点" % (Tfield.min(), Tfield.max(), len(uniq)))
        else:
            roe = ps / (D["gamma"] - 1.0) + 0.5 * ro * u2
        f["/VALUE/ro"][:] = ro.astype(np.float32)
        f["/VALUE/roUx"][:] = (ro * ux).astype(np.float32)
        f["/VALUE/roUy"][:] = (ro * uy).astype(np.float32)
        f["/VALUE/roUz"][:] = (ro * uz).astype(np.float32)
        f["/VALUE/roe"][:] = roe.astype(np.float32)
        for name, v in (("roK", ro * kk), ("roOmega", ro * om)):
            ds = "/VALUE/" + name
            if ds in f:
                f[ds][:] = v.astype(np.float32)
            else:
                f.create_dataset(ds, data=v.astype(np.float32))
        print("IC: n=%d  cavity nodes=%d  wall nodes=%d" % (n, int(incav.sum()), int(wall.sum())))


def index_copy(res, mesh):
    """同一メッシュの段間引き継ぎ: res の VALUE を mesh.h5 に index コピー (wall_dist は残す)。"""
    with h5py.File(res, "r") as s, h5py.File(mesh, "r+") as d:
        for k in s["VALUE"]:
            if k == "wall_dist":
                continue
            arr = np.array(s["VALUE/" + k])
            ds = "VALUE/" + k
            if ds in d:
                if d[ds].shape != arr.shape:
                    continue
                d[ds][...] = arr.astype(d[ds].dtype)
            else:
                d.create_dataset(ds, data=arr)


def run_forge(rd):
    r = subprocess.run([str(TOOLS / "run_case.sh"), str(rd)], env=ENV, capture_output=True, text=True)
    (rd / "run_case_stdout.log").write_text(r.stdout + r.stderr)
    return r.returncode


def stage(rd, tag, cfgtext, bctext, nsteps, keep=False):
    (rd / "solverConfig.yaml").write_text(cfgtext)
    (rd / "bcondConfig.yaml").write_text(bctext)
    rc = run_forge(rd)
    res = sorted(rd.glob("res_[0-9]*.h5"), key=lambda f: int(f.stem.split("_")[1]))
    # 完走判定は residual_history の最終 step で行う (出力間隔と nStepOuter が割り切れないと
    # 最後の res が出ないため)。段間引き継ぎに使う最終場が要るので outint=nsteps にしてある。
    last_step = -1
    rh = rd / "residual_history.csv"
    if rh.exists():
        tail = rh.read_text().strip().splitlines()
        if len(tail) > 1:
            last_step = int(float(tail[-1].split(",")[0]))
    ok = res and int(res[-1].stem.split("_")[1]) >= nsteps and last_step >= nsteps - 2
    print("  [%s] rc=%d last_res=%s last_step=%d" %
          (tag, rc, res[-1].name if res else "-", last_step), flush=True)
    if not ok:
        raise SystemExit("stage %s failed (rc=%d)" % (tag, rc))
    index_copy(res[-1], rd / "mesh.h5")
    if not keep:
        for f in rd.glob("res_*"):
            f.unlink()
    for nm in ("residual_history.csv", "residual_history.png", "CONVERGENCE_VERDICT.txt", "forge_run.log"):
        if (rd / nm).exists():
            shutil.move(str(rd / nm), str(rd / ("%s_%s%s" % (Path(nm).stem, tag, Path(nm).suffix))))


def cmd_run(a):
    rd = HERE / a.run
    if rd.exists():
        raise SystemExit("%s exists" % rd)
    rd.mkdir()
    shutil.copy(HERE / a.mesh if not os.path.isabs(a.mesh) else a.mesh, rd / "mesh.h5")
    shutil.copy(a.inlet_csv, rd / ("inlet_profile_%d.csv" % PID["inlet"]))
    gas = (a.gas or D["gas"]).upper()
    if gas == "TP":
        sdb = HERE / "species_db.yaml"
        if not sdb.exists():
            raise SystemExit("species_db.yaml が無い (setup で生成: mixture_pseudo_species)")
        shutil.copy(sdb, rd / "species_db.yaml")
    print("  gas =", gas)
    (rd / "GEN_ARGS").write_text(" ".join(sys.argv[1:]) + "\n")
    (rd / "probe.yaml").write_text("outStepInterval: 100\noutStepStart: 0\npoints:\nsurfaces:\n")
    patch_ic(rd / "mesh.h5", a.inlet_csv, gas=gas)
    if a.dry:
        (rd / "solverConfig.yaml").write_text(solver_cfg(a.main_steps, a.cfl, gas=gas))
        (rd / "bcondConfig.yaml").write_text(bcond("isothermal", inlet_profile=True))
        return
    ip = dict(inlet_profile=True)
    # S0: 全壁 slip・層流・1 次 (キャビティ内圧の平衡化)
    stage(rd, "S0_slip", solver_cfg(2000, 0.5, conv=0, lim=0, ninner=10, outint=2000, model="none", gas=gas),
          bcond("slip", **ip), 2000)
    # S1: no-slip 断熱 (層流)
    stage(rd, "S1_lam", solver_cfg(2000, 0.3, conv=0, lim=0, ninner=10, outint=2000, model="none", gas=gas),
          bcond("adiabatic", **ip), 2000)
    # S2: キャビティ等温壁 (層流)
    stage(rd, "S2_iso", solver_cfg(2000, 0.5, conv=0, lim=0, ninner=10, outint=2000, model="none", gas=gas),
          bcond("isothermal", **ip), 2000)
    # S3/S4: SST soft -> mid (1 次)
    stage(rd, "S3_sst_soft", solver_cfg(3000, 0.3, conv=0, lim=0, ninner=10, outint=3000, gas=gas), bcond("isothermal", **ip), 3000)
    stage(rd, "S4_sst_mid", solver_cfg(3000, 1.0, conv=0, lim=0, ninner=10, outint=3000, gas=gas), bcond("isothermal", **ip), 3000)
    # S5: 2 次ランプ
    for i, cv in enumerate([float(v) for v in a.ramp.split(",") if v]):
        stage(rd, "S5_ramp%d_cfl%g" % (i, cv), solver_cfg(2000, cv, outint=2000, gas=gas), bcond("isothermal", **ip), 2000)
    # S6: 本段
    (rd / "solverConfig.yaml").write_text(solver_cfg(a.main_steps, a.cfl, outint=a.out_int, gas=gas))
    (rd / "bcondConfig.yaml").write_text(bcond("isothermal", **ip))
    rc = run_forge(rd)
    print("main rc", rc)
    print((rd / "CONVERGENCE_VERDICT.txt").read_text()[-900:])


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("convert")
    c.add_argument("--msh", required=True)
    c.add_argument("--out", required=True)
    r = sub.add_parser("run")
    r.add_argument("--run", required=True)
    r.add_argument("--mesh", required=True)
    r.add_argument("--inlet-csv", required=True)
    r.add_argument("--main-steps", type=int, default=20000)
    r.add_argument("--cfl", type=float, default=2.0)
    r.add_argument("--out-int", type=int, default=2000)
    r.add_argument("--ramp", default="0.5,1,2")
    r.add_argument("--gas", default=None, choices=["CPG", "TP"], help="既定は case.json の gas")
    r.add_argument("--dry", action="store_true")
    a = ap.parse_args()
    (cmd_convert if a.cmd == "convert" else cmd_run)(a)


if __name__ == "__main__":
    main()
