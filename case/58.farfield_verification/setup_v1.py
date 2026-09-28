#!/usr/bin/env python3
"""V1 自由流の保持 (plan boundary-node-farfield-characteristic §6 V1) の run を作る (AWS で実行)。
直方体 hex (1 × 0.5 × 0.5 m、24 × 12 × 12)、6 面すべて farfield、初期値 = 自由流 (境界の floats と同じ状態)。
  python3 setup_v1.py RUN_DIR VARIANT [--steps N] [--bin DIR]
VARIANT: a = CPG M 0.5 / b = TP 外気 M 6 (SERN の EXH/AIR、Y_EXH 0) / c = b を xy 面内で 30° 傾ける / d = b + SST
TP の species_db.yaml / species_meta.yaml は SERN の run (SPECIES_SRC) から複製する。
"""
import math, os, shutil, subprocess, sys
import h5py, numpy as np, yaml

HERE = os.path.dirname(os.path.abspath(__file__))
SPECIES_SRC = os.path.join(HERE, "..", "46.sern_design", "run_0986_r4d_g3_base")
RU = 8.31446261815324


def nasa9_h_mass(e, T):
    a = e["nasa9_low"] if T < e["Tmid"] else e["nasa9_high"]
    hRT = -a[0] / T ** 2 + a[1] * math.log(T) / T + a[2] + a[3] * T / 2 + a[4] * T ** 2 / 3 + a[5] * T ** 3 / 4 + a[6] * T ** 4 / 5 + a[7] / T
    return hRT * RU * T / e["MW"]


def main():
    run, var = sys.argv[1], sys.argv[2]
    steps = int(sys.argv[sys.argv.index("--steps") + 1]) if "--steps" in sys.argv else 2000
    bdir = sys.argv[sys.argv.index("--bin") + 1] if "--bin" in sys.argv else os.path.expanduser("~/forge-pgrad-new/solver_density_cuda/build-ff")
    os.makedirs(run)
    shutil.copy(os.path.join(SPECIES_SRC, "probe.yaml"), run)
    subprocess.run([sys.executable, os.path.join(HERE, "make_box_msh.py"), os.path.join(run, "box.msh"), "24", "12", "12", "1.0", "0.5", "0.5"], check=True)
    tp = var in ("b", "c", "d")
    sst = var == "d"
    if tp:
        for f in ("species_db.yaml", "species_meta.yaml"):
            shutil.copy(os.path.join(SPECIES_SRC, f), run)
        db = yaml.safe_load(open(os.path.join(run, "species_db.yaml")))
        P, T = 2851.0, 220.0
        air = db["AIR"]; R = RU / air["MW"]; ro = P / (R * T)
        # γ(T) (frozen) で音速
        cpv = (nasa9_h_mass(air, T + 0.01) - nasa9_h_mass(air, T - 0.01)) / 0.02
        a = math.sqrt(cpv / (cpv - R) * R * T); M = 6.0
        e = (nasa9_h_mass(air, T) - nasa9_h_mass(air, 298.15)) - R * T
        phys = ('{thermalMethod: 2, viscMethod: 0, visc: 0.0, thermCond: 0.0, cp: 1004.5, gamma: 1.4, pMin: 20.0, '
                'species: ["EXH", "AIR"], speciesDBFile: "species_db.yaml", thermoHrefTemp: 298.15}')
    else:
        P, T, gam, cp = 101325.0, 300.0, 1.4, 1004.5
        R = cp - cp / gam; ro = P / (R * T); a = math.sqrt(gam * R * T); M = 0.5
        e = P / ((gam - 1) * ro)
        phys = "{thermalMethod: 0, viscMethod: 0, visc: 0.0, thermCond: 0.0, cp: 1004.5, gamma: 1.4}"
    ang = math.radians(30.0) if var == "c" else 0.0
    U = M * a; ux, uy, uz = U * math.cos(ang), U * math.sin(ang), 0.0
    k, om = (479.653, 119844.0) if sst else (0.0, 0.0)
    turb = ('{model: "sst", scalarDiffusion: 1, dilatationCorrection: 0, katoLaunder: 0, wallTreatmentSST: 0}' if sst else '{model: "none"}')
    open(os.path.join(run, "solverConfig.yaml"), "w").write(f"""mesh: {{discretization: "node", nodeWallDirichlet: 1, meshFileName: "box.h5", valueFileName: "box.h5"}}
gpu: 1
solver: "SLAU"
physProp: {phys}
time:
  unsteady: 0
  dualTime: 0
  last: {{nStepOuter: {steps}}}
  deltaT: {{control: 1, dt: 1e-8, cfl: 2.0, cfl_pseudo: 2.0, dt_min: 1e-12, dt_max: 1.0, blockDPLUR: 1, lowMachPrecond: 0, detectNaN: 1}}
  outStepStart: 0
  outStepInterval: {max(steps // 4, 1)}
  timeIntegration: 11
  nStepInner: 5
space: {{convMethod: 1, limiter: 2, pRef: {P}}}
turbulence: {turb}
initial: "uniform_p101325_u10"
""")
    fl = f"ro: {ro!r}, Ux: {ux!r}, Uy: {uy!r}, Uz: {uz!r}, Ps: {P!r}"
    if sst:
        fl += f", k: {k!r}, omega: {om!r}"
    if tp:
        fl += ", Y0: 0.0, Y1: 1.0"
    names = {1: "xmin", 2: "xmax", 3: "ymin", 4: "ymax", 5: "zmin", 6: "zmax"}
    open(os.path.join(run, "bcondConfig.yaml"), "w").write("".join(
        f"{names[p]}: {{physID: {p}, kind: farfield, outputHDFflg: 0, ints: , floats: {{{fl}}}}}\n" for p in range(1, 7)))
    # 変換器は書き終えた後に GPUassert で非ゼロ終了することがある (既知)。h5 ができたかで判定する
    subprocess.run([os.path.join(bdir, "convertGmshToForge"), "box.msh", "box.h5"], cwd=run,
                   stdout=open(os.path.join(run, "convert.log"), "w"), stderr=subprocess.STDOUT)
    if not os.path.exists(os.path.join(run, "box.h5")) or "writeInputH5: wrote /VIZMESH" not in open(os.path.join(run, "convert.log")).read():
        raise SystemExit(f"変換に失敗: {run}/convert.log")
    with h5py.File(os.path.join(run, "box.h5"), "r+") as f:
        V = f["VALUE"]; n = V["ro"].shape[0]
        vals = {"ro": ro, "roUx": ro * ux, "roUy": ro * uy, "roUz": ro * uz, "roe": ro * (e + 0.5 * U * U)}
        if sst:
            vals.update({"roK": ro * k, "roOmega": ro * om})
        if tp:
            vals.update({"roY0": 0.0, "roY1": ro})
        for kk, v in vals.items():
            if kk in V:
                V[kk][...] = np.full(n, v, dtype=V[kk].dtype)
            else:
                V.create_dataset(kk, data=np.full(n, v, dtype=V["ro"].dtype))
    open(os.path.join(run, "IC_FROM.txt"), "w").write(
        f"V1 自由流の保持 variant {var}: 一様場 ρ {ro:.6g}、u ({ux:.6g}, {uy:.6g}, {uz:.6g})、P {P}、T {T}、M {M}"
        + (f"、k {k}、ω {om}" if sst else "") + (" (TP, Y_EXH 0)" if tp else " (CPG)") + "\n")
    print(f"prepared {run}: M {M} U {U:.4g} ρ {ro:.6g}")


if __name__ == "__main__":
    main()
