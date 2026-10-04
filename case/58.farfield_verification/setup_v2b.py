#!/usr/bin/env python3
"""V2b 保存収支・V2f 局所逆流 (plan boundary-node-farfield-characteristic §6 V2b/V2f、配置は 2026-09-29 確定) の run を作る (AWS で実行)。
V1 と同じ箱 (1 × 0.5 × 0.5 m、24 × 12 × 12)、全 6 面 farfield、TP (EXH/AIR) + SST、自由流 = SERN 外気 (220 K、2851 Pa、Y_EXH 0)。
  V2b: python3 setup_v2b.py RUN --dir x|oblique [--ek 0|1] [--sfr 0|2] [--steps 3000]
       自由流 M 0.5。x = +x 方向 (流出配置)、oblique = (cos30 cos20, sin30 cos20, sin20) (流入配置)。
       初期値 = 内部全域 P∞・u∞ のまま T 600 K・Y_EXH 0.13・k・ω 10 倍。
  V2f: python3 setup_v2b.py RUN --v2f [--steps 3000]
       自由流 +x M 2 (xmax で Q_n = 2a∞)。初期値 = 自由流、ただし x ≥ 0.8 m の帯だけ U_x = −0.95 a・Y_EXH 0.13。
代表量 (判定の規格化) を v2b_meta.json に書く。
"""
import argparse, json, math, os, shutil, subprocess, sys
import h5py, numpy as np, yaml

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from setup_v1 import SPECIES_SRC, RU, nasa9_h_mass  # noqa: E402


def gas(db, Y0, T):
    """(R, e, c) : 混合の気体定数・内部エネルギー (datum 298.15 K の顕熱、setup_v1 と同じ規約)・凍結音速"""
    Ys = {"EXH": Y0, "AIR": 1.0 - Y0}
    R = sum(Y / db[s]["MW"] for s, Y in Ys.items()) * RU
    h = lambda T_: sum(Y * (nasa9_h_mass(db[s], T_) - nasa9_h_mass(db[s], 298.15)) for s, Y in Ys.items())
    cpv = (h(T + 0.01) - h(T - 0.01)) / 0.02
    return R, h(T) - R * T, math.sqrt(cpv / (cpv - R) * R * T)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run"); ap.add_argument("--dir", choices=("x", "oblique"), default="x")
    ap.add_argument("--ek", type=int, default=0); ap.add_argument("--sfr", type=int, default=0)
    ap.add_argument("--v2f", action="store_true"); ap.add_argument("--steps", type=int, default=3000)
    ap.add_argument("--cfl", type=float, default=2.0)
    ap.add_argument("--bin", default=os.path.expanduser("~/forge-pgrad-new/solver_density_cuda/build-ff"))
    a = ap.parse_args()
    run = a.run
    os.makedirs(run)
    for f in ("probe.yaml", "species_db.yaml", "species_meta.yaml"):
        shutil.copy(os.path.join(SPECIES_SRC, f), run)
    subprocess.run([sys.executable, os.path.join(HERE, "make_box_msh.py"), os.path.join(run, "box.msh"), "24", "12", "12", "1.0", "0.5", "0.5"], check=True)
    db = yaml.safe_load(open(os.path.join(run, "species_db.yaml")))
    P, T = 2851.0, 220.0
    R, e, c = gas(db, 0.0, T); ro = P / (R * T)
    M = 2.0 if a.v2f else 0.5
    if a.dir == "oblique" and not a.v2f:
        t30, t20 = math.radians(30.0), math.radians(20.0)
        d = np.array([math.cos(t30) * math.cos(t20), math.sin(t30) * math.cos(t20), math.sin(t20)])
    else:
        d = np.array([1.0, 0.0, 0.0])
    U = M * c; u = [float(v) for v in U * d]
    k, om = 479.653, 119844.0
    sfr = "" if a.sfr == 0 else f", speciesFaceReconstruction: {a.sfr}, speciesImplicitCoupling: 1"
    open(os.path.join(run, "solverConfig.yaml"), "w").write(f"""mesh: {{discretization: "node", nodeWallDirichlet: 1, meshFileName: "box.h5", valueFileName: "box.h5"}}
gpu: 1
solver: "SLAU"
physProp: {{thermalMethod: 2, viscMethod: 0, visc: 0.0, thermCond: 0.0, cp: 1004.5, gamma: 1.4, pMin: 20.0, species: ["EXH", "AIR"], speciesDBFile: "species_db.yaml", thermoHrefTemp: 298.15}}
time:
  unsteady: 0
  dualTime: 0
  last: {{nStepOuter: {a.steps}}}
  deltaT: {{control: 1, dt: 1e-8, cfl: {a.cfl}, cfl_pseudo: {a.cfl}, dt_min: 1e-12, dt_max: 1.0, blockDPLUR: 1, lowMachPrecond: 0, detectNaN: 1{sfr}}}
  outStepStart: 0
  outStepInterval: {max(a.steps // 6, 1)}
  timeIntegration: 11
  nStepInner: 5
space: {{convMethod: 1, limiter: 2, pRef: {P}}}
turbulence: {{model: "sst", scalarDiffusion: 1, dilatationCorrection: 0, katoLaunder: 0, wallTreatmentSST: 0, sstEnergyIncludesK: {a.ek}}}
initial: "uniform_p101325_u10"
""")
    fl = f"ro: {ro!r}, Ux: {u[0]!r}, Uy: {u[1]!r}, Uz: {u[2]!r}, Ps: {P!r}, k: {k!r}, omega: {om!r}, Y0: 0.0, Y1: 1.0"
    names = {1: "xmin", 2: "xmax", 3: "ymin", 4: "ymax", 5: "zmin", 6: "zmax"}
    open(os.path.join(run, "bcondConfig.yaml"), "w").write("".join(
        f"{names[p]}: {{physID: {p}, kind: farfield, outputHDFflg: 0, ints: , floats: {{{fl}}}}}\n" for p in range(1, 7)))
    subprocess.run([os.path.join(a.bin, "convertGmshToForge"), "box.msh", "box.h5"], cwd=run,
                   stdout=open(os.path.join(run, "convert.log"), "w"), stderr=subprocess.STDOUT)
    if not os.path.exists(os.path.join(run, "box.h5")) or "writeInputH5: wrote /VIZMESH" not in open(os.path.join(run, "convert.log")).read():
        raise SystemExit(f"変換に失敗: {run}/convert.log")
    ekf = 1.0 if a.ek else 0.0
    with h5py.File(os.path.join(run, "box.h5"), "r+") as f:
        V = f["VALUE"]; n = V["ro"].shape[0]
        x = np.array(f["MESH/COORD"]).reshape(-1, 3)[:, 0]
        if a.v2f:
            band = x >= 0.8 - 1e-9
            Ri, ei, ci = gas(db, 0.13, T)
            Yn = np.where(band, 0.13, 0.0)
            rn = np.where(band, P / (Ri * T), ro)
            en = np.where(band, ei, e)
            ux = np.where(band, -0.95 * ci, u[0])
            vel = np.stack([ux, np.zeros(n), np.zeros(n)], 1)
            kn, omn = np.full(n, k), np.full(n, om)
        else:
            Ri, ei, ci = gas(db, 0.13, 600.0)
            Yn = np.full(n, 0.13); rn = np.full(n, P / (Ri * 600.0)); en = np.full(n, ei)
            vel = np.tile(u, (n, 1)); kn, omn = np.full(n, 10 * k), np.full(n, 10 * om)
        ke = 0.5 * np.sum(vel ** 2, 1)
        vals = {"ro": rn, "roUx": rn * vel[:, 0], "roUy": rn * vel[:, 1], "roUz": rn * vel[:, 2],
                "roe": rn * (en + ke + ekf * kn), "roK": rn * kn, "roOmega": rn * omn, "roY0": rn * Yn, "roY1": rn * (1.0 - Yn)}
        for kk, v in vals.items():
            if kk in V:
                V[kk][...] = v.astype(V[kk].dtype)
            else:
                V.create_dataset(kk, data=v.astype(V["ro"].dtype))
    Hinf = e + P / ro + 0.5 * U * U + ekf * k
    meta = {"ro": ro, "U": U, "u": list(u), "P": P, "T": T, "e": e, "H": Hinf, "k": k, "omega": om, "M": M, "ek": a.ek, "sfr": a.sfr,
            "dir": "v2f" if a.v2f else a.dir}
    json.dump(meta, open(os.path.join(run, "v2b_meta.json"), "w"), indent=1)
    open(os.path.join(run, "IC_FROM.txt"), "w").write(
        ("V2f: 自由流 +x M 2、x ≥ 0.8 m の帯を U_x −0.95a・Y_EXH 0.13 に" if a.v2f else
         f"V2b: 自由流 M 0.5 方向 {a.dir}、内部全域 T 600 K・Y_EXH 0.13・k・ω 10 倍 (P∞・u∞ のまま)")
        + f"、sstEnergyIncludesK {a.ek}、speciesFaceReconstruction {a.sfr}\n")
    print(f"prepared {run}: M {M} U {U:.4g} ρ∞ {ro:.6g} c∞ {c:.5g}")


if __name__ == "__main__":
    main()
