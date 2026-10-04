#!/usr/bin/env python3
"""V2d TP での境界の独立参照 (plan boundary-node-farfield-characteristic §6 V2d) の forge run を作る (AWS で実行)。
V2a と同じ薄板チャネル (Δx 5 mm、y・z 各 3 層 slip、左右端 farfield)、dual-time。自由流 = SERN 外気 (220 K、2851 Pa、Y_EXH 0)。
  --gas cpg | tp1 (単成分 AIR、NASA-9) | tp2 (EXH/AIR)
  --test contact : V2d-1。P・u 一様 (P∞、u = ±0.5 c∞)、中心 0.5 m (長領域も同じ x) にガウス (半値全幅 0.1 m) の T 600 K・Y_EXH 0.13 の塊。
                   塊は右端 (u > 0) から出ていく。評価点 = 右端から 5 セル内側 (x = 0.975 m)。場を 10 step ごとに書く (T・Y の時系列)。
  --test acoustic: V2d-2。内部 = 一様の T 600 K・Y_EXH 0.13 (P∞、u = 0.3 c_i)、中心 0.4 m に右向きの圧力パルス (振幅 1e-3 P∞、FWHM 40 mm)。
                   farfield の自由流は外気 (内外で γ・ρc が違う)。評価点 x = 0.8 m。--nopulse でパルスなし対照。
リミッタ基準値は自由流で固定 (V2a と同じ理由)。
  python3 setup_v2d.py RUN --gas G --test T --len 1|3 [--dir 1|-1] [--nopulse] [--cfl-ac 0.113] [--nsub 20]
"""
import argparse, math, os, shutil, subprocess, sys
import h5py, numpy as np, yaml

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from setup_v1 import SPECIES_SRC, RU, nasa9_h_mass  # noqa: E402

P0, T0 = 2851.0, 220.0


class Gas:
    def __init__(self, kind, db):
        self.kind, self.db = kind, db

    def R(self, Y):
        if self.kind == "cpg":
            return 1004.5 - 1004.5 / 1.4
        if self.kind == "tp1":
            return RU / self.db["AIR"]["MW"]
        return (Y / self.db["EXH"]["MW"] + (1 - Y) / self.db["AIR"]["MW"]) * RU

    def h(self, Y, T):   # 顕熱 (datum 298.15 K、setup_v1 と同じ規約)
        if self.kind == "cpg":
            return 1004.5 * T
        sp = {"AIR": 1.0} if self.kind == "tp1" else {"EXH": Y, "AIR": 1 - Y}
        return sum(w * (nasa9_h_mass(self.db[s], T) - nasa9_h_mass(self.db[s], 298.15)) for s, w in sp.items())

    def e(self, Y, T):
        return self.h(Y, T) - self.R(Y) * T if self.kind != "cpg" else 1004.5 / 1.4 * T

    def c(self, Y, T):
        if self.kind == "cpg":
            return math.sqrt(1.4 * self.R(Y) * T)
        cp = (self.h(Y, T + 0.01) - self.h(Y, T - 0.01)) / 0.02
        R = self.R(Y)
        return math.sqrt(cp / (cp - R) * R * T)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run"); ap.add_argument("--gas", choices=("cpg", "tp1", "tp2"), required=True)
    ap.add_argument("--test", choices=("contact", "acoustic"), required=True); ap.add_argument("--len", type=float, required=True)
    ap.add_argument("--dir", type=int, default=1); ap.add_argument("--nopulse", action="store_true")
    ap.add_argument("--left-hot", action="store_true", help="診断: 左端 farfield の自由流を内部と同じ高温状態にする (左からの接触面を作らない)")
    ap.add_argument("--cfl-ac", type=float, default=0.113); ap.add_argument("--nsub", type=int, default=20)
    ap.add_argument("--bin", default=os.path.expanduser("~/forge-pgrad-new/solver_density_cuda/build-ff"))
    a = ap.parse_args()
    run = a.run
    os.makedirs(run)
    for f in ("species_db.yaml", "species_meta.yaml"):
        shutil.copy(os.path.join(SPECIES_SRC, f), run)
    db = yaml.safe_load(open(os.path.join(run, "species_db.yaml")))
    g = Gas(a.gas, db)
    dx, W, nyz = 0.005, 0.015, 3
    nx = int(round(a.len / dx))
    subprocess.run([sys.executable, os.path.join(HERE, "make_box_msh.py"), os.path.join(run, "chan.msh"), str(nx), str(nyz), str(nyz),
                    repr(a.len), repr(W), repr(W)], check=True)
    ro0 = P0 / (g.R(0.0) * T0); c0 = g.c(0.0, T0)
    Ti, Yi = 600.0, (0.13 if a.gas == "tp2" else 0.0)
    ci = g.c(Yi, Ti)
    if a.test == "contact":
        U = a.dir * 0.5 * c0; umax = abs(U) + max(c0, ci); tend = (0.475 + 0.15) / abs(U); probes = {"eval": 0.975 if a.dir > 0 else 0.025}
    else:
        # 終了時刻 = 反射到達窓の末尾 (t_ref + 3σ/(c−u)) + 余裕 0.02 m/(c−u)。旧式 (0.4/(c+u) + 0.35/(c−u)) は窓を覆わなかった (codex diagnose 3 回目)
        U = 0.3 * ci; umax = U + ci; sig_a = 0.04 / (2.0 * math.sqrt(2.0 * math.log(2.0)))
        tend = 0.6 / (ci + U) + 0.2 / (ci - U) + (3.0 * sig_a + 0.02) / (ci - U); probes = {"eval": 0.8, "mid": 0.4}
    dt = a.cfl_ac * dx / umax
    steps = int(math.ceil(tend / dt))
    if a.gas == "cpg":
        phys = "{thermalMethod: 0, viscMethod: 0, visc: 0.0, thermCond: 0.0, cp: 1004.5, gamma: 1.4}"
    else:
        spl = '["AIR"]' if a.gas == "tp1" else '["EXH", "AIR"]'
        phys = f'{{thermalMethod: 2, viscMethod: 0, visc: 0.0, thermCond: 0.0, cp: 1004.5, gamma: 1.4, pMin: 20.0, species: {spl}, speciesDBFile: "species_db.yaml", thermoHrefTemp: 298.15}}'
    outint = 10 if a.test == "contact" else steps
    open(os.path.join(run, "solverConfig.yaml"), "w").write(f"""mesh: {{discretization: "node", nodeWallDirichlet: 1, meshFileName: "chan.h5", valueFileName: "chan.h5"}}
gpu: 1
solver: "SLAU"
physProp: {phys}
time:
  unsteady: 1
  dualTime: 1
  nSubIterDualTime: {a.nsub}
  last: {{nStepOuter: {steps}}}
  deltaT: {{control: 0, dt: {dt!r}, cfl: 1.0, cfl_pseudo: 12.0, dt_min: 1e-12, dt_max: 1.0, blockDPLUR: 1, lowMachPrecond: 0, detectNaN: 1}}
  outStepStart: 0
  outStepInterval: {outint}
  timeIntegration: 11
  nStepInner: 5
space: {{convMethod: 1, limiter: 2, pRef: {P0}, limiterRefLength: 1.0, limiterRoRef: {ro0!r}, limiterPRef: {P0!r}, limiterARef: {c0!r}}}
turbulence: {{model: "none"}}
initial: "uniform_p101325_u10"
""")
    fl = f"ro: {ro0!r}, Ux: {U!r}, Uy: 0.0, Uz: 0.0, Ps: {P0!r}" + (", Y0: 0.0, Y1: 1.0" if a.gas == "tp2" else "")
    kinds = {1: ("xmin", "farfield"), 2: ("xmax", "farfield"), 3: ("ymin", "slip"), 4: ("ymax", "slip"), 5: ("zmin", "slip"), 6: ("zmax", "slip")}
    flh = fl
    if a.left_hot:
        roh = P0 / (g.R(Yi) * Ti)
        flh = f"ro: {roh!r}, Ux: {U!r}, Uy: 0.0, Uz: 0.0, Ps: {P0!r}" + (f", Y0: {Yi!r}, Y1: {1 - Yi!r}" if a.gas == "tp2" else "")
    open(os.path.join(run, "bcondConfig.yaml"), "w").write("".join(
        f"{nm}: {{physID: {p}, kind: {k}, outputHDFflg: 0, ints: , floats: {{{(flh if p == 1 else fl) if k == 'farfield' else ''}}}}}\n" for p, (nm, k) in kinds.items()))
    yp = W / 3.0 + 0.001
    open(os.path.join(run, "probe.yaml"), "w").write("outStepInterval: 1\noutStepStart: 0\npoints:\n" + "".join(
        f"  {k}: {{x: {v + 0.001}, y: {yp}, z: {yp}}}\n" for k, v in probes.items()) + "surfaces:\n")
    subprocess.run([os.path.join(a.bin, "convertGmshToForge"), "chan.msh", "chan.h5"], cwd=run,
                   stdout=open(os.path.join(run, "convert.log"), "w"), stderr=subprocess.STDOUT)
    if not os.path.exists(os.path.join(run, "chan.h5")) or "writeInputH5: wrote /VIZMESH" not in open(os.path.join(run, "convert.log")).read():
        raise SystemExit(f"変換に失敗: {run}/convert.log")
    with h5py.File(os.path.join(run, "chan.h5"), "r+") as f:
        V = f["VALUE"]; n = V["ro"].shape[0]
        x = np.array(f["MESH/COORD"]).reshape(-1, 3)[:, 0].astype(float)
        if a.test == "contact":
            sig = 0.1 / (2.0 * math.sqrt(2.0 * math.log(2.0)))
            w = np.exp(-0.5 * ((x - 0.5) / sig) ** 2)
            T = T0 + (Ti - T0) * w; Y = Yi * w; P = np.full(n, P0); u = np.full(n, U)
        else:
            sig = 0.04 / (2.0 * math.sqrt(2.0 * math.log(2.0)))
            roi = P0 / (g.R(Yi) * Ti)
            pp = 0.0 if a.nopulse else 1e-3 * P0 * np.exp(-0.5 * ((x - 0.4) / sig) ** 2)
            P = P0 + pp * np.ones(n); u = U + pp / (roi * ci); Y = np.full(n, Yi)
            # 等エントロピー擾乱: T を P から (小振幅、線形)
            gi = ci * ci * roi / P0
            T = Ti * (P / P0) ** ((gi - 1.0) / gi)
        R = np.array([g.R(y) for y in Y]); ro = P / (R * T)
        e = np.array([g.e(y, t) for y, t in zip(Y, T)])
        vals = {"ro": ro, "roUx": ro * u, "roUy": 0.0 * ro, "roUz": 0.0 * ro, "roe": ro * (e + 0.5 * u * u)}
        if a.gas == "tp2":
            vals.update({"roY0": ro * Y, "roY1": ro * (1 - Y)})
        for kk, v in vals.items():
            if kk in V:
                V[kk][...] = v.astype(V[kk].dtype)
            else:
                V.create_dataset(kk, data=v.astype(V["ro"].dtype))
    open(os.path.join(run, "IC_FROM.txt"), "w").write(
        f"V2d {a.test} gas {a.gas}: left_hot {a.left_hot}、L {a.len} m、U {U:.5g} (dir {a.dir})、nopulse {a.nopulse}、dt {dt:.4e}、steps {steps}、"
        f"自由流 ρ {ro0:.6g} c {c0:.5g}、内部 T {Ti} Y {Yi} c {ci:.5g}\n")
    print(f"prepared {run}: steps {steps}, dt {dt:.3e}")


if __name__ == "__main__":
    main()
