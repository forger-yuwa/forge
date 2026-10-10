"""一様な静止場 (plan architecture-float-state-double-geometry §6.3 の 3: 軸対称の closure の非劣化) を run の場のファイルに書く。
等温壁 (300 K) の節点から組成・気体定数 R = P/(ρT)・内部エネルギー e = ρE/ρ (u = 0) を取り、全節点を P0・T0・u = 0 にする。
入口・出口の境界の値も静止場と矛盾しない値に書き換える (入口 Pt = P0・Tt = T0、出口 Ps = Pt = P0)。
usage: python3 still_field.py <run_dir> [--P0 1.0e5]"""
import argparse, re, h5py, numpy as np
from pathlib import Path
ap = argparse.ArgumentParser(); ap.add_argument("run"); ap.add_argument("--P0", type=float, default=1.0e5); a = ap.parse_args()
run = Path(a.run); f = run / "nozzle.h5"
with h5py.File(f, "r+") as h:
    V = h["VALUE"]; T = np.asarray(V["T"][:]); P = np.asarray(V["P"][:]); ro = np.asarray(V["ro"][:]); roe = np.asarray(V["roe"][:])
    U2 = np.asarray(V["Ux"][:]) ** 2 + np.asarray(V["Uy"][:]) ** 2
    w = np.flatnonzero((np.abs(T - 300.0) < 1e-6) & (U2 < 1e-20))
    if w.size == 0: raise SystemExit("T = 300 K・u = 0 の節点 (等温壁) が無い")
    i = w[0]; T0 = 300.0; R = P[i] / (ro[i] * T[i]); e0 = roe[i] / ro[i]
    rho0 = a.P0 / (R * T0); n = ro.size
    def put(k, v):
        if k in V: V[k][...] = np.full(n, v, dtype=V[k].dtype)
    put("ro", rho0); put("P", a.P0); put("T", T0); put("roe", rho0 * e0)
    for k in ("Ux", "Uy", "Uz", "roUx", "roUy", "roUz"): put(k, 0.0)
    ys = sorted(k for k in V if re.fullmatch(r"Y\d+", k))
    for k in ys:
        yi = float(np.asarray(V[k][i])); put(k, yi); put("ro" + k, rho0 * yi)
    if "k" in V: put("k", 1.0e-6); put("roK", rho0 * 1.0e-6)
    if "omega" in V: put("omega", 10.0); put("roOmega", rho0 * 10.0)
    if "h0" in V: put("h0", e0 + R * T0)
print(f"[still_field] {run}: P0 {a.P0}・T0 300・ρ0 {rho0:.6g}・R {R:.6g}・e(300 K) {e0:.6g} (等温壁の節点 {i} から)")
b = run / "bcondConfig.yaml"; s = b.read_text()
s = re.sub(r"(Pt:\s*)[0-9.eE+-]+", lambda m: m.group(1) + f"{a.P0}", s)
s = re.sub(r"(Ps:\s*)[0-9.eE+-]+", lambda m: m.group(1) + f"{a.P0}", s)
s = re.sub(r"(Tt:\s*)[0-9.eE+-]+", lambda m: m.group(1) + "300.0", s)
b.write_text(s); print(f"[still_field] {b}: Pt・Ps = {a.P0}、Tt = 300.0 に書き換えた")
