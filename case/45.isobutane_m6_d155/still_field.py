"""一様な静止場 (plan architecture-float-state-double-geometry §6.3 の 3: 軸対称の closure の非劣化) を run の場のファイルに書く。
等温壁 (300 K) の節点から組成・気体定数 R = P/(ρT)・内部エネルギー e = ρE/ρ (u = 0) を取り、全節点を P0・T0・u = 0 にする。
入口・出口の境界の値も静止場と矛盾しない値に書き換える (入口 Pt = P0・Tt = T0、出口 Ps = Pt = P0)。
run の場のファイルが保存量だけ (restart の最小、P・T を持たない) のときは、組成・R・e を --state-from の res (P・T を持つ出力) の等温壁の節点から取る。
usage: python3 still_field.py <run_dir> [--P0 1.0e5] [--state-from <res.h5>]"""
import argparse, re, h5py, numpy as np
from pathlib import Path
ap = argparse.ArgumentParser(); ap.add_argument("run"); ap.add_argument("--P0", type=float, default=1.0e5); ap.add_argument("--state-from", default=None); a = ap.parse_args()
run = Path(a.run); f = run / "nozzle.h5"
with h5py.File(f, "r+") as h:
    V = h["VALUE"]
    S = V if "T" in V else (h5py.File(a.state_from, "r")["VALUE"] if a.state_from else None)
    if S is None: raise SystemExit(f"{f} に T が無い (保存量だけ)。--state-from に P・T を持つ res を渡す")
    if S["ro"].shape != V["ro"].shape: raise SystemExit("--state-from の節点数が run の場と違う")
    T = np.asarray(S["T"][:]); P = np.asarray(S["P"][:]); ro = np.asarray(S["ro"][:]); roe = np.asarray(S["roe"][:])
    U2 = np.asarray(S["Ux"][:]) ** 2 + np.asarray(S["Uy"][:]) ** 2
    w = np.flatnonzero((np.abs(T - 300.0) < 1e-6) & (U2 < 1e-20))
    if w.size == 0: raise SystemExit("T = 300 K・u = 0 の節点 (等温壁) が無い")
    i = w[0]; T0 = 300.0; R = P[i] / (ro[i] * T[i]); e0 = roe[i] / ro[i]
    rho0 = a.P0 / (R * T0); n = ro.size
    def put(k, v):
        if k in V: V[k][...] = np.full(n, v, dtype=V[k].dtype)
    put("ro", rho0); put("P", a.P0); put("T", T0); put("roe", rho0 * e0)
    for k in ("Ux", "Uy", "Uz", "roUx", "roUy", "roUz"): put(k, 0.0)
    ys = sorted(k for k in S if re.fullmatch(r"Y\d+", k))
    for k in ys:
        yi = float(np.asarray(S[k][i])); put(k, yi); put("ro" + k, rho0 * yi)
    put("k", 1.0e-6); put("roK", rho0 * 1.0e-6)
    put("omega", 10.0); put("roOmega", rho0 * 10.0)
    if "h0" in V: put("h0", e0 + R * T0)
print(f"[still_field] {run}: P0 {a.P0}・T0 300・ρ0 {rho0:.6g}・R {R:.6g}・e(300 K) {e0:.6g} (等温壁の節点 {i} から)")
b = run / "bcondConfig.yaml"; s = b.read_text()
s = re.sub(r"(Pt:\s*)[0-9.eE+-]+", lambda m: m.group(1) + f"{a.P0}", s)
s = re.sub(r"(Ps:\s*)[0-9.eE+-]+", lambda m: m.group(1) + f"{a.P0}", s)
s = re.sub(r"(Tt:\s*)[0-9.eE+-]+", lambda m: m.group(1) + "300.0", s)
b.write_text(s); print(f"[still_field] {b}: Pt・Ps = {a.P0}、Tt = 300.0 に書き換えた")
