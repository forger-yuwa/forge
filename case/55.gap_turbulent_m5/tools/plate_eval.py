#!/usr/bin/env python3
"""case/55 平板評価 — q_w(x) と T_w(x) を場から出す (断熱 run では T_w = T_aw)。

q_w = λ_w (dT/dn)_w。λ は forge の出力 `thermCond` (viscMethod 2) をそのまま使う。
dT/dn は壁ノードから内側 2 点の 2 次片側差分。case/50 の抽出器と同じ式。
"""
import argparse, json, sys
from pathlib import Path
import numpy as np
import h5py

CASE = Path(__file__).resolve().parents[1]


def dTdn_wall(y, T):
    h1, h2 = y[1] - y[0], y[2] - y[0]
    a = (T[1] - T[0]) / h1
    b = ((T[2] - T[0]) / h2 - a) / (h2 - h1)
    return a - b * h1


def evaluate(res, x_plate_end=0.450):
    with h5py.File(res, "r") as f:
        c = f["/MESH/COORD"][:].reshape(-1, 3)
        T = f["/VALUE/T"][:].astype(float)
        lam = f["/VALUE/thermCond"][:].astype(float)
        U = f["/VALUE/Ux"][:].astype(float)
        P = f["/VALUE/P"][:].astype(float)
    xs = np.unique(np.round(c[:, 0], 9))
    out = []
    for x in xs:
        if x < 1e-9 or x > x_plate_end + 1e-9:
            continue
        m = (np.abs(c[:, 0] - x) < 1e-9) & (c[:, 1] >= -1e-12)
        i = np.where(m)[0]
        if len(i) < 3:
            continue
        i = i[np.argsort(c[i, 1])]
        y = c[i, 1]
        if y[0] > 1e-9:
            continue
        q = lam[i[0]] * dTdn_wall(y[:3], T[i[:3]])
        # δ99 (0.99 U∞ の最初の y)
        u = U[i]
        Uinf = float(np.max(u))
        j = np.argmax(u >= 0.99 * Uinf)
        out.append((x, q, T[i[0]], lam[i[0]], P[i[0]], y[j] if j > 0 else np.nan))
    return np.array(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--res", default=None)
    ap.add_argument("--x-ref", type=float, default=0.34875)
    a = ap.parse_args()
    rd = CASE / a.run
    if list(rd.glob("res_nan_*.h5")):
        raise SystemExit(f"REFUSED: {a.run} は発散している")
    files = sorted(rd.glob("res_[0-9]*.h5"), key=lambda f: int(f.stem.split("_")[1]))
    res = rd / a.res if a.res else (files[-1] if files else rd / "mesh.h5")
    A = evaluate(res)
    np.savetxt(rd / "wall_q.csv", A, delimiter=",",
               header="x_m,q_w_W_m2,Tw_K,lambda_w,P_w,delta99_m", comments="")
    i = int(np.argmin(np.abs(A[:, 0] - a.x_ref)))
    setup = json.loads((rd / "case_setup.json").read_text(encoding="utf-8"))
    print(f"run {a.run}  res {res.name}  壁 {len(A)} 点 (x {A[0,0]*1e3:.1f}–{A[-1,0]*1e3:.1f} mm)"
          f"  平板 {setup.get('plate_thermal','?')}")
    print(f"  x = {A[i,0]*1e3:.2f} mm :  q_w = {A[i,1]*1e-3:8.3f} kW/m²   T_w = {A[i,2]:8.2f} K"
          f"   δ99 = {A[i,5]*1e3:6.3f} mm ({A[i,5]/2.5e-3:.2f} W)")
    if setup.get("plate_thermal") == "adiabatic":
        print(f"  → T_aw (実測) = {A[i,2]:.2f} K   (case/49 derived の TP 推定 1126.33 K と "
              f"{100*(A[i,2]/1126.33-1):+.2f} %)")
    print(f"  → {rd/'wall_q.csv'}")


if __name__ == "__main__":
    main()
