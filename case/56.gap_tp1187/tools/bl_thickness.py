#!/usr/bin/env python3
"""forge の境界層積分厚さを Table II の実測 $\\delta^*$ と比べる (分母とは独立の検算)。

TP-1187 Table II は run ごとに $\\delta^*$ を位置 I / II で載せている。これは $q_{FP}$ とは
**別の測定量**なので、forge の平板が正しい BL を作れているかの独立な確認になる。

定義は圧縮性 (質量流束): $\\delta^*=\\int(1-\\rho u/\\rho_e U_e)\\,dy$。
(case/51 で TN D-8233 の $\\theta$ が圧縮性定義と両立しなかった件があるので、ここでは
運動学定義も併記する。)
"""
import argparse, json, sys
from pathlib import Path
import numpy as np, h5py

CASE = Path(__file__).resolve().parents[1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--res", default=None)
    ap.add_argument("--x-trip", type=float, default=0.13)
    a = ap.parse_args()
    rd = CASE / a.run
    if list(rd.glob("res_nan_*.h5")):
        raise SystemExit(f"REFUSED: {a.run} は発散している")
    files = sorted(rd.glob("res_[0-9]*.h5"), key=lambda f: int(f.stem.split("_")[1]))
    res = rd / a.res if a.res else files[-1]
    setup = json.loads((rd / "case_setup.json").read_text(encoding="utf-8"))
    cond = json.loads((CASE / "conditions.json").read_text(encoding="utf-8"))
    ser = [s for s in cond["series"] if s["run"] == setup["tp1187_run"]][0]
    trip = 0.0 if setup.get("laminar") else a.x_trip

    with h5py.File(res, "r") as f:
        c = f["/MESH/COORD"][:].reshape(-1, 3)
        ro = f["/VALUE/ro"][:].astype(float)
        u = f["/VALUE/Ux"][:].astype(float)
    print(f"run {a.run}  res {res.name}  TP-1187 run {setup['tp1187_run']}")
    print(f"\n{'位置':>4} {'x 実':>7} {'x CFD':>7} {'δ* 圧縮性':>10} {'δ* 実測':>9} {'比':>6} "
          f"{'δ* 運動学':>10} {'δ99':>8} {'δ*/δ99':>8}")
    for loc, x_phys in setup["x_eval_m"].items():
        ds_m = ser.get(f"dstar_loc_{loc}_cm")
        x = x_phys - trip
        m = np.abs(c[:, 0] - x) < 2e-3
        if ds_m is None or not m.any():
            continue
        i = np.where(m)[0]
        x0 = c[i, 0][np.argmin(np.abs(c[i, 0] - x))]
        i = i[np.abs(c[i, 0] - x0) < 1e-9]
        i = i[np.argsort(c[i, 1])]
        y, uu, rr = c[i, 1], u[i], ro[i]
        # **縁は領域上端ではなく局所 BL 端**。鋭利前縁の弱い衝撃で BL 外側が ~7 % 圧縮されて
        # いるので、自由流の rho を使うと積分が打ち消し合って delta* が 0 になる (2026-09-20 に踏んだ)。
        d99_0 = y[int(np.argmax(uu >= 0.99 * np.max(uu)))]
        j = int(np.argmax(y >= 2.0 * d99_0))
        Ue, roe = float(uu[j]), float(rr[j])
        d99 = y[int(np.argmax(uu[:j + 1] >= 0.99 * Ue))]
        ds_c = np.trapz(1.0 - rr[:j + 1] * uu[:j + 1] / (roe * Ue), y[:j + 1])
        ds_k = np.trapz(1.0 - uu[:j + 1] / Ue, y[:j + 1])
        print(f"{loc:>4} {x_phys:7.2f} {x0:7.2f} {ds_c*1e2:10.3f} {ds_m:9.2f} "
              f"{ds_c*1e2/ds_m:6.3f} {ds_k*1e2:10.3f} {d99*1e2:8.3f} {ds_c/d99:8.3f}")
    print("\n単位は cm。δ* 実測は Table II。縁は局所 BL 端 (2 x δ99) で取っている。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
