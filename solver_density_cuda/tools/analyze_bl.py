#!/usr/bin/env python3
"""Boundary-layer diagnostics for the laminar channel / pipe verification cases.

Extracts wall-normal velocity profiles Ux(y) at several streamwise stations
from a forge result HDF5, plots them together with a Ux contour, and prints
bulk quantities (mean velocity, Re_D, Mach) so we can tell a real boundary
layer (Ux -> 0 at the wall) from the reported "wall acceleration" artefact.

Usage:
  analyze_bl.py --mesh pipe.h5 --result res_10.h5 --label "axisym pipe" \
                --half 0.005 --axisym 1 --out bl
"""
import argparse
from pathlib import Path

import h5py
import numpy as np


def _trapz(y, x):
    """台形積分。**numpy 2 で `np.trapz` が削除された**ので互換に包む (AWS の numpy で落ちた)。"""
    f = getattr(np, "trapezoid", None) or getattr(np, "trapz")
    return f(y, x)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load(mesh_path, result_path):
    with h5py.File(mesh_path, "r") as m:
        c = m["CELLS/centCoords"][:].reshape(-1, 3).astype(np.float64)
    with h5py.File(result_path, "r") as r:
        ux = r["VALUE/Ux"][:].astype(np.float64)
        uy = r["VALUE/Uy"][:].astype(np.float64)
        ro = r["VALUE/ro"][:].astype(np.float64)
        P = r["VALUE/P"][:].astype(np.float64)
    return c[:, 0], c[:, 1], ux, uy, ro, P


def nearest_column(x, x_target):
    xs = np.unique(np.round(x, 9))
    return xs[np.argmin(np.abs(xs - x_target))]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mesh", required=True)
    ap.add_argument("--result", required=True)
    ap.add_argument("--label", default="case")
    ap.add_argument("--half", type=float, required=True, help="wall-to-symmetry distance [m]")
    ap.add_argument("--length", type=float, default=0.1)
    ap.add_argument("--axisym", type=int, default=0)
    ap.add_argument("--cp", type=float, default=1004.5)
    ap.add_argument("--gamma", type=float, default=1.4)
    ap.add_argument("--visc", type=float, default=4.08e-3)
    ap.add_argument("--out", default="bl")
    args = ap.parse_args()

    x, y, ux, uy, ro, P = load(args.mesh, args.result)
    R = args.cp * (args.gamma - 1.0) / args.gamma
    T = P / (ro * R)
    a = np.sqrt(args.gamma * R * T)

    stations = [0.20, 0.50, 0.90]  # fractions of length
    fig, (axp, axc) = plt.subplots(1, 2, figsize=(12, 5))

    for frac in stations:
        xt = nearest_column(x, frac * args.length)
        m = np.isclose(x, xt)
        order = np.argsort(y[m])
        yy = y[m][order]
        uu = ux[m][order]
        axp.plot(uu, yy / args.half, "-o", ms=3, label=f"x/L={frac:.2f}")

    axp.axhline(1.0, color="k", lw=0.8, ls="--")
    axp.set_xlabel("Ux [m/s]")
    axp.set_ylabel("y / (wall distance)   (0=axis/center-ish, 1=wall)")
    axp.set_title(f"{args.label}: wall-normal Ux profiles")
    axp.legend()
    axp.grid(alpha=0.3)

    # contour from scattered centroids
    sc = axc.scatter(x, y, c=ux, s=4, cmap="viridis")
    axc.set_xlabel("x [m]")
    axc.set_ylabel("y [m]")
    axc.set_title(f"{args.label}: Ux field")
    fig.colorbar(sc, ax=axc, label="Ux [m/s]")

    fig.tight_layout()
    out_png = Path(f"{args.out}_profiles.png")
    fig.savefig(out_png, dpi=130)

    # bulk quantities at the outlet column
    xt = nearest_column(x, 0.95 * args.length)
    m = np.isclose(x, xt)
    yy = y[m]
    uu = ux[m]
    rr = ro[m]
    order = np.argsort(yy)
    yy, uu, rr = yy[order], uu[order], rr[order]
    if args.axisym:
        # area-weighted (annular) mean: integral u*2*pi*r dr / (pi R^2)
        num = _trapz(uu * rr * yy, yy)
        den = _trapz(rr * yy, yy)
        umean = num / den if den != 0 else float("nan")
        D = 2.0 * args.half
    else:
        num = _trapz(uu * rr, yy)
        den = _trapz(rr, yy)
        umean = num / den if den != 0 else float("nan")
        D = 2.0 * args.half
    romean = float(np.mean(rr))
    Re_D = romean * umean * D / args.visc
    ucl = float(uu[np.argmin(np.abs(yy - (0.0 if args.axisym else args.half)))])
    uwall = float(uu[np.argmax(yy)])  # cell nearest the wall (largest y)
    Mmean = umean / float(np.mean(a[m]))

    print(f"=== {args.label} (x/L~0.95) ===")
    print(f"  mean Ux        : {umean:.3f} m/s")
    print(f"  centre/axis Ux : {ucl:.3f} m/s")
    print(f"  near-wall Ux   : {uwall:.3f} m/s   (should be ~0 for a real BL)")
    print(f"  ro mean        : {romean:.4f} kg/m3")
    print(f"  Re_D           : {Re_D:.1f}")
    print(f"  mean Mach      : {Mmean:.4f}")
    print(f"  wrote          : {out_png}")


if __name__ == "__main__":
    main()
