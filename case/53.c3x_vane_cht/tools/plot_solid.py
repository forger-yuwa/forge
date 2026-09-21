#!/usr/bin/env python3
r"""連成 CHT の固体内部の温度場と固体メッシュを描く (報告の「The solid side」の図)。

界面壁温は連成 run の最終反復の `wall_profile_<phys>.csv`、内部は `Fem2DOperator.recover_interior` で復元する
(連成ループが解いたのと同じ作用素)。カラーマップは skill `forge-contour` に従い turbo。

usage: python3 case/53.c3x_vane_cht/tools/plot_solid.py <cht_run_dir> --solid solid_smooth.json --out solid.png
"""
import argparse, glob, json, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "solver_density_cuda/tools"))
from cht_loop import build_fem2d          # noqa: E402


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("run_dir"); ap.add_argument("--solid", required=True)
    ap.add_argument("--out", required=True); ap.add_argument("--phys", default="5")
    a = ap.parse_args()
    f = sorted(glob.glob(str(Path(a.run_dir) / f"it_*/wall_profile_{a.phys}.csv")))[-1]
    d = np.loadtxt(f, skiprows=1)
    op, perm = build_fem2d(json.load(open(a.solid)), d[:, :3])
    Tif = d[perm, 3]
    for _ in range(8):                      # k_s(T) は前回の内部温度で評価されるので数回まわして固定点に
        op.assemble(Tif); T = op.recover_interior(Tif)
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    from matplotlib.tri import Triangulation
    X = np.asarray(op.xy)[:, :2] * 100.0
    tri = Triangulation(X[:, 0], X[:, 1], np.asarray(op.tris))
    fig, ax = plt.subplots(1, 2, figsize=(12, 6.2))
    cs = ax[0].tricontourf(tri, T, levels=40, cmap="turbo")
    cb = fig.colorbar(cs, ax=ax[0], orientation="horizontal", location="bottom", fraction=0.05, pad=0.08)
    cb.set_label("metal temperature [K]")
    ax[0].set_title(f"{Path(a.run_dir).name}: {T.min():.0f}–{T.max():.0f} K", fontsize=10)
    ax[1].triplot(tri, lw=0.25, color="0.45")
    ax[1].set_title(f"conduction mesh: {len(X)} nodes, {len(op.tris)} triangles", fontsize=10)
    for x in ax: x.set_aspect("equal"); x.set_xlabel("x [cm]"); x.tick_params(labelsize=8)
    ax[0].set_ylabel("y [cm]")
    fig.tight_layout(); fig.savefig(a.out, dpi=110)
    print(f"[plot_solid] {f} -> {a.out}   T {T.min():.1f}-{T.max():.1f} K")


if __name__ == "__main__":
    main()
