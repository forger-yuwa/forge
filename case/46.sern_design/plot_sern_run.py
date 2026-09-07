#!/usr/bin/env python3
"""SERN run の可視化: Mach コンター + 設計壁 + 壁圧 (CFD vs MOC)。使い方: plot_sern_run.py <run_dir> <problem.yaml>"""
import sys, json, subprocess
from pathlib import Path
import numpy as np, h5py, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "design"))
from forge_design.evaluate import runner_sern as R
from forge_design.probdef import load_problem
from forge_design.geometry.rao_planar import p_over_pt
from matplotlib import font_manager
for _f in Path.home().glob(".fonts/NotoSansCJKjp-Regular.otf"):
    font_manager.fontManager.addfont(str(_f)); plt.rcParams["font.family"] = "Noto Sans CJK JP"
run = Path(sys.argv[1]); p = load_problem(sys.argv[2]); info = json.loads((run / "prepare_info.json").read_text())
H = info["H_m"]; g = p.gamma; ex = info["states"]["exhaust"]; en = info["states"]["ext"]
res = sorted(run.glob("res_[0-9]*.h5"), key=lambda f: int("".join(c for c in f.stem if c.isdigit())))[-1]
step = int("".join(c for c in res.stem if c.isdigit()))
with h5py.File(run / "sern.h5") as f:
    cc = f["CELLS/centCoords"][:].reshape(-1, 3)
with h5py.File(res) as f:
    ro = f["VALUE/ro"][:]; u = f["VALUE/roUx"][:] / ro; v = f["VALUE/roUy"][:] / ro; roe = f["VALUE/roe"][:]
P = (g - 1) * (roe - 0.5 * ro * (u * u + v * v)); T = P / (ro * p.R_gas); M = np.sqrt(u * u + v * v) / np.sqrt(g * p.R_gas * T)
x, y = cc[:, 0] / H, cc[:, 1] / H
tri = Triangulation(x, y)
fig, ax = plt.subplots(figsize=(12, 6))
cf = ax.tricontourf(tri, np.clip(M, 0, 8), levels=np.linspace(0, 8, 41), cmap="viridis")
plt.colorbar(cf, ax=ax, label=f"Mach (forge {info['model']} {info['discretization']}, step {step})")
ramp = np.loadtxt(run / "ramp_contour.csv", delimiter=",", skiprows=1) / H; cowl = np.loadtxt(run / "cowl_contour.csv", delimiter=",", skiprows=1) / H
ax.plot(ramp[:, 0], ramp[:, 1], "r-", lw=2); ax.plot(cowl[:, 0], cowl[:, 1], "r-", lw=2)
ax.set_aspect("equal"); ax.set_xlabel("x / H"); ax.set_ylabel("y / H"); ax.set_title(run.name)
fig.tight_layout(); fig.savefig(run / "mach_field.png", dpi=140); plt.close(fig)
# 壁圧: CFD (壁出力) vs MOC (逆設計の壁状態)
k, d, fr, _ = R.design_from_problem(p)
fig, ax = plt.subplots(figsize=(10, 4)); p0 = float(p_over_pt(ex["M"], g))
ax.plot(d.ramp_xy[:, 0], p_over_pt(d.ramp_M, g) / p0, "r--", label="MOC ramp"); ax.plot(d.cowl_xy[:, 0], p_over_pt(d.cowl_M, g) / p0, "b--", label="MOC cowl")
for name, pid, col in (("ramp", 4, "r"), ("cowl_in", 5, "b"), ("cowl_out", 6, "g")):
    fs = list(run.glob(f"res_{name}_{pid}_{step}.h5"))
    if not fs: continue
    with h5py.File(fs[0]) as f:
        xy = f["MESH/COORD"][:].reshape(-1, 3); ps = f["VALUE/Ps"][:]; cn = f["MESH/CONNE"][:].reshape(-1, 4)[:, 2:4]
    xw = xy[:, 0] if len(ps) == len(xy) else 0.5 * (xy[cn[:, 0], 0] + xy[cn[:, 1], 0])
    o = np.argsort(xw); ax.plot(xw[o] / H, ps[o] / ex["P"], col + "-", lw=1, label=f"forge {name}")
ax.axhline(en["P"] / ex["P"], color="k", ls=":", label="p_ext"); ax.set_yscale("log"); ax.set_xlabel("x / H"); ax.set_ylabel("p / p_in")
ax.legend(fontsize=8); ax.grid(alpha=.3); ax.set_title("wall pressure: forge vs MOC"); fig.tight_layout(); fig.savefig(run / "wall_pressure_cfd_vs_moc.png", dpi=140)
print("saved", run / "mach_field.png", run / "wall_pressure_cfd_vs_moc.png")
