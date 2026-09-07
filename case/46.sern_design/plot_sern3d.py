#!/usr/bin/env python3
"""3D SERN run の可視化: ランプ面圧力マップ (x–z) と、z=0 / z=W/2 直外 / 対称面の Mach 断面。使い方: plot_sern3d.py <run_dir> <problem.yaml>"""
import sys, json, glob
from pathlib import Path
import numpy as np, h5py, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib import font_manager
for f in Path.home().glob(".fonts/NotoSansCJKjp-Regular.otf"): font_manager.fontManager.addfont(str(f)); plt.rcParams["font.family"] = "Noto Sans CJK JP"
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "design"))
from forge_design.probdef import load_problem
run = Path(sys.argv[1]); p = load_problem(sys.argv[2]); info = json.loads((run / "prepare_info.json").read_text()); H = info["H_m"]; g = p.gamma
ex = info["states"]["exhaust"]; en = info["states"]["ext"]; hw = info["half_W_m"]
res = sorted(run.glob("res_[0-9]*.h5"), key=lambda f: int("".join(c for c in f.stem if c.isdigit())))[-1]; step = int("".join(c for c in res.stem if c.isdigit()))
# --- ランプ面圧力 (壁面出力) ---
wf = list(run.glob(f"res_ramp_4_{step}.h5"))[0]
with h5py.File(wf) as f:
    xyz = f["MESH/COORD"][:].reshape(-1, 3); ps = f["VALUE/Ps"][:]
fig, axs = plt.subplots(2, 1, figsize=(11, 7))
ax = axs[0]; sc = ax.scatter(xyz[:, 0] / H, xyz[:, 2] / H, c=ps / ex["P"], s=3, cmap="viridis", vmin=0, vmax=0.6)
ax.axhline(hw / H, color="w", ls="--", lw=1); ax.set_xlabel("x / H"); ax.set_ylabel("z / H"); ax.set_title(f"ランプ面圧力 p/p_in ({run.name}, step {step}; 破線 = ノズル幅端 z=W/2)")
plt.colorbar(sc, ax=ax, label="p / p_in")
# --- Mach 断面 (y–x) at z=0 (対称面) と z=W/2+ ---
with h5py.File(run / "sern.h5") as f: cc = f["CELLS/centCoords"][:].reshape(-1, 3)
with h5py.File(res) as f:
    ro = f["VALUE/ro"][:]; u = f["VALUE/roUx"][:] / ro; v = f["VALUE/roUy"][:] / ro; w = f["VALUE/roUz"][:] / ro; roe = f["VALUE/roe"][:]
P = (g - 1) * (roe - 0.5 * ro * (u * u + v * v + w * w)); T = P / (ro * p.R_gas); M = np.sqrt(u * u + v * v + w * w) / np.sqrt(g * p.R_gas * T)
ax = axs[1]
for zsel, col, lab in ((0.0, None, "z=0 (対称面)"),):
    thr = np.sort(cc[:, 2])[int(len(cc) / info["mesh"]["nz"] * 0.9)]; m = cc[:, 2] <= thr   # 最初の z 層 (対称面側)
    from matplotlib.tri import Triangulation
    tri = Triangulation(cc[m, 0] / H, cc[m, 1] / H); cf = ax.tricontourf(tri, np.clip(M[m], 0, 8), levels=np.linspace(0, 8, 33), cmap="viridis")
plt.colorbar(cf, ax=ax, label="Mach (z≈0)"); ax.set_aspect("equal"); ax.set_xlabel("x / H"); ax.set_ylabel("y / H"); ax.set_ylim(-4, 4)
fig.tight_layout(); fig.savefig(run / "ramp_pressure_map.png", dpi=140); print("saved", run / "ramp_pressure_map.png")
