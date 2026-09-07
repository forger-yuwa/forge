#!/usr/bin/env python3
# Taylor-Green: 運動エネルギー K(t) と総エントロピー S(t) の時間履歴を res_*.h5 から計算してプロット。
#
# 旧 plot_taylorGreen.py は log の "rok"/"ros"/"Total Time" 行 (旧 nagare ソルバ出力) を読む作りで、
# 現 forge は同行を出力しないため動かない。本スクリプトはその現行版で、保存場 res_*.h5 から
# K=Σ½ρ|u|²V, S=Σρ·cv·ln(P/ρ^γ)·V を直接積分する。
import glob, os, re
import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def _cell_volume(f, path):
    """セル体積: res には既定 (output.level 1) で volume が無いので、mesh h5 (solverConfig の meshFileName) の CELLS/volume を読む。"""
    import os, yaml
    if "VALUE/volume" in f:
        return f["VALUE/volume"][:].astype(np.float64)
    rd = os.path.dirname(os.path.abspath(path))
    mesh = yaml.safe_load(open(os.path.join(rd, "solverConfig.yaml")))["mesh"]["meshFileName"]
    with h5py.File(os.path.join(rd, mesh), "r") as m:
        return m["CELLS/volume"][:].astype(np.float64)


gamma = 1.4; cp = 0.4; cv = cp / gamma; dt = 0.007

def step_of(p):
    m = re.search(r"res_(\d+)\.h5$", os.path.basename(p)); return int(m.group(1)) if m else -1

L = 2.0 * np.pi  # 周期箱の辺長

def cv_weight(f, nCV):
    # node (median-dual) では周期マージ双対セルが res に重複格納される (VALUE 長 == MESH/COORD ノード数)。
    # 各境界 CV を多重度 m=2^(境界方向数) で割って一意 CV 1 個分に補正する。cell では VALUE 長 != ノード数
    # かつ重心が境界上に乗らないので weight=1 (no-op)。
    coord = f["MESH/COORD"][:].astype(np.float64)
    if coord.size != 3 * nCV:
        return np.ones(nCV)
    C = coord.reshape(-1, 3); tol = 1e-5
    m = np.ones(nCV)
    for k in range(3):
        onb = (np.abs(C[:, k]) < tol) | (np.abs(C[:, k] - L) < tol)
        m = m * np.where(onb, 2.0, 1.0)
    return 1.0 / m

def totals(path):
    with h5py.File(path, "r") as f:
        V = _cell_volume(f, path)
        ro = f["VALUE/ro"][:].astype(np.float64)
        P  = f["VALUE/P"][:].astype(np.float64)
        ux = f["VALUE/Ux"][:].astype(np.float64); uy = f["VALUE/Uy"][:].astype(np.float64); uz = f["VALUE/Uz"][:].astype(np.float64)
        w  = cv_weight(f, V.size)
    Vw = V * w
    K = np.sum(0.5 * ro * (ux*ux + uy*uy + uz*uz) * Vw)
    S = np.sum(ro * (cv * np.log(P / np.power(ro, gamma))) * Vw)
    return K, S

def series(run_dir):
    fs = sorted([p for p in glob.glob(os.path.join(run_dir, "res_*.h5"))
                 if step_of(p) >= 0 and "nan" not in p], key=step_of)
    t = []; K = []; S = []
    for p in fs:
        k, s = totals(p); t.append(step_of(p)*dt); K.append(k); S.append(s)
    return np.array(t), np.array(K), np.array(S)

RUNS = [
    ("node pure-KEEP visc (Re160)", "run_0007_node_keep_pure_visc", "C0", "-"),
    ("cell pure-KEEP visc (Re160)", "run_0008_cell_pure_visc",      "C1", "--"),
    ("node pure-KEEP inviscid",     "run_0009_node_keep_pure_eul",  "C2", "-"),
    ("cell pure-KEEP inviscid",     "run_0010_cell_pure_eul",       "C3", "--"),
]

base = os.path.dirname(os.path.abspath(__file__))
fig, ax = plt.subplots(1, 2, figsize=(13, 5))
for nm, d, c, ls in RUNS:
    t, K, S = series(os.path.join(base, d))
    if len(t) == 0:
        continue
    ax[0].plot(t, K/K[0], ls, color=c, lw=1.6, marker="o", ms=2.5, label=nm)
    ax[1].plot(t, (S-S[0])/abs(S[0]), ls, color=c, lw=1.6, marker="o", ms=2.5, label=nm)
ax[0].axhline(1.0, ls=":", c="k", lw=0.8)
ax[0].set_xlabel("time [s]"); ax[0].set_ylabel("K(t) / K(0)")
ax[0].set_title("Total kinetic energy"); ax[0].grid(alpha=0.3); ax[0].legend(fontsize=8)
ax[1].axhline(0.0, ls=":", c="k", lw=0.8)
ax[1].set_xlabel("time [s]"); ax[1].set_ylabel("(S(t) - S(0)) / |S(0)|")
ax[1].set_title("Total entropy change"); ax[1].grid(alpha=0.3); ax[1].legend(fontsize=8)
fig.suptitle("Taylor-Green 32^3 — KE & entropy time history (explicit RK4)")
fig.tight_layout()
out = os.path.join(base, "ke_entropy_history.png")
fig.savefig(out, dpi=120)
print("saved", out)

# 数値サマリ
print(f"\n{'run':32s}{'t_end':>7}{'K/K0_end':>10}{'(S-S0)/|S0|_end':>17}")
for nm, d, c, ls in RUNS:
    t, K, S = series(os.path.join(base, d))
    if len(t) == 0: continue
    print(f"{nm:32s}{t[-1]:7.2f}{K[-1]/K[0]:10.4f}{(S[-1]-S[0])/abs(S[0]):17.3e}")
