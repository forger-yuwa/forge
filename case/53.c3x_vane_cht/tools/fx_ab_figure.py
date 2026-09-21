#!/usr/bin/env python3
r"""面補間重み fx の A/B (幾何 fx 対 fx=0.5) を、壁熱流束の負圧面分布で SU2 と並べて描く。

上段: 壁熱流束 [kW/m²]。下段: 41 点移動平均を引いた残り (うねり)。凡例は図の外。
plan boundary-conjugate-heat-transfer §5.1 #59/#60。

usage: python3 case/53.c3x_vane_cht/tools/fx_ab_figure.py --ctrl RUN --trial RUN --out fig.png
"""
import argparse, sys
from pathlib import Path
import numpy as np, h5py
sys.path.insert(0, str(Path(__file__).parent))
from compare_h import arc_map
from resid_split import smooth
from report_figures import read_su2, C3, SU2_CTRL


def wall(rd):
    fs = sorted(Path(rd).glob("res_wall_5_[0-9]*.h5"), key=lambda p: int(p.stem.split("_")[-1]))
    with h5py.File(fs[-1]) as f:
        C = np.array(f["MESH/COORD"], float).reshape(-1, 3)[:, :2]
        return C, np.array(f["VALUE/iface_q_eff"], float), np.array(f["VALUE/iface_q_compact"], float)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--ctrl", required=True); ap.add_argument("--trial", required=True)
    ap.add_argument("--out", required=True); a = ap.parse_args()
    from scipy.spatial import cKDTree
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    C, qe0, qc0 = wall(a.ctrl); C1, qe1, qc1 = wall(a.trial); assert np.abs(C - C1).max() == 0
    XY, QS = read_su2(C3 / SU2_CTRL["su2"]); d, idx = cKDTree(XY).query(C); assert d.max() == 0.0
    s, ss = arc_map(C); o = np.argsort(s[ss]); sS = s[ss][o]; g = lambda y: y[ss][o] / 1e3
    cur = [("forge, geometric face weight (current)", g(qe0), "tab:red", 1.4),
           ("forge, face weight 0.5 (trial)", g(qe1), "k", 1.8),
           ("SU2 8.5, same mesh (wall-gradient flux)", g(QS[idx]), "tab:green", 1.4)]
    fig, ax = plt.subplots(2, 1, figsize=(10.5, 7.2), sharex=True)
    m = (sS >= 0.45) & (sS <= 0.95)
    for lab, y, c, lw in cur:
        w = y - smooth(y, 41)
        ax[0].plot(sS, y, color=c, lw=lw, label=lab)
        ax[1].plot(sS, w, color=c, lw=lw, label=f"{lab}\nrms {np.sqrt((w[m]**2).mean()):.2f} kW/m$^2$ on 0.45–0.95")
    ax[0].set_ylabel("wall heat flux [kW/m$^2$]"); ax[1].set_ylabel("flux minus 41-node moving mean [kW/m$^2$]")
    ax[1].set_xlabel("s/S, suction side"); ax[0].set_xlim(0.2, 1.0); ax[0].set_ylim(140, 205); ax[1].set_ylim(-6, 6)
    for x in ax: x.grid(alpha=.3); x.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=8, frameon=False)
    ax[0].set_title(f"C3X, uniform $T_w$ = 566 K — {Path(a.ctrl).name} vs {Path(a.trial).name}", fontsize=10)
    fig.tight_layout(); fig.savefig(a.out, dpi=110)
    print(f"[fx_ab_figure] -> {a.out}")


if __name__ == "__main__":
    main()
