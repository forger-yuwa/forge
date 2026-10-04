#!/usr/bin/env python3
r"""冷却孔の配置 (図 6/7 の (U,V) 系の定義どおり) を翼外形に重ねて描く。旧配置 (最小肉厚の最適化) を破線で併記。

usage: python3 case/53.c3x_vane_cht/tools/plot_holes.py --out holes.png
"""
import argparse, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).parent))
from check_geometry import _read_xy
ROOT = Path(__file__).resolve().parents[3]
V = {"C3X": ROOT / "case/53.c3x_vane_cht/ref", "Mark II": ROOT / "case/54.markii_vane_cht/ref"}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out", required=True); a = ap.parse_args()
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11, 7.5))
    for x, (name, ref) in zip(ax, V.items()):
        key = "c3x" if name == "C3X" else "markii"
        P = _read_xy(ref / f"vane_{key}_smooth.csv")
        H = np.array([[float(t) for t in l.strip().split(",")] for l in open(ref / f"cooling_holes_{key}.csv") if l[0].isdigit()])
        x.plot(*np.vstack([P, P[:1]]).T, "k", lw=1.2)
        th = np.linspace(0, 2 * np.pi, 60)
        for h in H:
            x.plot(h[5] + h[3] / 2 * np.cos(th), h[6] + h[3] / 2 * np.sin(th), color="tab:red", lw=1.4)
            x.plot(h[7] + h[3] / 2 * np.cos(th), h[8] + h[3] / 2 * np.sin(th), color="0.55", lw=0.9, ls="--")
            x.text(h[5], h[6], f"{int(h[0])}", ha="center", va="center", fontsize=7, color="tab:red")
        sh = np.hypot(H[:, 5] - H[:, 7], H[:, 6] - H[:, 8]).mean()
        x.set_aspect("equal"); x.set_title(f"{name}: holes moved {sh*10:.1f} mm along the chord", fontsize=10)
        x.set_xlabel("x [cm]"); x.grid(alpha=.3)
    ax[0].set_ylabel("y [cm]")
    ax[0].plot([], [], color="tab:red", lw=1.4, label="placed by the report's (U,V) frame"); ax[0].plot([], [], color="0.55", ls="--", label="earlier placement")
    fig.legend(loc="lower center", ncol=2, fontsize=9, frameon=False); fig.tight_layout(rect=[0, 0.04, 1, 1]); fig.savefig(a.out, dpi=110)
    print(f"[plot_holes] -> {a.out}")


if __name__ == "__main__":
    main()
