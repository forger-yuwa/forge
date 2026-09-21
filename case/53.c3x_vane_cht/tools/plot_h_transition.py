#!/usr/bin/env python3
r"""翼面の熱伝達率 $h(s)$ を、乱流モデルの違う run どうしで実測に重ねる (遷移モデル ON/OFF の比較図)。

横軸は前縁からの弧長 $s/S$ (負圧面を正、正圧面を負に取る)。$h=q/(T_g-T_w)$、$q$ は `iface_q_eff`。規約は `compare_h.py` と同じ。

usage: plot_h_transition.py --run run108|run42 --out fig.png RUN_DIR:label [RUN_DIR:label ...] [--gamma RUN_DIR]
"""
import argparse, sys
from pathlib import Path
import numpy as np, h5py
sys.path.insert(0, str(Path(__file__).parent))
import run_data as rd_mod
from compare_h import arc_map, TABLES, H0, TG


def load(rd, phys=5, flux="iface_q_eff"):
    fs = sorted(Path(rd).glob(f"res_wall_{phys}_[0-9]*.h5"), key=lambda p: int(p.stem.split("_")[-1]))
    with h5py.File(fs[-1], "r") as f:
        C = np.array(f["MESH/COORD"]).reshape(-1, 3); q = np.array(f["VALUE"][flux]); Tw = np.array(f["VALUE/Ts"]); ok = np.array(f["VALUE/iface_ok"]) > 0.5
    s, ss = arc_map(C[:, :2]); x = np.where(ss, s, -s); o = np.argsort(x)
    return x[o], q[o], Tw[o], ok[o], fs[-1]


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("runs", nargs="+"); ap.add_argument("--run", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--title", default=""); ap.add_argument("--ymax", type=float, default=1500.0)
    a = ap.parse_args()
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    from matplotlib import font_manager
    for fp in Path.home().glob(".fonts/NotoSansCJK*"): font_manager.fontManager.addfont(str(fp))
    if any("Noto Sans CJK JP" in f.name for f in font_manager.fontManager.ttflist): plt.rcParams["font.family"] = "Noto Sans CJK JP"
    coarse = {i for i, d in getattr(rd_mod, "RUN42_PARTIAL", {}).items() if a.run == "run42" and d["tol"] > 0.01}
    rows = [r for i, r in enumerate(TABLES[a.run]["rows"]) if r[3] is not None and i not in coarse]
    sd = np.array([r[0] for r in rows]); hd = np.array([r[3] for r in rows]) * H0; i0 = int(np.argmin(sd))
    xe = np.concatenate([-sd[:i0 + 1], sd[i0:]]); he = np.concatenate([hd[:i0 + 1], hd[i0:]])
    fig, ax = plt.subplots(figsize=(10, 5.2))
    ax.plot(xe, he, "ko", ms=5, mfc="none", label="実測 (NASA CR-168015)")
    for spec in a.runs:
        rd, _, lab = spec.partition(":"); x, q, Tw, ok, fn = load(rd); h = q / (TG[a.run] - Tw)
        ax.plot(x[ok], h[ok], lw=1.6, label=lab or Path(rd).name); print(f"{lab:<28} {fn}")
    ax.axvline(0, color="0.6", lw=0.8); ax.set_xlim(-1, 1); ax.set_ylim(0, a.ymax)
    ax.set_xlabel("前縁からの弧長 s/S   (← 正圧面 | 負圧面 →)"); ax.set_ylabel("熱伝達率 h [W/m²K]"); ax.grid(alpha=0.3)
    if a.title: ax.set_title(a.title)
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False); fig.tight_layout(); fig.savefig(a.out, dpi=140); print("wrote", a.out)


if __name__ == "__main__":
    main()
