#!/usr/bin/env python3
r"""連成 CHT の壁温を実測と比べる図 (V5 段 (c))。

usage: python3 case/53.c3x_vane_cht/tools/plot_tw.py <cht_run_dir> --run run108|run42 [--band DIR ...]
"""
import argparse, glob, sys
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from run_data import TABLES          # noqa: E402
from compare_h import arc_map, CASE  # noqa: E402


def load(run_dir):
    f = sorted(glob.glob(str(Path(run_dir) / "it_*/wall_profile_5.csv")))[-1]
    d = np.loadtxt(f, skiprows=1)
    s, ss = arc_map(d[:, :2])
    return s, ss, d[:, 3], f


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("--run", default="run108")
    ap.add_argument("--band", nargs="*", default=[], help="帯の run ディレクトリ (Tc / k_s 振り)")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    t = TABLES[a.run]
    rows = [r for r in t["rows"] if r[2] is not None]
    sd = np.array([r[0] for r in rows]); td = np.array([r[2] for r in rows]) * t["Tref"]
    i = int(np.argmin(sd))
    exp = {"PS": (sd[:i + 1][::-1], td[:i + 1][::-1]), "SS": (sd[i:], td[i:])}

    fig, ax = plt.subplots(figsize=(10, 5.8))
    for side, sgn, col in (("SS", +1, "tab:red"), ("PS", -1, "tab:blue")):
        se, te = exp[side]
        ax.plot(sgn * se, te, "o", color=col, ms=4, label="exp " + side, zorder=5)
    for rd, col, lw, lab in [(a.run_dir, "k", 2.0, "forge CHT")] + \
                            [(b, "0.55", 1.0, None) for b in a.band]:
        s, ss, Tw, _ = load(rd)
        for side, sgn in (("SS", +1), ("PS", -1)):
            sel = ss if side == "SS" else ~ss
            o = np.argsort(s[sel])
            ax.plot(sgn * s[sel][o], Tw[sel][o], "-", color=col, lw=lw,
                    label=(lab if side == "SS" else None), alpha=0.9)
    s, ss, Tw, f = load(a.run_dir)
    Te = np.array([np.interp(s[k], *exp["SS" if ss[k] else "PS"]) for k in range(len(s))])
    m = s <= 0.87; e = Tw[m] - Te[m]
    ax.axvspan(0, 0.25, color="0.9", zorder=0)
    ax.axvspan(0.87, 1.0, color="#ffe8e8", zorder=0); ax.axvspan(-1.0, -0.87, color="#ffe8e8", zorder=0)
    ax.set_xlabel("$-s/S$ (PS)   |   $+s/S$ (SS)"); ax.set_ylabel("$T_w$ [K]")
    ax.set_title(f"{t['vane']} {a.run} — coupled CHT wall temperature "
                 f"(bias {e.mean():+.1f} K, rms {np.sqrt((e**2).mean()):.1f} K over s/S<=0.87)",
                 fontsize=10)
    ax.grid(alpha=0.3); ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=9, frameon=False)
    out = a.out or str(Path(a.run_dir) / "Tw_compare.png")
    fig.tight_layout(); fig.savefig(out, dpi=110)
    print(f"[plot_tw] {f} -> {out}")
    # 領域別の壁温偏差 (compare_h と同じ区分。報告の表はここから転記する)
    for name, mm in (("PS", ~ss & (s <= 0.87)), ("SS laminar (s/S<0.25)", ss & (s < 0.25)),
                     ("SS post-transition", ss & (s >= 0.25) & (s <= 0.87)), ("all", s <= 0.87)):
        d = Tw[mm] - Te[mm]
        print(f"    {name:<24} n={mm.sum():3d}  bias {d.mean():+6.1f} K  rms {np.sqrt((d**2).mean()):5.1f} K  "
              f"max|err| {np.abs(d).max():5.1f} K")
    print(f"    Tw range {Tw.min():.1f}-{Tw.max():.1f} K, mean {Tw.mean():.1f} K")


if __name__ == "__main__":
    main()
