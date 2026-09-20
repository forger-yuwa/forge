#!/usr/bin/env python3
"""CFD (cavity_rear.csv) と W70 Fig 6(a) の実測 (ref/) を、両方の分母で並べて表と図にする。

usage: python3 tools/compare_ref.py RUN [--t0 run_0003_T0_A1_main] [--out compare.png]
"""
import argparse, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
CASE = HERE.parent


def jp_font():
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import font_manager, rcParams
    for p in Path.home().glob(".fonts/*CJK*"):
        try:
            font_manager.fontManager.addfont(str(p))
        except Exception:
            pass
    for name in ("Noto Sans CJK JP", "Noto Sans JP", "IPAGothic", "DejaVu Sans"):
        if any(f.name == name for f in font_manager.fontManager.ttflist):
            rcParams["font.family"] = name
            break
    rcParams["axes.unicode_minus"] = False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    rd = CASE / a.run
    ev = json.loads((rd / "cavity_eval.json").read_text(encoding="utf-8"))
    cfd = np.loadtxt(rd / "cavity_rear.csv", delimiter=",", skiprows=1)
    rows = [l.split(",") for l in (CASE / "ref/w70_fig6a_wd0063_rear.csv").read_text().splitlines()
            if l and not l.startswith("#") and not l.startswith("x_over_d")]
    xr = np.array([float(r[0]) for r in rows]); qr = np.array([float(r[1]) for r in rows])
    o = np.argsort(cfd[:, 0]); xc, qc = cfd[o, 0], cfd[o, 1]

    from burggraf import qs_over_qfp
    w, d = 1.270e-3, 20.32e-3
    print(f"run {a.run}  分母: 文献 {ev['q_lit']*1e-3:.2f} / forge smooth {ev['q_cfd']*1e-3:.2f} kW/m²")
    print(f"{'x/d':>7} {'実験':>8} {'理論(B4)':>10} {'CFD/文献':>10} {'実験/理論':>10} {'CFD/理論':>10}")
    for x, qe in zip(xr, qr):
        qi = float(np.interp(x, xc, qc)) / ev["q_lit"]
        qt = float(qs_over_qfp(x * d, w, d))
        print(f"{x:7.3f} {qe:8.3f} {qt:10.4f} {qi:10.5f} {qe/qt:10.2f} {qi/qt:10.3f}")
    print(f"\n開口面積平均 q̄_c/q_fp: CFD {ev['qbar_over_lit']:.3f} (文献分母) / "
          f"{ev['qbar_over_cfd']:.3f} (forge 分母)   vs 実験 1.07 (Fig 12 試読)")

    jp_font()
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7.2, 4.6), dpi=140)
    ax.semilogy(xr, qr, "o", ms=7, mfc="none", mew=1.6, color="#C4502A", label="実験 W70 Fig 6(a)")
    ax.semilogy(xc, qc / ev["q_lit"], "-", lw=2, color="#2E6F9E", label="forge (分母 = 文献 q_fp)")
    ax.semilogy(xc, qc / ev["q_cfd"], "--", lw=1.6, color="#2E6F9E", alpha=.7,
                label="forge (分母 = forge smooth)")
    xt = np.linspace(0.004, 1.0, 400)
    ax.semilogy(xt, qs_over_qfp(xt * d, w, d), "-.", lw=1.8, color="#4C7A34",
                label="Burggraf 非粘性コア理論 (B4)")
    ax.axhspan(1e-4, 0.02, color="0.5", alpha=.18)
    ax.text(0.62, 0.012, "計測限界帯 (q/q_fp ≲ 0.02)", fontsize=8.5, color="0.35")
    ax.axvline(1.270 / 20.32, color="0.4", ls=":", lw=1)
    ax.text(1.270 / 20.32 * 1.05, 0.5, "1 すきま幅", fontsize=8.5, color="0.35", rotation=90, va="top")
    ax.set_xlabel("x / d  (後壁を上端から下向き)")
    ax.set_ylabel("$q_s / q_{fp}$")
    ax.set_xlim(0, 0.6); ax.set_ylim(1e-4, 3)
    ax.grid(alpha=.25, which="both")
    ax.legend(fontsize=9, loc="upper right")
    ax.set_title("深キャビティ後壁の熱流束分布 (w/d = 0.063, M 6.9, 層流)", fontsize=11)
    fig.tight_layout()
    out = Path(a.out) if a.out else (rd / "compare_rear.png")
    fig.savefig(out)
    print(f"→ {out}")


if __name__ == "__main__":
    main()
