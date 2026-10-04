#!/usr/bin/env python3
r"""2D 計算結果の**標準コンタ図** (マッハ数・静圧・乱流粘性係数比・静温) を 1 枚に描く。

規約 (skill `forge-contour` が正本。ここはその実体化):

- **カラーマップは全量 `turbo` に統一する** (ParaView の Turbo と同じもの)。量ごとに色を変えない。
  色の違いが「量の違い」なのか「値の違い」なのかを読者に推測させないため。
- 周期翼列は `--pitch-cm` を与えると **1 ピッチの解を −1/0/+1 ピッチに並べて**描く
  (1 通路だけだと翼列として読めない)。
- $\mu_t/\mu$ は後流が境界層の 10 倍以上あるので **`--mut-clip` (既定 50) で頭打ち**にし、
  タイトルに明記する。線形の全域目盛では境界層が潰れ、対数では全面が一色になる。
- 他の量は 0.2–99.8 パーセンタイルで目盛を切る (1 点の外れ値で全体が潰れるのを防ぐ)。
- 凡例・カラーバーを**図の上に重ねない** (カラーバーは各パネルの下)。

`MESH/CONNE` は「**自分を含む語数**, 節点…」の可変長レコード (3 角 = 4 語、4 角 = 5 語)。
固定長だと思って reshape すると 4 角セルの開始位置から静かにずれる (2026-09-21 に実際に踏んだ)。

usage:
  python3 solver_density_cuda/tools/plot_field_contours.py <res_N.h5> --out fig.png \
      [--pitch-cm 11.773] [--ylim -1 24] [--title "..."] [--mut-clip 50]
"""
from __future__ import annotations

import argparse
import sys

import numpy as np
import h5py

CMAP = "turbo"          # **全量共通**。変えない (skill forge-contour)


def read(res):
    with h5py.File(res, "r") as f:
        V = f["VALUE"]
        X = np.array(f["MESH/COORD"]).reshape(-1, 3)[:, :2]
        a = np.array(f["MESH/CONNE"])
        d = {k: np.array(V[k]) for k in ("P", "T", "Ux", "Uy", "Uz", "sonic", "vis_lam", "vis_turb")
             if k in V}
    tri, i = [], 0
    while i < len(a):                       # 先頭は「自分を含む語数」
        n = int(a[i])
        idx = a[i + 1:i + n]
        if n == 4:
            tri.append(idx)
        elif n == 5:                        # 4 角は 2 枚の 3 角に割る
            tri.append(idx[[0, 1, 2]]); tri.append(idx[[0, 2, 3]])
        else:
            sys.exit(f"[plot_field_contours] 未対応のセル (語数 {n}) at {i}")
        i += n
    U = np.sqrt(d["Ux"] ** 2 + d["Uy"] ** 2 + d["Uz"] ** 2)
    F = {"Mach": U / d["sonic"], "P": d["P"] / 1000.0, "T": d["T"]}
    if "vis_turb" in d:
        F["mut"] = d["vis_turb"] / np.maximum(d["vis_lam"], 1e-30)
    return X, np.array(tri), F


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("res")
    ap.add_argument("--out", required=True)
    ap.add_argument("--pitch-cm", type=float, default=None, help="周期ピッチ [cm]。与えると 3 枚並べる")
    ap.add_argument("--ylim", type=float, nargs=2, default=None)
    ap.add_argument("--xlim", type=float, nargs=2, default=None)
    ap.add_argument("--title", default=None)
    ap.add_argument("--mut-clip", type=float, default=50.0)
    a = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.tri import Triangulation
    # 日本語タイトルが豆腐になるのを防ぐ (既定の DejaVu には CJK が無い)。
    # ~/.fonts の Noto Sans CJK JP を登録する。無ければ黙って既定のまま。
    try:
        from matplotlib import font_manager
        import pathlib as _pl
        for _f in _pl.Path.home().glob(".fonts/NotoSansCJK*"):
            font_manager.fontManager.addfont(str(_f))
        if any("Noto Sans CJK JP" == f.name for f in font_manager.fontManager.ttflist):
            plt.rcParams["font.family"] = "Noto Sans CJK JP"
            plt.rcParams["axes.unicode_minus"] = False
    except Exception:
        pass

    X, tri, F = read(a.res)
    Xc = X * 100.0
    if a.pitch_cm:
        n, reps = len(Xc), (-1, 0, 1)
        Xc = np.vstack([Xc + np.array([0.0, j * a.pitch_cm]) for j in reps])
        tri = np.vstack([tri + k * n for k in range(len(reps))])
        F = {k: np.tile(v, len(reps)) for k, v in F.items()}
    T = Triangulation(Xc[:, 0], Xc[:, 1], tri)

    fields = [("Mach", "Mach number"), ("P", "static pressure [kPa]")]
    if "mut" in F:
        fields.append(("mut", f"turbulent / laminar viscosity\n(clipped at {a.mut_clip:g}; "
                              f"max {F['mut'].max():.0f})"))
    fields.append(("T", "static temperature [K]"))

    fig, ax = plt.subplots(1, len(fields), figsize=(4.25 * len(fields), 6.2))
    for k_ax, (key, lab) in zip(np.atleast_1d(ax), fields):
        v = F[key]
        if key == "mut":
            vv = np.clip(v, 0.0, a.mut_clip)
            cs = k_ax.tricontourf(T, vv, levels=np.linspace(0, a.mut_clip, 41), cmap=CMAP, extend="max")
        else:
            lo, hi = np.percentile(v, 0.2), np.percentile(v, 99.8)
            vv = np.clip(v, lo, hi)
            cs = k_ax.tricontourf(T, vv, levels=40, cmap=CMAP)
        k_ax.tricontour(T, vv, levels=12, colors="k", linewidths=0.25, alpha=.35)
        cb = fig.colorbar(cs, ax=k_ax, fraction=0.055, pad=0.02, orientation="horizontal",
                          location="bottom")
        cb.ax.tick_params(labelsize=8)
        k_ax.set_aspect("equal"); k_ax.set_title(lab, fontsize=10); k_ax.tick_params(labelsize=8)
        if key == "Mach":
            k_ax.set_ylabel("y [cm]", fontsize=8)
        else:
            k_ax.set_yticklabels([])
        if a.ylim: k_ax.set_ylim(*a.ylim)
        if a.xlim: k_ax.set_xlim(*a.xlim)
    if a.title:
        fig.suptitle(a.title, fontsize=11)
    fig.tight_layout(rect=[0, 0.02, 1, 0.94])
    fig.savefig(a.out, dpi=110)
    print(f"[plot_field_contours] -> {a.out}  Mach max {F['Mach'].max():.3f}  "
          f"T {F['T'].min():.0f}-{F['T'].max():.0f} K"
          + (f"  mut/mu max {F['mut'].max():.0f}" if "mut" in F else ""))


if __name__ == "__main__":
    main()
