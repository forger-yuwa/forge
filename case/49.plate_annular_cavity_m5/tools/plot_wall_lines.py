#!/usr/bin/env python3
"""壁面の熱流束・熱伝達率のライン図 (ユーザ指定 2026-09-19)。

  側壁 (外筒内壁 / 円柱側面) … **深さ方向** の分布 (縦軸 z)
  底面                        … **周方向** の分布 (横軸 θ; θ=0 が上流, 180 度が下流)

熱伝達率は 2 定義を重ねる:
  (A) すきま中央面の局所 T0 基準  … 局所で最良だが FEM の BC にはできない (plan §4.7.3)
  (D) 断熱壁温 Taw 基準           … FEM に渡すのはこちら

側壁は押し出しメッシュの z 層そのもので集計する (等間隔 bin は中央の層が粗く空 bin が並ぶ)。

usage: python3 tools/plot_wall_lines.py RUN [RUN ...] [-o out.png] [--title ...]
"""
import argparse
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
CASE = HERE.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(CASE))
sys.path.insert(0, str(CASE.parents[1] / "solver_density_cuda" / "tools"))
import geom_common as gc  # noqa: E402
import cavity_eval as ce  # noqa: E402

SIDE = [("cav_outer", "外筒内壁"), ("cyl_side", "円柱側面")]
LS = ["-", "--", ":"]


def gather(run):
    run = Path(run)
    man = gc.load_manifest(run=run)
    D = ce.run_conditions(run)
    snaps = ce.snapshots(run)
    step = int(snaps[-1].stem.split("_")[1])
    c, v = ce.read(snaps[-1])
    from total_quantities import total_state
    T0 = np.asarray(total_state(str(run), str(snaps[-1]))["T0"], float)
    wh = ce.wall_href(ce.wall_heat(run, step, man, D), man, D, c, v, T0)
    Tw, Taw = float(D["wall_T"]), float(ce.taw_of(D))
    out = dict(run=run.name, Tw=Tw, Taw=Taw, side={}, floor=None)
    for g, _ in SIDE:
        if g not in wh:
            continue
        d = wh[g]
        zl, ib = np.unique(np.round(d["_z_node"].astype(np.float64), 9), return_inverse=True)
        den = np.bincount(ib, weights=d["_w_node"], minlength=len(zl))
        def av(a):
            return np.bincount(ib, weights=a * d["_w_node"], minlength=len(zl)) / np.maximum(den, 1e-30)
        q = av(d["_qin_node"])
        out["side"][g] = dict(z=zl, q=q, hA=av(np.nan_to_num(d["_href_node"], nan=0.0)),
                              hD=q / max(Taw - Tw, 1e-30))
    if "cav_floor" in wh:
        d = wh["cav_floor"]
        # θ=0 を上流 (-x 方向) に取る。半割なので 0..180 度
        th = np.degrees(np.arctan2(np.abs(d["_y_node"]), -d["_x_node"]))
        nb = 60
        edges = np.linspace(0.0, 180.0, nb + 1)
        ib = np.clip(np.digitize(th, edges) - 1, 0, nb - 1)
        den = np.bincount(ib, weights=d["_w_node"], minlength=nb)
        def av(a):
            return np.bincount(ib, weights=a * d["_w_node"], minlength=nb) / np.maximum(den, 1e-30)
        q = av(d["_qin_node"])
        q[den <= 0] = np.nan
        out["floor"] = dict(th=0.5 * (edges[1:] + edges[:-1]), q=q,
                            hA=av(np.nan_to_num(d["_href_node"], nan=0.0)),
                            hD=q / max(Taw - Tw, 1e-30))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    ap.add_argument("-o", "--out", default="wall_lines.png")
    ap.add_argument("--title", default="")
    a = ap.parse_args()
    R = [gather(r) for r in a.runs]

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    for f in Path.home().glob(".fonts/NotoSansCJK*"):
        try:
            font_manager.fontManager.addfont(str(f))
        except Exception:
            pass
    plt.rcParams["font.family"] = "Noto Sans CJK JP"

    fig, ax = plt.subplots(2, 3, figsize=(16.5, 9.4))
    for j, (g, jp) in enumerate(SIDE):
        b = ax[0][j]
        for m, r in enumerate(R):
            if g not in r["side"]:
                continue
            s = r["side"][g]
            b.plot(s["q"] / 1e3, s["z"] * 1e3, LS[m % 3], color="k", lw=1.7,
                   label="%.0f℃" % (r["Tw"] - 273.15))
        b.set_xscale("log")
        b.set_title("%s — 壁面熱流束 q''" % jp, fontsize=12)
        b.set_xlabel("q'' [kW/m²]")
        b.set_ylabel("z [mm]  (0 = 開口, 負が深い)")
        b.legend(fontsize=9)
        b.grid(alpha=0.3, which="both")

        b = ax[1][j]
        for m, r in enumerate(R):
            if g not in r["side"]:
                continue
            s = r["side"][g]
            b.plot(np.where(s["hA"] > 0, s["hA"], np.nan), s["z"] * 1e3, LS[m % 3],
                   color="tab:blue", lw=1.7,
                   label="(A) 局所T0 %.0f℃" % (r["Tw"] - 273.15))
            b.plot(np.where(s["hD"] > 0, s["hD"], np.nan), s["z"] * 1e3, LS[m % 3],
                   color="tab:green", lw=1.7,
                   label="(D) Taw %.0f℃" % (r["Tw"] - 273.15))
        b.set_xscale("log")
        b.set_title("%s — 熱伝達率 h" % jp, fontsize=12)
        b.set_xlabel("h [W/m²K]")
        b.set_ylabel("z [mm]")
        b.legend(fontsize=7.5, ncol=2)
        b.grid(alpha=0.3, which="both")

    b = ax[0][2]
    for m, r in enumerate(R):
        if r["floor"] is None:
            continue
        b.plot(r["floor"]["th"], r["floor"]["q"] / 1e3, LS[m % 3], color="k", lw=1.7,
               label="%.0f℃" % (r["Tw"] - 273.15))
    b.set_yscale("log")
    b.set_title("底面 — 壁面熱流束 q'' (周方向)", fontsize=12)
    b.set_xlabel("θ [deg]   (0 = 上流, 180 = 下流)")
    b.set_ylabel("q'' [kW/m²]")
    b.set_xticks([0, 45, 90, 135, 180])
    b.legend(fontsize=9)
    b.grid(alpha=0.3, which="both")

    b = ax[1][2]
    for m, r in enumerate(R):
        if r["floor"] is None:
            continue
        f = r["floor"]
        b.plot(f["th"], np.where(f["hA"] > 0, f["hA"], np.nan), LS[m % 3],
               color="tab:blue", lw=1.7, label="(A) 局所T0 %.0f℃" % (r["Tw"] - 273.15))
        b.plot(f["th"], np.where(f["hD"] > 0, f["hD"], np.nan), LS[m % 3],
               color="tab:green", lw=1.7, label="(D) Taw %.0f℃" % (r["Tw"] - 273.15))
    b.set_yscale("log")
    b.set_title("底面 — 熱伝達率 h (周方向)", fontsize=12)
    b.set_xlabel("θ [deg]   (0 = 上流, 180 = 下流)")
    b.set_ylabel("h [W/m²K]")
    b.set_xticks([0, 45, 90, 135, 180])
    b.legend(fontsize=7.5, ncol=2)
    b.grid(alpha=0.3, which="both")

    fig.suptitle(a.title or "case/49 — 壁面の熱流束・熱伝達率", fontsize=13)
    fig.tight_layout()
    fig.savefig(a.out, dpi=130)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
