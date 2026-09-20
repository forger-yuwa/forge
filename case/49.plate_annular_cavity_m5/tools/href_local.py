#!/usr/bin/env python3
"""基準温度 4 定義の **深さ方向の局所** 比較 (plan §4.14.1)。

面積平均 (`href_defs.py`) は「開口近傍の外部流駆動」と「深部の循環駆動」を混ぜてしまう。
ここでは同じ 4 定義を**深さ z ごと**に出し、各深さで壁温 20/500/1000 degC に対する
h の振れ幅 (max/min - 1) を測る。振れ幅が小さい深さ帯 = その定義が局所的にも使える帯。

  (A) すきま中央面の局所 T0   … その深さの中央面 T0。`cavity_eval.wall_href` の `_dT_node`
  (B) キャビティ体積平均 T    … 全体で 1 つのスカラー
  (C) 主流総温 Tt_inf        … 全体で 1 つのスカラー
  (D) 断熱壁温 Taw           … 全体で 1 つのスカラー

usage: python3 tools/href_local.py RUN RUN RUN [--nz 30] [--plot]
"""
import argparse
import json
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
from href_defs import cavity_mean_T  # noqa: E402

WALLS = [("cav_outer", "外筒内壁"), ("cyl_side", "円柱側面")]
DEFS = ["(A) すきま中央面 局所T0", "(B) キャビティ体積平均 T", "(C) 主流総温 Tt", "(D) 断熱壁温 Taw"]
COLS = ["tab:blue", "tab:orange", "tab:red", "tab:green"]


def profile(run, nz):
    """1 run から、壁ごとに深さ binned の q'' と 4 定義の ΔT を返す。"""
    run = Path(run)
    man = gc.load_manifest(run=run)
    D = ce.run_conditions(run)
    snaps = ce.snapshots(run)
    step = int(snaps[-1].stem.split("_")[1])
    c, v = ce.read(snaps[-1])
    from total_quantities import total_state
    T0 = np.asarray(total_state(str(run), str(snaps[-1]))["T0"], float)
    wh = ce.wall_heat(run, step, man, D)
    wh = ce.wall_href(wh, man, D, c, v, T0)
    Tcav = cavity_mean_T(snaps[-1], man)[0]
    Tw = float(D["wall_T"])
    Taw = float(ce.taw_of(D))
    Tt = float(D["Tt_tp"] if D.get("gas_used") == "TP" else D["Tt_cpg"])
    # **等間隔 bin でなく、押し出しメッシュの z 層そのもの**で集計する。
    # 層厚は開口・床の近くで 8 um、中央で 2.2 mm と 2 桁以上変わるので、等間隔 bin を切ると
    # 中央に空 bin が並び「データが無い」ように見える (実体はただの層の粗さ)。
    out = dict(run=run.name, Tw=Tw, Tcav=Tcav, Taw=Taw, Tt=Tt, walls={})
    for g, _ in WALLS:
        if g not in wh:
            continue
        d = wh[g]
        zn, qn, wt = d["_z_node"], d["_qin_node"], d["_w_node"]
        zl, ib = np.unique(np.round(zn.astype(np.float64), 9), return_inverse=True)
        n = len(zl)
        den = np.bincount(ib, weights=wt, minlength=n)
        qpp = np.bincount(ib, weights=qn * wt, minlength=n) / np.maximum(den, 1e-30)
        dTA = np.bincount(ib, weights=d["_dT_node"] * wt, minlength=n) / np.maximum(den, 1e-30)
        out["walls"][g] = dict(z=zl, qpp=qpp, dT=[dTA,
                                                  np.full(n, Tcav - Tw),
                                                  np.full(n, Tt - Tw),
                                                  np.full(n, Taw - Tw)])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--nz", type=int, default=25)
    ap.add_argument("--plot", action="store_true")
    a = ap.parse_args()

    R = [profile(r, a.nz) for r in a.runs]

    for g, jp in WALLS:
        if g not in R[0]["walls"]:
            continue
        z = R[0]["walls"][g]["z"]
        if any(len(r["walls"][g]["z"]) != len(z) for r in R):
            raise SystemExit("run 間で壁の z 層数が違う (同一メッシュでのみ比較できる)")
        print("\n=== %s — 深さごとの h [W/m2K] と、壁温に対する振れ幅 ===" % jp)
        print("  z[mm]  " + "".join("%22s" % d.split(")")[0][1:] for d in DEFS))
        print("         " + "".join("%22s" % ("20/500/1000degC  振れ") for _ in DEFS))
        for i in range(0, len(z), max(1, len(z) // a.nz)):
            row = "%7.1f  " % (z[i] * 1e3)
            for k in range(len(DEFS)):
                hs = [r["walls"][g]["qpp"][i] / r["walls"][g]["dT"][k][i] for r in R]
                if any(not np.isfinite(h) or h <= 0 for h in hs):
                    row += "%22s" % "  -"
                else:
                    row += "%16s%6s" % ("/".join("%.0f" % h for h in hs),
                                        "%+.0f%%" % (100 * (max(hs) / min(hs) - 1)))
            print(row)

    if a.plot:
        plot(R)


def plot(R):
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
    ws = [(g, jp) for g, jp in WALLS if g in R[0]["walls"]]
    fig, ax = plt.subplots(2, len(ws) * 2, figsize=(7.2 * len(ws), 9.0),
                           squeeze=False)
    ls = ["-", "--", ":"]
    for j, (g, jp) in enumerate(ws):
        z = R[0]["walls"][g]["z"]
        # 左: h(z) を定義ごと (壁温は線種)
        b = ax[0][2 * j]
        for k, lab in enumerate(DEFS):
            for m, r in enumerate(R):
                h = r["walls"][g]["qpp"] / r["walls"][g]["dT"][k]
                b.plot(np.where(h > 0, h, np.nan), z * 1e3, ls[m], color=COLS[k],
                       label=lab if m == 0 else None, lw=1.6)
        b.set_xscale("log")
        b.set_title("%s — h(深さ)\n実線20℃ / 破線500℃ / 点線1000℃" % jp, fontsize=11)
        b.set_xlabel("h [W/m²K]")
        b.set_ylabel("z [mm] (0 = 開口)")
        b.legend(fontsize=8)
        b.grid(alpha=0.3, which="both")
        # 右: 振れ幅(z)
        b = ax[0][2 * j + 1]
        for k, lab in enumerate(DEFS):
            sp = []
            for i in range(len(z)):
                hs = [r["walls"][g]["qpp"][i] / r["walls"][g]["dT"][k][i] for r in R]
                sp.append(100 * (max(hs) / min(hs) - 1)
                          if all(np.isfinite(h) and h > 0 for h in hs) else np.nan)
            b.plot(sp, z * 1e3, "-", color=COLS[k], label=lab, lw=1.8)
        b.axvline(20, color="0.5", ls="--", lw=1.0)
        b.set_xscale("log")
        b.set_xlim(1, 1000)
        b.set_title("%s — 壁温 20→1000℃ での h の振れ幅\n(小さいほど外挿できる; 破線 20%%)"
                    % jp, fontsize=11)
        b.set_xlabel("max/min − 1 [%]")
        b.set_ylabel("z [mm]")
        b.grid(alpha=0.3, which="both")
        # 下段: q'' と 局所ΔT
        b = ax[1][2 * j]
        for m, r in enumerate(R):
            b.plot(r["walls"][g]["qpp"] / 1e3, z * 1e3, ls[m], color="k",
                   label="%.0f℃" % (r["Tw"] - 273.15), lw=1.6)
        b.set_xscale("log")
        b.set_title("%s — 壁面熱流束 q''" % jp, fontsize=11)
        b.set_xlabel("q'' [kW/m²]")
        b.set_ylabel("z [mm]")
        b.legend(fontsize=9)
        b.grid(alpha=0.3, which="both")
        b = ax[1][2 * j + 1]
        for m, r in enumerate(R):
            b.plot(r["walls"][g]["dT"][0], z * 1e3, ls[m], color="tab:blue",
                   label="(A) 局所 %.0f℃" % (r["Tw"] - 273.15), lw=1.6)
        b.set_title("%s — (A) の局所温度差 T0_mid − Tw" % jp, fontsize=11)
        b.set_xlabel("ΔT [K]")
        b.set_ylabel("z [mm]")
        b.legend(fontsize=8)
        b.grid(alpha=0.3)
    fig.suptitle("case/49 形状2 同心 M9 — 基準温度 4 定義の**深さ方向**の比較", fontsize=13)
    fig.tight_layout()
    out = Path("href_local.png")
    fig.savefig(out, dpi=130)
    print("\nwrote", out)


if __name__ == "__main__":
    main()
