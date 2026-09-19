#!/usr/bin/env python3
"""すきま内部の流れ (ベクトル) と壁熱流束 (コンター) を描く。

usage: python3 tools/plot_cavity_fields.py RUN [--out fields.png] [--step N]

出す図:
  (a)(b) 外筒壁 / 内円柱側面の **壁熱流束 q'' コンター** — 周方向 θ × 深さ z に展開
         (θ=0° が下流側, 180° が上流側。壁面ダンプのノード値をそのまま三角形分割して描く)
  (c)    すきま中央面の **速度ベクトル** (u_θ, u_z) + 速度大きさコンター (同じ展開図)
  (d)(e) θ=0° (下流) / θ=180° (上流) の **半径-深さ断面**: (u_r, u_z) ベクトル + 温度コンター
"""
import argparse
import sys
from pathlib import Path

import h5py
import numpy as np
from scipy.interpolate import LinearNDInterpolator

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager  # noqa: E402

HERE = Path(__file__).resolve().parent
CASE = HERE.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(CASE))
import geom_common as gc  # noqa: E402
from cavity_eval import read, snapshots, wall_dump  # noqa: E402
from setup import load as load_conditions  # noqa: E402

for _f in Path.home().joinpath(".fonts").glob("NotoSansCJKjp-Regular.otf"):
    font_manager.fontManager.addfont(str(_f))
    matplotlib.rcParams["font.family"] = font_manager.FontProperties(fname=str(_f)).get_name()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--out", default=None)
    ap.add_argument("--step", type=int, default=None)
    a = ap.parse_args()
    man = gc.load_manifest()
    D = load_conditions()
    G = man["geometry"]
    Tw, Ro, Ri, off, dep = D["wall_T"], G["Ro"], G["Ri"], G["x_off"], G["depth"]

    snaps = snapshots(a.run)
    res = snaps[-1] if a.step is None else [s for s in snaps if int(s.stem.split("_")[1]) == a.step][0]
    step = int(res.stem.split("_")[1])
    c, v = read(res)

    # キャビティ + 開口すぐ上のノードだけで線形補間器を作る
    m = gc.cavity_mask(c[:, 0], c[:, 1], c[:, 2], man) | (
        (c[:, 2] < 0.004) & (c[:, 2] > -1e-9) & (np.hypot(c[:, 0], c[:, 1]) < Ro * 1.3))
    pts = c[m]
    L = {k: LinearNDInterpolator(pts, v[k][m]) for k in ("Ux", "Uy", "Uz", "T")}

    fig = plt.figure(figsize=(16.5, 9.4))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 1.15], hspace=0.30, wspace=0.26)

    # ---------- (a)(b) 壁熱流束コンター ----------
    for i, (grp, ttl) in enumerate((("cav_outer", "外筒壁 (r=%.0f mm)" % (Ro * 1e3)),
                                    ("cyl_side", "内円柱側面"))):
        ax = fig.add_subplot(gs[0, i])
        w = wall_dump(a.run, man["phys_id"][grp], step, name=grp)
        if w is None or "qwall" not in w:
            ax.text(.5, .5, "壁ダンプ無し", ha="center"); continue
        xyz = w["xyz"]
        cx = off if grp == "cyl_side" else 0.0
        th = np.degrees(np.arctan2(xyz[:, 1], xyz[:, 0] - cx))
        z = xyz[:, 2] * 1e3
        q = -np.asarray(w["qwall"], float) * 1e-3          # 壁に入る側を正 [kW/m2]
        sel = (z < -1e-6) & (z > -dep * 1e3 + 1e-6)
        lv = np.linspace(0, max(np.percentile(q[sel], 99.5), 1e-6), 21)
        cf = ax.tricontourf(th[sel], z[sel], np.clip(q[sel], lv[0], lv[-1]), levels=lv, cmap="inferno")
        fig.colorbar(cf, ax=ax, label="q'' [kW/m²]")
        ax.set_xlabel("周方向 θ [deg]  (0=下流, 180=上流)")
        ax.set_ylabel("深さ z [mm]")
        ax.set_title("(%s) %s の壁熱流束" % ("ab"[i], ttl), fontsize=12, loc="left")
        ax.set_xlim(0, 180); ax.set_ylim(-dep * 1e3, 0)

    # ---------- (c) すきま中央面の速度ベクトル ----------
    ax = fig.add_subplot(gs[0, 2])
    nth, nz = 61, 60
    tg = np.radians(np.linspace(0, 180, nth))
    zg = np.linspace(-dep * 0.999, -dep * 0.002, nz)
    rc = gc.gap_center_radius(tg, man)
    TH, ZZ = np.meshgrid(tg, zg)
    RC = np.tile(rc, (nz, 1))
    P = np.stack([RC * np.cos(TH), RC * np.sin(TH), ZZ], axis=-1).reshape(-1, 3)
    ux, uy, uz = (L["Ux"](P), L["Uy"](P), L["Uz"](P))
    ut = (-ux * np.sin(TH.ravel()) + uy * np.cos(TH.ravel())).reshape(nz, nth)
    uz = uz.reshape(nz, nth)
    sp = np.hypot(ut, uz)
    good = np.isfinite(sp)
    sp = np.where(good, sp, 0.0); ut = np.where(good, ut, 0.0); uz = np.where(good, uz, 0.0)
    cf = ax.contourf(np.degrees(TH), ZZ * 1e3, sp, levels=20, cmap="viridis")
    fig.colorbar(cf, ax=ax, label="|u| [m/s]")
    s = 3
    ax.quiver(np.degrees(TH)[::s, ::s], (ZZ * 1e3)[::s, ::s], ut[::s, ::s], uz[::s, ::s],
              color="w", scale_units="width", scale=np.nanmax(sp) * 14 + 1e-9, width=.0035)
    ax.set_xlabel("周方向 θ [deg]"); ax.set_ylabel("深さ z [mm]")
    ax.set_title("(c) すきま中央面の流れ (u_θ, u_z)", fontsize=12, loc="left")

    # ---------- (d)(e) 半径-深さ断面 ----------
    for i, (thd, ttl) in enumerate(((0.0, "θ=0° (下流側)"), (180.0, "θ=180° (上流側)"))):
        ax = fig.add_subplot(gs[1, i])
        t = np.radians(thd)
        ri = gc.inner_radius_at(np.array([t]), man)[0]
        rg = np.linspace(ri + 1e-5, Ro - 1e-5, 44)
        zg2 = np.linspace(-dep * 0.999, -dep * 0.0005, 180)
        RG, ZG = np.meshgrid(rg, zg2)
        P2 = np.stack([RG * np.cos(t), RG * np.sin(t), ZG], axis=-1).reshape(-1, 3)
        T2 = L["T"](P2).reshape(ZG.shape)
        ur = (L["Ux"](P2) * np.cos(t) + L["Uy"](P2) * np.sin(t)).reshape(ZG.shape)
        uz2 = L["Uz"](P2).reshape(ZG.shape)
        T2 = np.where(np.isfinite(T2), T2, Tw)
        ur = np.where(np.isfinite(ur), ur, 0.0); uz2 = np.where(np.isfinite(uz2), uz2, 0.0)
        cf = ax.contourf((RG - ri) * 1e3, ZG * 1e3, T2, levels=24, cmap="coolwarm")
        fig.colorbar(cf, ax=ax, label="T [K]")
        sz, sr = 6, 3
        ax.quiver(((RG - ri) * 1e3)[::sz, ::sr], (ZG * 1e3)[::sz, ::sr],
                  ur[::sz, ::sr], uz2[::sz, ::sr], color="k",
                  scale_units="height", scale=np.nanmax(np.hypot(ur, uz2)) * 9 + 1e-9, width=.004)
        ax.set_xlabel("内円柱側からの距離 [mm]  (幅 %.2f mm)" % ((Ro - ri) * 1e3))
        ax.set_ylabel("深さ z [mm]")
        ax.set_title("(%s) %s の断面 (u_r, u_z) と温度" % ("de"[i], ttl), fontsize=12, loc="left")

    # ---------- 説明 ----------
    ax = fig.add_subplot(gs[1, 2]); ax.axis("off")
    wh_txt = []
    for grp in ("cav_outer", "cyl_side", "cav_floor"):
        w = wall_dump(a.run, man["phys_id"][grp], step, name=grp)
        if w and "qwall" in w:
            q = -np.asarray(w["qwall"], float)
            Q = float(np.sum(q * w["w"]))
            wh_txt.append("  %-10s  Q = %6.2f W (半割),  q'' 平均 %5.0f / 最大 %6.0f W/m²"
                          % (grp, Q, Q / w["area"], q.max()))
    ax.text(0, 1, "run: %s  (step %d)\nすきま %.2f mm × 深さ %.0f mm, 壁 %g K, M%.2g\n\n%s\n\n"
            "θ=0° が下流側 (流れは +x)、θ=180° が上流側。\n"
            "(d)(e) の横軸は内円柱表面からの距離。"
            % (a.run, step, G["gap_nom"] * 1e3, dep * 1e3, Tw, D["mach"], "\n".join(wh_txt)),
            va="top", ha="left", fontsize=10.5, transform=ax.transAxes)

    out = Path(a.out or (Path(a.run) / "cavity_fields.png"))
    fig.suptitle("case/49  すきま内部の流れと壁熱流束  —  %s" % a.run, fontsize=13)
    fig.savefig(out, dpi=125, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
