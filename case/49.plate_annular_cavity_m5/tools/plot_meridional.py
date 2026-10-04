#!/usr/bin/env python3
"""対称面 (y=0) の子午断面をキャビティ**と上方の外部流れまで**まとめて描く。

`plot_cavity_fields.py` はキャビティ内部 (すきま中央面・θ 断面) しか描かないので、
開口より上で何が起きているか (境界層・開口上の剪断層・外部流の偏向) が見えない。
本ツールは z = -depth から z_top までを 1 枚に収める。

usage:
  python3 tools/plot_meridional.py <run_dir> [--step N] [--out FILE] [--zoom]

描くもの (半割モデルなので y=0 面の節点をそのまま使う):
  (a) 速度の大きさ + 流線 (キャビティ + 外部)
  (b) 温度
  (c) マッハ数 (外部の衝撃・膨張が見える)
  (d) 開口近傍のズーム (速度ベクトル)
"""
import argparse
import sys
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import geom_common as gc                      # noqa: E402

import matplotlib                             # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt               # noqa: E402
from matplotlib import font_manager           # noqa: E402
from scipy.interpolate import griddata        # noqa: E402

for fp in Path.home().joinpath(".fonts").glob("NotoSansCJKjp-Regular.otf"):
    font_manager.fontManager.addfont(str(fp))
    matplotlib.rcParams["font.family"] = font_manager.FontProperties(fname=str(fp)).get_name()


def snapshots(run):
    return sorted(Path(run).glob("res_[0-9]*.h5"), key=lambda f: int(f.stem.split("_")[1]))


def fluid_mask(X, Z, G):
    """対称面 (y=0) で**流体**の領域か。固体を塗らないために要る。

    z>0 は全て流体 (平板の上)。z<=0 は外筒 (|x|<Ro) の内側かつ内円柱 (|x-x_off|<Ri) の外側、
    深さ `depth` まで。これを外すと内円柱や平板の下まで補間で塗られ、図が読めなくなる。
    """
    Ro, Ri, off, dep = (G["Ro"] * 1e3, G["Ri"] * 1e3, G["x_off"] * 1e3, G["depth"] * 1e3)
    above = Z > 0.0
    slit = (np.abs(X) < Ro) & (np.abs(X - off) > Ri) & (Z > -dep)
    return above | slit


def grid(x, z, v, xlim, zlim, G, n=420, aspect=True):
    xi = np.linspace(*xlim, n)
    nz = max(40, int(n * (zlim[1] - zlim[0]) / (xlim[1] - xlim[0]))) if aspect else n
    zi = np.linspace(*zlim, nz)
    X, Z = np.meshgrid(xi, zi)
    V = griddata(np.stack([x, z], 1), v, (X, Z), method="linear")
    V = np.where(fluid_mask(X, Z, G), V, np.nan)      # **固体は塗らない**
    return X, Z, V


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--step", type=int, default=None)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    man = gc.load_manifest(run=a.run)
    G = man["geometry"]
    snaps = snapshots(a.run)
    if not snaps:
        raise SystemExit("res_*.h5 が無い: %s" % a.run)
    res = snaps[-1] if a.step is None else [s for s in snaps if int(s.stem.split("_")[1]) == a.step][0]
    step = int(res.stem.split("_")[1])
    with h5py.File(res, "r") as f:
        c = np.array(f["MESH/COORD"]).reshape(-1, 3)
        V = {k: np.array(f["VALUE/" + k]) for k in ("Ux", "Uy", "Uz", "T", "P", "ro")
             if "VALUE/" + k in f}
    # **対称面 (y=0) の節点**。半割モデルなので節点がちょうど乗っている。
    m = np.abs(c[:, 1]) < 1e-9
    if m.sum() < 500:
        m = np.abs(c[:, 1]) < 0.2e-3
    x, z = c[m, 0] * 1e3, c[m, 2] * 1e3
    ux, uz = V["Ux"][m], V["Uz"][m]
    umag = np.sqrt(ux ** 2 + V["Uy"][m] ** 2 + uz ** 2)
    T = V["T"][m]
    gam, cp = 1.4, 1004.5
    mach = umag / np.sqrt(np.maximum((gam - 1.0) * cp * T, 1e-9))

    Ro, dep, ztop = G["Ro"] * 1e3, G["depth"] * 1e3, G["z_top"] * 1e3
    xlim = (-2.2 * Ro, 2.6 * Ro)
    zlim = (-dep * 1.02, min(ztop, 2.0 * dep + 10.0))

    off = G["x_off"] * 1e3
    Ri = G["Ri"] * 1e3
    fig = plt.figure(figsize=(16.2, 9.6))
    gs = fig.add_gridspec(2, 4, width_ratios=[3.0, 3.0, 1.05, 1.05], hspace=0.30, wspace=0.55)
    ax_u, ax_T = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    ax_M, ax_z = fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])
    ax_up, ax_dn = fig.add_subplot(gs[:, 2]), fig.add_subplot(gs[:, 3])

    lims = {}
    for ttl, val, cmap, b in (("速度の大きさ |u| [m/s]", umag, "viridis", ax_u),
                              ("温度 T [K]", T, "inferno", ax_T),
                              ("マッハ数", mach, "coolwarm", ax_M)):
        X, Z, Vg = grid(x, z, val, xlim, zlim, G)
        f = np.isfinite(Vg)
        lo, hi = np.nanpercentile(Vg[f], [0.5, 99.5])
        if not (hi > lo):
            hi = lo + max(abs(lo), 1.0) * 1e-6
        cf = b.contourf(X, Z, np.clip(Vg, lo, hi), levels=np.linspace(lo, hi, 25), cmap=cmap)
        fig.colorbar(cf, ax=b, label=ttl, fraction=0.030, pad=0.02)
        lims[b] = (xlim, zlim)
        b.set_title(ttl, fontsize=11.5, loc="left")

    # 開口近傍 (キャビティ + 上方) を拡大し、ベクトルを重ねる
    xz = (-1.5 * Ro, 1.5 * Ro)
    zz = (-dep * 1.02, 0.7 * dep)
    X, Z, U = grid(x, z, ux, xz, zz, G, n=300)
    _, _, W = grid(x, z, uz, xz, zz, G, n=300)
    _, _, M = grid(x, z, umag, xz, zz, G, n=300)
    f = np.isfinite(M)
    lo, hi = np.nanpercentile(M[f], [0.5, 99.0])
    cf = ax_z.contourf(X, Z, np.clip(M, lo, hi), levels=np.linspace(lo, hi, 25), cmap="viridis")
    fig.colorbar(cf, ax=ax_z, label="|u| [m/s]", fraction=0.030, pad=0.02)
    lims[ax_z] = (xz, zz)
    st = 9
    Us, Ws = U[::st, ::st], W[::st, ::st]
    sp = np.hypot(Us, Ws)
    ax_z.quiver(X[::st, ::st], Z[::st, ::st], Us / np.maximum(sp, 1e-9), Ws / np.maximum(sp, 1e-9),
                color="w", scale=34, width=0.003)
    ax_z.set_title("開口近傍 (キャビティ + 上方)", fontsize=11.5, loc="left")

    # 上流側 / 下流側のすきまを縦長に拡大 (すきまは数 mm しかないので専用パネルが要る)
    gapmm = Ro - Ri
    for b, sgn, lab in ((ax_up, -1, "上流側すきま"), (ax_dn, +1, "下流側すきま")):
        if sgn < 0:
            xs = (-Ro - 0.35 * gapmm, -(Ri - off) + 0.35 * gapmm)
        else:
            xs = ((Ri + off) - 0.35 * gapmm, Ro + 0.35 * gapmm)
        zs = (-dep * 1.02, 0.5 * dep)
        X2, Z2, U2 = grid(x, z, ux, xs, zs, G, n=90, aspect=False)
        _, _, W2 = grid(x, z, uz, xs, zs, G, n=90, aspect=False)
        _, _, M2 = grid(x, z, umag, xs, zs, G, n=90, aspect=False)
        # 色は**キャビティ内部 (z<0)** で張る。外部流で張ると内部が真っ黒で読めない
        inside = np.isfinite(M2) & (Z2 < 0.0)
        f2 = inside if inside.any() else np.isfinite(M2)
        if f2.any():
            lo2, hi2 = np.nanpercentile(M2[f2], [1.0, 99.0])
            if not (hi2 > lo2):
                hi2 = lo2 + max(abs(lo2), 1.0) * 1e-6
            cf = b.contourf(X2, Z2, np.clip(M2, lo2, hi2), levels=np.linspace(lo2, hi2, 21),
                            cmap="viridis")
            fig.colorbar(cf, ax=b, label="|u| [m/s]", fraction=0.075, pad=0.04)
            s2 = 5
            Ua, Wa = U2[::s2, ::s2], W2[::s2, ::s2]
            spa = np.hypot(Ua, Wa)
            b.quiver(X2[::s2, ::s2], Z2[::s2, ::s2], Ua / np.maximum(spa, 1e-9),
                     Wa / np.maximum(spa, 1e-9), color="w", scale=26, width=0.012)
        b.set_title(lab, fontsize=11, loc="left")
        lims[b] = (xs, zs)

    for b in (ax_u, ax_T, ax_M, ax_z):
        b.set_aspect("equal")
    for b in (ax_u, ax_T, ax_M, ax_z, ax_up, ax_dn):
        b.set_xlabel("x [mm]  (流れは +x)")
        b.set_ylabel("z [mm]  (0 = 平板面)")
        b.axhline(0.0, color="0.35", lw=0.8, ls="--")
        for xx in (-Ro, Ro):
            b.plot([xx, xx], [-dep, 0], color="0.35", lw=0.9)
        for xx in (off - Ri, off + Ri):
            b.plot([xx, xx], [-dep, 0], color="0.55", lw=0.9)
        b.plot([-Ro, Ro], [-dep, -dep], color="0.35", lw=0.9)
        if b in lims:                       # ガイド線で軸が広がるので最後に固定する
            b.set_xlim(*lims[b][0]); b.set_ylim(*lims[b][1])

    gap = (G["Ro"] - G["Ri"]) * 1e3
    twall = float("nan")
    cj = Path(a.run) / "conditions.json"
    if cj.exists():
        import json
        twall = float(json.loads(cj.read_text()).get("wall_T", float("nan")))
    fig.suptitle("case/49  対称面 (y=0) の子午断面 — %s  (step %d;  外径 %.0f / 深さ %.0f / すきま %.1f mm"
                 ", 偏心 %.1f mm, 壁温 %.1f K = %.0f °C)"
                 % (a.run, step, 2 * Ro, dep, gap, off, twall, twall - 273.15),
                 fontsize=12.5)
    out = a.out or str(Path(a.run) / "meridional.png")
    fig.savefig(out, dpi=125, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
