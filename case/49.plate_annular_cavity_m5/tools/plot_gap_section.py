#!/usr/bin/env python3
"""すきま断面の速度ベクトル (ユーザ指定 2026-09-19)。

2 種類を出す。

  上段: **水平断面** (z = 一定) の環状すきまを上から見た図。面内 (u_x, u_y) をベクトルで、
        面内速度の大きさを色で。押し出しメッシュの z 層をそのまま使う (補間しない)。
        半割計算なので y<0 側は y>0 を鏡像にして全周を描く。
  下段: **すきまを周方向に展開した面** (θ-z)。(u_θ, u_z) をベクトルで。
        キャビティ内のノードを (θ, z) で平均して作る (中央面の半径を決め打ちしないので
        偏心でも使える)。循環が周方向のどこで上がり・下がりしているかが読める。

usage: python3 tools/plot_gap_section.py RUN [--z -1,-5,-12,-22] [-o out.png]
"""
import argparse
import sys
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
CASE = HERE.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(CASE))
import geom_common as gc  # noqa: E402


def load(run):
    run = Path(run)
    man = gc.load_manifest(run=run)
    snaps = sorted(run.glob("res_[0-9]*.h5"), key=lambda p: int(p.stem.split("_")[1]))
    with h5py.File(snaps[-1], "r") as f:
        c = np.array(f["MESH/COORD"]).reshape(-1, 3).astype(np.float64)
        u = np.stack([np.array(f["VALUE/U" + k]).astype(np.float64) for k in "xyz"], 1)
    m = gc.cavity_mask(c[:, 0], c[:, 1], c[:, 2], man)
    return man, c[m], u[m], int(snaps[-1].stem.split("_")[1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--z", default="-1,-5,-12.5,-22", help="水平断面の深さ [mm]")
    ap.add_argument("-o", "--out", default="gap_section.png")
    a = ap.parse_args()
    man, c, u, step = load(a.run)
    G = man["geometry"]
    zw = [float(x) for x in a.z.split(",")]

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

    # ---- 共通座標: θ (0=上流, 180=下流) と すきま内の無次元半径 fr (0=内円柱, 1=外筒) ----
    r_all = np.hypot(c[:, 0], c[:, 1])
    erx, ery = c[:, 0] / np.maximum(r_all, 1e-12), c[:, 1] / np.maximum(r_all, 1e-12)
    ur_all = u[:, 0] * erx + u[:, 1] * ery                  # 半径方向 (外向き正)
    utg = -u[:, 0] * ery + u[:, 1] * erx                    # 幾何 CCW
    thc = np.arctan2(np.abs(c[:, 1]), -c[:, 0])             # 0 = 上流
    ut_all = -np.where(c[:, 1] >= 0, utg, -utg)             # e_θ = -e_φ。θ 増加(下流向き)を正
    ri_all = gc.inner_radius_at(np.clip(thc, 0.0, np.pi), man)
    fr_all = np.clip((r_all - ri_all) / np.maximum(G["Ro"] - ri_all, 1e-12), 0.0, 1.0)
    thd_all = np.degrees(thc)
    zl = np.unique(np.round(c[:, 2], 9))

    ntb = 40
    teb = np.linspace(0.0, 180.0, ntb + 1)
    # 半径方向は**分位**で切る (境界層で壁際に密集しており等間隔だと中央が空になる)
    feb = np.unique(np.quantile(fr_all, np.linspace(0.0, 1.0, 23)))
    feb[0], feb[-1] = 0.0, 1.0
    nfb = len(feb) - 1
    TB, FB = np.meshgrid(0.5 * (teb[1:] + teb[:-1]), 0.5 * (feb[1:] + feb[:-1]))

    def slice_maps(sel):
        it = np.clip(np.digitize(thd_all[sel], teb) - 1, 0, ntb - 1)
        jf = np.clip(np.digitize(fr_all[sel], feb) - 1, 0, nfb - 1)
        ib = it * nfb + jf
        cnt = np.bincount(ib, minlength=ntb * nfb).astype(float)
        def av(v):
            return (np.bincount(ib, weights=v, minlength=ntb * nfb)
                    / np.maximum(cnt, 1)).reshape(ntb, nfb).T
        UZ, UT, UR = av(u[sel, 2]), av(ut_all[sel]), av(ur_all[sel])
        empty = cnt.reshape(ntb, nfb).T == 0
        for m in (UZ, UT, UR):
            m[empty] = np.nan
        return UZ, UT, UR

    ks = [int(np.argmin(np.abs(zl - z0 * 1e-3))) for z0 in zw]
    maps = [slice_maps(np.abs(c[:, 2] - zl[k]) < 1e-9) for k in ks]
    vlim = float(np.nanpercentile(np.abs(np.concatenate([m[0].ravel() for m in maps])), 98))

    n = len(zw)
    fig = plt.figure(figsize=(4.6 * n, 13.0))
    gs = fig.add_gridspec(3, n, height_ratios=[1.25, 1.0, 0.85], hspace=0.38, wspace=0.28)

    for j, (k, (UZ, UT, UR)) in enumerate(zip(ks, maps)):
        b = fig.add_subplot(gs[0, j])
        pm = b.pcolormesh(teb, feb, np.clip(UZ, -vlim, vlim), cmap="RdBu_r",
                          vmin=-vlim, vmax=vlim, shading="flat")
        sp = np.maximum(np.hypot(UT, UR), 1e-9)
        b.quiver(TB[::2, ::3], FB[::2, ::3], (UT / sp)[::2, ::3], (UR / sp)[::2, ::3],
                 color="k", scale=22, width=0.008)
        b.set_title("z = %.2f mm" % (zl[k] * 1e3), fontsize=12)
        b.set_xlabel("θ [deg]   0=上流 / 180=下流")
        b.set_xticks([0, 90, 180])
        if j == 0:
            b.set_ylabel("すきま内の位置\n0 = 内円柱壁 / 1 = 外筒壁\n(壁際を引き伸ばした目盛)")
        else:
            b.set_yticklabels([])
        if j == n - 1:
            fig.colorbar(pm, ax=b, label="鉛直速度 u_z [m/s]  赤=上昇 / 青=下降",
                         fraction=0.06, pad=0.03)

    # ---- 中段: 展開図 (θ × z)。すきま幅方向に平均した「全体の循環」 ----
    b = fig.add_subplot(gs[1, :])
    nt = 46
    te = np.linspace(0, 180, nt + 1)
    zlv, iz = np.unique(np.round(c[:, 2] * 1e3, 6), return_inverse=True)
    nz = len(zlv)
    it = np.clip(np.digitize(thd_all, te) - 1, 0, nt - 1)
    ib = it * nz + iz
    cnt = np.bincount(ib, minlength=nt * nz).astype(float)
    def av2(v):
        return (np.bincount(ib, weights=v, minlength=nt * nz)
                / np.maximum(cnt, 1)).reshape(nt, nz).T
    UT, UZ = av2(ut_all), av2(u[:, 2])
    SP = np.hypot(UT, UZ)
    SP[cnt.reshape(nt, nz).T == 0] = np.nan
    T, Z = np.meshgrid(0.5 * (te[1:] + te[:-1]), zlv)
    lo, hi = np.nanpercentile(SP, [1, 99])
    cf = b.contourf(T, Z, np.clip(SP, lo, hi), levels=np.linspace(lo, max(hi, lo + 1e-9), 25),
                    cmap="viridis")
    fig.colorbar(cf, ax=b, label="すきま面内の速さ [m/s]", fraction=0.025, pad=0.01)
    sp = np.maximum(SP, 1e-9)
    sz = max(1, nz // 22)
    b.quiver(T[::sz, ::2], Z[::sz, ::2], (UT / sp)[::sz, ::2], (UZ / sp)[::sz, ::2],
             color="w", scale=30, width=0.0024)
    b.set_title("すきまを周方向に展開 (θ × 深さ) — **すきま幅方向に平均**した全体の循環。"
                "矢印は向きのみ", fontsize=12)
    b.set_xlabel("θ [deg]   (0 = 上流, 180 = 下流)")
    b.set_ylabel("z [mm]  (0 = 開口)")
    b.set_xticks([0, 45, 90, 135, 180])

    # ---- 下段: 開口面の正味出入り (θ ごと) ----
    import cavity_eval as ce
    from scipy.interpolate import LinearNDInterpolator
    cc, vv = ce.read(sorted(Path(a.run).glob("res_[0-9]*.h5"),
                            key=lambda q: int(q.stem.split("_")[1]))[-1])
    zf = man["eval"].get("flux_depth_m") or man["eval"]["flux_depth_frac"] * G["depth"]
    zla = np.unique(np.round(cc[:, 2], 9))
    kk = int(np.argmin(np.abs(zla + zf)))
    sl = np.abs(cc[:, 2] - zla[kk]) < 1e-9
    li = {nm: LinearNDInterpolator(cc[sl][:, :2], vv[nm][sl]) for nm in ("ro", "Uz")}
    nth2, nr2 = 181, 60
    tg = np.linspace(0.0, np.pi, nth2)
    ria = gc.inner_radius_at(tg, man)
    fra = (np.arange(nr2) + 0.5) / nr2
    RR = ria[:, None] + fra[None, :] * (G["Ro"] - ria[:, None])
    dA = (G["Ro"] - ria[:, None]) / nr2 * RR * (np.pi / (nth2 - 1))
    DX, DY = gc.ray_dir(np.repeat(tg[:, None], nr2, axis=1))
    pq = np.stack([RR * DX, RR * DY], axis=-1).reshape(-1, 2)
    fx = (li["ro"](pq) * li["Uz"](pq)).reshape(nth2, nr2) * dA
    td = np.degrees(tg)
    b = fig.add_subplot(gs[2, :])
    b.fill_between(td, 0, fx.sum(1) * 1e6, where=(fx.sum(1) > 0), color="tab:red", alpha=0.25)
    b.fill_between(td, 0, fx.sum(1) * 1e6, where=(fx.sum(1) <= 0), color="tab:blue", alpha=0.25)
    b.plot(td, fx.sum(1) * 1e6, "k-", lw=2.0, label="正味 (赤帯=出 / 青帯=入)")
    b.plot(td, np.where(fx > 0, fx, 0).sum(1) * 1e6, "--", color="tab:red", lw=1.2, label="出 (上昇) の合計")
    b.plot(td, -np.where(fx < 0, fx, 0).sum(1) * 1e6, "--", color="tab:blue", lw=1.2, label="入 (下降) の合計")
    b.axhline(0, color="0.4", lw=0.9)
    b.set_title("開口面 (z = %.2f mm) をどこで出入りしているか — θ ごとに半径方向へ積分"
                % (zla[kk] * 1e3), fontsize=12)
    b.set_xlabel("θ [deg]   (0 = 上流, 180 = 下流)")
    b.set_ylabel("質量流束 [1e-6 kg/s]")
    b.set_xticks([0, 45, 90, 135, 180])
    b.legend(fontsize=9, ncol=3)
    b.grid(alpha=0.3)

    fig.suptitle("case/49  すきま断面の流れ — %s  (step %d;  偏心 %.1f mm)\n"
                 "上段: 深さごとの断面 (色=上下方向の速度、矢印=面内の向き)  / "
                 "中段: 幅方向に平均した循環  / 下段: 開口の出入り"
                 % (Path(a.run).name, step, G["x_off"] * 1e3), fontsize=13)
    fig.savefig(a.out, dpi=130, bbox_inches="tight")
    print("wrote", a.out)


if __name__ == "__main__":
    main()
