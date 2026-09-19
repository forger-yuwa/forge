#!/usr/bin/env python3
"""すきま内部の流れ (ベクトル) と壁熱流束 (コンター) を描く。

usage: python3 tools/plot_cavity_fields.py RUN [--out fields.png] [--step N]

出す図:
  (a)(b) 外筒壁 / 内円柱側面の **壁熱流束 q'' コンター** — 周方向 θ × 深さ z に展開
         (**θ=0° が上流**, 180° が下流。壁面ダンプのノード値をそのまま三角形分割して描く)
  (c)    すきま中央面の **速度ベクトル** (u_θ, u_z) + 速度大きさコンター (同じ展開図)
  (d)(e) θ=0° (上流) / θ=180° (下流) の **半径-深さ断面**: (u_r, u_z) ベクトル + 温度コンター
ベクトルは**長さを正規化して向きだけを示し**、速度の大きさは背景コンター (対数) で出す
(すきま内は速度が 3 桁以上変わるので、長さで表すと深部が見えない)。
"""
import argparse
import sys
from pathlib import Path

import h5py
import numpy as np
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import Delaunay

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager  # noqa: E402

HERE = Path(__file__).resolve().parent
CASE = HERE.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(CASE))
import geom_common as gc  # noqa: E402
from cavity_eval import read, snapshots, wall_dump, wall_heat, wall_href  # noqa: E402
from setup import load as load_conditions  # noqa: E402

for _f in Path.home().joinpath(".fonts").glob("NotoSansCJKjp-Regular.otf"):
    font_manager.fontManager.addfont(str(_f))
    matplotlib.rcParams["font.family"] = font_manager.FontProperties(fname=str(_f)).get_name()


def grid_mean(u, v, val, w, nu=120, nv=48, ulim=None, vlim=None):
    r"""散布値を (u,v) の規則格子に面積重み平均する (tricontourf の縞アーチファクト回避)。

    **2 種類の「値が無い」を区別する** (2026-09-19, 外とう壁の h_ref が階段状に見えた件):

    - **面が無いビン** (細い面の端で、そのビンに節点が 1 つも落ちない) → 近傍で埋めてよい。
    - **値が定義できないビン** (節点はあるが値が NaN。h_ref は $\Delta T\le$ 閾値の領域で
      定義できない) → **埋めてはいけない**。埋めると「定義できない領域」が近傍の値で
      塗られ、閾値の等高線がそのまま**階段状の不連続**として現れる。空白のまま返す。

    旧実装は `val * w` に NaN が 1 つでも入ったビンを NaN にし、両者をまとめて近傍で
    埋めていたため、$\Delta T=1$ K の等高線が段々の境界として描かれていた。
    """
    ulo, uhi = ulim if ulim else (np.min(u), np.max(u))
    vlo, vhi = vlim if vlim else (np.min(v), np.max(v))
    ue = np.linspace(ulo, uhi, nu + 1)
    ve = np.linspace(vlo, vhi, nv + 1)
    iu = np.clip(np.digitize(u, ue) - 1, 0, nu - 1)
    iv = np.clip(np.digitize(v, ve) - 1, 0, nv - 1)
    k = iv * nu + iu
    val = np.asarray(val, float)
    fin = np.isfinite(val)
    num = np.bincount(k[fin], weights=(val[fin] * w[fin]), minlength=nu * nv).reshape(nv, nu)
    den = np.bincount(k[fin], weights=w[fin], minlength=nu * nv).reshape(nv, nu)
    cov = np.bincount(k, weights=w, minlength=nu * nv).reshape(nv, nu)   # 面が有るか
    g = np.where(den > 0, num / np.maximum(den, 1e-30), np.nan)
    empty = (cov <= 0)                       # 面が無い = 埋めてよい
    if empty.any() and np.isfinite(g).any():
        from scipy.ndimage import distance_transform_edt
        src = np.where(np.isfinite(g), g, np.nan)
        idx = distance_transform_edt(~np.isfinite(src), return_distances=False,
                                     return_indices=True)
        filled = src[tuple(idx)]
        g = np.where(empty, filled, g)       # **未定義ビン (cov>0, den=0) は NaN のまま**
    UC = 0.5 * (ue[1:] + ue[:-1]); VC = 0.5 * (ve[1:] + ve[:-1])
    # numpy 2 系では meshgrid が tuple を返すので list 化してから連結する
    return list(np.meshgrid(UC, VC)) + [g]


def wall_uv(d, man, grp):
    """壁ノードを (θ[deg], 面内座標) に展開。側壁は深さ z[mm]、床はすきま横断距離[mm]。"""
    x, y, z = d["_x_node"], d["_y_node"], d["_z_node"]
    cx = man["geometry"]["x_off"] if grp == "cyl_side" else 0.0
    th = np.degrees(np.arctan2(y, -(x - cx)))
    if grp == "cav_floor":
        thg = np.degrees(np.arctan2(y, -x))
        ri = gc.inner_radius_at(np.radians(thg), man)
        return thg, (np.hypot(x, y) - ri) * 1e3, "内円柱側からの距離 [mm]"
    return th, z * 1e3, "深さ z [mm]"


def plot_htc(run, step, man, D, wh, out):
    """**熱伝達率 h_ref のコンター** (キャビティ 3 壁) + 基準温度のコンター。"""
    G = man["geometry"]
    grps = [g for g in ("cav_outer", "cyl_side", "cav_floor")
            if g in wh and "_href_node" in wh[g]]
    if not grps:
        print("  h_ref が無い -> 熱伝達率コンターはスキップ")
        return
    fig, ax = plt.subplots(2, len(grps), figsize=(5.4 * len(grps), 8.4), squeeze=False)
    names = {"cav_outer": "外筒壁 (r=%.0f mm)" % (G["Ro"] * 1e3),
             "cyl_side": "内円柱側面", "cav_floor": "底面"}
    for i, g in enumerate(grps):
        d = wh[g]
        u, v, vlab = wall_uv(d, man, g)
        for r, (val, lab, cmap) in enumerate(((d["_href_node"], "h_ref [W/m²K]", "viridis"),
                                              (d["_Tref_node"], "基準温度 T0_ref [K]", "coolwarm"))):
            U, V, Z = grid_mean(u, v, val, d["_w_node"], ulim=(0, 180))
            # **上限はリップ帯を除いた領域から決める** (熱流束の図と同じ理由)。リップは
            # 幾何的特異点で q'' が発散するので、そこを含めて正規化すると深部が真っ黒になる。
            lipmm = man["eval"].get("lip_band_m", 1.0e-3) * 1e3
            deep = np.isfinite(Z)
            if g in ("cav_outer", "cyl_side"):
                deep = deep & (V < -lipmm)
            lo = np.nanpercentile(Z[np.isfinite(Z)], 0.5) if r else 0.0
            hi = (np.nanpercentile(Z[deep], 99.0) if deep.any()
                  else np.nanpercentile(Z[np.isfinite(Z)], 99.5))
            zmax = float(np.nanmax(Z[np.isfinite(Z)])) if np.isfinite(Z).any() else float("nan")
            # **未定義ビンは塗らない** (灰色のハッチで示す)。h_ref は dT_ref <= 閾値の領域で
            # 定義できず、そこを近傍で埋めると閾値の等高線が階段状の不連続として現れる。
            und = ~np.isfinite(Z)
            ax[r][i].set_facecolor("0.82")
            cf = ax[r][i].contourf(U, V, np.clip(Z, lo, hi), levels=np.linspace(lo, hi, 21),
                                   cmap=cmap)
            fig.colorbar(cf, ax=ax[r][i], label=lab)
            if r == 0:
                frac = 100.0 * und.sum() / und.size
                note = "上限 %.4g で飽和 (最大 %.4g)" % (hi, zmax)
                if und.any():
                    note += "\n灰色 = 係数が定義できない領域 (dT_ref ≤ %.3g K) %.0f %%" % (
                        d.get("href_dT_min_K", 0.05), frac)
                ax[r][i].text(0.985, 0.02, note, transform=ax[r][i].transAxes,
                              fontsize=8.5, color="w", ha="right", va="bottom")
            ax[r][i].set_xlabel("周方向 θ [deg]  (0=上流, 180=下流)")
            ax[r][i].set_ylabel(vlab)
            ax[r][i].set_title("(%s) %s の %s" % ("abcdef"[r * len(grps) + i], names[g],
                                                  "熱伝達率" if r == 0 else "基準温度"),
                               fontsize=12, loc="left")
    hm = sum(wh[g]["h_ref"] * wh[g]["area_m2"] for g in grps) / sum(wh[g]["area_m2"] for g in grps)
    fig.suptitle("case/49  キャビティ壁の熱伝達率 h_ref = q''/(T0_ref − T_w) と基準温度  —  %s"
                 "  (面積平均 h_ref = %.1f W/m²K)" % (run, hm), fontsize=12.5)
    fig.tight_layout()
    fig.savefig(out, dpi=125, bbox_inches="tight")
    print("wrote", out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--out", default=None)
    ap.add_argument("--max-points", type=int, default=120000,
                    help="場の線形補間に使うノード数の上限 (3D Delaunay が重いので間引く)")
    ap.add_argument("--step", type=int, default=None)
    a = ap.parse_args()
    man = gc.load_manifest(run=a.run)
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
    idx = np.flatnonzero(m)
    # **3D Delaunay は点数に対して急激に重くなる** (1.00M 節点の run で RSS 4.3 GB・10 分超に
    # なり、同時に走っていた 2.74M 節点の計算を OOM で落としかけた)。図のサンプル格子は
    # 181x200 / 44x180 しかないので、雲を間引いても絵は変わらない。
    if idx.size > a.max_points:
        rng = np.random.default_rng(0)                    # 図の再現性のため固定 seed
        idx = np.sort(rng.choice(idx, a.max_points, replace=False))
        print("  補間点を %d -> %d に間引き (--max-points)" % (int(m.sum()), idx.size), flush=True)
    pts = c[idx]
    # **三角形分割は 1 回だけ作って 4 変数で共有する** (変数ごとに張ると 4 倍かかる)
    tri = Delaunay(pts)
    L = {k: LinearNDInterpolator(tri, v[k][idx]) for k in ("Ux", "Uy", "Uz", "T")}

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
        th = np.degrees(np.arctan2(xyz[:, 1], -(xyz[:, 0] - cx)))   # θ=0 が上流
        z = xyz[:, 2] * 1e3
        q = -np.asarray(w["qwall"], float) * 1e-3          # 壁に入る側を正 [kW/m2]
        sel = (z < -1e-6) & (z > -dep * 1e3 + 1e-6)
        # **カラースケールはリップ帯を除いた領域から決める**。開口リップは 90 度の鋭角で
        # q'' が h^-1/2 で発散する幾何的特異点なので (plan §4.4.2)、そこを含めて正規化すると
        # 図が真っ黒になり、深部の分布がまったく見えない (実測: 最大 418 kW/m2 対 平均 3 kW/m2)。
        # リップ帯は飽和させ、上限をタイトルに書く。
        lipmm = man["eval"].get("lip_band_m", 1.0e-3) * 1e3
        deep = sel & (z < -lipmm)
        hi = max(np.percentile(q[deep], 99.0), 1e-6) if deep.any() else 1e-6
        lv = np.linspace(0, hi, 21)
        qmax_all = float(np.max(q[sel])) if sel.any() else 0.0
        cf = ax.tricontourf(th[sel], z[sel], np.clip(q[sel], lv[0], lv[-1]), levels=lv, cmap="inferno")
        fig.colorbar(cf, ax=ax, label="q'' [kW/m²]")
        ax.set_xlabel("周方向 θ [deg]  (0=上流, 180=下流)")
        ax.set_ylabel("深さ z [mm]")
        ax.set_title("(%s) %s の壁熱流束" % ("ab"[i], ttl), fontsize=12, loc="left")
        # 飽和の断り書きは軸内に小さく置く (タイトルに入れると隣のパネルと重なる)
        ax.text(0.015, 0.015, "上限 %.3g kW/m² で飽和\nリップ最大 %.0f kW/m²" % (hi, qmax_all),
                transform=ax.transAxes, fontsize=8.5, color="w", va="bottom", ha="left")
        ax.set_xlim(0, 180); ax.set_ylim(-dep * 1e3, 0)

    # ---------- (c) すきま中央面の速度ベクトル ----------
    ax = fig.add_subplot(gs[0, 2])
    nth, nz = 61, 60
    tg = np.radians(np.linspace(0, 180, nth))
    zg = np.linspace(-dep * 0.999, -dep * 0.002, nz)
    rc = gc.gap_center_radius(tg, man)
    TH, ZZ = np.meshgrid(tg, zg)
    RC = np.tile(rc, (nz, 1))
    DX, DY = gc.ray_dir(TH)
    P = np.stack([RC * DX, RC * DY, ZZ], axis=-1).reshape(-1, 3)
    ux, uy, uz = (L["Ux"](P), L["Uy"](P), L["Uz"](P))
    # θ 方向の単位ベクトル = d/dθ(-cosθ, sinθ) = (sinθ, cosθ)
    ut = (ux * np.sin(TH.ravel()) + uy * np.cos(TH.ravel())).reshape(nz, nth)
    uz = uz.reshape(nz, nth)
    sp = np.hypot(ut, uz)
    good = np.isfinite(sp)
    sp = np.where(good, sp, 0.0); ut = np.where(good, ut, 0.0); uz = np.where(good, uz, 0.0)
    lo = max(np.percentile(sp[sp > 0], 2) if np.any(sp > 0) else 1e-3, 1e-3)
    hi = max(sp.max(), lo * 10)
    cf = ax.contourf(np.degrees(TH), ZZ * 1e3, np.clip(sp, lo, hi),
                     levels=np.geomspace(lo, hi, 24), cmap="viridis",
                     norm=matplotlib.colors.LogNorm(lo, hi))
    import matplotlib.ticker as mtick
    cb = fig.colorbar(cf, ax=ax, label="|u| [m/s] (対数)",
                      ticks=mtick.LogLocator(base=10, subs=(1.0, 3.0)))
    cb.ax.yaxis.set_major_formatter(mtick.FuncFormatter(
        lambda v, p: ("%g" % v) if v >= 1 else ("%.2g" % v)))
    cb.ax.minorticks_off()
    st, sz_ = 2, 3
    un = np.maximum(sp, 1e-12)
    ax.quiver(np.degrees(TH)[::sz_, ::st], (ZZ * 1e3)[::sz_, ::st],
              (ut / un)[::sz_, ::st], (uz / un)[::sz_, ::st], color="w",
              angles="xy", pivot="mid", scale=26, width=.0045, headwidth=3.6, headlength=4.2)
    ax.set_xlabel("周方向 θ [deg]  (0=上流, 180=下流)"); ax.set_ylabel("深さ z [mm]")
    ax.set_title("(c) すきま中央面の流れ (向き=矢印, 速さ=色)", fontsize=12, loc="left")

    # ---------- (d)(e) 半径-深さ断面 ----------
    for i, (thd, ttl) in enumerate(((0.0, "θ=0° (上流側)"), (180.0, "θ=180° (下流側)"))):
        ax = fig.add_subplot(gs[1, i])
        t = np.radians(thd)
        ri = gc.inner_radius_at(np.array([t]), man)[0]
        rg = np.linspace(ri + 1e-5, Ro - 1e-5, 44)
        zg2 = np.linspace(-dep * 0.999, -dep * 0.0005, 180)
        RG, ZG = np.meshgrid(rg, zg2)
        dxt, dyt = gc.ray_dir(t)
        P2 = np.stack([RG * dxt, RG * dyt, ZG], axis=-1).reshape(-1, 3)
        T2 = L["T"](P2).reshape(ZG.shape)
        ur = (L["Ux"](P2) * dxt + L["Uy"](P2) * dyt).reshape(ZG.shape)
        uz2 = L["Uz"](P2).reshape(ZG.shape)
        T2 = np.where(np.isfinite(T2), T2, Tw)
        ur = np.where(np.isfinite(ur), ur, 0.0); uz2 = np.where(np.isfinite(uz2), uz2, 0.0)
        cf = ax.contourf((RG - ri) * 1e3, ZG * 1e3, T2, levels=24, cmap="coolwarm")
        fig.colorbar(cf, ax=ax, label="T [K]")
        sz, sr = 7, 4
        spd = np.maximum(np.hypot(ur, uz2), 1e-12)
        ax.quiver(((RG - ri) * 1e3)[::sz, ::sr], (ZG * 1e3)[::sz, ::sr],
                  (ur / spd)[::sz, ::sr], (uz2 / spd)[::sz, ::sr], color="k",
                  angles="xy", pivot="mid", scale=22, width=.005, headwidth=3.6, headlength=4.2)
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
            "**θ=0° が上流側** (流れは +x なので θ=180° が下流)。\n"
            "(c)(d)(e) の矢印は**長さを正規化**して向きだけを示す (速さは色)。\n"
            "(d)(e) の横軸は内円柱表面からの距離。"
            % (a.run, step, G["gap_nom"] * 1e3, dep * 1e3, Tw, D["mach"], "\n".join(wh_txt)),
            va="top", ha="left", fontsize=10.5, transform=ax.transAxes)

    out = Path(a.out or (Path(a.run) / "cavity_fields.png"))
    fig.suptitle("case/49  すきま内部の流れと壁熱流束  —  %s" % a.run, fontsize=13)
    fig.savefig(out, dpi=125, bbox_inches="tight")
    print("wrote", out)

    # ---------------------------------------------------------------- 熱伝達率・底面
    wh = wall_heat(a.run, step, man, D)
    try:
        sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "solver_density_cuda" / "tools"))
        from total_quantities import total_state
        T0 = np.asarray(total_state(a.run, str(res))["T0"], float)
        wall_href(wh, man, D, c, v, T0)
    except Exception as e:                       # noqa: BLE001
        print("  WARNING: h_ref 不可 (%s)" % e)
    plot_htc(a.run, step, man, D, wh, Path(a.run) / "cavity_htc.png")
    plot_floor(a.run, step, man, D, wh, Path(a.run) / "cavity_floor.png")


def plot_floor(run, step, man, D, wh, out):
    """底面 (cav_floor) の熱流束・熱伝達率を θ × すきま横断位置に展開して描く。"""
    if "cav_floor" not in wh:
        print("  cav_floor の壁ダンプが無い -> 底面図はスキップ")
        return
    d = wh["cav_floor"]
    u, vv, vlab = wall_uv(d, man, "cav_floor")
    hr = d.get("_href_node")
    npan = 3 if hr is not None else 2
    fig, ax = plt.subplots(1, npan, figsize=(5.4 * npan, 4.6))
    U, V, Q = grid_mean(u, vv, d["_qin_node"] * 1e-3, d["_w_node"], nu=120, nv=30, ulim=(0, 180))
    hi = np.nanpercentile(Q, 99.5)
    cf = ax[0].contourf(U, V, np.clip(Q, 0, hi), levels=np.linspace(0, hi, 21), cmap="inferno")
    fig.colorbar(cf, ax=ax[0], label="q'' [kW/m²]")
    ax[0].set_title("(a) 底面の熱流束", fontsize=12, loc="left")
    if hr is not None:
        U, V, H = grid_mean(u, vv, hr, d["_w_node"], nu=120, nv=30, ulim=(0, 180))
        hi2 = np.nanpercentile(H, 99.5)
        cf = ax[1].contourf(U, V, np.clip(H, 0, hi2), levels=np.linspace(0, hi2, 21), cmap="viridis")
        fig.colorbar(cf, ax=ax[1], label="h_ref [W/m²K]")
        ax[1].set_title("(b) 底面の熱伝達率 (基準温度基準)", fontsize=12, loc="left")
    for b in ax[:npan - 1]:
        b.set_xlabel("周方向 θ [deg]  (0=上流, 180=下流)")
        b.set_ylabel(vlab)
        b.set_xlim(0, 180)
    th = u
    # 周方向分布 (すきま横断方向に面積平均)
    b = ax[npan - 1]
    nb = 60
    edges = np.linspace(0, 180, nb + 1)
    ib = np.clip(np.digitize(th, edges) - 1, 0, nb - 1)
    w = d["_w_node"]
    den = np.bincount(ib, weights=w, minlength=nb)
    qb = np.bincount(ib, weights=d["_qin_node"] * w, minlength=nb) / np.maximum(den, 1e-30)
    tc = 0.5 * (edges[1:] + edges[:-1])
    b.plot(tc, qb * 1e-3, lw=2.2, color="#dc2626", label="q'' [kW/m²]")
    b.set_xlabel("周方向 θ [deg]  (0=上流, 180=下流)")
    b.set_ylabel("q'' [kW/m²]", color="#dc2626")
    b.grid(alpha=.3)
    if hr is not None:
        hb = np.bincount(ib, weights=hr * w, minlength=nb) / np.maximum(den, 1e-30)
        b2 = b.twinx()
        b2.plot(tc, hb, lw=2.2, color="#1d4ed8", label="h_ref")
        b2.set_ylabel("h_ref [W/m²K]", color="#1d4ed8")
        Tb = np.bincount(ib, weights=d["_Tref_node"] * w, minlength=nb) / np.maximum(den, 1e-30)
        b.set_title("(c) 周方向分布 (横断平均)  基準温度 %.0f–%.0f K"
                    % (Tb.min(), Tb.max()), fontsize=12, loc="left")
    else:
        b.set_title("(c) 周方向分布 (横断平均)", fontsize=12, loc="left")
    fig.suptitle("case/49  キャビティ底面の熱流束・熱伝達率  —  %s  (Q=%.3f W 半割, q'' 平均 %.0f W/m²)"
                 % (run, d["Q_W"], d["qpp_mean"]), fontsize=12.5)
    fig.tight_layout()
    fig.savefig(out, dpi=130, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
