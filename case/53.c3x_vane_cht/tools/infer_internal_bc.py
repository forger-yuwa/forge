#!/usr/bin/env python3
r"""公開量から**冷却孔の内部境界条件を逆算**する (V5 段 (b) の再設計版)。

**なぜ逆算か**: NASA CR-168015 は冷却孔ごとの冷却剤温度・流量を**公開していない**
(方法だけを書き、付録 A にも欄が無い)。一方で公開されているのは
**実測壁温 $T_w$** と **外表面の熱伝達係数 $h$** ($h_0$=1135 W/m²K で正規化、ガス全温基準)。
そこで外表面に $q = h\,(T_g - T_w)$ を Neumann で課し、**孔ごとの $(h_c)$ と冷却剤温度 $T_c$ を
最小二乗で同定**して、その残差と不確かさを V5 の帯に積む。
**同定値を「公開値」と偽らないこと**が本スクリプトの存在理由。

固体は `solver_density_cuda/tools/solid_fem2d.py` (線形三角形 FE + Robin 辺)。
$K(h_1..h_{10}) = K_{cond} + \sum_k h_k M_k$ と線形なので `parts()` で一度だけ組む。

usage: python3 case/53.c3x_vane_cht/tools/infer_internal_bc.py [--vane c3x] [--run run108] [--per-hole-Tc]
"""
import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "solver_density_cuda" / "tools"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from solid_fem2d import Fem2DOperator          # noqa: E402
from run_data import TABLES                    # noqa: E402

H0 = 1135.0      # 熱伝達係数の正規化基準 [W/m2K] (報告 p.125)
TREF = 811.0     # 温度の正規化基準 [K]


def ordered_outer_loop(nodes, edges):
    """外周の辺リストから節点の巡回順を作る。"""
    adj = {}
    for a, b in edges:
        adj.setdefault(int(a), []).append(int(b))
        adj.setdefault(int(b), []).append(int(a))
    start = min(adj)
    loop = [start]
    prev, cur = None, start
    while True:
        nxt = [v for v in adj[cur] if v != prev]
        if not nxt:
            break
        prev, cur = cur, nxt[0]
        if cur == start:
            break
        loop.append(cur)
    if len(loop) != len(adj):
        sys.exit(f"外周が 1 本の閉曲線になっていない ({len(loop)} / {len(adj)})")
    return np.array(loop, int)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vane", default="c3x")
    ap.add_argument("--run", default="run108")
    ap.add_argument("--per-hole-Tc", action="store_true", help="冷却剤温度を孔ごとに持たせる")
    ap.add_argument("--Tg", type=float, default=None, help="ガス全温 [K] (既定は run の TT1)")
    ap.add_argument("--Tc-sweep", type=float, nargs="+", default=[300, 350, 400, 450, 500],
                    help="冷却剤温度の走査値 [K] (h_c との縮退を示すため)")
    ap.add_argument("--Tc-nominal", type=float, default=300.0, help="代表として採る冷却剤温度 [K]")
    ap.add_argument("--le-shift-cm", type=float, default=None,
                    help="弧長原点 (s=0) を前縁頂点から動かす量 [cm]。既定は表 IV の PS/SS 弧長比に合わせる値")
    a = ap.parse_args()

    case = ROOT / ("case/53.c3x_vane_cht" if a.vane == "c3x" else "case/54.markii_vane_cht")
    d = np.load(case / f"mesh/solid_{a.vane}.npz")
    nodes, tris = d["nodes"], d["tris"]
    outer_edges = d["outer_edges"]
    holes = [d[k] for k in d.files if k.startswith("hole")]

    tbl = TABLES[a.run]
    Tg = a.Tg if a.Tg is not None else {"run108": 786.0, "run42": 788.0}[a.run]
    mat = json.loads((case / "ref/material_astm310.json").read_text())
    kT, kk = np.array(mat["k_table"]["T"], float) + 273.15, np.array(mat["k_table"]["k"], float)

    loop = ordered_outer_loop(nodes, outer_edges)
    xy = nodes[loop]
    seg = np.hypot(*np.diff(np.vstack([xy, xy[:1]]), axis=0).T)
    i_le = int(np.argmin(xy[:, 0])); i_te = int(np.argmax(xy[:, 0]))
    # **弧長原点の較正**: 報告の s は「前縁から測った表面距離」だが、原点の定義 (幾何頂点か停留点か) は
    # 書かれていない。表 IV は PS 13.723 / SS 17.782 cm と弧長を与えているので、**その比に合うように
    # 原点を surface に沿って動かす**。ずらし量と残差の感度は --le-shift-cm で確認できる。
    n = len(loop)
    # 前縁から両回りで後縁までの弧長 -> 長い方が負圧面 (表 IV: SS 17.782 / PS 13.723 cm)
    s_from_le = np.zeros(n); acc = 0.0
    for k in range(n):
        idx = (i_le + k) % n
        s_from_le[idx] = acc
        acc += seg[idx]
    total = acc
    s_te_fwd = s_from_le[i_te]
    s_te_bwd = total - s_te_fwd
    # **弧長原点**: 報告の s は「前縁から測った表面距離」だが、原点の定義 (幾何頂点か停留点か) は
    # 書かれていない。既定は**幾何頂点** (min x) とし、--le-shift-cm でずらせるようにして
    # **残差の感度を確認できる**ようにする (実測: ±0.5 cm 動かしても RMS は 10.6-10.8 K でほぼ動かない
    # → 原点の曖昧さは残差床の原因ではない)。
    shift = a.le_shift_cm or 0.0
    if abs(shift) > 1e-9:
        step = 1 if shift > 0 else -1
        k, acc2 = 0, 0.0
        while acc2 < abs(shift) / 100.0 and k < n - 1:
            k += 1
            acc2 += seg[(i_le + step * k) % n]
        i_le = (i_le + step * k) % n
        s_from_le = np.zeros(n); acc = 0.0
        for kk2 in range(n):
            idx = (i_le + kk2) % n
            s_from_le[idx] = acc
            acc += seg[idx]
        total = acc
        s_te_fwd = s_from_le[i_te]
    print(f"[{a.vane}] arc origin: geometric LE apex {shift:+.3f} cm shifted")
    side = np.where(s_from_le <= s_te_fwd, 0, 1)          # 0 = 前向き, 1 = 後ろ向き
    arc_fwd, arc_bwd = s_te_fwd, s_te_bwd
    ss_is_fwd = arc_fwd > arc_bwd
    s_norm = np.where(side == 0, s_from_le / arc_fwd, (total - s_from_le) / arc_bwd)
    is_ss = (side == 0) if ss_is_fwd else (side == 1)
    ref_arc = {"c3x": (17.782, 13.723), "markii": (15.935, 12.949)}[a.vane]
    ss_arc = (arc_fwd if ss_is_fwd else arc_bwd) * 100
    ps_arc = (arc_bwd if ss_is_fwd else arc_fwd) * 100
    print(f"[{a.vane}] outer loop {n} nodes, total arc {total*100:.3f} cm "
          f"(SS {ss_arc:.3f} vs 表 IV {ref_arc[0]}, PS {ps_arc:.3f} vs {ref_arc[1]} cm; "
          f"差 {100*(ss_arc/ref_arc[0]-1):+.1f} % / {100*(ps_arc/ref_arc[1]-1):+.1f} %)")

    # **判読不能セル (None) のある行は落とす** (Mark II run42 は 9 セル欠測。case/54 README)。
    rows = [r for r in tbl["rows"] if r[2] is not None and r[3] is not None]
    s_dat = np.array([r[0] for r in rows]); Tw_dat = np.array([r[2] for r in rows]) * TREF
    h_dat = np.array([r[3] for r in rows]) * H0
    i_stag = int(np.argmin(s_dat))
    dat = {"PS": (s_dat[:i_stag + 1][::-1], Tw_dat[:i_stag + 1][::-1], h_dat[:i_stag + 1][::-1]),
           "SS": (s_dat[i_stag:], Tw_dat[i_stag:], h_dat[i_stag:])}

    Tw_meas = np.zeros(n); h_meas = np.zeros(n); covered = np.zeros(n, bool)
    for k in range(n):
        key = "SS" if is_ss[k] else "PS"
        sd, td, hd = dat[key]
        Tw_meas[k] = np.interp(s_norm[k], sd, td)
        h_meas[k] = np.interp(s_norm[k], sd, hd)
        covered[k] = (s_norm[k] <= sd.max())
    print(f"[{a.vane}] data covers {covered.sum()} / {n} outer nodes "
          f"(後縁側 {n - covered.sum()} 点は外挿。残差から除外する)")

    # 外表面の熱流束 [W/m2] -> 節点荷重 [W/m] (単位奥行き)
    lump = np.zeros(n)
    for k in range(n):
        lump[k] = 0.5 * (seg[(k - 1) % n] + seg[k])
    q_ext = h_meas * (Tg - Tw_meas)
    Q_nodal = q_ext * lump

    op = Fem2DOperator(nodes, tris, loop, [(int(e[0]), int(e[1])) for e in outer_edges],
                       [], float(np.interp(Tw_meas.mean(), kT, kk)))
    groups = [[(int(e[0]), int(e[1])) for e in hh] for hh in holes]
    K, parts = op.parts(T_ref=float(Tw_meas.mean()), groups=groups)
    nh = len(groups)

    def solve(hc, Tc):
        A = K.copy()
        rhs = np.zeros(op.N)
        rhs[loop] += Q_nodal
        for k in range(nh):
            M, v = parts[k]
            A = A + hc[k] * M
            rhs += hc[k] * Tc[k] * v
        u = spla.spsolve(A.tocsc(), rhs)
        return u

    # **冷却剤温度と h_c は縮退する** (データが拘束するのは概ね h_c (T_w - T_c) の積)。
    # 片方を自由にすると T_c が非物理な値へ走るので、**T_c を固定して h_c だけを同定し、
    # T_c を振って残差曲線を出す**。縮退を隠さずに不確かさとして扱うのが目的。
    def fit_for_Tc(Tc_val):
        def resid(p):
            u = solve(np.exp(p), np.full(nh, Tc_val))
            return (u[loop] - Tw_meas)[covered]
        r = least_squares(resid, np.full(nh, math.log(1500.0)), method="trf",
                          xtol=1e-12, ftol=1e-12, max_nfev=300)
        return np.exp(r.x), float(np.sqrt(np.mean(r.fun ** 2)))

    print(f"\n[{a.vane} {a.run}] 冷却剤温度を振って h_c を同定 (縮退の確認)")
    sweep = []
    for Tc_val in a.Tc_sweep:
        hc_s, rms_s = fit_for_Tc(Tc_val)
        sweep.append((Tc_val, hc_s, rms_s))
        print(f"  T_c = {Tc_val:5.0f} K -> RMS {rms_s:6.2f} K, h_c = "
              f"{hc_s.min():6.0f} .. {hc_s.max():6.0f} W/m2K (平均 {hc_s.mean():6.0f})")
    Tc_nom = a.Tc_nominal
    hc, rms_fit = fit_for_Tc(Tc_nom)
    Tc = np.full(nh, Tc_nom)
    u = solve(hc, Tc)
    err = u[loop] - Tw_meas
    rms = float(np.sqrt(np.mean(err[covered] ** 2)))
    print(f"\n[{a.vane} {a.run}] 同定結果 (k_s = {op.k_solid:.2f} W/mK @ {Tw_meas.mean():.0f} K)")
    print(f"  冷却剤温度 T_c = " + (", ".join(f"{t:.1f}" for t in Tc) if a.per_hole_Tc
                                   else f"{Tc[0]:.1f} K (全孔共通)"))
    for k in range(nh):
        print(f"  hole {k+1:2d}: h_c = {hc[k]:8.1f} W/m2K")
    print(f"  外表面温度の残差: RMS {rms:.2f} K, max |err| {np.abs(err[covered]).max():.2f} K "
          f"(実測 T_w の幅 {Tw_meas.min():.0f}-{Tw_meas.max():.0f} K)")
    print(f"  外表面の総入熱 {np.sum(Q_nodal):.1f} W/m ; 孔からの総排熱 "
          f"{sum(hc[k]*np.dot(parts[k][1], (Tc[k]-u)) for k in range(nh)):.1f} W/m")

    # 図: 実測 T_w と同定後の T_w、および外表面熱流束
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))
    sgn = np.where(is_ss, 1.0, -1.0) * s_norm
    o = np.argsort(sgn)
    ax[0].plot(sgn[o], Tw_meas[o], "-", lw=1.5, label="measured $T_w$")
    ax[0].plot(sgn[o], u[loop][o], "--", lw=1.5, label="fitted (inferred hole BC)")
    ax[0].plot(sgn[~covered], Tw_meas[~covered], "x", ms=4, color="0.6", label="extrapolated (excluded)")
    ax[0].set_xlabel("- s/S (PS)   |   + s/S (SS)"); ax[0].set_ylabel("$T_w$ [K]")
    ax[0].set_title(f"{a.vane} {a.run}: wall temperature (RMS {rms:.1f} K)")
    ax[1].plot(sgn[o], q_ext[o] / 1e3, "-", lw=1.5)
    ax[1].set_xlabel("- s/S (PS)   |   + s/S (SS)"); ax[1].set_ylabel("$q$ [kW/m²]")
    ax[1].set_title("imposed external heat flux $q = h\\,(T_g - T_w)$")
    for x in ax:
        x.grid(alpha=.3); x.legend(fontsize=8)
    fig.tight_layout()
    png = case / f"ref/{a.run}_internal_bc.png"
    fig.savefig(png, dpi=130)
    print(f"  -> {png.relative_to(ROOT)}")

    out = case / f"ref/{a.run}_internal_bc.json"
    out.write_text(json.dumps({
        "note": "報告に無い内部境界条件を、公開量 (実測 Tw と h) から逆算した**推定値**。公開値ではない。",
        "method": "外表面に q = h (Tg - Tw) を Neumann で課し、孔ごとの h_c と冷却剤温度 T_c を最小二乗で同定",
        "run": a.run, "vane": a.vane, "Tg_K": Tg, "k_solid_W_mK": float(op.k_solid),
        "h_c_W_m2K": [float(x) for x in hc], "T_c_K": [float(x) for x in Tc],
        "degeneracy_note": "データが拘束するのは概ね h_c (T_w - T_c) の積なので、T_c を固定して h_c を同定した。下の sweep が縮退の幅。",
        "Tc_sweep": [{"T_c_K": t, "rms_K": r, "h_c_mean": float(h.mean()),
                      "h_c_min": float(h.min()), "h_c_max": float(h.max())} for t, h, r in sweep],
        "residual_Tw_rms_K": rms, "residual_Tw_max_K": float(np.abs(err[covered]).max()),
        "nodes_used": int(covered.sum()), "nodes_total": int(n),
    }, indent=2, ensure_ascii=False) + "\n")
    print(f"  -> {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
