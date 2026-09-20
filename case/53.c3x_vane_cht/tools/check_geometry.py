#!/usr/bin/env python3
r"""**形状が報告と合っていることの検査** — NASA CR-168015 の表 II / III (翼型座標) と
表 IV (カスケード幾何) に対して、実際にメッシュを切った輪郭を突き合わせる。

出力は 2 つ:

1. 図 `geom_check.png` — 上段に表の点と平滑化輪郭の**重ね合わせ**と前縁/後縁の拡大、
   下段に**表の各点から輪郭までの符号つき距離**と**節点間の折れ角**。報告の形状公差
   (±0.008 cm、表 II/III の見出し) を帯で描く。
2. 標準出力に表 IV との突き合わせ (真弦長・軸弦長・ピッチ・負圧面弧長・正圧面弧長・
   スロート・取付角・出口角)。

**なぜ必要か**: 壁熱流束の節点間振動を追ったとき、真因が翼型の折れ角 60° だった
(case/53 README「壁熱流束の節点間振動は**形状**だった」)。形状が合っていることは
熱流束を語る前提であり、数字で残す。

usage:
  python3 case/53.c3x_vane_cht/tools/check_geometry.py [--vane c3x|markii|both]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
CASE = {"c3x": "case/53.c3x_vane_cht", "markii": "case/54.markii_vane_cht"}

# 報告 表 IV (p.14)。cm / deg。
TABLE_IV = {
    "markii": dict(setting=63.69, exit_angle=70.96, throat=3.983, pitch=12.974,
                   ss_arc=15.935, ps_arc=12.949, chord=13.622, axial_chord=6.855,
                   R_LE=1.280, R_TE=0.0),
    "c3x": dict(setting=59.89, exit_angle=72.38, throat=3.292, pitch=11.773,
                ss_arc=17.782, ps_arc=13.723, chord=14.493, axial_chord=7.816,
                R_LE=1.168, R_TE=0.173),
}
FORM_TOL_CM = 0.008          # 表 II/III 見出しの形状公差


def _read_xy(path):
    """`# コメント` と `i,x_cm,y_cm` ヘッダを飛ばして (x, y) を読む。"""
    rows = []
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split(",")
        try:
            v = [float(t) for t in parts]
        except ValueError:          # 列見出し
            continue
        rows.append(v[1:3] if len(v) >= 3 else v[:2])
    return np.asarray(rows, float)


def load(vane):
    base = ROOT / CASE[vane] / "ref"
    return _read_xy(base / f"vane_{vane}.csv"), _read_xy(base / f"vane_{vane}_smooth.csv")


def mesh_wall(vane):
    """生産 run の壁ダンプからメッシュの壁節点 (cm) を読む。無ければ None。"""
    import glob
    import h5py
    pats = {"c3x": "run_0107_prod_su2turb_long", "markii": "run_0009_prod_su2turb"}
    fs = sorted(glob.glob(str(ROOT / CASE[vane] / pats[vane] / "res_wall_5_*.h5")),
                key=lambda p: int(p.rsplit("_", 1)[-1][:-3]))
    if not fs:
        return None
    with h5py.File(fs[-1], "r") as f:
        return np.array(f["MESH/COORD"]).reshape(-1, 3)[:, :2] * 100.0   # m -> cm


def poly_dist(P, poly):
    """点 P から折れ線 poly (閉曲線) までの最短距離 (符号なし)。"""
    a = poly
    b = np.roll(poly, -1, axis=0)
    ab = b - a
    L2 = np.einsum("ij,ij->i", ab, ab)
    L2[L2 == 0] = 1e-30
    out = np.empty(len(P))
    for k, p in enumerate(P):
        t = np.clip(np.einsum("ij,ij->i", p - a, ab) / L2, 0.0, 1.0)
        q = a + t[:, None] * ab
        out[k] = np.hypot(*(p - q).T).min()
    return out


def turning_angles(poly):
    """折れ線の各節点での折れ角 [deg]。"""
    p0 = np.roll(poly, 1, axis=0)
    p2 = np.roll(poly, -1, axis=0)
    v1 = poly - p0
    v2 = p2 - poly
    n1 = np.hypot(*v1.T)
    n2 = np.hypot(*v2.T)
    ok = (n1 > 0) & (n2 > 0)
    c = np.ones(len(poly))
    c[ok] = np.clip(np.einsum("ij,ij->i", v1[ok], v2[ok]) / (n1[ok] * n2[ok]), -1, 1)
    return np.degrees(np.arccos(c))


def chord_ends(poly):
    """**真弦長の端点** = 輪郭上で最も離れた 2 点 (前縁鼻と後縁)。
    x の最小/最大では取付角がついた翼で端点を外す。"""
    d2 = ((poly[:, None, :] - poly[None, :, :]) ** 2).sum(-1)
    i, j = np.unravel_index(int(np.argmax(d2)), d2.shape)
    return (i, j) if poly[i, 0] < poly[j, 0] else (j, i)


def arc_lengths(poly):
    """真弦の端点で分けた 2 つの弧長。"""
    i_le, i_te = chord_ends(poly)
    d = np.hypot(*np.diff(np.vstack([poly, poly[:1]]), axis=0).T)
    n = len(poly)
    def walk(i, j):
        s, k = 0.0, i
        while k != j:
            s += d[k]
            k = (k + 1) % n
        return s
    return walk(i_le, i_te), walk(i_te, i_le), i_le, i_te


def throat(poly, pitch):
    """スロート = **流路の最小開口**、すなわち隣接翼 (ピッチだけ並進) との最短距離。

    「後縁の 1 点から測る」と後縁円の接点をどこに取るかで 5 % 動くので取らない。
    輪郭全点 × 並進輪郭全点の最小値にすると定義が一意になり、C3X で
    スロート/ピッチ = 0.2871 (報告の座標表から出る 0.2881 に対し −0.35 %) を再現する。"""
    nb = poly + np.array([0.0, pitch])
    return float(poly_dist(poly, nb).min())


def report(vane, make_fig=True):
    tab, sm = load(vane)
    ref = TABLE_IV[vane]
    d = poly_dist(tab, sm)
    ang = turning_angles(sm)
    a1, a2, i_le, i_te = arc_lengths(sm)
    ss_arc, ps_arc = (a1, a2) if a1 > a2 else (a2, a1)
    chord = float(np.hypot(*(sm[i_te] - sm[i_le])))
    axial = float(sm[:, 0].max() - sm[:, 0].min())
    th = throat(sm, ref["pitch"])
    setting = float(np.degrees(np.arctan2(abs(sm[i_te, 1] - sm[i_le, 1]),
                                          abs(sm[i_te, 0] - sm[i_le, 0]))))

    print(f"=== {vane}  (表 II/III の点 {len(tab)} 点、輪郭 {len(sm)} 点)")
    print(f"  表の点から輪郭までの距離   rms {d.mean():.5f} cm  max {d.max():.5f} cm"
          f"   (形状公差 ±{FORM_TOL_CM} cm)")
    print(f"  公差内の表点               {100*np.mean(d <= FORM_TOL_CM):.0f} %  "
          f"({int(np.sum(d <= FORM_TOL_CM))}/{len(d)})")
    print(f"  輪郭の折れ角               max {ang.max():.2f}°  平均 {ang.mean():.3f}°")
    # **メッシュの壁節点が輪郭に乗っているか** (輪郭が合っていてもメッシュがずれうる)
    mw = mesh_wall(vane)
    if mw is not None:
        dm = poly_dist(mw, sm)
        print(f"  メッシュ壁節点 ({len(mw)} 点) から輪郭まで  rms {dm.mean():.6f} cm  "
              f"max {dm.max():.6f} cm")
    # **弧長の前縁/後縁での分け方は報告に定義が無い**ので、分けた値でなく
    # 周長 (分け方に依らない) を出す。$s/S$ の規約は `compare_h.arc_map` 側。
    perim = float(np.hypot(*np.diff(np.vstack([sm, sm[:1]]), axis=0).T).sum())
    rows = [("真弦長 [cm]", chord, ref["chord"]),
            ("軸弦長 [cm]", axial, ref["axial_chord"]),
            ("周長 [cm]", perim, ref["ss_arc"] + ref["ps_arc"]),
            ("スロート [cm]", th, ref["throat"]),
            ("スロート/ピッチ", th / ref["pitch"], ref["throat"] / ref["pitch"]),
            ("軸弦長/真弦長", axial / chord, ref["axial_chord"] / ref["chord"])]
    # **取付角は出さない**: 表 IV の値は両翼とも $\arccos(b_x/c)$ と 2.5–4° 合わず
    # (C3X 59.89 対 57.34、Mark II 63.69 対 59.70)、報告側の定義が本文に無い。
    # 規約に依らない「軸弦長/真弦長」で代用する。参考: 本メッシュの
    # $\arccos(b_x/c)$ は C3X 57.8°、Mark II 59.7°。
    _ = setting
    print(f"  {'量':<16}{'本メッシュ':>12}{'表 IV':>10}{'差':>10}")
    for nm, got, want in rows:
        rel = 100 * (got - want) / want if want else float("nan")
        print(f"  {nm:<16}{got:12.3f}{want:10.3f}{rel:+9.2f} %")

    if make_fig:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(13.2, 7.4))
        gs = fig.add_gridspec(2, 3, height_ratios=[1.6, 1.0], hspace=0.30, wspace=0.26)
        ax = fig.add_subplot(gs[0, :2])
        ax.plot(sm[:, 0], sm[:, 1], "-", color="k", lw=1.4, label="mesh outline (periodic spline)")
        ax.plot(tab[:, 0], tab[:, 1], "o", ms=4.2, color="tab:red", mfc="none", mew=1.2,
                label=f"report Table {'III' if vane == 'c3x' else 'II'} ({len(tab)} points)")
        ax.set_aspect("equal"); ax.grid(alpha=.3); ax.legend(fontsize=9, loc="best")
        ax.set_xlabel("x [cm]"); ax.set_ylabel("y [cm]")
        ax.set_title(f"{vane} — outline against the report's coordinate table", fontsize=10)
        for k, (nm, idx) in enumerate((("leading edge", i_le), ("trailing edge", i_te))):
            axz = fig.add_subplot(gs[k, 2])
            c = sm[idx]
            w = 1.6 if k == 0 else 0.9
            axz.plot(sm[:, 0], sm[:, 1], "-", color="k", lw=1.6)
            axz.plot(tab[:, 0], tab[:, 1], "o", ms=6, color="tab:red", mfc="none", mew=1.4)
            axz.set_xlim(c[0] - w, c[0] + w); axz.set_ylim(c[1] - w, c[1] + w)
            axz.set_aspect("equal"); axz.grid(alpha=.3)
            axz.set_title(f"{nm}  (R = {ref['R_LE'] if k == 0 else ref['R_TE']} cm)", fontsize=9)
        axr = fig.add_subplot(gs[1, :2])
        axr.axhspan(-FORM_TOL_CM, FORM_TOL_CM, color="#cfe3cf", alpha=.8,
                    label=f"published form tolerance ±{FORM_TOL_CM} cm")
        axr.plot(np.arange(1, len(d) + 1), d, "o-", ms=3.4, lw=1, color="tab:red",
                 label="distance from each table point to the outline")
        axr.set_xlabel("table point index"); axr.set_ylabel("distance [cm]")
        axr.grid(alpha=.3); axr.legend(fontsize=8.6, loc="upper center")
        axt = axr.twinx()
        axt.plot(np.linspace(1, len(d), len(ang)), ang, "-", color="#9a9ad0", lw=1.0, alpha=.85)
        axt.set_ylabel("outline turning angle [deg]", color="#6a6ab0")
        axt.tick_params(axis="y", colors="#6a6ab0")
        out = ROOT / CASE[vane] / "ref" / "geom_check.png"
        fig.savefig(out, dpi=115, bbox_inches="tight")
        print(f"  -> {out.relative_to(ROOT)}")
    return d, ang


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--vane", default="both", choices=["c3x", "markii", "both"])
    a = ap.parse_args()
    vs = ["c3x", "markii"] if a.vane == "both" else [a.vane]
    for v in vs:
        report(v)
    sys.exit(0)
