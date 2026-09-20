#!/usr/bin/env python3
r"""翼型輪郭を**滑らかな曲線**にして書き出す (メッシュ生成はこちらを使う)。

**なぜ要るか** (2026-09-20, case/53 run_0002 で実測): 表の点 (78 点 / 周長 32 cm = 4〜6 mm 間隔) を
そのまま折れ線にすると、(a) 負圧面前縁の高曲率域で 1 弦 6 mm が **20°** も曲がり、(b) 復元した
前縁・後縁円弧と表の点の継ぎ目で接線が **最大 60°** 飛ぶ。その角で流れが局所加速し、
壁熱流束が隣り合う節点間で 52 → 302 kW/m² と振れて $h$ の分布が読めなくなった
(折れ角と $|q-\mathrm{smooth}(q)|$ の相関 0.57)。**形状の問題であってソルバの市松ではない**。

**やること**: 全点 (表 + 円弧) を通る**周期平滑化スプライン**を弦長で張る。重みは
報告の外形不確かさ **±0.008 cm** (表) と、構成した円弧の厳密さ (0.001 cm) から与え、
平滑化量 $s=m$ (残差が概ね $1\sigma$) にする。**点を動かして形を作らない**ための検査を付ける:

  1. 表の点からのずれ: RMS ≤ 0.008 cm、最大 ≤ 0.024 cm (報告の不確かさの 3 倍)
  2. 再標本化後の折れ角: 最大 ≤ `--max-turn` (既定 5°)
  3. 前縁・後縁の曲率半径が $R_{LE}$ / $R_{TE}$ の ±20 % 以内
  4. 自己交差が無い

**どれか落ちたら書き出さない** (メッシュを細かくして隠さない。[[geometry-must-not-follow-mesh-spacing]])。

usage: python3 case/53.c3x_vane_cht/tools/smooth_profile.py [--vane c3x] [--n 640] [--max-turn 5]
"""
import argparse
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import splprep, splev
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from gen_solid_mesh import REF                         # noqa: E402

RADII = {"c3x": dict(LE=1.168, TE=0.173), "markii": dict(LE=1.280, TE=0.0)}
SIG_TABLE = 0.008      # cm、報告 p.24 の外形プロファイル不確かさ
SIG_ARC = 0.001        # cm、円弧は構成値なのでほぼ厳密



ARC_TOL = 0.003        # cm、円弧区間に採る許容 (座標不確かさ 0.008 cm より厳しく取る)


def fit_arc_span(tab, i0, R, tol=ARC_TOL):
    """頂点 `i0` まわりで**半径 R の円に乗る最大区間**を探し、(中心, a, b) を返す。

    中心は R を固定した最小二乗で決める。**中心を「頂点 ± R を x 方向に取る」と仮定しない**:
    C3X は取付角 60° なので鼻の向きが x 軸ではなく、その仮定だと中心が 0.06 cm ずれ、
    区間端で表点が円から 0.05 cm 外れる。円弧をそこに繋ぐと**半径方向に段差ができ、
    折れ角 60° の角**になる (実測 2026-09-20: その角で壁熱流束が 52 → 302 kW/m² と振れた)。
    """
    n = len(tab)

    def fit(idx):
        Q = tab[idx]
        c0 = Q.mean(axis=0)
        r = least_squares(lambda c: np.hypot(*(Q - c).T) - R, c0)
        return r.x, float(np.max(np.abs(np.hypot(*(Q - r.x).T) - R)))

    best = None
    a, b = -1, 1
    while True:
        grown = False
        for da, db in ((a - 1, b), (a, b + 1)):
            if db - da + 1 > n - 4:
                continue
            idx = [(i0 + k) % n for k in range(da, db + 1)]
            c, dev = fit(idx)
            if dev <= tol:
                a, b, best, grown = da, db, (c, da, db, dev), True
                break
        if not grown:
            break
    if best is None or best[2] - best[1] + 1 < 3:
        return None
    return best


def arc_points(c, pa, pb, R, ds):
    """中心 c・半径 R の円上を、pa の方位から pb の方位へ (短い側) 刻む。端点も円上に置く。"""
    th_a = np.arctan2(pa[1] - c[1], pa[0] - c[0])
    th_b = np.arctan2(pb[1] - c[1], pb[0] - c[0])
    dth = (th_b - th_a) % (2 * np.pi)
    if dth > np.pi:
        dth -= 2 * np.pi
    m = max(3, int(abs(dth) * R / ds) + 1)
    th = th_a + dth * np.linspace(0.0, 1.0, m)
    return np.column_stack([c[0] + R * np.cos(th), c[1] + R * np.sin(th)])


def build_profile_arcs(vane, tab, ds=0.02):
    """表点 + **接線が繋がる**前縁/後縁円弧の閉輪郭を作り、(点列, 円弧マスク) を返す。"""
    pts = [np.array(p) for p in tab]
    is_arc = [False] * len(pts)
    for edge in ("TE", "LE"):
        R = RADII[vane][edge]
        if R <= 0:
            continue
        cur = np.array(pts)
        i0 = int(np.argmin(cur[:, 0])) if edge == "LE" else int(np.argmax(cur[:, 0]))
        got = fit_arc_span(cur, i0, R)
        if got is None:
            print(f"  [{edge}] 円弧区間が見つからない (tol {ARC_TOL} cm) — 表点のまま")
            continue
        c, a, b, dev = got
        n = len(cur)
        idx = [(i0 + k) % n for k in range(a, b + 1)]
        pa = c + (cur[idx[0]] - c) / np.linalg.norm(cur[idx[0]] - c) * R     # 端点を円上へ投影
        pb = c + (cur[idx[-1]] - c) / np.linalg.norm(cur[idx[-1]] - c) * R
        arc = arc_points(c, pa, pb, R, ds)
        print(f"  [{edge}] R={R} cm, center=({c[0]:.4f},{c[1]:.4f}), span {len(idx)} table points, "
              f"fit dev {dev:.5f} cm -> arc {len(arc)} points")
        keep = [k for k in range(n) if k not in set(idx)]
        # idx は巡回区間なので、残す側を順序どおりに繋ぐ
        start = (idx[-1] + 1) % n
        order = [(start + k) % n for k in range(n - len(idx))]
        new_pts = [cur[k] for k in order] + [p for p in arc]
        new_arc = [is_arc[k] for k in order] + [True] * len(arc)
        pts, is_arc = new_pts, new_arc
    return np.array(pts), np.array(is_arc)


def load_table(vane):
    f = REF[vane] / f"vane_{vane}.csv"
    pts = []
    for line in open(f):
        if line.startswith("#") or line.lstrip().startswith("i,") or not line.strip():
            continue
        p = line.split(",")
        pts.append((float(p[1]), float(p[2])))
    return np.array(pts)


def turn_angles(P):
    t = np.diff(np.vstack([P, P[:1]]), axis=0)
    L = np.hypot(*t.T)
    t = t / L[:, None]
    return np.degrees(np.arccos(np.clip(np.sum(t * np.roll(t, 1, 0), axis=1), -1, 1))), L


def self_intersects(P):
    """閉多角形の自己交差 (隣接辺は除く)。点数が数百なので素朴な O(n²) で十分。"""
    n = len(P)
    A = P
    B = np.roll(P, -1, axis=0)

    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    for i in range(n):
        for j in range(i + 2, n):
            if i == 0 and j == n - 1:
                continue
            d1 = cross(A[i], B[i], A[j]); d2 = cross(A[i], B[i], B[j])
            d3 = cross(A[j], B[j], A[i]); d4 = cross(A[j], B[j], B[i])
            if ((d1 > 0) != (d2 > 0)) and ((d3 > 0) != (d4 > 0)):
                return (i, j)
    return None


def curvature_radius(P, i, half=6):
    """節点 i まわりの曲率半径 [cm] (前後 half 点の円あてはめ)。"""
    n = len(P)
    idx = [(i + k) % n for k in range(-half, half + 1)]
    Q = P[idx]
    A = np.column_stack([2 * Q[:, 0], 2 * Q[:, 1], np.ones(len(Q))])
    b = Q[:, 0] ** 2 + Q[:, 1] ** 2
    c, *_ = np.linalg.lstsq(A, b, rcond=None)
    return float(np.sqrt(max(c[2] + c[0] ** 2 + c[1] ** 2, 0.0)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vane", default="c3x", choices=["c3x", "markii"])
    ap.add_argument("--n", type=int, default=3200, help="書き出す点数 (メッシュ側でさらに間引く)")
    ap.add_argument("--max-turn", type=float, default=5.0, help="再標本化後の折れ角の上限 [deg]")
    ap.add_argument("--smooth-scale", type=float, default=0.05, help="平滑化量 s = scale * m")
    a = ap.parse_args()

    tab = load_table(a.vane)           # 表の点だけ
    P0, is_arc = build_profile_arcs(a.vane, tab)
    is_tab = ~is_arc
    # **検査に使うのは輪郭に残っている表点だけ**。円弧区間に飲まれた点を基準にすると、
    # 円弧を採用した判断そのものを不合格にしてしまう。
    tab_kept = P0[is_tab]
    w = np.where(is_tab, 1.0 / SIG_TABLE, 1.0 / SIG_ARC)
    m = len(P0)
    ang0, L0 = turn_angles(P0)
    print(f"[{a.vane}] input {m} points ({is_tab.sum()} table / {m-is_tab.sum()} arc), "
          f"turn max {ang0.max():.1f}° p90 {np.percentile(ang0,90):.1f}°, "
          f"seg {L0.min()*10:.2f}–{L0.max()*10:.2f} mm")

    tck, u = splprep([P0[:, 0], P0[:, 1]], w=w, s=a.smooth_scale * m, per=True, k=3)
    # 弦長等間隔に近づけるため、細かく評価してから弧長で再標本化
    uu = np.linspace(0, 1, 20000, endpoint=False)
    X, Y = splev(uu, tck)
    Q = np.column_stack([X, Y])
    seg = np.hypot(*np.diff(np.vstack([Q, Q[:1]]), axis=0).T)
    s_cum = np.concatenate([[0.0], np.cumsum(seg)])
    total = s_cum[-1]
    tgt = np.linspace(0, total, a.n, endpoint=False)
    P = np.column_stack([np.interp(tgt, s_cum[:-1], Q[:, 0]),
                         np.interp(tgt, s_cum[:-1], Q[:, 1])])

    # ---- 検査 ----
    dev = np.array([np.min(np.hypot(*(Q - p).T)) for p in tab_kept])
    ang, L = turn_angles(P)
    i_le = int(np.argmin(P[:, 0])); i_te = int(np.argmax(P[:, 0]))
    # 円あてはめの窓は**期待半径に合わせる** (窓が半径より広いと丸みを過大に読む)
    ds = total / a.n
    R_le = curvature_radius(P, i_le, half=max(3, int(0.5 * RADII[a.vane]["LE"] / ds)))
    R_te = curvature_radius(P, i_te, half=max(3, int(0.5 * max(RADII[a.vane]["TE"], ds * 6) / ds)))
    xsec = self_intersects(P)
    fails = []
    print(f"[{a.vane}] smoothed: arc {total:.3f} cm -> {a.n} points ({total/a.n*10:.3f} mm)")
    print(f"  deviation from kept table points ({len(tab_kept)}): rms {dev.mean():.5f} cm, max {dev.max():.5f} cm "
          f"(uncertainty {SIG_TABLE} cm)")
    print(f"  turn angle: mean {ang.mean():.2f}° p90 {np.percentile(ang,90):.2f}° max {ang.max():.2f}°")
    print(f"  curvature radius: LE {R_le:.3f} cm (ref {RADII[a.vane]['LE']}), "
          f"TE {R_te:.3f} cm (ref {RADII[a.vane]['TE']})")
    if np.sqrt((dev ** 2).mean()) > SIG_TABLE:
        fails.append(f"deviation rms {np.sqrt((dev**2).mean()):.5f} > {SIG_TABLE} cm")
    if dev.max() > 3 * SIG_TABLE:
        fails.append(f"deviation max {dev.max():.5f} > {3*SIG_TABLE} cm")
    if ang.max() > a.max_turn:
        fails.append(f"turn angle max {ang.max():.2f}° > {a.max_turn}°")
    for nm, R, Rr in (("LE", R_le, RADII[a.vane]["LE"]), ("TE", R_te, RADII[a.vane]["TE"])):
        if Rr > 0 and abs(R - Rr) > 0.20 * Rr:
            fails.append(f"{nm} radius {R:.3f} cm vs {Rr} cm (>20 %)")
    if xsec:
        fails.append(f"self-intersection between segments {xsec}")
    if fails:
        print("VERDICT: FAIL")
        for f in fails:
            print("  -", f)
        sys.exit(1)

    out = REF[a.vane] / f"vane_{a.vane}_smooth.csv"
    with open(out, "w") as f:
        f.write(f"# {a.vane} vane — 平滑化輪郭 (周期スプライン, s={a.smooth_scale}*m)\n"
                f"# 入力: vane_{a.vane}_profile.csv (表の点 + 復元円弧)\n"
                f"# 検査: 表点からの max {dev.max():.5f} cm / 折れ角 max {ang.max():.2f}° / "
                f"R_LE {R_le:.3f} R_TE {R_te:.3f} cm\n")
        f.write("i,x_cm,y_cm\n")
        for i, (x, y) in enumerate(P, 1):
            f.write(f"{i},{x:.5f},{y:.5f}\n")
    print(f"VERDICT: PASS  -> {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
