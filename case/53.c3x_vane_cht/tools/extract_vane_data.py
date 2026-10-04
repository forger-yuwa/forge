#!/usr/bin/env python3
r"""NASA CR-168015 (Hylton et al. 1983) から**翼型座標と試験条件**を抽出する。

一次資料: `papers/cht/NASA-CR-168015_Hylton_1983.pdf` (NTRS 19830020105, **git 追跡外**)。
スキャン PDF の OCR なので、**数値は必ず検証してから書き出す**:

1. 表 II (Mark II) / III (C3X) は **cm と inch の両方**を載せているので、
   `cm == 2.54 * in` (相対 0.5 %) を全点で課す。**通らない点は捨てて REPAIRS から補う**。
2. `REPAIRS` は**ページ画像を目視して読んだ値**で、各行を同じ 2.54 則で検算してから入れてある。
3. 点数が揃わなければ**書き出さずに落ちる** (欠測のまま形状を作らない)。

出力:
  case/53.c3x_vane_cht/ref/vane_c3x.csv       78 点 (x_cm, y_cm)
  case/54.markii_vane_cht/ref/vane_markii.csv 60 点
  case/53.c3x_vane_cht/ref/test_conditions.csv 表 VIII/IX の主要行

usage: python3 case/53.c3x_vane_cht/tools/extract_vane_data.py [--pdf <path>] [--plot]
"""
import argparse
import math
import os
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]

# ---- 目視で読んだ補修値 (OCR が壊した行だけ)。各行 cm == 2.54*in を満たすことを確認済み ----
REPAIRS = {
    "markii": {   # page 11 (PDF 17)
        7:  (2.9812, 11.0833),   # (1.1737) (4.3635)
        22: (5.7269, 3.9444),    # (2.2547) (1.5529)
        24: (6.0295, 2.9741),    # (2.3738) (1.1709)
        27: (6.4554, 1.5128),    # (2.5415) (0.5956)
        57: (1.6337, 8.9891),    # (0.6432) (3.5390)
        59: (1.0208, 9.5456),    # (0.4019) (3.7581)
    },
    "c3x": {      # page 12 (PDF 18)
        # **報告自身が不整合**: 29 行は "0.4115 (0.0162)" と印字されているが 0.0162 in = 0.0411 cm。
        # **後縁円で決まる**: 点 30 (7.8161,-0.0053) が後縁頂点なので円中心は (7.8161-R_TE, -0.0053)
        # = (7.6431,-0.0053)。中心からの距離は (7.8115, 0.0411) で **0.174 cm = R_TE 0.173 に一致**、
        # (7.8115, 0.4115) では 0.449 cm で円から外れる。→ **inch 値 0.0162 が正、cm 側が誤植**。
        # (当初は点間隔から cm 側を採ったが、後縁円の検査で誤りと判明した。2026-09-20)
        29: (7.8115, 0.0411),
        11: (4.2885, 10.5766),   # (1.6884) (4.1640)
        65: (1.5519, 9.3932),    # (0.6110) (3.6981)
        70: (0.2822, 10.4094),   # (0.1111) (4.0982)
        75: (0.0064, 11.0802),   # (0.0025) (4.3623)
    },
}
NPTS = {"markii": 60, "c3x": 78}
PAGES = {"markii": 17, "c3x": 18}


def norm(t):
    t = t.replace('"', '.').replace('·', '.')
    t = re.sub(r'(?<=\d),(?=\d)', '.', t)
    t = re.sub(r'\bi0\b', '10', t)
    t = t.replace('II.', '11.').replace('(I.', '(1.').replace('(i.', '(1.')
    return t


def extract_points(doc, name):
    body = norm(doc[PAGES[name]].get_text())
    # **表ヘッダを捨てる**: "R_LE = 1.168 cm (0.460 in.) R_TE = 0.173 cm (0.068 in.)" は
    # cm == 2.54*in を満たすので、そのままだとデータ行として通ってしまう (実際に点 3 を潰した)。
    # 列見出し ("Position number" が 2 回) の後ろだけを読む。
    k = body.rfind("number")
    if k > 0:
        body = body[k:]
    nums = re.findall(r'-?\d+\.\d+|-?\d+', body)
    # 位置番号は 1..N で重複なし、かつ cm == 2.54*in を課す。番号の逐次チェックまでやると
    # OCR が 1 行壊したときに同期を失うので、ここは**重複禁止**までに留め、欠測は REPAIRS で埋める。
    pts, i = {}, 0
    while i + 4 < len(nums):
        try:
            pos = int(float(nums[i])); xc = float(nums[i + 1]); xi = float(nums[i + 2])
            yc = float(nums[i + 3]); yi = float(nums[i + 4])
        except ValueError:
            i += 1; continue
        okx = abs(xc - 2.54 * xi) <= max(0.01, 0.005 * abs(xc))
        oky = abs(yc - 2.54 * yi) <= max(0.01, 0.005 * abs(yc))
        if 1 <= pos <= NPTS[name] and okx and oky and pos not in pts:
            pts[pos] = (xc, yc); i += 5
        else:
            i += 1
    n_ocr = len(pts)
    for pos, xy in REPAIRS[name].items():
        pts.setdefault(pos, xy)
    missing = sorted(set(range(1, NPTS[name] + 1)) - set(pts))
    print(f"[{name}] OCR validated {n_ocr}/{NPTS[name]}, repaired {len(REPAIRS[name])}, missing {missing}")
    if missing:
        sys.exit(f"[{name}] {len(missing)} points missing — ページ画像を読んで REPAIRS に足すこと (推測で埋めない)")
    return [pts[i] for i in range(1, NPTS[name] + 1)]


RADII = {"c3x": {"LE": 1.168, "TE": 0.173}, "markii": {"LE": 1.280, "TE": 0.0}}


def check_round_edges(name, pts, tol=0.06):
    """前縁・後縁の丸みが表の点と整合するかを検査する (OCR 事故の最後の砦)。

    最も左 (前縁) / 最も右 (後縁) の点を頂点とみなし、内側へ R だけ入った点を円中心とする。
    R_LE / R_TE は表の見出しの値。**前縁と後縁で見方を変える**:

    - **後縁 (小半径)**: 近傍 (2R 以内) の点は**すべて円上**にあるはず。点列は後縁の丸みを
      刻んでいるので、外れる点は誤読。→ 実際にこれで C3X 点 29 の cm 側の誤植を捕まえた。
    - **前縁 (大半径)**: 翼面は円弧からすぐ離れるので「近傍なら円上」は課せない。代わりに
      **円の内側に入る点が無いこと** (翼面が鼻の円を食い込むことはない) と、
      **円上に乗る点が 3 点以上あること** (円弧が点列と繋がっている) を課す。
    """
    xs = [q[0] for q in pts]; ys = [q[1] for q in pts]
    bad_all = []
    for edge, sgn in (("LE", +1), ("TE", -1)):
        R = RADII[name][edge]
        if R <= 0:
            print(f"[{name}] {edge}: R=0 (blunt) — 円の検査は省略")
            continue
        i = (xs.index(min(xs)) if edge == "LE" else xs.index(max(xs)))
        cx, cy = xs[i] + sgn * R, ys[i]
        d = [math.hypot(x - cx, y - cy) for x, y in zip(xs, ys)]
        if edge == "TE":
            near = [(k, d[k]) for k in range(len(pts)) if d[k] < 2.0 * R]
            bad = [(k + 1, round(v, 4)) for k, v in near if abs(v - R) > tol * R]
            print(f"[{name}] TE circle R={R} cm at ({cx:.4f},{cy:.4f}): {len(near)} points nearby, "
                  f"{len(bad)} off the circle {bad if bad else ''}")
            bad_all += bad
        else:
            inside = [(k + 1, round(d[k], 4)) for k in range(len(pts)) if d[k] < (1.0 - tol) * R]
            on = [k for k in range(len(pts)) if abs(d[k] - R) <= tol * R]
            print(f"[{name}] LE circle R={R} cm at ({cx:.4f},{cy:.4f}): {len(on)} points on the circle, "
                  f"{len(inside)} inside it {inside if inside else ''}")
            bad_all += inside
            if len(on) < 3:
                bad_all.append(("LE", f"only {len(on)} points on the circle"))
    if bad_all:
        sys.exit(f"[{name}] 前縁/後縁の円と合わない点がある: {bad_all} — OCR か補修値を疑うこと")



# ---- 冷却孔 (図 6 = Mark II / 図 7 = C3X)。**図なので OCR テキストが無く、ページ画像から転記した**。
# 各行は (U_cm, U_in, V_cm, V_in, D_cm, D_in, Cr) で、書き出し時に cm == 2.54*in を機械検査する
# (転記ミスをそこで捕まえる)。U,V は **翼弦に沿った FE モデルの座標系** (図の U は厚み方向、V は前縁→後縁)。
# Cr は冷却孔 Nu の**助走区間補正係数** (報告 p.21 の手順)。
COOLING_HOLES = {
    "markii": [
        (1.448, 0.570,  0.711, 0.280, 0.630, 0.248, 1.118),
        (1.016, 0.400,  1.930, 0.760, 0.630, 0.248, 1.118),
        (2.083, 0.820,  1.854, 0.730, 0.630, 0.248, 1.118),
        (1.676, 0.660,  3.556, 1.400, 0.630, 0.248, 1.118),
        (1.524, 0.600,  5.182, 2.040, 0.630, 0.248, 1.118),
        (1.397, 0.550,  6.807, 2.680, 0.630, 0.248, 1.118),
        (1.143, 0.450,  8.433, 3.320, 0.630, 0.248, 1.118),
        (0.864, 0.340,  9.957, 3.920, 0.310, 0.122, 1.056),
        (0.635, 0.250, 11.303, 4.450, 0.310, 0.122, 1.056),
        (0.381, 0.150, 12.497, 4.920, 0.198, 0.078, 1.025),
    ],
    "c3x": [
        (2.377, 0.936,  1.311, 0.516, 0.630, 0.248, 1.118),
        (1.057, 0.416,  1.534, 0.604, 0.630, 0.248, 1.118),
        (1.981, 0.780,  3.119, 1.228, 0.630, 0.248, 1.118),
        (1.981, 0.780,  4.674, 1.840, 0.630, 0.248, 1.118),
        (1.869, 0.736,  6.182, 2.434, 0.630, 0.248, 1.118),
        (1.666, 0.656,  7.747, 3.050, 0.630, 0.248, 1.118),
        (1.412, 0.556,  9.235, 3.636, 0.630, 0.248, 1.118),
        (1.087, 0.428, 10.759, 4.236, 0.310, 0.122, 1.056),
        (0.737, 0.290, 12.253, 4.824, 0.310, 0.122, 1.056),
        (0.345, 0.136, 13.757, 5.416, 0.198, 0.078, 1.025),
    ],
}


def place_holes(name, profile, holes):
    r"""冷却孔の (U,V) を翼型の (x,y) 座標系へ写す。

    図 6/7 の (U,V) は FE モデルの座標系で、**変換行列は報告に書かれていない**。そこで
    **剛体変換 (回転 θ + 平行移動 + 鏡映の 2 通り) を、全孔が翼型内に入り最小肉厚が最大になるように
    最適化して決める**。初期値は「孔 1→孔 10」の向きを「前縁頂点→後縁頂点」に合わせたもの。
    最小肉厚が非正なら**落ちる** (推測で置かない)。肉厚は報告の不確かさ評価でも鍵になる量なので、
    求めた値をログと CSV に残す。
    """
    from scipy.optimize import minimize

    P = np.array(profile, float)
    UV = np.array([[h[0], h[2]] for h in holes], float)      # (U_cm, V_cm)
    rad = np.array([0.5 * h[4] for h in holes], float)       # 半径 [cm]
    le = P[int(np.argmin(P[:, 0]))]
    te = P[int(np.argmax(P[:, 0]))]

    seg_a = P
    seg_b = np.roll(P, -1, axis=0)

    def clearance(q):
        """点 q から輪郭までの最短距離 (線分距離)。"""
        d = seg_b - seg_a
        t = np.clip(((q - seg_a) * d).sum(1) / np.maximum((d * d).sum(1), 1e-30), 0.0, 1.0)
        proj = seg_a + t[:, None] * d
        return float(np.hypot(*(q - proj).T).min())

    def inside(q):
        x, y = q; c = False
        for i in range(len(P)):
            x1, y1 = P[i]; x2, y2 = P[(i + 1) % len(P)]
            if ((y1 > y) != (y2 > y)) and (x < (x2 - x1) * (y - y1) / (y2 - y1) + x1):
                c = not c
        return c

    def place(par, mirror):
        th, tx, ty = par
        R = np.array([[math.cos(th), -math.sin(th)], [math.sin(th), math.cos(th)]])
        M = np.array([[mirror, 0.0], [0.0, 1.0]])
        return (R @ (M @ UV.T)).T + np.array([tx, ty])

    def neg_min_wall(par, mirror):
        pos = place(par, mirror)
        w = [clearance(q) - r if inside(q) else -(clearance(q) + r) for q, r in zip(pos, rad)]
        return -min(w)

    # 初期値: 孔 1 -> 孔 10 の向きを 前縁 -> 後縁 に合わせる
    v_uv = UV[-1] - UV[0]
    v_xy = te - le
    th0 = math.atan2(v_xy[1], v_xy[0]) - math.atan2(v_uv[1], v_uv[0])
    best = None
    for mirror in (+1, -1):
        for dth in (0.0, 0.1, -0.1, 0.25, -0.25):
            R0 = np.array([[math.cos(th0 + dth), -math.sin(th0 + dth)],
                           [math.sin(th0 + dth), math.cos(th0 + dth)]])
            M0 = np.array([[mirror, 0.0], [0.0, 1.0]])
            t0 = le - (R0 @ (M0 @ UV[0]))
            res = minimize(neg_min_wall, [th0 + dth, t0[0], t0[1]], args=(mirror,),
                           method="Nelder-Mead",
                           options={"xatol": 1e-6, "fatol": 1e-8, "maxiter": 4000, "maxfev": 6000})
            if best is None or res.fun < best[0]:
                best = (res.fun, res.x, mirror)
    wall = -best[0]
    pos = place(best[1], best[2])
    th = math.degrees(best[1][0]) % 360.0
    print(f"[{name}] cooling holes placed: rotation {th:.2f} deg, mirror {best[2]:+d}, "
          f"**min wall thickness {wall:.4f} cm** ({wall/2.54:.4f} in)")
    if wall <= 0.02:
        sys.exit(f"[{name}] 最小肉厚が {wall:.4f} cm と小さすぎる — (U,V) の解釈か転記を見直すこと")
    return pos, wall


def write_holes(name, path, profile=None):
    """冷却孔表を検査して書き出す (cm == 2.54*in を全項目に課す)。

    profile を渡すと (U,V) を翼型座標 (x,y) に写した列も足す (place_holes)。
    """
    rows = COOLING_HOLES[name]
    for k, (uc, ui, vc, vi, dc, di, cr) in enumerate(rows, 1):
        for cm, inch, what in ((uc, ui, "U"), (vc, vi, "V"), (dc, di, "D")):
            if abs(cm - 2.54 * inch) > max(0.002, 0.005 * abs(cm)):
                sys.exit(f"[{name}] cooling hole {k} {what}: {cm} cm vs {inch} in — 転記ミス")
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        xy, clear = (place_holes(name, profile, rows) if profile is not None else (None, None))
        f.write(f"# {name} cooling holes — NASA CR-168015 Figure {'6' if name=='markii' else '7'} "
                f"(report p.{'15' if name=='markii' else '16'}), ページ画像から転記\n"
                "# U,V は FE モデルの翼弦座標系 [cm] (U=厚み方向, V=前縁→後縁)。D=直径 [cm]。\n"
                "# Cr = 冷却孔 Nusselt 数の助走区間補正係数 (報告の処理手順)\n"
                + ("# x_cm,y_cm は前縁頂点を原点・前縁→後縁を V 軸とした剛体変換で置いたもの "
                   f"(最小肉厚 {clear:.4f} cm)\n" if xy is not None else "")
                + ("hole,U_cm,V_cm,D_cm,Cr,x_cm,y_cm\n" if xy is not None else "hole,U_cm,V_cm,D_cm,Cr\n"))
        for k, (uc, ui, vc, vi, dc, di, cr) in enumerate(rows, 1):
            if xy is not None:
                f.write(f"{k},{uc:.3f},{vc:.3f},{dc:.3f},{cr:.3f},{xy[k-1][0]:.4f},{xy[k-1][1]:.4f}\n")
            else:
                f.write(f"{k},{uc:.3f},{vc:.3f},{dc:.3f},{cr:.3f}\n")
    print(f"  -> {path.relative_to(ROOT)} (10 holes, cm/in checked)")


def build_profile(name, pts, n_arc=60):
    """表の点列に**前縁・後縁の円弧を復元**して閉じた輪郭を作る。

    表は丸み部分を 2〜3 点しか刻まないので、そのまま結ぶと**鼻と後縁が角張る**
    (Mark II 前縁は 1 区間で 1.63 cm 飛ぶ)。R_LE / R_TE は見出しにあるので、
    **円上に乗っている点の連続区間**を見つけ、その両端の間を円弧で置き換える。
    メッシュ生成はこの輪郭を使う。
    """
    cur = [(float(q[0]), float(q[1])) for q in pts]
    for edge, sgn in (("TE", -1), ("LE", +1)):
        R = RADII[name][edge]
        if R <= 0:
            continue
        xs = np.array([q[0] for q in cur]); ys = np.array([q[1] for q in cur])
        i = int(np.argmin(xs)) if edge == "LE" else int(np.argmax(xs))
        cx, cy = xs[i] + sgn * R, ys[i]
        d = np.hypot(xs - cx, ys - cy)
        on = np.abs(d - R) <= 0.06 * R
        n = len(cur)
        if on.sum() < 2:
            continue
        # 円上の点の**巡回的に最長の連続区間** [a..b] を取る
        best = (0, -1, -1)
        for a in range(n):
            if not on[a] or (on[a - 1] and a > 0):      # 区間の先頭だけ見る
                continue
            L, b = 1, a
            while L < n and on[(a + L) % n]:
                b = (a + L) % n; L += 1
            if L > best[0]:
                best = (L, a, b)
        L, a, b = best
        if L < 2:
            continue
        th_a = math.atan2(ys[a] - cy, xs[a] - cx)
        th_b = math.atan2(ys[b] - cy, xs[b] - cx)
        dth = (th_b - th_a) % (2 * math.pi)
        if dth > math.pi:
            dth -= 2 * math.pi                          # 短い側の弧
        arc = [(cx + R * math.cos(th_a + dth * t), cy + R * math.sin(th_a + dth * t))
               for t in np.linspace(0.0, 1.0, n_arc)[1:-1]]
        if a <= b:
            cur = cur[:a + 1] + arc + cur[b:]
        else:                                           # 末尾→先頭をまたぐ区間
            cur = cur[b:a + 1] + arc
    return cur


def write_csv(path, pts, header):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        f.write(header)
        f.write("i,x_cm,y_cm\n")
        for i, (x, y) in enumerate(pts, 1):
            f.write(f"{i},{x:.4f},{y:.4f}\n")
    print(f"  -> {path.relative_to(ROOT)} ({len(pts)} points)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdf", default=str(ROOT / "papers/cht/NASA-CR-168015_Hylton_1983.pdf"))
    ap.add_argument("--plot", action="store_true")
    a = ap.parse_args()
    import fitz
    doc = fitz.open(a.pdf)

    meta = {
        "c3x": ("# C3X vane coordinates — NASA CR-168015 Table III (report p.12)\n"
                "# R_LE = 1.168 cm (0.460 in), R_TE = 0.173 cm (0.068 in)\n"
                "# 検証: 全点で cm == 2.54*in (相対 0.5%)。点 29 は報告の inch 値が誤植 (抽出スクリプト参照)\n"),
        "markii": ("# Mark II vane coordinates — NASA CR-168015 Table II (report p.11)\n"
                   "# R_LE = 1.280 cm (0.504 in), R_TE = 0.000 (blunt)\n"
                   "# 検証: 全点で cm == 2.54*in (相対 0.5%)\n"),
    }
    out = {"c3x": ROOT / "case/53.c3x_vane_cht/ref/vane_c3x.csv",
           "markii": ROOT / "case/54.markii_vane_cht/ref/vane_markii.csv"}
    pts = {}
    profiles = {}
    for name in ("c3x", "markii"):
        pts[name] = extract_points(doc, name)
        check_round_edges(name, pts[name])
        write_csv(out[name], pts[name], meta[name])
        prof = build_profile(name, pts[name])
        pp = out[name].with_name(out[name].stem + "_profile.csv")
        with open(pp, "w") as f:
            f.write(meta[name].replace("coordinates", "closed profile")
                    + "# 表の点 + R_LE / R_TE から復元した円弧 (メッシュ生成はこちらを使う)\n")
            f.write("i,x_cm,y_cm\n")
            for k, (x, y) in enumerate(prof, 1):
                f.write(f"{k},{x:.5f},{y:.5f}\n")
        print(f"  -> {pp.relative_to(ROOT)} ({len(prof)} points, arcs reconstructed)")
        profiles[name] = prof

    # 試験条件 (表 VIII/IX, report p.30 = PDF 36) — **ページ画像から読んだ値**。
    # 注意: 表の SI 列 ("PT1--Pa") は psia 列と 51.7 倍ずれており内部矛盾している。
    #       psia 側が文献の引用値 (C3X 4411 で ~3.2 atm) と合うので **psia を正**とする。
    write_holes("c3x", ROOT / "case/53.c3x_vane_cht/ref/cooling_holes_c3x.csv", profiles["c3x"])
    write_holes("markii", ROOT / "case/54.markii_vane_cht/ref/cooling_holes_markii.csv", profiles["markii"])

    cond = ROOT / "case/53.c3x_vane_cht/ref/test_conditions.csv"
    with open(cond, "w") as f:
        f.write("# NASA CR-168015 Table VIII (Mark II) / IX (C3X), report p.30 — ページ画像から読取\n"
                "# PT1 は psia 列を正とする (表の SI 列は psia と 51.7 倍不整合)\n"
                "vane,code,run,PT1_psia,PT1_kPa,TT1_K,M1,Re1_e6,M2,Re2_e6,Tu_pct,Tw_over_Tg\n"
                "markii,5411,42,48.89,337.1,788,0.19,0.56,1.04,2.01,6.5,0.68\n"
                "markii,4411,43,49.64,342.3,784,0.18,0.57,0.89,1.98,6.5,0.69\n"
                "c3x,4411,108,46.34,319.5,786,0.17,0.52,0.90,1.99,6.5,0.73\n"
                "c3x,5411,107,45.24,311.9,798,0.17,0.51,1.05,1.97,6.5,0.72\n")
    print(f"  -> {cond.relative_to(ROOT)}")

    if a.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(11, 5))
        for k, name in enumerate(("c3x", "markii")):
            p = pts[name]
            pr = profiles[name]
            ax[k].plot([q[0] for q in pr] + [pr[0][0]], [q[1] for q in pr] + [pr[0][1]], "-",
                       lw=1.2, color="0.4", label="profile (arcs reconstructed)")
            ax[k].plot([q[0] for q in p] + [p[0][0]], [q[1] for q in p] + [p[0][1]], "-o", ms=2.5,
                       label="table points")
            # **前縁/後縁の円弧は表に点が無い** (Mark II 前縁は 1 点 1.63 cm 飛ぶ) ので、
            # R_LE / R_TE から描き足す。メッシュ生成でも同じ扱いにすること。
            import numpy as _np
            xs = [q[0] for q in p]; ys = [q[1] for q in p]
            for edge, sgn in (("LE", +1), ("TE", -1)):
                R = RADII[name][edge]
                if R <= 0:
                    continue
                i = (xs.index(min(xs)) if edge == "LE" else xs.index(max(xs)))
                cx, cy = xs[i] + sgn * R, ys[i]
                th = _np.linspace(0, 2 * math.pi, 200)
                ax[k].plot(cx + R * _np.cos(th), cy + R * _np.sin(th), "--", lw=0.8,
                           label=f"{edge} circle R={R} cm")
            # 冷却孔 (place_holes で置いた位置と直径) も描いて、肉厚を目で確認できるようにする
            hp = ROOT / (f"case/53.c3x_vane_cht/ref/cooling_holes_c3x.csv" if name == "c3x"
                         else "case/54.markii_vane_cht/ref/cooling_holes_markii.csv")
            if hp.exists():
                import csv as _csv
                rows = [r for r in _csv.DictReader(l for l in open(hp) if not l.startswith("#"))]
                if rows and "x_cm" in rows[0]:
                    for r in rows:
                        cx0, cy0, D = float(r["x_cm"]), float(r["y_cm"]), float(r["D_cm"])
                        th = np.linspace(0, 2 * math.pi, 60)
                        ax[k].plot(cx0 + 0.5 * D * np.cos(th), cy0 + 0.5 * D * np.sin(th), "-",
                                   color="tab:red", lw=0.9)
                    ax[k].plot([], [], "-", color="tab:red", lw=0.9, label="cooling holes")
            ax[k].legend(fontsize=7)
            ax[k].set_title(f"{name} ({len(p)} pts)"); ax[k].set_aspect("equal"); ax[k].grid(alpha=.3)
            ax[k].set_xlabel("x [cm]"); ax[k].set_ylabel("y [cm]")
        fig.tight_layout()
        png = ROOT / "case/53.c3x_vane_cht/ref/vane_shapes.png"
        fig.savefig(png, dpi=130)
        print(f"  -> {png.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
