#!/usr/bin/env python3
"""case/49 の**全ヘキサ (構造化) メッシュ**を gmsh Python API で作る。

usage (mesh venv):
  CASE49_MANIFEST=manifest.json python3 cad/build_hex_mesh.py [--out cad/hex.msh] [--show-sizes]

考え方: **z=0 面を 2D 四角形でブロッキングし、上下に押し出す**だけで全域ヘキサになる。

    領域                        ブロック                        押し出し
    すきま (環状スリット)        4 セクタ × 半径 2 層 = 8 面      下向き (深さ方向、床側と開口側を細かく)
    内円柱の上面側               バタフライ (コア + 3 面)         上向き
    開口まわり                   円→角の O グリッド 4 面          上向き
    遠方                         矩形 3 面                        上向き

環状面は上下の押し出しで共有されるので開口は**内部面**になり、自動的に整合する。

tet+prism 版に対する利点 (2026-09-19 ユーザ質問への回答):
  - **レイヤーと空間メッシュの継ぎ目が無い** (半径方向の分布を直接指定する)
  - **「VL 総厚 ≤ 接線セルサイズ」の制約が消える** (層が自己交差しない)。内面をいくらでも細かくできる
  - すきま横断に両壁から自由に配れる (第一セルを両側で指定)
  - 同解像度での節点数が減る

偏心 (`x_off`) は、外筒軸からのレイが内円と交わる半径 r_i(φ) を使うだけでそのまま扱える
(トポロジは変わらないので、すきま量・深さ・マッハ数と同じくパラメータ)。
"""
import argparse
import json
import math
import os
import sys
from pathlib import Path

import gmsh

HERE = Path(__file__).resolve().parent
CASE = HERE.parent


# ---------------------------------------------------------------- 分布ユーティリティ
def progression_ratio(L, n, h1):
    """長さ L を n セル、第一セル h1 の等比で刻むときの比 p。"""
    if n <= 1:
        return 1.0
    if abs(h1 * n - L) < 1e-15 * L:
        return 1.0
    lo, hi = 1.0 + 1e-9, 10.0
    for _ in range(200):
        p = 0.5 * (lo + hi)
        s = h1 * (p ** n - 1.0) / (p - 1.0)
        if s < L:
            lo = p
        else:
            hi = p
    return 0.5 * (lo + hi)


def sizes_one_sided(L, n, h1):
    """片側 (始点側) を細かくする等比セル長の列。"""
    p = progression_ratio(L, n, h1)
    s = [h1 * p ** i for i in range(n)]
    f = L / sum(s)
    return [x * f for x in s]


def sizes_two_sided(L, n, h1, h2=None):
    """両端を細かくする列 (中央で最大)。n は偶数に丸める。"""
    h2 = h1 if h2 is None else h2
    n = max(2, n - (n % 2))
    a = sizes_one_sided(0.5 * L, n // 2, h1)
    b = sizes_one_sided(0.5 * L, n // 2, h2)
    return a + b[::-1]


def layers_from_sizes(sizes):
    """gmsh の Extrude Layers 用 (要素数リスト, 累積比リスト)。"""
    tot = sum(sizes)
    cum, s = [], 0.0
    for h in sizes:
        s += h
        cum.append(s / tot)
    return [1] * len(sizes), cum


# ---------------------------------------------------------------- 本体
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(HERE / "hex.msh"))
    ap.add_argument("--manifest", default=os.environ.get("CASE49_MANIFEST", "manifest.json"))
    ap.add_argument("--show-sizes", action="store_true")
    ap.add_argument("--no-mesh", action="store_true", help="形状だけ作って終わる (確認用)")
    ap.add_argument("--scale", type=float, default=1.0,
                    help="全方向の分割数を一律に掛ける (系統的な格子細分列を作る。第一セルは 1/scale)")
    a = ap.parse_args()

    MAN = json.loads((CASE / a.manifest).read_text())
    G, M, PID = MAN["geometry"], MAN["mesh"], MAN["phys_id"]
    if not G.get("half_model", True):
        raise SystemExit("いまは半割のみ対応 (全周は phi を 0..2pi に広げれば同じ構成で作れる)")
    if G.get("plug_cavity", False):
        raise SystemExit("塞ぎ形状は未対応")

    Ro, Ri, off, dep = G["Ro"], G["Ri"], G["x_off"], G["depth"]
    x0, x1, ym, zt = G["x_in"], G["x_out"], G["y_max"], G["z_top"]
    H = dict(M.get("hex", {}))
    # --- 系統細分: **全方向の分割数を一律 scale 倍、第一セルを 1/scale 倍** ---
    # これができるのがヘキサ (構造化) の利点。tet+prism では VL 総厚と接線サイズが
    # 結合していて独立に振れず、A〜D の 4 格子が系統列にならなかった (2026-09-19)。
    if abs(a.scale - 1.0) > 1e-9:
        for k in ("n_theta_per_45", "n_gap_half", "n_depth", "n_up", "n_ogrid",
                  "n_disc", "n_left", "n_right", "n_top"):
            H[k] = max(2, int(round(H.get(k, 10) * a.scale)))
        for k in ("gap_first_um", "up_first_um", "mouth_first_um", "floor_first_um",
                  "ogrid_first_mm", "disc_first_mm"):
            if k in H:
                H[k] = H[k] / a.scale
        print("scale %.3f -> %s" % (a.scale, {k: H[k] for k in sorted(H) if not k.startswith("_")}))
    n_az4 = int(H.get("n_theta_per_45", 30))          # 45° あたりの周方向セル数
    n_gap_half = int(H.get("n_gap_half", 12))         # すきま半分の半径方向セル数
    gap_first = H.get("gap_first_um", 10.0) * 1e-6
    n_dep = int(H.get("n_depth", 100))
    # 深さ方向は両端を細かくするが、**開口側と床側で別々に指定**する。
    # 床は死水域で勾配が緩く q'' も小さい一方、そこを詰めすぎると双対 CV が小さくなり
    # float32 のメッシュ量で CV 閉性が悪化する (2026-09-19 実測: 床の 246 CV が 2e-5)。
    mouth_first = H.get("mouth_first_um", H.get("depth_first_um", 20.0)) * 1e-6
    floor_first = H.get("floor_first_um", H.get("depth_first_um", 80.0)) * 1e-6
    n_up = int(H.get("n_up", 55))
    up_first = H.get("up_first_um", 10.0) * 1e-6
    Rs = H.get("Rs_mm", 45.0) * 1e-3                  # O グリッド外側の正方形半幅
    n_og = int(H.get("n_ogrid", 20))
    og_first = H.get("ogrid_first_mm", 0.3) * 1e-3    # 開口外周のリップ側 第一セル
    n_disc = int(H.get("n_disc", 14))
    disc_first = H.get("disc_first_mm", 0.3) * 1e-3
    n_left = int(H.get("n_left", 22))
    n_right = int(H.get("n_right", 38))
    n_top = int(H.get("n_top", 14))

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("cavity_hex")
    geo = gmsh.model.geo

    def P(x, y, z=0.0):
        return geo.addPoint(x, y, z)

    phis = [0.0, 0.25 * math.pi, 0.5 * math.pi, 0.75 * math.pi, math.pi]

    def r_inner(phi):
        """外筒軸からのレイ (cos φ, sin φ) が内円 (中心 (off,0), 半径 Ri) と交わる半径。"""
        c, s = math.cos(phi), math.sin(phi)
        return off * c + math.sqrt(max(Ri ** 2 - (off * s) ** 2, 0.0))

    # ---- 点 ----
    cen_o = P(0, 0)                      # 外筒軸
    cen_i = P(off, 0)                    # 内円柱軸
    rin = [r_inner(p) for p in phis]
    rmid = [0.5 * (rin[k] + Ro) for k in range(5)]
    I = [P(rin[k] * math.cos(phis[k]), rin[k] * math.sin(phis[k])) for k in range(5)]
    Mid = [P(rmid[k] * math.cos(phis[k]), rmid[k] * math.sin(phis[k])) for k in range(5)]
    O = [P(Ro * math.cos(phis[k]), Ro * math.sin(phis[k])) for k in range(5)]
    Sq = [P(Rs, 0), P(Rs, Rs), P(0, Rs), P(-Rs, Rs), P(-Rs, 0)]

    # 中間円 (すきま中央) は「円」ではないので、中点は直線でつなぐ代わりに
    # 内円と外円の中点を通る円弧近似にする: 中心 (off/2, 0)・半径 (Ri+Ro)/2 の円
    cen_m = P(0.5 * off, 0)

    # ---- 曲線 ----
    def arc(p1, c, p2):
        return geo.addCircleArc(p1, c, p2)

    def ln(p1, p2):
        return geo.addLine(p1, p2)

    aI = [arc(I[k], cen_i, I[k + 1]) for k in range(4)]          # 内円 (壁 cyl_side)
    aM = [arc(Mid[k], cen_m, Mid[k + 1]) for k in range(4)]      # すきま中央
    aO = [arc(O[k], cen_o, O[k + 1]) for k in range(4)]          # 外円 (壁 cav_outer)
    rIM = [ln(I[k], Mid[k]) for k in range(5)]                   # 内壁→中央
    rMO = [ln(Mid[k], O[k]) for k in range(5)]                   # 中央→外壁
    rOS = [ln(O[k], Sq[k]) for k in range(5)]                    # 外円→角
    eS = [ln(Sq[k], Sq[k + 1]) for k in range(4)]                # 角の辺

    # ---- 面: すきま (8) ----
    s_gap_in, s_gap_out = [], []
    for k in range(4):
        cl = geo.addCurveLoop([rIM[k], aM[k], -rIM[k + 1], -aI[k]])
        s_gap_in.append(geo.addPlaneSurface([cl]))
        cl = geo.addCurveLoop([rMO[k], aO[k], -rMO[k + 1], -aM[k]])
        s_gap_out.append(geo.addPlaneSurface([cl]))
    s_gap = s_gap_in + s_gap_out

    # ---- 面: O グリッド (4) ----
    s_og = []
    for k in range(4):
        cl = geo.addCurveLoop([rOS[k], eS[k], -rOS[k + 1], -aO[k]])
        s_og.append(geo.addPlaneSurface([cl]))

    # ---- 面: 内円柱上面のバタフライ (4) ----
    aa = 0.55 * Ri
    bb = 0.55 * Ri
    C0, C1, C2, C3 = P(off + aa, 0), P(off + aa, bb), P(off - aa, bb), P(off - aa, 0)
    cCore = [ln(C3, C0), ln(C0, C1), ln(C1, C2), ln(C2, C3)]
    lR0, lR1 = ln(C0, I[0]), ln(C1, I[1])
    lL1, lL0 = ln(C2, I[3]), ln(C3, I[4])
    s_disc = []
    s_disc.append(geo.addPlaneSurface([geo.addCurveLoop(cCore)]))                    # コア
    s_disc.append(geo.addPlaneSurface([geo.addCurveLoop([lR0, aI[0], -lR1, -cCore[1]])]))   # 右
    s_disc.append(geo.addPlaneSurface([geo.addCurveLoop([lR1, aI[1], aI[2], -lL1, -cCore[2]])]))  # 上
    s_disc.append(geo.addPlaneSurface([geo.addCurveLoop([lL1, aI[3], -lL0, -cCore[3]])]))   # 左

    # ---- 面: 遠方 3 ----
    D = [P(x0, 0), P(x0, Rs), P(x0, ym), P(x1, ym), P(x1, Rs), P(x1, 0)]
    lLb = ln(D[0], Sq[4])                 # x_in..-Rs (y=0)
    lLl = ln(D[0], D[1])                  # x=x_in, 0..Rs
    lLt = ln(D[1], Sq[3])                 # y=Rs, x_in..-Rs
    lRb = ln(Sq[0], D[5])                 # Rs..x_out (y=0)
    lRr = ln(D[5], D[4])                  # x=x_out, 0..Rs
    lRt = ln(Sq[1], D[4])                 # y=Rs, Rs..x_out
    lTl = ln(D[1], D[2])                  # x=x_in, Rs..ym
    lTt = ln(D[2], D[3])                  # y=ym
    lTr = ln(D[4], D[3])                  # x=x_out, Rs..ym
    s_far = []
    s_far.append(geo.addPlaneSurface([geo.addCurveLoop([lLb, -eS[3], -lLt, -lLl])]))         # 左
    s_far.append(geo.addPlaneSurface([geo.addCurveLoop([lRb, lRr, -lRt, -eS[0]])]))          # 右
    s_far.append(geo.addPlaneSurface([geo.addCurveLoop([lLt, -eS[2], -eS[1], lRt, lTr, -lTt, -lTl])]))  # 上

    geo.synchronize()

    # ---- 分割数 ----
    def tf(curve, n, kind=None, coef=1.0):
        gmsh.model.geo.mesh.setTransfiniteCurve(abs(curve), n + 1,
                                                kind or "Progression", coef)

    p_gap = progression_ratio(0.5 * (Ro - Ri), n_gap_half, gap_first)
    p_og = progression_ratio(Rs - Ro, n_og, og_first)
    p_disc = progression_ratio(Ri - aa, n_disc, disc_first)
    for k in range(4):
        for c in (aI[k], aM[k], aO[k]):
            tf(c, n_az4)
    for k in range(5):
        tf(rIM[k], n_gap_half, coef=p_gap)      # 内壁から細かく
        tf(rMO[k], n_gap_half, coef=1.0 / p_gap)  # 外壁側を細かく (終点側)
        tf(rOS[k], n_og, coef=p_og)
    for k in range(4):
        tf(eS[k], n_az4)
    # バタフライ
    tf(cCore[0], 2 * n_az4); tf(cCore[2], 2 * n_az4)
    tf(cCore[1], n_az4); tf(cCore[3], n_az4)
    for c in (lR0, lR1, lL1, lL0):
        tf(c, n_disc, coef=1.0 / p_disc)
    # 遠方
    tf(lLb, n_left); tf(lLt, n_left); tf(lRb, n_right); tf(lRt, n_right)
    tf(lLl, n_az4); tf(lRr, n_az4)
    tf(lTl, n_top); tf(lTr, n_top)
    tf(lTt, n_left + 2 * n_az4 + n_right)

    for s in s_gap + s_og + s_disc + s_far:
        gmsh.model.geo.mesh.setTransfiniteSurface(s)
        gmsh.model.geo.mesh.setRecombine(2, s)
    # コーナー指定が要る面 (5 辺以上)
    gmsh.model.geo.mesh.setTransfiniteSurface(s_disc[2], "Left", [C1, I[1], I[3], C2])
    gmsh.model.geo.mesh.setTransfiniteSurface(s_far[2], "Left", [D[1], D[2], D[3], D[4]])
    geo.synchronize()

    # ---- 押し出し ----
    up_sizes = sizes_one_sided(zt, n_up, up_first)
    nl, cl_ = layers_from_sizes(up_sizes)
    up = geo.extrude([(2, s) for s in s_gap + s_og + s_disc + s_far],
                     0, 0, zt, nl, cl_, True)
    dn_sizes = sizes_two_sided(dep, n_dep, mouth_first, floor_first)
    nl2, cl2 = layers_from_sizes(dn_sizes)
    dn = geo.extrude([(2, s) for s in s_gap], 0, 0, -dep, nl2, cl2, True)
    geo.synchronize()

    if a.show_sizes:
        print("すきま横断: 片側 %d セル 第一 %.3g m (比 %.4f), 計 %d セル"
              % (n_gap_half, gap_first, p_gap, 2 * n_gap_half))
        print("深さ: %d セル (開口側 %.3g m / 床側 %.3g m), 上向き: %d セル 第一 %.3g m"
              % (n_dep, mouth_first, floor_first, n_up, up_first))
        print("周方向: %d セル/半周" % (4 * n_az4))

    vols = [t for (d, t) in up + dn if d == 3]
    gmsh.model.addPhysicalGroup(3, vols, MAN["fluid_id"], "fluid")

    # ---- 境界面の分類 (位置で) ----
    bnd = gmsh.model.getBoundary([(3, v) for v in vols], combined=True, oriented=False)
    groups = {}
    tol = 1e-6
    for (d, t) in bnd:
        if d != 2:
            continue
        xa, ya, za, xb, yb, zb = gmsh.model.getBoundingBox(2, t)
        cx, cy, cz = 0.5 * (xa + xb), 0.5 * (ya + yb), 0.5 * (za + zb)
        dz, dy_, dx_ = zb - za, yb - ya, xb - xa
        if dx_ < tol and abs(cx - x0) < tol:
            nm = "inlet"
        elif dx_ < tol and abs(cx - x1) < tol:
            nm = "outlet"
        elif dz < tol and abs(cz - zt) < tol:
            nm = "top"
        elif dy_ < tol and abs(cy - ym) < tol:
            nm = "side"
        elif dy_ < tol and abs(cy) < tol:
            nm = "sym"
        elif dz < tol and abs(cz + dep) < tol:
            nm = "cav_floor"
        elif dz < tol and abs(cz) < tol:
            if math.hypot(cx - off, cy) < Ri:
                nm = "cyl_top"
            else:  # 開口まわり (O グリッド域) は plate_in、その外は plate
                rr = max(abs(xa), abs(xb), abs(ya), abs(yb))
                nm = "plate_in" if rr <= Rs + tol else "plate"
        elif za < -tol:
            # すきまの側壁: 外筒は bbox が ±Ro、内円柱は ±Ri (軸 off)
            rr = max(abs(xa), abs(xb), abs(ya), abs(yb))
            nm = "cav_outer" if abs(rr - Ro) < abs(rr - (abs(off) + Ri)) else "cyl_side"
        else:
            nm = "UNKNOWN(%.4f,%.4f,%.4f)" % (cx, cy, cz)
        groups.setdefault(nm, []).append(t)
    bad = [k for k in groups if k.startswith("UNKNOWN")]
    if bad:
        raise SystemExit("未分類の境界面: %s" % bad)
    for nm, tags in groups.items():
        if nm in PID:
            gmsh.model.addPhysicalGroup(2, tags, PID[nm], nm)
        else:
            print("WARNING: manifest に physID の無いグループ %s (面 %d)" % (nm, len(tags)))
    print("境界グループ:", {k: len(v) for k, v in sorted(groups.items())})
    miss = set(PID) - set(groups) - {"runup"}
    if miss:
        raise SystemExit("境界グループが足りない: %s" % sorted(miss))

    if a.no_mesh:
        gmsh.write(a.out.replace(".msh", ".geo_unrolled"))
        gmsh.finalize()
        return

    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.removeDuplicateNodes()
    ne = gmsh.model.mesh.getElements(3)
    NPE = {4: 4, 5: 8, 6: 6, 7: 5}          # tet/hex/prism/pyramid の節点数
    cnt = {t: len(tags) for t, tags in zip(ne[0], ne[1])}
    nn = len(gmsh.model.mesh.getNodes()[0])
    print("節点 %d,  要素: %s  (hex 率 %.1f %%)"
          % (nn, {"hex" if t == 5 else ("tet" if t == 4 else ("prism" if t == 6 else str(t))): c
                  for t, c in cnt.items()},
             100.0 * cnt.get(5, 0) / max(sum(cnt.values()), 1)))
    gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
    gmsh.write(a.out)
    print("wrote", a.out)
    gmsh.finalize()


if __name__ == "__main__":
    main()
