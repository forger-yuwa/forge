# build_geom.py -- M5 平板 + 環状深キャビティ の流体ドメイン (半割モデル) を生成する [単位: mm]
#
# 形状: 断熱平板 (z=0) に 半径 Ro / 深さ depth の円柱キャビティがあり、その中に
#       半径 Ri / 同深さの円柱が上面面一 (z=0) で埋まっている。内側円柱は x 方向に
#       x_off だけ偏心できる (下流ずらし) -> すきまが周方向に不均一になる。
#       流体側 = 「底閉じ環状スリット」+ 平板上の外部流領域。
#       **x 方向の偏心は y=0 対称を保つ** ので半割 (y>=0) のまま扱える。
#
# 実行:
#   QT_QPA_PLATFORM=offscreen /home/sano/opt/squashfs-root/usr/bin/freecadcmd build_geom.py
#   cwd の geom_config.json で P を上書き。使用値と派生量 (すきま min/max, グループ面積) を
#   geom_used.json に保存し、mesh_salome.py がそれを読んでサイズを決める (単一ソース)。
#
# 構築 (ブーリアンを壊さないコツ):
#   環状リングを箱内へ protrude だけ突出させて箱と fuse (同一平面 z=0 の coincident-face
#   fuse を避ける) -> y<0 を cut。**removeSplitter は呼ばない** (箱を x_plate で 2 分割した
#   seam と 細分パッチ円柱の seam を z=0 面に残し、runup / plate_in / plate を別グループに
#   できるようにするため)。
#
# 検証 (満たすまで作り直す; codex 2026-09-19 plan レビュー M4 の要求):
#   - 単一閉ソリッド / bbox / 総表面積が解析値と一致 (内部面の残留を検出)
#   - **グループ別面積が解析値と一致** (誤タグを検出。総面積一致では検出できない)
#   - 未分類 0 / 重複 0 / 全 face が必ず 1 グループ
#   - 面分類は重心半径ではなく **曲面種別 (Part.Cylinder/Plane) の軸・半径**で行う
#     (半円筒面の重心半径は 2R/pi で円筒上に無いため重心では判定できない)

import json
import math
import os
import traceback
from pathlib import Path

import FreeCAD as App  # noqa: F401  (freecadcmd が要求)
import Part
from FreeCAD import Vector

# ===== パラメータ: manifest.json (setup.py --resolve) から読む =====
# 形状の正本は ../case.json -> ../manifest.json。ここには既定値を持たない (転記事故の防止)。
MAN = json.loads((Path(__file__).resolve().parents[1] / "manifest.json").read_text())
_G = MAN["geometry"]
_LEN = ("Ro", "Ri", "x_off", "depth", "x_in", "x_plate", "x_out", "y_max", "z_top",
        "r_patch", "protrude")
P = {k: (_G[k] * 1e3 if k in _LEN else _G[k]) for k in _LEN + ("fuzzy", "plug_cavity")}
EXP_AREA = _G["group_area_mm2"]          # グループ別の期待面積 [mm^2] (正本は setup.py)
OUT = "cavity_plug_half.step" if P["plug_cavity"] else "cavity_fluid_half.step"

if os.path.exists("geom_config.json"):   # 実験用の一時上書き (常用しない)
    with open("geom_config.json") as f:
        _u = json.load(f)
    for k, v in _u.items():
        if k == "out":
            OUT = v
        elif not k.startswith("_"):
            P[k] = v
    print("WARNING: geom_config.json で manifest を上書き:", _u, flush=True)


def derived(P):
    """派生量: すきま min/max と、manifest のグループ別解析面積 [mm^2]"""
    Ro, Ri, off = P['Ro'], P['Ri'], P['x_off']
    a = dict(EXP_AREA)
    return dict(gap_min=Ro - Ri - abs(off), gap_max=Ro - Ri + abs(off), gap_nom=Ro - Ri,
                group_area=a, total_area=sum(a.values()),
                has_runup=P['x_plate'] > P['x_in'] + 1e-9,
                has_patch=P.get('r_patch', 0.0) > Ro,
                plug=bool(P.get('plug_cavity', False)))


def assemble(P):
    """流体 = 箱 (+助走分割) + 環状リング(突出) (+細分パッチ円柱) を fuse し、y<0 を cut"""
    y0, ym, zt = -P['y_max'], P['y_max'], P['z_top']
    parts = []
    if P['x_plate'] > P['x_in'] + 1e-9:
        base = Part.makeBox(P['x_plate'] - P['x_in'], ym - y0, zt, Vector(P['x_in'], y0, 0.0))
        parts.append(Part.makeBox(P['x_out'] - P['x_plate'], ym - y0, zt,
                                  Vector(P['x_plate'], y0, 0.0)))
    else:
        base = Part.makeBox(P['x_out'] - P['x_in'], ym - y0, zt, Vector(P['x_in'], y0, 0.0))
    # 環状リング: z = -depth .. +protrude (上端は箱の内部へ)。内側円柱は x_off 偏心。
    h = P['depth'] + P['protrude']
    outer = Part.makeCylinder(P['Ro'], h, Vector(0, 0, -P['depth']), Vector(0, 0, 1))
    inner = Part.makeCylinder(P['Ri'], h, Vector(P['x_off'], 0, -P['depth']), Vector(0, 0, 1))
    if not P.get('plug_cavity', False):
        parts.append(outer.cut(inner))
    # 細分パッチ: 箱の内部に収まる円柱を fuse (体積・表面積は不変、z=0 面に seam が残る)
    if P.get('r_patch', 0.0) > P['Ro']:
        parts.append(Part.makeCylinder(P['r_patch'], zt, Vector(0, 0, 0.0), Vector(0, 0, 1)))

    f = base.fuse(parts, P['fuzzy'])
    big = 4.0 * max(abs(P['x_in']), P['x_out'], P['y_max'], zt)
    return f.cut(Part.makeBox(big, big, big, Vector(-big / 2, -big, -big / 2)))  # y<0 を削る


def classify(fc, P, D, tol=1.0e-4):
    """face -> グループ名。円筒は Surface の軸・半径、平面は支持平面と範囲で判定する。"""
    s = fc.Surface
    kind = type(s).__name__
    if kind == "Cylinder":
        R, cx = s.Radius, s.Center.x
        if abs(R - P['Ro']) < tol and abs(cx) < tol:
            return "cav_outer"
        if abs(R - P['Ri']) < tol and abs(cx - P['x_off']) < tol:
            return "cyl_side"
        return "UNKNOWN_cyl_R%.4f_cx%.4f" % (R, cx)
    if kind != "Plane":
        return "UNKNOWN_surf_%s" % kind
    n = s.Axis
    c = fc.CenterOfMass
    if abs(abs(n.x) - 1.0) < 1e-6:                     # x 一定面
        if abs(c.x - P['x_in']) < tol:
            return "inlet"
        if abs(c.x - P['x_out']) < tol:
            return "outlet"
        return "UNKNOWN_planeX_%.4f" % c.x
    if abs(abs(n.y) - 1.0) < 1e-6:                     # y 一定面
        if abs(c.y - P['y_max']) < tol:
            return "side"
        if abs(c.y) < tol:
            return "sym"
        return "UNKNOWN_planeY_%.4f" % c.y
    if abs(abs(n.z) - 1.0) < 1e-6:                     # z 一定面
        if abs(c.z - P['z_top']) < tol:
            return "top"
        if abs(c.z + P['depth']) < tol:
            return "cav_floor"
        if abs(c.z) < tol:
            # 重心半径は使えない (面が曲線境界を持つと重心は境界円の内側に寄る)。
            # 頂点の最大半径で「その円内に収まる面か」を判定する。
            vs = [(v.Point.x, v.Point.y) for v in fc.Vertexes]
            rmax_o = max(math.hypot(vx, vy) for vx, vy in vs)          # 外筒軸まわり
            rmax_i = max(math.hypot(vx - P['x_off'], vy) for vx, vy in vs)  # 内円柱軸まわり
            xmax = max(vx for vx, _ in vs)
            if rmax_i < P['Ri'] + 1e-6:
                return "cyl_top"
            if D['has_runup'] and xmax < P['x_plate'] + 1e-6:
                return "runup"
            if D['has_patch'] and rmax_o < P['r_patch'] + 1e-6:
                return "plate_in"
            return "plate"
        return "UNKNOWN_planeZ_%.4f" % c.z
    return "UNKNOWN_plane_n(%.3f,%.3f,%.3f)" % (n.x, n.y, n.z)


def main():
    if abs(P['x_off']) >= P['Ro'] - P['Ri'] - 1e-9:
        raise SystemExit("ERROR: |x_off|=%g は すきま Ro-Ri=%g 以上 (内側円柱が外筒に接触)"
                         % (abs(P['x_off']), P['Ro'] - P['Ri']))
    D = derived(P)
    print("gap: nom %.3f  min %.3f  max %.3f mm" % (D['gap_nom'], D['gap_min'], D['gap_max']), flush=True)
    fluid = None
    for attempt in range(6):
        fluid = assemble(P)
        b = fluid.BoundBox
        closed = bool(fluid.Shells) and fluid.Shells[0].isClosed()
        aerr = abs(fluid.Area - D['total_area']) / D['total_area']
        ok = (fluid.isValid() and len(fluid.Solids) == 1 and closed
              and abs(b.XMin - P['x_in']) < 1e-3 and abs(b.XMax - P['x_out']) < 1e-3
              and abs(b.YMin) < 1e-3 and abs(b.YMax - P['y_max']) < 1e-3
              and abs(b.ZMin + (0.0 if D['plug'] else P['depth'])) < 1e-3 and abs(b.ZMax - P['z_top']) < 1e-3
              and aerr < 1e-6)
        print("attempt %d: valid=%s solids=%d closed=%s faces=%d area=%.4f (exp %.4f, err %.2e) "
              "bbox x[%.2f,%.2f] y[%.3f,%.2f] z[%.2f,%.2f]"
              % (attempt, fluid.isValid(), len(fluid.Solids), closed, len(fluid.Faces),
                 fluid.Area, D['total_area'], aerr,
                 b.XMin, b.XMax, b.YMin, b.YMax, b.ZMin, b.ZMax))
        if ok:
            break
    else:
        raise SystemExit("ERROR: 期待形状に到達しなかった")

    # --- 面分類と グループ別面積の照合 ---
    got = {}
    for fc in fluid.Faces:
        n = classify(fc, P, D)
        got.setdefault(n, [0, 0.0])
        got[n][0] += 1
        got[n][1] += fc.Area
    bad = [k for k in got if k.startswith("UNKNOWN")]
    if bad:
        raise SystemExit("ERROR: 未分類の面 %s" % {k: got[k] for k in bad})
    exp = D['group_area']
    miss = set(exp) - set(got)
    extra = set(got) - set(exp)
    if miss or extra:
        raise SystemExit("ERROR: グループ不一致 missing=%s extra=%s" % (miss, extra))
    print("%-10s %5s %14s %14s %10s" % ("group", "faces", "area", "expected", "rel.err"), flush=True)
    worst = 0.0
    for k in sorted(exp):
        e = exp[k]
        r = abs(got[k][1] - e) / max(e, 1e-12)
        worst = max(worst, r)
        print("%-10s %5d %14.4f %14.4f %10.2e" % (k, got[k][0], got[k][1], e, r), flush=True)
    print("worst group area rel.err = %.2e" % worst, flush=True)
    if worst > 1e-6:
        raise SystemExit("ERROR: グループ別面積が解析値と一致しない (誤タグの疑い)")

    fluid.exportStep(OUT)
    used = dict(P)
    used.update(out=OUT, gap_min=D['gap_min'], gap_max=D['gap_max'], gap_nom=D['gap_nom'],
                group_area=D['group_area'], total_area_expected=D['total_area'],
                total_area=fluid.Area, n_faces=len(fluid.Faces),
                group_faces={k: v[0] for k, v in got.items()})
    with open("geom_used.json", "w") as f:
        json.dump(used, f, indent=2)
    print("wrote", OUT, "and geom_used.json", flush=True)


try:
    main()
except SystemExit as e:
    print("FAILED:", e, flush=True)
    raise
except Exception:
    traceback.print_exc()
    raise
