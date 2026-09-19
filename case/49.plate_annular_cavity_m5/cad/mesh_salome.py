# mesh_salome.py -- STEP -> NETGEN tet + SMESH ViscousLayers (prism BL) -> MED  [case/49]
#
# 実行 (ヘッドレス):
#   export LD_LIBRARY_PATH=/home/sano/opt/salome/syslibs/usr/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH
#   /home/sano/opt/salome/SALOME-9.14.0-native-UB24.04-SRC/salome --keep-paths -t cad/mesh_salome.py \
#       args:<step>,<out.med>
#   (無引数なら manifest から決める)
#
# 設定はすべて ../manifest.json (setup.py --resolve)。このファイルに寸法・サイズの既定値は持たない。
#
# 面分類 (plan §4.1, codex 2026-09-19 plan-1 M4):
#   **重心半径は使わない** (半円筒面の重心半径は 2R/pi で円筒上に無い)。
#   円筒面は KindOfShape の軸・半径、平面は法線と支持位置、z=0 面の内外は頂点の最大半径で判定する。
#   分類後に **グループ別面積を manifest の解析値と照合**する (誤タグは総面積では見えない)。
import json
import math
import os
import sys

import salome

salome.salome_init()
import GEOM  # noqa: E402,F401
import SMESH  # noqa: E402
from salome.geom import geomBuilder  # noqa: E402
from salome.smesh import smeshBuilder  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
MAN = json.load(open(os.path.join(HERE, "..", "manifest.json")))
G, M = MAN["geometry"], MAN["mesh"]

argv = sys.argv[1:]
STEP = argv[0] if len(argv) > 0 and argv[0] else os.path.join(
    HERE, "cavity_plug_half.step" if G["plug_cavity"] else "cavity_fluid_half.step")
OUT = argv[1] if len(argv) > 1 and argv[1] else os.path.join(
    HERE, "cavity_plug.med" if G["plug_cavity"] else "cavity_salome.med")

geompy = geomBuilder.New()
smesh = smeshBuilder.New()

shape = geompy.ImportSTEP(STEP, False, True)
faces = geompy.ExtractShapes(shape, geompy.ShapeType["FACE"], True)

# --- 単位検出: STEP は mm、ImportSTEP が m に直すはず。manifest (m) の x 全長で確かめる ---
bb = geompy.BoundingBox(shape)                       # (xmin,xmax,ymin,ymax,zmin,zmax)
Lx_m = G["x_out"] - G["x_in"]
scale = (bb[1] - bb[0]) / Lx_m                       # 1.0 なら GEOM 単位 = m
print("bbox=", [round(b, 6) for b in bb], " unit scale (GEOM/m) =", round(scale, 6), flush=True)
S = scale                                            # m -> GEOM 単位
TOL = 1.0e-6 * max(1.0, S)


def _vertices(f):
    out = []
    for vtx in geompy.ExtractShapes(f, geompy.ShapeType["VERTEX"], True):
        out.append(geompy.PointCoordinates(vtx))
    return out


def classify(f):
    """face -> グループ名 (build_geom.py の classify と同じ規則)"""
    kind = geompy.KindOfShape(f)
    kname = str(kind[0])
    cdg = geompy.PointCoordinates(geompy.MakeCDG(f))
    Ro, Ri, off = G["Ro"] * S, G["Ri"] * S, G["x_off"] * S
    if "CYLINDER" in kname:
        # ['CYLINDER2D', xc,yc,zc, dx,dy,dz, R, H]
        cx, R = float(kind[1]), float(kind[7])
        if abs(R - Ro) < TOL and abs(cx) < TOL:
            return "cav_outer"
        if abs(R - Ri) < TOL and abs(cx - off) < TOL:
            return "cyl_side"
        return "UNKNOWN_cyl_R%.6f_cx%.6f" % (R, cx)
    # 平面: KindOfShape は PLANAR/POLYGON/RECTANGLE/DISK 等を返し名前が一定でないので、
    # **bbox の潰れている軸**で「どの座標が一定の平面か」を決める (軸平行面しか無い形状)。
    fb = geompy.BoundingBox(f)                      # (xmin,xmax,ymin,ymax,zmin,zmax)
    dx_, dy_, dz_ = fb[1] - fb[0], fb[3] - fb[2], fb[5] - fb[4]
    flat = min(dx_, dy_, dz_)
    if flat > TOL * 10:
        return "UNKNOWN_kind_%s_bbox(%.4g,%.4g,%.4g)" % (kname, dx_, dy_, dz_)
    x, y, z = cdg
    if flat == dx_:
        if abs(x - G["x_in"] * S) < TOL:
            return "inlet"
        if abs(x - G["x_out"] * S) < TOL:
            return "outlet"
        return "UNKNOWN_planeX_%.6f" % x
    if flat == dy_:
        if abs(y - G["y_max"] * S) < TOL:
            return "side"
        if abs(y) < TOL:
            return "sym"
        return "UNKNOWN_planeY_%.6f" % y
    if flat == dz_:
        if abs(z - G["z_top"] * S) < TOL:
            return "top"
        if abs(z + G["depth"] * S) < TOL:
            return "cav_floor"
        if abs(z) < TOL:
            vs = _vertices(f)
            rmax_o = max(math.hypot(p[0], p[1]) for p in vs)
            rmax_i = max(math.hypot(p[0] - off, p[1]) for p in vs)
            xmax = max(p[0] for p in vs)
            if (not G["plug_cavity"]) and rmax_i < Ri + 1e-6 * S:
                return "cyl_top"
            if G["has_runup"] and xmax < G["x_plate"] * S + 1e-6 * S:
                return "runup"
            if G["has_patch"] and rmax_o < G["r_patch"] * S + 1e-6 * S:
                return "plate_in"
            return "plate"
        return "UNKNOWN_planeZ_%.6f" % z
    return "UNKNOWN_plane_%s" % kname


bins, area = {}, {}
for f in faces:
    n = classify(f)
    bins.setdefault(n, []).append(f)
    area[n] = area.get(n, 0.0) + geompy.BasicProperties(f)[1]
bad = [k for k in bins if k.startswith("UNKNOWN")]
if bad:
    raise SystemExit("ERROR: 未分類の面 %s" % bad)
exp = G["group_area_mm2"]                      # [mm^2]
miss = set(exp) - set(bins)
extra = set(bins) - set(exp)
if miss or extra:
    raise SystemExit("ERROR: グループ不一致 missing=%s extra=%s" % (miss, extra))
print("%-10s %5s %14s %14s %10s" % ("group", "faces", "area[mm2]", "expected", "rel.err"), flush=True)
worst = 0.0
for k in sorted(exp):
    a_mm2 = area[k] / (S * S) * 1.0e6          # GEOM 単位 -> m^2 -> mm^2
    r = abs(a_mm2 - exp[k]) / max(exp[k], 1e-30)
    worst = max(worst, r)
    print("%-10s %5d %14.4f %14.4f %10.2e" % (k, len(bins[k]), a_mm2, exp[k], r), flush=True)
print("worst group area rel.err = %.2e" % worst, flush=True)
if worst > 1e-5:
    raise SystemExit("ERROR: グループ別面積が manifest と一致しない (誤タグの疑い)")

ggrp = {}
for name, fl in bins.items():
    g = geompy.CreateGroup(shape, geompy.ShapeType["FACE"])
    geompy.UnionList(g, fl)
    g.SetName(name)
    geompy.addToStudyInFather(shape, g, name)
    ggrp[name] = g

# ---------------------------------------------------------------- メッシュ
mesh = smesh.Mesh(shape, "cavity")
algo = mesh.Tetrahedron(smeshBuilder.NETGEN_1D2D3D)
par = algo.Parameters()
par.SetMaxSize(M["maxh_m"] * S)
par.SetMinSize(M["minh_m"] * S)
par.SetFineness(M["fineness"])
# 体積要素の成長率 (NETGEN 既定 0.3)。小さいほど壁から離れるときの粗大化が緩やかになり、
# **prism 層 -> tet の継ぎ目の段差**が小さくなる (ユーザ要件 2026-09-19)。
try:
    par.SetGrowthRate(M.get("growth_rate", 0.3))
    print("growth rate =", M.get("growth_rate", 0.3), flush=True)
except Exception as _e:                       # noqa: BLE001
    print("WARNING: SetGrowthRate 不可:", _e, flush=True)
par.SetSecondOrder(0)                       # forge は線形要素のみ
par.SetOptimize(1)

SIZE = {"plate": "size_plate_m", "plate_in": "size_plate_in_m",
        "cav_outer": "size_cav_m", "cyl_side": "size_cav_m",
        "cav_floor": "size_floor_m", "cyl_top": "size_cyl_top_m"}
n_loc = 0
for grp, key in SIZE.items():
    for f in bins.get(grp, []):
        geompy.addToStudyInFather(shape, f, "loc_%s_%d" % (grp, n_loc))
        par.SetLocalSizeOnShape(f, M[key] * S)
        n_loc += 1
print("local size on %d faces" % n_loc, flush=True)

# 開口の円エッジ (r=Ro, Ri @ z=0) を細分: せん断層が入る所
if not G["plug_cavity"]:
    n_e = 0
    esz = M["size_edge_opening_m"] * S
    for grp in ("cav_outer", "cyl_side", "cyl_top", "plate_in"):
        for f in bins.get(grp, []):
            for e in geompy.ExtractShapes(f, geompy.ShapeType["EDGE"], True):
                ex, ey, ez = geompy.PointCoordinates(geompy.MakeCDG(e))
                if abs(ez) < TOL:            # z=0 の縁だけ
                    geompy.addToStudyInFather(shape, e, "opening_edge_%d" % n_e)
                    par.SetLocalSizeOnShape(e, esz)
                    n_e += 1
    print("opening edge size %.4g on %d edges" % (esz, n_e), flush=True)

# Viscous layers: 壁グループにだけ張る
NO_LAYER = ["inlet", "outlet", "top", "side", "sym", "runup"]
ignore = []
for n in NO_LAYER:
    for f in bins.get(n, []):
        ignore.append(geompy.GetSubShapeID(shape, f))
algo.ViscousLayers(M["vl_total_m"] * S, M["vl_nlayers"], M["vl_stretch"], ignore, True)
print("VL: total %.6g m x %d layers (stretch %.3f, first %.3g m)"
      % (M["vl_total_m"], M["vl_nlayers"], M["vl_stretch"], M["vl_first_m"]), flush=True)

ok = mesh.Compute()
print("Compute:", ok, flush=True)
if not ok:
    raise SystemExit("ERROR: mesh.Compute() failed")

for name, g in ggrp.items():
    mesh.GroupOnGeom(g, name, SMESH.FACE)

print("nodes: %d  tet %d  prism %d  pyram %d  hex %d"
      % (mesh.NbNodes(), mesh.NbTetras(), mesh.NbPrisms(), mesh.NbPyramids(), mesh.NbHexas()),
      flush=True)
if mesh.NbNodes() > M["node_budget"]:
    print("WARNING: 節点数 %d が予算 %d を超過" % (mesh.NbNodes(), M["node_budget"]), flush=True)
mesh.ExportMED(OUT)
print("wrote", OUT, flush=True)

used = {"step": STEP, "out": OUT, "unit_scale": scale, "mesh": M,
        "nodes": mesh.NbNodes(), "tet": mesh.NbTetras(), "prism": mesh.NbPrisms(),
        "pyram": mesh.NbPyramids(), "group_faces": {k: len(v) for k, v in bins.items()}}
with open(os.path.join(os.path.dirname(os.path.abspath(OUT)), "mesh_settings_used.json"), "w") as f:
    json.dump(used, f, indent=2)
print("wrote mesh_settings_used.json", flush=True)
