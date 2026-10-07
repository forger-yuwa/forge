"""case/45 の生産の壁の CAD 向けの成果物と、STEP の独立の読み直し (2026-10-08、ユーザ「cad ある？ いろいろすすめてほしい」)。

  ref    (.venv-opt)  : 保存した壁 (_band_ab/prod_confirm/prep/wall_repr.json) から参照点と座標表を作る
                        → step/ref_points.csv (u = x [mm] と r [mm]、STEP の照合用)、step/wall_coordinates_mm.csv (加工向けの座標表)
  gmsh   (.venv-mesh) : gmsh (OpenCascade、FreeCAD と別の読み込み経路) で wall_physical.step を読み直し、曲線の本数・助変数の範囲・
                        参照点との差・x 軸まわりの回転面の面積を測る → step/gmsh_check.json
  surface (freecadcmd): FreeCAD で曲線を x 軸まわりに 360° 回して内面を作り、3D の STEP (wall_surface.step) と STL (wall_surface.stl)
                        を書き出す → step/surface_check.json

usage (case dir):
  /home/sano/work/forge/design/.venv-opt/bin/python cad_deliverables.py ref
  /home/sano/work/forge/.venv-mesh/bin/python cad_deliverables.py gmsh
  CAD_JOB=surface /home/sano/opt/squashfs-root/usr/bin/freecadcmd cad_deliverables.py
"""
import json
import os
import sys
from pathlib import Path

HERE = Path(os.environ.get("CAD_CASE_DIR", str(Path(__file__).resolve().parent)))
P = HERE / "_band_ab" / "prod_confirm"
STEP_DIR = P / "step"


def ref():
    import numpy as np
    sys.path.insert(0, str(HERE.parents[1] / "design"))
    from forge_design.geometry.wall_axismach import load_wall_file
    W = load_wall_file(P / "prep")
    s_mm = float(W["scale_m"]) * 1000.0
    x0, x1 = (float(v) for v in W["domain"])
    phys = W["physical"]
    # 照合用: 助変数 u = x [mm] の 20001 点
    xr = np.linspace(x0, x1, 20001)
    np.savetxt(STEP_DIR / "ref_points.csv", np.c_[xr * s_mm, phys.r(xr) * s_mm], delimiter=",", header="x_mm,r_mm", comments="")
    # 加工向けの座標表: 1 mm ごと + スロート付近 (|x| ≤ 50 mm) は 0.1 mm ごと。壁角 [deg] と曲率半径 [mm] も付ける
    xm = np.unique(np.r_[np.arange(np.ceil(x0 * s_mm), np.floor(x1 * s_mm) + 1, 1.0), np.arange(-50.0, 50.0001, 0.1), x0 * s_mm, x1 * s_mm])
    xm = xm[(xm >= x0 * s_mm) & (xm <= x1 * s_mm)]
    xt = xm / s_mm
    r = phys.r(xt) * s_mm
    r1 = phys.r(xt, 1)
    r2 = phys.r(xt, 2) / s_mm
    ang = np.degrees(np.arctan(r1))
    kap = r2 / (1.0 + r1 ** 2) ** 1.5
    rc = np.where(np.abs(kap) > 1e-12, 1.0 / np.where(np.abs(kap) > 1e-12, kap, 1.0), np.inf)
    hdr = ("x_mm,r_mm,wall_angle_deg,curvature_radius_mm  "
           "# case/45 生産の物理壁 (problem_d155_ns_prod.yaml、全域 1 本の 5 次 B-spline から評価)。原点 = 設計スロート、x = 流れ方向、"
           "r = 半径。曲率半径は符号付き (正 = 流路側から見て凸の向き、直管では inf)")
    np.savetxt(STEP_DIR / "wall_coordinates_mm.csv", np.c_[xm, r, ang, rc], delimiter=",", header=hdr, comments="", fmt="%.9f")
    out = {"n_ref": len(xr), "n_table": len(xm), "domain_mm": [x0 * s_mm, x1 * s_mm], "throat": W["throat"], "scale_mm_per_rt": s_mm}
    (STEP_DIR / "ref_summary.json").write_text(json.dumps(out, indent=1, ensure_ascii=False, default=float))
    print(json.dumps(out, indent=1, ensure_ascii=False, default=float))


def gmsh_check():
    import gmsh
    import numpy as np
    ref = np.loadtxt(STEP_DIR / "ref_points.csv", delimiter=",", skiprows=1)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("wall")
    ents = gmsh.model.occ.importShapes(str(STEP_DIR / "wall_physical.step"))
    gmsh.model.occ.synchronize()
    curves = gmsh.model.getEntities(1)
    out = {"gmsh": gmsh.__version__, "imported": ents, "n_curves": len(curves), "n_points": len(gmsh.model.getEntities(0))}
    tag = curves[0][1]
    lo, hi = gmsh.model.getParametrizationBounds(1, tag)
    out["param_range"] = [float(lo[0]), float(hi[0])]
    out["curve_type"] = gmsh.model.getType(1, tag)
    bb = gmsh.model.getBoundingBox(1, tag)
    out["bbox_mm"] = [float(v) for v in bb]
    # u = x [mm] の点で評価して参照と比べる (曲線の端は定義域の中に入れる)
    u = np.clip(ref[:, 0], lo[0], hi[0])
    xyz = np.asarray(gmsh.model.getValue(1, tag, u.tolist())).reshape(-1, 3)
    out["max_abs_dx_mm"] = float(np.abs(xyz[:, 0] - ref[:, 0]).max())
    out["max_abs_dr_mm"] = float(np.abs(xyz[:, 1] - ref[:, 1]).max())
    out["max_abs_z_mm"] = float(np.abs(xyz[:, 2]).max())
    out["r_min_mm"] = float(xyz[:, 1].min())
    out["x_at_r_min_mm"] = float(xyz[np.argmin(xyz[:, 1]), 0])
    # x 軸まわりに 360° 回す (gmsh は 2π を 1 回では作れないので π を 2 回)
    s1 = gmsh.model.occ.revolve([(1, tag)], 0, 0, 0, 1, 0, 0, np.pi)
    s2 = gmsh.model.occ.revolve([(1, tag)], 0, 0, 0, 1, 0, 0, -np.pi)
    gmsh.model.occ.synchronize()
    surfs = [e for e in s1 + s2 if e[0] == 2]
    out["revolve_surfaces"] = len(surfs)
    out["revolve_area_mm2"] = float(sum(gmsh.model.occ.getMass(2, t) for _, t in surfs))
    gmsh.finalize()
    (STEP_DIR / "gmsh_check.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
    print(json.dumps(out, indent=1, ensure_ascii=False))


def surface():
    import FreeCAD  # noqa: F401
    import Mesh
    import MeshPart
    import Part
    from FreeCAD import Vector
    s = Part.Shape()
    s.read(str(STEP_DIR / "wall_physical.step"))
    e = s.Edges[0]
    face = e.revolve(Vector(0, 0, 0), Vector(1, 0, 0), 360)
    shell = Part.Shell([face]) if face.ShapeType == "Face" else face
    shell.exportStep(str(STEP_DIR / "wall_surface.step"))
    m = MeshPart.meshFromShape(Shape=shell, LinearDeflection=0.05, AngularDeflection=0.05)
    m.write(str(STEP_DIR / "wall_surface.stl"))
    rd = Part.Shape()
    rd.read(str(STEP_DIR / "wall_surface.step"))
    bb = rd.BoundBox
    out = {"freecad": list(FreeCAD.Version()[:3]), "shape_type": shell.ShapeType, "is_valid": bool(shell.isValid()),
           "area_mm2": float(shell.Area), "readback_faces": len(rd.Faces), "readback_valid": bool(rd.isValid()),
           "readback_area_mm2": float(rd.Area), "bbox_mm": [bb.XMin, bb.YMin, bb.ZMin, bb.XMax, bb.YMax, bb.ZMax],
           "stl_triangles": int(m.CountFacets), "stl_linear_deflection_mm": 0.05}
    (STEP_DIR / "surface_check.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
    print(json.dumps(out, indent=1, ensure_ascii=False))


# freecadcmd の中では __name__ が "__main__" にならないので、CAD_JOB があれば名前によらず動く
if __name__ == "__main__" or os.environ.get("CAD_JOB"):
    job = os.environ.get("CAD_JOB") or (sys.argv[1] if len(sys.argv) > 1 else "")
    {"ref": ref, "gmsh": gmsh_check, "surface": surface}[job]()
