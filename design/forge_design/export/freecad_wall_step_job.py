"""`freecadcmd` の中で動く STEP の書き出し・読み直しの作業スクリプト (forge_design からは import しない)。

`forge_design.export.wall_step` が作業 JSON を書き、`WALL_STEP_JOB=<作業 JSON> freecadcmd <このファイル>` で呼ぶ
(スクリプトはファイルで渡す。`-c` は使わない)。結果は作業 JSON の `result` に書いたパスへ JSON で返す。

mode = write: 次数・制御点 (x, y) [mm]・異なるノット [mm]・重複度から `Part.BSplineCurve.buildFromPolesMultsKnots`
              (非周期・重み 1) で曲線を作り、辺 1 本だけを `exportStep` で書く。
mode = read:  `Part.Shape().read()` で読み直し、辺の数・曲線の型・次数・ノット・重複度・制御点・重み・辺の助変数の範囲と
              頂点の座標 (実際の端点) と、指定した助変数 u での位置・1 階・2 階微分 (`getD2`) を返す。`revolve` が真なら x 軸まわりに 360° 回して回転面 (内面) を作り、
              型・妥当性 (`isValid`)・面積を返す。
mode = surface: 2D の STEP (`step2d`、辺 1 本) を読み、辺を x 軸まわりに 360° 回した面 (内面) を殻にして `exportStep` で
              3D の STEP (`step3d`) に書き、`MeshPart.meshFromShape(LinearDeflection, AngularDeflection)` で STL (`stl`、mm) を書く。
              書いた 3D の STEP を `Part.Shape().read()` で読み直し、面の数・妥当性・面積・bbox を返す。面積は `Shape.Area` と、
              母線のノット区間ごとに切り出した面の和 (`readback_area_knot_split_mm2`、検査にはこちらを使う; `_area_knot_split`)
              (plans/accepted/tooling-nozzle-wt-deliverables.md §4.2 の cad/、case/45 の `cad_deliverables.py` の surface を一般化)。
"""
import json
import os
import traceback

import FreeCAD  # noqa: F401  (freecadcmd の中でだけ import できる)
import Part
from FreeCAD import Vector


def _write(job):
    c = Part.BSplineCurve()
    poles = [Vector(float(x), float(y), 0.0) for x, y in job["poles"]]
    c.buildFromPolesMultsKnots(poles, [int(m) for m in job["mults"]], [float(u) for u in job["knots"]],
                               False, int(job["degree"]))
    e = c.toShape()
    e.exportStep(job["step"])
    return {"degree": c.Degree, "n_poles": c.NbPoles, "is_rational": bool(c.isRational()),
            "first": c.FirstParameter, "last": c.LastParameter}


def _read(job):
    s = Part.Shape()
    s.read(job["step"])
    out = {"n_edges": len(s.Edges), "n_faces": len(s.Faces), "n_solids": len(s.Solids)}
    e = s.Edges[0]
    cc = e.Curve
    out["curve_type"] = type(cc).__name__
    if out["curve_type"] != "BSplineCurve":
        return out
    out.update({"degree": cc.Degree, "n_poles": cc.NbPoles, "is_rational": bool(cc.isRational()),
                "is_periodic": bool(cc.isPeriodic()), "knots": [float(u) for u in cc.getKnots()],
                "mults": [int(m) for m in cc.getMultiplicities()],
                "poles": [[p.x, p.y, p.z] for p in cc.getPoles()], "weights": [float(w) for w in cc.getWeights()],
                "first": cc.FirstParameter, "last": cc.LastParameter,
                "edge_first": e.FirstParameter, "edge_last": e.LastParameter,
                "vertices": [[v.Point.x, v.Point.y, v.Point.z] for v in e.Vertexes]})
    ev = []
    for u in job.get("u", []):
        p0, d1, d2 = cc.getD2(float(u))
        ev.append([[p0.x, p0.y, p0.z], [d1.x, d1.y, d1.z], [d2.x, d2.y, d2.z]])
    out["eval"] = ev
    if job.get("revolve"):
        rv = e.revolve(Vector(0, 0, 0), Vector(1, 0, 0), 360)
        out["revolve"] = {"shape_type": rv.ShapeType, "is_valid": bool(rv.isValid()), "area_mm2": float(rv.Area),
                          "n_faces": len(rv.Faces), "surface_type": (type(rv.Faces[0].Surface).__name__ if rv.Faces else None)}
    return out


def _area_knot_split(shape):
    """読み直した回転面の面積を、母線 (`BasisCurve`) のノット区間ごとに面を切り出して足す。
    `Shape.Area` (OCC の GProp) は長い B-spline を母線にした 1 枚の回転面では求積が粗い (case/45 の生産の壁で
    相対 −4.1e-4、`toNurbs()` してからでも +1.6e-4。ノット区間ごとに分けた和は ∫2πr ds と 1e-14 で一致; 2026-10-13 測定)。
    回転面でない面があれば None (呼び出し側で不合格にする)。"""
    tot = 0.0
    for f in shape.Faces:
        srf = f.Surface
        bc = getattr(srf, "BasisCurve", None)
        if bc is None or not hasattr(bc, "getKnots"):
            return None
        u0, u1, v0, v1 = f.ParameterRange               # 回転面: U = 回転角、V = 母線の助変数
        vs = [v0] + [float(k) for k in bc.getKnots() if v0 < k < v1] + [v1]
        tot += sum(srf.toShape(u0, u1, a, b).Area for a, b in zip(vs[:-1], vs[1:]))
    return float(tot)


def _surface(job):
    import MeshPart  # noqa: PLC0415  (freecadcmd の中でだけ import できる)
    s = Part.Shape()
    s.read(job["step2d"])
    if len(s.Edges) != 1:
        raise ValueError(f"2D の STEP の辺が 1 本でない ({len(s.Edges)} 本): {job['step2d']}")
    e = s.Edges[0]
    face = e.revolve(Vector(0, 0, 0), Vector(1, 0, 0), 360)
    shell = Part.Shell([face]) if face.ShapeType == "Face" else face
    shell.exportStep(job["step3d"])
    lin, ang = float(job["linear_deflection_mm"]), float(job["angular_deflection"])
    m = MeshPart.meshFromShape(Shape=shell, LinearDeflection=lin, AngularDeflection=ang)
    m.write(job["stl"])
    rd = Part.Shape()
    rd.read(job["step3d"])
    bb = rd.BoundBox
    return {"shape_type": shell.ShapeType, "is_valid": bool(shell.isValid()), "area_mm2": float(shell.Area),
            "surface_type": (type(shell.Faces[0].Surface).__name__ if shell.Faces else None),
            "readback_faces": len(rd.Faces), "readback_valid": bool(rd.isValid()), "readback_area_mm2": float(rd.Area),
            "readback_area_knot_split_mm2": _area_knot_split(rd),
            "readback_bbox_mm": [bb.XMin, bb.YMin, bb.ZMin, bb.XMax, bb.YMax, bb.ZMax],
            "stl_triangles": int(m.CountFacets), "stl_points": int(m.CountPoints), "stl_area_mm2": float(m.Area),
            "stl_linear_deflection_mm": lin, "stl_angular_deflection_rad": ang}


def main():
    job = json.load(open(os.environ["WALL_STEP_JOB"]))
    try:
        res = {"ok": True, "freecad_version": list(FreeCAD.Version()[:3])}
        # write・read の振り分けは従来どおり (surface 以外の mode は従来の分岐をそのまま通る)
        res.update(_write(job) if job["mode"] == "write" else (_surface(job) if job["mode"] == "surface" else _read(job)))
    except Exception:  # noqa: BLE001 — 失敗は呼び出し側に文字列で返す
        res = {"ok": False, "error": traceback.format_exc()}
    with open(job["result"], "w") as f:
        json.dump(res, f)


main()
