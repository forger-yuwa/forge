"""`freecadcmd` の中で動く STEP の書き出し・読み直しの作業スクリプト (forge_design からは import しない)。

`forge_design.export.wall_step` が作業 JSON を書き、`WALL_STEP_JOB=<作業 JSON> freecadcmd <このファイル>` で呼ぶ
(スクリプトはファイルで渡す。`-c` は使わない)。結果は作業 JSON の `result` に書いたパスへ JSON で返す。

mode = write: 次数・制御点 (x, y) [mm]・異なるノット [mm]・重複度から `Part.BSplineCurve.buildFromPolesMultsKnots`
              (非周期・重み 1) で曲線を作り、辺 1 本だけを `exportStep` で書く。
mode = read:  `Part.Shape().read()` で読み直し、辺の数・曲線の型・次数・ノット・重複度・制御点・重み・辺の助変数の範囲と
              頂点の座標 (実際の端点) と、指定した助変数 u での位置・1 階・2 階微分 (`getD2`) を返す。`revolve` が真なら x 軸まわりに 360° 回して回転面 (内面) を作り、
              型・妥当性 (`isValid`)・面積を返す。
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


def main():
    job = json.load(open(os.environ["WALL_STEP_JOB"]))
    try:
        res = {"ok": True, "freecad_version": list(FreeCAD.Version()[:3])}
        res.update(_write(job) if job["mode"] == "write" else _read(job))
    except Exception:  # noqa: BLE001 — 失敗は呼び出し側に文字列で返す
        res = {"ok": False, "error": traceback.format_exc()}
    with open(job["result"], "w") as f:
        json.dump(res, f)


main()
