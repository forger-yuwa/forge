"""3D の STEP (回転内面) を gmsh (OpenCascade) で読み直し、面の面積の合計と bbox を JSON で出す単体スクリプト
(plans/active/tooling-nozzle-wt-deliverables.md §6 V2 の独立の読み直し。FreeCAD と別の読み込み経路)。

forge_design に依存しない (`.venv-mesh` の python で動かす。design の venv に gmsh は無い)。単位は STEP の単位 (forge の STEP は mm)。

- **面積** `area_mm2`: gmsh が読んだ面の助変数 $(u, v)$ の上で、面要素 $|X_u \times X_v|$ (`getDerivative`) を
  合成 4 点 Gauss で積分する。区間数は $u$・$v$ の向きごとに倍にして変化の大きい向きを細かくし、相対変化が `--tol` 未満で止める
  (回転面の向きを仮定しない)。
  OCC の `getMass(2, tag)` (GProp) は `area_occ_mass_mm2` として参考に出すだけ: 長い B-spline を母線にした 1 枚の回転面では
  求積が粗く、case/45 の生産の壁で相対 −4.1e-4 ずれた (ノット区間ごとに分けた和と本スクリプトの求積は ∫2πr ds と 1e-14 で一致;
  2026-10-13 測定)。
- **bbox** `bbox_mm`: 面ごとの `occ.getBoundingBox` の和集合 (形状の許容差だけ広がりうる。B-spline 面は制御点の包)。
  呼び出し側は x の範囲を曲線の定義域と比べ、半径方向は参考にする。求積の点から取った範囲も `bbox_sampled_mm` に出す。

usage:
  /home/sano/work/forge/.venv-mesh/bin/python cad_check_gmsh.py WALL_SURFACE.step [--out check.json] [--tol 1e-12]
    → 標準出力に 1 行の JSON: {"gmsh", "n_surfaces", "surface_types", "area_mm2", "areas_mm2", "area_occ_mass_mm2", "bbox_mm", ...}
      (--out には同じ内容を字下げして書く)
"""
import argparse
import json
import sys
from pathlib import Path

GAUSS_N = 4


def _surface_area(gmsh, tag, tol: float, max_panels: int = 1 << 17, chunk: int = 400000):
    """面 tag の面積を助変数空間の合成 Gauss で求める。戻り: (面積, 情報 dict)。"""
    import numpy as np
    g, w = np.polynomial.legendre.leggauss(GAUSS_N)
    lo, hi = gmsh.model.getParametrizationBounds(2, tag)

    def nodes(a, b, n):
        e = np.linspace(a, b, n + 1)
        m, h = 0.5 * (e[:-1] + e[1:]), 0.5 * (e[1:] - e[:-1])
        return (m[:, None] + h[:, None] * g[None, :]).ravel(), (h[:, None] * w[None, :]).ravel()

    box = [np.inf, np.inf, np.inf, -np.inf, -np.inf, -np.inf]

    def area(nu, nv):
        uu, wu = nodes(lo[0], hi[0], nu)
        vv, wv = nodes(lo[1], hi[1], nv)
        tot = 0.0
        U, V = np.meshgrid(uu, vv, indexing="ij")
        W = (wu[:, None] * wv[None, :]).ravel()
        P = np.c_[U.ravel(), V.ravel()]
        for i in range(0, len(P), chunk):
            p = P[i:i + chunk]
            d = np.asarray(gmsh.model.getDerivative(2, tag, p.ravel().tolist()), dtype=float).reshape(-1, 6)
            tot += float((W[i:i + chunk] * np.linalg.norm(np.cross(d[:, :3], d[:, 3:]), axis=1)).sum())
            x = np.asarray(gmsh.model.getValue(2, tag, p.ravel().tolist()), dtype=float).reshape(-1, 3)
            for k in range(3):
                box[k] = min(box[k], float(x[:, k].min()))
                box[k + 3] = max(box[k + 3], float(x[:, k].max()))
        return tot

    nu, nv = 4, 64
    a = area(nu, nv)
    hist = []
    rel = float("inf")
    while True:
        au, av = area(2 * nu, nv), area(nu, 2 * nv)
        du, dv = abs(au - a) / max(abs(a), 1e-300), abs(av - a) / max(abs(a), 1e-300)
        hist.append({"nu": nu, "nv": nv, "area": a, "rel_change_u": du, "rel_change_v": dv})
        rel = max(du, dv)
        if rel < tol:
            a = av if dv >= du else au
            nu, nv = (nu, 2 * nv) if dv >= du else (2 * nu, nv)
            break
        if dv >= du:
            nv, a = 2 * nv, av
        else:
            nu, a = 2 * nu, au
        if max(nu, nv) > max_panels:
            break
    return a, {"param_bounds": [list(map(float, lo)), list(map(float, hi))], "gauss_points_per_panel": GAUSS_N,
               "panels_u": nu, "panels_v": nv, "rel_change_last": rel, "converged": bool(rel < tol), "tol": tol,
               "history": hist[-3:], "bbox_sampled": box}


def check(step_path: str, tol: float = 1e-12) -> dict:
    import gmsh
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("surface_check")
        ents = gmsh.model.occ.importShapes(str(step_path))
        gmsh.model.occ.synchronize()
        surfs = gmsh.model.getEntities(2)
        areas, quad, mass, bbs, samp = [], [], [], [], []
        for _, t in surfs:
            a, info = _surface_area(gmsh, t, tol)
            samp.append(info.pop("bbox_sampled"))
            areas.append(float(a))
            quad.append(info)
            mass.append(float(gmsh.model.occ.getMass(2, t)))
            bbs.append([float(v) for v in gmsh.model.occ.getBoundingBox(2, t)])

        def union(bl):
            if not bl:
                return None
            return [min(b[0] for b in bl), min(b[1] for b in bl), min(b[2] for b in bl),
                    max(b[3] for b in bl), max(b[4] for b in bl), max(b[5] for b in bl)]

        return {"gmsh": gmsh.__version__, "step": str(Path(step_path).resolve()),
                "imported": [list(e) for e in ents], "n_surfaces": len(surfs),
                "n_curves": len(gmsh.model.getEntities(1)), "n_volumes": len(gmsh.model.getEntities(3)),
                "surface_types": [gmsh.model.getType(2, t) for _, t in surfs],
                "areas_mm2": areas, "area_mm2": float(sum(areas)), "area_quadrature": quad,
                "area_quadrature_converged": bool(quad and all(q["converged"] for q in quad)),
                "areas_occ_mass_mm2": mass, "area_occ_mass_mm2": float(sum(mass)),
                "bbox_mm": union(bbs), "bbox_per_surface_mm": bbs, "bbox_sampled_mm": union(samp),
                "note": "area_mm2 = 助変数空間の合成 Gauss 求積 (区間数を倍にして収束)。area_occ_mass_mm2 = OCC の getMass "
                        "(1 枚の長い回転面では粗い、参考)。bbox_mm = occ.getBoundingBox の和集合 (許容差だけ広がりうる)"}
    finally:
        gmsh.finalize()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("step")
    ap.add_argument("--out", default=None, help="JSON の書き出し先 (既定は標準出力だけ)")
    ap.add_argument("--tol", type=float, default=1e-12, help="面積の求積を止める相対変化 (既定 1e-12)")
    a = ap.parse_args(argv)
    if not Path(a.step).is_file():
        print(json.dumps({"error": f"STEP が無い: {a.step}"}, ensure_ascii=False))
        return 2
    res = check(a.step, a.tol)
    if a.out:
        Path(a.out).write_text(json.dumps(res, indent=1, ensure_ascii=False))
    print(json.dumps(res, ensure_ascii=False))          # 標準出力は 1 行の JSON (呼び出し側が最後の行を読めるように)
    ok = res["n_surfaces"] >= 1 and res["area_mm2"] > 0.0 and res["area_quadrature_converged"]
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
