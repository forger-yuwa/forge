#!/usr/bin/env python3
"""case/56 横すきま断面 (2D) — **壁沿い境界層 + 四角形化**。

構造化 (縦線) ブロッキングは円弧の下端で破綻する: そこで壁が鉛直になり TFI の格子線も
鉛直なのでセルが潰れ、スキューが 0.99 に達した (2026-09-20 実測)。壁に沿った境界層を
gmsh の BoundaryLayer field で立て、残りを四角形化する。

幾何: 上面 y=0、すきま幅 W のスリット、両肩にエッジ半径 r (TP-1187 は r/W=1.39 なので
半径が支配的)。深さ D。
"""
import argparse, math, os, shutil, subprocess, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
MESH = HERE / "mesh"
ROOT = HERE.parents[1]
BUILD = ROOT / "solver_density_cuda" / "build"
TOOLS = ROOT / "solver_density_cuda" / "tools"
ENV = dict(os.environ, LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu/hdf5/serial:"
           + os.environ.get("LD_LIBRARY_PATH", ""))
sys.path.insert(0, str(HERE))
from gen_mesh_gap import CONV_CFG, CONV_BC                       # noqa: E402


def build(a):
    import gmsh
    W, r, D, H = a.w, a.r, a.depth, a.H
    xu, xd = -0.5 * W, 0.5 * W
    xvu, xvd = xu - r, xd + r
    gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("gap2d")
    o = gmsh.model.geo
    P = lambda x, y, lc: o.addPoint(x, y, 0.0, lc)

    lc_far, lc_wall, lc_gap = a.lc_far, a.lc_wall, a.lc_gap
    p = {}
    p["in0"] = P(a.x_in, 0.0, lc_far);      p["in1"] = P(a.x_in, H, lc_far)
    p["vu"] = P(xvu, 0.0, lc_gap);          p["vd"] = P(xvd, 0.0, lc_gap)
    p["su"] = P(xu, -r, lc_gap);            p["sd"] = P(xd, -r, lc_gap)
    p["cu"] = P(xvu, -r, lc_gap);           p["cd"] = P(xvd, -r, lc_gap)
    p["bu"] = P(xu, -D, lc_gap);            p["bd"] = P(xd, -D, lc_gap)
    p["pe0"] = P(a.x_plate_end, 0.0, lc_far); p["pe1"] = P(a.x_plate_end, H, lc_far)
    p["ou0"] = P(a.x_out, 0.0, lc_far);     p["ou1"] = P(a.x_out, H, lc_far)

    L = {}
    L["pl_in"] = o.addLine(p["in0"], p["vu"])
    L["arc_u"] = o.addCircleArc(p["vu"], p["cu"], p["su"])
    L["w_u"] = o.addLine(p["su"], p["bu"])
    L["floor"] = o.addLine(p["bu"], p["bd"])
    L["w_d"] = o.addLine(p["bd"], p["sd"])
    L["arc_d"] = o.addCircleArc(p["sd"], p["cd"], p["vd"])
    L["pl_dn"] = o.addLine(p["vd"], p["pe0"])
    L["slip"] = o.addLine(p["pe0"], p["ou0"])
    L["out"] = o.addLine(p["ou0"], p["ou1"])
    L["top"] = o.addLine(p["ou1"], p["in1"])
    L["inl"] = o.addLine(p["in1"], p["in0"])
    loop = [L[k] for k in ("pl_in", "arc_u", "w_u", "floor", "w_d", "arc_d",
                           "pl_dn", "slip", "out", "top", "inl")]
    surf = o.addPlaneSurface([o.addCurveLoop(loop)])
    o.synchronize()

    walls = [L[k] for k in ("pl_in", "arc_u", "w_u", "floor", "w_d", "arc_d", "pl_dn")]
    bl = gmsh.model.mesh.field.add("BoundaryLayer")
    gmsh.model.mesh.field.setNumbers(bl, "CurvesList", walls)
    gmsh.model.mesh.field.setNumber(bl, "Size", a.y1 * 1e-6)
    gmsh.model.mesh.field.setNumber(bl, "Ratio", a.ratio)
    gmsh.model.mesh.field.setNumber(bl, "Thickness", a.bl_thick)
    gmsh.model.mesh.field.setNumber(bl, "Quads", 1)
    gmsh.model.mesh.field.setNumbers(bl, "FanPointsList", [p["vu"], p["vd"]])
    gmsh.model.mesh.field.setAsBoundaryLayer(bl)

    gmsh.option.setNumber("Mesh.Algorithm", 8)          # Frontal-Delaunay for quads
    gmsh.option.setNumber("Mesh.RecombineAll", 1)
    gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 2)
    gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 1)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", a.curv)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 1)

    gmsh.model.addPhysicalGroup(1, [L["inl"]], 1, "inlet")
    gmsh.model.addPhysicalGroup(1, [L["out"]], 2, "outlet")
    gmsh.model.addPhysicalGroup(1, [L["top"]], 3, "top")
    gmsh.model.addPhysicalGroup(1, [L["pl_in"], L["pl_dn"]], 4, "plate")
    gmsh.model.addPhysicalGroup(1, [L["slip"]], 5, "slip")
    gmsh.model.addPhysicalGroup(1, [L[k] for k in ("arc_u", "w_u", "floor", "w_d", "arc_d")],
                                6, "gap")
    gmsh.model.addPhysicalGroup(2, [surf], 8, "fluid")
    gmsh.model.mesh.generate(2)
    MESH.mkdir(exist_ok=True)
    out = MESH / f"{a.tag}.msh"
    gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
    gmsh.write(str(out))
    n = len(gmsh.model.mesh.getNodes()[0])
    et, _, _ = gmsh.model.mesh.getElements(2)
    kinds = {3: "quad", 2: "tri"}
    print(f"[{a.tag}] W={W*1e3:.2f} r={r*1e3:.2f} (r/W={r/W:.2f}) D={D*1e2:.2f}cm  節点 {n}")
    print("        要素: " + ", ".join(f"{kinds.get(int(e),int(e))}" for e in et))
    gmsh.finalize()
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--w", type=float, default=0.18e-2)
    ap.add_argument("--r", type=float, default=0.25e-2)
    ap.add_argument("--depth", type=float, default=6.35e-2)
    ap.add_argument("--H", type=float, default=0.10)
    ap.add_argument("--x-in", type=float, default=-0.06)
    ap.add_argument("--x-plate-end", type=float, default=0.06)
    ap.add_argument("--x-out", type=float, default=0.09)
    ap.add_argument("--y1", type=float, default=8.0)
    ap.add_argument("--ratio", type=float, default=1.15)
    ap.add_argument("--bl-thick", type=float, default=6.0e-4)
    ap.add_argument("--lc-far", type=float, default=4.0e-3)
    ap.add_argument("--lc-wall", type=float, default=3.0e-4)
    ap.add_argument("--lc-gap", type=float, default=1.2e-4)
    ap.add_argument("--curv", type=float, default=30.0)
    ap.add_argument("--tag", default="gap2d")
    ap.add_argument("--no-convert", action="store_true")
    a = ap.parse_args()
    msh = build(a)
    if a.no_convert:
        return
    conv = MESH / "_conv"; conv.mkdir(exist_ok=True)
    (conv / "solverConfig.yaml").write_text(CONV_CFG)
    (conv / "bcondConfig.yaml").write_text(CONV_BC)
    r2 = subprocess.run([str(BUILD / "convertGmshToForge"), str(msh), "m.h5"],
                        cwd=conv, env=ENV, capture_output=True, text=True)
    (conv / f"convert_{a.tag}.log").write_text(r2.stdout + r2.stderr)
    if not (conv / "m.h5").exists():
        print(r2.stdout[-2500:], r2.stderr[-1500:]); sys.exit("convert failed")
    shutil.move(str(conv / "m.h5"), str(MESH / f"{a.tag}.h5"))
    q = subprocess.run([sys.executable, str(TOOLS / "check_mesh_quality.py"),
                        str(MESH / f"{a.tag}.h5")], capture_output=True, text=True)
    print("\n".join(l for l in q.stdout.splitlines() if any(
        k in l for k in ("aspect", "skew", "VERDICT", "cells"))))


if __name__ == "__main__":
    main()
