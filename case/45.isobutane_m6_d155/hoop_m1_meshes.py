#!/usr/bin/env python3
"""plan axisymmetric-freestream-hoop-gauge §4.11 (m1) の小さな試験の格子を作る (gmsh)。
    python hoop_m1_meshes.py <出力ディレクトリ>   → T1_tri_curved.msh・T2_mixed.msh・T3_step_multicorner.msh
物理グループ: inlet 1・outlet 2・wall 3・axis 4・fluid 5・(T3 だけ) wall2 6。軸は y = 0。"""
import math, sys
from pathlib import Path
import gmsh

OUT = Path(sys.argv[1]); OUT.mkdir(parents=True, exist_ok=True)
r_wall = lambda x: 1.0 - 0.4 * math.exp(-((x - 1.5) / 0.5) ** 2)
H = 0.05


def curved_duct(name, split_mixed):
    gmsh.model.add(name); g = gmsh.model.geo
    xs = [3.0 * i / 60 for i in range(61)]
    wp = [g.addPoint(x, r_wall(x), 0, H) for x in xs]
    a0 = g.addPoint(0.0, 0.0, 0, H); a1 = g.addPoint(3.0, 0.0, 0, H)
    inlet = g.addLine(a0, wp[0]); outlet = g.addLine(a1, wp[-1])
    if not split_mixed:
        wall = [g.addSpline(wp)]; axis = [g.addLine(a0, a1)]
        loop = g.addCurveLoop([axis[0], outlet, -wall[0], -inlet]); surfs = [g.addPlaneSurface([loop])]
    else:
        am = g.addPoint(1.5, 0.0, 0, H); im = 30   # x = 1.5 の壁の点
        wall = [g.addSpline(wp[:im + 1]), g.addSpline(wp[im:])]; axis = [g.addLine(a0, am), g.addLine(am, a1)]
        mid = g.addLine(am, wp[im])
        l1 = g.addCurveLoop([axis[0], mid, -wall[0], -inlet]); l2 = g.addCurveLoop([axis[1], outlet, -wall[1], -mid])
        surfs = [g.addPlaneSurface([l1]), g.addPlaneSurface([l2])]
    g.synchronize()
    if split_mixed:
        gmsh.model.mesh.setRecombine(2, surfs[0])   # 上流だけ四角形
    for tag, ents, nm in ((1, [inlet], "inlet"), (2, [outlet], "outlet"), (3, wall, "wall"), (4, axis, "axis")):
        gmsh.model.addPhysicalGroup(1, ents, tag, nm)
    gmsh.model.addPhysicalGroup(2, surfs, 5, "fluid")


def step_multicorner(name):
    gmsh.model.add(name); g = gmsh.model.geo
    P = [(0, 0), (3, 0), (3, 1.5), (1, 1.5), (1, 1.0), (0, 1.0)]
    h = [H, H, 0.02, 0.02, 0.02, H]
    p = [g.addPoint(x, y, 0, hh) for (x, y), hh in zip(P, h)]
    axis = g.addLine(p[0], p[1]); outlet = g.addLine(p[1], p[2]); wall2 = g.addLine(p[2], p[3])
    step = g.addLine(p[3], p[4]); wall_up = g.addLine(p[4], p[5]); inlet = g.addLine(p[5], p[0])
    loop = g.addCurveLoop([axis, outlet, wall2, step, wall_up, inlet]); s = g.addPlaneSurface([loop])
    g.synchronize()
    for tag, ents, nm in ((1, [inlet], "inlet"), (2, [outlet], "outlet"), (3, [step, wall_up], "wall"), (4, [axis], "axis"), (6, [wall2], "wall2")):
        gmsh.model.addPhysicalGroup(1, ents, tag, nm)
    gmsh.model.addPhysicalGroup(2, [s], 5, "fluid")


gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0); gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
for name, make in (("T1_tri_curved", lambda n: curved_duct(n, False)), ("T2_mixed", lambda n: curved_duct(n, True)), ("T3_step_multicorner", step_multicorner)):
    make(name)
    gmsh.model.mesh.generate(2)
    types, tags, _ = gmsh.model.mesh.getElements(2)
    cnt = {gmsh.model.mesh.getElementProperties(t)[0]: len(tg) for t, tg in zip(types, tags)}
    nn = len(gmsh.model.mesh.getNodes()[0])
    gmsh.write(str(OUT / f"{name}.msh"))
    print(f"{name}: 節点 {nn}、要素 {cnt}")
    gmsh.model.remove()
gmsh.finalize()
