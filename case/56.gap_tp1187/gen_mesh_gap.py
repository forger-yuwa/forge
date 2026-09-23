#!/usr/bin/env python3
"""case/56 横すきま断面 (2D, 平面) — **エッジ半径込み**。

TP-1187 のタイルはエッジ半径 r=0.25 cm で、すきま幅 W=0.10-0.41 cm より**大きい**
(run 8 の W=0.18 cm で r/W=1.39)。鋭いリップの溝ではなく「幅 6.8 mm の丸い谷の底に
1.8 mm のスリット」なので、**半径を落とすと実測ピーク (TC 93, 半径上, 4.26 q_FP) の
位置そのものが消える**。

断面ブロッキング (x-y, 上面 y=0):
    B1 上流外部  x[x_in, -W/2-r] y[0,H]
    B2a/b/c 谷の上 (3 分割: 円弧上 / スリット上 / 円弧上)
    B3 下流外部  x[+W/2+r, x_out] y[0,H]
    B4 谷        2 本の円弧で挟まれた曲線四辺形 y[-r, 0]
    B5 スリット  x[-W/2, +W/2] y[-D, -r]
    B6 出口バッファ (slip)

**このブロッキングは 3D (z 押し出し) でもそのまま使う。**
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

CONV_CFG = """mesh: {discretization: "node", nodeWallDirichlet: 1, nodeInletCornerWall: 1, meshFileName: "m.h5", valueFileName: "m.h5"}
gpu: 1
solver: "SLAU"
physProp: {thermalMethod: 0, viscMethod: 1, visc: 1.8e-5, thermCond: 0.0257, thermCondMethod: 1, prandtlLam: 0.72, cp: 1004.5, gamma: 1.4}
time:
  unsteady: 0
  dualTime: 0
  last: {nStepOuter: 10}
  deltaT: {control: 1, dt: 1e-8, cfl: 0.5, cfl_pseudo: 0.5, blockDPLUR: 1, dt_min: 1e-9, dt_max: 1.0, detectNaN: 1}
  outStepStart: 0
  outStepInterval: 10
  timeIntegration: 11
  nStepInner: 5
space: {convMethod: 0, limiter: 0}
turbulence: {model: "sst", scalarDiffusion: 1, wallTreatmentSST: 0, kInit: 100.0, omegaInit: 50000.0}
initial: "uniform_p101325_u10"
"""
CONV_BC = """inlet:  {physID: 1, kind: inlet_uniformVelocity, outputHDFflg: 0, ints: , floats: {ro: 0.0284, Ux: 2039.8, Uy: 0.0, Uz: 0.0, Ps: 1738.2, k: 100.0, omega: 50000.0}}
outlet: {physID: 2, kind: outlet_statPress, outputHDFflg: 0, ints: , floats: {Ps: 1738.2, Pt: 1738.2, Tt: 206.9}}
top:    {physID: 3, kind: slip, outputHDFflg: 0, ints: , floats: }
plate:  {physID: 4, kind: wall, outputHDFflg: 1, ints: , floats: {Ux: 0.0, Uy: 0.0, Uz: 0.0}}
slip:   {physID: 5, kind: slip, outputHDFflg: 0, ints: , floats: }
gap:    {physID: 6, kind: wall, outputHDFflg: 1, ints: , floats: {Ux: 0.0, Uy: 0.0, Uz: 0.0}}
"""


def ny_for(y1, H, r):
    return int(math.ceil(math.log(1.0 + H * (r - 1.0) / y1) / math.log(r))) + 1


def solve_r(y1, H, ny):
    lo, hi = 1.0001, 1.5
    for _ in range(200):
        r = 0.5 * (lo + hi)
        if H * (r - 1.0) / (r ** (ny - 1) - 1.0) > y1:
            lo = r
        else:
            hi = r
    return r


def build(a):
    import gmsh
    W, r, D, H = a.w, a.r, a.depth, a.H
    xu, xd = -0.5 * W, 0.5 * W                 # スリット側壁
    xvu, xvd = xu - r, xd + r                  # 谷の肩 (上面と半径の接点)
    gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("gap2d")
    g = gmsh.model.geo
    P = lambda x, y: g.addPoint(x, y, 0.0, 1.0)

    p_in0, p_in1 = P(a.x_in, 0.0), P(a.x_in, H)
    p_vu0, p_vu1 = P(xvu, 0.0), P(xvu, H)
    p_vd0, p_vd1 = P(xvd, 0.0), P(xvd, H)
    p_pe0, p_pe1 = P(a.x_plate_end, 0.0), P(a.x_plate_end, H)
    p_ou0, p_ou1 = P(a.x_out, 0.0), P(a.x_out, H)
    p_su, p_sd = P(xu, -r), P(xd, -r)          # 円弧の下端 = スリット上端
    # 谷を 3 分割するための上端点。円弧とスリット天の接合は接線が 90 度折れるので、
    # そこを 1 つのブロック内に含めるとスキューが 0.99 に達する (2026-09-20 実測)。
    p_tu, p_td = P(xu, H), P(xd, H)
    c_u, c_d = P(xvu, -r), P(xvd, -r)          # 円弧中心
    p_bu, p_bd = P(xu, -D), P(xd, -D)          # スリット床

    L = {}
    L["top_in"] = g.addLine(p_in1, p_vu1); L["top_v"] = g.addLine(p_vu1, p_vd1)
    L["top_dn"] = g.addLine(p_vd1, p_pe1);  L["top_bu"] = g.addLine(p_pe1, p_ou1)
    L["pl_in"] = g.addLine(p_in0, p_vu0);   L["pl_dn"] = g.addLine(p_vd0, p_pe0)
    L["bu"] = g.addLine(p_pe0, p_ou0)
    L["arc_u"] = g.addCircleArc(p_vu0, c_u, p_su)
    L["arc_d"] = g.addCircleArc(p_sd, c_d, p_vd0)
    L["slit_top"] = g.addLine(p_su, p_sd)
    L["w_u"] = g.addLine(p_su, p_bu); L["w_d"] = g.addLine(p_sd, p_bd)
    L["floor"] = g.addLine(p_bu, p_bd)
    L["v_su"] = g.addLine(p_su, p_tu); L["v_sd"] = g.addLine(p_sd, p_td)
    L["top_va"] = g.addLine(p_vu1, p_tu); L["top_vb"] = g.addLine(p_tu, p_td)
    L["top_vc"] = g.addLine(p_td, p_vd1)
    L["v_in"] = g.addLine(p_in0, p_in1); L["v_vu"] = g.addLine(p_vu0, p_vu1)
    L["v_vd"] = g.addLine(p_vd0, p_vd1); L["v_pe"] = g.addLine(p_pe0, p_pe1)
    L["v_ou"] = g.addLine(p_ou0, p_ou1)

    S = {}
    S["b1"] = g.addPlaneSurface([g.addCurveLoop([L["pl_in"], L["v_vu"], -L["top_in"], -L["v_in"]])])
    S["b2a"] = g.addPlaneSurface([g.addCurveLoop([L["arc_u"], L["v_su"], -L["top_va"], -L["v_vu"]])])
    S["b2b"] = g.addPlaneSurface([g.addCurveLoop([L["slit_top"], L["v_sd"], -L["top_vb"], -L["v_su"]])])
    S["b2c"] = g.addPlaneSurface([g.addCurveLoop([L["arc_d"], L["v_vd"], -L["top_vc"], -L["v_sd"]])])
    S["b3"] = g.addPlaneSurface([g.addCurveLoop([L["pl_dn"], L["v_pe"], -L["top_dn"], -L["v_vd"]])])
    S["b4"] = g.addPlaneSurface([g.addCurveLoop([L["bu"], L["v_ou"], -L["top_bu"], -L["v_pe"]])])
    S["b5"] = g.addPlaneSurface([g.addCurveLoop([L["w_u"], L["floor"], -L["w_d"], -L["slit_top"]])])
    g.synchronize()

    ny = ny_for(a.y1 * 1e-6, H, 1.08)
    ry = solve_r(a.y1 * 1e-6, H, ny)
    n_slit = a.n_slit
    n_arc = a.n_arc
    n_dep = a.n_depth
    rd_ = solve_r(a.y1 * 1e-6, D - r, n_dep)

    def tc(name, n, **kw):
        g.mesh.setTransfiniteCurve(L[name], n, **kw)
    for nm in ("v_in", "v_vu", "v_vd", "v_pe", "v_ou"):
        tc(nm, ny, meshType="Progression", coef=ry)
    tc("pl_in", a.n_up, meshType="Progression", coef=-1.02)
    tc("top_in", a.n_up, meshType="Progression", coef=-1.02)
    tc("pl_dn", a.n_dn, meshType="Progression", coef=1.02)
    tc("top_dn", a.n_dn, meshType="Progression", coef=1.02)
    tc("bu", a.n_buf, meshType="Progression", coef=1.05)
    tc("top_bu", a.n_buf, meshType="Progression", coef=1.05)
    for nm in ("arc_u", "arc_d"):
        tc(nm, n_arc, meshType="Progression", coef=1.0)
    for nm in ("slit_top", "floor"):
        tc(nm, n_slit)
    tc("top_va", n_arc); tc("top_vc", n_arc); tc("top_vb", n_slit)
    # 谷の縦線はスリット口 (y=-r) 側を細かく
    for nm in ("v_su", "v_sd"):
        tc(nm, ny, meshType="Progression", coef=ry)
    for nm in ("w_u", "w_d"):
        tc(nm, n_dep, meshType="Progression", coef=rd_)
    for k, s in S.items():
        g.mesh.setTransfiniteSurface(s)
        g.mesh.setRecombine(2, s)
    g.synchronize()

    gmsh.model.addPhysicalGroup(1, [L["v_in"]], 1, "inlet")
    gmsh.model.addPhysicalGroup(1, [L["v_ou"]], 2, "outlet")
    gmsh.model.addPhysicalGroup(1, [L[k] for k in ("top_in", "top_va", "top_vb", "top_vc",
                                                   "top_dn", "top_bu")], 3, "top")
    gmsh.model.addPhysicalGroup(1, [L["pl_in"], L["pl_dn"]], 4, "plate")
    gmsh.model.addPhysicalGroup(1, [L["bu"]], 5, "slip")
    gmsh.model.addPhysicalGroup(1, [L[k] for k in ("arc_u", "arc_d", "w_u", "w_d", "floor")], 6, "gap")
    gmsh.model.addPhysicalGroup(2, list(S.values()), 8, "fluid")
    gmsh.model.mesh.generate(2)
    out = MESH / f"{a.tag}.msh"
    gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
    gmsh.write(str(out))
    n = len(gmsh.model.mesh.getNodes()[0])
    gmsh.finalize()
    print(f"[{a.tag}] W={W*1e3:.2f} mm r={r*1e3:.2f} mm (r/W={r/W:.2f}) D={D*1e2:.2f} cm  節点 {n}")
    print(f"        ny={ny} (y1={a.y1} µm, r={ry:.5f})  スリット横断 {n_slit-1} セル  "
          f"深さ {n_dep-1}  半径上 {n_arc-1}")
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
    ap.add_argument("--n-up", type=int, default=241)
    ap.add_argument("--n-dn", type=int, default=241)
    ap.add_argument("--n-buf", type=int, default=41)
    ap.add_argument("--n-slit", type=int, default=41)
    ap.add_argument("--n-arc", type=int, default=61)
    ap.add_argument("--n-depth", type=int, default=241)
    ap.add_argument("--tag", default="gap2d")
    ap.add_argument("--no-convert", action="store_true")
    a = ap.parse_args()
    MESH.mkdir(exist_ok=True)
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
    print([l for l in q.stdout.splitlines() if "VERDICT" in l] or q.stdout[-300:])


if __name__ == "__main__":
    main()
