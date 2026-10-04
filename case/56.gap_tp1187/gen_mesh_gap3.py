#!/usr/bin/env python3
"""case/56 横すきま断面 (2D 構造化, 全四角形) — **谷を「壁が側辺」になるように切る**。

前の 2 方式が失敗した理由 (2026-09-20):
  1. 縦線ブロッキング: 円弧の下端で**壁が鉛直、TFI の格子線も鉛直**になりセルが潰れる
     (skew 0.99)。折れ点をブロック角に移しても直らない — 幾何的な退化なので。
  2. gmsh BoundaryLayer + 四角形化: メッシュはできるが非整合 (hanging node) が残り、
     forge の median-dual が閉じない (`dual faces not closed`, 496 節点)。

解: **谷ブロック V1 の"側辺"を円弧そのものにする**。壁が側辺なら、その格子線は壁に沿う
ので退化しない。壁法線方向は「スリット横断」方向になり、スリット本体 V2 と同じ分布を
共有できる。

    T1  上部        x[x_in,x_out] y[h1,H]        (底辺 = P1 + V1 + P2 の上辺、すべて y=h1)
    P1  上流壁襟    x[x_in,xvu]   y[0,h1]
    P2a 下流壁襟    x[xvd,xpe]    y[0,h1]
    P2b 出口襟      x[xpe,x_out]  y[0,h1]        (slip)
    V1a 谷の口      側辺 = **円弧そのもの**、底辺 = y=-r の 1.8mm、上辺 = y=0 の 6.8mm
    V1b 谷の上      x[xvu,xvd] y[0,h1]           (矩形)
    V2  スリット    x[xu,xd] y[-D,-r]            (完全な矩形, skew 0)

**この断面は 3D (z 押し出し) でもそのまま使う。**
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
from gen_mesh_gap import CONV_CFG, CONV_BC, solve_r              # noqa: E402


def build(a):
    import gmsh
    W, r, D, H, h1 = a.w, a.r, a.depth, a.H, a.h1
    xu, xd = -0.5 * W, 0.5 * W
    xvu, xvd = xu - r, xd + r
    gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("gap2d3")
    o = gmsh.model.geo
    P = lambda x, y: o.addPoint(x, y, 0.0, 1.0)

    q = {}
    q["in0"], q["in1"], q["inH"] = P(a.x_in, 0.0), P(a.x_in, h1), P(a.x_in, H)
    q["vu0"], q["vu1"] = P(xvu, 0.0), P(xvu, h1)
    q["vd0"], q["vd1"] = P(xvd, 0.0), P(xvd, h1)
    q["pe0"], q["pe1"] = P(a.x_plate_end, 0.0), P(a.x_plate_end, h1)
    q["ou0"], q["ou1"], q["ouH"] = P(a.x_out, 0.0), P(a.x_out, h1), P(a.x_out, H)
    q["su"], q["sd"] = P(xu, -r), P(xd, -r)
    q["cu"], q["cd"] = P(xvu, -r), P(xvd, -r)
    q["bu"], q["bd"] = P(xu, -D), P(xd, -D)

    L = {}
    L["pl_u"] = o.addLine(q["in0"], q["vu0"])
    L["pl_d"] = o.addLine(q["vd0"], q["pe0"])
    L["slip"] = o.addLine(q["pe0"], q["ou0"])
    L["off_u"] = o.addLine(q["in1"], q["vu1"])
    L["off_d"] = o.addLine(q["vd1"], q["pe1"])
    L["off_s"] = o.addLine(q["pe1"], q["ou1"])
    L["v_in"] = o.addLine(q["in0"], q["in1"]);  L["v_inH"] = o.addLine(q["in1"], q["inH"])
    L["v_vu"] = o.addLine(q["vu0"], q["vu1"]);  L["v_vd"] = o.addLine(q["vd0"], q["vd1"])
    L["v_pe"] = o.addLine(q["pe0"], q["pe1"])
    L["v_ou"] = o.addLine(q["ou0"], q["ou1"]);  L["v_ouH"] = o.addLine(q["ou1"], q["ouH"])
    L["arc_u"] = o.addCircleArc(q["vu0"], q["cu"], q["su"])
    L["arc_d"] = o.addCircleArc(q["sd"], q["cd"], q["vd0"])
    L["mouth"] = o.addLine(q["su"], q["sd"])
    L["w_u"] = o.addLine(q["su"], q["bu"]); L["w_d"] = o.addLine(q["sd"], q["bd"])
    L["floor"] = o.addLine(q["bu"], q["bd"])
    L["topH"] = o.addLine(q["inH"], q["ouH"])
    L["t_vu"] = o.addLine(q["vu1"], q["vd1"])   # V1b の上辺 (y=h1, 谷をまたぐ)
    # V1 を y=0 で 2 段に切る。円弧と縦線の接合 (xvu,0) は接線が 90 度折れるので、
    # 1 ブロックの側辺の途中に置くとそこだけ skew 0.99 になる (2026-09-20 実測 226 セル)。
    L["gap_top"] = o.addLine(q["vu0"], q["vd0"])

    S = {}
    S["P1"] = o.addPlaneSurface([o.addCurveLoop([L["pl_u"], L["v_vu"], -L["off_u"], -L["v_in"]])])
    S["P2a"] = o.addPlaneSurface([o.addCurveLoop([L["pl_d"], L["v_pe"], -L["off_d"], -L["v_vd"]])])
    S["P2b"] = o.addPlaneSurface([o.addCurveLoop([L["slip"], L["v_ou"], -L["off_s"], -L["v_pe"]])])
    # V1a: 側辺 = 円弧そのもの (折れ無し)。V1b: 完全な矩形。
    S["V1a"] = o.addPlaneSurface([o.addCurveLoop(
        [L["mouth"], L["arc_d"], -L["gap_top"], L["arc_u"]])])
    S["V1b"] = o.addPlaneSurface([o.addCurveLoop(
        [L["gap_top"], L["v_vd"], -L["t_vu"], -L["v_vu"]])])
    S["V2"] = o.addPlaneSurface([o.addCurveLoop([L["w_u"], L["floor"], -L["w_d"], -L["mouth"]])])
    S["T1"] = o.addPlaneSurface([o.addCurveLoop(
        [L["off_u"], L["t_vu"], L["off_d"], L["off_s"], L["v_ouH"], -L["topH"], -L["v_inH"]])])
    o.synchronize()

    y1 = a.y1 * 1e-6
    r_bl = solve_r(y1, h1, a.n_bl)
    r_top = solve_r(a.h2 * 1e-3, H - h1, a.n_top)
    r_dep = solve_r(y1, D - r, a.n_depth)

    def tc(name, n, **kw):
        o.mesh.setTransfiniteCurve(L[name], n, **kw)

    for nm in ("v_in", "v_vu", "v_vd", "v_pe", "v_ou"):
        tc(nm, a.n_bl, meshType="Progression", coef=r_bl)
    for nm in ("v_inH", "v_ouH"):
        tc(nm, a.n_top, meshType="Progression", coef=r_top)
    tc("pl_u", a.n_up, meshType="Progression", coef=-1.02)
    tc("off_u", a.n_up, meshType="Progression", coef=-1.02)
    tc("pl_d", a.n_dn, meshType="Progression", coef=1.02)
    tc("off_d", a.n_dn, meshType="Progression", coef=1.02)
    tc("slip", a.n_buf, meshType="Progression", coef=1.05)
    tc("off_s", a.n_buf, meshType="Progression", coef=1.05)
    for nm in ("mouth", "floor", "t_vu", "gap_top"):
        tc(nm, a.n_slit, meshType="Bump", coef=a.bump_slit)
    for nm in ("arc_u", "arc_d"):
        tc(nm, a.n_arc, meshType="Progression", coef=1.0)
    for nm in ("w_u", "w_d"):
        tc(nm, a.n_depth, meshType="Progression", coef=r_dep)
    tc("topH", a.n_up + a.n_slit + a.n_dn + a.n_buf - 3)

    o.mesh.setTransfiniteSurface(S["P1"]); o.mesh.setTransfiniteSurface(S["P2a"])
    o.mesh.setTransfiniteSurface(S["P2b"]); o.mesh.setTransfiniteSurface(S["V2"])
    o.mesh.setTransfiniteSurface(S["V1a"], "Left", [q["su"], q["sd"], q["vd0"], q["vu0"]])
    o.mesh.setTransfiniteSurface(S["V1b"])
    o.mesh.setTransfiniteSurface(S["T1"], "Left", [q["in1"], q["ou1"], q["ouH"], q["inH"]])
    for s in S.values():
        o.mesh.setRecombine(2, s)
    o.synchronize()

    gmsh.model.addPhysicalGroup(1, [L["v_in"], L["v_inH"]], 1, "inlet")
    gmsh.model.addPhysicalGroup(1, [L["v_ou"], L["v_ouH"]], 2, "outlet")
    gmsh.model.addPhysicalGroup(1, [L["topH"]], 3, "top")
    gmsh.model.addPhysicalGroup(1, [L["pl_u"], L["pl_d"]], 4, "plate")
    gmsh.model.addPhysicalGroup(1, [L["slip"]], 5, "slip")
    gmsh.model.addPhysicalGroup(1, [L[k] for k in ("arc_u", "arc_d", "w_u", "w_d", "floor")],
                                6, "gap")
    gmsh.model.addPhysicalGroup(2, list(S.values()), 8, "fluid")
    gmsh.model.mesh.generate(2)
    MESH.mkdir(exist_ok=True)
    out = MESH / f"{a.tag}.msh"
    gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
    gmsh.write(str(out))
    n = len(gmsh.model.mesh.getNodes()[0])
    et, tg, _ = gmsh.model.mesh.getElements(2)
    comp = ", ".join(f"{'quad' if int(e)==3 else ('tri' if int(e)==2 else int(e))}:{len(t)}"
                     for e, t in zip(et, tg))
    gmsh.finalize()
    print(f"[{a.tag}] W={W*1e3:.2f} r={r*1e3:.2f} (r/W={r/W:.2f}) D={D*1e2:.2f}cm h1={h1*1e3:.2f}mm")
    print(f"        節点 {n}, 要素 {comp}")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--w", type=float, default=0.18e-2)
    ap.add_argument("--r", type=float, default=0.25e-2)
    ap.add_argument("--depth", type=float, default=6.35e-2)
    ap.add_argument("--H", type=float, default=0.10)
    ap.add_argument("--h1", type=float, default=1.0e-3, help="壁襟の厚さ [m]")
    ap.add_argument("--h2", type=float, default=0.05, help="襟の外側 第一セル [mm]")
    ap.add_argument("--x-in", type=float, default=-0.06)
    ap.add_argument("--x-plate-end", type=float, default=0.06)
    ap.add_argument("--x-out", type=float, default=0.09)
    ap.add_argument("--y1", type=float, default=8.0)
    ap.add_argument("--n-bl", type=int, default=49)
    ap.add_argument("--n-top", type=int, default=49)
    ap.add_argument("--n-up", type=int, default=201)
    ap.add_argument("--n-dn", type=int, default=201)
    ap.add_argument("--n-buf", type=int, default=41)
    ap.add_argument("--n-slit", type=int, default=49)
    ap.add_argument("--bump-slit", type=float, default=0.02)
    ap.add_argument("--n-arc", type=int, default=61)
    ap.add_argument("--n-depth", type=int, default=241)
    ap.add_argument("--tag", default="gap2d3")
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
        print((r2.stdout + r2.stderr)[-1500:]); sys.exit("convert failed")
    shutil.move(str(conv / "m.h5"), str(MESH / f"{a.tag}.h5"))
    qq = subprocess.run([sys.executable, str(TOOLS / "check_mesh_quality.py"),
                         str(MESH / f"{a.tag}.h5")], capture_output=True, text=True)
    print("\n".join(l for l in qq.stdout.splitlines()
                    if any(k in l for k in ("aspect", "skew", "VERDICT", "cells:"))))


if __name__ == "__main__":
    main()
