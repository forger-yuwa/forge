#!/usr/bin/env python3
"""case/56 横すきま断面 (2D 構造化, 全四角形) — **接点にブロック角を置かない**。

3 方式目までの失敗と原因 (2026-09-20):
  1. 縦線ブロッキング → 円弧下端で壁も格子線も鉛直になり退化 (skew 0.99)
  2. gmsh BoundaryLayer + 四角形化 → 非整合が残り median-dual が閉じない
  3. 谷の側辺を円弧にする → **エッジ半径は平板に接する**ので、接点 (xvu,0) に
     ブロック角を置くと内角 180 度になり必ず skew→1 (226 セル)。y=0 で切っても同じ。

解: **襟 (collar) を壁に沿って連続させ、接点をブロック角にしない**。
K1 の下辺は「平板 + 円弧」を 1 本の辺として扱うので接点に角が来ない。
襟の外側の角 (xvu,h1) はオフセット円弧 (水平接線) と鉛直線の 90 度なので退化しない。

    K1  上流襟   下辺 pl_u+arc_u | 上辺 off_u+offarc_u | 左 入口(0..h1) | 右 radial_su
    K2  上流壁襟 下辺 w_u        | 上辺 offw_u
    K3  床襟     下辺 floor      | 上辺 offfloor
    K4  下流壁襟 (K2 の鏡)
    K5  下流襟   (K1 の鏡、slip まで含む)
    CORE スリット芯  x[xu+h1,xd-h1] y[-D+h1,-r]
    VOUT 谷の外側    下辺 offarc_u + y=-r + offarc_d、左右 鉛直 (h1..H)、上 y=H
    TU/TD 上部       x[x_in,xvu] / x[xvd,x_out], y[h1,H]
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
    W, r, D, H, h = a.w, a.r, a.depth, a.H, a.h1
    xu, xd = -0.5 * W, 0.5 * W
    xvu, xvd = xu - r, xd + r
    gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("gap2d4")
    o = gmsh.model.geo
    P = lambda x, y: o.addPoint(x, y, 0.0, 1.0)
    q = {}
    # 壁 (内側) 輪郭
    q["in0"] = P(a.x_in, 0.0); q["vu0"] = P(xvu, 0.0); q["su"] = P(xu, -r)
    q["bu"] = P(xu, -D); q["bd"] = P(xd, -D); q["sd"] = P(xd, -r)
    q["vd0"] = P(xvd, 0.0); q["pe0"] = P(a.x_plate_end, 0.0); q["ou0"] = P(a.x_out, 0.0)
    q["cu"] = P(xvu, -r); q["cd"] = P(xvd, -r)
    # オフセット輪郭 (襟の外側)
    q["in1"] = P(a.x_in, h); q["vu1"] = P(xvu, h); q["su1"] = P(xu + h, -r)
    q["bu1"] = P(xu + h, -D + h); q["bd1"] = P(xd - h, -D + h); q["sd1"] = P(xd - h, -r)
    q["vd1"] = P(xvd, h); q["pe1"] = P(a.x_plate_end, h); q["ou1"] = P(a.x_out, h)
    # 上部
    q["inH"] = P(a.x_in, H); q["vuH"] = P(xvu, H); q["vdH"] = P(xvd, H); q["ouH"] = P(a.x_out, H)

    L = {}
    L["pl_u"] = o.addLine(q["in0"], q["vu0"]); L["arc_u"] = o.addCircleArc(q["vu0"], q["cu"], q["su"])
    L["w_u"] = o.addLine(q["su"], q["bu"]); L["floor"] = o.addLine(q["bu"], q["bd"])
    L["w_d"] = o.addLine(q["bd"], q["sd"]); L["arc_d"] = o.addCircleArc(q["sd"], q["cd"], q["vd0"])
    L["pl_d"] = o.addLine(q["vd0"], q["pe0"]); L["slip"] = o.addLine(q["pe0"], q["ou0"])
    L["off_u"] = o.addLine(q["in1"], q["vu1"])
    L["oarc_u"] = o.addCircleArc(q["vu1"], q["cu"], q["su1"])
    L["ow_u"] = o.addLine(q["su1"], q["bu1"]); L["ofloor"] = o.addLine(q["bu1"], q["bd1"])
    L["ow_d"] = o.addLine(q["bd1"], q["sd1"])
    L["oarc_d"] = o.addCircleArc(q["sd1"], q["cd"], q["vd1"])
    L["off_d"] = o.addLine(q["vd1"], q["pe1"]); L["off_s"] = o.addLine(q["pe1"], q["ou1"])
    L["r_in"] = o.addLine(q["in0"], q["in1"]); L["r_su"] = o.addLine(q["su"], q["su1"])
    L["r_bu"] = o.addLine(q["bu"], q["bu1"]); L["r_bd"] = o.addLine(q["bd"], q["bd1"])
    L["r_sd"] = o.addLine(q["sd"], q["sd1"]); L["r_ou"] = o.addLine(q["ou0"], q["ou1"])
    L["core_top"] = o.addLine(q["su1"], q["sd1"])
    L["v_vu"] = o.addLine(q["vu1"], q["vuH"]); L["v_vd"] = o.addLine(q["vd1"], q["vdH"])
    L["v_inH"] = o.addLine(q["in1"], q["inH"]); L["v_ouH"] = o.addLine(q["ou1"], q["ouH"])
    L["tU"] = o.addLine(q["inH"], q["vuH"]); L["tV"] = o.addLine(q["vuH"], q["vdH"])
    L["tD"] = o.addLine(q["vdH"], q["ouH"])

    S = {}
    S["K1"] = o.addPlaneSurface([o.addCurveLoop(
        [L["pl_u"], L["arc_u"], L["r_su"], -L["oarc_u"], -L["off_u"], -L["r_in"]])])
    S["K2"] = o.addPlaneSurface([o.addCurveLoop([L["w_u"], L["r_bu"], -L["ow_u"], -L["r_su"]])])
    S["K3"] = o.addPlaneSurface([o.addCurveLoop([L["floor"], L["r_bd"], -L["ofloor"], -L["r_bu"]])])
    S["K4"] = o.addPlaneSurface([o.addCurveLoop([L["w_d"], L["r_sd"], -L["ow_d"], -L["r_bd"]])])
    S["K5"] = o.addPlaneSurface([o.addCurveLoop(
        [L["arc_d"], L["pl_d"], L["slip"], L["r_ou"], -L["off_s"], -L["off_d"],
         -L["oarc_d"], -L["r_sd"]])])
    S["CORE"] = o.addPlaneSurface([o.addCurveLoop(
        [L["ow_u"], L["ofloor"], L["ow_d"], -L["core_top"]])])
    S["VOUT"] = o.addPlaneSurface([o.addCurveLoop(
        [L["oarc_u"], L["core_top"], L["oarc_d"], L["v_vd"], -L["tV"], -L["v_vu"]])])
    S["TU"] = o.addPlaneSurface([o.addCurveLoop([L["off_u"], L["v_vu"], -L["tU"], -L["v_inH"]])])
    S["TD"] = o.addPlaneSurface([o.addCurveLoop(
        [L["off_d"], L["off_s"], L["v_ouH"], -L["tD"], -L["v_vd"]])])
    o.synchronize()

    y1 = a.y1 * 1e-6
    r_bl = solve_r(y1, h, a.n_bl)
    r_top = solve_r(a.h2 * 1e-3, H - h, a.n_top)
    n_side = a.n_up + a.n_arc - 1
    n_side_d = a.n_arc + a.n_dn + a.n_buf - 2

    def tc(name, n, **kw):
        o.mesh.setTransfiniteCurve(L[name], n, **kw)
    for nm in ("r_in", "r_su", "r_bu", "r_bd", "r_sd", "r_ou"):
        tc(nm, a.n_bl, meshType="Progression", coef=r_bl)
    for nm in ("pl_u", "off_u"):
        tc(nm, a.n_up, meshType="Progression", coef=-1.03)
    for nm in ("pl_d", "off_d"):
        tc(nm, a.n_dn, meshType="Progression", coef=1.03)
    for nm in ("slip", "off_s"):
        tc(nm, a.n_buf, meshType="Progression", coef=1.05)
    for nm in ("arc_u", "arc_d", "oarc_u", "oarc_d"):
        tc(nm, a.n_arc)
    for nm in ("w_u", "w_d", "ow_u", "ow_d"):
        tc(nm, a.n_depth, meshType="Progression", coef=a.dep_coef)
    for nm in ("floor", "ofloor", "core_top"):
        tc(nm, a.n_core, meshType="Bump", coef=0.35)
    for nm in ("v_vu", "v_vd", "v_inH", "v_ouH"):
        tc(nm, a.n_top, meshType="Progression", coef=r_top)
    tc("tU", a.n_up); tc("tD", a.n_dn + a.n_buf - 1); tc("tV", a.n_core)

    o.mesh.setTransfiniteSurface(S["K1"], "Left", [q["in0"], q["su"], q["su1"], q["in1"]])
    o.mesh.setTransfiniteSurface(S["K5"], "Left", [q["sd"], q["ou0"], q["ou1"], q["sd1"]])
    o.mesh.setTransfiniteSurface(S["VOUT"], "Left", [q["vu1"], q["vd1"], q["vdH"], q["vuH"]])
    o.mesh.setTransfiniteSurface(S["TD"], "Left", [q["vd1"], q["ou1"], q["ouH"], q["vdH"]])
    for k in ("K2", "K3", "K4", "CORE", "TU"):
        o.mesh.setTransfiniteSurface(S[k])
    for s in S.values():
        o.mesh.setRecombine(2, s)
    o.synchronize()

    gmsh.model.addPhysicalGroup(1, [L["r_in"], L["v_inH"]], 1, "inlet")
    gmsh.model.addPhysicalGroup(1, [L["r_ou"], L["v_ouH"]], 2, "outlet")
    gmsh.model.addPhysicalGroup(1, [L["tU"], L["tV"], L["tD"]], 3, "top")
    gmsh.model.addPhysicalGroup(1, [L["pl_u"], L["pl_d"]], 4, "plate")
    gmsh.model.addPhysicalGroup(1, [L["slip"]], 5, "slip")
    gmsh.model.addPhysicalGroup(1, [L[k] for k in ("arc_u", "w_u", "floor", "w_d", "arc_d")],
                                6, "gap")
    gmsh.model.addPhysicalGroup(2, list(S.values()), 8, "fluid")
    gmsh.model.mesh.generate(2)
    MESH.mkdir(exist_ok=True)
    out = MESH / f"{a.tag}.msh"
    gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
    gmsh.write(str(out))
    # **同一の gmsh モデルから SU2 用も書き出す** (procedures/su2-cross-check.md:
    # 「メッシュ違い」を交絡させないため同一 .geo から両方を出す)
    gmsh.write(str(MESH / f"{a.tag}.su2"))
    n = len(gmsh.model.mesh.getNodes()[0])
    et, tg, _ = gmsh.model.mesh.getElements(2)
    comp = ", ".join(f"{'quad' if int(e)==3 else ('tri' if int(e)==2 else int(e))}:{len(t)}"
                     for e, t in zip(et, tg))
    gmsh.finalize()
    print(f"[{a.tag}] W={W*1e3:.2f} r={r*1e3:.2f} (r/W={r/W:.2f}) h1={h*1e3:.2f}mm  節点 {n}, {comp}")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--w", type=float, default=0.18e-2)
    ap.add_argument("--r", type=float, default=0.25e-2)
    ap.add_argument("--depth", type=float, default=6.35e-2)
    ap.add_argument("--H", type=float, default=0.10)
    ap.add_argument("--h1", type=float, default=6.0e-4)
    ap.add_argument("--h2", type=float, default=0.06)
    ap.add_argument("--x-in", type=float, default=-0.06)
    ap.add_argument("--x-plate-end", type=float, default=0.06)
    ap.add_argument("--x-out", type=float, default=0.09)
    ap.add_argument("--y1", type=float, default=8.0)
    ap.add_argument("--n-bl", type=int, default=41)
    ap.add_argument("--n-top", type=int, default=61)
    ap.add_argument("--n-up", type=int, default=181)
    ap.add_argument("--n-dn", type=int, default=161)
    ap.add_argument("--n-buf", type=int, default=41)
    ap.add_argument("--n-arc", type=int, default=41)
    ap.add_argument("--n-core", type=int, default=25)
    ap.add_argument("--n-depth", type=int, default=201)
    ap.add_argument("--dep-coef", type=float, default=1.0)
    ap.add_argument("--tag", default="gap2d4")
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
        print((r2.stdout + r2.stderr)[-1200:]); sys.exit("convert failed")
    shutil.move(str(conv / "m.h5"), str(MESH / f"{a.tag}.h5"))
    qq = subprocess.run([sys.executable, str(TOOLS / "check_mesh_quality.py"),
                         str(MESH / f"{a.tag}.h5")], capture_output=True, text=True)
    print("\n".join(l for l in qq.stdout.splitlines()
                    if any(k in l for k in ("aspect", "skew", "VERDICT", "cells:"))))


if __name__ == "__main__":
    main()
