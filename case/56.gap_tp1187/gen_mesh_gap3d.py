#!/usr/bin/env python3
"""case/56 縦すきま × 横すきま交差の 3D メッシュ (全ヘキサ・構造化)。

TP-1187 の配列 (Fig 2(a)) は**半ピッチずらした 2 列**で、上流列のタイル間にある
流れ方向の縦すきま (長さ L=15.24 cm, 幅 W) が横すきまと交差し、**下流列のタイルが
その正面を塞ぐ** (前向き壁 θ=90°)。下流列の縦すきまは半ピッチ (7.62 cm) 横にある。
前向き壁の熱電対 TC 87-92 は上流縦すきまの中心線から 6 mm 以内にしか無く、
2D では衝突域を作れない (親 plan §5.1 #47/#27)。

座標: x = 流れ方向 (横すきまのスリット中心が x=0)、y = 壁法線 (タイル上面 y=0)、
z = スパン (z=0 = **上流**縦すきま中心線、z=Zh=7.62 cm = **下流**縦すきま中心線)。
両端とも対称面。足跡は

    上流縦すきま   z∈[0, W/2],      x∈[-L, -W/2]      (下流端は前向き壁で塞がれる)
    横すきま       z∈[0, Zh],       x∈[-W/2, W/2]
    下流縦すきま   z∈[Zh-W/2, Zh],  x∈[W/2, x_dn_end]

## 断面のブロッキング — `gen_mesh_gap.py` の B1-B6 を踏襲する

**襟方式 (`gen_mesh_gap4.py`) は使わない**。襟の外側輪郭はスリット口で鉛直になり、
その上のブロックの格子線も鉛直なので第一層が楔に潰れる (実測 skew 0.92、2D の
`gap2d_w018` は 0.99 で FAIL)。生産 2D run が使っている `gap2d_v4` = B1-B6 側は
skew 外れ 4 セル (0.01 %) の SOFT-PASS で、**谷を x=±W/2 の鉛直線で 3 列に割り、
円弧とスリット天の 90° の折れを同一ブロックに含めない**のが要点。

    b1a/b1b 上流上部       b2a 谷の上流肩   b2b スリット口の上   b2c 谷の下流肩
    b3a/b3b 下流上部       b4  出口バッファ (slip)             b5  スリット
    t1/t2   上流トレンチ (zone A のみ)      t3/t4 下流トレンチ (zone C のみ)

## z 方向 — 断面を 3 段に押し出す

    段 1  z∈[0, W/2]        S ∪ t1,t2    (上流縦すきま)
    段 2  z∈[W/2, Zh-W/2]   S            (タイル。縦すきま無し)
    段 3  z∈[Zh-W/2, Zh]    S' ∪ t3,t4   (下流縦すきま)

段 2/3 は前段の**天面コピーから押し出す**ので z=W/2 と z=Zh-W/2 の節点は自動的に
一致する。t1/t2 の z=W/2 側の面、t3/t4 の z=Zh-W/2 側の面がそのまま縦すきまの側壁。

## 既知の簡略 (照合量への影響を README/plan に書くこと)

1. **縦すきまの口にエッジ半径を付けていない** (z=W/2, z=Zh-W/2 は鋭いエッジ)。
   横すきま側の半径 r=0.25 cm は入れてある (TC 93/94 がその上に乗る)。縦すきま側を
   丸めるには断面が z 依存で変形するので別設計が要る。影響は縦すきまの取り込み流量を
   **過小**に見積もる方向 → 衝突加熱は下振れ。
2. **上流縦すきまの上流端 (x=-L) を無滑り壁で塞いでいる** (`--upstream wall`)。実際は
   そこに 1 本上流の横すきまがある。Table III(c)/(f) の L=15.24 / 30.48 の比較で
   衝突加熱が 26-43 % 動くので、**この端条件は A/B で確かめてから生産に使う**。
3. 下流縦すきまは `--l-dn` で打ち切って無滑り壁で塞ぐ (実長 15.24 cm)。
4. すきま床の y 分割は口側が細かく床側が粗い (2D の B1-B6 と同じ)。床は測定面でない。
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
sym:    {physID: 7, kind: slip, outputHDFflg: 0, ints: , floats: }
"""


def solve_r(y1, H, ny):
    """先頭間隔 y1・全長 H・**節点数** ny の等比公比。"""
    lo, hi = 1.0000001, 2.0
    for _ in range(300):
        rr = 0.5 * (lo + hi)
        if H * (rr - 1.0) / (rr ** (ny - 1) - 1.0) > y1:
            lo = rr
        else:
            hi = rr
    return rr


def ny_for(y1, H, ratio):
    return int(math.ceil(math.log(1.0 + H * (ratio - 1.0) / y1) / math.log(ratio))) + 1


def graded(total, first, ratio, dmax):
    """先頭 first から比 ratio で伸ばし dmax で頭打ちにした層厚 (合計 total)。"""
    d, acc = [], 0.0
    while acc < total - 1e-15:
        w = min(first * ratio ** len(d), dmax)
        d.append(w); acc += w
    k = total / acc
    return [v * k for v in d]


def cumfrac(d):
    s, c, acc = sum(d), [], 0.0
    for v in d:
        acc += v / s
        c.append(acc)
    c[-1] = 1.0
    return c


def sc(n, s):
    return max(3, int(round((n - 1) * s)) + 1)


def build(a):
    import gmsh
    import numpy as np
    W, r, D, H, h = a.w, a.r, a.depth, a.H, a.h1
    L, Zh = a.l_long, a.zhalf
    xu, xd = -0.5 * W, 0.5 * W                 # 横すきまのスリット側壁
    xvu, xvd = xu - r, xd + r                  # 谷の肩 (上面と半径の接点)
    xg, xe, xpe = -L, xd + a.l_dn, a.x_plate_end
    OPEN = (a.upstream == "open")          # 縦すきま上流端を上流横すきまへ開く (A/B の腕 B)
    x2d, x2u = xg, xg - W                  # 上流横すきまの下流壁 / 上流壁
    x2vu, x2vd = x2u - r, x2d + r          # その肩
    z1, z2 = 0.5 * W, Zh - 0.5 * W
    y1 = a.y1 * 1e-6 / a.scale

    for cond, msg in [(a.x_in < xg < xvu, "x_in < -L < xvu"),
                      (xvd < xe < xpe < a.x_out, "xvd < x_dn_end < x_plate_end < x_out"),
                      (z1 < z2, "W/2 < Zh-W/2"),
                      (h < 0.5 * W and h < r, "h1 < W/2 かつ h1 < r"),
                      (not OPEN or a.x_in < x2vu, "open: x_in < x2vu"),
                      (not OPEN or x2vd < xvu, "open: x2vd < xvu")]:
        if not cond:
            sys.exit(f"幾何が成立しない: {msg}")

    n_bl, n_arc = sc(a.n_bl, a.scale), sc(a.n_arc, a.scale)
    n_up0, n_up1 = sc(a.n_up0, a.scale), sc(a.n_up1, a.scale)
    n_dn0, n_dn1 = sc(a.n_dn0, a.scale), sc(a.n_dn1, a.scale)
    n_buf = sc(a.n_buf, a.scale)
    n_ch = max(2, (sc(a.n_core, a.scale) + 1) // 2)
    n_core = 2 * n_ch - 1
    n_dep, n_top = sc(a.n_depth, a.scale), sc(a.n_top, a.scale)
    # 谷の天面 (tV) の点数。VOUT は非構造なので**対辺に合わせる必要はなく**、
    # 細かくすると非構造メッシャが薄い三角形を量産する (scale 1.0 で prism の 85 % が
    # skew>0.9 になった)。2D の gap4 は n_core 相当の粗さで skew max 0.855。
    n_vtop = (2 * n_arc + n_core - 2 if a.n_vtop == 0
              else (sc(a.n_vtop, a.scale) if a.n_vtop > 0 else n_core))

    gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("gap3d")
    g = gmsh.model.geo
    P = lambda x, y: g.addPoint(x, y, 0.0, 1.0)
    q = {}
    # 壁 (内側) 輪郭
    q["in0"] = P(a.x_in, 0.0); q["g0"] = P(xg, 0.0); q["vu0"] = P(xvu, 0.0)
    q["su"] = P(xu, -r); q["bu"] = P(xu, -D); q["bd"] = P(xd, -D); q["sd"] = P(xd, -r)
    q["vd0"] = P(xvd, 0.0); q["e0"] = P(xe, 0.0); q["pe0"] = P(xpe, 0.0)
    q["ou0"] = P(a.x_out, 0.0)
    q["cu"] = P(xvu, -r); q["cd"] = P(xvd, -r)          # 円弧中心
    # オフセット輪郭 (襟の外側)
    q["in1"] = P(a.x_in, h); q["g1"] = P(xg, h); q["vu1"] = P(xvu, h)
    q["su1"] = P(xu + h, -r); q["bu1"] = P(xu + h, -D + h)
    q["bd1"] = P(xd - h, -D + h); q["sd1"] = P(xd - h, -r)
    q["vd1"] = P(xvd, h); q["e1"] = P(xe, h); q["pe1"] = P(xpe, h)
    q["ou1"] = P(a.x_out, h)
    # 上部
    q["inH"] = P(a.x_in, H); q["gH"] = P(xg, H); q["vuH"] = P(xvu, H)
    q["vdH"] = P(xvd, H); q["ouH"] = P(a.x_out, H)
    # トレンチ
    rc = 0.5 * W - h
    q["cc"] = P(0.0, -r); q["cm"] = P(0.0, -r - rc)
    q["gD"] = P(xg, -D)
    if OPEN:                                # 上流横すきま (腕 B)
        q["v2u0"] = P(x2vu, 0.0); q["s2u"] = P(x2u, -r); q["b2u"] = P(x2u, -D)
        q["b2d"] = P(x2d, -D); q["s2d"] = P(x2d, -r); q["v2d0"] = P(x2vd, 0.0)
        q["c2u"] = P(x2vu, -r); q["c2d"] = P(x2vd, -r)
        q["v2u1"] = P(x2vu, h); q["s2u1"] = P(x2u + h, -r); q["b2u1"] = P(x2u + h, -D + h)
        q["b2d1"] = P(x2d - h, -D + h); q["s2d1"] = P(x2d - h, -r); q["v2d1"] = P(x2vd, h)
        q["v2uH"] = P(x2vu, H); q["v2dH"] = P(x2vd, H)

    C = {}
    if OPEN:
        # in0 → v2u0 (平板) → 上流横すきま → v2d0 → vu0 (縦すきまの上の帯)
        C["pl_u0"] = g.addLine(q["in0"], q["v2u0"])
        C["arc2_u"] = g.addCircleArc(q["v2u0"], q["c2u"], q["s2u"])
        C["w2_u"] = g.addLine(q["s2u"], q["b2u"]); C["floor2"] = g.addLine(q["b2u"], q["b2d"])
        C["w2_d"] = g.addLine(q["b2d"], q["s2d"])
        C["arc2_d"] = g.addCircleArc(q["s2d"], q["c2d"], q["v2d0"])
        C["pl_u1"] = g.addLine(q["v2d0"], q["vu0"])
        C["off_u0"] = g.addLine(q["in1"], q["v2u1"])
        C["oarc2_u"] = g.addCircleArc(q["v2u1"], q["c2u"], q["s2u1"])
        C["ow2_u"] = g.addLine(q["s2u1"], q["b2u1"])
        C["ofloor2"] = g.addLine(q["b2u1"], q["b2d1"])
        C["ow2_d"] = g.addLine(q["b2d1"], q["s2d1"])
        C["oarc2_d"] = g.addCircleArc(q["s2d1"], q["c2d"], q["v2d1"])
        C["off_u1"] = g.addLine(q["v2d1"], q["vu1"])
        C["core_top2"] = g.addLine(q["s2u1"], q["s2d1"])
        C["r_s2u"] = g.addLine(q["s2u"], q["s2u1"]); C["r_b2u"] = g.addLine(q["b2u"], q["b2u1"])
        C["r_b2d"] = g.addLine(q["b2d"], q["b2d1"]); C["r_s2d"] = g.addLine(q["s2d"], q["s2d1"])
        C["v_v2u"] = g.addLine(q["v2u1"], q["v2uH"]); C["v_v2d"] = g.addLine(q["v2d1"], q["v2dH"])
        C["tU0"] = g.addLine(q["inH"], q["v2uH"]); C["tV2"] = g.addLine(q["v2uH"], q["v2dH"])
        C["tU1"] = g.addLine(q["v2dH"], q["vuH"])
        C["tr_f"] = g.addLine(q["b2d"], q["bu"])       # トレンチ床 x∈[x2d, xu]
    else:
        C["pl_u0"] = g.addLine(q["in0"], q["g0"])
        C["pl_u1"] = g.addLine(q["g0"], q["vu0"])
    C["arc_u"] = g.addCircleArc(q["vu0"], q["cu"], q["su"])
    C["w_u"] = g.addLine(q["su"], q["bu"]); C["floor"] = g.addLine(q["bu"], q["bd"])
    C["w_d"] = g.addLine(q["bd"], q["sd"])
    C["arc_d"] = g.addCircleArc(q["sd"], q["cd"], q["vd0"])
    C["pl_d0"] = g.addLine(q["vd0"], q["e0"]); C["pl_d1"] = g.addLine(q["e0"], q["pe0"])
    C["slip"] = g.addLine(q["pe0"], q["ou0"])
    C["oarc_u"] = g.addCircleArc(q["vu1"], q["cu"], q["su1"])
    C["ow_u"] = g.addLine(q["su1"], q["bu1"]); C["ofloor"] = g.addLine(q["bu1"], q["bd1"])
    C["ow_d"] = g.addLine(q["bd1"], q["sd1"])
    C["oarc_d"] = g.addCircleArc(q["sd1"], q["cd"], q["vd1"])
    C["off_d0"] = g.addLine(q["vd1"], q["e1"]); C["off_d1"] = g.addLine(q["e1"], q["pe1"])
    C["off_s"] = g.addLine(q["pe1"], q["ou1"])
    # 谷の底: 襟の円弧は su1/sd1 で**鉛直接線**なので、水平線でつなぐと 90° の折れが
    # でき、その上の非構造メッシュが薄い三角形を作る。半径 rc の 1/4 円 2 本で
    # 接線連続につなぐ。
    if a.valley_arc:
        C["ct1"] = g.addCircleArc(q["su1"], q["cc"], q["cm"])
        C["ct2"] = g.addCircleArc(q["cm"], q["cc"], q["sd1"])
        vbot = [C["ct1"], C["ct2"]]
    else:
        C["core_top"] = g.addLine(q["su1"], q["sd1"])
        vbot = [C["core_top"]]
    for k, (p0, p1) in {"r_in": ("in0", "in1"), "r_su": ("su", "su1"),
                        "r_bu": ("bu", "bu1"), "r_bd": ("bd", "bd1"), "r_sd": ("sd", "sd1"),
                        "r_e": ("e0", "e1"), "r_pe": ("pe0", "pe1"),
                        "r_ou": ("ou0", "ou1")}.items():
        C[k] = g.addLine(q[p0], q[p1])
    C["v_vu"] = g.addLine(q["vu1"], q["vuH"])
    C["v_vd"] = g.addLine(q["vd1"], q["vdH"]); C["v_inH"] = g.addLine(q["in1"], q["inH"])
    C["v_ouH"] = g.addLine(q["ou1"], q["ouH"])
    C["tV"] = g.addLine(q["vuH"], q["vdH"]); C["tD"] = g.addLine(q["vdH"], q["ouH"])
    if not OPEN:
        C["off_u0"] = g.addLine(q["in1"], q["g1"])
        C["off_u1"] = g.addLine(q["g1"], q["vu1"])
        C["v_g"] = g.addLine(q["g1"], q["gH"])
        C["tU0"] = g.addLine(q["inH"], q["gH"]); C["tU1"] = g.addLine(q["gH"], q["vuH"])
        C["r_g"] = g.addLine(q["g0"], q["g1"])
        C["tr_l"] = g.addLine(q["gD"], q["g0"]); C["tr_f"] = g.addLine(q["gD"], q["bu"])

    S = {}
    def surf(nm, loop):
        S[nm] = g.addPlaneSurface([g.addCurveLoop(loop)])
    if OPEN:
        # 上流横すきま: 主すきまと同じ襟ブロッキングを x2 station に作る
        surf("K1a", [C["pl_u0"], C["arc2_u"], C["r_s2u"], -C["oarc2_u"], -C["off_u0"],
                     -C["r_in"]])
        surf("K2b", [C["w2_u"], C["r_b2u"], -C["ow2_u"], -C["r_s2u"]])
        surf("K3b", [C["floor2"], C["r_b2d"], -C["ofloor2"], -C["r_b2u"]])
        surf("K4b", [C["w2_d"], C["r_s2d"], -C["ow2_d"], -C["r_b2d"]])
        surf("COREb", [C["ow2_u"], C["ofloor2"], C["ow2_d"], -C["core_top2"]])
        surf("VOUTb", [C["oarc2_u"], C["core_top2"], C["oarc2_d"], C["v_v2d"], -C["tV2"],
                       -C["v_v2u"]])
        surf("TU0", [C["off_u0"], C["v_v2u"], -C["tU0"], -C["v_inH"]])
        # 縦すきまの上の帯 + 主すきまの上流円弧をまとめた襟
        surf("K1b", [C["arc2_d"], C["pl_u1"], C["arc_u"], C["r_su"], -C["oarc_u"],
                     -C["off_u1"], -C["oarc2_d"], -C["r_s2d"]])
        surf("TU1", [C["off_u1"], C["v_vu"], -C["tU1"], -C["v_v2d"]])
    else:
        surf("K1a", [C["pl_u0"], C["r_g"], -C["off_u0"], -C["r_in"]])
        surf("K1b", [C["pl_u1"], C["arc_u"], C["r_su"], -C["oarc_u"], -C["off_u1"], -C["r_g"]])
        surf("TU0", [C["off_u0"], C["v_g"], -C["tU0"], -C["v_inH"]])
        surf("TU1", [C["off_u1"], C["v_vu"], -C["tU1"], -C["v_g"]])
    surf("K2", [C["w_u"], C["r_bu"], -C["ow_u"], -C["r_su"]])
    surf("K3", [C["floor"], C["r_bd"], -C["ofloor"], -C["r_bu"]])
    surf("K4", [C["w_d"], C["r_sd"], -C["ow_d"], -C["r_bd"]])
    surf("K5a", [C["arc_d"], C["pl_d0"], C["r_e"], -C["off_d0"], -C["oarc_d"], -C["r_sd"]])
    surf("K5b", [C["pl_d1"], C["r_pe"], -C["off_d1"], -C["r_e"]])
    surf("K5c", [C["slip"], C["r_ou"], -C["off_s"], -C["r_pe"]])
    surf("CORE", [C["ow_u"], C["ofloor"], C["ow_d"]] + [-c for c in vbot[::-1]])
    surf("VOUT", [C["oarc_u"]] + vbot + [C["oarc_d"], C["v_vd"], -C["tV"], -C["v_vu"]])
    surf("TD", [C["off_d0"], C["off_d1"], C["off_s"], C["v_ouH"], -C["tD"], -C["v_vd"]])
    # トレンチは 1 ブロック。上面と円弧の接点 vu0、円弧とスリット壁の接点 su の
    # うち、**接線連続な vu0 は辺の途中**に置く (角にすると内角 180° で潰れる)。
    if OPEN:
        surf("t1", [C["tr_f"], -C["w_u"], -C["arc_u"], -C["pl_u1"], -C["arc2_d"], -C["w2_d"]])
    else:
        surf("t1", [C["tr_f"], -C["w_u"], -C["arc_u"], -C["pl_u1"], -C["tr_l"]])
    g.synchronize()

    r_bl = solve_r(y1, h, n_bl)
    # 襟の外側 (VOUT) の第一層は**円弧の分割幅より十分小さく**する。同程度だと
    # 谷の上のセルが楔に潰れる (2D の gap4 は h2=0.06mm / 円弧幅 0.098mm = 0.6)。
    d_arc = r * 0.5 * math.pi / (n_arc - 1)
    h2 = a.h2_frac * d_arc
    r_top = solve_r(h2, H - h, n_top)
    rt = solve_r(y1, D, n_dep)

    def tc(nm, n, **kw):
        g.mesh.setTransfiniteCurve(C[nm], n, **kw)
    rad = ["r_in", "r_su", "r_bu", "r_bd", "r_sd", "r_e", "r_pe", "r_ou"]
    rad += (["r_s2u", "r_b2u", "r_b2d", "r_s2d"] if OPEN else ["r_g"])
    for nm in rad:
        tc(nm, n_bl, meshType="Progression", coef=r_bl)
    for nm in ("pl_u0", "off_u0", "tU0"):
        tc(nm, n_up0, meshType="Progression", coef=-1.05)
    for nm in ("pl_u1", "off_u1", "tU1"):
        tc(nm, n_up1, meshType="Progression", coef=-1.02)
    n_tr = n_up1 + (2 * n_arc - 2 if OPEN else n_arc - 1)   # トレンチ床 = 上辺の点数
    tc("tr_f", n_tr, meshType="Progression", coef=(1.0 if OPEN else -1.02))
    if OPEN:
        for nm in ("arc2_u", "arc2_d", "oarc2_u", "oarc2_d"):
            tc(nm, n_arc)
        for nm in ("w2_u", "w2_d", "ow2_u", "ow2_d"):
            tc(nm, n_dep, meshType="Progression", coef=a.dep_coef)
        for nm in ("floor2", "ofloor2", "core_top2"):
            tc(nm, n_core, meshType="Bump", coef=0.35)
        for nm in ("v_v2u", "v_v2d"):
            tc(nm, n_top, meshType="Progression", coef=r_top)
        tc("tV2", n_vtop)
    for nm in ("pl_d0", "off_d0"):
        tc(nm, n_dn0, meshType="Progression", coef=1.02)
    for nm in ("pl_d1", "off_d1"):
        tc(nm, n_dn1, meshType="Progression", coef=1.03)
    for nm in ("slip", "off_s"):
        tc(nm, n_buf, meshType="Progression", coef=1.05)
    for nm in ("arc_u", "arc_d", "oarc_u", "oarc_d"):
        tc(nm, n_arc)
    for nm in ("w_u", "w_d", "ow_u", "ow_d"):
        tc(nm, n_dep, meshType="Progression", coef=a.dep_coef)
    if not OPEN:
        tc("tr_l", n_dep, meshType="Progression", coef=a.dep_coef)
    for nm in ("floor", "ofloor"):
        tc(nm, n_core, meshType="Bump", coef=0.35)
    if a.valley_arc:
        for nm in ("ct1", "ct2"):
            tc(nm, n_ch)
    else:
        tc("core_top", n_core, meshType="Bump", coef=0.35)
    for nm in (["v_vu", "v_vd", "v_inH", "v_ouH"] + ([] if OPEN else ["v_g"])):
        tc(nm, n_top, meshType="Progression", coef=r_top)
    tc("tV", n_vtop); tc("tD", n_dn0 + n_dn1 + n_buf - 2)

    ms = g.mesh
    if OPEN:
        ms.setTransfiniteSurface(S["K1a"], "Left",
                                 [q["in0"], q["s2u"], q["s2u1"], q["in1"]])
        ms.setTransfiniteSurface(S["K1b"], "Left",
                                 [q["s2d"], q["su"], q["su1"], q["s2d1"]])
        # VOUTb も VOUT と同じ理由で構造化しない
        if a.vout_transfinite:
            ms.setTransfiniteSurface(S["VOUTb"], "Left",
                                     [q["v2u1"], q["v2d1"], q["v2dH"], q["v2uH"]])
        for k in ("K2b", "K3b", "K4b"):
            ms.setTransfiniteSurface(S[k])
        ms.setTransfiniteSurface(S["COREb"], "Left",
                                 [q["s2u1"], q["b2u1"], q["b2d1"], q["s2d1"]])
    else:
        ms.setTransfiniteSurface(S["K1b"], "Left", [q["g0"], q["su"], q["su1"], q["g1"]])
    ms.setTransfiniteSurface(S["K5a"], "Left", [q["sd"], q["e0"], q["e1"], q["sd1"]])
    # VOUT (谷の上) は transfinite で張らない。底辺が su1/sd1 で 90° 折れるので
    # 構造格子だとそこの第一層が楔に潰れる (実測 skew 0.98)。gmsh の非構造 +
    # recombine に任せると四角形主体で skew が収まる (2D の gap2d_v4 が実際そう)。
    if a.vout_transfinite:
        ms.setTransfiniteSurface(S["VOUT"], "Left", [q["vu1"], q["vd1"], q["vdH"], q["vuH"]])
    ms.setTransfiniteSurface(S["TD"], "Left", [q["vd1"], q["ou1"], q["ouH"], q["vdH"]])
    ms.setTransfiniteSurface(S["t1"], "Left",
                             ([q["b2d"], q["bu"], q["su"], q["s2d"]] if OPEN
                              else [q["gD"], q["bu"], q["su"], q["g0"]]))
    ms.setTransfiniteSurface(S["CORE"], "Left",
                             [q["su1"], q["bu1"], q["bd1"], q["sd1"]])
    for k in (["K2", "K3", "K4", "K5b", "K5c", "TU0", "TU1"] + ([] if OPEN else ["K1a"])):
        ms.setTransfiniteSurface(S[k])
    for s in S.values():
        ms.setRecombine(2, s)
    g.synchronize()

    for nm, s1, s2 in [("K1b", n_up1 + n_arc - 1 + (n_arc - 1 if OPEN else 0),
                        n_up1 + n_arc - 1 + (n_arc - 1 if OPEN else 0)),
                       ("t1", n_tr, n_tr),
                       ("K5a", n_arc + n_dn0 - 1, n_arc + n_dn0 - 1),
                       ("TD", n_dn0 + n_dn1 + n_buf - 2, n_dn0 + n_dn1 + n_buf - 2),
                       ("t1", n_up1 + n_arc - 1, n_up1 + n_arc - 1)]:
        if s1 != s2:
            sys.exit(f"transfinite 対辺点数の不一致: {nm} {s1} != {s2}")

    # --- z 方向の層 ---

    lay1 = graded(z1, y1, a.z_ratio, a.dz_gap * 1e-3 / a.scale)[::-1]
    half = graded(0.5 * (z2 - z1), lay1[-1], a.z_ratio, a.dz_far * 1e-3 / a.scale)
    lay2 = half + half[::-1]
    lay3 = graded(Zh - z2, y1, a.z_ratio, a.dz_gap * 1e-3 / a.scale)

    def ex(dts, dz, lay):
        return g.extrude(dts, 0, 0, dz, [1] * len(lay), cumfrac(lay), True)

    def tops(ret, tags):
        out, i = {}, 0
        for k, t in tags:
            nlat = len(gmsh.model.getBoundary([(2, t)], combined=False, oriented=False))
            assert ret[i][0] == 2 and ret[i + 1][0] == 3, (k, ret[i:i + 2])
            out[k] = ret[i][1]; i += 2 + nlat
        assert i == len(ret), (i, len(ret))
        return out

    keys1 = list(S)
    tags1 = [(k, S[k]) for k in keys1]
    top1 = tops(ex([(2, t) for _, t in tags1], z1, lay1), tags1)
    g.synchronize()
    keys2 = [k for k in keys1 if k != "t1"]
    tags2 = [(k, top1[k]) for k in keys2]
    top2 = tops(ex([(2, t) for _, t in tags2], z2 - z1, lay2), tags2)
    g.synchronize()

    # 下流トレンチを z=z2 のコピー曲線の上に建てる
    def pt_at(x, y, z, tol=1e-9):
        for _, t in gmsh.model.getEntities(0):
            b = gmsh.model.getBoundingBox(0, t)
            if abs(b[0] - x) < tol and abs(b[1] - y) < tol and abs(b[2] - z) < tol:
                return t
        sys.exit(f"点が見つからない ({x},{y},{z})")

    def cu_at(p0, p1):
        for _, t in gmsh.model.getEntities(1):
            if {abs(t2) for _, t2 in gmsh.model.getBoundary([(1, t)], oriented=False)} == {p0, p1}:
                return t
        sys.exit(f"曲線が見つからない {p0}-{p1}")

    q_bd, q_sd = pt_at(xd, -D, z2), pt_at(xd, -r, z2)
    q_vd, q_e = pt_at(xvd, 0.0, z2), pt_at(xe, 0.0, z2)
    c_wd, c_ad, c_pd = cu_at(q_bd, q_sd), cu_at(q_sd, q_vd), cu_at(q_vd, q_e)
    q_eD = g.addPoint(xe, -D, z2, 1.0)
    d_f, d_r = g.addLine(q_bd, q_eD), g.addLine(q_eD, q_e)
    s_t3 = g.addPlaneSurface([g.addCurveLoop([d_f, d_r, -c_pd, -c_ad, -c_wd])])
    g.synchronize()
    g.mesh.setTransfiniteCurve(d_f, n_arc + n_dn0 - 1, meshType="Progression", coef=1.02)
    g.mesh.setTransfiniteCurve(d_r, n_dep, meshType="Progression", coef=a.dep_coef)
    g.mesh.setTransfiniteSurface(s_t3, "Left", [q_bd, q_eD, q_e, q_sd])
    g.mesh.setRecombine(2, s_t3)
    g.synchronize()
    ex([(2, top2[k]) for k in keys2] + [(2, s_t3)], Zh - z2, lay3)
    g.synchronize()
    gmsh.option.setNumber("Mesh.Algorithm", a.algo2d)
    gmsh.option.setNumber("Mesh.RecombinationAlgorithm", a.recomb)
    gmsh.option.setNumber("Mesh.RecombineOptimizeTopology", 5)
    gmsh.model.mesh.generate(3)

    # --- 境界面の幾何分類 ---
    eps, groups = 1e-9, {k: [] for k in range(1, 8)}
    for _, tag in gmsh.model.getEntities(2):
        up, _ = gmsh.model.getAdjacencies(2, tag)
        if len(up) != 1:
            continue
        nt, nc, _ = gmsh.model.mesh.getNodes(2, tag, includeBoundary=True)
        if len(nt) == 0:
            continue
        nc = np.asarray(nc, dtype=float)
        xs, ys, zs = nc[0::3], nc[1::3], nc[2::3]
        flat = lambda v, c: abs(v.max() - c) < eps and abs(v.min() - c) < eps
        if flat(zs, 0.0) or flat(zs, Zh):
            k = 7
        elif flat(xs, a.x_in):
            k = 1
        elif flat(xs, a.x_out):
            k = 2
        elif flat(ys, H):
            k = 3
        elif ys.mean() > -eps and xs.mean() > xpe - eps:
            k = 5
        elif ys.mean() > -eps:
            k = 4
        else:
            k = 6
        groups[k].append(tag)
    names = {1: "inlet", 2: "outlet", 3: "top", 4: "plate", 5: "slip", 6: "gap", 7: "sym"}
    for k, tags in groups.items():
        if not tags:
            sys.exit(f"境界 {names[k]} に面が 1 つも分類されなかった")
        gmsh.model.addPhysicalGroup(2, tags, k, names[k])
    gmsh.model.addPhysicalGroup(3, [t for _, t in gmsh.model.getEntities(3)], 8, "fluid")

    gmsh.model.mesh.renumberNodes(); gmsh.model.mesh.renumberElements()
    MESH.mkdir(exist_ok=True)
    out = MESH / f"{a.tag}.msh"
    gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
    gmsh.write(str(out))
    nn = len(gmsh.model.mesh.getNodes()[0])
    et, tg, _ = gmsh.model.mesh.getElements(3)
    comp = {int(e): len(t) for e, t in zip(et, tg)}
    print(f"[{a.tag}] W={W*1e3:.2f}mm r={r*1e3:.2f}mm (r/W={r/W:.2f}) L_up={L*1e2:.2f}cm "
          f"L_dn={a.l_dn*1e2:.2f}cm Zh={Zh*1e2:.2f}cm scale={a.scale}")
    print(f"  z 層: [0,{z1*1e3:.2f}]={len(lay1)} [{z1*1e3:.2f},{z2*1e3:.2f}]={len(lay2)} "
          f"[{z2*1e3:.2f},{Zh*1e3:.2f}]={len(lay3)}   y1={y1*1e6:.2f}um h2={h2*1e6:.1f}um (円弧幅 {d_arc*1e6:.1f}um)")
    print(f"  n_bl={n_bl} n_arc={n_arc} n_core={n_core} n_depth={n_dep} "
          f"n_top={n_top} n_up1={n_up1} n_dn0={n_dn0} n_vtop={n_vtop}")
    print(f"  節点 {nn:,}  要素 " + ", ".join(
        f"{'hex' if k == 5 else ('prism' if k == 6 else k)}:{v:,}" for k, v in comp.items()))
    if set(comp) - {5}:
        print("  警告: ヘキサ以外が残っている")
    print(f"  -> {out}")
    gmsh.finalize()
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--w", type=float, default=0.18e-2, help="すきま幅 [m]")
    ap.add_argument("--r", type=float, default=0.25e-2, help="エッジ半径 [m]")
    ap.add_argument("--depth", type=float, default=6.35e-2)
    ap.add_argument("--l-long", type=float, default=15.24e-2, help="上流縦すきま長 [m]")
    ap.add_argument("--l-dn", type=float, default=5.0e-2, help="下流縦すきま長 (打ち切り) [m]")
    ap.add_argument("--upstream", choices=("wall", "open"), default="wall",
                    help="縦すきま上流端: wall=無滑り壁で塞ぐ / open=上流横すきまへ開く")
    ap.add_argument("--zhalf", type=float, default=7.62e-2, help="pitch/2 [m]")
    ap.add_argument("--H", type=float, default=0.08)
    ap.add_argument("--dz-gap", type=float, default=0.08, help="縦すきま内の dz 上限 [mm]")
    ap.add_argument("--dz-far", type=float, default=3.0, help="タイル側の dz 上限 [mm]")
    ap.add_argument("--z-ratio", type=float, default=1.15)
    ap.add_argument("--x-in", type=float, default=-0.1724)
    ap.add_argument("--x-plate-end", type=float, default=0.06)
    ap.add_argument("--x-out", type=float, default=0.09)
    ap.add_argument("--y1", type=float, default=8.0, help="壁第一間隔 [um]")
    ap.add_argument("--h1", type=float, default=6.0e-4, help="襟の厚さ [m]")
    ap.add_argument("--h2-frac", type=float, default=0.5,
                    help="襟の外側の先頭間隔 / 円弧の分割幅")
    ap.add_argument("--n-bl", type=int, default=29, help="襟の層数")
    ap.add_argument("--n-top", type=int, default=45)
    ap.add_argument("--n-core", type=int, default=21)
    ap.add_argument("--dep-coef", type=float, default=1.0)
    ap.add_argument("--n-up0", type=int, default=13)
    ap.add_argument("--n-up1", type=int, default=101)
    ap.add_argument("--n-arc", type=int, default=41)
    ap.add_argument("--n-dn0", type=int, default=71)
    ap.add_argument("--n-dn1", type=int, default=21)
    ap.add_argument("--n-buf", type=int, default=25)
    ap.add_argument("--n-depth", type=int, default=121)
    ap.add_argument("--n-vtop", type=int, default=-1,
                    help="谷の天面の点数 (-1=n_core, 0=対辺に合わせる, >0=その値)")
    ap.add_argument("--valley-arc", action="store_true", help="谷底を接線連続にする")
    ap.add_argument("--algo2d", type=int, default=6)
    ap.add_argument("--recomb", type=int, default=1)
    ap.add_argument("--vout-transfinite", action="store_true",
                    help="谷の上も構造格子にする (既定は非構造+recombine)")
    ap.add_argument("--scale", type=float, default=1.0)
    ap.add_argument("--tag", default="gap3d")
    ap.add_argument("--no-convert", action="store_true")
    a = ap.parse_args()
    msh = build(a)
    if a.no_convert:
        return
    conv = MESH / "_conv3d"; conv.mkdir(parents=True, exist_ok=True)
    (conv / "solverConfig.yaml").write_text(CONV_CFG)
    (conv / "bcondConfig.yaml").write_text(CONV_BC)
    r2 = subprocess.run([str(BUILD / "convertGmshToForge"), str(msh), "m.h5"],
                        cwd=conv, env=ENV, capture_output=True, text=True)
    (conv / f"convert_{a.tag}.log").write_text(r2.stdout + r2.stderr)
    if not (conv / "m.h5").exists():
        print((r2.stdout + r2.stderr)[-2000:]); sys.exit("convert failed")
    shutil.move(str(conv / "m.h5"), str(MESH / f"{a.tag}.h5"))
    qq = subprocess.run([sys.executable, str(TOOLS / "check_mesh_quality.py"),
                         str(MESH / f"{a.tag}.h5")], capture_output=True, text=True)
    print(qq.stdout[-1200:])


if __name__ == "__main__":
    main()
