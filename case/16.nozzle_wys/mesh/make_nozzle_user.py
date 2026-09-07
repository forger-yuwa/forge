#!/usr/bin/env python3
"""ユーザ指定 (2026-08-19) の Wyslouzil Fig.3 ノズル形状 (`nozzle_user_profile.py`, 単位 mm) の .geo を生成する。

上壁 (半高) は 5 区間 (直管 / 直線収縮 / 3 次収縮 / 3 次膨張ブレンド / 直線膨張) を別カーブにして
x=-38 mm のキンク (0 → -0.33 の勾配不連続) をスプラインでなまさない。区間ごとに Transfinite で
点数を与え、Transfinite Surface (4 コーナー, 複合辺) で構造化 quad を作る。

usage:
  make_nozzle_user.py cell        -> nozzle_user_2d.geo         (1 層押し出し pseudo-2D, cell 用, 壁クラスタ)
  make_nozzle_user.py planar      -> nozzle_user_2d_planar.geo  (平面 2D, node 用, 壁クラスタ)
  make_nozzle_user.py planar_inv  -> nozzle_user_2d_planar_inv.geo (平面 2D, node 非粘性用, 一様)
  make_nozzle_user.py 3d          -> nozzle_user_3d.geo (12.7 mm 押し出し, 側壁も no-slip 壁, z 両側幾何集中; node/cell 共用)
単位 mm (Mesh.ScalingFactor 0.001)。throat x=0。
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nozzle_user_profile import y_user_mm, X_STR, X_CONV, X_CUB, X_TH, X_BLEND, X_EXIT

HERE = os.path.dirname(os.path.abspath(__file__))

# 区間ごとの streamwise 節点数 (区間端を共有)。目安 Δx ≈ 0.42 mm (収縮〜膨張)、直管は 0.63 mm
SEG = [  # (x0, x1, n_nodes, kind)
    (X_STR,   X_CONV,  36,  "line"),
    (X_CONV,  X_CUB,   40,  "line"),
    (X_CUB,   X_TH,    50,  "spline"),
    (X_TH,    X_BLEND, 20,  "spline"),
    (X_BLEND, X_EXIT,  205, "line"),
]
NY = 120
BUMP = float(os.environ.get("NOZZLE_BUMP", "0.004"))   # 壁集中 (両端): 既存 nozzle_fig3_2d と同じ。env NOZZLE_BUMP で上書き (壁厚感度試験用)
SPAN_MM = 12.7      # z 押し出し (cell 用 1 層 / 3d モードの全幅 = Wyslouzil 断面幅)
DENSE_DX = 0.25     # スプライン制御点間隔 [mm]
# 3d モード: z を両側幾何集中で分割 (側壁 no-slip を解像)。1 層 = 1 要素、片側 NHALF_Z 層、第一層 Z1_MM。
NHALF_Z = 16        # 片側層数 (全 2*NHALF_Z 層 = 33 節点)
Z1_MM = 0.004       # 側壁第一層厚 [mm] (= 4 µm; 2D 壁法線 y1≈0.5–5 µm と同程度)


def z_layer_heights():
    """両側幾何集中の累積正規化高さ列 (Extrude Layers 用)。公比 r は等比和 = 半幅 で二分法。"""
    a = Z1_MM / SPAN_MM
    def ssum(r):
        return a * NHALF_Z if abs(r - 1) < 1e-12 else a * (r ** NHALF_Z - 1) / (r - 1)
    lo, hi = 1.0, 10.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if ssum(mid) < 0.5 else (lo, mid)
    r = 0.5 * (lo + hi)
    pos = [0.0]
    for k in range(NHALF_Z):
        pos.append(pos[-1] + a * r ** k)
    pos[-1] = 0.5
    h = pos[1:] + [1.0 - pos[NHALF_Z - k] for k in range(1, NHALF_Z + 1)]
    h[-1] = 1.0
    return r, h


def build_geo(mode):
    L = []
    L.append("Geometry.LineNumbers = 1;")
    L.append("Mesh.ScalingFactor = 0.001; // mm -> m")
    L.append("// Wyslouzil Fig.3 nozzle, user profile 2026-08-19 (nozzle_user_profile.py):")
    L.append(f"//   throat half 2.5 mm @x=0, exit half {y_user_mm(X_EXIT):.4f} mm (A/A*={y_user_mm(X_EXIT)/2.5:.4f}), inlet half 12.7 mm, x=[{X_STR},{X_EXIT}] mm")
    pid = 1
    cid = 1

    def add_pt(x, y):
        nonlocal pid
        L.append(f"Point({pid}) = {{{x:.7f}, {y:.7f}, 0.0, 1.0}};")
        pid += 1
        return pid - 1

    def wall_curves(sign):
        """sign=+1 上壁 / -1 下壁。(curve ids [inlet->outlet], first point id, last point id)"""
        nonlocal cid
        curves = []
        first = None
        prev = None
        for (x0, x1, n, kind) in SEG:
            if prev is None:
                prev = add_pt(x0, sign * y_user_mm(x0)); first = prev
            if kind == "line":
                p1 = add_pt(x1, sign * y_user_mm(x1))
                L.append(f"Line({cid}) = {{{prev}, {p1}}};")
            else:
                nd = max(int(round((x1 - x0) / DENSE_DX)), 4)
                xs = np.linspace(x0, x1, nd + 1)[1:]
                ids = [prev]
                for x in xs:
                    ids.append(add_pt(x, sign * y_user_mm(x)))
                p1 = ids[-1]
                L.append(f"Spline({cid}) = {{{','.join(map(str, ids))}}};")
            L.append(f"Transfinite Line {{{cid}}} = {n};")
            curves.append(cid); cid += 1; prev = p1
        return curves, first, prev

    top, p_in_top, p_out_top = wall_curves(+1.0)
    bot, p_in_bot, p_out_bot = wall_curves(-1.0)
    l_out = cid; L.append(f"Line({l_out}) = {{{p_out_top}, {p_out_bot}}};"); cid += 1
    l_in = cid;  L.append(f"Line({l_in}) = {{{p_in_top}, {p_in_bot}}};"); cid += 1
    if mode == "planar_inv":
        L.append(f"Transfinite Line {{{l_in}, {l_out}}} = {NY};   // 非粘性用: 壁クラスタなし")
    else:
        L.append(f"Transfinite Line {{{l_in}, {l_out}}} = {NY} Using Bump {BUMP};")
    loop = ", ".join([str(c) for c in top] + [str(l_out)] + [str(-c) for c in reversed(bot)] + [str(-l_in)])
    L.append(f"Curve Loop(1) = {{{loop}}};")
    L.append("Plane Surface(1) = {1};")
    L.append(f"Transfinite Surface {{1}} = {{{p_in_top}, {p_out_top}, {p_out_bot}, {p_in_bot}}};")
    L.append("Recombine Surface(1);")
    ntop = sum(s[2] for s in SEG) - (len(SEG) - 1)
    L.append(f"// streamwise nodes = {ntop}, NY = {NY}")
    if mode == "cell":
        L.append("")
        L.append("// pseudo-2D (1 層押し出し), frontback=slip 用に別 physID (cell 用)")
        L.append(f"e[] = Extrude {{0, 0, {SPAN_MM}}} {{ Surface{{1}}; Layers{{1}}; Recombine; }};")
        # e[0]=back(top surface), e[1]=volume, e[2..]=lateral in loop order:
        #   top curves (len(top)) -> outlet -> bottom curves (len(bot)) -> inlet
        nt, nb = len(top), len(bot)
        top_faces = [f"e[{2+i}]" for i in range(nt)]
        out_face = f"e[{2+nt}]"
        bot_faces = [f"e[{2+nt+1+i}]" for i in range(nb)]
        in_face = f"e[{2+nt+1+nb}]"
        L.append(f'Physical Surface("inlet", 1)  = {{{in_face}}};')
        L.append(f'Physical Surface("outlet", 2) = {{{out_face}}};')
        L.append(f'Physical Surface("wall", 3)   = {{{", ".join(top_faces + bot_faces)}}};')
        L.append('Physical Surface("frontback", 4) = {1, e[0]};')
        L.append('Physical Volume("fluid", 5) = {e[1]};')
    elif mode == "3d":
        r, h = z_layer_heights()
        nl = len(h)
        L.append("")
        L.append(f"// 3D: z 押し出し {SPAN_MM} mm を両側幾何集中 {nl} 層 (第一層 {Z1_MM*1e3:.1f} µm, 公比 {r:.3f})。4 壁すべて no-slip 壁 (physID 3)")
        L.append(f"e[] = Extrude {{0, 0, {SPAN_MM}}} {{ Surface{{1}}; Layers{{ {{{', '.join(['1'] * nl)}}}, {{{', '.join(f'{v:.8f}' for v in h)}}} }}; Recombine; }};")
        nt, nb = len(top), len(bot)
        top_faces = [f"e[{2+i}]" for i in range(nt)]
        out_face = f"e[{2+nt}]"
        bot_faces = [f"e[{2+nt+1+i}]" for i in range(nb)]
        in_face = f"e[{2+nt+1+nb}]"
        L.append(f'Physical Surface("inlet", 1)  = {{{in_face}}};')
        L.append(f'Physical Surface("outlet", 2) = {{{out_face}}};')
        L.append(f'Physical Surface("wall", 3)   = {{{", ".join(top_faces + bot_faces)}, 1, e[0]}};   // 輪郭壁 + 側壁 (front/back)')
        L.append('Physical Volume("fluid", 5) = {e[1]};')
    else:
        L.append("")
        L.append("// planar 2D (node 用): 押し出しなし。physID は forge 規約 (inlet 1 / outlet 2 / wall 3 / fluid 5)")
        L.append(f'Physical Curve("inlet", 1)  = {{{l_in}}};')
        L.append(f'Physical Curve("outlet", 2) = {{{l_out}}};')
        L.append(f'Physical Curve("wall", 3)   = {{{", ".join(map(str, top + bot))}}};')
        L.append('Physical Surface("fluid", 5) = {1};')
    return "\n".join(L) + "\n"


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "cell"
    assert mode in ("cell", "planar", "planar_inv", "3d")
    name = {"cell": "nozzle_user_2d", "planar": "nozzle_user_2d_planar", "planar_inv": "nozzle_user_2d_planar_inv",
            "3d": "nozzle_user_3d"}[mode]
    suffix = os.environ.get("NOZZLE_SUFFIX", "")
    out = os.path.join(HERE, name + suffix + ".geo")
    with open(out, "w") as f:
        f.write(build_geo(mode))
    print("wrote", out)


if __name__ == "__main__":
    main()
