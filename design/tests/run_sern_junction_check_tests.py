#!/usr/bin/env python3
"""⑤ SERN 接続模型 (case/46.sern_design/cad/hex_junction_model.py) の壁第一層検査 `first_layer_check` の単体テスト
(plan tooling-sern-mesh-blocking §5.1 B4-0。python3 で実行、gmsh 不要)。

山場: 旧検査は対向節点が別の壁上の節点でも除外せず、層数も設定の記録だけだったので、全 6 面が壁で内部節点 0 の立方体でも
6 壁 ok になった (2026-10-06 plan レビュー M1)。人工の格子を渡して
(a) 全 6 面が壁の単一立方体 -> 全壁 FAIL、(b) 内部点の無い / 層が足りない薄い層 -> FAIL、(c) 正常な壁層格子 (凹の稜つき) -> PASS、
(d) 成長率 > 1.2 の層・期待タグの欠落 -> FAIL を確かめる。
"""
import sys
import types
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "case/46.sern_design/cad"))
try:
    import gmsh  # noqa: F401
except ImportError:                       # 検査の中核は gmsh を使わない。import だけ通す
    sys.modules["gmsh"] = types.ModuleType("gmsh")
from hex_junction_model import first_layer_check  # noqa: E402

FAIL = 0


def check(name, cond, detail=""):
    global FAIL
    print(("ok   " if cond else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
    if not cond:
        FAIL += 1


def box(xs, ys, zs):
    """直交格子 (gmsh のヘキサ節点順) と 6 面の境界四角形 (xmin, xmax, ymin, ymax, zmin, zmax)"""
    nx, ny, nz = len(xs), len(ys), len(zs)
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij"); xyz = np.stack([X.ravel(), Y.ravel(), Z.ravel()], 1)
    nid = lambda i, j, k: (i * ny + j) * nz + k
    hx = []
    for i in range(nx - 1):
        for j in range(ny - 1):
            for k in range(nz - 1):
                b = [nid(i, j, k), nid(i + 1, j, k), nid(i + 1, j + 1, k), nid(i, j + 1, k)]
                hx.append(b + [n + 1 for n in b])
    hx = np.array(hx, np.int64)
    def face(ax, end):
        idx = [range(nx), range(ny), range(nz)]; n = [nx, ny, nz]; idx[ax] = [0 if end == 0 else n[ax] - 1]
        a, b = [d for d in range(3) if d != ax]; q = []
        for u in range(n[a] - 1):
            for v in range(n[b] - 1):
                c = [0, 0, 0]; c[ax] = idx[ax][0]
                pts = []
                for du, dv in ((0, 0), (1, 0), (1, 1), (0, 1)):
                    c[a], c[b] = u + du, v + dv; pts.append(nid(*c))
                q.append(pts)
        return np.array(q, np.int64)
    F = {nm: face(ax, e) for nm, (ax, e) in zip(("xmin", "xmax", "ymin", "ymax", "zmin", "zmax"),
                                               [(0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1)])}
    return xyz, hx, F


def prog(h, g, n, tail=0, htail=None):
    """第一間隔 h・公比 g の n 区間、続けて一様 htail を tail 区間"""
    d = [h * g ** k for k in range(n)] + [htail or h * g ** n] * tail
    return np.r_[0.0, np.cumsum(d)]


H1, NL = 1.0e-3, 6
DR = float(prog(H1, 1.15, NL)[-1])
run = lambda xyz, hx, bq, **kw: first_layer_check(xyz, hx, bq, H1, H1, NL, NL, DR, **kw)

# (a) 全 6 面が壁の単一立方体 (内部節点 0)
xyz, hx, F = box([0.0, H1], [0.0, H1], [0.0, H1])
tags = ("ramp", "vehicle", "sidewall_in", "sidewall_out", "cowl_in", "cowl_out")
fl = run(xyz, hx, {t: F[f] for t, f in zip(tags, ("ymin", "ymax", "xmin", "xmax", "zmin", "zmax"))})
check("(a) 単一立方体: 6 壁すべて FAIL", all(not fl[t]["ok"] for t in tags), str({t: fl[t]["ok"] for t in tags}))
check("(a) 単一立方体: 全節点が評価不能", all(fl[t]["unevaluable_nodes"] == fl[t]["nodes"] == 4 for t in tags),
      str({t: fl[t]["unevaluable_nodes"] for t in tags}))

# (b1) 壁 2 面に挟まれた 1 セル厚の板 (内部点なし)
xs = np.linspace(0, 10 * H1, 6); zs = np.linspace(0, 10 * H1, 6)
xyz, hx, F = box(xs, [0.0, H1], zs)
fl = run(xyz, hx, {"ramp": F["ymin"], "vehicle": F["ymax"], "symmetry": F["xmin"]})
check("(b1) 1 セル厚の板: 両壁 FAIL・全節点評価不能",
      not fl["ramp"]["ok"] and not fl["vehicle"]["ok"] and fl["ramp"]["unevaluable_nodes"] == fl["ramp"]["nodes"],
      f"ramp {fl['ramp']['unevaluable_nodes']}/{fl['ramp']['nodes']}")
# (b2) 2 セル厚 (第一内部点はあるが層が 1 つしか無い)
xyz, hx, F = box(xs, [0.0, H1, 2 * H1], zs)
fl = run(xyz, hx, {"ramp": F["ymin"], "vehicle": F["ymax"]})
check("(b2) 2 セル厚: 第一層の比は 1 でも層数不足で FAIL",
      not fl["ramp"]["ok"] and fl["ramp"]["unevaluable_nodes"] == 0 and fl["ramp"]["layers_min"] == 1 and abs(fl["ramp"]["ratio_max"] - 1) < 1e-12,
      f"layers_min {fl['ramp']['layers_min']} / 期待 {fl['ramp']['layers_expected']}, short {fl['ramp']['layers_short']}")
# (b3) 両壁に層はあるが期待層数より薄い (壁 -> 壁の間に NL - 2 層ずつ)
ys = prog(H1, 1.15, NL - 2); ys = np.r_[ys, 2 * ys[-1] - ys[-2::-1]]
xyz, hx, F = box(xs, ys, zs)
fl = run(xyz, hx, {"ramp": F["ymin"], "vehicle": F["ymax"]})
check("(b3) 層の足りない流路: 両壁 FAIL (layers_min < 期待)",
      not fl["ramp"]["ok"] and not fl["vehicle"]["ok"] and fl["ramp"]["layers_min"] < NL,
      f"layers_min {fl['ramp']['layers_min']} / {NL}")

# (c) 正常: y = 0 (ramp) と x = 0 (sidewall_in) が凹の稜で接し、両壁とも第一層 H1・公比 1.15 で NL + 3 層。他の面は壁でない
xs = prog(H1, 1.15, NL + 3, tail=2); ys = prog(H1, 1.15, NL + 3, tail=2); zs = np.linspace(0, 20 * H1, 5)
xyz, hx, F = box(xs, ys, zs)
bq = {"ramp": F["ymin"], "sidewall_in": F["xmin"], "symmetry": F["zmin"], "far_side": F["xmax"]}
fl = run(xyz, hx, bq, expected_tags=("ramp", "sidewall_in"))
for t in ("ramp", "sidewall_in"):
    f = fl[t]
    check(f"(c) 正常格子 {t}: PASS", f["ok"],
          f"評価 {f['evaluated']}/{f['nodes']} 比 {f['ratio_min']:.4f}–{f['ratio_max']:.4f} 層 {f['layers_min']}≥{f['layers_expected']} "
          f"隣接比 {f['adjacent_ratio_max']:.4f} 稜の対角 {f['layers_skipped_diag']}")
check("(c) 稜の節点は対角の内部点で評価 (評価不能 0)", fl["ramp"]["unevaluable_nodes"] == 0 and fl["ramp"]["layers_skipped_diag"] == len(zs))

# (d1) 成長率 1.3 の層 -> 隣接間隔比で FAIL (第一層と層数は合格する)
ys = prog(H1, 1.3, NL + 3, tail=2)
xyz, hx, F = box(np.linspace(0, 10 * H1, 4), ys, np.linspace(0, 10 * H1, 4))
fl = run(xyz, hx, {"ramp": F["ymin"]})
check("(d1) 公比 1.3: FAIL (adjacent_ratio > 1.2)", not fl["ramp"]["ok"] and fl["ramp"]["adjacent_ratio_gt"] > 0 and fl["ramp"]["layers_short"] == 0,
      f"adjacent_ratio_max {fl['ramp']['adjacent_ratio_max']:.4f}")
# (d2) 公比 1.19 -> PASS (1.2 ちょうどは浮動小数の丸めで上下するので境界の内側で見る)
ys = prog(H1, 1.19, NL + 3, tail=2)
xyz, hx, F = box(np.linspace(0, 10 * H1, 4), ys, np.linspace(0, 10 * H1, 4))
fl = run(xyz, hx, {"ramp": F["ymin"]})
check("(d2) 公比 1.19: PASS", fl["ramp"]["ok"], f"adjacent_ratio_max {fl['ramp']['adjacent_ratio_max']:.6f}")
# (d3) 期待した壁タグが無い -> FAIL
fl = run(xyz, hx, {"ramp": F["ymin"]}, expected_tags=("ramp", "cowl_base"))
check("(d3) 期待タグ cowl_base の欠落: FAIL", fl["ramp"]["ok"] and not fl["cowl_base"]["ok"] and fl["cowl_base"].get("missing"))
# (d4) 端面タグは h1e と端面の期待層数で見る
fl = first_layer_check(xyz, hx, {"cowl_base": F["ymin"]}, 9 * H1, H1, 99, NL, DR)
check("(d4) 端面 (cowl_base) は h1e・nl_end で判定", fl["cowl_base"]["ok"] and fl["cowl_base"]["layers_expected"] == NL)

print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAILED'}")
sys.exit(1 if FAIL else 0)
