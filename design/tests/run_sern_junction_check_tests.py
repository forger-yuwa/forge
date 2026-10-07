#!/usr/bin/env python3
"""⑤ SERN 接続模型 (case/46.sern_design/cad/hex_junction_model.py) の検査の単体テスト
(plan tooling-sern-mesh-blocking §5.1 B4-0 / B4-1。python3 で実行、gmsh 不要)。

B4-0: 壁第一層検査 `first_layer_check`。B4-1: 出力実座標の隣接間隔比 `adjacent_spacing_check`、後流 Δx `wake_dx_check`、
形状ゲート `shape_check` (と MOC 輪郭 `Contour` / カウル外壁のオフセット `Geom.yo`) — 末尾の (e)(f)(g)。
B4-2: 継ぎ目の分布則 (片側等比 `prog_nodes`・継ぎ目の間隔の連動 `match_spacing`・z の節点数の規則) と `cowl_side` の一般壁分類 — (h)(i)。
B4-4: リング対角辺の非対称分布 `ring_diag_spacing` (壁側端の保持・コア側端の比・全長・区間数・内部隣接比) — (j)。
B4-4 (§6.5): SW 端間隔の選定則 `sw_end_spacing` (A = match_spacing・B = s/2^(1/4)) — (k)、全長診断の評価位置 `section_targets`
(要求 x = 評価 x、station に丸めない、物理 station 上は両側の領域) — (l)。
B4b-4: 1D 分布の桁あふれ (prog_r・geo_sum・march・prog_spacing) と安全な二分区間、旧と同じビット、遠方帯の実間隔列での判定 far_n — (m)。
B4b-5: 壁別の第一内部点距離 (wall_h1・プリセット・cowl_side = sidewall_out)、隅の対角の連動と解けない組合せの検出 (diag_first・corner_ratios・
station_conflicts)、SW 列などの両端分布 (two_end_*)、first_layer_check の壁別目標、事前見積もり count_mesh — (n)。

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
from hex_junction_model import (first_layer_check, adjacent_spacing_check, wake_dx_check, shape_check, Geom, P0,  # noqa: E402
                                prog_nodes, prog_r, prog_n, match_spacing, ring_diag_spacing, sw_end_spacing, section_targets,
                                END_TAGS, WALL_TAGS,
                                prog_spacing, geo_sum, far_n, march,                                           # B4b-4
                                wall_h1, PRESETS, WALL_GENERAL, WALL_ENDS, diag_first, corner_ratios, station_conflicts,  # B4b-5
                                CORNER_TOL, two_end_spacing, two_end_n_range, two_end_law, count_mesh, edge_counts)

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

# (d5) 端面の縁の節点: 端面と同一平面の内部節点が座標の丸めで n.d = +1e-18 になっても第一内部点に選ばない (B4-1 の偽 FAIL の再現)。
#      壁 (x = 0, z >= zm) の上流 z > zm は固体、z < zm は流体。下流は y が x に比例してずれる (傾いた壁) ので、
#      +x の隣より同一平面の隣 (Δz = 0.1 H1) の方が法線からの横ずれが小さい
xs = np.r_[-3 * H1, -2 * H1, -H1, prog(H1, 1.1, 8)]; ys = np.array([0.0, H1, 2 * H1]); zs = np.array([0.0, 0.5, 0.9, 1.0, 1.5, 2.0]) * H1
xyz, hx, F = box(xs, ys, zs); xyz = xyz.copy(); ny_, nz_ = len(ys), len(zs); i0, k0 = 3, 3
cen = xyz[hx].mean(1); hx = hx[~((cen[:, 0] < 0) & (cen[:, 2] > zs[k0]))]
nid = lambda i, j, k: (i * ny_ + j) * nz_ + k
q = np.array([[nid(i0, j, k), nid(i0, j + 1, k), nid(i0, j + 1, k + 1), nid(i0, j, k + 1)] for j in range(ny_ - 1) for k in range(k0, nz_ - 1)])
inplane = (xyz[:, 0] == 0) & (xyz[:, 2] < zs[k0]); xyz[inplane, 0] = 1e-18
xyz[:, 1] += 0.5 * np.maximum(xyz[:, 0], 0.0)
fl = first_layer_check(xyz, hx, {"cowl_base": q}, 9 * H1, H1, 99, 3, DR)
check("(d5) 端面の縁: 同一平面の隣 (n.d = 1e-18) を選ばず +x の隣で評価 (PASS)", fl["cowl_base"]["ok"] and fl["cowl_base"]["ratio_min"] > 0.99,
      f"ratio {fl['cowl_base']['ratio_min']:.3g}–{fl['cowl_base']['ratio_max']:.3g}")

# ---------------------------------------------------------------- B4-1
# (e) 出力実座標の隣接間隔比: 面を共有する 2 ヘキサの、面に直交する辺の長さの比 (ブロック内・継ぎ目を跨ぐものとも)
xs = prog(H1, 1.19, 8); ys = prog(H1, 1.15, 6); zs = np.linspace(0, 4 * H1, 3)
xyz, hx, F = box(xs, ys, zs)
a = adjacent_spacing_check(xyz, hx)
check("(e1) 公比 1.19 (x)・1.15 (y): PASS", a["ok"] and abs(a["x"]["ratio_max"] - 1.19) < 1e-9 and abs(a["section"]["ratio_max"] - 1.15) < 1e-9,
      f"x {a['x']['ratio_max']:.6f} section {a['section']['ratio_max']:.6f}")
# 継ぎ目で x 間隔が 2 倍に飛ぶ: 2 つのブロック (labels 0/1) をまたぐ面で検出し、ブロックの組で数える
xs2 = np.r_[np.linspace(0, 4 * H1, 5), 4 * H1 + 2 * H1 * np.arange(1, 4)]
xyz, hx, F = box(xs2, ys, zs); lab = (xyz[hx].mean(1)[:, 0] > 4 * H1).astype(int)
a = adjacent_spacing_check(xyz, hx, labels=lab, label_names=["L", "R"])
check("(e2) 継ぎ目で x 間隔 2 倍: FAIL (x 方向・比 2・組 L|R)", not a["ok"] and abs(a["x"]["ratio_max"] - 2) < 1e-9 and a["section"]["n_gt"] == 0
      and list(a.get("gt_by_block_pair", {}))[:1] == ["L|R|x"], f"x {a['x']['ratio_max']:.4f} 組 {list(a.get('gt_by_block_pair', {}))}")
# 断面内で公比 1.3 -> FAIL (section)
xyz, hx, F = box(np.linspace(0, 3 * H1, 4), prog(H1, 1.3, 5), zs)
a = adjacent_spacing_check(xyz, hx)
check("(e3) 断面内の公比 1.3: FAIL (section)", not a["ok"] and a["section"]["n_gt"] > 0 and a["x"]["n_gt"] == 0, f"section {a['section']['ratio_max']:.4f}")
# 傾いた (15°) 格子でも x 方向の分類が保たれ、比は辺長で測る
th = np.radians(15.0); xyz, hx, F = box(prog(H1, 1.1, 6), np.linspace(0, 3 * H1, 4), zs)
xyz = xyz.copy(); xyz[:, 1] += np.tan(th) * xyz[:, 0]
a = adjacent_spacing_check(xyz, hx)
check("(e4) 15° に傾けた x 方向の公比 1.1: PASS・x に分類", a["ok"] and abs(a["x"]["ratio_max"] - 1.1) < 1e-9 and a["x"]["pairs"] > 0, f"x {a['x']['ratio_max']:.6f}")

# (f) 後流 Δx: [L, L + 0.02] の全 x 辺で Δx <= cap
Hm, cap = 0.1, 1.0e-4                      # H = 0.1 m、cap = 1e-3 H
xw = np.r_[np.linspace(0.0, 0.002, 21), 0.002 + np.cumsum(1e-4 * 1.15 ** np.arange(1, 6))]
xyz, hx, F = box(xw, [0.0, 1e-3, 2e-3], [0.0, 1e-3])
w = wake_dx_check(xyz, hx, {"cowl": (0.0, 0.002)}, cap, unit=Hm)
check("(f1) 後流 20 × 1e-3 H: PASS", w["ok"] and w["cowl"]["edges"] > 0 and abs(w["cowl"]["dx_max"] - 1e-3) < 1e-12, str(w["cowl"]))
xw2 = np.r_[np.linspace(0.0, 0.0018, 13), 0.002]          # 1.5e-3 H の間隔 + 2e-3 H の間隔
xyz, hx, F = box(xw2, [0.0, 1e-3], [0.0, 1e-3])
w = wake_dx_check(xyz, hx, {"cowl": (0.0, 0.002)}, cap, unit=Hm)
check("(f2) 後流に 1.5e-3 H 超の間隔: FAIL", not w["ok"] and w["cowl"]["dx_max"] > 1e-3, f"dx_max {w['cowl']['dx_max']:.3g} H")
w = wake_dx_check(xyz, hx, {"cowl": (0.5, 0.6)}, cap, unit=Hm)
check("(f3) 区間に辺が無い: 評価不能で FAIL", not w["ok"] and w["cowl"]["edges"] == 0)

# (g) 形状ゲート: 直線の MOC 輪郭 (ランプ 15°・カウル -5°) で Geom を作り、節点を参照輪郭の上 / 少しずらして置く
H = 0.1; P = dict(P0); P["TC"] = 0.005
xr = np.linspace(-0.002, 0.2, 60); xc = np.linspace(0.0, P["LCOWL"] * H, 50)
G = Geom(P, np.stack([xr, H + np.tan(np.radians(15)) * xr], 1), np.stack([xc, -np.tan(np.radians(5)) * xc], 1))
t = P["TC"] * H; sc = -np.tan(np.radians(5)); L = P["LCOWL"] * H
check("(g0) 外壁と後縁平面の交点 y = y_c(L) - t/cos θ (オフセット曲線を延長して平面で切る)",
      abs(G.yo_te - (sc * L - t * np.sqrt(1 + sc * sc))) < 1e-15 and abs(G.yo(0.5 * L) - (sc * 0.5 * L - t * np.sqrt(1 + sc * sc))) < 1e-15,
      f"yo_te {G.yo_te:.12g} 期待 {sc * L - t * np.sqrt(1 + sc * sc):.12g}")
xq = np.linspace(0.0, L, 41); z0 = np.zeros_like(xq)
pin = np.stack([xq, G.yc(xq), z0], 1); pr = np.stack([xq, G.yr(xq), z0], 1); pout = np.stack([xq, G.yo(xq), z0], 1)
zb = np.linspace(0, 0.1, 5); yb = np.linspace(G.yo_te, float(G.yc(L)), 4)
pbase = np.array([(L, y_, z_) for y_ in yb for z_ in zb]); pend = np.array([(G.LSW, y_, z_) for y_ in (0.0, 0.05) for z_ in zb])
def gates(pin, pr, pout, pbase, pend):
    xyz = np.vstack([pin, pr, pout, pbase, pend]); o = np.cumsum([0, len(pin), len(pr), len(pout), len(pbase), len(pend)])
    bq = {t_: np.arange(o[i], o[i + 1]) for i, t_ in enumerate(("cowl_in", "ramp", "cowl_out", "cowl_base", "sidewall_end"))}
    return shape_check(xyz, bq, G)
s_ = gates(pin, pr, pout, pbase, pend)
check("(g1) 参照輪郭の上の節点・法線オフセットの外壁・平面の端面: PASS", s_["ok"],
      f"cowl_in {s_['cowl_in']['dist_max']:.2e} ramp {s_['ramp']['dist_max']:.2e} 厚さ {s_['cowl_out_thickness']['err_max']:.2e}")
s_ = gates(pin + [0, 2e-6 * H, 0], pr, pout, pbase, pend)
check("(g2) cowl_in を 2e-6 H ずらす: FAIL (cowl_in のみ)", not s_["ok"] and not s_["cowl_in"]["ok"] and s_["ramp"]["ok"], f"{s_['cowl_in']['dist_max']:.2e}")
s_ = gates(pin, pr + [0, -2e-6 * H, 0], pout, pbase, pend)
check("(g3) ramp を 2e-6 H ずらす: FAIL (ramp のみ)", not s_["ok"] and not s_["ramp"]["ok"] and s_["cowl_in"]["ok"], f"{s_['ramp']['dist_max']:.2e}")
pout_v = np.stack([xq, G.yc(xq) - t, z0], 1)              # y 方向に t ずらしただけ (法線オフセットでない) -> 法線厚 t cos θ
s_ = gates(pin, pr, pout_v, pbase, pend)
check("(g4) 外壁を y 方向オフセットにする: FAIL (厚さ)", not s_["ok"] and not s_["cowl_out_thickness"]["ok"], f"err {s_['cowl_out_thickness']['err_min']:.2e} H")
s_ = gates(pin, pr, pout, pbase + [1e-8 * H, 0, 0], pend)
check("(g5) cowl_base を x に 1e-8 H ずらす: FAIL", not s_["ok"] and not s_["cowl_base_x"]["ok"], f"{s_['cowl_base_x']['dx_max']:.2e}")
pb2 = pbase.copy(); pb2[pb2[:, 1] == pb2[:, 1].min(), 1] = float(G.yc(L)) - t          # 外側交点を y 方向オフセットの位置に置く
s_ = gates(pin, pr, pout, pb2, pend)
check("(g6) cowl_base の外側交点の y 違い: FAIL", not s_["ok"] and not s_["cowl_base_y"]["ok"], f"err_out {s_['cowl_base_y']['err_out']:.2e}")
s_ = shape_check(np.vstack([pin, pr]), {"cowl_in": np.arange(len(pin)), "ramp": len(pin) + np.arange(len(pr))}, G)
check("(g7) 必要なタグ (cowl_out ほか) の欠落: FAIL", not s_["ok"] and s_.get("cowl_out", {}).get("missing"))

# ---------------------------------------------------------------- B4-2
# (h) 継ぎ目の分布則 (1D)。単位 H
h1, g = 1.6e-4, 1.2
z = prog_nodes(1.0, 48, h1, at_end=True); dz = np.diff(z)
check("(h1) §6.3 B′: 48 区間・終端 (継ぎ目側) の間隔 = h1、公比 1.153867、最大間隔 0.133488 H",
      abs(dz[-1] - h1) < 1e-15 and abs(z[-1] - 1.0) < 1e-14 and abs(dz[0] / dz[1] - 1.153867) < 5e-7 and abs(dz.max() - 0.133488) < 5e-7,
      f"終端 {dz[-1]:.6e} 公比 {dz[0] / dz[1]:.6f} 最大 {dz.max():.6f}")
check("(h2) 片側等比の隣接比は全区間で一定 (<= 1.2)", np.allclose(dz[:-1] / dz[1:], dz[0] / dz[1], rtol=1e-12) and (dz[:-1] / dz[1:]).max() <= g)
g_t = 1.0 + 0.75 * (g - 1.0); nz = prog_n(1.0, h1, g_t) + 1
z2 = prog_nodes(1.0, nz - 1, h1, at_end=True); d2 = np.diff(z2)
check("(h3) z の節点数の規則 NZ = prog_n(W/2, h1, 1.15) + 1 = 50: 第一間隔 >= h1・長さ <= W/2 の全辺で公比 <= 1.15",
      nz == 50 and (d2[:-1] / d2[1:]).max() <= g_t and prog_r(0.92, nz - 1, 1.2 * h1) <= g_t, f"NZ {nz} 公比 {(d2[:-1] / d2[1:]).max():.6f}")
# 継ぎ目の連動: 隣の間隔の全部と比 g 以内、両側の比が等しい (幾何平均)
v = match_spacing([1.0, 1.4392], g)
check("(h4) match_spacing: 隣 1.0・1.4392 の両方と比 <= 1.2 で、両側の比が等しい", max(1.4392 / v, v / 1.0) <= g and abs(1.4392 / v - v / 1.0) < 1e-12, f"{v:.6f}")
try:
    match_spacing([1.0, 1.5], g); bad = False
except ValueError:
    bad = True
check("(h5) match_spacing: 比 1.5 > 1.2² の隣は連動できず ValueError (分布だけでは解けない継ぎ目で止める)", bad)
# 生産形状 (ランプ 15°・カウル 5°) の隅の連動の鎖: 対角 e2/e6 の第一間隔 s2/s6、e5 の端 hend、側壁外面帯 kv の間で全部比 <= 1.2
c, cc = 1 / np.cos(np.radians(15.0)), 1 / np.cos(np.radians(5.0)); s2, s6 = np.sqrt(2) * (1 + c) / 2, np.sqrt(2) * (1 + cc) / 2
hend = 0.5 * (c + cc) * P0["E5K"]; kv = P0["KV"]
F = match_spacing([hend, s2, s6, (s2 + s6) / (2 * np.sqrt(2))], g); a3 = match_spacing([s6, F, kv], g); a4 = match_spacing([s2, kv], g)
v21 = match_spacing([s6, hend, cc], g); b21 = match_spacing([F, 1.0], g)
links = [(F, hend), (F, s2), (F, s6), (F, (s2 + s6) / (2 * np.sqrt(2))), (a3, s6), (a3, F), (a3, kv), (a4, s2), (a4, kv), (v21, s6), (v21, hend), (v21, cc), (b21, F), (b21, 1.0)]
rmax = max(max(p / q, q / p) for p, q in links)
check("(h6) 隅の連動の鎖 (F・a3・a4・v21・b21) が全部比 <= 1.2", rmax <= g, f"最大比 {rmax:.4f}  F {F:.4f} a3 {a3:.4f} a4 {a4:.4f} v21 {v21:.4f} b21 {b21:.4f} (/h1)")
w_kv1 = g * 1.0 / (s2 / g) - 1
check("(h7) KV = 1 では a4 の窓 [s2/1.2, 1.2 h1] の幅が 0.1 % 未満 (KV = 1.02 を置く理由)", 0 <= w_kv1 < 1e-3, f"窓の相対幅 {w_kv1:.2e}")

# (i) cowl_side は一般壁 (h1・NL)。端面は sidewall_end・cowl_base だけ (2026-10-06 仕様訂正)
check("(i1) END_TAGS = (sidewall_end, cowl_base)、cowl_side は壁タグだが端面でない",
      set(END_TAGS) == {"sidewall_end", "cowl_base"} and "cowl_side" in WALL_TAGS and "cowl_side" not in END_TAGS)
ys = prog(H1, 1.15, NL + 3, tail=2)
xyz, hx, F_ = box(np.linspace(0, 10 * H1, 4), ys, np.linspace(0, 10 * H1, 4))
fl = first_layer_check(xyz, hx, {"cowl_side": F_["ymin"]}, H1, 9 * H1, NL, 3, DR)
check("(i2) cowl_side は h1 (端面の h1e でない) と一般壁の層数 NL で判定: PASS", fl["cowl_side"]["ok"] and fl["cowl_side"]["layers_expected"] == NL
      and abs(fl["cowl_side"]["ratio_min"] - 1) < 1e-12, f"比 {fl['cowl_side']['ratio_min']:.4f} 期待層数 {fl['cowl_side']['layers_expected']}")
fl = first_layer_check(xyz, hx, {"cowl_side": F_["ymin"]}, 9 * H1, H1, NL, 3, DR)
check("(i3) cowl_side の第一層が h1 の 1/9 (h1e 相当の設定で切った格子) なら FAIL", not fl["cowl_side"]["ok"], f"比 {fl['cowl_side']['ratio_min']:.4f}")

# (j) リング対角辺の非対称分布 (plan §6.4 の B)。対称の基準分布 (両端細分、gmsh の Bump の代わりに対数放物の間隔列) から作る
n_r = 52; tt = (np.arange(n_r) + 0.5) / n_r; base = np.exp(2.2 * (1 - (2 * tt - 1) ** 2)) ** -1
base = base / base.sum() * 0.08 * np.sqrt(2)                       # 全長 = √2 dr (dr 0.08 H)
qb = np.maximum(base[1:] / base[:-1], base[:-1] / base[1:]).max()
d = ring_diag_spacing(base, 1 / np.sqrt(2), g); qd = np.maximum(d[1:] / d[:-1], d[:-1] / d[1:]).max()
check("(j1) 非対称分布: 区間数・全長・壁側端を保持し、コア側端 / 壁側端 = 1/√2、内部隣接比 <= 1.2",
      len(d) == n_r and abs(d.sum() - base.sum()) < 1e-15 and abs(d[0] / base[0] - 1) < 1e-12 and abs(d[-1] / d[0] - 1 / np.sqrt(2)) < 1e-12 and qd <= g,
      f"全長差 {d.sum() - base.sum():.1e} 壁側 {d[0] / base[0]:.12f} 端比 {d[-1] / d[0]:.12f} 内部比 {qd:.4f} (基準 {qb:.4f})")
d1 = ring_diag_spacing(base, base[-1] / base[0], g)
check("(j2) 比 = 基準の端比 (対称なら 1) では基準分布そのもの (A と同一)", np.allclose(d1, base, rtol=1e-12, atol=0), f"最大相対差 {np.abs(d1 / base - 1).max():.1e}")
try:
    ring_diag_spacing(base[:8] / base[:8].sum(), 0.05, g); bad = False
except ValueError:
    bad = True
check("(j3) 少ない区間数で端比 0.05 は内部隣接比 > 1.2 になり ValueError (分布だけでは満たせないときは止める)", bad)

# (k) SW の端間隔の選定則 (plan §6.5)。値は §6.4 B の x/H = 1.5 断面 (/h1): s2 1.439157589、s6 1.416914595、F・kv は P0 の規則
s2k, s6k, kvk = 1.439157589, 1.416914595, P0["KV"]; Fk = 1.2054973751476779; R0 = 2 ** 0.25
a3A, a4A = sw_end_spacing(s2k, s6k, Fk, kvk, g, 0)
check("(k1) 規則 0 (A) は従来の match_spacing と同一", a3A == match_spacing([s6k, Fk, kvk], g) and a4A == match_spacing([s2k, kvk], g), f"a3 {a3A:.6f} a4 {a4A:.6f}")
a3B, a4B = sw_end_spacing(s2k, s6k, Fk, kvk, g, 1); tN = (s2k + s6k) / (2 * np.sqrt(2))
check("(k2) 規則 1 (B): 両隅 s6/a3・s2/a4 と e5 中央の平均間隔 / t_N がどれも 2^(1/4) (codex 値 a3 1.191478、a4 1.210182)",
      abs(s6k / a3B - R0) < 1e-14 and abs(s2k / a4B - R0) < 1e-14 and abs(0.5 * (a3B + a4B) / tN - R0) < 1e-14 and abs(a3B - 1.191478404) < 1e-9 and abs(a4B - 1.210182457) < 1e-9,
      f"a3 {a3B:.9f} a4 {a4B:.9f} 中央比 {0.5 * (a3B + a4B) / tN:.9f}")
rA = 0.5 * (a3A + a4A) / tN
check("(k3) A の e5 中央の平均間隔比は 1.1952 (§6.4 B の N_SIDE|SW を再現)、B はそれより小さい", abs(rA - 1.195204443) < 1e-8 and 0.5 * (a3B + a4B) / tN < rA, f"A {rA:.9f}")
rk = max(max(a3B / kvk, kvk / a3B), max(a4B / kvk, kvk / a4B), max(a3B / Fk, Fk / a3B))
check("(k4) B の端間隔は外側帯 kv・e8 の F とも比 <= 1.2 (codex: kv に対して最大 1.186)", rk <= g, f"最大比 {rk:.6f}")
try:
    sw_end_spacing(s2k, s6k, Fk, kvk, g, 2); bad = False
except ValueError:
    bad = True
check("(k5) 未定義の規則番号は ValueError", bad)

# (l) 全長診断の評価位置: 要求 x/H をそのまま評価位置にする (最近傍 station への丸めなし)
Hl = 0.1; ivl = [(0.0, 0.08, "A"), (0.08, 0.12, "B"), (0.12, 0.15, "C"), (0.15, 1.0091642885701801, "C")]
xq = [0.0, 0.8, 1.0, 1.2, 1.37, 3.0, 8.0, 10.0, 10.0916428857018]
tg = section_targets(xq, ivl, Hl)
check("(l1) 評価 x = 要求 x/H × H (全要求点で一致、丸めなし)", [t[0] for t in tg] == xq and all(t[1] == t[0] * Hl for t in tg),
      " ".join(f"{t[0]:g}->{t[1] / Hl:.15g}" for t in tg))
rgs = {t[0]: t[2] for t in tg}
check("(l2) 物理 station 上 (x/H 0.8・1.2) は両側の領域、区間内部は 1 つ、x/H 1.5 の station は C 同士",
      rgs[0.8] == {"A", "B"} and rgs[1.2] == {"B", "C"} and rgs[1.0] == {"B"} and rgs[1.37] == {"C"} and rgs[10.0] == {"C"} and rgs[0.0] == {"A"}, str({k: sorted(v) for k, v in rgs.items()}))
try:
    section_targets([10.2], ivl, Hl); bad = False
except ValueError:
    bad = True
check("(l3) 診断範囲 (参照輪郭の終端) の外は ValueError", bad)

# ---------------------------------------------------------------- B4b-4 (plan §5.1、§4.15 Minor・§4.16 (6)): 1D 分布の桁あふれと最大間隔の判定
def prog_r_legacy(L, n, h):
    """B4b-4 より前の prog_r (比較用の写し)。r**n があふれると OverflowError、区間 [1+1e-12, 50] / [1e-3, 1-1e-12] の外の根には端を返す"""
    if abs(n * h - L) < 1e-12 * L: return 1.0
    lo, hi = (1.0 + 1e-12, 50.0) if n * h < L else (1e-3, 1.0 - 1e-12)
    for _ in range(200):
        r = 0.5 * (lo + hi)
        lo, hi = (r, hi) if h * (r ** n - 1) / (r - 1) < L else (lo, r)
    return 0.5 * (lo + hi)


try:
    prog_r_legacy(1.224, 220, 0.001); ovf = False
except OverflowError:
    ovf = True
r220 = prog_r(1.224, 220, 0.001); s220 = sum(0.001 * r220 ** i for i in range(220))
check("(m1) prog_r(1.224, 220, 0.001): 旧は OverflowError、新は解けて Σ h r^i = L (相対 1e-12)", ovf and r220 > 1 and abs(s220 / 1.224 - 1) < 1e-12,
      f"r {r220:.12f} Σ {s220:.15g}")
r2 = prog_r(1.0, 2, 0.001)
check("(m2) 安全な二分区間: prog_r(1.0, 2, 0.001) = 999 (旧は区間の端 50 を返していた)", abs(r2 - 999.0) < 1e-9 and prog_r_legacy(1.0, 2, 0.001) == 50.0, f"{r2:.12g}")
r3 = prog_r(1.0005, 5, 1.0); s3 = sum(r3 ** i for i in range(5))
check("(m3) r < 1e-3 の根 (旧の下端 1e-3 の外): prog_r(1.0005, 5, 1.0) で Σ = L", r3 < 1e-3 and abs(s3 - 1.0005) < 1e-12, f"r {r3:.6g}")
bad = 0
for args in ((1.0, 1, 0.5), (1.0, 3, 1.5)):
    try:
        prog_r(*args)
    except ValueError:
        bad += 1
check("(m4) 解の無い入力 (n = 1 で h != L、h >= L) は ValueError (旧は区間の端を黙って返した)", bad == 2 and prog_r(0.5, 1, 0.5) == 1.0)
# 既定の生成で使う範囲 (L: 壁帯 DR・z 幅 W/2・遠方帯・端面区間、h: 2 µm〜1 mm の /H 相当) で、旧が桁あふれせず根が旧区間内なら旧と同じビット
nd = ne = 0
for L in (0.008, 0.0915, 0.1, 0.03, 0.002, 1.0e-3):
    for h in (2e-6, 7e-6, 1.6e-5, 6.4e-5, 1e-4, 3e-4):
        for n in range(2, 90):
            try:
                old = prog_r_legacy(L, n, h)
            except OverflowError:
                continue
            if old in (50.0,) or old <= 1.0005e-3: continue
            ne += 1; nd += prog_r(L, n, h) != old
check("(m5) 旧が桁あふれせず根が旧区間内の入力では旧と同じビット (既定の生成結果のビット同一の前提)", nd == 0 and ne > 1000, f"{ne} 組中 不一致 {nd}")
d = prog_spacing(1.224, 220, 0.001)
check("(m6) prog_spacing(1.224, 220, 0.001): 有限・第一間隔 h・全長 L", np.isfinite(d).all() and abs(d[0] / 0.001 - 1) < 1e-12 and abs(d.sum() / 1.224 - 1) < 1e-14,
      f"末尾 {d[-1]:.6g} 最大 {d.max():.6g}")
check("(m7) geo_sum: q^n があふれる所は inf (旧式は OverflowError)、それ以外は旧式 (q^n - 1)/(q - 1) と同じ値",
      geo_sum(1.2, 5000) == float("inf") and geo_sum(1.15, 40) == (1.15 ** 40 - 1) / (1.15 - 1))
# 遠方帯の区間数: 旧は 1 組の概算と末尾の近似式 L (1 - 1/r) で判定 (既定の接続模型で近似 0.0991 H、実際は 0.110 H > HFAR 0.10 H)。新は実間隔列の最大
Lf, h0, hmax, gt = 0.915, 0.012, 0.10, 1.15                     # /H。帯長 0.915 H・第一間隔 0.012 H (壁帯の外の間隔の桁)
n_old = prog_n(Lf, h0, gt) + 1
while Lf * (1 - 1 / prog_r_legacy(Lf, n_old, h0)) > hmax and n_old < 400: n_old += 1
n_new, mx = far_n([(Lf, h0)], hmax, gt)
check("(m8) far_n は実間隔列の最大で判定: 旧の近似判定の区間数では実最大が上限を超え、新は上限以内",
      prog_spacing(Lf, n_old, h0).max() > hmax and mx <= hmax and n_new > n_old and abs(mx - prog_spacing(Lf, n_new, h0).max()) < 1e-15,
      f"旧 n {n_old} (実最大 {prog_spacing(Lf, n_old, h0).max():.4f}) → 新 n {n_new} (実最大 {mx:.4f})")
n2_, mx2 = far_n([(Lf, h0), (0.5, 0.012), (0.915, 0.03)], hmax, gt)
check("(m9) far_n は全組の最大で判定 (帯長・第一間隔の違う組を同じ区間数で)", mx2 <= hmax and all(prog_spacing(L_, n2_, h_).max() <= hmax for L_, h_ in ((Lf, h0), (0.5, 0.012), (0.915, 0.03))),
      f"n {n2_} 最大 {mx2:.4f}")
try:
    res = march([0.002, 8.6], 1e-4, 1.15, 1.2, 1e-3); ovf_m = False
except OverflowError:
    ovf_m = True
check("(m10) march: 長い区間 (8.6、上限 1e-3 で 1 万区間級) でも OverflowError にならない (旧の (q^n - 1)/(q - 1) はあふれる)", not ovf_m and res[1][0] > 3890,
      f"区間数 {res[1][0] if not ovf_m else '-'}")

# ---------------------------------------------------------------- B4b-5 (plan §5.1、§4.16): 壁別の第一内部点距離
P_ = dict(P0); hw0 = wall_h1(P_)
check("(n1) 既定: 一般壁は H1、端面 (sidewall_end・cowl_base) は H1_END", all(hw0[t] == P0["H1"] for t in WALL_GENERAL) and all(hw0[t] == P0["H1_END"] for t in WALL_ENDS))
P_ = dict(P0, H1_sidewall_out=6.4e-4); hw1 = wall_h1(P_)
check("(n2) cowl_side を指定しなければ sidewall_out に従う (同一平面 z = z_o で z 分布を共有)", hw1["cowl_side"] == 6.4e-4 and hw1["sidewall_out"] == 6.4e-4)
bad = 0
for extra in (dict(H1_sidewall_out=6.4e-4, H1_cowl_side=1.6e-4), dict(H1_vehicle_top=1e-4), dict(H1_ramp=1.6e-4, H1_vehicle=6.4e-4)):
    try:
        wall_h1(dict(P0, **extra))
    except ValueError:
        bad += 1
check("(n3) cowl_side と sidewall_out・vehicle と ramp (同一平面) に違う値 / 未知の壁タグの指定は ValueError", bad == 3)
_wv = wall_h1(dict(P0, H1_ramp=1.6e-4))
check("(n3b) vehicle は指定しなければ ramp に従う (plan §4.16 の仕様訂正)", _wv["vehicle"] == _wv["ramp"] == 1.6e-4)
um = lambda d, t: round(d["H1_" + t] * 0.1 * 1e6, 9)          # /H -> µm (H = 0.1 m)
pd, pf = PRESETS["design"], PRESETS["final"]
check("(n4) プリセット design (plan §4.16 の表): 細かくする壁 16 µm、緩める壁 64 µm (帯 vehicle は ramp に従う)、端面 0.1 mm、後流 Δx <= 0.1 mm を 0.02 H、HX 0.10・HX_FAR 0.16 H",
      all(um(pd, t) == 16 for t in ("ramp", "cowl_in", "cowl_out", "sidewall_in")) and all(um(pd, t) == 64 for t in ("sidewall_out", "cowl_side")) and "H1_vehicle" not in pd
      and wall_h1(dict(P0, **pd))["vehicle"] == pd["H1_ramp"]
      and all(um(pd, t) == 100 for t in WALL_ENDS) and pd["LWAKE"] == 0.02 and pd["DXW"] == 1e-3 and pd["HX"] == 0.10 and pd["HX_FAR"] == 0.16 and pd["TC"] == 0.005)
check("(n5) プリセット final (暫定): cowl_in 3・ramp 7・sidewall_in 2・cowl_out 16 µm、緩める壁・端面・後流は design と同じ、HX 0.05・HX_FAR 0.08 H",
      (um(pf, "cowl_in"), um(pf, "ramp"), um(pf, "sidewall_in"), um(pf, "cowl_out")) == (3, 7, 2, 16)
      and all(pf[k] == pd[k] for k in pd if k not in ("H1_ramp", "H1_cowl_in", "H1_cowl_out", "H1_sidewall_in", "HX", "HX_FAR")) and (pf["HX"], pf["HX_FAR"]) == (0.05, 0.08))
# 隅の対角の壁側端: 全壁が等しければ旧式 h1 0.5 (1 + c) k elen/dr と同じビット、違えば 2 壁の目標の算術平均
h_, c_, el_, dr_ = 1.6e-5, 1.0352761804106076, 0.008 * np.sqrt(2), 0.008
check("(n6) diag_first: 全壁が等しいとき旧式と同じビット", diag_first(h_, h_, c_, 1.0, el_, dr_) == h_ * 0.5 * (1 + c_) * 1.0 * el_ / dr_)
s_ = diag_first(2e-6, 7e-6, c_, 1.0, el_, dr_); rs, rw = corner_ratios(s_, el_, dr_, 2e-6, 7e-6, c_)
check("(n7) 隅の比: 対角の第一節点の法線成分 = 2 壁の目標の算術平均 (sidewall_in 2・ramp 7 µm では比 2.31 / 0.64 で ±5 % の外)",
      abs(s_ * dr_ / el_ - 0.5 * (2e-6 + c_ * 7e-6)) < 1e-18 and rs > 1 + CORNER_TOL and rw < 1 - CORNER_TOL, f"sidewall_in {rs:.4f} ramp {rw:.4f}")
rs, rw = corner_ratios(diag_first(h_, h_, c_, 1.0, el_, dr_), el_, dr_, h_, h_, c_)
check("(n8) 全壁 16 µm (既定) の隅の比は ±5 % 以内 (ランプ傾き c = 1.035 で sidewall_in 1.018・ramp 0.983)", abs(rs - 1) <= CORNER_TOL and abs(rw - 1) <= CORNER_TOL, f"{rs:.4f} {rw:.4f}")
el = dict(e2=el_, e6=el_)
# 同一平面の継ぎ目の検出そのものを試す: vehicle を ramp に従わせる前の組 (ramp 16・vehicle 64 µm) を直接与える (wall_h1 は今は止める、n3)
hwd = {k: v * 0.1 for k, v in wall_h1(dict(P0, **PRESETS["design"])).items()}; hwd["vehicle"] = 6.4e-4 * 0.1
sd = dict(e2=diag_first(hwd["sidewall_in"], hwd["ramp"], c_, 1.0, el_, dr_), e6=diag_first(hwd["sidewall_in"], hwd["cowl_in"], 1.0038, 1.0, el_, dr_))
at, av = (0.5 * hwd[t] * (c_ + 1.0038) * P0["E5K"] for t in ("ramp", "vehicle"))
iss = station_conflicts(1.0, c_, 1.0038, sd, el, dr_, hwd, at, av, True)
check("(n9) design (ramp 16・vehicle 64 µm): x >= L_sw の同一平面の継ぎ目 ramp | vehicle を解けない組合せとして止める (隅は解ける)",
      [i["kind"] for i in iss] == ["coplanar"] and iss[0]["ratio_ramp"] > 1.1 and iss[0]["ratio_vehicle"] < 0.9, str([(i["kind"], round(i.get("ratio_ramp", 0), 3), round(i.get("ratio_vehicle", 0), 3)) for i in iss]))
check("(n10) 同じ値の組で SW の無い区間 (x < L_sw) は止めない (側壁が ramp と vehicle を隔てる)", station_conflicts(0.5, c_, 1.0038, sd, el, dr_, hwd, at, av, False) == [])
hwf = {k: v * 0.1 for k, v in wall_h1(dict(P0, **PRESETS["final"])).items()}
sf = dict(e2=diag_first(hwf["sidewall_in"], hwf["ramp"], c_, 1.0, el_, dr_), e6=diag_first(hwf["sidewall_in"], hwf["cowl_in"], 1.0038, 1.0, el_, dr_))
iss = station_conflicts(0.5, c_, 1.0038, sf, el, dr_, hwf, 0.5 * hwf["ramp"] * 2.04, 0.5 * hwf["vehicle"] * 2.04, False)
check("(n11) final (cowl_in 3・ramp 7・sidewall_in 2 µm): 隅 WU・WL の両方を止める", sorted(i["where"][:2] for i in iss) == ["WL", "WU"], str([i["where"] for i in iss]))
# 両端分布 (SW 列: z = W/2 側 a4 ≈ 19 µm、z = z_o 側 kv = 1.02 × 64 µm、長さ t_sw 0.5 mm)
a4_, kv_, Ls, gtw = 1.9364e-5, 1.02 * 6.4e-5, 5.0e-4, 1.19
rng_ = two_end_n_range(Ls, a4_, kv_, gtw); d = two_end_spacing(Ls, rng_[0], a4_, kv_, gtw); q = np.maximum(d[1:] / d[:-1], d[:-1] / d[1:]).max()
check("(n12) SW 列の両端分布: 端の間隔が a4・kv (相対 1e-12)、全長、隣接比 <= 1.19", abs(d[0] / a4_ - 1) < 1e-12 and abs(d[-1] / kv_ - 1) < 1e-12 and abs(d.sum() / Ls - 1) < 1e-14 and q <= gtw * (1 + 1e-12),
      f"区間数 {rng_} 最大比 {q:.5f} 最大間隔 {d.max() * 1e6:.2f} µm")
ok_rng = all(two_end_spacing(Ls, n_, a4_, kv_, gtw) is not None for n_ in range(rng_[0], rng_[1] + 1))
bad = 0
for n_ in (rng_[0] - 1, rng_[1] + 1):
    try:
        two_end_spacing(Ls, n_, a4_, kv_, gtw)
    except ValueError:
        bad += 1
check("(n13) two_end_n_range の範囲内は全部解け、範囲の外 (1 少ない / 1 多い) は ValueError", ok_rng and bad == 2, f"範囲 {rng_}")
d = two_end_spacing(0.09, 90, 1.66e-5, 6.6e-5, gtw); q = np.maximum(d[1:] / d[:-1], d[:-1] / d[1:]).max()
check("(n14) 長い辺 (側壁帯の y、90 区間) は中央を一定の M で頭打ちにして全長に合わせる (端 a・b、隣接比 <= 1.19)",
      abs(d[0] / 1.66e-5 - 1) < 1e-12 and abs(d[-1] / 6.6e-5 - 1) < 1e-12 and q <= gtw * (1 + 1e-12) and np.isclose(d.max(), d[45]), f"M {d.max() * 1e3:.3f} mm")
law = two_end_law(Ls, rng_[0] + 1, a4_, kv_, gtw)
check("(n15) two_end_law: 節点数 n、相対位置 0..1 の単調列 (apply_node_laws で直線辺に置く)", law[0] == rng_[0] + 1 and law[1] == "Nodes" and law[2][0] == 0.0 and abs(law[2][-1] - 1) < 1e-15
      and all(b > a for a, b in zip(law[2][:-1], law[2][1:])))
# first_layer_check の壁別目標: 2 壁 (y = 0 の ramp は第一層 H1、x = 0 の sidewall_in は 2 H1) の箱
xs_ = prog(2 * H1, 1.15, NL + 3, tail=2); ys_ = prog(H1, 1.15, NL + 3, tail=2)
xyz, hx, F_ = box(xs_, ys_, np.linspace(0, 20 * H1, 5))
bq = {"ramp": F_["ymin"], "sidewall_in": F_["xmin"]}
fl = first_layer_check(xyz, hx, bq, {"ramp": H1, "sidewall_in": 2 * H1}, None, {"ramp": NL, "sidewall_in": NL}, None, DR, expected_tags=("ramp", "sidewall_in"))
check("(n16) 壁別の目標 (ramp H1・sidewall_in 2 H1) で両壁 PASS、report に目標を記録", fl["ramp"]["ok"] and fl["sidewall_in"]["ok"] and fl["sidewall_in"]["target_h1"] == 2 * H1,
      f"ramp {fl['ramp']['ratio_min']:.3f}–{fl['ramp']['ratio_max']:.3f} sidewall_in {fl['sidewall_in']['ratio_min']:.3f}–{fl['sidewall_in']['ratio_max']:.3f}")
fl = first_layer_check(xyz, hx, bq, H1, H1, NL, NL, DR, expected_tags=("ramp", "sidewall_in"))
check("(n17) 同じ格子を共通の h1 で判定すると sidewall_in は FAIL (目標が壁別に効いている)", fl["ramp"]["ok"] and not fl["sidewall_in"]["ok"])
# 事前見積もり: 節点数の数え方 (count_mesh) を記録済みの接続模型 (B4a 板厚 0.005 H、scale 1: 4,349,372 節点・4,259,394 ヘキサ) で
cnt = dict(NZ=50, NR=62, NY=85, NSW=19, NTC=19, NL_co=30, NL_so=30, NFY=20, NFZ=20)
regs = [rg for rg, n_ in zip("AABBBCC", [11, 28, 20, 23, 25, 20, 26]) for _ in range(n_)]
check("(n18) count_mesh: 記録済み模型 (scale 1・板厚 0.005 H) の節点数・ヘキサ数を再現", count_mesh(regs, edge_counts(cnt)) == (4349372, 4259394), str(count_mesh(regs, edge_counts(cnt))))

print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAILED'}")
sys.exit(1 if FAIL else 0)
