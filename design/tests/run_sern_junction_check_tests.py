#!/usr/bin/env python3
"""⑤ SERN 接続模型 (case/46.sern_design/cad/hex_junction_model.py) の検査の単体テスト
(plan tooling-sern-mesh-blocking §5.1 B4-0 / B4-1。python3 で実行、gmsh 不要)。

B4-0: 壁第一層検査 `first_layer_check`。B4-1: 出力実座標の隣接間隔比 `adjacent_spacing_check`、後流 Δx `wake_dx_check`、
形状ゲート `shape_check` (と MOC 輪郭 `Contour` / カウル外壁のオフセット `Geom.yo`) — 末尾の (e)(f)(g)。
B4-2: 継ぎ目の分布則 (片側等比 `prog_nodes`・継ぎ目の間隔の連動 `match_spacing`・z の節点数の規則) と `cowl_side` の一般壁分類 — (h)(i)。

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
                                prog_nodes, prog_r, prog_n, match_spacing, END_TAGS, WALL_TAGS)

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

print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAILED'}")
sys.exit(1 if FAIL else 0)
