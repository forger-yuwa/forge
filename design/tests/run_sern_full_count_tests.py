#!/usr/bin/env python3
"""⑤ SERN 全体の節点集計 (case/46.sern_design/cad/count_full_domain.py) の単体テスト
(plan tooling-sern-mesh-blocking §5.1 B4b-6。python3 で実行、gmsh 不要)。

(a) 汎用の数え上げ count_complex が、接続模型の範囲で生成器の count_mesh と一致する (記録済み 4 模型。全体の表を接続模型の行に限っても同じ)、Euler 標数 1
(b) 区間ごとの存在ブロックの規則 (U0/A 26・B 27・C 29・D 31、固体 = SW・CW1・CW2・V1・V2 の出入り) と全体の表の位相
    (全流体の断面で各辺を参照するブロックは 1 か 2、1 のものが外周 = 断面の境界ちょうど)
(c) 重複除去: 手で数えられる 2 ブロック × 3 区間の例 (共有辺・区間の境での和集合・一意辺の数・Euler 標数)
(d) x の station: 接続模型の手順 (junction_layout) が記録済みの区間数 (design / default) を再現。全体の手順 (full_layout) を合成形状で:
    物理 station を含む・単調・設計値の隣接比 <= 1.2・後流区間とフィレット区間の Δx 上限・マイタ区間 1 区間・区間名
(e) 形状 (FullGeom): x < 0 の平坦内壁と外壁のマイタ (頂点 x_m = t (1 - 1/cos θ)/tan θ、θ 5°・t 0.005 H で −0.000218305 H)、機体上面線 b(x) の端と下限
(f) 行 v の両端分布の区間数: 端 16 µm / 64 µm で長さ 0.02〜1.75 H は解けない (検出する)、上端 16 µm なら解ける
(g) メモリのモデル (600 万節点で変換 8.94 GiB・本体 7.78 GiB = codex 2026-10-07 の計算値) と仮想の遠方則 far_n_capped
同値類の節点数 (section_classes) の build との一致は Bump の較正に gmsh が要るので、ツールの --selfcheck で確かめる (mesh venv)。
"""
import math
import sys
import types
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "case/46.sern_design/cad"))
try:
    import gmsh  # noqa: F401
except ImportError:                       # 数え上げは gmsh を使わない。import だけ通す
    sys.modules["gmsh"] = types.ModuleType("gmsh")
import hex_junction_model as hjm  # noqa: E402
import count_full_domain as C  # noqa: E402

FAIL = 0


def check(name, cond, detail=""):
    global FAIL
    print(("ok   " if cond else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
    if not cond:
        FAIL += 1


# ---------------------------------------------------------------- (a) count_mesh との一致
for r in C.selfcheck_pure():
    check(f"(a1) {r['model']}: count_complex (接続模型の表) = count_mesh = 記録値", r["generic"] == r["recorded"] == r["count_mesh"], str(r["generic"]))
    check(f"(a2) {r['model']}: 全体の表を接続模型の行に限っても同じ、Euler 標数 1", r["full_table"] == r["recorded"] and r["euler"] == [1, 1], str(r["euler"]))

# ---------------------------------------------------------------- (b) 存在ブロックの規則と表の位相
edge, blk = C.full_tables(); present = C.full_present(blk)
check("(b1) 区間ごとの流体ブロック数 U0/A 26・B 27・C 29・D 31 (下書き §2)",
      [len(present[r]) for r in C.REGIONS] == [26, 26, 27, 29, 31], str([len(present[r]) for r in C.REGIONS]))
check("(b2) 固体: SW は B から、CW1・CW2 は C から、V1・V2 は D から流体 (上流では数えない)",
      "SW" not in present["A"] and "SW" in present["B"] and "CW1" not in present["B"] and {"CW1", "CW2"} <= present["C"]
      and "V1" not in present["C"] and {"V1", "V2"} <= present["D"] and present["U0"] == present["A"])
check("(b3) 全ブロックの 4 辺が辺の表にあり、V1 の下辺は N_TOP の上辺 e1、V2 の下辺は SW の上辺 hz41",
      all(e in edge for b in blk for e, _ in blk[b]) and blk["V1"][0][0] == "e1" and blk["N_TOP"][0][0] == "e1" and blk["V2"][0][0] == "hz41" and blk["SW"][2][0] == "hz41")
refs = {}
for b in present["D"]:
    for e, _ in blk[b]: refs[e] = refs.get(e, 0) + 1
outer = {e for e, n in refs.items() if n == 1}
exp_outer = ({"hz0%d" % k for k in range(4)} | {"hz7%d" % k for k in range(4)} | {"vy%d4" % j for j in range(7)}
             | {"vy00", "vy10", "vy20", "e4", "e11", "e9", "vy40", "vy50", "vy60"})
check("(b4) 全流体 (D) の断面: 各辺の参照は 1 か 2、参照 1 の辺 = 外周 (y_bot・上端・z = 0・z_far) ちょうど",
      set(refs.values()) <= {1, 2} and outer == exp_outer, str(sorted(outer ^ exp_outer)))
check("(b5) 帰属の順 BLOCK_ORDER が全 31 ブロックをちょうど 1 回含み、(行, 列) が決まる",
      sorted(C.BLOCK_ORDER) == sorted(blk) and len(set(C.BLOCK_ORDER)) == 31 and all(C.block_cell(b)[0] in C.ROWS for b in blk))
cnt_ = dict(NZ=50, NR=62, NY=109, NSW=12, NTC=19, NL_co=30, NL_so=22, NFY=289, NFZ=34, NV=73, NL_vt=22, NT=46)
ec = C.full_edge_counts(cnt_)
check("(b6) 新しい同値類: 行 v = NV、行 tb = NL_vt + 1、行 tc = NT + 1、y 線 5–7 の z は列ごとに NZ・NSW・NL_so + 1・NFZ + 1",
      ec["vy42"] == 73 and ec["vy50"] == 23 and ec["vy64"] == 47 and ec["hz70"] == 50 and ec["hz51"] == 12 and ec["hz62"] == 23 and ec["hz53"] == 35
      and ec["e1"] == 50 and ec["vy04"] == 290)
res_full = C.count_complex(["U0"] * 3 + ["A"] * 4 + ["B"] * 2 + ["C"] * 5 + ["D"] * 3, present, blk, edge, ec, order=C.BLOCK_ORDER)
check("(b7) 全体の表で Euler 標数 1 (固体は入口面・対称面に付くので流体は可縮)、帰属の和 = 節点数", res_full["euler"] == 1
      and sum(res_full["nodes_by_block"].values()) == res_full["nodes"] and sum(res_full["hexes_by_block"].values()) == res_full["hexes"], str(res_full["euler"]))

# ---------------------------------------------------------------- (c) 手で数えられる例
# 断面: B0 (y 0..1、z 3 節点 × y 4 節点) の上に B1 (y 1..2、z 3 × y 5) が辺 h1 を共有。区間 0 は B0 だけ、区間 1・2 は両方
E_ = {"h0": ("A", "B"), "h1": ("C", "D"), "h2": ("E", "F"), "v0": ("A", "C"), "v1": ("B", "D"), "v2": ("C", "E"), "v3": ("D", "F")}
B_ = {"B0": [("h0", 1), ("v1", 1), ("h1", -1), ("v0", -1)], "B1": [("h1", 1), ("v3", 1), ("h2", -1), ("v2", -1)]}
EC_ = {"h0": 3, "h1": 3, "h2": 3, "v0": 4, "v1": 4, "v2": 5, "v3": 5}
s2 = C.section_stats({"B0", "B1"}, B_, E_, EC_)
check("(c1) 断面: 3x4 + 3x5 - 共有 3 = 24 節点、辺 17 + 22 - 2 = 37、四角形 6 + 8", (s2["nodes"], s2["edges"], s2["quads"]) == (24, 37, 14), str((s2["nodes"], s2["edges"], s2["quads"])))
rc = C.count_complex(["X", "Y", "Y"], {"X": {"B0"}, "Y": {"B0", "B1"}}, B_, E_, EC_)
check("(c2) 3 区間: 節点 12 + 24 × 3 = 84 (区間の境の station は前後の和集合)、ヘキサ 6 + 14 × 2 = 34", (rc["nodes"], rc["hexes"]) == (84, 34), str((rc["nodes"], rc["hexes"])))
rb = C.count_complex(["X", "X"], {"X": {"B0"}}, B_, E_, EC_)
check("(c3) 単一ブロックの直方体 (3 x 3 x 4 節点) の一意辺 = 2·3·4 + 3·2·4 + 3·3·3 = 75、Euler 1",
      (rb["nodes"], rb["edges"], rb["euler"]) == (36, 75, 1), str((rb["nodes"], rb["edges"], rb["euler"])))
check("(c4) L 字の 2 ブロック × 3 区間も Euler 1、station の帰属は下流側の区間", rc["euler"] == 1 and rc["stations_by_region"] == {"X": 1, "Y": 3})

# ---------------------------------------------------------------- (d) x の station
for pre, nx_ref, ns_ref in (("design", [8, 28, 20, 23, 25, 20, 26], 151), (None, [11, 28, 20, 23, 25, 20, 26], 154)):
    P = dict(hjm.P0)
    if pre: P.update(hjm.PRESETS[pre])
    lay = C.junction_layout(P)
    check(f"(d1) junction_layout ({pre or 'default'}) = 生成器の記録 (区間数 {nx_ref}・station {ns_ref})",
          [n for n, _ in lay["xlaw"]] == nx_ref and len(lay["xs"]) == ns_ref, str([n for n, _ in lay["xlaw"]]))


def synth_contours(H=0.1, Lr=10.0, th_c=5.0, xf1=-0.0131652, xf2=0.0127167):
    """合成の輪郭 (m): ランプ = x_f1..0 で平坦 y = 1、0..L_r で 1 + 1.7 (1 - (1 - x/L_r)^2) (端の傾き 0)。カウル = 傾き -θ_c の直線 0..1.2"""
    xr = np.r_[np.linspace(xf1, 0.0, 20, endpoint=False), np.linspace(0.0, Lr, 2001)]
    yr = np.where(xr < 0, 1.0, 1.0 + 1.7 * (1.0 - (1.0 - xr / Lr) ** 2))
    xc = np.linspace(0.0, 1.2, 201); yc = -math.tan(math.radians(th_c)) * xc
    return dict(ramp=np.stack([xr, yr], 1) * H, cowl=np.stack([xc, yc], 1) * H), dict(x_f1=xf1, x_f2=xf2)


cont, meta = synth_contours()
P = C.make_P("design"); P["MEET_SCAN"] = 120
G = C.FullGeom(P, cont["ramp"], cont["cowl"], meta); H = G.H
lay = C.full_layout(P, G); xs = np.array(lay["xs"]) / H; dx = np.diff(xs)
phys = [-0.5, meta["x_f1"], G.x_m / H, 0.0, meta["x_f2"], 0.5, 0.8, 0.82, 1.0, 1.2, 1.22, 1.5, G.x0 / H, G.L_ramp / H - 0.3, G.L_ramp / H, G.L_ramp / H + 0.02,
        G.L_ramp / H + 0.3, G.XEND / H]
check("(d2) full_layout: 物理 station を全部含む (丸めなし)・単調増加", all(np.min(np.abs(xs - p)) < 1e-9 for p in phys) and np.all(dx > 0),
      str([p for p in phys if np.min(np.abs(xs - p)) >= 1e-9]))
check("(d3) 設計値の x 隣接比 <= 1.2 (出会う位置を間隔列全体で選ぶ)", lay["x_adjacent_ratio_design"] <= 1.2 * (1 + 1e-9), "%.4f at %.4f" % (lay["x_adjacent_ratio_design"], lay["x_adjacent_ratio_at"]))


def dx_in(a, b):
    m = (xs[:-1] >= a - 1e-12) & (xs[1:] <= b + 1e-12); return dx[m]


cap_f = (meta["x_f2"] - meta["x_f1"]) / P["FILLET_NMIN"]
check("(d4) 後流区間の Δx 上限 (側壁・カウル 0.001 H、ベース 0.004 H) とフィレット区間の Δx <= (x_f2 - x_f1)/12",
      dx_in(0.8, 0.82).max() <= 1e-3 * (1 + 1e-9) and dx_in(1.2, 1.22).max() <= 1e-3 * (1 + 1e-9)
      and dx_in(G.L_ramp / H, G.L_ramp / H + 0.02).max() <= 4e-3 * (1 + 1e-9) and dx_in(meta["x_f1"], meta["x_f2"]).max() <= cap_f * (1 + 1e-9)
      and len(dx_in(meta["x_f1"], meta["x_f2"])) >= 12)
check("(d5) マイタ区間 [x_m, 0] は 1 区間、端面の第一間隔 (ベース 17 µm・側壁/カウル端 0.1 mm)",
      len(dx_in(G.x_m / H, 0.0)) == 1 and abs(dx_in(G.L_ramp / H, G.L_ramp / H + 0.02)[0] - 1.7e-4) < 1e-12
      and abs(dx_in(0.8, 0.82)[0] - 1e-3) < 1e-12 and abs(dx_in(0.5, 0.8)[-1] - 1e-3) < 1e-12)
regs = lay["regions"]; xm_ = 0.5 * (xs[:-1] + xs[1:])
reg_ok = all((r == "U0") == (x < 0) and (r == "D") == (x > G.L_ramp / H) and (r == "B") == (0.8 < x < 1.2) for r, x in zip(regs, xm_))
check("(d6) station 区間の区間名: x < 0 = U0、(L_sw, L_cowl) = B、x > L_ramp = D", reg_ok)
P2 = dict(P); P2["MITER_STATION"] = 0
lay2 = C.full_layout(P2, G)
check("(d7) マイタ頂点を station にしないと station が減る (折れ点の局所細分のコスト)", lay2["n_stations"] < lay["n_stations"],
      "%d -> %d" % (lay["n_stations"], lay2["n_stations"]))

# ---------------------------------------------------------------- (e) 形状
th = math.radians(5.0); t = 0.005
xm_ref = t * (1 - 1 / math.cos(th)) / math.tan(th)
check("(e1) 外壁のマイタ頂点 x_m = t (1 - 1/cos θ)/tan θ (θ 5°・t 0.005 H で -0.000218305 H、§4.15 (4))",
      abs(G.x_m / H - xm_ref) < 1e-9 and abs(xm_ref + 0.000218305) < 1e-9, "%.9f" % (G.x_m / H))
xx = np.array([-0.3, -0.01, G.x_m / H - 1e-10, G.x_m / H + 1e-10, 0.5]) * H        # ±1e-10 H (傾き 5° の段差 1.7e-11 H)
check("(e2) x < 0: 内壁 y = 0・傾き 0、外壁は x_m まで y = -t、x_m で連続 (段差 < 1e-9 H)",
      np.all(G.yc(xx[:2]) == 0) and np.all(G.dyc(xx[:2]) == 0) and np.allclose(G.yo(xx[:3]), -t * H) and abs(G.yo(xx[3]) - G.yo(xx[2])) < 1e-9 * H)
check("(e3) 0 <= x では生成器の形状と同じ (内壁・外壁・ランプ)", abs(G.yc(xx[4]) - hjm.Geom.yc(G, xx[4])) == 0 and abs(G.yo(xx[4]) - hjm.Geom.yo(G, xx[4])) == 0
      and abs(G.yr(xx[4]) - hjm.Geom.yr(G, xx[4])) == 0)
xb = np.linspace(-0.5, G.XEND / H, 4001) * H
check("(e4) 機体上面線: x < x0 で y_veh、L_ramp で y_e + t_base、どこでも b >= y_r + t_base、下流は y_r + t_base (ランプは θ_e の直線)",
      abs(G.b(0.0) - G.y_veh) < 1e-15 and abs(G.b(G.L_ramp) - (G.ramp.yb + G.tb)) < 1e-12 and np.all(G.b(xb) >= G.yr(xb) + G.tb - 1e-15)
      and abs(G.b(G.L_ramp + 0.1 * H) - (G.ramp.yb + G.th_e * 0.1 * H + G.tb)) < 1e-12)
check("(e5) 機体上面線の傾き: テーパ始点 x0 で 0 (C¹)、L_ramp の上流側で min(θ_e, 0) - tan 3°",
      abs(G.db(G.x0 + 1e-9)) < 1e-6 and abs(G.db(G.L_ramp - 1e-12) - G.m1) < 1e-6)

# ---------------------------------------------------------------- (f) 行 v
xs_m = np.r_[np.linspace(-0.5, G.L_ramp / H, 300), G.L_ramp / H + 0.01, G.XEND / H] * H
svs = [dict(a_v=1.6e-4 * H * 1.03)] * len(xs_m)
top, tinfo = C.top_classes(P, G, list(xs_m), svs)
check("(f1) design (上端 64 µm・下端 16.5 µm、長さ 0.02〜1.75 H): 共通部分が空 = 解けない組合せとして返す (区間数は長い辺の下限)",
      tinfo["conflict"] is not None and tinfo["NV_range"][0] > tinfo["NV_range"][1] and top["NV"] == tinfo["conflict"]["n_lo"] + 1,
      str(tinfo["NV_range"]))
P3 = dict(P); P3["H1_vehicle_top"] = 1.6e-4
top3, tinfo3 = C.top_classes(P3, G, list(xs_m), svs)
check("(f2) 上端 16 µm なら解ける (下限 <= 上限)・行 tb は prog_n(DR, h_vt, 1.2) + 4 区間", tinfo3["conflict"] is None and tinfo3["NV_range"][0] <= tinfo3["NV_range"][1]
      and top3["NL_vt"] == hjm.prog_n(0.08 * H, 1.6e-4 * H, 1.2) + 4 and top["NL_vt"] == hjm.prog_n(0.08 * H, 6.4e-4 * H, 1.2) + 4, str(tinfo3["NV_range"]))
try:
    C.wall_h1_full(dict(P, H1_vehicle_side=1.6e-4)); bad = False
except ValueError:
    bad = True
check("(f3) vehicle_side と sidewall_out (同一平面 z = z_o) の第一層が違えば ValueError", bad)

# ---------------------------------------------------------------- (g) メモリ・仮想の遠方則
m6 = C.memory(6.0e6)
check("(g1) 600 万節点: 変換 8.94 GiB・本体 7.78 GiB (codex 2026-10-07 の計算値)", abs(m6["convert"] - 8.94) < 0.01 and abs(m6["solver_host"] - 7.78) < 0.01,
      "%.3f / %.3f" % (m6["convert"], m6["solver_host"]))
check("(g2) far_n_capped: L 10・第一間隔 1・上限 2・公比 2 -> 1 + ceil(9/2) = 6 区間", C.far_n_capped([(10.0, 1.0)], 2.0, 2.0)[0] == 6)
n_p, _ = hjm.far_n([(12.24, 0.0177)], 0.1, 1.15); n_c, _ = C.far_n_capped([(12.24, 0.0177)], 0.1, 1.15)
check("(g3) 行 a 相当 (12.24 H・第一間隔 0.018 H・上限 0.1 H): 生成器の等比 1 本は capped の約 2 倍の区間数", n_p > 1.8 * n_c, "%d vs %d" % (n_p, n_c))

print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAILED'}")
sys.exit(1 if FAIL else 0)
