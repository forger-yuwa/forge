#!/usr/bin/env python3
"""SERN 全体 (全ヘキサ・ブロッキング) の節点集計 — gmsh で格子を張らずに、重複除去した節点・ヘキサ・一意辺の数とメモリを出す
(plan tooling-sern-mesh-blocking §5.1 B4b-6。§4.15・§4.16、codex diagnose 2026-10-07 Major 3「全体を生成する前に区間別の実分割数と
共有節点を正確に集計」)。

接続は下書き notes/investigations/2026-10-06-sern-b4a1/connection-table-draft.md (区間 U0/A/B/C/D、行 a/b/c/de/v/tb/tc × 列 0–3、
追加 12 ブロック) を §4.15 の採用で読む: 機体幅 V1 (外面 z = z_o)、側壁跡上の帯は vehicle (= ramp と同一平面、ramp の第一層に従う)、
下側 L1 (水平 y_bot、跡線は後縁接線の直線)、x < 0 のカウル内壁は平坦 + 外壁のマイタ (両方の折れ点を station に)、機体上面線は CAD 側で固定、
t_base 0.02 H、ランプフィレットは局所 x 細分 (旧の計 12 区間を下限)。

接続模型の生成器 hex_junction_model.py の関数 (march・prog_n・prog_r・far_n・Bump・match_spacing・diag_first・station_conflicts・
two_end_n_range・ring_diag_spacing・edge_counts ほか) をそのまま使う。build() の中の「断面の分布 (節点数は同値類ごとに 1 つ)」の規則は
関数の外から呼べないので、同じ演算を section_classes に写した (--selfcheck で build の ESTIMATE と同値類の節点数が一致することを確かめる)。

  x の station (full_layout): 物理 station の区間ごとに生成器の march (端面から離れる向き、後流区間の間隔上限、最終間隔 <= 1.15 hmax)。
      物理 station: -L_up、フィレット始点 x_f1、外壁のマイタ頂点 x_m、内壁の折れ点 0、フィレット終点 x_f2、L_sw - ZONE、L_sw、L_sw + LWAKE、
      (L_sw + L_cowl)/2、L_cowl、L_cowl + LWAKE、L_cowl + ZONE、テーパ始点 x0、L_ramp ± ZONE、L_ramp、L_ramp + LWAKE_B、x_out。
      [x_m, 0] は 1 区間 (21.8 µm)、そこから両向きに等間隔で始めてフィレット区間は Δx <= (x_f2 - x_f1)/FILLET_NMIN。
      細かい源が向かい合う区間 (フィレット | sidewall_end、cowl_base | vehicle_base) は両側から march し、出会う位置 (格子の station、
      物理入力ではない) を両側の最終間隔の比が最小になる所に置く。B は生成器と同じ中点 (L_sw + L_cowl)/2。
  断面 (行 y × 列 z の同値類): 接続模型の同値類は section_classes (build と同じ規則を全体の station 列で)、新しい行は top_classes
      (行 tb = 機体上面の壁帯 [prog_n(DR, h_vt, 1.2) + 4 区間]、行 tc = 遠方帯 [far_n]、行 v = 両端分布 [下端 = 行 de の上端、上端 = 行 tb の第一間隔]
      の区間数の共通部分)。列 0–3 の z 分布は行 v/tb/tc にもそのまま伝わる (内部面)。
  節点 (count_complex): 区間ごとに存在するブロック (固体は数えない) の面の和集合を station ごとに数える (端面 station は前後の区間の和集合)。
      節点 = Σ_station (点 + 辺の内部 + 面の内部)、ヘキサ = Σ_区間 Σ_ブロック (ni-1)(nj-1)、一意辺 = Σ_station 断面の辺 + Σ_区間 x 方向の辺。
      面と Euler 標数 (V - E + F - C = 1) も数えて自己検査にする。接続模型の表に限れば hex_junction_model.count_mesh と一致する。

usage (mesh venv。Bump の較正だけ gmsh の 1D を使う。全体の格子は張らない):
  .venv-mesh/bin/python case/46.sern_design/cad/count_full_domain.py --contours DIR [--preset design|final] [--set k=v ...]
      [--sens] [--final-ref] [--selfcheck] [--out OUT.json]
  DIR は export_contours.py の出力 (ramp_contour.csv / cowl_contour.csv / contours_meta.json)
"""
import sys, math, json, argparse
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import hex_junction_model as hjm  # noqa: E402

gmsh = hjm.gmsh

# ---------------------------------------------------------------- 全体の入力 (/H)。生成器の P0 + プリセットに足す
# 領域寸法は生産 YAML (problem_3d_prod_m6on_wallres_lswx08.yaml) と旧メッシャ mesh_sern3d の値。y_bot は L1 (§4.15 (3))、z_far は絶対 2.50 H (§4.14-9)
FULL_P0 = dict(YBOT=-12.320694, ZFAR=2.5,
               L_UP=0.5, X_OUT_EXTRA=2.0, TOP_DEPTH=2.0,          # 上流長・出口 = L_ramp + 2 H・機体上の外部流の厚さ (旧 mesh.L_up / x_out_extra / top_depth)
               CLEAR=0.06, T_BASE=0.02, TAPER=0.35, WEDGE_DEG=3.0,  # 機体上面線 (旧 vehicle_clearance / t_base / vehicle_taper / vehicle_wedge_deg)
               LWAKE_B=0.02, DXW_B=0.004,                          # 機体ベースの後流: 区間長 (仮定: カウルの LWAKE と同じ)・Δx 上限 = t_base/5
               FILLET_NMIN=12,                                     # ランプフィレット [x_f1, x_f2] の区間数の下限 (旧の計 12、§4.15 (8))
               MITER_STATION=1,                                    # 外壁のマイタ頂点 x_m を station にする (§4.15 (4))。0 は感度用
               RAMP_ZONE=1,                                        # L_ramp ± ZONE の物理 station (端面の上下流に ZONE、生成器の L_sw - ZONE・L_cowl + ZONE と同じ扱い)
               FAR_LAW="progression",                              # 遠方帯の分布則: progression = 生成器 (等比 1 本)、capped = 仮想の「公比 g_t で HFAR まで増やして一様」
               NY_CAP_SW_ONLY=0,                                   # 1 = NY の AR 上限 (AR_TAN × 側壁の h1) を側壁のある x <= L_sw の station だけに掛ける (感度用。
                                                                   #     生成器は全 station の最大の高さに掛ける。側壁は L_sw で終わるので下流の高い断面には壁層セルが無い)
               MEET_SCAN=400)                                      # 出会う位置の探索の分割数
# 全体だけにある壁 (生成器の wall_h1 は未知の H1_* を拒むので、ここで別に持つ)
FULL_WALLS = ("vehicle_top", "vehicle_side", "vehicle_base")
FULL_KEYS = (set(FULL_P0) - {"YBOT", "ZFAR"}) | {"H1_" + t for t in FULL_WALLS} | {"HFAR_Y", "HFAR_Z", "HFAR_T"}   # YBOT・ZFAR は Geom が読む
# plan §4.16 の暫定の試験値。design: vehicle_top 64 µm・vehicle_side (= sidewall_out = cowl_side、同一平面) 64 µm・vehicle_base 17 µm
FULL_PRESETS = {
    "design": dict(hjm.PRESETS["design"], H1_vehicle_top=6.4e-4, H1_vehicle_side=6.4e-4, H1_vehicle_base=1.7e-4),
    "final": dict(hjm.PRESETS["final"], H1_vehicle_top=6.4e-4, H1_vehicle_side=6.4e-4, H1_vehicle_base=1.7e-4),
}
# 参考: 隅に接する 3 壁 (ramp・cowl_in・sidewall_in) を共通 2 µm にそろえた final (§4.16: 隅の比が ±5 % の外だと今の構成で作れないため)
FINAL_COMMON2 = dict(FULL_PRESETS["final"], H1_ramp=2.0e-5, H1_cowl_in=2.0e-5, H1_sidewall_in=2.0e-5)

# メモリのモデル (B/節点 と切片)。生成+検査は接続模型の実測 (約 1.45 kB/節点)、変換器 1579 B/節点 + 117 MiB (B4-5)、
# 本体 1371 B/節点 + 120 MiB (plan architecture-solver-host-memory R1–R3)、GPU 約 1.4 kB/節点 (B4b-3)
MEM_MODEL = dict(generate_check=(1450.0, 0.0), convert=(1579.0, 117.0), solver_host=(1371.0, 120.0), solver_gpu=(1400.0, 0.0))
BUDGET_NODES = 6.0e6
GIB = 1024.0 ** 3


def memory(nodes):
    """節点数 -> 工程別のメモリ (GiB)。モデルは MEM_MODEL (見積もりで実行保証ではない)"""
    return {k: (a * nodes + b * 1024.0 ** 2) / GIB for k, (a, b) in MEM_MODEL.items()}


def split_P(P):
    """全体の P から、生成器 (wall_h1 ほか) に渡す接続模型の P を作る (全体だけのキーを落とす)"""
    return {k: v for k, v in P.items() if k not in FULL_KEYS}


def wall_h1_full(P):
    """壁タグ -> 第一内部点距離 (/H)。接続模型の 9 タグは生成器の wall_h1、全体の 3 壁を足す。
    vehicle_side は sidewall_out と同一平面 (z = z_o) で列 2 の z 分布 (hz?2) を共有するので同じ値 (違えば ValueError、cowl_side と同じ規則)"""
    out = hjm.wall_h1(split_P(P))
    out["vehicle_top"] = float(P.get("H1_vehicle_top", P["H1"]))
    out["vehicle_base"] = float(P.get("H1_vehicle_base", P.get("H1_END", P["H1"])))
    vs = float(P.get("H1_vehicle_side", out["sidewall_out"]))
    if not math.isclose(vs, out["sidewall_out"], rel_tol=1e-12, abs_tol=0.0):
        raise ValueError("vehicle_side (%.6g H) と sidewall_out (%.6g H) の第一層が違う: 同一平面 z = z_o で z 分布 hz?2 を共有する" % (vs, out["sidewall_out"]))
    out["vehicle_side"] = vs
    return out


# ---------------------------------------------------------------- 全体の接続表 (下書き §3: 行 j × 列 k、N はバタフライ)
ROWS = ("a", "b", "c", "de", "v", "tb", "tc")                 # 行 j = 0..6 (y 線 j と j+1 の間)
# y 線: 0 YBOT / 1 y_o - DR / 2 y_o / 3 y_c / 4 y_r / 5 b / 6 b + DR / 7 b + top_depth。z 線: 0 / z_w / z_o / z_o + DR / z_far
NEW_BLOCKS = {"V1": (4, 0), "V2": (4, 1), "Sb_v": (4, 2), "Sc_v": (4, 3),
              "T1b": (5, 0), "T2b": (5, 1), "Sb_tb": (5, 2), "Sc_tb": (5, 3),
              "T1c": (6, 0), "T2c": (6, 1), "Sb_tc": (6, 2), "Sc_tc": (6, 3)}
SOLID = {"U0": {"CW1", "CW2", "SW", "V1", "V2"}, "A": {"CW1", "CW2", "SW", "V1", "V2"}, "B": {"CW1", "CW2", "V1", "V2"},
         "C": {"V1", "V2"}, "D": set()}
REGIONS = ("U0", "A", "B", "C", "D")


def _gp(j, k):
    return hjm.GNAME.get((j, k), "G%d%d" % (j, k))


def full_tables():
    """全体の (辺 -> 端点, ブロック -> 4 辺の (辺名, 向き))。接続模型の表 (hex_junction_model の EDGE・BLK) をそのまま含み、
    y 線 5–7 の水平辺 hz5k–hz7k、行 4–6 の鉛直辺 vy4k–vy6k、追加 12 ブロックを足す。V1 の下辺 hz40 は N_TOP の上辺 e1 (PT0 -> WU) の別名"""
    edge = dict(hjm.EDGE)
    alias = dict(hjm.ALIAS); alias["hz40"] = ("e1", 1)
    for j in range(5, 8):
        for k in range(4): edge["hz%d%d" % (j, k)] = (_gp(j, k), _gp(j, k + 1))
    for j in range(4, 7):
        for k in range(5): edge["vy%d%d" % (j, k)] = (_gp(j, k), _gp(j + 1, k))

    def ref(name, sign=1):
        if name.startswith("-"): return ref(name[1:], -sign)
        if name in alias: b, s_ = alias[name]; return b, s_ * sign
        return name, sign
    blk = dict(hjm.BLK)
    for b, (j, k) in NEW_BLOCKS.items(): blk[b] = [ref(e) for e in hjm.grid_block(j, k)]
    return edge, blk


def full_present(blk):
    """区間 -> 流体ブロックの集合 (固体は数えない)。U0/A 26・B 27 (+SW)・C 29 (+CW1, CW2)・D 31 (+V1, V2)"""
    return {r: set(blk) - SOLID[r] for r in REGIONS}


def block_cell(b):
    """ブロック -> (行, 列)。内訳の集計用"""
    if b.startswith("N_"): return "de", 0
    fixed = dict(SW=("de", 1), CW1=("c", 0), CW2=("c", 1), U1a=("a", 0), U1b=("b", 0), U2a=("a", 1), U2b=("b", 1),
                 V1=("v", 0), V2=("v", 1), T1b=("tb", 0), T2b=("tb", 1), T1c=("tc", 0), T2c=("tc", 1))
    if b in fixed: return fixed[b]
    if b[0] == "S": return b.split("_")[1], (2 if b[1] == "b" else 3)
    raise KeyError(b)


# 節点の帰属の順 (下の行から、行の中は列 0 -> 3。共有する点・辺は最初に現れる存在ブロックへ)
BLOCK_ORDER = ["U1a", "U2a", "Sb_a", "Sc_a", "U1b", "U2b", "Sb_b", "Sc_b", "CW1", "CW2", "Sb_c", "Sc_c",
               "N_BOT", "N_CORE", "N_SIDE", "N_TOP", "SW", "Sb_de", "Sc_de", "V1", "V2", "Sb_v", "Sc_v",
               "T1b", "T2b", "Sb_tb", "Sc_tb", "T1c", "T2c", "Sb_tc", "Sc_tc"]


def full_edge_counts(cnt):
    """辺 -> 節点数。接続模型の同値類は生成器の edge_counts、全体の新しい同値類:
    列 0–3 の z (NZ・NSW・NL_so + 1・NFZ + 1) は y 線 5–7 にもそのまま、行 v は NV (節点数)、行 tb は NL_vt + 1、行 tc は NT + 1"""
    ec = hjm.edge_counts(cnt)
    for j in (5, 6, 7):
        ec["hz%d0" % j] = cnt["NZ"]; ec["hz%d1" % j] = cnt["NSW"]; ec["hz%d2" % j] = cnt["NL_so"] + 1; ec["hz%d3" % j] = cnt["NFZ"] + 1
    for k in range(5):
        ec["vy4%d" % k] = cnt["NV"]; ec["vy5%d" % k] = cnt["NL_vt"] + 1; ec["vy6%d" % k] = cnt["NT"] + 1
    return ec


# ---------------------------------------------------------------- 重複除去した数え上げ (gmsh 非依存の純関数)
def section_stats(nb, blk, edge, ec, order=None):
    """断面 1 枚 (ブロック集合 nb の面の和集合) の節点・辺・四角形の数と、節点の帰属 (ブロック -> 数)。
    節点 = 点 (重複除去) + 辺の内部 (辺は重複除去) + 面の内部。辺 = 辺ごとの (n - 1) + 面の内部の格子線。四角形 = Σ (ni - 1)(nj - 1)。
    共有する点・辺は order で最初に現れる存在ブロックに帰属させる (合計は帰属の決め方によらない)"""
    blocks = [b for b in (order or sorted(nb)) if b in nb]
    if len(blocks) != len(nb): raise ValueError("order に無いブロック: %s" % sorted(set(nb) - set(blocks)))
    eds, pts = {}, {}
    for b in blocks:
        for e, _ in blk[b]:
            eds.setdefault(e, b)
            for p in edge[e]: pts.setdefault(p, b)
    own = {}; n_edges = n_quads = 0
    for b in blocks:
        ni, nj = ec[blk[b][0][0]], ec[blk[b][1][0]]
        own[b] = (ni - 2) * (nj - 2)
        n_edges += (ni - 1) * (nj - 2) + (nj - 1) * (ni - 2)
        n_quads += (ni - 1) * (nj - 1)
    for e, b in eds.items(): own[b] += ec[e] - 2; n_edges += ec[e] - 1
    for p, b in pts.items(): own[b] += 1
    return dict(nodes=sum(own.values()), edges=n_edges, quads=n_quads, own=own)


def count_complex(regions, present, blk, edge, ec, order=None):
    """区間ごとの存在ブロックから、重複除去した節点・ヘキサ・一意辺・四角形と Euler 標数を数える (gmsh 非依存の純関数)。
    regions は station 区間ごとの区間名の列 (station 数 - 1 個)、present は 区間名 -> ブロック集合。
    station s の断面は前後の区間の存在ブロックの和集合 (端面 station では片側にしか無いブロックの面が端面になる)。
    3D の数: 節点 V = Σ_s N2(nb_s)、辺 E = Σ_s E2(nb_s) + Σ_i N2(P_i) (x 方向の辺)、四角形 F = Σ_s F2(nb_s) + Σ_i E2(P_i)、ヘキサ C = Σ_i F2(P_i)。
    流体領域は可縮 (固体は入口面・対称面に付いている) なので V - E + F - C = 1 になる (数え方の自己検査)"""
    ns = len(regions) + 1; cache = {}

    def st(nb):
        k = frozenset(nb)
        if k not in cache: cache[k] = section_stats(nb, blk, edge, ec, order)
        return cache[k]
    V = E = F = C = 0
    own_tot, hex_blk = {}, {}
    nodes_rg = {r: 0 for r in dict.fromkeys(regions)}; hex_rg = dict.fromkeys(nodes_rg, 0); st_rg = dict.fromkeys(nodes_rg, 0)
    for s in range(ns):
        nb = (present[regions[s - 1]] if s > 0 else set()) | (present[regions[s]] if s < ns - 1 else set())
        q = st(nb); V += q["nodes"]; E += q["edges"]; F += q["quads"]
        rg = regions[min(s, ns - 2)]; nodes_rg[rg] += q["nodes"]; st_rg[rg] += 1          # station s は下流側の区間に数える (最後の station は最後の区間)
        for b, v in q["own"].items(): own_tot[b] = own_tot.get(b, 0) + v
    for rg in regions:
        q = st(present[rg]); E += q["nodes"]; F += q["edges"]; C += q["quads"]; hex_rg[rg] += q["quads"]
        for b in present[rg]:
            ni, nj = ec[blk[b][0][0]], ec[blk[b][1][0]]; hex_blk[b] = hex_blk.get(b, 0) + (ni - 1) * (nj - 1)
    return dict(nodes=V, hexes=C, edges=E, quads=F, euler=V - E + F - C, nodes_by_block=own_tot, hexes_by_block=hex_blk,
                nodes_by_region=nodes_rg, hexes_by_region=hex_rg, stations_by_region=st_rg, intervals=len(regions))


def breakdown_rows_cols(res):
    """節点・ヘキサの (行, 列) 別の表 (count_complex の帰属から)"""
    out = {}
    for key, src in (("nodes", res["nodes_by_block"]), ("hexes", res["hexes_by_block"])):
        tab = {r: [0, 0, 0, 0] for r in ROWS}
        for b, v in src.items():
            r, c = block_cell(b); tab[r][c] += v
        out[key] = tab
    return out


# ---------------------------------------------------------------- 遠方帯の区間数
def far_n_capped(pairs, hmax, g_t):
    """仮想の分布則 (感度用。生成器には無い): 第一間隔 h0 から公比 g_t で hmax まで増やし、その先は一様 (<= hmax)。
    pairs の全組で全長を覆う最小の区間数 (の最大) と、その最大間隔の上限 hmax を返す"""
    nmax = 0
    for L, h0 in pairs:
        k, S, d = 0, 0.0, h0
        while d < hmax and S + d < L: S += d; d *= g_t; k += 1
        nmax = max(nmax, k + max(int(math.ceil((L - S) / hmax - 1e-12)), 0 if S >= L else 1))
    return nmax, hmax


FAR_FN = {"progression": hjm.far_n, "capped": far_n_capped}


# ---------------------------------------------------------------- 形状
def _scal(x, v):
    return float(v) if np.ndim(x) == 0 else v


class FullGeom(hjm.Geom):
    """SERN 全体の形状 (m 単位)。生成器の Geom (MOC 輪郭の B-spline、カウル外壁のオフセット、断面の点) を引き継ぎ、全体で違う所だけ上書きする:
      - ランプ: x < x_f1 (フィレット始点) は平坦 (輪郭の始点の y)、L_ramp より下流は θ_e (輪郭の最終区間の傾き) の直線 (下書き §5 E1)
      - カウル内壁: x < 0 は平坦 y = 0、外壁は x <= x_m で平坦 y = -t、x_m (平坦の外壁と傾いた外壁の交点 = マイタ頂点) より下流は生成器のオフセット
        (§4.15 (4)。x_m は外壁のオフセット曲線から二分法で解く。カウルの跡は生成器の後縁接線の延長のまま = L1)
      - 機体上面線 b(x) (CAD 固定、§4.15 (6)): y_veh = 輪郭の最大 y + CLEAR、テーパ始点 x0 = L_ramp (1 - TAPER)、
        エルミートで y_e + t_base に着地 (端の傾き min(θ_e, 0) - tan WEDGE)、どこでも b >= y_r + t_base、L_ramp より下流は y_r + t_base
      - 領域: YBOT (L1)、ZFAR (絶対)、x_out = L_ramp + X_OUT_EXTRA。meta は contours_meta.json (x_f1・x_f2)"""
    def __init__(s, P, ramp_xy, cowl_xy, meta):
        H = P["H"]; xe = float(ramp_xy[-1, 0])
        Pj = split_P(P); Pj["XEND"] = xe / H
        while Pj["XEND"] * H > xe: Pj["XEND"] = float(np.nextafter(Pj["XEND"], 0.0))
        super().__init__(Pj, ramp_xy, cowl_xy)
        s.L_ramp = xe; s.XEND = xe + P["X_OUT_EXTRA"] * H
        s.th_e = float((ramp_xy[-1, 1] - ramp_xy[-2, 1]) / (ramp_xy[-1, 0] - ramp_xy[-2, 0]))
        s.x_f1, s.x_f2 = float(meta["x_f1"]) * H, float(meta["x_f2"]) * H
        if abs(s.x_f1 - s.ramp.xa) > 1e-9 * H: raise ValueError("ランプ輪郭の始点 %.9g m がフィレット始点 x_f1 %.9g m と一致しない" % (s.ramp.xa, s.x_f1))
        lo, hi = -2.0 * s.TC, 0.0                                       # マイタ頂点: 傾いた外壁 (生成器の yo) が y = -t を横切る x
        f = lambda x_: float(hjm.Geom.yo(s, x_)) + s.TC
        if not (f(lo) > 0.0 > f(hi)): raise ValueError("外壁のマイタ頂点が [-2t, 0] に無い")
        for _ in range(200):
            m = 0.5 * (lo + hi); lo, hi = (m, hi) if f(m) > 0.0 else (lo, m)
        s.x_m = 0.5 * (lo + hi)
        s.y_veh = float(ramp_xy[:, 1].max()) + P["CLEAR"] * H; s.tb = P["T_BASE"] * H
        s.taper_len = P["TAPER"] * s.L_ramp; s.x0 = s.L_ramp - s.taper_len
        s.y_end = s.ramp.yb + s.tb; s.m1 = min(s.th_e, 0.0) - math.tan(math.radians(P["WEDGE_DEG"]))
        s.TOPD = P["TOP_DEPTH"] * H

    def yr(s, x):
        xv = np.asarray(x, float); v = s.ramp.y(np.clip(xv, s.ramp.xa, s.L_ramp))
        return _scal(x, np.where(xv > s.L_ramp, s.ramp.yb + s.th_e * (xv - s.L_ramp), v))

    def dyr(s, x):
        xv = np.asarray(x, float); v = s.ramp.dy(np.clip(xv, s.ramp.xa, s.L_ramp))
        return _scal(x, np.where(xv > s.L_ramp, s.th_e, np.where(xv < s.ramp.xa, 0.0, v)))

    def yc(s, x):
        xv = np.asarray(x, float); return _scal(x, np.where(xv < 0.0, 0.0, s.cowl.y(xv)))

    def dyc(s, x):
        xv = np.asarray(x, float); return _scal(x, np.where(xv < 0.0, 0.0, s.cowl.dy(xv)))

    def yo(s, x):
        xv = np.asarray(x, float); return _scal(x, np.where(xv <= s.x_m, -s.TC, hjm.Geom.yo(s, np.maximum(xv, s.x_m))))

    def dyo(s, x):
        xv = np.asarray(x, float); return _scal(x, np.where(xv <= s.x_m, 0.0, hjm.Geom.dyo(s, np.maximum(xv, s.x_m))))

    def _herm(s, xv):
        u = np.clip((xv - s.x0) / s.taper_len, 0.0, 1.0); tl = s.taper_len
        y = (2 * u ** 3 - 3 * u ** 2 + 1) * s.y_veh + (-2 * u ** 3 + 3 * u ** 2) * s.y_end + (u ** 3 - u ** 2) * tl * s.m1
        dy = ((6 * u ** 2 - 6 * u) * s.y_veh + (-6 * u ** 2 + 6 * u) * s.y_end) / tl + (3 * u ** 2 - 2 * u) * s.m1
        return np.where(xv < s.x0, s.y_veh, y), np.where(xv < s.x0, 0.0, dy)

    def b(s, x):
        """機体上面線 (行 v の上端・行 tb の下端)"""
        xv = np.asarray(x, float); h, _ = s._herm(xv); lo = s.yr(xv) + s.tb
        return _scal(x, np.where(xv >= s.L_ramp, lo, np.maximum(h, lo)))

    def db(s, x):
        xv = np.asarray(x, float); h, dh = s._herm(xv); lo = s.yr(xv) + s.tb
        return _scal(x, np.where((xv >= s.L_ramp) | (lo > h), s.dyr(xv), dh))


# ---------------------------------------------------------------- x の station
def _end_spacing(L, n, r):
    """march の 1 区間 (区間数 n・公比 r、march の向き) の最終間隔"""
    w = abs(r) ** np.arange(n); return L * float(w[-1] / w.sum())


def run_march(segs, h0, hmax, caps, P):
    """源から離れる向きに並べた区間 segs ((始, 終) の列、向きは問わない) を生成器の march で分割する。戻り値: ((n, r) の列, 最終間隔)"""
    g = P["G"]; g_t = 1.0 + 0.75 * (g - 1.0); Ls = [abs(b_ - a_) for a_, b_ in segs]
    res = hjm.march(Ls, h0, g_t, g, hmax, caps=caps)
    return res, _end_spacing(Ls[-1], *res[-1])


def stations_from_laws(xp, xlaw):
    """物理 station xp と区間ごとの (n, r) (r < 0 は終端側が細かい) から格子の全 x 位置を作る (生成器の build と同じ式)"""
    xs = [xp[0]]
    for i in range(len(xp) - 1):
        xa, xb = xp[i], xp[i + 1]; n, r = xlaw[i]; q = abs(r) if r > 0 else 1.0 / abs(r)
        w = q ** np.arange(n); cw = np.cumsum(w) / w.sum()
        xs += [xa + (xb - xa) * float(c_) for c_ in cw[:-1]] + [xb]
    return xs


def junction_layout(P):
    """接続模型の station 列 (生成器 build の「物理 station」「x 分布」「格子の全 x 位置」と同じ手順。断面隅フィレット 0 = PROFILE 1 点のみ)。
    戻り値: dict(xp, ivp_regions, xlaw, xs, regions)。regions は station 区間ごとの A/B/C"""
    if len(hjm.PROFILE) != 1: raise NotImplementedError("断面隅フィレット (PROFILE 2 点以上) の station は未対応")
    H = P["H"]; LSW, LCOWL, XEND, ZONE, LWAKE = (P[k] * H for k in ("LSW", "LCOWL", "XEND", "ZONE", "LWAKE"))
    hw = {t: v * H for t, v in hjm.wall_h1(split_P(P)).items()}
    xm = 0.5 * (LSW + LCOWL)
    xs = {p[0] * H for p in hjm.PROFILE} | {LSW, LCOWL, XEND, xm, LSW - ZONE, LCOWL + ZONE, LSW + LWAKE, LCOWL + LWAKE}
    xp = sorted(x for x in xs if -1e-12 <= x <= XEND + 1e-12)
    xp = [x for i, x in enumerate(xp) if i == 0 or x - xp[i - 1] > 1e-9]
    reg = lambda xa, xb: "A" if xb <= LSW + 1e-12 else ("B" if xb <= LCOWL + 1e-12 else "C")
    ivp = [(xp[i], xp[i + 1], reg(xp[i], xp[i + 1])) for i in range(len(xp) - 1)]
    hx = P["HX"] * H; g = P["G"]; g_t = 1.0 + 0.75 * (g - 1.0); dxw = P["DXW"] * H; h1e_sw, h1e_cb = hw["sidewall_end"], hw["cowl_base"]
    near = lambda a, b: abs(a - b) <= 1e-9 * H
    capf = lambda i: dxw if any(near(ivp[i][0], L_) and near(ivp[i][1], L_ + LWAKE) for L_ in (LSW, LCOWL)) else None
    idx = {r_: [i for i, v in enumerate(ivp) if v[2] == r_] for r_ in "ABC"}
    xlaw = {}

    def run(ids, toward_plus, hmax, h0):
        res = hjm.march([ivp[i][1] - ivp[i][0] for i in ids], h0, g_t, g, hmax, caps=[capf(i) for i in ids])
        for i, (n, r) in zip(ids, res): xlaw[i] = (n, r if toward_plus else -r)
    run(idx["A"][::-1], False, hx, h1e_sw)
    run([i for i in idx["B"] if ivp[i][1] <= xm + 1e-12], True, hx, h1e_sw); run([i for i in idx["B"] if ivp[i][0] >= xm - 1e-12][::-1], False, hx, h1e_cb)
    run(idx["C"], True, P["HX_FAR"] * H, h1e_cb)
    xs = stations_from_laws(xp, xlaw)
    regions = [reg(xs[i], xs[i + 1]) for i in range(len(xs) - 1)]
    return dict(xp=xp, ivp_regions=[v[2] for v in ivp], xlaw=[xlaw[i] for i in range(len(ivp))], xs=xs, regions=regions)


def _spacings(segs, res):
    """march の結果 ((n, r) の列) を、源から離れる向きの間隔列にする"""
    out = []
    for (a_, b_), (n, r) in zip(segs, res):
        w = abs(r) ** np.arange(n); out += list(abs(b_ - a_) * w / w.sum())
    return out


def _meet(left, right, xl, xr, inner, P, scan):
    """両側に細かい源がある区間 [xl, xr] の出会う位置 m を決める。left/right = (固定区間の (n, r) を作る march の入力: segs, h0, hmax, caps)。
    左の run は left["segs"] + [xl, m] (inner の物理 station で分割)、右の run は right["segs"] + [m, xr] を右から。
    m を (xl, xr) の scan 点で試し、左の源から右の源までの間隔列 (左の run + 右の run の逆順) の隣接比の最大が g 以内で区間数が最小の位置を選ぶ
    (g 以内が無ければ隣接比の最大が最小の位置)。march は短い区間で公比が合わないと一様に落とす (生成器の規則) ので、出会う点だけでなく列全体で見る"""
    g = P["G"]; best = None
    for k in range(1, scan + 1):
        m = xl + (xr - xl) * k / (scan + 1.0)
        if any(abs(m - p) < 1e-6 * P["H"] for p in inner): continue
        lp = [xl] + [p for p in inner if p < m] + [m]; rp = [m] + [p for p in inner if p > m] + [xr]
        lsegs = left["segs"] + [(lp[i], lp[i + 1]) for i in range(len(lp) - 1)]
        rsegs = right["segs"] + [(rp[i + 1], rp[i]) for i in range(len(rp) - 1)][::-1]
        try:
            lres, le = run_march(lsegs, left["h0"], left["hmax"], left["caps"] + [None] * (len(lsegs) - len(left["caps"])), P)
            rres, re_ = run_march(rsegs, right["h0"], right["hmax"], right["caps"] + [None] * (len(rsegs) - len(right["caps"])), P)
        except ValueError:
            continue
        sp = np.array(_spacings(lsegs, lres) + _spacings(rsegs, rres)[::-1])
        rmax = float(np.maximum(sp[1:] / sp[:-1], sp[:-1] / sp[1:]).max()); n_tot = len(sp)
        ok = rmax <= g * (1 + 1e-9)
        key = (0 if ok else 1, n_tot if ok else rmax, rmax)
        if best is None or key < best[0]: best = (key, m, lsegs, lres, rsegs, rres, max(le / re_, re_ / le), rmax)
    if best is None: raise ValueError("出会う位置が見つからない [%.6g, %.6g]" % (xl, xr))
    return dict(m=best[1], left=list(zip(best[2], best[3])), right=list(zip(best[4], best[5])), ratio=best[6], ratio_run_max=best[7])


def full_layout(P, G):
    """全体の station 列。物理 station・区間名・区間ごとの (n, r)、出会う位置、設計値の x 隣接比を返す (モジュール冒頭の説明)"""
    H = P["H"]; ZONE, LWAKE = P["ZONE"] * H, P["LWAKE"] * H
    hw = {t: v * H for t, v in wall_h1_full(P).items()}
    hx, hxf = P["HX"] * H, P["HX_FAR"] * H; dxw, dxwb, lwb = P["DXW"] * H, P["DXW_B"] * H, P["LWAKE_B"] * H
    LUP, LSW, LCOWL, LR, XOUT = P["L_UP"] * H, G.LSW, G.LCOWL, G.L_ramp, G.XEND
    xf1, xf2, xm_, x0 = G.x_f1, G.x_f2, G.x_m, G.x0
    use_miter = bool(P.get("MITER_STATION", 1)); ramp_zone = bool(P.get("RAMP_ZONE", 1))
    cap_f = (xf2 - xf1) / float(P["FILLET_NMIN"])
    kink = (xm_, 0.0) if use_miter else None
    h0k = (0.0 - xm_) if use_miter else cap_f                       # 折れ点の隣の第一間隔 (マイタ区間と等間隔で始める。マイタ無しはフィレットの上限)
    # 上流 [-L_up, 0]: マイタ区間 (1 区間) から上流へ (フィレット区間は Δx 上限、U0 は HX)
    up = run_march([(xm_ if use_miter else 0.0, xf1), (xf1, -LUP)], h0k, hx, [cap_f, None], P)
    # フィレット | sidewall_end: 左 = 0 から下流 ([0, x_f2] は Δx 上限)、右 = L_sw から上流 ([L_sw - ZONE, L_sw])
    A = _meet(dict(segs=[(0.0, xf2)], h0=h0k, hmax=hx, caps=[cap_f]),
              dict(segs=[(LSW, LSW - ZONE)], h0=hw["sidewall_end"], hmax=hx, caps=[None]), xf2, LSW - ZONE, [], P, int(P["MEET_SCAN"]))
    # B: 生成器と同じ (中点で出会う)
    xmB = 0.5 * (LSW + LCOWL)
    B1 = run_march([(LSW, LSW + LWAKE), (LSW + LWAKE, xmB)], hw["sidewall_end"], hx, [dxw, None], P)
    B2 = run_march([(LCOWL, xmB)], hw["cowl_base"], hx, [None], P)
    # cowl_base | vehicle_base: 左 = L_cowl から下流 (後流区間の Δx 上限・L_cowl + ZONE)、右 = L_ramp から上流 (L_ramp - ZONE)。x0 は途中の物理 station
    cl = [(LCOWL, LCOWL + LWAKE), (LCOWL + LWAKE, LCOWL + ZONE)]
    cr = [(LR, LR - ZONE)] if ramp_zone else []
    C = _meet(dict(segs=cl, h0=hw["cowl_base"], hmax=hxf, caps=[dxw, None]),
              dict(segs=cr, h0=hw["vehicle_base"], hmax=hxf, caps=[None] * len(cr)),
              LCOWL + ZONE, (LR - ZONE) if ramp_zone else LR, [x0], P, int(P["MEET_SCAN"]))
    # D: L_ramp から下流 (ベース後流の Δx 上限・L_ramp + ZONE・x_out)
    dsegs = [(LR, LR + lwb)] + ([(LR + lwb, LR + ZONE), (LR + ZONE, XOUT)] if ramp_zone else [(LR + lwb, XOUT)])
    D = run_march(dsegs, hw["vehicle_base"], hxf, [dxwb] + [None] * (len(dsegs) - 1), P)
    # 区間 (始, 終) -> (n, 向きつきの r)
    laws = {}

    def put(seg, nr):
        a_, b_ = seg; n, r = nr
        laws[(min(a_, b_), max(a_, b_))] = (n, r if b_ > a_ else -r)
    for seg, nr in zip([(xm_ if use_miter else 0.0, xf1), (xf1, -LUP)], up[0]): put(seg, nr)
    for side in ("left", "right"):
        for seg, nr in A[side]: put(seg, nr)
        for seg, nr in C[side]: put(seg, nr)
    for seg, nr in zip([(LSW, LSW + LWAKE), (LSW + LWAKE, xmB)], B1[0]): put(seg, nr)
    for seg, nr in zip([(LCOWL, xmB)], B2[0]): put(seg, nr)
    for seg, nr in zip(dsegs, D[0]): put(seg, nr)
    if kink: laws[kink] = (1, 1.0)
    xp = sorted({p for k_ in laws for p in k_})
    xp = [x for i, x in enumerate(xp) if i == 0 or x - xp[i - 1] > 1e-12 * H]
    xlaw = []
    for i in range(len(xp) - 1):
        key = next((k_ for k_ in laws if abs(k_[0] - xp[i]) <= 1e-12 * H and abs(k_[1] - xp[i + 1]) <= 1e-12 * H), None)
        if key is None: raise RuntimeError("物理区間 [%.9g, %.9g] の分布が無い" % (xp[i], xp[i + 1]))
        xlaw.append(laws[key])

    def region(xa, xb):
        if xb <= 1e-12 * H: return "U0"
        if xb <= LSW + 1e-12 * H: return "A"
        if xb <= LCOWL + 1e-12 * H: return "B"
        if xb <= LR + 1e-12 * H: return "C"
        return "D"
    xs = stations_from_laws(xp, xlaw)
    regions = [region(xs[i], xs[i + 1]) for i in range(len(xs) - 1)]
    dxs = np.diff(xs); xr_ = np.maximum(dxs[1:] / dxs[:-1], dxs[:-1] / dxs[1:]); iw = int(np.argmax(xr_))
    phys = [dict(x_a=xp[i] / H, x_b=xp[i + 1] / H, region=region(xp[i], xp[i + 1]), n=int(xlaw[i][0]), r=float(xlaw[i][1])) for i in range(len(xp) - 1)]
    far = np.array([r_ in "CD" for r_ in regions])
    return dict(xs=xs, regions=regions, phys=phys, n_stations=len(xs), dx_min=float(dxs.min()) / H, dx_max=float(dxs.max()) / H,
                dx_max_far=float(dxs[far].max()) / H if far.any() else 0.0,
                x_adjacent_ratio_design=float(xr_.max()), x_adjacent_ratio_at=float(xs[iw + 1]) / H,
                x_adjacent_ratio_gt=[[float(xs[i + 1]) / H, float(xr_[i])] for i in np.flatnonzero(xr_ > P["G"] * (1 + 1e-9))[:12]],
                meet=dict(A=dict(x=A["m"] / H, ratio=A["ratio"], ratio_run_max=A["ratio_run_max"]),
                          C=dict(x=C["m"] / H, ratio=C["ratio"], ratio_run_max=C["ratio_run_max"])),
                x_m=xm_ / H, x0=x0 / H, L_ramp=LR / H, x_out=XOUT / H, fillet_cap=cap_f / H, kink_dx=(0.0 - xm_) / H if use_miter else None)


# ---------------------------------------------------------------- 断面の同値類 (生成器 build の規則)
class GBump(hjm.Bump):
    """生成器の Bump (gmsh の 1D で較正)。gmsh が閉じていたら開き直す (build の ESTIMATE が finalize するため)"""
    def spacing(self, n, coef):
        if not gmsh.isInitialized():
            gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0)
        return super().spacing(n, coef)


def section_classes(P, G, xs, sw_flags, bump, far_fn=None):
    """station 列 xs の全域から、接続模型の同値類の節点数を決める。build() の「断面の分布 (節点数は同値類ごとに 1 つ)」と
    「B4b-5: 壁別 h1 の事前検査」と同じ規則・同じ演算 (断面隅フィレット 0 の経路だけ。tilt は 1)。sw_flags[s] は station s の断面に SW があるか。
    far_fn は遠方帯の区間数の則 (既定 = 生成器の far_n)。P の HFAR_Y・HFAR_Z (無ければ HFAR) で遠方帯ごとの上限を変えられる (感度用)。
    戻り値: (counts, info)。counts は build の counts と同じ鍵、info は解けない組合せ・両端分布の範囲・遠方の最大間隔"""
    far_fn = far_fn or hjm.far_n
    Pj = split_P(P); H = G.H; g = P["G"]; DR = G.DR; g_t = 1.0 + 0.75 * (g - 1.0)
    if any(r_ != 0 for p in G.prof for r_ in p[1:]): raise NotImplementedError("断面隅フィレット (r > 0) は未対応")
    hw = {t: v * H for t, v in hjm.wall_h1(Pj).items()}
    h_r, h_v, h_si, h_so, h_ci, h_co = (hw[t] for t in ("ramp", "vehicle", "sidewall_in", "sidewall_out", "cowl_in", "cowl_out"))
    h_ring = min(h_r, h_si, h_ci)
    NLw = {t: hjm.prog_n(DR, hw[t], g) + 4 for t in hw}
    NL = hjm.prog_n(DR, h_ring, g) + 4
    NL_co, NL_so = NLw["cowl_out"], NLw["sidewall_out"]
    xa_ = np.array(xs); Lmax = float(np.max(G.yr(xa_) - G.yc(xa_)))
    ccs = np.sqrt(1 + G.dyc(xa_) ** 2); LTmax = float(np.max(G.yc(xa_) - G.yo(xa_)))
    tan_max = P["AR_TAN"] * min(h_si, h_so)
    cap_sw = bool(P.get("NY_CAP_SW_ONLY", 0)); in_sw = xa_ <= G.LSW + 1e-12 * H
    if not cap_sw:                                      # 生成器の規則: 全 station の最大の高さに AR の上限を掛ける
        NY = bump.fit_n(min(h_r, h_ci, h_v) / Lmax, g, hmax_rel=tan_max / Lmax)
    else:                                               # 感度: 上限は側壁のある station だけ、公比の条件は全 station の最大の高さで
        Lsw = float(np.max(G.yr(xa_[in_sw]) - G.yc(xa_[in_sw])))
        NY = max(bump.fit_n(min(h_r, h_ci, h_v) / Lsw, g, hmax_rel=tan_max / Lsw), bump.fit_n(min(h_r, h_ci, h_v) / Lmax, g))
    NSW = bump.fit_n(min(h_si, h_so) / G.TSW, g)
    NTC = bump.fit_n(min(h_ci, h_co) * float(ccs.min()) / LTmax, g)
    NZ = max(int(P["NZ"]), hjm.prog_n(G.ZW, h_ring, g_t) + 1)
    cmax = float(np.max(np.sqrt(1 + G.dyr(xa_) ** 2))); s2_rel = h_ring * 0.5 * (1 + cmax) / DR
    NR = 2 * NL
    while True:
        cf_ = bump.coef(NR + 1, s2_rel); sp_ = bump.spacing(NR + 1, cf_)
        if np.maximum(sp_[1:] / sp_[:-1], sp_[:-1] / sp_[1:]).max() <= g_t: break
        NR += 2
    kv = P["KV"] * h_so

    def band_out(hh, nl):
        hh = min(hh, 0.999 * DR / nl); return hh * hjm.prog_r(DR, nl, hh) ** nl
    NFY, far_y = far_fn([(float(G.yo(x_)) - DR - G.YBOT, band_out(h_co * math.sqrt(1 + G.dyo(x_) ** 2), NL_co)) for x_ in xs],
                        P.get("HFAR_Y", P["HFAR"]) * H, g_t)
    NFZ, far_z = far_fn([(G.ZFAR - G.ZO - DR, band_out(kv if j >= 3 else h_so, NL_so)) for j in range(5)], P.get("HFAR_Z", P["HFAR"]) * H, g_t)
    asym = dict(SW=not (h_so == h_si == h_r == h_ci), TC=h_ci != h_co, Y=not (h_r == h_v == h_ci))
    g_two = 1.0 + 0.95 * (g - 1.0)
    rc = float(P.get("RING_CORE_RATIO", 1.0))

    def bl(n, Led, h):
        cf = bump.coef(n, h / Led); return (n, "Progression", 1.0) if cf is None else (n, "Bump", cf)
    conflicts, svs = [], []
    for x, swp in zip(xs, sw_flags):                   # build の seam_vals (断面隅フィレット 0)
        pt, _, (ru, rl) = G.section(x)
        elen = {"e2": math.dist(pt["WU"], pt["CUR"]), "e6": math.dist(pt["WL"], pt["CLR"])}
        c = math.sqrt(1 + G.dyr(x) ** 2); cc = math.sqrt(1 + G.dyc(x) ** 2); co = math.sqrt(1 + G.dyo(x) ** 2)
        k2 = k6 = 1.0
        s = dict(e2=hjm.diag_first(h_si, h_r, c, k2, elen["e2"], DR), e6=hjm.diag_first(h_si, h_ci, cc, k6, elen["e6"], DR), e4=h_r * c, e9=h_ci * cc)
        a_top, a_bot, a_v = (0.5 * h_ * (c + cc) * P["E5K"] for h_ in (h_r, h_ci, h_v))
        sv = dict(x=x, c=c, cc=cc, co=co, s=s, a_top=a_top, a_bot=a_bot, a_v=a_v, LT=G.yc(x) - G.yo(x))
        conflicts += hjm.station_conflicts(x / H, c, cc, s, elen, DR, hw, a_top, a_v, swp)
        sc = dict(e2=s["e2"], e6=s["e6"])
        if abs(rc - 1.0) > 1e-12:
            for e in ("e2", "e6"):
                n_, typ_, cf_ = bl(NR + 1, elen[e], s[e])
                base = bump.spacing(n_, cf_) * elen[e] if typ_ == "Bump" else np.full(n_ - 1, elen[e] / (n_ - 1))
                sc[e] = float(hjm.ring_diag_spacing(base, rc, g)[-1])
        try:
            F = hjm.match_spacing([a_top, a_bot, s["e2"], s["e6"], (s["e2"] + s["e6"]) / (2 * math.sqrt(2))], g)
            if abs(rc - 1.0) > 1e-12: hjm.match_spacing([sc["e2"], sc["e6"], (sc["e2"] + sc["e6"]) / (2 * math.sqrt(2))], g)
            hjm.match_spacing([s["e4"], s["e9"]], g)
            a3, a4 = hjm.sw_end_spacing(s["e2"], s["e6"], F, kv, g, int(P.get("SW_END_RULE", 0)))
            if not asym["SW"]: hjm.match_spacing([F, h_so], g)
            if not asym["TC"]: hjm.match_spacing([s["e6"], a_bot, h_co * co], g)
            v21t = hjm.match_spacing([s["e6"], a_bot], g) if asym["TC"] else None
        except ValueError as e_:
            conflicts.append(dict(kind="seam", where="継ぎ目の間隔の連動 (match_spacing)", x=x / H, reason=str(e_))); svs.append(sv); continue
        sv.update(F=F, a3=a3, a4=a4, v21t=v21t); svs.append(sv)
    asym_ranges = {}

    def class_n(cls, items, n_min):
        rng = [hjm.two_end_n_range(L_, a_, b_, g_two) for L_, a_, b_ in items]
        if any(r_ is None for r_ in rng) or not rng:
            conflicts.append(dict(kind="count", where=cls, x=float("nan"), reason="両端の間隔が違う分布の区間数が取れない")); return n_min
        lo, hi = max(max(r_[0] for r_ in rng), n_min), min(r_[1] for r_ in rng)
        asym_ranges[cls] = (lo, hi)
        if lo > hi: conflicts.append(dict(kind="count", where=cls, x=float("nan"), reason="両端の間隔が違う分布の区間数の範囲が全 station で交わらない (下限 %d > 上限 %d)" % (lo, hi)))
        return lo
    ok_sv = [v for v in svs if "F" in v]
    if asym["SW"] and ok_sv:
        NSW = class_n("SW 列 (hz?1)", [(G.TSW, v["F"], h_so) for v in ok_sv] + [(G.TSW, v[k], kv) for v in ok_sv for k in ("a3", "a4")], 1) + 1
    if asym["TC"] and ok_sv:
        NTC = class_n("CW 列 (vy2k)", [(v["LT"], h_co * v["co"], h_ci * v["cc"]) for v in ok_sv] + [(v["LT"], h_co * v["co"], v["v21t"]) for v in ok_sv], 1) + 1
    if asym["Y"] and ok_sv:                             # build の y_class と同じ (側壁帯の y: Bump を曲げて上端だけ変える)
        def y_class(Ly, v, n):
            n_, typ_, cf_ = bl(n, Ly, v["a_bot"])
            base = bump.spacing(n_, cf_) * Ly if typ_ == "Bump" else np.full(n_ - 1, Ly / (n_ - 1))
            top32 = v["a_top"] if v["a_top"] == v["a_v"] else math.sqrt(v["a_top"] * v["a_v"])
            return [hjm.ring_diag_spacing(base, top32 / base[0], g), hjm.ring_diag_spacing(base, v["a_v"] / base[0], g)]
        ly = lambda x_: float(G.yr(x_) - G.yc(x_)); ny0 = NY
        capped = (lambda v: v["x"] <= G.LSW + 1e-12 * H) if cap_sw else (lambda v: True)
        while True:
            try:
                dmax = [float(d_.max()) if capped(v) else 0.0 for v in ok_sv for d_ in y_class(ly(v["x"]), v, NY)]     # 全 station で解けること (公比) + 上限
                if max(dmax) <= tan_max: break
            except ValueError:
                pass
            NY += 1
            if NY > 3 * ny0:
                conflicts.append(dict(kind="count", where="側壁帯の y (e5・vy3k)", x=float("nan"), reason="両端の間隔が違う分布が NY %d..%d で解けない" % (ny0, NY))); break
        asym_ranges["側壁帯の y (e5・vy3k)"] = (ny0 - 1, NY - 1)
    counts = dict(NL=NL, NL_co=NL_co, NL_so=NL_so, NL_END={t: NLw[t] for t in hjm.WALL_ENDS}, NR=NR, NY=NY, NZ=NZ, NSW=NSW, NTC=NTC, NFY=NFY, NFZ=NFZ)
    info = dict(conflicts=conflicts, asym=asym, asym_ranges={k: list(v) for k, v in asym_ranges.items()}, Lmax=Lmax / H, LTmax=LTmax / H,
                far_spacing_max=dict(y=far_y / H, z=far_z / H), wall_layers=NLw, tan_max=tan_max / H, svs=svs)
    return counts, info


def top_classes(P, G, xs, svs, far_fn=None):
    """全体の新しい行の同値類 (機体上の外部流): 行 tb (機体上面の壁帯 b..b+DR、NL_vt 区間 = prog_n(DR, h_vt, 1.2) + 4、生成器の壁帯と同じ則)、
    行 tc (遠方帯 b+DR..b+top_depth、far_fn [既定 far_n]、上限 HFAR_T)、行 v (y_r..b、両端分布の区間数)。
    行 v の端: 下端 = 行 de の上端 (vy33・vy34 の Bump の端 a_v、station ごと)、上端 = 行 tb の第一間隔 h_vt √(1 + b'²)。
    区間数は全 station の two_end_n_range の共通部分の最小、かつ最大間隔 <= AR_TAN × h_vehicle_side (機体側面 z = z_o の壁層セルの接線幅、生成器の NY と同じ考え)。
    共通部分が空なら解けない組合せとして返し、区間数は下限 (長い辺の要求) を使う (呼び出し側で仮定として明示する)"""
    far_fn = far_fn or hjm.far_n
    H = G.H; g = P["G"]; DR = G.DR; g_t = 1.0 + 0.75 * (g - 1.0); g_two = 1.0 + 0.95 * (g - 1.0)
    hwf = {t: v * H for t, v in wall_h1_full(P).items()}
    h_vt, h_vs = hwf["vehicle_top"], hwf["vehicle_side"]
    NL_vt = hjm.prog_n(DR, h_vt, g) + 4
    xa_ = np.array(xs); cb = np.sqrt(1 + G.db(xa_) ** 2)

    def band_out(hh, nl):
        hh = min(hh, 0.999 * DR / nl); return hh * hjm.prog_r(DR, nl, hh) ** nl
    NT, far_t = far_fn([(G.TOPD - DR, band_out(h_vt * float(c_), NL_vt)) for c_ in cb], P.get("HFAR_T", P["HFAR"]) * H, g_t)
    Lv = G.b(xa_) - G.yr(xa_); top = h_vt * cb; bot = np.array([v["a_v"] for v in svs])
    rng = [hjm.two_end_n_range(float(L_), float(a_), float(t_), g_two) for L_, a_, t_ in zip(Lv, bot, top)]
    if any(r_ is None for r_ in rng): raise ValueError("行 v: 区間数が取れない station がある")
    lo_i = int(np.argmax([r_[0] for r_ in rng])); hi_i = int(np.argmin([r_[1] for r_ in rng]))
    n_lo, n_hi = rng[lo_i][0], rng[hi_i][1]
    tan_v = P["AR_TAN"] * h_vs; n = n_lo
    longest = np.argsort(-Lv)[:8]
    while max(float(hjm.two_end_spacing(float(Lv[i]), n, float(bot[i]), float(top[i]), g_two).max()) for i in longest) > tan_v:
        n += 1
    conflict = None
    if n > n_hi:
        conflict = dict(kind="count", where="行 v (vy4k: y_r..b)", n_lo=int(n), n_hi=int(n_hi),
                        x_long=float(xs[lo_i]) / H, L_long=float(Lv[lo_i]) / H, x_short=float(xs[hi_i]) / H, L_short=float(Lv[hi_i]) / H,
                        ends_um=[float(bot[hi_i]) * 1e6, float(top[hi_i]) * 1e6],
                        reason="長い辺 (上流の機体高さ) が要る区間数が、短い辺 (ベース近くの t_base) に入る区間数の上限を超える (端の間隔と公比 <= %.3g)" % g_two)
    return dict(NV=int(n) + 1, NL_vt=int(NL_vt), NT=int(NT)), dict(conflict=conflict, NV_range=[int(n_lo), int(n_hi)], far_t=far_t / H,
                                                                  Lv_max=float(Lv.max()) / H, Lv_min=float(Lv.min()) / H, tan_v=tan_v / H)


# ---------------------------------------------------------------- 全体の集計
def full_count(P, contours, meta, bump):
    """全体を数える。戻り値は report の dict (節点・ヘキサ・一意辺、区間別・行列別の内訳、同値類、メモリ、予算判定、仮定に効く値)"""
    G = FullGeom(P, contours["ramp"], contours["cowl"], meta); H = G.H
    far_fn = FAR_FN[P["FAR_LAW"]]
    lay = full_layout(P, G); xs, regions = lay["xs"], lay["regions"]
    edge, blk = full_tables(); present = full_present(blk)
    ns = len(xs)
    sw_flags = [("SW" in ((present[regions[s - 1]] if s > 0 else set()) | (present[regions[s]] if s < ns - 1 else set()))) for s in range(ns)]
    cnt, info = section_classes(P, G, xs, sw_flags, bump, far_fn)
    top, tinfo = top_classes(P, G, xs, info["svs"], far_fn)
    cnt.update(top)
    ec = full_edge_counts(cnt)
    res = count_complex(regions, present, blk, edge, ec, order=BLOCK_ORDER)
    hwf = wall_h1_full(P)
    # 遠方の AR の概算: 遠方の最大間隔 (C/D の x、遠方帯の y・z) / 全域に伝わる最小の間隔 (細かくする壁の第一層。壁帯の分布は行・列ごとに全幅へ伝わる)
    s_min = min(hwf[t] for t in ("ramp", "cowl_in", "cowl_out", "sidewall_in"))          # /H
    far_max = max(lay["dx_max_far"], info["far_spacing_max"]["y"], info["far_spacing_max"]["z"], tinfo["far_t"])     # /H
    rep = dict(nodes=res["nodes"], hexes=res["hexes"], unique_edges=res["edges"], quads=res["quads"], euler=res["euler"],
               n_stations=lay["n_stations"], counts=cnt, wall_h1_um={t: v * H * 1e6 for t, v in hwf.items()},
               nodes_by_region=res["nodes_by_region"], hexes_by_region=res["hexes_by_region"], stations_by_region=res["stations_by_region"],
               rows_cols=breakdown_rows_cols(res), nodes_by_block=res["nodes_by_block"],
               layout={k: v for k, v in lay.items() if k not in ("xs", "regions")},
               conflicts=[{k: v for k, v in c_.items()} for c_ in info["conflicts"]] + ([tinfo["conflict"]] if tinfo["conflict"] else []),
               asym_ranges=info["asym_ranges"], NV_range=tinfo["NV_range"], far_spacing_max=dict(info["far_spacing_max"], t=tinfo["far_t"]),
               far_law=P["FAR_LAW"], Lmax=info["Lmax"], Lv=[tinfo["Lv_min"], tinfo["Lv_max"]],
               ar_far_estimate=far_max / s_min, ar_far_parts=dict(dx_max_far=lay["dx_max_far"], far_max=far_max, s_min_um=s_min * H * 1e6),
               memory_gib=memory(res["nodes"]), budget_nodes=BUDGET_NODES, budget_ratio=res["nodes"] / BUDGET_NODES,
               within_budget=res["nodes"] <= BUDGET_NODES,
               geometry=dict(L_ramp=G.L_ramp / H, x_out=G.XEND / H, x0=G.x0 / H, x_m=G.x_m / H, x_f1=G.x_f1 / H, x_f2=G.x_f2 / H,
                             y_veh=G.y_veh / H, theta_e_deg=math.degrees(math.atan(G.th_e)), y_bot=G.YBOT / H, z_far=G.ZFAR / H,
                             z_o=G.ZO / H, top_depth=P["TOP_DEPTH"], t_base=P["T_BASE"]))
    return rep


def make_P(preset="design", sets=()):
    """全体の P (生成器の P0 + プリセット + FULL_P0 + --set)"""
    P = dict(hjm.P0); P.update(FULL_P0)
    if preset == "final-common2": P.update(FINAL_COMMON2)
    else: P.update(FULL_PRESETS[preset])
    P["PRESET"] = preset
    for kv in sets:
        k, v = kv.split("=")
        P[k] = v if k == "FAR_LAW" else float(v)
    return P


def load_meta(d):
    return json.load(open(Path(d) / "contours_meta.json"))


SENS = [  # (名前, 説明, 変更) — 既定の design から 1 つずつ
    ("hx_far_0.24", "HX_FAR 0.16 -> 0.24 H (C・D の x)", dict(HX_FAR=0.24)),
    ("hx_0.15", "HX 0.10 -> 0.15 H (U0・A・B の x)", dict(HX=0.15)),
    ("hfar_0.20", "HFAR 0.10 -> 0.20 H (遠方帯の全部: 行 a・列 3・行 tc)", dict(HFAR=0.20)),
    ("hfar_0.50", "HFAR 0.10 -> 0.50 H (遠方帯の全部。codex 2026-10-06 は AR で却下した値)", dict(HFAR=0.50)),
    ("hfar_a_0.50", "行 a の遠方だけ HFAR 0.50 H (プルームが横切る帯を粗くする。他は 0.10)", dict(HFAR_Y=0.50)),
    ("far_capped", "遠方帯の分布則を「公比 g_t で HFAR まで増やして一様」に (生成器の等比 1 本の則の代わり。仮想)", dict(FAR_LAW="capped")),
    ("far_capped_a_0.50", "遠方帯を capped かつ行 a だけ HFAR 0.50 H", dict(FAR_LAW="capped", HFAR_Y=0.50)),
    ("vehicle_top_16um", "vehicle_top 64 -> 16 µm (行 v の両端分布が解ける組合せ)", dict(H1_vehicle_top=1.6e-4)),
    ("no_miter_station", "外壁のマイタ頂点を station にしない (フィレット 12 区間・内壁の折れ点 0 だけ)", dict(MITER_STATION=0)),
    ("no_ramp_zone", "L_ramp ± ZONE の物理 station を置かない (march の等比 1 本が長い区間に掛かる)", dict(RAMP_ZONE=0)),
    ("lwake_b_0.08", "ベース後流の Δx 上限区間 0.02 -> 0.08 H (= 4 t_base)", dict(LWAKE_B=0.08)),
    ("ny_cap_sw_only", "NY の AR 上限 (900 h1) を側壁のある x <= L_sw だけに掛ける (生成器の規則の変更)", dict(NY_CAP_SW_ONLY=1)),
]


def summarize(rep):
    """感度表の 1 行"""
    c = rep["counts"]
    return dict(nodes=rep["nodes"], hexes=rep["hexes"], unique_edges=rep["unique_edges"], n_stations=rep["n_stations"],
                budget_ratio=rep["budget_ratio"], within_budget=rep["within_budget"],
                NY=c["NY"], NZ=c["NZ"], NFY=c["NFY"], NFZ=c["NFZ"], NV=c["NV"], NT=c["NT"], NL_vt=c["NL_vt"], NSW=c["NSW"], NTC=c["NTC"], NR=c["NR"],
                ar_far_estimate=rep["ar_far_estimate"], x_adjacent_ratio_design=rep["layout"]["x_adjacent_ratio_design"],
                n_conflicts=len(rep["conflicts"]), memory_gib=rep["memory_gib"])


# ---------------------------------------------------------------- 自己検査
# 記録済みの接続模型 (同値類の節点数と区間ごとの x 区間数 -> count_mesh の値)。B4a 板厚 0.005 H scale 1 は gmsh の実数 (design_tests n18)、
# design / default / final は B4b-5 の ESTIMATE (design は gmsh 生成の PASS 格子と同数、notes/investigations/2026-10-07-sern-b4b5/)
RECORDED = [
    ("b4a_tc0.005_s1", dict(NZ=50, NR=62, NY=85, NSW=19, NTC=19, NL_co=30, NL_so=30, NFY=20, NFZ=20), [11, 28, 20, 23, 25, 20, 26], (4349372, 4259394)),
    ("design", dict(NZ=50, NR=62, NY=85, NSW=12, NTC=19, NL_co=30, NL_so=22, NFY=22, NFZ=22), [8, 28, 20, 23, 25, 20, 26], (4034091, 3949416)),
    ("default", dict(NZ=50, NR=62, NY=85, NSW=19, NTC=19, NL_co=30, NL_so=30, NFY=22, NFZ=22), [11, 28, 20, 23, 25, 20, 26], (4433456, 4342320)),
    ("final", dict(NZ=65, NR=92, NY=337, NSW=43, NTC=39, NL_co=30, NL_so=22, NFY=22, NFZ=22), [11, 28, 20, 23, 25, 20, 26], (15559061, 15326600)),
]


def junction_regions(nx, labels="AABBBCC"):
    return [rg for rg, n_ in zip(labels, nx) for _ in range(n_)]


def selfcheck_pure():
    """gmsh 不要の自己検査: 汎用の数え上げが接続模型の範囲で count_mesh と一致、全体の表を接続模型の行に限っても一致、Euler 標数 1"""
    out = []
    edge, blk = full_tables()
    for name, cnt, nx, ref in RECORDED:
        regs = junction_regions(nx); ec = hjm.edge_counts(cnt)
        cm = hjm.count_mesh(regs, ec)
        a = count_complex(regs, hjm.PRESENT, hjm.BLK, hjm.EDGE, ec)
        b_ = count_complex(regs, hjm.PRESENT, blk, edge, ec)        # 全体の表 (新しい行は存在しない区間名だけ)
        out.append(dict(model=name, count_mesh=list(cm), recorded=list(ref), generic=[a["nodes"], a["hexes"]], full_table=[b_["nodes"], b_["hexes"]],
                        euler=[a["euler"], b_["euler"]],
                        ok=tuple(cm) == tuple(ref) and (a["nodes"], a["hexes"]) == tuple(ref) and (b_["nodes"], b_["hexes"]) == tuple(ref) and a["euler"] == b_["euler"] == 1))
    return out


def selfcheck_generator(contours, bump):
    """gmsh を使う自己検査: build(ESTIMATE=1) と、station 列 (junction_layout) と同値類 (section_classes) が一致すること。
    接続模型 3 プリセット (default/design/final) と、design を参照ランプ輪郭の終端まで延ばした全長 (SECX と同じ XEND)"""
    out = []
    cases = [("default", {}), ("design", hjm.PRESETS["design"]), ("final", hjm.PRESETS["final"]), ("design_full_length", dict(hjm.PRESETS["design"], SECX=[1.0]))]
    for name, pre in cases:
        P = dict(hjm.P0); P.update(pre); P["ESTIMATE"] = 1
        est = hjm.build("/dev/null", dict(P), contours, quiet=True)
        if P.get("SECX") is not None:                     # build と同じ XEND の延長
            xe = float(contours["ramp"][-1, 0]) / P["H"]
            while xe * P["H"] > float(contours["ramp"][-1, 0]): xe = float(np.nextafter(xe, 0.0))
            P["XEND"] = xe
        lay = junction_layout(P)
        G = hjm.Geom(P, contours["ramp"], contours["cowl"])
        ns = len(lay["xs"])
        sw = [("SW" in ((hjm.PRESENT[lay["regions"][s - 1]] if s > 0 else set()) | (hjm.PRESENT[lay["regions"][s]] if s < ns - 1 else set()))) for s in range(ns)]
        cnt, info = section_classes(P, G, lay["xs"], sw, bump)
        nx = [n for n, _ in lay["xlaw"]]
        cm = hjm.count_mesh(lay["regions"], hjm.edge_counts(cnt))
        same_counts = {k: (cnt[k], est["counts"][k]) for k in est["counts"] if cnt[k] != est["counts"][k]}
        n_conf = (len(info["conflicts"]), len(est["conflicts"]))
        out.append(dict(case=name, nx=nx, nx_build=est["nx"], stations=(ns, est["stations"]), nodes=(cm[0], est["nodes"]), hexes=(cm[1], est["hexes"]),
                        counts_diff=same_counts, conflicts=n_conf,
                        ok=nx == est["nx"] and ns == est["stations"] and cm == (est["nodes"], est["hexes"]) and not same_counts and n_conf[0] == n_conf[1]))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--contours", required=True)
    ap.add_argument("--preset", default="design", choices=sorted(FULL_PRESETS) + ["final-common2"])
    ap.add_argument("--set", nargs="*", default=[])
    ap.add_argument("--sens", action="store_true", help="感度 (SENS の 1 因子ずつ)")
    ap.add_argument("--final-ref", action="store_true", help="参考: 隅を共通 2 µm にそろえた final")
    ap.add_argument("--selfcheck", action="store_true")
    ap.add_argument("--out")
    a = ap.parse_args()
    contours = hjm.load_contours(a.contours); meta = load_meta(a.contours); bump = GBump()
    report = dict(contours={k: contours[k + "_src"] for k in ("ramp", "cowl")}, meta=meta, mem_model=MEM_MODEL)
    if a.selfcheck:
        report["selfcheck_pure"] = selfcheck_pure()
        report["selfcheck_generator"] = selfcheck_generator(contours, bump)
    P = make_P(a.preset, a.set)
    base = full_count(P, contours, meta, bump)
    report["params"] = {k: P[k] for k in sorted(P) if not isinstance(P[k], (list, dict))}
    report["base"] = base
    if a.sens:
        report["sens"] = {}
        for name, desc, ch in SENS:
            Pv = dict(P); Pv.update(ch)
            report["sens"][name] = dict(desc=desc, change=ch, **summarize(full_count(Pv, contours, meta, bump)))
    if a.final_ref:
        Pf = make_P("final-common2", a.set)
        fr = full_count(Pf, contours, meta, bump)
        report["final_common2"] = dict(desc="隅 (ramp・cowl_in・sidewall_in) を共通 2 µm、cowl_out 16 µm、緩める壁は design と同じ、HX 0.05・HX_FAR 0.08 H",
                                       **summarize(fr), wall_h1_um=fr["wall_h1_um"], conflicts=fr["conflicts"], NV_range=fr["NV_range"])
        Pf["NY_CAP_SW_ONLY"] = 1
        fr2 = full_count(Pf, contours, meta, bump)
        report["final_common2_nycap_sw"] = dict(desc="同上 + NY の AR 上限を x <= L_sw だけに (生成器の規則の変更)", **summarize(fr2), conflicts=fr2["conflicts"])
    if gmsh.isInitialized(): gmsh.finalize()
    if a.out: json.dump(report, open(a.out, "w"), indent=1, ensure_ascii=False, default=float)
    s = summarize(base)
    print(json.dumps(dict(nodes=s["nodes"], hexes=s["hexes"], unique_edges=s["unique_edges"], n_stations=s["n_stations"],
                          budget_ratio=round(s["budget_ratio"], 3), counts=base["counts"], euler=base["euler"],
                          memory_gib={k: round(v, 2) for k, v in base["memory_gib"].items()}, conflicts=len(base["conflicts"])), indent=1, ensure_ascii=False))
    if a.selfcheck:
        print("selfcheck_pure", all(r["ok"] for r in report["selfcheck_pure"]), "selfcheck_generator", all(r["ok"] for r in report["selfcheck_generator"]))


if __name__ == "__main__":
    main()
