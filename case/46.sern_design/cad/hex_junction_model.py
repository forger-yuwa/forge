#!/usr/bin/env python3
"""接続模型 (全ヘキサ): 側壁終端・カウル終端を含む 3D 接続を gmsh の transfinite ブロックで切る。
設計は plans/active/tooling-sern-mesh-blocking.md §4.10 (19 ブロック/断面 × 3 区間)。

  N  : ダクト内のバタフライ (TOP/SIDE/BOT/CORE)。壁曲線は弧長一様の 1 本スプライン (§4.7)
  SW : 側壁の跡 (x >= L_sw)          CW1/CW2 : カウルの跡 (x >= L_cowl)
  U  : カウル下 (2 列 × 2 帯)         S : 側壁の外 (2 列 × 4 帯)

第一層は「辺長」でなく壁面の 3D 法線に対する距離で決める (§4.10 (e))。端面 (sidewall_end / cowl_base) にも壁層を置く。
境界タグは「参照する体積が 1 つの面」を機械的に拾い、(面の名前, 区間) の表で決める。

usage (mesh venv):  .venv-mesh/bin/python case/46.sern_design/cad/hex_junction_model.py OUT_PREFIX [--scale S] [--h1 H1/H] [--set k=v ...]
出力: OUT.msh / OUT_report.json / OUT_sections.npz (描画用)
"""
import sys, math, json, argparse
import numpy as np
import gmsh

# ---------------------------------------------------------------- 物理入力 (すべて /H)
P0 = dict(H=0.1, ZW=1.0, TSW=0.05, TC=0.02, LSW=1.2, LCOWL=1.6, XEND=2.4, DR=0.08, YBOT=-1.0, ZFAR=2.0,
          RAMP_A=0.30, RAMP_B=-0.02,            # y_r/H = 1 + A x/H + B (x/H)^2
          ZONE=0.3,                             # 端面の上下流に置く物理 station までの距離
          H1=4.0e-5, G=1.2, NZ=41, HX=0.03, HX_FAR=0.03, M_SPL=129,   # HX: 45° の隅セルも AR <= 1000 に入る x 間隔 (h1 の 750 倍)
          CORE_MILD=1.0, MILD1=5.0, MILD2=25.0, E5_BLEND=1.1,
          HFAR=0.10, AR_TAN=900.0)                            # 外部遠方の最大格子幅 (端面層が断面全体に伝播するので、遠方の幅が AR を決める)
# (x/H, 上隅 r/H, 下隅 r/H)。区間内は smoothstep (両端で r' = 0)
# 半径を変える区間の長さは 2 r 以上にする: 短いと r'(x) で壁法線が x に傾き、円弧上の第一層が薄くなる (長さ r で -19 %、2 r で -5 %)
PROFILE = [(0.00, 0.00, 0.06), (0.20, 0.10, 0.06), (0.60, 0.10, 0.06), (0.90, 0.10, 0.06), (1.10, 0.00, 0.00), (1.20, 0.00, 0.00)]


def smooth(t):
    t = min(max(t, 0.0), 1.0); return t * t * (3.0 - 2.0 * t)


# ---------------------------------------------------------------- 1D 分布
def prog_r(L, n, h):
    """n 区間・第一間隔 h・全長 L の公比 (r < 1 も許す)"""
    if abs(n * h - L) < 1e-12 * L: return 1.0
    lo, hi = (1.0 + 1e-12, 50.0) if n * h < L else (1e-3, 1.0 - 1e-12)
    for _ in range(200):
        r = 0.5 * (lo + hi)
        lo, hi = (r, hi) if h * (r ** n - 1) / (r - 1) < L else (lo, r)
    return 0.5 * (lo + hi)


def prog_n(L, h, g):
    """公比 <= g で覆える最小の区間数"""
    if h >= L: return 1
    return max(int(math.ceil(math.log(1.0 + L * (g - 1.0) / h) / math.log(g) - 1e-9)), 1)


class Bump:
    """gmsh の Bump 分布を単位長の直線で実測して較正する (間隔は弧長の放物線、端/中央 = coef)"""
    def __init__(self): self.c = {}

    def spacing(self, n, coef):
        k = (n, float("%.10e" % coef))
        if k not in self.c:
            gmsh.model.add("cal"); g = gmsh.model.geo
            l = g.addLine(g.addPoint(0, 0, 0), g.addPoint(1, 0, 0)); g.mesh.setTransfiniteCurve(l, n, "Bump", k[1]); g.synchronize()
            gmsh.model.mesh.generate(1); _, x, _ = gmsh.model.mesh.getNodes(1, l, includeBoundary=True)
            self.c[k] = np.diff(np.sort(x.reshape(-1, 3)[:, 0])); gmsh.model.remove()
        return self.c[k]

    def coef(self, n, h_rel):
        """端の間隔/全長 = h_rel になる係数。一様以上なら None (一様)"""
        if h_rel >= 0.98 / (n - 1): return None
        mk = ("c", n, float("%.6e" % h_rel))
        if mk in self.c: return self.c[mk]
        lo, hi = 1e-8, 0.999
        for _ in range(60):
            c = math.sqrt(lo * hi); lo, hi = (c, hi) if self.spacing(n, c)[0] < h_rel else (lo, c)
        self.c[mk] = math.sqrt(lo * hi); return self.c[mk]

    def fit_n(self, h_rel, g, n0=11, hmax_rel=1.0):
        n = n0
        while True:
            c = self.coef(n, h_rel)
            if c is None: return n
            h = self.spacing(n, c)
            if np.maximum(h[1:] / h[:-1], h[:-1] / h[1:]).max() <= g and h.max() <= hmax_rel: return n
            n += 2 if n < 60 else 4


def march(lengths, h0, g_t, g, hmax):
    """端面から離れる向きに区間を順にたどり、区間ごとの (区間数, 公比) を決める。最初の区間は第一間隔 = h0。
    以降は「第一間隔 = 前の最終間隔 × 公比」で継ぎ目の間隔比も公比に揃える。なるべく粗く (公比 <= g_t、最終間隔 <= 1.15 hmax)"""
    out, last = [], None
    for L in lengths:
        found = None
        for lim in (g_t, g):
            for n in range(1, int(L / (h0 if last is None else last / g)) + 4):
                if last is None:
                    r = prog_r(L, n, h0); first = h0
                elif n == 1:
                    r = L / last; first = L
                else:
                    f = lambda q: last * q * ((q ** n - 1) / (q - 1) if abs(q - 1) > 1e-12 else n)
                    if f(1.0 / g) > L or f(g) < L: continue
                    lo, hi = 1.0 / g, g
                    for _ in range(100):
                        q = 0.5 * (lo + hi); lo, hi = (q, hi) if f(q) < L else (lo, q)
                    r = 0.5 * (lo + hi); first = last * r
                if not (1.0 / g - 1e-9 <= r <= lim + 1e-9): continue
                end = first * (1.0 if n == 1 else r ** (n - 1))
                if end > 1.15 * hmax: continue
                found = (n, 1.0 if n == 1 else r, end); break
            if found: break
        if not found:                       # 物理 station の間隔が格子間隔より狭い等。一様で入れ、継ぎ目の比は報告側で検査する
            n = max(int(round(L / min(last, hmax))), 1); found = (n, 1.0, L / n)
        out.append(found[:2]); last = found[2]
    return out


# ---------------------------------------------------------------- 形状
class Geom:
    def __init__(s, P):
        s.P = P; H = P["H"]; s.H = H
        for k in ("ZW", "TSW", "TC", "LSW", "LCOWL", "XEND", "DR", "YBOT", "ZFAR", "ZONE"): setattr(s, k, P[k] * H)
        s.ZO = s.ZW + s.TSW; s.YC = 0.0; s.YCL = -s.TC
        s.prof = [(a * H, b * H, c * H) for a, b, c in PROFILE]

    def yr(s, x): xh = x / s.H; return s.H * (1.0 + s.P["RAMP_A"] * xh + s.P["RAMP_B"] * xh * xh)
    def dyr(s, x): return s.P["RAMP_A"] + 2.0 * s.P["RAMP_B"] * x / s.H

    def r(s, x):
        pr = s.prof
        if x <= pr[0][0]: return pr[0][1], pr[0][2]
        for (xa, ua, la), (xb, ub, lb) in zip(pr[:-1], pr[1:]):
            if x <= xb:
                w = smooth((x - xa) / (xb - xa)); return ua + w * (ub - ua), la + w * (lb - la)
        return pr[-1][1], pr[-1][2]

    def section(s, x):
        """断面の点 (y,z) と 3 本の壁曲線 (密な折れ線)。角度は z 軸から y 軸へ測る"""
        ru, rl = s.r(x); yr, yc, ZW, DR = s.yr(x), s.YC, s.ZW, s.DR
        arc = lambda cy, cz, r, a0, a1: np.stack([cy + r * np.sin(np.linspace(a0, a1, 200)), cz + r * np.cos(np.linspace(a0, a1, 200))], 1)
        if ru > 0: up1, up2 = arc(yr - ru, ZW - ru, ru, math.pi / 2, math.pi / 4), arc(yr - ru, ZW - ru, ru, math.pi / 4, 0.0)
        else: up1 = up2 = np.array([(yr, ZW)])
        if rl > 0: lo1, lo2 = arc(yc + rl, ZW - rl, rl, 0.0, -math.pi / 4), arc(yc + rl, ZW - rl, rl, -math.pi / 4, -math.pi / 2)
        else: lo1 = lo2 = np.array([(yc, ZW)])
        walls = dict(e1=np.vstack([[(yr, 0.0)], [(yr, ZW - ru)], up1]),
                     e5=np.vstack([up2, [(yr - ru, ZW)], [(yc + rl, ZW)], lo1]),
                     e8=np.vstack([lo2, [(yc, ZW - rl)], [(yc, 0.0)]]))
        # 内側線は壁を dr だけ内側へずらした線 (r > dr なら同心の円弧 r - dr、r <= dr なら角)。リングの法線方向の厚みが円弧上でも dr になる
        iu, il = max(ru - DR, 0.0), max(rl - DR, 0.0)
        if iu > 0: iu1, iu2 = arc(yr - ru, ZW - ru, iu, math.pi / 2, math.pi / 4), arc(yr - ru, ZW - ru, iu, math.pi / 4, 0.0)
        else: iu1 = iu2 = np.array([(yr - DR, ZW - DR)])
        if il > 0: il1, il2 = arc(yc + rl, ZW - rl, il, 0.0, -math.pi / 4), arc(yc + rl, ZW - rl, il, -math.pi / 4, -math.pi / 2)
        else: il1 = il2 = np.array([(yc + DR, ZW - DR)])
        walls.update(e3=np.vstack([[(yr - DR, 0.0)], [(yr - DR, ZW - max(ru, DR))], iu1])[::-1],
                     e7=np.vstack([iu2, [(yr - max(ru, DR), ZW - DR)], [(yc + max(rl, DR), ZW - DR)], il1])[::-1],
                     e10=np.vstack([il2, [(yc + DR, ZW - max(rl, DR))], [(yc + DR, 0.0)]])[::-1])
        pt = dict(PT0=(yr, 0.0), WU=tuple(up1[-1]), CUR=tuple(iu1[-1]), CUL=(yr - DR, 0.0),
                  WL=tuple(lo1[-1]), CLR=tuple(il1[-1]), CLL=(yc + DR, 0.0), PB0=(yc, 0.0))
        ys = [s.YBOT, s.YCL - DR, s.YCL, yc, yr]; zs = [0.0, ZW, s.ZO, s.ZO + DR, s.ZFAR]
        for j in range(5):
            for k in range(5):
                if (j, k) in ((3, 0), (3, 1), (4, 1), (4, 0)): continue
                pt["G%d%d" % (j, k)] = (ys[j], zs[k])
        return pt, walls, (ru, rl)


# 断面の辺: 名前 -> (始点, 終点)。hz30 = -e8、vy31 = -e5 は別名
N_EDGE = dict(e1=("PT0", "WU"), e2=("WU", "CUR"), e3=("CUR", "CUL"), e4=("CUL", "PT0"), e5=("WU", "WL"), e6=("WL", "CLR"),
              e7=("CLR", "CUR"), e8=("WL", "PB0"), e9=("PB0", "CLL"), e10=("CLL", "CLR"), e11=("CUL", "CLL"))
ALIAS = {"hz30": ("e8", -1), "vy31": ("e5", -1)}
GNAME = {(3, 0): "PB0", (3, 1): "WL", (4, 1): "WU", (4, 0): "PT0"}
gp = lambda j, k: GNAME.get((j, k), "G%d%d" % (j, k))
EDGE = dict(N_EDGE)
for j in range(5):
    for k in range(4):
        if (j, k) != (4, 0) and "hz%d%d" % (j, k) not in ALIAS: EDGE["hz%d%d" % (j, k)] = (gp(j, k), gp(j, k + 1))
for j in range(4):
    for k in range(5):
        if (j, k) != (3, 0) and "vy%d%d" % (j, k) not in ALIAS: EDGE["vy%d%d" % (j, k)] = (gp(j, k), gp(j + 1, k))


def ref(name, sign=1):
    if name.startswith("-"): return ref(name[1:], -sign)
    if name in ALIAS: b, s_ = ALIAS[name]; return b, s_ * sign
    return name, sign


def grid_block(j, k): return ["hz%d%d" % (j, k), "vy%d%d" % (j, k + 1), "-hz%d%d" % (j + 1, k), "-vy%d%d" % (j, k)]


BLK = {"N_TOP": ["e1", "e2", "e3", "e4"], "N_SIDE": ["e5", "e6", "e7", "-e2"], "N_BOT": ["e8", "e9", "e10", "-e6"], "N_CORE": ["-e3", "-e7", "-e10", "-e11"],
       "SW": grid_block(3, 1), "CW1": grid_block(2, 0), "CW2": grid_block(2, 1),
       "U1a": grid_block(0, 0), "U1b": grid_block(1, 0), "U2a": grid_block(0, 1), "U2b": grid_block(1, 1)}
for k, c in ((2, "b"), (3, "c")):
    for j, b in ((0, "a"), (1, "b"), (2, "c"), (3, "de")): BLK["S%s_%s" % (c, b)] = grid_block(j, k)
BLK = {b: [ref(e) for e in loop] for b, loop in BLK.items()}
ALWAYS = [b for b in BLK if b not in ("SW", "CW1", "CW2")]
PRESENT = {"A": set(ALWAYS), "B": set(ALWAYS) | {"SW"}, "C": set(ALWAYS) | {"SW", "CW1", "CW2"}}
WALL_TAGS = ("ramp", "vehicle", "sidewall_in", "sidewall_out", "sidewall_end", "cowl_in", "cowl_out", "cowl_side", "cowl_base")


def xface_tag(e, reg):
    """参照数 1 の x 面 (辺 e を x に掃引した面) のタグ。表に無ければ None = 失敗"""
    if e == "e1" or e == "hz41": return "ramp"
    if e in ("hz42", "hz43"): return "vehicle"
    if e == "e5": return "sidewall_in" if reg == "A" else None
    if e == "vy32": return "sidewall_out" if reg == "A" else None
    if e == "e8": return "cowl_in" if reg in "AB" else None
    if e == "hz31": return "cowl_in" if reg == "B" else None
    if e in ("hz20", "hz21"): return "cowl_out" if reg in "AB" else None
    if e == "vy22": return "cowl_side" if reg in "AB" else None
    if e in ("e4", "e11", "e9", "vy00", "vy10") or (e == "vy20" and reg == "C"): return "symmetry"
    if e.startswith("hz0"): return "far_bottom"
    if e.startswith("vy") and e.endswith("4"): return "far_side"
    return None


# ---------------------------------------------------------------- 本体
def build(out, P, quiet=False):
    G = Geom(P); H = G.H; h1 = P["H1"] * H; g = P["G"]; DR = G.DR
    gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0)
    bump = Bump()
    # ---- station (すべて物理位置) と区間
    xs = {p[0] for p in G.prof} | {G.LSW, G.LCOWL, G.XEND, 0.5 * (G.LSW + G.LCOWL), G.LSW - G.ZONE, G.LCOWL + G.ZONE}
    for (xa, ua, la), (xb, ub, lb) in zip(G.prof[:-1], G.prof[1:]):
        if (ua, la) != (ub, lb):
            for ra, rb in ((ua, ub), (la, lb)):      # r が dr を横切る位置: 内側線の角が円弧に変わり CUR(x) が折れるので station にする
                if (ra - G.DR) * (rb - G.DR) < 0:
                    lo, hi = xa, xb
                    for _ in range(80):
                        xm_ = 0.5 * (lo + hi); rm = ra + smooth((xm_ - xa) / (xb - xa)) * (rb - ra)
                        lo, hi = (xm_, hi) if (rm - G.DR) * (ra - G.DR) > 0 else (lo, xm_)
                    xs.add(0.5 * (lo + hi))
            xs.add(0.5 * (xa + xb))          # 対角辺の第一間隔 s1 = h1 L/dr は r について連続なので、r = 0 の手前に特別な station は要らない
    xs = sorted(x for x in xs if -1e-12 <= x <= G.XEND + 1e-12)
    xs = [x for i, x in enumerate(xs) if i == 0 or x - xs[i - 1] > 1e-9]
    reg = lambda xa, xb: "A" if xb <= G.LSW + 1e-12 else ("B" if xb <= G.LCOWL + 1e-12 else "C")
    ivs = [(xs[i], xs[i + 1], reg(xs[i], xs[i + 1])) for i in range(len(xs) - 1)]
    # ---- 事前判定
    assert 0 < G.LSW < G.LCOWL < G.XEND and G.TC > 0 and G.TSW > 0 and G.YBOT < G.YCL - DR and G.ZO + DR < G.ZFAR
    for x in xs + [0.5 * (a + b) for a, b, _ in ivs]:
        pt, _, (ru, rl) = G.section(x)
        if x >= G.LSW - 1e-12 and (ru > 0 or rl > 0): raise ValueError("x >= L_sw で r != 0")
        if ru + rl + 2 * DR >= G.yr(x) - G.YC or max(ru, rl) + DR >= G.ZW: raise ValueError("フィレットが入らない x=%.4g" % x)
        for e in ("e2", "e6"):
            if math.dist(pt[N_EDGE[e][0]], pt[N_EDGE[e][1]]) < 0.5 * DR: raise ValueError("コアが壁に近すぎる (%s, x=%.4g)" % (e, x))
    # ---- x 分布: 端面から離れる向きに march。A は L_sw から上流へ、B は両端から中央へ、C は L_cowl から下流へ
    h1e = P.get("H1_END", P["H1"]) * H; hx = P["HX"] * H; g_t = 1.0 + 0.75 * (g - 1.0)
    idx = {r_: [i for i, v in enumerate(ivs) if v[2] == r_] for r_ in "ABC"}
    xlaw = {}
    def run(ids, toward_plus, hmax):
        res = march([ivs[i][1] - ivs[i][0] for i in ids], h1e, g_t, g, hmax)
        for i, (n, r) in zip(ids, res): xlaw[i] = (n, r if toward_plus else -r)       # gmsh: 負 = 逆向き
    run(idx["A"][::-1], False, hx)
    xm = 0.5 * (G.LSW + G.LCOWL)
    run([i for i in idx["B"] if ivs[i][1] <= xm + 1e-12], True, hx); run([i for i in idx["B"] if ivs[i][0] >= xm - 1e-12][::-1], False, hx)
    run(idx["C"], True, P["HX_FAR"] * H)
    dxs = []
    for i, (xa, xb, _) in enumerate(ivs):
        n, r = xlaw[i]; q = abs(r) if r > 0 else 1.0 / abs(r); w = q ** np.arange(n); dxs += list((xb - xa) * w / w.sum())
    dxs = np.array(dxs); x_ratio_max = float(np.maximum(dxs[1:] / dxs[:-1], dxs[:-1] / dxs[1:]).max())
    # ---- 断面の分布 (節点数は同値類ごとに 1 つ)
    NL = prog_n(DR, h1, g) + 4                                           # 壁法線の区間数
    Lmax = max(G.yr(x) - G.YC for x in xs); cmax = math.sqrt(1 + max(abs(G.dyr(x)) for x in xs) ** 2)
    tan_max = P["AR_TAN"] * h1            # 壁層セルの接線方向の幅の上限 (リングは最大 45° 傾くので、壁層でも AR <= 1000 に収める)
    NY = bump.fit_n(h1 / Lmax, g, hmax_rel=tan_max / Lmax); NSW = bump.fit_n(h1 / G.TSW, g); NTC = bump.fit_n(h1 / G.TC, g); NZ = int(P["NZ"])
    s_out = DR * (1 - 1 / prog_r(DR, NL, h1))                            # 壁帯の最外間隔 (概算)
    def far_n(L):
        n = prog_n(L, s_out * g_t, g_t) + 1
        while L * (1 - 1 / prog_r(L, n, s_out * g_t)) > P["HFAR"] * H and n < 400: n += 1      # 最終間隔 = L (r-1)/r
        return n
    NFY = far_n(G.YCL - DR - G.YBOT); NFZ = far_n(G.ZFAR - G.ZO - DR)
    m1, m2, mc = P["MILD1"] * h1, P["MILD2"] * h1, P["CORE_MILD"] * h1

    def laws(x, pt, elen):
        """station x での各辺の (節点数, 種別, 係数)"""
        c = math.sqrt(1 + G.dyr(x) ** 2); ru, rl = G.r(x); L = {}
        # gmsh の transfinite 補間は、層 1 を「辺長に対する比 eta = s1/L」で置く (X = (1-eta) 壁 + eta 内側線、eta は両側辺の間で線形)。
        # 内側線は平面壁から dr だけ離れた平行線なので、**リングの全側辺で eta = h1/dr に揃えれば平面部の法線距離は厳密に h1**。
        # 対角辺 e2/e6 は長さが sqrt2 dr - (sqrt2-1) r なので s1 = h1 L/dr (r=0 で sqrt2 h1 = 接する両壁の条件と同じ)。
        # フィレット円弧上は内側線までの法線距離が dr より長いぶん厚くなる (45° 点で L_e2/dr 倍。薄くはならない)
        def tilt(e, wp):                      # 45° 点の壁面法線の x 方向の傾き a_n = d(壁点)/dx . d (r'(x) とランプ勾配で生じる)。法線距離は 1/sqrt(1+a_n^2) 倍に縮む
            a_, b_ = N_EDGE[e]; d = np.array(pt[b_]) - np.array(pt[a_]); d /= np.linalg.norm(d); dx = 1e-5 * H
            x0, x1 = max(x - dx, 0.0), min(x + dx, G.XEND); w_ = lambda xx: np.array(G.section(xx)[0][wp])
            return math.sqrt(1 + (float((w_(x1) - w_(x0)) @ d) / (x1 - x0)) ** 2)
        # 傾きの補正は 3 割だけ掛ける: 全部掛けると e2 側の eta が増えて平面部が厚くなる (円弧上 -9 % -> -6 %、平面部 +3 % 以内)
        k2 = (1 + 0.3 * (tilt("e2", "WU") - 1)) if ru > 0 else 1.0; k6 = (1 + 0.3 * (tilt("e6", "WL") - 1)) if rl > 0 else 1.0
        s = dict(e2=h1 * 0.5 * (1 + c) * k2 * elen["e2"] / DR, e6=h1 * k6 * elen["e6"] / DR, e4=h1 * c, e9=h1)
        for e, sg in (("e2", 1), ("e6", 1), ("e9", 1), ("e4", -1)): L[e] = (NL + 1, "Progression", sg * prog_r(elen[e], NL, s[e]))
        for e in ("e1", "e3", "e8", "e10", "hz00", "hz10", "hz20"): L[e] = (NZ, "Progression", 1.0)
        hend = 0.5 * h1 * (1 + c)                                        # 両端細分: ランプ側とカウル側の平均
        w = smooth((x - (x_r0 - P["E5_BLEND"] * H)) / (P["E5_BLEND"] * H))   # e5 は上流で一様、フィレットが消える位置 x_r0 までに h1 へ
        huni = elen["e5"] / (NY - 1); h5 = (1 - w) * huni + w * hend      # 間隔で線形に混ぜる (対数で混ぜると節点の移動が上流側に偏り、壁面の格子線が x-y 面で 45° 傾く)
        def bl(n, Led, h):
            cf = bump.coef(n, h / Led); return (n, "Progression", 1.0) if cf is None else (n, "Bump", cf)
        L["e5"] = bl(NY, elen["e5"], h5)
        for e in ("e7", "e11"): L[e] = bl(NY, elen[e], max(mc, h5 * elen[e] / elen["e5"]))
        for k in (2, 3, 4): L["vy3%d" % k] = bl(NY, elen["vy32"], hend)
        for j, hh in ((0, m2), (1, h1), (2, h1), (3, h1), (4, h1)): L["hz%d1" % j] = bl(NSW, G.TSW, hh)
        cap = lambda hh: min(hh, 0.999 * DR / NL)
        for k, hh in ((0, h1), (1, h1), (2, h1), (3, h1), (4, m2)):
            L["vy2%d" % k] = bl(NTC, G.TC, hh)
            hh = cap(hh); r_ = prog_r(DR, NL, hh); L["vy1%d" % k] = (NL + 1, "Progression", -r_)
            L["vy0%d" % k] = (NFY + 1, "Progression", -prog_r(G.YCL - DR - G.YBOT, NFY, hh * r_ ** NL))
        for j, hh in ((0, m2), (1, h1), (2, h1), (3, h1), (4, h1)):
            hh = cap(hh); r_ = prog_r(DR, NL, hh); L["hz%d2" % j] = (NL + 1, "Progression", r_)
            L["hz%d3" % j] = (NFZ + 1, "Progression", prog_r(G.ZFAR - G.ZO - DR, NFZ, hh * r_ ** NL))
        return L

    # ---- 断面の前計算 (gmsh の較正モデルは本体モデルを作る前に使い切る)
    x_r0 = min([xb for (xa, ua, la), (xb, ub, lb) in zip(G.prof[:-1], G.prof[1:]) if ub == 0 and lb == 0 and all(q[1] == 0 and q[2] == 0 for q in G.prof if q[0] >= xb)] or [G.LSW])
    need_blk = lambda si: (PRESENT[ivs[si - 1][2]] if si > 0 else set()) | (PRESENT[ivs[si][2]] if si < len(ivs) else set())
    SEC = []
    for si, x in enumerate(xs):
        pt, walls, _ = G.section(x); nb = need_blk(si)
        ed = sorted({e for b in nb for e, _ in BLK[b]}); elen, wr = {}, {}
        for e in ed:
            a_, b_ = EDGE[e]
            if e in walls:
                W = np.asarray(walls[e]); dd = np.r_[0, np.cumsum(np.linalg.norm(np.diff(W, axis=0), axis=1))]
                W = W[np.r_[True, np.diff(dd) > 1e-14]]; dd = np.r_[0, np.cumsum(np.linalg.norm(np.diff(W, axis=0), axis=1))]
                sv = np.linspace(0, dd[-1], int(P["M_SPL"])); wr[e] = np.stack([np.interp(sv, dd, W[:, q]) for q in (0, 1)], 1); elen[e] = dd[-1]
            else: elen[e] = math.dist(pt[a_], pt[b_])
        for e in ("e5", "e7", "e11", "vy32"): elen.setdefault(e, math.dist(pt[EDGE[e][0]], pt[EDGE[e][1]]))
        SEC.append(dict(pt=pt, nb=nb, ed=ed, elen=elen, wr=wr, law=laws(x, pt, elen)))
    # ---- gmsh モデル
    gmsh.model.add("junction"); ge = gmsh.model.geo; gm = ge.mesh
    S = []
    for si, x in enumerate(xs):
        Q = SEC[si]; pt, ed = Q["pt"], Q["ed"]; pn = sorted({p for e in ed for p in EDGE[e]})
        pid = {p: ge.addPoint(x, pt[p][0], pt[p][1]) for p in pn}; eid = {}
        for e in ed:
            a_, b_ = EDGE[e]
            if e in Q["wr"]: eid[e] = ge.addSpline([pid[a_]] + [ge.addPoint(x, p_[0], p_[1]) for p_ in Q["wr"][e][1:-1]] + [pid[b_]])
            else: eid[e] = ge.addLine(pid[a_], pid[b_])
        face = {b: ge.addSurfaceFilling([ge.addCurveLoop([sg * eid[e] for e, sg in BLK[b]])]) for b in sorted(Q["nb"])}
        S.append(dict(x=x, pt=pt, pid=pid, eid=eid, face=face, law=Q["law"]))
    vols, xfaces, refs = {}, {}, {}
    for ii, (xa, xb, rg) in enumerate(ivs):
        A, B = S[ii], S[ii + 1]; pm = G.section(0.5 * (xa + xb))[0]; blks = sorted(PRESENT[rg])
        ed = sorted({e for b in blks for e, _ in BLK[b]}); pn = sorted({p for e in ed for p in EDGE[e]}); xl = {}
        for p in pn:
            ya, yb, ym = np.array(A["pt"][p]), np.array(B["pt"][p]), np.array(pm[p])
            if np.abs(0.5 * (ya + yb) - ym).max() < 1e-13: xl[p] = ge.addLine(A["pid"][p], B["pid"][p])
            else: xl[p] = ge.addSpline([A["pid"][p], ge.addPoint(0.5 * (xa + xb), ym[0], ym[1]), B["pid"][p]])
            n, r = xlaw[ii]; gm.setTransfiniteCurve(xl[p], n + 1, "Progression", r)
        for e in ed:
            a, b = EDGE[e]
            xfaces[(ii, e)] = ge.addSurfaceFilling([ge.addCurveLoop([A["eid"][e], xl[b], -B["eid"][e], -xl[a]])])
        for b in blks:
            fs = [("s", ii, b), ("s", ii + 1, b)] + [("x", ii, e) for e, _ in BLK[b]]
            for f in fs: refs[f] = refs.get(f, 0) + 1
            sl = ge.addSurfaceLoop([A["face"][b], B["face"][b]] + [xfaces[(ii, e)] for e, _ in BLK[b]]); vols[(ii, b)] = ge.addVolume([sl])
    for St in S:
        for e, tg in St["eid"].items():
            n, typ, cf = St["law"][e]; gm.setTransfiniteCurve(tg, n, typ, cf)
        for f in St["face"].values(): gm.setTransfiniteSurface(f); gm.setRecombine(2, f)
    for f in xfaces.values(): gm.setTransfiniteSurface(f); gm.setRecombine(2, f)
    for v in vols.values(): gm.setTransfiniteVolume(v)
    # ---- 境界タグ: 参照数 1 の面だけ
    groups, untag = {}, []
    for f, c in refs.items():
        assert c in (1, 2), (f, c)
        if c == 2: continue
        if f[0] == "s":
            _, si, b = f; x = xs[si]
            if si == 0: t = "inlet" if b.startswith("N_") else "ext_in"
            elif si == len(xs) - 1: t = "outlet"
            elif b == "SW" and abs(x - G.LSW) < 1e-12: t = "sidewall_end"
            elif b in ("CW1", "CW2") and abs(x - G.LCOWL) < 1e-12: t = "cowl_base"
            else: t = None
            tg = S[si]["face"][b]
        else:
            _, ii, e = f; t = xface_tag(e, ivs[ii][2]); tg = xfaces[(ii, e)]
        if t is None: untag.append(f)
        else: groups.setdefault(t, []).append(tg)
    # 表にあるのに共有面だった、の逆向きの取りこぼしも見る (参照数 2 の面にタグは付かない)
    ge.synchronize()
    for t in sorted(groups): gmsh.model.addPhysicalGroup(2, groups[t], name=t)
    gmsh.model.addPhysicalGroup(3, list(vols.values()), name="fluid")
    if P.get("DIM2"):                       # 断面だけ切って 2D の歪みをブロック別・station 別に見る (分布の調整用。速い)
        gmsh.model.mesh.generate(2); ntag, xyz, _ = gmsh.model.mesh.getNodes(); xyz = xyz.reshape(-1, 3)
        idx = np.zeros(int(ntag.max()) + 1, np.int64); idx[ntag.astype(np.int64)] = np.arange(len(ntag)); res = {}
        for si, St in enumerate(S):
            for b, f in St["face"].items():
                q = xyz[idx[np.asarray(gmsh.model.mesh.getElements(2, f)[2][0], np.int64)].reshape(-1, 4)]; sk = np.zeros(len(q))
                for c_ in range(4):
                    u, v = q[:, (c_ - 1) % 4] - q[:, c_], q[:, (c_ + 1) % 4] - q[:, c_]
                    th = np.degrees(np.arccos(np.clip(np.einsum('ij,ij->i', u, v) / (np.linalg.norm(u, axis=1) * np.linalg.norm(v, axis=1)), -1, 1)))
                    sk = np.maximum(sk, np.abs(th - 90) / 90)
                res.setdefault(b, []).append((round(xs[si] / H, 3), round(float(sk.max()), 3), int((sk > 0.7).sum())))
        for b in res: print("%-7s" % b, " ".join("%.3f:%.2f(%d)" % t for t in res[b] if t[1] > 0.3))
        print("counts", dict(NL=NL, NY=NY, NZ=NZ, NSW=NSW, NTC=NTC, NFY=NFY, NFZ=NFZ), "nx", [xlaw[i][0] for i in range(len(ivs))], "sum", sum(xlaw[i][0] for i in range(len(ivs))))
        gmsh.finalize(); return dict(VERDICT="DIM2")
    gmsh.model.mesh.generate(3)
    rep = dict(params={k: P[k] for k in sorted(P)}, stations=[x / H for x in xs], regions=[v[2] for v in ivs], nx=[abs(xlaw[i][0]) for i in range(len(ivs))],
               x_ratio=[abs(xlaw[i][1]) for i in range(len(ivs))], x_adjacent_ratio_max=x_ratio_max, dx_min=float(dxs.min()), dx_max=float(dxs.max()), counts=dict(NL=NL, NY=NY, NZ=NZ, NSW=NSW, NTC=NTC, NFY=NFY, NFZ=NFZ),
               n_blocks=len(vols), blocks_per_region={r_: len(PRESENT[r_]) for r_ in "ABC"}, untagged_boundary_faces=[str(u) for u in untag],
               tags={t: len(v) for t, v in groups.items()})
    rep.update(check(G, P, vols, groups, h1, h1e, out))
    ok = (not untag and x_ratio_max <= g * 1.02 and rep["ar_gt5000"] == 0 and rep["ar1000_skew_max"] <= 0.30 and rep["other_elems"] == {} and rep["neg_jac_f32"] == 0 and rep["mixed_sign_volumes"] == 0 and rep["skew_max"] <= 0.90
          and rep["topology"]["ok"] and all(v["ok"] for v in rep["first_layer"].values()))
    rep["VERDICT"] = "PASS" if ok else "FAIL"
    json.dump(rep, open(out + "_report.json", "w"), indent=1, ensure_ascii=False)
    if not quiet: print(json.dumps({k: v for k, v in rep.items() if k != "params"}, indent=1, ensure_ascii=False))
    gmsh.write(out + ".msh"); gmsh.finalize(); return rep


# ---------------------------------------------------------------- 検査
NBR = {0: (1, 3, 4), 1: (2, 0, 5), 2: (3, 1, 6), 3: (0, 2, 7), 4: (7, 5, 0), 5: (4, 6, 1), 6: (5, 7, 2), 7: (6, 4, 3)}
HF = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
HE = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7)]


def cornerJ(Pn):
    return np.stack([np.einsum('ij,ij->i', np.cross(Pn[:, a] - Pn[:, c], Pn[:, b] - Pn[:, c]), Pn[:, d] - Pn[:, c]) for c, (a, b, d) in NBR.items()], 1)


def check(G, P, vols, groups, h1, h1e, out):
    from scipy.spatial import cKDTree
    M = gmsh.model.mesh; ntag, xyz, _ = M.getNodes(); xyz = xyz.reshape(-1, 3)
    idx = np.zeros(int(ntag.max()) + 1, np.int64); idx[ntag.astype(np.int64)] = np.arange(len(ntag))
    other, hx_l, mixed = {}, [], 0
    for v in vols.values():
        et, _, cn = M.getElements(3, v)
        for t, c in zip(et, cn):
            if int(t) != 5: other[int(t)] = other.get(int(t), 0) + len(c); continue
            hv = idx[np.asarray(c, np.int64)].reshape(-1, 8); hx_l.append(hv)
    et2, _, cn2 = M.getElements(2)
    for t, c in zip(et2, cn2):
        if int(t) != 3: other[int(t)] = other.get(int(t), 0) + len(c)
    hx = np.vstack(hx_l); used = np.unique(hx); x32 = xyz.astype(np.float32).astype(float)
    J = np.vstack([cornerJ(xyz[hx[i:i + 400000]]) for i in range(0, len(hx), 400000)])
    J32 = np.vstack([cornerJ(x32[hx[i:i + 400000]]) for i in range(0, len(hx), 400000)])
    sg = np.sign(np.median(J)); o = 0
    for hv in hx_l:
        s_ = np.sign(J[o:o + len(hv)]); mixed += int(not (np.all(s_ > 0) or np.all(s_ < 0)) or np.sign(s_.flat[0]) != sg); o += len(hv)
    sk = np.zeros(len(hx))
    for f in HF:
        Pf = [xyz[hx[:, q]] for q in f]
        for c in range(4):
            a, b = Pf[(c - 1) % 4] - Pf[c], Pf[(c + 1) % 4] - Pf[c]
            th = np.degrees(np.arccos(np.clip(np.einsum('ij,ij->i', a, b) / (np.linalg.norm(a, axis=1) * np.linalg.norm(b, axis=1)), -1, 1)))
            sk = np.maximum(sk, np.abs(th - 90) / 90)
    o = 0; per = {}
    Led = np.stack([np.linalg.norm(xyz[hx[:, a]] - xyz[hx[:, b]], axis=1) for a, b in HE], 1); AR = Led.max(1) / Led.min(1)
    for (ii, b), hv in zip(vols.keys(), hx_l):
        q = per.setdefault(b, dict(skew_max=0.0, ar_max=0.0, n=0)); m = slice(o, o + len(hv)); o += len(hv)
        if sk[m].max() > q["skew_max"]: q["skew_max"] = float(sk[m].max()); q["skew_at_interval"] = ii
        if AR[m].max() > q["ar_max"]: q["ar_max"] = float(AR[m].max()); q["ar_at_interval"] = ii
        q["n"] += len(hv)
        mm = (AR[m] > 1000) & (sk[m] > 0.30)
        if mm.any():
            q["ar1000_skewed"] = q.get("ar1000_skewed", 0) + int(mm.sum()); i_ = int(np.flatnonzero(mm)[np.argmax(AR[m][mm])])
            cand = (float(AR[m][i_]), float(sk[m][i_]), [float(v) for v in xyz[hv[i_]].mean(0) / G.H], [float(v) for v in np.ptp(xyz[hv[i_]], axis=0)])
            if cand[0] > q.get("ar1000_skewed_worst", (0,))[0]: q["ar1000_skewed_worst"] = cand
    # ---- 位相: 全ヘキサの面を数え、1 回だけ現れる面 = 境界面がタグ付き四角形とちょうど一致すること
    def keys(Q):                                   # 面 (4 節点) -> 並べ替えに使う 2 本の 64 bit キー
        Q = np.sort(Q.astype(np.int64), axis=1); return (Q[:, 0] << 32) | Q[:, 1], (Q[:, 2] << 32) | Q[:, 3]
    h32 = hx.astype(np.int32); k1, k2 = keys(np.vstack([h32[:, f] for f in HF])); od = np.lexsort((k2, k1)); k1, k2 = k1[od], k2[od]; del od
    new = np.r_[True, (k1[1:] != k1[:-1]) | (k2[1:] != k2[:-1])]; st = np.flatnonzero(new); cnt = np.diff(np.r_[st, len(k1)])
    b1, b2 = k1[st[cnt == 1]], k2[st[cnt == 1]]; del k1, k2
    bq = {}
    for t, surfs in groups.items():
        q = []
        for s_ in surfs:
            _, _, c2 = M.getElements(2, s_); q.append(idx[np.asarray(c2[0], np.int64)].reshape(-1, 4))
        bq[t] = np.vstack(q)
    a1, a2 = keys(np.vstack(list(bq.values()))); od = np.lexsort((a2, a1)); a1, a2 = a1[od], a2[od]
    dup = int(((a1[1:] == a1[:-1]) & (a2[1:] == a2[:-1])).sum())
    same = len(b1) == len(a1) and bool(np.all(b1 == a1) and np.all(b2 == a2))
    topo = dict(faces_gt2=int((cnt > 2).sum()), boundary_faces=int(len(b1)), tagged_quads=int(len(a1)), quads_in_two_tags=dup,
                unused_nodes=int(len(xyz) - len(used)), ok=bool(same and (cnt > 2).sum() == 0 and dup == 0))
    # ---- 壁第一層: 壁タグ別・節点別の |n . dx| / h1
    wallnode = np.zeros(len(xyz), bool)
    for t in WALL_TAGS:
        if t in bq: wallnode[bq[t].ravel()] = True
    cand = np.where(wallnode[hx].sum(1) >= 4)[0]; key = {}
    for hi in cand:
        for f in HF:
            if wallnode[hx[hi, list(f)]].all(): key[tuple(sorted(hx[hi, list(f)]))] = (hi, f)
    multi = np.zeros(len(xyz), int)
    for t in WALL_TAGS:
        if t in bq: multi[np.unique(bq[t])] += 1
    ridge = xyz[multi >= 2]; tree = cKDTree(ridge) if len(ridge) else None
    fl, acc = {}, {}; nrm_all = np.zeros((len(xyz), 3))
    for t in WALL_TAGS:                      # 1 巡目: 壁タグ別の節点法線 (隣接する壁四角形の面積ベクトルの和) と対向節点
        if t not in bq: continue
        q = bq[t]; nrm = np.zeros((len(xyz), 3)); opp = -np.ones(len(xyz), np.int64); uneval = 0
        for row in q:
            hit = key.get(tuple(sorted(row)))
            if hit is None: uneval += 1; continue
            hi, f = hit; loc = list(f); Pq = xyz[hx[hi, loc]]; n_ = np.cross(Pq[2] - Pq[0], Pq[3] - Pq[1])
            cen_in = xyz[hx[hi]].mean(0) - Pq.mean(0); n_ = n_ * np.sign(n_ @ cen_in)
            for l in loc:
                o_ = [m for m in NBR[l] if m not in loc][0]; nrm[hx[hi, l]] += n_; opp[hx[hi, l]] = hx[hi, o_]
        acc[t] = (nrm, opp, uneval)
    for t in WALL_TAGS:
        if t not in bq: continue
        q = bq[t]; nrm, opp, uneval = acc[t]
        nd = np.unique(q); nd = nd[opp[nd] >= 0]; n_ = nrm[nd] / np.linalg.norm(nrm[nd], axis=1)[:, None]
        ns = nrm[nd].copy()                                                   # タグの境が滑らか (フィレットの 45° 点など) なら隣のタグの面も法線に入れる。稜 (50° 超) は入れない
        for t2 in acc:
            if t2 == t: continue
            v2 = acc[t2][0][nd]; l2 = np.linalg.norm(v2, axis=1); okk = l2 > 0
            cs = np.zeros(len(nd)); cs[okk] = np.einsum('ij,ij->i', n_[okk], v2[okk]) / l2[okk]
            ns += np.where((cs > math.cos(math.radians(50.0)))[:, None], v2, 0.0)
        n_ = ns / np.linalg.norm(ns, axis=1)[:, None]
        target = h1e if t in ("sidewall_end", "cowl_base") else h1
        ratio = np.abs(np.einsum('ij,ij->i', n_, xyz[opp[nd]] - xyz[nd])) / target
        # 区分: (i) フィレット円弧の上 (+ 余白 dr/4) は [-10 %, +45 %] (バタフライ位相では 45° 点で最大 sqrt2 倍厚い。薄くはしない)、
        #       (ii) 稜 (複数の壁タグに属する節点) から 2 dr 以内は ±10 %、(iii) それ以外の平面部は ±5 %
        arc = np.zeros(len(nd), bool); edge = np.zeros(len(nd), bool)
        if tree is not None: edge |= tree.query(xyz[nd])[0] < 2 * G.DR
        if t in ("ramp", "sidewall_in", "cowl_in"):
            dl = 0.05 * G.H      # r -> 0 の端点 (フィレットの生え際) も円弧の区分に入れる
            rr = np.array([np.max([G.r(x_ + q) for q in (-dl, 0.0, dl)], axis=0) for x_ in xyz[nd, 0]]); yrn = np.array([G.yr(x_) for x_ in xyz[nd, 0]]); mg = 0.5 * G.DR
            arc |= (rr[:, 0] > 0) & (xyz[nd, 2] > G.ZW - rr[:, 0] - mg) & (xyz[nd, 1] > yrn - rr[:, 0] - mg)
            arc |= (rr[:, 1] > 0) & (xyz[nd, 2] > G.ZW - rr[:, 1] - mg) & (xyz[nd, 1] < G.YC + rr[:, 1] + mg)
            for (xa, ua, la), (xb, ub, lb) in zip(G.prof[:-1], G.prof[1:]):       # 半径が変わる区間 (フィレット遷移) はダクト内の壁全周を ±10 % の区分に入れる (plan §4.10 (e))
                if (ua, la) != (ub, lb): edge |= (xyz[nd, 0] >= xa - 1e-12) & (xyz[nd, 0] <= xb + 1e-12)
        edge &= ~arc; core = ~(edge | arc); dev = ratio - 1
        bad = (core & (np.abs(dev) > 0.05)) | (edge & (np.abs(dev) > 0.10)) | (arc & ((dev < -0.10) | (dev > 0.45)))
        fl[t] = dict(nodes=int(len(nd)), unevaluable_quads=int(uneval), ratio_min=float(ratio.min()), ratio_max=float(ratio.max()))
        for nm, mk in (("core", core), ("edge", edge), ("arc", arc)):
            fl[t][nm] = dict(nodes=int(mk.sum()), out=int((bad & mk).sum()), min=float(ratio[mk].min()) if mk.any() else 1.0, max=float(ratio[mk].max()) if mk.any() else 1.0)
        fl[t]["ok"] = bool(uneval == 0 and not bad.any()); dev = np.abs(dev)
        fl[t]["_bad"] = [[float(v) for v in xyz[nd[i]] / G.H] + [float(ratio[i]), "core" if core[i] else ("edge" if edge[i] else "arc")] for i in np.flatnonzero(bad)[:8]]
        fl[t]["_worst"] = [[float(v) for v in xyz[nd[i]]] + [float(ratio[i])] for i in np.argsort(-dev)[:5]]
    # ---- 描画用
    dump = {"wall_" + t: xyz[bq[t]] for t in bq if t in WALL_TAGS}
    for xq in P.get("DUMP_X", [0.0, 0.6, 1.2, 1.4, 1.6, 2.0]):
        xv = xq * G.H; tol = 1e-9
        for f in HF:
            Pf = xyz[hx[:, f]]; m = np.all(np.abs(Pf[:, :, 0] - xv) < tol, axis=1)
            if m.any(): dump.setdefault("sec_%.3f" % xq, []).append(Pf[m])
    dump = {k: (np.vstack(v) if isinstance(v, list) else v) for k, v in dump.items()}
    np.savez_compressed(out + "_sections.npz", **dump)
    return dict(nodes=int(len(used)), hexes=int(len(hx)), other_elems=other, neg_jac=int((J * sg <= 0).sum()), neg_jac_f32=int((J32 * sg <= 0).sum()),
                mixed_sign_volumes=mixed, skew_max=float(sk.max()), skew_p99=float(np.percentile(sk, 99)), skew_gt090=int((sk > 0.9).sum()),
                skew_gt070=int((sk > 0.7).sum()), ar_max=float(AR.max()), ar_gt1000=int((AR > 1000).sum()), ar_gt5000=int((AR > 5000).sum()), ar1000_skew_max=float(sk[AR > 1000].max()) if (AR > 1000).any() else 0.0,
                min_edge=float(Led.min()), per_block=per, topology=topo, first_layer=fl)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("out"); ap.add_argument("--scale", type=float, default=1.0)
    ap.add_argument("--set", nargs="*", default=[]); a = ap.parse_args(); P = dict(P0)
    for kv in a.set: k, v = kv.split("="); P[k] = float(v)
    if a.scale != 1.0:          # 真の細分列: h1 を 1/scale、成長率を 1/scale 乗、接線分割を scale 倍
        P["H1"] /= a.scale; P["G"] = P["G"] ** (1.0 / a.scale); P["NZ"] = int(round((P["NZ"] - 1) * a.scale)) + 1
        P["HX"] /= a.scale; P["HX_FAR"] /= a.scale
    r = build(a.out, P, quiet=False); sys.exit(0 if r["VERDICT"] == "PASS" else 1)
