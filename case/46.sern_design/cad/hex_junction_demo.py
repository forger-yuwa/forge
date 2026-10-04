#!/usr/bin/env python3
"""お試し形状 (全ヘキサ): SERN ダクトの半スパン断面をバタフライ型 C リング (上・側・下) + コアで切り、
x 方向に transfinite 体積でつなぐ。隅フィレットは上 (ランプ×側壁) と下 (カウル×側壁) で**独立な半径プロファイル r(x)**。
壁断面曲線は**弧長一様に標本化した 1 本のスプライン** (r → 0 でも位相が変わらない。plan tooling-sern-mesh-blocking §4.7)。

usage (mesh venv):  .venv-mesh/bin/python case/46.sern_design/cad/hex_junction_demo.py OUT_PREFIX
出力: OUT_PREFIX.msh / OUT_PREFIX.vtk / OUT_PREFIX_sections.npz (断面メッシュの描画用)
"""
import sys, math, json
import numpy as np
import gmsh

H = 0.1                      # [m]
ZW = 1.0 * H                 # 側壁内面 z = W/2
DR = 0.08 * H                # リング厚 (物理入力)
H1 = 4.0e-5 * H              # 壁第一層 = 4 µm
NL, NZ, NY = 40, 41, 41      # 壁法線層数 / z 方向 / y 方向の節点数
M_SPL = 129                  # 壁断面スプラインの標本点数
ramp = lambda x: H * (1.0 + 0.30 * (x / H) - 0.02 * (x / H) ** 2)      # ランプ輪郭 y_r(x)
cowl = lambda x: 0.0                                                   # カウル上面

# (x/H, 上隅 r/H, 下隅 r/H)。上 = 徐変で始まり、終端の手前 1 r_f で消す。下 = 入口から一定、同じく終端で消す
PROFILE = [(0.00, 0.00, 0.06), (0.15, 0.10, 0.06), (0.60, 0.10, 0.06), (1.00, 0.10, 0.06), (1.10, 0.00, 0.00), (1.20, 0.00, 0.00)]
NX_PER_H = 40                # x 方向の分割 (1 H あたり)


def prog(L, n, h):
    lo, hi = 1.0 + 1e-9, 3.0
    for _ in range(200):
        r = 0.5 * (lo + hi)
        lo, hi = (r, hi) if h * (r ** n - 1) / (r - 1) < L else (lo, r)
    return 0.5 * (lo + hi)


def resample(P, m):
    """折れ線 P を弧長一様に m 点へ"""
    P = np.asarray(P); d = np.r_[0, np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))]
    s = np.linspace(0, d[-1], m)
    return np.stack([np.interp(s, d, P[:, k]) for k in range(P.shape[1])], 1)


def arc(cy, cz, r, a0, a1, n=200):
    a = np.linspace(a0, a1, n); return np.stack([cy + r * np.sin(a), cz + r * np.cos(a)], 1)


def section(x, ru, rl):
    """断面の 8 頂点 (y,z) と 3 本の壁曲線 (密な折れ線)。角度 a は z 軸から y 軸へ測る"""
    yr, yc = ramp(x), cowl(x)
    Su, Sl = (yr, ZW), (yc, ZW)
    # 上隅: 中心 (yr-ru, ZW-ru)。ランプ接点 a=90°、中点 45°、側壁接点 0°
    if ru > 0:
        up1 = arc(yr - ru, ZW - ru, ru, math.pi / 2, math.pi / 4); up2 = arc(yr - ru, ZW - ru, ru, math.pi / 4, 0.0)
    else:
        up1 = np.array([Su]); up2 = np.array([Su])
    # 下隅: 中心 (yc+rl, ZW-rl)。側壁接点 a=0°、中点 -45°、カウル接点 -90°
    if rl > 0:
        lo1 = arc(yc + rl, ZW - rl, rl, 0.0, -math.pi / 4); lo2 = arc(yc + rl, ZW - rl, rl, -math.pi / 4, -math.pi / 2)
    else:
        lo1 = np.array([Sl]); lo2 = np.array([Sl])
    WU, WL = up1[-1], lo1[-1]
    wall_top = np.vstack([[(yr, 0.0)], [(yr, ZW - ru)], up1])
    wall_side = np.vstack([up2, [(yr - ru, ZW)], [(yc + rl, ZW)], lo1])
    wall_bot = np.vstack([lo2, [(yc, ZW - rl)], [(yc, 0.0)]])
    pts = dict(PT0=(yr, 0.0), WU=tuple(WU), CUR=(yr - DR, ZW - DR), CUL=(yr - DR, 0.0),
               WL=tuple(WL), CLR=(yc + DR, ZW - DR), CLL=(yc + DR, 0.0), PB0=(yc, 0.0))
    return pts, dict(top=wall_top, side=wall_side, bot=wall_bot)


def main(out):
    gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0)
    g = gmsh.model.geo; gm = g.mesh
    names = ["PT0", "WU", "CUR", "CUL", "WL", "CLR", "CLL", "PB0"]
    # 断面辺: (名前, 始点, 終点, 種別)。種別 wall = スプライン
    EDG = [("e1", "PT0", "WU", "wall:top"), ("e2", "WU", "CUR", "norm+"), ("e3", "CUR", "CUL", "z"), ("e4", "CUL", "PT0", "norm-"),
           ("e5", "WU", "WL", "wall:side"), ("e6", "WL", "CLR", "norm+"), ("e7", "CLR", "CUR", "y"),
           ("e8", "WL", "PB0", "wall:bot"), ("e9", "PB0", "CLL", "norm+"), ("e10", "CLL", "CLR", "z"), ("e11", "CUL", "CLL", "y")]
    BLK = {"TOP": ["e1", "e2", "e3", "e4"], "SIDE": ["e5", "e6", "e7", "-e2"], "BOT": ["e8", "e9", "e10", "-e6"],
           "CORE": ["-e3", "-e7", "-e10", "-e11"]}
    secs = []
    for xh, ru, rl in PROFILE:
        x = xh * H; pts, walls = section(x, ru * H, rl * H)
        pid = {n: g.addPoint(x, pts[n][0], pts[n][1]) for n in names}
        eid, elen = {}, {}
        for nm, a, b, kind in EDG:
            if kind.startswith("wall"):
                W = resample(walls[kind.split(":")[1]], M_SPL)
                mid = [g.addPoint(x, p[0], p[1]) for p in W[1:-1]]
                eid[nm] = g.addSpline([pid[a]] + mid + [pid[b]])
            else:
                eid[nm] = g.addLine(pid[a], pid[b])
            elen[nm] = math.dist(pts[a], pts[b])
        face = {}
        for bn, loop in BLK.items():
            cl = g.addCurveLoop([(-eid[e[1:]] if e.startswith("-") else eid[e]) for e in loop])
            face[bn] = g.addSurfaceFilling([cl])
        secs.append(dict(x=x, pid=pid, eid=eid, elen=elen, face=face, pts=pts))
    vols, wall_faces, sec_faces = [], [], []
    for i in range(len(secs) - 1):
        A, B = secs[i], secs[i + 1]; xm = 0.5 * (A["x"] + B["x"])
        rum = 0.5 * (PROFILE[i][1] + PROFILE[i + 1][1]) * H; rlm = 0.5 * (PROFILE[i][2] + PROFILE[i + 1][2]) * H
        pm, _ = section(xm, rum, rlm)
        xl = {n: g.addSpline([A["pid"][n], g.addPoint(xm, pm[n][0], pm[n][1]), B["pid"][n]]) for n in names}
        xf = {}
        for nm, a, b, kind in EDG:
            cl = g.addCurveLoop([A["eid"][nm], xl[b], -B["eid"][nm], -xl[a]])
            xf[nm] = g.addSurfaceFilling([cl])
            if kind.startswith("wall"): wall_faces.append(xf[nm])
        nx = max(int(round(NX_PER_H * (B["x"] - A["x"]) / H)), 4) + 1
        for n in names: gm.setTransfiniteCurve(xl[n], nx)
        for bn, loop in BLK.items():
            sl = g.addSurfaceLoop([A["face"][bn], B["face"][bn]] + [xf[e.lstrip("-")] for e in loop])
            vols.append(g.addVolume([sl]))
        for s_ in xf.values(): gm.setTransfiniteSurface(s_); gm.setRecombine(2, s_)
    for S in secs:
        for nm, a, b, kind in EDG:
            if kind.startswith("norm"):
                r = prog(S["elen"][nm], NL, H1)
                gm.setTransfiniteCurve(S["eid"][nm], NL + 1, "Progression", r if kind == "norm+" else -r)
            else:
                gm.setTransfiniteCurve(S["eid"][nm], NZ if kind in ("z", "wall:top", "wall:bot") else NY)
        for f in S["face"].values(): gm.setTransfiniteSurface(f); gm.setRecombine(2, f)
    for v in vols: gm.setTransfiniteVolume(v)
    g.synchronize(); gmsh.model.mesh.generate(3); gmsh.model.mesh.removeDuplicateNodes()
    # ---- 検査 ----
    et, tags, conn = gmsh.model.mesh.getElements(3)
    ntag, xyz, _ = gmsh.model.mesh.getNodes(); xyz = xyz.reshape(-1, 3); idx = np.zeros(int(ntag.max()) + 1, int); idx[ntag.astype(int)] = np.arange(len(ntag))
    cnt = {int(t): len(tg) for t, tg in zip(et, tags)}
    hx = idx[np.asarray(conn[list(et).index(5)], int)].reshape(-1, 8)
    nb = {0: (1, 3, 4), 1: (2, 0, 5), 2: (3, 1, 6), 3: (0, 2, 7), 4: (7, 5, 0), 5: (4, 6, 1), 6: (5, 7, 2), 7: (6, 4, 3)}
    def cornerJ(P): return np.stack([np.linalg.det(np.stack([P[:, a] - P[:, c], P[:, b] - P[:, c], P[:, d] - P[:, c]], 1)) for c, (a, b, d) in nb.items()], 1)
    J = cornerJ(xyz[hx]); sg = np.sign(np.median(J)); J32 = cornerJ(xyz.astype(np.float32).astype(float)[hx])
    HF = [(0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]; sk = np.zeros(len(hx))
    for f in HF:
        Pf = [xyz[hx[:, q]] for q in f]
        for c in range(4):
            a = Pf[(c - 1) % 4] - Pf[c]; b = Pf[(c + 1) % 4] - Pf[c]
            th = np.degrees(np.arccos(np.clip(np.einsum('ij,ij->i', a, b) / (np.linalg.norm(a, axis=1) * np.linalg.norm(b, axis=1)), -1, 1)))
            sk = np.maximum(sk, np.maximum((th - 90) / 90, (90 - th) / 90))
    E = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7)]
    L = np.stack([np.linalg.norm(xyz[hx[:, a]] - xyz[hx[:, b]], axis=1) for a, b in E], 1); AR = L.max(1) / L.min(1)
    rep = dict(nodes=int(len(xyz)), hexes=int(len(hx)), other_elems={k: v for k, v in cnt.items() if k != 5},
               neg_jac=int((J * sg <= 0).sum()), neg_jac_f32=int((J32 * sg <= 0).sum()),
               skew_max=float(sk.max()), skew_p99=float(np.percentile(sk, 99)), skew_mean=float(sk.mean()), skew_gt090=int((sk > 0.9).sum()),
               ar_max=float(AR.max()), ar_p99=float(np.percentile(AR, 99)), ar_gt5000=int((AR > 5000).sum()), min_edge=float(L.min()), h1=H1)
    print(json.dumps(rep, indent=1, ensure_ascii=False))
    json.dump(rep, open(out + "_report.json", "w"), indent=1)
    # ---- 描画用データ: 各断面のメッシュ (quad) と壁面メッシュ ----
    dump = {}
    for k, S in enumerate(secs):
        qs = []
        for f in S["face"].values():
            t2, tg2, c2 = gmsh.model.mesh.getElements(2, f)
            qs.append(idx[np.asarray(c2[0], int)].reshape(-1, 4))
        dump["sec%d" % k] = xyz[np.vstack(qs)]; dump["secx%d" % k] = S["x"]
    wq = []
    for f in wall_faces:
        t2, tg2, c2 = gmsh.model.mesh.getElements(2, f); wq.append(idx[np.asarray(c2[0], int)].reshape(-1, 4))
    dump["wall"] = xyz[np.vstack(wq)]
    np.savez_compressed(out + "_sections.npz", **dump)
    gmsh.option.setNumber("Mesh.SaveAll", 1); gmsh.write(out + ".msh"); gmsh.write(out + ".vtk"); gmsh.finalize()


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "hex_junction_demo")
