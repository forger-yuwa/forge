#!/usr/bin/env python3
"""G0 / G1 / G2 を 1 本の GPU run (1 step) で測る (plan boundary-node-periodic-gradient-fix §6)。

変種 (`python3 g_suite.py <variant> [--scratch DIR]`):
  tgv          TGV 32^3 三重周期 (一様直交。2/4/8 member = 面・辺・角)
  tgv_mirror   同じメッシュを点対称 (x→x0+xL−x, y, z も) に写す = 各 group の root (最小 index) が継ぎ目の反対側へ移る (root 交換)
  tgv_bcswap   bcondConfig の記述順を逆にする (root は index 最小で決まるので不変のはず → tgv とビット比較)
  tgv_repeat   tgv と同一設定の再実行 (atomicAdd 順序による run 間差の床)
  tgv_shift    座標原点を (100,100,100) へ移動 (float32 座標の丸めが差分に乗る)
  jitter32     32^3 の節点ジッタ ±0.2h (継ぎ目の両側で stencil が非対称)
  channel      x・z 周期 + y 両壁 (x 等比格子、壁∩継ぎ目)

場 (ρ=1 で焼く。w = 継ぎ目中心の局所座標、gharness.wrapped_coords):
  LSQ (G0/G2): Ux = 0.5 + w·(0,1,0)、Uy = −0.3 + w·(1,0,0)、Uz = 0.2 + w·(0.3,−0.7,0.5) (斜め)、ρ = 1 (非零定数)
  GG  (G1)   : k = 2 + 0.1 w·(1,2,−1)、ω = 100 + w·(2,−1,3)、ξ (受動トレーサ、化学種 Y と同じ GG カーネル) = 0.5 + 0.02 w·(1,−1,2)
判定: 状態は res_0.h5 (step 1 の残差組立はこの状態で行う)、勾配は res_1.h5 (k/ω は res_0 では未計算)。
"""
import argparse
import os
import shutil
import sys

import h5py
import numpy as np
import yaml

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gharness as G          # noqa: E402
import mkmesh                 # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
SCR_DEFAULT = "/tmp/claude-1000/-home-sano-work-forge/4b0c8643-66fb-4fde-8878-7c8a5104c060/scratchpad"
TGV_SRC = os.path.join(SCR_DEFAULT, "lsqseam_m1", "Taylor-Green.h5")

LSQ_FIELDS = {"Ux": (0.5, (0.0, 1.0, 0.0)), "Uy": (-0.3, (1.0, 0.0, 0.0)), "Uz": (0.2, (0.3, -0.7, 0.5))}
GG_FIELDS = {"k": (2.0, (0.1, 0.2, -0.1), "K", "divide"), "omega": (100.0, (2.0, -1.0, 3.0), "Omega", "divide"),
             "Xi": (0.5, (0.02, -0.02, 0.04), "Xi", "mulinv")}
CLASSES = [(1, "内部"), (2, "面 (2 member)"), (4, "辺 (4 member)"), (8, "角 (8 member)")]


def transform_mesh(h5, mode, offset=None):
    """座標系の変換 (ノード座標・双対重心・面重心・面ベクトル)。mode: 'mirror' | 'shift'。"""
    with h5py.File(h5, "r+") as f:
        X = f["MESH/COORD"][()].reshape(-1, 3).astype(np.float64)
        lo, hi = X.min(0), X.max(0)
        for name in ("MESH/COORD", "CELLS/centCoords", "PLANES/centCoords"):
            a = f[name][()].reshape(-1, 3).astype(np.float64)
            a = (lo + hi) - a if mode == "mirror" else a + np.asarray(offset)
            f[name][...] = a.astype(np.float32).reshape(-1)
        if mode == "mirror":
            s = f["PLANES/surfVect"][()]
            f["PLANES/surfVect"][...] = -s


def bcond_variant(src, mode):
    d = yaml.safe_load(open(src))
    if mode == "mirror":
        for v in d.values():
            if v.get("kind") == "periodic":
                for k in ("dx", "dy", "dz"):
                    v["floats"][k] = -float(v["floats"][k])
    elif mode == "swap":
        # 記述順を逆にする。physID は h5 の BCONDS/<id> (面の集合) に結び付いているので入れ替えられない
        # (入れ替えると並進量が面と合わず setPeriodicPartner が停止する: 2026-09-26 実測)。
        # root は union-find で「group 内の最小 index」に決まる (mesh.cpp buildPeriodicNodeGroups の unite) ので、
        # bcond の記述順・対の向きでは root は交換できない → ビット同一になることを確認する。root の実交換は tgv_mirror。
        d = dict(list(d.items())[::-1])
    return d


def prepare(variant, scr):
    run = os.path.join(scr, "g_harness_" + variant)
    tgv_bc = os.path.join(HERE, "bcondConfig.yaml")
    visc = 0.0
    if variant.startswith("tgv"):
        src, bc = TGV_SRC, tgv_bc          # tgv_repeat は tgv と同一設定 (run 間の atomicAdd 順序の床)
        if variant == "tgv_mirror":
            bc = bcond_variant(tgv_bc, "mirror")
        elif variant == "tgv_bcswap":
            bc = bcond_variant(tgv_bc, "swap")
    elif variant == "jitter32":
        src, bc, q = mkmesh.make_box_h5(os.path.join(scr, "g_harness_mesh", "box32_j"), 32, 0.2)
        print("mesh quality:", q)
    elif variant == "channel":
        src, bc, q = mkmesh.make_channel_h5(os.path.join(scr, "g_harness_mesh", "channel"))
        print("mesh quality:", q)
        visc = 1.8e-5
    else:
        raise SystemExit("unknown variant " + variant)
    cfg = G.base_solver_cfg(sst=True, tracer=True, visc=visc)
    h5 = G.make_run(run, src, cfg, bc)
    if variant == "tgv_mirror":
        transform_mesh(h5, "mirror")
    elif variant == "tgv_shift":
        transform_mesh(h5, "shift", (100.0, 100.0, 100.0))
    bcp = os.path.join(run, "bcondConfig.yaml")
    msh = G.Mesh(h5, bcp)
    w, clean, axes = G.wrapped_coords(msh)
    fields = {}
    for name, (c0, a) in list(LSQ_FIELDS.items()):
        fields[name] = c0 + w @ np.asarray(a)
    for name, (c0, a, _, _) in GG_FIELDS.items():
        fields[name] = c0 + w @ np.asarray(a)
    ro = 1.0
    extra = {}
    if variant == "channel":
        # 場の丸め床 ε·|ρ|/(h·|a|) を 1e-5 より十分下げるため勾配を大きく取る (h_min≈0.028、ρ∈[0.25,1.75])
        ro_a = np.array([0.8, -0.4, 0.6])
        ro = 1.2 + w @ ro_a
        fields["ro"] = ro
        extra["ro_a"] = ro_a
    G.bake(h5, G.prim_state_fields(msh, ro, fields["Ux"], fields["Uy"], fields["Uz"], 101325.0,
                                   k=fields["k"], om=fields["omega"], xi=fields["Xi"]))
    np.savez(os.path.join(run, "baked.npz"), w=w, clean=clean, **fields, **extra)
    return run


def evaluate(variant, run, ref_run=None):
    h5 = os.path.join(run, "mesh.h5")
    msh = G.Mesh(h5, os.path.join(run, "bcondConfig.yaml"))
    bk = np.load(os.path.join(run, "baked.npz"))
    s0 = G.read_res(os.path.join(run, "res_0.h5"), ["Ux", "Uy", "Uz", "ro", "P", "k", "omega", "Xi"])
    gn = []
    for v in list(LSQ_FIELDS) + ["ro", "P"]:
        gn += [f"d{v}d{c}" for c in "xyz"]
    for v in GG_FIELDS.values():
        gn += [f"d{v[2]}d{c}" for c in "xyz"]
    g1 = G.read_res(os.path.join(run, "res_1.h5"), gn)
    grad = lambda v: np.stack([g1[f"d{v}d{c}"] for c in "xyz"], axis=1).astype(np.float64)
    out = []
    P = out.append
    P(f"# {variant}: plan boundary-node-periodic-gradient-fix §6 G0/G1/G2")
    P(f"run: {run}")
    P(G.provenance(run).rstrip())
    P(f"harness revision: {G.git_rev()}")
    P(f"CV {msh.nCells}, 周期 group {len(msh.groups)}, member 数の内訳 " +
      ", ".join(f"{n}:{(msh.nmember == n).sum()}" for n in (1, 2, 4, 8) if (msh.nmember == n).any()) +
      f", 周期対の最大残差 {msh.partner_resid:.2e}, 同値類 (重複 >1) {msh.alpha() is not None and msh.n_cls_merged}")
    V13 = np.cbrt(msh.vmerged64)
    # 状態の確認: res_0 の値 = 焼いた値 (f32)
    unchanged = {}
    for v in list(LSQ_FIELDS) + list(GG_FIELDS) + (["ro"] if "ro_a" in bk.files else []):
        ref32 = bk[v].astype(np.float32)
        unchanged[v] = np.abs(s0[v] - ref32) <= 4 * np.spacing(np.abs(ref32))   # ρ≠1 では原始量 = 保存量/ρ の丸めを許す
        if not unchanged[v].all():
            P(f"  注: res_0 の {v} は {(~unchanged[v]).sum()} 節点で焼いた値と異なる (境界条件の上書き)。G0 はその stencil を除外")
    verdict = {}

    # ---------------- G0 (LSQ、解析値) ----------------
    P("\n## G0: 合併 stencil LSQ の線形場 (解析勾配との差)")
    P("閾値: 非零成分 |g_c−a_c| ≤ 1e-5·|a|、ゼロ成分 |g_c| ≤ 1e-6·max|φ|/h (h = 合併体積^(1/3))。"
      "対象 = 折返し不連続にも境界条件の上書きにも触れない節点")
    P("| 場 | a | 区分 | 節点数 | 非零成分 最大誤差/|a| | ゼロ成分 最大 |g|·h/max|φ| | 判定 |")
    P("| --- | --- | --- | --- | --- | --- | --- |")
    g0ok = True
    for v, (c0, a) in LSQ_FIELDS.items():
        a = np.asarray(a)
        ok_st = unchanged[v].copy()
        bad = ~ok_st[msh.inc_j] | ~ok_st[msh.inc_m]
        badg = np.zeros(msh.nCells, bool); badg[msh.root[msh.inc_m[bad]]] = True
        sel0 = bk["clean"] & ~badg[msh.root] & ok_st
        g = grad(v)
        phimax = np.abs(s0[v][sel0]).max()
        for nm, lab in CLASSES:
            sel = sel0 & (msh.nmember == nm)
            if not sel.any():
                continue
            nz = a != 0
            e_nz = (np.abs(g[sel][:, nz] - a[nz]).max() / np.linalg.norm(a)) if nz.any() else 0.0
            e_z = (np.abs(g[sel][:, ~nz]) * V13[sel][:, None]).max() / phimax if (~nz).any() else 0.0
            ok = e_nz <= 1e-5 and e_z <= 1e-6
            g0ok &= ok
            P(f"| {v} | {tuple(a)} | {lab} | {sel.sum()} | {e_nz:.2e} | {e_z:.2e} | {'ok' if ok else 'NG'} |")
    if "ro_a" in bk.files:          # channel: ρ を線形にして壁∩継ぎ目の節点も G0 に入れる (壁は ρ を上書きしない)
        a = bk["ro_a"]; v = "ro"
        ok_st = unchanged[v]
        bad = ~ok_st[msh.inc_j] | ~ok_st[msh.inc_m]
        badg = np.zeros(msh.nCells, bool); badg[msh.root[msh.inc_m[bad]]] = True
        sel0 = bk["clean"] & ~badg[msh.root] & ok_st
        wall = np.zeros(msh.nCells, bool); wall[msh.pc[msh.bnd_planes, 0]] = True; wall = wall[msh.root] | wall
        g = grad("ro")
        phimax = np.abs(s0["ro"][sel0]).max()
        for nm, lab in CLASSES:
            for wl, wlab in ((False, ""), (True, "∩壁")):
                sel = sel0 & (msh.nmember == nm) & (wall == wl)
                if not sel.any():
                    continue
                e_nz = np.abs(g[sel] - a).max() / np.linalg.norm(a)
                ok = e_nz <= 1e-5
                g0ok &= ok
                flo = G.EPS32 * phimax / (V13[sel].min() * np.linalg.norm(a))
                P(f"| ρ (線形) | {tuple(np.round(a, 3))} | {lab}{wlab} | {sel.sum()} | {e_nz:.2e} | - | {'ok' if ok else 'NG'} "
                  f"(場の丸め ε·max|ρ|/(h·|a|) = {flo:.1e}) |")
        verdict["G0"] = g0ok
        P(f"VERDICT G0 ({variant}): {'PASS' if g0ok else 'FAIL'}")
    else:
      # 非零定数場 ρ=1
      gro = grad("ro")
      e_ro = (np.abs(gro) * V13[:, None]).max() / 1.0
      P(f"| ρ (定数 1) | (0,0,0) | 全節点 | {msh.nCells} | - | {e_ro:.2e} | {'ok' if e_ro <= 1e-6 else 'NG'} |")
      g0ok &= e_ro <= 1e-6
      verdict["G0"] = g0ok
      P(f"VERDICT G0 ({variant}): {'PASS' if g0ok else 'FAIL'}")

    # ---------------- G2 (LSQ、CPU double 参照) ----------------
    P("\n## G2: GPU LSQ と CPU double の合併 stencil LSQ (格納 float32 座標・res_0 の場) の差")
    P("閾値: max_i |∇φ_gpu − ∇φ_ref| ≤ 1e-5·S、S = max_i |∇φ_ref| (全節点)")
    P("S は折返し不連続の節点で大きくなるので、不連続に触れない節点だけの S_clean と差も併記する (判定は両方)")
    P("| 場 | S | 最大差 | 最大差/S | S_clean | 最大差 (clean) / S_clean | 最大差 (ε·max|φ|/h 単位) | 判定 |")
    P("| --- | --- | --- | --- | --- | --- | --- | --- |")
    g2ok = True
    floor_rows = {}
    for v in list(LSQ_FIELDS) + ["P"]:
        ref = G.lsq_merged_ref(msh, s0[v])
        g = grad(v)
        S = np.linalg.norm(ref, axis=1).max()
        dif = np.abs(g - ref).max(axis=1)
        unit = G.EPS32 * np.abs(s0[v]).max() / V13
        e = dif.max()
        cl = bk["clean"]
        Sc = np.linalg.norm(ref[cl], axis=1).max()
        ec = dif[cl].max()
        ok = (e <= 1e-5 * S and ec <= 1e-5 * Sc) if v != "P" else True
        if v != "P":
            g2ok &= ok
        P(f"| {v} | {S:.4g} | {e:.2e} | {e / max(S, 1e-300):.2e} | {Sc:.4g} | {ec / max(Sc, 1e-300):.2e} | "
          f"{(dif / unit).max():.1f} | {'ok' if ok else 'NG'}{' (参考: 定数場)' if v == 'P' else ''} |")
        if v == "Ux":
            floor_rows["dUx (LSQ)"] = (dif, unit)
    P(f"  (打ち切り group 数 {msh.lsq_ndegen})")
    verdict["G2"] = g2ok
    P(f"VERDICT G2 ({variant}): {'PASS' if g2ok else 'FAIL'}")

    # ---------------- G1 (GG) ----------------
    P("\n## G1: node GG (k, ω は ransGradient、ξ は受動種 = 化学種と同じ species_gradient_d)")
    P("G1-a: GPU − CPU float32 再現 ≤ 4ε·max|φ|/h。G1-b: GPU − CPU double ≤ max(1e-5·S, 2·N_max·ε·max|φ|/h) "
      "(継ぎ目で超えたら 2·N_max·n_member)。h = 合併体積^(1/3) (区分内の最小)、N_max = 区分内の節点の面数の最大、"
      "S = 折返し不連続にも境界条件の上書きにも触れない節点での max|∇φ_ref|。CPU 再現は plan §4.2 どおり周期半割面を除外")
    P("診断列 (判定外): 「半割面込み」= 周期半割面も面値 φ[ic0] で積算する CPU float32 再現との差 "
      "(GPU の plane_cells では周期半割面の ic1 が ghost なので `excludePeriodic` が効かない仮説の検査)")
    P("| 場 | 区分 | 節点数 | N_max | G1-a 最大差 [ε·max|φ|/h] | G1-b 最大差 [ε·max|φ|/h] | G1-b 閾値 [同] | 判定 a / b | 半割面込み再現との差 [同] |")
    P("| --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    fc = msh.face_count()
    g1a_ok = g1b_ok = True
    for v, (c0, a, gname, div) in GG_FIELDS.items():
        g = grad(gname)
        c32 = G.gg_merged(msh, s0[v], "f32", div).astype(np.float64)
        c64 = G.gg_merged(msh, s0[v], "f64", div)
        # S: 折返し不連続にも境界条件の上書き (壁の k/ω ピン) にも触れない節点での最大 (閾値を緩めない)
        okv = unchanged[v]
        badv = ~okv[msh.inc_j] | ~okv[msh.inc_m]
        bg = np.zeros(msh.nCells, bool); bg[msh.root[msh.inc_m[badv]]] = True
        S = np.linalg.norm(c64[bk["clean"] & ~bg[msh.root] & okv], axis=1).max()
        phimax = np.abs(s0[v]).max()
        c32p = G.gg_merged(msh, s0[v], "f32", div, include_periodic=True).astype(np.float64)
        da = np.abs(g - c32).max(axis=1); db = np.abs(g - c64).max(axis=1); dp = np.abs(g - c32p).max(axis=1)
        for nm, lab in CLASSES:
            sel = msh.nmember == nm
            if not sel.any():
                continue
            h = V13[sel].min()
            unit = G.EPS32 * phimax / h
            Nmax = int(fc[sel].max())
            ea, eb = da[sel].max() / unit, db[sel].max() / unit
            thr_b = max(1e-5 * S / unit, 2.0 * Nmax)
            note = ""
            if eb > thr_b and nm > 1:
                thr_b = max(1e-5 * S / unit, 2.0 * Nmax * nm); note = " (×n_member)"
            oka, okb = ea <= 4.0, eb <= thr_b
            g1a_ok &= oka; g1b_ok &= okb
            P(f"| {v} | {lab} | {sel.sum()} | {Nmax} | {ea:.2f} | {eb:.2f} | {thr_b:.1f}{note} | "
              f"{'ok' if oka else 'NG'} / {'ok' if okb else 'NG'} | {dp[sel].max() / unit:.2f} |")
        floor_rows[{"k": "dK (GG)", "omega": "dΩ (GG)", "Xi": "dξ (GG、Y の代理)"}[v]] = (db, G.EPS32 * phimax / V13)
    # 対の周期半割面の面ベクトル不一致 |S_a+S_b|/|S| と、定数場で継ぎ目に出る GG 勾配の予測 (半割面込み・double)
    from scipy.spatial import cKDTree
    mis = 0.0
    for pid in sorted(msh.periodic_ids):
        q = int(msh.bcfg[pid]["ints"]["partnerBCID"]); fl = msh.bcfg[pid]["floats"]
        off = np.array([fl.get("dx", 0.0), fl.get("dy", 0.0), fl.get("dz", 0.0)])
        pa, pb = msh.bconds[pid]["iPlanes"], msh.bconds[q]["iPlanes"]
        kk = cKDTree(msh.pcent[pb].astype(np.float64)).query(msh.pcent[pa].astype(np.float64) + off)[1]
        ss = msh.sv[pa].astype(np.float64) + msh.sv[pb[kk]].astype(np.float64)
        mis = max(mis, (np.linalg.norm(ss, axis=1) / np.linalg.norm(msh.sv[pa].astype(np.float64), axis=1)).max())
    cst = G.gg_merged(msh, np.ones(msh.nCells, np.float32), "f64", "divide", include_periodic=True)
    P(f"  対の周期半割面の面ベクトル不一致 max|S_a+S_b|/|S_a| = {mis:.2e}。"
      f"定数場 φ≡1 の GG 勾配 (半割面込み、double): 継ぎ目 max {np.abs(cst[msh.nmember > 1]).max() * V13.min():.2e}·φ/h、"
      f"内部 max {np.abs(cst[msh.nmember == 1]).max() * V13.min():.2e}·φ/h")
    verdict["G1a"], verdict["G1b"] = g1a_ok, g1b_ok
    P(f"VERDICT G1-a ({variant}): {'PASS' if g1a_ok else 'FAIL'}")
    P(f"VERDICT G1-b ({variant}): {'PASS' if g1b_ok else 'FAIL'}")

    # ---------------- 床の表 ----------------
    P("\n## 4 量の床 (GPU − CPU double、ε·max|φ|/h 単位、h = 各節点の合併体積^(1/3))")
    P("| 量 | 内部 最大 | 内部 中央値 | 継ぎ目 (member≥2) 最大 | 継ぎ目 中央値 |")
    P("| --- | --- | --- | --- | --- |")
    seam = msh.nmember > 1
    for nm, (dif, unit) in floor_rows.items():
        r = dif / unit
        P(f"| {nm} | {r[~seam].max():.1f} | {np.median(r[~seam]):.2f} | {r[seam].max():.1f} | {np.median(r[seam]):.2f} |")

    if ref_run is not None:
        P(f"\n## root 不変性: {ref_run} との比較")
        P("atomicAdd の順序で run ごとに ulp 級の差が出る (同一設定の再実行 tgv_repeat で床を測る)。判定: 全勾配配列で差 ≤ 4ε·max|φ|/h")
        P("| 配列 | 差のある節点 | 最大差 [ε·max|φ|/h] |")
        P("| --- | --- | --- |")
        mref = G.Mesh(os.path.join(ref_run, "mesh.h5"), os.path.join(ref_run, "bcondConfig.yaml"))
        rootsame = bool((mref.root == msh.root).all())
        worst = 0.0
        with h5py.File(os.path.join(ref_run, "res_1.h5"), "r") as fa, h5py.File(os.path.join(run, "res_1.h5"), "r") as fb:
            for k in gn:
                a_, b_ = fa["VALUE/" + k][()].astype(np.float64), fb["VALUE/" + k][()].astype(np.float64)
                base = {"K": "k", "Omega": "omega", "Xi": "Xi"}.get(k[1:-2], k[1:-2])
                unit = G.EPS32 * np.abs(s0[base]).max() / V13.min() if base in s0 else 1.0
                e = np.abs(a_ - b_).max() / unit
                worst = max(worst, e)
                P(f"| {k} | {int((a_ != b_).sum())} | {e:.2f} |")
        P(f"  CPU で再構成した root 配列が同一: {rootsame}")
        verdict["rootorder"] = rootsame and worst <= 4.0
        P(f"VERDICT root-order ({variant}): {'PASS' if verdict['rootorder'] else 'FAIL'} (最大 {worst:.2f} ε·max|φ|/h)")
    return "\n".join(out) + "\n", verdict


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("variant")
    ap.add_argument("--scratch", default=SCR_DEFAULT)
    ap.add_argument("--eval-only", action="store_true")
    a = ap.parse_args()
    run = os.path.join(a.scratch, "g_harness_" + a.variant)
    if not a.eval_only:
        run = prepare(a.variant, a.scratch)
        G.run_forge(run)
    ref = os.path.join(a.scratch, "g_harness_tgv") if a.variant in ("tgv_bcswap", "tgv_repeat") else None
    txt, v = evaluate(a.variant, run, ref)
    outp = os.path.join(HERE, f"G_{a.variant}.txt")
    open(outp, "w").write(txt)
    print(txt)
    print("->", outp)


if __name__ == "__main__":
    main()
