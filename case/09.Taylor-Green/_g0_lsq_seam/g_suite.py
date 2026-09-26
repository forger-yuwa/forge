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
G1-b (2026-09-26 codex result M3 → §5.1 #6c): 上限は 2·N_max·ε·max|φ|/h のみ (n_member 倍・1e-5·S との max の自動緩和は削除)、
  継ぎ目/内部の誤差比 ≤ 2 も機械判定に入れる。
GPU 定数場 (同 #6c、2026-09-26 #8c で 3 区分に): 別 run `g_harness_<variant>_const` で k・ω・ξ = const を焼いて 1 step。
  壁なし継ぎ目 ≤ 4ε|φ|/h、壁∩継ぎ目は φ·ΣS_f/V (格納 float32 面ベクトルの double 閉包) との差 ≤ 4ε|φ|/h、非継ぎ目壁は報告。
scratch は `--scratch` または環境変数 G_HARNESS_SCRATCH。tgv 系の入力は <scratch>/lsqseam_m1/Taylor-Green.h5。
"""
import argparse
import os
import shutil
import subprocess
import sys

import h5py
import numpy as np
import yaml

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gharness as G          # noqa: E402
import mkmesh                 # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
SCR_DEFAULT = os.environ.get("G_HARNESS_SCRATCH",
                             "/tmp/claude-1000/-home-sano-work-forge/4b0c8643-66fb-4fde-8878-7c8a5104c060/scratchpad")
TGV_REL = os.path.join("lsqseam_m1", "Taylor-Green.h5")      # scratch からの相対 (AWS でも同じ配置にコピーして使う)

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


# GPU 定数場 (plan §6 G1「定数場: 継ぎ目 ≤ 4ε|φ|/h」、§5.1 #6c): k・ω・ξ を定数で焼いて 1 step
CONST_FIELDS = {"k": (2.0, "K"), "omega": (100.0, "Omega"), "Xi": (0.5, "Xi")}


def prepare(variant, scr, const=False):
    run = os.path.join(scr, "g_harness_" + variant + ("_const" if const else ""))
    TGV_SRC = os.path.join(scr, TGV_REL)
    tgv_bc = os.path.join(HERE, "bcondConfig.yaml")
    visc = 0.0
    if variant.startswith("tgv"):
        src, bc = TGV_SRC, tgv_bc          # tgv_repeat は tgv と同一設定 (run 間の atomicAdd 順序の床)
        if variant == "tgv_mirror":
            bc = bcond_variant(tgv_bc, "mirror")
        elif variant == "tgv_bcswap":
            bc = bcond_variant(tgv_bc, "swap")
    elif variant in ("jitter32", "channel"):
        # 既にあるメッシュは作り直さない (本 run と定数場 run、再実行で同じ入力を使う)
        md = os.path.join(scr, "g_harness_mesh", "box32_j" if variant == "jitter32" else "channel")
        if os.path.exists(os.path.join(md, "mesh.h5")) and os.path.exists(os.path.join(md, "quality.txt")):
            src, bc = os.path.join(md, "mesh.h5"), os.path.join(md, "bcondConfig.yaml")
            q = [l for l in open(os.path.join(md, "quality.txt")).read().splitlines() if "VERDICT" in l][-1]
        elif variant == "jitter32":
            src, bc, q = mkmesh.make_box_h5(md, 32, 0.2)
        else:
            src, bc, q = mkmesh.make_channel_h5(md)
        print("mesh quality:", q)
        if variant == "channel":
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
        fields[name] = (CONST_FIELDS[name][0] + 0.0 * w[:, 0]) if const else (c0 + w @ np.asarray(a))
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
    P("G1-a: GPU − CPU float32 再現 ≤ 4ε·max|φ|/h。G1-b (plan §6 G1、2026-09-26 codex result M3 で自動緩和を削除): "
      "GPU − CPU double ≤ 2·N_max·ε·max|φ|/h **かつ** 継ぎ目/内部の誤差比 ≤ 2 (下の「G1-b 誤差比」表。"
      "内部が厳密 0 の成分は継ぎ目 ≤ 4ε·max|φ|/h)。h = 合併体積^(1/3) (区分内の最小)、N_max = 区分内の節点の面数の最大。"
      "CPU 再現は plan §4.2 どおり周期半割面を除外")
    P("診断列 (判定外): 「半割面込み」= 周期半割面も面値 φ[ic0] で積算する CPU float32 再現との差 "
      "(GPU の plane_cells では周期半割面の ic1 が ghost なので `excludePeriodic` が効かない仮説の検査)")
    P("| 場 | 区分 | 節点数 | N_max | G1-a 最大差 [ε·max|φ|/h] | G1-b 最大差 [ε·max|φ|/h] | G1-b 閾値 [同] | 判定 a / b | 半割面込み再現との差 [同] |")
    P("| --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    fc = msh.face_count()
    g1a_ok = g1b_ok = True
    ratio_rows = {}
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
            thr_b = 2.0 * Nmax                      # 面数から導く上限のみ (1e-5·S との max・n_member 倍の緩和は削除)
            oka, okb = ea <= 4.0, eb <= thr_b
            g1a_ok &= oka; g1b_ok &= okb
            P(f"| {v} | {lab} | {sel.sum()} | {Nmax} | {ea:.2f} | {eb:.2f} | {thr_b:.1f} | "
              f"{'ok' if oka else 'NG'} / {'ok' if okb else 'NG'} | {dp[sel].max() / unit:.2f} |")
        floor_rows[{"k": "dK (GG)", "omega": "dΩ (GG)", "Xi": "dξ (GG、Y の代理)"}[v]] = (db, G.EPS32 * phimax / V13)
        ratio_rows[v] = (np.abs(g - c64) / (G.EPS32 * phimax / V13)[:, None])   # 節点ごとの h で正規化した成分誤差
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
    # G1-b 誤差比 (機械判定): 継ぎ目 (member≥2) / 内部 (member 1) の最大誤差の比。誤差は G1-b 表と同じ量
    # (GPU − CPU double の全成分の最大絶対差) を節点ごとの ε·max|φ|/h (h = 合併体積^(1/3)) で割ったもの
    # (上の「4 量の床」と同じ正規化)。判定は場ごと: 比 ≤ 2。内部の最大誤差が厳密 0 の場は継ぎ目 ≤ 4 (ε·max|φ|/h)。
    # 成分ごとの比は参考 (判定外) として併記する。
    P("\n### G1-b 誤差比 (継ぎ目/内部。判定は場ごと (全成分の max) ≤ 2、内部が厳密 0 なら継ぎ目 ≤ 4ε·max|φ|/h。成分ごとは参考)")
    P("| 場 | 成分 | 内部 最大 [ε·max|φ|/h] | 継ぎ目 最大 [同] | 継ぎ目/内部 | 判定 |")
    P("| --- | --- | --- | --- | --- | --- |")
    seam_ = msh.nmember > 1
    ratio_ok = True
    for v, r in ratio_rows.items():
        ei, es = r[~seam_].max(), r[seam_].max()
        if ei > 0:
            rr = es / ei; ok = rr <= 2.0; rs = f"{rr:.2f}"
        else:
            ok = es <= 4.0; rs = "内部 0 → 継ぎ目 ≤ 4"
        ratio_ok &= ok
        P(f"| {v} | 全成分 (判定) | {ei:.2f} | {es:.2f} | {rs} | {'ok' if ok else 'NG'} |")
        for ci, cn in enumerate("xyz"):
            ci_, cs_ = r[~seam_, ci].max(), r[seam_, ci].max()
            P(f"| {v} | {cn} (参考) | {ci_:.2f} | {cs_:.2f} | {cs_ / ci_ if ci_ > 0 else float('nan'):.2f} | - |")
    verdict["G1b_ratio"] = ratio_ok
    g1b_ok = g1b_ok and ratio_ok
    P(f"G1-b 誤差比: {'ok' if ratio_ok else 'NG'}")
    verdict["G1a"], verdict["G1b"] = g1a_ok, g1b_ok
    P(f"VERDICT G1-a ({variant}): {'PASS' if g1a_ok else 'FAIL'}")
    P(f"VERDICT G1-b ({variant}): {'PASS' if g1b_ok else 'FAIL'} (上限 2·N_max・誤差比 ≤ 2 の両方)")

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


def evaluate_const(variant, run):
    """GPU 定数場: k・ω・ξ = const で 1 step、res_1 の勾配を 3 区分で機械判定する (plan §5.1 #8c、codex result-2 M3)。

    区分 (節点の group 単位。壁 = group のいずれかの member が非周期境界面 (壁・slip 等) を持つ):
      壁なし継ぎ目 (member≥2 ∧ 非壁): max_c |∂φ/∂x_c| ≤ 4ε|φ|/h                         → 判定
      壁∩継ぎ目   (member≥2 ∧ 壁)  : max_c |∂φ/∂x_c − φ·(ΣS_f)_c/V| ≤ 4ε|φ|/h            → 判定
      非継ぎ目壁   (member 1 ∧ 壁)   : 同じ参照との差 (と参照そのものの大きさ)                  → 報告 (判定外)
      内部         (member 1 ∧ 非壁) : max_c |∂φ/∂x_c|                                          → 報告 (判定外)
    参照 φ·ΣS_f/V: ΣS_f = GPU が積算する面 (内部双対面の incidence を符号付き + 非周期境界面、周期半割面は除外) の
    **格納 float32 面ベクトル** (PLANES/surfVect) の double 和を group 全 member で合算、V = 合併体積 (double)。
    = `gharness.gg_merged(φ≡1, 'f64', include_periodic=False)` × φ。壁半割面込みの GG は定数場でも閉包 ΣS_f ≠ 0 の分だけ
    勾配を持つ (float32 の面ベクトルの閉包誤差、継ぎ目に依らない既知制約) ので、壁の節点はこの参照からのずれで判定する。
    h = 節点の合併体積^(1/3)、単位 ε|φ|/h。
    除外 (現行どおり): 境界条件の上書き (壁の k/ω ピン等) で res_0 の値が焼いた定数と違う節点と、その節点に stencil が
    触れる group。k/ω は壁でピンされるので壁近傍の定数場が壊れる (定数場の試験として意味を持たない) ため。
    ξ は壁で上書きされないので除外 0。除外数は表に併記する。"""
    h5 = os.path.join(run, "mesh.h5")
    msh = G.Mesh(h5, os.path.join(run, "bcondConfig.yaml"))
    s0 = G.read_res(os.path.join(run, "res_0.h5"), list(CONST_FIELDS))
    V13 = np.cbrt(msh.vmerged64)
    out = []
    P = out.append
    P(f"\n## G1 定数場 (GPU、k・ω・ξ = const を焼いて 1 step、plan §6 G1・§5.1 #8c の 3 区分)")
    P(f"run: {run}")
    P(G.provenance(run).rstrip())
    P("判定: 壁なし継ぎ目 max_c |∂φ| ≤ 4ε|φ|/h、壁∩継ぎ目 max_c |∂φ − φ·ΣS_f/V| ≤ 4ε|φ|/h "
      "(ΣS_f = 格納 float32 面ベクトルの double 和、周期半割面を除く・非周期境界面を含む、V = 合併体積)。"
      "非継ぎ目壁は同じ参照で報告、内部は参考。単位 ε|φ|/h (h = 節点の合併体積^(1/3))")
    P("除外 (現行どおり): 境界条件の上書きで res_0 が焼いた定数と違う節点と、それに stencil が触れる group "
      "(k/ω の壁ピン: 定数場が BC で壊れるので試験にならない)。ξ は上書きされないので除外 0")
    # 壁 (非周期境界面を持つ group)
    wall = np.zeros(msh.nCells, bool); wall[msh.pc[msh.bnd_planes, 0]] = True
    wg = np.zeros(msh.nCells, bool); wg[msh.root[wall]] = True; wall = wg[msh.root]
    seam = msh.nmember > 1
    classes = [("壁なし継ぎ目", seam & ~wall, "abs", True), ("壁∩継ぎ目", seam & wall, "ref", True),
               ("非継ぎ目壁", ~seam & wall, "ref", False), ("内部", ~seam & ~wall, "abs", False)]
    P("区分の節点数 (除外前): " + "、".join(f"{lab} {int(m.sum())}" for lab, m, _, _ in classes))
    clos = G.gg_merged(msh, np.ones(msh.nCells, np.float32), "f64", "divide", include_periodic=False)   # ΣS_f/V (double)
    P("| 場 | φ | 区分 | 節点 (除外) | 比較 | 最大 [ε|φ|/h] | 参照 |φ·ΣS_f/V| の最大 [同] | 判定 |")
    P("| --- | --- | --- | --- | --- | --- | --- | --- |")
    names = [f"d{g}d{c}" for _, g in CONST_FIELDS.values() for c in "xyz"]
    g1 = G.read_res(os.path.join(run, "res_1.h5"), names)
    ok_all = True
    for v, (c0, gname) in CONST_FIELDS.items():
        c32 = np.float32(c0)
        okv = s0[v] == c32
        badv = ~okv[msh.inc_j] | ~okv[msh.inc_m]
        bg = np.zeros(msh.nCells, bool); bg[msh.root[msh.inc_m[badv]]] = True
        sel = okv & ~bg[msh.root]
        g = np.stack([g1[f"d{gname}d{c}"] for c in "xyz"], 1).astype(np.float64)
        unit = G.EPS32 * abs(c0) / V13
        ref = c0 * clos
        r_abs = np.abs(g).max(axis=1) / unit
        r_ref = np.abs(g - ref).max(axis=1) / unit
        r_clo = np.abs(ref).max(axis=1) / unit
        for lab, m, kind, judged in classes:
            mm = m & sel
            nex = int((m & ~sel).sum())
            if not mm.any():
                st = ("対象 0 (除外)" if m.any() else "対象 0") if judged else "-"
                if judged and not m.any() and lab == "壁なし継ぎ目":
                    ok_all = False; st = "NG (対象 0)"
                P(f"| {v} | {c0:g} | {lab} | 0 ({nex}) | - | - | - | {st} |")
                continue
            r = r_abs if kind == "abs" else r_ref
            e = r[mm].max()
            ok = e <= 4.0
            if judged:
                ok_all &= ok
            P(f"| {v} | {c0:g} | {lab} | {int(mm.sum())} ({nex}) | {'|∂φ|' if kind == 'abs' else '|∂φ − φΣS/V|'} | {e:.2f} | "
              f"{r_clo[mm].max():.2f} | {('ok' if ok else 'NG') if judged else '報告'} |")
    P(f"VERDICT G1 定数場 ({variant}, 3 区分: 壁なし継ぎ目 ≤ 4ε、壁∩継ぎ目 |∂φ − φΣS_f/V| ≤ 4ε): {'PASS' if ok_all else 'FAIL'}")
    return "\n".join(out) + "\n", ok_all


# ====================================================================================================
# S0 (plan gradient-scalar-lsq-unification §6 S0-a〜e): mesh.scalarGradient: lsq のスカラー勾配作用素
# ====================================================================================================
# 使い方 (AWS):
#   python3 g_suite.py s0 <variant> [--scratch DIR] [--kinds sin,const,same] [--eval-only]
#   python3 g_suite.py s0e_case39 --scratch DIR --src RUN_DIR [--eval-only]
# 変種: tgv, tgv_mirror, tgv_bcswap, tgv_repeat, tgv_shift, jitter32, channel (既存 7 本)、
#       box_slip (周期なし・全面 slip の 32³ ジッタ箱 = S0-c の直接比較用)、
#       axi / axi_m1 (軸対称 × 並進周期 = #6a の r6a_prep 61×41、SST + ξ + 2 成分。axi_m1 は mesh.axisymMethod: 1)。
# 場の種類 (kind。run = <scratch>/s0_<variant>_<kind>):
#   sin   : k・ω・ξ を滑らかな周期 sin 場 (gharness.periodic_sin) で焼く → S0-a (double 参照 ≤ 1e-5·S)、S0-d (jitter32)、S0-e (channel)
#   const : k=2・ω=100・ξ=0.5 の定数場 → S0-b
#   same  : ρ = k = ω = ξ = q (q = m/2048, m∈[1024,2047] の量子化 sin 場: q² が float32 で厳密なので ρk/ρ = q が厳密)
#           → S0-c (NS の dρ とスカラー勾配のビット一致)。軸対称 (TP の状態を作り直せない) は対象外
# 化学種 Y の勾配 (dY{s}d*) は出力変数に無い (variables.hpp output_cellValNames) ので S0-a の Y (5 種)・S0-c の dY は測れない。
S0_VARIANTS = ("tgv", "tgv_mirror", "tgv_bcswap", "tgv_repeat", "tgv_shift", "jitter32", "channel", "box_slip", "axi", "axi_m1")
S0_SIN = {"k": (2.0, 0.5, (0.3, 1.1, -0.4)), "omega": (100.0, 20.0, (1.7, -0.6, 0.9)), "Xi": (0.5, 0.2, (-0.8, 0.4, 2.1))}
S0_GN = {"k": "K", "omega": "Omega", "Xi": "Xi"}
S0_CONST = {"k": 2.0, "omega": 100.0, "Xi": 0.5}


def _s0_source(variant, scr):
    """(src_h5, bcond (path or dict), solver cfg dict, 追加ファイル, 変換 mode, 合併を期待するか)。"""
    extra = []
    if variant.startswith("tgv"):
        src, bc = os.path.join(scr, TGV_REL), os.path.join(HERE, "bcondConfig.yaml")
        if variant == "tgv_mirror":
            bc = bcond_variant(bc, "mirror")
        elif variant == "tgv_bcswap":
            bc = bcond_variant(bc, "swap")
        cfg = G.base_solver_cfg(sst=True, tracer=True)
        tr = {"tgv_mirror": "mirror", "tgv_shift": "shift"}.get(variant)
        return src, bc, cfg, extra, tr, True
    if variant in ("jitter32", "channel", "box_slip"):
        sub = {"jitter32": "box32_j", "channel": "channel", "box_slip": "box32_slip"}[variant]
        md = os.path.join(scr, "g_harness_mesh", sub)
        if os.path.exists(os.path.join(md, "mesh.h5")) and os.path.exists(os.path.join(md, "quality.txt")):
            src, bc = os.path.join(md, "mesh.h5"), os.path.join(md, "bcondConfig.yaml")
        elif variant == "jitter32":
            src, bc, _ = mkmesh.make_box_h5(md, 32, 0.2)
        elif variant == "box_slip":
            src, bc, _ = mkmesh.make_box_h5(md, 32, 0.2, mkmesh.slip_box_bcond())
        else:
            src, bc, _ = mkmesh.make_channel_h5(md)
        q = [l for l in open(os.path.join(md, "quality.txt")).read().splitlines() if "VERDICT" in l]
        print("mesh quality:", q[-1] if q else "?")
        cfg = G.base_solver_cfg(sst=True, tracer=True, visc=1.8e-5 if variant == "channel" else 0.0)
        return src, bc, cfg, extra, None, variant != "box_slip"
    if variant in ("axi", "axi_m1"):
        pd = os.path.join(scr, "r6a_prep")
        cfg = yaml.safe_load(open(os.path.join(pd, "solverConfig.yaml")))
        cfg["turbulence"] = {"model": "sst"}
        if variant == "axi_m1":
            cfg["mesh"]["axisymMethod"] = 1
        extra = [os.path.join(pd, "species_db.yaml")]
        return os.path.join(pd, "axi.h5"), os.path.join(pd, "bcondConfig.yaml"), cfg, extra, None, False
    raise SystemExit("unknown S0 variant " + variant)


def _s0_quantized(f):
    """[0.5, 1) の量子化 q = m/2048 (q² が float32 で厳密)。"""
    m = np.clip(np.rint(f * 2048.0), 1024, 2047)
    return m / 2048.0


def s0_prepare(variant, scr, kind):
    run = os.path.join(scr, f"s0_{variant}_{kind}")
    src, bc, cfg, extra, tr, merge = _s0_source(variant, scr)
    cfg = dict(cfg)
    cfg["mesh"] = dict(cfg.get("mesh", {}), scalarGradient="lsq")
    h5name = "axi.h5" if variant.startswith("axi") else "mesh.h5"
    h5 = G.make_run(run, src, cfg, bc, h5name=h5name)
    for f in extra:
        shutil.copy(f, run)
    if tr == "mirror":
        transform_mesh(h5, "mirror")
    elif tr == "shift":
        transform_mesh(h5, "shift", (100.0, 100.0, 100.0))
    msh = G.Mesh(h5, os.path.join(run, "bcondConfig.yaml"))
    w, clean, axes = G.wrapped_coords(msh)
    x = msh.xyz.astype(np.float64)
    ext = x.max(0) - x.min(0)
    if kind == "sin":
        phi = {v: G.periodic_sin(w, ext, c0, amp, ph) for v, (c0, amp, ph) in S0_SIN.items()}
    elif kind == "const":
        phi = {v: np.full(msh.nCells, c) for v, c in S0_CONST.items()}
    elif kind == "same":
        if variant.startswith("axi"):
            raise SystemExit("S0-c (same) は軸対称では作らない")
        q = _s0_quantized(G.periodic_sin(w, ext, 0.75, 0.2, (0.4, -0.3, 1.2)))
        phi = {"k": q, "omega": q, "Xi": q}
    else:
        raise SystemExit("unknown kind " + kind)
    if variant.startswith("axi"):
        # TP 2 成分の状態 (#6a で焼いた ρ・u・e・Y) はそのまま使い、k・ω・ξ だけ焼く
        with h5py.File(h5, "r") as f:
            ro = np.asarray(f["VALUE/ro"], dtype=np.float64)
        G.bake(h5, {"roK": ro * phi["k"], "roOmega": ro * phi["omega"], "roXi": ro * phi["Xi"]})
    else:
        ro = phi["k"] if kind == "same" else 1.0
        fl = {n: c0 + w @ np.asarray(a) for n, (c0, a) in LSQ_FIELDS.items()}
        G.bake(h5, G.prim_state_fields(msh, ro, fl["Ux"], fl["Uy"], fl["Uz"], 101325.0,
                                       k=phi["k"], om=phi["omega"], xi=phi["Xi"]))
    np.savez(os.path.join(run, "baked.npz"), w=w, clean=clean, merge=merge, **phi)
    return run, merge


def _s0_load(run):
    h5 = os.path.join(run, "axi.h5" if os.path.exists(os.path.join(run, "axi.h5")) else "mesh.h5")
    msh = G.Mesh(h5, os.path.join(run, "bcondConfig.yaml"))
    bk = np.load(os.path.join(run, "baked.npz"))
    s0 = G.read_res(os.path.join(run, "res_0.h5"), ["ro", "k", "omega", "Xi"])
    gn = [f"d{g}d{c}" for g in list(S0_GN.values()) + ["ro"] for c in "xyz"]
    g1 = G.read_res(os.path.join(run, "res_1.h5"), gn)
    grad = lambda g: np.stack([g1[f"d{g}d{c}"] for c in "xyz"], axis=1)
    return msh, bk, s0, grad


def _s0_head(run, title):
    return [f"\n## {title}", f"run: {run}", G.provenance(run).rstrip(),
            "起動エコー: " + "; ".join(l.strip() for l in open(os.path.join(run, "forge_run.log"))
                                   if l.startswith("'scalarGradient' effective") or "回転周期" in l)]


def _nan_lines(run):
    bad = []
    for st in (0, 1):
        with h5py.File(os.path.join(run, f"res_{st}.h5"), "r") as f:
            for k in f["VALUE"]:
                a = f["VALUE/" + k][()]
                if a.dtype.kind == "f" and not np.isfinite(a).all():
                    bad.append(f"res_{st}:{k}")
    return bad


def s0_eval_sin(variant, run):
    """S0-a (double 参照 ≤ 1e-5·S)、S0-d (jitter32: GG 参照との差 > 1e-5·S)。"""
    msh, bk, s0, grad = _s0_load(run)
    merged = bool(bk["merge"])
    out = _s0_head(run, f"S0-a 純作用素 ({variant}, sin 場, 参照 = {'合併' if merged else '非合併 (root = 自分)'} LSQ double)")
    P = out.append
    nan = _nan_lines(run)
    P(f"NaN/Inf: {nan if nan else 'なし'}")
    P("判定: max_i |∇φ_gpu − ∇φ_ref| ≤ 1e-5·S、S = max_i |∇φ_ref| (全節点)。参照の入力は res_0 の場 (境界条件の上書き後)")
    P("| 場 | S | 最大差 | 最大差/S | 区分別 最大差/S (member 1/2/4/8) | 判定 |")
    P("| --- | --- | --- | --- | --- | --- |")
    ok = True
    vd = {}
    for v, gname in S0_GN.items():
        ref = G.lsq_merged_ref(msh, s0[v], merged=merged)
        g = grad(gname).astype(np.float64)
        S = np.linalg.norm(ref, axis=1).max()
        dif = np.abs(g - ref).max(axis=1)
        e = dif.max()
        cls = ", ".join(f"{nm}:{dif[msh.nmember == nm].max() / S:.1e}" for nm in (1, 2, 4, 8) if (msh.nmember == nm).any())
        o = e <= 1e-5 * S and np.isfinite(e)
        ok &= o
        P(f"| {v} | {S:.4g} | {e:.3e} | {e / S:.2e} | {cls} | {'ok' if o else 'NG'} |")
        vd[v] = (ref, g, S)
    P(f"  (打ち切り group/節点数 {msh.lsq_ndegen})")
    P("化学種 Y: dY{s}d* は出力変数に無いので測れない (保留)")
    P(f"VERDICT S0-a ({variant}): {'PASS' if ok and not nan else 'FAIL'}")
    res = {"S0-a": ok and not nan}
    if variant == "jitter32":
        P("\n## S0-d 検出力 (負の対照): GPU LSQ と CPU double GG (合併) の差が 1e-5·S を超えること")
        P("| 場 | S (LSQ 参照) | max|g_lsq − g_GG| | /S | 判定 |")
        P("| --- | --- | --- | --- | --- |")
        okd = True
        for v, (ref, g, S) in vd.items():
            gg = G.gg_merged(msh, s0[v], "f64", "divide")
            e = np.abs(g - gg).max()
            o = e > 1e-5 * S
            okd &= o
            P(f"| {v} | {S:.4g} | {e:.3e} | {e / S:.2e} | {'ok' if o else 'NG'} |")
        P(f"VERDICT S0-d ({variant}): {'PASS' if okd else 'FAIL'}")
        res["S0-d"] = okd
    return "\n".join(out) + "\n", res


def s0_eval_const(variant, run):
    """S0-b: ξ は全節点 (壁込み) で勾配 == 0。k/ω は stencil が定数だけの節点で == 0、ピン値に触れる節点は double 参照と ≤ 4ε。"""
    msh, bk, s0, grad = _s0_load(run)
    merged = bool(bk["merge"])
    out = _s0_head(run, f"S0-b 定数場 ({variant}, k=2・ω=100・ξ=0.5)")
    P = out.append
    nan = _nan_lines(run)
    P(f"NaN/Inf: {nan if nan else 'なし'}")
    P("判定: 定数だけの stencil の節点 (group 単位) は ∇φ == 0 (厳密)。境界条件で値が変わった節点に stencil が触れる group は "
      "|∇φ_gpu − ∇φ_ref| ≤ 4ε·max_stencil|φ|/h (参照 = res_0 の場の LSQ double、h = 節点の合併体積^(1/3))。ξ は全節点で == 0")
    root = msh.root if merged else np.arange(msh.nCells)
    V13 = np.cbrt(msh.vmerged64 if merged else msh.vol.astype(np.float64))
    wall = np.zeros(msh.nCells, bool); wall[msh.pc[msh.bnd_planes, 0]] = True
    P("| 場 | 区分 | 節点 | 非零の節点 | 最大 |∇φ| or 差 [ε·max|φ|/h] | 判定 |")
    P("| --- | --- | --- | --- | --- | --- |")
    ok = True
    for v, gname in S0_GN.items():
        c = np.float32(S0_CONST[v])
        g = grad(gname).astype(np.float64)
        okv = s0[v] == c
        badinc = ~okv[msh.inc_j] | ~okv[msh.inc_m]
        touch = np.zeros(msh.nCells, bool); touch[root[msh.inc_m[badinc]]] = True
        touch = touch[root]
        pure = ~touch
        nz = np.any(g[pure] != 0.0, axis=1)
        o1 = not nz.any()
        ok &= o1
        wl = "(壁込み)" if (wall & pure).any() else ""
        P(f"| {v} | stencil が定数のみ {wl} | {int(pure.sum())} | {int(nz.sum())} | {np.abs(g[pure]).max():.3e} (絶対値) | {'ok' if o1 else 'NG'} |")
        if v == "Xi" and touch.any():
            ok = False
            P(f"| {v} | 値が変わった節点に触れる | {int(touch.sum())} | - | - | NG (ξ は全節点で定数のはず) |")
        elif touch.any():
            ref = G.lsq_merged_ref(msh, s0[v], merged=merged)
            # stencil の max|φ|
            mx = np.abs(s0[v]).astype(np.float64)
            smax = mx.copy()
            np.maximum.at(smax, root[msh.inc_m], mx[msh.inc_j])
            smax = smax[root]
            unit = G.EPS32 * smax / V13
            r = (np.abs(g - ref).max(axis=1) / unit)[touch]
            o2 = r.max() <= 4.0
            ok &= o2
            P(f"| {v} | ピン値に触れる | {int(touch.sum())} | - | {r.max():.2f} | {'ok' if o2 else 'NG'} |")
    P(f"VERDICT S0-b ({variant}): {'PASS' if ok and not nan else 'FAIL'}")
    return "\n".join(out) + "\n", {"S0-b": ok and not nan}


def s0_eval_same(variant, run):
    """S0-c: ρ = k = ω = ξ = q の場で dρ (NS LSQ) とスカラー勾配を比べる。周期なし (box_slip) は全節点ビット一致。
    周期あり: member 1・2 はビット一致、3 以上は |差| ≤ 4ε·Σ_m|p_m| (p_m = 合併係数の member ごとの部分和、double)。"""
    msh, bk, s0, grad = _s0_load(run)
    out = _s0_head(run, f"S0-c NS との一致 ({variant}, ρ = k = ω = ξ = q、q = m/2048)")
    P = out.append
    nan = _nan_lines(run)
    P(f"NaN/Inf: {nan if nan else 'なし'}")
    gro = grad("ro")
    per = bool(bk["merge"])
    part = np.abs(G.lsq_partials(msh, s0["ro"])) if per else None
    sabs = np.zeros((msh.nCells, 3))
    if per:
        np.add.at(sabs, msh.root, part); sabs = sabs[msh.root]
    P("比較対象: res_0 でスカラーの値が ρ とビット一致する節点だけで stencil (group 単位) が閉じる節点 (境界条件で k・ω が変わる節点に触れる group は除外)")
    P("| 場 | 区分 | 節点 (除外) | ビット不一致 | 最大 |差| | 最大 |差|/(ε·Σ|p_m|) | 判定 |")
    P("| --- | --- | --- | --- | --- | --- | --- |")
    ok = True
    for v, gname in S0_GN.items():
        same = s0[v].view(np.uint32) == s0["ro"].view(np.uint32)
        badinc = ~same[msh.inc_j] | ~same[msh.inc_m]
        touch = np.zeros(msh.nCells, bool); touch[msh.root[msh.inc_m[badinc]]] = True
        sel0 = same & ~touch[msh.root]
        g = grad(gname)
        mism = np.any(g.view(np.uint32) != gro.view(np.uint32), axis=1)
        dif = np.abs(g.astype(np.float64) - gro.astype(np.float64))
        for nm in sorted(set(msh.nmember.tolist())):
            m = msh.nmember == nm
            sel = sel0 & m
            nex = int((m & ~sel0).sum())
            if not sel.any():
                P(f"| {v} | member {nm} | 0 ({nex}) | - | - | - | - |")
                continue
            if nm <= 2:
                o = not mism[sel].any()
                rr = "-"
            else:
                lim = 4.0 * G.EPS32 * sabs[sel]
                o = bool(np.all(dif[sel] <= lim))
                rr = f"{(dif[sel] / np.maximum(G.EPS32 * sabs[sel], 1e-300)).max():.2f}"
            ok &= o
            P(f"| {v} | member {nm}{' (ビット一致を要求)' if nm <= 2 else ''} | {int(sel.sum())} ({nex}) | {int(mism[sel].sum())} | "
              f"{dif[sel].max():.3e} | {rr} | {'ok' if o else 'NG'} |")
    P("化学種 dY と dξ のビット一致 (5 種、チャンク 4+1): dY{s}d* が出力変数に無いので測れない (保留)")
    P(f"VERDICT S0-c ({variant}): {'PASS' if ok and not nan else 'FAIL'}")
    return "\n".join(out) + "\n", {"S0-c": ok and not nan}


def s0_eval_e(run, label, h5name=None):
    """S0-e (診断): res_0 (初期化の applyBconds 直後の状態) と res_1 の k・ω の周期 group 内差 (member − root)。
    wall_y_eff は出力変数に無いので測れない (保留)。"""
    h5 = os.path.join(run, h5name) if h5name else os.path.join(run, "mesh.h5")
    msh = G.Mesh(h5, os.path.join(run, "bcondConfig.yaml"))
    out = [f"\n## S0-e BC 後の周期 group 内の同値性 (診断、{label})", f"run: {run}", G.provenance(run).rstrip()]
    P = out.append
    wall = np.zeros(msh.nCells, bool); wall[msh.pc[msh.bnd_planes, 0]] = True
    wg = np.zeros(msh.nCells, bool); wg[msh.root[wall]] = True; wallg = wg[msh.root]
    mem = msh.root != np.arange(msh.nCells)
    P(f"周期 group {len(msh.groups)} (壁∩継ぎ目の group {len({int(msh.root[c]) for c in np.nonzero(mem & wallg)[0]})})")
    P("| res | 場 | 区分 | member 節点 | group 内差が非零の member | 最大 |φ_m − φ_root| | 最大 相対 |")
    P("| --- | --- | --- | --- | --- | --- | --- |")
    nz_any = False
    for st in (0, 1):
        s = G.read_res(os.path.join(run, f"res_{st}.h5"), ["k", "omega", "roK", "roOmega", "wall_y_eff"])
        for v in ("k", "omega", "roK", "roOmega", "wall_y_eff"):
            if v not in s:
                continue
            a = s[v].astype(np.float64)
            d = np.abs(a - a[msh.root])
            for lab, m in (("壁∩継ぎ目", mem & wallg), ("壁なし継ぎ目", mem & ~wallg)):
                if not m.any():
                    continue
                nzm = int((d[m] != 0).sum())
                nz_any |= nzm > 0
                rel = (d[m] / np.maximum(np.abs(a[msh.root][m]), 1e-300)).max()
                P(f"| {st} | {v} | {lab} | {int(m.sum())} | {nzm} | {d[m].max():.3e} | {rel:.2e} |")
    if "wall_y_eff" not in s:
        P("wall_y_eff: res に無い (extraFields 未指定か、#4a 以前のバイナリ)")
    P(f"S0-e 記録 ({label}): k・ω{'・wall_y_eff' if 'wall_y_eff' in s else ''} の group 内差 {'非零あり' if nz_any else 'すべて 0'} (判定ではなく記録)")
    return "\n".join(out) + "\n", {"S0-e_nonzero": nz_any}


def s0e_case39_prepare(scr, src_run):
    """case/39 run_0039_r1_gradfix_new_ext の res_800000 を restart_field.py で入力へ写し、lsq で 1 step。"""
    run = os.path.join(scr, "s0e_case39")
    if os.path.exists(run):
        shutil.rmtree(run)
    os.makedirs(run)
    cfg = yaml.safe_load(open(os.path.join(src_run, "solverConfig.yaml")))
    for f in (cfg["mesh"]["meshFileName"], cfg["mesh"]["valueFileName"], "bcondConfig.yaml", "probe.yaml"):
        if os.path.exists(os.path.join(src_run, f)):
            shutil.copy(os.path.join(src_run, f), run)
    cfg["mesh"]["scalarGradient"] = "lsq"
    out_cfg = cfg.setdefault("output", {})
    out_cfg["extraFields"] = list(out_cfg.get("extraFields") or []) + ["wall_y_eff"]   # #4a の extraFields 登録 (post-gather 値)
    cfg["time"]["last"]["nStepOuter"] = 1
    cfg["time"]["outStepStart"] = 0
    cfg["time"]["outStepInterval"] = 1
    with open(os.path.join(run, "solverConfig.yaml"), "w") as fp:
        yaml.safe_dump(cfg, fp, sort_keys=False)
    rf = os.path.join(G.REPO, "solver_density_cuda", "tools", "restart_field.py")
    r = subprocess.run(["python3", rf, os.path.join(src_run, "res_800000.h5"), os.path.join(run, cfg["mesh"]["valueFileName"])],
                       capture_output=True, text=True)
    open(os.path.join(run, "restart_field.log"), "w").write(r.stdout + r.stderr)
    if r.returncode != 0:
        raise SystemExit("restart_field failed: " + r.stdout + r.stderr)
    return run, cfg["mesh"]["meshFileName"]


def s0_main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("variant")
    ap.add_argument("--scratch", default=SCR_DEFAULT)
    ap.add_argument("--kinds", default="sin,const,same")
    ap.add_argument("--eval-only", action="store_true")
    ap.add_argument("--src", default=None, help="s0e_case39: 起点 run ディレクトリ")
    a = ap.parse_args(argv)
    if a.variant == "s0e_case39":
        run = os.path.join(a.scratch, "s0e_case39")
        cfg = yaml.safe_load(open(os.path.join(a.src, "solverConfig.yaml")))
        mname = cfg["mesh"]["meshFileName"]
        if not a.eval_only:
            run, mname = s0e_case39_prepare(a.scratch, a.src)
            G.run_forge(run)
        txt, _ = s0_eval_e(run, "case/39 run_0039_r1_gradfix_new_ext res_800000 起点、lsq 1 step", h5name=mname)
        txt = f"# S0-e case/39 (plan gradient-scalar-lsq-unification §6 S0-e)\nharness revision: {G.git_rev()}\n" + txt
        outp = os.path.join(HERE, "S0e_case39.txt")
        open(outp, "w").write(txt); print(txt); print("->", outp)
        return
    if a.variant not in S0_VARIANTS:
        raise SystemExit("unknown S0 variant " + a.variant)
    kinds = [k for k in a.kinds.split(",") if k]
    if a.variant.startswith("axi"):
        kinds = [k for k in kinds if k != "same"]
    txt = f"# S0 {a.variant} (plan gradient-scalar-lsq-unification §6 S0、mesh.scalarGradient: lsq)\nharness revision: {G.git_rev()}\n"
    verdict = {}
    for kind in kinds:
        run = os.path.join(a.scratch, f"s0_{a.variant}_{kind}")
        if not a.eval_only:
            run, merge = s0_prepare(a.variant, a.scratch, kind)
            G.run_forge(run, expect_merge=merge)
        if kind == "sin":
            t, v = s0_eval_sin(a.variant, run)
            if a.variant == "channel":
                te, ve = s0_eval_e(run, "channel sin 場、lsq 1 step")
                t += te; v.update(ve)
        elif kind == "const":
            t, v = s0_eval_const(a.variant, run)
        else:
            t, v = s0_eval_same(a.variant, run)
        txt += t
        verdict.update(v)
    txt += "\n## まとめ\n" + "\n".join(f"{k}: {('PASS' if x else 'FAIL') if not k.endswith('nonzero') else ('非零あり' if x else '0')}"
                                       for k, x in verdict.items()) + "\n"
    outp = os.path.join(HERE, f"S0_{a.variant}.txt")
    open(outp, "w").write(txt)
    print(txt)
    print("->", outp)


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "s0":
        return s0_main(sys.argv[2:])
    ap = argparse.ArgumentParser()
    ap.add_argument("variant")
    ap.add_argument("--scratch", default=SCR_DEFAULT)
    ap.add_argument("--eval-only", action="store_true")
    ap.add_argument("--no-const", action="store_true", help="GPU 定数場の run を省く (旧形式の再評価用)")
    a = ap.parse_args()
    run = os.path.join(a.scratch, "g_harness_" + a.variant)
    runc = os.path.join(a.scratch, "g_harness_" + a.variant + "_const")
    if not a.eval_only:
        run = prepare(a.variant, a.scratch)
        G.run_forge(run)
        if not a.no_const:
            runc = prepare(a.variant, a.scratch, const=True)
            G.run_forge(runc)
    ref = os.path.join(a.scratch, "g_harness_tgv") if a.variant in ("tgv_bcswap", "tgv_repeat") else None
    txt, v = evaluate(a.variant, run, ref)
    if not a.no_const:
        tc, vc = evaluate_const(a.variant, runc)
        txt += tc
    outp = os.path.join(HERE, f"G_{a.variant}.txt")
    open(outp, "w").write(txt)
    print(txt)
    print("->", outp)


if __name__ == "__main__":
    main()
