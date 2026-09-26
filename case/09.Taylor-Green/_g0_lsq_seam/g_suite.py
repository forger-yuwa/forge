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


def main():
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
