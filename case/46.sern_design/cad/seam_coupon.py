#!/usr/bin/env python3
"""継ぎ目分布の局所 A′/B′ (plan tooling-sern-mesh-blocking §6.3、§5.1 B4-2 (1))。

生産板厚 (TC 0.005 H) の x/H = 1.0 断面から `U1b|U2b` (カウル外壁の壁帯、z = W/2 の継ぎ目) だけを切り出した接続試験片を作り、
本体 (hex_junction_model.py) と同じ検査関数で判定する。変えるのは `U1b` の z 分布則 (hz10・hz20) だけ:
  A′: z 一様 49 節点 (継ぎ目比の予測 ZW/48/h1 = 130.21)
  B′: 継ぎ目側 (z = W/2) の第一間隔 h1 の片側等比 49 節点 (48 区間で 1 H、公比 1.153867、最大間隔 0.133488 H)
形状・接続・壁法線分布 (vy1k: カウル外壁へ NL 層)・`NSW` (U2b の z、両端 h1 の Bump) は本体の laws と同じ式で固定。
x は x/H = 1.0 から HX (0.05 H) 間隔の短い押出し (3 区間・4 station、全節点が station 面上)。CFD はしない。

usage (mesh venv):  .venv-mesh/bin/python case/46.sern_design/cad/seam_coupon.py OUT_DIR --contours DIR
出力: OUT_DIR/coupon_{A,B}.msh・coupon_{A,B}_report.json・coupon_report.json (分岐判定つき)
"""
import sys, json, math, argparse
from pathlib import Path
import numpy as np
import gmsh

sys.path.insert(0, str(Path(__file__).resolve().parent))
from hex_junction_model import (P0, Geom, Bump, prog_r, prog_n, load_contours, first_layer_check, adjacent_spacing_check,  # noqa: E402
                                hex_quality)

# 断面の点 (y 段 j, z 段 k) と辺。U1b = (j 1..2, k 0..1)、U2b = (j 1..2, k 1..2)
BLOCKS = {"U1b": ("hz10", "vy11", "-hz20", "-vy10"), "U2b": ("hz11", "vy12", "-hz21", "-vy11")}
EDGES = {"hz10": ((1, 0), (1, 1)), "hz20": ((2, 0), (2, 1)), "hz11": ((1, 1), (1, 2)), "hz21": ((2, 1), (2, 2)),
         "vy10": ((1, 0), (2, 0)), "vy11": ((1, 1), (2, 1)), "vy12": ((1, 2), (2, 2))}


def build_coupon(out, P, contours, variant, x0=1.0, nx=3, nzu=49):
    G = Geom(P, contours["ramp"], contours["cowl"]); H = G.H; h1 = P["H1"] * H; g = P["G"]; DR = G.DR
    if not (G.LSW < x0 * H and (x0 + nx * P["HX"]) * H < G.LCOWL): raise ValueError("試験片は区間 B (L_sw < x < L_cowl) に置く")
    gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0)
    bump = Bump(); NSW = bump.fit_n(h1 / G.TSW, g); NL = prog_n(DR, h1, g) + 4
    gmsh.model.add("coupon"); ge = gmsh.model.geo; gm = ge.mesh
    xs = [(x0 + i * P["HX"]) * H for i in range(nx + 1)]
    law = {}
    # U1b の z 分布 (hz10・hz20 は同値類なので同じ則)。gmsh の Progression は始点 (z = 0) から係数倍で並ぶので、継ぎ目側を細かくするには 1/r
    if variant == "A":
        law["z_u1"] = (nzu, "Progression", 1.0); pred = G.ZW / (nzu - 1) / h1
    else:
        r = prog_r(G.ZW, nzu - 1, h1); law["z_u1"] = (nzu, "Progression", 1.0 / r); pred = 1.0
    cf = bump.coef(NSW, h1 / G.TSW); law["z_u2"] = (NSW, "Progression", 1.0) if cf is None else (NSW, "Bump", cf)
    S = []
    for x in xs:
        co = math.sqrt(1 + G.dyo(x) ** 2); yo = G.yo(x); ys = {1: yo - DR, 2: yo}; zs = {0: 0.0, 1: G.ZW, 2: G.ZO}
        pid = {(j, k): ge.addPoint(x, ys[j], zs[k]) for j in (1, 2) for k in (0, 1, 2)}
        eid = {e: ge.addLine(pid[a], pid[b]) for e, (a, b) in EDGES.items()}
        hh = min(h1 * co, 0.999 * DR / NL); rv = prog_r(DR, NL, hh)          # 本体 laws の vy1k と同じ (壁 = 終点側 j = 2 へ細かく)
        L = {"hz10": law["z_u1"], "hz20": law["z_u1"], "hz11": law["z_u2"], "hz21": law["z_u2"]}
        for e in ("vy10", "vy11", "vy12"): L[e] = (NL + 1, "Progression", -rv)
        loop = lambda b: ge.addCurveLoop([(-1 if s.startswith("-") else 1) * eid[s.lstrip("-")] for s in BLOCKS[b]])
        face = {b: ge.addSurfaceFilling([loop(b)]) for b in BLOCKS}
        S.append(dict(pid=pid, eid=eid, face=face, law=L))
    vols, xf, tagq = {}, {}, []
    for i in range(nx):
        A, B = S[i], S[i + 1]; xl = {}
        for p in A["pid"]:
            xl[p] = ge.addLine(A["pid"][p], B["pid"][p]); gm.setTransfiniteCurve(xl[p], 2, "Progression", 1.0)
        for e, (a, b) in EDGES.items():
            xf[(i, e)] = ge.addSurfaceFilling([ge.addCurveLoop([A["eid"][e], xl[b], -B["eid"][e], -xl[a]])])
        for b, loop_ in BLOCKS.items():
            sl = ge.addSurfaceLoop([A["face"][b], B["face"][b]] + [xf[(i, s.lstrip("-"))] for s in loop_]); vols[(i, b)] = ge.addVolume([sl])
        tagq += [xf[(i, "hz20")], xf[(i, "hz21")]]                     # カウル外壁 (cowl_out)
    for St in S:
        for e, tg in St["eid"].items(): n, typ, c_ = St["law"][e]; gm.setTransfiniteCurve(tg, n, typ, c_)
        for f in St["face"].values(): gm.setTransfiniteSurface(f); gm.setRecombine(2, f)
    for f in xf.values(): gm.setTransfiniteSurface(f); gm.setRecombine(2, f)
    for v in vols.values(): gm.setTransfiniteVolume(v)
    ge.synchronize()
    gmsh.model.addPhysicalGroup(2, tagq, name="cowl_out"); gmsh.model.addPhysicalGroup(3, list(vols.values()), name="fluid")
    gmsh.model.mesh.generate(3)
    M = gmsh.model.mesh; ntag, xyz, _ = M.getNodes(); xyz = xyz.reshape(-1, 3)
    idx = np.zeros(int(ntag.max()) + 1, np.int64); idx[ntag.astype(np.int64)] = np.arange(len(ntag))
    hx_l, lab = [], []
    names = sorted(BLOCKS)
    for (i, b), v in vols.items():
        et, _, cn = M.getElements(3, v); assert list(et) == [5], et
        hv = idx[np.asarray(cn[0], np.int64)].reshape(-1, 8); hx_l.append(hv); lab.append(np.full(len(hv), names.index(b), np.int16))
    hx = np.vstack(hx_l); lab = np.concatenate(lab)
    q = [idx[np.asarray(M.getElements(2, s_)[2][0], np.int64)].reshape(-1, 4) for s_ in tagq]
    bq = {"cowl_out": np.vstack(q)}
    # ---- 本体と同じ検査関数
    J, J32, sk, AR, _ = hex_quality(xyz, hx); sg = np.sign(np.median(J))
    fl = first_layer_check(xyz, hx, bq, h1, h1, NL, NL, DR, wall_tags=("cowl_out",), end_tags=(), expected_tags=("cowl_out",), gmax=g, unit=H)
    adj = adjacent_spacing_check(xyz, hx, gmax=g, labels=lab, label_names=names, unit=H)
    # z 分布の実測 (U1b の壁面 = 継ぎ目に向かう間隔列、x = x0 の station)
    w = xyz[np.unique(bq["cowl_out"])]; w = w[np.abs(w[:, 0] - xs[0]) < 1e-12]; zz = np.sort(w[:, 2]); dz = np.diff(zz)
    u1 = zz[1:] <= G.ZW + 1e-12; dz1 = dz[u1]
    rep = dict(variant=variant, x0=x0, stations=[x / H for x in xs], nodes=int(len(xyz)), hexes=int(len(hx)), NL=NL, NSW=NSW, nz_u1b=nzu,
               z_law_u1b=law["z_u1"], seam_ratio_predicted=pred,
               seam_ratio_measured=adj.get("max_by_block_pair", {}).get("U1b|U2b|section"),
               in_block_section_ratio={k: v for k, v in adj.get("max_by_block_pair", {}).items() if k.split("|")[0] == k.split("|")[1]},
               u1b_dz_seam=float(dz1[-1]) / H, u1b_dz_max=float(dz1.max()) / H, u1b_ratio_max=float(np.max(dz1[:-1] / dz1[1:])) if len(dz1) > 1 else None,
               neg_jac=int((J * sg <= 0).sum()), neg_jac_f32=int((J32 * sg <= 0).sum()), skew_max=float(sk.max()), ar_max=float(AR.max()),
               ar_gt1000=int((AR > 1000).sum()), ar_gt5000=int((AR > 5000).sum()),
               ar1000_skew_max=float(sk[AR > 1000].max()) if (AR > 1000).any() else 0.0,
               first_layer={k: v for k, v in fl["cowl_out"].items() if not k.startswith("_")}, adjacent_spacing={k: v for k, v in adj.items() if k != "worst"})
    gates = dict(adjacent_spacing=adj["ok"], first_layer=fl["cowl_out"]["ok"], neg_jac_f32=rep["neg_jac_f32"] == 0, skew=rep["skew_max"] <= 0.90,
                 ar=rep["ar_gt5000"] == 0 and rep["ar1000_skew_max"] <= 0.30)
    rep["gates"] = gates; rep["VERDICT"] = "PASS" if all(gates.values()) else "FAIL"
    json.dump(rep, open(out + "_report.json", "w"), indent=1, ensure_ascii=False)
    gmsh.write(out + ".msh"); gmsh.finalize(); return rep


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("outdir"); ap.add_argument("--contours", required=True)
    a = ap.parse_args(); od = Path(a.outdir); od.mkdir(parents=True, exist_ok=True)
    P = dict(P0); P["TC"] = 0.005; C = load_contours(a.contours); res = {}
    for v in ("A", "B"):
        res[v] = build_coupon(str(od / f"coupon_{v}"), P, C, v)
        r = res[v]
        print(f"{v}′: 継ぎ目比 実測 {r['seam_ratio_measured']:.4f} (予測 {r['seam_ratio_predicted']:.4f})  U1b dz 継ぎ目 {r['u1b_dz_seam']:.4e} H 最大 {r['u1b_dz_max']:.6f} H "
              f"公比 {r['u1b_ratio_max']:.6f}  ブロック内 {r['in_block_section_ratio']}  第一層 {r['first_layer']['ratio_min']:.4f}–{r['first_layer']['ratio_max']:.4f} "
              f"層 {r['first_layer']['layers_min']}/{r['first_layer']['layers_expected']}  負J32 {r['neg_jac_f32']} skew {r['skew_max']:.4f} AR {r['ar_max']:.1f}  "
              f"隣接比 max {r['adjacent_spacing']['ratio_max']:.4f}  {r['VERDICT']} {r['gates']}")
    A, B = res["A"], res["B"]
    a_rep = abs(A["seam_ratio_measured"] / A["seam_ratio_predicted"] - 1) < 1e-3
    b_ok = B["VERDICT"] == "PASS" and B["seam_ratio_measured"] <= 1.2
    branch = ("A′ で予測の跳びを再現し B′ で比 <= 1.2 かつ品質・壁層合格 -> 全継ぎ目へ (B4-2)" if a_rep and b_ok else
              ("B′ でも不合格 -> 片側分布だけで修復できるという仮説を棄却 (止める)" if not b_ok else "A′ が予測の跳びを再現しない (止める)"))
    json.dump(dict(A=A, B=B, A_reproduces=a_rep, B_pass=b_ok, branch=branch), open(od / "coupon_report.json", "w"), indent=1, ensure_ascii=False)
    print("分岐:", branch); sys.exit(0 if a_rep and b_ok else 1)
