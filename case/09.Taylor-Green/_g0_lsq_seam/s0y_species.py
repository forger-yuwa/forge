#!/usr/bin/env python3
"""#4c (plan gradient-scalar-lsq-unification §5.1 #4c): 化学種 Y の未測定項目を 5 種 TP で測る (mesh.scalarGradient: lsq、1 step)。

5 種 = case/16 `run_0471` の [H2O, N2, O2, AR, CO2] (thermalMethod 2、thermoHrefTemp 298.15、species_db は run_0471 のもの)。
場: ρ = 1 (Y = ρY/ρ・ξ = ρξ/ρ が厳密)、Y0 = Y4 = ξ = q (q = m/2048 の量子化 sin 場、q∈[0.13,0.37] < 0.5)、
    Y1..Y3 は正の sin 重みで Σ = 1 − 2q を double で作り、Y1・Y2 を 2^-20 格子に量子化して Y3 = 1 − 2q − Y1 − Y2 (Σ が厳密に 1)。
    k・ω は S0 と同じ sin 場、速度は S0 と同じ線形場、T = 300 K の内部エネルギー (NASA9、sensible datum 298.15 K) で roe を焼く。
測る項目:
  S0-a (Y): dY{s} (res_1、extraFields で出力、周期 gather 後) が res_0 の Y{s} を入れた double LSQ 参照 (合併) と ≤ 1e-5·S
  S0-c (dY と dξ): gather 前 (FORGE_DUMP_PREGATHER の species_lsq.loop1 / passive_lsq.loop1) で dY0 = dξ、dY4 = dξ がビット一致
       (Y0 は先頭チャンク NV=4、Y4 はチャンク境界の孤立変数 NV=1、ξ は受動種 NV=1)。gather 後 (res_1) は member ≤ 2 でビット一致、
       3 以上は各自の gather 前配列の float32 順列和の集合に含まれる。周期なし (box_slip) は res_1 同士をビット比較
  Σ dY = 0: res_1 の全節点・全成分で |Σ_s dY_s| ≤ 4ε Σ_s |dY_s|
  S0-e (wall_y_eff): channel の res_0 / res_1 の wall_y_eff (extraFields) の周期 group 内差 (記録)

使い方 (AWS):
  python3 s0y_species.py <variant> --scratch DIR --species-db PATH [--eval-only]
  変種: tgv, jitter32, channel, box_slip
"""
import argparse
import os
import shutil
import sys

import h5py
import numpy as np
import yaml

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gharness as G  # noqa: E402
import g_suite as GS  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
SPECIES = ["H2O", "N2", "O2", "AR", "CO2"]
NS = len(SPECIES)
RU = 8.314462618
T0 = 300.0
TREF = 298.15


def nasa9_h(c, T):
    """NASA9 の h [J/mol] (係数 a1..a7, b1, b2)。"""
    a = c
    hRT = -a[0] * T ** -2 + a[1] * np.log(T) / T + a[2] + a[3] * T / 2 + a[4] * T ** 2 / 3 + a[5] * T ** 3 / 4 \
        + a[6] * T ** 4 / 5 + a[7] / T
    return hRT * RU * T


def e_sensible(db, Y, T):
    """sensible datum (h_s(TREF) = 0) の内部エネルギー [J/kg]: Σ Y_s (h_s(T) − h_s(TREF))/MW_s − R_mix T。"""
    e = np.zeros_like(Y[0])
    Rm = np.zeros_like(Y[0])
    for s, nm in enumerate(SPECIES):
        d = db[nm]
        c = d["nasa9_low"] if T < d["Tmid"] else d["nasa9_high"]
        cr = d["nasa9_low"] if TREF < d["Tmid"] else d["nasa9_high"]
        dh = (nasa9_h(c, T) - nasa9_h(cr, TREF)) / d["MW"]
        e += Y[s] * dh
        Rm += Y[s] * RU / d["MW"]
    return e - Rm * T


def prepare(variant, scr, species_db):
    run = os.path.join(scr, f"s0y_{variant}")
    src, bc, cfg, extra, tr, merge = GS._s0_source(variant, scr)
    cfg = dict(cfg)
    cfg["mesh"] = dict(cfg.get("mesh", {}), scalarGradient="lsq")
    cfg["physProp"] = dict(cfg["physProp"], thermalMethod=2, cp=1039.0, gamma=1.4, species=SPECIES,
                           speciesDBFile="species_db.yaml", thermoHrefTemp=TREF)
    ex = ["dt_local", "volume", "wall_y_eff"] + [f"dY{s}d{c}" for s in range(NS) for c in "xyz"]
    cfg["output"] = {"level": 2, "extraFields": ex}
    h5 = G.make_run(run, src, cfg, bc)
    shutil.copy(species_db, os.path.join(run, "species_db.yaml"))
    msh = G.Mesh(h5, os.path.join(run, "bcondConfig.yaml"))
    w, clean, axes = G.wrapped_coords(msh)
    x = msh.xyz.astype(np.float64)
    ext = x.max(0) - x.min(0)
    q = np.rint((0.25 + 0.12 * (G.periodic_sin(w, ext, 0.0, 1.0, (0.4, -0.3, 1.2)) / 1.7)) * 2048.0) / 2048.0
    R = 1.0 - 2.0 * q
    wts = [G.periodic_sin(w, ext, 1.0, 0.3, ph) for ph in ((0.1, 0.9, -0.5), (1.3, -0.2, 0.7), (-0.9, 0.5, 2.0))]
    sw = wts[0] + wts[1] + wts[2]
    g20 = 2.0 ** 20
    Y1 = np.rint(R * wts[0] / sw * g20) / g20
    Y2 = np.rint(R * wts[1] / sw * g20) / g20
    Y3 = R - Y1 - Y2
    Y = [q, Y1, Y2, Y3, q]
    assert (np.array(Y) > 0).all() and np.all(np.abs(sum(Y) - 1.0) == 0.0)
    k = G.periodic_sin(w, ext, *GS.S0_SIN["k"])
    om = G.periodic_sin(w, ext, *GS.S0_SIN["omega"])
    db = yaml.safe_load(open(species_db))
    fl = {n: c0 + w @ np.asarray(a) for n, (c0, a) in GS.LSQ_FIELDS.items()}
    ro = np.ones(msh.nCells)
    e = e_sensible(db, Y, T0)
    out = {"ro": ro, "roUx": ro * fl["Ux"], "roUy": ro * fl["Uy"], "roUz": ro * fl["Uz"],
           "roe": ro * (e + 0.5 * (fl["Ux"] ** 2 + fl["Uy"] ** 2 + fl["Uz"] ** 2)),
           "roK": ro * k, "roOmega": ro * om, "roXi": ro * q}
    for s in range(NS):
        out[f"roY{s}"] = ro * Y[s]
    G.bake(h5, out)
    np.savez(os.path.join(run, "baked.npz"), merge=merge, q=q)
    return run, merge


def evaluate(variant, run):
    msh = G.Mesh(os.path.join(run, "mesh.h5"), os.path.join(run, "bcondConfig.yaml"))
    bk = np.load(os.path.join(run, "baked.npz"))
    merged = bool(bk["merge"])
    out = [f"# #4c 化学種 5 種 ({', '.join(SPECIES)}) の Y 項目 ({variant}、lsq 1 step)",
           f"harness revision: {G.git_rev()}", f"run: {run}", G.provenance(run).rstrip(),
           "起動エコー: " + "; ".join(l.strip() for l in open(os.path.join(run, "forge_run.log"))
                                   if l.startswith("'scalarGradient' effective") or "[FORGE_DUMP_PREGATHER]" in l)]
    P = out.append
    nan = GS._nan_lines(run)
    P(f"NaN/Inf: {nan if nan else 'なし'}")
    s0 = G.read_res(os.path.join(run, "res_0.h5"), [f"Y{s}" for s in range(NS)] + ["Xi", "T"])
    P(f"res_0 の T: min {s0['T'].min():.2f} / max {s0['T'].max():.2f} K、ΣY − 1: max |.| "
      f"{np.abs(sum(s0[f'Y{s}'].astype(np.float64) for s in range(NS)) - 1.0).max():.2e}")
    g1 = G.read_res(os.path.join(run, "res_1.h5"), [f"dY{s}d{c}" for s in range(NS) for c in "xyz"] + [f"dXid{c}" for c in "xyz"])
    gY = [np.stack([g1[f"dY{s}d{c}"] for c in "xyz"], 1) for s in range(NS)]
    gX = np.stack([g1[f"dXid{c}"] for c in "xyz"], 1)
    res = {}
    # S0-a (Y)
    P(f"\n## S0-a (Y): res_1 の dY{{s}} (gather 後) vs {'合併' if merged else '非合併'} LSQ double (入力 res_0 の Y{{s}})、≤ 1e-5·S")
    P("| 種 | S | 最大差 | 最大差/S | 判定 |")
    P("| --- | --- | --- | --- | --- |")
    ok = True
    for s in range(NS):
        ref = G.lsq_merged_ref(msh, s0[f"Y{s}"], merged=merged)
        S = np.linalg.norm(ref, axis=1).max()
        e = np.abs(gY[s].astype(np.float64) - ref).max()
        o = e <= 1e-5 * S
        ok &= o
        P(f"| Y{s} ({SPECIES[s]}) | {S:.4g} | {e:.3e} | {e / S:.2e} | {'ok' if o else 'NG'} |")
    P(f"VERDICT S0-a Y ({variant}): {'PASS' if ok and not nan else 'FAIL'}")
    res["S0-a(Y)"] = ok and not nan
    # S0-c (dY と dξ)
    P("\n## S0-c (dY と dξ): 前提 res_0 で Y0 = Y4 = ξ (ビット)")
    same0 = (s0["Y0"].view(np.uint32) == s0["Xi"].view(np.uint32)); same4 = (s0["Y4"].view(np.uint32) == s0["Xi"].view(np.uint32))
    P(f"res_0 で Y0≠ξ の節点 {int((~same0).sum())}、Y4≠ξ の節点 {int((~same4).sum())} (全 {msh.nCells})")
    okc = bool(same0.all() and same4.all())
    dump = os.path.join(run, "pregather")
    have = os.path.exists(dump + ".species_lsq.loop1") and os.path.exists(dump + ".passive_lsq.loop1")
    P("| 段 | 比較 | 対象節点 | ビット不一致 | 判定 |")
    P("| --- | --- | --- | --- | --- |")
    if have:
        pY = G.read_pregather(dump + ".species_lsq.loop1")
        pX = G.read_pregather(dump + ".passive_lsq.loop1")["Xi"]
        for s in (0, 4):
            m = np.any(pY[f"Y{s}"].view(np.uint32) != pX.view(np.uint32), axis=1)
            o = not m.any()
            okc &= o
            P(f"| gather 前 (ダンプ) | dY{s} vs dξ | {msh.nCells} | {int(m.sum())} | {'ok' if o else 'NG'} |")
    else:
        okc = False
        P("| gather 前 | ダンプが無い | - | - | NG |")
    for s in (0, 4):
        m = np.any(gY[s].view(np.uint32) != gX.view(np.uint32), axis=1)
        le2 = msh.nmember <= 2
        o = not m[le2].any()
        okc &= o
        P(f"| gather 後 (res_1) | dY{s} vs dξ、member ≤ 2 | {int(le2.sum())} | {int(m[le2].sum())} | {'ok' if o else 'NG'} |")
        if (~le2).any() and have:
            bad = 0
            for name, pre, post in ((f"Y{s}", pY[f"Y{s}"], gY[s]), ("Xi", pX, gX)):
                sets = G.gather_perm_sums(msh, pre)
                u = post.view(np.uint32)
                for rt, s3 in sets.items():
                    if len(msh.groups[rt]) <= 2:
                        continue
                    for node in msh.groups[rt]:
                        for c in range(3):
                            bad += u[node, c].item() not in s3[c]
            o = bad == 0
            okc &= o
            P(f"| gather 後 (res_1) | dY{s}・dξ、member ≥ 3 は順列和の集合に含まれる | {int((~le2).sum())} | 集合外 {bad} "
              f"(参考: dY{s}≠dξ の節点 {int(m[~le2].sum())}) | {'ok' if o else 'NG'} |")
    P(f"VERDICT S0-c dY・dξ ({variant}): {'PASS' if okc and not nan else 'FAIL'}")
    res["S0-c(dY)"] = okc and not nan
    # Σ dY = 0
    P("\n## Σ_s dY_s = 0: 全節点・全成分で |Σ dY| ≤ 4ε Σ|dY_s| (res_1、gather 後)")
    sm = sum(g.astype(np.float64) for g in gY)
    sa = sum(np.abs(g.astype(np.float64)) for g in gY)
    ratio = np.abs(sm) / np.maximum(G.EPS32 * sa, 1e-300)
    viol = np.abs(sm) > 4 * G.EPS32 * sa
    P(f"| 最大 |Σ dY| | 最大 |Σ dY|/(ε Σ|dY|) | 超過 (節点×成分) | 超過の member 数 | 判定 |")
    P("| --- | --- | --- | --- | --- |")
    vm = sorted(set(msh.nmember[np.any(viol, axis=1)].tolist()))
    o = not viol.any()
    P(f"| {np.abs(sm).max():.3e} | {ratio.max():.2f} | {int(viol.sum())} | {vm} | {'ok' if o else 'NG'} |")
    P(f"VERDICT ΣdY ({variant}): {'PASS' if o and not nan else 'FAIL'}")
    res["SumdY"] = o and not nan
    # S0-e wall_y_eff
    if variant == "channel":
        P("\n## S0-e wall_y_eff (記録): res_0 / res_1 の周期 group 内差 (member − root)")
        wall = np.zeros(msh.nCells, bool); wall[msh.pc[msh.bnd_planes, 0]] = True
        wg = np.zeros(msh.nCells, bool); wg[msh.root[wall]] = True; wallg = wg[msh.root]
        mem = msh.root != np.arange(msh.nCells)
        P("| res | 区分 | member | 非零 | 最大 |差| | 最大 相対 |")
        P("| --- | --- | --- | --- | --- | --- |")
        nz = False
        for st in (0, 1):
            a = G.read_res(os.path.join(run, f"res_{st}.h5"), ["wall_y_eff"]).get("wall_y_eff")
            if a is None:
                P(f"| {st} | wall_y_eff が res に無い | - | - | - | - |"); continue
            a = a.astype(np.float64)
            d = np.abs(a - a[msh.root])
            for lab, m in (("壁∩継ぎ目", mem & wallg), ("壁なし継ぎ目", mem & ~wallg)):
                nzm = int((d[m] != 0).sum()); nz |= nzm > 0
                P(f"| {st} | {lab} | {int(m.sum())} | {nzm} | {d[m].max():.3e} | "
                  f"{(d[m] / np.maximum(np.abs(a[msh.root][m]), 1e-300)).max():.2e} |")
        P(f"S0-e wall_y_eff 記録 (channel): {'非零あり' if nz else 'すべて 0'}")
        res["S0-e(wall_y_eff)_nonzero"] = nz
    txt = "\n".join(out) + "\n"
    return txt, res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("variant", choices=("tgv", "jitter32", "channel", "box_slip"))
    ap.add_argument("--scratch", required=True)
    ap.add_argument("--species-db")
    ap.add_argument("--eval-only", action="store_true")
    a = ap.parse_args()
    run = os.path.join(a.scratch, f"s0y_{a.variant}")
    if not a.eval_only:
        run, merge = prepare(a.variant, a.scratch, a.species_db)
        os.environ["FORGE_DUMP_PREGATHER"] = os.path.join(run, "pregather")
        G.run_forge(run, expect_merge=merge)
        os.environ.pop("FORGE_DUMP_PREGATHER", None)
    txt, res = evaluate(a.variant, run)
    txt += "\n## まとめ\n" + "\n".join(f"{k}: {v}" for k, v in res.items()) + "\n"
    outp = os.path.join(HERE, f"S0Y_{a.variant}.txt")
    open(outp, "w").write(txt)
    print(txt)
    print("->", outp)


if __name__ == "__main__":
    main()
