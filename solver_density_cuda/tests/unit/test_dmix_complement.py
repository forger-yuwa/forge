#!/usr/bin/env python3
"""混合平均拡散係数の分子の補数形 (plans/active/condensation-two-phase-transport.md §5.1 #3b) の判別 A/B (CFD 0 step)。
独立 CUDA 係数ハーネス test_dmix_complement.cu で、同じ入力に A = 現行分子 1−X_i / B = 補数形 Σ_{j≠i} X_j を評価し、
独立 double 参照 (test_transport_gas_phase.py の d_binary / d_mix と同じ Chapman–Enskog + Neufeld Ω(1,1)) と比べる。

  python3 solver_density_cuda/tests/unit/test_dmix_complement.py --forge BIN [--base-forge OLD_BIN] [--seed-run DIR] [--keep]

入力:
  (G2)  test_transport_gas_phase.py の G2 と同じ状態 (seed の外部 DB の [MIXDRY, H2O]、TP carrier 凝縮、3 組成 × g/Y_w 6 水準 ×
        T 11 水準、ρ 0.37)。forge の FORGE_TRANSPORT_PROBE (0 step) で読み戻したセル・ghost の float T・P・ρ・ρY・ρg を
        そのままハーネスに渡す (重複状態は 1 つにまとめる)。
  (bin) 二成分 [MIXDRY, H2O]: 各種を微量 X = 1e-2〜1e-8 にした組成 × T × P。
  (tri) 三成分 [MIXDRY, H2O, CO2 (LJ 3.941 Å / 195.2 K)] (異なる D_ij): 微量 2 種の組合せと中程度の混合。
  (pure) 純成分端点 (二成分・三成分) と n = 1。
合格 (事前固定; plan §5.1 #3b):
  B が G2 全点で独立 double 参照と相対 ≤1e-5、NaN/Inf・非正 0 件。二成分の非純成分点で両種 D が D₁₂ (double) と ≤1e-5。
  三成分も独立 double 参照と ≤1e-5。n = 1 は 0、純成分端点の純種は自己拡散フォールバック (A/B とも同じ値) で参照と ≤1e-5。
  判別: A は G2 で FAIL すること (FAIL しなければ元の FAIL を再現できていない)。
  ハーネスの忠実度: --base-forge (#3b 以前) の probe D が A と、--forge の probe D が B とビット一致 (情報として件数を出す)。
規約: [PASS]/[FAIL]、失敗があれば非ゼロ終了。
"""
import argparse, math, os, shutil, subprocess, sys, tempfile

import numpy as np
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
SOLVER = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)
import test_transport_gpu as tg  # noqa: E402
import test_transport_gas_phase as gp  # noqa: E402
from transport_reference import Reference  # noqa: E402

FAIL = 0
TOL = 1e-5


def check(ok, what):
    global FAIL
    print(("[PASS] " if ok else "[FAIL] ") + what, flush=True)
    if not ok:
        FAIL += 1


def rel(a, b):
    a, b = np.asarray(a, dtype=np.float64), np.asarray(b, dtype=np.float64)
    return np.abs(a - b) / np.maximum(np.maximum(np.abs(a), np.abs(b)), 1e-300)


def build(root, arch):
    exe = os.path.join(root, "test_dmix_complement")
    cmd = ["nvcc", "-std=c++17", "-O2", f"-arch={arch}", "--expt-relaxed-constexpr", "-I", SOLVER, "-o", exe,
           os.path.join(HERE, "test_dmix_complement.cu")]
    r = subprocess.run(["nice", "-n", "19"] + cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout + r.stderr)
        return None
    return exe


def harness(exe, root, species, states, name):
    """species: [(MW, sigma, eps)]、states: [(idx list, iw, T, P, ro, rog, roY list)]。戻り値 [(X, DA, DB)] (float32 配列)。"""
    pin, pout = os.path.join(root, name + ".in"), os.path.join(root, name + ".out")
    with open(pin, "w") as f:
        f.write(f"{len(species)} {len(states)}\n")
        for s in species:
            f.write("%.17g %.17g %.17g\n" % s)
        for idx, iw, T, P, ro, rog, roY in states:
            f.write(" ".join([str(len(idx))] + [str(i) for i in idx] + [str(iw)] + ["%.17g" % v for v in (T, P, ro, rog)]
                             + ["%.17g" % v for v in roY]) + "\n")
    r = subprocess.run([exe, pin, pout], capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(r.stdout + r.stderr)
    out = []
    for (idx, *_), line in zip(states, open(pout)):
        v = np.array([np.float32(t) for t in line.split()], dtype=np.float32)
        n = len(idx)
        out.append((v[:n], v[n:2 * n], v[2 * n:3 * n]))
    assert len(out) == len(states)
    return out


def X_of_Y(mw, Y):
    w = [y / m for y, m in zip(Y, mw)]
    return [v / sum(w) for v in w]


def summary(name, eA, eB, badA, badB):
    print(f"[INFO] {name}: points {eA.size}, A (1-X_i) max rel {eA.max():.3e}, B (sum_j X_j) max rel {eB.max():.3e}; "
          f"NaN/Inf/non-positive A {badA} B {badB}", flush=True)


def nbad(v):
    v = np.asarray(v, dtype=np.float64)
    return int(np.sum(~np.isfinite(v)) + np.sum(v <= 0))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--forge", required=True)
    ap.add_argument("--base-forge", default="")
    ap.add_argument("--seed-run", default="")
    ap.add_argument("--blocksize", type=int, default=256)
    ap.add_argument("--arch", default="sm_86")
    ap.add_argument("--keep", action="store_true")
    a = ap.parse_args()
    a.forge = os.path.abspath(a.forge)
    if a.base_forge:
        a.base_forge = os.path.abspath(a.base_forge)
    if not a.seed_run:
        for c in (os.path.join(tg.REPO, tg.SEED_REL), os.path.join(os.path.dirname(tg.REPO), "forge", tg.SEED_REL)):
            if os.path.exists(os.path.join(c, "nozzle.h5")):
                a.seed_run = c
                break
    if not a.seed_run:
        print("[FAIL] seed run with nozzle.h5 not found (--seed-run)")
        return 1
    root = tempfile.mkdtemp(prefix="dmix_complement_")
    print(f"work dir: {root}\nseed: {a.seed_run}\nforge: {a.forge}\nbase: {a.base_forge or '(none)'}", flush=True)
    exe = build(root, a.arch)
    if exe is None:
        check(False, "harness build")
        return 1
    seed_db = yaml.safe_load(open(os.path.join(a.seed_run, "species_db.yaml")))
    spA = ["MIXDRY", "H2O"]
    trA = {"MIXDRY": "kinetic", "H2O": "custom:h2o_iapws_cea_v1"}
    refA = Reference(spA, trA, seed_db)
    dsp = [{"MW": float(seed_db[n]["MW"]), "sigma": float(seed_db[n]["LJ_sigma"]), "eps": float(seed_db[n]["LJ_eps_kB"])} for n in spA]
    dsp.append({"MW": 0.0440095, "sigma": 3.941, "eps": 195.2})   # 三成分の 3 種目 (CO2 の LJ; D_ij が他の 2 対と異なる)
    hsp = [(d["MW"], d["sigma"], d["eps"]) for d in dsp]
    mw = [d["MW"] for d in dsp]
    iw = 1

    # ------------------------------------------------------------------ (G2) 既存 G2 の状態
    R = tg.Runner(a, root)
    comps = [[0.99, 0.01], [0.94, 0.06], [0.6, 0.4]]
    fr = [0.0, 0.3, 0.9, 0.99, 1.0, 1.2]
    states = [(T, 0.37, gp.y_of_x(refA, X), f * gp.y_of_x(refA, X)[iw]) for X in comps for f in fr for T in gp.T_LIST]
    probes = {}
    for which, exe_f in (("new", a.forge), ("old", a.base_forge)):
        if not exe_f:
            continue
        d = gp.make(R, f"g2_{which}", spA, trA, True, keep_db=True)
        out, _, err = gp.probe(R, exe_f, d, states, 2, True)
        if out is None:
            check(False, f"G2 {which}: probe failed: {err}")
            return 1
        probes[which] = out
    out = probes["new"]
    uniq, cell_key = {}, []
    for p, q in enumerate(out["passes"]):
        for i in range(out["h"]["nCells_all"]):
            key = (np.float32(q["T"][i]).tobytes(), np.float32(q["P"][i]).tobytes(), np.float32(q["ro"][i]).tobytes(),
                   np.float32(q["rog"][i]).tobytes(), np.asarray(q["roY"][:, i], np.float32).tobytes())
            if key not in uniq:
                uniq[key] = (float(q["T"][i]), float(q["P"][i]), float(q["ro"][i]), float(q["rog"][i]),
                             [float(v) for v in q["roY"][:, i]])
            cell_key.append((p, i, key))
    keys = list(uniq)
    hin = [([0, 1], iw, *uniq[k][:4], uniq[k][4]) for k in keys]
    res = harness(exe, root, hsp, hin, "g2")
    eA, eB, Dref_all, DA_all, DB_all, Xg_all = [], [], [], [], [], []
    for k, (X, DA, DB) in zip(keys, res):
        T, P, ro, rog, roY = uniq[k]
        Xg = refA.X_from_Y(gp.gas_Y(roY, rog, ro, iw))
        Dr = gp.d_mix(dsp[:2], Xg, T, P)
        eA.append(rel(DA, Dr)); eB.append(rel(DB, Dr))
        Dref_all.append(Dr); DA_all.append(DA); DB_all.append(DB); Xg_all.append(Xg)
    eA, eB = np.array(eA), np.array(eB)
    DA_all, DB_all, Xg_all = np.array(DA_all), np.array(DB_all), np.array(Xg_all)
    badA, badB = nbad(DA_all), nbad(DB_all)
    summary(f"G2 ({len(keys)} unique states from {len(cell_key)} cells+ghosts x passes)", eA, eB, badA, badB)
    worst = np.unravel_index(np.argmax(eA), eA.shape)
    kw = keys[worst[0]]
    print(f"[INFO] G2 A worst point: species {worst[1]}, 1-X_s {1.0 - Xg_all[worst]:.3e}, T {uniq[kw][0]:.6g} K, P {uniq[kw][1]:.6g} Pa, "
          f"A rel {eA[worst]:.3e}, B rel {eB[worst]:.3e}", flush=True)
    check(eB.max() <= TOL and badB == 0,
          f"G2 B (complement): max rel vs independent double reference {eB.max():.3e} (<= {TOL:g}); NaN/Inf/non-positive {badB}")
    check(eA.max() > TOL,
          f"G2 A (current 1-X_i) discrimination: max rel {eA.max():.3e} (> {TOL:g}: the original G2 FAIL is reproduced)")
    # ハーネスの忠実度と旧バイナリとの比 (情報)
    kidx = {k: n for n, k in enumerate(keys)}
    nA_ = out["h"]["nCells_all"]
    for which, arr in (("new", DB_all), ("old", DA_all)):
        if which not in probes:
            continue
        diff = 0
        for (p, i, k) in cell_key:
            if not np.array_equal(probes[which]["passes"][p]["D"][:, i], arr[kidx[k]]):
                diff += 1
        print(f"[INFO] harness fidelity: forge probe D ({which} binary) vs harness {'B' if which == 'new' else 'A'} "
              f"bit-identical in {len(cell_key) - diff}/{len(cell_key)} cells+ghosts", flush=True)
    if "old" in probes:
        r_new_old = []
        for (p, i, k) in cell_key:
            r_new_old.append(probes["new"]["passes"][p]["D"][:, i].astype(np.float64) / probes["old"]["passes"][p]["D"][:, i])
        r_new_old = np.array(r_new_old)
        dev = np.abs(r_new_old - 1.0)
        print(f"[INFO] new/old forge probe D over G2 cells+ghosts: max |D_new/D_old - 1| {dev.max():.3e}, "
              f"bit-identical {int(np.sum(dev == 0))}/{dev.size}", flush=True)
    # 代表組成 (気相 X_H2O) ごとの B/A (旧バイナリとの比) — 気相組成で束ねる
    print("[INFO] D_B/D_A by gas-phase X_H2O (G2 states; per species max |ratio-1|):", flush=True)
    xw = Xg_all[:, 1]
    for lo, hi in ((0.0, 1e-8), (1e-8, 1e-6), (1e-6, 1e-4), (1e-4, 1e-2), (1e-2, 0.1), (0.1, 1.0)):
        m = (xw > lo) & (xw <= hi) if lo > 0 else (xw <= hi)
        if m.any():
            rr = np.abs(DB_all[m].astype(np.float64) / DA_all[m] - 1.0)
            print(f"         X_H2O in ({lo:g}, {hi:g}]: n {int(m.sum())}, MIXDRY {rr[:, 0].max():.3e}, H2O {rr[:, 1].max():.3e}", flush=True)

    # ------------------------------------------------------------------ (bin) 二成分の微量成分 1e-2〜1e-8
    TP = [(T, P) for T in (200.0, 300.0, 1000.0, 3000.0) for P in (1.0e3, 1.0e5, 1.0e6)]
    hin, meta = [], []
    for minor in (0, 1):
        for e in range(2, 9):
            xm = 10.0 ** (-e)
            X = [1.0 - xm, xm] if minor == 1 else [xm, 1.0 - xm]
            Y = [float(np.float32(v)) for v in X_of_Y([1 / m for m in mw[:2]], X)]   # Y = X M / Σ X M
            for T, P in TP:
                hin.append(([0, 1], -1, T, P, 1.0, 0.0, Y))
                meta.append((minor, xm, T, P))
    res = harness(exe, root, hsp, hin, "bin")
    eA = np.array([rel(r[1], [gp.d_binary(dsp[0], dsp[1], T, P)] * 2) for r, (_, _, T, P) in zip(res, meta)])
    eB = np.array([rel(r[2], [gp.d_binary(dsp[0], dsp[1], T, P)] * 2) for r, (_, _, T, P) in zip(res, meta)])
    badA, badB = nbad([r[1] for r in res]), nbad([r[2] for r in res])
    summary(f"binary trace (x_minor 1e-2..1e-8, both species minor, {len(TP)} (T,P))", eA, eB, badA, badB)
    for e in range(2, 9):
        m = np.array([abs(mt[1] - 10.0 ** (-e)) < 1e-30 for mt in meta])
        print(f"         x_minor 1e-{e}: A max rel {eA[m].max():.3e}, B max rel {eB[m].max():.3e}", flush=True)
    check(eB.max() <= TOL and badB == 0,
          f"binary B: both species D vs D12 (double) max rel {eB.max():.3e} (<= {TOL:g}) at x_minor 1e-2..1e-8; NaN/Inf/non-positive {badB}")

    # ------------------------------------------------------------------ (tri) 三成分・異なる D_ij
    trace = [1e-2, 1e-4, 1e-6, 1e-8]
    comps3 = []
    for dom in range(3):
        for t1 in trace:
            for t2 in trace:
                X = [0.0, 0.0, 0.0]
                o = [s for s in range(3) if s != dom]
                X[o[0]], X[o[1]] = t1, t2
                X[dom] = 1.0 - t1 - t2
                comps3.append(X)
    comps3 += [[0.5, 0.3, 0.2], [0.2, 0.5, 0.3], [1 / 3, 1 / 3, 1 / 3], [0.7, 0.29, 0.01]]
    hin, meta = [], []
    for X in comps3:
        Y = [float(np.float32(v)) for v in X_of_Y([1 / m for m in mw], X)]
        for T, P in TP:
            hin.append(([0, 1, 2], -1, T, P, 1.0, 0.0, Y))
            meta.append((X_of_Y(mw, Y), T, P))
    res = harness(exe, root, hsp, hin, "tri")
    eA = np.array([rel(r[1], gp.d_mix(dsp, Xr, T, P)) for r, (Xr, T, P) in zip(res, meta)])
    eB = np.array([rel(r[2], gp.d_mix(dsp, Xr, T, P)) for r, (Xr, T, P) in zip(res, meta)])
    badA, badB = nbad([r[1] for r in res]), nbad([r[2] for r in res])
    summary(f"ternary (distinct D_ij, {len(comps3)} compositions x {len(TP)} (T,P))", eA, eB, badA, badB)
    check(eB.max() <= TOL and badB == 0,
          f"ternary B: vs independent double reference max rel {eB.max():.3e} (<= {TOL:g}); NaN/Inf/non-positive {badB}")

    # ------------------------------------------------------------------ (pure) 純成分端点と n = 1
    hin, meta = [], []
    for idx in ([0, 1], [0, 1, 2]):
        for p in range(len(idx)):
            Y = [1.0 if s == p else 0.0 for s in range(len(idx))]
            for T, P in TP:
                hin.append((idx, -1, T, P, 1.0, 0.0, Y))
                meta.append((idx, p, T, P))
    for s in range(3):
        for T, P in TP:
            hin.append(([s], -1, T, P, 1.0, 0.0, [1.0]))
            meta.append(([s], 0, T, P))
    res = harness(exe, root, hsp, hin, "pure")
    worst, worst_self, n1_bad, same_bad, bad, dAB = 0.0, 0.0, 0, 0, 0, 0.0
    for r, (idx, p, T, P) in zip(res, meta):
        DA, DB = r[1], r[2]
        if len(idx) == 1:
            n1_bad += int(DA[0] != 0.0 or DB[0] != 0.0)
            continue
        sps = [dsp[i] for i in idx]
        Xr = [1.0 if s == p else 0.0 for s in range(len(idx))]
        worst = max(worst, float(rel(DB, gp.d_mix(sps, Xr, T, P)).max()))
        worst_self = max(worst_self, float(rel(DB[p], gp.d_binary(sps[p], sps[p], T, P))))
        same_bad += int(not np.array_equal(DA, DB))
        dAB = max(dAB, float(rel(DA, DB).max()))
        bad += nbad(DB)
    check(worst <= TOL and worst_self <= TOL and n1_bad == 0 and bad == 0,
          f"pure endpoints: B vs reference max rel {worst:.3e}, pure species = self-diffusion fallback max rel {worst_self:.3e} (<= {TOL:g}); "
          f"n==1 nonzero {n1_bad}; NaN/Inf/non-positive {bad}")
    # 純種は A/B 同じ自己拡散フォールバック。他の種の分子は A = 1−0 = 1、B = X_pure (float の X_from_Y で 1 から 1 ulp ずれうる)
    print(f"[INFO] pure endpoints (n >= 2): A and B differ in {same_bad} of {len(meta) - 3 * len(TP)} states, "
          f"max rel {dAB:.3e} (non-pure species: A numerator 1-0 = 1, B numerator = float X_pure)", flush=True)

    if not a.keep:
        shutil.rmtree(root, ignore_errors=True)
    print("ALL PASS" if FAIL == 0 else f"FAILED: {FAIL}")
    return 0 if FAIL == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
