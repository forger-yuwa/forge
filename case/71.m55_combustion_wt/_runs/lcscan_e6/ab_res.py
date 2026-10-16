"""V9 (b′) 解像度の A/B: カード B の入力で n_axis_inv だけ 2400 / 4800 を変え、7 点の L_c で x_F・終端の Mach 線との整合・切り口・壁の検査を見る。
plan tooling-nozzle-card-local-mesh §6 V9 (b′) (事前登録 2026-10-16)。"""
import os, sys, json, time
for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[k] = "1"
import numpy as np

PROB, N, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
LCS = [37.170, 37.1885, 37.2070, 37.2075, 37.2255, 37.2440, 37.260]

def one(Lc):
    from forge_design.probdef import load_problem
    from forge_design.evaluate import runner_axismach as RA
    p = load_problem(PROB)
    p.geometry["Lc_mode"] = "explicit"
    p.geometry["n_axis_inv"] = N
    p.dv["L_c"] = {"value": float(Lc)}
    t0 = time.time()
    try:
        d = RA.design_chain(p)
    except Exception as e:
        return {"L_c": Lc, "err": repr(e)[:400], "t": time.time() - t0}
    ex, law = d["exit"], d["law"]
    x_F, r_F, x_E = float(ex["x_F"]), float(ex["r_F"]), float(law.x_E)
    M_d = float(law(x_E))
    wi = d["wall_inv"]
    return {"L_c": Lc, "err": None, "t": time.time() - t0, "x_F": x_F, "r_F": r_F, "x_E": x_E, "M_d": M_d,
            "term_line": x_F - x_E - r_F * np.sqrt(M_d * M_d - 1.0),
            "cut_dx": x_F - float(wi[-1, 0]), "cut_dr": r_F - float(wi[-1, 1]),
            "qa_violations": list(d["qa"]["violations"]), "wall_validate": list(d["wall"].validate()),
            "x_E_bracket": ex.get("x_E_bracket"), "x_E_weight": ex.get("x_E_weight"), "x_F_nodes": ex.get("x_F_nodes")}

if __name__ == "__main__":
    nproc = int(sys.argv[4]) if len(sys.argv) > 4 else 1
    if nproc > 1:
        from multiprocessing import Pool
        with Pool(nproc) as pool:
            rows = pool.map(one, LCS)
    else:
        rows = [one(L) for L in LCS]
    json.dump({"problem": PROB, "N": N, "rows": rows}, open(out, "w"), indent=1)
    for r in rows:
        print(json.dumps({k: r.get(k) for k in ("L_c", "x_F", "term_line", "cut_dx", "qa_violations", "err", "t")}))
