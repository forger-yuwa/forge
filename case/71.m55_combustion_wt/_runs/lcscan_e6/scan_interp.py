"""#6e の試作: 終端 C⁺ を x_E を挟む 2 本の網の C⁺ の x_F の線形補間にする (MOC の網・壁は変えない)。
forge-7d の scan.py と同じ形で x_F(L_c) を掃く。INTERP=0 なら現行 (最寄りの節点)。"""
import os, sys, json, time
for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[k] = "1"
import numpy as np
from multiprocessing import Pool

PROB = sys.argv[1]; N = int(sys.argv[2]); lo, hi, n = float(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5])
nproc = int(sys.argv[6]); out = sys.argv[7]

def one(Lc):
    from forge_design.probdef import load_problem
    from forge_design.evaluate import runner_axismach as RA
    p = load_problem(PROB)
    p.geometry["Lc_mode"] = "explicit"
    p.geometry["n_axis_inv"] = N
    p.dv["L_c"] = {"value": float(Lc)}
    L_tgt = float(p.dv["L_total"]["value"])
    t0 = time.time()
    try:
        d = RA.design_chain(p)
        xF = float(d["exit"]["x_F"])
        return {"L_c": Lc, "x_F": xF, "f": xF - L_tgt, "x_E": float(d["law"].x_E), "t": time.time() - t0, "err": None,
                "exit": {k: (float(v) if isinstance(v, (int, float, np.floating)) else None) for k, v in d["exit"].items()
                         if isinstance(v, (int, float, np.floating))}}
    except Exception as e:
        return {"L_c": Lc, "x_F": None, "f": None, "t": time.time() - t0, "err": repr(e)[:300]}

if __name__ == "__main__":
    import forge_design.evaluate.runner_axismach as RA
    print("module", RA.__file__, flush=True)
    Ls = list(np.linspace(lo, hi, n))
    rows = []
    with Pool(nproc) as pool:
        for r in pool.imap_unordered(one, Ls):
            rows.append(r)
    rows.sort(key=lambda r: r["L_c"])
    json.dump({"problem": PROB, "N": N, "rows": rows}, open(out, "w"), indent=1)
    print("done", len(rows), sum(1 for r in rows if r["err"]), "errors")
