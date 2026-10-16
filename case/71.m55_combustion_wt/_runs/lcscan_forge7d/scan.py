"""from_length の x_F(L_c) の形を見る掃引 (カード B、M6 の feature/card-local-mesh 2871a489 のコード)。
L_c を固定 (Lc_mode explicit) して design_chain を回し、x_F − L_total を記録する。"""
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
    if os.environ.get("ALIGN_XE"):
        # 診断の A/B: x_E に最も近い軸の節点を x_E に乗せる (上流の節点 x_m〜x_end を 2 区間の線形の写像で伸縮;
        # 間隔の変化は ±(半間隔)/(x_E − x_m) ≈ ±0.15 %)。x_A = x0 (カード B の実測 0.3387028221366673)
        from forge_design.geometry import moc_inverse as MI
        orig = MI._axis_grid
        xE = 0.3387028221366673 + float(Lc)
        def aligned(x0, x_end, n, dx0=None, *a, **k):
            g = np.array(orig(x0, x_end, n, dx0, *a, **k), dtype=float)
            j = int(np.argmin(np.abs(g - xE)))
            m = int(np.argmin(np.abs(g - 0.5 * xE)))
            xm, xj = g[m], g[j]
            lo = (g >= xm) & (g <= xj); hi = g > xj
            g2 = g.copy()
            g2[lo] = xm + (g[lo] - xm) * (xE - xm) / (xj - xm)
            g2[hi] = xE + (g[hi] - xj) * (x_end - xE) / (x_end - xj)
            g2[j] = xE
            return g2
        MI._axis_grid = aligned
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
            print(json.dumps({k: r[k] for k in ("L_c", "f", "t", "err")}), flush=True)
    rows.sort(key=lambda r: r["L_c"])
    json.dump({"problem": PROB, "N": N, "rows": rows}, open(out, "w"), indent=1)
