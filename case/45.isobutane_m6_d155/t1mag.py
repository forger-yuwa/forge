import h5py, numpy as np, csv
C45 = "/home/ubuntu/forge-wallfit/case/45.isobutane_m6_d155"; C56 = "/home/ubuntu/forge-wallfit/case/56.gap_tp1187"
def vals(r, n):
    with h5py.File(f"{r}/res_{n}.h5") as h: return {k: np.asarray(h["VALUE/" + k][:], dtype=np.float64) for k in h["VALUE"].keys()}
def hist(r):
    rows = [x for x in csv.DictReader(open(f"{r}/residual_history.csv")) if x["phase"] == "outer_end"]
    cols = [c for c in rows[0] if c.startswith("rms_") and not c.startswith("rms_dq")]
    return np.array([[float(x[c]) for c in cols] for x in rows]), cols
for tag, base, n, (o, nw, o2) in (("FP64 case/45", C45, 20, ("run_0390_fg_id64_old", "run_0391_fg_id64_new", "run_0392_fg_id64_old2")),
                                  ("float case/56", C56, 20, ("run_0079_fg_id32_old", "run_0080_fg_id32_new", "run_0081_fg_id32_old2"))):
    A, B, A2 = vals(f"{base}/{o}", n), vals(f"{base}/{nw}", n), vals(f"{base}/{o2}", n)
    print(f"== {tag}: 20 step 後の場の差 (最大 |差| / 最大 |値|)")
    for k in ("ro", "roUx", "roUy", "roe", "roK", "roOmega", "P", "T"):
        if k not in A: continue
        s = np.max(np.abs(A[k]))
        print(f"   {k:8s} 新−旧 {np.max(np.abs(B[k]-A[k]))/s:.2e}   旧−旧 (再実行) {np.max(np.abs(A2[k]-A[k]))/s:.2e}")
    ha, cols = hist(f"{base}/{o}"); hb, _ = hist(f"{base}/{nw}"); ha2, _ = hist(f"{base}/{o2}")
    rn = np.max(np.abs(hb - ha) / np.maximum(np.abs(ha), 1e-300)); rr = np.max(np.abs(ha2 - ha) / np.maximum(np.abs(ha), 1e-300))
    print(f"   残差履歴 (全列・全 step) の相対差の最大: 新−旧 {rn:.2e}   旧−旧 {rr:.2e}")
