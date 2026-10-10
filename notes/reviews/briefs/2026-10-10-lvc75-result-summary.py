# マスク 7/5 の腕の観測の要約 (記録用)
import csv, glob, h5py, numpy as np
runs = ["run_0570_lvc75_m7_a1", "run_0571_lvc75_m5_b1", "run_0572_lvc75_m7_a2", "run_0573_lvc75_m5_b2"]
for r in runs:
    rows = [x for x in csv.DictReader(open(r + "/residual_history.csv")) if x[list(x.keys())[2]] == "outer_begin"]
    c0 = list(rows[0].keys())[0]; d = {int(x[c0]): x for x in rows}
    print(r, "rms_ro", " ".join("%d:%.2e" % (s, float(d[s]["rms_ro"])) for s in (0, 20, 50, 100, 200, 300, 500, 1000, 1500, 1999) if s in d))
    print("   rms_roe", " ".join("%d:%.2e" % (s, float(d[s]["rms_roe"])) for s in (0, 20, 50, 100, 200, 300, 500, 1000, 1500, 1999) if s in d))
    s_ = sorted(d); ro = np.array([float(d[s]["rms_ro"]) for s in s_])
    m = int(np.argmax(np.where(np.isfinite(ro), ro, -1)))
    print(f"   rms_ro の最大 {ro[m]:.2e} (step {s_[m]})、最初に 10 倍を超えた step {next((s for s in s_ if float(d[s]['rms_ro']) > 10 * ro[0]), None)}")
with h5py.File(runs[0] + "/res_0.h5", "r") as h:
    X = h["MESH"]["COORD"][...].reshape(-1, 3); wd = h["VALUE"]["wall_dist"][...]
    T0, r0 = h["VALUE"]["T"][...], h["VALUE"]["ro"][...]
def top(r, step):
    with h5py.File(f"{r}/res_{step}.h5", "r") as h:
        T, ro = h["VALUE"]["T"][...], h["VALUE"]["ro"][...]
    rt = np.abs(T - T0) / T0; rr = np.abs(ro - r0) / r0
    it = np.argsort(-rt)[:5]
    cols = sorted(set((np.where(rt > 0.1)[0] // 121).tolist()))
    return (f"T の相対変化 >0.1: {int(np.sum(rt > 0.1))} 節点、最大 {rt.max():.2f}、列の範囲 {cols[:3]}…{cols[-3:]} (x {X[cols[0]*121,0]:.3f}〜{X[cols[-1]*121,0]:.3f})" if cols else f"T の相対変化 >0.1: 0、最大 {rt.max():.3f}") + \
           f"; 上位 {[(int(i), round(float(X[i,0]),3), round(float(X[i,1]),4), '%.1e' % wd[i], round(float(rt[i]),2)) for i in it[:3]]}"
for r in runs:
    for st in (100, 500, 1000, 2000):
        try: print(r, st, top(r, st))
        except Exception as e: pass
    nf = sorted(glob.glob(r + "/res_nan_*.h5"))
    if nf:
        with h5py.File(nf[0], "r") as h:
            bad = np.unique(np.concatenate([np.where(~np.isfinite(h["VALUE"][k][...]))[0] for k in h["VALUE"].keys() if h["VALUE"][k].shape == wd.shape]))
        cols = np.unique(bad // 121)
        print(r, "非有限:", bad.size, "節点、列", cols.tolist(), "x", [round(float(X[c*121, 0]), 3) for c in cols[[0, -1]]] if cols.size else None)
