import json, sys
import numpy as np
for fn in sys.argv[1:]:
    d = json.load(open(fn))
    rows = [r for r in d["rows"] if r["f"] is not None]
    nerr = len(d["rows"]) - len(rows)
    L = np.array([r["L_c"] for r in rows]); f = np.array([r["f"] for r in rows])
    A = np.vstack([L - L.mean(), np.ones_like(L)]).T
    (s, c), *_ = np.linalg.lstsq(A, f, rcond=None)
    res = f - (s * (L - L.mean()) + c)
    dL = np.diff(L); df = np.diff(f); loc = df / dL
    jump = df - s * dL
    i0 = np.argmin(np.abs(f))
    print(f"== {fn}: N {d['N']}  点 {len(rows)} (失敗 {nerr})  L_c {L.min():.4f}〜{L.max():.4f}  刻み {np.median(dL):.4g}")
    print(f"  直線の当てはめ: 傾き dx_F/dL_c = {s:.4f}、直線からの残差 rms {res.std():.3g}・最大 {np.abs(res).max():.3g} r_t")
    print(f"  隣どうしの局所の傾き: 中央値 {np.median(loc):.4f}、5〜95 % {np.percentile(loc,5):.3f}〜{np.percentile(loc,95):.3f}、最小 {loc.min():.3f}・最大 {loc.max():.3f}")
    print(f"  段 (隣どうしの差 − 傾き×刻み): 最大 |段| {np.abs(jump).max():.3g} r_t、|段| > 1e-3 の数 {int((np.abs(jump)>1e-3).sum())}、> 2.4e-3 (公差の半分) {int((np.abs(jump)>2.38e-3).sum())}")
    print(f"  |f| の最小 {abs(f[i0]):.3g} r_t (L_c {L[i0]:.5f})、f の符号が変わる区間 {int((np.sign(f[:-1])!=np.sign(f[1:])).sum())} 個、|f| ≤ 0.00476 の点 {int((np.abs(f)<=0.004762).sum())}")
