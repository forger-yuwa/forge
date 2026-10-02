# limiter-inlet-column-oscillation §5.1 #4 (H-R): face-level recomputation of the scaled Venkatakrishnan psi at the inlet columns
# from the every-step window runs. Reproduces limiter_d.cu (scaled=1, node, convMethod 1): delta_m = grad.(x_j-x_i)/2,
# delta+- = qmax/qmin over {i} U edge neighbours - q_i, all divided by q_ref; eps2 = (K h / L)^3, h = sqrt(volume) (planar 2D).
# Neighbours = primal quad edges (= internal dual faces). Validation: recomputed psi vs the solver's limiter_* output.
import h5py, numpy as np, glob, sys
C = "/home/ubuntu/forge-species-tp/case/16.nozzle_wys/"
run = sys.argv[1]
K, L = 0.05, 0.1570673715
REF = {"ro": 0.5279511755, "Ux": 330.8270384, "Uy": 330.8270384, "P": 43128.33253}
GR = {"ro": "drod", "Ux": "dUxd", "Uy": "dUyd", "P": "dPd"}
fs = sorted(glob.glob(C + run + "/res_*.h5"), key=lambda s: int(s.split("_")[-1][:-3]))
fs = [f for f in fs if not f.endswith("/res_0.h5")]
f0 = h5py.File(fs[0]); xyz = np.array(f0["MESH/COORD"][...]).reshape(-1, 3).astype(np.float64)
cn = np.array(f0["MESH/CONNE"][...]); edges = set(); k = 0
while k < cn.size:
    t = cn[k]; n = {5: 4, 4: 3}[t]; nd = cn[k + 1:k + 1 + n]
    for a in range(n):
        u, v = nd[a], nd[(a + 1) % n]; edges.add((min(u, v), max(u, v)))
    k += 1 + n
E = np.array(sorted(edges)); N = xyz.shape[0]
nb = [[] for _ in range(N)]
for u, v in E: nb[u].append(v); nb[v].append(u)
x = xyz[:, 0] * 1e3
S = np.where((np.abs(x + 60) < 0.05) | (np.abs(x + 59.37) < 0.05))[0]
print(run, "snapshots", len(fs), "edges", len(E), "S nodes", S.size)
stats = {q: {"match": [], "bind": [], "dp_ulp": [], "dm": [], "dp": [], "e2": [], "psi": []} for q in REF}
for f in fs:
    V = h5py.File(f)["VALUE"]
    vol = np.array(V["volume"][...]).ravel().astype(np.float64); h = np.sqrt(vol); e2 = (K * h / L) ** 3
    for q in REF:
        Q = np.array(V[q][...]).ravel().astype(np.float32)
        gx = np.array(V[GR[q] + "x"][...]).ravel(); gy = np.array(V[GR[q] + "y"][...]).ravel()
        lim = np.array(V["limiter_" + q][...]).ravel()
        ps, bd, dpu, dms, dps = [], [], [], [], []
        for i in S:
            js = np.array(nb[i]); qi = Q[i]
            qmax = max(qi, Q[js].max()); qmin = min(qi, Q[js].min())
            dx = 0.5 * (xyz[js, 0] - xyz[i, 0]); dy = 0.5 * (xyz[js, 1] - xyz[i, 1])
            dm = (gx[i] * dx + gy[i] * dy) / REF[q]
            dp = np.where(dm > 0, (qmax - qi) / REF[q], (qmin - qi) / REF[q])
            num = dp * dp + e2[i] + 2 * dm * dp; den = dp * dp + 2 * dm * dm + dp * dm + e2[i]
            psi_f = np.where(den > 0, num / den, 1.0); a = np.argmin(psi_f)
            ps.append(min(max(psi_f[a], 0), 1)); bd.append(js[a]); dms.append(dm[a]); dps.append(dp[a])
            dpu.append(abs(dp[a]) / (np.spacing(np.float32(max(abs(qmax), abs(qmin)))) / REF[q]))
        ps = np.array(ps)
        stats[q]["match"].append(np.abs(ps - lim[S])); stats[q]["bind"].append(bd); stats[q]["dp_ulp"].append(dpu)
        stats[q]["dm"].append(dms); stats[q]["dp"].append(dps); stats[q]["e2"].append(e2[S]); stats[q]["psi"].append(ps)
for q in REF:
    s = {k2: np.array(v) for k2, v in stats[q].items()}
    sw = np.mean(s["bind"][1:] != s["bind"][:-1], axis=0)        # 各節点の制限面の入れ替わり率 (step ごと)
    psd = s["psi"].std(0); w = np.argsort(psd)[::-1][:5]
    print(f"[{q}] recomputed vs solver psi: max|d| {s['match'].max():.2e} p99 {np.quantile(s['match'],0.99):.2e} | "
          f"binding-face switch rate p50 {np.median(sw):.3f} max {sw.max():.3f} | psi std p50 {np.median(psd):.3e} max {psd.max():.3e}")
    for j in w:
        i = S[j]
        print(f"   node x {x[i]:7.2f} y {xyz[i,1]*1e3:7.3f}: psi mean {s['psi'][:,j].mean():.3f} std {psd[j]:.3f} switch {sw[j]:.2f} | "
              f"|dm| p50 {np.median(np.abs(s['dm'][:,j])):.2e} |dp| p50 {np.median(np.abs(s['dp'][:,j])):.2e} min {np.abs(s['dp'][:,j]).min():.2e} "
              f"eps2 {s['e2'][0,j]:.2e} (eps {np.sqrt(s['e2'][0,j]):.2e}) | dp in ulps p50 {np.median(s['dp_ulp'][:,j]):.1f} min {s['dp_ulp'][:,j].min():.1f}")
