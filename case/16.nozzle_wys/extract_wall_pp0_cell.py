import sys, h5py, numpy as np, yaml, os
run, res, out = sys.argv[1], sys.argv[2], sys.argv[3]; p0 = 59070.0
mesh = yaml.safe_load(open(run+'/solverConfig.yaml'))['mesh']['meshFileName']
m = h5py.File(run+'/'+mesh, 'r'); r = h5py.File(run+'/'+res, 'r')
P = np.array(r['VALUE/P']); wd = np.array(r['VALUE/wall_dist']); n = len(P)
c = np.array(m['CELLS/centCoords']).reshape(-1, 3)[:n]
G = np.array(r['VALUE/g_0']) if 'VALUE/g_0' in r else None
xs = np.round(c[:, 0], 6); ux = np.unique(xs); xx, pp, gg = [], [], []
for xv in ux:
    cand = np.where((xs == xv) & (c[:, 1] > 0))[0]
    if len(cand) == 0: continue
    j = cand[np.argmin(wd[cand])]; xx.append(c[j, 0]*1e3); pp.append(P[j]/p0); gg.append(G[j] if G is not None else 0.0)
cols = {'x_mm': np.array(xx), 'contour_wall': np.array(pp)}
if G is not None: cols['contour_g'] = np.array(gg)
np.savetxt(out, np.column_stack(list(cols.values())), delimiter=',', header=','.join(cols), comments='')
print('wrote', out, len(xx), 'stations (cell: first cell off the contour wall)')
