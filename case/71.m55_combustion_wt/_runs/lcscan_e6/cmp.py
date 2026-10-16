import json, sys
import numpy as np
def load(f):
    r = [x for x in json.load(open(f))["rows"] if x["err"] is None]
    return np.array([x["L_c"] for x in r]), np.array([x["f"] for x in r])
for name, f in zip(("現行 (最寄り)", "補間"), sys.argv[1:3]):
    L, F = load(f)
    dL = np.diff(L); s = np.diff(F) / dL
    print(f"{name}: 点 {len(L)}  f {F.min():+.5f}〜{F.max():+.5f}  局所の傾き 中央値 {np.median(s):.4f} 最小 {s.min():.4f} 最大 {s.max():.4f}  "
          f"隣の差の最大 {np.abs(np.diff(F)).max():.5f} r_t (L_c の刻み {dL[0]:.4g})")
    fit = np.polyfit(L, F, 1); res = F - np.polyval(fit, L)
    print(f"   直線の当てはめ: 傾き {fit[0]:.5f}  残差 {np.abs(res).max():.2e} r_t (最大)、{res.std():.2e} (標準偏差)")
