#!/usr/bin/env python3
r"""C3X 一様壁: forge と SU2 の**遷移モデルどうし**を同一メッシュ・同一節点で比べる (plan turbulence-transition-lm2009 §5.1 #5b)。

壁熱流束は forge が `iface_q_eff`、SU2 が `Heat_Flux` (vol_solution.vtu の壁節点)。比べる量:
  * 領域別 (正圧面 / 負圧面 s/S<0.25 / 負圧面 0.25–0.87) の平均 $h$ と forge/SU2 の比
  * 実測に対する平均偏差 (`compare_h.py` と同じ規約: 熱電対位置でなく全壁節点、$s/S\le0.87$)
  * 負圧面の遷移位置: $s/S>0.15$ で $h$ が最小になる位置 (遷移直前) と、その下流で (最小+最大)/2 を上向きに横切る位置

usage: compare_su2_lm.py FORGE_RUN SU2_DIR [--Tw 566] [--run run108]
"""
import argparse, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).parent))
import report_figures as rf
from compare_h import TABLES, H0
from scipy.spatial import cKDTree

ap = argparse.ArgumentParser(); ap.add_argument("forge"); ap.add_argument("su2"); ap.add_argument("--Tw", type=float, default=566.0); ap.add_argument("--run", default="run108")
a = ap.parse_args()
s, ss, V, C, step = rf.wall(Path(a.forge).resolve())
XY, QS = rf.read_su2(Path(a.su2) / "vol_solution.vtu")
d, idx = cKDTree(XY).query(C); assert d.max() == 0.0, d.max()
Tg = 786.0; hf = np.array(V["iface_q_eff"]) / (Tg - a.Tw); hs = QS[idx] / (Tg - a.Tw)
rows = [r for r in TABLES[a.run]["rows"] if r[3] is not None]; sd = np.array([r[0] for r in rows]); hd = np.array([r[3] for r in rows]) * H0; i0 = int(np.argmin(sd))
exp = {"PS": (sd[:i0 + 1][::-1], hd[:i0 + 1][::-1]), "SS": (sd[i0:], hd[i0:])}
he = np.array([np.interp(s[k], *exp["SS" if ss[k] else "PS"]) for k in range(len(s))])
print(f"forge {a.forge} (step {step})  vs  SU2 {a.su2}   [uniform Tw {a.Tw} K, {len(s)} coincident wall nodes]")
print(f"{'region':<26}{'h forge':>10}{'h SU2':>10}{'forge/SU2-1':>13}{'forge vs exp':>14}{'SU2 vs exp':>12}")
for name, m in (("pressure side", (~ss) & (s <= 0.87)), ("suction s/S<0.25", ss & (s < 0.25)), ("suction 0.25-0.87", ss & (s >= 0.25) & (s <= 0.87)), ("all (s/S<=0.87)", s <= 0.87)):
    print(f"{name:<26}{hf[m].mean():10.1f}{hs[m].mean():10.1f}{hf[m].mean()/hs[m].mean()-1:+13.1%}{((hf[m]-he[m])/he[m]).mean():+14.1%}{((hs[m]-he[m])/he[m]).mean():+12.1%}")
for lab, h in (("forge", hf), ("SU2", hs)):
    m = ss & (s > 0.15) & (s < 0.6); o = np.argsort(s[m]); x, y = s[m][o], h[m][o]
    k = 41; ys = np.convolve(np.pad(y, k // 2, mode="edge"), np.ones(k) / k, "valid") if len(y) > k else y
    i = int(np.argmin(ys)); mid = 0.5 * (ys[i] + ys[i:].max()); j = i + int(np.argmax(ys[i:] >= mid))
    print(f"  {lab:<6} suction side: h minimum at s/S = {x[i]:.3f} (h/h0 {ys[i]/H0:.3f}), crosses half-rise at s/S = {x[j]:.3f}, plateau max h/h0 {ys[i:].max()/H0:.3f}")
