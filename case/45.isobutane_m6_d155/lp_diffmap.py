"""§6.47 の記録 (判定ではない): pRef の下で float と FP64 の差が x のどこに残るか。A で回す。
各節点の「自分の step 0 からの変化」の (float − FP64) を FP64 の step 0 で割る。float は Bp1・Bp2 (pRef)、F0 (pRef なし)、FP64 は Ap (pRef)。
step 20,000・25,000・30,000。領域ごと (入口 x<−8・収縮 −8〜−2・スロート −2〜2・超音速 2〜20・超音速 >20) に、列の区分 (中心 j<100・壁側 j≥100) の
|ΔT|/T の中央値と 90 % 点、|ΔU|/U の中央値、θ_r の差 (x の分布を区分で平均) を出す。"""
import sys
from pathlib import Path
import numpy as np, h5py
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cold_xcheck as XC  # noqa: E402
with h5py.File(HERE / "run_0487_tr_fp64/nozzle.h5", "r") as h: xy = np.asarray(h["MESH/COORD"][:], float).reshape(-1, 3)
ni, nj, S = XC.mesh_info(len(xy)); X = (xy[:, 0] / S).reshape(ni, nj); J = np.tile(np.arange(nj), (ni, 1))
yb_x, yb = XC.common_yb()
def load(run, s):
    xy_, ro, ux, uy, T, k = XC.load_forge(HERE / run, s)
    o = XC.reduce_fields(xy_, ro, ux, uy, T, k, False, ni, nj, S, yb_x, yb)
    return dict(T=T.reshape(ni, nj), U=np.hypot(ux, uy).reshape(ni, nj), th=o["theta_r"], x=o["x"])
R = {"Ap": "run_0712_lp_ap", "Bp1": "run_0710_lp_bp1", "Bp2": "run_0711_lp_bp2", "F0": "run_0709_lp_f0"}
D0 = {k: load(r, 0) for k, r in R.items()}
REG = [("入口 x<-8", X < -8), ("収縮 -8..-2", (X >= -8) & (X < -2)), ("スロート -2..2", (X >= -2) & (X < 2)), ("超音速 2..20", (X >= 2) & (X < 20)), ("超音速 >20", X >= 20)]
for s in (20000, 30000):
    Ds = {k: load(r, s) for k, r in R.items()}
    print(f"== step {s} (float − FP64 + pRef の、自分の step 0 からの変化の差)")
    for f in ("Bp1", "Bp2", "F0"):
        dT = ((Ds[f]["T"] - D0[f]["T"]) - (Ds["Ap"]["T"] - D0["Ap"]["T"])) / D0["Ap"]["T"]
        dU = np.abs(((Ds[f]["U"] - D0[f]["U"]) - (Ds["Ap"]["U"] - D0["Ap"]["U"]))) / np.maximum(D0["Ap"]["U"], 1e-3)
        dth = ((Ds[f]["th"] - D0[f]["th"]) - (Ds["Ap"]["th"] - D0["Ap"]["th"])) / np.abs(D0["Ap"]["th"])
        row = []
        for nm, m in REG:
            core, wall = m & (J < 100) & (J >= 1), m & (J >= 100) & (J <= nj - 2)
            xm = np.array([np.any(m[i]) for i in range(ni)])
            th = dth[xm & np.isfinite(dth)]
            row.append(f"{nm}: |ΔT|/T 中心 {np.median(np.abs(dT[core])):.1e}/{np.percentile(np.abs(dT[core]), 90):.1e} 壁側 {np.median(np.abs(dT[wall])):.1e}/{np.percentile(np.abs(dT[wall]), 90):.1e}、|ΔU|/U 中心 {np.median(dU[core]):.1e}、Δθ_r 平均 {np.mean(th) * 100 if th.size else float('nan'):+.2f} %")
        print(f"  {f}:\n    " + "\n    ".join(row))
