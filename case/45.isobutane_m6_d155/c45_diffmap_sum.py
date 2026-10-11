"""c45_diffmap.npz (B で作る大きい中間) から、x ごとの要約を作る (記録。判定ではない)。→ _band_ab/cold_pair/c45_diffmap_sum.npz
- θ_A0(x) (A の step 0 の θ_r、r_t 単位)、各 step の float − FP64 の θ_r の差 (絶対値、r_t 単位) と相対値
- q_w(x_w) の float − FP64 の差の相対値
- 場の ρ・|U|・T・P の float − FP64 の差の相対値の、列 (j) の区分ごとの |·| の最大: 軸 (j 0〜8)・中心 (j 9〜99)・壁側 (j 100〜119)・壁 (j 120)
- A の step 0 のマッハ数: 軸 (j 0)・中心 (j 60)・壁の 1 つ内側 (j 119)"""
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
z = np.load(HERE / "_band_ab/cold_pair/c45_diffmap.npz")
steps = z["steps"]; out = {"steps": steps, "x": z["x"], "x_w": z["x_w"], "xg0": z["xg"][:, 0], "rw": z["rg"][:, -1]}
nj = z["xg"].shape[1]
J = {"axis": slice(0, 9), "core": slice(9, 100), "nearwall": slice(100, nj - 1), "wall": slice(nj - 1, nj)}
if "M0" in z:
    for nm, j in (("axis", 0), ("mid", 60), ("wall1", nj - 2)): out[f"M0_{nm}"] = z["M0"][:, j]
for a in ("B", "Bp"):
    for s in steps:
        out[f"dth_rel_{a}_{s}"] = z[f"d_theta_{a}_{s}"]
        out[f"dqw_rel_{a}_{s}"] = z[f"d_qw_{a}_{s}"]
        for q in ("ro", "U", "T", "P"):
            k = f"f_{q}_{a}_{s}"
            if k not in z: continue
            f = np.abs(z[k])
            for nm, sl in J.items(): out[f"fmax_{q}_{nm}_{a}_{s}"] = np.nanmax(f[:, sl], axis=1)
import sys; sys.path.insert(0, str(HERE)); import cold_xcheck as XC  # noqa: E402
yb_x, yb = XC.common_yb(); xy, ro, ux, uy, T, k = XC.load_forge(HERE / "run_0487_tr_fp64", 0); ni, nj2, S = XC.mesh_info(len(ro))
o = XC.reduce_fields(xy, ro, ux, uy, T, k, False, ni, nj2, S, yb_x, yb)
out["theta_A0"] = o["theta_r"]; out["qw_A0"] = o["q_w"]; out["yb_i"] = np.interp(o["x"], yb_x, yb)
np.savez_compressed(HERE / "_band_ab/cold_pair/c45_diffmap_sum.npz", **out); print("→ _band_ab/cold_pair/c45_diffmap_sum.npz", len(out), "配列")
