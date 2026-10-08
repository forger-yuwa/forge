"""SU2 の restart (ASCII CSV) を forge の CPG の run の初期値から作る (plan tooling-nozzle-isothermal-wall-chain §5.1 #26)。
forge と SU2 を同じ初期値から始める。節点の順番は msh と同じ (座標で照合して確かめる)。
SU2 の Energy は k を含む (CNSVariable::SetPrimVar で E − ½u² − k が静的エネルギー) ので、forge (k を含まない) の P を保って
Energy = roe_forge + ρk とする。Turb_Kin_Energy・Omega は原始変数 (k・ω)。
usage: python3 cold_su2_ic.py <ic_npz> <nozzle.su2 (座標の照合用)> <出力 csv>
"""
import sys
import numpy as np

ic = np.load(sys.argv[1])
# 座標の照合は SU2 の格子ファイル (nozzle.su2) の節点で行う
with open(sys.argv[2]) as fh:
    for ln, line in enumerate(fh):
        if line.startswith("NPOIN="):
            npt = int(line.split("=")[1].split()[0]); break
ref = np.loadtxt(sys.argv[2], skiprows=ln + 1, max_rows=npt, usecols=(0, 1))
xy = ic["coord"]
if xy.shape != ref.shape or float(np.max(np.abs(xy - ref))) > 1e-14:
    raise SystemExit(f"節点の座標が SU2 と一致しない (最大差 {float(np.max(np.abs(xy - ref))) if xy.shape == ref.shape else 'shape'})")
ro = ic["ro"]; k = ic["roK"] / ro; om = ic["roOmega"] / ro
E = ic["roe"] + ic["roK"]
n = len(ro)
A = np.column_stack([np.arange(n), xy[:, 0], xy[:, 1], ro, ic["roUx"], ic["roUy"], E, k, om])
np.savetxt(sys.argv[3], A, delimiter=", ", fmt=["%d"] + ["%.15e"] * 8, comments="",
           header='"PointID","x","y","Density","Momentum_x","Momentum_y","Energy","Turb_Kin_Energy","Omega"')
print(f"{sys.argv[3]}: {n} 点、座標の最大差 {float(np.max(np.abs(xy - ref))):.1e}、ρ {ro.min():.3g}〜{ro.max():.3g}、k {k.min():.3g}〜{k.max():.3g}、ω {om.min():.3g}〜{om.max():.3g}")
