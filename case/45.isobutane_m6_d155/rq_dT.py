"""§6.40 の記録 (判定ではない): 各 run の step 2,000 の温度を、FP64 の C1・C2 の平均と K で比べる。対象は x/r_t < −2 の内点。"""
import os, sys
import numpy as np, h5py
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cold_xcheck as XC  # noqa: E402
os.chdir(os.path.dirname(os.path.abspath(__file__)))
with h5py.File("run_0487_tr_fp64/nozzle.h5", "r") as h: xy0 = np.asarray(h["MESH/COORD"][:], float)
ni, nj, Sc = XC.mesh_info(xy0.size // 3); X = (xy0.reshape(-1, 3)[:, 0] / Sc).reshape(ni, nj)
I, J = np.meshgrid(np.arange(ni), np.arange(nj), indexing="ij")
DOM = ((X < -2) & (I >= 1) & (J >= 1) & (J <= nj - 2)).ravel()
HOT = ((np.abs(X + 7.1) < 0.15) & (J >= 100) & (J <= 106)).ravel()
def T(r, s): return np.asarray(h5py.File(f"{r}/res_{s}.h5", "r")["VALUE/T"][:], float)
ALL = ("run_0595_rq_c1", "run_0596_rq_c2", "run_0597_rq_c3", "run_0598_rq_e1", "run_0599_rq_e2", "run_0620_rq_e3", "run_0621_rq_f1", "run_0622_rq_f2", "run_0623_rq_f3")
runs = [r for r in ALL if os.path.exists(f"{r}/RUN_RC")]
ref = np.mean([T(r, 2000) for r in ("run_0595_rq_c1", "run_0596_rq_c2")], axis=0)
print(f"基準 = FP64 の C1・C2 の平均 (step 2,000)。対象 = 入口〜収縮部の内点 {int(DOM.sum())} 節点、差の最大の所 (x/r_t≈−7, j 100〜106) {int(HOT.sum())} 節点 (温度 {np.mean(ref[HOT]):.1f} K)")
print("run                 RMS ΔT [K]   最大 |ΔT| [K]   |ΔT|>1K の割合   差の最大の所の平均 ΔT [K]")
for r in runs + ["run_0487_tr_fp64", "run_0488_tr_f32", "run_0489_tr_f32b"]:
    d = T(r, 2000) - ref
    print(f"{r:18s} {np.sqrt(np.mean(d[DOM] ** 2)):11.4f} {np.max(np.abs(d[DOM])):14.3f} {np.mean(np.abs(d[DOM]) > 1) * 100:14.2f} % {np.mean(d[HOT]):+20.3f}")
