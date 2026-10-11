"""float と FP64 の差が x のどこから育つかを見る (記録。判定ではない)。float plan §6.27 の軌跡の run を使う (新しい run は回さない)。
A = run_0487_tr_fp64 (FP64)、B = run_0488_tr_f32・B' = run_0489_tr_f32b (float)。共通の Q32 の起点からの各 step で:
  - θ_r(x) と壁の q_w(x)・τ_w(x) の「自分の step 0 からの変化」の float − FP64 (A の step 0 の値で割る)
  - 場の ρ・|U|・T の節点ごとの「step 0 からの変化」の float − FP64 (A の step 0 の値で割る)。壁法線の列 j ごとに x の分布として保存
  - A の step 0 のマッハ数 (軸 j=0、帯の外縁付近、壁の 1 つ内側) の x 分布
使い方 (B で): python3 c45_diffmap.py → _band_ab/cold_pair/c45_diffmap.npz (作図は手元で c45_diffmap_plot.py)。"""
import sys
from pathlib import Path
import numpy as np, h5py
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cold_xcheck as XC  # noqa: E402
RUNS = {"A": "run_0487_tr_fp64", "B": "run_0488_tr_f32", "Bp": "run_0489_tr_f32b"}
STEPS = (0, 1000, 2000, 5000, 10000, 20000, 30000)
yb_x, yb = XC.common_yb()
def load(run, st):
    xy, ro, ux, uy, T, k = XC.load_forge(HERE / run, st)
    with h5py.File(HERE / run / f"res_{st}.h5", "r") as h:
        c = np.asarray(h["VALUE/sonic"][:], float) if "VALUE/sonic" in h else None
        P = np.asarray(h["VALUE/P"][:], float) if "VALUE/P" in h else None
    ni, nj, S = XC.mesh_info(len(ro))
    o = XC.reduce_fields(xy, ro, ux, uy, T, k, False, ni, nj, S, yb_x, yb)
    U = np.hypot(ux, uy)
    return dict(x=o["x"], x_w=o["x_w"], theta=o["theta_r"], qw=o["q_w"], tauw=o["tau_w"],
                ro=ro.reshape(ni, nj), U=U.reshape(ni, nj), T=T.reshape(ni, nj), P=None if P is None else P.reshape(ni, nj),
                M=None if c is None else (U / c).reshape(ni, nj), xg=(xy[:, 0] / S).reshape(ni, nj), rg=(xy[:, 1] / S).reshape(ni, nj))
D = {a: {s: load(r, s) for s in STEPS} for a, r in RUNS.items()}
out = {"steps": np.array(STEPS), "x": D["A"][0]["x"], "x_w": D["A"][0]["x_w"], "xg": D["A"][0]["xg"], "rg": D["A"][0]["rg"]}
if D["A"][0]["M"] is not None: out["M0"] = D["A"][0]["M"]
for a in ("B", "Bp"):
    for s in STEPS:
        for q in ("theta", "qw", "tauw"):
            ref = D["A"][0][q]
            out[f"d_{q}_{a}_{s}"] = ((D[a][s][q] - D[a][0][q]) - (D["A"][s][q] - D["A"][0][q])) / np.where(np.abs(ref) > 0, np.abs(ref), np.nan)
        for q in ("ro", "U", "T", "P"):
            if D["A"][0][q] is None: continue
            ref = np.abs(D["A"][0][q])
            out[f"f_{q}_{a}_{s}"] = ((D[a][s][q] - D[a][0][q]) - (D["A"][s][q] - D["A"][0][q])) / np.where(ref > 0, ref, np.nan)
# 簡単な要約: 各 step で |Δθ_r| が 0.1 % を超える最初の x (B と B')
for a in ("B", "Bp"):
    for s in STEPS[1:]:
        d = np.abs(out[f"d_theta_{a}_{s}"]); ok = np.isfinite(d) & (d > 1e-3)
        print(f"{a} step {s:6d}: |Δθ_r| > 0.1 % の最初の x/r_t = {out['x'][ok][0] if ok.any() else None}  最大 {np.nanmax(d) * 100:.3f} % (x/r_t {out['x'][np.nanargmax(d)]:.1f})")
np.savez_compressed(HERE / "_band_ab/cold_pair/c45_diffmap.npz", **out); print("→ _band_ab/cold_pair/c45_diffmap.npz")
