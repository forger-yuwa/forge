"""c45_diffmap.npz から、差の大きい節点の場所と大きさを見る (記録)。step 1000・10000・30000、float B (run_0488)。
- |U| と T の差 (相対) の上位の節点の x/r_t・j・その点の |U|・マッハ数
- 収縮部 (x/r_t ∈ [−8, −2)) と入口 (x < −8) で、相対差が 1e-3 を超える節点の割合 (列の区分ごと)"""
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
z = np.load(HERE / "_band_ab/cold_pair/c45_diffmap.npz")
xg = z["xg"]; nj = xg.shape[1]; M0 = z["M0"] if "M0" in z else None
for s in (1000, 10000, 30000):
    for q in ("U", "T"):
        f = np.abs(z[f"f_{q}_B_{s}"]); ff = np.where(np.isfinite(f), f, -1)
        idx = np.argsort(ff, axis=None)[::-1][:6]
        top = [(round(float(xg.flat[k]), 2), int(k % nj), f"{ff.flat[k]:.2e}", f"M {M0.flat[k]:.4f}" if M0 is not None else "") for k in idx]
        print(f"step {s:5d} {q}: 上位 (x/r_t, j, 相対差, M) {top}")
    for nm, m in (("入口 x<-8", xg < -8), ("収縮 -8..-2", (xg >= -8) & (xg < -2)), ("スロート -2..2", (xg >= -2) & (xg < 2)), ("超音速 >2", xg >= 2)):
        fu = np.abs(z[f"f_U_B_{s}"]); ft = np.abs(z[f"f_T_B_{s}"])
        print(f"   {nm:14s}: |ΔU|/U > 1e-3 の節点の割合 {np.mean(fu[m] > 1e-3) * 100:5.1f} %、|ΔT|/T > 1e-4 の割合 {np.mean(ft[m] > 1e-4) * 100:5.1f} %、ΔT/T の中央値 {np.nanmedian(ft[m]):.1e}")
