"""plan condensation-two-phase-default #4sr: 改訂尺度での最大比の面を分解し、上界の前提 (逆数・ctg が正規化数) を全評価面で確認する。
  python3 decompose_max.py tp_faces.h5 [A.bin faces.bin.idx.npy]
"""
import sys, numpy as np, h5py
f = h5py.File(sys.argv[1]); g = lambda k: np.array(f[k][...])
ev = g("face/skip") == 0
rho0 = g("on_in/rho0").astype(np.float32); rho1 = g("on_in/rho1").astype(np.float32)
ct = g("on_in/ct").astype(np.float32); geo = g("on_in/geo").astype(np.float32)
rg0 = g("on_in/rg0").astype(np.float64); rg1 = g("on_in/rg1").astype(np.float64)
tiny = np.finfo(np.float32).tiny
ir0 = (np.float32(1) / rho0); ir1 = (np.float32(1) / rho1); ctg_f = ct * geo
print("precondition (all evaluated faces): ir0/ir1 normal:", bool(np.all(np.abs(ir0[ev]) >= tiny)), bool(np.all(np.abs(ir1[ev]) >= tiny)),
      "| ctg_f normal or zero:", bool(np.all((np.abs(ctg_f[ev]) >= tiny) | (ctg_f[ev] == 0))), "| ctg range", np.abs(ctg_f[ev]).min(), np.abs(ctg_f[ev]).max())
Jf = g("on_f/Jl").astype(np.float64); Jd = g("on_d/Jl")
u = 2.0 ** -24; eta = 2.0 ** -150; EPS32 = 2.0 ** -23
ctg = ct.astype(np.float64) * geo.astype(np.float64)
Al = np.abs(ctg) * (np.abs(rg0 / rho0) + np.abs(rg1 / rho1))
eabs = (3.0 * (1 + u) ** 3 * np.abs(ctg) + 1.0) * eta
rel = 8 * EPS32 * Al
err = np.abs(Jf - Jd)
with np.errstate(divide="ignore", invalid="ignore"):
    r = np.where(ev & ((rel + eabs) > 0), err / (rel + eabs), 0.0)
i = int(np.nanargmax(r))
p0 = rg0[i] * ir0[i]; p1 = rg1[i] * ir1[i]
print(f"max ratio {r[i]:.4f} at face {i}: Jf {Jf[i]:.6e} Jd {Jd[i]:.6e} err {err[i]:.3e} (= {err[i]/eta:.2f} eta)")
print(f"  rel term 8eps*A_l {rel[i]:.3e} | abs term E_abs {eabs[i]:.3e} | ctg {ctg[i]:.3e} | rg0 {rg0[i]:.3e} rg1 {rg1[i]:.3e} | p0 {p0:.3e} p1 {p1:.3e} | d {p1-p0:.3e}")
print("faces with ratio > 0.9:", int(np.sum(r > 0.9)), "| > 1:", int(np.sum(r > 1)))
