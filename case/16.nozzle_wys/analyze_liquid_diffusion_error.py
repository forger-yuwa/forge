"""凝縮域の拡散の近似誤差を既存 run の場から見積もる (流束の大きさの比のみ。温度への影響の予測ではない:
二相 EOS では液の乱流輸送の潜熱項が相殺し、温度に効くのは R_w T 程度 — codex 2026-09-27 諮問 notes/reviews/2026-09-27-condensation-diffusion-error-diagnose.md) (plan condensation-two-phase-transport §5.1, 2026-09-27)。

現行実装 (コード読み取り):
  - 総水分 roY_w は Fick 拡散 J_w = -rho (D_w + D_t) grad(Y_w) (D_w: kinetic 混合平均, D_t = mu_t/(rho Sc_t)) で、
    エネルギーに h_v J_w を加える (speciesTransport_d.cu species_diffusion_d)。
  - 液 rog (凝縮モーメント) は拡散しない (passive_diffusion_d はトレーサだけが呼ぶ)。
  => 拡散で動くのは実質「蒸気」(g はそのまま) で、エネルギー勘定は蒸気として整合している。

目標モデル (plan §4.2): 蒸気は -rho D_w grad(Y_v) - rho D_t grad(Y_v)、液は -rho D_t grad(g) (分子拡散なし) で h_l を運ぶ。
差 (現行 - 目標):
  (i)  余計な蒸気の分子流束      E1 = rho D_w |grad g|            (蒸気が液の勾配で押される)
  (ii) 乱流で液を蒸気として運ぶ  E2 = rho D_t |grad g|  (質量は同じ, 相が違う) -> エネルギー差 L * E2
比較の基準:
  - 蒸気の分子流束               F1 = rho D_w |grad Y_v|
  - 熱流束 (層流 + 乱流)         Q  = (lambda + cp mu_t / Pr_t) |grad T|
出力: g > g_min のノードでの比 E1/F1, L*E2/Q の分位点と、液の乱流混合の相当温度 L*|dg|/cp のスケール。

usage: python3 analyze_liquid_diffusion_error.py RUN_DIR [--gmin 1e-4] [--sct 0.7] [--prt 0.9] [--cp 1040]
"""
import argparse, glob, math, os
import h5py, numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("run"); ap.add_argument("--gmin", type=float, default=1e-4)
ap.add_argument("--sct", type=float, default=0.7); ap.add_argument("--prt", type=float, default=0.9)
ap.add_argument("--cp", type=float, default=1040.0)
a = ap.parse_args()

res = sorted(glob.glob(os.path.join(a.run, "res_[0-9]*.h5")), key=lambda s: int(os.path.basename(s)[4:-3]))[-1]
mesh_h5 = [f for f in glob.glob(os.path.join(a.run, "*.h5")) if not os.path.basename(f).startswith("res_")][0]
m = h5py.File(mesh_h5)
xy = m["MESH/COORD"][:].reshape(-1, 3)[:, :2]
# VIZMESH/CONNE は XDMF mixed 形式 [型, 節点...] (4 = 三角形 3 節点, 5 = 四角形 4 節点)
vc = m["VIZMESH/CONNE"][:]; tris = []; i = 0
while i < len(vc):
    ty = int(vc[i])
    if ty == 5: q4 = vc[i+1:i+5]; tris += [[q4[0], q4[1], q4[2]], [q4[0], q4[2], q4[3]]]; i += 5
    elif ty == 4: tris.append(list(vc[i+1:i+4])); i += 4
    else: raise ValueError(f"未対応の要素型 {ty}")
tri = np.asarray(tris, dtype=np.int64)

v = h5py.File(res)["VALUE"]
ro, T, P, g = v["ro"][:], v["T"][:], v["P"][:], v["g_0"][:]
Yw = v["roY1"][:] / ro
Yv = np.maximum(Yw - g, 0.0)
mut, lam = v["vis_turb"][:], v["thermCond"][:]

def grad(f):
    # 三角形ごとの線形勾配をノードへ面積重み平均 (P1 の勾配回復)
    p0, p1, p2 = xy[tri[:, 0]], xy[tri[:, 1]], xy[tri[:, 2]]
    d1, d2 = p1 - p0, p2 - p0
    det = d1[:, 0]*d2[:, 1] - d1[:, 1]*d2[:, 0]
    ok = np.abs(det) > 0
    f0, f1, f2 = f[tri[:, 0]], f[tri[:, 1]], f[tri[:, 2]]
    df1, df2 = f1 - f0, f2 - f0
    gx = np.where(ok, (df1*d2[:, 1] - df2*d1[:, 1]) / np.where(ok, det, 1), 0.0)
    gy = np.where(ok, (d1[:, 0]*df2 - d2[:, 0]*df1) / np.where(ok, det, 1), 0.0)
    A = 0.5*np.abs(det)
    sx = np.zeros(len(xy)); sy = np.zeros(len(xy)); sa = np.zeros(len(xy))
    for k in range(3):
        np.add.at(sx, tri[:, k], A*gx); np.add.at(sy, tri[:, k], A*gy); np.add.at(sa, tri[:, k], A)
    sa = np.where(sa > 0, sa, np.nan)
    return np.hypot(sx/sa, sy/sa)

gT, gg, gYv = grad(T), grad(g), grad(Yv)

# H2O–N2 二元拡散 (Chapman–Enskog, Neufeld Ω11) を D_w の代わりに使う (希薄 H2O なので混合平均 ≈ 二元)
def om11(Ts): return 1.06036/Ts**0.15610 + 0.19300/np.exp(0.47635*Ts) + 1.03587/np.exp(1.52996*Ts) + 1.76474/np.exp(3.89411*Ts)
Mi, Mj, s, e = 18.0153, 28.0134, 0.5*(2.605+3.621), math.sqrt(572.4*97.53)
Dw = 1.8583e-3*np.sqrt(T**3*(1/Mi+1/Mj)) / ((P/101325.0)*s*s*om11(T/e)) * 1e-4   # [m^2/s]
Dt = mut/(ro*a.sct)
L = 2.5e6

E1 = ro*Dw*gg; F1 = ro*Dw*gYv
E2 = ro*Dt*gg; Q = (lam + a.cp*mut/a.prt)*gT
mask = (g > a.gmin) & np.isfinite(gg) & np.isfinite(gT)
def q(x): return np.percentile(x[mask], [50, 95, 99, 100])
r1 = E1/np.maximum(F1, 1e-300); r2 = L*E2/np.maximum(Q, 1e-300)
turb = mask & (mut > 10*v["vis_lam"][:])
print(f"run: {a.run}  res: {os.path.basename(res)}  nodes g>{a.gmin}: {mask.sum()}  (うち mu_t/mu>10: {turb.sum()})")
print("  (i)  余計な蒸気分子流束 / 蒸気分子流束  E1/F1  p50/p95/p99/max = " + " / ".join(f"{x:.3g}" for x in q(r1)))
print("  (ii) 液を蒸気として乱流輸送したエネルギー差 / 熱流束  L*E2/Q  p50/p95/p99/max = " + " / ".join(f"{x:.3g}" for x in q(r2)))
if turb.any():
    print("       (mu_t/mu>10 の乱流域だけ) p50/p95/max = " + " / ".join(f"{x:.3g}" for x in np.percentile(r2[turb], [50, 95, 100])))
i = np.argmax(np.where(mask, r2, -1))
print(f"  (ii) の最大点: x={xy[i,0]*1e3:.2f} mm y={xy[i,1]*1e3:.3f} mm wall_dist={v['wall_dist'][i]*1e3:.4f} mm g={g[i]:.4g} T={T[i]:.1f} mu_t/mu={mut[i]/v['vis_lam'][i]:.3g}")
print(f"  D_w (H2O-N2) p50 on mask = {np.percentile(Dw[mask],50):.3g} m^2/s, D_t p50 = {np.percentile(Dt[mask],50):.3g} m^2/s")
print(f"  液の乱流混合の相当温度スケール L*g_max/cp = {L*g.max()/a.cp:.1f} K (g_max {g.max():.4g})")
