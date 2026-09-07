#!/usr/bin/env python3
"""⑤ SERN 平面 MOC + key point 逆設計の単体テスト (plan §6 S1)。pytest 非依存。"""
import sys, time
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from forge_design.geometry.moc_sern import PlanarMOC, SernKernelSpec, wall_forces, _solve_nu_minus_mu  # noqa: E402
from forge_design.geometry.moc_kernel import pm_nu, pm_mach, mu_vec  # noqa: E402
from forge_design.geometry.rao_planar import theta_lip, p_over_pt, lip_residual  # noqa: E402

FAIL = 0
def check(name, cond, detail=""):
    global FAIL
    print(("ok   " if cond else "FAIL ") + name + ("" if not detail else f"  [{detail}]"))
    if not cond:
        FAIL += 1

g = 1.4
def area_ratio(M):
    return (1.0 / M) * ((2.0 / (g + 1.0)) * (1.0 + 0.5 * (g - 1.0) * M * M)) ** ((g + 1.0) / (2.0 * (g - 1.0)))

# --- 0. ν−μ 逆解 ------------------------------------------------------------
M = np.array([1.2, 2.0, 4.0, 8.0])
Mr = _solve_nu_minus_mu(pm_nu(M, g) - mu_vec(M), g)
check("ν−μ 二分法の逆解 1e-10", np.max(np.abs(Mr - M)) < 1e-9, f"{np.max(np.abs(Mr-M)):.1e}")

# --- 1. 一様流保持 -------------------------------------------------------------
t0 = time.time()
k = PlanarMOC(SernKernelSpec(M_in=2.0, theta_r0=0.0, theta_c0=0.0, L_cowl=3.0, x_max=3.0, nj=101, dx=5e-3)).march()
check("一様流 (θ_r0=θ_c0=0) が機械精度で保持", np.max(np.abs(k.M - 2.0)) < 1e-12 and np.max(np.abs(k.TH)) < 1e-14,
      f"dM {np.max(np.abs(k.M-2.0)):.1e}, {time.time()-t0:.2f}s")

# --- 2. ランプ角部扇の解析比較 (カウル直線, 反射到達前の領域) ----------------------
M_in, th_r = 2.0, np.deg2rad(10.0)
t0 = time.time()
k = PlanarMOC(SernKernelSpec(M_in=M_in, theta_r0=th_r, theta_c0=0.0, L_cowl=6.0, x_max=6.0, nj=301, dx=2e-3)).march()
tm = time.time() - t0
nu_in = float(pm_nu(M_in, g)); mu_in = float(mu_vec(M_in))
errs = []
for (x, yfrac) in [(0.5, 0.95), (1.0, 0.9), (1.5, 0.97), (2.0, 0.85)]:
    y = yfrac * k.y_up(x) + (1 - yfrac) * 0.0
    phi = np.arctan2(y - 1.0, x)
    th, nu, Mk = k.state_at(x, y)
    if phi <= -mu_in:
        Ma, tha = M_in, 0.0
    else:
        Mfan = _solve_nu_minus_mu(np.array([phi + nu_in]), g)[0]
        Mr_ = pm_mach(nu_in + th_r, g)
        if phi >= th_r - float(mu_vec(Mr_)):
            Ma, tha = Mr_, th_r
        else:
            Ma, tha = Mfan, float(pm_nu(Mfan, g)) - nu_in
    errs.append((abs(Mk - Ma), abs(th - tha)))
errs = np.array(errs)
check("ランプ扇: 単純波域で M 誤差 < 1e-6, θ 誤差 < 1e-6", errs[:, 0].max() < 1e-6 and errs[:, 1].max() < 1e-6,
      f"dM {errs[:,0].max():.1e} dθ {errs[:,1].max():.1e}, march {tm:.2f}s")
# 扇のカウル反射後: 壁 (カウル) 上で θ=0、下壁 M は反射後に ν = ν_in + 2θ_r (単純波が壁で反射)
Mw = k.M[-1, 0]
check("カウル反射後の下壁 M = M(ν_in + 2θ_r)", abs(Mw - pm_mach(nu_in + 2 * th_r, g)) < 2e-4,
      f"{Mw:.5f} vs {pm_mach(nu_in + 2*th_r, g):.5f}")

# --- 3. 対称 MLN 極限 (鋭角スロート、両壁 θ_max=(ν_e−ν_in)/2) -------------------------
M_in, M_e = 1.5, 3.0
th_max = 0.5 * (float(pm_nu(M_e, g)) - float(pm_nu(M_in, g)))
t0 = time.time()
k = PlanarMOC(SernKernelSpec(M_in=M_in, theta_r0=th_max, theta_c0=th_max, L_cowl=8.0, x_max=8.0, nj=401, dx=2e-3)).march()
d = k.design_ramp(M_c=M_e, f=0.5, ds=5e-4)
tm = time.time() - t0
xc, yc, Mc, thc = d.key_point
check("MLN: key point が対称線上 (|y_c−0.5| < 2e-3)", abs(yc - 0.5) < 2e-3, f"y_c={yc:.5f}")
check("MLN: θ_c ≈ 0 (< 0.05°)", abs(np.rad2deg(thc)) < 0.05, f"θ_c={np.rad2deg(thc):.4f}°")
check("MLN: 壁終端角 θ_e ≈ 0 (< 0.05°)", abs(np.rad2deg(d.info['theta_e'])) < 0.05, f"{np.rad2deg(d.info['theta_e']):.4f}°")
H_e = 2.0 * (d.lip_e[1] - 0.5)
AR = area_ratio(M_e) / area_ratio(M_in)
check("MLN: 出口高さ = 等エントロピー面積比 (0.2%)", abs(H_e / AR - 1) < 2e-3, f"H_e={H_e:.5f} vs {AR:.5f}")
check("MLN: c–e を横切る流量 = f (0.2%)", abs(d.mass_fraction_check / 0.5 - 1) < 2e-3, f"{d.mass_fraction_check:.5f}")
check("MLN: a–c 上で K⁻ = θ+ν 一定 (< 1e-4)", d.info["K_minus_spread"] < 1e-4, f"{d.info['K_minus_spread']:.1e}")
# 推力: p_a = p_e のとき総推力 = 出口運動量 γ M_e² p_e/p_in H_e (一様出口)
pe = float(p_over_pt(M_e, g) / p_over_pt(M_in, g))
# 対称なのでランプ側だけ (上半分): 壁積分は上壁のみ → 総推力 = 入口流れ推力(上半分) + 上壁推力
fr = wall_forces(d, M_in, g, pa_over_pin=pe)
F_exit_half = g * M_e ** 2 * pe * (H_e / 2)
# wall_forces は入口全高 (H=1) の流れ推力 + ランプ壁 + カウル壁 (カウルは直線 θ_c0 の下壁) を足す。
# 対称 MLN ではカウル壁 (直線) の推力は上壁の kernel 直線部と同じ。ここでは上半分だけで検算する:
inlet_half = 0.5 * (g * M_in ** 2 + (1 - pe))
F_half = inlet_half + fr["T_ramp"]
check("MLN: 上半分の総推力 = 出口運動量流束 (0.3%)", abs(F_half / F_exit_half - 1) < 3e-3,
      f"{F_half:.5f} vs {F_exit_half:.5f}, {tm:.1f}s")
print(f"     MLN: L_ramp={d.L_ramp:.4f}, x_a={d.foot_a[0]:.4f}, rays={d.info['n_rays']}")

# --- 4. カウル付き SERN (自由境界) の逆設計 + 格子収束 -------------------------------
def run_sern(nj, dx, ds=1e-3):
    k = PlanarMOC(SernKernelSpec(M_in=2.5, theta_r0=np.deg2rad(15.0), theta_c0=np.deg2rad(5.0),
                                 L_cowl=1.0, x_max=9.0, nj=nj, dx=dx, p_ext_over_p_in=0.05)).march()
    d = k.design_ramp(M_c=3.8, f=0.6, ds=ds)
    return k, d, wall_forces(d, 2.5, g)
t0 = time.time()
k1, d1, f1 = run_sern(301, 2e-3)
t1 = time.time() - t0
k2, d2, f2 = run_sern(601, 1e-3, 5e-4)
check("SERN: 逆設計が完走 (警告なし)", len(d1.info["warnings"]) == 0, f"{d1.info['warnings']}")
check("SERN: a は直線ランプ上 (0 < x_a < x_c)", 0.0 < d1.foot_a[0] < d1.key_point[0], f"x_a={d1.foot_a[0]:.4f} x_c={d1.key_point[0]:.4f}")
check("SERN: c–e 流量 = f (0.3%)", abs(d1.mass_fraction_check / 0.6 - 1) < 3e-3, f"{d1.mass_fraction_check:.5f}")
check("SERN: 壁終端角 = θ_c", abs(d1.info["theta_e"] - d1.key_point[3]) < 1e-6)
check("SERN: 格子倍で C_T 変化 < 0.3%", abs(f1["C_T"] / f2["C_T"] - 1) < 3e-3, f"{f1['C_T']:.5f} vs {f2['C_T']:.5f}")
check("SERN: 格子倍で L_ramp 変化 < 0.5%", abs(d1.L_ramp / d2.L_ramp - 1) < 5e-3, f"{d1.L_ramp:.4f} vs {d2.L_ramp:.4f}")
print(f"     SERN: p_TE/p_in={k1.p_te_over_p_in:.4f}, θ_c={np.rad2deg(d1.key_point[3]):.2f}°, L_ramp={d1.L_ramp:.3f}, "
      f"C_T={f1['C_T']:.4f} (wall {f1['C_T_wall']:.4f}) C_L={f1['C_L']:.4f} C_M={f1['C_M']:.4f}, {t1:.1f}s")
# 縁条件 (Rao) の残差診断
pc = float(p_over_pt(d1.key_point[2], g) / p_over_pt(2.5, g))
print(f"     SERN: 縁条件残差 θ_c − θ_lip = {np.rad2deg(lip_residual(d1.key_point[2], d1.key_point[3], pc / f1['pa_over_pin'], g)):.3f}°")

# --- 5. 縁条件 ------------------------------------------------------------------
check("θ_lip(pe=pa) = 0", theta_lip(3.0, 1.0, g) == 0.0)
check("θ_lip は pe/pa に単調増加", theta_lip(3.0, 1.5, g) < theta_lip(3.0, 2.0, g) < theta_lip(3.0, 3.0, g))

# --- 6. 最適性の数値掃引 (診断: 等長の候補間で C_T 最大点の縁条件残差) ---------------
print("--- 最適性掃引 (診断、assert なし): kernel = test 4, p_a = p_ext ---")
rows = []
for f in (0.45, 0.55, 0.65, 0.75, 0.85):
    for Mc in (3.0, 3.3, 3.6, 3.9, 4.2, 4.5):
        try:
            d = k1.design_ramp(M_c=Mc, f=f)
        except ValueError as e:
            print(f"     skip f={f} M_c={Mc}: {str(e)[:70]}")
            continue
        fr = wall_forces(d, 2.5, g)
        pc = float(p_over_pt(d.key_point[2], g) / p_over_pt(2.5, g))
        rows.append((f, Mc, d.L_ramp, fr["C_T"], np.rad2deg(d.key_point[3]),
                     np.rad2deg(lip_residual(d.key_point[2], d.key_point[3], pc / fr["pa_over_pin"], g))))
rows = np.array(rows)
print("     f     M_c    L_ramp   C_T     θ_c[°]  θ_c−θ_lip[°]")
for r in rows[np.argsort(rows[:, 2])]:
    print("     " + "  ".join(f"{v:7.3f}" for v in r))
np.savetxt(Path(__file__).with_name("sern_optimality_sweep.csv"), rows, delimiter=",",
           header="f,M_c,L_ramp,C_T,theta_c_deg,lip_residual_deg", comments="")

print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAIL'}")
sys.exit(1 if FAIL else 0)
