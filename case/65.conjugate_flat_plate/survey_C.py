#!/usr/bin/env python3
r"""C (共役平板) の条件選定 (plan §5.1 #3b)。Blasius 流速・源項なしの参照解 (case/64 の conjugate_ref.py、平面) で、
板の厚さ比 b/L と伝導比 k_s/k_f ごとに、固体の厚さ方向温度差 ΔT_s と板断面の軸方向熱量 Q_ax を出す。forge の結果は見ない。"""
import sys
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "64.conjugate_pipe_wall"))
import conjugate_ref as cr

def blasius():
    """f''' + f f''/2 = 0, f(0)=f'(0)=0, f'(∞)=1。f''(0) を射撃法で決める。"""
    def rhs(e, y): return [y[1], y[2], -0.5 * y[0] * y[2]]
    lo, hi = 0.2, 0.5
    for _ in range(60):
        s = 0.5 * (lo + hi); sol = solve_ivp(rhs, [0, 12], [0, 0, s], rtol=1e-11, atol=1e-12)
        lo, hi = (s, hi) if sol.y[1, -1] < 1 else (lo, s)
    sol = solve_ivp(rhs, [0, 12], [0, 0, s], rtol=1e-11, atol=1e-12, dense_output=True)
    return s, sol

fpp0, BL = blasius()
k_f, cp, rho, mu = 5.700284e-2, 1004.5, 1.176829, 4.085818e-5     # case/63 と同じ定数物性 (Pr 0.72)
U = 0.1 * 347.1887                                                 # M 0.1
nu = mu / rho

def run(bL, ks_ratio, ReL=1e4, h_back=1e8, Th=1.0, nx_scale=1.0, ny_f=48, ny_s=16):
    L = ReL * nu / U; b = bL * L; k_s = ks_ratio * k_f
    dx0 = 0.004 * L / nx_scale
    xs = np.concatenate([cr.graded(-0.5 * L, 0.0, int(30 * nx_scale), dx0, toward_a=False),
                         cr.graded(0.0, L, int(120 * nx_scale), dx0)[1:], cr.graded(L, 1.5 * L, int(20 * nx_scale), dx0 * 5)[1:]])
    delta = 5 * L / np.sqrt(ReL); H = 6 * delta
    yf = cr.graded(0.0, H, ny_f, H / ny_f / 8)
    ys = np.concatenate([np.linspace(-b, 0, ny_s + 1)[:-1], yf])
    js = ny_s                                                       # 界面の行
    mat = np.full((len(xs) - 1, len(ys) - 1), cr.FLUID)
    for i in range(len(xs) - 1):
        xm = 0.5 * (xs[i] + xs[i + 1])
        mat[i, :js] = cr.SOLID if 0 <= xm <= L else cr.VOID
    nx, ny = len(xs), len(ys)
    u = np.zeros((nx, ny)); v = np.zeros((nx, ny))
    for i, x in enumerate(xs):
        for j in range(js, ny):
            y = ys[j]
            if x <= 0 or x > L:
                xe = max(x, 1e-12) if x > L else None
            if x <= 0:
                u[i, j] = U; continue
            eta = y * np.sqrt(U / (nu * x)); eta = min(eta, 12)
            f, fp, _ = BL.sol(eta)
            u[i, j] = U * fp; v[i, j] = 0.5 * np.sqrt(nu * U / x) * (eta * fp - f)
    rob = lambda xm, y: (h_back, Th) if (0 <= xm <= L and abs(y + b) < 1e-15) else None
    P = cr.Problem(xs, ys, mat, k_f, k_s, cp, False, 0.0, ro=np.full((nx, ny), rho), u=u, v=v, robin=rob)
    T = P.solve()
    plate = (xs >= 0) & (xs <= L)
    dTs = T[plate, 0] - T[plate, js]                                # 裏面 − 界面
    # 板の中の節点だけで x 微分を取る (板の外の未使用節点を拾うと端で桁違いになる — 初版の誤り)
    ip = np.where(plate)[0]
    Tx = np.gradient(T[ip], xs[ip], axis=0)
    trap = getattr(np, "trapezoid", None) or np.trapz
    Qax_p = np.array([-trap(k_s * Tx[k, :js + 1], ys[:js + 1]) for k in range(len(ip))])
    Qtot = P.robin_heat()
    win = (xs >= 0.2 * L) & (xs <= 0.9 * L)
    return dict(L=L, b=b, k_s=k_s, Ti_win=(T[win, js].min(), T[win, js].max()), dTs_win=(dTs[(xs[plate] >= 0.2 * L) & (xs[plate] <= 0.9 * L)].min(), dTs.max()),
                Qax_max=np.abs(Qax_p).max() / Qtot, Qtot=Qtot)

if __name__ == "__main__":
    print(f"Blasius f''(0) = {fpp0:.6f} (文献 0.332057)")
    print("  b/L  ks/kf   L[mm]  窓内の界面温度 θ_i=(T_i−T∞)/(T_h−T∞)   窓内の厚さ方向差 min / 全体 max   |Q_ax|max/Q_tot")
    for bL in (0.05, 0.2):
        for ksr in (1, 10, 100):
            d = run(bL, ksr)
            print(f" {bL:4.2f}  {ksr:5d}  {d['L']*1e3:6.3f}   {d['Ti_win'][0]:.3f} … {d['Ti_win'][1]:.3f}                    {d['dTs_win'][0]:.3f} / {d['dTs_win'][1]:.3f}            {d['Qax_max']:.4f}")
