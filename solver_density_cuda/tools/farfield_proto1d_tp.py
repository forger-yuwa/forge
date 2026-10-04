#!/usr/bin/env python3
"""遠方境界 (farfield) の TP ホストゲート V0h (plan boundary-node-farfield-characteristic §6)。

`farfield_proto1d.py` (CPG) の採用方式を thermally-perfect の 2 擬似種 (SERN の EXH / AIR、NASA-9、`FrozenGas`) に広げる。
node 型 1D 格子 (両端の境界節点は Δx/2 の半 CV)、内部面は TP SLAU + MUSCL (ρ, u, P, Y、van Leer)、境界半割面は 1 次で、
外側状態 (内部エントロピーの TRRS [frozen γ_i]、自由流の密度・組成・エントロピー、真空置換、原始変数の超音速ブレンド) + TP HLLC (Davis 波速)。
時間は SSP-RK3 (時間精度、CFL 0.25)。保存量 (ρ, ρu, ρE, ρY_EXH)、E = e_sens(T) + u²/2 (forge と同じ href 298.15 K)。

試験 (合否は plan §6 V0h):
  V0p  保存形の混合による圧力誤差 (境界と無関係、記録のみ)
  T1   一様な外気 (T 220 K、P 2851 Pa、AIR) の音響反射 M 0 / 0.3
  T1b  高温の排気 (T 600 K、Y_EXH 0.13) の中を出ていく音響、自由流は外気 (γ・R・ρc が違う)。パルスなしの対照も
  T3   接触波 (温度・組成の塊が出ていく): 短領域 − 長領域 (境界が加える誤差) と、長領域自身の誤差を分けて記録
  T2   超音速流出の極限 (外気を掃引しても F = F(U_i))、T4 連続性、T6 正値性・真空置換・HLL 退避の計数

  cd <repo>/design && python3 ../solver_density_cuda/tools/farfield_proto1d_tp.py    (数分、CPU のみ)
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "design"))
from forge_design.gas.frozen import FrozenGas  # noqa: E402

EXH_MASS = {"N2": 0.732638, "H2O": 0.241109, "H2": 0.00103232, "AR": 0.0125424, "OH": 0.00420498, "O2": 0.00521324,
            "NO": 0.00253027, "H": 5.15746e-05, "O": 0.000242324, "CO2": 0.000378318, "CO": 5.73291e-05}
GE = FrozenGas(EXH_MASS, "EXH"); GA = FrozenGas.air()
TT = np.arange(100.0, 4000.0, 0.5)
HE, HA = GE.h_sens(TT), GA.h_sens(TT)
CPE, CPA = GE.cp_mass(TT), GA.cp_mass(TT)
RE, RA = GE.R, GA.R
CNT = {"hll": 0, "vac": 0}


def mix(Y):
    return Y * RE + (1 - Y) * RA


def h_of(T, Y):
    return Y * np.interp(T, TT, HE) + (1 - Y) * np.interp(T, TT, HA)


def cp_of(T, Y):
    return Y * np.interp(T, TT, CPE) + (1 - Y) * np.interp(T, TT, CPA)


def e_of(T, Y):
    return h_of(T, Y) - mix(Y) * T


def T_from_e(e, Y, T0=None):
    T = np.full_like(e, 400.0) if T0 is None else T0.copy()
    for _ in range(30):
        f = e_of(T, Y) - e; cv = cp_of(T, Y) - mix(Y)
        T = np.clip(T - f / cv, 101.0, 3999.0)
        if np.max(np.abs(f / cv)) < 1e-9 * np.max(T):
            break
    return T


def gam(T, Y):
    cp = cp_of(T, Y); return cp / (cp - mix(Y))


def snd(T, Y):
    return np.sqrt(gam(T, Y) * mix(Y) * T)


def prim_cons(r, u, p, Y):
    T = p / (r * mix(Y)); E = e_of(T, Y) + 0.5 * u * u
    return np.stack([r, r * u, r * E, r * Y])


def cons_prim(U, T0=None):
    r = U[0]; u = U[1] / r; Y = np.clip(U[3] / r, 0.0, 1.0)
    e = U[2] / r - 0.5 * u * u
    T = T_from_e(e, Y, T0); p = r * mix(Y) * T
    return r, u, p, Y, T


def phys_flux(r, u, p, Y):
    T = p / (r * mix(Y)); H = h_of(T, Y) + 0.5 * u * u
    m = r * u
    return np.stack([m, m * u + p, m * H, m * Y])


def slau_tp(rL, uL, pL, YL, rR, uR, pR, YR):
    TL = pL / (rL * mix(YL)); TR = pR / (rR * mix(YR))
    cL = snd(TL, YL); cR = snd(TR, YR); ch = 0.5 * (cL + cR)
    Mp, Mm = uL / ch, uR / ch
    bp = np.where(np.abs(Mp) >= 1, 0.5 * (1 + np.sign(Mp)), 0.25 * (Mp + 1) ** 2 * (2 - Mp))
    bm = np.where(np.abs(Mm) >= 1, 0.5 * (1 - np.sign(Mm)), 0.25 * (Mm - 1) ** 2 * (2 + Mm))
    g = -np.maximum(np.minimum(Mp, 0), -1) * np.minimum(np.maximum(Mm, 0), 1)
    vh = (rL * np.abs(uL) + rR * np.abs(uR)) / (rL + rR)
    vhp = (1 - g) * vh + g * np.abs(uL); vhm = (1 - g) * vh + g * np.abs(uR)
    Mh = np.minimum(1.0, np.sqrt(0.5 * (uL * uL + uR * uR)) / ch); chi = (1 - Mh) ** 2
    md = 0.5 * (rL * (uL + vhp) + rR * (uR - vhm) - chi / ch * (pR - pL))
    pt = 0.5 * (pL + pR) + 0.5 * (bp - bm) * (pL - pR) + (1 - chi) * (bp + bm - 1) * 0.5 * (pL + pR)
    HL = h_of(TL, YL) + 0.5 * uL * uL; HR = h_of(TR, YR) + 0.5 * uR * uR
    up = md >= 0
    return np.stack([md, md * np.where(up, uL, uR) + pt, md * np.where(up, HL, HR), md * np.where(up, YL, YR)])


def hllc_tp(rL, uL, pL, YL, rR, uR, pR, YR):
    TL = pL / (rL * mix(YL)); TR = pR / (rR * mix(YR))
    cL = snd(TL, YL); cR = snd(TR, YR)
    SL = np.minimum(uL - cL, uR - cR); SR = np.maximum(uL + cL, uR + cR)
    Sm = (pR - pL + rL * uL * (SL - uL) - rR * uR * (SR - uR)) / (rL * (SL - uL) - rR * (SR - uR))
    rsL = rL * (SL - uL) / (SL - Sm); rsR = rR * (SR - uR) / (SR - Sm)
    bad = ~((SL <= Sm) & (Sm <= SR) & (rsL > 0) & (rsR > 0) & np.isfinite(Sm))
    CNT["hll"] += int(np.sum(bad))
    UL = prim_cons(rL, uL, pL, YL); UR = prim_cons(rR, uR, pR, YR)
    FL = phys_flux(rL, uL, pL, YL); FR = phys_flux(rR, uR, pR, YR)

    def star(r, u, p, Y, S, U):
        f = r * (S - u) / (S - Sm)
        return np.stack([f, f * Sm, f * (U[2] / r + (Sm - u) * (Sm + p / (r * (S - u)))), f * Y])
    Fc = np.where(SL >= 0, FL, np.where(Sm >= 0, FL + SL * (star(rL, uL, pL, YL, SL, UL) - UL),
                  np.where(SR > 0, FR + SR * (star(rR, uR, pR, YR, SR, UR) - UR), FR)))
    Fh = (SR * FL - SL * FR + SL * SR * (UR - UL)) / (SR - SL)
    Fh = np.where(SL >= 0, FL, np.where(SR <= 0, FR, Fh))
    return np.where(bad, Fh, Fc)


def smooth01(x):
    t = np.clip(x, 0.0, 1.0); return t * t * (3 - 2 * t)


def outer_state(ri, ui, pi, Yi, rinf, uinf, pinf, Yinf, band=0.1):
    """採用方式の外側状態 (plan §4.2) の TP 版。frozen γ: 内部側 γ_i、自由流側 γ_∞。"""
    Ti = pi / (ri * mix(Yi)); gi = gam(Ti, Yi); ci = np.sqrt(gi * pi / ri)
    Tinf = pinf / (rinf * mix(Yinf)); ginf = gam(Tinf, Yinf); ainf = np.sqrt(ginf * pinf / rinf)
    z = (gi - 1) / (2 * gi)
    cpo = ci * (pinf / pi) ** z
    B = ci + cpo - 0.5 * (gi - 1) * (uinf - ui)
    vac = B <= 1e-3 * (ci + cpo)
    Bs = np.where(vac, ci + cpo, B)
    ps = (Bs / (ci / pi ** z + cpo / pinf ** z)) ** (1 / z)
    us = ui + 2 * ci / (gi - 1) * (1 - (ps / pi) ** z)
    rb = rinf * (ps / pinf) ** (1 / ginf); ub, pb = us, ps; yb = Yinf * np.ones_like(ri)
    wq = smooth01(((-1 + band) - uinf / ainf) / band)
    rb = (1 - wq) * rb + wq * rinf; ub = (1 - wq) * ub + wq * uinf; pb = (1 - wq) * pb + wq * pinf; yb = (1 - wq) * yb + wq * Yinf
    wi = smooth01((ui / ci - (1 - band)) / band)
    rb = (1 - wi) * rb + wi * ri; ub = (1 - wi) * ub + wi * ui; pb = (1 - wi) * pb + wi * pi; yb = (1 - wi) * yb + wi * Yi
    bad = vac | ~(np.isfinite(rb) & np.isfinite(pb) & (rb > 0) & (pb > 0))
    CNT["vac"] += int(np.sum(bad))
    return (np.where(bad, ri, rb), np.where(bad, ui, ub), np.where(bad, pi, pb), np.where(bad, Yi, yb))


def bflux(kind, ri, ui, pi, Yi, inf):
    rinf, uinf, pinf, Yinf = (np.full_like(ri, v) for v in inf)
    if kind == "ghost_hllc":
        return hllc_tp(ri, ui, pi, Yi, rinf, uinf, pinf, Yinf)
    if kind == "ghost_slau":
        return slau_tp(ri, ui, pi, Yi, rinf, uinf, pinf, Yinf)
    if kind == "farfield":
        rb, ub, pb, yb = outer_state(ri, ui, pi, Yi, rinf, uinf, pinf, Yinf)
        return hllc_tp(ri, ui, pi, Yi, rb, ub, pb, yb)
    raise ValueError(kind)


def vanleer(a, b):
    return np.where(a * b > 0, 2 * a * b / (a + b + 1e-300), 0.0)


def run(kind, x0, x1, dx, init, t_end, probe, inf, cfl=0.25):
    n = int(round((x1 - x0) / dx)) + 1
    x = np.linspace(x0, x1, n); vol = np.full(n, dx); vol[0] = vol[-1] = 0.5 * dx
    U = prim_cons(*init(x)); ip = int(np.argmin(np.abs(x - probe)))
    Tc = [None]
    inf_l = (inf[0], -inf[1], inf[2], inf[3])

    def rhs(U):
        r, u, p, Y, T = cons_prim(U, Tc[0]); Tc[0] = T
        W = np.stack([r, u, p, Y]); d = np.diff(W, axis=1)
        s = np.zeros_like(W); s[:, 1:-1] = vanleer(d[:, :-1], d[:, 1:])
        WL = W[:, :-1] + 0.5 * s[:, :-1]; WR = W[:, 1:] - 0.5 * s[:, 1:]
        F = slau_tp(*WL, *WR)
        R = np.zeros_like(U); R[:, :-1] -= F; R[:, 1:] += F
        R[:, -1:] -= bflux(kind, r[-1:], u[-1:], p[-1:], Y[-1:], inf)
        Fl = bflux(kind, r[:1], -u[:1], p[:1], Y[:1], inf_l) * np.array([[1.0], [-1.0], [1.0], [1.0]])
        R[:, :1] -= Fl
        return R / vol
    t = 0.0; ts = [0.0]; ps = [cons_prim(U)[2][ip]]
    while t < t_end - 1e-15:
        r, u, p, Y, T = cons_prim(U, Tc[0])
        dt = min(cfl * dx / np.max(np.abs(u) + snd(T, Y)), t_end - t)
        U1 = U + dt * rhs(U); U2 = 0.75 * U + 0.25 * (U1 + dt * rhs(U1)); U = U / 3 + 2 / 3 * (U2 + dt * rhs(U2))
        t += dt; ts.append(t); ps.append(cons_prim(U, Tc[0])[2][ip])
    return np.array(ts), np.array(ps)


P_INF, T_INF = 2851.0, 220.0
R_INF = P_INF / (RA * T_INF); A_INF = float(snd(np.array([T_INF]), np.array([0.0]))[0])


def acoustic(kind, dx, M, Tin=T_INF, Yin=0.0, amp=1e-3, pulse=True):
    u0 = M * A_INF
    rin = P_INF / (mix(Yin) * Tin); cin = float(snd(np.array([Tin]), np.array([Yin]))[0])
    sig = 0.04 / 2.3548; inf = (R_INF, u0, P_INF, 0.0)

    def init(x):
        pp = (amp * P_INF * np.exp(-0.5 * ((x - 0.5) / sig) ** 2)) if pulse else 0 * x
        return rin + pp / cin ** 2, u0 + pp / (rin * cin), P_INF + pp, np.full_like(x, Yin)
    t_end = 0.5 / (u0 + cin) + 0.4 / max(cin - u0, 1.0) + 0.1 / cin
    ts, pa = run(kind, -1.0, 1.0, dx, init, t_end, 0.8, inf)
    if not pulse:
        return np.max(np.abs(pa - P_INF)) / (amp * P_INF)
    tl, pl = run(kind, -1.0, 3.0, dx, init, t_end, 0.8, inf)
    return np.max(np.abs(pa - np.interp(ts, tl, pl))) / np.max(np.abs(pl - P_INF))


def contact(kind, dx):
    u0 = 0.3 * A_INF; inf = (R_INF, u0, P_INF, 0.0); sig = 0.1 / 2.3548

    def init(x):
        w = np.exp(-0.5 * ((x - 0.4) / sig) ** 2)
        T = T_INF + (600.0 - T_INF) * w; Y = 0.13 * w
        return P_INF / (mix(Y) * T), np.full_like(x, u0), np.full_like(x, P_INF), Y
    t_end = 1.2 / u0 * 0.5
    ts, pa = run(kind, -1.0, 1.0, dx, init, t_end, 0.9, inf)
    tl, pl = run(kind, -1.0, 3.0, dx, init, t_end, 0.9, inf)
    return np.max(np.abs(pa - np.interp(ts, tl, pl))) / P_INF, np.max(np.abs(pl - P_INF)) / P_INF


def main():
    print(f"自由流 (外気): T {T_INF} K, P {P_INF} Pa, ρ {R_INF:.5g}, a {A_INF:.4g} m/s, γ {float(gam(np.array([T_INF]), np.array([0.0]))[0]):.4f}")
    # V0p
    rH = P_INF / (mix(0.13) * 600.0); UH = prim_cons(np.array([rH]), np.array([0.0]), np.array([P_INF]), np.array([0.13]))
    UA = prim_cons(np.array([R_INF]), np.array([0.0]), np.array([P_INF]), np.array([0.0]))
    r, u, p, Y, T = cons_prim(0.9 * UH + 0.1 * UA)
    print(f"\n== V0p 保存形の混合 (0.9 高温排気 + 0.1 外気): T {T[0]:.2f} K, Y {Y[0]:.5f}, P 誤差 {(p[0] / P_INF - 1) * 100:+.4f} % (境界と無関係、記録のみ)")
    print("\n== T2 超音速流出の極限: 内部 (外気, u = 1.1 a) 固定、自由流の速度 −2〜2 a・圧力 0.2〜5 倍を掃引 → max|F/F(U_i) − 1|")
    ua = np.linspace(-2, 2, 41) * A_INF; pa = np.array([0.2, 1.0, 5.0]) * P_INF
    UU, PP = np.meshgrid(ua, pa, indexing="ij"); n = UU.size
    ri = np.full(n, R_INF); ui = np.full(n, 1.1 * A_INF); pi = np.full(n, P_INF); Yi = np.zeros(n)
    Fi = phys_flux(ri[:1], ui[:1], pi[:1], Yi[:1])[:, 0]
    for k in ("ghost_slau", "ghost_hllc", "farfield"):
        err = []
        for j in range(n):
            F = bflux(k, ri[j:j + 1], ui[j:j + 1], pi[j:j + 1], Yi[j:j + 1], (R_INF, UU.ravel()[j], PP.ravel()[j], 0.0))[:, 0]
            err.append(np.max(np.abs(F[:3] / Fi[:3] - 1)))
        print(f"  {k:11s} {max(err):.2e}")
    print("\n== T4 連続性 (内部の速度を −1.2〜1.2 a で 1e-4 a 刻み、自由流は平行、内部は高温排気): 隣接差の最大 / 規模")
    uu = np.arange(-1.2, 1.2, 1e-4) * 486.0; o = np.ones_like(uu)
    rH_ = P_INF / (mix(0.13) * 600.0)
    for k in ("ghost_hllc", "farfield"):
        CNT["hll"] = CNT["vac"] = 0
        F = bflux(k, rH_ * o, uu, P_INF * o, 0.13 * o, (R_INF, 0.0, P_INF, 0.0))
        jump = np.max(np.abs(np.diff(F, axis=1)), axis=1) / (np.max(np.abs(F), axis=1) + 1e-30)
        print(f"  {k:11s} " + "  ".join(f"{j:.1e}" for j in jump) + f"   HLL 退避 {CNT['hll']}  真空置換 {CNT['vac']}")
    print("\n== T6 正値性・置換の計数 (内部: T 200–2000 K・Y 0/0.13/1・P 0.1–10 倍・u −3〜3 a、自由流は外気で u −3〜10 a)")
    CNT["hll"] = CNT["vac"] = 0
    Ts = np.array([200.0, 600.0, 2000.0]); Ys = np.array([0.0, 0.13, 1.0]); Ps = np.array([0.1, 1.0, 10.0]) * P_INF
    us = np.linspace(-3, 3, 13) * A_INF; uo = np.linspace(-3, 10, 27) * A_INF
    G = np.meshgrid(Ts, Ys, Ps, us, uo, indexing="ij"); T_, Y_, P_, U_, UO = (g.ravel() for g in G)
    r_ = P_ / (mix(Y_) * T_)
    nonfin = 0
    for j in range(0, T_.size, 2000):
        s = slice(j, j + 2000)
        for uo_val in np.unique(UO[s]):
            m = UO[s] == uo_val
            F = bflux("farfield", r_[s][m], U_[s][m], P_[s][m], Y_[s][m], (R_INF, uo_val, P_INF, 0.0))
            nonfin += int(np.sum(~np.isfinite(F)))
    print(f"  {T_.size} 点: 非有限 {nonfin}, 真空置換 {CNT['vac']}, HLL 退避 {CNT['hll']}")
    print("\n== T1 一様な外気の音響反射 (短 − 長 / 入射)")
    for M in (0.0, 0.3):
        for dx in (5e-3, 2.5e-3):
            CNT["hll"] = CNT["vac"] = 0
            row = "  ".join(f"{k}:{acoustic(k, dx, M) * 100:6.2f}%" for k in ("ghost_slau", "ghost_hllc", "farfield"))
            print(f"  M {M:.1f} Δx {dx:.1e}  {row}   (退避 {CNT['hll']}, 真空 {CNT['vac']})")
    print("\n== T1b 高温の排気 (T 600 K, Y 0.13) を出ていく音響 (u 0.5 a∞)、自由流は外気")
    for dx in (5e-3, 2.5e-3):
        CNT["hll"] = CNT["vac"] = 0
        row = "  ".join(f"{k}:{acoustic(k, dx, 0.5, 600.0, 0.13) * 100:6.2f}%" for k in ("ghost_hllc", "farfield"))
        nop = acoustic("farfield", dx, 0.5, 600.0, 0.13, pulse=False)
        print(f"  Δx {dx:.1e}  {row}   パルスなし対照 (farfield): |P−P∞| / パルス振幅 = {nop:.2e}   (退避 {CNT['hll']}, 真空 {CNT['vac']})")
    print("\n== T3 接触波 (T 600 K・Y 0.13 の塊が出ていく、M 0.3): (a) 境界が加える誤差 |短−長|/P、(b) 長領域自身の |P−P∞|/P")
    for dx in (5e-3, 2.5e-3):
        row = []
        for k in ("ghost_hllc", "farfield"):
            a_, b_ = contact(k, dx); row.append(f"{k}: (a) {a_:.2e} (b) {b_:.2e}")
        print(f"  Δx {dx:.1e}  " + "  ".join(row))


if __name__ == "__main__":
    main()
