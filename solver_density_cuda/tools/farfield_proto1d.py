#!/usr/bin/env python3
"""遠方境界 (farfield) の境界流束の候補を 1 次元で比較するホスト試作 (plan boundary-node-farfield-characteristic §5.1 #1g)。

node (median-dual) 方式を 1 次元で模す: 節点 x_0..x_N、内部 CV 幅 Δx、両端の境界 CV は Δx/2。
内部面は SLAU (forge の convectiveFlux_slau_d と同じ式) + MUSCL (van Leer、原始変数)、境界半割面は 1 次 (境界節点の値) で候補の流束。
時間積分は SSP-RK3、CFL 0.25 (音速 + |u| 基準)。CPG γ = 1.4、無次元 (ρ∞ = 1, P∞ = 1/γ, c∞ = 1)。

候補 (右端 = 遠方境界。左端は十分遠い同じ遠方境界):
  slau_ghost : F_SLAU(U_i, U_∞)
  hllc_ghost : F_HLLC(U_i, U_∞)
  roe_ghost  : F_Roe(U_i, U_∞) (Harten の entropy fix)
  char_hllc  : 局所線形化の特性振幅 (w± = P ± ρ_i c_i u_n) で U_b を組み立て、F_HLLC(U_i, U_b)
  char_slau  : 同上で F_SLAU(U_i, U_b)

試験:
  T1 音響反射: 右向きだけの純音波 (ガウス、振幅 1e-3 P∞) を M 0 / 0.3 で。短領域 [-1, 1] と長領域 [-1, 3] の観測点 (x = 0.8) の圧力差の最大 / 入射振幅。Δx 3 水準。
  T2 超音速流出の極限: U_i = (1, 2, 1/γ)、U_∞ = (1, 3, 1/γ) で面流束が内部の物理流束に一致するか。
  T3 接触波: M 0.3、P 一様・u 一様で、温度 2.7 倍 (密度 1/2.7) のガウス塊が右端から出ていく。観測点の |P − P∞|/P∞ の最大。
  T4 連続性: U_∞ の法線速度 0 (平行) で U_i の速度を −1.2〜1.2 まで 1e-4 刻みに掃引し、面流束の隣接差の最大 / 流束規模。
  T5 対照: M 0 で右端を slip (鏡像) にしたときの反射率 (≈ 1 になること = 試験が反射を検出できる)。

  python3 solver_density_cuda/tools/farfield_proto1d.py   (数十秒、CPU のみ)
"""
import numpy as np

G = 1.4
P0 = 1.0 / G


def prim2cons(r, u, p):
    return np.stack([r, r * u, p / (G - 1) + 0.5 * r * u * u])


def cons2prim(U):
    r = U[0]; u = U[1] / r; p = (G - 1) * (U[2] - 0.5 * r * u * u)
    return r, u, p


def phys_flux(r, u, p):
    E = p / (G - 1) + 0.5 * r * u * u
    return np.stack([r * u, r * u * u + p, (E + p) * u])


def slau(rL, uL, pL, rR, uR, pR):
    """forge convectiveFlux_slau_d.inc.cuh と同じ式 (1D、接線速度 0、SLAU1)。"""
    cL = np.sqrt(G * pL / rL); cR = np.sqrt(G * pR / rR); ch = 0.5 * (cL + cR)
    Mp = uL / ch; Mm = uR / ch
    bp = np.where(np.abs(Mp) >= 1, 0.5 * (1 + np.sign(Mp)), 0.25 * (Mp + 1) ** 2 * (2 - Mp))
    bm = np.where(np.abs(Mm) >= 1, 0.5 * (1 - np.sign(Mm)), 0.25 * (Mm - 1) ** 2 * (2 + Mm))
    g = -np.maximum(np.minimum(Mp, 0), -1) * np.minimum(np.maximum(Mm, 0), 1)
    vh = (rL * np.abs(uL) + rR * np.abs(uR)) / (rL + rR)
    vhp = (1 - g) * vh + g * np.abs(uL); vhm = (1 - g) * vh + g * np.abs(uR)
    Mh = np.minimum(1.0, np.sqrt(0.5 * (uL * uL + uR * uR)) / ch); chi = (1 - Mh) ** 2
    md = 0.5 * (rL * (uL + vhp) + rR * (uR - vhm) - chi / ch * (pR - pL))
    pt = 0.5 * (pL + pR) + 0.5 * (bp - bm) * (pL - pR) + (1 - chi) * (bp + bm - 1) * 0.5 * (pL + pR)
    HL = G / (G - 1) * pL / rL + 0.5 * uL * uL; HR = G / (G - 1) * pR / rR + 0.5 * uR * uR
    up = md >= 0
    return np.stack([md, md * np.where(up, uL, uR) + pt, md * np.where(up, HL, HR)])


def roe_avg(rL, uL, pL, rR, uR, pR):
    sL, sR = np.sqrt(rL), np.sqrt(rR)
    HL = G / (G - 1) * pL / rL + 0.5 * uL * uL; HR = G / (G - 1) * pR / rR + 0.5 * uR * uR
    u = (sL * uL + sR * uR) / (sL + sR); H = (sL * HL + sR * HR) / (sL + sR)
    c = np.sqrt(np.maximum((G - 1) * (H - 0.5 * u * u), 1e-12))
    return u, c, H, sL * sR


def hllc(rL, uL, pL, rR, uR, pR):
    cL = np.sqrt(G * pL / rL); cR = np.sqrt(G * pR / rR)
    ut, ct, _, _ = roe_avg(rL, uL, pL, rR, uR, pR)
    SL = np.minimum(uL - cL, ut - ct); SR = np.maximum(uR + cR, ut + ct)
    Sm = (pR - pL + rL * uL * (SL - uL) - rR * uR * (SR - uR)) / (rL * (SL - uL) - rR * (SR - uR))
    FL = phys_flux(rL, uL, pL); FR = phys_flux(rR, uR, pR)
    UL = prim2cons(rL, uL, pL); UR = prim2cons(rR, uR, pR)

    def star(r, u, p, S, U):
        f = r * (S - u) / (S - Sm)
        E = U[2]
        return np.stack([f, f * Sm, f * (E / r + (Sm - u) * (Sm + p / (r * (S - u))))])
    UsL = star(rL, uL, pL, SL, UL); UsR = star(rR, uR, pR, SR, UR)
    F = np.where(SL >= 0, FL, np.where(Sm >= 0, FL + SL * (UsL - UL), np.where(SR > 0, FR + SR * (UsR - UR), FR)))
    return F


def roe(rL, uL, pL, rR, uR, pR):
    u, c, H, rt = roe_avg(rL, uL, pL, rR, uR, pR)
    dr, du, dp = rR - rL, uR - uL, pR - pL
    a1 = (dp - rt * c * du) / (2 * c * c); a2 = dr - dp / (c * c); a3 = (dp + rt * c * du) / (2 * c * c)
    l1, l2, l3 = np.abs(u - c), np.abs(u), np.abs(u + c)
    d = 0.1 * c
    fix = lambda l: np.where(l < d, (l * l + d * d) / (2 * d), l)
    l1, l3 = fix(l1), fix(l3)
    K1 = np.stack([np.ones_like(u), u - c, H - u * c]); K2 = np.stack([np.ones_like(u), u, 0.5 * u * u]); K3 = np.stack([np.ones_like(u), u + c, H + u * c])
    return 0.5 * (phys_flux(rL, uL, pL) + phys_flux(rR, uR, pR)) - 0.5 * (l1 * a1 * K1 + l2 * a2 * K2 + l3 * a3 * K3)


def char_state(ri, ui, pi, rinf, uinf, pinf):
    """局所線形化の特性振幅 (Whitfield–Janus): w+ = P + Z u は内部、w- = P - Z u は外、エントロピーは流向の側。"""
    ci = np.sqrt(G * pi / ri); Z = ri * ci
    sup_out = ui / ci >= 1; sup_in = ui / ci <= -1
    pb = 0.5 * (pi + pinf) + 0.5 * Z * (ui - uinf)
    ub = 0.5 * (ui + uinf) + (pi - pinf) / (2 * Z)
    rs = np.where(ub > 0, ri, rinf); ps = np.where(ub > 0, pi, pinf)
    rb = rs + (pb - ps) / (ci * ci)
    rb = np.where(sup_out, ri, np.where(sup_in, rinf, rb)); ub = np.where(sup_out, ui, np.where(sup_in, uinf, ub)); pb = np.where(sup_out, pi, np.where(sup_in, pinf, pb))
    return rb, ub, pb


FALLBACK = {"n": 0}


def hll_flux(rL, uL, pL, rR, uR, pR, SL, SR):
    FL = phys_flux(rL, uL, pL); FR = phys_flux(rR, uR, pR)
    UL = prim2cons(rL, uL, pL); UR = prim2cons(rR, uR, pR)
    Fm = (SR * FL - SL * FR + SL * SR * (UR - UL)) / (SR - SL)
    return np.where(SL >= 0, FL, np.where(SR <= 0, FR, Fm))


def hllc_davis(rL, uL, pL, rR, uR, pR):
    """HLLC、波速は各側の u ± c から直接 (Davis、EOS に依らない)。S* ∉ [SL, SR] か星密度 ≤ 0 なら HLL に退避して数える。"""
    cL = np.sqrt(G * pL / rL); cR = np.sqrt(G * pR / rR)
    SL = np.minimum(uL - cL, uR - cR); SR = np.maximum(uL + cL, uR + cR)
    Sm = (pR - pL + rL * uL * (SL - uL) - rR * uR * (SR - uR)) / (rL * (SL - uL) - rR * (SR - uR))
    rsL = rL * (SL - uL) / (SL - Sm); rsR = rR * (SR - uR) / (SR - Sm)
    bad = ~((SL <= Sm) & (Sm <= SR) & (rsL > 0) & (rsR > 0) & np.isfinite(Sm))
    FALLBACK["n"] += int(np.sum(bad))
    FL = phys_flux(rL, uL, pL); FR = phys_flux(rR, uR, pR)
    UL = prim2cons(rL, uL, pL); UR = prim2cons(rR, uR, pR)

    def star(r, u, p, S, U):
        f = r * (S - u) / (S - Sm)
        return np.stack([f, f * Sm, f * (U[2] / r + (Sm - u) * (Sm + p / (r * (S - u))))])
    Fc = np.where(SL >= 0, FL, np.where(Sm >= 0, FL + SL * (star(rL, uL, pL, SL, UL) - UL),
                  np.where(SR > 0, FR + SR * (star(rR, uR, pR, SR, UR) - UR), FR)))
    return np.where(bad, hll_flux(rL, uL, pL, rR, uR, pR, SL, SR), Fc)


def charghost_state(ri, ui, pi, rinf, uinf, pinf):
    """外側の状態: 圧力・法線速度は境界節点の Z = ρ_i c_i で線形化した特性量 (外向き w+ は内部、内向き w- は外気)、
    密度 (エントロピー)・接線速度・組成は常に外気側。超音速流入面 (Q_n ≤ -a∞) は外気そのもの。"""
    ci = np.sqrt(G * pi / ri); Z = ri * ci; ainf = np.sqrt(G * pinf / rinf)
    pb = 0.5 * (pi + pinf) + 0.5 * Z * (ui - uinf)
    ub = 0.5 * (ui + uinf) + (pi - pinf) / (2 * Z)
    rb = rinf * np.maximum(pb / pinf, 1e-12) ** (1.0 / G)       # 外気のエントロピーで
    supin = uinf <= -ainf
    return np.where(supin, rinf, rb), np.where(supin, uinf, ub), np.where(supin, pinf, pb)


def bflux(kind, ri, ui, pi, rinf, uinf, pinf):
    if kind == "slau_ghost": return slau(ri, ui, pi, rinf, uinf, pinf)
    if kind == "hllc_ghost": return hllc(ri, ui, pi, rinf, uinf, pinf)
    if kind == "roe_ghost": return roe(ri, ui, pi, rinf, uinf, pinf)
    if kind in ("char_hllc", "char_slau"):
        rb, ub, pb = char_state(ri, ui, pi, rinf, uinf, pinf)
        return (hllc if kind == "char_hllc" else slau)(ri, ui, pi, rb, ub, pb)
    if kind == "hllcd_ghost": return hllc_davis(ri, ui, pi, rinf, uinf, pinf)
    if kind == "charghost_hllcd":
        rb, ub, pb = charghost_state(ri, ui, pi, rinf, uinf, pinf)
        return hllc_davis(ri, ui, pi, rb, ub, pb)
    if kind == "slip":
        return slau(ri, ui, pi, ri, -ui, pi)
    raise ValueError(kind)


def vanleer(a, b):
    return np.where(a * b > 0, 2 * a * b / (a + b + 1e-300), 0.0)


def run(kind, x0, x1, dx, init, t_end, probe, uinf, cfl=0.25, left_kind=None):
    n = int(round((x1 - x0) / dx)) + 1
    x = np.linspace(x0, x1, n); vol = np.full(n, dx); vol[0] = vol[-1] = 0.5 * dx
    r, u, p = init(x)
    U = prim2cons(r, u, p)
    rinf, pinf = 1.0, P0
    ip = int(np.argmin(np.abs(x - probe)))
    lk = left_kind or ("hllc_ghost" if kind == "slip" else kind)

    def rhs(U):
        r, u, p = cons2prim(U)
        W = np.stack([r, u, p])
        d = np.diff(W, axis=1)
        s = np.zeros_like(W); s[:, 1:-1] = vanleer(d[:, :-1], d[:, 1:])
        WL = W[:, :-1] + 0.5 * s[:, :-1]; WR = W[:, 1:] - 0.5 * s[:, 1:]
        F = slau(WL[0], WL[1], WL[2], WR[0], WR[1], WR[2])
        R = np.zeros_like(U)
        R[:, :-1] -= F; R[:, 1:] += F
        Fr = bflux(kind, r[-1:], u[-1:], p[-1:], np.array([rinf]), np.array([uinf]), np.array([pinf]))
        R[:, -1:] -= Fr
        Fl = bflux(lk, r[:1], -u[:1], p[:1], np.array([rinf]), np.array([-uinf]), np.array([pinf]))  # 左端は外向き法線 −x
        Fl = Fl * np.array([[1.0], [-1.0], [1.0]])
        R[:, :1] -= Fl
        return R / vol
    t = 0.0; ts, ps_ = [0.0], [p[ip]]
    while t < t_end - 1e-12:
        r, u, p = cons2prim(U)
        dt = min(cfl * dx / np.max(np.abs(u) + np.sqrt(G * p / r)), t_end - t)
        U1 = U + dt * rhs(U); U2 = 0.75 * U + 0.25 * (U1 + dt * rhs(U1)); U = U / 3 + 2 / 3 * (U2 + dt * rhs(U2))
        t += dt
        ts.append(t); ps_.append(cons2prim(U)[2][ip])
    return np.array(ts), np.array(ps_)


def acoustic(kind, M, dx, amp=1e-3, left_kind=None):
    uinf = M
    sig = 0.04 / 2.3548

    def init(x):
        pp = amp * P0 * np.exp(-0.5 * ((x - 0.5) / sig) ** 2)
        return 1.0 + pp, uinf + pp, P0 + pp          # 右向き純音波: ρ' = p'/c², u' = p'/(ρc) (ρ=c=1)
    t_end = 0.5 / (1 + M) + 0.4 / max(1 - M, 1e-3) + 0.1
    ts, pa = run(kind, -1.0, 1.0, dx, init, t_end, 0.8, uinf, left_kind=left_kind)
    tl, pl = run("hllc_ghost" if kind == "slip" else kind, -1.0, 3.0, dx, init, t_end, 0.8, uinf, left_kind=left_kind)
    pl_i = np.interp(ts, tl, pl)
    A = np.max(np.abs(pl - P0))
    return np.max(np.abs(pa - pl_i)) / A


def acoustic_hot(kind, dx, amp=1e-3):
    """codex plan-7 M1 の設定: 内部 ρ = 220/600 (温度 600/220 倍)、P = P∞、u = 0.5 (内外)、内部音速に整合する右向き純音波。"""
    uinf = 0.5; rh = 220.0 / 600.0; ch = np.sqrt(G * P0 / rh)
    sig = 0.04 / 2.3548

    def init(x):
        pp = amp * P0 * np.exp(-0.5 * ((x - 0.5) / sig) ** 2)
        return rh + pp / ch ** 2, uinf + pp / (rh * ch), P0 + pp
    t_end = 0.5 / (0.5 + ch) + 0.4 / max(ch - 0.5, 1e-3) + 0.1
    ts, pa = run(kind, -1.0, 1.0, dx, init, t_end, 0.8, uinf)
    tl, pl = run(kind, -1.0, 3.0, dx, init, t_end, 0.8, uinf)
    pl_i = np.interp(ts, tl, pl)
    return np.max(np.abs(pa - pl_i)) / np.max(np.abs(pl - P0))


def contact(kind, dx):
    uinf = 0.3
    sig = 0.1 / 2.3548

    def init(x):
        ratio = 1.0 + 1.7 * np.exp(-0.5 * ((x - 0.4) / sig) ** 2)   # 温度比 1 → 2.7
        return 1.0 / ratio, np.full_like(x, uinf), np.full_like(x, P0)
    ts, pa = run(kind, -1.0, 1.0, dx, init, 2.5, 0.9, uinf)
    return np.max(np.abs(pa - P0)) / P0


def main():
    kinds = ["slau_ghost", "hllc_ghost", "roe_ghost", "char_hllc", "char_slau", "hllcd_ghost", "charghost_hllcd"]
    print("== T2 超音速流出の極限 (U_i u=2, U_inf u=3): 面流束 / 内部の物理流束 − 1")
    one = np.array([1.0]); Fi = phys_flux(one, 2 * one, P0 * one)[:, 0]
    for k in kinds:
        F = bflux(k, one, 2 * one, P0 * one, one, 3 * one, P0 * one)[:, 0]
        print(f"  {k:11s} " + "  ".join(f"{(F[j] / Fi[j] - 1):+.2e}" for j in range(3)))
    print("\n== T4 連続性 (U_inf の u = 0、U_i の u を −1.2..1.2 で 1e-4 刻み): 隣接差の最大 / 流束規模")
    uu = np.arange(-1.2, 1.2, 1e-4); onev = np.ones_like(uu)
    for k in kinds:
        F = bflux(k, onev, uu, P0 * onev, onev, 0 * onev, P0 * onev)
        jump = np.max(np.abs(np.diff(F, axis=1)), axis=1) / (np.max(np.abs(F), axis=1) + 1e-30)
        print(f"  {k:11s} " + "  ".join(f"{j:.1e}" for j in jump) + ("   (刻み 1e-4 に対し滑らか)" if np.all(jump < 1e-3) else "   ← 跳びあり"))
    print("\n== T1 音響反射率 (短領域 − 長領域 の最大 / 入射振幅)")
    for M in (0.0, 0.3):
        for dx in (5e-3, 2.5e-3, 1.25e-3):
            row = "  ".join(f"{k}:{acoustic(k, M, dx) * 100:6.2f}%" for k in kinds)
            print(f"  M {M:.1f} Δx {dx:.2e}  {row}")
    print(f"  対照 T5 (M 0、右端 slip): {acoustic('slip', 0.0, 2.5e-3) * 100:.1f}%")
    print("\n== T1b 高温の内部を出ていく音響 (codex plan-7 M1: T 600/220、u 0.5)")
    for dx in (5e-3, 2.5e-3, 1.25e-3):
        print(f"  Δx {dx:.2e}  " + "  ".join(f"{k}:{acoustic_hot(k, dx) * 100:6.2f}%" for k in ("hllc_ghost", "hllcd_ghost", "charghost_hllcd")))
    print("\n== T4b 共通速度を加えた掃引 (内 T 10 倍・P 1e4 倍の極端な比、codex plan-7 M2 型) : 隣接差の最大 / 規模、HLL 退避回数")
    for k in ("hllcd_ghost", "charghost_hllcd"):
        FALLBACK["n"] = 0
        uu = np.arange(-3.0, 3.0, 1e-4); o = np.ones_like(uu)
        F = bflux(k, 0.1 * o, uu, 1e4 * P0 * o, 1e-4 * o, uu, 1e-3 * P0 * o)
        jump = np.max(np.abs(np.diff(F, axis=1)), axis=1) / (np.max(np.abs(F), axis=1) + 1e-30)
        print(f"  {k:16s} " + "  ".join(f"{j:.1e}" for j in jump) + f"  退避 {FALLBACK['n']}  有限 {bool(np.all(np.isfinite(F)))}")
    print("\n== T3 接触波 (温度 2.7 倍の塊が右端から流出、M 0.3): 観測点 |P − P∞|/P∞ の最大")
    for dx in (5e-3, 2.5e-3):
        print(f"  Δx {dx:.2e}  " + "  ".join(f"{k}:{contact(k, dx):.2e}" for k in kinds))


if __name__ == "__main__":
    main()
