#!/usr/bin/env python3
"""V2d の独立参照解 (plan boundary-node-farfield-characteristic §6 V2d)。forge と独立の 1 次元 Euler (セル中心 FV)。

forge の node 離散化・境界とは独立に書く: セル中心、内部面は HLLC (Davis 波速)、原始変数 (ρ, u, P, Y) の MUSCL (MC 制限)、
SSP-RK3、両端は外挿 (領域を十分長く取り、評価時間内に端の影響が評価点へ届かないようにする)。熱物性だけは
`farfield_proto1d_tp.py` の `FrozenGas` の表 (EXH / AIR、NASA-9、href 298.15 K) と CPG (γ 1.4、cp 1004.5) を使う。

  --case contact  : V2d-1。P∞ 2851 Pa、u = 0.5 c∞、中心 0.5 m の T 600 K・Y_EXH 0.13 (tp2) のガウス塊 (FWHM 0.1 m)。評価点 x 0.975
  --case acoustic : V2d-2 (隔離配置)。内部 = 一様 T 600 K・Y 0.13 (tp2)、u = 0.3 c_i、中心 0.4 m の右向きパルス (1e-3 P∞、FWHM 0.04 m)。評価点 x 0.8
  --gas cpg | tp1 | tp2、--dx (既定 1.25e-3。forge の Δx 5 mm では参照自身が収束しない)。出力 = 評価点の時系列 CSV (t, P, T, Y) と、--dx /4 との差 (参照自身の収束確認) を --check で。

  cd <repo>/design && python3 ../solver_density_cuda/tools/ref1d_euler_tp.py --case contact --gas tp2 --check --out ref_contact_tp2.csv
"""
import argparse, math, sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import farfield_proto1d_tp as TP  # noqa: E402  (熱物性の表だけを使う)

P0, T0 = 2851.0, 220.0


class Gas:
    def __init__(self, kind):
        self.kind = kind

    def R(self, Y):
        return (1004.5 - 1004.5 / 1.4) + 0 * Y if self.kind == "cpg" else TP.mix(Y if self.kind == "tp2" else 0 * Y)

    def e(self, T, Y):
        return 1004.5 / 1.4 * T if self.kind == "cpg" else TP.e_of(T, Y if self.kind == "tp2" else 0 * Y)

    def T_from_e(self, e, Y, T0=None):
        return e / (1004.5 / 1.4) if self.kind == "cpg" else TP.T_from_e(e, Y if self.kind == "tp2" else 0 * Y, T0)

    def c(self, T, Y):
        return np.sqrt(1.4 * self.R(Y) * T) if self.kind == "cpg" else TP.snd(T, Y if self.kind == "tp2" else 0 * Y)


def run(g, case, dx, cfl=0.25):
    x = np.arange(-0.5 + 0.5 * dx, 2.0, dx)   # 評価時間内に端の影響が評価点へ届かない長さ
    Tinf = T0; ro_inf = P0 / (g.R(np.array([0.0])) * Tinf)[0]; c_inf = g.c(np.array([Tinf]), np.array([0.0]))[0]
    Yi = 0.13 if g.kind == "tp2" else 0.0
    if case == "contact":
        sig = 0.1 / (2 * math.sqrt(2 * math.log(2))); w = np.exp(-0.5 * ((x - 0.5) / sig) ** 2)
        T = Tinf + (600.0 - Tinf) * w; Y = Yi * w; P = np.full_like(x, P0); U = 0.5 * c_inf; u = np.full_like(x, U)
        tend = (0.475 + 0.15) / U; xe = 0.975
    else:
        ci = g.c(np.array([600.0]), np.array([Yi]))[0]; roi = P0 / (g.R(np.array([Yi])) * 600.0)[0]
        U = 0.3 * ci; sig = 0.04 / (2 * math.sqrt(2 * math.log(2)))
        pp = 1e-3 * P0 * np.exp(-0.5 * ((x - 0.4) / sig) ** 2)
        P = P0 + pp; u = U + pp / (roi * ci); Y = np.full_like(x, Yi)
        gi = ci * ci * roi / P0; T = 600.0 * (P / P0) ** ((gi - 1) / gi)
        tend = 0.6 / (ci + U) + 0.2 / (ci - U) + (3 * sig + 0.02) / (ci - U); xe = 0.8
    r = P / (g.R(Y) * T)
    Q = np.stack([r, r * u, r * (g.e(T, Y) + 0.5 * u * u), r * Y])
    ie = int(np.argmin(np.abs(x - xe)))
    Tc = T.copy()

    def prim(Q):
        nonlocal Tc
        r = Q[0]; u = Q[1] / r; Y = np.clip(Q[3] / r, 0.0, 1.0); e = Q[2] / r - 0.5 * u * u
        Tc = g.T_from_e(e, Y, Tc); P = r * g.R(Y) * Tc
        return r, u, P, Y, Tc

    def flux(rL, uL, pL, YL, rR, uR, pR, YR):
        TL = pL / (rL * g.R(YL)); TR = pR / (rR * g.R(YR)); cL = g.c(TL, YL); cR = g.c(TR, YR)
        EL = g.e(TL, YL) + 0.5 * uL * uL; ER = g.e(TR, YR) + 0.5 * uR * uR
        SL = np.minimum(uL - cL, uR - cR); SR = np.maximum(uL + cL, uR + cR)
        Sm = (pR - pL + rL * uL * (SL - uL) - rR * uR * (SR - uR)) / (rL * (SL - uL) - rR * (SR - uR))
        FL = np.stack([rL * uL, rL * uL * uL + pL, (rL * EL + pL) * uL, rL * uL * YL])
        FR = np.stack([rR * uR, rR * uR * uR + pR, (rR * ER + pR) * uR, rR * uR * YR])
        UL = np.stack([rL, rL * uL, rL * EL, rL * YL]); UR = np.stack([rR, rR * uR, rR * ER, rR * YR])

        def star(r, u, p, E, Y, S):
            f = r * (S - u) / (S - Sm)
            return np.stack([f, f * Sm, f * (E + (Sm - u) * (Sm + p / (r * (S - u)))), f * Y])
        return np.where(SL >= 0, FL, np.where(Sm >= 0, FL + SL * (star(rL, uL, pL, EL, YL, SL) - UL),
                        np.where(SR > 0, FR + SR * (star(rR, uR, pR, ER, YR, SR) - UR), FR)))

    def mm(a, b):   # MC (monotonized central): minmod は Δx 5 mm で散逸が大きく参照の収束条件を満たさなかった
        return np.where(a * b > 0, np.sign(a) * np.minimum(np.minimum(2 * np.abs(a), 2 * np.abs(b)), 0.5 * np.abs(a + b)), 0.0)

    def rhs(Q):
        W = np.stack(prim(Q)[:4])
        We = np.concatenate([W[:, :1], W[:, :1], W, W[:, -1:], W[:, -1:]], 1)   # 外挿 (2 ゴースト)
        d = mm(We[:, 1:-1] - We[:, :-2], We[:, 2:] - We[:, 1:-1])               # 各セル (ゴースト 1 つ込み) の勾配
        L = We[:, 1:-2] + 0.5 * d[:, :-1]; R = We[:, 2:-1] - 0.5 * d[:, 1:]    # 面 i-1/2 .. i+1/2
        F = flux(L[0], L[1], L[2], L[3], R[0], R[1], R[2], R[3])
        return -(F[:, 1:] - F[:, :-1]) / dx

    t = 0.0; out = []
    while t < tend:
        r_, u_, P_, Y_, T_ = prim(Q)
        dt = min(cfl * dx / np.max(np.abs(u_) + g.c(T_, Y_)), tend - t)
        out.append((t, P_[ie], T_[ie], Y_[ie]))
        Q1 = Q + dt * rhs(Q)
        Q2 = 0.75 * Q + 0.25 * (Q1 + dt * rhs(Q1))
        Q = Q / 3.0 + 2.0 / 3.0 * (Q2 + dt * rhs(Q2))
        t += dt
    r_, u_, P_, Y_, T_ = prim(Q); out.append((t, P_[ie], T_[ie], Y_[ie]))
    return np.array(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", choices=("contact", "acoustic"), required=True); ap.add_argument("--gas", choices=("cpg", "tp1", "tp2"), required=True)
    ap.add_argument("--dx", type=float, default=1.25e-3); ap.add_argument("--check", action="store_true"); ap.add_argument("--out", default=None)
    a = ap.parse_args()
    g = Gas(a.gas)
    s = run(g, a.case, a.dx)
    if a.out:
        np.savetxt(a.out, s, delimiter=",", header="t,P,T,Y", comments="")
    msg = f"{a.case} {a.gas} dx {a.dx}: 評価点 max|P−P∞| {np.max(np.abs(s[:, 1] - P0)):.4g} Pa、max T {np.max(s[:, 2]):.2f} K、max Y {np.max(s[:, 3]):.4f}"
    if a.check:
        f = run(g, a.case, a.dx / 4)
        dP = np.max(np.abs(s[:, 1] - np.interp(s[:, 0], f[:, 0], f[:, 1])))
        dT = np.max(np.abs(s[:, 2] - np.interp(s[:, 0], f[:, 0], f[:, 2])))
        dY = np.max(np.abs(s[:, 3] - np.interp(s[:, 0], f[:, 0], f[:, 3])))
        if a.case == "contact":
            tol = {"P": 1e-3 * P0, "T": 0.02 * 380.0, "Y": 0.002}
        else:
            tol = {"P": 0.05 * 1e-3 * P0, "T": np.inf, "Y": np.inf}
        ok = dP <= tol["P"] / 5 and dT <= tol["T"] / 5 and dY <= tol["Y"] / 5
        msg += (f"\n  Δx/4 との差: P {dP:.3g} Pa (許容の 1/5 = {tol['P'] / 5:.3g})、T {dT:.3g} K (1/5 = {tol['T'] / 5:.3g})、Y {dY:.2e} (1/5 = {tol['Y'] / 5:.3g})"
                f" → 参照の収束 {'OK' if ok else 'NG'}")
    print(msg)


if __name__ == "__main__":
    main()
