#!/usr/bin/env python3
"""lump を含む化学種拡散の縮約規約の参照実装 (double) と判定 (plan thermophysics-solver-owned-species-db §4.4 確定版, #7a, §6 V4a/V4c/V4d)。

ソルバに依存しない独立参照。式は forge の現行演算に合わせる:
  D_rq: Chapman–Enskog + Neufeld Ω(1,1) (thermo_d.cuh thermo_Dbinary; σ_rq = (σ_r+σ_q)/2, ε_rq = √(ε_r ε_q))
  混合平均: D_r = Σ_{q≠r} X_q / Σ_{q≠r} X_q/D_rq (補数形), 純成分は自己拡散 D_rr
  流束: J_s = −ρ D_s ∇Y_s,  補正 J_s* = J_s − Y_s Σ J,  エネルギー q = Σ h_s J_s*
規約 (B, 採用): 輸送種 → 実種モル分率 X_r = Σ_s X_s E_sr、非 lump は D_{r(i)}、lump は D_L = Σ_r a_{r|L} D_r (a = lump 内質量分率)。
対照: full (全実種を輸送種として同じ演算)、A (対称 Blanc: ラベル対の係数 1/D_ij = Σ_{p∈i} Σ_{q∈j} x_p x_q / D_pq)、
      旧 (lump を質量分率平均 LJ の擬似分子として扱う; speciesDB.cpp の PROVISIONAL)。
使い方: python3 test_lump_diffusion_reduction.py  (終了コード 0 = 全判定 PASS)
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "design"))

FAIL = []


def check(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  {detail}" if detail else ""))
    if not ok:
        FAIL.append(name)


# ---------------------------------------------------------------- 物性 (forge と同じ式、double)
def omega11(Ts):
    Ts = min(max(Ts, 0.3), 100.0)
    return (1.06036 * Ts ** -0.15610 + 0.19300 * math.exp(-0.47635 * Ts)
            + 1.03587 * math.exp(-1.52996 * Ts) + 1.76474 * math.exp(-3.89411 * Ts))


def dbinary(a, b, T, P):
    """a, b = (MW [kg/mol], σ [Å], ε/kB [K])。[m²/s]"""
    Mi, Mj = a[0] * 1000.0, b[0] * 1000.0
    sig = 0.5 * (a[1] + b[1])
    eps = math.sqrt(max(a[2] * b[2], 1e-30))
    Dcm2 = 1.8583e-3 * math.sqrt(T ** 3 * (1.0 / Mi + 1.0 / Mj)) / ((P / 101325.0) * sig * sig * omega11(T / eps))
    return Dcm2 * 1e-4


def dmix(Dmat, X, r):
    """forge の補数形混合平均。Dmat[r][q] は二元係数。"""
    num = den = 0.0
    for q in range(len(X)):
        if q == r:
            continue
        num += X[q]
        den += X[q] / max(Dmat[r][q], 1e-300)
    return Dmat[r][r] if den < 1e-300 else num / den


# ---------------------------------------------------------------- モデル
class Mixture:
    """実種 (MW, σ, ε, h) とラベル (lump は実種のモル分率 dict、非 lump は {実種: 1})。"""

    def __init__(self, reals, labels, Dmat=None):
        self.rn = list(reals.keys())
        self.R = [reals[k] for k in self.rn]          # (MW, σ, ε, h)
        self.labels = labels                           # [(name, {real: x})]
        self.Dmat_given = Dmat
        nr = len(self.rn)
        self.E = []
        for _, comp in labels:
            s = sum(comp.values())
            self.E.append([comp.get(k, 0.0) / s for k in self.rn])
        self.ML = [sum(e[r] * self.R[r][0] for r in range(nr)) for e in self.E]
        # lump 内質量分率 a_{r|L}
        self.A = [[self.E[L][r] * self.R[r][0] / self.ML[L] for r in range(nr)] for L in range(len(labels))]
        self.hL = [sum(self.A[L][r] * self.R[r][3] for r in range(nr)) for L in range(len(labels))]

    def Dreal(self, T, P):
        if self.Dmat_given is not None:
            return self.Dmat_given
        n = len(self.R)
        return [[dbinary(self.R[i], self.R[j], T, P) for j in range(n)] for i in range(n)]

    def X_label(self, Y):
        w = [Y[s] / self.ML[s] for s in range(len(Y))]
        t = sum(w)
        return [v / t for v in w]

    def X_real(self, Y):
        XL = self.X_label(Y)
        return [sum(XL[s] * self.E[s][r] for s in range(len(Y))) for r in range(len(self.R))]

    # 規約 B
    def D_B(self, Y, T, P):
        Dm = self.Dreal(T, P)
        Xr = self.X_real(Y)
        Dr = [dmix(Dm, Xr, r) for r in range(len(self.R))]
        return [sum(self.A[s][r] * Dr[r] for r in range(len(self.R))) for s in range(len(Y))]

    # 対照 A: 対称 Blanc のラベル二元係数で混合平均
    def D_A(self, Y, T, P):
        Dm = self.Dreal(T, P)
        n = len(Y)
        DL = [[1.0 / sum(self.E[i][p] * self.E[j][q] / Dm[p][q] for p in range(len(self.R)) for q in range(len(self.R)))
               for j in range(n)] for i in range(n)]
        XL = self.X_label(Y)
        return [dmix(DL, XL, i) for i in range(n)]

    # 対照 旧: 質量分率平均 LJ の擬似分子
    def D_old(self, Y, T, P):
        nr = len(self.R)
        pseudo = []
        for L in range(len(Y)):
            sig = sum(self.A[L][r] * self.R[r][1] for r in range(nr))
            eps = sum(self.A[L][r] * self.R[r][2] for r in range(nr))
            pseudo.append((self.ML[L], sig, eps))
        Dm = [[dbinary(pseudo[i], pseudo[j], T, P) for j in range(len(Y))] for i in range(len(Y))]
        XL = self.X_label(Y)
        return [dmix(Dm, XL, i) for i in range(len(Y))]

    def label_fluxes(self, D, Y, gY, rho=1.0):
        J = [-rho * D[s] * gY[s] for s in range(len(Y))]
        S = sum(J)
        Js = [J[s] - Y[s] * S for s in range(len(Y))]
        q = sum(self.hL[s] * Js[s] for s in range(len(Y)))
        Jreal = [sum(self.A[s][r] * Js[s] for s in range(len(Y))) for r in range(len(self.R))]
        return Js, q, Jreal

    def full_fluxes(self, Y, gY, T, P, rho=1.0):
        nr = len(self.R)
        Yr = [sum(self.A[s][r] * Y[s] for s in range(len(Y))) for r in range(nr)]
        gYr = [sum(self.A[s][r] * gY[s] for s in range(len(Y))) for r in range(nr)]
        Dm = self.Dreal(T, P)
        Xr = self.X_real(Y)
        Dr = [dmix(Dm, Xr, r) for r in range(nr)]
        J = [-rho * Dr[r] * gYr[r] for r in range(nr)]
        S = sum(J)
        Js = [J[r] - Yr[r] * S for r in range(nr)]
        q = sum(self.R[r][3] * Js[r] for r in range(nr))
        return Js, q


def relerr(a, b, scale):
    return abs(a - b) / max(abs(scale), 1e-300)


def norm(v):
    return math.sqrt(sum(x * x for x in v))


# ---------------------------------------------------------------- V4a: codex の合成例 (非重複)
def v4a_synthetic():
    print("V4a (合成例: MW (2,4,8), X (0.4,0.4,0.2), D12 1 / D13 2 / D23 4, lump = 種 1・2 等モル, h (1,5,0))")
    reals = {"r1": (0.002, 1, 1, 1.0), "r2": (0.004, 1, 1, 5.0), "r3": (0.008, 1, 1, 0.0)}
    Dm = [[1.0, 1.0, 2.0], [1.0, 1.0, 4.0], [2.0, 4.0, 1.0]]
    mix = Mixture(reals, [("L12", {"r1": 1.0, "r2": 1.0}), ("r3", {"r3": 1.0})], Dmat=Dm)
    Y = [0.6, 0.4]
    gY = [1.0, -1.0]
    Xr = mix.X_real(Y)
    check("実種モル分率 = (0.4, 0.4, 0.2)", max(abs(a - b) for a, b in zip(Xr, [0.4, 0.4, 0.2])) < 1e-12, str([round(x, 12) for x in Xr]))
    Jf, qf = mix.full_fluxes(Y, gY, 0, 0)
    JB, qB, JBr = mix.label_fluxes(mix.D_B(Y, 0, 0), Y, gY)
    JA, qA, _ = mix.label_fluxes(mix.D_A(Y, 0, 0), Y, gY)
    scale = norm(Jf)
    print(f"    外部種 J3*: full {Jf[2]:.9f}  B {JB[1]:.9f}  A {JA[1]:.9f}   q: full {qf:.9f}  B {qB:.9f}  A {qA:.9f}")
    check("B: 外部種の補正後流束 = full (≤1e-12)", relerr(JB[1], Jf[2], scale) <= 1e-12, f"{relerr(JB[1], Jf[2], scale):.2e}")
    check("B: lump の補正後流束 = 構成実種の和 (≤1e-12)", relerr(JB[0], Jf[0] + Jf[1], scale) <= 1e-12, f"{relerr(JB[0], Jf[0] + Jf[1], scale):.2e}")
    check("A: 外部種の補正後流束は full と不一致 (判別)", relerr(JA[1], Jf[2], scale) > 1e-6, f"{relerr(JA[1], Jf[2], scale):.2e}")
    # codex の独立計算 (notes/reviews/2026-10-05-lump-blanc-diffusion-diagnose.md) の規格化値と照合: ρ=1・D0=1・∇Y_L=1
    check("codex の値の再現 (full J3* 2.115555556, B q −7.757037037, A J3* 2.666666667)",
          abs(Jf[2] - 2.115555556) < 5e-9 and abs(qB + 7.757037037) < 5e-9 and abs(JA[1] - 2.666666667) < 5e-9 and abs(qf + 7.875555556) < 5e-9,
          f"full J3* {Jf[2]:.9f} q {qf:.9f}, B q {qB:.9f}, A J3* {JA[1]:.9f}")
    print(f"    V4c 記録: エネルギーの共分散項 q_B − q_full = {qB - qf:+.9f} (規格化)")


# ---------------------------------------------------------------- 実在種のデータ (設計側の共通データ; 既定 ljSource)
def real_db(names):
    from forge_design.gas import semiperfect as sp
    lj = sp.lj_params()
    out = {}
    for k in names:
        MW = sp.SPECIES_NASA9[k]["MW"]
        sig, eps = lj[k]
        out[k] = [MW, sig, eps, 0.0]
    return out


def h_mass(name, T):
    from forge_design.gas import semiperfect as sp
    try:
        return sp.h_mass(name, T)
    except Exception:
        return None


def v4a_va3():
    print("V4a (case/44 va3: MIXDRY = N2/O2/Ar/CO2 + H2O、非重複、T 300/1000 K、1 atm)")
    comp = {"N2": 6.64860e-1, "O2": 2.16072e-1, "AR": 7.97588e-3, "CO2": 4.90034e-2}
    db = real_db(["N2", "O2", "AR", "CO2", "H2O"])
    for k in db:
        db[k][3] = (db[k][0] * 1e3) * 0.37 + len(k)   # 共分散項を出すための任意の h (判定は恒等式のみ)
    mix = Mixture({k: tuple(v) for k, v in db.items()}, [("MIXDRY", comp), ("H2O", {"H2O": 1.0})])
    for T in (300.0, 1000.0):
        for Yw in (0.04, 0.2):
            Y = [1 - Yw, Yw]
            gY = [-1.0, 1.0]
            Jf, qf = mix.full_fluxes(Y, gY, T, 101325.0)
            JB, qB, _ = mix.label_fluxes(mix.D_B(Y, T, 101325.0), Y, gY)
            JO, _, _ = mix.label_fluxes(mix.D_old(Y, T, 101325.0), Y, gY)
            scale = norm(Jf)
            e_ext = relerr(JB[1], Jf[4], scale)
            e_lump = relerr(JB[0], sum(Jf[:4]), scale)
            check(f"T {T:.0f} Y_H2O {Yw}: B 外部種 = full ({e_ext:.1e})、lump = 構成和 ({e_lump:.1e})", e_ext <= 1e-12 and e_lump <= 1e-12)
            DB = mix.D_B(Y, T, 101325.0)[1]
            DO = mix.D_old(Y, T, 101325.0)[1]
            print(f"      V4g 記録: H2O の D 新/旧 − 1 = {DB / DO - 1:+.4%}   補正後 H2O 流束 旧 vs full {relerr(JO[1], Jf[4], abs(Jf[4])):+.3%}")


# ---------------------------------------------------------------- V4d: SERN (重複) — 新規約は旧平均 LJ より悪くない
def v4d_sern():
    print("V4d (SERN m6_on: EXH (CEA 凍結 11 種) と AMB (乾燥空気 4 種) が N2/O2/Ar/CO2 を共有、∇Y_EXH = 1)")
    exh = {'N2': 0.6389, 'H2O': 0.32695, 'H2': 0.01251, 'AR': 0.00767, 'OH': 0.00604, 'O2': 0.00398, 'NO': 0.00206,
           'H': 0.00125, 'O': 0.00037, 'CO2': 0.00021, 'CO': 0.00005}
    from forge_design.gas.frozen import AIR_MOLE
    amb = {k.upper(): v for k, v in AIR_MOLE.items()}
    names = sorted(set(exh) | set(amb))
    db = real_db(names)
    mix = Mixture({k: tuple(v) for k, v in db.items()}, [("EXH", exh), ("AMB", amb)])
    rows = []
    allok = True
    for T in (300.0, 1000.0, 2000.0):
        for Ye in (0.1, 0.3, 0.5, 0.7, 0.9):
            Y = [Ye, 1 - Ye]
            gY = [1.0, -1.0]
            Jf, _ = mix.full_fluxes(Y, gY, T, 101325.0)
            _, _, JBr = mix.label_fluxes(mix.D_B(Y, T, 101325.0), Y, gY)
            _, _, JOr = mix.label_fluxes(mix.D_old(Y, T, 101325.0), Y, gY)
            nf = norm(Jf)
            eB = norm([a - b for a, b in zip(JBr, Jf)]) / nf
            eO = norm([a - b for a, b in zip(JOr, Jf)]) / nf
            DB = mix.D_B(Y, T, 101325.0)
            DO = mix.D_old(Y, T, 101325.0)
            rows.append((T, Ye, eB, eO, DB[0] / DO[0] - 1, DB[1] / DO[1] - 1))
            allok &= eB <= eO
    print("      T [K]  Y_EXH   ε_新 (B)    ε_旧 (平均 LJ)   D_EXH 新/旧−1   D_AMB 新/旧−1")
    for T, Ye, eB, eO, d0, d1 in rows:
        print(f"      {T:6.0f}  {Ye:4.1f}   {eB:9.4f}   {eO:9.4f}        {d0:+8.3%}       {d1:+8.3%}")
    check("全 15 状態で ε_新 ≤ ε_旧 (事前登録 V4d)", allok)


if __name__ == "__main__":
    v4a_synthetic()
    v4a_va3()
    v4d_sern()
    print("\nALL PASS" if not FAIL else f"\nFAIL: {len(FAIL)} 件: " + "; ".join(FAIL))
    sys.exit(0 if not FAIL else 1)
