#!/usr/bin/env python3
"""二相拡散 (plans/active/condensation-two-phase-transport.md §4.2) のホスト参照実装と §6 単体合格条件の判定 (#4 の諮問材料; CFD 0 step)。

cuda_forge のカーネルは使わない。§4.2 の離散作用素 (二点差分) と更新 (蒸気/液を変数にした点対角の前処理つき固定点反復) を
numpy で参照実装し、§6 の事前固定の数値で判定する。カーネル実装の前に、方針どおりの作用素と更新で合格条件が満たせるか
(と、満たせない場合はどこか) を上位に諮るための材料。

  python3 solver_density_cuda/tests/unit/test_twophase_diffusion_harness.py [--quick]

記号 (§4.2): ρ_g = ρ(1−g)、気相内の質量分率 z_k = ρY_k/ρ_g (k ≠ 水)、z_v = (ρY_w − ρg)/ρ_g (Σz で正規化)。
  分子 (気相のみ): j_k⁰ = −ρ_g,f D_k,f (z_k,R − z_k,L) a_f (a_f = A/d)、補正 j_k = j_k⁰ − z_k,f Σ_j j_j⁰。
  乱流 (全相共通): J_t(φ) = −Γ_f (φ_R − φ_L) a_f (Γ = μ_t/Sc_t、φ = Y_k, Y_w, g, Q_n, 蒸気 v = Y_w − g)。
  組み立て: J_k = j_k + J_t(Y_k)、J_v = j_v + J_t(v)、J_l = J_t(g)、J_w = J_v + J_l (蒸気を先に作る; 旧順 J_v = J_w − J_l も比較)。
  エネルギー: q = Σ_{k≠水} h_k(T_f) J_k + h_v(T_f) J_w − L(T_f) J_l。
  EOS (forge の carrier 二相 EOS, cuda_forge/condensationEOS_d.cuh cond_twophase_resid と同じ式):
    e = Σ_{k≠水} Y_k (h_k − R_k T) + Y_w (h_v − R_w T) + g (R_w T − L(T))   (総水分を蒸気として数えた e_gas + g(R_wT − L))。
  更新 (契約 = 保存形後退 Euler の固定点): (V/Δt)(ρφ − ρφⁿ) = R(φ)。点対角 (蒸気: 分子 + 乱流、液・Q: 乱流) は反復の前処理。
    増分は蒸気 ρv と液 ρg を変数に作り Δ(ρY_w) = Δρv + Δρg で戻す (writeback 'increment')。比較: ρY_w = ρg' + ρv' ('rebuild')。
  再正規化の係数は ρg と Q にも掛ける。実現可能性クランプ (0 ≤ ρg ≤ ρY_w, Q ≥ 0, ρY ≥ 0) は保険として最後に掛け、補正量を数える。

試験 (plan §6 の事前固定の数値; 判定語 [PASS]/[FAIL]、[INFO] は判定外の記録):
  S1 分子流束の構造 (i)  三成分 (等分子量, z = [0.2,0.3,0.5], 二元 D_12/D_13/D_23 = [1,2,3]e-5 → 混合平均) の気相一様・∇g = 0.1 で
                          |j_k| ≤ 1e-6·ρ_g·max D·|∇g| (float64 / float32)。旧案 (−ρD_k∇Y_k + 気相補正) の偽流束 ~1.3e-7 m/s を併記。
  S2 面恒等式 (ii)       ランダム面で |Σ_気相 j_k| ≤ 8ε₃₂ Σ|j_k⁰|、|J_w − J_v − J_l| ≤ 8ε₃₂(|J_w|+|J_v|+|J_l|) (float32 で計算、float64 で評価)。
  S3 3 セル判別         閉じた直列 3 セル、ρ=V=Δt=1、乱流係数 1、分子・移流・相変化なし、g=[0.1,0.2,0.2]、Y_w=g+0.1。
                          A = 点対角 1 回更新、B = 解き切った後退 Euler (BE 残差 ≤1e-7 相対)。B: 総液量・各種総量の相対変化 ≤1e-6
                          (1・1000 更新、再正規化・クランプ前、float64 集計)、全セル ρv ≥ 0、|ΣρY − ρ| ≤ 8ε₃₂ρ。規則: A FAIL・B PASS で契約を支持。
  S4 非負               2 セル例 (蒸気 0、Y_w = g = 0.1/0.2、面係数・V/Δτ = 1) の解き切った更新で ρv ≥ 0 (クランプ 0)。
                          蒸気 0 近傍のランダム 1D 100 セル × 10⁴ 擬似反復でクランプ補正量 0。
  S5 保存               閉じた箱・拡散単独の解き切った物理時間更新 1000 回で Σ ρY_k V・Σ ρg V・Σ ρQ V 相対 ≤1e-6、Σ ρE V は初期 Σ|ρE|V
                          分母で ≤1e-6、初期総量 0 の量は絶対 ≤1e-12。
  S6 エネルギー         分子拡散 0・総組成一様・乱流のみ・g のみ勾配の 1D で、1 更新の温度変化の EOS 接線投影が −R_wT r_g Δt/(ρ c_v,eff) と ≤1 %
                          (有限更新の温度比較は +1e-3 K)。判別: 潜熱流束なし / エネルギーだけ潜熱 (液を拡散しない) は FAIL すること。
"""
import argparse, math, os, sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TOOLS = os.path.normpath(os.path.join(HERE, "..", "..", "tools"))
sys.path.insert(0, HERE)
sys.path.insert(0, TOOLS)
from total_quantities import _TPGas      # noqa: E402  forge と同じ NASA-9 評価 (範囲外は定 cp 外挿)
from test_tpgas_lowT import DB as NASA_DB  # noqa: E402  N2 / H2O の NASA-9 係数

RU = 8.314462618
EPS32 = float(np.finfo(np.float32).eps)
g_fail = 0


def verdict(ok, what):
    global g_fail
    print(("[PASS] " if ok else "[FAIL] ") + what)
    if not ok:
        g_fail += 1
    return ok


def info(what):
    print("[INFO] " + what)


# ===================================================================== 熱物性
class Species:
    """気相種 1 つ。h(T) は質量基準 (datum 込み), cp(T), R。NASA-9 (_TPGas) か定 cp。"""

    def __init__(self, name, R, h, cp):
        self.name, self.R, self._h, self._cp = name, R, h, cp

    def h(self, T):
        T = np.asarray(T, dtype=np.float64)
        return np.asarray(self._h(np.atleast_1d(T)), dtype=np.float64).reshape(T.shape)

    def cp(self, T):
        T = np.asarray(T, dtype=np.float64)
        return np.asarray(self._cp(np.atleast_1d(T)), dtype=np.float64).reshape(T.shape)


def nasa_species(name):
    gas = _TPGas(NASA_DB, [name], 0.0)   # Tref=0: datum 込み (forge thermo_h_mass と同じ)
    R = RU / NASA_DB[name]["MW"]
    return Species(name, R, lambda T: gas.h([1.0], T), lambda T: gas.cp([1.0], T))


def cpg_species(name, MW, cp, h298):
    R = RU / MW
    return Species(name, R, lambda T: h298 + cp * (T - 298.15), lambda T: cp + 0.0 * T)


class Latent:
    """L(T) = h_v(T) − h_l(T)、h_l = h_v(T0) − L0 + c_l (T − T0) (forge の h2o_latent_pair と同じ「h_v − h_l」の構造。
    液の h は簡易形で代用: 接線投影の検証には L の具体形は効かず、L と L' が EOS と流束で同じ関数であることだけが要る)。"""

    def __init__(self, vapor, L0=2.501e6, T0=273.15, cl=4186.0):
        self.v, self.cl, self.T0 = vapor, cl, T0
        self.hl0 = float(vapor.h(T0)) - L0

    def L(self, T):
        return self.v.h(T) - (self.hl0 + self.cl * (np.asarray(T, dtype=np.float64) - self.T0))

    def dL(self, T):
        return self.v.cp(T) - self.cl


class Thermo:
    def __init__(self, nonwater, vapor, latent):
        self.nw, self.v, self.lat = nonwater, vapor, latent

    def e_mix(self, T, Yk, Yw, g):
        e = Yw * (self.v.h(T) - self.v.R * T) + g * (self.v.R * T - self.lat.L(T))
        for k, sp in enumerate(self.nw):
            e = e + Yk[k] * (sp.h(T) - sp.R * T)
        return e

    def cv_eff(self, T, Yk, Yw, g):
        c = Yw * (self.v.cp(T) - self.v.R) + g * (self.v.R - self.lat.dL(T))
        for k, sp in enumerate(self.nw):
            c = c + Yk[k] * (sp.cp(T) - sp.R)
        return c

    def T_from_state(self, st, T0):
        """double Newton (forge の cond_twophase_polish と同じく残差 1e-9|e| + 0.05 J/kg まで)。"""
        rho = st["rho"].astype(np.float64)
        Yk = st["rY"].astype(np.float64) / rho
        Yw = st["rYw"].astype(np.float64) / rho
        g = st["rg"].astype(np.float64) / rho
        e = st["rE"].astype(np.float64) / rho
        T = np.array(T0, dtype=np.float64, copy=True)
        for _ in range(50):
            G = self.e_mix(T, Yk, Yw, g) - e
            if np.all(np.abs(G) <= 1e-12 * np.abs(e) + 1e-6):
                break
            T = T - G / self.cv_eff(T, Yk, Yw, g)
        return T


# ===================================================================== 離散作用素
def mixavg_D(X, Dbin):
    """混合平均 (補数形, #3b): D_k = Σ_{j≠k} X_j / Σ_{j≠k} X_j/D_kj。X: (n, nf)。純成分は二元 D の平均で代用。"""
    n = X.shape[0]
    D = np.empty_like(X)
    for k in range(n):
        num = np.zeros_like(X[0])
        den = np.zeros_like(X[0])
        for j in range(n):
            if j == k:
                continue
            num = num + X[j]
            den = den + X[j] / X.dtype.type(Dbin[k, j])
        fb = X.dtype.type(np.mean([Dbin[k, j] for j in range(n) if j != k]))
        D[k] = np.where(den > 1e-30, num / np.where(den > 1e-30, den, 1), fb)
    return D


class Problem:
    """閉じた 1D 列 (N セル, 内部面 N−1)。face 係数 a (=A/d), Γ (=μ_t/Sc_t)、分子は二元 D 行列 × スケール molc (=0 で分子なし)。
    気相種は [非水 0..K-1, 蒸気] の順。MW は分子量 (X の換算)。"""

    def __init__(self, N, a, Gam, Dbin, MW, V, molc=1.0, thermo=None):
        self.N, self.a, self.Gam = N, np.asarray(a, float), np.asarray(Gam, float)
        self.Dbin, self.MW, self.V, self.molc, self.thermo = np.asarray(Dbin, float), np.asarray(MW, float), np.asarray(V, float), molc, thermo
        self.L = np.arange(N - 1)
        self.R = self.L + 1


def make_state(rho, rY, rYw, rg, rQ, rE, dt):
    return {"rho": np.asarray(rho, dt), "rY": np.atleast_2d(np.asarray(rY, dt)), "rYw": np.asarray(rYw, dt),
            "rg": np.asarray(rg, dt), "rQ": np.atleast_2d(np.asarray(rQ, dt)), "rE": np.asarray(rE, dt)}


def copy_state(st):
    return {k: v.copy() for k, v in st.items()}


def cast_state(st, dt):
    return {k: v.astype(dt) for k, v in st.items()}


def gas_composition(st, dt):
    rhog = st["rho"] - st["rg"]
    rv = st["rYw"] - st["rg"]
    z = np.vstack([st["rY"], rv[None, :]]) / rhog
    z = z / z.sum(axis=0)
    return z, rhog, rv


def fluxes(st, T, pb, dt, corr="mean", order="vapor", design="new", energy="full"):
    """面流束 (向き L→R が正)。design 'new' = §4.2、'old' = 改訂前案 (混合物基準 −ρD_k∇Y_k [水は蒸気 Y_w − g] + 気相補正)。
    energy: 'full' = §4.2、'no_latent' = −L J_l を入れない、'latent_only' = 液を拡散させずエネルギーにだけ L·Γ∇g (禁止案)。"""
    L, R = pb.L, pb.R
    a = pb.a.astype(dt)
    G = pb.Gam.astype(dt)
    half = dt(0.5)
    rho = st["rho"]
    z, rhog, rv = gas_composition(st, dt)
    nG = z.shape[0]
    K = nG - 1
    # 分子 (気相内)
    zf = half * (z[:, L] + z[:, R])
    Xf = (zf / pb.MW.astype(dt)[:, None])
    Xf = Xf / Xf.sum(axis=0)
    Df = mixavg_D(Xf, pb.Dbin * pb.molc) if pb.molc > 0 else np.zeros_like(Xf)
    if design == "new":
        rhog_f = half * (rhog[L] + rhog[R])
        j0 = -rhog_f * Df * (z[:, R] - z[:, L]) * a
    else:
        # 旧案: 気相種の混合物基準分率 Y_k (水は蒸気 Y_w − g) を ρ_f で拡散
        Ymix = np.vstack([st["rY"], rv[None, :]]) / rho
        rho_f = half * (rho[L] + rho[R])
        j0 = -rho_f * Df * (Ymix[:, R] - Ymix[:, L]) * a
    S = j0.sum(axis=0)
    if corr == "mean":
        zc = zf
    else:   # 'upwind': 補正の質量流束 −Σj⁰ (L→R が正) が出ていくセルの z
        zc = np.where(S <= 0, z[:, L], z[:, R])
    j = j0 - zc * S
    # 乱流 (全相共通)
    def Jt(phi):
        return -G * (phi[..., R] - phi[..., L]) * a
    Yk = st["rY"] / rho
    Yw = st["rYw"] / rho
    g = st["rg"] / rho
    v = rv / rho
    Q = st["rQ"] / rho
    Jk = j[:K] + Jt(Yk)
    Jl = Jt(g)
    if energy == "latent_only":
        # 禁止案: 液は拡散させず (総水分は現行どおり ∇Y_w)、エネルギーにだけ潜熱流束 −L Γ∇g を足す
        Jw = j[K] + Jt(Yw)
        Jv = Jw
        Jl_mass = np.zeros_like(Jl)
        JQ = np.zeros_like(Q[..., :len(L)])
    else:
        if order == "vapor":
            Jv = j[K] + Jt(v)
            Jw = Jv + Jl
        else:
            Jw = j[K] + Jt(Yw)
            Jv = Jw - Jl
        Jl_mass = Jl
        JQ = Jt(Q)
    out = {"j0": j0, "j": j, "S": S, "Jk": Jk, "Jv": Jv, "Jw": Jw, "Jl": Jl_mass, "JQ": JQ, "Df": Df}
    if pb.thermo is not None and T is not None:
        th = pb.thermo
        Tf = 0.5 * (T[L] + T[R])
        q = th.v.h(Tf) * Jw.astype(np.float64)
        for k, sp in enumerate(th.nw):
            q = q + sp.h(Tf) * Jk[k].astype(np.float64)
        if energy in ("full", "latent_only"):
            q = q - th.lat.L(Tf) * Jl.astype(np.float64)
        out["q"] = q.astype(dt)
    return out


def residuals(F, pb, dt):
    """R_i = Σ (面からの流入)。面流束 F (L→R) は L から −F、R へ +F。"""
    N = pb.N

    def div(Fx):
        Fx = np.atleast_2d(Fx)
        r = np.zeros((Fx.shape[0], N), dtype=dt)
        np.add.at(r, (slice(None), pb.L), -Fx)
        np.add.at(r, (slice(None), pb.R), Fx)
        return r
    out = {"k": div(F["Jk"]), "v": div(F["Jv"])[0], "l": div(F["Jl"])[0], "w": div(F["Jw"])[0], "Q": div(F["JQ"])}
    if "q" in F:
        out["E"] = div(F["q"])[0]
    return out


def diagonals(st, F, pb, dt, corr="mean", design="new"):
    """点対角 (前処理): ∂(流出)/∂(ρφ_i)。蒸気・非水: 分子 a ρ_g,f D_f/ρ_g,i + 補正の流出 + 乱流 aΓ/ρ_i、液・Q: 乱流 aΓ/ρ_i。"""
    N = pb.N
    L, R = pb.L, pb.R
    a = pb.a.astype(dt)
    G = pb.Gam.astype(dt)
    half = dt(0.5)
    rho = st["rho"]
    rhog = rho - st["rg"]
    rhog_f = half * (rhog[L] + rhog[R])
    turb = np.zeros(N, dt)
    np.add.at(turb, L, a * G / rho[L])
    np.add.at(turb, R, a * G / rho[R])
    nG = F["Df"].shape[0]
    mol = np.zeros((nG, N), dt)
    S = F["S"]
    rho_f = half * (rho[L] + rho[R])
    for k in range(nG):
        if design == "new":
            cL = a * rhog_f * F["Df"][k] / rhog[L]
            cR = a * rhog_f * F["Df"][k] / rhog[R]
        else:   # 改訂前: 混合物基準 ρ_f D/ρ_i
            cL = a * rho_f * F["Df"][k] / rho[L]
            cR = a * rho_f * F["Df"][k] / rho[R]
        # 補正 −z_c S の流出 (mean: 自セル分 1/2、upwind: 流出側の風上セル)。改訂前案は現行 species_diffusion_d と同じく補正を対角に入れない
        if design != "new":
            pass
        elif corr == "mean":   # 補正流束 −z_f S (L→R) の自セル分: L からの流出は −S > 0
            cL = cL + half * np.maximum(-S, 0) / rhog[L]
            cR = cR + half * np.maximum(S, 0) / rhog[R]
        else:
            cL = cL + np.maximum(-S, 0) / rhog[L]
            cR = cR + np.maximum(S, 0) / rhog[R]
        np.add.at(mol[k], L, cL)
        np.add.at(mol[k], R, cR)
    return {"k": mol[:-1] + turb, "v": mol[-1] + turb, "w": mol[-1] + turb, "l": turb}


def jacobi_update(st, st_n, Rr, d, M, dt, vars_mode="vl", writeback="increment", has_energy=False):
    """BE 残差 r = R − M(ρφ − ρφⁿ) を点対角 (M + d) で割った 1 回の更新 (定常の擬似反復なら st_n = st)。"""
    new = copy_state(st)
    rv = st["rYw"] - st["rg"]
    rv_n = st_n["rYw"] - st_n["rg"]
    r_k = Rr["k"] - M * (st["rY"] - st_n["rY"])
    r_l = Rr["l"] - M * (st["rg"] - st_n["rg"])
    r_Q = Rr["Q"] - M * (st["rQ"] - st_n["rQ"])
    new["rY"] = st["rY"] + r_k / (M + d["k"])
    dl = r_l / (M + d["l"])
    new["rg"] = st["rg"] + dl
    new["rQ"] = st["rQ"] + r_Q / (M + d["l"])
    if vars_mode == "vl":
        r_v = Rr["v"] - M * (rv - rv_n)
        dv = r_v / (M + d["v"])
        if writeback == "increment":
            new["rYw"] = st["rYw"] + (dv + dl)
        else:
            new["rYw"] = new["rg"] + (rv + dv)
    else:   # 'wl': 総水分と液をそれぞれの点対角で (改訂前の更新)
        r_w = Rr["w"] - M * (st["rYw"] - st_n["rYw"])
        new["rYw"] = st["rYw"] + r_w / (M + d["w"])
    if has_energy:
        r_E = Rr["E"] - M * (st["rE"] - st_n["rE"])
        new["rE"] = st["rE"] + r_E / M
    return new


def line_update(st, st_n, Rr, d, M, pb, dt, writeback="increment", has_energy=False):
    """前処理を「点対角 (分子・補正) + 乱流の三重対角そのもの」にした更新 (1D の参照用。固定点は jacobi_update と同じ BE 解)。
    乱流は線形で係数が状態に依らないので、分子が 0 の問題では 1 回で質量側の BE 解に着く。"""
    N = pb.N
    rho = st["rho"].astype(np.float64)
    aG = (pb.a * pb.Gam).astype(np.float64)
    T = np.zeros((N, N))
    np.add.at(T, (pb.L, pb.L), aG / rho[pb.L])
    np.add.at(T, (pb.R, pb.R), aG / rho[pb.R])
    np.add.at(T, (pb.L, pb.R), -aG / rho[pb.R])
    np.add.at(T, (pb.R, pb.L), -aG / rho[pb.L])
    turb = np.diag(T).copy()
    Mf = np.broadcast_to(np.asarray(M, dtype=np.float64), (N,))

    def solve(r, dmol):
        A = T + np.diag(Mf + dmol)
        return np.linalg.solve(A, r.astype(np.float64)).astype(dt)
    new = copy_state(st)
    rv = st["rYw"] - st["rg"]
    rv_n = st_n["rYw"] - st_n["rg"]
    r_k = Rr["k"] - M * (st["rY"] - st_n["rY"])
    r_l = Rr["l"] - M * (st["rg"] - st_n["rg"])
    r_Q = Rr["Q"] - M * (st["rQ"] - st_n["rQ"])
    r_v = Rr["v"] - M * (rv - rv_n)
    for k in range(st["rY"].shape[0]):
        new["rY"][k] = st["rY"][k] + solve(r_k[k], d["k"][k].astype(np.float64) - turb)
    dl = solve(r_l, np.zeros(N))
    new["rg"] = st["rg"] + dl
    for k in range(st["rQ"].shape[0]):
        new["rQ"][k] = st["rQ"][k] + solve(r_Q[k], np.zeros(N))
    dv = solve(r_v, d["v"].astype(np.float64) - turb)
    if writeback == "increment":
        new["rYw"] = st["rYw"] + (dv + dl)
    else:
        new["rYw"] = new["rg"] + (rv + dv)
    if has_energy:
        r_E = Rr["E"] - M * (st["rE"] - st_n["rE"])
        new["rE"] = st["rE"] + r_E / M
    return new


def be_residual_norms(st, st_n, Rr, M, has_energy):
    rv = st["rYw"] - st["rg"]
    rv_n = st_n["rYw"] - st_n["rg"]
    r = {"k": Rr["k"] - M * (st["rY"] - st_n["rY"]), "v": Rr["v"] - M * (rv - rv_n),
         "l": Rr["l"] - M * (st["rg"] - st_n["rg"]), "Q": Rr["Q"] - M * (st["rQ"] - st_n["rQ"])}
    if has_energy:
        r["E"] = Rr["E"] - M * (st["rE"] - st_n["rE"])
    return {k: float(np.max(np.abs(v.astype(np.float64)))) for k, v in r.items()}


def be_solve(st_n, pb, M, dt, T0=None, tol=1e-7, maxit=20000, corr="mean", order="vapor",
             vars_mode="vl", writeback="increment", energy="full", stats=None, floor_eps=16.0, precond="point"):
    """解き切った後退 Euler: 点対角を前処理にした固定点反復を BE 残差 ≤ tol × (反復開始時の残差) まで (成分ごと)。
    成分ごとに丸め床 16ε·max|M ρφ| (ε はその精度の machine eps) 以下も収束とみなす (開始時の残差自体が丸め程度の成分
    — 例: 一様な蒸気の ρY_w − ρg — は相対基準が意味を持たないため)。床で止めた回数は stats['floor']。"""
    has_energy = pb.thermo is not None and T0 is not None
    st = copy_state(st_n)
    T = T0
    r0 = None
    for it in range(maxit):
        if has_energy:
            T = pb.thermo.T_from_state(st, T)
        F = fluxes(st, T, pb, dt, corr=corr, order=order, energy=energy)
        Rr = residuals(F, pb, dt)
        nr = be_residual_norms(st, st_n, Rr, M, has_energy)
        if r0 is None:
            r0 = nr
            fl = {k: floor_eps * float(np.finfo(dt).eps) * float(np.max(np.abs((M * v).astype(np.float64))))
                  for k, v in (("k", st["rY"]), ("v", st["rYw"] - st["rg"]), ("l", st["rg"]), ("Q", st["rQ"]), ("E", st["rE"]))}
        rel = max((nr[k] / r0[k]) if r0[k] > 0 else (0.0 if nr[k] == 0 else math.inf) for k in nr)
        if rel <= tol:
            return st, T, it, rel
        if all(nr[k] <= max(tol * r0[k], fl[k]) for k in nr):
            if stats is not None:
                stats["floor"] = stats.get("floor", 0) + 1
            return st, T, it, rel
        d = diagonals(st, F, pb, dt, corr=corr)
        if precond == "line":
            st = line_update(st, st_n, Rr, d, M, pb, dt, writeback, has_energy)
        else:
            st = jacobi_update(st, st_n, Rr, d, M, dt, vars_mode, writeback, has_energy)
    return st, T, maxit, rel


def renormalize(st, dt):
    """forge species_renormalize_d と同じ (ρY ≥ 0 化 → Σ = ρ へ比例) で、係数を ρg と Q にも掛ける (§4.2)。補正量を返す。"""
    s = copy_state(st)
    neg = float(np.sum(np.minimum(s["rY"], 0).astype(np.float64))) + float(np.sum(np.minimum(s["rYw"], 0).astype(np.float64)))
    s["rY"] = np.maximum(s["rY"], 0)
    s["rYw"] = np.maximum(s["rYw"], 0)
    tot = s["rY"].astype(np.float64).sum(axis=0) + s["rYw"].astype(np.float64)
    fac = s["rho"].astype(np.float64) / tot
    s["rY"] = (s["rY"].astype(np.float64) * fac).astype(dt)
    s["rYw"] = (s["rYw"].astype(np.float64) * fac).astype(dt)
    s["rg"] = (s["rg"].astype(np.float64) * fac).astype(dt)
    s["rQ"] = (s["rQ"].astype(np.float64) * fac).astype(dt)
    return s, abs(neg), float(np.max(np.abs(fac - 1.0)))


def realizability_clamp(st, V):
    """0 ≤ ρg ≤ ρY_w、Q ≥ 0 (forge cond_realizability_clamp_d の上下限部分)。補正量 Σ|Δ|V (double) と作動数を返す。"""
    s = copy_state(st)
    r0 = s["rg"].copy()
    s["rg"] = np.minimum(np.maximum(s["rg"], 0), s["rYw"])
    dg = np.abs(s["rg"].astype(np.float64) - r0.astype(np.float64))
    q0 = s["rQ"].copy()
    s["rQ"] = np.maximum(s["rQ"], 0)
    dq = np.abs(s["rQ"].astype(np.float64) - q0.astype(np.float64))
    return s, float(np.sum(dg * V)), int(np.count_nonzero(dg)) + int(np.count_nonzero(dq)), float(np.sum(dq * V))


def totals(st, V):
    V = V.astype(np.float64)
    t = {"rYw": float(np.sum(st["rYw"].astype(np.float64) * V)), "rg": float(np.sum(st["rg"].astype(np.float64) * V)),
         "rv": float(np.sum((st["rYw"].astype(np.float64) - st["rg"].astype(np.float64)) * V)),
         "rE": float(np.sum(st["rE"].astype(np.float64) * V)), "absE": float(np.sum(np.abs(st["rE"].astype(np.float64)) * V))}
    for k in range(st["rY"].shape[0]):
        t["rY%d" % k] = float(np.sum(st["rY"][k].astype(np.float64) * V))
    for k in range(st["rQ"].shape[0]):
        t["rQ%d" % k] = float(np.sum(st["rQ"][k].astype(np.float64) * V))
    return t


# ===================================================================== S1 分子流束の構造 (i)
def test_structure(quick):
    print("\n=== S1 分子流束の構造 (i): 三成分・気相一様・g のみ勾配 ===")
    Dbin = np.array([[0, 1, 2], [1, 0, 3], [2, 3, 0]], float) * 1e-5   # D_12, D_13, D_23
    MW = np.array([0.028, 0.028, 0.028])
    zg = np.array([0.2, 0.3, 0.5])     # [非水 0, 非水 1, 蒸気]
    rho = 1.0
    res = {}
    for dx, gL in [(1.0, 0.1), (1.0e-3, 0.1)]:
        gR = gL + 0.1 * dx    # ∇g = 0.1 m⁻¹
        for dtn, dt in [("float64", np.float64), ("float32", np.float32)]:
            g = np.array([gL, gR])
            rY = np.vstack([(1 - g) * zg[0], (1 - g) * zg[1]]) * rho
            rv = (1 - g) * zg[2] * rho
            st = make_state([rho, rho], rY, rv + g * rho, g * rho, np.zeros((1, 2)), [0, 0], dt)
            pb = Problem(2, [1.0 / dx], [0.0], Dbin, MW, [dx, dx])
            for design in ("new", "old"):
                F = fluxes(st, None, pb, dt, design=design)
                X = (zg / MW) / np.sum(zg / MW)
                Dk = mixavg_D(X[:, None], Dbin)[:, 0]
                rhog = rho * (1 - 0.5 * (gL + gR))
                lim = 1e-6 * rhog * Dk.max() * 0.1
                jmax = float(np.max(np.abs(F["j"].astype(np.float64))))
                res[(dx, dtn, design)] = (jmax, lim, F["j"].astype(np.float64)[:, 0] / rho)
    # 判定は plan の幾何 (dx = 1 m, Δg = 0.1 を 1 面で) で
    for dtn in ("float64", "float32"):
        jmax, lim, j = res[(1.0, dtn, "new")]
        verdict(jmax <= lim, f"S1(i) §4.2 気相基準 [{dtn}, dx=1 m]: max|j_k| = {jmax:.3e} ≤ 1e-6·ρ_g·maxD·|∇g| = {lim:.3e}")
        jo, _, jv = res[(1.0, dtn, "old")]
        info(f"S1(i) 旧案 [{dtn}, dx=1 m]: j/ρ = [{jv[0]:.4e}, {jv[1]:.4e}, {jv[2]:.4e}] m/s (codex 2026-09-27: [-1.3182e-7, -6.1364e-8, 1.9318e-7]); "
             f"許容の {jo/lim:.2e} 倍 → 判別 {'可' if jo > lim else '不可'}")
        verdict(jo > 100 * lim, f"S1(i) 判別: 旧案は同じ判定で FAIL すること [{dtn}] (max|j| = {jo:.3e} > 許容)")
    for dtn in ("float64", "float32"):
        jmax, lim, _ = res[(1.0e-3, dtn, "new")]
        info(f"S1(i) 尺度依存 [{dtn}, dx=1e-3 m, Δg=1e-4/面]: max|j| = {jmax:.3e}、許容 {lim:.3e} (比 {jmax/lim:.2e})"
             + (" — float32 の z 丸め ~ε₃₂ z/Δg が効く。許容は ∇g 基準なので面あたり Δg が小さいと float32 で満たせない" if dtn == "float32" else ""))


# ===================================================================== S2 面恒等式 (ii)
def test_face_identities(quick):
    print("\n=== S2 面恒等式 (ii): float32 で計算、格納値から float64 で評価 ===")
    rng = np.random.default_rng(20261002)
    nF = 20000 if quick else 200000
    worst = {}
    for K in (1, 2):
        nG = K + 1
        Dbin = np.array([[0, 1, 2], [1, 0, 3], [2, 3, 0]], float) * 1e-5 if nG == 3 else np.array([[0, 2e-5], [2e-5, 0]])
        MW = np.array([0.028, 0.032, 0.018]) if nG == 3 else np.array([0.028, 0.018])
        # 面ごとに独立な 2 セル (2i, 2i+1)。3 割は両側がほぼ同じ状態 (勾配が小さい面)
        rho = rng.uniform(0.05, 2.0, size=(nF, 2))
        g = rng.uniform(0.0, 0.3, size=(nF, 2)) * (rng.random((nF, 2)) < 0.8)
        vfrac = rng.choice([0.0, 1e-8, 1e-4, 1e-2], size=(nF, 2)) * rng.random((nF, 2))
        zr = np.stack([rng.dirichlet(np.ones(K), size=nF) for _ in range(2)], axis=1) if K > 1 else np.ones((nF, 2, 1))
        same = rng.random(nF) < 0.3
        pert = 1 + 1e-5 * rng.standard_normal((same.sum(),))
        g[same, 1] = g[same, 0] * pert
        vfrac[same, 1] = vfrac[same, 0] * pert
        rho[same, 1] = rho[same, 0]
        zr[same, 1] = zr[same, 0]
        rhog = rho * (1 - g)
        rv = rhog * vfrac
        rY = (zr * (rhog - rv)[..., None]).reshape(2 * nF, K).T
        a = rng.uniform(0.1, 10, nF)
        Gam = rng.uniform(0, 1e-3, nF) * (rng.random(nF) < 0.9)
        pb = Problem(2 * nF, np.zeros(2 * nF - 1), np.zeros(2 * nF - 1), Dbin, MW, np.ones(2 * nF))
        pb.L = np.arange(0, 2 * nF, 2)
        pb.R = pb.L + 1
        pb.a, pb.Gam = a, Gam
        st = make_state(rho.reshape(-1), rY, (rv + g * rho).reshape(-1), (g * rho).reshape(-1),
                        np.zeros((1, 2 * nF)), np.zeros(2 * nF), np.float32)
        for corr in ("mean", "upwind"):
            for order in ("vapor", "total"):
                F = fluxes(st, None, pb, np.float32, corr=corr, order=order)
                j = F["j"].astype(np.float64)
                sabs = np.sum(np.abs(F["j0"].astype(np.float64)), axis=0)
                m = sabs > 0
                ws = float(np.max(np.abs(j.sum(axis=0))[m] / (8 * EPS32 * sabs[m]))) if m.any() else 0.0
                Jw, Jv, Jl = (F[k].astype(np.float64) for k in ("Jw", "Jv", "Jl"))
                den = np.abs(Jw) + np.abs(Jv) + np.abs(Jl)
                m = den > 0
                wi = float(np.max(np.abs(Jw - Jv - Jl)[m] / (8 * EPS32 * den[m]))) if m.any() else 0.0
                worst[(K, corr, order)] = (ws, wi)
    for (K, corr, order), (ws, wi) in worst.items():
        tag = f"気相 {K+1} 種, 補正の面 z = {corr}, 順序 = {order}"
        default = (corr == "mean" and order == "vapor")
        msg1 = f"S2(ii) |Σ j_k| / (8ε₃₂Σ|j⁰|) 最大 {ws:.3f} [{tag}, {nF} 面]"
        msg2 = f"S2(ii) |J_w−J_v−J_l| / (8ε₃₂Σ|J|) 最大 {wi:.3f} [{tag}, {nF} 面]"
        if default:
            verdict(ws <= 1.0, msg1)
            verdict(wi <= 1.0, msg2)
        else:
            info(msg1 + (" (≤1)" if ws <= 1 else " (>1)"))
            info(msg2 + (" (≤1)" if wi <= 1 else " (>1)"))


# ===================================================================== S3 3 セル判別
def three_cell_state(dt):
    g = np.array([0.1, 0.2, 0.2])
    Yw = g + 0.1
    rQ = np.vstack([g * 3.0, g * 2.0, g * 1.0])   # Q も質量当たり g に比例 (値は任意)
    return make_state(np.ones(3), (1 - Yw)[None, :], Yw, g, rQ, np.zeros(3), dt)


def test_three_cell(quick):
    print("\n=== S3 3 セル判別 (A = 点対角 1 回更新 / B = 解き切った後退 Euler) ===")
    out = {}
    for dtn, dt, fe in [("float64", np.float64, 16.0), ("float32", np.float32, 16.0), ("float32", np.float32, 4.0)]:
        pb = Problem(3, [1, 1], [1, 1], np.array([[0, 1e-5], [1e-5, 0]]), np.array([0.028, 0.018]), np.ones(3), molc=0.0)
        M = dt(1.0)
        st0 = three_cell_state(dt)
        t0 = totals(st0, pb.V)
        for scheme in (("A", "B") if fe == 16.0 else ("B",)):
            st = copy_state(st0)
            rec = {}
            nit_max, stats = 0, {}
            for n in range(1, 1001):
                if scheme == "A":
                    F = fluxes(st, None, pb, dt)
                    Rr = residuals(F, pb, dt)
                    d = diagonals(st, F, pb, dt)
                    st = jacobi_update(st, st, Rr, d, M, dt)
                else:
                    st, _, nit, rel = be_solve(st, pb, M, dt, stats=stats, floor_eps=fe, maxit=5000)
                    nit_max = max(nit_max, nit)
                if n in (1, 1000):
                    t = totals(st, pb.V)
                    rels = {k: (t[k] - t0[k]) / abs(t0[k]) for k in t0 if k not in ("rE", "absE") and t0[k] != 0}
                    rv = (st["rYw"].astype(np.float64) - st["rg"].astype(np.float64))
                    sumY = st["rY"].astype(np.float64).sum(axis=0) + st["rYw"].astype(np.float64)
                    rec[n] = (rels, float(rv.min()), float(np.max(np.abs(sumY - st["rho"].astype(np.float64)) / st["rho"].astype(np.float64))))
            out[(dtn, scheme, fe)] = (rec, nit_max, stats.get("floor", 0))
    for (dtn, scheme, fe), (rec, nit, nfl) in out.items():
        for n in (1, 1000):
            rels, vmin, sY = rec[n]
            worst = max(abs(v) for v in rels.values())
            line = (f"S3 {scheme} [{dtn}] {n} 更新: 総液量 {rels['rg']:+.7e}、総水分 {rels['rYw']:+.3e}、非水 {rels['rY0']:+.3e}、"
                    f"Q0 {rels['rQ0']:+.3e} (最大 |相対| {worst:.3e})、min ρv {vmin:.3e}、max|ΣρY−ρ|/ρ {sY:.2e}")
            if scheme == "B":
                line += f" (反復 ≤{nit}; BE 残差 ≤1e-7 相対か丸め床 {fe:g}ε 以下で停止、床で停止 {nfl} 回)"
            if scheme == "B" and dtn == "float64":
                verdict(worst <= 1e-6 and vmin >= 0 and sY <= 8 * EPS32, line)
            elif scheme == "A" and dtn == "float64":
                print("[A]    " + line)
            else:
                info(line)
    out = {(k[0], k[1]): (v[0],) for k, v in out.items() if k[2] == 16.0}
    a_fail = max(abs(v) for v in out[("float64", "A")][0][1][0].values()) > 1e-6
    info("S3 float32 の B は停止基準 (丸め床) に依存する: 残した BE 残差の和がそのまま総量の誤差になる (保存は固定点でだけ成り立つ)")
    b_ok = all(max(abs(v) for v in out[("float64", "B")][0][n][0].values()) <= 1e-6 for n in (1, 1000))
    if a_fail and b_ok:
        verdict(True, "S3 判別規則: A が FAIL・B が PASS → 点対角更新が保存破れの原因、B の更新契約 (保存形 BE の固定点) を支持")
    elif not b_ok:
        verdict(False, "S3 判別規則: B も FAIL → #4 不合格 (面交換・反復残差・commit の整合を調べる)")
    else:
        verdict(False, "S3 判別規則: A も合格 → 実装が縮約式と異なる (保存を成立させている追加処理を特定する)")
    a1 = out[("float64", "A")][0][1][0]["rg"]
    a1000 = out[("float64", "A")][0][1000][0]["rg"]
    info(f"S3 A の総液量: 1 更新 {a1*100:+.7f} %、1000 更新 {a1000*100:+.7f} % (上位検算 +3.3333346 % / +2.8571441 %)")


# ===================================================================== S4 非負
def test_nonnegativity(quick):
    print("\n=== S4 非負 ===")
    # 2 セル例: 蒸気 0、Y_w = g = 0.1/0.2、分子・乱流の面係数 1、V/Δτ = 1、ρ = 1
    for dtn, dt in [("float64", np.float64), ("float32", np.float32)]:
        g = np.array([0.1, 0.2])
        st0 = make_state([1, 1], (1 - g)[None, :], g, g, (g * 2)[None, :], [0, 0], dt)
        pb = Problem(2, [1.0], [1.0], np.array([[0, 1.0], [1.0, 0]]), np.array([0.028, 0.018]), np.ones(2), molc=1.0)
        M = dt(1.0)
        stB, _, nit, rel = be_solve(st0, pb, M, dt)
        rvB = stB["rYw"].astype(np.float64) - stB["rg"].astype(np.float64)
        _, corr, nact, _ = realizability_clamp(stB, pb.V)
        line = f"S4 2 セル [{dtn}] 解き切った擬似時間更新 (新, 蒸気/液変数): ρv = [{rvB[0]:.3e}, {rvB[1]:.3e}]、クランプ補正 {corr:.3e} ({nact} 件; 反復 {nit})"
        if dtn == "float64":
            verdict(rvB.min() >= 0 and corr == 0, line)
        else:
            info(line)
        # 判別: 改訂前 (総水分と液をそれぞれの点対角で 1 回)
        F = fluxes(st0, None, pb, dt, design="old")
        Rr = residuals(F, pb, dt)
        d = diagonals(st0, F, pb, dt, design="old")
        stO = jacobi_update(st0, st0, Rr, d, M, dt, vars_mode="wl")
        rvO = stO["rYw"].astype(np.float64) - stO["rg"].astype(np.float64)
        info(f"S4 2 セル [{dtn}] 改訂前 (総水分/液の点対角 1 回): Y_w' = [{float(stO['rYw'][0]):.8f}, {float(stO['rYw'][1]):.8f}]、"
             f"g' = [{float(stO['rg'][0]):.8f}, {float(stO['rg'][1]):.8f}]、ρv = [{rvO[0]:.4e}, {rvO[1]:.4e}] (codex: 0.13333334 / 0.15 / −0.016666666)")
    # 補正の面 z (plan §4.2 は面値を定めていない) の非負性: 拡散係数の差が大きい三成分で、蒸気 0 のセルから補正流束が蒸気を
    # 持ち出す例 (情報; mean は −z_v,f S の隣セル分が非対角に負で入るので M 行列でなくなる)。
    #   L: (A 0, B 1, 蒸気 0)、R: (A 0.99, B 0, 蒸気 0.01)、g = 0、二元 D_AB = D_Av = 10、D_Bv = 0.1 (相対値)、面係数・V/Δτ = 1
    for corr in ("mean", "upwind"):
        st0 = make_state([1, 1], np.array([[0.0, 0.99], [1.0, 0.0]]), [0.0, 0.01], [0.0, 0.0], np.zeros((1, 2)), [0, 0], np.float64)
        pb = Problem(2, [1.0], [0.0], np.array([[0, 10, 10], [10, 0, 0.1], [10, 0.1, 0]]), np.array([0.028, 0.028, 0.028]), np.ones(2))
        stB, _, nitB, _ = be_solve(st0, pb, 1.0, np.float64, corr=corr, maxit=5000)
        F = fluxes(st0, None, pb, np.float64, corr=corr)
        Rr = residuals(F, pb, np.float64)
        d = diagonals(st0, F, pb, np.float64, corr=corr)
        stJ = jacobi_update(st0, st0, Rr, d, 1.0, np.float64)
        vB = float(stB["rYw"][0] - stB["rg"][0])
        vJ = float(stJ["rYw"][0] - stJ["rg"][0])
        info(f"S4 補正の面 z = {corr} の反例 (D 比 100, 組成の段差 1): 蒸気 0 のセル L の ρv — 点対角 1 回 {vJ:+.3e}、解き切った更新 {vB:+.3e} "
             f"(反復 {nitB}{' = 上限: 固定点反復が収束しない' if nitB >= 5000 else ''}; Σj⁰ = {float(F['S'][0]):+.3f}, j_v⁰ = {float(F['j0'][2][0]):+.4f}) "
             f"→ {'負になる' if min(vB, vJ) < 0 else '非負'}")
    # ランダム 1D 100 セル × 10⁴ 擬似反復
    nit = 2000 if quick else 10000
    rng0 = np.random.default_rng(7)
    N = 100
    base = {}
    for K in (1, 2):
        rng = np.random.default_rng(100 + K)
        rho = rng.uniform(0.1, 1.0, N)
        g = rng.uniform(0.0, 0.3, N)
        cat = rng.integers(0, 3, N)
        vfrac = np.where(cat == 0, 0.0, np.where(cat == 1, rng.uniform(0, 1e-7, N), rng.uniform(0, 1e-3, N)))
        rhog = rho * (1 - g)
        rv = rhog * vfrac
        zrest = rng.dirichlet(np.ones(K), size=N).T if K > 1 else np.ones((1, N))
        rY = zrest * (rhog - rv)
        rQ = np.vstack([g * rho * 1e3, g * rho * 1e6, g * rho * 1e9]) * rng.uniform(0.5, 2, (3, N))
        a = rng.uniform(0.0, 2.0, N - 1)
        Gam = rng.uniform(0.0, 1.0, N - 1)
        nG = K + 1
        Dbin = np.array([[0, 1, 2], [1, 0, 3], [2, 3, 0]], float)[:nG, :nG] if nG == 3 else np.array([[0, 1.0], [1.0, 0]])
        Dbin = Dbin * rng.uniform(0.5, 1.5)
        MW = np.array([0.028, 0.032, 0.018]) if nG == 3 else np.array([0.028, 0.018])
        M = 10 ** rng.uniform(-1, 1, N)
        base[K] = (rho, rY, rv + g * rho, g * rho, rQ, a, Gam, Dbin, MW, M)
    del rng0
    variants = [(K, corr, dtn, wb) for K in (1, 2) for corr in ("mean", "upwind") for dtn in ("float64", "float32")
                for wb in ("increment", "rebuild")]
    for (K, corr, dtn, wb) in variants:
        dt = np.float64 if dtn == "float64" else np.float32
        rho, rY, rYw, rg, rQ, a, Gam, Dbin, MW, M = base[K]
        pb = Problem(N, a, Gam, Dbin, MW, np.ones(N), molc=1.0)
        st = make_state(rho, rY, rYw, rg, rQ, np.zeros(N), dt)
        Md = M.astype(dt)
        tot_corr, n_act, vmin, tot_rn = 0.0, 0, math.inf, 0.0
        for it in range(nit):
            F = fluxes(st, None, pb, dt, corr=corr)
            Rr = residuals(F, pb, dt)
            d = diagonals(st, F, pb, dt, corr=corr)
            st = jacobi_update(st, st, Rr, d, Md, dt, writeback=wb)
            rvv = st["rYw"].astype(np.float64) - st["rg"].astype(np.float64)
            vmin = min(vmin, float(rvv.min()))
            st, rn_neg, _ = renormalize(st, dt)
            tot_rn += rn_neg
            st, c, na, cq = realizability_clamp(st, pb.V)
            tot_corr += c + cq
            n_act += na
        line = (f"S4 ランダム 100 セル × {nit} 擬似反復 [気相 {K+1} 種, 補正 z = {corr}, {dtn}, ρY_w の戻し = {wb}]: "
                f"クランプ補正 Σ|Δ|V = {tot_corr:.3e} ({n_act} 件)、ρY<0 の再正規化補正 {tot_rn:.3e}、min ρv (クランプ前) {vmin:.3e}")
        if dtn == "float64" and wb == "increment" and corr == "mean":
            verdict(tot_corr == 0 and n_act == 0 and tot_rn == 0, line)
        else:
            info(line)


# ===================================================================== S5 保存
def test_conservation(quick):
    print("\n=== S5 保存: 閉じた箱・拡散単独・解き切った物理時間更新 ===")
    nstep = 200 if quick else 1000
    N = 24
    rng = np.random.default_rng(31)
    th = Thermo([nasa_species("N2"), cpg_species("AR", 0.039948, 520.3, 0.0), cpg_species("HE", 0.0040026, 5193.0, 0.0)],
                nasa_species("H2O"), None)
    th.lat = Latent(th.v)
    rho = rng.uniform(0.2, 0.6, N)
    g = rng.uniform(0.0, 0.02, N)
    vfrac = rng.uniform(0.0, 0.01, N) * (rng.random(N) < 0.7)
    rhog = rho * (1 - g)
    rv = rhog * vfrac
    zr = rng.dirichlet([3.0, 1.0], size=N).T
    rY = np.vstack([zr[0] * (rhog - rv), zr[1] * (rhog - rv), np.zeros(N)])   # HE は初期総量 0
    rQ = np.vstack([g * rho * 1e3, g * rho * 1e6, np.zeros(N)])              # ρQ0 も初期総量 0 の量として
    T0 = rng.uniform(205.0, 260.0, N)
    Yk = rY / rho
    rE = rho * th.e_mix(T0, Yk, (rv + g * rho) / rho, g)
    # D は 1e-4 m²/s 級、面 a = 1/dx (dx = 1 mm)、Γ = 1e-4 kg/m/s 級; Δt は拡散数 ~1
    Dbin = np.array([[0, 2.0, 6.0, 2.2], [2.0, 0, 5.5, 1.8], [6.0, 5.5, 0, 7.0], [2.2, 1.8, 7.0, 0]]) * 1e-5
    MW = np.array([0.028, 0.0399, 0.0040, 0.018])
    dx = 1e-3
    a = np.full(N - 1, 1.0 / dx)
    Gam = rng.uniform(0.0, 2e-4, N - 1)
    V = np.full(N, dx)
    dtime = 2e-3   # 乱流の拡散数 Γ/ρ Δt/dx² ≈ 1、分子 ≈ 0.2〜1.4
    pb = Problem(N, a, Gam, Dbin, MW, V, molc=1.0, thermo=th)
    for dtn, dt, fe in [("float64", np.float64, 16.0), ("float32", np.float32, 16.0), ("float32", np.float32, 4.0)]:
        st = make_state(rho, rY, rv + g * rho, g * rho, rQ, rE, dt)
        t0 = totals(st, V)
        M = (V / dtime).astype(dt)
        T = T0.copy()
        nit_max, rel_max, vmin, sYmax = 0, 0.0, math.inf, 0.0
        for n in range(nstep):
            st, T, nit, rel = be_solve(st, pb, M, dt, T0=T, maxit=4000, precond="line", floor_eps=fe)
            nit_max, rel_max = max(nit_max, nit), max(rel_max, rel)
            rvv = st["rYw"].astype(np.float64) - st["rg"].astype(np.float64)
            vmin = min(vmin, float(rvv.min()))
            sY = st["rY"].astype(np.float64).sum(axis=0) + st["rYw"].astype(np.float64)
            sYmax = max(sYmax, float(np.max(np.abs(sY - st["rho"].astype(np.float64)) / st["rho"].astype(np.float64))))
        t = totals(st, V)
        rel_nz = {k: (t[k] - t0[k]) / abs(t0[k]) for k in t0 if k not in ("rE", "absE", "rv") and t0[k] != 0}
        abs_z = {k: abs(t[k] - t0[k]) for k in t0 if k not in ("rE", "absE", "rv") and t0[k] == 0}
        eE = abs(t["rE"] - t0["rE"]) / t0["absE"]
        worst = max(abs(v) for v in rel_nz.values())
        zworst = max(abs_z.values()) if abs_z else 0.0
        line = (f"S5 [{dtn}, 丸め床 {fe:g}ε] {nstep} 更新 (反復 ≤{nit_max}): 最大 |相対| {worst:.3e} "
                f"({', '.join(f'{k} {v:+.1e}' for k, v in rel_nz.items())})、ρE {eE:.3e}、初期 0 の量 {zworst:.1e}、"
                f"min ρv {vmin:.2e}、max|ΣρY−ρ|/ρ {sYmax:.1e}")
        if dtn == "float64":
            verdict(worst <= 1e-6 and eE <= 1e-6 and zworst <= 1e-12 and vmin >= 0, line)
        else:
            info(line)
        T_end = T
    info(f"S5 温度範囲 初期 {T0.min():.1f}–{T0.max():.1f} K → 最終 {T_end.min():.2f}–{T_end.max():.2f} K (熱伝導なし: 種のエンタルピー流束だけで動く)")


# ===================================================================== S6 エネルギー
def test_energy(quick):
    print("\n=== S6 エネルギー: 乱流のみ・総組成一様・g のみ勾配 (分子 0) ===")
    th = Thermo([nasa_species("N2")], nasa_species("H2O"), None)
    th.lat = Latent(th.v)
    N = 200
    Lx = 0.1
    dx = Lx / N
    x = (np.arange(N) + 0.5) * dx
    rho0, Yw0, T0 = 0.2, 0.011, 210.0
    g = 0.006 + 0.004 * np.cos(math.pi * x / Lx)
    Gam = 1e-3
    pb = Problem(N, np.full(N - 1, 1 / dx), np.full(N - 1, Gam), np.array([[0, 1e-5], [1e-5, 0]]),
                 np.array([0.028, 0.018]), np.full(N, dx), molc=0.0, thermo=th)
    rho = np.full(N, rho0)
    Tn = np.full(N, T0)
    Yk = np.array([[1 - Yw0] * N])
    rE = rho * th.e_mix(Tn, Yk, np.full(N, Yw0), g)
    k = math.pi / Lx
    dtime = 0.25 / ((Gam / rho0) * k * k)
    st0 = make_state(rho, (1 - Yw0) * rho[None, :], Yw0 * rho, g * rho, (g * rho)[None, :], rE, np.float64)
    M = pb.V / dtime
    Rw = th.v.R
    cv = th.cv_eff(Tn, Yk, np.full(N, Yw0), g)
    Lc = th.lat.L(Tn)
    results = {}
    for variant in ("full", "no_latent", "latent_only"):
        st1, T1, nit, rel = be_solve(st0, pb, M, np.float64, T0=Tn, energy=variant, precond="line", maxit=200)
        drg = st1["rg"] - st0["rg"]
        drE = st1["rE"] - st0["rE"]
        if variant == "latent_only":
            # 液は動かないので、比較する r_g は「動くはずだった」液の変化 (full の値) を使う
            drg_ref = results["full"][3]
        else:
            drg_ref = drg
        dT_ana = -Rw * Tn * drg_ref / (rho * cv)               # = −R_wT r_g Δt/(ρ c_v,eff)
        # EOS 接線投影: 状態 n のまわりの線形化。エネルギー流束の h・L は Tⁿ の面値で評価し (流束の質量側は解いた n+1)、
        # δT = (δ(ρe) − (R_wT − L) δ(ρg)) / (ρ c_v,eff) (ρ, ρY は不変)。L(T_f^{n+1}) で評価すると J_l·∇L (2 次) が入り、
        # cos の節の近く (r_g → 0) で相対差が大きく見える — 物理の項で誤りではないので、有限更新の比較 (+1e-3 K) 側で見る。
        Ft = fluxes(st1, Tn, pb, np.float64, energy=variant)
        drE_tan = residuals(Ft, pb, np.float64)["E"] / M
        dT_proj = (drE_tan - (Rw * Tn - Lc) * drg) / (rho * cv)
        dT_fin = T1 - Tn
        sel = np.abs(dT_ana) > 0.05 * np.max(np.abs(dT_ana))
        e_proj = float(np.max(np.abs(dT_proj[sel] - dT_ana[sel]) / np.abs(dT_ana[sel])))
        e_projT1 = float(np.max(np.abs(((drE - (Rw * Tn - Lc) * drg) / (rho * cv))[sel] - dT_ana[sel]) / np.abs(dT_ana[sel])))
        e_fin = float(np.max(np.abs(dT_fin[sel] - dT_ana[sel]) - (0.01 * np.abs(dT_ana[sel]) + 1e-3)))
        results[variant] = (e_proj, e_fin, float(np.max(np.abs(dT_ana))), drg, nit, e_projT1)
    e_proj, e_fin, dTmax, drg, nit, e_projT1 = results["full"]
    verdict(e_proj <= 0.01, f"S6 §4.2 (−L J_l あり): 接線投影と解析値の最大相対差 {e_proj:.2e} ≤ 1 % (|ΔT| 最大 {dTmax:.3f} K, Δg 最大 {np.max(np.abs(drg))/rho0:.2e}, 反復 {nit})")
    info(f"S6 参考: 面の h・L を T^(n+1) で評価した δ(ρe) で投影すると最大相対差 {e_projT1:.2e} (J_l·∇L の 2 次項; |r_g| が最大の 5 % 以上のセル)")
    verdict(e_fin <= 0, f"S6 §4.2 有限更新 (EOS Newton 反転): |ΔT − 解析| ≤ 1 %·|解析| + 1e-3 K (超過量の最大 {e_fin:+.2e} K)")
    for variant, label in (("no_latent", "潜熱流束なし (h_v J_w だけ)"), ("latent_only", "エネルギーだけ潜熱流束・液は拡散しない (禁止案)")):
        ep = results[variant][0]
        verdict(ep > 0.01, f"S6 判別: {label} は FAIL すること (接線投影と解析値の相対差 {ep:.2e})")
    # 連続解との照合 (情報): r_g の解析値 Γ ∂²g/∂x² を BE の 1 更新と比べる
    rg_cont = -Gam * k * k * 0.004 * np.cos(math.pi * x / Lx) * dtime
    info(f"S6 参考: BE 1 更新の Δ(ρg) / 連続の r_gΔt (中央 1/2 区間) = {np.median(drg[N//4:3*N//4] / rg_cont[N//4:3*N//4]):.4f} "
         f"(BE の減衰 1/(1+Dk²Δt) = {1/(1+0.25):.4f})")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quick", action="store_true", help="反復数を減らす (開発用; 判定は正式でない)")
    ap.add_argument("--only", default="", help="S1..S6 をカンマ区切りで")
    args = ap.parse_args()
    tests = {"S1": test_structure, "S2": test_face_identities, "S3": test_three_cell,
             "S4": test_nonnegativity, "S5": test_conservation, "S6": test_energy}
    sel = [s.strip() for s in args.only.split(",") if s.strip()] or list(tests)
    for name in sel:
        tests[name](args.quick)
    print(f"\n{'ALL PASS' if g_fail == 0 else f'{g_fail} FAIL'}" + (" (--quick: 正式判定ではない)" if args.quick else ""))
    return 0 if g_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
