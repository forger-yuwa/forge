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
  S1 分子流束の構造 (i)  三成分 (等分子量, z = [0.2,0.3,0.5]) の気相一様・∇g = 0.1。S1a = 作用素 (D_k = [1,2,3]e-5 を固定、plan §6 の字義)、
                          S1b = 係数 (二元 D_12/D_13/D_23 = [1,2,3]e-5 の混合平均; 手計算と照合) とその係数での作用素。各々 3 つの許容に分ける
                          (§5.1 #4a (3)): 厳密入力 float64 は 1e-6·ρ_g·max D·|∇g|、入力量子化 (格納 float32 を float64 で評価) は ε₃₂ の誤差伝播
                          の許容 + 物理の許容、演算誤差 (float32 − 同じ格納入力の float64) は演算の誤差伝播の許容 (s1_tolerances)。旧案は FAIL すること。
  S2 面恒等式 (ii)       ランダム面で |Σ_気相 j_k| ≤ 8ε₃₂ Σ|j_k⁰|、|J_w − J_v − J_l| ≤ 8ε₃₂(|J_w|+|J_v|+|J_l|) (float32 で計算、float64 で評価)。
  S3 3 セル判別         閉じた直列 3 セル、ρ=V=Δt=1、乱流係数 1、分子・移流・相変化なし、g=[0.1,0.2,0.2]、Y_w=g+0.1。
                          A = 点対角 1 回更新、B = 解き切った後退 Euler (BE 残差 ≤1e-7 相対)。B: 総液量・各種総量の相対変化 ≤1e-6
                          (1・1000 更新、再正規化・クランプ前、float64 集計)、全セル ρv ≥ 0、|ΣρY − ρ| ≤ 8ε₃₂ρ。規則: A FAIL・B PASS で契約を支持。
                          float32 の B は停止則ごとに累積 ≤1e-6・**反復上限到達なし** で判定 (§5.1 #4a (2))。
  S4 非負               2 セル例 (蒸気 0、Y_w = g = 0.1/0.2、面係数・V/Δτ = 1) の解き切った更新で ρv ≥ 0 (クランプ 0)。
                          蒸気 0 近傍のランダム 1D 100 セル × 10⁴ 擬似反復でクランプ補正量 0。補正の面 z の反例 (両セル g = 0, D 比 100) で
                          風上 z の非負と反復の収束を記録 (§5.1 #4a (5))。
  S5 保存               閉じた箱・拡散単独の解き切った物理時間更新 1000 回で Σ ρY_k V・Σ ρg V・Σ ρQ V 相対 ≤1e-6、Σ ρE V は初期 Σ|ρE|V
                          分母で ≤1e-6、初期総量 0 の量は絶対 ≤1e-12。float32 は停止則ごと (上限到達は FAIL)。
  S6 エネルギー         分子拡散 0・総組成一様・乱流のみ・g のみ勾配の 1D で、1 更新の温度変化の EOS 接線投影が −R_wT r_g Δt/(ρ c_v,eff) と ≤1 %
                          (有限更新の温度比較は +1e-3 K)。判別: 潜熱流束なし / エネルギーだけ潜熱 (液を拡散しない) は FAIL すること。
  S7 固定点判別         ρ=V=1、ρY_w=0.3、ρg=0.1、R_wᵗʳ=0、R_gᵗʳ=−0.02、S=+0.02、P_v=3、P_g=2、Q 残差 0 で 1 擬似更新 (float32/float64)。
                          A = P1 (輸送とソースを別の前処理で分割)、B = 非分割の全残差 R_v = R_w − R_g。同じ limiter・commit (vl_limit_commit)。
                          合格: A が ΔρY_w ≈ −0.003333 動き、B は全増分・補正が厳密に 0 (§5.1 #4a (1))。
  S8 A/B (§5.1 #4b)       S5 の更新経路だけを替える: A = line_update / B = 全残差変換 (r_v = r_w − r_l) + vl_limit_commit (B_LIMITS 固定)。
                          float32・upwind・6ε・持ち越しなし・同じ前処理 (line)・同じ初期値・1000 更新・上限 3000。A・B は同じ停止成分 (k, v, w, l, Q, E)。
                          停止判定の残差を float32 で見る版と float64 で再評価する版の両方。見る量: 再正規化前の全量誤差・原方程式残差 (float64)・
                          上限・最小蒸気・θ・補正量。判定規則: A 合格・B 不合格 → 写像、両者合格 → 除外、両者不合格 → 停止則。
  S9 相変化ソースの反復 (§5.1 #4b) B に GrowthSource (簡易成長則) と移流 ṁ (ρ の更新 = 基点 φⁿρ_new) とソース Jacobian (蒸気・液 2×2) を入れて
                          収束・保存 (Σρ, ΣρY_w, ΣρY_k, ΣρE)・非負を見る。
  S10 θ=0 の停止           蒸気 0 で凝縮を要求する状態: 状態は動かないが原方程式残差は 0 でなく、停止判定は上限到達 (FAIL) にする。
  S11 反例を B で           両セル g=0・D 比 100 を新経路 B で、緩和 1 / 0.5 (緩和は前処理の後・limiter の前)。
  補正の面 z の既定は upwind (codex diagnose 2026-10-02 ① で採用)。mean は現行 species_diffusion_d と同じ比較用。
  be_solve が反復上限で返った更新は「解き切っていない」= FAIL (codex diagnose 2026-10-02 ④)。
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

    def __init__(self, N, a, Gam, Dbin, MW, V, molc=1.0, thermo=None, Dfixed=None):
        self.Dfixed = None if Dfixed is None else np.asarray(Dfixed, float)   # 種ごとの D_k を直接与える (混合平均しない作用素試験)
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


def fluxes(st, T, pb, dt, corr="upwind", order="vapor", design="new", energy="full"):
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
    if pb.Dfixed is not None:
        Df = np.broadcast_to(pb.Dfixed.astype(dt)[:, None], Xf.shape).copy()
    else:
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
        # 丸めの尺度 (停止の床用): 生成エンタルピー込みの h は大きいので、q の打ち消し前の大きさ
        qa = np.abs(th.v.h(Tf)) * (np.abs(Jw) + np.abs(Jv) + np.abs(Jl)).astype(np.float64) + np.abs(th.lat.L(Tf)) * np.abs(Jl).astype(np.float64)
        for k, sp in enumerate(th.nw):
            qa = qa + np.abs(sp.h(Tf)) * np.abs(Jk[k]).astype(np.float64)
        # 流束の被演算子の丸め (例: 一様な Y_w の J_w = J_v + J_l は丸め程度だが、h_v ~ 1e7 J/kg を掛けると効く)
        Ga = (G * a).astype(np.float64)
        ops = np.abs(Yw[L]) + np.abs(Yw[R]) + np.abs(g[L]) + np.abs(g[R])
        qa = qa + np.abs(th.v.h(Tf)) * Ga * ops.astype(np.float64)
        for k, sp in enumerate(th.nw):
            qa = qa + np.abs(sp.h(Tf)) * Ga * (np.abs(Yk[k][L]) + np.abs(Yk[k][R])).astype(np.float64)
        out["q_abs"] = qa
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


def diagonals(st, F, pb, dt, corr="upwind", design="new"):
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
    if "w" in Rr:
        r["w"] = Rr["w"] - M * (st["rYw"] - st_n["rYw"])
    if has_energy:
        r["E"] = Rr["E"] - M * (st["rE"] - st_n["rE"])
    return {k: float(np.max(np.abs(v.astype(np.float64)))) for k, v in r.items()}


def be_residual_sums(st, st_n, pb, T, M, corr, order, energy, has_energy):
    """成分ごと (各種・蒸気・液・各 Q・E) の (|Σ_i r_i|, Σ_i |M ρφ_i|, sqrt(Σ_i (M ρφ_i)²))。格納値を float64 にして評価する。"""
    s64, n64 = cast_state(st, np.float64), cast_state(st_n, np.float64)
    M64 = np.asarray(M, np.float64)
    R = residuals(fluxes(s64, T, pb, np.float64, corr=corr, order=order, energy=energy), pb, np.float64)
    rv, rvn = s64["rYw"] - s64["rg"], n64["rYw"] - n64["rg"]
    out = []

    def add(r, x):
        out.append((abs(float(np.sum(r))), float(np.sum(np.abs(M64 * x))), float(np.sqrt(np.sum((M64 * x) ** 2)))))
    for k in range(s64["rY"].shape[0]):
        add(R["k"][k] - M64 * (s64["rY"][k] - n64["rY"][k]), s64["rY"][k])
    add(R["v"] - M64 * (rv - rvn), rv)
    add(R["l"] - M64 * (s64["rg"] - n64["rg"]), s64["rg"])
    for k in range(s64["rQ"].shape[0]):
        add(R["Q"][k] - M64 * (s64["rQ"][k] - n64["rQ"][k]), s64["rQ"][k])
    if has_energy:
        add(R["E"] - M64 * (s64["rE"] - n64["rE"]), s64["rE"])
    return out


def be_solve(st_n, pb, M, dt, T0=None, tol=1e-7, maxit=20000, corr="upwind", order="vapor",
             vars_mode="vl", writeback="increment", energy="full", stats=None, floor_eps=16.0, precond="point",
             check="state", damp=1.0, sum_tol=None, sum_q=None, carry=None, stop_keys=("k", "v", "l", "Q", "E")):
    """解き切った後退 Euler: 点対角を前処理にした固定点反復を BE 残差 ≤ tol × (反復開始時の残差) まで (成分ごと)。
    成分ごとに丸め床 floor_eps·ε·max|M ρφ| 以下も収束とみなす (開始時の残差自体が丸め程度の成分 — 例: 一様な蒸気の
    ρY_w − ρg — は相対基準が意味を持たないため)。check='f64' は停止判定の残差を格納値から float64 で評価し直す
    (更新は dt のまま)。damp は増分の緩和係数 (1 = 緩和なし)。
    sum_tol (保存を見た停止): 停止にさらに |Σ_i r_q,i| ≤ sum_tol·Σ_i |M ρφ_q,i| (成分ごと, 格納値から float64 で評価) を課す。
    BE 残差の和はそのステップの総量変化そのもの (Σ_i R_i = 0) なので、これが 1 ステップの保存誤差を直接抑える。
    sum_q (量子化尺度の保存停止): |Σ_i r_q,i| ≤ sum_q·ε·sqrt(Σ_i (M ρφ_q,i)²) — 残した残差の和を格納値の丸めの酔歩と同じ大きさまで下げる
    (それより下は格納の丸めで決まり反復では下げられない)。
    carry (残差の持ち越し; 停止則ではなく離散式の変更): 前ステップで解き残した BE 残差 (成分ごとの float64 配列) を今ステップの右辺に足す。
    累積の保存誤差が最後のステップの残差だけになる (telescoping)。終了時の残差は stats['r_final'] に返す。
    **反復上限で返った解は「解き切っていない」**: stats['cap'] を数え、info['capped'] = True を返す。呼び出し側はそれを FAIL として扱う
    (codex diagnose 2026-10-02 ④)。stats: 'conv' (相対基準で停止), 'floor' (丸め床で停止), 'cap' (上限)。"""
    if stats is None:
        stats = {}
    has_energy = pb.thermo is not None and T0 is not None
    st = copy_state(st_n)
    T = T0
    r0 = None
    hist = []
    for it in range(maxit):
        if has_energy:
            T = pb.thermo.T_from_state(st, T)
        F = fluxes(st, T, pb, dt, corr=corr, order=order, energy=energy)
        Rr = residuals(F, pb, dt)
        if carry is not None:
            for k_ in carry:
                Rr[k_] = Rr[k_] + carry[k_].astype(dt)
        if check == "f64" or carry is not None:
            st64, st64n = cast_state(st, np.float64), cast_state(st_n, np.float64)
            R64 = residuals(fluxes(st64, T, pb, np.float64, corr=corr, order=order, energy=energy), pb, np.float64)
            if carry is not None:
                for k_ in carry:
                    R64[k_] = R64[k_] + carry[k_]
            nr = be_residual_norms(st64, st64n, R64, np.asarray(M, np.float64), has_energy)
            M64 = np.asarray(M, np.float64)
            rf = {"k": R64["k"] - M64 * (st64["rY"] - st64n["rY"]), "l": R64["l"] - M64 * (st64["rg"] - st64n["rg"]),
                  "v": R64["v"] - M64 * ((st64["rYw"] - st64["rg"]) - (st64n["rYw"] - st64n["rg"])),
                  "Q": R64["Q"] - M64 * (st64["rQ"] - st64n["rQ"])}
            if has_energy:
                rf["E"] = R64["E"] - M64 * (st64["rE"] - st64n["rE"])
            stats["r_final"] = rf
        else:
            nr = be_residual_norms(st, st_n, Rr, M, has_energy)
        nr = {k: v for k, v in nr.items() if k in stop_keys}   # 停止判定に使う成分 (#4b: 原方程式の変数 k,w,l,Q,E に揃える選択肢)
        hist.append(nr)
        if r0 is None:
            r0 = nr
            fl = {k: floor_eps * float(np.finfo(dt).eps) * float(np.max(np.abs((M * v).astype(np.float64))))
                  for k, v in (("k", st["rY"]), ("v", st["rYw"] - st["rg"]), ("w", st["rYw"]), ("l", st["rg"]), ("Q", st["rQ"]), ("E", st["rE"]))}
            if has_energy and "q_abs" in F:   # エネルギーは面のエンタルピー流束の打ち消し前の大きさも尺度に入れる
                qa = np.zeros(pb.N)
                np.add.at(qa, pb.L, F["q_abs"])
                np.add.at(qa, pb.R, F["q_abs"])
                fl["E"] = fl["E"] + floor_eps * float(np.finfo(dt).eps) * float(np.max(qa))
        rel = max((nr[k] / r0[k]) if r0[k] > 0 else (0.0 if nr[k] == 0 else math.inf) for k in nr)
        sum_ok = True
        if sum_tol is not None or sum_q is not None:
            sums = be_residual_sums(st, st_n, pb, T, M, corr, order, energy, has_energy)
            if sum_tol is not None:
                sum_ok = all(a <= sum_tol * b for a, b, _ in sums)
            if sum_q is not None:
                sum_ok = sum_ok and all(a <= sum_q * float(np.finfo(dt).eps) * q2 for a, _, q2 in sums)
        if rel <= tol and sum_ok:
            stats["conv"] = stats.get("conv", 0) + 1
            return st, T, it, rel
        if sum_ok and all(nr[k] <= max(tol * r0[k], fl[k]) for k in nr):
            stats["floor"] = stats.get("floor", 0) + 1
            return st, T, it, rel
        d = diagonals(st, F, pb, dt, corr=corr)
        if precond == "line":
            new = line_update(st, st_n, Rr, d, M, pb, dt, writeback, has_energy)
        else:
            new = jacobi_update(st, st_n, Rr, d, M, dt, vars_mode, writeback, has_energy)
        if damp != 1.0:
            w = dt(damp)
            for k in ("rY", "rYw", "rg", "rQ", "rE"):
                new[k] = st[k] + w * (new[k] - st[k])
        st = new
    stats["cap"] = stats.get("cap", 0) + 1
    stats["last_hist"] = hist
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
def s1_tolerances(st32, pb, rho_f_g, Dk):
    """入力量子化と演算誤差の許容 (面ごと・種ごと; 1 面問題)。z_k = ρY_k/(ρ−ρg) を Σz で正規化した値の 1 次の誤差伝播:
      入力量子化 (格納 float32 の各入力 x に |δx| ≤ ε/2·|x|): |δz_k| ≤ ε z_k c_q,  c_q = 1 + (ρ+ρg)/ρ_g  (分子・分母と Σz の正規化)
      演算 (同じ格納入力を float32 で計算): |δz_k| ≤ ε z_k c_a,  c_a = 3 + n_G  (減算・除算・n_G 項の和・正規化の除算)
      流束: A_k = ρ_g,f D_k a (|δz_k,L| + |δz_k,R|)、補正込みで tol_k = A_k + z_f,k Σ_j A_j (補正前の j⁰ の誤差が補正項に入る分)。
    物理の許容 (1e-6 ρ_g maxD |∇g|) とは別に足す (codex diagnose 2026-10-02 ④: 入力量子化と演算誤差を分ける)。"""
    s64 = cast_state(st32, np.float64)
    z, rhog, _ = gas_composition(s64, np.float64)
    nG = z.shape[0]
    a = float(pb.a[0])
    cq = 1.0 + (s64["rho"] + s64["rg"]) / rhog
    out = {}
    for name, c in (("q", cq), ("a", np.full(2, 3.0 + nG))):
        A = rho_f_g * Dk * a * EPS32 * (z[:, 0] * c[0] + z[:, 1] * c[1])
        zf = 0.5 * (z[:, 0] + z[:, 1])
        out[name] = A + zf * A.sum()
    return out


def s1_case(dx, gL, Dbin, Dfixed, zg, MW):
    """1 面 (dx, ∇g = 0.1) で 新 (§4.2) / 旧案 の流束を、厳密入力 (float64)・格納 float32 入力を float64 で・float32 で計算。"""
    rho = 1.0
    gR = gL + 0.1 * dx
    g = np.array([gL, gR])
    rY = np.vstack([(1 - g) * zg[0], (1 - g) * zg[1]]) * rho
    rv = (1 - g) * zg[2] * rho
    pb = Problem(2, [1.0 / dx], [0.0], Dbin, MW, [dx, dx], Dfixed=Dfixed)
    st64 = make_state([rho, rho], rY, rv + g * rho, g * rho, np.zeros((1, 2)), [0, 0], np.float64)
    st32 = cast_state(st64, np.float32)
    st32in64 = cast_state(st32, np.float64)
    r = {}
    for design in ("new", "old"):
        r[(design, "exact")] = fluxes(st64, None, pb, np.float64, design=design)["j"][:, 0]
        r[(design, "q")] = fluxes(st32in64, None, pb, np.float64, design=design)["j"][:, 0]
        r[(design, "f32")] = fluxes(st32, None, pb, np.float32, design=design)["j"].astype(np.float64)[:, 0]
    Df = fluxes(st64, None, pb, np.float64)["Df"][:, 0]
    rhog_f = rho * (1 - 0.5 * (gL + gR))
    lim = 1e-6 * rhog_f * Df.max() * 0.1
    tol = s1_tolerances(st32, pb, rhog_f, Df)
    return r, lim, tol, Df


def test_structure(quick):
    print("\n=== S1 分子流束の構造 (i): 三成分・気相一様・g のみ勾配 (∇g = 0.1 m⁻¹) ===")
    zg = np.array([0.2, 0.3, 0.5])     # [非水 0, 非水 1, 蒸気]
    MW = np.array([0.028, 0.028, 0.028])
    Dbin = np.array([[0, 1, 2], [1, 0, 3], [2, 3, 0]], float) * 1e-5   # D_12, D_13, D_23
    # (S1a) 作用素試験: plan §6 の字義どおり D_k = [1,2,3]e-5 を種の拡散係数として固定 (混合平均しない)
    # (S1b) 係数試験: 二元 D を混合平均した D_k (補数形) — 係数の値と、その係数での作用素
    Dmix_ref = np.array([0.8 / 0.55, 0.7 / (0.2 + 0.5 / 3), 0.5 / 0.2]) * 1e-5   # 手計算 (等分子量で X = z)
    Dmix64 = mixavg_D(zg[:, None], Dbin)[:, 0]
    Dmix32 = mixavg_D(zg.astype(np.float32)[:, None], Dbin.astype(np.float32))[:, 0].astype(np.float64)
    e64 = float(np.max(np.abs(Dmix64 / Dmix_ref - 1)))
    e32 = float(np.max(np.abs(Dmix32 / Dmix_ref - 1)))
    verdict(e64 <= 1e-12 and e32 <= 8 * EPS32,
            f"S1b 係数: 混合平均 D_k = [{Dmix64[0]:.6e}, {Dmix64[1]:.6e}, {Dmix64[2]:.6e}] と手計算の相対差 float64 {e64:.1e} (≤1e-12)、float32 {e32:.1e} (≤8ε₃₂)")
    for tag, Dfixed in (("S1a 作用素 (D_k 固定 [1,2,3]e-5)", np.array([1.0, 2.0, 3.0]) * 1e-5), ("S1b 作用素 (混合平均 D_k)", None)):
        for dx, gL in ((1.0, 0.1), (1.0e-3, 0.1)):
            r, lim, tol, Df = s1_case(dx, gL, Dbin, Dfixed, zg, MW)
            geo = f"dx={dx:g} m, Δg={0.1*dx:g}/面"
            je = np.abs(r[("new", "exact")])
            jq = np.abs(r[("new", "q")])
            ja = np.abs(r[("new", "f32")] - r[("new", "q")])
            verdict(je.max() <= lim, f"{tag} [{geo}] 厳密入力 float64: max|j_k| = {je.max():.3e} ≤ 1e-6·ρ_g·maxD·|∇g| = {lim:.3e}")
            verdict(np.all(jq <= tol["q"] + lim), f"{tag} [{geo}] 入力量子化 (格納 float32 を float64 で): max|j_k| = {jq.max():.3e}、"
                    f"許容 (量子化 + 物理) 最小 {np.min(tol['q'] + lim):.3e} (比の最大 {np.max(jq / (tol['q'] + lim)):.3f})")
            verdict(np.all(ja <= tol["a"]), f"{tag} [{geo}] 演算誤差 (float32 − 同じ入力の float64): max = {ja.max():.3e}、"
                    f"許容最小 {np.min(tol['a']):.3e} (比の最大 {np.max(ja / tol['a']):.3f})")
            jo = np.abs(r[("old", "exact")])
            oldv = r[("old", "exact")]
            sep = f" 旧案 j/ρ = [{oldv[0]:.4e}, {oldv[1]:.4e}, {oldv[2]:.4e}] m/s" + \
                  (" (codex 2026-09-27: [-1.3182e-7, -6.1364e-8, 1.9318e-7])" if Dfixed is None else "")
            verdict(np.max(jo) > np.max(tol["q"]) + lim,
                    f"{tag} [{geo}] 判別: 旧案の偽流束 max {jo.max():.3e} は許容 (量子化 + 物理) の {jo.max()/(np.max(tol['q'])+lim):.1e} 倍で FAIL すること —{sep}")


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
        default = (corr == "upwind" and order == "vapor")   # 採用 (codex diagnose 2026-10-02 ①)
        msg1 = f"S2(ii) |Σ j_k| / (8ε₃₂Σ|j⁰|) 最大 {ws:.3f} [{tag}, {nF} 面]"
        msg2 = f"S2(ii) |J_w−J_v−J_l| / (8ε₃₂Σ|J|) 最大 {wi:.3f} [{tag}, {nF} 面]"
        if default:
            verdict(ws <= 1.0, msg1)
            verdict(wi <= 1.0, msg2)
        else:
            info(msg1 + (" (≤1)" if ws <= 1 else " (>1)"))
            info(msg2 + (" (≤1)" if wi <= 1 else " (>1)"))


# ===================================================================== S3 3 セル判別
# float32 の「解き切り」の停止則 (§5.1 #4a (2))。反復上限で返った更新は FAIL (be_solve の stats['cap'])。
#   floor kε: 成分ごとに BE 残差 ≤ 1e-7 相対 または ≤ kε·max|Mρφ|。carry: 解き残した BE 残差を次ステップの右辺へ持ち越す (停止則ではなく
#   離散式の変更 — 累積誤差が最後のステップの残差だけになる; 未承認の候補)。
STOP_RULES_F32 = [("床 16ε", 16.0, False), ("床 8ε", 8.0, False), ("床 6ε", 6.0, False), ("床 4ε", 4.0, False), ("床 16ε + 残差持ち越し", 16.0, True)]


def run_be_series(st0, pb, M, dt, nstep, fe, carry, maxit, T0=None, rec_at=(), **kw):
    """解き切った BE 更新を nstep 回。戻り値: (最終状態, 最終 T, {n: 状態}, stats, 最大反復数)。"""
    st, T, stats, nit_max, car, recs = copy_state(st0), T0, {}, 0, None, {}
    for n in range(1, nstep + 1):
        st, T, nit, _ = be_solve(st, pb, M, dt, T0=T, maxit=maxit, floor_eps=fe, stats=stats,
                                 check=("f64" if carry else "state"), carry=car, **kw)
        nit_max = max(nit_max, nit)
        if carry:
            car = stats["r_final"]
        if n in rec_at:
            recs[n] = copy_state(st)
    return st, T, recs, stats, nit_max


def three_cell_state(dt):
    g = np.array([0.1, 0.2, 0.2])
    Yw = g + 0.1
    rQ = np.vstack([g * 3.0, g * 2.0, g * 1.0])   # Q も質量当たり g に比例 (値は任意)
    return make_state(np.ones(3), (1 - Yw)[None, :], Yw, g, rQ, np.zeros(3), dt)


def test_three_cell(quick):
    print("\n=== S3 3 セル判別 (A = 点対角 1 回更新 / B = 解き切った後退 Euler) ===")
    out = {}
    for dtn, dt, fe in [("float64", np.float64, 16.0), ("float32", np.float32, 16.0)]:
        pb = Problem(3, [1, 1], [1, 1], np.array([[0, 1e-5], [1e-5, 0]]), np.array([0.028, 0.018]), np.ones(3), molc=0.0)
        M = dt(1.0)
        st0 = three_cell_state(dt)
        t0 = totals(st0, pb.V)
        for scheme in (("A", "B") if dtn == "float64" else ("A",)):
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
            out[(dtn, scheme, fe)] = (rec, nit_max, stats.get("floor", 0), stats.get("cap", 0))
    for (dtn, scheme, fe), (rec, nit, nfl, ncap) in out.items():
        for n in (1, 1000):
            rels, vmin, sY = rec[n]
            worst = max(abs(v) for v in rels.values())
            line = (f"S3 {scheme} [{dtn}] {n} 更新: 総液量 {rels['rg']:+.7e}、総水分 {rels['rYw']:+.3e}、非水 {rels['rY0']:+.3e}、"
                    f"Q0 {rels['rQ0']:+.3e} (最大 |相対| {worst:.3e})、min ρv {vmin:.3e}、max|ΣρY−ρ|/ρ {sY:.2e}")
            if scheme == "B":
                line += f" (反復 ≤{nit}; BE 残差 ≤1e-7 相対か丸め床 {fe:g}ε 以下で停止、床で停止 {nfl} 回、上限到達 {ncap} 回)"
            if scheme == "B" and dtn == "float64":
                verdict(worst <= 1e-6 and vmin >= 0 and sY <= 8 * EPS32 and ncap == 0, line)
            elif scheme == "A" and dtn == "float64":
                print("[A]    " + line)
            else:
                info(line)
    out = {(k[0], k[1]): (v[0],) for k, v in out.items() if k[2] == 16.0}
    a_fail = max(abs(v) for v in out[("float64", "A")][0][1][0].values()) > 1e-6
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
    # float32 の B: 停止則ごとに 1000 更新の累積保存と上限到達 (上限 = FAIL)
    pb = Problem(3, [1, 1], [1, 1], np.array([[0, 1e-5], [1e-5, 0]]), np.array([0.028, 0.018]), np.ones(3), molc=0.0)
    st0 = three_cell_state(np.float32)
    t0 = totals(st0, pb.V)
    info("S3 は分子拡散 0 なので補正項 (z_c Σj⁰) は 0 で、補正の面 z (upwind = 採用 / mean = 現行) は結果に効かない (下の 6ε で両方を示す)")
    for name, fe, carry, corr in [r_ + ("upwind",) for r_ in STOP_RULES_F32] + [("床 6ε", 6.0, False, "mean")]:
        st, _, _, stats, nit = run_be_series(st0, pb, np.float32(1.0), np.float32, 1000, fe, carry, 5000, corr=corr)
        name = f"z={corr}, {name}"
        t = totals(st, pb.V)
        rels = {k: (t[k] - t0[k]) / abs(t0[k]) for k in t0 if k not in ("rE", "absE") and t0[k] != 0}
        worst = max(abs(v) for v in rels.values())
        cap = stats.get("cap", 0)
        ok = worst <= 1e-6 and cap == 0
        line = (f"S3 B [float32, {name}] 1000 更新: 最大 |相対| {worst:.3e} (総液量 {rels['rg']:+.2e}, 総水分 {rels['rYw']:+.2e})、"
                f"反復上限到達 {cap} 回、最大反復 {nit} → {'合格' if ok else '不合格'}")
        if fe == 6.0 and not carry and corr == "upwind":
            verdict(ok, line + " [#4b の固定条件]")
        else:
            info(line)


# ===================================================================== S4 非負
def test_nonnegativity(quick):
    print("\n=== S4 非負 ===")
    # 2 セル例: 蒸気 0、Y_w = g = 0.1/0.2、分子・乱流の面係数 1、V/Δτ = 1、ρ = 1
    for dtn, dt in [("float64", np.float64), ("float32", np.float32)]:
        g = np.array([0.1, 0.2])
        st0 = make_state([1, 1], (1 - g)[None, :], g, g, (g * 2)[None, :], [0, 0], dt)
        pb = Problem(2, [1.0], [1.0], np.array([[0, 1.0], [1.0, 0]]), np.array([0.028, 0.018]), np.ones(2), molc=1.0)
        M = dt(1.0)
        st4 = {}
        stB, _, nit, rel = be_solve(st0, pb, M, dt, stats=st4)
        rvB = stB["rYw"].astype(np.float64) - stB["rg"].astype(np.float64)
        _, corr, nact, _ = realizability_clamp(stB, pb.V)
        line = (f"S4 2 セル [{dtn}] 解き切った擬似時間更新 (新, 蒸気/液変数): ρv = [{rvB[0]:.3e}, {rvB[1]:.3e}]、クランプ補正 {corr:.3e} "
                f"({nact} 件; 反復 {nit}, 上限到達 {st4.get('cap', 0)})")
        if dtn == "float64":
            verdict(rvB.min() >= 0 and corr == 0 and st4.get("cap", 0) == 0, line)
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
    # 補正の面 z の反例 (両セル g = 0; codex diagnose 2026-10-02 ①⑤): 拡散係数の差が大きい三成分で、蒸気 0 のセルから補正流束が蒸気を
    # 持ち出す例。mean は −z_v,f S の隣セル分が非対角に負で入るので M 行列でなくなる。風上 (補正流束 −Σj⁰ の流出側の z) を採用。
    #   L: (A 0, B 1, 蒸気 0)、R: (A 0.99, B 0, 蒸気 0.01)、g = 0、二元 D_AB = D_Av = 10、D_Bv = 0.1 (相対値; D 比 100)、面係数・V/Δτ = 1
    def cex(dt):
        st0 = make_state([1, 1], np.array([[0.0, 0.99], [1.0, 0.0]]), [0.0, 0.01], [0.0, 0.0], np.zeros((1, 2)), [0, 0], dt)
        pb = Problem(2, [1.0], [0.0], np.array([[0, 10, 10], [10, 0, 0.1], [10, 0.1, 0]]), np.array([0.028, 0.028, 0.028]), np.ones(2))
        return st0, pb
    for corr in ("upwind", "mean"):
        for dtn, dt in (("float64", np.float64), ("float32", np.float32)):
            st0, pb = cex(dt)
            M = dt(1.0)
            F = fluxes(st0, None, pb, dt, corr=corr)
            Rr = residuals(F, pb, dt)
            d = diagonals(st0, F, pb, dt, corr=corr)
            vJ = float((jacobi_update(st0, st0, Rr, d, M, dt)["rYw"][0]))
            res = {}
            for damp in (1.0, 0.5):
                stats = {}
                with np.errstate(all="ignore"):   # mean の緩和なしは発散してオーバーフローする (記録だけ)
                    stB, _, nit, rel = be_solve(st0, pb, M, dt, corr=corr, maxit=5000, stats=stats, damp=damp)
                hist = stats.get("last_hist", [])
                res[damp] = (stB, nit, stats, hist)
            sol = None
            if dtn == "float64":
                try:
                    from scipy.optimize import fsolve

                    def Fn(x):
                        s_ = copy_state(st0)
                        s_["rY"] = x[:4].reshape(2, 2)
                        s_["rYw"] = x[4:6]
                        R_ = residuals(fluxes(s_, None, pb, np.float64, corr=corr), pb, np.float64)
                        return np.concatenate([(R_["k"] - (s_["rY"] - st0["rY"])).ravel(), R_["w"] - (s_["rYw"] - st0["rYw"])])
                    x, _, ier, _ = fsolve(Fn, np.concatenate([st0["rY"].ravel(), st0["rYw"]]), full_output=True, xtol=1e-14)
                    sol = (x[4:6], float(np.max(np.abs(Fn(x)))), ier)
                except ImportError:
                    sol = None
            st1, n1, s1, h1 = res[1.0]
            st5, n5, s5, h5 = res[0.5]
            v1 = (st1["rYw"] - st1["rg"]).astype(np.float64)
            v5 = (st5["rYw"] - st5["rg"]).astype(np.float64)
            cap1, cap5 = s1.get("cap", 0) > 0, s5.get("cap", 0) > 0
            hs = ""
            if cap1 and h1:
                pick = [0, 1, 2, 3, 5, 10, 20, 50, 100, 500, 1000, 4999]
                hs = "; 残差履歴 (max|r_k|, max|r_v|) " + ", ".join(f"#{i}: ({h1[i]['k']:.2e}, {h1[i]['v']:.2e})" for i in pick if i < len(h1))
            tag = f"S4 反例 [補正 z = {corr}, {dtn}]"
            if corr == "upwind":
                verdict(vJ >= 0, f"{tag} 点対角 1 回: 蒸気 0 のセル L の ρv {vJ:+.4e} ≥ 0")
                verdict(not cap1, f"{tag} 緩和なし点対角の固定点反復: {'上限 5000 に到達 (収束しない)' if cap1 else f'{n1} 回で収束'}、"
                        f"最終 ρv = [{v1[0]:+.3e}, {v1[1]:+.3e}]{hs}")
                verdict(not cap5 and v5.min() >= 0, f"{tag} 緩和 0.5 の固定点反復: {'上限に到達' if cap5 else f'{n5} 回で収束'}、ρv = [{v5[0]:+.6e}, {v5[1]:+.6e}]")
            else:
                info(f"{tag} 点対角 1 回 ρv_L {vJ:+.4e}; 緩和なし {'上限到達' if cap1 else f'{n1} 回'} / 緩和 0.5 {'上限到達' if cap5 else f'{n5} 回'}、"
                     f"ρv = [{v5[0]:+.6e}, {v5[1]:+.6e}]")
            if sol is not None:
                info(f"{tag} 参照 (Newton, scipy fsolve): BE 解は存在し ρv = [{sol[0][0]:+.6e}, {sol[0][1]:+.6e}] (|F| {sol[1]:.1e}, ier {sol[2]})")


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
        if dtn == "float64" and wb == "increment" and corr == "upwind":
            verdict(tot_corr == 0 and n_act == 0 and tot_rn == 0, line)
        else:
            info(line)


# ===================================================================== S7 固定点判別 (P1 / 非分割)
def vl_limit_commit(rho, rYw, rg, rQ, dv, dg, dQ, dt, dg_max=math.inf, dT_per_dg=0.0, dT_max=math.inf, frac=1.0):
    """蒸気・液の制限と commit を一体で定義する (§5.1 #4a; codex diagnose 2026-10-02 ②)。入力の増分は前処理済み (δρv, δρg, δρQ_n)。
    セル配列でもスカラーでもよい (rQ, dQ は (nQ, N) か (nQ,))。
      1. 閾値 θ_thr = min(1, dg_max/|δg/ρ|, dT_max/(|δg/ρ|·dT_per_dg))  — 全増分 (輸送 + ソース + BDF) に対して。
      2. 対 (蒸気, 液) の非負 θ_vg = min(1, ρv/(−δρv) [δρv<0], ρg/(−δρg) [δρg<0])  — 蒸気は**更新前の** ρv = ρY_w − ρg で見る
         (現行の液用 θ は更新済み総水分の avail を見て、avail = 0 で枯渇制限が掛からない; condensationUpdateLimiter_d.cuh:78-80)。
      3. θ = min(θ_thr, θ_vg) を δρv・δρg・δρQ の**全部に共通**に掛ける (蒸気と液を別々に切らない — 片方だけ切ると総水分が増減する)。
      4. Q は成分ごとの非負化 (θδQ_n < −ρQ_n なら −ρQ_n に切る; 切った量を補正として数える)。
      5. commit: ρg' = ρg + θδρg、ρY_w' = ρY_w + fl(θδρv + θδρg) (増分形: 増分 0 なら格納値は不変)。
         丸めで fl(ρY_w' − ρg') < 0 になったら ρY_w' = ρg' (補正として数える; 正確な算術では θ_vg が蒸気 ≥ 0 を保証する)。
    frac (< 1 は #4b の診断用): 非負の θ_vg を境界までの距離の frac 倍に止める (1 = 境界ちょうど; 既定)。
    残差が 0 なら θ に依らず増分 0 → 固定点は動かない。
    戻り値: (ρY_w', ρg', ρQ', θ (float64 配列), 補正量 dict — limiter_v/g は θ で保留した増分 (質量の補正ではない)、
    Q_cut・v_round は状態を書き換えた補正 (非保存))。"""
    f64 = lambda x: np.asarray(x, dtype=np.float64)
    rv = f64(rYw) - f64(rg)
    adg = np.abs(f64(dg)) / f64(rho)
    th = np.ones_like(adg)
    with np.errstate(divide="ignore", invalid="ignore"):
        th = np.where(adg > 0, np.minimum(th, dg_max / adg), th)
        dTp = np.broadcast_to(f64(dT_per_dg), adg.shape)
        th = np.where((adg > 0) & (dTp > 0), np.minimum(th, dT_max / (adg * dTp)), th)
        th = np.where(f64(dv) < 0, np.minimum(th, frac * rv / (-f64(dv))), th)
        th = np.where(f64(dg) < 0, np.minimum(th, frac * f64(rg) / (-f64(dg))), th)
    th = np.maximum(th, 0.0)
    T_ = th.astype(dt)
    one = dt(1.0)
    corr = {"limiter_v": float(np.sum(f64((one - T_) * np.abs(dv)))), "limiter_g": float(np.sum(f64((one - T_) * np.abs(dg))))}
    dQn = T_ * dQ
    cut = np.minimum(dQn + rQ, 0)
    corr["Q_cut"] = float(np.sum(np.abs(f64(cut))))
    rQn = rQ + (dQn - cut)
    rg_n = rg + T_ * dg
    rYw_n = rYw + (T_ * dv + T_ * dg)
    neg = (rYw_n - rg_n) < 0
    corr["v_round"] = float(np.sum(np.where(neg, f64(rg_n) - f64(rYw_n), 0.0)))
    rYw_n = np.where(neg, rg_n, rYw_n).astype(dt)
    if np.ndim(rYw_n) == 0:
        rYw_n, rg_n, th = dt(rYw_n), dt(rg_n), float(th)
    return rYw_n, rg_n, rQn, th, corr


def vl_increments(scheme, Rw_tr, Rg_tr, S, Pv, Pg):
    """A = P1 (輸送とソースを別の前処理で分割): δρg = (R_gᵗʳ + S)/P_g、δρv = (R_wᵗʳ − R_gᵗʳ)/P_v − S/P_g。
    B = 非分割の全残差: R_g = R_gᵗʳ + S (+ BDF)、R_v = R_w − R_g、δρv = R_v/P_v、δρg = R_g/P_g。"""
    if scheme == "A":
        return (Rw_tr - Rg_tr) / Pv - S / Pg, (Rg_tr + S) / Pg
    Rg = Rg_tr + S
    return (Rw_tr - Rg) / Pv, Rg / Pg


# ===================================================================== 新経路 B (全残差変換 + vl_limit_commit) の BE 反復 (§5.1 #4b)
B_LIMITS = {"dg_max": 5.0e-3, "dT_max": 1.0}   # forge 既定 condDgMaxStep / condDTmaxStep (input/solverConfig.hpp:647-648)。#4b では実行前にこの値で固定


class GrowthSource:
    """相変化ソースの簡易成長則 (#4b (3) の反復試験用; forge の CNT/Gyarmathy ではない)。体積当たり [kg/m³/s]:
        S = (ρv − ρv_sat(T))/τ                    (ρv ≥ ρv_sat: 凝縮)
        S = (ρv − ρv_sat(T))/τ · ρg/(ρg + ρg_ref)  (ρv < ρv_sat: 蒸発; 液 0 で 0)
        ρv_sat = p_sat(T)/(R_w T)、p_sat = p0 exp((L0/R_w)(1/T0 − 1/T)) (定 L0 の Clausius–Clapeyron)。
    S は液の残差に +S V、総水分には入らない (蒸気 → 液の移し替え)。エネルギーにもソースは無く、潜熱は EOS e = e_gas + g(R_wT − L)
    経由で温度に出る (forge の carrier 二相 EOS と同じ扱い)。前処理用の Jacobian (≥ 0):
        J_g = −∂S/∂(ρg) ≈ f (dρv_sat/dT)(∂T/∂ρg)/τ − ((ρv − ρv_sat)/τ)·∂f/∂ρg、∂T/∂ρg|_{ρe,ρY_w} = (L − R_wT)/(ρ c_v,eff)、
        ∂f/∂ρg = ρg_ref/(ρg + ρg_ref)² (蒸発側)、J_v = ∂S/∂(ρv) ≈ f/τ
        (f = 1 凝縮 / ρg/(ρg+ρg_ref) 蒸発。組成変化による T の依存は無視 — 前処理にしか使わない)。"""

    def __init__(self, tau, rg_ref, Rw, p0=611.2, T0=273.15, L0=2.501e6):
        self.tau, self.rg_ref, self.Rw, self.p0, self.T0, self.L0 = tau, rg_ref, Rw, p0, T0, L0

    def __call__(self, st64, T, thermo):
        rho, rg = st64["rho"], st64["rg"]
        rv = st64["rYw"] - rg
        ps = self.p0 * np.exp((self.L0 / self.Rw) * (1 / self.T0 - 1 / T))
        rvs = ps / (self.Rw * T)
        d = rv - rvs
        f = np.where(d >= 0, 1.0, np.maximum(rg, 0) / (np.maximum(rg, 0) + self.rg_ref))
        S = d / self.tau * f
        Yk = st64["rY"] / rho
        cv = thermo.cv_eff(T, Yk, st64["rYw"] / rho, rg / rho)
        L = thermo.lat.L(T)
        dTdg = (L - self.Rw * T) / (rho * cv)
        drvs = rvs * (self.L0 / (self.Rw * T * T) - 1 / T)
        dfdg = np.where(d >= 0, 0.0, self.rg_ref / (np.maximum(rg, 0) + self.rg_ref) ** 2)
        Jg = np.maximum(f * drvs * dTdg / self.tau - d / self.tau * dfdg, 0.0)   # 蒸発側の ∂f/∂ρg (液 → 0 で支配的) を含む
        Jv = f / self.tau
        return S, Jg, Jv


def div_faces(Fx, pb, dt):
    Fx = np.atleast_2d(Fx)
    r = np.zeros((Fx.shape[0], pb.N), dtype=dt)
    np.add.at(r, (slice(None), pb.L), -Fx)
    np.add.at(r, (slice(None), pb.R), Fx)
    return r


def full_residuals(st, st_n, pb, T, M, dt, corr="upwind", T64=None):
    """原方程式 (保存形 BE) の残差 r = R(φ) − M(ρφ − ρφⁿ) を (k, w, l, v_dir, Q, E) で。R = 拡散 (+ 移流 pb.mdot) (+ 相変化 pb.source)。
    w = 総水分 (ソース無し)、l = 液 (+S V)、v_dir = 蒸気の直接の流束残差 (A の line_update が使う形; 比較用)。精度 dt で評価。"""
    F = fluxes(st, T, pb, dt, corr=corr)
    R = residuals(F, pb, dt)
    rho = st["rho"]
    extra_diag = np.zeros(pb.N, dt)
    if getattr(pb, "mdot", None) is not None:
        m = pb.mdot.astype(dt)
        up = lambda phi: np.where(m >= 0, phi[..., pb.L], phi[..., pb.R])
        Yk, g = st["rY"] / rho, st["rg"] / rho
        v = (st["rYw"] - st["rg"]) / rho
        Q, e = st["rQ"] / rho, st["rE"] / rho
        Fv, Fl = m * up(v), m * up(g)
        R["k"] = R["k"] + div_faces(m * up(Yk), pb, dt)       # 移流流束 F = ṁ φ_up (L→R): L から出て R へ入る
        R["v"] = R["v"] + div_faces(Fv, pb, dt)[0]
        R["l"] = R["l"] + div_faces(Fl, pb, dt)[0]
        R["w"] = R["w"] + div_faces(Fv + Fl, pb, dt)[0]
        R["Q"] = R["Q"] + div_faces(m * up(Q), pb, dt)
        if "E" in R:
            R["E"] = R["E"] + div_faces(m * up(e[None, :]), pb, dt)[0]
        np.add.at(extra_diag, pb.L, np.maximum(m, 0) / rho[pb.L])
        np.add.at(extra_diag, pb.R, np.maximum(-m, 0) / rho[pb.R])
    Jg = Jv = np.zeros(pb.N, dt)
    if getattr(pb, "source", None) is not None:
        S, Jg64, Jv64 = pb.source(cast_state(st, np.float64), np.asarray(T if T64 is None else T64, np.float64), pb.thermo)
        SV = (S * pb.V).astype(dt)
        R["l"] = R["l"] + SV
        R["v"] = R["v"] - SV
        Jg, Jv = (Jg64 * pb.V).astype(dt), (Jv64 * pb.V).astype(dt)
    r = {"k": R["k"] - M * (st["rY"] - st_n["rY"]), "w": R["w"] - M * (st["rYw"] - st_n["rYw"]),
         "l": R["l"] - M * (st["rg"] - st_n["rg"]),
         "v_dir": R["v"] - M * ((st["rYw"] - st["rg"]) - (st_n["rYw"] - st_n["rg"])),
         "Q": R["Q"] - M * (st["rQ"] - st_n["rQ"])}
    if "E" in R:
        r["E"] = R["E"] - M * (st["rE"] - st_n["rE"])
    return r, F, extra_diag, Jg, Jv


ORIG_KEYS = ("k", "w", "l", "Q", "E")


def orig_scales(st_n, pb, T, M, floor_eps, corr="upwind"):
    """原方程式の残差の判定尺度 (ステップ開始時): r0_q (= ステップ開始時の残差、float64) と丸め床 floor_eps·ε₃₂·max|M ρφ_q|。"""
    s64 = cast_state(st_n, np.float64)
    M64 = np.asarray(M, np.float64)
    r0, F0, _, _, _ = full_residuals(s64, s64, pb, T, M64, np.float64, corr=corr)
    fl = {"k": s64["rY"], "w": s64["rYw"], "l": s64["rg"], "Q": s64["rQ"], "E": s64["rE"]}
    fl = {k: floor_eps * EPS32 * float(np.max(np.abs(M64 * v))) for k, v in fl.items()}
    if "q_abs" in F0:   # E は面のエンタルピー流束の大きさの床も (be_solve / be_solve_B の停止と同じ定義)
        qa = np.zeros(pb.N)
        np.add.at(qa, pb.L, F0["q_abs"])
        np.add.at(qa, pb.R, F0["q_abs"])
        fl["E"] += floor_eps * EPS32 * float(np.max(qa))
    return {k: float(np.max(np.abs(v))) for k, v in r0.items() if k in ORIG_KEYS}, fl


def orig_residual_ratio(st, st_n, pb, T, M, r0, fl, tol=1e-7, corr="upwind"):
    """原方程式の残差を格納値から float64 で評価し、成分ごとに max|r_q| / max(tol·r0_q, 床_q) の最大を返す (≤1 で残差合格)。
    更新量ゼロ (θ=0 の停止) を残差合格の代用にしないための判定。"""
    s64, n64 = cast_state(st, np.float64), cast_state(st_n, np.float64)
    r, _, _, _, _ = full_residuals(s64, n64, pb, T, np.asarray(M, np.float64), np.float64, corr=corr)
    out = {}
    for k in ORIG_KEYS:
        if k not in r:
            continue
        den = max(tol * r0.get(k, 0.0), fl[k])
        out[k] = float(np.max(np.abs(r[k]))) / den if den > 0 else (0.0 if np.max(np.abs(r[k])) == 0 else math.inf)
    return max(out.values()), out


def be_solve_B(st_n, pb, M, dt, T0=None, tol=1e-7, maxit=3000, corr="upwind", floor_eps=6.0, precond="line",
               relax=1.0, limits=None, stats=None, src_block=True, frac=1.0, check64=False, stop_keys=("k", "v", "w", "l", "Q", "E"),
               floor_diag=False):
    """新経路 B の解き切った BE (§5.1 #4b)。各反復:
      (0) ステップ開始時に ρ を連続の式で更新 (pb.mdot; BE で厳密) し、全輸送量の基点を φⁿ ρ_new に移す (= forge の φ_N δρ 項)。
      (1) 原方程式の残差 r_k, r_w, r_l (+S V), r_Q, r_E を精度 dt で組む。
      (2) 全残差変換 r_v = r_w − r_l (float dt の減算; ソース・BDF・移流・拡散を全部含む)。
      (3) 前処理 (A と同じ: precond='line' は乱流の三重対角 + 点対角 [分子・補正の流出・移流の流出・ソース Jacobian]、'point' は点対角だけ):
          δρv = P_v⁻¹ r_v、δρg = P_g⁻¹ r_l、δρQ = P_g⁻¹ r_Q、δρY_k = P_k⁻¹ r_k、δρE = r_E/M。
          src_block=True (相変化ソースがあるとき): 蒸気と液を 2N の連成系として解く — 点ごとの 2×2 ブロック
            P_vv = M + P_vᵗʳ + V S_v、P_vg = V S_g、P_gv = −V S_v、P_gg = M + P_gᵗʳ − V S_g  (S_v = ∂S/∂ρv ≈ J_v, S_g = ∂S/∂ρg ≈ −J_g)
            (列の和 P_vv+P_gv, P_vg+P_gg はソースを含まない = 総水分の前処理はソースに反応しない)。False は対角だけ (P_v + V J_v, P_g + V J_g)。
      (4) **緩和 ω (relax) をここで全増分に掛ける** — 前処理の後・limiter の前 (limiter が見るのは緩和後の増分)。
      (5) vl_limit_commit (B_LIMITS の dg_max・dT_max、dT_per_dg = L/c_v,eff) で蒸気・液・Q を commit、非水種と ρE は増分をそのまま足す。
    停止: 原方程式の残差 (k, w, l, Q, E) が成分ごとに ≤ tol·r0 か ≤ 床 (floor_eps·ε·max|Mρφ|; E は流束の大きさの床も)。上限は FAIL。
    stats: 'cap'/'conv'/'floor'、'theta_min'、'theta_lt1' (θ<1 のセル·反復数)、'theta0' (最終反復で θ=0 のセル数)、補正量の累積。"""
    if stats is None:
        stats = {}
    lim = dict(B_LIMITS if limits is None else limits)
    has_energy = pb.thermo is not None and T0 is not None
    st = copy_state(st_n)
    if getattr(pb, "mdot", None) is not None:
        Rrho = div_faces(pb.mdot.astype(np.float64), pb, np.float64)[0]
        rho_new = (st_n["rho"].astype(np.float64) + Rrho / np.asarray(M, np.float64)).astype(dt)
        fac = rho_new / st_n["rho"]
        for k in ("rY", "rYw", "rg", "rQ", "rE"):
            st[k] = (st_n[k] * fac).astype(dt)
        st["rho"] = rho_new
    T = T0
    r0 = None
    w = dt(relax)
    for it in range(maxit):
        if has_energy:
            T = pb.thermo.T_from_state(st, T)
        r, F, xdiag, Jg, Jv = full_residuals(st, st_n, pb, T, M, dt, corr=corr)
        if check64:   # 停止判定の残差を格納値から float64 で評価 (原方程式の残差そのもの)
            r64, _, _, _, _ = full_residuals(cast_state(st, np.float64), cast_state(st_n, np.float64), pb, T, np.asarray(M, np.float64), np.float64, corr=corr)
            r64["v"] = r64["w"] - r64["l"]
            nr = {k: float(np.max(np.abs(r64[k]))) for k in stop_keys if k in r64}
        else:
            rr_ = dict(r)
            rr_["v"] = r["w"] - r["l"]
            nr = {k: float(np.max(np.abs(rr_[k].astype(np.float64)))) for k in stop_keys if k in rr_}
        if r0 is None:
            r0 = nr
            M64 = np.asarray(M, np.float64)
            fl = {"k": st["rY"], "w": st["rYw"], "v": st["rYw"] - st["rg"], "l": st["rg"], "Q": st["rQ"], "E": st["rE"]}
            fl = {k: floor_eps * float(np.finfo(dt).eps) * float(np.max(np.abs(M64 * v.astype(np.float64)))) for k, v in fl.items()}
            if has_energy and "q_abs" in F:
                qa = np.zeros(pb.N)
                np.add.at(qa, pb.L, F["q_abs"])
                np.add.at(qa, pb.R, F["q_abs"])
                fl["E"] += floor_eps * float(np.finfo(dt).eps) * float(np.max(qa))
            if floor_diag:
                # 床の尺度を M だけでなく前処理の対角全体 (M + 輸送 + 移流の流出 + V·ソース Jacobian) で取る (#4b 診断):
                # 格納値の丸め δ(ρφ) が残差に出る大きさは (M + 対角)·ε|ρφ| なので、硬いソースでは M だけの床が丸めより下になる。
                d0 = diagonals(st, F, pb, dt, corr=corr)
                xx = xdiag.astype(np.float64)
                eps_ = floor_eps * float(np.finfo(dt).eps)
                Mf0 = np.broadcast_to(np.asarray(M, np.float64), (pb.N,))
                rvv = np.abs((st["rYw"] - st["rg"]).astype(np.float64))
                rgg = np.abs(st["rg"].astype(np.float64))
                Jv_, Jg_ = Jv.astype(np.float64), Jg.astype(np.float64)
                # 残差の各成分が、各変数の格納丸め ε|x_j| から受ける変化 Σ_j |∂r/∂x_j| ε|x_j| (対角近似)
                fl["v"] = max(fl["v"], eps_ * float(np.max((Mf0 + d0["v"].astype(np.float64) + xx + Jv_) * rvv + Jg_ * rgg)))
                fl["l"] = max(fl["l"], eps_ * float(np.max((Mf0 + d0["l"].astype(np.float64) + xx + Jg_) * rgg + Jv_ * rvv)))
                fl["w"] = max(fl["w"], eps_ * float(np.max(np.abs((Mf0 + d0["v"].astype(np.float64) + xx) * st["rYw"].astype(np.float64)))))
                fl["k"] = max(fl["k"], eps_ * float(np.max(np.abs((Mf0 + d0["k"].astype(np.float64).max(axis=0) + xx) * st["rY"].astype(np.float64).max(axis=0)))))
        rel = max((nr[k] / r0[k]) if r0[k] > 0 else (0.0 if nr[k] == 0 else math.inf) for k in nr)
        if rel <= tol:
            stats["conv"] = stats.get("conv", 0) + 1
            stats["theta0_last"] = stats.get("_t0", 0)
            return st, T, it, rel
        if all(nr[k] <= max(tol * r0[k], fl[k]) for k in nr):
            stats["floor"] = stats.get("floor", 0) + 1
            stats["theta0_last"] = stats.get("_t0", 0)
            return st, T, it, rel
        d = diagonals(st, F, pb, dt, corr=corr)
        r_v = r["w"] - r["l"]                                    # (2) 全残差変換
        rho = st["rho"]
        aG = (pb.a * pb.Gam).astype(np.float64)
        turb = np.zeros(pb.N)
        np.add.at(turb, pb.L, aG / rho[pb.L].astype(np.float64))
        np.add.at(turb, pb.R, aG / rho[pb.R].astype(np.float64))
        Mf = np.broadcast_to(np.asarray(M, np.float64), (pb.N,))
        Tm = np.diag(turb)
        if precond == "line":
            Tm[pb.L, pb.R] -= aG / rho[pb.R].astype(np.float64)
            Tm[pb.R, pb.L] -= aG / rho[pb.L].astype(np.float64)

            def solve(rr, dpt):   # dpt: 点対角のうち乱流以外
                return np.linalg.solve(Tm + np.diag(Mf + dpt), rr.astype(np.float64)).astype(dt)
        else:
            def solve(rr, dpt):
                return (rr.astype(np.float64) / (Mf + turb + dpt)).astype(dt)
        x = xdiag.astype(np.float64)
        dk = np.vstack([solve(r["k"][k], d["k"][k].astype(np.float64) - turb + x) for k in range(st["rY"].shape[0])])
        if src_block and getattr(pb, "source", None) is not None:
            N_ = pb.N
            Tv = Tm + np.diag(Mf + d["v"].astype(np.float64) - turb + x)
            Tg = Tm + np.diag(Mf + x)
            Sv, Sg = Jv.astype(np.float64), -Jg.astype(np.float64)          # V·∂S/∂ρv, V·∂S/∂ρg
            A = np.zeros((2 * N_, 2 * N_))
            A[:N_, :N_] = Tv + np.diag(Sv)
            A[:N_, N_:] = np.diag(Sg)
            A[N_:, :N_] = -np.diag(Sv)
            A[N_:, N_:] = Tg - np.diag(Sg)
            sol = np.linalg.solve(A, np.concatenate([r_v.astype(np.float64), r["l"].astype(np.float64)]))
            dv, dg = sol[:N_].astype(dt), sol[N_:].astype(dt)
        else:
            dv = solve(r_v, d["v"].astype(np.float64) - turb + x + Jv.astype(np.float64))
            dg = solve(r["l"], x + Jg.astype(np.float64))
        dQ = np.vstack([solve(r["Q"][k], x) for k in range(st["rQ"].shape[0])])
        dE = (r["E"] / M) if "E" in r else None
        if relax != 1.0:                                          # (4) 緩和の位置
            dk, dv, dg, dQ = w * dk, w * dv, w * dg, w * dQ
            if dE is not None:
                dE = w * dE
        dTp = 0.0
        if has_energy:
            T64 = np.asarray(T, np.float64)
            s64 = cast_state(st, np.float64)
            dTp = pb.thermo.lat.L(T64) / pb.thermo.cv_eff(T64, s64["rY"] / s64["rho"], s64["rYw"] / s64["rho"], s64["rg"] / s64["rho"])
        rYw_n, rg_n, rQ_n, th, corr_ = vl_limit_commit(st["rho"], st["rYw"], st["rg"], st["rQ"], dv, dg, dQ, dt,
                                                       dg_max=lim["dg_max"], dT_per_dg=dTp, dT_max=lim["dT_max"], frac=frac)
        th = np.atleast_1d(th)
        stats["theta_min"] = min(stats.get("theta_min", 1.0), float(th.min()))
        stats["theta_lt1"] = stats.get("theta_lt1", 0) + int(np.count_nonzero(th < 1))
        stats["_t0"] = int(np.count_nonzero(th == 0))
        stats["theta0_cells"] = np.flatnonzero(th == 0)
        for k_, v_ in corr_.items():
            stats[k_] = stats.get(k_, 0.0) + v_
        new = copy_state(st)
        new["rY"] = st["rY"] + dk
        new["rYw"], new["rg"], new["rQ"] = rYw_n, rg_n, rQ_n
        if dE is not None:
            new["rE"] = st["rE"] + dE
        st = new
    stats["cap"] = stats.get("cap", 0) + 1
    stats["theta0_last"] = stats.get("_t0", 0)
    return st, T, maxit, rel


def test_fixed_point(quick):
    print("\n=== S7 固定点判別: 輸送とソースが釣り合う状態で 1 擬似更新 (A = P1 / B = 非分割の全残差) ===")
    for dtn, dt in (("float32", np.float32), ("float64", np.float64)):
        c = lambda x: dt(x)
        rho, rYw, rg = c(1.0), c(0.3), c(0.1)
        rQ = np.array([0.3, 0.2, 0.1], dt)
        dQ = np.zeros(3, dt)                 # Q の全残差 0
        Rw_tr, Rg_tr, S, Pv, Pg = c(0.0), c(-0.02), c(0.02), c(3.0), c(2.0)
        out = {}
        for sch in ("A", "B"):
            dv, dg = vl_increments(sch, Rw_tr, Rg_tr, S, Pv, Pg)
            rYw_n, rg_n, rQn, th, corr = vl_limit_commit(rho, rYw, rg, rQ, dv, dg, dQ, dt)
            out[sch] = (float(dv), float(dg), float(rYw_n) - float(rYw), float(rg_n) - float(rg),
                        float(np.max(np.abs(rQn.astype(np.float64) - rQ.astype(np.float64)))), th, corr)
        a, b = out["A"], out["B"]
        line = lambda t, o: (f"{t} [{dtn}]: δρv {o[0]:+.7e}, δρg {o[1]:+.7e}, ΔρY_w {o[2]:+.7e}, Δρg {o[3]:+.3e}, max|ΔρQ| {o[4]:.1e}, "
                             f"θ {o[5]:.3f}, 補正 {', '.join(f'{k} {v:.1e}' for k, v in o[6].items())}")
        verdict(abs(a[2] + 0.02 / 2 - 0.02 / 3) <= 1e-6, line("S7 A = P1", a) + " — ΔρY_w ≈ −0.003333 動くこと")
        verdict(all(x == 0 for x in b[:5]) and all(v == 0 for v in b[6].values()), line("S7 B = 非分割", b) + " — 全増分・補正が厳密に 0")
    # 情報: 非固定点での制限と commit の振る舞い
    for desc, (rYw, rg, Rw_tr, Rg_tr, S) in (("蒸気 0 で凝縮を要求 (ρY_w = ρg = 0.1, S = +0.01)", (0.1, 0.1, 0.0, 0.0, 0.01)),
                                           ("蒸気が足りない凝縮 (ρv = 0.002, S = +0.02)", (0.102, 0.1, 0.0, 0.0, 0.02)),
                                           ("輸送で蒸気流入・液流出 + 蒸発 (S = −0.01)", (0.3, 0.1, 0.01, -0.005, -0.01))):
        dt = np.float64
        dv, dg = vl_increments("B", Rw_tr, Rg_tr, S, 3.0, 2.0)
        rYw_n, rg_n, rQn, th, corr = vl_limit_commit(1.0, rYw, rg, np.array([0.3, 0.2, 0.1]), dv, dg, np.zeros(3), dt)
        info(f"S7 B の制限 [{desc}]: δρv {dv:+.4e}, δρg {dg:+.4e} → θ {th:.4f}, ρv' {rYw_n - rg_n:.4e}, ρg' {rg_n:.4e}, "
             f"ΔρY_w {rYw_n - rYw:+.4e} (= θ(δρv+δρg) {th*(dv+dg):+.4e})")


# ===================================================================== S5 保存
def s5_problem():
    """S5 の問題 (閉じた箱・拡散単独、24 セル、気相 3 種 [N2, AR, HE(初期 0)] + 蒸気 + 液 + Q)。S5・S8 で共有。"""
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
    return pb, (rho, rY, rv + g * rho, g * rho, rQ, rE), T0, V, dtime


def test_conservation(quick):
    print("\n=== S5 保存: 閉じた箱・拡散単独・解き切った物理時間更新 (補正の面 z: upwind = 採用 / mean = 現行 species_diffusion_d) ===")
    nstep = 200 if quick else 1000
    pb, args, T0, V, dtime = s5_problem()
    rho, rY, rYw0, rg0, rQ, rE = args
    rv = rYw0 - rg0
    g = rg0 / rho
    runs = ([("upwind", "float64", np.float64, "床 16ε", 16.0, False)]
            + [("upwind", "float32", np.float32, n_, f_, c_) for n_, f_, c_ in STOP_RULES_F32]
            + [("mean", "float64", np.float64, "床 16ε", 16.0, False), ("mean", "float32", np.float32, "床 6ε", 6.0, False)])
    for corr, dtn, dt, name, fe, carry in runs:
        st = make_state(rho, rY, rv + g * rho, g * rho, rQ, rE, dt)
        t0 = totals(st, V)
        M = (V / dtime).astype(dt)
        st, T, _, stats, nit_max = run_be_series(st, pb, M, dt, nstep, fe, carry, 3000, T0=T0.copy(), precond="line", corr=corr)
        cap = stats.get("cap", 0)
        t = totals(st, V)
        rvv = st["rYw"].astype(np.float64) - st["rg"].astype(np.float64)
        vmin = float(rvv.min())
        sY = st["rY"].astype(np.float64).sum(axis=0) + st["rYw"].astype(np.float64)
        sYmax = float(np.max(np.abs(sY - st["rho"].astype(np.float64)) / st["rho"].astype(np.float64)))
        rel_nz = {k: (t[k] - t0[k]) / abs(t0[k]) for k in t0 if k not in ("rE", "absE", "rv") and t0[k] != 0}
        abs_z = {k: abs(t[k] - t0[k]) for k in t0 if k not in ("rE", "absE", "rv") and t0[k] == 0}
        eE = abs(t["rE"] - t0["rE"]) / t0["absE"]
        worst = max(abs(v) for v in rel_nz.values())
        zworst = max(abs_z.values()) if abs_z else 0.0
        ok = worst <= 1e-6 and eE <= 1e-6 and zworst <= 1e-12 and vmin >= 0 and cap == 0
        line = (f"S5 [{dtn}, z={corr}, {name}] {nstep} 更新 (最大反復 {nit_max}, 上限到達 {cap} 回): 最大 |相対| {worst:.3e} "
                f"({', '.join(f'{k} {v:+.1e}' for k, v in rel_nz.items())})、ρE {eE:.3e}、初期 0 の量 {zworst:.1e}、"
                f"min ρv {vmin:.2e}、max|ΣρY−ρ|/ρ (最終) {sYmax:.1e} → {'合格' if ok else '不合格'}")
        if dtn == "float64" and corr == "upwind":
            verdict(ok, line)
        elif corr == "upwind" and fe == 6.0 and not carry:
            verdict(ok, line + " [#4b の固定条件]")
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
        st6 = {}
        st1, T1, nit, rel = be_solve(st0, pb, M, np.float64, T0=Tn, energy=variant, precond="line", maxit=200, stats=st6)
        if st6.get("cap", 0):
            verdict(False, f"S6 [{variant}] 反復上限 200 に到達 (解き切っていない)")
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


# ===================================================================== S8–S11 (§5.1 #4b)
def series_with_residual(path, st0, pb, M, dt, nstep, T0, fe=6.0, maxit=3000, corr="upwind", relax=1.0, mdot_fn=None, precond="line", **kwB):
    """path 'A' (be_solve + line_update) / 'B' (be_solve_B) で nstep 回の解き切った更新。各ステップの終状態で原方程式の残差を float64 で判定し、
    θ=0 で止まったセルがあればその残差を記録する。戻り値 dict。"""
    st, T = copy_state(st0), (None if T0 is None else np.asarray(T0, np.float64).copy())
    stats, nit_max, ratio_max, n_bad, vmin, gmin, t0_records = {}, 0, 0.0, 0, math.inf, math.inf, []
    worst_comp = {}
    for n in range(nstep):
        if mdot_fn is not None:
            pb.mdot = mdot_fn(n)
        st_n = copy_state(st)
        Tn = T
        if path == "A":
            r0, fl = orig_scales(st_n, pb, Tn, M, fe, corr=corr)
            st, T, nit, _ = be_solve(st_n, pb, M, dt, T0=Tn, maxit=maxit, floor_eps=fe, stats=stats, precond="line", corr=corr,
                                     check=("f64" if kwB.get("check64") else "state"), stop_keys=kwB.get("stop_keys_A", ("k", "v", "l", "Q", "E")))
        else:
            stB = {}
            kwB_ = {k_: v_ for k_, v_ in kwB.items() if k_ != "stop_keys_A"}
            with np.errstate(all="ignore"):
                st, T, nit, _ = be_solve_B(st_n, pb, M, dt, T0=Tn, maxit=maxit, floor_eps=fe, stats=stB, corr=corr, relax=relax, precond=precond, **kwB_)
            for k_, v_ in stB.items():
                if k_ in ("theta_min",):
                    stats[k_] = min(stats.get(k_, 1.0), v_)
                elif k_ in ("theta0_cells", "_t0"):
                    continue
                else:
                    stats[k_] = stats.get(k_, 0) + v_
            # 尺度は ρ 更新後の基点で取る (B は反復開始前に ρ と基点を更新する)
            stn2 = copy_state(st_n)
            r0, fl = orig_scales(stn2, pb, Tn, M, fe, corr=corr)
        nit_max = max(nit_max, nit)
        ratio, comp = orig_residual_ratio(st, st_n, pb, T, M, r0, fl, corr=corr)
        for k_, v_ in comp.items():
            worst_comp[k_] = max(worst_comp.get(k_, 0.0), v_)
        ratio_max = max(ratio_max, ratio)
        n_bad += int(ratio > 1.0)
        if path == "B" and stB.get("theta0_last", 0) > 0:
            s64, n64 = cast_state(st, np.float64), cast_state(st_n, np.float64)
            rr, _, _, _, _ = full_residuals(s64, n64, pb, T, np.asarray(M, np.float64), np.float64, corr=corr)
            cells = stB["theta0_cells"]
            t0_records.append((n, cells.tolist(), {k: float(np.max(np.abs(np.atleast_2d(rr[k])[..., cells]))) for k in ("w", "l")}, ratio))
        rv = st["rYw"].astype(np.float64) - st["rg"].astype(np.float64)
        vmin = min(vmin, float(rv.min()))
        gmin = min(gmin, float(st["rg"].astype(np.float64).min()))
    return {"st": st, "T": T, "stats": stats, "nit_max": nit_max, "ratio_max": ratio_max, "n_bad": n_bad, "vmin": vmin, "gmin": gmin,
            "theta0": t0_records, "comp": worst_comp}


def conservation_report(st, t0, V):
    t = totals(st, V)
    rel = {k: (t[k] - t0[k]) / abs(t0[k]) for k in t0 if k not in ("rE", "absE", "rv", "rg") and t0[k] != 0}
    zero = {k: abs(t[k] - t0[k]) for k in t0 if k not in ("rE", "absE", "rv", "rg") and t0[k] == 0}
    eE = abs(t["rE"] - t0["rE"]) / t0["absE"] if t0["absE"] > 0 else 0.0
    return rel, zero, eE, t


def test_ab_update_map(quick):
    print("\n=== S8 A/B: S5 の更新経路だけを替える (A = line_update / B = 全残差変換 + vl_limit_commit) ===")
    nstep = 200 if quick else 1000
    pb, args, T0, V, dtime = s5_problem()
    info(f"S8 固定条件: float32・補正 z = upwind・停止 = BE 残差 ≤1e-7 相対 または 床 6ε (E は流束の大きさの床も)・持ち越しなし・前処理 line (A と同じ)・"
         f"{nstep} 物理更新・各更新の上限 3000 回 (到達 = FAIL)。B の制限 dg_max = {B_LIMITS['dg_max']:g}、dT_max = {B_LIMITS['dT_max']:g} K "
         f"(forge 既定 condDgMaxStep/condDTmaxStep、実行前に固定)、緩和なし。再正規化は掛けない (全量は再正規化前)")
    for chk in (False, True):
        tag = "停止判定の残差 = float64 で再評価 (原方程式)" if chk else "停止判定の残差 = float32 (#4a と同じ)"
        print(f"--- S8 [{tag}]; 判定対象の A・B は同じ停止成分 (k, v, w, l, Q, E) で止める ---")
        res = {}
        U = ("k", "v", "w", "l", "Q", "E")
        for path, lim, keys in (("A", None, U), ("B", B_LIMITS, U), ("B (閾値なし; 診断)", {"dg_max": math.inf, "dT_max": math.inf}, U),
                                ("A (停止成分 k,v,l,Q,E = #4a; 診断)", None, ("k", "v", "l", "Q", "E")),
                                ("A (停止成分 k,w,l,Q,E; 診断)", None, ORIG_KEYS), ("B (停止成分 k,w,l,Q,E; 診断)", B_LIMITS, ORIG_KEYS)):
            st0 = make_state(*args, np.float32)
            t0 = totals(st0, V)
            M = (V / dtime).astype(np.float32)
            kw = {"check64": chk}
            if path.startswith("B"):
                kw["limits"] = lim
                kw["stop_keys"] = keys
            else:
                kw["stop_keys_A"] = keys
            out = series_with_residual(path[0], st0, pb, M, np.float32, nstep, T0, fe=6.0, **kw)
            rel, zero, eE, t = conservation_report(out["st"], t0, V)
            worst = max(abs(v) for v in rel.values())
            zw = max(zero.values()) if zero else 0.0
            st = out["st"]
            sY = st["rY"].astype(np.float64).sum(axis=0) + st["rYw"].astype(np.float64)
            sYm = float(np.max(np.abs(sY - st["rho"].astype(np.float64)) / st["rho"].astype(np.float64)))
            cap = out["stats"].get("cap", 0)
            ok = worst <= 1e-6 and eE <= 1e-6 and zw <= 1e-12 and cap == 0 and out["vmin"] >= 0 and out["n_bad"] == 0
            line = (f"S8 {path} [{'f64 停止' if chk else 'f32 停止'}]: 保存 最大 |相対| {worst:.3e} ({', '.join(f'{k} {v:+.1e}' for k, v in rel.items())})、ρE {eE:.2e}、"
                    f"初期 0 の量 {zw:.1e}、原方程式残差 (float64) の最大比 {out['ratio_max']:.3f} (>1 のステップ {out['n_bad']}; 成分別 "
                    f"{', '.join(f'{k} {v:.2f}' for k, v in out['comp'].items())})、上限到達 {cap}、最大反復 {out['nit_max']}、min ρv {out['vmin']:.3e}、"
                    f"|ΣρY−ρ|/ρ (最終) {sYm:.1e}")
            if path.startswith("B"):
                s_ = out["stats"]
                line += (f"、θ 最小 {s_.get('theta_min', 1.0):.3f}・θ<1 のセル·反復 {s_.get('theta_lt1', 0)}、保留した増分 Σ (1−θ)|δ| 蒸気 {s_.get('limiter_v', 0):.2e} / "
                         f"液 {s_.get('limiter_g', 0):.2e}、状態の補正 Q_cut {s_.get('Q_cut', 0):.2e} / v_round {s_.get('v_round', 0):.2e}、θ=0 で止まったステップ {len(out['theta0'])}")
                for rec in out["theta0"][:3]:
                    line += f" [θ=0 ステップ {rec[0]} セル {rec[1]}: |r_w| {rec[2]['w']:.1e} |r_l| {rec[2]['l']:.1e}]"
            if path in ("A", "B"):
                verdict(ok, line)
                res[path] = ok
            else:
                info(line)
        a, b = res["A"], res["B"]
        if a and b:
            verdict(True, f"S8 判定規則 [{tag}]: A・B とも合格 → この試験では更新写像による保存上の懸念を除外")
        elif a and not b:
            verdict(False, f"S8 判定規則 [{tag}]: A 合格・B 不合格 → 更新写像 (全残差変換・limiter・commit) の問題")
        elif not a and not b:
            verdict(False, f"S8 判定規則 [{tag}]: A・B とも不合格 → 停止則の問題へ戻る")
        else:
            verdict(False, f"S8 判定規則 [{tag}]: A 不合格・B 合格 (規則に無い組合せ) → 上位へ")


def source_problem(N=12):
    th = Thermo([nasa_species("N2")], nasa_species("H2O"), None)
    th.lat = Latent(th.v)
    dx = 1e-3
    x = (np.arange(N) + 0.5) * dx
    rho = np.full(N, 0.3)
    T0 = np.where(x < N * dx / 2, 230.0, 260.0)
    left = x < N * dx / 2
    g = np.where(left, 0.0, 1.0e-3)
    yv = np.where(left, 1.0e-3, 5.0e-4)
    rYw = rho * (yv + g)
    rg = rho * g
    rY = (rho - rYw)[None, :]
    rQ = np.vstack([rg * 1e3, rg * 1e6, rg * 1e9])
    rE = rho * th.e_mix(T0, rY / rho, rYw / rho, g)
    pb = Problem(N, np.full(N - 1, 1 / dx), np.full(N - 1, 1e-4), np.array([[0, 2.2e-5], [2.2e-5, 0]]), np.array([0.028, 0.018]),
                 np.full(N, dx), molc=1.0, thermo=th)
    pb.source = GrowthSource(tau=1e-3, rg_ref=1e-6, Rw=th.v.R)
    xf = (np.arange(N - 1) + 1) * dx
    m0 = 0.0072
    pb.mdot = None
    mdot_fn = lambda n: m0 * np.sin(math.pi * xf / (N * dx)) * math.cos(2 * math.pi * n / 50)
    return pb, (rho, rY, rYw, rg, rQ, rE), T0, mdot_fn, 2e-3


def test_source_iteration(quick):
    print("\n=== S9 B の反復試験: 実際の相変化ソース (成長則 GrowthSource) + ρ の更新 (移流 ṁ, φ_N δρ) + ソース Jacobian ===")
    nstep = 50 if quick else 200
    info("S9 条件: 12 セル閉じた箱、N2 + H2O (NASA-9)、左半分 230 K 過飽和 (Y_v 1e-3, g 0)・右半分 260 K 未飽和 (g 1e-3, Y_v 5e-4)、"
         "乱流 Γ 1e-4・分子 D 2.2e-5、ṁ_f = 7.2e-3·sin(πx/L)·cos(2πn/50) kg/m²/s (ρ が ±1 割程度往復)、τ = 1e-3 s、Δt = 2e-3 s、"
         f"{nstep} 物理更新、停止 6ε、上限 3000、B_LIMITS {B_LIMITS}。Q はソースなし (移流・拡散だけ)")
    variants = [(dtn, dt, relax, 1.0, True, nstep, True, False) for dtn, dt in (("float64", np.float64), ("float32", np.float32)) for relax in (1.0, 0.5)]
    variants += [("float32", np.float32, relax, 1.0, True, nstep, False, True) for relax in (1.0, 0.5)]
    variants += [("float32", np.float32, 1.0, 0.5, True, nstep, False, False), ("float64", np.float64, 1.0, 1.0, False, 3, False, False)]
    for dtn, dt, relax, frac, sblock, nst, judged, fdiag in variants:
        if True:
            pb, args, T0, mdot_fn, dtime = source_problem()
            st0 = make_state(*args, dt)
            t0 = totals(st0, pb.V)
            rho0 = float(np.sum(st0["rho"].astype(np.float64) * pb.V))
            M = (pb.V / dtime).astype(dt)
            out = series_with_residual("B", st0, pb, M, dt, nst, T0, fe=6.0, relax=relax, mdot_fn=mdot_fn, frac=frac, src_block=sblock, floor_diag=fdiag)
            st = out["st"]
            t = totals(st, pb.V)
            cons = {k: (t[k] - t0[k]) / abs(t0[k]) for k in ("rYw", "rY0", "rE")}
            cons["rho"] = (float(np.sum(st["rho"].astype(np.float64) * pb.V)) - rho0) / rho0
            worst = max(abs(v) for v in cons.values())
            s = out["stats"]
            cap = s.get("cap", 0)
            ok = cap == 0 and worst <= 1e-6 and out["vmin"] >= 0 and out["gmin"] >= 0 and out["n_bad"] == 0
            g_end = st["rg"].astype(np.float64) / st["rho"].astype(np.float64)
            line = (f"S9 [{dtn}, 緩和 {relax}, 前処理 {'蒸気・液 2×2 連成' if sblock else '対角のみ'}, θ_vg の距離 ×{frac}, 床 {'6ε·(M+対角)' if fdiag else '6ε·M'}, {nst} 更新] 上限到達 {cap}、最大反復 {out['nit_max']}、原方程式残差の最大比 {out['ratio_max']:.3f} (>1 のステップ {out['n_bad']})、"
                    f"保存 {', '.join(f'{k} {v:+.1e}' for k, v in cons.items())}、min ρv {out['vmin']:.3e}、min ρg {out['gmin']:.3e}、"
                    f"θ 最小 {s.get('theta_min', 1.0):.3f}・θ<1 のセル·反復 {s.get('theta_lt1', 0)}、Q_cut {s.get('Q_cut', 0):.1e}、v_round {s.get('v_round', 0):.1e}、"
                    f"θ=0 のステップ {len(out['theta0'])}、終状態 g {np.nanmin(g_end):.2e}–{np.nanmax(g_end):.2e}・T {np.nanmin(out['T']):.1f}–{np.nanmax(out['T']):.1f} K")
            if judged:
                verdict(ok, line)
            else:
                info(line + " (診断; 原方程式残差の比は 6ε·M の床に対して)")
            for rec in out["theta0"][:3]:
                info(f"S9 [{dtn}, 緩和 {relax}, ×{frac}] θ=0 で止まったセル (ステップ {rec[0]}, セル {rec[1]}): 原方程式残差 |r_w| {rec[2]['w']:.2e}・|r_l| {rec[2]['l']:.2e}、全体の比 {rec[3]:.3f}")


class ConstSource:
    """θ=0 停止の記録用: 一定の凝縮要求 S [kg/m³/s] (Jacobian 0)。"""

    def __init__(self, S):
        self.S = S

    def __call__(self, st64, T, thermo):
        n = st64["rho"].shape[0]
        return np.full(n, self.S), np.zeros(n), np.zeros(n)


def test_theta0_record(quick):
    print("\n=== S10 θ=0 で止まった状態の原方程式残差 (固定点との区別) ===")
    for dtn, dt in (("float64", np.float64), ("float32", np.float32)):
        pb = Problem(2, [1.0], [0.0], np.array([[0, 1e-5], [1e-5, 0]]), np.array([0.028, 0.018]), np.ones(2), molc=0.0)
        pb.source = ConstSource(0.01)
        st0 = make_state([1, 1], [[0.9, 0.9]], [0.1, 0.1], [0.1, 0.1], [[0.3, 0.3], [0.2, 0.2], [0.1, 0.1]], [0, 0], dt)
        stats = {}
        st, _, nit, rel = be_solve_B(st0, pb, dt(1.0), dt, maxit=200, stats=stats, precond="point")
        r, _, _, _, _ = full_residuals(cast_state(st, np.float64), cast_state(st0, np.float64), pb, None, 1.0, np.float64)
        moved = float(np.max(np.abs(st["rg"].astype(np.float64) - st0["rg"].astype(np.float64))))
        verdict(stats.get("cap", 0) > 0 and moved == 0 and np.max(np.abs(r["l"])) > 0,
                f"S10 [{dtn}] 蒸気 0・一定の凝縮要求 S = +0.01: θ = 0 で状態は動かない (最大変化 {moved:.1e}) が、原方程式残差 |r_l| = "
                f"{np.max(np.abs(r['l'])):.3e}・|r_w − r_l| = {np.max(np.abs(r['w'] - r['l'])):.3e} は 0 でない → 上限到達 {stats.get('cap', 0)} (反復 {nit}) "
                f"で FAIL 扱い。**動かない状態 ≠ 固定点** を停止判定が区別することの確認")


def test_counterexample_B(quick):
    print("\n=== S11 反例 (両セル g=0, D 比 100, 風上 z) を新経路 B で: 緩和 1 / 0.5 (緩和は前処理の後・limiter の前) ===")
    for dtn, dt in (("float64", np.float64), ("float32", np.float32)):
        for relax in (1.0, 0.5):
            st0 = make_state([1, 1], np.array([[0.0, 0.99], [1.0, 0.0]]), [0.0, 0.01], [0.0, 0.0], np.zeros((1, 2)), [0, 0], dt)
            pb = Problem(2, [1.0], [0.0], np.array([[0, 10, 10], [10, 0, 0.1], [10, 0.1, 0]]), np.array([0.028, 0.028, 0.028]), np.ones(2))
            stats = {}
            with np.errstate(all="ignore"):
                st, _, nit, rel = be_solve_B(st0, pb, dt(1.0), dt, maxit=5000, stats=stats, precond="point", relax=relax, floor_eps=6.0)
            rv = (st["rYw"] - st["rg"]).astype(np.float64)
            cap = stats.get("cap", 0)
            line = (f"S11 [{dtn}, 緩和 {relax}] {'上限 5000 に到達 (収束しない)' if cap else f'{nit} 回で収束'}、ρv = [{rv[0]:+.6e}, {rv[1]:+.6e}] "
                    f"(Newton 解 [+1.195754e-03, +8.804246e-03])、θ 最小 {stats.get('theta_min', 1.0):.3f}・θ<1 {stats.get('theta_lt1', 0)}")
            if relax == 0.5:
                verdict(cap == 0 and rv.min() >= 0 and abs(rv[0] - 1.195754e-3) < 1e-6, line)
            else:
                info(line)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quick", action="store_true", help="反復数を減らす (開発用; 判定は正式でない)")
    ap.add_argument("--only", default="", help="S1..S6 をカンマ区切りで")
    args = ap.parse_args()
    tests = {"S1": test_structure, "S2": test_face_identities, "S3": test_three_cell,
             "S4": test_nonnegativity, "S5": test_conservation, "S6": test_energy, "S7": test_fixed_point,
             "S8": test_ab_update_map, "S9": test_source_iteration, "S10": test_theta0_record, "S11": test_counterexample_B}
    sel = [s.strip() for s in args.only.split(",") if s.strip()] or list(tests)
    for name in sel:
        tests[name](args.quick)
    print(f"\n{'ALL PASS' if g_fail == 0 else f'{g_fail} FAIL'}" + (" (--quick: 正式判定ではない)" if args.quick else ""))
    return 0 if g_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
