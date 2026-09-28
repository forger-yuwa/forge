#!/usr/bin/env python3
r"""case/63 Graetz の温度再現 A/B (plan `boundary-cht-axisymmetric-graetz.md` §5.1 #6g、事前登録 2026-09-29)。

対照 run (ΔT=0) の保存場の ρ・u・v・p を固定し、forge と独立に定常温度場を解く:

    ρ c_p (u T_x + v T_r) = k [T_xx + (1/r)(r T_r)_r] + α (Φ + u p_x + v p_r)
    Φ = τ_xx u_x + τ_rr v_r + τ_θθ v/r + τ_xr (u_r + v_x)   (τ = μ(∇u+∇uᵀ) − (2/3)μ(∇·u)I、軸対称)

A: α=0 / B: α=1。熱 BC は共通: 入口列 = 保存場の T (Dirichlet)、加熱区間の壁節点 = 保存場の壁温 (Dirichlet)、
上流・下流の壁は断熱 (∂T/∂r=0)、軸は対称、出口列は零勾配。**x=0 の温度は与えない**。
離散化: 構造格子の節点で有限差分。対流は x 方向 2 次風上 (入口の次の列は 1 次)、r 方向は中心差分。拡散は 2 次。

    python3 temp_reproduce.py selftest                 # 合成 2 試験 (実データの前に不確かさを確認)
    python3 temp_reproduce.py run <run> [--last 25]    # A/B と判定
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import h5py
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

import graetz_common as gc

TOL = 0.005          # 登録: 分布・幅の再現 [K]
UNC = TOL / 3        # 登録: 後処理離散化の不確かさ [K]


class Grid:
    def __init__(self, xyz):
        x = np.round(xyz[:, 0], 10)
        self.xs = np.unique(x)
        cols = [np.where(x == xv)[0] for xv in self.xs]
        self.idx = np.array([c[np.argsort(xyz[c, 1])] for c in cols])     # (nx, nr+1) 節点番号
        self.r = xyz[self.idx[0], 1]
        self.nx, self.ny = self.idx.shape
        self.h = self.r[1] - self.r[0]


def grads(G: Grid, f):
    """節点場 f (節点番号順) の ∂/∂x, ∂/∂r を (nx, ny) で返す (2 次)。"""
    F = f[G.idx]
    fx = np.gradient(F, G.xs, axis=0, edge_order=2)
    fr = np.gradient(F, G.r, axis=1, edge_order=2)
    return F, fx, fr


def source(G: Grid, ro, u, v, p):
    mu = gc.MU
    U, ux, ur = grads(G, u)
    V, vx, vr = grads(G, v)
    _, px, pr = grads(G, p)
    R = np.broadcast_to(G.r, U.shape)
    v_r = np.where(R > 0, V / np.where(R > 0, R, 1.0), vr)                  # 軸上は極限 ∂v/∂r
    div = ux + vr + v_r
    txx = 2 * mu * ux - 2 / 3 * mu * div
    trr = 2 * mu * vr - 2 / 3 * mu * div
    ttt = 2 * mu * v_r - 2 / 3 * mu * div
    txr = mu * (ur + vx)
    phi = txx * ux + trr * vr + ttt * v_r + txr * (ur + vx)
    work = U * px + V * pr
    return phi, work


def solve_T(G: Grid, ro, u, v, S, T_inlet, wall_dir):
    """定常温度場を直接解く。wall_dir[i] が有限なら列 i の壁節点を Dirichlet。戻り値 (nx, ny)。"""
    nx, ny, h = G.nx, G.ny, G.h
    RO, U, V = ro[G.idx], u[G.idx], v[G.idx]
    k, cp = gc.K_F, gc.CP
    N = nx * ny
    idx = lambda i, j: i * ny + j
    rows, cols, vals = [], [], []
    b = np.zeros(N)

    def add(p, q, w):
        rows.append(p); cols.append(q); vals.append(w)

    xs = G.xs
    for i in range(nx):
        for j in range(ny):
            p = idx(i, j)
            if i == 0:
                add(p, p, 1.0); b[p] = T_inlet[j]; continue
            if j == ny - 1 and np.isfinite(wall_dir[i]):
                add(p, p, 1.0); b[p] = wall_dir[i]; continue
            if i == nx - 1:
                add(p, p, 1.0); add(p, idx(i - 1, j), -1.0); continue
            a = {}
            # 対流 ρ c_p u T_x (2 次風上、u ≥ 0 を仮定。入口の次の列は 1 次)
            c = RO[i, j] * cp * U[i, j]
            if i >= 2:
                d1, d2 = xs[i] - xs[i - 1], xs[i] - xs[i - 2]
                w0 = (d1 + d2) / (d1 * d2); w1 = -d2 / (d1 * (d2 - d1)); w2 = d1 / (d2 * (d2 - d1))
                for ii, w in ((i, w0), (i - 1, w1), (i - 2, w2)):
                    a[idx(ii, j)] = a.get(idx(ii, j), 0.0) + c * w
            else:
                dx = xs[i] - xs[i - 1]
                a[p] = a.get(p, 0.0) + c / dx; a[idx(i - 1, j)] = a.get(idx(i - 1, j), 0.0) - c / dx
            # 拡散 −k T_xx (非等間隔 3 点)
            dm, dp = xs[i] - xs[i - 1], xs[i + 1] - xs[i]
            cm, cpp = 2 / (dm * (dm + dp)), 2 / (dp * (dm + dp))
            a[idx(i - 1, j)] = a.get(idx(i - 1, j), 0.0) - k * cm
            a[idx(i + 1, j)] = a.get(idx(i + 1, j), 0.0) - k * cpp
            a[p] = a.get(p, 0.0) + k * (cm + cpp)
            # 半径方向: 対流 ρ c_p v T_r (中心) と拡散 −k (T_rr + T_r/r)
            if j == 0:                                   # 軸: (1/r)(r T_r)_r → 2 T_rr、T_r = 0
                a[p] = a.get(p, 0.0) + k * 4 / h ** 2
                a[idx(i, 1)] = a.get(idx(i, 1), 0.0) - k * 4 / h ** 2
            elif j == ny - 1:                            # 断熱壁: ゴースト T_{N+1} = T_{N−1}、T_r = 0
                a[p] = a.get(p, 0.0) + k * 2 / h ** 2
                a[idx(i, j - 1)] = a.get(idx(i, j - 1), 0.0) - k * 2 / h ** 2
            else:
                rj = G.r[j]
                cv = RO[i, j] * cp * V[i, j] / (2 * h)
                a[idx(i, j + 1)] = a.get(idx(i, j + 1), 0.0) + cv - k * (1 / h ** 2 + 1 / (2 * h * rj))
                a[idx(i, j - 1)] = a.get(idx(i, j - 1), 0.0) - cv - k * (1 / h ** 2 - 1 / (2 * h * rj))
                a[p] = a.get(p, 0.0) + k * 2 / h ** 2
            for q, w in a.items():
                add(p, q, w)
            b[p] = S[i, j]
    A = sp.csr_matrix((vals, (rows, cols)), shape=(N, N))
    return spla.spsolve(A.tocsc(), b).reshape(nx, ny)


# ------------------------------------------------------------------ 合成試験
def selftest():
    here = Path(__file__).resolve().parent
    with h5py.File(here / "mesh" / "graetz_r32.h5", "r") as m:
        xyz = np.asarray(m["MESH/COORD"][:], float).reshape(-1, 3)
    G = Grid(xyz)
    n = len(xyz)
    r = xyz[:, 1]
    umax = 2 * gc.U_M
    u = umax * np.clip(1 - (r / gc.R) ** 2, 0, None)
    v = np.zeros(n); ro = np.full(n, gc.RHO); p = np.full(n, gc.P_OUT)
    ok = True
    # (1) 散逸だけの十分発達解: T − T_w = μ U_max² (1 − (r/R)⁴)/(4k)、壁は全域 Dirichlet
    phi, work = source(G, ro, u, v, p)
    Tex_r = gc.T_IN + gc.MU * umax ** 2 * (1 - (G.r / gc.R) ** 4) / (4 * gc.K_F)
    T = solve_T(G, ro, u, v, phi, Tex_r, np.full(G.nx, gc.T_IN))
    err1 = np.abs(T - Tex_r[None, :]).max()
    amp = Tex_r.max() - gc.T_IN
    print(f"(1) 散逸の十分発達解 (振幅 {amp:.4f} K): max|T − T_ex| = {err1:.2e} K  (許容 {UNC:.2e})  {'ok' if err1 <= UNC else 'NG'}")
    ok &= err1 <= UNC
    # (2) Graetz (α=0、放物速度、上流断熱・加熱区間 T_w = T_in + dT、下流断熱) を ellip (Pe 720) と比較
    import graetz_ref
    dT = 0.3
    Pe = gc.PE
    xu, xd = gc.L_UP / (gc.D * Pe), gc.L_DOWN / (gc.D * Pe)
    xe, nue, tbe, xie, the = graetz_ref.solve_ellip(120, Pe, xu, 0.12, xd, 2400, 400, 400, return_field=True)
    wall = np.where((G.xs >= -1e-9) & (G.xs <= gc.L_HEAT + 1e-7), gc.T_IN + dT, np.nan)
    T2 = solve_T(G, ro, u, v, np.zeros((G.nx, G.ny)), np.full(G.ny, gc.T_IN), wall)
    worst = 0.0
    # x=0 は壁温の段差の角 (特異点) で、基準の ellip 自身が格子を倍にして θ が 0.693→0.702 と動く (未収束) ので
    # 判定に使わない (2026-09-29、実データの前に変更: plan §5.1 #6g)。参考として出すだけ
    i0 = int(np.argmin(np.abs(G.xs)))
    ie0 = int(np.argmin(np.abs(xe)))
    e0 = np.abs(T2[i0, :-1] - (gc.T_IN + dT - dT * np.interp(G.r / gc.R, xie, the[ie0])[:-1])).max()
    print(f"(2) 参考 x=0 (壁温段差 {dT} K の角、基準未収束): max|T − T_ref| = {e0:.2e} K → 段差あたり {e0/dT:.3f}")
    for xmm in (1.0, 4.3, 14.4):
        i = int(np.argmin(np.abs(G.xs - xmm * 1e-3)))
        xp = G.xs[i] / (gc.D * Pe)
        ie = int(np.argmin(np.abs(xe - xp)))
        Tref = gc.T_IN + dT - dT * np.interp(G.r / gc.R, xie, the[ie])
        e = np.abs(T2[i, :-1] - Tref[:-1]).max()
        worst = max(worst, e)
        print(f"(2) Graetz x={G.xs[i]*1e3:7.3f} mm (ellip x⁺ {xe[ie]:.2e} vs {xp:.2e}): max|T − T_ref| (壁除く) = {e:.2e} K")
    print(f"    最大 {worst:.2e} K (ΔT {dT} K、許容 {UNC:.2e})  {'ok' if worst <= UNC else 'NG'}")
    ok &= worst <= UNC
    print("VERDICT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


# ------------------------------------------------------------------ 実データ
def run_ab(run: Path, last: int):
    with h5py.File(run / "mesh.h5", "r") as m:
        xyz = np.asarray(m["MESH/COORD"][:], float).reshape(-1, 3)
    G = Grid(xyz)
    steps = sorted(int(p.name[4:-3]) for p in run.glob("res_[0-9]*.h5"))[-last:]
    rows = []
    for st in steps:
        with h5py.File(run / f"res_{st}.h5", "r") as h:
            V = h["VALUE"]
            ro, u, v, p, T = (np.asarray(V[k][:], float) for k in ("ro", "Ux", "Uy", "P", "T"))
        phi, work = source(G, ro, u, v, p)
        Ts = T[G.idx]
        heat = (G.xs >= -1e-7) & (G.xs <= gc.L_HEAT + 1e-7)
        wall = np.where(heat, Ts[:, -1], np.nan)
        TA = solve_T(G, ro, u, v, np.zeros_like(phi), Ts[0], wall)
        TB = solve_T(G, ro, u, v, phi + work, Ts[0], wall)
        i0 = int(np.argmin(np.abs(G.xs)))
        s, a, b = Ts[i0, :-1], TA[i0, :-1], TB[i0, :-1]
        step_corner = abs(Ts[i0 - 1, -1] - Ts[i0, -1])     # 上流の断熱壁と加熱壁の角の壁温段差
        rows.append([st, np.ptp(s), np.ptp(a), np.ptp(b), np.abs(a - s).max(), np.abs(b - s).max(),
                     s[0], b[0], a[0], s[-1], b[-1], a[-1], step_corner])
    R = np.array(rows)
    out = run / "temp_reproduce_series.csv"
    np.savetxt(out, R, delimiter=",", comments="", fmt=["%d"] + ["%.10e"] * 12,
               header="step,width_saved,width_A,width_B,maxdiff_A,maxdiff_B,Taxis_saved,Taxis_B,Taxis_A,T1_saved,T1_B,T1_A,corner_step")
    L = R[-1]
    print(f"=== 温度再現 A/B {run}  (step {int(L[0])}、x=0 の内部節点 {G.ny-1} 点)")
    print(f"  保存場: 幅 {L[1]:.5f} K、軸 {L[6]:.4f} K、壁の隣 {L[9]:.4f} K")
    print(f"  A (α=0): 幅 {L[2]:.5f} K、分布の最大差 {L[4]:.5f} K (軸 {L[8]:.4f}、壁の隣 {L[11]:.4f})")
    print(f"  B (α=1): 幅 {L[3]:.5f} K、分布の最大差 {L[5]:.5f} K (軸 {L[7]:.4f}、壁の隣 {L[10]:.4f})")
    okB = L[5] <= TOL and abs(L[3] - L[1]) <= TOL
    okA = L[4] <= TOL and abs(L[2] - L[1]) <= TOL
    print(f"  B 再現 (分布・幅とも ≤ {TOL} K): {okB}   A 再現: {okA}")
    print(f"  角の壁温段差 {L[12]:.4f} K → 合成試験の比 0.048 K/K から角の効果の見積もり ≈ {0.048*L[12]:.1e} K")
    print(f"  参考: 末尾 {len(R)} 枚で B の最大差の変動 {np.ptp(R[:,5]):.2e} K、保存場の幅の変動 {np.ptp(R[:,1]):.2e} K")
    if okB and not okA:
        v = "第 1 仮説を支持 (A で再現せず B で再現)"
    elif not okB:
        v = "第 1 仮説を棄却 (B でも再現しない) — 離散化・境界・抽出を候補に戻す"
    else:
        v = "散逸・圧力仕事を主因とする説明は支持しない (A でも再現)"
    print(f"VERDICT: {v}   系列 → {out}")
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mode", choices=["selftest", "run"])
    ap.add_argument("run", nargs="?")
    ap.add_argument("--last", type=int, default=25)
    a = ap.parse_args()
    if a.mode == "selftest":
        return selftest()
    return run_ab(Path(a.run), a.last)


if __name__ == "__main__":
    sys.exit(main())
