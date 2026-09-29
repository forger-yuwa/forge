#!/usr/bin/env python3
r"""流速場を与えた定常の共役伝熱を、流体と固体を一つの系として直接解く独立参照解 (forge と独立の実装)。

plan [`boundary-cht-conjugate-benchmarks.md`](../../plans/active/boundary-cht-conjugate-benchmarks.md) §4.1。
case/64 (A: 厚肉管、軸対称) と case/65 (C: 共役平板、平面) で共用する。

離散化: テンソル格子 (x_i, y_j) の**節点中心の有限体積** (双対セル)。各主セル (i,j)–(i+1,j+1) に材料 (流体 / 固体 / 空) を割り当て、
熱伝導の面コンダクタンスは辺に接する 2 つの主セルの寄与の和 (各主セルが半分ずつ: k_cell × 相手方向の半幅 / 辺長)。
軸対称では面積・体積に r を掛ける (y = r)。界面の節点は双対セルの半分が流体・半分が固体になり、温度と熱流束の連続が自然に入る。

    流体:  ρ c_p (u T_x + v T_y) = ∇·(k_f ∇T) + α (Φ + u p_x + v p_y)
    固体:  0 = ∇·(k_s ∇T)

対流は保存形 + 発散補正 (双対セルの面の質量流束 × 面の温度 − 節点の温度)。面の温度は x 面が 2 次風上、y 面が中心。
境界: 入口列の流体だけの節点は Dirichlet T_in。出口列は零勾配 (T_N = T_{N−1})。軸 (y=0、軸対称) と空セルに面する辺は断熱。
固体外面の Robin は `robin(x_mid, y) -> (h, T_c)` を返す関数で外周の辺ごとに与える (None なら断熱)。
壁 Dirichlet (自己検査用) は `wall_dirichlet(x) -> T or None`。

    python3 conjugate_ref.py selftest
"""
from __future__ import annotations

import sys

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

FLUID, SOLID, VOID = 0, 1, 2


class Problem:
    def __init__(self, xs, ys, mat, k_f, k_s, cp, axisym, T_in, ro=None, u=None, v=None, source=None,
                 robin=None, wall_dirichlet=None, j_wall=None, top_dirichlet=None, source_solid=None,
                 solid_axial=1.0, exact_boundary=None):
        self.xs, self.ys = np.asarray(xs, float), np.asarray(ys, float)
        self.nx, self.ny = len(xs), len(ys)
        self.mat = np.asarray(mat)                          # (nx-1, ny-1)
        self.k_f, self.k_s, self.cp, self.axisym, self.T_in = k_f, k_s, cp, axisym, T_in
        z = np.zeros((self.nx, self.ny))
        self.ro = z + 1 if ro is None else ro
        self.u = z if u is None else u
        self.v = z if v is None else v
        self.S = z if source is None else source
        self.robin, self.wall_dirichlet, self.j_wall = robin, wall_dirichlet, j_wall
        self.top_dirichlet = top_dirichlet
        self.Ss = z if source_solid is None else source_solid       # 固体の体積源 (製造解の試験用)
        self.solid_axial = solid_axial                              # 固体の x 方向伝導の倍率 (0 で軸方向伝導なし: 診断用)
        self.exact_boundary = exact_boundary                        # f(x, y): 外周の全節点を Dirichlet (製造解の試験用)

    def kcell(self, i, j):
        if i < 0 or j < 0 or i >= self.nx - 1 or j >= self.ny - 1:
            return None
        m = self.mat[i, j]
        return None if m == VOID else (self.k_f if m == FLUID else self.k_s)

    def rw(self, y):
        return y if self.axisym else 1.0

    def solve(self):
        nx, ny, xs, ys = self.nx, self.ny, self.xs, self.ys
        idx = lambda i, j: i * ny + j
        N = nx * ny
        rows, cols, vals = [], [], []
        b = np.zeros(N)
        active = np.zeros((nx, ny), bool)
        fluidnode = np.zeros((nx, ny), bool)
        for i in range(nx):
            for j in range(ny):
                for di in (-1, 0):
                    for dj in (-1, 0):
                        c = self.kcell(i + di, j + dj)
                        if c is not None:
                            active[i, j] = True
                            if self.mat[i + di, j + dj] == FLUID:
                                fluidnode[i, j] = True
        self.active, self.fluidnode = active, fluidnode

        def add(p, q, w):
            rows.append(p); cols.append(q); vals.append(w)

        for i in range(nx):
            for j in range(ny):
                p = idx(i, j)
                if not active[i, j]:
                    add(p, p, 1.0); b[p] = self.T_in; continue
                # 入口: 流体だけに囲まれた節点は Dirichlet
                if self.exact_boundary is not None and (i in (0, nx - 1) or j == ny - 1 or (j == 0 and not self.axisym)):
                    add(p, p, 1.0); b[p] = self.exact_boundary(xs[i], ys[j]); continue
                if i == 0 and fluidnode[i, j] and all(
                        self.kcell(i + di, j + dj) in (None, self.k_f) and
                        (self.kcell(i + di, j + dj) is None or self.mat[i + di, j + dj] == FLUID)
                        for di in (-1, 0) for dj in (-1, 0)):
                    add(p, p, 1.0); b[p] = self.T_in; continue
                if i == nx - 1:
                    add(p, p, 1.0); add(p, idx(i - 1, j), -1.0); continue
                if self.wall_dirichlet is not None and j == self.j_wall:
                    tw = self.wall_dirichlet(xs[i])
                    if tw is not None:
                        add(p, p, 1.0); b[p] = tw; continue
                if self.top_dirichlet is not None and j == ny - 1:
                    add(p, p, 1.0); b[p] = self.top_dirichlet; continue
                a = {}
                # 伝導: 4 方向の辺。辺に接する主セルの寄与を半分ずつ足す
                for (di, dj) in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                    ii, jj = i + di, j + dj
                    if ii < 0 or jj < 0 or ii >= nx or jj >= ny:
                        continue
                    G = 0.0
                    if di != 0:                                  # x 方向の辺: 主セルは (min(i,ii), j-1) と (min(i,ii), j)
                        ci = min(i, ii); dx = abs(xs[ii] - xs[i])
                        for cj, ylo, yhi in ((j - 1, None, None), (j, None, None)):
                            k = self.kcell(ci, cj)
                            if k is None:
                                continue
                            y0, y1 = (0.5 * (ys[j - 1] + ys[j]), ys[j]) if cj == j - 1 else (ys[j], 0.5 * (ys[j] + ys[j + 1]))
                            if self.mat[ci, cj] == SOLID:
                                k = k * self.solid_axial
                            G += k * abs(y1 - y0) * self.rw(0.5 * (y0 + y1)) / dx
                    else:                                        # y 方向の辺
                        cj = min(j, jj); dy = abs(ys[jj] - ys[j]); ymid = 0.5 * (ys[j] + ys[jj])
                        for ci in (i - 1, i):
                            k = self.kcell(ci, cj)
                            if k is None:
                                continue
                            x0, x1 = (0.5 * (xs[i - 1] + xs[i]), xs[i]) if ci == i - 1 else (xs[i], 0.5 * (xs[i] + xs[i + 1]))
                            G += k * abs(x1 - x0) * self.rw(ymid) / dy
                    if G > 0:
                        a[idx(ii, j + dj if di == 0 else j)] = a.get(idx(ii, j + dj if di == 0 else j), 0.0) - G
                        a[p] = a.get(p, 0.0) + G
                # 双対セルの流体部分の体積 (対流・源項用)
                Vf = 0.0
                for di in (-1, 0):
                    for dj in (-1, 0):
                        ci, cj = i + di, j + dj
                        if 0 <= ci < nx - 1 and 0 <= cj < ny - 1 and self.mat[ci, cj] == FLUID:
                            xa, xb = (0.5 * (xs[i - 1] + xs[i]), xs[i]) if di == -1 else (xs[i], 0.5 * (xs[i] + xs[i + 1]))
                            ya, yb = (0.5 * (ys[j - 1] + ys[j]), ys[j]) if dj == -1 else (ys[j], 0.5 * (ys[j] + ys[j + 1]))
                            Vf += abs(xb - xa) * abs(yb - ya) * self.rw(0.5 * (ya + yb))
                if Vf > 0:
                    # 対流 (保存形 + 発散補正): Σ_面 c_p m_f (T_f − T_p)。m_f は双対セルの面を通る質量流束 (流体部分)、
                    # T_f は x 面が 2 次風上 (直線外挿)、y 面が中心。発散 0 の流れでは Σ c_p m_f T_f と同値で、熱収支が離散的に閉じる。
                    cp = self.cp
                    def afluid(ci_list):
                        """x 面の流体部分の面積 ∫ r dr (軸対称) / 高さ (平面)。軸の双対セルでも有限 (Δr²/8)。
                        初版は高さ × r(y_j) で、軸上 (r=0) の軸方向対流が消えていた (2026-09-29 codex plan レビュー M1)。"""
                        As = 0.0
                        for cj in (j - 1, j):
                            if 0 <= cj < ny - 1 and any(0 <= ci < nx - 1 and self.mat[ci, cj] == FLUID for ci in ci_list):
                                y0, y1 = (0.5 * (ys[j - 1] + ys[j]), ys[j]) if cj == j - 1 else (ys[j], 0.5 * (ys[j] + ys[j + 1]))
                                As += abs(y1 - y0) * self.rw(0.5 * (y0 + y1))
                        return As
                    for side in (1, -1):                               # x 面 (i+1/2 と i−1/2)
                        ii = i + side
                        if not (0 <= ii < nx):
                            continue
                        ci = min(i, ii)
                        Af = afluid([ci])
                        if Af == 0.0:
                            continue
                        m = 0.5 * (self.ro[i, j] * self.u[i, j] + self.ro[ii, j] * self.u[ii, j]) * Af * side
                        if m == 0.0:
                            continue
                        xf = 0.5 * (xs[i] + xs[ii])
                        up, dn = (i, ii) if m > 0 else (ii, i)             # 風上側の節点
                        upup = up - (dn - up)                              # さらに風上
                        wts = {up: 1.0}
                        if 0 <= upup < nx:
                            g = (xf - xs[up]) / (xs[up] - xs[upup])
                            wts = {up: 1.0 + g, upup: -g}
                        for q, w in wts.items():
                            a[idx(q, j)] = a.get(idx(q, j), 0.0) + cp * m * w
                        a[p] = a.get(p, 0.0) - cp * m
                    for side in (1, -1):                               # y 面 (中心)
                        jj = j + side
                        if not (0 <= jj < ny):
                            continue
                        cj = min(j, jj); W = 0.0
                        for ci in (i - 1, i):
                            if 0 <= ci < nx - 1 and 0 <= cj < ny - 1 and self.mat[ci, cj] == FLUID:
                                x0, x1 = (0.5 * (xs[i - 1] + xs[i]), xs[i]) if ci == i - 1 else (xs[i], 0.5 * (xs[i] + xs[i + 1]))
                                W += abs(x1 - x0)
                        if W == 0.0:
                            continue
                        m = 0.5 * (self.ro[i, j] * self.v[i, j] + self.ro[i, jj] * self.v[i, jj]) * W * self.rw(0.5 * (ys[j] + ys[jj])) * side
                        if m == 0.0:
                            continue
                        a[idx(i, jj)] = a.get(idx(i, jj), 0.0) + 0.5 * cp * m
                        a[p] = a.get(p, 0.0) + 0.5 * cp * m - cp * m
                    b[p] += self.S[i, j] * Vf
                Vs = 0.0
                for di in (-1, 0):
                    for dj in (-1, 0):
                        ci, cj = i + di, j + dj
                        if 0 <= ci < nx - 1 and 0 <= cj < ny - 1 and self.mat[ci, cj] == SOLID:
                            xa, xb = (0.5 * (xs[i - 1] + xs[i]), xs[i]) if di == -1 else (xs[i], 0.5 * (xs[i] + xs[i + 1]))
                            ya, yb = (0.5 * (ys[j - 1] + ys[j]), ys[j]) if dj == -1 else (ys[j], 0.5 * (ys[j] + ys[j + 1]))
                            Vs += abs(xb - xa) * abs(yb - ya) * self.rw(0.5 * (ya + yb))
                if Vs > 0:
                    b[p] += self.Ss[i, j] * Vs
                # Robin: 外周の辺 (空セルか領域外に面する辺) の固体側
                if self.robin is not None:
                    for dj in (1, -1):                          # 上下の外周 (y 方向の境界)
                        for ci in (i - 1, i):
                            if not (0 <= ci < nx - 1):
                                continue
                            inside = self.kcell(ci, j - 1 if dj == 1 else j)
                            outside = self.kcell(ci, j if dj == 1 else j - 1)
                            if inside is not None and self.mat[ci, j - 1 if dj == 1 else j] == SOLID and outside is None:
                                x0, x1 = (0.5 * (xs[i - 1] + xs[i]), xs[i]) if ci == i - 1 else (xs[i], 0.5 * (xs[i] + xs[i + 1]))
                                hr = self.robin(0.5 * (x0 + x1), ys[j])
                                if hr is not None:
                                    h, Tc = hr
                                    A = abs(x1 - x0) * self.rw(ys[j])
                                    a[p] = a.get(p, 0.0) + h * A; b[p] += h * A * Tc
                for q, w in a.items():
                    add(p, q, w)
        A = sp.csr_matrix((vals, (rows, cols)), shape=(N, N))
        self.T = spla.spsolve(A.tocsc(), b).reshape(nx, ny)
        return self.T

    # ---------------------------------------------------------------- 診断
    def robin_heat(self):
        """固体外面の Robin から入った熱 (per rad または per m)。"""
        xs, ys, T = self.xs, self.ys, self.T
        Q = 0.0
        for i in range(self.nx):
            for j in range(self.ny):
                for dj in (1, -1):
                    for ci in (i - 1, i):
                        if not (0 <= ci < self.nx - 1):
                            continue
                        cin = j - 1 if dj == 1 else j
                        cout = j if dj == 1 else j - 1
                        if not (0 <= cin < self.ny - 1) or self.mat[ci, cin] != SOLID:
                            continue
                        if self.kcell(ci, cout) is not None:
                            continue
                        x0, x1 = (0.5 * (xs[i - 1] + xs[i]), xs[i]) if ci == i - 1 else (xs[i], 0.5 * (xs[i] + xs[i + 1]))
                        hr = self.robin(0.5 * (x0 + x1), ys[j]) if self.robin else None
                        if hr is None:
                            continue
                        h, Tc = hr
                        Q += h * abs(x1 - x0) * self.rw(ys[j]) * (Tc - T[i, j])
        return Q

    def enthalpy_flux(self, i):
        """列 i の対流エンタルピー流束 ∫ρ u c_p (T − T_in) y dy (流体節点、台形)。"""
        y = self.ys; w = np.array([self.rw(t) for t in y])
        f = self.ro[i] * self.u[i] * self.cp * (self.T[i] - self.T_in) * w
        return np.trapezoid(f, y) if hasattr(np, "trapezoid") else np.trapz(f, y)


# ------------------------------------------------------------------ 自己検査
def graded(a, b, n, h0, toward_a=True):
    """[a,b] を n 区間、端 a (toward_a) で幅 h0 から等比で広げる。"""
    L = b - a
    lo, hi = 1.0 + 1e-12, 3.0
    for _ in range(200):
        q = 0.5 * (lo + hi)
        s = h0 * (q ** n - 1) / (q - 1)
        lo, hi = (q, hi) if s < L else (lo, q)
    w = h0 * q ** np.arange(n)
    w *= L / w.sum()
    x = a + np.concatenate([[0.0], np.cumsum(w)])
    return x if toward_a else (a + b - x)[::-1]


def pipe_problem(nr_f, nr_s, nx_scale, R, r_o, k_f, k_s, rho, cp, Um, L_u, L_h, L_d, heat_robin, wall_T=None,
                 solid_everywhere=True):
    """A 型: 管 (流体 0..R) + 固体殻 (R..r_o)。Poiseuille、源項なし。"""
    h0 = 0.02e-3 / nx_scale
    xu = graded(-L_u, 0.0, int(40 * nx_scale), h0, toward_a=False)
    xh = graded(0.0, L_h, int(160 * nx_scale), h0)
    xd = graded(L_h, L_h + L_d, int(20 * nx_scale), max(h0, (xh[-1] - xh[-2])))
    xs = np.concatenate([xu, xh[1:], xd[1:]])
    yf = np.linspace(0, R, nr_f + 1)
    ys = np.concatenate([yf, np.linspace(R, r_o, nr_s + 1)[1:]]) if nr_s > 0 else yf
    mat = np.full((len(xs) - 1, len(ys) - 1), FLUID)
    mat[:, nr_f:] = SOLID
    nx, ny = len(xs), len(ys)
    u = np.zeros((nx, ny)); ro = np.full((nx, ny), rho)
    u[:, :nr_f + 1] = 2 * Um * (1 - (yf / R) ** 2)[None, :]
    rob = None
    if heat_robin is not None:
        h, Tc = heat_robin
        rob = lambda xm, y: (h, Tc) if (0.0 <= xm <= L_h and abs(y - r_o) < 1e-12) else None
    wd = None
    if wall_T is not None:
        wd = lambda x: wall_T if (-1e-12 <= x <= L_h + 1e-12) else None
    return Problem(xs, ys, mat, k_f, k_s, cp, True, 0.0, ro=ro, u=u, robin=rob, wall_dirichlet=wd, j_wall=nr_f), xs, ys


def selftest():
    sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[1] / "63.graetz_cht"))
    import graetz_ref
    R, k_f, cp, rho = 1e-3, 5.700284e-2, 1004.5, 1.176829
    Um = 17.35944; D = 2 * R
    Pe = rho * cp * Um * D / k_f
    L_u, L_h, L_d = 10e-3, 0.1728, 10e-3
    ok = True
    # (i) 固体なし・壁 Dirichlet → graetz_ref ellip と一致
    P, xs, ys = pipe_problem(60, 0, 1.0, R, R, k_f, k_f, rho, cp, Um, L_u, L_h, L_d, None, wall_T=1.0)
    P.T_in = 0.0
    T = P.solve()
    xu, xd = L_u / (D * Pe), L_d / (D * Pe)
    xe, nue, tbe, xie, the = graetz_ref.solve_ellip(120, Pe, xu, 0.12, xd, 1600, 300, 300, return_field=True)
    worst = 0.0
    for xmm in (4.3, 14.4, 43.0, 144.0):
        i = int(np.argmin(np.abs(xs - xmm * 1e-3))); xp = xs[i] / (D * Pe); ie = int(np.argmin(np.abs(xe - xp)))
        ref = 1.0 - np.interp(ys / R, xie, the[ie])                     # θ = (T−T_w)/(T_in−T_w)、T_in=0・T_w=1
        e = np.abs(T[i] - ref).max(); worst = max(worst, e)
    print(f"(i) 固体なし・壁温一様 vs graetz_ref ellip: max |ΔT|/ΔT = {worst:.2e} (許容 3e-3)  {'ok' if worst < 3e-3 else 'NG'}")
    ok &= worst < 3e-3
    # (ii) 熱収支: Robin から入った熱 = 出口のエンタルピー流束
    P, xs, ys = pipe_problem(40, 20, 1.0, R, 2 * R, k_f, 10 * k_f, rho, cp, Um, 20e-3, 0.05, 10e-3, (100.0, 10.0))
    P.T_in = 0.0; P.solve()
    Qin = P.robin_heat(); Qout = P.enthalpy_flux(P.nx - 1)
    e2 = abs(Qin - Qout) / abs(Qin)
    print(f"(ii) 熱収支: Robin 入熱 {Qin:.6e} / 出口エンタルピー {Qout:.6e} → 相対差 {e2:.2e} (許容 1e-3)  {'ok' if e2 < 1e-3 else 'NG'}")
    ok &= e2 < 1e-3
    # (iii) k_s→大・h 大 → 界面温度が T_c に近づく (等温壁の極限)
    P, xs, ys = pipe_problem(40, 10, 1.0, R, 1.2 * R, k_f, 1e6 * k_f, rho, cp, Um, 20e-3, 0.05, 10e-3, (1e8, 1.0))
    P.T_in = 0.0; T = P.solve()
    heat = (xs >= 0) & (xs <= 0.05)
    e3 = np.abs(T[heat, 40] - 1.0).max()
    print(f"(iii) k_s→大: 加熱区間の界面温度と T_c の差 max = {e3:.2e} (許容 1e-3)  {'ok' if e3 < 1e-3 else 'NG'}")
    ok &= e3 < 1e-3
    ok &= mms(True) & mms(False)
    print("VERDICT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def mms(axisym):
    """製造解 (2 材料・有限の固体抵抗・v ≠ 0・体積源あり・軸を含む)。誤差が格子の倍増で約 1/4 (2 次) に減ること。
    流体 0 ≤ y ≤ 1、固体 1 ≤ y ≤ 1.5、x ∈ [0, 1]。T = F(x) g(y)、F = 1 + x − x²/2、
    g_f = 1 + a y² (軸で対称)、g_s = g_f(1) + (k_f/k_s) g_f'(1) (y − 1) (界面で温度・熱流束が連続)。
    u = 1 − y²、v = 0.3 y (1 − y) (流体)、ρ = c_p = 1、k_f 0.1、k_s 0.5。境界 (x 両端・上端、平面では下端も) は厳密解の Dirichlet。"""
    kf, ks, a = 0.1, 0.5, 0.8
    F = lambda x: 1 + x - x * x / 2; Fx = lambda x: 1 - x; Fxx = -1.0
    gf = lambda y: 1 + a * y * y; gfy = lambda y: 2 * a * y
    gs = lambda y: gf(1.0) + (kf / ks) * gfy(1.0) * (y - 1.0)
    def Tex(x, y): return F(x) * (gf(y) if y <= 1.0 + 1e-14 else gs(y))
    errs = []
    for n in (8, 16, 32):
        xs = np.linspace(0, 1, 2 * n + 1)
        ys = np.concatenate([np.linspace(0, 1, n + 1), np.linspace(1, 1.5, n // 2 + 1)[1:]])
        mat = np.full((len(xs) - 1, len(ys) - 1), FLUID); mat[:, n:] = SOLID
        X, Y = np.meshgrid(xs, ys, indexing="ij")
        fl = Y <= 1.0 + 1e-14
        u = np.where(fl, 1 - Y ** 2, 0.0); v = np.where(fl, 0.3 * Y * (1 - Y), 0.0)
        rw = Y if axisym else np.ones_like(Y)
        # 流体の源: u T_x + v T_y − k_f ∇²T  (軸対称 ∇² = T_xx + T_yy + T_y / y、y→0 は 2 T_yy)
        Tf_x = Fx(X) * gf(Y); Tf_y = F(X) * gfy(Y); lap_f = Fxx * gf(Y) + F(X) * 2 * a + (F(X) * 2 * a if axisym else 0.0)
        Sf = u * Tf_x + v * Tf_y - kf * lap_f
        gsyy = 0.0; gsy = (kf / ks) * gfy(1.0)
        lap_s = Fxx * np.vectorize(gs)(Y) + F(X) * gsyy + (F(X) * gsy / np.where(Y > 0, Y, 1.0) if axisym else 0.0)
        Ss = -ks * lap_s
        P = Problem(xs, ys, mat, kf, ks, 1.0, axisym, 0.0, ro=np.ones_like(X), u=u, v=v, source=np.where(fl, Sf, 0.0),
                    source_solid=np.where(Y >= 1.0 - 1e-14, Ss, 0.0), exact_boundary=Tex)   # 界面の行にも固体側の源を与える (双対セルの固体半分)
        T = P.solve()
        E = np.abs(T - np.vectorize(Tex)(X, Y)).max()
        errs.append(E)
    r1, r2 = errs[0] / errs[1], errs[1] / errs[2]
    good = r2 > 3.0
    print(f"(mms {'軸対称' if axisym else '平面'}) 2 材料・v≠0・源あり: max 誤差 {errs[0]:.2e} / {errs[1]:.2e} / {errs[2]:.2e}、比 {r1:.2f} / {r2:.2f} (最後の比 > 3 で 2 次)  {'ok' if good else 'NG'}")
    return good


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "selftest":
        sys.exit(selftest())
    print(__doc__)
