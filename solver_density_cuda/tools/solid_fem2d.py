#!/usr/bin/env python3
r"""CHT の固体側: **一般 2 次元領域**の伝導 (`fem2d`)。冷却孔は Robin 境界。

薄肉シェル (`solid_shell.py`) では表せない固体 — 内部に冷却孔を持つ厚いタービン翼など — 用。
plan §4.4c/§4.4d、methods/boundary.md「共役熱伝達 (CHT)」。

## 定式化

固体全節点温度を $u$、界面抽出を $E$ (界面節点 ← 固体節点) として

$$\left(K_s + E^{\mathsf T} D_f E\right)u^{k+1} = b_s + E^{\mathsf T}\!\left[Q_f(Eu^{k}) + D_f\,Eu^{k}\right]$$

$K_s$ は線形三角形 FE の剛性 + Robin 辺の寄与、$b_s$ は Robin 辺の荷重。
$Q_f$ は**積分済みの節点荷重** [W] (平面 2D は W/m) で、$E^{\mathsf T}$ で載せるときに**面積を再乗算しない**。

`cht_loop.py` / `FixedPointDriver` からは `ShellOperator` と同じ顔 (`n`, `coords`, `area`,
`assemble`, `driver`) に見えるよう、**内部節点を消去した Schur 補元**を界面作用素として返す:

$$S = K_{ii} - K_{io}K_{oo}^{-1}K_{oi},\qquad
  b^{\rm iface} = b_i - K_{io}K_{oo}^{-1}b_o$$

これで固定点も残差もシェルと同じ言葉になる (単位も [W/K] と [W])。

自己検証: python3 solver_density_cuda/tools/test_solid_fem2d.py
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from solid_shell import FixedPointDriver, ShellOperator


class Fem2DOperator:
    """平面 2 次元固体 (単位奥行き)。荷重・熱量の単位は W/m。

    nodes       : (N,2) or (N,3) 節点座標 [m]
    tris        : (M,3) 三角形
    iface_nodes : 界面 (ガス側表面) の節点 index。**流体側の壁節点と 1 対 1** に対応させる
    iface_edges : 界面の辺 [(n0,n1), ...] (集中長さ = 荷重変換に使う)
    robin_edges : [(n0, n1, h, T_c), ...] 冷却孔など。$h$ [W/m2K], $T_c$ [K]
    k_solid     : 定数 or (T[], k[]) テーブル
    """

    def __init__(self, nodes, tris, iface_nodes, iface_edges, robin_edges, k_solid):
        self.xy = np.asarray(nodes, float)[:, :2]
        self.tris = np.asarray(tris, int)
        self.iface = np.asarray(iface_nodes, int)
        self.iface_edges = [tuple(e) for e in iface_edges]
        self.robin_edges = list(robin_edges)
        self.k_solid = k_solid
        self.N = len(self.xy)
        self.n = len(self.iface)
        self.coords = np.zeros((self.n, 3))
        self.coords[:, :2] = self.xy[self.iface]
        self.area = self._iface_lumped()
        self.u = None                      # 直近の固体全温度 (k_s(T) の Picard 用)
        self._other = np.setdiff1d(np.arange(self.N), self.iface)
        if self.area.min() <= 0:
            raise ValueError("interface lumped length must be positive (check iface_edges)")

    # ---- 幾何 ----
    def _iface_lumped(self):
        a = np.zeros(self.n)
        pos = {int(g): i for i, g in enumerate(self.iface)}
        for (n0, n1) in self.iface_edges:
            L = float(np.linalg.norm(self.xy[n1] - self.xy[n0]))
            for g in (n0, n1):
                if int(g) in pos:
                    a[pos[int(g)]] += 0.5 * L
        return a

    def k_of(self, T):
        if isinstance(self.k_solid, tuple):
            Tt, kt = self.k_solid
            return np.interp(np.asarray(T, float), np.asarray(Tt, float), np.asarray(kt, float))
        return np.full(np.shape(T), float(self.k_solid))

    # ---- 全系の組立て ----
    def assemble_full(self, u):
        k = self.k_of(u)
        rows, cols, vals = [], [], []
        for t in self.tris:
            p = self.xy[t]
            b = np.array([p[1, 1] - p[2, 1], p[2, 1] - p[0, 1], p[0, 1] - p[1, 1]])
            c = np.array([p[2, 0] - p[1, 0], p[0, 0] - p[2, 0], p[1, 0] - p[0, 0]])
            det = (p[1, 0] - p[0, 0]) * (p[2, 1] - p[0, 1]) - (p[2, 0] - p[0, 0]) * (p[1, 1] - p[0, 1])
            area = 0.5 * abs(det)
            if area <= 0:
                continue
            ke = float(np.mean(k[t])) / (4.0 * area) * (np.outer(b, b) + np.outer(c, c))
            for i in range(3):
                for j in range(3):
                    rows.append(t[i]); cols.append(t[j]); vals.append(ke[i, j])
        K = sp.csr_matrix((vals, (rows, cols)), shape=(self.N, self.N))
        b = np.zeros(self.N)
        # Robin 辺: ∫ h (T - T_c) v ds を線形辺の consistent 行列で入れる
        for (n0, n1, h, Tc) in self.robin_edges:
            L = float(np.linalg.norm(self.xy[n1] - self.xy[n0]))
            if L <= 0:
                continue
            Me = h * L / 6.0 * np.array([[2.0, 1.0], [1.0, 2.0]])
            idx = [int(n0), int(n1)]
            for i in range(2):
                for j in range(2):
                    K[idx[i], idx[j]] += Me[i, j]
                b[idx[i]] += h * Tc * L / 2.0
        return K.tocsr(), b

    # ---- 逆問題・繰り返し解析用: 作用素を部品に分けて先に組む ----
    def parts(self, T_ref=None, groups=None):
        r"""伝導剛性と **Robin 群ごとの行列・荷重**を一度だけ組んで返す。

        $K(h_1..h_m) = K_{cond} + \sum_k h_k M_k$, $b = \sum_k h_k T_{c,k} v_k$ なので、
        $h_k, T_{c,k}$ を振る解析 (内部条件の同定など) は**毎回組み直さずに済む**。

        groups: 辺リストの配列 (省略時は self.robin_edges を 1 群として扱う)
        戻り値: (K_cond, [(M_k, v_k), ...])
        """
        u_ref = np.full(self.N, float(T_ref) if T_ref is not None else 300.0)
        saved, self.robin_edges = self.robin_edges, []      # 伝導だけ組む
        K_cond, _ = self.assemble_full(u_ref)
        self.robin_edges = saved
        if groups is None:
            groups = [[(e[0], e[1]) for e in self.robin_edges]]
        out = []
        for g in groups:
            rows, cols, vals = [], [], []
            v = np.zeros(self.N)
            for e in g:
                n0, n1 = int(e[0]), int(e[1])
                L = float(np.linalg.norm(self.xy[n1] - self.xy[n0]))
                if L <= 0:
                    continue
                Me = L / 6.0 * np.array([[2.0, 1.0], [1.0, 2.0]])
                idx = [n0, n1]
                for i in range(2):
                    for j in range(2):
                        rows.append(idx[i]); cols.append(idx[j]); vals.append(Me[i, j])
                    v[idx[i]] += L / 2.0
            out.append((sp.csr_matrix((vals, (rows, cols)), shape=(self.N, self.N)), v))
        return K_cond, out

    # ---- 界面への縮約 (Schur 補元) ----
    def assemble(self, T_iface):
        u = self.u if self.u is not None else np.full(self.N, float(np.mean(T_iface)))
        K, b = self.assemble_full(u)
        i, o = self.iface, self._other
        Kii = K[i][:, i].toarray() if self.n <= 400 else K[i][:, i]
        Kio = K[i][:, o]
        Koi = K[o][:, i]
        Koo = K[o][:, o].tocsc()
        lu = spla.splu(Koo)
        W = np.column_stack([lu.solve(np.asarray(Koi[:, j].todense()).ravel()) for j in range(self.n)])
        S = (Kii if isinstance(Kii, np.ndarray) else Kii.toarray()) - Kio @ W
        bi = b[i] - Kio @ lu.solve(b[o])
        self._lu, self._Kio, self._b = lu, Kio, b
        return sp.csr_matrix(S), np.asarray(bi).ravel()

    def recover_interior(self, T_iface):
        """界面温度から内部節点を復元する (診断・可視化用)。"""
        i, o = self.iface, self._other
        Koi = self._Kio.T
        uo = self._lu.solve(self._b[o] - np.asarray((Koi @ T_iface)).ravel())
        u = np.zeros(self.N)
        u[i] = T_iface
        u[o] = uo
        self.u = u
        return u

    def interior_residual(self, T_iface):
        r"""**内部節点の残差** $\max|K_{oo}u_o + K_{oi}T_i - b_o|$ [W/m] (G-if の 1 項目)。

        界面へ縮約した Schur 作用素だけを見ていると、内部の復元 (`recover_interior`) と
        物性の局所評価がずれていても気づけない。同じ評価温度で全系を組み直して残差を測る。
        """
        u = self.recover_interior(np.asarray(T_iface, float))
        K, b = self.assemble_full(u)
        res = np.asarray(K.dot(u)).ravel() - np.asarray(b).ravel()
        o = self._other
        return float(np.max(np.abs(res[o]))) if len(o) else 0.0

    def solve(self, Qf, T0=None, iters=15, tol=1e-10):
        """節点熱荷重 Qf [W/m] を界面に与えて解く (k_s(T) は Picard)。"""
        T = np.full(self.n, 300.0) if T0 is None else np.asarray(T0, float).copy()
        for _ in range(iters):
            S, bi = self.assemble(T)
            Tn = spla.spsolve(S.tocsc(), bi + np.asarray(Qf, float))
            self.recover_interior(Tn)
            if np.max(np.abs(Tn - T)) < tol * max(1.0, np.max(np.abs(Tn))):
                return Tn
            T = Tn
        return T

    def driver(self, T0, Df0=None, anderson=5, Df_cap=1e12):
        return FixedPointDriver(self, T0, Df0, anderson, Df_cap)

    # FixedPointDriver が使う静的ヘルパはシェルと共有する
    _merit_M = staticmethod(ShellOperator._merit_M)
    _linsolve = staticmethod(ShellOperator._linsolve)
    _anderson = staticmethod(ShellOperator._anderson)
