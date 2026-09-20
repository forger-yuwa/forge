#!/usr/bin/env python3
r"""CHT の固体側モデル (薄肉シェル)。**未知数はガス側表面温度 $T_w$**。

仕様の正本は methods/boundary.md「共役熱伝達 (CHT)」、設計判断は
plans/active/boundary-conjugate-heat-transfer.md。ここはその実装。

## モデル

$$\nabla_s\cdot(k_s t\,\nabla_s T_w) + q_{gas} - \frac{T_w - T_b}{R_{tot}} = 0,\qquad
  R_{tot} = \frac{t}{k_s} + R_{back}$$

- **背面等温でも $t/k_s$ を落とさない** (落とすと $T_w=T_b$ に退化する)。
- `mode`:
  - `local1d` — 面内伝導を落とした極限 ($k_s t\to0$)。点ごとに $T_w = T_b + q\,R_{tot}$。
  - `shell2d` — 面内伝導つき。境界面 (primal facet) の線形 FE で $\nabla_s$ を離散化する。
- **断熱・孤立系** ($R_{back}=\infty$ かつ端部断熱) は定数零空間を持ち、正味入熱が非零なら定常解が無い。
  `assemble` は零空間を検出し、適合条件 ($\sum Q \approx 0$) を満たさない入力を例外にする。

## 連成反復 (固定点を保存する形)

$$\left(A_s + D_f\right)T^{k+1} = b_s + Q_f(T^k) + D_f T^k$$

$D_f$ が何であっても固定点 $A_s T^* = b_s + Q_f(T^*)$ は変わらない ($D_f$ は収束速度だけを決める)。
受理は**メリット関数** $\Phi(r)=r^\mathsf{T}(A_s+D_f)^{-1}r$ の降下で判定する
(**残差最大ノルムの単調減少を受理条件にしない**: 収束する反復を棄却する反例がある)。

単位: 熱荷重 $Q_f$ は**積分済みの節点荷重 [W]** (平面 2D は W/m)。面積を再乗算しないこと。

使い方:
  from solid_shell import SolidModel, ShellOperator
  model = SolidModel.from_json("solid.json")
  op    = ShellOperator(coords, faces, model)      # shell2d
  Tw    = op.solve(Qf)                             # 片方向 (与えられた熱荷重で解く)
  Tw, info = op.couple(Qf_callable, T0)            # 連成反復 (受理判定つき)

自己検証: python3 solver_density_cuda/tools/test_solid_shell.py
"""
from __future__ import annotations

import json
import math
from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


# --------------------------------------------------------------------------- model
@dataclass
class SolidModel:
    """固体の物性と背面条件。thickness/k は定数またはノード配列 (長さ n)。"""
    mode: str = "shell2d"                 # local1d | shell2d
    thickness: float | np.ndarray = 1.0e-3
    k_solid: float | np.ndarray = 15.0    # [W/mK] 定数 or ノード配列
    k_table: tuple | None = None          # (T[], k[]) を与えると k_s(T) に温度依存を入れる
    back_kind: str = "isothermal"         # isothermal | coolant | layers | adiabatic
    T_b: float = 300.0
    h_c: float = 0.0                      # coolant: [W/m2K]
    layers: list = field(default_factory=list)   # layers: [{"t":..,"k":..}, ...]

    @staticmethod
    def from_json(path: str) -> "SolidModel":
        with open(path) as f:
            d = json.load(f)
        back = d.get("back", {})
        kt = d.get("k_table")
        return SolidModel(
            mode=d.get("mode", "shell2d"),
            thickness=np.asarray(d["thickness"], float) if isinstance(d.get("thickness"), list) else float(d.get("thickness", 1.0e-3)),
            k_solid=np.asarray(d["k_solid"], float) if isinstance(d.get("k_solid"), list) else float(d.get("k_solid", 15.0)),
            k_table=(np.asarray(kt["T"], float), np.asarray(kt["k"], float)) if kt else None,
            back_kind=back.get("kind", "isothermal"),
            T_b=float(back.get("T_b", 300.0)),
            h_c=float(back.get("h", 0.0)),
            layers=back.get("layers", []),
        )

    # ---- 物性 ----
    def k_of(self, T) -> np.ndarray:
        if self.k_table is not None:
            Tt, kt = self.k_table
            return np.interp(np.asarray(T, float), Tt, kt)
        return np.asarray(self.k_solid, float) * np.ones_like(np.asarray(T, float))

    def t_of(self, n: int) -> np.ndarray:
        t = np.asarray(self.thickness, float)
        return t if t.ndim else np.full(n, float(t))

    def R_back(self, n: int) -> np.ndarray:
        """背面環境までの抵抗 (表面→背面の t/k_s は含めない)。"""
        if self.back_kind == "isothermal":
            return np.zeros(n)
        if self.back_kind == "coolant":
            if not self.h_c > 0.0:
                raise ValueError("back.kind=coolant requires h > 0")
            return np.full(n, 1.0 / self.h_c)
        if self.back_kind == "layers":
            R = sum(float(L["t"]) / float(L["k"]) for L in self.layers)
            return np.full(n, R)
        if self.back_kind == "adiabatic":
            return np.full(n, np.inf)
        raise ValueError(f"unknown back.kind={self.back_kind}")

    def R_tot(self, T) -> np.ndarray:
        """ガス側表面から背面環境までの**全抵抗** [m2K/W]。t/k_s を必ず含む。"""
        T = np.asarray(T, float)
        n = T.size
        return self.t_of(n) / self.k_of(T) + self.R_back(n)


# --------------------------------------------------------------------------- geometry
def _tri_stiffness(p0, p1, p2, kt):
    """三角形の線形 FE 剛性 (曲面: 三角形平面上で評価)。kt = k_s t [W/K]."""
    v1, v2 = p1 - p0, p2 - p0
    nrm = np.cross(v1, v2)
    area2 = np.linalg.norm(nrm)
    if area2 <= 0:
        return None, 0.0
    area = 0.5 * area2
    # 三角形平面の 2D 座標へ
    e1 = v1 / np.linalg.norm(v1)
    e2 = np.cross(nrm / area2, e1)
    q = np.array([[0.0, 0.0],
                  [np.dot(v1, e1), np.dot(v1, e2)],
                  [np.dot(v2, e1), np.dot(v2, e2)]])
    b = np.array([q[1, 1] - q[2, 1], q[2, 1] - q[0, 1], q[0, 1] - q[1, 1]])
    c = np.array([q[2, 0] - q[1, 0], q[0, 0] - q[2, 0], q[1, 0] - q[0, 0]])
    K = (kt / (4.0 * area)) * (np.outer(b, b) + np.outer(c, c))
    return K, area


def _faces_to_tris(faces):
    """面を三角形に割り、(三角形, 重み) を返す。

    **四角形は 2 通りの対角線分割を 1/2 ずつ使う**。片方の対角線だけで割ると、
    集中面積 (質量集中) が対角の 2 頂点と残り 2 頂点で非対称になり、
    **境界の節点で O(1) の荷重不均衡**が残って全体の収束次数が 2 次 -> 1 次に落ちる
    (実測: 帯フィン問題で rate 1.00)。両分割の平均なら対称性が戻る。
    """
    for f in faces:
        f = list(f)
        if len(f) == 2:
            continue                      # 線要素は stiffness/_lumped_area が直接扱う
        if len(f) == 3:
            yield f, 1.0
        elif len(f) == 4:
            yield [f[0], f[1], f[2]], 0.5
            yield [f[0], f[2], f[3]], 0.5
            yield [f[1], f[2], f[3]], 0.5
            yield [f[1], f[3], f[0]], 0.5
        elif len(f) > 4:
            for i in range(1, len(f) - 1):
                yield [f[0], f[i], f[i + 1]], 1.0


# --------------------------------------------------------------------------- operator
class ShellOperator:
    """境界面メッシュ上の薄肉シェル作用素。

    coords : (n,3) ノード座標 [m]
    faces  : primal 境界面の節点リスト (三角/四角/多角)
    model  : SolidModel
    axisym : True なら軸対称 (y = 半径) として面内伝導・面積に r 重みを掛ける
    """

    def __init__(self, coords, faces, model: SolidModel, axisym: bool = False):
        self.coords = np.asarray(coords, float)
        self.n = len(self.coords)
        self.faces = [list(f) for f in faces]
        self.model = model
        self.axisym = bool(axisym)
        self.area = self._lumped_area()          # 節点集中面積 [m2] (軸対称は 2πr 重み)
        if self.area.min() <= 0:
            raise ValueError("lumped nodal area must be positive (check the face list)")

    # ---- 幾何 ----
    def _rw(self, idx):
        if not self.axisym:
            return 1.0
        return max(float(np.mean(self.coords[idx, 1])), 0.0) * 2.0 * math.pi

    def _lumped_area(self):
        a = np.zeros(self.n)
        # 線要素 (2 節点) = 2D 平面ケースの壁。単位奥行きあたりで扱う (荷重は W/m)。
        for f in self.faces:
            if len(f) != 2:
                continue
            p = self.coords[f]
            Le = float(np.linalg.norm(p[1] - p[0]))
            a[f[0]] += 0.5 * Le * self._rw(f)
            a[f[1]] += 0.5 * Le * self._rw(f)
        for tri, wgt in _faces_to_tris(self.faces):
            p = self.coords[tri]
            ar = 0.5 * np.linalg.norm(np.cross(p[1] - p[0], p[2] - p[0]))
            a[tri] += wgt * ar * self._rw(tri) / 3.0
        return a

    def stiffness(self, T):
        """面内伝導の剛性 A_in [W/K] (shell2d)。local1d ではゼロ行列。"""
        if self.model.mode == "local1d":
            return sp.csr_matrix((self.n, self.n))
        k = self.model.k_of(T)
        t = self.model.t_of(self.n)
        rows, cols, vals = [], [], []
        for f in self.faces:
            if len(f) != 2:
                continue
            p = self.coords[f]
            Le = float(np.linalg.norm(p[1] - p[0]))
            if not Le > 0:
                continue
            kt = float(np.mean(k[f] * t[f])) * self._rw(f) / Le
            for a_ in range(2):
                for b_ in range(2):
                    rows.append(f[a_]); cols.append(f[b_]); vals.append(kt * (1.0 if a_ == b_ else -1.0))
        for tri, wgt in _faces_to_tris(self.faces):
            kt = float(np.mean(k[tri] * t[tri])) * self._rw(tri)
            K, _ = _tri_stiffness(*self.coords[tri], kt)
            if K is None:
                continue
            for a in range(3):
                for b in range(3):
                    rows.append(tri[a]); cols.append(tri[b]); vals.append(wgt * K[a, b])
        return sp.csr_matrix((vals, (rows, cols)), shape=(self.n, self.n))

    def assemble(self, T):
        """A_s T = b_s + Q_f の A_s [W/K] と b_s [W] (背面項込み)。"""
        R = self.model.R_tot(T)
        with np.errstate(divide="ignore"):
            g = np.where(np.isfinite(R) & (R > 0), self.area / np.where(R > 0, R, 1.0), 0.0)
        A = self.stiffness(T) + sp.diags(g)
        b = g * self.model.T_b
        return A.tocsr(), b

    # ---- 片方向 (熱荷重を与えて解く) ----
    def solve(self, Qf, T0=None, iters=20, tol=1e-10):
        """節点熱荷重 Qf [W] (固体向き正) を与えて T_w を解く。k_s(T) は Picard。"""
        Qf = np.asarray(Qf, float)
        T = np.full(self.n, self.model.T_b, float) if T0 is None else np.asarray(T0, float).copy()
        for _ in range(iters):
            A, b = self.assemble(T)
            self._check_solvable(A, b + Qf)
            Tn = self._linsolve(A, b + Qf, T)
            if np.max(np.abs(Tn - T)) < tol * max(1.0, np.max(np.abs(Tn))):
                T = Tn
                break
            T = Tn
        return T

    def _check_solvable(self, A, rhs):
        """断熱・孤立系 (零空間あり) を検出して、適合条件を満たさない入力を弾く。"""
        s = np.asarray(A.sum(axis=1)).ravel()
        if np.max(np.abs(s)) > 1e-10 * max(1.0, np.max(np.abs(A.diagonal()))):
            return                        # Robin 項があるので正則
        net = float(np.sum(rhs))
        scale = float(np.sum(np.abs(rhs))) or 1.0
        raise ValueError(
            "solid operator is singular (adiabatic/isolated: constant null space). "
            f"net load = {net:.6g} W (|load| sum {scale:.6g}); "
            "give a Robin/isothermal back condition or pin the constant."
        )

    @staticmethod
    def _linsolve(A, rhs, T0):
        n = A.shape[0]
        if n <= 2000:
            return spla.spsolve(A.tocsc(), rhs)
        M = spla.LinearOperator(A.shape, matvec=lambda x: x / A.diagonal())
        x, info = spla.cg(A, rhs, x0=T0, rtol=1e-12, atol=0.0, maxiter=5000, M=M)
        if info != 0:
            raise RuntimeError(f"CG failed (info={info})")
        return x

    # ---- 連成反復 ----
    def couple(self, Qf_of_T, T0, Df0=None, max_iter=200, tol_K=1e-6, tol_rel=1e-6,
               anderson=5, Df_cap=1e12, verbose=False):
        r"""$(A_s+D_f)T^{k+1} = b_s + Q_f(T^k) + D_f T^k$ を回す (固定点は $D_f$ に依らない)。

        Qf_of_T : T_w (n,) -> 節点熱荷重 [W] (固体向き正)
        Df0     : 界面コンダクタンス [W/K] の初期推定。**上界ではない** (plan §4.2 の反例)。
        anderson: Anderson 加速の深さ (0 で素の固定点反復)。スカラー $D_f$ では扱えない
                  **非対角な流体応答**はこれで潰す。加速候補は残差ノルムが減らなければ棄却する。
        受理     : メリット関数 $\Phi=r^\mathsf{T}(A_s+D_f)^{-1}r$ の降下 (**$D_f$ を固定して比較**)。
                  降下しなければ $D_f$ を倍にして退避・再試行する。
                  **残差最大ノルムの単調減少は受理条件にしない** (収束する反復を棄却する反例がある)。
        収束     : max|ΔT| < tol_K [K] かつ max|r| < tol_rel * max(|Q_f|, |b_s|)
        """
        T = np.asarray(T0, float).copy()
        Df = np.full(self.n, 1.0) if Df0 is None else np.asarray(Df0, float).copy()
        hist = []
        Ts, Gs = [], []                      # Anderson の履歴 (T_k, G(T_k))

        def residual(Tx):
            Qx = Qf_of_T(Tx)
            Ax, bx = self.assemble(Tx)
            return Ax @ Tx - bx - Qx, Qx, Ax, bx

        r, Q, A, b = residual(T)
        for it in range(max_iter):
            accepted = False
            for _ in range(12):
                M = (A + sp.diags(Df)).tocsr()
                phi0 = self._merit_M(M, r)
                G = self._linsolve(M, b + Q + Df * T, T)          # 素の固定点写像
                cands = [("plain", G)]
                if anderson > 0 and len(Ts) >= 1:
                    acc = self._anderson(Ts + [T], Gs + [G], anderson)
                    if acc is not None:
                        cands.append(("anderson", acc))
                best = None
                for name, Tc in cands:
                    rc, Qc, Ac, bc = residual(Tc)
                    phic = self._merit_M((Ac + sp.diags(Df)).tocsr(), rc)
                    if best is None or phic < best[0]:
                        best = (phic, name, Tc, rc, Qc, Ac, bc)
                if best is not None and best[0] < phi0:
                    phic, name, Tc, rc, Qc, Ac, bc = best
                    dT = float(np.max(np.abs(Tc - T)))
                    Ts.append(T.copy()); Gs.append(G.copy())
                    if len(Ts) > anderson + 1:
                        Ts.pop(0); Gs.pop(0)
                    T, r, Q, A, b = Tc, rc, Qc, Ac, bc
                    accepted, used = True, name
                    break
                Df = np.minimum(Df * 2.0, Df_cap)                 # 降下しない -> D_f を上げて退避
                Ts.clear(); Gs.clear()                            # 重みが変わったので履歴は捨てる
            scale = max(float(np.max(np.abs(Q))), float(np.max(np.abs(b))), 1e-30)
            res_rel = float(np.max(np.abs(r))) / scale
            hist.append({"phi": best[0] if accepted else None, "accepted": accepted,
                         "used": used if accepted else None, "Df_mean": float(np.mean(Df)),
                         "dT": dT if accepted else None, "res_rel": res_rel})
            if verbose:
                print(f"  [couple] it={it:3d} {used if accepted else 'REJECT':8s} dT={dT:.3e} K "
                      f"res_rel={res_rel:.3e} <Df>={np.mean(Df):.4g}")
            if not accepted:
                return T, {"converged": False, "reason": "no descent (retry limit)", "history": hist}
            if dT < tol_K and res_rel < tol_rel:
                return T, {"converged": True, "iters": it + 1, "dT": dT, "res_rel": res_rel, "history": hist}
        return T, {"converged": False, "reason": "max_iter", "dT": dT, "res_rel": res_rel, "history": hist}

    @staticmethod
    def _anderson(Ts, Gs, depth):
        """Anderson 加速 (type-II)。残差 F_k = G_k - T_k の最小二乗で混合係数を決める。"""
        m = min(depth, len(Ts) - 1)
        if m < 1:
            return None
        F = [g - t for t, g in zip(Ts, Gs)]
        dF = np.array([F[-1] - F[-1 - i] for i in range(1, m + 1)]).T      # (n, m)
        dG = np.array([Gs[-1] - Gs[-1 - i] for i in range(1, m + 1)]).T
        if not np.all(np.isfinite(dF)) or np.linalg.norm(dF) <= 0:
            return None
        try:
            gamma, *_ = np.linalg.lstsq(dF, F[-1], rcond=1e-12)
        except np.linalg.LinAlgError:
            return None
        if not np.all(np.isfinite(gamma)) or np.max(np.abs(gamma)) > 1e6:
            return None
        return Gs[-1] - dG @ gamma

    def driver(self, T0, Df0=None, anderson=5, Df_cap=1e12):
        """`FixedPointDriver` を作る (CFD のように $Q_f$ の評価が高価な場合)。

        `couple()` は候補ごとに $Q_f$ を評価する (安い場合向け)。CFD 連成では 1 反復に CFD 1 回しか
        使えないので、**棄却は「最後に良かった状態へ退避して $D_f$ を上げる」形**で行う。
        """
        return FixedPointDriver(self, T0, Df0, anderson, Df_cap)

    @staticmethod
    def _merit_M(M, r):
        """メリット関数 Φ = r^T M^{-1} r。**比較の間は M (= A_s + D_f) を固定する**。"""
        return float(r @ spla.spsolve(M.tocsc(), r))


    def driver(self, T0, Df0=None, anderson=5, Df_cap=1e12):
        """`FixedPointDriver` を作る。CFD のように $Q_f$ の評価が高価な場合はこちらを使う。

        `couple()` は候補ごとに $Q_f$ を評価する (安い場合向け)。CFD 連成では 1 反復に CFD 1 回しか
        使えないので、**棄却は「次の候補を素の反復に戻して $D_f$ を上げる」形**で行う。
        """
        return FixedPointDriver(self, T0, Df0, anderson, Df_cap)


class FixedPointDriver:
    r"""外部連成ループ用の固定点ドライバ (CFD 1 回/反復)。

    使い方:
        drv = op.driver(T0, Df0=...)
        Tw = drv.T
        while True:
            Qf = run_cfd_and_extract(Tw)        # 節点熱荷重 [W] (固体向き正)
            Tw, info = drv.advance(Qf)
            if info["converged"]: break

    受理: 固定重み $M=A_s+D_f$ での $\Phi=r^\mathsf{T}M^{-1}r$ が増えたら、その反復を**退避**して
    (最後に良かった $T$ に戻し) $D_f$ を倍にし、Anderson 履歴を捨てる。
    収束: max|ΔT| < tol_K かつ max|r| < tol_rel * スケール が `n_consec` 回連続。
    """

    def __init__(self, op: ShellOperator, T0, Df0=None, anderson=5, Df_cap=1e12):
        self.op = op
        self.T = np.asarray(T0, float).copy()
        self.Df = np.full(op.n, 1.0) if Df0 is None else np.asarray(Df0, float).copy()
        self.anderson = int(anderson)
        self.Df_cap = float(Df_cap)
        self.Ts, self.Gs = [], []
        self.best = None            # (phi, T)
        self.history = []
        self.it = 0
        self._ok_streak = 0

    def _assemble(self, T):
        """固体作用素を組む。`fem2d` のように内部節点を持つ作用素では、**内部温度を復元してから
        組み直す** (codex result M4, 2026-09-20)。これをしないと $k_s(T)$ が界面平均温度で
        全域一様に評価され、`solve()` の局所 Picard とは別の方程式を解くことになる
        (温度依存円環で 2.07 K ずれた状態を `converged` にしていた)。"""
        A, b = self.op.assemble(T)
        if hasattr(self.op, "recover_interior"):
            self.op.recover_interior(T)      # self.u を現在の T に整合させる
            A, b = self.op.assemble(T)       # 局所 k_s(T) で組み直す
        return A, b

    def advance(self, Qf, tol_K=1e-6, tol_rel=1e-6, n_consec=2,
                tol_abs_W=None, tol_solid=None):
        """最新の $Q_f(T_k)$ を受け取り、次の $T_{k+1}$ を返す。"""
        Qf = np.asarray(Qf, float)
        A, b = self._assemble(self.T)
        M = (A + sp.diags(self.Df)).tocsr()
        r = A @ self.T - b - Qf
        phi = ShellOperator._merit_M(M, r)
        # **規格化は $\max|Q_f|$ のみ** (codex result M5): 背面温度を含む大きな $b$ を混ぜると、
        # 物理的な不釣合いが 100 % でも res_rel が 1e-6 に見えて合格してしまう。
        scale = max(float(np.max(np.abs(Qf))), 1e-30)
        res_rel = float(np.max(np.abs(r))) / scale

        rejected = False
        if self.best is not None and phi > self.best[0]:
            # 退避: **最後に良かった状態の $T$ と $Q_f$ を組で戻す** (codex result M2)。
            # 旧実装は $T$ だけ戻して $Q_f$ は棄却された $T$ のものを使っており、
            # 異なる評価点を混ぜた更新になっていた (反例: 正 297.33 に対し 315.67 を返す)。
            self.T = self.best[1].copy()
            Qf = self.best[2].copy()
            self.Df = np.minimum(self.Df * 2.0, self.Df_cap)
            self.Ts.clear(); self.Gs.clear()
            self._ok_streak = 0
            rejected = True
            A, b = self._assemble(self.T)
            M = (A + sp.diags(self.Df)).tocsr()
            # **重みを変えたら基準メリットも新しい重みで測り直す** (混在比較の禁止)
            r_best = A @ self.T - b - Qf
            self.best = (ShellOperator._merit_M(M, r_best), self.T.copy(), Qf.copy())
        else:
            self.best = (phi, self.T.copy(), Qf.copy())

        G = ShellOperator._linsolve(M, b + Qf + self.Df * self.T, self.T)
        Tn = G
        used = "plain"
        if self.anderson > 0 and not rejected and len(self.Ts) >= 1:
            acc = ShellOperator._anderson(self.Ts + [self.T], self.Gs + [G], self.anderson)
            if acc is not None and np.all(np.isfinite(acc)):
                Tn, used = acc, "anderson"
        if not rejected:
            self.Ts.append(self.T.copy()); self.Gs.append(G.copy())
            if len(self.Ts) > self.anderson + 1:
                self.Ts.pop(0); self.Gs.pop(0)

        dT = float(np.max(np.abs(Tn - self.T)))
        # ---- 界面ゲート G-if (plan boundary-conjugate-heat-transfer §6、codex result M5) ----
        # **独立に満たすものを全部見る**。旧実装は dT と res_rel だけで、$D_f$ を倍にして
        # 更新が小さくなっただけの状態を「収束」にできた (反例: A_s=1000, b=3e5, Q_f=1 で
        # res_rel 3.3e-6 に見えて物理的な不釣合いは 100 %)。
        #   res_abs   : 界面残差の絶対値 [W] (面積あたりでなく節点荷重の単位)
        #   res_rel   : $\max|r| / \max|Q_f|$ (規格化は $Q_f$ のみ。b は混ぜない)
        #   res_solid : 固体内部の残差 (界面へ縮約した作用素と内部復元の整合)
        #   dT        : 温度更新の絶対値
        #   Df_grown  : この反復で $D_f$ を上げた (= 退避した) なら収束と認めない
        res_abs = float(np.max(np.abs(r)))
        res_solid = float("nan")
        if hasattr(self.op, "interior_residual"):
            try:
                res_solid = float(self.op.interior_residual(self.T))
            except Exception:
                res_solid = float("nan")
        ok_abs = True if tol_abs_W is None else (res_abs < tol_abs_W)
        ok_solid = True
        if tol_solid is not None and np.isfinite(res_solid):
            ok_solid = (res_solid < tol_solid)
        conv_now = ((dT < tol_K) and (res_rel < tol_rel) and ok_abs and ok_solid
                    and not rejected)
        self._ok_streak = self._ok_streak + 1 if conv_now else 0
        info = {"iter": self.it, "phi": phi, "res_rel": res_rel, "dT": dT, "used": used,
                "res_abs": res_abs, "res_solid": res_solid,
                "rejected": rejected, "Df_mean": float(np.mean(self.Df)),
                "converged": self._ok_streak >= n_consec}
        self.history.append(info)
        self.T = Tn
        self.it += 1
        return self.T, info


