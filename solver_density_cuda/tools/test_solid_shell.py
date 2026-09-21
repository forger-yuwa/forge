#!/usr/bin/env python3
r"""solid_shell.py の単体検証 (解析解との照合)。

AGENTS.md / plans/active/boundary-conjugate-heat-transfer.md §6 の
「格子収束で連続解析解への接近を調べる試験」と「離散残差の機械精度検査」を**分けて**実施する。

  T1 直列抵抗 (機械精度)   : 一様 q・背面等温 → T_w = T_b + q (t/k_s + R_back)。local1d / shell2d 両方。
  T2 フィン (格子収束)     : (k t) T'' - (T-T_b)/R = -q(x), q = q0 cos(pi x/L), 端部断熱
                             → T = T_b + (q0/(k t)) cos(pi x/L) / (1/(k t R) + (pi/L)^2)。O(h^2) を確認。
  T3 零空間の検出          : 背面断熱 + 正味入熱 != 0 は**解が無い**ので例外にする。
  T4 温度依存 k_s          : k_s(T) テーブルで T_w = T_b + q R_tot(T_w) の非線形解に一致。
  T5 連成反復の受理判定    : 非対角な流体応答 H (SPD) に対し、メリット関数 + line search で収束する。
                             **最大ノルムの単調減少を要求すると棄却してしまう反例**を含む。

usage: python3 solver_density_cuda/tools/test_solid_shell.py
"""
import math
import os
import sys

import numpy as np
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from solid_shell import ShellOperator, SolidModel, FixedPointDriver   # noqa: E402

FAILS = []


def check(name, ok, detail=""):
    print(("PASS  " if ok else "FAIL  ") + name + ("   " + detail if detail else ""))
    if not ok:
        FAILS.append(name)


def strip_mesh(nx, L=0.1, w=0.01):
    """x 方向 nx+1 節点・y 方向 2 節点の帯 (1 列の四角形)。"""
    xs = np.linspace(0.0, L, nx + 1)
    coords = np.array([[x, y, 0.0] for x in xs for y in (0.0, w)])
    faces = [[2 * i, 2 * i + 2, 2 * i + 3, 2 * i + 1] for i in range(nx)]
    return coords, faces, xs


def t1_series_resistance():
    coords, faces, xs = strip_mesh(8)
    for mode in ("local1d", "shell2d"):
        for back in ({"kind": "isothermal"}, {"kind": "coolant", "h": 5000.0}):
            m = SolidModel(mode=mode, thickness=4.0e-4, k_solid=15.0,
                           back_kind=back["kind"], T_b=300.0, h_c=back.get("h", 0.0))
            op = ShellOperator(coords, faces, m)
            q = 1.0e5                                    # [W/m2] 一様
            T = op.solve(q * op.area)                    # 節点荷重 [W]
            R = 4.0e-4 / 15.0 + (1.0 / 5000.0 if back["kind"] == "coolant" else 0.0)
            exact = 300.0 + q * R
            err = float(np.max(np.abs(T - exact)))
            check(f"T1 series resistance [{mode}, {back['kind']}]",
                  err < 1e-9 * exact, f"max|T-exact| = {err:.3e} K (exact {exact:.4f} K)")


def t2_fin_convergence():
    L, w = 0.1, 0.01
    t, ks, hc, Tb, q0 = 4.0e-4, 15.0, 20.0, 300.0, 1.0e5   # h_c は緩く (m h << 1 で 2 次を見る)
    R = t / ks + 1.0 / hc
    kt = ks * t
    denom = 1.0 / (kt * R) + (math.pi / L) ** 2

    def exact(x):
        return Tb + (q0 / kt) * np.cos(math.pi * x / L) / denom

    errs = []
    for nx in (20, 40, 80, 160):
        coords, faces, xs = strip_mesh(nx, L, w)
        m = SolidModel(mode="shell2d", thickness=t, k_solid=ks, back_kind="coolant", T_b=Tb, h_c=hc)
        op = ShellOperator(coords, faces, m)
        qn = q0 * np.cos(math.pi * coords[:, 0] / L) * op.area      # 節点荷重 [W]
        T = op.solve(qn)
        e = float(np.max(np.abs(T - exact(coords[:, 0]))))
        errs.append(e)
    rates = [math.log2(errs[i] / errs[i + 1]) for i in range(len(errs) - 1)]
    check("T2 fin: grid convergence 2nd order",
          min(rates) > 1.8, "errors " + " ".join(f"{e:.3e}" for e in errs)
          + " | rates " + " ".join(f"{r:.2f}" for r in rates))


def t3_nullspace():
    coords, faces, _ = strip_mesh(8)
    m = SolidModel(mode="shell2d", thickness=4.0e-4, k_solid=15.0, back_kind="adiabatic")
    op = ShellOperator(coords, faces, m)
    try:
        op.solve(1.0e5 * op.area)
        check("T3 adiabatic isolated solid is rejected", False, "no exception raised")
    except ValueError as e:
        check("T3 adiabatic isolated solid is rejected", "singular" in str(e), str(e)[:70])


def t4_temperature_dependent_k():
    coords, faces, _ = strip_mesh(4)
    Tt = np.array([300.0, 1000.0]); kt = np.array([10.0, 30.0])
    m = SolidModel(mode="local1d", thickness=1.0e-3, k_solid=10.0, k_table=(Tt, kt),
                   back_kind="isothermal", T_b=300.0)
    op = ShellOperator(coords, faces, m)
    q = 2.0e5
    T = op.solve(q * op.area)
    # 非線形方程式 T = Tb + q (t/k(T)) を独立に解く (固定点)
    Tx = 300.0
    for _ in range(200):
        Tx = 300.0 + q * (1.0e-3 / float(np.interp(Tx, Tt, kt)))
    err = float(np.max(np.abs(T - Tx)))
    check("T4 temperature dependent k_s", err < 1e-6 * Tx, f"max|T-ref| = {err:.3e} K (ref {Tx:.4f} K)")


def t5_coupling_acceptance():
    """非対角な流体応答に対する受理判定。codex 3 巡目 #1 の反例を含む。"""
    coords, faces, _ = strip_mesh(6)
    m = SolidModel(mode="shell2d", thickness=1.0e-3, k_solid=20.0, back_kind="coolant", h_c=1000.0, T_b=300.0)
    op = ShellOperator(coords, faces, m)
    n = op.n
    # 流体応答: Q_f(T) = H (T_aw - T), H は非対角 SPD (隣接節点を結ぶ伝導網)
    Taw = 1200.0
    main = np.full(n, 2.0); off = np.full(n - 1, -1.0)
    H = sp.diags([off, main, off], [-1, 0, 1]).tocsr() * 50.0
    Qf = lambda T: H @ (Taw - np.asarray(T, float))      # noqa: E731
    # 厳密解: (A_s + H) T = b_s + H Taw   (A_s, b_s は温度非依存)
    A, b = op.assemble(np.full(n, 300.0))
    import scipy.sparse.linalg as spla
    Texact = spla.spsolve((A + H).tocsc(), b + H @ np.full(n, Taw))
    T, info = op.couple(Qf, np.full(n, 300.0), Df0=np.full(n, 1.0),
                        max_iter=200, tol_K=1e-7, tol_rel=1e-7, anderson=5)
    err = float(np.max(np.abs(T - Texact)))
    check("T5 coupling reaches the exact fixed point (non-diagonal fluid response)",
          info.get("converged", False) and err < 1e-5 * np.max(np.abs(Texact - 300.0)),
          f"converged={info.get('converged')} iters={info.get('iters')} max|T-exact|={err:.3e} K")
    # 素の固定点反復 (Anderson なし) ではスカラー D_f が非対角応答を捌けず止まる = 加速が必要なことの記録
    T0_, info0 = op.couple(Qf, np.full(n, 300.0), Df0=np.full(n, 1.0),
                           max_iter=200, tol_K=1e-7, tol_rel=1e-7, anderson=0)
    check("T5c plain fixed point is the slow one (Anderson is needed, not cosmetic)",
          not info0.get("converged", False) or info0.get("iters", 1e9) > info.get("iters", 0),
          f"plain: converged={info0.get('converged')} iters={info0.get('iters')} "
          f"max|T-exact|={np.max(np.abs(T0_-Texact)):.3e} K")
    # 反例: 最大ノルムの単調減少を受理条件にすると全候補が棄却される構成でも、
    #       メリット関数なら降下する (plan §4.2)。
    As = sp.diags([0.1, 0.1]).tocsr()
    Hs = sp.csr_matrix(np.array([[2.0, -1.0], [-1.0, 2.0]]))
    Df = np.array([1.0, 10.0])
    r = np.array([1.0, 1.0])
    M = (As + sp.diags(Df)).tocsc()
    rtrial = r - (As + Hs) @ spla.spsolve(M, r)
    phi0 = float(r @ spla.spsolve(M, r)); phi1 = float(rtrial @ spla.spsolve(M, rtrial))
    check("T5b max-norm rejects a convergent step, merit function accepts it",
          np.max(np.abs(rtrial)) > np.max(np.abs(r)) and phi1 < phi0,
          f"max-norm {np.max(np.abs(r)):.3f} -> {np.max(np.abs(rtrial)):.3f}, Phi {phi0:.6f} -> {phi1:.6f}")


def test_gif_normalisation():
    r"""**規格化は $\max|Q_f|$ のみ** (codex result M5 の反例, 2026-09-20)。

    $A_s=1000$ W/K, $b=300000$ W, $T=300.002$ K, $Q_f=1$ W とすると残差は $|r|=1$ W。
    旧実装は $\max(|Q_f|,|b|)$ で割っていたので `res_rel` が **3.33e-6** に見え、
    物理的な不釣合いが **100 %** でも 2 反復で `converged=True` になっていた。
    """
    class Op1:
        n = 1
        def assemble(self, T):
            return sp.csr_matrix(np.array([[1000.0]])), np.array([300000.0])
    drv = FixedPointDriver(Op1(), np.array([300.002]), Df0=np.array([1.0e9]), anderson=0)
    _, info = drv.advance(np.array([1.0]), tol_K=1e-3, tol_rel=1e-3, n_consec=1)
    ok = (info["res_rel"] > 0.5) and (not info["converged"])
    check("T7 G-if normalises by max|Qf| only (not b)", ok,
          f"res_rel={info['res_rel']:.4g} (旧実装 3.33e-06), res_abs={info['res_abs']:.4g} W, "
          f"converged={info['converged']}")
    return ok


def test_reject_keeps_state_pair():
    """**棄却時に $T$ と $Q_f$ を組で戻す** (codex result M2, 2026-09-20 の反例)。

    $A_s=1,\\ b_s=0,\\ Q(T)=3300-10T,\\ T_0=301,\\ D=1$。旧実装は $T$ だけ `best` に戻し、
    $Q_f$ は棄却された $T$ で評価したものを使っていたため 315.667 K を返していた。
    同じ評価点の組で更新すれば 297.333 K になる。
    """
    class Op1:
        n = 1
        def assemble(self, T):
            return sp.csr_matrix(np.array([[1.0]])), np.array([0.0])
    op = Op1()
    Q = lambda T: 3300.0 - 10.0 * T
    drv = FixedPointDriver(op, np.array([301.0]), Df0=np.array([1.0]), anderson=0)
    T1, _ = drv.advance(np.array([Q(301.0)]))
    T2, info = drv.advance(np.array([Q(float(T1[0]))]))
    ok = abs(float(T2[0]) - 297.3333333) < 1e-6
    check("T6 rejection restores (T, Qf) as one state", ok,
          f"first={float(T1[0]):.4f} -> after reject {float(T2[0]):.4f} K "
          f"(正 297.3333, 旧実装 315.6667; rejected={info['rejected']})")
    return ok


def test_noise_aware_acceptance():
    r"""**ノイズを知った受理** (plan boundary-conjugate-heat-transfer §5.1 #63)。

    固体は節点ごとに独立な抵抗 ($A_s$=50 W/K, 背面 300 K)、流体は $Q(T)=h(T_g-T)$ ($h$=200 W/K, $T_g$=800 K)
    に反復ごとの乱数 (σ=40 W) を足す。不動点は $T^*=(50\cdot300+200\cdot800)/250=700$ K。
    旧ロジック (σ を渡さない) は不動点の手前でノイズに負けて棄却を続け、$D_f$ が倍々に増えて凍る。
    新ロジックはノイズ幅 (σ/(A_s+h)=0.16 K の数倍) 以内に着き、`converged` を返す。
    """
    n, As, Tb, h, Tg, sig = 40, 50.0, 300.0, 200.0, 800.0, 40.0
    Df0 = 20.0 * h          # 安全側に大きい感度 = 1 反復の縮小率 0.94。ノイズに負けやすい実機の状況
    class OpN:
        def __init__(self): self.n = n
        def assemble(self, T):
            return sp.identity(n, format="csr") * As, np.full(n, As * Tb)
    Tstar = (As * Tb + h * Tg) / (As + h)
    out = {}
    for mode in ("old", "new", "understated"):
        rng = np.random.default_rng(7)
        drv = FixedPointDriver(OpN(), np.full(n, 400.0), Df0=np.full(n, Df0), anderson=0)
        T = drv.T; conv_it = None
        for it in range(400):
            Qf = h * (Tg - T) + rng.normal(0.0, sig, n)
            # "understated": 呼び出し側の σ が実際の 1/20 (実機と同じ状況)。回し直しのずれから学習できること
            kw = {} if mode == "old" else {"Qf_sigma": np.full(n, sig if mode == "new" else sig / 20.0)}
            T, info = drv.advance(Qf, tol_K=(1e-3 if mode == "old" else 5e-3), tol_rel=1e-3, n_consec=(3 if mode == "old" else 8), **kw)
            if info["converged"]:
                conv_it = it; break
        out[mode] = (float(np.abs(T - Tstar).max()), float(np.mean(drv.Df)), conv_it)
    err_old, df_old, _ = out["old"]; err_new, df_new, it_new = out["new"]
    ok = (err_new < 5 * sig / (As + h)) and (it_new is not None) and (df_new == Df0) and (df_old > 100 * Df0) and (err_old > 3 * err_new)
    err_u, df_u, it_u = out["understated"]
    ok_u = (err_u < 5 * sig / (As + h)) and (df_u < 64 * Df0)
    check("T9 understated sigma: reproducibility learned from re-evaluation, no freeze", ok_u,
          f"max|T-T*|={err_u:.3f} K, D_f grown x{df_u/Df0:.0f} (旧ロジックは x{df_old/Df0:.0f}), converged at it {it_u}")
    check("T8 noise-aware acceptance reaches the fixed point; the old rule freezes", ok,
          f"new: max|T-T*|={err_new:.3f} K, D_f unchanged, converged at it {it_new} | "
          f"old: max|T-T*|={err_old:.3f} K with D_f grown x{df_old/Df0:.0f}")
    return ok and ok_u


if __name__ == "__main__":
    t1_series_resistance()
    t2_fin_convergence()
    t3_nullspace()
    t4_temperature_dependent_k()
    t5_coupling_acceptance()
    test_reject_keeps_state_pair()
    test_gif_normalisation()
    test_noise_aware_acceptance()
    print()
    if FAILS:
        print("VERDICT: FAIL (%d)" % len(FAILS))
        for f in FAILS:
            print("  -", f)
        sys.exit(1)
    print("VERDICT: PASS (all)")


