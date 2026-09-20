#!/usr/bin/env python3
r"""solid_fem2d.py の単体検証 (解析解との照合)。plan §4.4d の「単体検証」項目。

円環 (内半径 a に Robin の冷却孔、外半径 b が CHT 界面) の 1 次元 (半径方向) 伝導:

  Q = q * 2 pi b                 … 外周から入る全熱量 (単位奥行き, W/m)
  T(a) = T_c + Q / (2 pi a h)    … 内周の Robin
  T(r) = T(a) + Q/(2 pi k) ln(r/a)

  T1 円環: T(b) と T(a) を解析解と照合 + 格子収束 (2 次)
  T2 Schur: 界面へ縮約した解と全系直接解が一致する (機械精度)
  T3 非一様荷重: 周方向に変わる q でも Schur 経由の解が全系解と一致する
  T4 連成: 非対角な流体応答に対して driver が厳密固定点へ収束する

usage: python3 solver_density_cuda/tools/test_solid_fem2d.py
"""
import math
import os
import sys

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from solid_fem2d import Fem2DOperator   # noqa: E402

FAILS = []


def check(name, ok, detail=""):
    print(("PASS  " if ok else "FAIL  ") + name + ("   " + detail if detail else ""))
    if not ok:
        FAILS.append(name)


def annulus(a, b, nr, nt):
    """円環の構造三角形メッシュ。戻り値: nodes, tris, outer_nodes, outer_edges, inner_edges"""
    rs = np.linspace(a, b, nr + 1)
    th = np.linspace(0.0, 2.0 * math.pi, nt, endpoint=False)
    nodes = np.array([[r * math.cos(t), r * math.sin(t)] for r in rs for t in th])
    idx = lambda i, j: i * nt + (j % nt)          # noqa: E731
    tris = []
    for i in range(nr):
        for j in range(nt):
            n00, n01, n10, n11 = idx(i, j), idx(i, j + 1), idx(i + 1, j), idx(i + 1, j + 1)
            tris.append([n00, n10, n11])
            tris.append([n00, n11, n01])
    outer = [idx(nr, j) for j in range(nt)]
    outer_edges = [(idx(nr, j), idx(nr, j + 1)) for j in range(nt)]
    inner_edges = [(idx(0, j), idx(0, j + 1)) for j in range(nt)]
    return nodes, np.array(tris), np.array(outer), outer_edges, inner_edges


def build(a, b, nr, nt, k, h, Tc):
    nodes, tris, outer, oe, ie = annulus(a, b, nr, nt)
    robin = [(n0, n1, h, Tc) for (n0, n1) in ie]
    return Fem2DOperator(nodes, tris, outer, oe, robin, k), nodes


def t1_annulus():
    a, b, k, h, Tc, q = 0.002, 0.010, 20.0, 5000.0, 400.0, 5.0e4
    Q = q * 2 * math.pi * b
    Ta = Tc + Q / (2 * math.pi * a * h)
    Tb = Ta + Q / (2 * math.pi * k) * math.log(b / a)
    errs = []
    for (nr, nt) in ((4, 24), (8, 48), (16, 96)):
        op, nodes = build(a, b, nr, nt, k, h, Tc)
        T = op.solve(q * op.area)
        u = op.recover_interior(T)
        r = np.linalg.norm(nodes, axis=1)
        inner = u[r < a * 1.0001]
        errs.append(max(abs(float(np.mean(T)) - Tb), abs(float(np.mean(inner)) - Ta)))
        last = (float(np.mean(T)), float(np.mean(inner)))
    rates = [math.log2(errs[i] / errs[i + 1]) for i in range(len(errs) - 1)]
    check("T1 annulus with a Robin hole vs analytic",
          errs[-1] < 0.02 * (Tb - Tc) and min(rates) > 1.6,
          f"T(b)={last[0]:.3f} (exact {Tb:.3f}) , T(a)={last[1]:.3f} (exact {Ta:.3f}) ; "
          f"errors " + " ".join(f"{e:.3e}" for e in errs) + " ; rates " + " ".join(f"{r:.2f}" for r in rates))


def t2_schur_equals_full():
    a, b, k, h, Tc, q = 0.002, 0.010, 20.0, 5000.0, 400.0, 5.0e4
    op, nodes = build(a, b, 8, 48, k, h, Tc)
    Qf = q * op.area
    T_schur = op.solve(Qf)
    K, bvec = op.assemble_full(np.full(op.N, 400.0))
    rhs = bvec.copy()
    rhs[op.iface] += Qf
    u_full = spla.spsolve(K.tocsc(), rhs)
    err = float(np.max(np.abs(T_schur - u_full[op.iface])))
    check("T2 Schur reduction reproduces the full solve",
          err < 1e-8 * max(1.0, float(np.max(np.abs(u_full)))), f"max|T_schur - T_full| = {err:.3e} K")


def t3_nonuniform_load():
    a, b, k, h, Tc = 0.002, 0.010, 20.0, 5000.0, 400.0
    op, nodes = build(a, b, 8, 48, k, h, Tc)
    th = np.arctan2(op.coords[:, 1], op.coords[:, 0])
    q = 5.0e4 * (1.0 + 0.7 * np.cos(2 * th))          # 周方向に変わる熱流束
    Qf = q * op.area
    T_schur = op.solve(Qf)
    K, bvec = op.assemble_full(np.full(op.N, 400.0))
    rhs = bvec.copy(); rhs[op.iface] += Qf
    u_full = spla.spsolve(K.tocsc(), rhs)
    err = float(np.max(np.abs(T_schur - u_full[op.iface])))
    check("T3 non-uniform interface load", err < 1e-8 * float(np.max(np.abs(u_full))),
          f"T(b) spread {T_schur.max()-T_schur.min():.2f} K ; max|T_schur - T_full| = {err:.3e} K")


def t4_coupling():
    a, b, k, h, Tc = 0.002, 0.010, 20.0, 5000.0, 400.0
    op, nodes = build(a, b, 8, 48, k, h, Tc)
    n = op.n
    Taw = 1600.0
    hf = 800.0                                        # 流体側の膜係数 [W/m2K]
    off = np.full(n - 1, -0.15 * hf)
    H = (sp.diags([off, np.full(n, hf), off], [-1, 0, 1]).tocsr()) * sp.diags(op.area).tocsr()
    Qf = lambda T: np.asarray(H @ (Taw - np.asarray(T, float))).ravel()   # noqa: E731
    S, bi = op.assemble(np.full(n, 600.0))
    Texact = spla.spsolve((S + H).tocsc(), bi + np.asarray(H @ np.full(n, Taw)).ravel())
    drv = op.driver(np.full(n, 600.0), Df0=np.full(n, hf) * op.area)
    T = drv.T
    info = {}
    for _ in range(200):
        T, info = drv.advance(Qf(T), tol_K=1e-7, tol_rel=1e-7)
        if info["converged"]:
            break
    err = float(np.max(np.abs(T - Texact)))
    check("T4 coupling on the 2D solid reaches the exact fixed point",
          info.get("converged", False) and err < 1e-5 * float(np.max(np.abs(Texact - Tc))),
          f"converged={info.get('converged')} iters={info.get('iter')} max|T-exact|={err:.3e} K "
          f"(T range {T.min():.1f}..{T.max():.1f} K)")


def t5_local_k_of_T():
    r"""**温度依存 $k_s(T)$ が局所で解かれているか** (codex result M4, 2026-09-20)。

    外部ドライバは `assemble()` しか呼ばず `recover_interior()` を呼ばなかったため、
    `self.u` が `None` のまま = 全固体節点が**界面平均温度**とみなされ、$k_s$ が全域一様に
    評価されていた。`solve()` の局所 Picard とは別の方程式であり、codex の再現では
    壁温が 2.07 K ずれた状態を `converged=True` にしていた。
    """
    a, b, h, Tc = 0.002, 0.010, 5000.0, 400.0
    kT = (np.array([300.0, 600.0, 900.0, 1200.0]), np.array([12.0, 18.0, 26.0, 36.0]))
    op, nodes = build(a, b, 8, 48, kT, h, Tc)
    n = op.n
    Qf = np.full(n, 0.0)
    # 界面に一様熱荷重 (単位奥行き)
    Qf = 9.0e4 * op.area
    T_full = op.solve(Qf, T0=np.full(n, 600.0), iters=60, tol=1e-13)
    # **ドライバ側は必ず新しい作用素で始める**。同じ op を使い回すと `solve()` が残した
    # `self.u` が初期値になり、「内部温度を復元していない」欠陥が隠れる。
    op2, _ = build(a, b, 8, 48, kT, h, Tc)
    drv = op2.driver(np.full(n, 600.0), Df0=np.full(n, h) * op2.area)
    T = drv.T; info = {}
    for _ in range(300):
        T, info = drv.advance(Qf, tol_K=1e-8, tol_rel=1e-8)
        if info["converged"]:
            break
    err = float(np.max(np.abs(T - T_full)))
    check("T5 driver solves the same k_s(T) equation as the full Picard solve",
          info.get("converged", False) and err < 1e-3,
          f"driver {T.mean():.4f} K vs full {T_full.mean():.4f} K, max|diff| {err:.3e} K "
          f"(修正前は 2.07 K ずれて converged=True だった)")



if __name__ == "__main__":
    t1_annulus(); t2_schur_equals_full(); t3_nonuniform_load(); t4_coupling()
    print()
    if FAILS:
        print("VERDICT: FAIL (%d)" % len(FAILS))
        for f in FAILS:
            print("  -", f)
        sys.exit(1)
    print("VERDICT: PASS (all)")
