#!/usr/bin/env python3
r"""軸対称の固体 FE (`fem2d`, $r=y$ の重み) を検証する。

plan boundary-cht-axisymmetric-fem2d §6 V-ax1 (事前登録の合格ライン)。C++ (`conjugate/solidFem2d.cpp`、
`solid_fem2d_tool --axisym`) と参照オラクル (`tools/solid_fem2d.py` の `Fem2DOperator(axisym=True)`) の両方を見る。

  (a) C++ ↔ Python: 軸対称で $K$・$b$ の相対差 ≤1e-12、同一荷重での解 ≤1e-9 K
  (b) 独立な辺積分: 半径が変わる辺 3 通り (半径方向・斜め・軸の近く) で Robin 行列・荷重・$A_i^r$・`q_hole` を
      **この評価器の中で独立に取った Gauss–Legendre 積分** ($\int h N_iN_j r\,ds$ 等、5 点 = 9 次まで厳密) と照合、
      相対 ≤1e-13。C++↔Python の一致だけでは共通の誤りを排除できないため
  (c) 厚肉円筒殻: 内面 $r_1$ に一様熱流束 $q_1$ (荷重 $Q_i=q_1A_i^r$)、外面 $r_2=2r_1$ に Robin。
      $T(r)=T_c+q_1r_1[1/(hr_2)+\ln(r_2/r)/k]$。半径方向 4/8/16 層で収束率 ≥1.8、16 層で内面温度誤差
      ≤0.1 % of 殻の温度降下。**平面のまま組むと FAIL すること** (検出力) も確認する。
      格子は軸方向も同率で細かくする (縦横比固定)。軸方向固定の結果は参考表として出す (判定外)
  (d) 軸対称円板: $x\in[0,t]$, $r\in[r_1,4r_1]$、界面 $x=0$ に一様熱流束、Robin $x=t$。解は $x$ の 1 次関数
      (線形要素で厳密) → 全節点温度誤差 ≤1e-9 K、`q_hole` 節点値 = consistent 積分の厳密値 (相対 ≤1e-12)
  (e) 収支: 界面入熱 = Robin 持ち去り (`q_hole` の総和) 相対 ≤1e-12、(c)(d) とも

`--base-exe` に変更前のツールを渡すと、**平面 (`--axisym` なし)** の matvec / solve 出力が
変更後のツールと**バイト一致**することも確認する (V-ax1 追加項目: 平面のビット同一)。

使い方:
  python3 solver_density_cuda/tools/test_solid_fem2d_axisym.py \
      --exe solver_density_cuda/build-ypls/solid_fem2d_tool \
      [--base-exe /path/to/solid_fem2d_tool_before] \
      [--planar-h5 case/58.conjugate_slot/mesh/solid_front_tc300_nl16.h5 ...]
"""
from __future__ import annotations

import argparse
import math
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from solid_fem2d import Fem2DOperator  # noqa: E402
from solid_mesh_to_h5 import write_solid_h5  # noqa: E402
from test_solid_mesh_to_h5 import solve_all  # noqa: E402

FAILS: list[str] = []
ENV = dict(os.environ)
ENV.setdefault("LD_LIBRARY_PATH", "/usr/lib/x86_64-linux-gnu/hdf5/serial")


def check(name: str, ok: bool, detail: str = ""):
    print(f"  {'PASS' if ok else 'FAIL'}  {name}" + (f"   {detail}" if detail else ""))
    if not ok:
        FAILS.append(name)


# ---------------------------------------------------------------------------
# ツール呼び出し
# ---------------------------------------------------------------------------
def tool(exe, mode, h5, inp=None, axisym=True, td=None, raw=False):
    """solid_fem2d_tool を回して出力ファイルのパス (raw) か読んだ値を返す。"""
    td = Path(td)
    fi, fo = td / f"in_{mode}.txt", td / f"out_{mode}.txt"
    args = [exe] + (["--axisym"] if axisym else []) + [mode, str(h5)]
    if mode == "lumped":
        args += [str(fo)]
    elif mode == "field":
        np.savetxt(fi, np.asarray(inp, float), fmt="%.17e")
        fo = td / "field"
        args += [str(fi), str(fo)]
    else:
        np.savetxt(fi, np.asarray(inp, float), fmt="%.17e")
        args += [str(fi), str(fo)]
    r = subprocess.run(args, env=ENV, capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"[test_solid_fem2d_axisym] {' '.join(args)} が落ちた:\n{r.stdout}\n{r.stderr}")
    if mode == "field":
        with h5py.File(str(fo) + ".h5", "r") as f:
            return np.asarray(f["VALUE/q_hole"][:], float)
    if raw:
        return fo
    if mode == "matrix":
        return read_matrix(fo)
    return np.atleast_1d(np.loadtxt(fo))


def read_matrix(path):
    """`matrix` 出力 (下三角の "i j a" と "b i b_i") を dict と配列にする。"""
    K, b = {}, {}
    for line in Path(path).read_text().splitlines():
        t = line.split()
        if t[0] == "b":
            b[int(t[1])] = float(t[2])
        else:
            K[(int(t[0]), int(t[1]))] = float(t[2])
    n = len(b)
    return K, np.array([b[i] for i in range(n)])


def op_from_h5(path, axisym):
    with h5py.File(path, "r") as f:
        nodes = np.asarray(f["MESH/COORD"][:], float)
        tris = np.asarray(f["MESH/TRIS"][:], int)
        iface = np.asarray(f["IFACE/NODES"][:], int)
        ie = np.asarray(f["IFACE/EDGES"][:], int)
        re_, rh, rtc = f["ROBIN/EDGES"][:], f["ROBIN/H"][:], f["ROBIN/TC"][:]
        kT, kV = np.asarray(f["SOLID/K_T"][:], float), np.asarray(f["SOLID/K_V"][:], float)
    robin = [(int(a), int(b), float(h), float(t)) for (a, b), h, t in zip(re_, rh, rtc)]
    k_solid = (kT, kV) if len(kT) > 1 else float(kV[0])
    return Fem2DOperator(nodes, tris, iface, [tuple(e) for e in ie], robin, k_solid, axisym=axisym)


def lower_dict(Ksp):
    K = Ksp.tocoo()
    out = {}
    for i, j, v in zip(K.row, K.col, K.data):
        if i >= j and v != 0.0:
            out[(int(i), int(j))] = out.get((int(i), int(j)), 0.0) + float(v)
    return out


def rel_dict(A, B):
    keys = set(A) | set(B)
    d = max(abs(A.get(k, 0.0) - B.get(k, 0.0)) for k in keys)
    s = max(abs(v) for v in B.values())
    return d / s


# ---------------------------------------------------------------------------
# メッシュ
# ---------------------------------------------------------------------------
def rect_mesh(x0, x1, y0, y1, nx, ny, jitter=0.0, seed=0, diag="alt"):
    """矩形 (x: nx 分割, y: ny 分割) の三角形分割。境界の辺の節点 index を返す。

    diag="alt" は対角の向きを市松に交互、"uniform" は全セル同じ向き。"""
    xs, ys = np.linspace(x0, x1, nx + 1), np.linspace(y0, y1, ny + 1)
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    nodes = np.column_stack([X.ravel(), Y.ravel()])
    idx = lambda i, j: i * (ny + 1) + j          # noqa: E731
    if jitter > 0:
        rng = np.random.default_rng(seed)
        for i in range(1, nx):
            for j in range(1, ny):
                nodes[idx(i, j), 0] += jitter * (x1 - x0) / nx * (rng.random() - 0.5)
                nodes[idx(i, j), 1] += jitter * (y1 - y0) / ny * (rng.random() - 0.5)
    tris = []
    for i in range(nx):
        for j in range(ny):
            a, b, c, d = idx(i, j), idx(i + 1, j), idx(i + 1, j + 1), idx(i, j + 1)
            if diag == "uniform" or (i + j) % 2 == 0:   # 既定は対角の向きを交互にして偏りを消す
                tris += [[a, b, c], [a, c, d]]
            else:
                tris += [[a, b, d], [b, c, d]]
    edges = {
        "y0": [(idx(i, 0), idx(i + 1, 0)) for i in range(nx)],
        "y1": [(idx(i, ny), idx(i + 1, ny)) for i in range(nx)],
        "x0": [(idx(0, j), idx(0, j + 1)) for j in range(ny)],
        "x1": [(idx(nx, j), idx(nx, j + 1)) for j in range(ny)],
    }
    return nodes, np.array(tris), edges


# ---------------------------------------------------------------------------
# (b) 独立な辺積分 (Gauss–Legendre 5 点)
# ---------------------------------------------------------------------------
GX, GW = np.polynomial.legendre.leggauss(5)


def edge_quad(pa, pb, fn):
    """辺 pa→pb 上で ∫ fn(Na, Nb, r) ds を Gauss–Legendre で取る (fn はベクトル/行列を返してよい)。"""
    L = math.hypot(pb[0] - pa[0], pb[1] - pa[1])
    acc = 0.0
    for xi, w in zip(GX, GW):
        s = 0.5 * (xi + 1.0)                  # 0..1
        Na, Nb = 1.0 - s, s
        r = Na * pa[1] + Nb * pb[1]           # r = y (座標を辺上で線形補間)
        acc = acc + 0.5 * w * L * np.asarray(fn(Na, Nb, r), float)
    return acc


def test_b(exe, td):
    print("\n(b) 独立な辺積分 (Gauss–Legendre 5 点との照合、相対 ≤1e-13)")
    h, Tc = 1234.5, 321.0
    cases = {
        "半径方向 r_a=1,r_b=2": ((0.0, 1.0), (0.0, 2.0)),
        "斜め (0,0.5)->(0.8,1.7)": ((0.0, 0.5), (0.8, 1.7)),
        "軸の近く (0,1e-4)->(0.3,0.02)": ((0.0, 1.0e-4), (0.3, 0.02)),
    }
    worst = 0.0
    for name, (pa, pb) in cases.items():
        # 1 三角形。Robin 辺 = 界面辺 = (a,b)。k=0 で剛性を消し、組んだ行列を Robin 辺そのものにする
        pc = (pa[0] + 0.7 * (pb[1] - pa[1]) + 0.1, 0.5 * (pa[1] + pb[1]) + 0.05)
        nodes = np.array([pa, pb, pc])
        h5p = Path(td) / "edge.h5"
        write_solid_h5(h5p, nodes, [[0, 1, 2]], [(0, 1)], [(0, 1)], [h], [Tc], [300.0], [0.0],
                       rcm=False, source="edge")
        Kc, bc = tool(exe, "matrix", h5p, np.full(3, 300.0), td=td)
        Ac = tool(exe, "lumped", h5p, td=td)
        u = np.array([450.0, 700.0, 500.0])
        qh = tool(exe, "field", h5p, u, td=td)

        M_ref = edge_quad(pa, pb, lambda Na, Nb, r: h * r * np.outer([Na, Nb], [Na, Nb]))
        b_ref = edge_quad(pa, pb, lambda Na, Nb, r: h * Tc * r * np.array([Na, Nb]))
        A_ref = edge_quad(pa, pb, lambda Na, Nb, r: r * np.array([Na, Nb]))
        q_ref = edge_quad(pa, pb, lambda Na, Nb, r: h * r * (Na * u[0] + Nb * u[1] - Tc) * np.array([Na, Nb]))

        M_cpp = np.array([[Kc.get((0, 0), 0.0), Kc.get((1, 0), 0.0)],
                          [Kc.get((1, 0), 0.0), Kc.get((1, 1), 0.0)]])
        # 界面の並び (IFACE/NODES = 辺の節点を昇順) は (0,1)
        rels = {
            "M": np.max(np.abs(M_cpp - M_ref) / np.abs(M_ref)),
            "b": np.max(np.abs(bc[:2] - b_ref) / np.abs(b_ref)),
            "A^r": np.max(np.abs(Ac - A_ref) / np.abs(A_ref)),
            "q_hole": np.max(np.abs(qh[:2] - q_ref) / np.abs(q_ref)),
        }
        extra_ok = (Kc.get((2, 2), 0.0) == 0.0 and bc[2] == 0.0 and qh[2] == 0.0)
        # Python オラクルも同じ参照と照合する
        op = op_from_h5(h5p, axisym=True)
        Kp, bp = op.assemble_full(np.full(3, 300.0))
        Kp = Kp.toarray()
        rels_py = {
            "M": np.max(np.abs(Kp[:2, :2] - M_ref) / np.abs(M_ref)),
            "b": np.max(np.abs(bp[:2] - b_ref) / np.abs(b_ref)),
            "A^r": np.max(np.abs(op.area - A_ref) / np.abs(A_ref)),
        }
        r_all = max(max(rels.values()), max(rels_py.values()))
        worst = max(worst, r_all)
        check(f"(b) {name}", r_all <= 1e-13 and extra_ok,
              "C++ " + " ".join(f"{k} {v:.1e}" for k, v in rels.items())
              + " | Py " + " ".join(f"{k} {v:.1e}" for k, v in rels_py.items()))
    return worst


# ---------------------------------------------------------------------------
# (a) C++ ↔ Python
# ---------------------------------------------------------------------------
def test_a(exe, td):
    print("\n(a) C++ ↔ Python (軸対称、K・b 相対 ≤1e-12、解 ≤1e-9 K)")
    # 軸に内部節点が乗る非一様メッシュ + 温度依存 k + 位置で変わる h/T_c
    nodes, tris, E = rect_mesh(0.0, 0.03, 0.0, 0.012, 12, 9, jitter=0.4, seed=3)
    # 軸 (y=0) は断熱、界面は上面 y1、Robin は x0/x1 (半径が変わる辺)
    robin = E["x0"] + E["x1"]
    rng = np.random.default_rng(7)
    hh = 2000.0 + 3000.0 * rng.random(len(robin))
    tc = 350.0 + 100.0 * rng.random(len(robin))
    h5p = Path(td) / "a.h5"
    write_solid_h5(h5p, nodes, tris, E["y1"], robin, hh, tc,
                   [300.0, 600.0, 900.0], [20.0, 25.0, 28.0], source="axisym a")
    op = op_from_h5(h5p, axisym=True)
    N, n = op.N, op.n
    x = op.xy[:, 0]
    u = 300.0 + 500.0 * (x - x.min()) / (x.max() - x.min()) + 50.0 * np.sin(300.0 * op.xy[:, 1])
    Kc, bc = tool(exe, "matrix", h5p, u, td=td)
    Kp, bp = op.assemble_full(u)
    rK = rel_dict(Kc, lower_dict(Kp))
    rb = np.max(np.abs(bc - bp)) / np.max(np.abs(bp))
    check("(a) 組立 K・b", rK <= 1e-12 and rb <= 1e-12, f"K 相対 {rK:.2e}  b 相対 {rb:.2e}")

    Ac = tool(exe, "lumped", h5p, td=td)
    rA = np.max(np.abs(Ac - op.area)) / np.max(op.area)
    check("(a) 界面の集中量 A_i^r", rA <= 1e-12, f"相対 {rA:.2e}")

    Qf = op.area * (2.0e5 + 1.0e5 * rng.random(n))
    T_iface = op.solve(Qf, iters=60, tol=1e-13)
    u_py = solve_all(op, T_iface)
    u_cpp = tool(exe, "solve", h5p, Qf, td=td)
    d = float(np.max(np.abs(u_cpp - u_py)))
    check("(a) 求解 (同一荷重)", d <= 1e-9, f"全 {N} 節点 max|ΔT| = {d:.2e} K "
          f"(T {u_cpp.min():.1f}..{u_cpp.max():.1f} K)")
    return rK, rb, d


def test_a_real(exe, td, h5src, shift):
    """実在の固体メッシュを y 方向に平行移動し (r > 0 にして) 軸対称として突き合わせる。"""
    h5p = Path(td) / "real_shift.h5"
    shutil.copy(h5src, h5p)
    with h5py.File(h5p, "r+") as f:
        c = f["MESH/COORD"][:]
        c[:, 1] += shift
        f["MESH/COORD"][...] = c
        ic = f["IFACE/COORD"][:]
        ic[:, 1] += shift
        f["IFACE/COORD"][...] = ic
    op = op_from_h5(h5p, axisym=True)
    x = op.xy[:, 0]
    u = 300.0 + 400.0 * (x - x.min()) / max(1e-30, x.max() - x.min())
    Kc, bc = tool(exe, "matrix", h5p, u, td=td)
    Kp, bp = op.assemble_full(u)
    rK = rel_dict(Kc, lower_dict(Kp))
    rb = np.max(np.abs(bc - bp)) / np.max(np.abs(bp))
    rng = np.random.default_rng(11)
    Qf = op.area * (3.0e5 + 2.0e5 * rng.random(op.n))
    T_iface = op.solve(Qf, iters=60, tol=1e-13)
    u_py = solve_all(op, T_iface)
    u_cpp = tool(exe, "solve", h5p, Qf, td=td)
    d = float(np.max(np.abs(u_cpp - u_py)))
    check(f"(a) 実形状 {Path(h5src).name} (y+{shift} m)", rK <= 1e-12 and rb <= 1e-12 and d <= 1e-9,
          f"N={op.N}  K 相対 {rK:.2e}  b 相対 {rb:.2e}  max|ΔT| {d:.2e} K")


# ---------------------------------------------------------------------------
# (c) 厚肉円筒殻  (e) 収支
# ---------------------------------------------------------------------------
def test_c(exe, td):
    print("\n(c) 厚肉円筒殻 (収束率 ≥1.8、16 層で内面誤差 ≤0.1 % of 温度降下) / (e) 収支")
    r1, r2, W = 0.01, 0.02, 0.004
    k, h, Tc, q1 = 15.0, 3000.0, 400.0, 2.0e5
    Tex = lambda r: Tc + q1 * r1 * (1.0 / (h * r2) + np.log(r2 / r) / k)   # noqa: E731
    drop = q1 * r1 * math.log(r2 / r1) / k                                   # 殻の温度降下
    res = {"axisym": [], "planar": []}
    bal_worst = 0.0
    # 格子は**全方向を同率で細かくする** (軸方向 nx = 3/4·nl、セルの縦横比を固定)。
    # 軸方向を固定して半径方向だけ細かくすると、双対性の評価 (L2 ≲ h_max·エネルギー誤差) で
    # h_max = Δx が縮まないので 2 次にならない (下の参考表。対角の向きが揃ったメッシュでは 1 次)。
    for nl in (4, 8, 16):
        nodes, tris, E = rect_mesh(0.0, W, r1, r2, 3 * nl // 4, nl)
        h5p = Path(td) / f"shell_{nl}.h5"
        write_solid_h5(h5p, nodes, tris, E["y0"], E["y1"], [h] * len(E["y1"]), [Tc] * len(E["y1"]),
                       [300.0], [k], source=f"shell nl={nl}")
        op = op_from_h5(h5p, axisym=True)
        Ar = tool(exe, "lumped", h5p, td=td)
        Qf = q1 * Ar                                          # 界面荷重 Q_i = q1 A_i^r [W/rad]
        for geom in ("axisym", "planar"):
            u = tool(exe, "solve", h5p, Qf, axisym=(geom == "axisym"), td=td)
            r = op.xy[:, 1]
            e_in = float(np.max(np.abs(u[op.iface] - Tex(r1))))
            e_all = float(np.max(np.abs(u - Tex(r))))
            res[geom].append((e_in, e_all))
            if geom == "axisym":
                qh = tool(exe, "field", h5p, u, td=td)
                bal = abs(qh.sum() - Qf.sum()) / abs(Qf.sum())
                bal_worst = max(bal_worst, bal)
                # Python オラクルでも同じ問題を解く
                T_i = op.solve(Qf, iters=30, tol=1e-13)
                u_py = solve_all(op, T_i)
                dpy = float(np.max(np.abs(u_py - u)))
                print(f"    nl={nl:2d}  内面誤差 {e_in:.3e} K  全節点 {e_all:.3e} K  "
                      f"収支 {bal:.1e}  |C++−Py| {dpy:.1e} K")
    e_in = [e[0] for e in res["axisym"]]
    e_all = [e[1] for e in res["axisym"]]
    rates = [math.log2(e_all[i] / e_all[i + 1]) for i in range(2)]
    rates_in = [math.log2(e_in[i] / e_in[i + 1]) for i in range(2)]
    check("(c) 収束率 (全節点 max 誤差)", min(rates) >= 1.8,
          "誤差 " + " ".join(f"{e:.3e}" for e in e_all) + " K ; rate " + " ".join(f"{x:.2f}" for x in rates)
          + " (内面 rate " + " ".join(f"{x:.2f}" for x in rates_in) + ")")
    check("(c) 16 層の内面誤差 ≤0.1 % of 温度降下", e_in[-1] <= 1e-3 * drop,
          f"{e_in[-1]:.3e} K ≤ {1e-3 * drop:.3e} K (温度降下 {drop:.3f} K, "
          f"{100 * e_in[-1] / drop:.4f} %)")
    p_in = res["planar"][-1][0]
    check("(c) 検出力: 平面のまま組むと FAIL する", p_in > 1e-3 * drop,
          f"平面で組んだ 16 層の内面誤差 {p_in:.3e} K = {100 * p_in / drop:.1f} % of 温度降下")
    check("(e) 収支 (c) 界面入熱 = Robin 持ち去り", bal_worst <= 1e-12, f"相対 最大 {bal_worst:.2e}")

    # 参考 (判定しない): 軸方向の分割数を固定して半径方向だけ細かくした場合
    print("    参考 (判定外): 軸方向 3 分割固定・半径方向 4/8/16/32 層の全節点 max 誤差")
    for diag in ("alt", "uniform"):
        es = []
        for nl in (4, 8, 16, 32):
            nodes, tris, E = rect_mesh(0.0, W, r1, r2, 3, nl, diag=diag)
            h5p = Path(td) / f"shellfix_{nl}.h5"
            write_solid_h5(h5p, nodes, tris, E["y0"], E["y1"], [h] * len(E["y1"]),
                           [Tc] * len(E["y1"]), [300.0], [k], source=f"shell fixed-x nl={nl}")
            Ar = tool(exe, "lumped", h5p, td=td)
            u = tool(exe, "solve", h5p, q1 * Ar, td=td)
            with h5py.File(h5p, "r") as f:
                r = np.asarray(f["MESH/COORD"][:, 1], float)
            es.append(float(np.max(np.abs(u - Tex(r)))))
        rr = [math.log2(es[i] / es[i + 1]) for i in range(len(es) - 1)]
        print(f"      対角 {diag:7s}: 誤差 " + " ".join(f"{e:.3e}" for e in es)
              + " K ; rate " + " ".join(f"{x:.2f}" for x in rr))


# ---------------------------------------------------------------------------
# (d) 軸対称円板  (e) 収支
# ---------------------------------------------------------------------------
def test_d(exe, td):
    print("\n(d) 軸対称円板 (全節点誤差 ≤1e-9 K、q_hole = consistent 厳密値 相対 ≤1e-12) / (e) 収支")
    t, r1 = 0.003, 0.005
    r2 = 4.0 * r1
    k, h, Tc, q = 12.0, 2500.0, 350.0, 3.0e5
    Tex = lambda x: Tc + q / h + q * (t - x) / k           # noqa: E731
    # 内部節点をずらした非一様メッシュ (境界は動かさない)
    nodes, tris, E = rect_mesh(0.0, t, r1, r2, 5, 13, jitter=0.5, seed=5)
    h5p = Path(td) / "disc.h5"
    write_solid_h5(h5p, nodes, tris, E["x0"], E["x1"], [h] * len(E["x1"]), [Tc] * len(E["x1"]),
                   [300.0], [k], source="disc")
    op = op_from_h5(h5p, axisym=True)
    Ar = tool(exe, "lumped", h5p, td=td)
    Qf = q * Ar
    u = tool(exe, "solve", h5p, Qf, td=td)
    err = float(np.max(np.abs(u - Tex(op.xy[:, 0]))))
    T_i = op.solve(Qf, iters=30, tol=1e-13)
    err_py = float(np.max(np.abs(solve_all(op, T_i) - Tex(op.xy[:, 0]))))
    check("(d) 全節点温度 (C++ / Python)", err <= 1e-9 and err_py <= 1e-9,
          f"max|T−T_ex| C++ {err:.2e} K / Py {err_py:.2e} K  (N={op.N})")

    # q_hole の厳密値: T−Tc = q/h (Robin 面で一定) → q_hole_i = q ∫ N_i r ds (Robin 辺の consistent 積分)
    with h5py.File(h5p, "r") as f:
        re_ = np.asarray(f["ROBIN/EDGES"][:], int)
    y = op.xy[:, 1]
    q_ex = np.zeros(op.N)
    for a, b in re_:
        L = abs(y[b] - y[a])                                # Robin 面は x=t (r 方向の辺)
        q_ex[a] += q * L / 6.0 * (2.0 * y[a] + y[b])
        q_ex[b] += q * L / 6.0 * (y[a] + 2.0 * y[b])
    qh = tool(exe, "field", h5p, u, td=td)
    rq = float(np.max(np.abs(qh - q_ex)) / np.max(np.abs(q_ex)))
    check("(d) q_hole 節点値 = consistent 積分の厳密値", rq <= 1e-12, f"相対 {rq:.2e}")
    bal = abs(qh.sum() - Qf.sum()) / abs(Qf.sum())
    check("(e) 収支 (d) 界面入熱 = Robin 持ち去り", bal <= 1e-12,
          f"{Qf.sum():.10e} vs {qh.sum():.10e} W/rad (相対 {bal:.2e})")


# ---------------------------------------------------------------------------
# 平面のビット同一 (変更前のツールと比べる)
# ---------------------------------------------------------------------------
def test_planar_bits(exe, base, h5s, td):
    print("\n(追加) 平面のビット同一 (変更前ツールとの matvec / solve 出力のバイト比較)")
    for h5 in h5s:
        with h5py.File(h5, "r") as f:
            xy = np.asarray(f["MESH/COORD"][:], float)
            nI = len(f["IFACE/NODES"])
        N = len(xy)
        x = xy[:, 0]
        rng = np.random.default_rng(2)
        u = 300.0 + 400.0 * (x - x.min()) / max(1e-30, x.max() - x.min()) + rng.random(N)
        Qf = 1.0 + rng.random(nI)
        ok, det = True, []
        for mode, inp in (("matvec", u), ("solve", Qf)):
            d1, d2 = Path(td) / "new", Path(td) / "old"
            d1.mkdir(exist_ok=True)
            d2.mkdir(exist_ok=True)
            o1 = tool(exe, mode, h5, inp, axisym=False, td=d1, raw=True).read_bytes()
            o2 = tool(base, mode, h5, inp, axisym=False, td=d2, raw=True).read_bytes()
            same = (o1 == o2)
            ok &= same
            det.append(f"{mode} {'一致' if same else '不一致'} ({len(o1)} B)")
        check(f"平面ビット同一 {h5}", ok, ", ".join(det))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--exe", default="solver_density_cuda/.build-native/release/solid_fem2d_tool")
    ap.add_argument("--base-exe", default=None, help="変更前の solid_fem2d_tool (平面のビット同一)")
    ap.add_argument("--planar-h5", nargs="*",
                    default=["case/58.conjugate_slot/mesh/solid_front_tc300_nl16.h5",
                             "case/53.c3x_vane_cht/mesh/solid_c3x.h5"],
                    help="平面のビット同一と実形状の (a) に使う固体 h5")
    a = ap.parse_args()
    if not Path(a.exe).exists():
        sys.exit(f"[test_solid_fem2d_axisym] {a.exe} が無い (make solid_fem2d_tool)")
    h5s = [p for p in a.planar_h5 if Path(p).exists()]

    with tempfile.TemporaryDirectory() as td:
        test_b(a.exe, td)
        test_a(a.exe, td)
        for h5 in h5s:
            with h5py.File(h5, "r") as f:
                ymin, ymax = float(f["MESH/COORD"][:, 1].min()), float(f["MESH/COORD"][:, 1].max())
            test_a_real(a.exe, td, h5, shift=(ymax - ymin) * 0.5 - ymin)
        test_c(a.exe, td)
        test_d(a.exe, td)
        if a.base_exe:
            test_planar_bits(a.exe, a.base_exe, h5s, td)
        else:
            print("\n(追加) 平面のビット同一: --base-exe 未指定のため省略")

    print(f"\nVERDICT: {'PASS (all)' if not FAILS else 'FAIL: ' + ', '.join(FAILS)}")
    sys.exit(1 if FAILS else 0)


if __name__ == "__main__":
    main()
