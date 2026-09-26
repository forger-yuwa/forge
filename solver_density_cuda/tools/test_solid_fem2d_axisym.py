#!/usr/bin/env python3
r"""軸対称の固体 FE (`fem2d`, $r=y$ の重み) を検証する。

plan boundary-cht-axisymmetric-fem2d §6 V-ax1 (事前登録の合格ライン)。C++ (`conjugate/solidFem2d.cpp`、
`solid_fem2d_tool --axisym`) と参照オラクル (`tools/solid_fem2d.py` の `Fem2DOperator(axisym=True)`) の両方を見る。

  (a) C++ ↔ Python: 軸対称で $K$・$b$ の相対差 ≤1e-12、同一荷重での解 ≤1e-9 K
  (b) 独立な辺積分: 半径が変わる辺 3 通り (半径方向・斜め・軸の近く) で Robin 行列・荷重・$A_i^r$・`q_hole` を
      **この評価器の中で独立に取った Gauss–Legendre 積分** ($\int h N_iN_j r\,ds$ 等、5 点 = 9 次まで厳密) と照合、
      相対 ≤1e-13。C++↔Python の一致だけでは共通の誤りを排除できないため
  旧 (c) 厚肉円筒殻 (旧パラメータ r1=0.01, r2=0.02, W=0.004, k=15, h=3000, Tc=400, q1=2e5):
      **旧登録の判定として出す (記録のみ、終了コードに使わない)**。最初の格子 (軸 3 固定・対角交互・
      半径 4/8/16 層) の次数 ≥1.8 と、改訂登録 (軸 3N/4・同方向対角) の 8→16・16→32 の次数 ≥1.8。
      縦横比固定・対角交互の結果は「事後選定の格子 (判定に使わない)」として表示のみ。
  (c‴) 置き換え試験 (新パラメータ r1=0.02, r2=0.06, W=0.01, k=40, h=500, Tc=350, q1=1e5、
      N=8/16/32/64、軸 nx=round(N W/(r2−r1))、C++ と Python の差 ≤1e-9 K を全対象で合否に含める):
      実行前の判別 A/B (旧パラメータ・N=8・軸 6・周期、A 同方向 = 0 / B 交互 = ±C a R d/[3(R²−d²/4)]、
      |差|/(q1 r1 Δx) ≤1e-11) が成立したときだけ走らせる。
      (c1) 周期端面・同方向対角: 2D = 独立 1D FE 解 ≤1e-8 K、残差 ‖Ku−b‖∞/‖b‖∞ ≤1e-11
           (C++ tool は周期端面を持たないので `matrix` の組立て行列を Python 側で周期自由度に集約して解く)
      (c2) 自然端面・交互対角: 次数 16→32・32→64 ≥1.8、N=64 で ≤0.05 % of 殻の温度降下
      (c3) 自然端面・同方向対角: N=64 で ≤0.1 %、誤差が N とともに単調減少 (次数は判定しない)
      (c4) 検出力: 平面のまま組むと (c2) の N=64 が FAIL
  (d) 軸対称円板: $x\in[0,t]$, $r\in[r_1,4r_1]$、界面 $x=0$ に一様熱流束、Robin $x=t$。解は $x$ の 1 次関数
      (線形要素で厳密) → 全節点温度誤差 ≤1e-9 K、`q_hole` 節点値 = consistent 積分の厳密値 (相対 ≤1e-12)
  (e) 収支: 界面入熱 = Robin 持ち去り (`q_hole` の総和) 相対 ≤1e-12、旧 (c)・(c‴)・(d)

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
import scipy.sparse as sp
import scipy.sparse.linalg as spla

sys.path.insert(0, str(Path(__file__).resolve().parent))
from solid_fem2d import Fem2DOperator  # noqa: E402
from solid_mesh_to_h5 import write_solid_h5  # noqa: E402
from test_solid_mesh_to_h5 import solve_all  # noqa: E402

FAILS: list[str] = []
OLD_FAILS: list[str] = []      # 旧 V-ax1(c) 登録の判定 (記録のみ。終了コードに使わない)
ENV = dict(os.environ)
ENV.setdefault("LD_LIBRARY_PATH", "/usr/lib/x86_64-linux-gnu/hdf5/serial")


def check(name: str, ok: bool, detail: str = "", bucket=None):
    print(f"  {'PASS' if ok else 'FAIL'}  {name}" + (f"   {detail}" if detail else ""))
    if not ok:
        (FAILS if bucket is None else bucket).append(name)


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
# 円筒殻の共通部品
# ---------------------------------------------------------------------------
def shell_h5(td, tag, r1, r2, W, nx, nr, h, Tc, k, diag):
    """円筒殻 (x∈[0,W], r∈[r1,r2]) の固体 h5 を書く。界面 = 内面 r1、Robin = 外面 r2、端面 x=0,W は自然 (断熱)。"""
    nodes, tris, E = rect_mesh(0.0, W, r1, r2, nx, nr, diag=diag)
    h5p = Path(td) / f"shell_{tag}.h5"
    write_solid_h5(h5p, nodes, tris, E["y0"], E["y1"], [h] * len(E["y1"]), [Tc] * len(E["y1"]),
                   [300.0], [k], source=f"shell {tag}")
    return h5p


def full_from_lower(Kl, n):
    """`matrix` 出力の下三角 dict を対称な疎行列にする。"""
    rows, cols, vals = [], [], []
    for (i, j), v in Kl.items():
        rows.append(i); cols.append(j); vals.append(v)
        if i != j:
            rows.append(j); cols.append(i); vals.append(v)
    return sp.csr_matrix((vals, (rows, cols)), shape=(n, n))


def iface_load(op, Ar, q1):
    """界面荷重 Q_i = q1 A_i^r を全節点ベクトルに散らす。"""
    Q = np.zeros(op.N)
    Q[op.iface] += q1 * np.asarray(Ar, float)
    return Q


def periodic_map(xy, W):
    """x=W の節点を同じ r の x=0 の節点に同一視する集約行列 P (N×N_d) を返す。"""
    tol = 1e-9 * W
    left = {round(float(y), 12): i for i, (x, y) in enumerate(xy) if abs(x) <= tol}
    dof = -np.ones(len(xy), int)
    nd = 0
    for i, (x, y) in enumerate(xy):
        if abs(x - W) > tol:
            dof[i] = nd
            nd += 1
    for i, (x, y) in enumerate(xy):
        if abs(x - W) <= tol:
            dof[i] = dof[left[round(float(y), 12)]]
    P = sp.csr_matrix((np.ones(len(xy)), (np.arange(len(xy)), dof)), shape=(len(xy), nd))
    return P


def fe1d(r, r1, r2, nr, k, h, Tc, q1):
    """独立な 1 次元 FE 解 ((c″) の式): T_N = T_c + q1 r1/(h r2)、
    T_j − T_{j+1} = q1 r1 Δr / [k (r_j + r_{j+1})/2]。節点の r で値を返す。"""
    rj = np.linspace(r1, r2, nr + 1)
    C = q1 * r1
    T = np.empty(nr + 1)
    T[nr] = Tc + C / (h * r2)
    for j in range(nr - 1, -1, -1):
        T[j] = T[j + 1] + C * (rj[j + 1] - rj[j]) / (k * 0.5 * (rj[j] + rj[j + 1]))
    jj = np.rint((np.asarray(r) - r1) / ((r2 - r1) / nr)).astype(int)
    return T[jj]


def linres(K, u, rhs):
    return float(np.max(np.abs(K @ u - rhs)) / np.max(np.abs(rhs)))


# ---------------------------------------------------------------------------
# 旧 (c) 厚肉円筒殻 — 旧登録の判定 (記録のみ。終了コードに使わない)  /  (e) 収支
# ---------------------------------------------------------------------------
def test_c_old(exe, td):
    print("\n旧 (c) 厚肉円筒殻 [旧登録の判定・記録のみ。終了コードには使わない。(e) 収支は終了コードに含める]")
    r1, r2, W = 0.01, 0.02, 0.004
    k, h, Tc, q1 = 15.0, 3000.0, 400.0, 2.0e5
    Tex = lambda r: Tc + q1 * r1 * (1.0 / (h * r2) + np.log(r2 / r) / k)   # noqa: E731
    drop = q1 * r1 * math.log(r2 / r1) / k                                   # 殻の温度降下

    def run(nx_of, nls, diag, tag, geoms=("axisym",), bal=False, py=False):
        out = {g: [] for g in geoms}
        bal_worst = 0.0
        for nl in nls:
            h5p = shell_h5(td, f"old_{tag}_{nl}", r1, r2, W, nx_of(nl), nl, h, Tc, k, diag)
            op = op_from_h5(h5p, axisym=True)
            Ar = tool(exe, "lumped", h5p, td=td)
            Qf = q1 * Ar
            for geom in geoms:
                u = tool(exe, "solve", h5p, Qf, axisym=(geom == "axisym"), td=td)
                r = op.xy[:, 1]
                out[geom].append((float(np.max(np.abs(u[op.iface] - Tex(r1)))),
                                  float(np.max(np.abs(u - Tex(r))))))
                if geom == "axisym" and bal:
                    qh = tool(exe, "field", h5p, u, td=td)
                    bal_worst = max(bal_worst, abs(qh.sum() - Qf.sum()) / abs(Qf.sum()))
                if geom == "axisym" and py:
                    T_i = op.solve(Qf, iters=30, tol=1e-13)
                    dpy = float(np.max(np.abs(solve_all(op, T_i) - u)))
                    print(f"    nl={nl:2d}  内面誤差 {out[geom][-1][0]:.3e} K  全節点 {out[geom][-1][1]:.3e} K"
                          f"  |C++−Py| {dpy:.1e} K")
        return out, bal_worst

    def rates(es):
        return [math.log2(es[i] / es[i + 1]) for i in range(len(es) - 1)]

    # 最初の格子 (初稿の判定): 軸方向 3 分割固定・対角交互・半径方向 4/8/16 層
    print("  最初の格子 (軸 3 固定・対角交互・半径 4/8/16 層) — 旧登録の判定")
    res, bal_worst = run(lambda nl: 3, (4, 8, 16), "alt", "first", geoms=("axisym", "planar"),
                         bal=True, py=True)
    e_all = [e[1] for e in res["axisym"]]
    e_in = [e[0] for e in res["axisym"]]
    rr = rates(e_all)
    check("旧(c) 最初の格子: 収束率 ≥1.8 (全節点 max 誤差)", min(rr) >= 1.8,
          "誤差 " + " ".join(f"{e:.3e}" for e in e_all) + " K ; rate " + " ".join(f"{x:.2f}" for x in rr),
          bucket=OLD_FAILS)
    check("旧(c) 最初の格子: 16 層の内面誤差 ≤0.1 % of 温度降下", e_in[-1] <= 1e-3 * drop,
          f"{e_in[-1]:.3e} K ({100 * e_in[-1] / drop:.4f} %, 温度降下 {drop:.3f} K)", bucket=OLD_FAILS)
    p_in = res["planar"][-1][0]
    check("旧(c) 最初の格子: 検出力 (平面のまま組むと FAIL)", p_in > 1e-3 * drop,
          f"平面で組んだ 16 層の内面誤差 {p_in:.3e} K = {100 * p_in / drop:.1f} %", bucket=OLD_FAILS)
    check("(e) 収支 旧(c) 最初の格子 界面入熱 = Robin 持ち去り", bal_worst <= 1e-12,
          f"相対 最大 {bal_worst:.2e}")

    # 改訂登録 (縦横比固定 軸 3N/4・全節点 max 誤差・両対角): 同方向対角は (c′) B で FAIL
    print("  改訂登録 (軸 3N/4・同方向対角・半径 4/8/16/32 層) — 8→16・16→32 の次数で判定")
    res, _ = run(lambda nl: 3 * nl // 4, (4, 8, 16, 32), "uniform", "rev_uni")
    e_all = [e[1] for e in res["axisym"]]
    rr = rates(e_all)
    check("旧(c) 改訂登録 同方向対角: 8→16・16→32 の次数 ≥1.8", min(rr[1:]) >= 1.8,
          "誤差 " + " ".join(f"{e:.3e}" for e in e_all) + " K ; rate " + " ".join(f"{x:.3f}" for x in rr),
          bucket=OLD_FAILS)

    # 事後選定の格子 (判定に使わない): 最初の格子の FAIL を見た後に選んだ 軸 3N/4・対角交互
    res, _ = run(lambda nl: 3 * nl // 4, (4, 8, 16), "alt", "posthoc")
    e_all = [e[1] for e in res["axisym"]]
    e_in = [e[0] for e in res["axisym"]]
    print("    事後選定の格子 (判定に使わない; 軸 3N/4・対角交互・4/8/16 層): 誤差 "
          + " ".join(f"{e:.3e}" for e in e_all) + " K ; rate " + " ".join(f"{x:.2f}" for x in rates(e_all))
          + f" ; 16 層の内面 {100 * e_in[-1] / drop:.4f} %")


# ---------------------------------------------------------------------------
# (c‴) 置き換え試験の実行前判別 A/B (旧パラメータ・N=8・軸 6 分割・周期)
# ---------------------------------------------------------------------------
def test_c3_precheck(exe, td):
    """周期自由度に集約した C++ 行列に独立 1D FE 場を代入し、内部節点残差を予測と比べる。

    予測 (面の向きが揃った節点ごとに要素寄与を足した閉形式): A (同方向) = 0、
    B (交互) = ±C a R d / [3(R² − d²/4)] (C=q1 r1, a=軸刻み, d=半径刻み, R=節点半径)。
    符号は残差 = K u − b の定義で、4 本の対角が集まる節点 (rect_mesh の (i+j) 偶数) で +、そうでない節点で −。
    判定: |残差 − 予測| / (q1 r1 a) ≤1e-11。戻り値 True なら (c‴) へ進んでよい。"""
    print("\n(c‴) 実行前の判別 A/B [旧パラメータ・N=8・軸 6 分割・周期端面・C++ 行列 (matrix) を Python で周期集約]")
    r1, r2, W = 0.01, 0.02, 0.004
    k, h, Tc, q1 = 15.0, 3000.0, 400.0, 2.0e5
    nr, nx = 8, 6
    a, d, C = W / nx, (r2 - r1) / nr, q1 * r1
    ok_all = True
    for lab, diag in (("A 同方向", "uniform"), ("B 交互", "alt")):
        h5p = shell_h5(td, f"pre_{diag}", r1, r2, W, nx, nr, h, Tc, k, diag)
        op = op_from_h5(h5p, axisym=True)
        Kl, _b = tool(exe, "matrix", h5p, np.full(op.N, 300.0), td=td)
        P = periodic_map(op.xy, W)
        Kd = (P.T @ full_from_lower(Kl, op.N) @ P).tocsr()
        # 集約後の各自由度の代表節点 (x<W の節点)
        rep = np.array([np.flatnonzero((P[:, j].toarray().ravel() > 0) & (op.xy[:, 0] < W - 1e-9 * W))[0]
                        for j in range(P.shape[1])])
        xr, rr = op.xy[rep, 0], op.xy[rep, 1]
        ud = fe1d(rr, r1, r2, nr, k, h, Tc, q1)
        res = Kd @ ud                                            # 内部節点は b = Q = 0
        p_i = np.rint(xr / a).astype(int)
        q_i = np.rint((rr - r1) / d).astype(int)
        inner = (q_i > 0) & (q_i < nr)
        if diag == "uniform":
            pred = np.zeros(len(rep))
        else:
            sgn = np.where((p_i + q_i) % 2 == 0, 1.0, -1.0)
            pred = sgn * C * a * rr * d / (3.0 * (rr ** 2 - d ** 2 / 4.0))
        dev = float(np.max(np.abs(res[inner] - pred[inner])) / (C * a))
        ok = dev <= 1e-11
        ok_all &= ok
        check(f"(c‴) 判別 {lab}: |残差 − 予測| / (q1 r1 Δx) ≤1e-11", ok,
              f"内部 {int(inner.sum())} 自由度  max|残差| {np.max(np.abs(res[inner])):.6e} W/rad  "
              f"max|予測| {np.max(np.abs(pred[inner])):.6e}  max|差| {dev * C * a:.3e} → 正規化 {dev:.3e}")
    return ok_all


# ---------------------------------------------------------------------------
# (c‴) 置き換え試験 (新パラメータ)  /  (e) 収支
# ---------------------------------------------------------------------------
def test_c3(exe, td):
    print("\n(c‴) 置き換え試験 [新パラメータ r1=0.02 r2=0.06 W=0.01 k=40 h=500 Tc=350 q1=1e5、"
          "N=8/16/32/64、軸 nx=round(N W/(r2−r1))]")
    r1, r2, W = 0.02, 0.06, 0.01
    k, h, Tc, q1 = 40.0, 500.0, 350.0, 1.0e5
    Tex = lambda r: Tc + q1 * r1 * (1.0 / (h * r2) + np.log(r2 / r) / k)   # noqa: E731
    drop = q1 * r1 * math.log(r2 / r1) / k
    Ns = (8, 16, 32, 64)
    nx_of = lambda N: int(round(N * W / (r2 - r1)))                          # noqa: E731
    print(f"    殻の温度降下 {drop:.6f} K、軸分割 " + " ".join(str(nx_of(N)) for N in Ns))
    dcp_worst, bal_worst = 0.0, 0.0

    def rates(es):
        return [math.log2(es[i] / es[i + 1]) for i in range(len(es) - 1)]

    # (c1) 周期端面・同方向対角: C++ は matrix で組んだ行列を Python で周期集約して直接求解
    print("  (c1) 周期端面・同方向対角 (C++ は `matrix` の組立て行列を Python 側で周期自由度に集約して求解)")
    c1_ok = True
    for N in Ns:
        h5p = shell_h5(td, f"c1_{N}", r1, r2, W, nx_of(N), N, h, Tc, k, "uniform")
        op = op_from_h5(h5p, axisym=True)
        P = periodic_map(op.xy, W)
        u1 = fe1d(op.xy[:, 1], r1, r2, N, k, h, Tc, q1)
        out = {}
        Ar_c = tool(exe, "lumped", h5p, td=td)
        Kl, bc = tool(exe, "matrix", h5p, np.full(op.N, 300.0), td=td)
        Kp, bp = op.assemble_full(np.full(op.N, 300.0))
        for who, K, b, Ar in (("C++", full_from_lower(Kl, op.N), bc, Ar_c), ("Py", Kp, bp, op.area)):
            Kd = (P.T @ K @ P).tocsc()
            rhs = P.T @ (b + iface_load(op, Ar, q1))
            ud = spla.spsolve(Kd, rhs)
            out[who] = (P @ ud, linres(Kd, ud, rhs))
        e1 = max(float(np.max(np.abs(out[w][0] - u1))) for w in out)
        rs = max(out[w][1] for w in out)
        dcp = float(np.max(np.abs(out["C++"][0] - out["Py"][0])))
        dcp_worst = max(dcp_worst, dcp)
        ok = e1 <= 1e-8 and rs <= 1e-11 and dcp <= 1e-9
        c1_ok &= ok
        check(f"(c1) N={N:2d}: |2D − 1D FE| ≤1e-8 K・残差 ≤1e-11・|C++−Py| ≤1e-9 K", ok,
              f"|2D−1D| {e1:.3e} K  残差 {rs:.2e}  |C++−Py| {dcp:.2e} K  (C++ 残差 {out['C++'][1]:.2e} / "
              f"Py {out['Py'][1]:.2e})")

    # (c2)(c3)(c4) 自然断熱端面: C++ は solve、Python はオラクルの組立てを直接求解
    def natural(diag, geoms=("axisym",), Nlist=Ns):
        nonlocal dcp_worst, bal_worst
        errs = {g: [] for g in geoms}
        for N in Nlist:
            h5p = shell_h5(td, f"nat_{diag}_{N}", r1, r2, W, nx_of(N), N, h, Tc, k, diag)
            op_ax = op_from_h5(h5p, axisym=True)
            Ar = tool(exe, "lumped", h5p, td=td)
            Qf = q1 * Ar                                        # 平面で組むときも同じ荷重 [W/rad] を与える
            r = op_ax.xy[:, 1]
            for geom in geoms:
                ax = geom == "axisym"
                op = op_ax if ax else op_from_h5(h5p, axisym=False)
                u_c = tool(exe, "solve", h5p, Qf, axisym=ax, td=td)
                Kl, bc = tool(exe, "matrix", h5p, u_c, axisym=ax, td=td)
                rhs_c = bc + iface_load(op, Ar, q1)
                rc = linres(full_from_lower(Kl, op.N), u_c, rhs_c)
                Kp, bp = op.assemble_full(np.full(op.N, 300.0))
                rhs_p = bp + iface_load(op, Ar, q1)
                u_p = spla.spsolve(Kp.tocsc(), rhs_p)
                rp = linres(Kp, u_p, rhs_p)
                dcp = float(np.max(np.abs(u_c - u_p)))
                dcp_worst = max(dcp_worst, dcp)
                e = float(np.max(np.abs(u_c - Tex(r))))
                ep = float(np.max(np.abs(u_p - Tex(r))))
                ipos = int(np.argmax(np.abs(u_c - Tex(r))))
                errs[geom].append(e)
                msg = ""
                if ax:
                    qh = tool(exe, "field", h5p, u_c, td=td)
                    bal = abs(qh.sum() - Qf.sum()) / abs(Qf.sum())
                    bal_worst = max(bal_worst, bal)
                    msg = f"  収支 {bal:.1e}"
                print(f"    {geom:6s} 対角 {diag:7s} N={N:2d} nx={nx_of(N):2d}: 全節点 max 誤差 C++ {e:.6e} K "
                      f"({100 * e / drop:.5f} %) / Py {ep:.6e} K  位置 (x={op.xy[ipos, 0]:.4g}, r={r[ipos]:.4g})"
                      f"  残差 C++ {rc:.1e} Py {rp:.1e}  |C++−Py| {dcp:.1e} K{msg}")
        return errs

    print("  (c2) 自然断熱端面・交互対角")
    e2 = natural("alt")["axisym"]
    r2s = rates(e2)
    check("(c2) 次数 16→32 と 32→64 がともに ≥1.8 (全節点 max 誤差)", r2s[1] >= 1.8 and r2s[2] >= 1.8,
          "誤差 " + " ".join(f"{e:.4e}" for e in e2) + " K ; 次数 " + " ".join(f"{x:.3f}" for x in r2s))
    check("(c2) N=64 で ≤0.05 % of 殻の温度降下", e2[-1] <= 5e-4 * drop,
          f"{e2[-1]:.4e} K = {100 * e2[-1] / drop:.5f} % (上限 {5e-4 * drop:.4e} K)")

    print("  (c3) 自然断熱端面・同方向対角 (次数は判定しない — 旧要件からの緩和。二次精度の検証ではない)")
    e3 = natural("uniform")["axisym"]
    mono = all(e3[i + 1] < e3[i] for i in range(len(e3) - 1))
    check("(c3) N=64 で ≤0.1 % of 殻の温度降下", e3[-1] <= 1e-3 * drop,
          f"{e3[-1]:.4e} K = {100 * e3[-1] / drop:.5f} % (上限 {1e-3 * drop:.4e} K)")
    check("(c3) 誤差が N とともに単調減少", mono,
          "誤差 " + " ".join(f"{e:.4e}" for e in e3) + " K ; 次数 (参考・判定しない) "
          + " ".join(f"{x:.3f}" for x in rates(e3)))

    print("  (c4) 検出力: (c2) の N=64 を平面のまま組む")
    e4 = natural("alt", geoms=("planar",), Nlist=(64,))["planar"][-1]
    check("(c4) 平面のまま組むと (c2) の N=64 が FAIL (> 0.05 %)", e4 > 5e-4 * drop,
          f"平面 {e4:.4e} K = {100 * e4 / drop:.3f} % of 温度降下")

    check("(c‴) 全対象の |C++−Py| ≤1e-9 K (まとめ)", dcp_worst <= 1e-9, f"最大 {dcp_worst:.2e} K")
    check("(e) 収支 (c‴) 自然端面 界面入熱 = Robin 持ち去り", bal_worst <= 1e-12, f"相対 最大 {bal_worst:.2e}")


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
    ap.add_argument("--precheck-only", action="store_true",
                    help="(c‴) の実行前判別 A/B だけを回す (登録: 置き換え試験より先に実行する)")
    a = ap.parse_args()
    if not Path(a.exe).exists():
        sys.exit(f"[test_solid_fem2d_axisym] {a.exe} が無い (make solid_fem2d_tool)")
    h5s = [p for p in a.planar_h5 if Path(p).exists()]

    with tempfile.TemporaryDirectory() as td:
        if a.precheck_only:
            ok = test_c3_precheck(a.exe, td)
            print(f"\nVERDICT (判別 A/B): {'PASS' if ok else 'FAIL — (c‴) を走らせず再調査'}")
            sys.exit(0 if ok else 1)
        test_b(a.exe, td)
        test_a(a.exe, td)
        for h5 in h5s:
            with h5py.File(h5, "r") as f:
                ymin, ymax = float(f["MESH/COORD"][:, 1].min()), float(f["MESH/COORD"][:, 1].max())
            test_a_real(a.exe, td, h5, shift=(ymax - ymin) * 0.5 - ymin)
        test_c_old(a.exe, td)
        # 登録: 判別 A/B が「A = 0・B = 予測式」にならなければ置き換え試験を走らせない
        if test_c3_precheck(a.exe, td):
            test_c3(a.exe, td)
        else:
            check("(c‴) 置き換え試験", False, "判別 A/B が予測どおりでないため実行しない (再調査)")
        test_d(a.exe, td)
        if a.base_exe:
            test_planar_bits(a.exe, a.base_exe, h5s, td)
        else:
            print("\n(追加) 平面のビット同一: --base-exe 未指定のため省略")

    print(f"\nVERDICT (旧 V-ax1(c) 登録・記録のみ、終了コードに使わない): "
          f"{'PASS' if not OLD_FAILS else 'FAIL: ' + ', '.join(OLD_FAILS)}")
    print(f"VERDICT ((a)(b)(c‴)(d)(e)): {'PASS (all)' if not FAILS else 'FAIL: ' + ', '.join(FAILS)}")
    sys.exit(1 if FAILS else 0)


if __name__ == "__main__":
    main()
