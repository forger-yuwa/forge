#!/usr/bin/env python3
r"""C++ の固体 FE (`conjugate/solidFem2d.cpp`) を **Python 実装を真値として**検証する。

plan boundary-conjugate-heat-transfer §6 V4b (a)。移植で問題が変わっていないことを、
**後段 (連成) では気づけない粒度**で閉じておく。合格ラインは §6 に事前登録したもの:

  a1 組立    : 同じ u に対する $K_su$ の相対差 ≤ 1e-12
  a2 求解    : 実荷重 (published run の `iface_q_eff` × 集中辺長) で全節点 max|ΔT| ≤ 1e-6 K、
               孔 Robin の持ち去り総量の相対差 ≤ 1e-9
  a4 内部残差: C++ の解で ≤ 1e-9·max|Q_f|、その解を Python の残差に通しても ≤ 1e-8·max|Q_f|

使い方:
  python3 solver_density_cuda/tools/test_solid_fem2d_cpp.py \
      --h5   case/53.c3x_vane_cht/mesh/solid_c3x.h5 \
      --wall /home/sano/work/forge/case/53.c3x_vane_cht/run_0144_cht_published/it_027/res_wall_4_4000.h5 \
      --exe  solver_density_cuda/.build-native/release/solid_fem2d_tool
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from test_solid_mesh_to_h5 import build_from_h5, solve_all  # noqa: E402
from solid_mesh_to_h5 import write_solid_h5  # noqa: E402

FAILS: list[str] = []


def check(name: str, ok: bool, detail: str = ""):
    print(f"  {'PASS' if ok else 'FAIL'}  {name}" + (f"   {detail}" if detail else ""))
    if not ok:
        FAILS.append(name)


def run_tool(exe: str, mode: str, h5: str, inp: np.ndarray | None, n_out: int) -> np.ndarray:
    env = dict(os.environ)
    env.setdefault("LD_LIBRARY_PATH", "/usr/lib/x86_64-linux-gnu/hdf5/serial")
    with tempfile.TemporaryDirectory() as td:
        fi, fo = Path(td) / "in.txt", Path(td) / "out.txt"
        np.savetxt(fi, np.asarray(inp, float), fmt="%.17e")
        r = subprocess.run([exe, mode, h5, str(fi), str(fo)], env=env,
                           capture_output=True, text=True)
        if r.returncode != 0:
            sys.exit(f"[test_solid_fem2d_cpp] {mode} が落ちた:\n{r.stdout}\n{r.stderr}")
        print(f"    (C++ {mode}) {r.stdout.strip()}")
        out = np.loadtxt(fo)
    if out.size != n_out:
        sys.exit(f"[test_solid_fem2d_cpp] {mode} の出力長 {out.size} != {n_out}")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--h5", required=True, help="solid_mesh_to_h5.py が作った固体 HDF5")
    ap.add_argument("--wall", default=None,
                    help="a2 の荷重に使う壁ダンプ (`iface_q_eff` を持つ res_wall_*.h5)。"
                         "省略すると合成荷重で a2 を回す")
    ap.add_argument("--exe", default="solver_density_cuda/.build-native/release/solid_fem2d_tool")
    a = ap.parse_args()

    f = h5py.File(a.h5, "r")
    op = build_from_h5(f)
    N, n = op.N, op.n
    print(f"[test_solid_fem2d_cpp] {a.h5}  N={N} iface={n}")

    # ---- a1 組立 ----
    x = op.xy[:, 0]
    u = 300.0 + 400.0 * (x - x.min()) / max(1e-30, (x.max() - x.min()))
    K, _ = op.assemble_full(u)
    ku_py = np.asarray(K.dot(u)).ravel()
    ku_cpp = run_tool(a.exe, "matvec", a.h5, u, N)
    rel = np.max(np.abs(ku_cpp - ku_py)) / max(1e-30, np.max(np.abs(ku_py)))
    check("a1 組立 K_s u", rel <= 1e-12, f"相対差 {rel:.3e}")

    # ---- a2 求解 ----
    if a.wall:
        w = h5py.File(a.wall, "r")
        if "iface_q_eff" not in w["VALUE"]:
            sys.exit(f"[test_solid_fem2d_cpp] {a.wall} に iface_q_eff が無い "
                     "(interfaceDiag: 1 で回した run の壁ダンプを渡すこと)")
        q = np.asarray(w["VALUE"]["iface_q_eff"][:], float)
        wc = np.asarray(w["MESH/COORD"][:], float).reshape(-1, 3)[:, :2]
        # 壁ダンプ順 -> 固体界面順 (座標一致。cht_loop.build_fem2d と同じ規約)
        idx = np.empty(n, int)
        for i, p in enumerate(op.coords[:, :2]):
            k = int(np.argmin(np.hypot(*(wc - p).T)))
            if np.hypot(*(wc[k] - p)) > 1e-7:
                sys.exit(f"[test_solid_fem2d_cpp] 界面節点 {i} に一致する壁節点が無い")
            idx[i] = k
        Qf = q[idx] * op.area                      # [W/m2] * 集中辺長 [m] = [W/m]
        src = Path(a.wall).name
    else:
        rng = np.random.default_rng(1)
        Qf = op.area * (3.0e5 + 2.0e5 * rng.random(n))
        src = "合成荷重"

    T_py_iface = op.solve(Qf, iters=50, tol=1e-12)
    u_py = solve_all(op, T_py_iface)
    u_cpp = run_tool(a.exe, "solve", a.h5, Qf, N)
    dmax = float(np.max(np.abs(u_cpp - u_py)))
    check(f"a2 求解 (荷重: {src})", dmax <= 1e-6, f"全 {N} 節点 max|ΔT| = {dmax:.3e} K")

    # 孔 Robin の持ち去り総量
    def hole_total(uu):
        s = 0.0
        for (n0, n1, h, Tc) in op.robin_edges:
            L = float(np.linalg.norm(op.xy[n1] - op.xy[n0]))
            s += h * L * (0.5 * (uu[n0] + uu[n1]) - Tc)
        return s
    h_py, h_cpp = hole_total(u_py), hole_total(u_cpp)
    rel_h = abs(h_cpp - h_py) / max(1e-30, abs(h_py))
    check("a2b 孔の持ち去り総量", rel_h <= 1e-9,
          f"{h_cpp:.6f} vs {h_py:.6f} W/m  (相対 {rel_h:.3e})")

    # ---- a4 内部残差 ----
    Kc, bc = op.assemble_full(u_cpp)
    r = np.asarray(Kc.dot(u_cpp)).ravel() - np.asarray(bc).ravel()
    r_int = float(np.max(np.abs(r[op._other])))
    qmax = float(np.max(np.abs(Qf)))
    check("a4 内部残差 (C++ の解を Python の残差に通す)", r_int <= 1e-8 * qmax,
          f"{r_int:.3e} <= {1e-8*qmax:.3e} W/m")

    # 収支: 界面が受け取った熱 = 孔が持ち去った熱
    q_if = float(np.sum(r[op.iface]))
    bal = abs(q_if - h_cpp) / max(1e-30, abs(q_if))
    check("a4b 収支 (界面入熱 = 孔の持ち去り)", bal <= 1e-9,
          f"{q_if:.6f} vs {h_cpp:.6f} W/m  (相対 {bal:.3e})")

    # ---- a3 円環解析解 (格子収束。Python の T1 を C++ で回す) ----
    import math
    ra, rb, kk, hh, Tc, qq = 0.002, 0.010, 20.0, 5000.0, 400.0, 5.0e4
    Qtot = qq * 2 * math.pi * rb
    Ta_ex = Tc + Qtot / (2 * math.pi * ra * hh)
    Tb_ex = Ta_ex + Qtot / (2 * math.pi * kk) * math.log(rb / ra)
    errs, last = [], (0.0, 0.0)
    with tempfile.TemporaryDirectory() as td:
        for (nr, nt) in ((4, 24), (8, 48), (16, 96)):
            rs = np.linspace(ra, rb, nr + 1)
            th = np.linspace(0.0, 2.0 * math.pi, nt, endpoint=False)
            nodes = np.array([[r * math.cos(t), r * math.sin(t)] for r in rs for t in th])
            ix = lambda i, j: i * nt + (j % nt)          # noqa: E731
            tris = [[ix(i, j), ix(i + 1, j), ix(i + 1, j + 1)] for i in range(nr) for j in range(nt)] \
                 + [[ix(i, j), ix(i + 1, j + 1), ix(i, j + 1)] for i in range(nr) for j in range(nt)]
            oe = [(ix(nr, j), ix(nr, j + 1)) for j in range(nt)]
            ie = [(ix(0, j), ix(0, j + 1)) for j in range(nt)]
            h5p = Path(td) / f"ann_{nr}.h5"
            write_solid_h5(h5p, nodes, tris, oe, ie, [hh] * len(ie), [Tc] * len(ie),
                           [300.0], [kk], source="annulus")
            fa = h5py.File(h5p, "r")
            opa = build_from_h5(fa)
            ua = run_tool(a.exe, "solve", str(h5p), qq * opa.area, opa.N)
            rr = np.hypot(opa.xy[:, 0], opa.xy[:, 1])
            Tb_num = float(np.mean(ua[opa.iface]))
            Ta_num = float(np.mean(ua[rr < ra * 1.0001]))
            errs.append(max(abs(Tb_num - Tb_ex), abs(Ta_num - Ta_ex)))
            last = (Tb_num, Ta_num)
    rates = [math.log2(errs[i] / errs[i + 1]) for i in range(len(errs) - 1)]
    check("a3 円環解析解 (格子収束)", errs[-1] < 0.02 * (Tb_ex - Tc) and min(rates) > 1.6,
          f"T(b)={last[0]:.3f} (厳密 {Tb_ex:.3f}) T(a)={last[1]:.3f} (厳密 {Ta_ex:.3f}) ; "
          + "誤差 " + " ".join(f"{e:.3e}" for e in errs)
          + " ; rate " + " ".join(f"{r:.2f}" for r in rates))

    print(f"\nVERDICT: {'PASS (all)' if not FAILS else 'FAIL: ' + ', '.join(FAILS)}")
    sys.exit(1 if FAILS else 0)


if __name__ == "__main__":
    main()
