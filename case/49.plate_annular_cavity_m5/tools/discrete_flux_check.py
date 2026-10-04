#!/usr/bin/env python3
r"""開口収支を **ソルバの双対面そのもの**で組み直す (plan §4.8.6.3 / 残作業 #28・#4)。

**なぜ要るか**: `cavity_eval.py` の開口流束は、節点層を $(\theta, r)$ 格子へ補間してから
求積している。壁入熱はソルバの離散量 (`qwall`) なので**土俵が違う**。実際 mixA では
深さ 10 mm 以深で 6 W 合わない (§4.8.6.3) 一方、再構成の対流流束は壁温分布にほとんど反応しない。

ここでは再構成を使わず、**`PLANES/STRUCT` の双対面リストと `PLANES/surfVect` の面積ベクトル**を
そのまま使う。深さ $z_c$ を跨ぐ面 (片側の DOF が $z<z_c$、もう片側が $z\ge z_c$) だけを集め、

$$\dot m_f = \rho_f\,(\mathbf u_f\cdot\mathbf S),\qquad
  H_f = \dot m_f\,h_{0,f},\qquad
  Q_{cond,f} = k_{eff,f}\,(T_b-T_a)\,\frac{\mathbf d\cdot\mathbf S}{|\mathbf d|^2}$$

を面ごとに積む。面値は**中心平均**と**風上** (質量流束の符号) の両方を出す
(ソルバは SLAU なので厳密には一致しないが、**幾何と面リストはソルバと同一**になる)。

定常なら「$z_c$ 以深の壁入熱 = 面を通る正味流入」が成り立つはずで、
再構成側の 6 W がここで消えるなら原因は求積・補間、残るなら物理か壁入熱側である。

usage:
  python3 tools/discrete_flux_check.py <run_dir> [--depth 1,2,5,10] [--step N]
"""
import argparse
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
CASE = HERE.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(CASE))
import h5py                       # noqa: E402
import cavity_eval as ce          # noqa: E402
import geom_common as gc          # noqa: E402

CAV_WALLS = ("cav_outer", "cyl_side", "cav_floor")


def dual_faces(mesh):
    """`PLANES/STRUCT` を舐めて、面ごとの (DOF a, DOF b) を返す。

    レコードは [nn, node..., nc, cell...] の繰り返し。node 方式では cell = 節点 index。
    面の数は `PLANES/surfArea` の長さ。
    """
    with h5py.File(mesh, "r") as f:
        st = np.array(f["PLANES/STRUCT"])
        npl = len(np.array(f["PLANES/surfArea"]))
    a = np.full(npl, -1, dtype=np.int64)
    b = np.full(npl, -1, dtype=np.int64)
    i = ip = 0
    n = len(st)
    while i < n and ip < npl:
        nn = int(st[i]); i += 1 + nn
        nc = int(st[i]); i += 1
        if nc >= 2:
            a[ip] = int(st[i]); b[ip] = int(st[i + 1])
        i += nc
        ip += 1
    return a, b


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--depth", default="1,2,5,10", help="評価深さ [mm] のリスト")
    ap.add_argument("--step", type=int, default=None)
    a = ap.parse_args()
    run = Path(a.run)
    man = gc.load_manifest(run=run)
    D = ce.run_conditions(run)
    step = a.step if a.step is not None else int(ce.snapshots(run)[-1].stem.split("_")[1])
    mesh = run / "mesh.h5"

    with h5py.File(mesh, "r") as f:
        xyz = np.array(f["MESH/COORD"]).reshape(-1, 3)
        S = np.array(f["PLANES/surfVect"]).reshape(-1, 3)
    fa, fb = dual_faces(str(mesh))
    ok = (fa >= 0) & (fb >= 0)
    fa, fb, S = fa[ok], fb[ok], S[ok]
    # surfVect の向きを a->b に揃える
    d = xyz[fb] - xyz[fa]
    sgn = np.sign(np.einsum("ij,ij->i", S, d))
    sgn[sgn == 0] = 1.0
    S = S * sgn[:, None]

    with h5py.File(run / ("res_%d.h5" % step), "r") as f:
        V = {k: np.array(f["VALUE/" + k]) for k in
             ("ro", "Ux", "Uy", "Uz", "h0", "T", "vis_lam", "vis_turb") if "VALUE/" + k in f}
    need = ("ro", "Ux", "Uy", "Uz", "h0")
    miss = [k for k in need if k not in V]
    if miss:
        print("res に %s が無い (output.level 1 以上が要る)" % ",".join(miss)); return 2
    prl = float(D.get("prandtl_lam", 0.72))
    prt = float(D.get("prandtl_turb", 0.9))
    mu = V.get("vis_lam"); mut = V.get("vis_turb")
    U = np.stack([V["Ux"], V["Uy"], V["Uz"]], axis=-1)

    wh = ce.wall_heat(run, step, man, D)
    z = xyz[:, 2]
    dn = np.linalg.norm(d, axis=1)

    print("=== 離散流束による収支: %s  step %d ===" % (run, step))
    print("  双対面 %d / 節点 %d" % (len(fa), len(xyz)))
    print("  %-7s %10s %12s %12s %12s %10s" %
          ("z [mm]", "壁入熱 [W]", "対流(中心)", "対流(風上)", "伝導 [W]", "面数"))
    for dv in [float(x) for x in a.depth.split(",") if x]:
        zc = -dv * 1e-3
        below = z < zc
        cross = below[fa] != below[fb]
        if not cross.any():
            print("  %-7.3g 跨ぐ面が無い" % dv); continue
        ia, ib, Sc, dc, dnc = fa[cross], fb[cross], S[cross], d[cross], dn[cross]
        rof = 0.5 * (V["ro"][ia] + V["ro"][ib])
        uf = 0.5 * (U[ia] + U[ib])
        mdot = rof * np.einsum("ij,ij->i", uf, Sc)          # a -> b が正
        h0c = 0.5 * (V["h0"][ia] + V["h0"][ib])
        h0u = np.where(mdot >= 0, V["h0"][ia], V["h0"][ib])  # 風上
        # 「下側」への流入を正にする: a が下なら a->b は流出
        sflow = np.where(below[ia], -1.0, 1.0)
        Hc = float(np.sum(sflow * mdot * h0c))
        Hu = float(np.sum(sflow * mdot * h0u))
        Qc = float("nan")
        if "T" in V and mu is not None:
            cpf = ce.cp_of_T(D, 0.5 * (V["T"][ia] + V["T"][ib]))
            muf = 0.5 * (mu[ia] + mu[ib])
            mutf = 0.5 * (mut[ia] + mut[ib]) if mut is not None else 0.0
            kf = cpf * (muf / prl + mutf / prt)
            gproj = np.einsum("ij,ij->i", dc, Sc) / np.maximum(dnc ** 2, 1e-30)
            qc = kf * (V["T"][ib] - V["T"][ia]) * gproj      # a -> b が正 (熱は高温から低温へ)
            Qc = float(np.sum(sflow * (-qc)))
        qs = ce.wall_Q_below(wh, CAV_WALLS, zc)
        print("  %-7.3g %10.4f %12.4f %12.4f %12.4f %10d"
              % (dv, qs, Hc, Hu, Qc, int(cross.sum())))
    return 0


if __name__ == "__main__":
    sys.exit(main())
