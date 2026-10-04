#!/usr/bin/env python3
r"""等温壁の入熱を **2 通り**で出して突き合わせる (plan §4.8.6.2 / 残作業 #41)。

**なぜ要るか**: node × 等温壁 × 強制形 (`nodeWallDirichlet: 1`, `nodeIsothermalEnergyBC` 非 1) では
壁ノードの `roe` 行が decouple されるので、`viscousFlux_d.cu` が壁ノードに足した `heatflux` は
**捨てられる**。壁ダンプの `qwall` はその捨てられた値 (壁ノードの片側勾配 $k_{eff}\nabla T\cdot S$) で、
「ソルバが実際に流体から抜いた熱」とは限らない。

そこで、弱形式 (`nodeIsothermalEnergyBC: 1`, `weakIsothermalWall_d.cu`) と**同じ式**を後処理で組む:

$$q_{in} = \frac{k_{eff}}{y_1}\,(T_I - T_w) \qquad [\mathrm{W/m^2}]$$

$T_I$ は**第一内部 DOF** の温度、$y_1$ はそこまでの壁法線距離 (`check_wall_resolution.py` と同じ
引き方)。両者が合えば `qwall` は信用してよく、系統的にずれるなら CV 収支の欠損はここに居る。

usage:
  python3 tools/wall_heat_crosscheck.py <run_dir> [--step N] [--groups cav_floor,cav_outer]
"""
import argparse
import os
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
CASE = HERE.parent
ROOT = CASE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(CASE))
sys.path.insert(0, str(ROOT / "solver_density_cuda" / "tools"))
import h5py                                        # noqa: E402
import cavity_eval as ce                           # noqa: E402
import geom_common as gc                           # noqa: E402
from check_wall_resolution import dof_neighbors    # noqa: E402

WALLS = ("cav_outer", "cyl_side", "cav_floor", "cyl_top")


def first_interior(mesh, pid, nb, coords, align_min=0.5):
    """壁 DOF -> (y1, 第一内部 DOF index, 評価できたか)。`check_wall_resolution` と同じ引き方で、
    **内部 DOF の index も返す** (温度を読むため)。"""
    with h5py.File(mesh, "r") as f:
        key = "BCONDS/%d" % pid
        if key not in f:
            return None
        icells = np.array(f[key + "/iCells"])
        iplanes = np.array(f[key + "/iPlanes"])
        svec = np.array(f["PLANES/surfVect"]).reshape(-1, 3)
    nrm = {}
    for c, p in zip(icells, iplanes):
        c = int(c)
        nrm[c] = nrm.get(c, np.zeros(3)) + svec[int(p)]
    out = {}
    for c, n in nrm.items():
        ln = float(np.linalg.norm(n))
        if ln <= 0.0:
            out[c] = (np.nan, -1, False)
            continue
        nh = n / ln
        best, bestal = None, -1.0
        for (j, _p) in nb.get(c, []):
            d = coords[j] - coords[c]
            dn = float(np.linalg.norm(d))
            if dn <= 0.0:
                continue
            al = abs(float(np.dot(d, nh))) / dn
            if al > bestal:
                bestal, best = al, (abs(float(np.dot(d, nh))), j)
        out[c] = (best[0], best[1], True) if (best is not None and bestal >= align_min) \
            else (np.nan, -1, False)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--step", type=int, default=None)
    ap.add_argument("--groups", default=None)
    ap.add_argument("--align-min", type=float, default=0.5)
    a = ap.parse_args()
    run = Path(a.run)
    man = gc.load_manifest(run=run)
    D = ce.run_conditions(run)
    PID = man["phys_id"]
    step = a.step if a.step is not None else int(ce.snapshots(run)[-1].stem.split("_")[1])
    mesh = run / "mesh.h5"
    res = run / ("res_%d.h5" % step)
    if not res.exists():
        print("res が無い: %s" % res); return 2

    with h5py.File(mesh, "r") as f:
        coords = np.array(f["MESH/COORD"]).reshape(-1, 3)
    with h5py.File(res, "r") as f:
        V = {k: np.array(f["VALUE/" + k]) for k in ("T", "vis_lam", "vis_turb") if "VALUE/" + k in f}
    if "T" not in V:
        print("res に VALUE/T が無い"); return 2
    prt = float(D.get("prandtl_turb", 0.9))
    prl = float(D.get("prandtl_lam", 0.72))
    nb = dof_neighbors(str(mesh))
    tree = cKDTree(coords)
    want = set(x.strip() for x in a.groups.split(",")) if a.groups else set(WALLS)

    print("=== 壁入熱クロスチェック: %s  step %d ===" % (run, step))
    print("  `qwall`   : 壁ダンプ (壁ノードの片側勾配。強制形では残差ごと捨てられる値)")
    print("  弱形式    : k_eff (T_I - T_w) / y1  (第一内部 DOF から。SU2 型と同じ式)")
    print()
    print("  %-11s %10s %10s %8s %10s %10s" % ("壁", "qwall [W]", "弱形式 [W]", "比", "y1 中央 [m]", "評価率"))
    tot_q = tot_w = 0.0
    for g in WALLS:
        if g not in want or g not in PID:
            continue
        w = ce.wall_dump(run, PID[g], step, name=g)
        if w is None or "qwall" not in w:
            continue
        qin = -np.asarray(w["qwall"], float)
        area = np.asarray(w["w"], float)
        Qq = float(np.sum(qin * area))
        fi = first_interior(str(mesh), PID[g], nb, coords, a.align_min)
        if fi is None:
            print("  %-11s BCONDS が無い" % g); continue
        d, idx = tree.query(np.asarray(w["xyz"], float), k=1)
        ok = d < 1e-9
        y1 = np.full(len(idx), np.nan)
        jj = np.full(len(idx), -1, dtype=np.int64)
        for n, c in enumerate(idx):
            v = fi.get(int(c))
            if v is not None and v[2]:
                y1[n], jj[n] = v[0], v[1]
        good = ok & np.isfinite(y1) & (jj >= 0) & (y1 > 0)
        if not good.any():
            print("  %-11s 第一内部 DOF が引けない" % g); continue
        Tw = V["T"][idx[good]]
        Ti = V["T"][jj[good]]
        mul = V["vis_lam"][idx[good]] if "vis_lam" in V else np.full(good.sum(), D["mu_inf"])
        mut = V["vis_turb"][idx[good]] if "vis_turb" in V else np.zeros(good.sum())
        cp = ce.cp_of_T(D, Tw)
        keff = mul * cp / prl + cp * mut / prt
        Qw = float(np.sum(keff * (Ti - Tw) / y1[good] * area[good]))
        frac = float(area[good].sum() / max(area.sum(), 1e-30))
        Qw_full = Qw / max(frac, 1e-30)            # 評価できた面積で外挿
        print("  %-11s %10.4f %10.4f %8.3f %10.3e %9.1f %%"
              % (g, Qq, Qw_full, Qw_full / Qq if abs(Qq) > 1e-9 else float("nan"),
                 float(np.median(y1[good])), 100 * frac))
        if g != "cyl_top":
            tot_q += Qq
            tot_w += Qw_full
    print()
    print("  キャビティ 3 壁 合計: qwall %.4f W  vs  弱形式 %.4f W  (比 %.3f, 差 %+.4f W)"
          % (tot_q, tot_w, tot_w / tot_q if abs(tot_q) > 1e-9 else float("nan"), tot_w - tot_q))
    return 0


if __name__ == "__main__":
    sys.exit(main())
