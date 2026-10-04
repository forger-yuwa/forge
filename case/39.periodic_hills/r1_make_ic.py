#!/usr/bin/env python3
"""R1 (plan boundary-node-periodic-gradient-fix §6 R1) の初期場と stage_manifest を作る。

旧 (1266aba1) / 新 (565959c7) のバイナリで**同一メッシュ・同一 IC・同一段構成**にするため、
IC はここで決定論的に作る (乱数なし)。段構成は README の run_0001 → run_0002 の実績に倣う:

  S1_spin : 静止 IC (u=0, P=101325 Pa, T=300 K, k=kInit, ω=omegaInit) から固定体積力でスピンアップ
  S2_dev  : S1 最終場の ρ・P・k・ω を保ったまま、速度だけ**成形 IC** (列ごとに質量流束一定の放物線,
            丘頂断面のバルク 52 m/s = run_0002 の u_max 78 と同じ) に張り替えて発達させる

使い方:
    python3 r1_make_ic.py static  MESH_H5 OUT_H5
    python3 r1_make_ic.py shaped  SRC_RES_H5 MESH_H5 OUT_H5 [--ub 52]
    python3 r1_make_ic.py manifest RUN_DIR     (solverConfig_S1_spin.yaml / _S2_dev.yaml から)

保存量の index コピーは `solver_density_cuda/tools/restart_field.py` と同じ考え方 (原始量から組み直さない)。
ただし mesh h5 に roK/roOmega が無いと restart_field は写さないので、ここでは作成して書く。
"""
import argparse
import os
import shutil
import sys

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "mesh"))
from make_hill_mesh import H, LX, LY_TOP, wall_y  # noqa: E402

# run_0012 と同じ物性 (CPG 空気)
CP, GAM = 1004.5, 1.4
R = CP * (GAM - 1.0) / GAM
P0, T0 = 101325.0, 300.0
K0, W0 = 1.0, 1000.0  # solverConfig の kInit / omegaInit と同値 (recommended-settings §1.2)


def node_xyz(f, n):
    """DOF (node) の座標。node 変換の h5 では CELLS/centCoords が DOF 数ぶんある。"""
    for key in ("CELLS/centCoords", "MESH/COORD"):
        if key in f:
            c = f[key][:].astype(np.float64).reshape(-1, 3)
            if len(c) == n:
                return c
    sys.exit("DOF 数 %d に合う座標データセットが無い" % n)


def put(v, name, arr):
    arr = arr.astype(v["ro"].dtype)
    if name in v:
        v[name][:] = arr
    else:
        v.create_dataset(name, data=arr)


def cmd_static(mesh, out):
    shutil.copyfile(mesh, out)
    with h5py.File(out, "r+") as f:
        v = f["VALUE"]
        n = v["ro"].shape[0]
        ro = np.full(n, P0 / (R * T0))
        put(v, "ro", ro)
        for c in ("roUx", "roUy", "roUz"):
            put(v, c, np.zeros(n))
        put(v, "roe", np.full(n, P0 / (GAM - 1.0)))
        put(v, "roK", ro * K0)
        put(v, "roOmega", ro * W0)
    print("static IC -> %s (ro=%.6f, k=%g, omega=%g)" % (out, P0 / (R * T0), K0, W0))


def cmd_shaped(src, mesh, out, ub):
    shutil.copyfile(mesh, out)
    with h5py.File(src, "r") as s, h5py.File(out, "r+") as f:
        sv, v = s["VALUE"], f["VALUE"]
        n = v["ro"].shape[0]
        for c in ("ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega"):
            if sv[c].shape[0] != n:
                sys.exit("SRC と MESH の DOF 数が違う (%s: %d vs %d)" % (c, sv[c].shape[0], n))
        ro = sv["ro"][:].astype(np.float64)
        ek_old = 0.5 * (sv["roUx"][:].astype(np.float64) ** 2 + sv["roUy"][:].astype(np.float64) ** 2
                        + sv["roUz"][:].astype(np.float64) ** 2) / ro
        ei = sv["roe"][:].astype(np.float64) - ek_old  # 内部エネルギー (S1 のまま保つ)
        xyz = node_xyz(f, n)
        x, y = xyz[:, 0], xyz[:, 1]
        xs = np.linspace(0.0, LX, 8001)
        yw = np.interp(np.clip(x, 0.0, LX), xs, [wall_y(xi) for xi in xs])
        hcol = LY_TOP - yw
        eta = np.clip((y - yw) / hcol, 0.0, 1.0)
        # 列ごとに ∫ρu dy が丘頂断面 (高さ LY_TOP-H) のバルク ub と同じになる放物線 (ρ は ro0 で近似)
        ucol = ub * (LY_TOP - H) / hcol
        u = 6.0 * ucol * eta * (1.0 - eta)
        roUx = ro * u
        # 壁ノードは u=0 (eta=0/1)。index コピーで保存量を書く
        put(v, "ro", ro)
        put(v, "roUx", roUx)
        put(v, "roUy", np.zeros(n))
        put(v, "roUz", np.zeros(n))
        put(v, "roe", ei + 0.5 * ro * u * u)
        for c in ("roK", "roOmega"):
            put(v, c, sv[c][:].astype(np.float64))
    print("shaped IC -> %s (src %s, U_b(crest)=%.1f m/s, u_max(crest)=%.1f)" % (out, src, ub, 1.5 * ub))


def cmd_manifest(run):
    sys.path.insert(0, os.path.join(HERE, "..", "..", "solver_density_cuda", "tools"))
    from stage_manifest import StageManifest
    bc = open(os.path.join(run, "bcondConfig.yaml")).read()
    sm = StageManifest(run)
    for tag in ("S1_spin", "S2_dev"):
        cfg = open(os.path.join(run, "solverConfig_%s.yaml" % tag)).read()
        sm.add(tag, cfg, bc, history="residual_history_%s.csv" % tag)
    print("wrote", sm.write())


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("static"); a.add_argument("mesh"); a.add_argument("out")
    b = sub.add_parser("shaped"); b.add_argument("src"); b.add_argument("mesh"); b.add_argument("out")
    b.add_argument("--ub", type=float, default=52.0)
    c = sub.add_parser("manifest"); c.add_argument("run")
    a_ = ap.parse_args()
    if a_.cmd == "static":
        cmd_static(a_.mesh, a_.out)
    elif a_.cmd == "shaped":
        cmd_shaped(a_.src, a_.mesh, a_.out, a_.ub)
    else:
        cmd_manifest(a_.run)


if __name__ == "__main__":
    main()
