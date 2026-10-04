"""`check_wall_resolution.py` の距離計算の回帰試験 (解析的に $y_1$ が分かる配置)。

壁面ごとの第一内部点を `PLANES/STRUCT` の接続から法線方向に引く部分を、**合成メッシュ**で検証する。
ここが壊れると y1+ が桁で変わる (2026-09-19 に実際に 8 倍ずれた)。

    python3 solver_density_cuda/tools/test_wall_resolution.py

試す配置:

- **一様直交**: 壁から 1e-5 m に第一内部点 -> y1 = 1e-5
- **非一様**: 壁ごとに第一層厚が違う -> 壁ごとに別の y1 が出る (全域共通値にしない)
- **斜交**: 法線から 60° 傾いた隣接しかない -> `--align-min 0.5` で**評価不能**にする
- **角**: 2 つの壁面を共有するノード -> 平均法線に最も沿う隣接を選ぶ
"""
import os
import sys

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import check_wall_resolution as W      # noqa: E402

ok = True


def chk(name, cond, detail=""):
    global ok
    ok = ok and bool(cond)
    print("  [%s] %-52s %s" % ("OK " if cond else "NG ", name, detail))


def build(tmp, coords, planes, bconds):
    """planes: [(cell0, cell1, normal)]、bconds: {physID: [wall dof]}"""
    st = []
    svec = []
    for c0, c1, n in planes:
        st += [3, 0, 0, 0]            # nNodes とダミー節点 (距離計算では使わない)
        st += [2, c0, c1]
        svec.append(n)
    with h5py.File(tmp, "w") as f:
        f.create_dataset("MESH/COORD", data=np.asarray(coords, float).ravel())
        f.create_dataset("CELLS/centCoords", data=np.asarray(coords, float).ravel())
        f.create_dataset("PLANES/STRUCT", data=np.asarray(st, np.int32))
        f.create_dataset("PLANES/surfVect", data=np.asarray(svec, float).ravel())
        f.create_dataset("PLANES/surfArea", data=np.ones(len(planes)))
        for pid, dofs in bconds.items():
            if pid == "_planes":
                continue
            f.create_dataset("BCONDS/%d/iCells" % pid, data=np.asarray(dofs, np.int32))
            f.create_dataset("BCONDS/%d/iPlanes" % pid,
                             data=np.asarray(bconds["_planes"][pid], np.int32))
    return tmp


def main():
    import tempfile
    tmp = os.path.join(tempfile.mkdtemp(), "m.h5")

    # --- 一様直交: dof0 が壁 (z=0)、dof1 が z=1e-5 ---
    coords = [[0, 0, 0], [0, 0, 1e-5], [0, 0, 3e-5]]
    planes = [(0, 1, [0, 0, -1.0]), (1, 2, [0, 0, -1.0])]
    build(tmp, coords, planes, {1: [0], "_planes": {1: [0]}})
    nb = W.dof_neighbors(tmp)
    wd = W.wall_first_distance(tmp, 1, nb, 0.5, coords=np.array(coords))
    chk("一様直交: y1 = 1e-5", abs(wd[0][0] - 1e-5) < 1e-12 and wd[0][2],
        "y1=%.4e ok=%s" % (wd[0][0], wd[0][2]))

    # --- 非一様: 2 つの壁 dof が別々の第一層厚を持つ ---
    coords = [[0, 0, 0], [0, 0, 1e-5], [1, 0, 0], [1, 0, 4e-5]]
    planes = [(0, 1, [0, 0, -1.0]), (2, 3, [0, 0, -1.0])]
    build(tmp, coords, planes, {1: [0, 2], "_planes": {1: [0, 1]}})
    nb = W.dof_neighbors(tmp)
    wd = W.wall_first_distance(tmp, 1, nb, 0.5, coords=np.array(coords))
    chk("非一様: 壁ごとに別の y1 (1e-5 と 4e-5)",
        abs(wd[0][0] - 1e-5) < 1e-12 and abs(wd[2][0] - 4e-5) < 1e-12,
        "y1=%.4e, %.4e" % (wd[0][0], wd[2][0]))

    # --- 斜交: 法線から 60 度 (cos=0.5 未満は評価不能) ---
    d = 1e-5
    coords = [[0, 0, 0], [d * np.sin(np.radians(75)), 0, d * np.cos(np.radians(75))]]
    planes = [(0, 1, [0, 0, -1.0])]
    build(tmp, coords, planes, {1: [0], "_planes": {1: [0]}})
    nb = W.dof_neighbors(tmp)
    wd = W.wall_first_distance(tmp, 1, nb, 0.5, coords=np.array(coords))
    chk("斜交 75 度 (沿い 0.26) -> 評価不能", not wd[0][2], "ok=%s" % wd[0][2])

    # --- 角: 2 面を共有 (法線が平均されて 45 度)、真下と真横の隣接から選ぶ ---
    coords = [[0, 0, 0], [0, 0, 2e-5], [2e-5, 0, 0]]
    planes = [(0, 1, [0, 0, -1.0]), (0, 2, [-1.0, 0, 0])]
    build(tmp, coords, planes, {1: [0], "_planes": {1: [0, 1]}})
    nb = W.dof_neighbors(tmp)
    wd = W.wall_first_distance(tmp, 1, nb, 0.5, coords=np.array(coords))
    chk("角 (2 面共有): 平均法線に沿う隣接を選び評価できる", wd[0][2],
        "y1=%.4e ok=%s" % (wd[0][0], wd[0][2]))

    print("\nVERDICT: %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
