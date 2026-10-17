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

    # --- 軸対称の判定 (CLI、--weight area): 旧キー physProp.isAxisymmetric だけの run も軸対称 (2πr の重み) ---
    # solver と同じ優先順 (mesh が優先、無ければ旧 physProp; plan tooling-generator-defaults-and-launch-gate §5.1 #21)。
    # 2026-10-15 までは mesh.isAxisymmetric だけを見て、旧キーだけの run の面積を平面 (単位奥行き) で数えていた
    import subprocess
    R, dx = 0.5, 1e-3
    for tag, mesh_key, pp_key, want_axi in (("mesh 1", ", isAxisymmetric: 1", "", True),
                                            ("旧 physProp だけ 1 (負例: 旧実装は平面)", "", ", isAxisymmetric: 1", True),
                                            ("mesh 0・physProp 1 (mesh が優先)", ", isAxisymmetric: 0", ", isAxisymmetric: 1", False),
                                            ("どちらも無し", "", "", False)):
        rd = tempfile.mkdtemp()
        coords = [[i * dx, R, 0] for i in range(3)] + [[i * dx, R - 1e-5, 0] for i in range(3)]
        planes = [(i, i + 3, [0, 1.0, 0]) for i in range(3)]
        build(os.path.join(rd, "mesh.h5"), coords, planes, {1: [0, 1, 2], "_planes": {1: [0, 1, 2]}})
        with h5py.File(os.path.join(rd, "res_wall_1_10.h5"), "w") as f:
            f.create_dataset("MESH/COORD", data=np.asarray(coords[:3], float).ravel())
            f.create_dataset("MESH/CONNE", data=np.array([2, 2, 0, 1, 2, 2, 1, 2], np.int32))
            for k, v in (("twall_x", 1.0), ("twall_y", 0.0), ("twall_z", 0.0), ("ro", 1.0), ("Ts", 300.0)):
                f.create_dataset("VALUE/" + k, data=np.full(3, v))
        with open(os.path.join(rd, "solverConfig.yaml"), "w") as f:
            f.write('mesh: {discretization: "node", meshFileName: "mesh.h5"%s}\nphysProp: {viscMethod: 0, visc: 1.8e-5%s}\n'
                    % (mesh_key, pp_key))
        p = subprocess.run([sys.executable, os.path.join(HERE, "check_wall_resolution.py"), rd, "--weight", "area",
                            "--profile-csv", os.path.join(rd, "prof.csv")], capture_output=True, text=True)
        try:
            wsum = float(np.loadtxt(os.path.join(rd, "prof_wall.csv"), delimiter=",", skiprows=1, ndmin=2)[:, -1].sum())
        except OSError:
            wsum = float("nan")
        want = 2.0 * dx * (2.0 * np.pi * R if want_axi else 1.0)
        chk("軸対称の判定 (%s): 面積の重みの和" % tag, p.returncode == 0 and abs(wsum / want - 1.0) < 1e-12,
            "和 %.6e (期待 %.6e、%s)、rc %d" % (wsum, want, "2πr" if want_axi else "単位奥行き", p.returncode))
        chk("  run_axisymmetric (%s)" % tag, W.run_axisymmetric(
            {"mesh": dict([("isAxisymmetric", int(mesh_key[-1]))] if mesh_key else []),
             "physProp": dict([("isAxisymmetric", int(pp_key[-1]))] if pp_key else [])}) is want_axi)

    print("\nVERDICT: %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
