"""冷却壁用の格子 (forge の nozzle.msh、Gmsh 4.1、座標 17 桁) を SU2 形式に変換する (plan tooling-nozzle-isothermal-wall-chain §5.1 #26)。
座標は 17 桁のまま書き、SU2 と forge が同じ節点を解くようにする (gmsh の su2 書き出しの桁数に頼らない)。
書いた後に読み直して、節点の座標が msh と一致すること (完全一致) と、マーカーの要素数を確かめる。
usage: python3 cold_su2_mesh.py <nozzle.msh> <nozzle.su2>
"""
import sys
import numpy as np


def read_msh41(path):
    L = open(path).read().split("\n")
    names = {}
    i = L.index("$PhysicalNames")
    for k in range(int(L[i + 1])):
        d, tag, nm = L[i + 2 + k].split(maxsplit=2)
        names[int(tag)] = nm.strip('"')
    # 曲線の実体 → 物理タグ
    i = L.index("$Entities"); n0, n1, n2, n3 = map(int, L[i + 1].split()); j = i + 2 + n0
    curve_phys = {}
    for k in range(n1):
        a = L[j + k].split(); tag = int(a[0]); nphys = int(a[7]); curve_phys[tag] = int(a[8]) if nphys else None
    i = L.index("$Nodes"); nb, nn = map(int, L[i + 1].split()[:2]); j = i + 2
    tags, X = [], []
    while not L[j].startswith("$EndNodes"):
        d, ent, par, n = map(int, L[j].split()); j += 1
        t = [int(L[j + k]) for k in range(n)]; j += n
        X += [list(map(float, L[j + k].split()[:3])) for k in range(n)]; j += n
        tags += t
    order = np.argsort(tags); tags = np.array(tags)[order]; X = np.array(X)[order]
    if not np.array_equal(tags, np.arange(1, len(tags) + 1)):
        raise SystemExit("節点のタグが 1..N の連番でない")
    i = L.index("$Elements"); j = i + 2
    quads, lines = [], {}
    while not L[j].startswith("$EndElements"):
        d, ent, et, n = map(int, L[j].split()); j += 1
        for k in range(n):
            a = list(map(int, L[j + k].split()))
            if et == 3:
                quads.append(a[1:5])
            elif et == 1:
                lines.setdefault(names[curve_phys[ent]], []).append(a[1:3])
        j += n
    return X, np.array(quads) - 1, {k: np.array(v) - 1 for k, v in lines.items()}


def write_su2(path, X, quads, lines):
    with open(path, "w") as f:
        f.write("NDIME= 2\n")
        f.write(f"NELEM= {len(quads)}\n")
        for k, q in enumerate(quads):
            f.write(f"9 {q[0]} {q[1]} {q[2]} {q[3]} {k}\n")
        f.write(f"NPOIN= {len(X)}\n")
        for k, x in enumerate(X):
            f.write(f"{x[0]:.17g} {x[1]:.17g} {k}\n")
        f.write(f"NMARK= {len(lines)}\n")
        for nm in ("inlet", "outlet", "wall", "axis"):
            e = lines[nm]
            f.write(f"MARKER_TAG= {nm}\nMARKER_ELEMS= {len(e)}\n")
            for a in e:
                f.write(f"3 {a[0]} {a[1]}\n")


def read_su2_points(path):
    L = open(path).read().split("\n")
    i = next(k for k, l in enumerate(L) if l.startswith("NPOIN="))
    n = int(L[i].split("=")[1].split()[0])
    return np.array([list(map(float, L[i + 1 + k].split()[:2])) for k in range(n)])


if __name__ == "__main__":
    X, quads, lines = read_msh41(sys.argv[1])
    # 四角形の向き (反時計回り) を確かめる
    P = X[quads][:, :, :2]
    area = 0.5 * np.sum(P[:, :, 0] * np.roll(P[:, :, 1], -1, axis=1) - np.roll(P[:, :, 0], -1, axis=1) * P[:, :, 1], axis=1)
    if not np.all(area > 0):
        raise SystemExit(f"反時計回りでない四角形が {int(np.sum(area <= 0))} 個")
    write_su2(sys.argv[2], X, quads, lines)
    Y = read_su2_points(sys.argv[2])
    same = np.array_equal(Y, X[:, :2])
    print(f"節点 {len(X)}、四角形 {len(quads)}、マーカー {{{', '.join(f'{k}: {len(v)}' for k, v in lines.items())}}}、"
          f"読み直しの座標が msh と完全一致 {same}、最小面積 {area.min():.3e} m²")
    if not same:
        raise SystemExit("座標が一致しない")
