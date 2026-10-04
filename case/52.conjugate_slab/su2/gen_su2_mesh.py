#!/usr/bin/env python3
r"""case/52 の共役スラブを **SU2 の 2 ゾーン (流体 + 固体) メッシュ**にする。

plan boundary-conjugate-heat-transfer §6 V2。forge の V1 と**同じ問題**を SU2 の multizone CHT
(`MARKER_CHT_INTERFACE`) で解き、解析解と forge に突き合わせる。

問題 (case/52 README と同一):

    流体層  H=0.01 m, k_f=0.0241 W/mK, 上面 350 K 固定
    固体層  t=0.002 m, k_s=0.01 W/mK  (= t/k_s = 0.2 m2K/W), 背面 300 K 固定
    界面    y=0 で連成
    解析解  q = (350-300)/(H/k_f + t/k_s) = 81.3090 W/m2,  T_w = 300 + 0.2 q = 316.2618 K

**ゾーンごとに別ファイル**にする (`MULTIZONE_MESH= NO`)。1 ファイルに `IZONE=` で詰める形式もあるが、
ゾーンごとのマーカ名の衝突を気にせずに済むこちらを採る。

界面の節点は**流体側と固体側で一致させる** (SU2 は内挿できるが、内挿誤差を V2 の差に混ぜない)。
側面は対称面にして 1 次元性を保つ。

使い方:
  python3 case/52.conjugate_slab/su2/gen_su2_mesh.py [--nx 20] [--nyf 16] [--nys 8]
"""
from __future__ import annotations

import argparse
from pathlib import Path


def write_su2(path: Path, x0, x1, y0, y1, nx, ny, tags):
    """構造 quad の 2D メッシュを SU2 native 形式で書く。

    tags = {"bottom": name, "top": name, "left": name, "right": name}
    """
    xs = [x0 + (x1 - x0) * i / nx for i in range(nx + 1)]
    ys = [y0 + (y1 - y0) * j / ny for j in range(ny + 1)]
    nid = lambda i, j: j * (nx + 1) + i          # noqa: E731

    lines = ["NDIME= 2", f"NELEM= {nx*ny}"]
    e = 0
    for j in range(ny):
        for i in range(nx):
            lines.append(f"9 {nid(i,j)} {nid(i+1,j)} {nid(i+1,j+1)} {nid(i,j+1)} {e}")
            e += 1
    lines.append(f"NPOIN= {(nx+1)*(ny+1)}")
    for j in range(ny + 1):
        for i in range(nx + 1):
            lines.append(f"{xs[i]:.16e} {ys[j]:.16e} {nid(i,j)}")

    marks = [
        (tags["bottom"], [(nid(i, 0), nid(i + 1, 0)) for i in range(nx)]),
        (tags["top"], [(nid(i, ny), nid(i + 1, ny)) for i in range(nx)]),
        (tags["left"], [(nid(0, j), nid(0, j + 1)) for j in range(ny)]),
        (tags["right"], [(nid(nx, j), nid(nx, j + 1)) for j in range(ny)]),
    ]
    lines.append(f"NMARK= {len(marks)}")
    for name, edges in marks:
        lines.append(f"MARKER_TAG= {name}")
        lines.append(f"MARKER_ELEMS= {len(edges)}")
        for a, b in edges:
            lines.append(f"3 {a} {b}")
    path.write_text("\n".join(lines) + "\n")
    return (nx + 1) * (ny + 1), nx * ny


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nx", type=int, default=20)
    ap.add_argument("--nyf", type=int, default=16, help="流体層の分割 (forge の V1 と同じ)")
    ap.add_argument("--nys", type=int, default=8, help="固体層の分割")
    ap.add_argument("--W", type=float, default=0.02)
    ap.add_argument("--H", type=float, default=0.01, help="流体層の厚さ [m]")
    ap.add_argument("--t", type=float, default=0.002, help="固体層の厚さ [m]")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    out = Path(a.out) if a.out else Path(__file__).resolve().parent
    out.mkdir(parents=True, exist_ok=True)

    # 流体: y = 0 (界面) .. H (上面 350 K)
    nf = write_su2(out / "fluid.su2", 0.0, a.W, 0.0, a.H, a.nx, a.nyf,
                   {"bottom": "fluid_interface", "top": "fluid_top",
                    "left": "fluid_left", "right": "fluid_right"})
    # 固体: y = -t (背面 300 K) .. 0 (界面)
    ns = write_su2(out / "solid.su2", 0.0, a.W, -a.t, 0.0, a.nx, a.nys,
                   {"bottom": "solid_back", "top": "solid_interface",
                    "left": "solid_left", "right": "solid_right"})

    kf, ks = 0.0241, 0.01
    Rf, Rs = a.H / kf, a.t / ks
    q = (350.0 - 300.0) / (Rf + Rs)
    Tw = 300.0 + q * Rs
    print(f"[gen_su2_mesh] fluid.su2  {nf[0]} nodes / {nf[1]} quads  (y 0 .. {a.H})")
    print(f"[gen_su2_mesh] solid.su2  {ns[0]} nodes / {ns[1]} quads  (y {-a.t} .. 0)")
    print(f"  R_f = H/k_f = {Rf:.6f},  R_s = t/k_s = {Rs:.6f}  [m2K/W]")
    print(f"  解析解: q = {q:.4f} W/m2,  T_w = {Tw:.4f} K")
    print(f"  界面の節点: 流体 {a.nx+1} / 固体 {a.nx+1} (x 方向に一致)")


if __name__ == "__main__":
    main()
