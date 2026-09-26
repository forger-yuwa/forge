#!/usr/bin/env python3
r"""V6′ — 深いスロットの両側壁を CHT にする 2D ケースのメッシュ生成 (自己完結)。

plan [`boundary-conjugate-heat-transfer.md`](../../plans/accepted/boundary-conjugate-heat-transfer.md)
§6 V6′。**他セッションの資産に依存しない**ことが要件なので、幾何もメッシュも本ファイルで作る。

## 何を試すケースか

深いスロット ($D/W\ge16$) の**前壁と後壁に別々の背面温度**を与えて `conjugate` で解く。
深部は $Pe=|u|W/\alpha\ll1$ で伝導支配になるので、**幅方向 1 次元伝導 + 固体の直列抵抗**の
解析解に漸近しなければならない:

    q*    = (T_b1 - T_b2) / (t1/k_s1 + W/k_f + t2/k_s2)
    T_w1* = T_b1 - q* t1/k_s1        (前壁の界面温度)
    T_w2* = T_b2 + q* t2/k_s2        (後壁の界面温度)

**一様壁温だと漸近解が $T\equiv T_w$・$q\equiv0$ の自明解 (null 検査) にしかならない**ので、
壁温差をつけるのが要点。

## 幾何 (既定)

    y=H  ┌────────────────────────────────┐ top (slip)
         │                                │
    y=0  ┤──────┬──┬──────────────────────┤ plate (等温)
                │  │  ← スロット 幅 W
                │  │     前壁 = slot_front (x=0, conjugate 背面 T_b1)
                │  │     後壁 = slot_back  (x=W, conjugate 背面 T_b2)
    y=-D        └──┘  slot_bottom (等温)

流れは左から右へ M=2 の層流。超音速なので出口に背圧を与えなくてよい。

## 使い方

    python3 case/58.conjugate_slot/gen_mesh.py            # .geo と .msh を作る
    python3 case/58.conjugate_slot/gen_mesh.py --nw 60    # スロット幅方向を細かく
"""
from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent


def build_geo(a) -> str:
    W, D, H = a.W, a.D, a.H
    x0, x1 = 0.0, W                       # スロットの前壁 / 後壁
    xin, xout = -a.Lup, W + a.Ldn
    L: list[str] = []
    P = L.append
    P("// case/58 — 深いスロット両側壁 CHT (V6′)。gen_mesh.py が生成。編集しないこと。")
    P("Geometry.PointNumbers = 0;  lc = 0.005;")
    P(f"xin = {xin:.9f}; x0 = {x0:.9f}; x1 = {x1:.9f}; xout = {xout:.9f};")
    P(f"H = {H:.9f}; D = {D:.9f};")

    # --- 上側チャネル: 3 ブロック (上流 / スロット直上 / 下流)
    for i, (x, tag) in enumerate([("xin", 1), ("x0", 2), ("x1", 3), ("xout", 4)]):
        P(f"Point({tag}) = {{{x}, 0.0, 0.0, lc}};")
        P(f"Point({tag+4}) = {{{x}, H, 0.0, lc}};")
    for i in range(3):                                  # 下辺 1,2,3
        P(f"Line({i+1}) = {{{i+1}, {i+2}}};")
    for i in range(3):                                  # 上辺 4,5,6
        P(f"Line({i+4}) = {{{i+5}, {i+6}}};")
    for i in range(4):                                  # 縦 7,8,9,10
        P(f"Line({i+7}) = {{{i+1}, {i+5}}};")
    # 法線方向 (y): 壁側を細かく
    P(f"Transfinite Line {{7, 8, 9, 10}} = {a.ny} Using Progression {a.ry:.6f};")
    P(f"Transfinite Line {{1, 4}} = {a.nup} Using Progression {1.0/a.rx:.6f};")   # 上流: 開口へ寄せる
    P(f"Transfinite Line {{2, 5}} = {a.nw};")                                     # 開口幅 = スロットと同じ
    P(f"Transfinite Line {{3, 6}} = {a.ndn} Using Progression {a.rx:.6f};")       # 下流
    for i, (b, v, t, u) in enumerate([(1, 8, 4, 7), (2, 9, 5, 8), (3, 10, 6, 9)]):
        P(f"Curve Loop({i+1}) = {{{b}, {v}, -{t}, -{u}}};")
        P(f"Plane Surface({i+1}) = {{{i+1}}};  Transfinite Surface {{{i+1}}};"
          f"  Recombine Surface {{{i+1}}};")

    # --- スロット (幅 W × 深さ D)。上辺は上側ブロック 2 の下辺 (Line 2) を共有
    P(f"Point(11) = {{x0, -D, 0.0, lc}};")
    P(f"Point(12) = {{x1, -D, 0.0, lc}};")
    P("Line(11) = {11, 12};")      # 底
    P("Line(12) = {2, 11};")       # 前壁 (x=x0)
    P("Line(13) = {3, 12};")       # 後壁 (x=x1)
    P(f"Transfinite Line {{11}} = {a.nw};")
    P(f"Transfinite Line {{12, 13}} = {a.nd} Using Progression {a.rd:.6f};")
    P("Curve Loop(4) = {12, 11, -13, -2};")
    P("Plane Surface(4) = {4};  Transfinite Surface {4};  Recombine Surface {4};")

    P('Physical Curve("inlet",       1) = {7};')
    P('Physical Curve("outlet",      2) = {10};')
    P('Physical Curve("top",         3) = {4, 5, 6};')
    P('Physical Curve("plate",       4) = {1, 3};')
    P('Physical Curve("slot_front",  5) = {12};')
    P('Physical Curve("slot_back",   6) = {13};')
    P('Physical Curve("slot_bottom", 7) = {11};')
    P("Physical Surface(\"fluid\", 100) = {1, 2, 3, 4};")
    P("Mesh.RecombineAll = 1;  Mesh.ElementOrder = 1;  Mesh.MshFileVersion = 4.1;")
    return "\n".join(L) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--W", type=float, default=1.0e-3, help="スロット幅 [m]")
    ap.add_argument("--D", type=float, default=20.0e-3, help="スロット深さ [m] (D/W>=16 が要件)")
    ap.add_argument("--H", type=float, default=8.0e-3, help="上側チャネル高さ [m]")
    ap.add_argument("--Lup", type=float, default=10.0e-3, help="開口より上流の長さ [m]")
    ap.add_argument("--Ldn", type=float, default=20.0e-3, help="開口より下流の長さ [m]")
    ap.add_argument("--nw", type=int, default=41, help="幅方向の節点数 (開口とスロットで共通)")
    ap.add_argument("--nd", type=int, default=321, help="深さ方向の節点数")
    ap.add_argument("--ny", type=int, default=81, help="チャネル法線方向の節点数")
    ap.add_argument("--nup", type=int, default=81, help="上流方向の節点数")
    ap.add_argument("--ndn", type=int, default=121, help="下流方向の節点数")
    ap.add_argument("--ry", type=float, default=1.06, help="法線方向の伸長比 (壁側を細かく)")
    ap.add_argument("--rx", type=float, default=1.02, help="流れ方向の伸長比")
    ap.add_argument("--rd", type=float, default=1.004, help="深さ方向の伸長比 (開口側を細かく)")
    ap.add_argument("--name", default="slot_wd20")
    a = ap.parse_args()

    if a.D / a.W < 16.0:
        raise SystemExit(f"D/W = {a.D/a.W:.2f} < 16 — V6′ の要件を満たさない")

    mesh = HERE / "mesh"
    mesh.mkdir(parents=True, exist_ok=True)
    geo = mesh / f"{a.name}.geo"
    geo.write_text(build_geo(a))
    print(f"[gen_mesh] {geo}")
    print(f"  W {a.W*1e3:.3f} mm / D {a.D*1e3:.3f} mm  (D/W = {a.D/a.W:.1f})")
    print(f"  幅方向 {a.nw} 節点 -> Δx {a.W/(a.nw-1)*1e6:.1f} µm")
    print(f"  深さ方向 {a.nd} 節点 -> 開口側 Δy {a.D/(a.nd-1)*1e6:.1f} µm (伸長 {a.rd})")

    msh = mesh / f"{a.name}.msh"
    r = subprocess.run(["gmsh", "-2", str(geo), "-o", str(msh)],
                       capture_output=True, text=True)
    tail = "\n".join(r.stdout.strip().split("\n")[-4:])
    print(tail)
    if r.returncode != 0:
        print(r.stderr[-2000:])
        raise SystemExit("gmsh が失敗した")
    print(f"[gen_mesh] {msh}")


if __name__ == "__main__":
    main()
