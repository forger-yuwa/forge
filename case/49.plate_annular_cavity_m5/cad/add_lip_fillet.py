#!/usr/bin/env python3
r"""開口リップに**フィレット**を入れる (plan §5.1 #17)。

リップは 90° の凸エッジ (固体が流体側へ尖っている) で、$q''\propto s^{-1/2}$ と発散するため
ピーク熱流束が格子収束しない (§4.4.2)。実物にはエッジ半径があるので、**形状側で特異点を除いて**
総入熱の半径依存を出す。

`build_hex_mesh.py` は z=0 の 2D ブロッキングを上下に押し出すので、$r$ を $z$ の関数に
できない = 形状生成の段階ではフィレットを作れない。そこで**生成済みメッシュの節点を変位**させる。

外筒リップ (r=Ro, z=0) の場合、固体は r>Ro かつ z<0 の四半分。その角を半径 rf で丸めると
壁面は中心 (Ro+rf, -rf) の円弧になる:

    r_w(z) = Ro + rf - sqrt(rf^2 - (z+rf)^2)      (-rf <= z <= 0)
    z_w(r) = -rf + sqrt(rf^2 - (r-Ro-rf)^2)       (Ro <= r <= Ro+rf)

壁がこれだけ動くので、内部節点も**減衰させながら**同じ向きに動かす (減衰長 L)。
L を rf より十分大きく取れば歪みは rf/L 程度に収まる。内円柱リップ (r=r_i, z=0) も同様で、
こちらは固体が r<r_i 側なので変位の向きが逆になる。

usage (mesh venv):
  python3 cad/add_lip_fillet.py --in cad/hex.msh --out cad/hex_rf0.1.msh --rf 0.1 [--blend 10]
"""
import argparse
import json
import math
import os
import sys
from pathlib import Path

import gmsh
import numpy as np

HERE = Path(__file__).resolve().parent
CASE = HERE.parent
sys.path.insert(0, str(CASE))
sys.path.insert(0, str(CASE / "tools"))


def smooth(t):
    """0<=t<=1 で 1 -> 0 に落ちる滑らかな重み (端で 1 階微分も 0)。"""
    t = np.clip(t, 0.0, 1.0)
    return 1.0 - (3.0 * t ** 2 - 2.0 * t ** 3)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="src", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--rf", type=float, required=True, help="フィレット半径 [mm]")
    ap.add_argument("--blend", type=float, default=10.0,
                    help="減衰長 / rf。大きいほど歪みが小さい (既定 10)")
    ap.add_argument("--manifest", default=os.environ.get("CASE49_MANIFEST", "manifest.json"))
    a = ap.parse_args()
    rf = a.rf * 1e-3
    L = a.blend * rf

    man = json.loads((CASE / a.manifest).read_text())
    G = man["geometry"]
    Ro, Ri, off = float(G["Ro"]), float(G["Ri"]), float(G["x_off"])
    print("フィレット rf = %.3g mm, 減衰長 L = %.3g mm  (Ro %.1f / Ri %.1f / 偏心 %.2f mm)"
          % (a.rf, L * 1e3, Ro * 1e3, Ri * 1e3, off * 1e3))

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.open(a.src)
    tags, coord, _ = gmsh.model.mesh.getNodes()
    P = np.array(coord).reshape(-1, 3)
    x, y, z = P[:, 0], P[:, 1], P[:, 2]

    r_out = np.hypot(x, y)                       # 外筒軸からの半径
    r_in = np.hypot(x - off, y)                  # 内円柱軸からの半径
    moved = 0

    # 角 (s=0, z=0) を中心 (rf, -rf) 半径 rf の円弧で置き換える。(s = 壁からの符号付き距離、
    # 固体側が正。z は鉛直)。円弧は下端 (0, -rf) = 角度 180 度 から 上端 (rf, 0) = 90 度。
    #
    # **素朴に「同じ z の円弧点へ」「同じ s の円弧点へ」と写すと写像が反転する**
    # (すきま側の上端が円弧の上端へ、プレート側の内端が下端へ飛び、両者が入れ替わる)。
    # 壁に沿った順序を保つため、**円弧上の角度を弧長パラメータで割り当てる**:
    #   すきま側 (0, z), t=(z+rf)/rf in [0,1]  ->  phi = 180 - 45 t   (t=0 で不動、t=1 で中点)
    #   プレート側 (s, 0), u=s/rf in [0,1]     ->  phi = 90 + 45 (1-u) (u=1 で不動、u=0 で中点)
    # 両端 (z=-rf, s=rf) が不動点で、角では双方が円弧の中点 (135 度) に集まる。
    def arc_point(phi):
        return rf + rf * np.cos(phi), -rf + rf * np.sin(phi)

    dx = np.zeros_like(x); dy = np.zeros_like(y); dz = np.zeros_like(z)

    for name, rr, sgn, cx in (("外筒", r_out, +1.0, 0.0), ("内円柱", r_in, -1.0, off)):
        R = Ro if sgn > 0 else Ri
        s = sgn * (rr - R)                       # 固体側が正
        er_x = (x - cx) / np.maximum(rr, 1e-12)
        er_y = y / np.maximum(rr, 1e-12)
        # (a) すきま側 (z<=0, 流体 s<0)
        m = (z <= 0.0) & (z > -rf) & (s <= 0.0) & (s > -L)
        n_a = int(m.sum())
        if n_a:
            t = (z[m] + rf) / rf
            phi = np.radians(180.0 - 45.0 * t)
            ps, pz = arc_point(phi)
            w = smooth(-s[m] / L)
            dx[m] += sgn * ps * w * er_x[m]
            dy[m] += sgn * ps * w * er_y[m]
            dz[m] += (pz - z[m]) * w
        # (b) プレート / 円柱上面側 (z>=0, 固体 s>0, s<rf)
        m2 = (z >= 0.0) & (z < L) & (s > 0.0) & (s <= rf)
        n_b = int(m2.sum())
        if n_b:
            u = s[m2] / rf
            phi = np.radians(90.0 + 45.0 * (1.0 - u))
            ps, pz = arc_point(phi)
            w2 = smooth(z[m2] / L)
            dx[m2] += sgn * (ps - s[m2]) * w2 * er_x[m2]
            dy[m2] += sgn * (ps - s[m2]) * w2 * er_y[m2]
            dz[m2] += pz * w2
        print("  %s: すきま側 %d 点 / 上面側 %d 点" % (name, n_a, n_b))

    # **局所セル厚を超える変位は出さない** (2026-09-20)。開口直下の第一層は 8 um しかなく、
    # 鉛直変位の最大 rf*(1-sin135deg) = 0.293*rf がこれを超えるとセルが反転し、
    # 変換器の双対面閉性 (median-dual) が破れて黙って使えないメッシュになる。
    slit = (np.abs(r_out - Ro) < 0.15 * (Ro - Ri)) & (z <= 0.0) & (z > -5.0 * rf)
    if slit.sum() > 10:
        zl = np.unique(np.round(z[slit], 9))[::-1]
        first = float(abs(zl[1] - zl[0])) if len(zl) > 1 else float("inf")
    else:
        first = float("inf")
    dz_max = rf * (1.0 - math.sin(math.radians(135.0)))
    print("  開口直下の第一層 %.2f um / 鉛直変位の最大 %.2f um" % (first * 1e6, dz_max * 1e6))
    if dz_max > 0.3 * first:
        gmsh.finalize()
        raise SystemExit(
            "REFUSED: rf = %.4g mm では鉛直変位 %.2f um が第一層 %.2f um の 30 %% を超える。\n"
            "  この押し出しメッシュでは rf <= %.4g mm までしか節点変位で作れない。\n"
            "  実物のエッジ半径 (0.05-0.5 mm) を入れるには、リップ周りを (r,z) 面の O グリッドで\n"
            "  包むようブロッキングを組み直す必要がある (build_hex_mesh.py の設計変更)。"
            % (a.rf, dz_max * 1e6, first * 1e6, 0.3 * first / (1.0 - math.sin(math.radians(135.0))) * 1e3))

    P2 = P + np.stack([dx, dy, dz], axis=1)
    # **動かした節点だけ書き戻す** (gmsh の API は `setNode` (単数) しか無い)。
    mv = (np.abs(dx) + np.abs(dy) + np.abs(dz)) > 0.0
    for t, q in zip(tags[mv], P2[mv]):
        gmsh.model.mesh.setNode(int(t), [float(q[0]), float(q[1]), float(q[2])], [])
    print("  最大変位 %.4g mm" % (np.max(np.linalg.norm(P2 - P, axis=1)) * 1e3))
    gmsh.write(a.out)
    gmsh.finalize()
    print("動かした節点 %d 点 -> %s" % (int(mv.sum()), a.out))


if __name__ == "__main__":
    main()
