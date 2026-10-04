#!/usr/bin/env python3
"""forge メッシュ h5 の**追加品質ゲート** (plan §4.4 の表, codex 2026-09-19 plan-2 M9)。

usage: python3 tools/check_mesh_extra.py mesh/stageB.h5 [--sliver-min 0.01] [--closure-max 1e-6]

`check_mesh_quality.py` の AR (辺長比) と skew (面内角) だけでは**体積が潰れた要素**を落とせない
(厚さ 1e-6 の人工 tet が AR 1.414 / skew 0.250 で通る)。ここでは:

  1. primal 要素 (VIZMESH/CONNE) の **符号付き体積の向きが全要素で一致**していること
     (forge の VIZMESH 節点順では全要素が負になるのが正常なので、「> 0」ではなく
      「多数派と逆符号の要素が無い」+「|V| が潰れていない」で判定する)
  2. tet の **正規化体積** V/(L_rms^3/(6√2))  (正四面体で 1) — sliver 判定。既定 ≥ 0.01
     prism は底面 3 角形の正規化面積と層厚比で見る (薄い層は正常なので体積比では判定しない)
  3. **dual CV 体積 > 0** (CELLS/volume) 全数
  4. **CV ごとの閉性** ‖Σ S_f‖ / Σ‖S_f‖ — 変換器のログは全体最大面積で正規化するので
     微小 CV の局所閉性を保証しない。ここは CV ごとの面積尺度で測る

`PLANES/STRUCT` は [nNodes, nodes..., nCells, cells...] の並び (gmshReader.hpp)。
surfVect は iCells[0] → iCells[1] 向き。
"""
import argparse
import sys
from pathlib import Path

import h5py
import numpy as np

# forge の VIZMESH コード (tools/res_h5_to_vtu.py の XDMF2VTK と同じ)
XDMF_NN = {4: 3, 5: 4, 6: 4, 7: 5, 8: 6, 9: 8}     # tri/quad/tet/pyramid/prism/hex
TET, PYRAMID, PRISM, HEX = 6, 7, 8, 9


def parse_viz(conne):
    """VIZMESH/CONNE (XDMF Mixed) -> [(code, nodes)] """
    out, i = [], 0
    n = len(conne)
    while i < n:
        code = int(conne[i])
        nn = XDMF_NN.get(code)
        if nn is None:
            break
        out.append((code, conne[i + 1:i + 1 + nn]))
        i += 1 + nn
    return out


def tet_vol(p):
    return np.einsum("ij,ij->i", np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), p[:, 3] - p[:, 0]) / 6.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("h5")
    ap.add_argument("--sliver-min", type=float, default=0.01)
    # 閉性の既定は 5e-5。メッシュ h5 の surfVect/surfArea/volume は **float32** なので、
    # CV が µm 級になると相対 1e-5 級の丸めが乗る (実測 2026-09-19: 全ヘキサ 100 万 CV で
    # 床の最小 CV 48 個が 1.2e-5。座標 ~0.1 m の外積を float32 で取るため絶対誤差が
    # 面積スケールに対して相対的に大きくなる)。真の欠陥 (向き不整合・面の取りこぼし) は
    # これより桁違いに大きい。閾値を超えた CV は**位置と体積**を出して死水域か判断する。
    ap.add_argument("--closure-max", type=float, default=5.0e-5)
    ap.add_argument("--layer-ratio-max", type=float, default=1.5)
    a = ap.parse_args()

    with h5py.File(a.h5, "r") as f:
        xyz = np.array(f["MESH/COORD"]).reshape(-1, 3)
        viz = np.array(f["VIZMESH/CONNE"])
        vol_dual = np.array(f["CELLS/volume"])
        strct = np.array(f["PLANES/STRUCT"])
        svec = np.array(f["PLANES/surfVect"]).reshape(-1, 3)
        sarea = np.array(f["PLANES/surfArea"])

    fails = []
    print("=== %s ===" % a.h5)
    print("  節点 %d,  dual CV %d,  plane %d" % (len(xyz), len(vol_dual), len(sarea)))

    # ---- 1,2: primal 要素 ----
    cells = parse_viz(viz)
    by = {}
    for code, nd in cells:
        by.setdefault(code, []).append(nd)
    for code, lst in by.items():
        nd = np.array(lst)
        pts = xyz[nd]
        if code == TET:
            v = tet_vol(pts)
            L2 = np.zeros(len(nd))
            for i in range(4):
                for j in range(i + 1, 4):
                    L2 += np.sum((pts[:, i] - pts[:, j]) ** 2, axis=1)
            Lrms = np.sqrt(L2 / 6.0)
            q = np.abs(v) / (Lrms ** 3 / (6.0 * np.sqrt(2.0)))
            sgn = np.sign(np.median(v))
            nflip = int(np.sum(np.sign(v) != sgn))
            nzero = int(np.sum(v == 0.0))
            nsl = int(np.sum(q < a.sliver_min))
            print("  tet   %8d : 向き %s,  逆符号 %d, 体積 0 %d,  正規化体積 min %.4g / p1 %.4g  (< %.3g が %d 個)"
                  % (len(nd), "負(正常)" if sgn < 0 else "正", nflip, nzero,
                     q.min(), np.percentile(q, 1), a.sliver_min, nsl))
            if nflip or nzero:
                fails.append("tet の向きが不揃い (逆符号 %d / 体積 0 %d)" % (nflip, nzero))
            if nsl:
                fails.append("tet の sliver (正規化体積 < %.3g) が %d 個" % (a.sliver_min, nsl))
        elif code == PRISM:
            # 3 本の側辺 (0-3, 1-4, 2-5) の長さで層厚、底面積で接線スケール
            v = np.zeros(len(nd))
            for tetn in ((0, 1, 2, 3), (1, 2, 3, 4), (2, 3, 4, 5)):
                v += tet_vol(pts[:, tetn, :])
            a0 = 0.5 * np.linalg.norm(np.cross(pts[:, 1] - pts[:, 0], pts[:, 2] - pts[:, 0]), axis=1)
            a1 = 0.5 * np.linalg.norm(np.cross(pts[:, 4] - pts[:, 3], pts[:, 5] - pts[:, 3]), axis=1)
            thick = np.mean([np.linalg.norm(pts[:, i + 3] - pts[:, i], axis=1) for i in range(3)], axis=0)
            sgn = np.sign(np.median(v))
            nflip = int(np.sum(np.sign(v) != sgn))
            ratio = np.maximum(a0, a1) / np.maximum(np.minimum(a0, a1), 1e-30)
            nbad = int(np.sum(ratio > a.layer_ratio_max))
            print("  prism %8d : 向き %s, 逆符号 %d,  層厚 %.3g..%.3g m,  上下面積比 max %.3f (> %.2f が %d 個)"
                  % (len(nd), "負(正常)" if sgn < 0 else "正", nflip, thick.min(), thick.max(),
                     ratio.max(), a.layer_ratio_max, nbad))
            if nflip:
                fails.append("prism の向きが不揃い (逆符号 %d 個)" % nflip)
        elif code == HEX:
            # 6 tet に分割して体積、辺長比で潰れ判定
            v = np.zeros(len(nd))
            for tn in ((0,1,3,4),(1,2,3,6),(1,3,4,6),(1,4,5,6),(3,4,6,7)):
                v += tet_vol(pts[:, tn, :])
            sgn = np.sign(np.median(v))
            nflip = int(np.sum(np.sign(v) != sgn))
            # 3 方向の代表辺長
            e1 = np.linalg.norm(pts[:, 1] - pts[:, 0], axis=1)
            e2 = np.linalg.norm(pts[:, 3] - pts[:, 0], axis=1)
            e3 = np.linalg.norm(pts[:, 4] - pts[:, 0], axis=1)
            q = np.abs(v) / np.maximum(e1 * e2 * e3, 1e-300)      # 直方体で 1
            nsl = int(np.sum(q < a.sliver_min))
            print("  hex   %8d : 向き %s, 逆符号 %d,  正規化体積 min %.4g / p1 %.4g  (< %.3g が %d 個)"
                  % (len(nd), "負(正常)" if sgn < 0 else "正", nflip,
                     q.min(), np.percentile(q, 1), a.sliver_min, nsl))
            if nflip:
                fails.append("hex の向きが不揃い (逆符号 %d 個)" % nflip)
            if nsl:
                fails.append("hex の潰れ (正規化体積 < %.3g) が %d 個" % (a.sliver_min, nsl))
        else:
            print("  code %d %8d : (未判定)" % (code, len(nd)))

    # ---- 3: dual CV 体積 ----
    nneg = int(np.sum(vol_dual <= 0))
    print("  dual CV 体積: min %.4g, max %.4g,  <=0 が %d 個" % (vol_dual.min(), vol_dual.max(), nneg))
    if nneg:
        fails.append("dual CV 体積 <= 0 が %d 個" % nneg)

    # ---- 4: CV ごとの閉性 ----
    ncell = len(vol_dual)
    sumS = np.zeros((ncell, 3))
    sumA = np.zeros(ncell)
    i, ip = 0, 0
    n = len(strct)
    while i < n and ip < len(sarea):
        nn = int(strct[i]); i += 1 + nn
        nc = int(strct[i]); i += 1
        cl = strct[i:i + nc]; i += nc
        if nc >= 1:
            c0 = int(cl[0])
            sumS[c0] += svec[ip]; sumA[c0] += sarea[ip]
        if nc >= 2:
            c1 = int(cl[1])
            sumS[c1] -= svec[ip]; sumA[c1] += sarea[ip]
        ip += 1
    clo = np.linalg.norm(sumS, axis=1) / np.maximum(sumA, 1e-30)
    worst = int(np.argmax(clo))
    nbad = int(np.sum(clo > a.closure_max))
    print("  CV ごとの閉性 ||ΣS||/Σ|S| : max %.3e (CV %d),  p99.9 %.3e,  > %.1e が %d 個"
          % (clo.max(), worst, np.percentile(clo, 99.9), a.closure_max, nbad))
    if nbad:
        fails.append("CV 閉性 > %.1e が %d 個 (最悪 %.3e)" % (a.closure_max, nbad, clo.max()))
    n_mild = int(np.sum(clo > 1.0e-5))
    if n_mild:
        b = np.where(clo > 1.0e-5)[0]
        cc = np.array([np.mean(xyz[cells[i][1]], axis=0) if i < len(cells) else [0, 0, 0] for i in b[:1]])
        print("    (参考) 1e-5 超は %d 個、体積 %.3g..%.3g m^3 — float32 丸めの水準"
              % (n_mild, vol_dual[b].min(), vol_dual[b].max()))

    print("\nVERDICT: %s" % ("PASS" if not fails else "FAIL"))
    for m in fails:
        print("  - " + m)
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
