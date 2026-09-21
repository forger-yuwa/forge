#!/usr/bin/env python3
"""node (median-dual) メッシュの双対 CV を **CV ごと**に検査する。

変換器のログ (`[buildMedianDual] closure max|sum dS|`) は全域の最大面積で正規化した 1 つの数で、双対体積は和しか出ない。
局所の反転や、1 つの CV だけ閉じていない不具合はそこに埋もれる。本ツールは変換後の h5 から

  * 各 CV の閉性  |Σ S_f| / Σ |S_f|   (内部面は owner に +S、neighbor に -S。境界面は owner に +S)
  * 各 CV の体積  CELLS/volume > 0 かつ有限

を出し、`VERDICT: PASS / FAIL` を返す。判定不能 (面が 1 枚も無い CV、データセットの欠落) は FAIL。

usage: check_dual_closure.py MESH.h5 [--tol 1e-5] [--worst 5]
"""
import sys, argparse
import numpy as np
import h5py


def parse_struct(st, nplane):
    """PLANES/STRUCT = [nNodes, nodes..., nCells, cells...] の繰り返し → (owner, neighbor)。neighbor = -1 は境界面"""
    own = np.full(nplane, -1, np.int64); nei = np.full(nplane, -1, np.int64)
    i = ip = 0; n = len(st)
    while i < n and ip < nplane:
        i += 1 + int(st[i]); nc = int(st[i]); i += 1
        own[ip] = st[i]
        if nc >= 2: nei[ip] = st[i + 1]
        i += nc; ip += 1
    return own, nei, ip


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("mesh"); ap.add_argument("--tol", type=float, default=1e-5)
    ap.add_argument("--worst", type=int, default=5); a = ap.parse_args()
    with h5py.File(a.mesh, "r") as f:
        for k in ("PLANES/STRUCT", "PLANES/surfVect", "CELLS/volume", "CELLS/centCoords"):
            if k not in f:
                print("VERDICT: FAIL (判定不能: %s が無い)" % k); return 1
        S = np.array(f["PLANES/surfVect"], dtype=np.float64).reshape(-1, 3); vol = np.array(f["CELLS/volume"], dtype=np.float64)
        cc = np.array(f["CELLS/centCoords"]).reshape(-1, 3); st = np.array(f["PLANES/STRUCT"])
    ncv = len(vol); own, nei, got = parse_struct(st, len(S))
    if got != len(S) or (own < 0).any() or own.max() >= ncv or nei.max() >= ncv:
        print("VERDICT: FAIL (判定不能: PLANES/STRUCT を %d / %d 面しか読めない、または CV 番号が範囲外)" % (got, len(S))); return 1
    inner = nei >= 0; A = np.linalg.norm(S, axis=1)
    sumS = np.zeros((ncv, 3)); sumA = np.zeros(ncv); nf = np.zeros(ncv, np.int64)
    for d in range(3):
        sumS[:, d] = np.bincount(own, S[:, d], ncv) - np.bincount(nei[inner], S[inner, d], ncv)
    sumA = np.bincount(own, A, ncv) + np.bincount(nei[inner], A[inner], ncv)
    nf = np.bincount(own, minlength=ncv) + np.bincount(nei[inner], minlength=ncv)
    noface = int((nf == 0).sum()); ok_f = nf > 0
    clo = np.full(ncv, np.inf); clo[ok_f] = np.linalg.norm(sumS[ok_f], axis=1) / sumA[ok_f]
    badv = ~np.isfinite(vol) | (vol <= 0)
    print("CV %d / 双対面 %d (内部 %d, 境界 %d)" % (ncv, len(S), int(inner.sum()), int((~inner).sum())))
    print("体積: min %.3e  max %.3e  非正または非有限 %d" % (np.nanmin(vol), np.nanmax(vol), int(badv.sum())))
    fin = clo[ok_f]
    print("閉性 |ΣS|/Σ|S|: median %.2e  p99 %.2e  max %.2e   (> %.0e: %d CV, 面の無い CV: %d)" % (
        np.median(fin), np.percentile(fin, 99), fin.max(), a.tol, int((fin > a.tol).sum()), noface))
    for i in np.argsort(-np.where(ok_f, clo, -1))[:a.worst]:
        print("   worst CV %d  closure %.2e  vol %.3e  faces %d  at (%.5f, %.5f, %.5f)" % (i, clo[i], vol[i], nf[i], *cc[i]))
    ok = noface == 0 and not badv.any() and fin.max() <= a.tol
    print("VERDICT: %s" % ("PASS" if ok else "FAIL")); return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
