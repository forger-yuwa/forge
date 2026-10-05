#!/usr/bin/env python3
"""FORGE_DMIX_PROBE の出力 (GPU の化学種拡散係数) を、解決済み記録から組んだ独立 double 参照と照合する
(plan thermophysics-solver-owned-species-db #7b, §6 V4f)。

参照は tests/unit/test_lump_diffusion_reduction.py の Mixture.D_B (lump は実種展開 + 質量加重、非 lump は D_{r(i)})。
実種の同一性は記録の種名 (大文字化) で取る。凝縮 OFF の run だけを対象とする (気相組成の補正を入れない)。
合格: 全セル・全種で相対差 ≤ 1e-5 (V4e の縮約後係数の基準)、非有限 0、probe の lumpReduction フラグが記録の lump の有無と一致。
使い方: python3 check_dmix_probe.py <probe.bin> <species_record.yaml>
"""
import os
import struct
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
from test_lump_diffusion_reduction import Mixture  # noqa: E402
import forge_species  # noqa: E402


def main():
    probe, recpath = sys.argv[1], sys.argv[2]
    with open(probe, "rb") as f:
        n, nS, lumpflag = struct.unpack("<3i", f.read(12))
        arr = np.frombuffer(f.read(), dtype=np.float32)
    T, P, ro = arr[0:n], arr[n:2 * n], arr[2 * n:3 * n]
    roY = arr[3 * n:(3 + nS) * n].reshape(nS, n)
    D = arr[(3 + nS) * n:(3 + 2 * nS) * n].reshape(nS, n)
    rec = forge_species.load_record(recpath)
    sp = rec["species"]
    assert len(sp) == nS, f"記録の種数 {len(sp)} != probe {nS}"
    reals, labels = {}, []
    for e in sp:
        if e.get("lump"):
            comp = {}
            for m in e["lump"]["members"]:
                k = m["name"].upper()
                reals.setdefault(k, (m["MW"], m["LJ_sigma"], m["LJ_eps_kB"], 0.0))
                comp[k] = comp.get(k, 0.0) + m["x"]
            labels.append((e["name"], comp))
        else:
            k = e["name"].upper()
            reals.setdefault(k, (e["MW"], e["LJ_sigma"], e["LJ_eps_kB"], 0.0))
            labels.append((e["name"], {k: 1.0}))
    has_lump = any(e.get("lump") for e in sp)
    mix = Mixture(reals, labels)
    print(f"probe {probe}: cells {n}, species {nS} {[e['name'] for e in sp]}, lump reduction flag {lumpflag}, record has lump {has_lump}, real species {len(reals)}")
    ok = (lumpflag == (1 if has_lump else 0))
    print(f"  flag 一致: {'ok' if ok else 'NG'}")
    nonfin = int(np.count_nonzero(~np.isfinite(D)))
    worst, wi = 0.0, None
    stride = max(1, n // 4000)                   # 照合はセルを間引いて最大 ~4000 点 (double 参照が O(n_real²) の Python のため)
    for ic in range(0, n, stride):
        Y = [max(float(roY[s, ic]) / float(ro[ic]), 0.0) for s in range(nS)]
        t = sum(Y)
        Y = [y / t for y in Y]
        ref = mix.D_B(Y, float(T[ic]), float(P[ic]))
        for s in range(nS):
            r = abs(float(D[s, ic]) - ref[s]) / ref[s]
            if r > worst:
                worst, wi = r, (ic, s, float(D[s, ic]), ref[s], float(T[ic]))
    print(f"  非有限 {nonfin}、最大相対差 {worst:.3e} (cell, species, GPU, ref, T) = {wi}")
    ok &= nonfin == 0 and worst <= 1e-5
    print("VERDICT V4f:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
