#!/usr/bin/env python3
"""plan §6.3 の追加記録 (追加実行なし): SERN g3 の twall_z の最悪差の位置と ULP、condClampCorrQ_0 が最大になった run の位置の値。

    python3 record_extra.py twall  --sern-root ~/forge-r8/case/46.sern_design
    python3 record_extra.py clamp  --root ~/forge-b4/case/66.hostmem_regression

出力は標準出力 (notes.txt に写す)。読むだけで何も書かない。
"""
import argparse
import itertools
import os

import h5py
import numpy as np

SERN = {"B1": "run_1072_hm_g3_base_r1", "B2": "run_1073_hm_g3_base_r2", "B3": "run_1074_hm_g3_base_r3",
        "N1": "run_1075_hm_g3_new_r1", "N2": "run_1076_hm_g3_new_r2", "N3": "run_1077_hm_g3_new_r3"}


def ukey(x):
    i = np.asarray(x, dtype=np.float32).view(np.int32).astype(np.int64)
    return np.where(i >= 0, i, -(i & 0x7FFFFFFF))


def ulp(a, b):
    return np.abs(ukey(a) - ukey(b))


def g(x):
    return "%.9g" % x


def twall(a):
    fn = "res_vehicle_base_18_100.h5"
    R = {k: os.path.join(os.path.expanduser(a.sern_root), v) for k, v in SERN.items()}
    V = {k: h5py.File(os.path.join(v, fn), "r")["VALUE/twall_z"][()] for k, v in R.items()}
    C = h5py.File(os.path.join(R["B1"], fn), "r")["MESH/COORD"][()].reshape(-1, 3)
    b, n = ["B1", "B2", "B3"], ["N1", "N2", "N3"]
    cross = [(x, y) for x in b for y in n]
    within = list(itertools.combinations(b, 2)) + list(itertools.combinations(n, 2))

    def dm(p):
        return float(np.max(np.abs(V[p[0]].astype(np.float64) - V[p[1]].astype(np.float64))))
    wd, ws = max(cross, key=dm), max(within, key=dm)
    print(f"SERN g3 {fn} VALUE/twall_z ({len(V['B1'])} 点、MESH/COORD {len(C)} 点 = 値は境界節点ごと)")
    print(f"  D_abs = {dm(wd):.6g} (組 {wd[0]}-{wd[1]})、S_abs = {dm(ws):.6g} (組 {ws[0]}-{ws[1]})")
    print(f"  値の範囲 (B1) {V['B1'].min():.6g} 〜 {V['B1'].max():.6g}、max|x| {np.max(np.abs(V['B1'])):.6g}")
    for name, p in (("D (ビルド間) の最悪位置", wd), ("S (同ビルド内) の最悪位置", ws)):
        d = np.abs(V[p[0]].astype(np.float64) - V[p[1]].astype(np.float64))
        i = int(np.argmax(d))
        print(f"  {name}: index {i}、座標 {C[i].tolist()}")
        print("    6 本の値: " + ", ".join(f"{k} {g(V[k][i])}" for k in SERN))
        print(f"    その値の float32 刻み (ulp) {float(np.spacing(np.float32(abs(V['B1'][i])))):.6g}")
        pairs = list(itertools.combinations(b, 2)) + list(itertools.combinations(n, 2)) + cross
        print("    ULP 距離: " + ", ".join(f"{x}-{y} {int(ulp(V[x][i], V[y][i]))}" for x, y in pairs))
        for q in ("twall_x", "twall_y", "utau", "ypls", "Ps", "Ts", "qwall"):
            vals = {k: h5py.File(os.path.join(R[k], fn), "r")["VALUE/" + q][()][i] for k in SERN}
            print(f"    {q:8s} " + ", ".join(f"{k} {g(v)}" for k, v in vals.items()))
    ub = np.max(np.stack([ulp(V[x], V[y]) for x, y in itertools.combinations(b, 2)]), axis=0)
    un = np.max(np.stack([ulp(V[x], V[y]) for x, y in itertools.combinations(n, 2)]), axis=0)
    ux = np.max(np.stack([ulp(V[x], V[y]) for x, y in cross]), axis=0)
    print(f"  全点の ULP (各点で組の最大): base 内 中央値 {int(np.median(ub))}・最大 {int(ub.max())}、new 内 中央値 {int(np.median(un))}・最大 {int(un.max())}、"
          f"ビルド間 中央値 {int(np.median(ux))}・最大 {int(ux.max())}")
    w = np.maximum(ub, un)
    print(f"  ビルド間の ULP が同ビルド内の ULP (両ビルドの大きい方) を超える点: {int(np.sum(ux > w))} / {len(ux)}"
          f" (2 倍を超える点 {int(np.sum(ux > 2 * w))})")


def clamp(a):
    root = os.path.expanduser(a.root)
    cases = [("c44dual_ckpt100", "run_0106_c44dual_ckpt100_new_r1",
              ["run_0029_c44dual_ckpt100_base_r1", "run_0038_c44dual_ckpt100_base_r2", "run_0068_c44dual_ckpt100_base_r3",
               "run_0106_c44dual_ckpt100_new_r1", "run_0136_c44dual_ckpt100_new_r2", "run_0166_c44dual_ckpt100_new_r3"]),
             ("c44dual_restart100", "run_0034_c44dual_restart100_base_r1",
              ["run_0034_c44dual_restart100_base_r1", "run_0039_c44dual_restart100_base_r2", "run_0069_c44dual_restart100_base_r3",
               "run_0107_c44dual_restart100_new_r1", "run_0137_c44dual_restart100_new_r2", "run_0167_c44dual_restart100_new_r3"])]
    names = ["condClampCorrQ_0", "condClampCorr_0", "condLim_0", "condR30_0", "rog_0", "roQ0_0", "roQ1_0", "roQ2_0",
             "g_0", "Q0_0", "Q1_0", "Q2_0", "ro", "T", "P", "condS_0", "condTheta_0", "condDrdt_0", "passiveFloorCorr_Q0_0",
             "passiveFloorCorr_Q1_0", "passiveFloorCorr_Q2_0", "passiveLimCorr_Q0_0", "passiveFloorCorr_g_0"]
    for cfg, top, runs in cases:
        f = h5py.File(os.path.join(root, top, "res_100.h5"), "r")
        c = f["VALUE/condClampCorrQ_0"][()]
        i = int(np.argmax(np.abs(c)))
        C = f["MESH/COORD"][()].reshape(-1, 3)
        nz = int(np.count_nonzero(c))
        print(f"{cfg}: condClampCorrQ_0 が最大の run {top}: index {i}、座標 {C[i].tolist()}、値 {g(c[i])} (非零 {nz} 節点)")
        order = np.argsort(-np.abs(c))[:5]
        print("  上位 5 節点: " + ", ".join(f"{int(j)}:{g(c[j])}" for j in order))
        print(f"  {'量':24s} " + " ".join(f"{r.split('_')[1] + '_' + r.split('_')[-1]:>16s}" for r in runs))
        for q in names:
            vals = []
            for r in runs:
                h = h5py.File(os.path.join(root, r, "res_100.h5"), "r")
                vals.append(g(h["VALUE/" + q][()][i]) if "VALUE/" + q in h else "-")
            print(f"  {q:24s} " + " ".join(f"{v:>16s}" for v in vals))
        print("  (出力は補正後の状態だけ。condClampCorrQ_0 は max_k |ΔQ_k| / max(|Q_k,before|, 1e-30) の最大 [condensationRealizability_d.cuh:313] で、"
              "補正前の Q_before は出力に無い)")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("twall")
    s.add_argument("--sern-root", required=True)
    s = sub.add_parser("clamp")
    s.add_argument("--root", required=True)
    a = ap.parse_args()
    {"twall": twall, "clamp": clamp}[a.cmd](a)


if __name__ == "__main__":
    main()
