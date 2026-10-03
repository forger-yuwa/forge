#!/usr/bin/env python3
"""V2d の forge run と独立参照 1D (`ref1d_euler_tp.py`、Δx 0.3125 mm、参照自身は Δx/4 との差が許容の 1/5 以内) の比較 (記録用、合否に入れない;
plan boundary-node-farfield-characteristic §5.1 #3e)。評価点の時系列で、forge の時刻に参照を線形補間して差を取る。
  接触波 (V2d-1): max_t |T_forge − T_ref| と max_t |P_forge − P_ref| (forge の short / long)
  音響 (V2d-2 隔離配置): max_t |P_forge − P_ref| / 入射振幅 (参照の max|P − P∞|)
  python3 compare_ref1d.py REF_DIR
"""
import sys
import numpy as np

REF = sys.argv[1]
P0 = 2851.0


def probe(run):
    d = np.genfromtxt(f"{run}/point_probe_0.out", delimiter=",", names=True)
    return np.asarray(d["TotalTime"], float), np.asarray(d["P"], float), np.asarray(d["T"], float)


def ref(case, gas):
    d = np.genfromtxt(f"{REF}/ref_{case}_{gas}_dx0p3125.csv", delimiter=",", names=True)
    return d["t"], d["P"], d["T"]


for gas, runs in (("cpg", ("run_0070_v2d_cpg_contact_short", "run_0071_v2d_cpg_contact_long")),
                  ("tp1", ("run_0075_v2d_tp1_contact_short", "run_0076_v2d_tp1_contact_long")),
                  ("tp2", ("run_0080_v2d_tp2_contact_short", "run_0081_v2d_tp2_contact_long"))):
    tr, Pr, Tr = ref("contact", gas)
    for r in runs:
        t, P, T = probe(r); k = t <= tr[-1]
        print(f"contact {gas} {r}: max|T − T_ref| {np.max(np.abs(T[k] - np.interp(t[k], tr, Tr))):.3f} K (ピーク T_ref {Tr.max():.2f}、forge {T[k].max():.2f})、"
              f"max|P − P_ref| {np.max(np.abs(P[k] - np.interp(t[k], tr, Pr))):.4g} Pa")
for gas, runs in (("cpg", ("run_0132_v2d2_cpg_short_dt2", "run_0133_v2d2_cpg_long_dt2")),
                  ("tp1", ("run_0136_v2d2_tp1_short_dt2", "run_0137_v2d2_tp1_long_dt2")),
                  ("tp2", ("run_0140_v2d2_tp2_short_dt2", "run_0141_v2d2_tp2_long_dt2"))):
    tr, Pr, Tr = ref("acoustic", gas); A = np.max(np.abs(Pr - P0))
    for r in runs:
        t, P, T = probe(r); k = t <= tr[-1]
        d = np.abs(P[k] - np.interp(t[k], tr, Pr))
        print(f"acoustic {gas} {r}: max|P − P_ref| {d.max():.4g} Pa = {d.max() / A:.2%} of 参照の入射振幅 {A:.4g} Pa (t {t[k][d.argmax()]:.4e})、forge の振幅 {np.max(np.abs(P[k] - P0)):.4g} Pa")
