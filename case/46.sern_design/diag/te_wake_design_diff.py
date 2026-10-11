#!/usr/bin/env python3
"""修正後の格子での L_sw 1.0 − 0.8 の設計差と加重目的 ΔJ = 0.5Δ6 + 0.3Δ10 + 0.2Δ4 (C_T_with_shear) を V3_EVAL.txt から出す (plan tooling-sern-te-wake-grid §6 #7・#8)。
使い方: te_wake_design_diff.py <m10 0.8> <m10 1.0> <m6 0.8> <m6 1.0> <m4 0.8> <m4 1.0>  (各 V3_EVAL.txt)
"""
import re, sys
def read(p):
    out = {}
    for line in open(p):
        m = re.match(r"\s+(C_T_with_shear|C_T|C_L|C_M)\s+平均\s+(\S+)\s+a\s+(\S+)\s+前窓差\s+(\S+)", line)
        if m: out[m.group(1)] = tuple(float(m.group(i)) for i in (2, 3, 4))
    return out
p = sys.argv[1:]; ops = ("m10_on", "m6_on", "m4_off"); W = {"m6_on": 0.5, "m10_on": 0.3, "m4_off": 0.2}
dJ = 0.0; uJ = 0.0
for k, op in enumerate(ops):
    a, b = read(p[2 * k]), read(p[2 * k + 1])
    print(f"== {op} (L_sw 1.0 − 0.8)")
    for q in ("C_T", "C_T_with_shear", "C_L", "C_M"):
        d = b[q][0] - a[q][0]; u = max(a[q][1], a[q][2], 1e-6) + max(b[q][1], b[q][2], 1e-6)
        print(f"   {q:15s} 0.8 {a[q][0]:+.7f}  1.0 {b[q][0]:+.7f}  差 {d:+.4e}  時間変動幅 {u:.2e}")
    d6 = b["C_T_with_shear"][0] - a["C_T_with_shear"][0]; dJ += W[op] * d6
    uJ += W[op] * (max(a["C_T_with_shear"][1], a["C_T_with_shear"][2], 1e-6) + max(b["C_T_with_shear"][1], b["C_T_with_shear"][2], 1e-6))
print(f"ΔJ (C_T_with_shear、重み 0.5/0.3/0.2) = {dJ:+.4e}   うち時間変動幅の寄与 {uJ:.2e} (E_J の格子・領域感度の幅は含まない)")
