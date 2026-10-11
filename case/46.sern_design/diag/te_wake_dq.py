#!/usr/bin/env python3
"""plan tooling-sern-te-wake-grid §6 #6: 処置の前後 (A = 0 格子 / B = 局所化) の D_q を v3_farfield_eval の出力 (V3_EVAL.txt) から計算する。
D_q = |平均差| + max(a_A, d_A, 1e-6) + max(a_B, d_B, 1e-6)、許容 C_T・C_T_with_shear・C_L ≤ 5e-5、C_M ≤ 5e-4。
使い方: te_wake_dq.py <A の V3_EVAL.txt> <B の V3_EVAL.txt>
"""
import re, sys
TOL = {"C_T": 5e-5, "C_T_with_shear": 5e-5, "C_L": 5e-5, "C_M": 5e-4}
def read(p):
    out = {}
    for line in open(p):
        m = re.match(r"\s+(C_T_with_shear|C_T|C_L|C_M)\s+平均\s+(\S+)\s+a\s+(\S+)\s+前窓差\s+(\S+)", line)
        if m: out[m.group(1)] = tuple(float(m.group(i)) for i in (2, 3, 4))
    return out
A, B = read(sys.argv[1]), read(sys.argv[2])
ok = True
for q in ("C_T", "C_T_with_shear", "C_L", "C_M"):
    (mA, aA, dA), (mB, aB, dB) = A[q], B[q]
    dq = abs(mB - mA) + max(aA, dA, 1e-6) + max(aB, dB, 1e-6)
    p = dq <= TOL[q]; ok &= p
    print(f"{q:15s} A {mA:+.7f}  B {mB:+.7f}  |diff| {abs(mB-mA):.3e}  D_q {dq:.3e}  tol {TOL[q]:.0e}  {'PASS' if p else 'FAIL'}")
print("D_Q VERDICT:", "CARRY_OVER" if ok else "RETAKE")
