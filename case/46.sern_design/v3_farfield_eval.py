#!/usr/bin/env python3
"""farfield plan §6 V3 の比較 (R4d と同じ規則)。force_history.csv (collect() が作る) から、
末尾 10000 step の平均・a = max|値 − 平均|・前窓差 (その前 10000 step の平均との差、≤ 0.1ε が窓条件)、
2 run の D = |Δ平均| + a₁ + a₂ を出す。ε: C_T・C_T_with_shear・C_L 5e-4、C_M 5e-3。
  python3 v3_farfield_eval.py RUN...                 各 run の統計
  python3 v3_farfield_eval.py --pair RUN_A RUN_B      D (判定)。両 run の窓条件 (全 4 量)・末尾窓の標本数 ≥ 20・有限性を
                                                     必要条件にし、満たさなければ D を判定せず「判定不能」(codex 2026-09-29 M)
  --tau X (--pair と併用): 全量の閾値を X × ε に置き換える (#4a の τ = 0.2)
"""
import csv, sys
import numpy as np

Q = ("C_T", "C_T_with_shear", "C_L", "C_M")
EPS = {"C_T": 5e-4, "C_T_with_shear": 5e-4, "C_L": 5e-4, "C_M": 5e-3}


def stats(run):
    rows = list(csv.DictReader(open(run + "/force_history.csv")))
    st = np.array([int(x["step"]) for x in rows])
    t = st > st.max() - 10000; p = (st > st.max() - 20000) & ~t
    out = {}
    for q in Q:
        v = np.array([float(x[q]) for x in rows])
        m = v[t].mean()
        out[q] = (m, float(np.max(np.abs(v[t] - m))), abs(m - v[p].mean()) if p.any() else float("nan"))
    return out, int(st.max()), int(t.sum())


def prereq(run):
    s, mx, n = stats(run)
    bad = []
    if n < 20:
        bad.append(f"末尾窓の標本 {n} < 20")
    for q in Q:
        m, a, pw = s[q]
        if not all(np.isfinite([m, a, pw])):
            bad.append(f"{q} 非有限")
        elif pw > 0.1 * EPS[q]:
            bad.append(f"{q} 窓条件 {pw:.2e} > {0.1 * EPS[q]:.1e}")
    return bad


if __name__ == "__main__":
    if sys.argv[1] == "--pair":
        a, b = sys.argv[2], sys.argv[3]
        tau = float(sys.argv[sys.argv.index("--tau") + 1]) if "--tau" in sys.argv else 1.0
        bad = [f"{r}: {x}" for r in (a, b) for x in prereq(r)]
        if bad:
            print(f"{a} vs {b}: 判定不能 (必要条件未達: {'; '.join(bad)})")
            sys.exit(2)
        sa, _, _ = stats(a); sb, _, _ = stats(b)
        ok = True
        parts = []
        for q in Q:
            d = abs(sa[q][0] - sb[q][0]) + sa[q][1] + sb[q][1]
            lim = tau * EPS[q]
            ok &= d <= lim
            parts.append(f"{q} D {d:.2e} Δ平均 {abs(sa[q][0] - sb[q][0]):.2e} ({'≤' if d <= lim else '>'}{tau:g}ε)")
        print(f"{a} vs {b}: " + "  ".join(parts) + f"  → {'全量 D ≤ ' + f'{tau:g}ε' if ok else 'D > ' + f'{tau:g}ε あり'}")
    else:
        for r in sys.argv[1:]:
            s, mx, n = stats(r)
            print(f"{r}: 最終 step {mx}、末尾窓 {n} 標本")
            for q in Q:
                m, a, pw = s[q]
                print(f"   {q:15s} 平均 {m:.7f}  a {a:.2e}  前窓差 {pw:.2e} (≤ 0.1ε: {'OK' if pw <= 0.1 * EPS[q] else 'NG'})")
