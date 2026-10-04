#!/usr/bin/env python3
"""R3 (plan boundary-node-periodic-gradient-fix §6) の判定: 非周期・軸対称で勾配・リミタ配列が旧新で変わらないこと。

§6 R3「判定の具体化 (2026-09-26、測定前)」をそのまま実装する:
  同一入力状態から 1 step、output.level 2 の res_1.h5 の全勾配 (d*d[xyz])・リミタ (limiter_*) 配列を旧 2 回・新 2 回で比較。
  - 旧同士がビット一致する配列 → 旧新 (old_a vs new_a) もビット一致を要求。
  - 旧同士でも一致しない配列 → 旧新の最大差 ≤ 旧同士の最大差 × 2 かつ 不一致数の桁 (floor(log10)) が同じ。
ケースの VERDICT は res_1 の全配列が PASS なら PASS。res_0 (step 前の勾配) も同じ規則で表に出すが VERDICT には入れない。
保存量 (ro, roU*, roe, roK, roOmega, roY*) の res_1 比較は参考 (判定なし)。NaN/Inf は全配列・全 run で数える。

入力は r3_prepare.py が作る <scratch>/r3_<case>_<old|new>_<a|b>/。

使い方: python3 r3_compare.py <scratch_dir> [--out R3.txt]
"""
import argparse
import math
import os
import re
import sys

import h5py
import numpy as np

CASES = ("case48", "case16", "axi")
RUNS = ("old_a", "old_b", "new_a", "new_b")
PAIRS = (("old_a", "old_b"), ("new_a", "new_b"), ("old_a", "new_a"))
GRAD_RX = re.compile(r"^d.+d[xyz]$")
STATE = ("ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega", "roY0", "roY1")


def read_values(path):
    with h5py.File(path, "r") as f:
        return {k: np.array(f["VALUE"][k]) for k in f["VALUE"].keys()}


def cmp(a, b):
    """(不一致数 = ビット比較, 最大絶対差 (NaN を除く))。"""
    if a.shape != b.shape:
        return None, None
    ai = a.view(np.uint32) if a.dtype == np.float32 else a.view(np.uint64)
    bi = b.view(np.uint32) if b.dtype == np.float32 else b.view(np.uint64)
    n = int(np.count_nonzero(ai != bi))
    d = np.abs(a.astype(np.float64) - b.astype(np.float64))
    d = d[np.isfinite(d)]
    return n, (float(d.max()) if d.size else 0.0)


def digits(n):
    return int(math.floor(math.log10(n))) if n > 0 else -1


def judge(oo, on):
    """§6 R3 の規則。oo = (n, max) 旧同士、on = 旧新。"""
    n_oo, m_oo = oo
    n_on, m_on = on
    if n_oo is None or n_on is None:
        return "判定不能 (形状不一致)"
    if n_oo == 0:
        return "PASS" if n_on == 0 else "FAIL (旧同士ビット一致・旧新不一致)"
    if n_on == 0:
        return "PASS (旧同士は不一致、旧新はビット一致)"
    ok_max = m_on <= 2.0 * m_oo
    ok_dig = digits(n_on) == digits(n_oo)
    if ok_max and ok_dig:
        return "PASS (ノイズ床)"
    why = []
    if not ok_max:
        why.append(f"最大差 {m_on:.3e} > 2×{m_oo:.3e}")
    if not ok_dig:
        why.append(f"不一致数の桁 {n_on} vs {n_oo}")
    return "FAIL (" + "; ".join(why) + ")"


def fmt(nm):
    n, m = nm
    return "shape!" if n is None else f"{n:>7d} {m:10.3e}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("scratch")
    ap.add_argument("--out", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "R3.txt"))
    args = ap.parse_args()
    lines = []
    P = lines.append
    P("R3: 非周期・軸対称の勾配・リミタ配列 旧新比較 (plan boundary-node-periodic-gradient-fix §6 R3、2026-09-26)")
    P(f"run root: {os.path.abspath(args.scratch)}/r3_<case>_<old|new>_<a|b>")
    for c in CASES:
        prov = os.path.join(args.scratch, f"r3_{c}_new_a", "RUN_PROVENANCE.txt")
        provo = os.path.join(args.scratch, f"r3_{c}_old_a", "RUN_PROVENANCE.txt")
        if os.path.exists(prov) and os.path.exists(provo):
            sha_n = [l for l in open(prov) if l.startswith("forge_sha256")][0].split(":")[1].strip()
            sha_o = [l for l in open(provo) if l.startswith("forge_sha256")][0].split(":")[1].strip()
            P(f"forge sha256: old {sha_o[:12]}…  new {sha_n[:12]}…")
            break
    P("列: 不一致数 (ビット比較) / 最大絶対差。判定は oo=old_a/old_b と on=old_a/new_a に §6 R3 の規則を当てる。")
    verdicts = {}
    for c in CASES:
        P("")
        P("=" * 118)
        P(f"case {c}")
        vals = {}
        for step in (0, 1):
            vals[step] = {r: read_values(os.path.join(args.scratch, f"r3_{c}_{r}", f"res_{step}.h5")) for r in RUNS}
        # NaN / Inf
        bad = []
        for step in (0, 1):
            for r in RUNS:
                for k, a in vals[step][r].items():
                    if a.dtype.kind == "f":
                        nb = int(np.count_nonzero(~np.isfinite(a)))
                        if nb:
                            bad.append(f"res_{step} {r} {k}: 非有限 {nb}")
        P("NaN/Inf: " + ("なし (全 run・res_0/res_1 の全 VALUE)" if not bad else ""))
        for b in bad:
            P("  " + b)
        ok_all = True
        for step in (1, 0):
            names = sorted(k for k in vals[step]["old_a"] if GRAD_RX.match(k) or k.startswith("limiter_"))
            P("")
            P(f"-- res_{step}.h5 勾配・リミタ ({len(names)} 配列)" + ("  [VERDICT 対象]" if step == 1 else "  [参考: step 前]"))
            P(f"{'array':<10} {'old_a/old_b':>18} {'new_a/new_b':>18} {'old_a/new_a':>18}  判定")
            for k in names:
                if any(k not in vals[step][r] for r in RUNS):
                    P(f"{k:<10} 一部の run に無い"); ok_all = ok_all and step != 1; continue
                res = {p: cmp(vals[step][p[0]][k], vals[step][p[1]][k]) for p in PAIRS}
                j = judge(res[PAIRS[0]], res[PAIRS[2]])
                if step == 1 and not j.startswith("PASS"):
                    ok_all = False
                amax = float(np.nanmax(np.abs(vals[step]["old_a"][k].astype(np.float64))))
                P(f"{k:<10} {fmt(res[PAIRS[0]])} {fmt(res[PAIRS[1]])} {fmt(res[PAIRS[2]])}  {j}   (max|値| {amax:.3e})")
        P("")
        P("-- res_1.h5 保存量 (参考、判定なし)")
        P(f"{'array':<10} {'old_a/old_b':>18} {'new_a/new_b':>18} {'old_a/new_a':>18}")
        for k in STATE:
            if all(k in vals[1][r] for r in RUNS):
                res = {p: cmp(vals[1][p[0]][k], vals[1][p[1]][k]) for p in PAIRS}
                P(f"{k:<10} {fmt(res[PAIRS[0]])} {fmt(res[PAIRS[1]])} {fmt(res[PAIRS[2]])}")
        v = "PASS" if (ok_all and not bad) else ("FAIL" if not bad else "FAIL (非有限値あり)")
        verdicts[c] = v
        P("")
        P(f"VERDICT {c}: {v}")
    P("")
    P("=" * 118)
    P("SUMMARY: " + "  ".join(f"{c}={v}" for c, v in verdicts.items()))
    P("VERDICT: " + ("PASS" if all(v == "PASS" for v in verdicts.values()) else "FAIL"))
    txt = "\n".join(lines) + "\n"
    with open(args.out, "w") as fp:
        fp.write(txt)
    sys.stdout.write(txt)


if __name__ == "__main__":
    main()
