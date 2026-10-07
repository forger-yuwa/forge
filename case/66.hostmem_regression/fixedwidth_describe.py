#!/usr/bin/env python3
"""固定幅の独立 A/B (plan §6.3) 段階 2 の記述 (判定には使わない)。AWS で回す。読むだけ。

    python3 fixedwidth_describe.py --plan fixedwidth_c44dual_ckpt100/plan.json

1. 各 run (新規 12 本と既存 base 3 本) の condClampCorrQ_0 の最大・位置 (節点 index・座標)、> 1e22 の状態に入ったか
   (ビルドごとの本数)、そのとき condClampCorrQ_0 × 1e-30 がその節点の roQ0_0/roQ1_0/roQ2_0 のどれかと一致するか
   (分母が床 1e-30 だったことと整合するか)。
2. fixedwidth_eval.py eval が T 超過とした量 (result/fixedwidth_result.txt から読む) について、run・値・T 比・最悪の既存 base・
   最大差の位置 (index・座標) とその位置の X と B_i の値。
"""
import argparse
import json
import os
import re
import sys

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fixedwidth_eval as fw  # noqa: E402

FLOOR = 1e-30


def g(x):
    return "%.6g" % x


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--plan", required=True)
    ap.add_argument("--threshold", type=float, default=1e22)
    a = ap.parse_args()
    plan = json.load(open(a.plan))
    pdir = os.path.dirname(os.path.abspath(a.plan))
    runs = [(j["run"], j["group"]) for j in plan["runs"]]
    base = plan["base_existing"]
    out = []
    out.append(f"1. condClampCorrQ_0 の最大と位置 (> {a.threshold:g} を「1.46e23 型の状態」として数える。記述であって判定には使わない)")
    out.append(f"   {'run':30s} {'build':5s} {'最大':>12s} {'index':>6s}  {'座標 (x, y)':22s} {'>1e22':5s} 床との整合 (max×1e-30 と roQk_0 の相対差が最小の k)")
    count = {"b": [0, 0], "n": [0, 0], "existing": [0, 0]}
    for r, grp in runs + [(b, "existing") for b in base]:
        f = h5py.File(os.path.join(HERE, r, "res_100.h5"), "r")
        c = f["VALUE/condClampCorrQ_0"][()]
        i = int(np.argmax(np.abs(c)))
        C = f["MESH/COORD"][()].reshape(-1, 3)
        big = abs(float(c[i])) > a.threshold
        count[grp][0] += 1
        count[grp][1] += 1 if big else 0
        best = None
        for k in range(3):
            q = float(f[f"VALUE/roQ{k}_0"][()][i])
            rel = abs(q - float(c[i]) * FLOOR) / max(abs(q), 1e-300)
            if best is None or rel < best[1]:
                best = (k, rel, q)
        cons = f"roQ{best[0]}_0 {g(best[2])} (相対差 {best[1]:.1e})" + (" → 整合" if best[1] < 1e-5 else "")
        out.append(f"   {r:30s} {grp:5s} {g(c[i]):>12s} {i:6d}  ({C[i][0]:.5f}, {C[i][1]:.5f})  {'yes' if big else 'no':5s} {cons}")
    out.append("   > 1e22 の本数: 新規 base " + f"{count['b'][1]}/{count['b'][0]}、新規 new {count['n'][1]}/{count['n'][0]}、"
               f"既存 base (T の算定に使った 3 本) {count['existing'][1]}/{count['existing'][0]}")
    # 2. T 超過の一覧
    res = os.path.join(pdir, "result", "fixedwidth_result.txt")
    trows = {(r["file"], r["name"]): r for r in fw.read_T(os.path.join(pdir, "T_frozen.tsv"))}
    out.append("")
    out.append("2. T を超えた量 (fixedwidth_eval.py の結果) の位置")
    out.append(f"   {'run':24s} {'量':32s} {'v':>11s} {'T':>11s} {'v/T':>8s} 最悪  index  座標 (x, y)              X の値        B_i の値")
    cur = None
    for line in open(res):
        m = re.match(r"(run_\S+) \((b|n)\):", line)
        if m:
            cur = m.group(1)
            continue
        m = re.match(r"\s+(\S+?):(\S+) \[(.+?)\] v (\S+) T (\S+) — (.+)$", line.rstrip("\n"))
        if not (m and cur):
            continue
        fn, name, mode, v, t, why = m.groups()
        wb = re.search(r"最悪 B(\d)", why)
        bi = int(wb.group(1)) - 1 if wb else 0
        x = np.asarray(fw.load_quantity(os.path.join(HERE, cur), fn, name), dtype=np.float64)
        b = np.asarray(fw.load_quantity(os.path.join(HERE, base[bi]), fn, name), dtype=np.float64)
        d = np.abs(x - b)
        i = int(np.argmax(d))
        C = h5py.File(os.path.join(HERE, cur, fn), "r")["MESH/COORD"][()].reshape(-1, 3) if fn.endswith(".h5") else None
        xy = f"({C[i][0]:.5f}, {C[i][1]:.5f})" if C is not None and i < len(C) else "-"
        vv = float(d[i])                       # 結果ファイルの丸めた値でなく、配列から求め直した値
        tt = trows[(fn, name)]["T"]
        vt = vv / tt if tt > 0 else float("inf")
        out.append(f"   {cur:24s} {name[:32]:32s} {vv:11.4e} {tt:11.4e} {vt:8.4g} B{bi + 1}  {i:5d}  {xy:24s} {g(x[i]):>13s} {g(b[i]):>13s}")
    print("\n".join(out))


if __name__ == "__main__":
    main()
