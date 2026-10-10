#!/usr/bin/env python3
"""後縁下流の局所変形の 0 step 比較: A (te_wake_blend_H 0) と B (1.0) の最終 node 格子で、対応するヘキサごとの
最小スケール済みヤコビアンと AR を比べる (plan tooling-sern-te-wake-grid §5.1 #4、codex diagnose 2026-10-11 g4-admission-jacobian)。
読むだけ。変形領域 = 動いた節点を 1 個以上含むヘキサ。
使い方: te_wake_hex_compare.py <A の run> <B の run> [--thr 0.65] [--json out.json]
"""
import argparse, json, sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import te_wake_grid_check as T

ap = argparse.ArgumentParser(); ap.add_argument("A"); ap.add_argument("B"); ap.add_argument("--thr", type=float, default=0.65); ap.add_argument("--json")
a = ap.parse_args()
GA, GB = T.load(a.A), T.load(a.B)
cA, cB = np.asarray(GA["coords"], float), np.asarray(GB["coords"], float)
hA, hB = np.asarray(GA["hexes"]), np.asarray(GB["hexes"])
H = GA["H"]
out = {"A": a.A, "B": a.B, "thr": a.thr}
out["same_hexes"] = bool(hA.shape == hB.shape and np.array_equal(hA, hB))
out["same_n_nodes"] = bool(cA.shape == cB.shape)
if not (out["same_hexes"] and out["same_n_nodes"]):
    print("UNDECIDABLE: 接続または節点数が違う"); print(json.dumps(out)); sys.exit(2)
moved = np.any(cA != cB, axis=1)
reg = moved[hA].any(axis=1)
def q(c, h):
    J, Js = T.corner_jacobians(c, h)
    sgn = np.sign(np.median(J))
    s = (Js * sgn).min(axis=1)
    X = c[h]
    edges = [(0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4),(0,4),(1,5),(2,6),(3,7)]
    L = np.stack([np.linalg.norm(X[:, i] - X[:, j], axis=1) for i, j in edges], axis=1)
    ar = L.max(axis=1) / np.maximum(L.min(axis=1), 1e-300)
    return s, ar
sA, arA = q(cA, hA); sB, arB = q(cB, hB)
thr = a.thr
new_low = (sA >= thr) & (sB < thr)
worse_low = (sA < thr) & (sB < sA)
outside_diff = (~reg) & ((sA != sB) | (arA != arB))
cen = cB[hB].mean(axis=1) / H
def loc(mask, key, n=5):
    idx = np.flatnonzero(mask)
    if idx.size == 0: return []
    idx = idx[np.argsort(key[idx])][:n]
    return [{"hex": int(k), "A": float(sA[k]), "B": float(sB[k]), "centroid_over_H": [round(float(v), 5) for v in cen[k]]} for k in idx]
out.update({
    "n_hex": int(hA.shape[0]), "n_moved_nodes": int(moved.sum()), "n_region_hex": int(reg.sum()),
    "min_A": float(sA.min()), "min_B": float(sB.min()),
    "n_lt_thr_A": int((sA < thr).sum()), "n_lt_thr_B": int((sB < thr).sum()),
    "region_n_lt_thr_A": int((reg & (sA < thr)).sum()), "region_n_lt_thr_B": int((reg & (sB < thr)).sum()),
    "region_min_A": float(sA[reg].min()) if reg.any() else None, "region_min_B": float(sB[reg].min()) if reg.any() else None,
    "new_lt_thr": int(new_low.sum()), "worsened_existing_lt_thr": int(worse_low.sum()),
    "outside_region_changed_hex": int(outside_diff.sum()),
    "ar_gt5000_A": int((arA > 5000).sum()), "ar_gt5000_B": int((arB > 5000).sum()), "ar_max_A": float(arA.max()), "ar_max_B": float(arB.max()),
    "examples_new_lt_thr": loc(new_low, sB), "examples_worsened": loc(worse_low, sB - sA),
})
# AR > 5000 の位置 (A)
idx = np.flatnonzero(arA > 5000)
if idx.size:
    cx, cy = cen[idx, 0], cen[idx, 1]
    out["ar_gt5000_A_xH_range"] = [float(cx.min()), float(cx.max())]; out["ar_gt5000_A_yH_range"] = [float(cy.min()), float(cy.max())]
    out["ar_gt5000_A_in_region"] = int(reg[idx].sum())
    hist, edges = np.histogram(cx, bins=[-1, 0, 0.5, 1.2, 1.5, 2.2, 5, 12])
    out["ar_gt5000_A_xH_hist"] = list(zip([float(e) for e in edges[:-1]], hist.tolist()))
for k, v in out.items():
    if not isinstance(v, list): print(f"{k}: {v}")
print("examples_new_lt_thr:", out["examples_new_lt_thr"][:3]); print("examples_worsened:", out["examples_worsened"][:3])
print("ar_gt5000_A_xH_hist:", out.get("ar_gt5000_A_xH_hist"))
v = "LOCAL_DEGRADATION" if (out["new_lt_thr"] or out["worsened_existing_lt_thr"]) else "NO_LOCAL_DEGRADATION"
print("HEX_COMPARE VERDICT:", v)
if a.json: json.dump(out, open(a.json, "w"), indent=1)
