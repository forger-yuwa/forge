#!/usr/bin/env python3
"""2 つの res_*.h5 の VALUE を全キーで比較: max|Δ| と max|Δ|/max|ref|。回帰 (ビット同等/反復ノイズ) の判定用。
usage: diff_res.py REF.h5 NEW.h5 [--keys k1,k2,...] [--tol REL] [--tolfile noise.json] [--factor 2]
  --keys 省略時は保存量+原始量+凝縮量の既定集合 (診断量 limiter_*/res_*/wf_* 等は除外)。欠落キー・形状不一致・非有限値は FAIL。
  --tol: 全キー共通の相対許容値。--tolfile: 変数別許容値 JSON ({key: rel}) を --factor 倍して使う (反復ノイズ床から作る)。
  --dump out.json で本比較の rel を書き出す (ノイズ床の記録用)。不合格は非ゼロ終了。"""
import sys, json, argparse, h5py, numpy as np
DEFAULT_KEYS = ["ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega", "roY0", "roY1", "rog_0", "roQ0_0", "roQ1_0", "roQ2_0",
                "P", "T", "Ux", "Uy", "Uz", "k", "omega", "Y0", "Y1", "g_0", "Q0_0", "Q1_0", "Q2_0", "sonic", "gamma", "cp", "Rmix", "vis_turb", "h0"]
ap = argparse.ArgumentParser(); ap.add_argument("ref"); ap.add_argument("new"); ap.add_argument("--keys", default=None)
ap.add_argument("--tol", type=float, default=None); ap.add_argument("--tolfile", default=None); ap.add_argument("--factor", type=float, default=2.0)
ap.add_argument("--dump", default=None)
a = ap.parse_args(); A = h5py.File(a.ref, "r")["VALUE"]; B = h5py.File(a.new, "r")["VALUE"]
keys = a.keys.split(",") if a.keys else [k for k in DEFAULT_KEYS if k in A]
tolmap = json.load(open(a.tolfile)) if a.tolfile else {}
worst = 0.0; rows = []; fail = False; rels = {}
for k in keys:
    if k not in A or k not in B: rows.append((k, "MISSING", float("inf"), 0, "FAIL")); fail = True; continue
    x = np.array(A[k], dtype=np.float64); y = np.array(B[k], dtype=np.float64)
    if x.shape != y.shape: rows.append((k, "SHAPE MISMATCH", float("inf"), 0, "FAIL")); fail = True; continue
    if not (np.all(np.isfinite(x)) and np.all(np.isfinite(y))): rows.append((k, "NON-FINITE", float("inf"), 0, "FAIL")); fail = True; continue
    d = np.abs(x - y); den = np.max(np.abs(x)); rel = (d.max()/den) if den > 0 else (0.0 if d.max() == 0 else np.inf)
    rels[k] = rel; nz = int((d > 0).sum()); worst = max(worst, rel)
    if a.tolfile and k not in tolmap:   # ノイズ床の記録に無いキー (例: 旧バイナリが出力しなかった gamma) は判定せず記録のみ
        lim = None; verdict = "n/a (no noise ref)"
    else:
        lim = tolmap[k]*a.factor if a.tolfile else a.tol
        verdict = "" if lim is None else ("ok" if rel <= lim else "FAIL")
    if verdict == "FAIL": fail = True
    rows.append((k, f"{d.max():.3e}", rel, nz, (f"{verdict} (tol {lim:.1e})" if lim is not None else "")))
print("| key | max|Δ| | max|Δ|/max|ref| | #diff | 判定 |")
print("|---|---|---|---|---|")
for k, dm, rel, nz, v in rows: print(f"| {k} | {dm} | {rel:.2e} | {nz} | {v} |")
print(f"WORST rel = {worst:.3e} over {len(rows)} keys -> {'FAIL' if fail else ('PASS' if (a.tol is not None or a.tolfile) else 'n/a')}")
if a.dump: json.dump(rels, open(a.dump, "w"), indent=1)
sys.exit(1 if fail else 0)
