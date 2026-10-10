#!/usr/bin/env python3
"""case/48 の精度の A/B (plan architecture-float-state-double-geometry §6.21) の判定量の時系列を、各 res_<N>.h5 から作る。

    python3 c48_prec_series.py <run_dir> [--steps 2000:48000:2000]   → <run_dir>/prec_series.csv

判定量は `cooled_plate_eval.py` の定義をそのまま使う (同じ関数を呼ぶ):
- x = 0.3 / 0.6 / 0.9 の θ・δ*・Cf・q_w (`station()`: 壁法線の列、ρu が 99.5 % に達する最初の点 je を縁とする台形積分、壁の量は 3 点の片側差分)
- CD = ∫Cf dx・HF = ∫q_w dx (`plate_integrals()`: x ∈ [0.002, 0.998] の 200 点)
記録として各断面の縁の点 je と縁の高さ y_edge も出す (je の飛びを見るため)。
指定した step の res が 1 つでも欠ける・判定量に非有限があると、CSV を書かずに終了コード 2 で止める (判定不能にするため)。"""
import argparse, csv, math, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import cooled_plate_eval as E  # noqa: E402

XS = (0.3, 0.6, 0.9)
QTY = ("theta", "dstar", "Cf", "qw")
ap = argparse.ArgumentParser(); ap.add_argument("run"); ap.add_argument("--steps", default="2000:48000:2000"); a = ap.parse_args()
s0, s1, ds = (int(v) for v in a.steps.split(":"))
run = Path(a.run); steps = list(range(s0, s1 + 1, ds))
miss = [n for n in steps if not (run / f"res_{n}.h5").exists()]
if miss: sys.exit(print(f"[c48_prec_series] {run.name}: res が無い step {miss[:5]}… ({len(miss)} 個)") or 2)
cols = ["step"] + [f"{q}_{x}" for x in XS for q in QTY] + ["CD", "HF"] + [f"je_{x}" for x in XS] + [f"yedge_{x}" for x in XS]
rows = []
for n in steps:
    D = E.load_forge(run / f"res_{n}.h5")
    r = {"step": n}
    for x in XS:
        st = E.station(D, x)
        for q in QTY: r[f"{q}_{x}"] = st[q]
        r[f"je_{x}"] = st["je"]; r[f"yedge_{x}"] = st["y_edge"]
    pi = E.plate_integrals(D); r["CD"] = pi["CD"]; r["HF"] = pi["HF"]
    bad = [k for k, v in r.items() if not math.isfinite(float(v))]
    if bad: sys.exit(print(f"[c48_prec_series] {run.name}: step {n} の {bad} が非有限") or 2)
    rows.append(r)
with open(run / "prec_series.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols); w.writeheader(); [w.writerow(r) for r in rows]
print(f"[c48_prec_series] {run.name}: {len(rows)} 点 → {run / 'prec_series.csv'}")
