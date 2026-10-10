#!/usr/bin/env python3
"""plan architecture-float-state-double-geometry §6.23 (事前登録) の je の後処理の A/B の主判定。c48_prec_series.py で
A = prec_series_dbl.csv (動的 je、共通の double の座標)、B = prec_series_dbl_je76.csv (x = 0.6 の je を 76 に固定) を作った後に回す。
主判定: run_0057_pv_f32a の δ*(0.6) の相対変動幅 (最大 − 最小)/平均 の比 B/A ≤ 0.5 → 「je の切り替えが変動幅の過半を説明する」を支持、> 0.5 → 棄却。
記録: 4 本の δ*(0.6)・θ(0.6) の A・B の変動幅と je の範囲。結果は case/48/c48_je_ab.json。"""
import csv, json
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parents[1]
RUNS = ("run_0057_pv_f32a", "run_0058_pv_f32b", "run_0059_pv_f64a", "run_0060_pv_f64b")
def rng(run, name, k):
    S = [x for x in csv.DictReader(open(HERE / run / name)) if 24000 <= int(x["step"]) <= 48000]
    v = np.array([float(x[k]) for x in S]); j = [int(float(x["je_0.6"])) for x in S]
    if len(v) != 13 or not np.all(np.isfinite(v)): raise SystemExit(f"判定不能: {run} {name} の窓が 13 点でないか非有限")
    return float(np.ptp(v) / abs(v.mean())), [min(j), max(j)]
out = {"runs": {}}
for r in RUNS:
    out["runs"][r] = {}
    for k in ("dstar_0.6", "theta_0.6"):
        (ra, ja), (rb, jb) = rng(r, "prec_series_dbl.csv", k), rng(r, "prec_series_dbl_je76.csv", k)
        out["runs"][r][k] = dict(A=ra, B=rb, ratio=rb / ra if ra > 0 else None, je_A=ja, je_B=jb)
        print(f"  {r} {k:10s} A 幅 {ra:.2e} (je {ja})  B 幅 {rb:.2e} (je {jb})  B/A {rb / ra if ra > 0 else float('nan'):.3f}")
m = out["runs"]["run_0057_pv_f32a"]["dstar_0.6"]
out["verdict"] = ("支持 (je の切り替えが変動幅の過半を説明する)" if m["ratio"] is not None and m["ratio"] <= 0.5
                  else "棄却 (je を固定しても変動幅の半分を超えて残る)")
print(f"== VERDICT §6.23 (je の A/B): {out['verdict']}  — run_0057 の δ*(0.6): B/A = {m['ratio']:.3f}")
(HERE / "c48_je_ab.json").write_text(json.dumps(out, indent=1, ensure_ascii=False)); print("→", HERE / "c48_je_ab.json")
