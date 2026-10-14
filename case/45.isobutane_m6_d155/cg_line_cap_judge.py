"""plan tooling-nozzle-core-grid §4.32 (事前登録、codex 諮問 2026-10-11): ライン陰解法の方向別 dt の上限 50 (A) 対 1 (B) の判定。
usage (AWS で、両 run の終了後): python3 cg_line_cap_judge.py <A の run> <B の run> → _band_ab/core_grid/line_cap_judge.json
  1. 実効の dt: res_1000 の dt_local を A と B で比べる (B が小さい節点がなければ、この A/B は仮説を検査できていない)。
  2. 残差: residual_history の outer_end の行で、窓 10k〜15k と 15k〜20k の各列の中央値。平均流の非ゼロ列 (rms_ro・roUx・roUy・roe・roY0・roY1)
     のすべてで両窓とも B/A ≤ 1/3 → 「dt の拡大の抑制は有望」。SST の列 (roK・roOmega) は併記し、B/A > 1 なら「平均流の停滞は緩和、連成系は未解決」。
     dt は縮んだが条件を満たさない → 「上限 1 への変更だけで 2 万 step 内に停滞は崩れない」。非有限 → 発散。列で割れる → 判定不能。延長しない。
  3. 量: 出力ごとの 8 量 (gref_series_<run>_<n>.csv をまとめる) を、Gref の point の止まった値 (GrefP3 の 180k〜200k) と B0 の窓 (GrefA の 40k〜60k) に対して出す
     (記録。「同じ結果」の識別用の閾値は相対 0.1 %・M 0.005 %pt で、この試験では到達を主張しない)。"""
import csv, glob, json, math, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; OUTD = HERE / "_band_ab" / "core_grid"
MEAN = ("rms_ro", "rms_roUx", "rms_roUy", "rms_roe", "rms_roY0", "rms_roY1"); SST = ("rms_roK", "rms_roOmega")
REQ8 = ("theta_r_w65_75", "theta_r_w84_94", "delta_loc_w65_75", "delta_loc_w84_94", "Q_w", "exit_M_line", "eta0.0_mean_pct", "eta0.1_mean_pct")
A, B = sys.argv[1:3]; out = {"registration": "plans/active/tooling-nozzle-core-grid.md §4.32", "A": A, "B": B}
import h5py
def dtl(run, n):
    with h5py.File(HERE / run / f"res_{n}.h5", "r") as h:
        return np.asarray(h["VALUE/dt_local"][:], float) if "VALUE/dt_local" in h else None
da, db = dtl(A, 1000), dtl(B, 1000)
if da is None or db is None:
    out["dt"] = {"ok": False, "note": "dt_local が無い"}
else:
    r = db / np.where(da > 0, da, np.nan); smaller = int(np.count_nonzero(r < 1 - 1e-9))
    out["dt"] = {"ok": smaller > 0, "nodes": int(r.size), "B_smaller_nodes": smaller, "ratio_min": float(np.nanmin(r)), "ratio_median": float(np.nanmedian(r)),
                 "ratio_max": float(np.nanmax(r)), "B_larger_nodes": int(np.count_nonzero(r > 1 + 1e-9))}
def res(run):
    rows = [x for x in csv.DictReader(open(HERE / run / "residual_history.csv")) if x["phase"] == "outer_end"]
    return {k: np.array([float(x[k]) for x in rows]) for k in ("step",) + MEAN + SST}
RA, RB = res(A), res(B); win = ((10000, 15000), (15000, 20000)); out["residual"] = {}
nonfinite = any(not np.all(np.isfinite(R[k])) for R in (RA, RB) for k in MEAN + SST)
crit_all = True; split = False
for lo, hi in win:
    w = {}
    for k in MEAN + SST:
        ma = float(np.median(RA[k][(RA["step"] >= lo) & (RA["step"] < hi)])); mb = float(np.median(RB[k][(RB["step"] >= lo) & (RB["step"] < hi)]))
        w[k] = {"A": ma, "B": mb, "B_over_A": mb / ma if ma > 0 else None}
    out["residual"][f"{lo}-{hi}"] = w
    oks = [w[k]["B_over_A"] is not None and w[k]["B_over_A"] <= 1 / 3 for k in MEAN]
    crit_all &= all(oks); split |= (any(oks) and not all(oks))
sst_up = any((out["residual"][f"{lo}-{hi}"][k]["B_over_A"] or 0) > 1 for lo, hi in win for k in SST)
# 量
def merged(run):
    rows = []
    for f in sorted(glob.glob(str(OUTD / f"gref_series_{run}_*.csv")), key=lambda p: int(p.rsplit("_", 1)[1][:-4])):
        rows += list(csv.DictReader(open(f)))
    with open(OUTD / f"gref_series_{run}.csv", "w") as fh:
        if rows:
            fh.write(",".join(rows[0].keys()) + "\n"); [fh.write(",".join(r.values()) + "\n") for r in rows]
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]} if rows else None
def tmean(D, k, lo, hi):
    m = (D["step"] >= lo) & (D["step"] <= hi); return float(np.trapezoid(D[k][m], D["step"][m]) / (hi - lo))
def load(nm):
    rows = list(csv.DictReader(open(OUTD / f"gref_series_{nm}.csv"))); return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}
SA, SB = merged(A), merged(B); P3 = load("GrefP3"); GA = load("GrefA")
out["quantities"] = {}
for k in REQ8:
    mp = tmean(P3, k, 180000, 200000); mb0 = tmean(GA, k, 40000, 60000); D0 = abs(mb0)
    f = (lambda v: v) if k.startswith("eta") else (lambda v, D0=D0: 100 * v / D0)
    out["quantities"][k] = {arm: {"steps": S["step"].tolist(), "minus_point_conv": [f(x - mp) for x in S[k]], "minus_B0": [f(x - mb0) for x in S[k]]}
                            for arm, S in (("A", SA), ("B", SB)) if S is not None}
    out["quantities"][k]["B0_minus_point_conv"] = f(mb0 - mp)
out["result"] = ("発散 (非有限)" if nonfinite else "検査できていない (B の dt が縮んでいない)" if not out["dt"]["ok"] else
                 ("dt の拡大の抑制は有望" + ("。ただし SST の残差は B の方が大きい (平均流の停滞は緩和、連成系は未解決)" if sst_up else "")) if crit_all else
                 "判定不能 (平均流の列で結果が割れた)" if split else "上限 1 への変更だけで 2 万 step 内に停滞は崩れない")
(OUTD / "line_cap_judge.json").write_text(json.dumps(out, ensure_ascii=False, indent=1))
print("dt:", out["dt"])
for wk, w in out["residual"].items():
    print(f"残差 {wk}:", "  ".join(f"{k[4:]} {v['B_over_A']:.3g}" for k, v in w.items()))
for k in REQ8:
    q = out["quantities"][k]
    for arm in ("A", "B"):
        if arm in q:
            print(f"  {k:16s} {arm} − point の止まった値:", " ".join(f"{x:+.3f}" for x in q[arm]["minus_point_conv"][::4]), f"(B0 − point {q['B0_minus_point_conv']:+.3f})")
print("結果:", out["result"])
