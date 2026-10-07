"""帯の選び方 E (band_select="edge") のオフライン検証 T1〜T6 (plan verification-m6-axis-wave-mesh-su2 §5.1 #7a)。

合格条件の出典: notes/reviews/2026-10-04-deltastar-extraction-robust-method-diagnose.md §3。
usage: design/.venv-opt/bin/python test_edge_band.py   → _band_ab/edge_tests/summary.json と標準出力の表
"""
import json
import os
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, "/home/sano/work/forge/design")
from forge_design.metrics.deltastar import (deltastar_from_core_matched_euler, band_local_deficit,  # noqa: E402
                                            edge_band_positions, pspline_uniform, _load_structured)
from forge_design.feedback.deltastar_loop import extract_and_merge  # noqa: E402

C = Path(__file__).resolve().parent
OUT = C / "_band_ab" / "edge_tests"
OUT.mkdir(parents=True, exist_ok=True)
X_LO, X_HI = 3.0, 93.0
RES = {}


def snapshot_dir(run, res_name):
    """_load_structured は最後の res を読むので、1 枚だけを見せる一時 run を作る (シンボリックリンク)。"""
    d = Path(tempfile.mkdtemp(prefix=f"{run}_{res_name}_", dir=OUT))
    for f in ("nozzle.h5", "prepare_info.json"):
        os.symlink(C / run / f, d / f)
    os.symlink(C / run / res_name, d / res_name)
    return d


FIELDS = {
    "F1": ("run_0022_ns_ib_pass0", "run_0001_euler_shortest_dry"),
    "F2": ("run_0038_ns_final_rt77p02", "run_0037_euler_rt77p02"),
    "F3": ("run_0042_ns_restart_ctrl", "run_0037_euler_rt77p02"),
}
F3P = {"res_8000": snapshot_dir("run_0042_ns_restart_ctrl", "res_8000.h5"),
       "res_11000": snapshot_dir("run_0042_ns_restart_ctrl", "res_11000.h5")}


def rho(d):
    return d["x"], d["delta_r_raw"] / d["delta_in"]


def window(x):
    return (x >= X_LO) & (x <= X_HI)


def step_metric(x, r):
    """S(x) = median(ρ′,(x,x+2]) − median(ρ′,[x−2,x))、ρ′ = ρ − P-spline(ノット 10 r_t, λ=1)。戻り: (x, S)。"""
    m = window(x) & np.isfinite(r)
    xx, rr = x[m], r[m]
    rp = rr - pspline_uniform(xx, rr, knot=10.0, lam=1.0)
    S = np.full(len(xx), np.nan)
    for i in range(len(xx)):
        a = (xx > xx[i]) & (xx <= xx[i] + 2.0); b = (xx >= xx[i] - 2.0) & (xx < xx[i])
        if a.sum() >= 3 and b.sum() >= 3:
            S[i] = np.median(rp[a]) - np.median(rp[b])
    return xx, S


def t1_vals(x, r):
    xx, S = step_metric(x, r)
    hi = (xx >= 8.0) & np.isfinite(S); lo = (xx < 8.0) & np.isfinite(S)
    k_hi = np.nanargmax(np.abs(np.where(hi, S, np.nan))); k_lo = np.nanargmax(np.abs(np.where(lo, S, np.nan)))
    return dict(maxS_ge8=float(abs(S[k_hi])), x_ge8=float(xx[k_hi]), maxS_lt8=float(abs(S[k_lo])), x_lt8=float(xx[k_lo]))


def ratio_diff(xa, ra, xb, rb):
    m = window(xa)
    rbi = np.interp(xa[m], xb, rb)
    D = ra[m] / rbi - 1.0
    return xa[m], D


# ------------------------------------------------------------------ 抽出 (E 既定・adaptive・位相・スナップショット)
EXT = {}
for F, (ns, eu) in FIELDS.items():
    EXT[(F, "edge")] = deltastar_from_core_matched_euler(C / ns, C / eu, band_select="edge", return_ladders=True)
    EXT[(F, "adaptive")] = deltastar_from_core_matched_euler(C / ns, C / eu)
    EXT[(F, "edge_phase")] = deltastar_from_core_matched_euler(C / ns, C / eu, band_select="edge", edge_y0_shift=1.25 ** 0.125)
for k, d in F3P.items():
    EXT[("F3", f"edge_{k}")] = deltastar_from_core_matched_euler(d, C / "run_0037_euler_rt77p02", band_select="edge")
    EXT[("F3", f"adaptive_{k}")] = deltastar_from_core_matched_euler(d, C / "run_0037_euler_rt77p02")

# ------------------------------------------------------------------ T1 段差
T1 = {}
for F in FIELDS:
    for meth in ("edge", "adaptive"):
        T1[f"{F}/{meth}"] = t1_vals(*rho(EXT[(F, meth)]))
RES["T1"] = dict(vals=T1, pass_=all(T1[f"{F}/edge"]["maxS_ge8"] <= 0.003 and T1[f"{F}/edge"]["maxS_lt8"] <= 0.006 for F in FIELDS))


# ------------------------------------------------------------------ T2 パラメータ感度 (梯子は共通、帯の位置だけ振る)
def remeasure(F, c, eps, win):
    """既定抽出の梯子から帯の位置を作り直し、同じ測定器で δ_r を測る。"""
    ns, eu = FIELDS[F]
    d = EXT[(F, "edge")]; det = d["edge_detail"]
    N = _load_structured(C / ns); E = _load_structured(C / eu)
    xE = E["x"][:, 0]; rwE_col = E["r"][:, -1]
    din = np.array([N["r"][i][-1] for i in det["i"]]) - np.interp(det["x"], xE, rwE_col)
    rwN = np.array([N["r"][i][-1] for i in det["i"]])
    pos = edge_band_positions(det["x"], det["ladders_y"], det["ladders_d"], din, rwN, eps=eps, c=c, window=win)
    out = np.full(len(det["i"]), np.nan)
    for k, i in enumerate(det["i"]):
        x = det["x"][k]
        j = int(np.clip(np.searchsorted(xE, x) - 1, 0, len(xE) - 2))
        w = float(np.clip((x - xE[j]) / max(xE[j + 1] - xE[j], 1e-30), 0.0, 1.0))
        qf = lambda r, j=j, w=w: (1 - w) * np.interp(r, E["r"][j], E["q"][j]) + w * np.interp(r, E["r"][j + 1], E["q"][j + 1])
        rw_e = (1 - w) * rwE_col[j] + w * rwE_col[j + 1]
        rr = band_local_deficit(N["r"][i], N["q"][i], qf, rw_e, y_b_fixed=float(pos["y_b"][k]))
        out[k] = rr["delta_r"] if rr else np.nan
    return det["x"], out / din


T2 = {}
consistency = {}
for F in FIELDS:
    xd, rd = remeasure(F, 1.25, 0.003, 3.0)
    x0, r0 = rho(EXT[(F, "edge")])
    consistency[F] = float(np.nanmax(np.abs(np.interp(x0, xd, rd) / r0 - 1)))
    var = {}
    for c in (1.0, 1.25, 1.56):
        for eps in (0.002, 0.003, 0.005):
            for win in (2.0, 3.0, 5.0):
                var[(c, eps, win)] = remeasure(F, c, eps, win)
    worst = dict(spread=0.0, median=0.0)
    keys = list(var)
    for a in range(len(keys)):
        for b in range(a + 1, len(keys)):
            xa, ra = var[keys[a]]; xb, rb = var[keys[b]]
            xx, D = ratio_diff(xa, ra, xb, rb)
            D = D[np.isfinite(D)]
            med = float(np.median(D)); spr = float(np.max(np.abs(D - med)))
            if spr > worst["spread"]:
                worst.update(spread=spr, spread_pair=[keys[a], keys[b]])
            if abs(med) > abs(worst["median"]):
                worst.update(median=med, median_pair=[keys[a], keys[b]])
    T2[F] = worst
RES["T2"] = dict(vals=T2, consistency_default_vs_function=consistency,
                 pass_=all(T2[F]["spread"] <= 0.003 and abs(T2[F]["median"]) <= 0.01 for F in FIELDS))

# ------------------------------------------------------------------ T3 梯子位相
T3 = {}
for F in FIELDS:
    xx, D = ratio_diff(*rho(EXT[(F, "edge_phase")]), *rho(EXT[(F, "edge")]))
    k = np.nanargmax(np.abs(D)); T3[F] = dict(max=float(abs(D[k])), x=float(xx[k]))
RES["T3"] = dict(vals=T3, pass_=all(T3[F]["max"] <= 0.001 for F in FIELDS))

# ------------------------------------------------------------------ T4 スナップショット
T4 = {}
for k in F3P:
    for meth in ("edge", "adaptive"):
        xx, D = ratio_diff(*rho(EXT[("F3", f"{meth}_{k}")]), *rho(EXT[("F3", meth)]))
        j = np.nanargmax(np.abs(D)); T4[f"{meth}_{k}_vs_res_12000"] = dict(max=float(abs(D[j])), x=float(xx[j]))
RES["T4"] = dict(vals=T4, pass_=all(T4[f"edge_{k}_vs_res_12000"]["max"] <= 0.002 for k in F3P))

# ------------------------------------------------------------------ T5 壁 (extract_and_merge、λ 昇格なし、c=1.0 vs 1.56 の高域差)
T5 = {}
for F, (ns, eu) in FIELDS.items():
    nx = {}
    lam_ok = True
    for c in (1.0, 1.25, 1.56):
        od = OUT / f"T5_{F}_c{c}"
        try:
            s = extract_and_merge(C / ns, C / eu, omega=1.0, out_dir=od, max_lam_factor=1.0, band_select="edge", edge_c=c)
            lam = s["smooth"]["lam"]
        except RuntimeError as e:
            lam = f"FAIL: {e}"; lam_ok = False
            continue
        t = np.loadtxt(od / "delta_r_next.csv", delimiter=",", skiprows=1)
        nx[c] = (t[:, 0], t[:, 1])
    hp = np.nan
    if 1.0 in nx and 1.56 in nx:
        x = nx[1.0][0]; m = window(x)
        d = nx[1.56][1][m] - nx[1.0][1][m]
        d_hp = d - pspline_uniform(x[m], d, knot=20.0, lam=1.0)
        hp = float(np.max(np.abs(d_hp)))
    T5[F] = dict(lam_not_escalated=lam_ok, hp_max_rt=hp)
RES["T5"] = dict(vals=T5, pass_=all(T5[F]["lam_not_escalated"] and T5[F]["hp_max_rt"] <= 3e-4 for F in FIELDS))

# ------------------------------------------------------------------ T6 合成 (縁が格子を横切る)
# design/tests/run_deltastar_tests.py の同名関数の写し (あちらは import すると試験一式を実行して exit する)
def wall_clustered_grid(rw, n=97, first_frac=4.5e-5):
    lo, hi = 1.0001, 1.5
    for _ in range(100):
        g = 0.5 * (lo + hi)
        tot = first_frac * (g ** (n - 1) - 1) / (g - 1)
        lo, hi = (g, hi) if tot < 1.0 else (lo, g)
    dy = first_frac * g ** np.arange(n - 1)
    y = np.concatenate([[0.0], np.cumsum(dy)]); y = y / y[-1]
    return (rw * (1.0 - y))[::-1]


def profile_smooth_bl(r, q_core_fn, rw_ns, delta99, n_pow=7.0):
    y = rw_ns - r
    f = np.where(y < delta99, (np.maximum(y, 0) / delta99) ** (1.0 / n_pow), 1.0)
    return q_core_fn(r) * f


def equiv_delta_exact(r_fine, q_ref, q_ns, rw):
    D = 2 * np.pi * np.trapezoid((q_ref - q_ns) * r_fine, r_fine)
    seg = 0.5 * (q_ref[1:] * r_fine[1:] + q_ref[:-1] * r_fine[:-1]) * np.diff(r_fine)
    F = 2 * np.pi * np.concatenate([[0.0], np.cumsum(seg[::-1])])
    return rw - float(np.interp(D, F, r_fine[::-1]))

xs6 = np.linspace(0.0, 99.0, 100)
rw6 = 5.0 + 3.0 * xs6 / 99.0
dstar6 = 0.02 + 0.6 * (xs6 / 99.0)                    # 真の δ* を滑らかに増やす
ratio99 = 8.0 - 6.0 * (0.5 - 0.5 * np.cos(np.pi * xs6 / 99.0))   # δ99/δ* を 8→2 に滑らかに
qE6 = lambda r_, rw: 1.0 + 0.2 * (r_ / rw) ** 2
truth, ladY, ladD, din6, rwN6, cols = [], [], [], [], [], []
for k, x in enumerate(xs6):
    rw = rw6[k]; r = wall_clustered_grid(rw)
    wave = lambda r_, rw=rw, x=x: qE6(r_, rw) * (1.0 + 0.01 * np.sin(2 * np.pi * r_ / rw + 0.3 * x))
    # δ99 を δ* の比で決める (べき乗則の指数 n は δ*/δ99 = 1/(n+1) から)
    n_pow = max(ratio99[k] - 1.0, 1.0); d99 = ratio99[k] * dstar6[k]
    rf = np.linspace(0, rw, 200001)
    tru = equiv_delta_exact(rf, wave(rf), profile_smooth_bl(rf, wave, rw, d99, n_pow=n_pow), rw)
    qN = profile_smooth_bl(r, wave, rw, d99, n_pow=n_pow)
    qf = lambda r_, rw=rw: qE6(r_, rw)
    y0 = max(1.2 * tru, 0.016 * rw); ys = []
    y = y0
    while y <= 0.5 * rw:
        ys.append(y); y *= 1.25 ** 0.25
    ds = []
    for yb in ys:
        rr = band_local_deficit(r, qN, qf, rw, delta_in=tru, y_b_fixed=yb)
        ds.append(rr["delta_r"] if rr else np.nan)
    truth.append(tru); ladY.append(np.array(ys)); ladD.append(np.array(ds)); din6.append(tru); rwN6.append(rw); cols.append((r, qN, qf, rw))
pos6 = edge_band_positions(xs6, ladY, ladD, np.array(din6), np.array(rwN6))
est6 = np.array([band_local_deficit(r, qN, qf, rw, delta_in=din6[k], y_b_fixed=float(pos6["y_b"][k]))["delta_r"]
                 for k, (r, qN, qf, rw) in enumerate(cols)])
adp6 = np.array([band_local_deficit(r, qN, qf, rw, delta_in=din6[k])["delta_r"] for k, (r, qN, qf, rw) in enumerate(cols)])
truth = np.array(truth)
r6 = est6 / truth; ra6 = adp6 / truth
def jumps(v):
    return float(np.max(np.abs(np.diff(v))))
xx6, S6 = step_metric(xs6, r6); _, Sa6 = step_metric(xs6, ra6)
RES["T6"] = dict(edge=dict(adj_jump_max=jumps(r6), maxS=float(np.nanmax(np.abs(S6))), bias_min=float(r6.min() - 1), bias_max=float(r6.max() - 1)),
                 adaptive=dict(adj_jump_max=jumps(ra6), maxS=float(np.nanmax(np.abs(Sa6))), bias_min=float(ra6.min() - 1), bias_max=float(ra6.max() - 1)),
                 note="真値 = 波込みコアを参照にした厳密等価排除厚。δ_in = 真値 (固定点)。x 間隔 1 なので S の窓 ±2 は各側 2 点 → step_metric は 3 点未満で NaN になるため、ここでは S を隣接跳びで代用")
RES["T6"]["pass_"] = (RES["T6"]["edge"]["adj_jump_max"] <= 0.001 and abs(RES["T6"]["edge"]["bias_min"]) <= 0.03 and abs(RES["T6"]["edge"]["bias_max"]) <= 0.03)

(OUT / "summary.json").write_text(json.dumps(RES, indent=1, default=str))
print(json.dumps(RES, indent=1, default=str))
