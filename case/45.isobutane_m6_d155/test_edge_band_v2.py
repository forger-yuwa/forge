"""帯の選び方 E のオフライン検証 v2 (plan verification-m6-axis-wave-mesh-su2 §6「事前登録 v2」)。

初回 (test_edge_band.py) からの変更は試験器の誤りの修正だけ (しきい値は据え置き):
T1′ 分母を積分法 (CONTUR) に、x<15 は絶対量 / T2′ c∈{1.12,1.25,1.56}・窓内のみ / T3 位相は新梯子の半刻み /
T5′ c=1.25 vs 1.56、判定 x∈[3,50] ([50,93] は併記) / T6′ 実機レンジの合成 (Δx=0.1 で 900 断面)。
usage: design/.venv-opt/bin/python test_edge_band_v2.py  → _band_ab/edge_tests_v2/summary.json
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
OUT = C / "_band_ab" / "edge_tests_v2"
OUT.mkdir(parents=True, exist_ok=True)
X_LO, X_HI = 3.0, 93.0
X_ABS = 15.0           # これより上流は絶対量 (|S|·δ ≤ 3e-4 r_t) で判定
ABS_TOL = 3e-4
REL_TOL = 0.003
RES = {}
CT = np.loadtxt(C / "run_0022_ns_ib_pass0/delta_r_initial.csv", delimiter=",", skiprows=1)
contur = lambda x: np.interp(x, CT[:, 0], CT[:, 1])
STEP = 1.25 ** 0.125


def snapshot_dir(run, res_name):
    d = Path(tempfile.mkdtemp(prefix=f"{run}_{res_name}_", dir=OUT))
    for f in ("nozzle.h5", "prepare_info.json"):
        os.symlink(C / run / f, d / f)
    os.symlink(C / run / res_name, d / res_name)
    return d


FIELDS = {"F1": ("run_0022_ns_ib_pass0", "run_0001_euler_shortest_dry"),
          "F2": ("run_0038_ns_final_rt77p02", "run_0037_euler_rt77p02"),
          "F3": ("run_0042_ns_restart_ctrl", "run_0037_euler_rt77p02")}
F3P = {"res_8000": snapshot_dir("run_0042_ns_restart_ctrl", "res_8000.h5"),
       "res_11000": snapshot_dir("run_0042_ns_restart_ctrl", "res_11000.h5")}


def win(x):
    return (x >= X_LO) & (x <= X_HI)


def step_metric(x, r):
    m = win(x) & np.isfinite(r)
    xx, rr = x[m], r[m]
    rp = rr - pspline_uniform(xx, rr, knot=10.0, lam=1.0)
    S = np.full(len(xx), np.nan)
    for i in range(len(xx)):
        a = (xx > xx[i]) & (xx <= xx[i] + 2.0); b = (xx >= xx[i] - 2.0) & (xx < xx[i])
        if a.sum() >= 3 and b.sum() >= 3:
            S[i] = np.median(rp[a]) - np.median(rp[b])
    return xx, S


def split_judge(xx, rel, delta_ref):
    """x<15: 絶対 |rel|·δ ≤ 3e-4 r_t、x≥15: 相対 ≤ 0.3 %。戻り: dict(最大値と位置, pass)。"""
    out = {}
    lo = (xx < X_ABS) & np.isfinite(rel); hi = (xx >= X_ABS) & np.isfinite(rel)
    a = np.abs(rel[lo]) * delta_ref[lo]; j = int(np.argmax(a))
    out["abs_lt15_rt"] = float(a[j]); out["x_abs"] = float(xx[lo][j])
    b = np.abs(rel[hi]); k = int(np.argmax(b))
    out["rel_ge15"] = float(b[k]); out["x_rel"] = float(xx[hi][k])
    out["pass"] = bool(out["abs_lt15_rt"] <= ABS_TOL and out["rel_ge15"] <= REL_TOL)
    return out


def rho_c(d):
    """T1′ の量: δ_r,raw / δ_CONTUR(x)。"""
    return d["x"], d["delta_r_raw"] / contur(d["x"])


def diff_on_window(xa, ra, xb, rb):
    m = win(xa)
    return xa[m], ra[m] / np.interp(xa[m], xb, rb) - 1.0


# ------------------------------------------------------------------ 抽出
EXT = {}
for F, (ns, eu) in FIELDS.items():
    EXT[(F, "edge")] = deltastar_from_core_matched_euler(C / ns, C / eu, band_select="edge", return_ladders=True)
    EXT[(F, "adaptive")] = deltastar_from_core_matched_euler(C / ns, C / eu)
    EXT[(F, "edge_phase")] = deltastar_from_core_matched_euler(C / ns, C / eu, band_select="edge", edge_y0_shift=STEP ** 0.5)
for k, d in F3P.items():
    EXT[("F3", f"edge_{k}")] = deltastar_from_core_matched_euler(d, C / "run_0037_euler_rt77p02", band_select="edge")
    EXT[("F3", f"adaptive_{k}")] = deltastar_from_core_matched_euler(d, C / "run_0037_euler_rt77p02")

# 縁推定の成否 (y* NaN の数、x 区間ごと)
RES["edge_ystar_nan"] = {}
for F in FIELDS:
    det = EXT[(F, "edge")]["edge_detail"]; x = det["x"]; ys = det["y_star"]
    RES["edge_ystar_nan"][F] = {f"[{a},{b})": int(np.isnan(ys[(x >= a) & (x < b)]).sum()) for a, b in ((-13, 0), (0, 3), (3, 15), (15, 50), (50, 96))}

# ------------------------------------------------------------------ T1′
T1 = {}
for F in FIELDS:
    for meth in ("edge", "adaptive"):
        xx, S = step_metric(*rho_c(EXT[(F, meth)]))
        T1[f"{F}/{meth}"] = split_judge(xx, S, contur(xx))
RES["T1p"] = dict(vals=T1, pass_=all(T1[f"{F}/edge"]["pass"] for F in FIELDS))


# ------------------------------------------------------------------ T2′
def remeasure(F, c, eps, wnd):
    ns, eu = FIELDS[F]
    det = EXT[(F, "edge")]["edge_detail"]
    N = _load_structured(C / ns); E = _load_structured(C / eu)
    xE = E["x"][:, 0]; rwE_col = E["r"][:, -1]
    rwN = np.array([N["r"][i][-1] for i in det["i"]])
    din = rwN - np.interp(det["x"], xE, rwE_col)
    pos = edge_band_positions(det["x"], det["ladders_y"], det["ladders_d"], din, rwN, eps=eps, c=c, window=wnd)
    out = np.full(len(det["i"]), np.nan)
    for k, i in enumerate(det["i"]):
        x = det["x"][k]
        if not (X_LO <= x <= X_HI):
            continue
        j = int(np.clip(np.searchsorted(xE, x) - 1, 0, len(xE) - 2))
        w = float(np.clip((x - xE[j]) / max(xE[j + 1] - xE[j], 1e-30), 0.0, 1.0))
        qf = lambda r, j=j, w=w: (1 - w) * np.interp(r, E["r"][j], E["q"][j]) + w * np.interp(r, E["r"][j + 1], E["q"][j + 1])
        rw_e = (1 - w) * rwE_col[j] + w * rwE_col[j + 1]
        rr = band_local_deficit(N["r"][i], N["q"][i], qf, rw_e, y_b_fixed=float(pos["y_b"][k]))
        out[k] = rr["delta_r"] if rr else np.nan
    return det["x"], out


T2 = {}; consistency = {}
for F in FIELDS:
    xd, dd = remeasure(F, 1.25, 0.003, 3.0)
    d0 = EXT[(F, "edge")]; m = win(d0["x"])
    consistency[F] = float(np.nanmax(np.abs(np.interp(d0["x"][m], xd, dd) / d0["delta_r_raw"][m] - 1)))
    var = {(c, e, w): remeasure(F, c, e, w) for c in (1.12, 1.25, 1.56) for e in (0.002, 0.003, 0.005) for w in (2.0, 3.0, 5.0)}
    keys = list(var); worst = dict(abs_lt15_rt=0.0, rel_ge15=0.0, median=0.0)
    for a in range(len(keys)):
        for b in range(a + 1, len(keys)):
            xx, D = diff_on_window(*var[keys[a]], *var[keys[b]])
            fin = np.isfinite(D); xx, D = xx[fin], D[fin]
            med = float(np.median(D)); j = split_judge(xx, D - med, contur(xx))
            if j["abs_lt15_rt"] > worst["abs_lt15_rt"]:
                worst.update(abs_lt15_rt=j["abs_lt15_rt"], abs_pair=[keys[a], keys[b]], x_abs=j["x_abs"])
            if j["rel_ge15"] > worst["rel_ge15"]:
                worst.update(rel_ge15=j["rel_ge15"], rel_pair=[keys[a], keys[b]], x_rel=j["x_rel"])
            if abs(med) > abs(worst["median"]):
                worst.update(median=med, median_pair=[keys[a], keys[b]])
    worst["pass"] = bool(worst["abs_lt15_rt"] <= ABS_TOL and worst["rel_ge15"] <= REL_TOL and abs(worst["median"]) <= 0.01)
    T2[F] = worst
RES["T2p"] = dict(vals=T2, consistency_default_vs_function=consistency, pass_=all(T2[F]["pass"] for F in FIELDS))

# ------------------------------------------------------------------ T3 位相
T3 = {}
for F in FIELDS:
    a = EXT[(F, "edge_phase")]; b = EXT[(F, "edge")]
    xx, D = diff_on_window(a["x"], a["delta_r_raw"], b["x"], b["delta_r_raw"])
    k = int(np.nanargmax(np.abs(D)))
    xs_, S = step_metric(xx, D)
    T3[F] = dict(max=float(abs(D[k])), x=float(xx[k]), D_phase_step_max=float(np.nanmax(np.abs(S))))
RES["T3"] = dict(vals=T3, pass_=all(T3[F]["max"] <= 0.001 for F in FIELDS))

# ------------------------------------------------------------------ T4
T4 = {}
for k in F3P:
    for meth in ("edge", "adaptive"):
        a = EXT[("F3", f"{meth}_{k}")]; b = EXT[("F3", meth)]
        xx, D = diff_on_window(a["x"], a["delta_r_raw"], b["x"], b["delta_r_raw"])
        j = int(np.nanargmax(np.abs(D))); T4[f"{meth}_{k}_vs_res_12000"] = dict(max=float(abs(D[j])), x=float(xx[j]))
RES["T4"] = dict(vals=T4, pass_=all(T4[f"edge_{k}_vs_res_12000"]["max"] <= 0.002 for k in F3P))

# ------------------------------------------------------------------ T5′
T5 = {}
for F, (ns, eu) in FIELDS.items():
    nx = {}; lam_ok = True
    for c in (1.25, 1.56):
        od = OUT / f"T5_{F}_c{c}"
        try:
            extract_and_merge(C / ns, C / eu, omega=1.0, out_dir=od, max_lam_factor=1.0, band_select="edge", edge_c=c)
        except RuntimeError:
            lam_ok = False; continue
        t = np.loadtxt(od / "delta_r_next.csv", delimiter=",", skiprows=1); nx[c] = (t[:, 0], t[:, 1])
    v = dict(lam_not_escalated=lam_ok)
    if len(nx) == 2:
        x = nx[1.25][0]; m = win(x)
        d = nx[1.56][1][m] - nx[1.25][1][m]; hp = d - pspline_uniform(x[m], d, knot=20.0, lam=1.0); xm = x[m]
        v.update(hp_3_50=float(np.abs(hp[xm < 50]).max()), hp_50_93_info=float(np.abs(hp[xm >= 50]).max()))
    v["pass"] = bool(lam_ok and v.get("hp_3_50", np.inf) <= ABS_TOL)
    T5[F] = v
RES["T5p"] = dict(vals=T5, pass_=all(T5[F]["pass"] for F in FIELDS))


# ------------------------------------------------------------------ T6′ 実機レンジの合成
def wall_clustered_grid(rw, n=97, first_frac=4.5e-5):
    lo, hi = 1.0001, 1.5
    for _ in range(100):
        g = 0.5 * (lo + hi)
        tot = first_frac * (g ** (n - 1) - 1) / (g - 1)
        lo, hi = (g, hi) if tot < 1.0 else (lo, g)
    dy = first_frac * g ** np.arange(n - 1)
    y = np.concatenate([[0.0], np.cumsum(dy)]); y = y / y[-1]
    return (rw * (1.0 - y))[::-1]


def profile_bl(r, q_core, rw, d99, n_pow):
    y = rw - r
    return q_core * np.where(y < d99, (np.maximum(y, 0) / d99) ** (1.0 / n_pow), 1.0)


def equiv_delta_exact(rf, q_ref, q_ns, rw):
    D = 2 * np.pi * np.trapezoid((q_ref - q_ns) * rf, rf)
    seg = 0.5 * (q_ref[1:] * rf[1:] + q_ref[:-1] * rf[:-1]) * np.diff(rf)
    Fc = 2 * np.pi * np.concatenate([[0.0], np.cumsum(seg[::-1])])
    return rw - float(np.interp(D, Fc, rf[::-1]))


W22 = np.loadtxt(C / "run_0022_ns_ib_pass0/wall_physical.csv", delimiter=",", skiprows=1) / 0.0771
xs6 = np.arange(0.0, 90.0, 0.1)
truth, ladY, ladD, din6, rwN6, cols = [], [], [], [], [], []
for x in xs6:
    rw = float(np.interp(x, W22[:, 0], W22[:, 1])); r = wall_clustered_grid(rw)
    dstar = float(contur(x)); ratio99 = 2.5 + 5.5 * np.exp(-x / 4.0)
    n_pow = max(ratio99 - 1.0, 1.0); d99 = ratio99 * dstar
    qE = lambda r_, rw=rw: 1.0 + 0.2 * (r_ / rw) ** 2
    core = lambda r_, rw=rw, x=x: 0.999 * qE(r_) * (1.0 + 0.002 * np.sin(2 * np.pi * r_ / (1.3 * rw) + 0.05 * x))
    rf = np.linspace(0, rw, 100001)
    tru = equiv_delta_exact(rf, core(rf), profile_bl(rf, core(rf), rw, d99, n_pow), rw)
    qN = profile_bl(r, core(r), rw, d99, n_pow)
    y0 = max(1.2 * tru, 0.016 * rw); ys = []; y = y0
    while y <= 0.5 * rw:
        ys.append(y); y *= STEP
    ds = []
    for yb in ys:
        rr = band_local_deficit(r, qN, qE, rw, delta_in=tru, y_b_fixed=yb)
        ds.append(rr["delta_r"] if rr else np.nan)
    truth.append(tru); ladY.append(np.array(ys)); ladD.append(np.array(ds)); din6.append(tru); rwN6.append(rw); cols.append((r, qN, qE, rw))
truth = np.array(truth)
pos6 = edge_band_positions(xs6, ladY, ladD, np.array(din6), np.array(rwN6))
est6 = np.array([band_local_deficit(r, qN, qE, rw, delta_in=din6[k], y_b_fixed=float(pos6["y_b"][k]))["delta_r"] for k, (r, qN, qE, rw) in enumerate(cols)])
adp6 = np.array([band_local_deficit(r, qN, qE, rw, delta_in=din6[k])["delta_r"] for k, (r, qN, qE, rw) in enumerate(cols)])
r6, ra6 = est6 / truth, adp6 / truth
m6 = xs6 >= X_LO


def t6(v):
    xx, S = step_metric(xs6, v)
    return dict(adj_jump_max=float(np.max(np.abs(np.diff(v[m6])))), maxS=float(np.nanmax(np.abs(S))),
                bias_min=float(v[m6].min() - 1), bias_max=float(v[m6].max() - 1),
                ystar_nan=int(np.isnan(pos6["y_star"][m6]).sum()))


RES["T6p"] = dict(edge=t6(r6), adaptive=t6(ra6))
RES["T6p"]["pass_"] = bool(RES["T6p"]["edge"]["adj_jump_max"] <= 0.001 and RES["T6p"]["edge"]["maxS"] <= 0.002
                           and abs(RES["T6p"]["edge"]["bias_min"]) <= 0.03 and abs(RES["T6p"]["edge"]["bias_max"]) <= 0.03)
np.savetxt(OUT / "T6p_series.csv", np.c_[xs6, truth, r6, ra6, pos6["y_star"] / truth, pos6["y_b"] / truth], delimiter=",",
           header="x,truth,edge_over_truth,adaptive_over_truth,ystar_over_delta,yb_over_delta", comments="")

RES["all_pass"] = all(RES[k]["pass_"] for k in ("T1p", "T2p", "T3", "T4", "T5p", "T6p"))
(OUT / "summary.json").write_text(json.dumps(RES, indent=1, default=str))
print(json.dumps(RES, indent=1, default=str))
