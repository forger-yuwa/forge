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
F3P_UNUSED = {"res_8000": snapshot_dir("run_0042_ns_restart_ctrl", "res_8000.h5"),
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



EXT = {}
for F, (ns, eu) in FIELDS.items():
    EXT[(F, "edge")] = deltastar_from_core_matched_euler(C / ns, C / eu, band_select="edge", return_ladders=True)
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
    var = {(c, e, w): remeasure(F, c, e, w) for c in (1.25, 1.56) for e in (0.002, 0.003, 0.005) for w in (2.0, 3.0, 5.0)}
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
RES["T2pp"] = dict(vals=T2, consistency_default_vs_function=consistency, pass_=all(T2[F]["pass"] for F in FIELDS))


(OUT / "summary_T2pp.json").write_text(json.dumps(RES, indent=1, default=str))
print(json.dumps(RES, indent=1, default=str))
