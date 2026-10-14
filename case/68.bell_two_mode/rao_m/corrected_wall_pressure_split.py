#!/usr/bin/env python3
r"""③ベル (M) L100: δ* 補正壁の壁圧差の分解 — plan tooling-nozzle-bell-two-mode §5.1 #7 の診断 (codex (diagnose) 2026-10-13 の A/B)。

壁 w と共通の点 s (判定の edM と同じ r_s = r_wall,NS(x) − y_b(x)、y_b は run_0019 の edge の帯の位置で凍結、一定 x の断面に補間)で
  e_w    = (p_Nw − p_Ew)/p_Ew                 (現行の壁圧の指標 pw の被積分)
  e_s    = (p_Ns − p_Es)/p_Ew                 (共通座標での圧力差)
  e_near = [(p_Nw − p_Ns) − (p_Ew − p_Es)]/p_Ew  (共通の点から各々の壁までの圧力増分の差)
  e_w = e_s + e_near
を x 1.17〜L・Δx 重み・RMS (`judge_corrected_wall` と同じ) で集計する。診断の線 (凍結した合否は変えない): A = e_w − e_near (= e_s) の RMS が
e_w の RMS の 30 % 以下かつ e_near の平均が負、B = それ以外。組: U (run_0019/run_0018)・V-ext (run_0026/run_0024)・V-zero (run_0027/run_0024)・
V-ext の狙い細格子 (run_0026/run_0025)。`--res-ns` で NS の別のスナップショット (時系列) を読む。出力 corrected_wall_pressure_split.json。"""
import argparse
import json
import sys
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import judge_corrected_wall as j  # noqa: E402

PAIRS = {"U": ("run_0019_tp_m_ns_l100", "run_0018_tp_m_euler_l100_pin"),
         "Vext": ("run_0026_tp_m_ns_l100_corr_vext", "run_0024_tp_m_euler_l100_inv"),
         "Vzero": ("run_0027_tp_m_ns_l100_corr_vzero", "run_0024_tp_m_euler_l100_inv"),
         "Vext_fine": ("run_0026_tp_m_ns_l100_corr_vext", "run_0025_tp_m_euler_l100_inv_fine")}
REGS = [(1.17, 2), (2, 4), (4, 8), (8, 13), (13, 18), (18, 99)]


def load(run: Path, res: str | None = None) -> dict:
    d = j.load(run)
    if res is not None:
        with h5py.File(run / res) as f:
            d["P"] = f["VALUE/P"][:].reshape(d["P"].shape)
        d["res"] = res
    xt = d["x"][:, 0]
    r = np.empty_like(d["r"]); P = np.empty_like(d["P"])
    for jj in range(d["x"].shape[1]):
        r[:, jj] = np.interp(xt, d["x"][:, jj], d["r"][:, jj]); P[:, jj] = np.interp(xt, d["x"][:, jj], d["P"][:, jj])
    d.update(xt=xt, rc=r, Pc=P)
    return d


def at_point(d, xs, rs):
    out = np.empty(len(xs))
    xt = d["xt"]
    for k, (xx, rr) in enumerate(zip(xs, rs)):
        i = int(np.clip(np.searchsorted(xt, xx), 1, len(xt) - 1)); a = (xx - xt[i - 1]) / (xt[i] - xt[i - 1])
        rc = (1 - a) * d["rc"][i - 1] + a * d["rc"][i]; Pc = (1 - a) * d["Pc"][i - 1] + a * d["Pc"][i]
        if rr > rc[-1]:
            raise ValueError(f"{d['run']}: 共通の点 r {rr} が壁 {rc[-1]} の外 (x {xx})")
        out[k] = np.interp(rr, rc, Pc)
    return out


def split(ns, eu, yb_x, yb):
    L = min(ns["x"][-1, -1], eu["x"][-1, -1])
    xg = np.linspace(0.0, L, j.NX); w = j.trap_w(xg); cal = xg >= j.X_CAL
    xs, ws = xg[cal], w[cal]
    pNw = np.interp(xs, ns["x"][:, -1], ns["P"][:, -1]); pEw = np.interp(xs, eu["x"][:, -1], eu["P"][:, -1])
    rs = np.interp(xs, ns["xt"], ns["rc"][:, -1]) - np.interp(xs, yb_x, yb)
    pNs, pEs = at_point(ns, xs, rs), at_point(eu, xs, rs)
    e_w = pNw / pEw - 1.0; e_s = (pNs - pEs) / pEw; e_near = ((pNw - pNs) - (pEw - pEs)) / pEw
    out = {"res_ns": ns["res"], "res_eu": eu["res"], "rms_e_w": j.wrms(e_w, ws), "rms_e_s": j.wrms(e_s, ws), "rms_e_near": j.wrms(e_near, ws),
           "mean_e_w": float(np.average(e_w, weights=ws)), "mean_e_s": float(np.average(e_s, weights=ws)),
           "mean_e_near": float(np.average(e_near, weights=ws)), "closure_max": float(np.max(np.abs(e_w - e_s - e_near)))}
    out["ratio_rms_s_over_w"] = out["rms_e_s"] / out["rms_e_w"]
    out["diag"] = "A" if (out["ratio_rms_s_over_w"] <= 0.30 and out["mean_e_near"] < 0) else "B"
    out["regions"] = []
    for lo, hi in REGS:
        m = (xs >= lo) & (xs < hi)
        out["regions"].append({"x": [lo, min(hi, float(L))], "mean_e_w": float(np.mean(e_w[m])), "mean_e_s": float(np.mean(e_s[m])),
                               "mean_e_near": float(np.mean(e_near[m]))})
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--res-ns", default=None, help="NS のスナップショット名 (既定は最終)")
    ap.add_argument("--pairs", nargs="*", default=list(PAIRS))
    ap.add_argument("--out", default=str(HERE / "corrected_wall_pressure_split.json"))
    a = ap.parse_args()
    ed = np.load(j.EDGE, allow_pickle=True)
    res = {"plan": "plans/active/tooling-nozzle-bell-two-mode.md §5.1 #7", "consult": "notes/reviews/2026-10-13-bell-m-corrected-wall-result-diagnose.md",
           "line": "A: rms(e_s) ≤ 0.30·rms(e_w) かつ mean(e_near) < 0", "pairs": {}}
    for n in a.pairs:
        nsr, eur = PAIRS[n]
        ns = load(j.CASE / nsr, a.res_ns); eu = load(j.CASE / eur)
        res["pairs"][n] = split(ns, eu, ed["x"], ed["band_y_b"])
        r = res["pairs"][n]
        print(f"{n:9s} rms e_w {r['rms_e_w']*100:.3f}%  e_s {r['rms_e_s']*100:.3f}%  e_near {r['rms_e_near']*100:.3f}%  "
              f"mean w/s/near {r['mean_e_w']*100:+.3f}/{r['mean_e_s']*100:+.3f}/{r['mean_e_near']*100:+.3f}%  s/w {r['ratio_rms_s_over_w']:.2f} → {r['diag']}"
              f"  (closure {r['closure_max']:.1e})")
        print("          " + "  ".join(f"[{g['x'][0]:g},{g['x'][1]:.4g}) {g['mean_e_w']*100:+.2f}={g['mean_e_s']*100:+.2f}{g['mean_e_near']*100:+.2f}" for g in r["regions"]))
    Path(a.out).write_text(json.dumps(res, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
