#!/usr/bin/env python3
r"""③ベル (M): δ* 補正の縁を MOC の壁 (+ 凍結源 Euler) に替えても kδ = 1.0289 を据え置けるか — plan tooling-nozzle-bell-two-mode §5.1 #7c (7c-0)
の事前登録の A/B (codex (diagnose) 2026-10-13 option-design)。変えるのは縁の M の供給元だけ (CFD 0 step)。

固定: 受理済みの NS 抽出値 d (edge の帯・一定 x、L100 = run_0019、L70 = run_0020)、区間 (L100 = 較正の区間 = 帯の安定 ∧ M3 STEADY、
L70 = 盲検の判定の区間 — `judge_l70_blind.stable_mask`・`steady_mask` と同じ)、重み (Δx の台形)、kδ (`kdelta_l100_calibration.json`)。
比べるもの: C_euler (較正・盲検で使った contur の npz) と C_moc (`check_edge_moc_vs_euler.moc_edge_C`、同じ x 格子) を抽出の x へ同じ線形補間。
誤差 e = kδ·C/d − 1 の出口・相対 RMS・局所の最大。判定 (L70): A = MOC の縁でも出口 ≤ 10 %・RMS ≤ 10 %・局所 ≤ 20 %、B = 外れる。
出力 check_edge_moc_transfer.json。"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import check_edge_moc_vs_euler as ce  # noqa: E402
import judge_l70_blind as jl  # noqa: E402

CASES = {"L100": dict(ns_npz="dstar_l100_ns_run0019_constx_edge.npz", series="dstar_series_l100_edge.npz",
                      contur="contur_l100_euler_run0018.npz", problem="problem_bell_tp_m_ns_l100.yaml",
                      wall="rao_m/rao_pin_parts/L100_pin_n321_wall.npz", eu="run_0018_tp_m_euler_l100_pin"),
         "L70": dict(ns_npz="dstar_l70_ns_run0020_constx_edge.npz", series="dstar_series_l70_edge.npz",
                     contur="contur_l70_euler_run0023.npz", problem="problem_bell_tp_m_ns_l70.yaml",
                     wall="rao_m/rao_pin_parts/L70_pin_n321_wall.npz", eu="run_0023_tp_m_euler_l70_pin")}


def errs(e, J, w):
    return {"exit_err": float(e[-1]), "rms": float(np.sqrt(np.average(e[J] ** 2, weights=w[J]))), "local_max": float(np.max(np.abs(e[J]))),
            "median": float(np.median(e[J]))}


def main():
    kd = float(json.loads((HERE / "kdelta_l100_calibration.json").read_text())["k_delta"])
    out = {"k_delta": kd, "lines_L70": jl.LINES, "cases": {}}
    for tag, c in CASES.items():
        ed = np.load(HERE / c["ns_npz"], allow_pickle=True)
        x, d = ed["x"], ed["delta_r_use"]
        stab, _ = jl.stable_mask(ed)
        steady, nsnap = jl.steady_mask(x, stab, HERE / c["series"])
        J = stab & steady
        w = np.zeros_like(x); w[1:-1] = 0.5 * (x[2:] - x[:-2]); w[0] = 0.5 * (x[1] - x[0]); w[-1] = 0.5 * (x[-1] - x[-2])
        co = np.load(HERE / c["contur"], allow_pickle=True)
        xm, Cm, Ce = ce.moc_edge_C(c["problem"], c["wall"], c["eu"])
        if not np.allclose(np.interp(xm, co["x"], co["delta_r_ann"]), Ce):
            raise SystemExit(f"{tag}: contur の npz と再計算の C_euler が一致しない")
        Ceu = np.interp(x, co["x"], co["delta_r_ann"]); Cmo = np.interp(x, xm, Cm)
        r = {"n_J": int(J.sum()), "n_snapshots": nsnap, "x_J_min": float(x[J].min()),
             "euler_edge": errs(kd * Ceu / d - 1.0, J, w), "moc_edge": errs(kd * Cmo / d - 1.0, J, w),
             "C_moc_over_C_euler_J": {"min": float(np.min(Cmo[J] / Ceu[J])), "max": float(np.max(Cmo[J] / Ceu[J]))}}
        if tag == "L70":
            m = r["moc_edge"]
            r["verdict"] = "A" if (abs(m["exit_err"]) <= jl.LINES["exit"] and m["rms"] <= jl.LINES["rms"] and m["local_max"] <= jl.LINES["local"]) else "B"
        out["cases"][tag] = r
        print(tag, json.dumps(r, ensure_ascii=False))
    (HERE / "check_edge_moc_transfer.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
