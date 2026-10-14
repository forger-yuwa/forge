#!/usr/bin/env python3
r"""③ベル (M): δ* の較正 (kδ・contur_v2) を L70 で盲検に検証する — plan tooling-nozzle-bell-two-mode §5.1 #7 の凍結 (`e69d56e7`) どおり。

凍結した設定: kδ = 1.0289 (幅 1.0280〜1.0289、`kdelta_l100_calibration.json`)、理論式 = contur_v2 (`contur_l70_euler_run0023.npz`)、
抽出 = edge の帯・一定 x (`dstar_l70_ns_run0020_constx_edge.npz`)、判定の区間 = 帯の安定 (s_band ≤ 2 %・y_b ≥ y*・梯子 [y_b, 2y_b] ≤ 2 %)
∧ M3 STEADY (`dstar_series_l70_edge.npz` に tail 0.4・drift 0.005・osc 0.01・min 20・漸近値の最終比 ≤ 1 %)、誤差 e = (kδ·C − d)/d
(分母の床 max(s_band, 1 %)·d で割った値も併記)、重み = Δx の台形。
L70 自身のゲート: 判定の区間が出口を含み、x 方向に連続した区間を持つこと (抽出のゲートに落ちたら移植性は判定不能、L70 で較正し直さない)。
既登録の線: 出口の誤差 ≤ 10 %・相対 RMS ≤ 10 %・局所 ≤ 20 %。kδ の幅の両端でも判定する。出力: l70_blind_judgement.json。"""
import json
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "solver_density_cuda" / "tools"))
import check_quasisteady as cq  # noqa: E402

LINES = {"exit": 0.10, "rms": 0.10, "local": 0.20}


def stable_mask(ed):
    x = ed["x"]; d = ed["delta_r_use"]; fac, rel = ed["yb_factors"], ed["rel_yb"]
    s_band = np.nanmax(np.abs(rel[:, [k for k, f in enumerate(fac) if f > 1]]), axis=1)
    yb, ys = ed["band_y_b"], ed["edge_y_star"]
    lad = np.full(len(x), np.nan)
    for i in range(len(x)):
        y, dd = ed["ladder_y"][i], ed["ladder_d"][i]
        ok = np.isfinite(dd) & (dd > 0) & (y >= yb[i]) & (y <= 2 * yb[i])
        if ok.sum() >= 3:
            lad[i] = np.nanmax(np.abs(dd[ok] / d[i] - 1))
    return (x >= 0) & (s_band <= 0.02) & (yb >= ys) & (np.nan_to_num(lad, nan=9) <= 0.02), s_band


def steady_mask(x, cand, series):
    z = np.load(series); steps = z["step"].astype(float); D = z["delta_r_raw"]
    if not np.allclose(z["x"], x):
        raise SystemExit("時系列の断面の x が抽出と違う")
    out = np.zeros(len(x), bool)
    for i in np.where(cand)[0]:
        v, det, _ = cq.classify_series(steps, D[:, i], 0.4, 0.005, 0.01, 20)
        m = re.search(r"最終比 ([+-]?[0-9.eE+-]+) %", det or "")
        a = float(m.group(1)) / 100 if m else 0.0
        out[i] = v == "STEADY" and abs(a) <= 0.01
    return out, len(steps)


def main():
    cal = json.loads((HERE / "kdelta_l100_calibration.json").read_text())
    ed = np.load(HERE / "dstar_l70_ns_run0020_constx_edge.npz", allow_pickle=True)
    co = np.load(HERE / "contur_l70_euler_run0023.npz", allow_pickle=True)
    x = ed["x"]; d = ed["delta_r_use"]
    stab, s_band = stable_mask(ed)
    steady, nsnap = steady_mask(x, stab, HERE / "dstar_series_l70_edge.npz")
    J = stab & steady
    C = np.interp(x, co["x"], co["delta_r_ann"])
    w = np.zeros_like(x); w[1:-1] = 0.5 * (x[2:] - x[:-2]); w[0] = 0.5 * (x[1] - x[0]); w[-1] = 0.5 * (x[-1] - x[-2])
    gate = {"n_x_ge0": int((x >= 0).sum()), "n_stable": int(stab.sum()), "n_judge": int(J.sum()), "n_snapshots": nsnap,
            "exit_in_judge": bool(J[-1]), "judge_x_min": float(x[J].min()) if J.any() else None}
    gate["ok"] = bool(gate["exit_in_judge"] and J.sum() >= 100 and nsnap >= 20)
    res = {"frozen": {"k_delta": cal["k_delta"], "k_delta_range": cal["k_delta_range"], "commit": "e69d56e7"}, "lines": LINES, "gate": gate}
    for name, kd in (("k_delta", cal["k_delta"]), ("k_lo", cal["k_delta_range"][0]), ("k_hi", cal["k_delta_range"][1])):
        e = kd * C / d - 1.0
        ef = (kd * C - d) / (np.maximum(np.nan_to_num(s_band, nan=1.0), 0.01) * d) * 0.01   # 床つき (参考)
        ex = float(e[-1]); rms = float(np.sqrt(np.average(e[J] ** 2, weights=w[J]))); loc = float(np.max(np.abs(e[J])))
        res[name] = {"k": kd, "exit_err": ex, "rms": rms, "local_max": loc, "x_at_local_max": float(x[J][np.argmax(np.abs(e[J]))]),
                     "pass": bool(abs(ex) <= LINES["exit"] and rms <= LINES["rms"] and loc <= LINES["local"]),
                     "median": float(np.median(e[J]))}
    regions = []
    for lo, hi in [(0, 1.17), (1.17, 2), (2, 4), (4, 8), (8, 16), (16, 99)]:
        sel = J & (x >= lo) & (x < hi)
        e = cal["k_delta"] * C / d - 1.0
        regions.append({"x": [lo, hi], "n": int(sel.sum()), "median": float(np.median(e[sel])) if sel.any() else None,
                        "maxabs": float(np.max(np.abs(e[sel]))) if sel.any() else None})
    res["regions"] = regions
    res["verdict"] = ("判定不能 (L70 の抽出・M3 のゲート)" if not gate["ok"] else
                      ("PASS (定数 1 個の移植性を支持)" if res["k_delta"]["pass"] and res["k_lo"]["pass"] and res["k_hi"]["pass"] else "FAIL (定数 1 個の移植性を棄却)"))
    (HERE / "l70_blind_judgement.json").write_text(json.dumps(res, indent=1, ensure_ascii=False))
    print(json.dumps(res, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
