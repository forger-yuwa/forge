#!/usr/bin/env python3
r"""③ベル (M) L100: δ_r の時系列 (M3) の判定 — plan tooling-nozzle-bell-two-mode §5.1 #7 の事前登録の線で。

入力: `dstar_series_l100_edge.npz` (B で `bell_deltastar.deltastar_series(..., band_select="edge")`、run_0019 の本段の 21 枚)。
登録の線 (codex 2026-10-12): `check_quasisteady` の classify_series に tail 0.4・drift 0.005・osc 0.01・min_snaps 20。対象は
代表 x (0.5・1・2・4・8・16)・出口・x 積分の L2 (Δx の台形の重み、x ≥ 0)・安定な全断面 (edge_stable_kdelta_provisional.json の定義を
最終場で再計算した集合)。全対象で STEADY。単調な系列は漸近値との差 (最終比) も 1 % 以内。
出力: dstar_series_l100_edge_m3.json。"""
import json
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "solver_density_cuda" / "tools"))
import check_quasisteady as cq  # noqa: E402

TAIL, DRIFT, OSC, MIN = 0.4, 0.005, 0.01, 20
ASYM_TOL = 0.01
ASYM_RE = re.compile(r"最終比 ([+-]?[0-9.eE+-]+) %")


def judge(steps, v):
    verdict, detail, _ = cq.classify_series(steps, v, TAIL, DRIFT, OSC, MIN)
    m = ASYM_RE.search(detail or "")
    asym = float(m.group(1)) / 100.0 if m else None
    ok = verdict == "STEADY" and (asym is None or abs(asym) <= ASYM_TOL)
    return {"verdict": verdict, "asym_final_ratio": asym, "ok": bool(ok), "detail": detail}


def main():
    z = np.load(HERE / "dstar_series_l100_edge.npz")
    steps, x, D = z["step"].astype(float), z["x"], z["delta_r_raw"]
    ed = np.load(HERE / "dstar_l100_ns_run0019_constx_edge.npz", allow_pickle=True)
    if not np.allclose(ed["x"], x):
        raise SystemExit("断面の x が抽出の最終場と違う")
    # 安定な断面 (edge_stable_kdelta_provisional.json と同じ定義)
    fac, rel = ed["yb_factors"], ed["rel_yb"]
    s_band = np.nanmax(np.abs(rel[:, [k for k, f in enumerate(fac) if f > 1]]), axis=1)
    d = ed["delta_r_use"]; yb, ys = ed["band_y_b"], ed["edge_y_star"]
    lad = np.full(len(x), np.nan)
    for i in range(len(x)):
        y, dd = ed["ladder_y"][i], ed["ladder_d"][i]
        ok = np.isfinite(dd) & (dd > 0) & (y >= yb[i]) & (y <= 2 * yb[i])
        if ok.sum() >= 3:
            lad[i] = np.nanmax(np.abs(dd[ok] / d[i] - 1))
    m0 = x >= 0
    stable = m0 & (s_band <= 0.02) & (yb >= ys) & (np.nan_to_num(lad, nan=9) <= 0.02)
    w = np.zeros_like(x); w[1:-1] = 0.5 * (x[2:] - x[:-2]); w[0] = 0.5 * (x[1] - x[0]); w[-1] = 0.5 * (x[-1] - x[-2])
    res = {"lines": {"tail": TAIL, "drift": DRIFT, "osc": OSC, "min_snaps": MIN, "asym_tol": ASYM_TOL},
           "n_snapshots": int(len(steps)), "steps": steps.tolist(), "n_stable": int(stable.sum())}
    targets = {}
    for s in (0.5, 1.0, 2.0, 4.0, 8.0, 16.0):
        i = int(np.argmin(abs(x - s)))
        targets[f"x{s:g}"] = (D[:, i], bool(stable[i]))
    targets["exit"] = (D[:, -1], bool(stable[-1]))
    targets["L2_x_ge0"] = (np.sqrt(np.sum(w[m0] * D[:, m0] ** 2, axis=1) / np.sum(w[m0])), True)
    targets["L2_stable"] = (np.sqrt(np.sum(w[stable] * D[:, stable] ** 2, axis=1) / np.sum(w[stable])), True)
    res["targets"] = {k: {**judge(steps, v), "stable_section": st} for k, (v, st) in targets.items()}
    per = []
    for i in np.where(stable)[0]:
        j = judge(steps, D[:, i]); per.append((float(x[i]), j["verdict"], j["asym_final_ratio"], j["ok"]))
    bad = [p for p in per if not p[3]]
    res["sections"] = {"n": len(per), "n_ok": len(per) - len(bad),
                       "verdicts": {v: sum(1 for p in per if p[1] == v) for v in sorted(set(p[1] for p in per))},
                       "not_ok_first": bad[:30],
                       "not_ok_x_ranges": _ranges([p[0] for p in bad], x)}
    # 最終 4 枚の相対変化 (参考)
    res["last4_rel_change_max_stable"] = float(np.nanmax(np.abs(D[-1, stable] / D[-4, stable] - 1)))
    res["all_ok"] = bool(all(t["ok"] for t in res["targets"].values()) and not bad)
    (HERE / "dstar_series_l100_edge_m3.json").write_text(json.dumps(res, indent=1, ensure_ascii=False, default=float))
    print(json.dumps({k: res[k] for k in ("n_snapshots", "n_stable", "all_ok", "last4_rel_change_max_stable")}, ensure_ascii=False))
    for k, t in res["targets"].items():
        print(f"{k:10s} {t['verdict']:12s} asym {t['asym_final_ratio']} ok {t['ok']}  | {t['detail'][:110]}")
    print("sections", res["sections"]["n"], "ok", res["sections"]["n_ok"], res["sections"]["verdicts"], "not ok ranges", res["sections"]["not_ok_x_ranges"][:10])


def _ranges(xs, xall):
    if not xs:
        return []
    xs = sorted(xs); out = [[xs[0], xs[0], 1]]
    dx = float(np.median(np.diff(xall)))
    for v in xs[1:]:
        if v - out[-1][1] <= 3 * dx:
            out[-1][1] = v; out[-1][2] += 1
        else:
            out.append([v, v, 1])
    return out


if __name__ == "__main__":
    main()
