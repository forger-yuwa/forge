#!/usr/bin/env python3
"""③ベル (M) モード: 円弧 kernel の TP の Rao 輪郭を L_k 5 本で解く (記録、合否の判定はしない) —
plan tooling-nozzle-bell-two-mode §5.1 #5 の (d) (2026-10-11 diagnostician の判断
notes/reviews/2026-10-11-bell-s1-generator-holes-diagnostician.md の「次の一手」)。

生成器 `geometry/rao_axisym.py` の `design_rao` (kernel = "arc") を、§5.1 #5 (a)(b) の修正 (DE を r パラメータに、
中の根の枝を η_D の連続性で追う、受理条件 θ_E ≤ θ_B) の後で回す。

設定: ガス = 基準ケースの凍結 TP (`RocketExhaust.moc_gas`、../exhaust_lox_lh2_pc100_of6.yaml)、R_t 1、ε 50 (固定)、
初期線 = Hall の n 161、kernel = 円弧 (曲率半径 R_t)、L_k = 等面積比 15° 円錐の長さ (√ε − 1)/tan15° = 22.6575 r_t
(`wall_bell_c2.cone15_length`、plan §8-1) の 60・70・80・90・100 %。DE の r 格子・Goursat の点数・走査点数は `design_rao` の
既定 (n 161 で dr 3.2/n、m = q = n/2、θ_B の走査 `N_SCAN_SB` 点、η_D の走査 `N_SCAN_D` 点)。

出力 (このディレクトリ):
- `rao_lk_n161.json`: L_k ごとの記録 (根の有無と理由、θ_B・η_D・θ_E・θ_max・x_E・r_E・κ_W(B⁺) − κ_B・推力の 2 会計の差、
  受理した候補・捨てた候補 [θ_E > θ_B]・lost、θ_B の走査の全根、枝) と来歴 (commit・ソースの sha256)。
- `rao_lk_n161.csv`: 同じものの 1 行 1 L_k の表。
- `rao_l60_scan.csv`: L_k 60 % の θ_B 36〜46° を細かく (`N_SCAN_SB` の倍の点数で) 走査した D の根の表
  (θ_B, η_D, x_E, θ_E と r_E・det_D・枝)。根の無い走査点は理由つきの 1 行。

実行: design/.venv-opt の Python で `python case/68.bell_two_mode/rao_m/run_rao_lk.py` (L_k 5 本と 60 % の走査を並列、
数十分)。`--lk 0.7,0.8` で一部だけ、`--no-scan60` で走査を省く (一部だけ回したときは出力を `*_partial` に書く)。
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
import sys
import time
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "design"))
from forge_design.gas.rocket_exhaust import RocketExhaust  # noqa: E402
from forge_design.geometry import rao_axisym as ra  # noqa: E402
from forge_design.geometry.wall_bell_c2 import cone15_length  # noqa: E402

EXHAUST = HERE.parent / "exhaust_lox_lh2_pc100_of6.yaml"
EPS = 50.0
R_T = 1.0
N_LINE = 161
L_FRACS = (0.6, 0.7, 0.8, 0.9, 1.0)
SCAN60 = {"L_frac": 0.6, "theta_B_deg": (36.0, 46.0), "n_scan": 2 * ra.N_SCAN_SB}
SOURCES = ("design/forge_design/geometry/rao_axisym.py", "design/forge_design/gas/isentrope.py",
           "design/forge_design/geometry/moc_forward.py", "design/forge_design/geometry/moc_kernel.py")


def git_info() -> dict:
    head = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "-C", str(ROOT), "status", "--porcelain", "design/forge_design", "design/tests"],
                           capture_output=True, text=True).stdout
    sha = {s: hashlib.sha256((ROOT / s).read_bytes()).hexdigest()[:16] for s in SOURCES}
    return {"commit": head, "design_tree_dirty": [ln for ln in dirty.splitlines() if ln.strip()], "source_sha256_16": sha}


def _gas_and_line():
    ig = ra.as_isentrope(RocketExhaust.from_file(EXHAUST).moc_gas())
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        line = ra._hall_line(R_T, ig, N_LINE)
    return ig, line


def _jsonable(v):
    if isinstance(v, dict):
        return {str(k): _jsonable(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [_jsonable(x) for x in v]
    if isinstance(v, (np.floating, float)):
        return None if not np.isfinite(v) else float(v)
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, np.ndarray):
        return _jsonable(v.tolist())
    return v


def run_lk(frac: float) -> dict:
    """L_k 1 本: `design_rao` (円弧 kernel) の結果から表に載せる量を抜く。"""
    t0 = time.time()
    ig, line = _gas_and_line()
    L = frac * cone15_length(EPS)
    res = ra.design_rao(ig, R_t=R_T, L=L, eps=EPS, kernel="arc", line=line)
    rec = {"L_frac": frac, "L_rt": L, "status": res["status"], "why": res.get("why"), "n": res["n"], "dr": res["dr"],
           "m": res["m"], "q": res["q"], "eta_tol": res.get("eta_tol"),
           "rejected": res.get("rejected", []), "lost": res.get("lost", []), "branches": res.get("branches", {}),
           "scan_mid": res.get("scan_mid", [])}
    if res["status"] == "ok":
        E, th = res["E"], res["thrust"]
        rec.update({
            "theta_B_deg": float(np.rad2deg(res["theta_B"])), "s_B": res["s_B"], "eta_D": res["eta_D"],
            "eta0_BD": res["eta0_BD"], "N_BD": res["N_BD"], "det_D": res["det_D"], "branch": res["branch"],
            "M_D": float(res["D"]["M"]), "theta_D_deg": float(np.rad2deg(res["D"]["th"])),
            "theta_E_deg": float(np.rad2deg(res["theta_E"])), "theta_max_deg": float(np.rad2deg(res["theta_max"])),
            "x_E": float(E["x"]), "r_E": float(E["r"]), "M_E": float(E["M"]), "x_B": res["x_B"], "r_B": res["r_B"],
            "kappa_B": res["kappa_B"], "kappaW_B": res["kappaW_B"], "kappaW_B_minus_kappa_B": res["kappaW_B"] - res["kappa_B"],
            "thrust_rel_diff": th["rel_diff"], "mdot_rel_diff": th["mdot_rel_diff"], "F_cs": th["F_cs"], "F_wall": th["F_wall"],
            "closure": res["closure"], "n_wall_missing": res["n_wall_missing"], "goursat_nfail": res["goursat_nfail"],
            "streamline_resid_max": res["streamline_resid_max"], "CF": res["CF"], "isp_ratio": res["isp_ratio"],
            "edge": res["edge"], "candidates": res["candidates"],
        })
    rec["time_s"] = time.time() - t0
    return _jsonable(rec)


def run_scan60() -> dict:
    """L_k 60 % の θ_B 36〜46° を細かく走査 (`solve_thetaB` に同じ範囲・倍の点数を渡す。kernel 場は design_rao と同じ引数)。"""
    t0 = time.time()
    ig, line = _gas_and_line()
    L = SCAN60["L_frac"] * cone15_length(EPS)
    rng = ra.sB_range(R_T, 0.0)
    net = ra.kernel_field(ig, R_T, 0.0, line, rng[1], L_march=ra.default_L_march(L))
    lo, hi = (np.deg2rad(v) * R_T for v in SCAN60["theta_B_deg"])           # 円弧: s_B = θ_B R_t
    hi = min(hi, net.s_need)
    mid = ra.solve_thetaB(net, L, eps=EPS, dr=3.2 / len(line["x"]), s_range=(lo, hi), n_scan=SCAN60["n_scan"])
    return _jsonable({"L_frac": SCAN60["L_frac"], "L_rt": L, "theta_B_range_deg": [float(np.rad2deg(lo)), float(np.rad2deg(hi))],
                      "n_scan": SCAN60["n_scan"], "scan": mid["scan"], "branches": mid["branches"], "eta_tol": mid["eta_tol"],
                      "roots": [ra.candidate_summary(c) for c in mid["roots"]],
                      "rejected": [ra.candidate_summary(c) for c in mid["rejected"]], "lost": mid["lost"], "why": mid["why"],
                      "time_s": time.time() - t0})


CSV_COLS = ("L_frac", "L_rt", "status", "theta_B_deg", "eta_D", "eta0_BD", "N_BD", "det_D", "theta_E_deg", "theta_max_deg",
            "x_E", "r_E", "M_E", "kappaW_B_minus_kappa_B", "thrust_rel_diff", "mdot_rel_diff", "closure_dx", "closure_dr",
            "n_wall_missing", "goursat_nfail", "n_candidates", "n_rejected", "rejected", "n_lost", "time_s", "why")


def lk_row(r: dict) -> dict:
    row = {k: r.get(k) for k in CSV_COLS}
    if r.get("closure"):
        row["closure_dx"], row["closure_dr"] = r["closure"]["dx"], r["closure"]["dr"]
    row["n_candidates"] = len(r.get("candidates") or [])
    row["n_rejected"] = len(r.get("rejected") or [])
    row["rejected"] = "; ".join(f"θ_B {c['theta_B_deg']:.3f}°・η_D {c['eta_D']:.2f}・x_E {c['x_E']:.3f}・θ_E {c['theta_E_deg']:.3f}°"
                                for c in r.get("rejected") or [])
    row["n_lost"] = len(r.get("lost") or [])
    return row


def scan_rows(sc: dict) -> list:
    rows = []
    for p in sc["scan"]:
        if not p["roots"]:
            rows.append({"theta_B_deg": p["theta_B_deg"], "s_B": p["s_B"], "eta0": p["eta0"], "N": p["N"], "eta_D": None,
                         "x_E": None, "r_E": None, "theta_E_deg": None, "det_D": None, "branch": None,
                         "theta_E_le_theta_B": None, "why": p["why"]})
        for q in p["roots"]:
            rows.append({"theta_B_deg": p["theta_B_deg"], "s_B": p["s_B"], "eta0": p["eta0"], "N": p["N"], "eta_D": q["eta_D"],
                         "x_E": q["x_E"], "r_E": q["r_E"], "theta_E_deg": q["theta_E_deg"], "det_D": q["det_D"],
                         "branch": q.get("branch"), "theta_E_le_theta_B": q["theta_E_deg"] <= p["theta_B_deg"], "why": None})
    return rows


def write_csv(path: Path, rows: list, cols) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(cols))
        w.writeheader()
        for r in rows:
            w.writerow({k: ("" if r.get(k) is None else (f"{r[k]:.10g}" if isinstance(r[k], float) else r[k])) for k in cols})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lk", default=",".join(str(f) for f in L_FRACS))
    ap.add_argument("--no-scan60", action="store_true")
    ap.add_argument("--workers", type=int, default=6)
    a = ap.parse_args()
    fracs = tuple(float(v) for v in a.lk.split(","))
    full = fracs == L_FRACS and not a.no_scan60
    suffix = "" if full else "_partial"
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        futs = {f: ex.submit(run_lk, f) for f in fracs}
        fut_scan = None if a.no_scan60 else ex.submit(run_scan60)
        lks = []
        for f in fracs:
            r = futs[f].result()
            print(f"L {f:.0%} ({r['L_rt']:.3f} r_t): {r['status']}"
                  + (f" θ_B {r['theta_B_deg']:.3f}° η_D {r['eta_D']:.2f} θ_E {r['theta_E_deg']:.3f}° θ_max {r['theta_max_deg']:.3f}° "
                     f"x_E {r['x_E']:.4f} r_E {r['r_E']:.4f} κ_W(B⁺)−κ_B {r['kappaW_B_minus_kappa_B']:+.4f} "
                     f"推力の 2 会計の差 {r['thrust_rel_diff']:+.2e}" if r["status"] == "ok" else f" — {r['why']}")
                  + f" (捨てた候補 {len(r['rejected'])}、{r['time_s']:.0f} s)", flush=True)
            lks.append(r)
        sc = fut_scan.result() if fut_scan is not None else None
    out = {"plan": "plans/active/tooling-nozzle-bell-two-mode.md §5.1 #5 (d)",
           "decision_record": "notes/reviews/2026-10-11-bell-s1-generator-holes-diagnostician.md",
           "setup": {"gas": "EXH (RocketExhaust.moc_gas)", "exhaust": str(EXHAUST.relative_to(ROOT)), "R_t": R_T, "eps": EPS,
                     "initial_line": f"Hall n {N_LINE}", "kernel": "arc", "L_cone15_rt": cone15_length(EPS),
                     "N_SCAN_SB": ra.N_SCAN_SB, "N_SCAN_D": ra.N_SCAN_D, "ETA_LINK_FRAC": ra.ETA_LINK_FRAC,
                     "THETA_B_RANGE_DEG": ra.THETA_B_RANGE_DEG},
           **git_info(), "L_k": lks, "total_s": time.time() - t0}
    if sc is not None:
        out["scan60_summary"] = {k: sc[k] for k in ("theta_B_range_deg", "n_scan", "roots", "rejected", "lost", "why", "time_s")}
    (HERE / f"rao_lk_n161{suffix}.json").write_text(json.dumps(_jsonable(out), ensure_ascii=False, indent=1))
    write_csv(HERE / f"rao_lk_n161{suffix}.csv", [lk_row(r) for r in lks], CSV_COLS)
    if sc is not None:
        write_csv(HERE / f"rao_l60_scan{suffix}.csv", scan_rows(sc),
                  ("theta_B_deg", "s_B", "eta0", "N", "eta_D", "x_E", "r_E", "theta_E_deg", "det_D", "branch",
                   "theta_E_le_theta_B", "why"))
        print(f"60 % の走査 (θ_B {sc['theta_B_range_deg'][0]:.2f}〜{sc['theta_B_range_deg'][1]:.2f}°、{sc['n_scan']} 点): "
              f"受理 {len(sc['roots'])}・捨てた {len(sc['rejected'])}・lost {len(sc['lost'])} ({sc['time_s']:.0f} s)", flush=True)
    print(f"total {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
