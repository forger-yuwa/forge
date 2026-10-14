#!/usr/bin/env python3
r"""③ベル (M): 同じ軸長で出口の面積比 ε を広げたときの、Rao (MOC、非粘性) と 15° 相当の円錐 (同じ長さ・同じ ε) の推力の比 —
研究目標「円錐より 2 %」(ユーザ 2026-10-15) がどこで超えるかの探索。plan tooling-nozzle-bell-two-mode §5.1 #8d。

Rao: `design_rao` (円弧 kernel、CFD ピンの初期線 = 凍結源 run_0015 の res_30000、n 321、`run_pin_line.py` と同じ)、固定 ε と固定 L。
η_Rao = F_cs / (π·F_1D,nd(ε)) (同じ MOC の気体の 1 次元理想推力、真空)。
円錐: 同じ L・ε の直線壁の半角 α = atan((√ε − 1)/L) (スロートの円弧は無視)、発散係数 λ = (1 + cos α)/2 を運動量の推力にだけ掛ける
(圧力の推力 p_e A_e は損なし)、η_cone = (λ·ṁu_e + p_e A_e)/(ṁu_e + p_e A_e)、出口の状態は同じ MOC の気体の 1 次元等エントロピー。
推力の比 = (η_Rao/MDOT_LINE)/η_cone − 1。η_Rao はスロートの流量の不足 (ピン線の ṁ/ṁ_1D = MDOT_LINE、どの ε でも同じスロート) を含み、
円錐の λ のモデルは含まないので、Rao を流量あたりの推力に直してそろえる (2026-10-15 訂正: 当初は η_Rao/η_cone − 1 で、上乗せを
約 0.74 ポイント過小に出していた)。同じスロートの円錐も同じ流量の不足を受けるので、これが同じ長さ・同じ出口径の推力の比。
粘性の損は両者で同程度と仮定 (含めない)。`--post` は既存の rao_eps_scan.json から csv をこの式で作り直す (生成器は回さない)。
B で並列 4 本 (1 本 25 分): `python run_rao_eps_scan.py [--cases L100:80 L100:100 ...] [--workers 4]`。出力 rao_eps_scan.json・.csv。
"""
from __future__ import annotations

import os

for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_k, "1")

import argparse  # noqa: E402
import csv  # noqa: E402
import json  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from concurrent.futures import ProcessPoolExecutor, as_completed  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2] / "design"))
sys.path.insert(0, str(HERE))
import run_pin_line as pl  # noqa: E402
from forge_design.feedback.cfd_initial_line import moc_initial_line  # noqa: E402
from forge_design.geometry import rao_axisym as ra  # noqa: E402

MDOT_LINE = 0.992622     # ピン線の ṁ/ṁ_1D (rao_pin_line.json、run_0018 の V2 の参照値)
DEFAULT = ("L100:80", "L100:100", "L100:120", "L100:150", "L80:80", "L80:100")


def cone(ig, eps: float, L: float) -> dict:
    i1 = ig.ideal_1d(eps, 0.0)
    a = float(np.arctan((np.sqrt(eps) - 1.0) / L))
    lam = 0.5 * (1.0 + np.cos(a))
    out = {"alpha_deg": float(np.degrees(a)), "lambda": float(lam), "ideal_keys": sorted(i1.keys())}
    if all(k in i1 for k in ("F_mom", "F_p")):
        Fm, Fp = float(i1["F_mom"]), float(i1["F_p"])
    else:
        # 1 次元の出口状態から: ṁu_e = F − p_e·ε (A* = 1)、p_e は等エントロピー
        F = float(i1["F"])
        pe = float(i1["p_e"]) if "p_e" in i1 else None
        if pe is None:
            raise KeyError(f"ideal_1d に F_mom/F_p も p_e も無い ({sorted(i1)})")
        Fp = pe * eps
        Fm = F - Fp
    out.update(F_mom=Fm, F_p=Fp, eta_cone=float((lam * Fm + Fp) / (Fm + Fp)))
    return out


def job(spec: str) -> dict:
    t0 = time.time()
    lk, e = spec.split(":")
    frac, eps = float(lk[1:]) / 100.0, float(e)
    L = frac * pl.L100
    ig = pl.gas()
    th, _ = pl.pinned_throat(pl.FREEZE, ig)
    line = moc_initial_line(th, ig.moc, pl.N_LINE)
    res = ra.design_rao(ig, R_t=pl.R_T, L=L, eps=eps, kernel="arc", line=line, p_a=pl.P_A)
    rec = {"case": spec, "L_frac": frac, "L_rt": L, "eps": eps, "status": res["status"], "why": res.get("why"), "time_s": None}
    c = cone(ig, eps, L)
    rec["cone"] = c
    if res["status"] == "ok":
        F1D = float(ig.ideal_1d(eps, 0.0)["F"])
        eta = float(res["thrust"]["F_cs"]) / (np.pi * F1D)
        rec.update(eta_rao=eta, theta_B_deg=float(np.degrees(res["theta_B"])), theta_E_deg=float(np.degrees(res["theta_E"])),
                   theta_max_deg=float(np.degrees(res["theta_max"])), thrust_rel_diff=res["thrust"]["rel_diff"],
                   eta_rao_specific=eta / MDOT_LINE, gain_vs_cone=(eta / MDOT_LINE) / c["eta_cone"] - 1.0)
    rec["time_s"] = time.time() - t0
    return rec


def write_csv(out: dict) -> None:
    L100 = pl.L100
    with open(HERE / "rao_eps_scan.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["case", "L_frac", "eps", "L_over_L15cone", "status", "eta_rao", "eta_rao_specific", "theta_E_deg", "cone_alpha_deg",
                    "eta_cone", "gain_vs_cone"])
        for r in sorted(out["rows"], key=lambda r: (r.get("eps", 0), -r.get("L_frac", 0))):
            c = r.get("cone", {})
            e = r.get("eta_rao")
            es = e / MDOT_LINE if e else None
            L15 = (np.sqrt(r["eps"]) - 1.0) / np.tan(np.radians(15.0))
            w.writerow([r["case"], r.get("L_frac"), r.get("eps"), r["L_frac"] * L100 / L15, r["status"], e, es, r.get("theta_E_deg"),
                        c.get("alpha_deg"), c.get("eta_cone"), (es / c["eta_cone"] - 1.0) if es else None])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", nargs="*", default=list(DEFAULT))
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--post", action="store_true", help="生成器を回さず rao_eps_scan.json から csv を作り直す")
    a = ap.parse_args()
    if a.post:
        out = json.loads((HERE / "rao_eps_scan.json").read_text())
        out["note_2026_10_15"] = "gain_vs_cone は csv の式 (流量あたりの推力でそろえた値) を使う — json の旧 gain_vs_cone は η_Rao/η_cone − 1 (過小)"
        (HERE / "rao_eps_scan.json").write_text(json.dumps(out, ensure_ascii=False, indent=1, default=float))
        write_csv(out)
        print("post done")
        return
    out = {"plan": "plans/active/tooling-nozzle-bell-two-mode.md §5.1 #8d", "freeze": pl.FREEZE, "n_line": pl.N_LINE, "rows": []}
    # ε 50 の参照 (既存の生成器の記録から — 同じ式で円錐も)
    ig = pl.gas()
    for k, f in (("L100", 1.0), ("L90", 0.9), ("L80", 0.8), ("L70", 0.7)):
        rp = HERE / "rao_pin_parts" / f"{k}_pin_n321.json"
        if rp.exists():
            r = json.loads(rp.read_text())
            eta = float(r["F_cs"]) / (np.pi * float(ig.ideal_1d(50.0, 0.0)["F"]))
            c = cone(ig, 50.0, f * pl.L100)
            out["rows"].append({"case": f"{k}:50", "L_frac": f, "eps": 50.0, "status": "ok (既存)", "eta_rao": eta,
                                "eta_rao_specific": eta / MDOT_LINE, "theta_E_deg": r.get("theta_E_deg"), "cone": c,
                                "gain_vs_cone": (eta / MDOT_LINE) / c["eta_cone"] - 1.0})
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        futs = {ex.submit(job, s): s for s in a.cases}
        for fu in as_completed(futs):
            try:
                r = fu.result()
            except Exception as e:  # noqa: BLE001
                r = {"case": futs[fu], "status": "error", "why": f"{type(e).__name__}: {e}"}
            out["rows"].append(r)
            print(json.dumps({k: r.get(k) for k in ("case", "status", "eta_rao", "theta_E_deg", "gain_vs_cone", "time_s")}, ensure_ascii=False,
                             default=float), "cone", json.dumps({k: r.get("cone", {}).get(k) for k in ("alpha_deg", "eta_cone")}), flush=True)
            (HERE / "rao_eps_scan.json").write_text(json.dumps(out, ensure_ascii=False, indent=1, default=float))
    write_csv(out)
    print("done")


if __name__ == "__main__":
    main()
