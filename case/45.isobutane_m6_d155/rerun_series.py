"""plan tooling-rerun-conditions §6 (事前登録): 検証 run の量の時系列と準定常判定。
usage: python3 rerun_series.py RUN [EULER_REF] → RUN/quantities_series.csv (step, exitM_A, mdot[, delta_E]) + check_quasisteady (§6 の量ごとの閾値)
exitM_A = 出口断面の η∈[0.05,0.7] 節点平均 (nozzle_report と同じ)、mdot = 2π∫ρU_x r dr (r_t 単位) の x∈(−2.5, x_E) の中央値
(nozzle_report の流量比と同じ積分)、delta_E = E 法の未緩和 δ_E(x_F) (EULER_REF を与えたときだけ)。res_0 は含めない。"""
import csv, os, shutil, subprocess, sys
from pathlib import Path
import numpy as np
C = Path(__file__).resolve().parent; sys.path.insert(0, str(C.parents[1] / "design"))
from forge_design.report.nozzle_report import _res_files, load_field  # noqa: E402
from forge_design.metrics.deltastar import _load_structured  # noqa: E402
run = C / sys.argv[1]; eu = (C / sys.argv[2]) if len(sys.argv) > 2 else None
TOOLS = C.parents[1] / "solver_density_cuda/tools"


def mdot_of(rd, res):
    tmp = Path(f"/tmp/rs_{rd.name}_{res}"); shutil.rmtree(tmp, ignore_errors=True); tmp.mkdir()
    for f in ("bcondConfig.yaml", "solverConfig.yaml", "prepare_info.json"): shutil.copy(rd / f, tmp / f)
    os.symlink((rd / "nozzle.h5").resolve(), tmp / "nozzle.h5"); os.symlink((rd / res).resolve(), tmp / res)
    d = _load_structured(tmp); xs = d["x"][:, 0]; m = 2 * np.pi * np.trapezoid(d["q"] * d["r"], d["r"], axis=1)
    out = float(np.median(m[(xs > -2.5) & (xs < d["info"]["x_E"])]))
    if eu is not None:
        from forge_design.feedback.deltastar_loop import extract_and_merge, read_delta_r_next
        extract_and_merge(tmp, eu, band_select="edge"); nx = read_delta_r_next(tmp / "delta_r_next.csv")
        G = load_field(rd, res); xF = float(G["X"][-1, 0]); dE = float(np.interp(xF, nx["x_rt"], nx["delta_E"]))
    else:
        dE = None
    shutil.rmtree(tmp); return out, dE


rows = []
for res in _res_files(run):
    G = load_field(run, res); eta = G["R"][-1] / G["R"][-1, -1]; core = (eta >= 0.05) & (eta <= 0.7)
    md, dE = mdot_of(run, res)
    r = dict(step=int(res.split("_")[1].split(".")[0]), exitM_A=float(G["V"]["M"][-1][core].mean()), mdot=md)
    if dE is not None: r["delta_E"] = dE
    rows.append(r)
with open(run / "quantities_series.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
t, p = rows[-5:], rows[-10:-5]
for k in rows[0]:
    if k == "step": continue
    a = np.array([r[k] for r in t]); b = np.array([r[k] for r in p]) if len(p) == 5 else None
    print(f"{k}: tail5 mean {a.mean():.7g} ptp {np.ptp(a):.3g} prev5diff {(a.mean() - b.mean()) if b is not None else float('nan'):.3g}")
thr = {"delta_E": ("5e-5", "1e-4"), "exitM_A": ("3e-6", "3e-6"), "mdot": ("1e-5", "2e-5")}
for k, (dr, osc) in thr.items():
    if k not in rows[0]: continue
    q = subprocess.run([sys.executable, str(TOOLS / "check_quasisteady.py"), str(run), "--series-csv", str(run / "quantities_series.csv"), "--series-cols", k,
                        "--tail", "0.4", "--min-snaps", "10", "--drift", dr, "--osc", osc], capture_output=True, text=True)
    print("\n".join(l for l in q.stdout.splitlines() if k in l))
