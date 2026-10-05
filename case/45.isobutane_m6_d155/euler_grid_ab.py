"""plan tooling-nozzle-cfd-pinned-initial-line §5.1 #11f 手順 2 (事前登録 2026-10-06): Euler G0/G1 の格子 A/B。
usage: python3 euler_grid_ab.py prep RUN {G0|G1}   (準備 + 段階起動 soft → 本段 2 次 cfl 2・relax 0.7・12000 step)
       python3 euler_grid_ab.py eval RUN_G0 RUN_G1  (出口コア M: 自格子平均と G1 帯内 η 節点の共通標本、時系列 CSV)
IC = 旧較正 Euler 場 run_0086 (G0 は同格子なので restart_field、G1 は interp_field)。設計壁が旧較正と同じことを wall_design.csv で検査。"""
import csv, json, subprocess, sys
from pathlib import Path
import numpy as np
C = Path(__file__).resolve().parent; sys.path.insert(0, str(C.parents[1] / "design"))
from forge_design.evaluate import runner_axismach as RA  # noqa: E402
OLD = C / "run_0086_euler_wallfit_pincal_r1_ext6k"
if sys.argv[1] == "prep":
    rd, arm = C / sys.argv[2], sys.argv[3]
    info = RA.prepare(C / f"problem_d155_euler_pin_{arm}.yaml", rd, nsteps=12000, ic_from=(OLD if arm == "G1" else None), cfl_main=2.0, implicit_relax=0.7)
    (rd / "prepare_info.json").write_text(json.dumps(info, indent=1, default=str))
    if arm == "G0":
        last = sorted(OLD.glob("res_[0-9]*.h5"), key=lambda f: int(f.stem.split("_")[1]))[-1]
        r = subprocess.run([sys.executable, str(C.parents[1] / "solver_density_cuda/tools/restart_field.py"), str(last), str(rd / "nozzle.h5")], capture_output=True, text=True)
        print(r.stdout.strip().splitlines()[-1]); assert r.returncode == 0
    a = np.genfromtxt(rd / "wall_design.csv", delimiter=",", names=True); b = np.genfromtxt(OLD / "wall_design.csv", delimiter=",", names=True)
    k = a.dtype.names
    print("wall_design vs run_0086: max|Δr| =", float(np.max(np.abs(np.interp(b[k[0]], a[k[0]], a[k[1]]) - b[k[1]]))))
    rc = RA.run_staged(rd, cfl_main=2.0, mid_stage=False, stages="soft"); print("forge exit", rc); raise SystemExit(rc)
from forge_design.report.nozzle_report import _res_files, load_field  # noqa: E402
r0, r1 = C / sys.argv[2], C / sys.argv[3]


def prof(run, res):
    G = load_field(run, res); eta = G["R"][-1] / G["R"][-1, -1]; m = (eta >= 0.05) & (eta <= 0.7); o = np.argsort(eta)
    return eta, G["V"]["M"][-1], m, o


e1, _, m1, _ = prof(r1, _res_files(r1)[-1]); etaS = e1[m1]
out = {}
for run in (r0, r1):
    rows = []
    for res in _res_files(run):
        eta, M, m, o = prof(run, res)
        rows.append(dict(step=int(res.split("_")[1].split(".")[0]), M_own=float(M[m].mean()), M_common=float(np.interp(etaS, eta[o], M[o]).mean())))
    with open(run / "exitM_series.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    t, p5 = rows[-5:], rows[-10:-5]
    out[run.name] = {k: dict(mean=float(np.mean([r[k] for r in t])), ptp=float(np.ptp([r[k] for r in t])),
                             prev5_diff=float(np.mean([r[k] for r in t]) - np.mean([r[k] for r in p5])) if len(p5) == 5 else None) for k in ("M_own", "M_common")}
    q = subprocess.run([sys.executable, str(C.parents[1] / "solver_density_cuda/tools/check_quasisteady.py"), str(run), "--series-csv", str(run / "exitM_series.csv"),
                        "--series-cols", "M_own,M_common"], capture_output=True, text=True)
    out[run.name]["quasisteady"] = [l for l in q.stdout.splitlines() if "M_own" in l or "M_common" in l or ("===" in l and "csv" in l)]
A, B = out[r0.name]["M_common"], out[r1.name]["M_common"]
D = B["mean"] - A["mean"]; u = (A["ptp"] + B["ptp"]) / 2
out["D_common"] = dict(D=D, interval=[D - u, D + u], verdict=("A: 負側の格子依存 (< −2e-4)" if D + u < -2e-4 else "B: ±2e-4 内 (持越し主因説棄却)" if (D - u > -2e-4 and D + u < 2e-4) else "保留 (またぐ/正側)"))
(C / "_band_ab" / "euler_grid_ab.json").write_text(json.dumps(out, indent=1)); print(json.dumps(out, indent=1))
