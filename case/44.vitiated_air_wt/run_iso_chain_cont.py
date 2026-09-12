#!/usr/bin/env python3
"""等温壁チェーンの継続 run (plan tooling-nozzle-isothermal-wall-chain §4.5):
  run_0114_va_ns_fine_iso300_samewall_cont : run_0109 (S4a-B) を同一壁・同一メッシュで +36000 step (熱場の静定: Q_w drift 9.6 % → 1 % 目標)
  run_0115_va_ns_fine_iso300_ib_pass1      : S4b pass 1 (run_0110 の場から帯局所抽出 → 符号付き δ_r, ω=1.0, warm start)
usage: design/.venv-opt/bin/python run_iso_chain_cont.py [--wait-for FILE:TOKEN]
"""
import argparse, json, sys, time
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]; sys.path.insert(0, str(ROOT / "design"))
CASE = Path(__file__).resolve().parent
EULER = CASE / "run_0005_va_R2_LU6_Lc8"; REF = CASE / "run_0107_va_R2_LU6_Lc8_ns_ib_pass0"
P_ISO = CASE / "problem_va_R2_LU6_Lc8_ns_fine_iso300.yaml"
RB = CASE / "run_0109_va_ns_fine_iso300_samewall"; RBC = CASE / "run_0114_va_ns_fine_iso300_samewall_cont"
RC = CASE / "run_0110_va_ns_fine_iso300_ib_pass0"; RC1 = CASE / "run_0115_va_ns_fine_iso300_ib_pass1"


def log(m): print(time.strftime("%H:%M:%S"), m, flush=True)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--wait-for", default=None); ap.add_argument("--steps", type=int, default=36000)
    a = ap.parse_args()
    if a.wait_for:
        f, tok = a.wait_for.split(":")
        while not (Path(f).exists() and tok in Path(f).read_text()): time.sleep(30)
        log("GPU free")
    from forge_design.evaluate.runner_axismach import prepare_ns, run_staged_ns, collect
    info = prepare_ns(P_ISO, RBC, nsteps=a.steps, ic_from=RB, delta_r_csv=str(REF / "delta_r_initial.csv"), offset="radial", euler_ref=EULER)
    (RBC / "prepare_info.json").write_text(json.dumps(info, indent=1, default=str)); log(f"prepared {RBC.name}")
    rc = run_staged_ns(RBC, stages="none"); log(f"{RBC.name} rc={rc}")
    m = collect(P_ISO, RBC); (RBC / "metrics.json").write_text(json.dumps(m, indent=1, default=str))
    from forge_design.feedback.deltastar_loop import run_pass
    run_pass(str(P_ISO), str(EULER), str(RC), str(RC1), omega=1.0, ic_from=str(RC), nsteps=a.steps, stages="none")
    log("pass1 done")


if __name__ == "__main__":
    main()
