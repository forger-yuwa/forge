#!/usr/bin/env python3
"""case/44 等温壁チェーン (plan tooling-nozzle-isothermal-wall-chain §4.5) の実行ドライバ。

S4a-A: run_0108 断熱・細分メッシュ・run_0107 と同じ物理壁 (delta_r_csv = run_0107/delta_r_initial.csv)、IC = run_0107 (cross-mesh)
S4a-B: run_0109 等温 300 K・同じ壁・同じメッシュ、IC = run_0108 (index コピー + soft/mid 段)
S4b  : run_0110 等温再設計 pass 0 (積分法初期壁 Tw 300, 符号付き δ*)、IC = run_0109
usage: design/.venv-opt/bin/python run_iso_chain.py [--stage a|b|c|all] [--wait-for FILE:TOKEN]
"""
import argparse, json, shutil, subprocess, sys, time
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]; sys.path.insert(0, str(ROOT / "design"))
CASE = Path(__file__).resolve().parent
EULER = CASE / "run_0005_va_R2_LU6_Lc8"
REF = CASE / "run_0107_va_R2_LU6_Lc8_ns_ib_pass0"
P_AD = CASE / "problem_va_R2_LU6_Lc8_ns_fine.yaml"
P_ISO = CASE / "problem_va_R2_LU6_Lc8_ns_fine_iso300.yaml"
RA = CASE / "run_0108_va_ns_fine_ad_samewall"; RB = CASE / "run_0109_va_ns_fine_iso300_samewall"; RC = CASE / "run_0110_va_ns_fine_iso300_ib_pass0"


def log(msg):
    print(time.strftime("%H:%M:%S"), msg, flush=True)


def run_ns(problem, run_dir, ic_from, stages="full", **kw):
    from forge_design.evaluate.runner_axismach import prepare_ns, run_staged_ns, collect
    info = prepare_ns(problem, run_dir, ic_from=ic_from, euler_ref=EULER, **kw)
    (run_dir / "prepare_info.json").write_text(json.dumps(info, indent=1, default=str))
    log(f"prepared {run_dir.name}: mesh {info['mesh']} wall_thermal {info.get('wall_thermal')} dstar_src {info['dstar_source'][:60]}")
    log(Path(run_dir, "MESH_QUALITY.txt").read_text().splitlines()[-1])
    t0 = time.time(); rc = run_staged_ns(run_dir, stages=stages)
    log(f"{run_dir.name} rc={rc} {time.time()-t0:.0f}s")
    m = collect(problem, run_dir); (run_dir / "metrics.json").write_text(json.dumps(m, indent=1, default=str))
    log(f"{run_dir.name} verdict: {m.get('convergence_verdict')}")
    return rc


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--stage", default="all"); ap.add_argument("--wait-for", default=None)
    a = ap.parse_args()
    if a.wait_for:
        f, tok = a.wait_for.split(":")
        while not (Path(f).exists() and tok in Path(f).read_text()):
            time.sleep(30)
        log("GPU queue free")
    if a.stage in ("a", "all"):
        run_ns(P_AD, RA, ic_from=REF, delta_r_csv=str(REF / "delta_r_initial.csv"), offset="radial")
    if a.stage in ("b", "all"):
        run_ns(P_ISO, RB, ic_from=RA, delta_r_csv=str(REF / "delta_r_initial.csv"), offset="radial")
    if a.stage in ("c", "all"):
        from forge_design.feedback.deltastar_loop import run_pass0_integral
        run_pass0_integral(str(P_ISO), str(EULER), str(RC), str(RB), initializer={"model": "contur"}, stages="full")


if __name__ == "__main__":
    main()
