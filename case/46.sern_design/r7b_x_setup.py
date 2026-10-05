#!/usr/bin/env python3
"""R7b 4-1 (plan tooling-nozzle-sern-chain §5.1 R7b ④、2026-10-05 事前登録): `L_sw_exact` 形状の m6_on・g3 対を作る (AWS で実行)。
  A = run_1050_r7b_x_lsw08 (problem_3d_prod_m6on_wallres_lswx08.yaml、初期場 run_1034 [A・g3] から節点番号で)
  B = run_1051_r7b_x_lsw10 (problem_3d_prod_m6on_wallres_lswx10.yaml、初期場 run_1040 [B・g3] から節点番号で)
格子と bcond は runner の prepare、solverConfig は元 run の時間・space を継承し physProp・turbulence だけ現行 YAML の生成値 (R7a と同じ)。
"""
import os, re, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "design"))
from forge_design.evaluate import runner_sern as R2, runner_sern3d as R3  # noqa: E402
from forge_design.probdef import load_problem  # noqa: E402

for dst, prob, src in (("run_1050_r7b_x_lsw08", "problem_3d_prod_m6on_wallres_lswx08.yaml", "run_1034_r7a_lsw08_1"),
                       ("run_1051_r7b_x_lsw10", "problem_3d_prod_m6on_wallres_lswx10.yaml", "run_1040_r7a_lsw10_1_c40k")):
    d = os.path.join(HERE, dst); s_ = os.path.join(HERE, src)
    r = subprocess.run([sys.executable, "-m", "forge_design.evaluate.runner_sern3d", os.path.join(HERE, prob), d, "--prepare-only"],
                       cwd=os.path.join(ROOT, "design"), capture_output=True, text=True)
    open(os.path.join(d, "PREPARE.log"), "w").write(r.stdout + r.stderr)
    if r.returncode or not os.path.exists(os.path.join(d, "sern.h5")):
        raise SystemExit(f"prepare 失敗 {dst}: {r.stderr[-600:]}")
    p = load_problem(os.path.join(HERE, prob)); R2.design_snapshot(p); R2.select_operating_point(p, None); st = R2.gas_states(p)
    gen = R3._solver_config(p, 20000, 500, 0.25, float(st["ext"]["P"]))
    cfg = open(os.path.join(s_, "solverConfig.yaml")).read()
    for k in ("physProp:", "turbulence:"):
        cfg = re.sub(rf"^{k}.*$", [l for l in gen.splitlines() if l.startswith(k)][0], cfg, count=1, flags=re.M)
    open(os.path.join(d, "solverConfig.yaml"), "w").write(cfg)
    r = subprocess.run([sys.executable, os.path.join(HERE, "r7b_index_restart.py"), os.path.join(s_, "res_20000.h5"), os.path.join(s_, "sern.h5"),
                        os.path.join(d, "sern.h5")], capture_output=True, text=True)
    open(os.path.join(d, "RESTART_FROM.txt"), "w").write(r.stdout + r.stderr)
    if r.returncode:
        raise SystemExit(f"r7b_index_restart 失敗 {dst}: {r.stdout[-600:]}{r.stderr[-600:]}")
    open(os.path.join(d, "IC_FROM.txt"), "w").write(f"R7b 4-1: r7b_index_restart.py {src} res_20000 -> sern.h5 (節点番号で、タグ・接続同一を検査)。config = {prob} 生成 + {src} の時間・space\n")
    print(dst, r.stdout.strip().splitlines()[0][:200])
