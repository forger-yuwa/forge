#!/usr/bin/env python3
"""R7a 判別試験 (plan tooling-nozzle-sern-chain §5.1 R7a、2026-10-05 事前登録) の run を作る (AWS で実行)。
  A = L_sw 0.8 (生産 problem_3d_prod_m6on_wallres.yaml、run_1017 と同一格子) × 3: run_1034–1036
  B = L_sw 1.0 (problem_3d_prod_m6on_wallres_lsw10.yaml、格子を作り直す) × 3: run_1037–1039
初期場は run_1017_ff4f_g3_2p50 の最終場 (A: restart_field、B: interp_field)。時間・space の行は run_1017 を引き継ぎ、
physProp・turbulence は現行 YAML から生成 (prod_restart_setup と同じ)。各群 1 本目を作り、2・3 本目はその入力の複製 (反復は同一入力)。
  python3 r7a_setup.py
"""
import json, os, re, shutil, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "design"))
from forge_design.evaluate import runner_sern as R2, runner_sern3d as R3  # noqa: E402
from forge_design.probdef import load_problem  # noqa: E402

SRC = os.path.join(HERE, os.environ.get("R7A_SRC", "run_1017_ff4f_g3_2p50")); SRC_RES = os.path.join(SRC, "res_20000.h5")
TOOLS = os.path.join(ROOT, "solver_density_cuda", "tools")
GROUPS = {"A": ("problem_3d_prod_m6on_wallres.yaml", ["run_1034_r7a_lsw08_1", "run_1035_r7a_lsw08_2", "run_1036_r7a_lsw08_3"]),
          "B": ("problem_3d_prod_m6on_wallres_lsw10.yaml", ["run_1037_r7a_lsw10_1", "run_1038_r7a_lsw10_2", "run_1039_r7a_lsw10_3"])}
if os.environ.get("R7A_G4"):   # g4 の対 (2026-10-05): 各群 1 本、B は双子の内外を保存して写す (r7a_twin_restart.py)
    GROUPS = {"A": ("problem_3d_prod_m6on_g4.yaml", ["run_1045_r7a_g4_lsw08"]),
              "B": ("problem_3d_prod_m6on_g4_lsw10.yaml", ["run_1046_r7a_g4_lsw10"])}
INPUTS = ("solverConfig.yaml", "bcondConfig.yaml", "sern.h5", "probe.yaml", "prepare_info.json", "cowl_contour.csv", "ramp_contour.csv",
          "MESH_QUALITY.txt", "species_meta.yaml", "IC_FROM.txt", "RESTART_FROM.txt")


def merged_config(p, prob_dst):
    """run_1017 の solverConfig の physProp・turbulence 行だけを現行 YAML の生成値に置き換える"""
    R2.design_snapshot(p); R2.select_operating_point(p, None); st = R2.gas_states(p)
    gen = R3._solver_config(p, 20000, 500, 0.25, float(st["ext"]["P"]))
    s = open(os.path.join(SRC, "solverConfig.yaml")).read()
    for k in ("physProp:", "turbulence:"):
        line = [l for l in gen.splitlines() if l.startswith(k)][0]
        s = re.sub(rf"^{k}.*$", line, s, count=1, flags=re.M)
    open(os.path.join(prob_dst, "solverConfig.yaml"), "w").write(s)
    open(os.path.join(prob_dst, "bcondConfig.yaml"), "w").write(R3._bcond_config(p, st))
    R2.write_species_db(p, prob_dst, R2.frozen_gases(p))


def first_A(dst, prob):
    os.makedirs(dst)
    for f in ("sern.h5", "probe.yaml", "prepare_info.json", "cowl_contour.csv", "ramp_contour.csv", "MESH_QUALITY.txt"):
        shutil.copy(os.path.join(SRC, f), dst)
    merged_config(load_problem(os.path.join(HERE, prob)), dst)
    r = subprocess.run([sys.executable, os.path.join(TOOLS, "restart_field.py"), SRC_RES, os.path.join(dst, "sern.h5"), "--force-species"],
                       capture_output=True, text=True)
    open(os.path.join(dst, "RESTART_FROM.txt"), "w").write(r.stdout + r.stderr)
    if r.returncode:
        raise SystemExit(f"restart_field 失敗: {r.stdout[-500:]}{r.stderr[-500:]}")
    open(os.path.join(dst, "IC_FROM.txt"), "w").write(f"R7a A (L_sw 0.8): restart_field.py {os.path.basename(SRC)} res_20000 -> sern.h5 (同一格子、--force-species)。config = {prob} 生成 + run_1017 の時間・space\n")


def first_B(dst, prob):
    # 格子と入力を runner で作る (段階起動はしない: prepare だけ)
    r = subprocess.run([sys.executable, "-m", "forge_design.evaluate.runner_sern3d", os.path.join(HERE, prob), dst, "--prepare-only"],
                       cwd=os.path.join(ROOT, "design"), capture_output=True, text=True)
    open(os.path.join(dst, "PREPARE.log"), "w").write(r.stdout + r.stderr)
    if r.returncode or not os.path.exists(os.path.join(dst, "sern.h5")):
        raise SystemExit(f"prepare 失敗: {r.stdout[-800:]}{r.stderr[-800:]}")
    merged_config(load_problem(os.path.join(HERE, prob)), dst)
    if os.environ.get("R7A_G4"):   # 双子の内外を保存する写像 (codex 2026-10-05; interp_field は厚さ 0 の板の反対側を拾う)
        r = subprocess.run([sys.executable, os.path.join(HERE, "r7a_twin_restart.py"), SRC_RES, os.path.join(SRC, "sern.h5"), os.path.join(dst, "sern.h5")],
                           capture_output=True, text=True)
    else:
        r = subprocess.run([sys.executable, os.path.join(TOOLS, "interp_field.py"), SRC_RES, os.path.join(dst, "sern.h5"), "--force-species"],
                           capture_output=True, text=True)
    open(os.path.join(dst, "RESTART_FROM.txt"), "w").write(r.stdout + r.stderr)
    if r.returncode:
        raise SystemExit(f"interp_field 失敗: {r.stdout[-800:]}{r.stderr[-800:]}")
    open(os.path.join(dst, "IC_FROM.txt"), "w").write(f"R7a B (L_sw 1.0): {'r7a_twin_restart.py' if os.environ.get('R7A_G4') else 'interp_field.py'} {os.path.basename(SRC)} res_20000 -> sern.h5 (格子を作り直し)。config = {prob} 生成 + run_1017 の時間・space\n")


for g, (prob, runs) in GROUPS.items():
    d0 = os.path.join(HERE, runs[0])
    (first_A if g == "A" else first_B)(d0, prob)
    for r in runs[1:]:
        d = os.path.join(HERE, r); os.makedirs(d)
        for f in INPUTS:
            if os.path.exists(os.path.join(d0, f)):
                shutil.copy(os.path.join(d0, f), d)
        with open(os.path.join(d, "IC_FROM.txt"), "a") as fh:
            fh.write(f"入力は {runs[0]} の複製 (反復)\n")
    # 側壁の実位置 (メッシャは既存 station に丸める): sidewall_in の節点の x 最大
    import h5py, numpy as np
    with h5py.File(os.path.join(d0, "sern.h5")) as f:
        xyz = np.array(f["MESH/COORD"]).reshape(-1, 3)
        bc = open(os.path.join(d0, "bcondConfig.yaml")).read()
        pid = int(re.search(r"^sidewall_in: \{physID: (\d+)", bc, re.M).group(1))
        xs = xyz[f[f"BCONDS/{pid}/iCells"][:], 0]
    print(g, prob, runs, f"側壁 sidewall_in の x 最大 = {xs.max():.6f} m")
