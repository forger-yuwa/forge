#!/usr/bin/env python3
"""生産 3D config (現行の problem YAML から runner_sern3d が生成するもの) で、既存 run の最終場から同一格子で継続する run を作る
(AWS で実行)。2026-10-01: farfield (side_far) + 種ごとの輸送物性 (R8) の組み合わせの初回確認用。
  - solverConfig.yaml / bcondConfig.yaml / species_meta.yaml / species_db_external.yaml は problem YAML から runner で生成
    (時間・出力・CFL・space の行は SRC run の solverConfig を引き継ぎ、physProp と turbulence の行だけ生成値に置き換える)
  - 格子と初期場は SRC run の sern.h5 と最終 res を restart_field.py で (種の記録が無い/名前が変わった場は --force-species で 1 回だけ移行)
  python3 prod_restart_setup.py DST SRC_RUN SRC_RES PROBLEM_YAML [--force-species]
"""
import os, re, shutil, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "design"))
from forge_design.evaluate import runner_sern as R2, runner_sern3d as R3  # noqa: E402
from forge_design.probdef import load_problem  # noqa: E402


def main():
    dst, src, src_res, prob = sys.argv[1:5]
    force = "--force-species" in sys.argv
    os.makedirs(dst)
    for f in ("sern.h5", "probe.yaml", "prepare_info.json", "cowl_contour.csv", "ramp_contour.csv", "MESH_QUALITY.txt"):
        if os.path.exists(os.path.join(src, f)):
            shutil.copy(os.path.join(src, f), dst)
    p = load_problem(prob)
    R2.design_snapshot(p); R2.select_operating_point(p, None); st = R2.gas_states(p)
    gen = R3._solver_config(p, 20000, 500, 0.25, float(st["ext"]["P"]))
    gl = {k: [l for l in gen.splitlines() if l.startswith(k)][0] for k in ("physProp:", "turbulence:")}
    s = open(os.path.join(src, "solverConfig.yaml")).read()
    for k, v in gl.items():
        s = re.sub(rf"^{k}.*$", v, s, count=1, flags=re.M)
    open(os.path.join(dst, "solverConfig.yaml"), "w").write(s)
    open(os.path.join(dst, "bcondConfig.yaml"), "w").write(R3._bcond_config(p, st))
    R2.write_species_db(p, dst, R2.frozen_gases(p))
    cmd = [sys.executable, os.path.join(ROOT, "solver_density_cuda", "tools", "restart_field.py"), src_res, os.path.join(dst, "sern.h5")]
    if force:
        cmd.append("--force-species")
    r = subprocess.run(cmd, capture_output=True, text=True)
    open(os.path.join(dst, "RESTART_FROM.txt"), "w").write(r.stdout + r.stderr)
    if r.returncode != 0:
        raise SystemExit(f"restart_field に失敗:\n{r.stdout[-800:]}{r.stderr[-800:]}")
    open(os.path.join(dst, "IC_FROM.txt"), "w").write(
        f"IC: restart_field.py {src_res} -> sern.h5{' (--force-species: 種の名前変更 AIR→AMB を 1 回だけ移行)' if force else ''}。"
        f"config は {os.path.basename(prob)} (現行の生産: side_far farfield + 種ごとの輸送物性) から生成、時間・space は {src} を引き継ぐ\n")
    print("prepared", dst)
    print(" ", gl["physProp:"][:160], "...")
    print("  side_far:", [l for l in open(os.path.join(dst, "bcondConfig.yaml")).read().splitlines() if l.startswith("side_far")][0][:90])


if __name__ == "__main__":
    main()
