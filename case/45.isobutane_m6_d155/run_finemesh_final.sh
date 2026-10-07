#!/bin/bash
# plan tooling-nozzle-cfd-pinned-initial-line §5.1 #11 ③ (codex diagnose 2026-10-05 pred-out: ③ 単独、④ は ③ の全ゲート判定後に別投入)。
# ② の結果 (c2pin_solve_fine.json: k_f 1.055734、r_t 76.6715 mm、problem_d155_ns_finemesh_pin_final.yaml) で最終 NS:
# IC = run_0103 res_60000 (interp_field) → 段階起動 full → 本段 2 次 cfl 1・60000 step・5000 ごと。報告は失敗を隠さない。
set -e
: "${FORGE_BIN:=$HOME/forge-wallfit-bin/solver_density_cuda/build/forge}"; : "${REAL_CONVERTER:=$HOME/forge-wallfit-bin/solver_density_cuda/build/convertGmshToForge}"
: "${FORGE_CUDA_BLOCKSIZE:=128}"; export FORGE_BIN REAL_CONVERTER FORGE_CUDA_BLOCKSIZE
export FORGE_CONVERTER="$(cd "$(dirname "$0")" && pwd)/conv_tolerant.sh"
cd "$(dirname "$0")"
F1=run_0103_ns_finemesh_pass_cfl1; F2=run_0107_ns_finemesh_final; EU=run_0086_euler_wallfit_pincal_r1_ext6k
KF=$(python3 -c "import json;print(json.load(open('c2pin_solve_fine.json'))['k_f'])")
python3 prep_c2pin.py $F2 $KF --problem problem_d155_ns_finemesh_pin_final.yaml --ic $F1 --stages full
sed -i -E "s/cfl: [0-9.]+, cfl_pseudo: [0-9.]+/cfl: 1.0, cfl_pseudo: 1.0/; s/nStepOuter: [0-9]+/nStepOuter: 60000/; s/outStepInterval: [0-9]+/outStepInterval: 5000/" $F2/solverConfig.yaml
python3 -c 'import sys; sys.path.insert(0,"../../design"); from forge_design.evaluate.runner_axismach import run_staged_ns; from pathlib import Path; rc=run_staged_ns(Path(sys.argv[1]), stages="full"); print("forge exit", rc); sys.exit(rc)' $F2
(cd ../../design && python3 -m forge_design.report.nozzle_report ../case/45.isobutane_m6_d155/$F2 --euler ../case/45.isobutane_m6_d155/$EU --no-pptx --wall-over-frac 5 > ../case/45.isobutane_m6_d155/$F2/report_stdout.log 2>&1); echo "report rc=$?"
echo "done $F2"; echo ALLDONE
