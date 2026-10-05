#!/bin/bash
# plan tooling-nozzle-cfd-pinned-initial-line §5.1 #11 ④ (ユーザ決定 2026-10-06「A」: ③ の波 η0.1 DRIFTING は未達のまま記録して ④ へ)。
# 最終 NS run_0117 の壁 (k_f 1.054129、r_t 76.6539 mm) で凝縮 ON: IC = run_0117 res_60000 を convert_species_field (conserve)、cfl 1・18000 step・1000 ごと → run_0118
set -e
: "${FORGE_BIN:=$HOME/forge-wallfit-bin/solver_density_cuda/build/forge}"; : "${REAL_CONVERTER:=$HOME/forge-wallfit-bin/solver_density_cuda/build/convertGmshToForge}"
: "${FORGE_CUDA_BLOCKSIZE:=128}"; export FORGE_BIN REAL_CONVERTER FORGE_CUDA_BLOCKSIZE
export FORGE_CONVERTER="$(cd "$(dirname "$0")" && pwd)/conv_tolerant.sh"
cd "$(dirname "$0")"
SRC=run_0117_ns_recal_final_ext; CD=run_0118_ns_recal_final_cond; EU=run_0114_euler_pin_G1_recal_ext6k
KF=$(python3 -c "import json;print(json.load(open('c2pin_solve_recal.json'))['k_f'])")
python3 - "$CD" "$KF" "$SRC" <<'PY'
import sys, json; sys.path.insert(0, "../../design")
from pathlib import Path
from forge_design.evaluate.runner_axismach import prepare_ns
cd, kf = Path(sys.argv[1]), float(sys.argv[2])
info = prepare_ns(Path("problem_d155_ns_finemesh_recal_final_cond.yaml"), cd, nsteps=18000, ic_from=None,
                  initializer={"model": "contur", "a_crocco": 1.0, "cf_scale": kf, "n_scale": 1.0}, cfl_main=1.0, implicit_relax=None)
info["stages"] = {"stages": "none", "ramp": None, "ramp_steps": 1000}; info["restart_from"] = sys.argv[3]
(cd / "prepare_info.json").write_text(json.dumps(info, indent=1, default=str))
PY
grep -E "cfl_pseudo|nStepOuter|outStepInterval|implicitRelax" $CD/solverConfig.yaml
python3 ../../solver_density_cuda/tools/convert_species_field.py "$SRC/res_60000.h5" "$CD/nozzle.h5" --meta "$CD/species_meta.yaml" --src-run "$SRC" --dst-run "$CD" > "$CD/convert_species_field.log" 2>&1
tail -1 "$CD/convert_species_field.log"
python3 -c 'import sys; sys.path.insert(0,"../../design"); from forge_design.evaluate.runner_axismach import run_staged_ns; from pathlib import Path; rc=run_staged_ns(Path(sys.argv[1]), stages="none"); print("forge exit", rc); sys.exit(rc)' $CD
(cd ../../design && python3 -m forge_design.report.nozzle_report ../case/45.isobutane_m6_d155/$CD --euler ../case/45.isobutane_m6_d155/$EU --no-pptx --wall-over-frac 5 > ../case/45.isobutane_m6_d155/$CD/report_stdout.log 2>&1); echo "report rc=$?"
python3 cond_series.py $CD
echo "done $CD"; echo ALLDONE
