#!/bin/bash
# 最終採用の条件 ⑤ (plan tooling-nozzle-cfd-pinned-initial-line §9 2026-10-05): run_0092 の壁で凝縮 ON 1 本 (run_0093)。IC=run_0092 を convert_species_field (conserve) で。
set -e
: "${FORGE_BIN:=$HOME/forge-integ/solver_density_cuda/build/forge}"; : "${FORGE_CUDA_BLOCKSIZE:=128}"; export FORGE_BIN FORGE_CUDA_BLOCKSIZE
export FORGE_CONVERTER="$(cd "$(dirname "$0")" && pwd)/conv_tolerant.sh"
cd "$(dirname "$0")"
FN=run_0092_ns_c2pin_pass2; CD=run_0093_ns_c2pin_pass2_cond
KF=$(python3 -c "import json;print(json.load(open('c2pin_solve_pass2.json'))['k_f'])")
python3 - "$CD" "$KF" <<'PY'
import sys, json; sys.path.insert(0, "../../design")
from pathlib import Path
from forge_design.evaluate.runner_axismach import prepare_ns
cd, kf = Path(sys.argv[1]), float(sys.argv[2])
info = prepare_ns(Path("problem_d155_ns_c2pin_final_cond.yaml"), cd, nsteps=12000, ic_from=None,
                  initializer={"model": "contur", "a_crocco": 1.0, "cf_scale": kf, "n_scale": 1.0}, cfl_main=1.0, implicit_relax=None)
info["stages"] = {"stages": "none", "ramp": None, "ramp_steps": 1000}; info["restart_from"] = "run_0092_ns_c2pin_pass2"
(cd / "prepare_info.json").write_text(json.dumps(info, indent=1, default=str))
PY
LAST=$(ls $FN/res_[0-9]*.h5 | sort -V | tail -1)
python3 ../../solver_density_cuda/tools/convert_species_field.py "$LAST" "$CD/nozzle.h5" --meta "$CD/species_meta.yaml" --src-run "$FN" --dst-run "$CD" > "$CD/convert_species_field.log" 2>&1
tail -1 "$CD/convert_species_field.log"
python3 -c 'import sys; sys.path.insert(0,"../../design"); from forge_design.evaluate.runner_axismach import run_staged_ns; from pathlib import Path; rc=run_staged_ns(Path(sys.argv[1]), stages="none"); print("forge exit", rc); sys.exit(rc)' $CD
(cd ../../design && python3 -m forge_design.report.nozzle_report ../case/45.isobutane_m6_d155/$CD --euler ../case/45.isobutane_m6_d155/run_0086_euler_wallfit_pincal_r1_ext6k --no-pptx > ../case/45.isobutane_m6_d155/$CD/report_stdout.log 2>&1)
echo "done $CD"; echo ALLDONE
