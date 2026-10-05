#!/bin/bash
# run_c2pin.sh の ④ (凝縮 ON) だけを回す (2026-10-05: restart_field が液相の種構成差で拒否したため convert_species_field に差し替えて再開)。
set -e
: "${FORGE_BIN:=$HOME/forge-integ/solver_density_cuda/build/forge}"; : "${FORGE_CUDA_BLOCKSIZE:=128}"; export FORGE_BIN FORGE_CUDA_BLOCKSIZE
cd "$(dirname "$0")"
FN=run_0090_ns_c2pin_final; CD=run_0091_ns_c2pin_final_cond
LAST=$(ls $FN/res_[0-9]*.h5 | sort -V | tail -1)
python3 ../../solver_density_cuda/tools/convert_species_field.py "$LAST" "$CD/nozzle.h5" --meta "$CD/species_meta.yaml" --src-run "$FN" --dst-run "$CD" > "$CD/convert_species_field.log" 2>&1
tail -3 "$CD/convert_species_field.log"
python3 -c 'import sys; sys.path.insert(0,"../../design"); from forge_design.evaluate.runner_axismach import run_staged_ns; from pathlib import Path; rc=run_staged_ns(Path(sys.argv[1]), stages="none"); print("forge exit", rc); sys.exit(rc)' $CD
echo "done $CD"; echo ALLDONE
