#!/bin/bash
# P5 (plan tooling-nozzle-cfd-pinned-initial-line §5.1 #7): CFD ピン + joint 壁 + 物理壁解析経路で C2 方式の最終設計を作り直す (AWS)。
# ① pass 1 NS (k_f 1.02573, 段階起動 full) ② E で出口 δ → r_t・k_f を解く ③ 最終 NS (IC=pass 1) ④ 凝縮 ON (restart_field) 。
set -e
: "${FORGE_BIN:=$HOME/forge-integ/solver_density_cuda/build/forge}"; : "${FORGE_CUDA_BLOCKSIZE:=128}"
export FORGE_BIN FORGE_CUDA_BLOCKSIZE
export FORGE_CONVERTER="$(cd "$(dirname "$0")" && pwd)/conv_tolerant.sh"
cd "$(dirname "$0")"
EU=run_0086_euler_wallfit_pincal_r1_ext6k
P1=run_0089_ns_c2pin_pass1; FN=run_0090_ns_c2pin_final; CD=run_0091_ns_c2pin_final_cond
RUNNS='import sys; sys.path.insert(0,"../../design"); from forge_design.evaluate.runner_axismach import run_staged_ns; import json; from pathlib import Path; rd=Path(sys.argv[1]); st=json.loads((rd/"prepare_info.json").read_text())["stages"]["stages"]; rc=run_staged_ns(rd, stages=st); print("forge exit", rc); sys.exit(rc)'
python3 prep_c2pin.py $P1 1.02573
python3 -c "$RUNNS" $P1
echo "done $P1"
python3 c2pin_solve.py $P1 $EU
KF=$(python3 -c "import json;print(json.load(open('c2pin_solve.json'))['k_f'])")
python3 prep_c2pin.py $FN $KF --problem problem_d155_ns_c2pin_final.yaml --ic $P1
python3 -c "$RUNNS" $FN
echo "done $FN"
python3 - "$FN" "$CD" "$KF" <<'PY'
import sys, json; sys.path.insert(0, "../../design")
from pathlib import Path
from forge_design.evaluate.runner_axismach import prepare_ns
fn, cd, kf = Path(sys.argv[1]), Path(sys.argv[2]), float(sys.argv[3])
info = prepare_ns(Path("problem_d155_ns_c2pin_final_cond.yaml"), cd, nsteps=12000, ic_from=None,
                  initializer={"model": "contur", "a_crocco": 1.0, "cf_scale": kf, "n_scale": 1.0}, cfl_main=1.0, implicit_relax=None)
info["stages"] = {"stages": "none", "ramp": None, "ramp_steps": 1000}; info["restart_from"] = str(fn)
(cd / "prepare_info.json").write_text(json.dumps(info, indent=1, default=str))
PY
LAST=$(ls $FN/res_[0-9]*.h5 | sort -V | tail -1)
python3 ../../solver_density_cuda/tools/restart_field.py "$LAST" "$CD/nozzle.h5" --dst-run "$CD" > "$CD/restart_field.log" 2>&1
tail -2 "$CD/restart_field.log"
python3 -c "$RUNNS" $CD
echo "done $CD"
echo ALLDONE
