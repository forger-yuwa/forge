#!/bin/bash
# §6 (1″): run_0012・run_0013 を wt_nozzle_v2 (区間 1 万・保持なし) で再生する (AWS B、道具は ~/forge-cardmesh-tools、run には書かない)
set -uo pipefail
C=$HOME/forge-cardmesh/case/71.m55_combustion_wt
T=$HOME/forge-cardmesh-tools
PY=$HOME/venv-cardmesh/bin/python
W=$HOME/m55_p5/v2_replay; mkdir -p $W
for spec in "run_0012_ns_free 120000" "run_0013_ns_free260 114000"; do
  set -- $spec; R=$1; TSMAX=$2
  ( cd $T/design && PYTHONPATH=$PWD $PY - "$C/$R" "$W/criteria_v2_$R.yaml" <<'PYC'
import sys, json, yaml, copy
from pathlib import Path
from forge_design.evaluate import settle_presets as SP
run, out = Path(sys.argv[1]), Path(sys.argv[2])
cur = yaml.safe_load((run / "settle_criteria.yaml").read_text())
ev = copy.deepcopy(cur["source"]["evaluate_settle"]); ev["preset"] = "wt_nozzle_v2"
info = json.loads((run / "prepare_info.json").read_text())
cfg = yaml.safe_load((run / "solverConfig.yaml").read_text()); bc = yaml.safe_load((run / "bcondConfig.yaml").read_text())
crit = SP.resolve(ev, info, cfg, bc, segment=SP.main_stage_segment())
out.write_text(yaml.safe_dump(crit, allow_unicode=True, sort_keys=False))
print("wrote", out, crit["source"].get("preset_version"), "W", crit["window"], "hold", crit["hold"])
PYC
  ) > $W/make_$R.log 2>&1 || { echo "make criteria $R failed"; cat $W/make_$R.log; continue; }
  cat $W/make_$R.log
  ( time $PY $T/solver_density_cuda/tools/settle_virtual_stop.py $C/$R --criteria $W/criteria_v2_$R.yaml --out $W/$R --ts-max $TSMAX ) > $W/replay_$R.log 2>&1
  echo "== $R"; cut -c1-320 $W/$R/VERDICT.txt
done
echo REPLAY_DONE
