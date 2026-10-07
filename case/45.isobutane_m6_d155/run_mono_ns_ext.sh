#!/bin/bash
# plan tooling-nozzle-throat-monotone-r2 §6 N (事前登録「未達の量があれば延長 1 回 (20000 step)」): run_0147 のオーバーシュート η0.1 が DRIFTING。
#  run_0149_ns_mono_final_ext: run_0147 の入力を複製し、res_60000 を restart_field (同一メッシュ、ビット一致) で継続。本段 cfl 1・20000 step・5000 ごと。
#  判定: run_0147 の時系列 (5000〜60000) に延長分を step + 60000 で連結し、末尾 5 枚 (60000〜80000) で throat_mono_ns_verdicts.py と同じ判定。
# usage (AWS, case dir): bash run_mono_ns_ext.sh   **実行中にこのファイルを編集しないこと**
set -euo pipefail
: "${FORGE_BIN:=$HOME/forge-wallfit-bin/solver_density_cuda/build/forge}"; : "${REAL_CONVERTER:=$HOME/forge-wallfit-bin/solver_density_cuda/build/convertGmshToForge}"
: "${FORGE_CUDA_BLOCKSIZE:=128}"; export FORGE_BIN REAL_CONVERTER FORGE_CUDA_BLOCKSIZE
export FORGE_CONVERTER="$(cd "$(dirname "$0")" && pwd)/conv_tolerant.sh"
cd "$(dirname "$0")"
TOOLS=../../solver_density_cuda/tools
B=run_0147_ns_mono_final; E=run_0149_ns_mono_final_ext; SRC=run_0117_ns_recal_final_ext; EU=run_0143_euler_wallfit_monoG1_r1
PN=problem_d155_ns_finemesh_recal_final_mono.yaml
[ -e "$E" ] && { echo "$E が既にある — 止める"; exit 2; }
[ -f "$B/res_60000.h5" ] || { echo "$B/res_60000.h5 が無い"; exit 2; }
mkdir "$E"
# 入力だけを複製する (res_*・ログ・報告・判定ファイルは持ち込まない)
for f in solverConfig.yaml bcondConfig.yaml nozzle.h5 prepare_info.json probe.yaml species_meta.yaml MESH_QUALITY.txt wall_design.csv delta_r_initial.csv delta_r_initial.json; do
  [ -e "$B/$f" ] && cp "$B/$f" "$E/"
done
cp $B/resolved_species_*.yaml "$E/" 2>/dev/null || true
ls "$E"
python3 $TOOLS/restart_field.py "$B/res_60000.h5" "$E/nozzle.h5" > "$E/restart_field.log" 2>&1; tail -1 "$E/restart_field.log"
grep -q "ビット一致" "$E/restart_field.log" || { echo "restart_field がビット一致を確認していない — 止める"; exit 1; }
sed -i -E "s/nStepOuter: [0-9]+/nStepOuter: 20000/; s/outStepInterval: [0-9]+/outStepInterval: 5000/" "$E/solverConfig.yaml"
grep -E "cfl_pseudo|nStepOuter|outStepInterval|implicitRelax" "$E/solverConfig.yaml"
python3 - "$E" "$B" <<'PY'
import json, sys
from pathlib import Path
e, b = Path(sys.argv[1]), sys.argv[2]
info = json.loads((e / "prepare_info.json").read_text()); info["stages"] = {"stages": "none", "ramp": None, "ramp_steps": 1000}
info["restart_from"] = f"{b}/res_60000.h5 (restart_field)"; info["extends"] = b
(e / "prepare_info.json").write_text(json.dumps(info, indent=1, default=str))
PY
python3 -c 'import sys; sys.path.insert(0,"../../design"); from forge_design.evaluate.runner_axismach import run_staged_ns; from pathlib import Path; rc=run_staged_ns(Path(sys.argv[1]), stages="none"); print("forge exit", rc); sys.exit(rc)' "$E"
python3 $TOOLS/check_convergence.py "$E" --segment > "$E/CONVERGENCE_VERDICT_segment.txt" 2>&1 || true
grep -m1 -E "^=== .*-> " "$E/CONVERGENCE_VERDICT_segment.txt" || true
(cd ../../design && python3 -m forge_design.report.nozzle_report ../case/45.isobutane_m6_d155/$E --euler ../case/45.isobutane_m6_d155/$EU --no-pptx --wall-over-frac 5 > ../case/45.isobutane_m6_d155/$E/report_stdout.log 2>&1) && echo "report $E rc=0" || echo "report $E rc=$?"
EULER_REF=$EU SOLVE_JSON=c2pin_solve_recal.json FINAL_PROBLEM=$PN python3 exitM_sampling_ab.py $E $SRC > $E/exitM_sampling_stdout.log 2>&1 && echo "series $E ok" || echo "series $E rc=$?"
# 連結 (延長分は step + 60000)
python3 - "$B" "$E" <<'PY'
import csv, sys
from pathlib import Path
b, e = Path(sys.argv[1]), Path(sys.argv[2])
rb = list(csv.DictReader(open(b / "quantities_series.csv"))); re_ = list(csv.DictReader(open(e / "quantities_series.csv")))
for r in re_:
    r["step"] = str(int(float(r["step"])) + 60000)
rows = rb + [r for r in re_ if int(r["step"]) > 60000]
with open(e / "quantities_series_joined.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rb[0].keys())); w.writeheader(); w.writerows(rows)
print("joined rows", len(rows), "last step", rows[-1]["step"])
PY
cp "$E/quantities_series.csv" "$E/quantities_series_extonly.csv"; cp "$E/quantities_series_joined.csv" "$E/quantities_series.csv"
python3 throat_mono_ns_verdicts.py "$E" dry
echo ALLDONE
