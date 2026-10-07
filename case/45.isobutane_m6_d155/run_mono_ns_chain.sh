#!/bin/bash
# plan tooling-nozzle-throat-monotone-r2 §6 N・K (ユーザ決定 2026-10-06「NS と凝縮計算を回し直して、結果まとめといて」)。
#  N: run_0147_ns_mono_final — 単調壁 (problem_d155_ns_finemesh_recal_final_mono.yaml)、r_t・k_f・Md_moc_offset は run_0117 と同じ。
#     IC = run_0117 res_60000 の保存量を ic_index_map.py (検証付き番号写像) で写す (壁際の第 1 セル 0.34 µm < 壁の移動のため最近傍は使わない)。
#     段階起動なし (IC は同じ格子位相の収束済み NS)、本段 2 次 cfl 1・relax 0.7・60000 step・5000 ごと (run_0117 と同じ本段設定)。
#  K: run_0148_ns_mono_final_cond — run_0118 と同じ手順 (IC = N の res_60000 を convert_species_field、cfl 1・18000 step・1000 ごと)。
#  後処理: nozzle_report (--no-pptx)、exitM_sampling_ab.py (quantities_series.csv)、check_convergence --segment、cond_series.py。
#  Euler 参照は単調壁の Euler run_0143 (同じ設計壁)。報告の pptx はローカルで作る (.venv-pptx)。
# usage (AWS, case dir): bash run_mono_ns_chain.sh   **実行中にこのファイルを編集しないこと**
set -euo pipefail
: "${FORGE_BIN:=$HOME/forge-wallfit-bin/solver_density_cuda/build/forge}"; : "${REAL_CONVERTER:=$HOME/forge-wallfit-bin/solver_density_cuda/build/convertGmshToForge}"
: "${FORGE_CUDA_BLOCKSIZE:=128}"; export FORGE_BIN REAL_CONVERTER FORGE_CUDA_BLOCKSIZE
export FORGE_CONVERTER="$(cd "$(dirname "$0")" && pwd)/conv_tolerant.sh"
cd "$(dirname "$0")"
TOOLS=../../solver_density_cuda/tools
SRC=run_0117_ns_recal_final_ext; EU=run_0143_euler_wallfit_monoG1_r1
N=run_0147_ns_mono_final; K=run_0148_ns_mono_final_cond
PN=problem_d155_ns_finemesh_recal_final_mono.yaml; PK=problem_d155_ns_finemesh_recal_final_mono_cond.yaml
for d in $N $K; do [ -e "$d" ] && { echo "$d が既にある — 止める"; exit 2; }; done
[ -f "$SRC/res_60000.h5" ] || { echo "$SRC/res_60000.h5 が無い"; exit 2; }
KF=$(python3 -c "import json;print(json.load(open('c2pin_solve_recal.json'))['k_f'])"); echo "k_f $KF"
RUNNS='import sys; sys.path.insert(0,"../../design"); from forge_design.evaluate.runner_axismach import run_staged_ns; from pathlib import Path; rc=run_staged_ns(Path(sys.argv[1]), stages="none"); print("forge exit", rc); sys.exit(rc)'
report() { (cd ../../design && python3 -m forge_design.report.nozzle_report ../case/45.isobutane_m6_d155/$1 --euler ../case/45.isobutane_m6_d155/$EU --no-pptx --wall-over-frac 5 > ../case/45.isobutane_m6_d155/$1/report_stdout.log 2>&1) && echo "report $1 rc=0" || echo "report $1 rc=$?"; }
# --- N
python3 prep_c2pin.py $N $KF --problem $PN --stages none
sed -i -E "s/cfl: [0-9.]+, cfl_pseudo: [0-9.]+/cfl: 1.0, cfl_pseudo: 1.0/; s/nStepOuter: [0-9]+/nStepOuter: 60000/; s/outStepInterval: [0-9]+/outStepInterval: 5000/" $N/solverConfig.yaml
grep -E "cfl_pseudo|nStepOuter|outStepInterval|implicitRelax|convMethod" $N/solverConfig.yaml
python3 ic_index_map.py $SRC/res_60000.h5 $N/nozzle.h5 --mode index --forge "$FORGE_BIN" 2>&1 | tail -4
python3 -c "import json,sys; d=json.load(open('$N/IC_MAP.json')); print('IC_MAP', d.get('VERDICT')); sys.exit(0 if d.get('VERDICT')=='OK' else 1)"
python3 -c "$RUNNS" $N
python3 $TOOLS/check_convergence.py $N --segment > $N/CONVERGENCE_VERDICT_segment.txt 2>&1 || true
grep -m1 -E "^=== .*-> " $N/CONVERGENCE_VERDICT_segment.txt || true
report $N
EULER_REF=$EU SOLVE_JSON=c2pin_solve_recal.json FINAL_PROBLEM=$PN python3 exitM_sampling_ab.py $N $SRC > $N/exitM_sampling_stdout.log 2>&1 && echo "series $N ok" || echo "series $N rc=$?"
echo "done $N"
# --- K
python3 - "$K" "$KF" "$N" "$PK" <<'PY'
import sys, json; sys.path.insert(0, "../../design")
from pathlib import Path
from forge_design.evaluate.runner_axismach import prepare_ns
cd, kf = Path(sys.argv[1]), float(sys.argv[2])
info = prepare_ns(Path(sys.argv[4]), cd, nsteps=18000, ic_from=None,
                  initializer={"model": "contur", "a_crocco": 1.0, "cf_scale": kf, "n_scale": 1.0}, cfl_main=1.0, implicit_relax=None)
info["stages"] = {"stages": "none", "ramp": None, "ramp_steps": 1000}; info["restart_from"] = sys.argv[3]
(cd / "prepare_info.json").write_text(json.dumps(info, indent=1, default=str))
PY
grep -E "cfl_pseudo|nStepOuter|outStepInterval|implicitRelax" $K/solverConfig.yaml
python3 $TOOLS/convert_species_field.py "$N/res_60000.h5" "$K/nozzle.h5" --meta "$K/species_meta.yaml" --src-run "$N" --dst-run "$K" > "$K/convert_species_field.log" 2>&1
tail -1 "$K/convert_species_field.log"
python3 -c "$RUNNS" $K
python3 $TOOLS/check_convergence.py $K --segment > $K/CONVERGENCE_VERDICT_segment.txt 2>&1 || true
grep -m1 -E "^=== .*-> " $K/CONVERGENCE_VERDICT_segment.txt || true
report $K
python3 cond_series.py $K > $K/cond_series_stdout.log 2>&1 && echo "cond_series ok" || echo "cond_series rc=$?"
echo "done $K"; echo ALLDONE
