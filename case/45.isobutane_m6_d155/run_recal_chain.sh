#!/bin/bash
# plan tooling-nozzle-cfd-pinned-initial-line §5.1 #11f の NS 連鎖 (事前登録 2026-10-06): 出口較正をやり直した壁 (Md_moc_offset +3.770e-4) で
#  ① 細分 NS pass (k_f 1.055734、r_t 76.6715 mm; IC = run_0109 を interp_field → 段階起動 full → 2 次 cfl 1・60000 step・5000 ごと) → run_0115
#  ② 未緩和 δ_E で k_f・r_t を一発で解く (Euler 参照は新しい固定参照 run_0114)。健全性チェック: r_t 76.67 ± 0.10 mm・k_f 1.056 ± 0.010 を外れたら停止 (諮問)
#  ③ 最終 NS (IC = ① を interp_field → 段階起動 → cfl 1・60000 step) → run_0116。④ (凝縮) は ③ のゲート判定後に別投入 (自動続行しない)
set -e
: "${FORGE_BIN:=$HOME/forge-wallfit-bin/solver_density_cuda/build/forge}"; : "${REAL_CONVERTER:=$HOME/forge-wallfit-bin/solver_density_cuda/build/convertGmshToForge}"
: "${FORGE_CUDA_BLOCKSIZE:=128}"; export FORGE_BIN REAL_CONVERTER FORGE_CUDA_BLOCKSIZE
export FORGE_CONVERTER="$(cd "$(dirname "$0")" && pwd)/conv_tolerant.sh"
cd "$(dirname "$0")"
EU=run_0114_euler_pin_G1_recal_ext6k; P0=run_0109_ns_finemesh_final_ext; F1=run_0115_ns_recal_pass; F2=run_0116_ns_recal_final
RUN='import sys; sys.path.insert(0,"../../design"); from forge_design.evaluate.runner_axismach import run_staged_ns; from pathlib import Path; rc=run_staged_ns(Path(sys.argv[1]), stages="full"); print("forge exit", rc); sys.exit(rc)'
cfl1() { sed -i -E "s/cfl: [0-9.]+, cfl_pseudo: [0-9.]+/cfl: 1.0, cfl_pseudo: 1.0/; s/nStepOuter: [0-9]+/nStepOuter: 60000/; s/outStepInterval: [0-9]+/outStepInterval: 5000/" $1/solverConfig.yaml; }
# ①
if [ ! -f $F1/res_60000.h5 ]; then
  KF0=$(python3 -c "import json;print(json.load(open('c2pin_solve_fine.json'))['k_f'])")
  python3 prep_c2pin.py $F1 $KF0 --problem problem_d155_ns_finemesh_recal.yaml --ic $P0 --stages full
  cfl1 $F1; python3 -c "$RUN" $F1
fi
echo "done $F1"
# ②
python3 c2pin_solve.py $F1 $EU --base problem_d155_ns_finemesh_recal.yaml --out problem_d155_ns_finemesh_recal_final --cond-steps 18000 --cond-out 1000
cp c2pin_solve.json c2pin_solve_recal.json
python3 - <<'PY'
import json, sys
d = json.load(open("c2pin_solve_recal.json")); kf = d["k_f"]; rt = d["solve_rt"]["r_t_m"] * 1e3
ok = abs(kf - 1.056) <= 0.010 and abs(rt - 76.67) <= 0.10
print(f"SANITY k_f {kf:.6f} (1.056±0.010)  r_t {rt:.4f} mm (76.67±0.10)  -> {'IN' if ok else 'OUT'}  hist {d['k_f_hist'][-3:]}")
sys.exit(0 if ok else 3)
PY
echo "done solve"
# ③
KF=$(python3 -c "import json;print(json.load(open('c2pin_solve_recal.json'))['k_f'])")
python3 prep_c2pin.py $F2 $KF --problem problem_d155_ns_finemesh_recal_final.yaml --ic $F1 --stages full
cfl1 $F2; python3 -c "$RUN" $F2
(cd ../../design && python3 -m forge_design.report.nozzle_report ../case/45.isobutane_m6_d155/$F2 --euler ../case/45.isobutane_m6_d155/$EU --no-pptx --wall-over-frac 5 > ../case/45.isobutane_m6_d155/$F2/report_stdout.log 2>&1); echo "report rc=$?"
echo "done $F2"; echo ALLDONE
