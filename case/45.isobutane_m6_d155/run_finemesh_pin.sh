#!/bin/bash
# plan tooling-nozzle-cfd-pinned-initial-line §5.1 #11 (事前登録 2026-10-05): 細分メッシュの生産 pass。
#  ① run_0092 の壁 (r_t 76.7531 mm, k_f = c2pin_solve_pass2.json) を細分メッシュで NS、IC = run_0092 (interp_field)、12000 step・1000 step ごと → run_0098 (段階起動)
#  ② 未緩和 δ_E で k_f・r_t を一発で解く (予測 k_f 1.037 ± 0.005、r_t 76.75 ± 0.03 mm; 外れたら止めて諮問)
#  ③ 最終 NS (IC = ① の場、interp_field)、12000 step・1000 step ごと → run_0099
#  ④ 凝縮 ON (IC = run_0099 を convert_species_field --mode conserve)、18000 step・1000 step ごと → run_0100
#  報告は --no-pptx (pptx はローカルで作る)。判定 (check_*) は終了後に別途。
set -e
: "${FORGE_BIN:=$HOME/forge-integ/solver_density_cuda/build/forge}"; : "${FORGE_CUDA_BLOCKSIZE:=128}"; export FORGE_BIN FORGE_CUDA_BLOCKSIZE
export FORGE_CONVERTER="$(cd "$(dirname "$0")" && pwd)/conv_tolerant.sh"
cd "$(dirname "$0")"
EU=run_0086_euler_wallfit_pincal_r1_ext6k; P0=run_0092_ns_c2pin_pass2
F1=run_0098_ns_finemesh_pass_staged; F2=run_0099_ns_finemesh_final; F3=run_0100_ns_finemesh_final_cond
# 2026-10-05: 1 回目 (run_0095、IC 補間 → 本段 cfl 5 直行) は step 20 で発散 (x/r_t 55〜70 の壁際の薄セルで T → 6000 K 上限)。
# 手順書 divergence-and-startup の段階起動 (soft 1 次 cfl 0.5 → mid 1 次 cfl 1 → 本段) を同じ IC に掛ける。番号 0096/0097 は欠番。
RUN='import sys; sys.path.insert(0,"../../design"); from forge_design.evaluate.runner_axismach import run_staged_ns; from pathlib import Path; rc=run_staged_ns(Path(sys.argv[1]), stages=(sys.argv[2] if len(sys.argv) > 2 else "none")); print("forge exit", rc); sys.exit(rc)'
report() { (cd ../../design && python3 -m forge_design.report.nozzle_report ../case/45.isobutane_m6_d155/$1 --euler ../case/45.isobutane_m6_d155/$EU --no-pptx > ../case/45.isobutane_m6_d155/$1/report_stdout.log 2>&1) || echo "report $1 failed (rc=$?)"; }

# ① (段階起動: 本段の残差履歴の最終 step で完了を判定)
if [ ! -f $F1/residual_history.csv ] || [ "$(tail -1 $F1/residual_history.csv | cut -d, -f1)" != "11999" ]; then
  KF0=$(python3 -c "import json;print(json.load(open('c2pin_solve_pass2.json'))['k_f'])")
  python3 prep_c2pin.py $F1 $KF0 --problem problem_d155_ns_finemesh_pin.yaml --ic $P0 --stages full
  python3 -c "$RUN" $F1 full
fi
echo "done $F1"

# ②
python3 c2pin_solve.py $F1 $EU --base problem_d155_ns_finemesh_pin.yaml --out problem_d155_ns_finemesh_pin_final --cond-steps 18000 --cond-out 1000
cp c2pin_solve.json c2pin_solve_fine.json
python3 - <<'PY'
import json, sys
d = json.load(open("c2pin_solve_fine.json")); kf = d["k_f"]; rt = d["solve_rt"]["r_t_m"] * 1e3
ok = abs(kf - 1.037) <= 0.005 and abs(rt - 76.75) <= 0.03
print(f"PRED k_f {kf:.5f} (1.037±0.005)  r_t {rt:.4f} mm (76.75±0.03)  -> {'IN' if ok else 'OUT'}")
sys.exit(0 if ok else 3)
PY
echo "done solve"

# ③
KF=$(python3 -c "import json;print(json.load(open('c2pin_solve_fine.json'))['k_f'])")
python3 prep_c2pin.py $F2 $KF --problem problem_d155_ns_finemesh_pin_final.yaml --ic $F1
python3 -c "$RUN" $F2
report $F2
echo "done $F2"

# ④
python3 - "$F3" "$KF" "$F2" <<'PY'
import sys, json; sys.path.insert(0, "../../design")
from pathlib import Path
from forge_design.evaluate.runner_axismach import prepare_ns
cd, kf = Path(sys.argv[1]), float(sys.argv[2])
info = prepare_ns(Path("problem_d155_ns_finemesh_pin_final_cond.yaml"), cd, nsteps=18000, ic_from=None,
                  initializer={"model": "contur", "a_crocco": 1.0, "cf_scale": kf, "n_scale": 1.0}, cfl_main=1.0, implicit_relax=None)
info["stages"] = {"stages": "none", "ramp": None, "ramp_steps": 1000}; info["restart_from"] = sys.argv[3]
(cd / "prepare_info.json").write_text(json.dumps(info, indent=1, default=str))
PY
LAST=$(ls $F2/res_[0-9]*.h5 | sort -V | tail -1)
python3 ../../solver_density_cuda/tools/convert_species_field.py "$LAST" "$F3/nozzle.h5" --meta "$F3/species_meta.yaml" --src-run "$F2" --dst-run "$F3" > "$F3/convert_species_field.log" 2>&1
tail -1 "$F3/convert_species_field.log"
python3 -c "$RUN" $F3
report $F3
echo "done $F3"; echo ALLDONE
