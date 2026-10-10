#!/bin/bash
# plan tooling-sern-te-wake-grid §5.1 #4 / #4f: g4 の局所化格子 (run_1082、prepare 済み) を領域別の一様な初期場から既存 3D の段階起動
# → 本段の最終場から同一設定の 20000 step の判定区間 (run_1085)。共通バイナリ (カウンタ入り)、evaluate.floor_events 1。
set -u
cd /home/ubuntu/forge-r8
export FORGE_ALLOW_UNVERIFIED_SPECIES=1 FORGE_CUDA_BLOCKSIZE=128
C=case/46.sern_design; P=problem_3d_prod_3op_wallres_lswx08_tewake_g4_L10.yaml
R=run_1082_tewake_g4_L10_m10; N=run_1085_tewake_g4_L10_cont20k
# 0) 双子の割り当て (prepare が書いた初期場 = 領域別の一様な初期場)
python3 $C/diag/twin_ic_check.py --g3 $C/run_1084_tewake_g3_L10_m10 --src-res $C/run_1084_tewake_g3_L10_m10/res_20000.h5 --g4 $C/$R --json $C/$R/TWIN_IC.json > $C/$R/TWIN_IC_VERDICT.txt 2>&1; trc=$?
tail -1 $C/$R/TWIN_IC_VERDICT.txt
[ $trc -ne 0 ] && { echo "TWIN NOT OK rc=$trc -> stop"; echo G4_LOCAL_DONE_STOP $(date); exit 0; }
# 1) 段階起動 (runner_sern3d.main の prepare 以降と同じ呼び出し)
(cd design && python3 - <<PY > ../$C/$R.staged.log 2>&1
from forge_design.evaluate.runner_sern3d import run_staged, collect
from forge_design.evaluate.runner_sern import load_problem
o = (load_problem("../$C/$P").raw.get("opt") or {})
rc = run_staged("../$C/$R", "full", int(o.get("soft_steps", 2000)), soft_cfl=float(o.get("soft_cfl", 0.5)), soft_conv=int(o.get("soft_conv", 0)),
                warm_lam_steps=int(o.get("warm_lam_steps", 0)), warm_lam_cfl=float(o.get("warm_lam_cfl", 0.2)), mid_steps=int(o.get("mid_steps", 0)))
out = collect("../$C/$P", "../$C/$R", rc=rc)
g = out["gates"]; print("STAGED rc", rc, "GATES", g["verdict"], g["fail_class"], g["reasons"])
PY
)
tail -2 $C/$R.staged.log
last=$(ls $C/$R | grep -E '^res_[0-9]+\.h5$' | sed 's/res_\([0-9]*\)\.h5/\1/' | sort -n | tail -1)
[ -z "$last" ] && { echo "no res in $R -> stop"; echo G4_LOCAL_DONE_STOP $(date); exit 0; }
# 2) 判定区間 20000 step (同一設定・restart_field)
[ -e $C/$N ] && { echo "EXISTS $N"; exit 1; }
mkdir $C/$N; for f in solverConfig.yaml bcondConfig.yaml sern.h5 probe.yaml prepare_info.json cowl_contour.csv ramp_contour.csv MESH_QUALITY.txt species_meta.yaml; do cp $C/$R/$f $C/$N/; done
sed -i "s/nStepOuter: [0-9]*/nStepOuter: 20000/" $C/$N/solverConfig.yaml
python3 solver_density_cuda/tools/restart_field.py $C/$R/res_$last.h5 $C/$N/sern.h5 > $C/$N/RESTART_FROM.txt 2>&1 || python3 solver_density_cuda/tools/restart_field.py $C/$R/res_$last.h5 $C/$N/sern.h5 --force-species >> $C/$N/RESTART_FROM.txt 2>&1
tail -1 $C/$N/RESTART_FROM.txt
echo "IC: restart_field $R res_$last (段階起動の本段の最終場、同一設定で判定区間 20000)。plan tooling-sern-te-wake-grid §5.1 #4" > $C/$N/IC_FROM.txt
mon() { r=$1
  while true; do
    done_flag=0; [ -f $C/$r/.done ] && done_flag=1
    files=$(ls $C/$r | grep -E '^res_[0-9]+\.h5$' | sed 's/res_\([0-9]*\)\.h5/\1/' | sort -n); lst=$(echo "$files" | tail -1)
    for n in $files; do
      [ "$n" = "$lst" ] && [ $done_flag = 0 ] && continue
      if [ "$n" = 0 ]; then python3 $C/diag/te_monitor.py $C/$r $C/$r/res_0.h5 --csv $C/$r/TE_MONITOR_step0.csv --lb-h 1.0 >> $C/$r/TE_MONITOR.log 2>&1; rm -f $C/$r/res_0.h5 $C/$r/res_0.xmf; continue; fi
      grep -q "^$n," $C/$r/TE_MONITOR.csv 2>/dev/null || python3 $C/diag/te_monitor.py $C/$r $C/$r/res_$n.h5 --lb-h 1.0 >> $C/$r/TE_MONITOR.log 2>&1
      [ "$n" != "$lst" ] && rm -f $C/$r/res_$n.h5 $C/$r/res_$n.xmf
    done
    [ $done_flag = 1 ] && break; sleep 30
  done; }
( mon $N ) &
solver_density_cuda/tools/run_case.sh $C/$N > $C/$N/run_case_stdout.log 2>&1 < /dev/null
rc=$?; echo "$N rc=$rc $(date)"; touch $C/$N/.done; wait
(cd design && python3 -c "
from forge_design.evaluate import runner_sern3d as R
out=R.collect('../$C/$P','../$C/$N',rc=$rc)
g=out['gates']; print('GATES',g['verdict'],g['fail_class'],g['reasons'])
" > ../$C/$N/COLLECT.txt 2>&1)
python3 solver_density_cuda/tools/check_floor_events.py $C/$N --window-steps 10000 > $C/$N/FLOOR_EVENTS_VERDICT.txt 2>&1; echo "floor rc=$?"
python3 solver_density_cuda/tools/check_quasisteady.py --series-csv $C/$N/force_history.csv --series-cols C_T,C_T_with_shear,C_L,C_M --tail 0.5 > $C/$N/QUASISTEADY_forces.txt 2>&1
python3 solver_density_cuda/tools/check_quasisteady.py --series-csv $C/$N/TE_MONITOR.csv --series-cols R_TE_Tmin,R_SE_Tmin,RET_Tmin,ALL_Tmin --tail 0.5 > $C/$N/QUASISTEADY_Tmin.txt 2>&1
(cd $C && python3 v3_farfield_eval.py $N > $N/V3_EVAL.txt 2>&1)
echo G4_LOCAL_DONE $(date)
