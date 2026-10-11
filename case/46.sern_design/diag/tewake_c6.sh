#!/bin/bash
# plan tooling-sern-te-wake-grid §5.1 #6a: 処置の前後 (A = te_wake 0 / B = 局所化) を m6_on・m4_off で各 20000 step (AWS-B)。
set -u
cd /home/ubuntu/forge-r8
export FORGE_ALLOW_UNVERIFIED_SPECIES=1 FORGE_CUDA_BLOCKSIZE=128
C=case/46.sern_design; PA=problem_3d_prod_3op_wallres_lswx08_tewake_g3_A0.yaml; PB=problem_3d_prod_3op_wallres_lswx08_tewake_g3_L10.yaml
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
prep() { prob=$1; r=$2; op=$3
  [ -e $C/$r ] && { echo "EXISTS $r"; exit 1; }
  (cd design && python3 -m forge_design.evaluate.runner_sern3d ../$C/$prob ../$C/$r --op $op --prepare-only > ../$C/$r.prepare.log 2>&1 < /dev/null); echo "$r prepare rc=$?"
  rm -f $C/$r/sern.msh
  cp $C/$r/solverConfig.yaml $C/$r/solverConfig.yaml.prepare; sed -i "s/nStepOuter: [0-9]*/nStepOuter: 20000/" $C/$r/solverConfig.yaml
  diff $C/$r/solverConfig.yaml.prepare $C/$r/solverConfig.yaml > $C/$r/solverConfig.yaml.diff_from_prepare; }
run1() { prob=$1; r=$2
  ( mon $r ) &
  solver_density_cuda/tools/run_case.sh $C/$r > $C/$r/run_case_stdout.log 2>&1 < /dev/null
  rc=$?; echo "$r rc=$rc $(date)"; touch $C/$r/.done; wait
  (cd design && python3 -c "
from forge_design.evaluate import runner_sern3d as R
out=R.collect('../$C/$prob','../$C/$r',rc=$rc)
g=out['gates']; print('GATES',g['verdict'],g['fail_class'],g['reasons'])
" > ../$C/$r/COLLECT.txt 2>&1)
  python3 solver_density_cuda/tools/check_floor_events.py $C/$r --window-steps 10000 > $C/$r/FLOOR_EVENTS_VERDICT.txt 2>&1
  python3 solver_density_cuda/tools/check_quasisteady.py --series-csv $C/$r/force_history.csv --series-cols C_T,C_T_with_shear,C_L,C_M --tail 0.5 > $C/$r/QUASISTEADY_forces.txt 2>&1
  python3 solver_density_cuda/tools/check_quasisteady.py --series-csv $C/$r/TE_MONITOR.csv --series-cols R_TE_Tmin,R_SE_Tmin,RET_Tmin,ALL_Tmin --tail 0.5 > $C/$r/QUASISTEADY_Tmin.txt 2>&1
  (cd $C && python3 v3_farfield_eval.py $r > $r/V3_EVAL.txt 2>&1)
  echo "$r post done $(date)"; }
do_op() { op=$1; src=$2; A=$3; B=$4
  prep $PA $A $op; prep $PB $B $op
  python3 solver_density_cuda/tools/restart_field.py $C/$src/res_20000.h5 $C/$A/sern.h5 > $C/$A/RESTART_FROM.txt 2>&1; echo "$A ic rc=$?"; tail -1 $C/$A/RESTART_FROM.txt
  echo "IC: restart_field $src res_20000 (旧評価の最終場)。plan tooling-sern-te-wake-grid §5.1 #6a の A (te_wake 0)" > $C/$A/IC_FROM.txt
  python3 solver_density_cuda/tools/restart_field_deformed.py $C/$src/res_20000.h5 $C/$B/sern.h5 --src-mesh $C/$A/sern.h5 --list-moved $C/$B/MOVED_NODES.csv > $C/$B/IC_TRANSFER.txt 2>&1; echo "$B ic rc=$?"; tail -1 $C/$B/IC_TRANSFER.txt
  echo "IC: restart_field_deformed $src res_20000 → 局所化格子、src-mesh = $A の格子。plan tooling-sern-te-wake-grid §5.1 #6a の B" > $C/$B/IC_FROM.txt
  run1 $PA $A & sleep 120; run1 $PB $B & wait; }
do_op m6_on run_1050_r7b_x_lsw08 run_1087_tewake_m6_A0 run_1088_tewake_m6_L10
do_op m4_off run_1058_r7b_m4_A_c run_1089_tewake_m4_A0 run_1090_tewake_m4_L10
echo C6_DONE $(date)
