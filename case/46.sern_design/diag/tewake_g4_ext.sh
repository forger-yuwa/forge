#!/bin/bash
# plan tooling-sern-te-wake-grid §6 (事前登録: 未定常なら +20000 を 1 回): run_1085 の最終場から同一設定で +20000 (run_1086)。
set -u
cd /home/ubuntu/forge-r8
export FORGE_ALLOW_UNVERIFIED_SPECIES=1 FORGE_CUDA_BLOCKSIZE=128
C=case/46.sern_design; P=problem_3d_prod_3op_wallres_lswx08_tewake_g4_L10.yaml
S=run_1085_tewake_g4_L10_cont20k; N=run_1086_tewake_g4_L10_cont40k
[ -e $C/$N ] && { echo "EXISTS $N"; exit 1; }
mkdir $C/$N; for f in solverConfig.yaml bcondConfig.yaml sern.h5 probe.yaml prepare_info.json cowl_contour.csv ramp_contour.csv MESH_QUALITY.txt species_meta.yaml; do cp $C/$S/$f $C/$N/; done
python3 solver_density_cuda/tools/restart_field.py $C/$S/res_20000.h5 $C/$N/sern.h5 > $C/$N/RESTART_FROM.txt 2>&1 || python3 solver_density_cuda/tools/restart_field.py $C/$S/res_20000.h5 $C/$N/sern.h5 --force-species >> $C/$N/RESTART_FROM.txt 2>&1
tail -1 $C/$N/RESTART_FROM.txt
echo "IC: restart_field $S res_20000 (同一設定の +20000 延長、事前登録の延長 1 回)。plan tooling-sern-te-wake-grid §6" > $C/$N/IC_FROM.txt
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
echo G4_EXT_DONE $(date)
