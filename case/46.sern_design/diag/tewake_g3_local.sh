#!/bin/bash
# plan tooling-sern-te-wake-grid §5.1 #4c/#4d: g3 の局所化格子 (run_1084) で冷点の消失を再確認する。
# 初期場 = 冷点を含む run_1055 の最終場 (restart_field_deformed、src-mesh = run_1083 の 0 格子 [座標は run_1055 と同じ署名])、20000 step。
set -u
cd /home/ubuntu/forge-r8
C=case/46.sern_design; R=run_1084_tewake_g3_L10_m10; SRC=$C/run_1055_r7b_m10_A_c; REF=$C/run_1083_tewake_g3_A0_ref
export FORGE_ALLOW_UNVERIFIED_SPECIES=1
cp $C/$R/solverConfig.yaml $C/$R/solverConfig.yaml.prepare
sed -i "s/nStepOuter: [0-9]*/nStepOuter: 20000/" $C/$R/solverConfig.yaml
diff $C/$R/solverConfig.yaml.prepare $C/$R/solverConfig.yaml > $C/$R/solverConfig.yaml.diff_from_prepare
python3 solver_density_cuda/tools/restart_field_deformed.py $SRC/res_20000.h5 $C/$R/sern.h5 --src-mesh $REF/sern.h5 --list-moved $C/$R/MOVED_NODES.csv > $C/$R/IC_TRANSFER.txt 2>&1; echo "ic rc=$?"; tail -3 $C/$R/IC_TRANSFER.txt
echo "IC: restart_field_deformed run_1055_r7b_m10_A_c res_20000 (冷点を含む床到達場、診断用の継続元) → run_1084 (局所化格子)、--src-mesh run_1083 (0 格子、座標署名 = run_1055)。plan tooling-sern-te-wake-grid §5.1 #4d" > $C/$R/IC_FROM.txt
mon() { r=$1
  while true; do
    done_flag=0; [ -f $C/$r/.done ] && done_flag=1
    files=$(ls $C/$r | grep -E '^res_[0-9]+\.h5$' | sed 's/res_\([0-9]*\)\.h5/\1/' | sort -n); last=$(echo "$files" | tail -1)
    for n in $files; do
      [ "$n" = "$last" ] && [ $done_flag = 0 ] && continue
      if [ "$n" = 0 ]; then python3 $C/diag/te_monitor.py $C/$r $C/$r/res_0.h5 --csv $C/$r/TE_MONITOR_step0.csv --lb-h 1.0 >> $C/$r/TE_MONITOR.log 2>&1; rm -f $C/$r/res_0.h5 $C/$r/res_0.xmf; continue; fi
      grep -q "^$n," $C/$r/TE_MONITOR.csv 2>/dev/null || python3 $C/diag/te_monitor.py $C/$r $C/$r/res_$n.h5 --lb-h 1.0 >> $C/$r/TE_MONITOR.log 2>&1
      [ "$n" != "$last" ] && rm -f $C/$r/res_$n.h5 $C/$r/res_$n.xmf
    done
    [ $done_flag = 1 ] && break; sleep 30
  done; }
( mon $R ) &
FORGE_CUDA_BLOCKSIZE=128 solver_density_cuda/tools/run_case.sh $C/$R > $C/$R/run_case_stdout.log 2>&1 < /dev/null
rc=$?; echo "$R rc=$rc $(date)"; touch $C/$R/.done; wait
(cd design && python3 -c "
from forge_design.evaluate import runner_sern3d as R
out=R.collect('../$C/problem_3d_prod_3op_wallres_lswx08_tewake_g3_L10.yaml','../$C/$R',rc=$rc)
g=out['gates']; print('GATES',g['verdict'],g['fail_class'],g['reasons'])
" > ../$C/$R/COLLECT.txt 2>&1)
python3 solver_density_cuda/tools/check_floor_events.py $C/$R --window-steps 10000 > $C/$R/FLOOR_EVENTS_VERDICT.txt 2>&1; echo "floor rc=$?"
python3 solver_density_cuda/tools/check_quasisteady.py --series-csv $C/$R/force_history.csv --series-cols C_T,C_T_with_shear,C_L,C_M --tail 0.5 > $C/$R/QUASISTEADY_forces.txt 2>&1
python3 solver_density_cuda/tools/check_quasisteady.py --series-csv $C/$R/TE_MONITOR.csv --series-cols R_TE_Tmin,R_SE_Tmin,RET_Tmin,ALL_Tmin --tail 0.5 > $C/$R/QUASISTEADY_Tmin.txt 2>&1
(cd $C && python3 v3_farfield_eval.py $R > $R/V3_EVAL.txt 2>&1)
echo G3_LOCAL_DONE $(date)
