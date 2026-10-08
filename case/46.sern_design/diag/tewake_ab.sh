#!/bin/bash
# plan convection-zero-thickness-edge-reconstruction §6.0: 格子の判別 A/B (A′ = te_wake_blend_H 0、B = 1.0)。
# 各 20000 step を同時に回す (共通バイナリ・共通設定、w と 1 面介入は無効)。500 step ごとの全場は te_monitor.py で
# 後縁まわりの温度を記録してから消す (最後の 1 枚は残す)。終了後に collect (力・ゲート)・check_quasisteady (力 4 量と最低温度)。
cd /home/ubuntu/forge-r8
C=case/46.sern_design
mon() { r=$1; lb=$2
  while true; do
    done_flag=0; [ -f $C/$r/.done ] && done_flag=1
    files=$(ls $C/$r | grep -E '^res_[0-9]+\.h5$' | sed 's/res_\([0-9]*\)\.h5/\1/' | sort -n)
    last=$(echo "$files" | tail -1)
    for n in $files; do
      [ "$n" = "$last" ] && [ $done_flag = 0 ] && continue
      if [ "$n" = 0 ]; then
        python3 $C/diag/te_monitor.py $C/$r $C/$r/res_0.h5 --csv $C/$r/TE_MONITOR_step0.csv --lb-h $lb >> $C/$r/TE_MONITOR.log 2>&1
        rm -f $C/$r/res_0.h5 $C/$r/res_0.xmf; continue
      fi
      grep -q "^$n," $C/$r/TE_MONITOR.csv 2>/dev/null || python3 $C/diag/te_monitor.py $C/$r $C/$r/res_$n.h5 --lb-h $lb >> $C/$r/TE_MONITOR.log 2>&1
      [ "$n" != "$last" ] && rm -f $C/$r/res_$n.h5 $C/$r/res_$n.xmf
    done
    [ $done_flag = 1 ] && break
    sleep 30
  done; }
run() { r=$1; lb=$2; prob=$3
  ( mon $r $lb ) &
  FORGE_CUDA_BLOCKSIZE=128 FORGE_ALLOW_UNVERIFIED_SPECIES=1 solver_density_cuda/tools/run_case.sh $C/$r > $C/$r/run_case_stdout.log 2>&1 < /dev/null
  rc=$?; echo "$r rc=$rc $(date)"; touch $C/$r/.done; wait
  (cd design && FORGE_ALLOW_UNVERIFIED_SPECIES=1 python3 -c "
from forge_design.evaluate import runner_sern3d as R
out=R.collect('../$C/$prob','../$C/$r',rc=$rc)
g=out['gates']; print('GATES',g['verdict'],g['fail_class'],g['reasons'])
" > ../$C/$r/COLLECT.txt 2>&1)
  python3 solver_density_cuda/tools/check_quasisteady.py --series-csv $C/$r/force_history.csv --series-cols C_T,C_T_with_shear,C_L,C_M --tail 0.5 > $C/$r/QUASISTEADY_forces.txt 2>&1
  python3 solver_density_cuda/tools/check_quasisteady.py --series-csv $C/$r/TE_MONITOR.csv --series-cols R_TE_Tmin,R_SE_Tmin,RET_Tmin,ALL_Tmin --tail 0.5 > $C/$r/QUASISTEADY_Tmin.txt 2>&1
  (cd $C && python3 v3_farfield_eval.py $r > $r/V3_EVAL.txt 2>&1)
  echo "$r post done $(date)"; }
run run_1078_tewake_A0_m10 1.0 problem_3d_prod_3op_wallres_lswx08_tewake_A0.yaml &
sleep 120
run run_1079_tewake_B10_m10 1.0 problem_3d_prod_3op_wallres_lswx08_tewake_B10.yaml &
wait
echo TEWAKE_AB_DONE $(date)
