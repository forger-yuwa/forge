#!/bin/bash
# condensation-two-phase-default #4g3a: clean build 8a9b673b (+ register usage), regression ON/OFF 200 steps old 6663cbc3 x3 vs new x2,
# FORGE_DIAG_TP_OPERATOR on the G2 finals (ON 0561-0563, OFF 0567-0569), projection mask from run_0592, judge per pair.
set -e; set +o pipefail
R=$HOME/forge-species-tp; C=$R/case/16.nozzle_wys; ST=$HOME/forge-species/g3a_status.txt; N=$R/notes/investigations/2026-10-04-twophase-g3
while pgrep -x make > /dev/null; do sleep 20; done
cd $R && git fetch -q origin feature/species-transport && git checkout -q 8a9b673b
cd solver_density_cuda && rm -rf build && cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_CUDA_ARCHITECTURES=86 -DCMAKE_CXX_FLAGS="-I/usr/local/cuda/include -I/usr/local/cuda/include/cccl" > $HOME/forge-species/build_g3a.log 2>&1 && make -C build -j4 >> $HOME/forge-species/build_g3a.log 2>&1 || { echo BUILDFAIL >> $ST; grep -n " error" $HOME/forge-species/build_g3a.log | head -20 >> $ST; exit 1; }
cp build/forge $HOME/forge-species/forge_8a9b673b; echo "BUILT 8a9b673b" >> $ST
cuobjdump --dump-resource-usage build/forge 2>/dev/null | grep -A1 -E "condensation_source_f_d|twophase_diffusion_d|species_diffusion_d|species_advection_faceY_d" | grep -oE "Function [^ :]+|REG:[0-9]+" | paste - - | head -12 >> $HOME/forge-species/g3a_regs.txt || true
export FORGE_CUDA_BLOCKSIZE=256 FORGE_ALLOW_UNVERIFIED_SPECIES=1
OLD=$HOME/forge-species/forge_6663cbc3; NEW=$HOME/forge-species/forge_8a9b673b; IC=$HOME/forge-species/src_run0482/res_48000.h5
QS=ro,roUx,roUy,roe,P,T,h0,roK,roOmega,roY0,roY1,rog_0,roQ0_0,roQ1_0,roQ2_0,g_0
n=594
for arm in on off; do
  src=$C/run_0561_g2_onprod_r1; [ $arm = off ] && src=$C/run_0567_g2_off_r1
  ro=""; rn=""
  for v in old1 old2 old3 new1 new2; do
    D=$C/run_0${n}_g3areg_${arm}_$v; mkdir $D
    cp $src/solverConfig.yaml $src/bcondConfig.yaml $src/probe.yaml $src/species_db.yaml $src/nozzle_user_2d.h5 $D/
    python3 $R/solver_density_cuda/tools/restart_field.py $IC $D/nozzle_user_2d.h5 --force-species > $D/restart.log 2>&1
    python3 - $D/solverConfig.yaml <<'PY'
import sys, yaml
p=sys.argv[1]; c=yaml.safe_load(open(p)); c["time"]["last"]["nStepOuter"]=200; c["time"]["outStepInterval"]=200
yaml.safe_dump(c, open(p,"w"), default_flow_style=None, sort_keys=False, width=200)
PY
    printf "plan condensation-two-phase-default #4g3a: 診断分岐の回帰 ($arm, $v), IC run_0482 res_48000, 200 step\n" > $D/IC_FROM.txt
    B=$OLD; [[ $v == new* ]] && B=$NEW
    FORGE_BIN=$B $R/solver_density_cuda/tools/run_case.sh $D > $D/run_case.out 2>&1 || true
    [[ $v == old* ]] && ro="$ro $D" || rn="$rn $D"; n=$((n+1))
  done
  python3 $R/solver_density_cuda/tools/check_field_regress.py --repeat $ro --candidate $rn --quantities $QS > $HOME/forge-species/g3areg_${arm}.txt 2>&1 || true
  echo "regress $arm $(grep VERDICT $HOME/forge-species/g3areg_${arm}.txt | tail -1)" >> $ST
done
python3 $N/extract_projection_mask.py $C/run_0592_pj_update_A/tp_update.h5 $HOME/forge-species/proj.npy > $HOME/forge-species/mask.log 2>&1 || echo "mask failed" >> $ST
for src in run_0561_g2_onprod_r1 run_0562_g2_onprod_r2 run_0563_g2_onprod_r3 run_0567_g2_off_r1 run_0568_g2_off_r2 run_0569_g2_off_r3; do
  D=$C/run_0${n}_g3a_op_${src#run_0}; mkdir $D; n=$((n+1))
  cp $C/$src/solverConfig.yaml $C/$src/bcondConfig.yaml $C/$src/probe.yaml $C/$src/species_db.yaml $C/$src/nozzle_user_2d.h5 $D/
  python3 $R/solver_density_cuda/tools/restart_field.py $C/$src/res_48000.h5 $D/nozzle_user_2d.h5 --force-species > $D/restart.log 2>&1
  printf "plan condensation-two-phase-default #4g3a G3-a: 作用素の収支診断 (更新なし)\nIC: $src/res_48000.h5\nbuild 8a9b673b\n診断 run — 収束ゲートの対象外\n" > $D/IC_FROM.txt
  cd $D && { FORGE_DIAG_TP_OPERATOR=$D/tp_operator.h5 $NEW > forge_run.log 2>&1; echo "diag $src rc=$?" >> $ST; } || true
  grep -E "\[tp-operator\].*(coverage|mismatch|refus)" forge_run.log | head -3 >> $ST || true
done
cd $C
ON=( $(ls -d run_06*_g3a_op_*onprod* | sort) ); OFF=( $(ls -d run_06*_g3a_op_*off* | sort) )
for k in 0 1 2; do
  python3 $N/g3a_judge.py ${ON[$k]}/tp_operator.h5 ${OFF[$k]}/tp_operator.h5 --mask $HOME/forge-species/proj.npy --mask-arm on --csv $HOME/forge-species/g3a_pair$k.csv > $HOME/forge-species/G3A_pair$k.txt 2>&1; echo "judge pair $k rc=$?" >> $ST
done
echo DONE >> $ST
