#!/bin/bash
# #7e 性能: GPU が自分だけのときに測る。測定中に他プロセスが出た回は捨てて測り直す (最大 3 回)。
set +e
H=$HOME; O=$H/integ/v7e; ST=$O/perf.txt; : > $ST
R=$H/forge-integ; C46=$R/case/46.sern_design; S=$H/forge-r8/case/46.sern_design/run_1033_r11_1
OLD=$H/integ/forge_integ; ARR=$H/integ/v7c/forge_tab; NOA=$O/forge_noarr
export FORGE_CUDA_BLOCKSIZE=128 FORGE_ALLOW_UNVERIFIED_SPECIES=1 FORGE_LUMPDIFF_TABLE=0
D=$C46/run_1919_lumpdiff_perf_clean; mkdir -p $D
for f in solverConfig.yaml bcondConfig.yaml probe.yaml species_meta.yaml sern.h5; do cp $S/$f $D/; done; cp $S/resolved_species_*.yaml $D/
python3 $R/solver_density_cuda/tools/restart_field.py $S/res_12000.h5 $D/sern.h5 > $D/restart.log 2>&1
python3 - $D/solverConfig.yaml <<'PY'
import sys, yaml
p=sys.argv[1]; c=yaml.safe_load(open(p)); c["time"]["last"]["nStepOuter"]=5000; c["time"]["outStepInterval"]=100000
yaml.safe_dump(c, open(p,"w"), default_flow_style=None, sort_keys=False, width=200)
PY
printf "plan thermophysics-solver-owned-species-db #7e 性能 (GPU 占有時のみ、5000 step、同じディレクトリで binary を替えて計測)\n" > $D/IC_FROM.txt
procs() { nvidia-smi --query-compute-apps=pid --format=csv,noheader | grep -c . ; }
cd $D
for round in 1 2; do
 for v in old arr noarr; do
  B=$NOA; [ $v = old ] && B=$OLD; [ $v = arr ] && B=$ARR
  for att in 1 2 3; do
    w=0; while [ "$(procs)" != "0" ] && [ $w -lt 720 ]; do sleep 30; w=$((w+1)); done   # 最大 6 時間待つ
    ( $B > forge_run_$v.log 2>&1 ) & fp=$!
    clean=1; while kill -0 $fp 2>/dev/null; do [ "$(procs)" -gt 1 ] && clean=0; sleep 2; done
    t=$(grep -h "wall, 5000 steps" forge_run_$v.log | tail -1)
    rm -f res_*.h5 res_*.xmf
    if [ $clean = 1 ]; then echo "round $round $v: $t" >> $ST; break; else echo "round $round $v attempt $att: contaminated ($t)" >> $ST; fi
  done
 done
done
rm -f sern.h5
echo DONE >> $ST
