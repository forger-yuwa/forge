#!/bin/bash
# plan thermophysics-solver-owned-species-db #7e: 局所配列を持たない縮約の性能とビット一致。
set +e
H=$HOME; O=$H/integ/v7e; ST=$O/status.txt; mkdir -p $O; : > $ST
R=$H/forge-integ
cd $R && git fetch -q origin +refs/heads/feature/species-transport:refs/remotes/origin/feature/species-transport; git checkout -q --detach 029dc631
[ "$(git -C $R rev-parse --short=8 HEAD)" = "029dc631" ] || { echo "WRONG HEAD" >> $ST; exit 1; }
cd $R/solver_density_cuda && rm -rf build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_CUDA_ARCHITECTURES=86 "-DCMAKE_CXX_FLAGS=-I/usr/local/cuda/include -I/usr/local/cuda/include/cccl" > $O/build.log 2>&1 && make -C build -j4 >> $O/build.log 2>&1 || { echo BUILDFAIL >> $ST; grep -n " error" $O/build.log | head -20 >> $ST; exit 1; }
cp build/forge $O/forge_noarr; echo "BUILT 029dc631" >> $ST; rm -rf build/cuda_forge/CMakeFiles
OLD=$H/integ/forge_integ; ARR=$H/integ/v7c/forge_tab; NOA=$O/forge_noarr
export FORGE_CUDA_BLOCKSIZE=128 FORGE_ALLOW_UNVERIFIED_SPECIES=1 FORGE_LUMPDIFF_TABLE=0
S=$H/forge-r8/case/46.sern_design/run_1033_r11_1; C46=$R/case/46.sern_design
mk46() { local D=$C46/$1; mkdir $D || { echo "EXISTS $D" >> $ST; return 1; }
  for f in solverConfig.yaml bcondConfig.yaml probe.yaml species_meta.yaml sern.h5; do cp $S/$f $D/; done
  cp $S/resolved_species_*.yaml $D/ 2>/dev/null
  python3 $R/solver_density_cuda/tools/restart_field.py $S/res_12000.h5 $D/sern.h5 > $D/restart.log 2>&1 || echo "restart fail $1" >> $ST
  printf "plan thermophysics-solver-owned-species-db #7e: $2\nIC: forge-r8 case/46 run_1033_r11_1 res_12000\n" > $D/IC_FROM.txt; }
mk46 run_1912_lumpdiff_noarr_probe "probe (配列なし)" && (cd $C46/run_1912_lumpdiff_noarr_probe && FORGE_DMIX_PROBE=$PWD/dmix_noarr.bin $NOA > forge_run.log 2>&1 && FORGE_DMIX_PROBE=$PWD/dmix_arr.bin $ARR > forge_run_arr.log 2>&1; cmp dmix_noarr.bin dmix_arr.bin && echo "probe bitwise IDENTICAL (noarr vs arr)" >> $ST || echo "probe DIFFERS" >> $ST)
setsteps() { python3 - $1 <<'PY'
import sys, yaml
p=sys.argv[1]; c=yaml.safe_load(open(p)); c["time"]["last"]["nStepOuter"]=10000; c["time"]["outStepInterval"]=10000
yaml.safe_dump(c, open(p,"w"), default_flow_style=None, sort_keys=False, width=200)
PY
}
n=1913
for v in old1 arr1 noarr1 old2 arr2 noarr2; do
  D=run_${n}_lumpdiff_perf_$v; mk46 $D "SERN 2D 10000 step ($v)" || { n=$((n+1)); continue; }; setsteps $C46/$D/solverConfig.yaml
  B=$NOA; [[ $v == old* ]] && B=$OLD; [[ $v == arr* ]] && B=$ARR
  other=$(nvidia-smi --query-compute-apps=pid --format=csv,noheader | wc -l)
  FORGE_BIN=$B $R/solver_density_cuda/tools/run_case.sh $C46/$D > $C46/$D/run_case.out 2>&1
  other2=$(nvidia-smi --query-compute-apps=pid --format=csv,noheader | wc -l)
  echo "$D $(grep -h "wall, 10000 steps" $C46/$D/forge_run.log | tail -1)  [GPU procs before/after: $other/$other2]" >> $ST
  rm -f $C46/$D/sern.h5 $C46/$D/res_0.h5
  n=$((n+1))
done
echo DONE >> $ST
