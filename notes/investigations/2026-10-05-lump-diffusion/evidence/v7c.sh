#!/bin/bash
# plan thermophysics-solver-owned-species-db #7c の AWS 検証: 表引き版の V4f と性能、表 vs 式の場の比較。
set +e
H=$HOME; O=$H/integ/v7c; ST=$O/status.txt; mkdir -p $O; : > $ST
R=$H/forge-integ
cd $R && git fetch -q origin +refs/heads/feature/species-transport:refs/remotes/origin/feature/species-transport; git checkout -q --detach 36cadb3b
[ "$(git -C $R rev-parse --short=8 HEAD)" = "36cadb3b" ] || { echo "WRONG HEAD" >> $ST; exit 1; }
cd $R/solver_density_cuda && rm -rf build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_CUDA_ARCHITECTURES=86 "-DCMAKE_CXX_FLAGS=-I/usr/local/cuda/include -I/usr/local/cuda/include/cccl" > $O/build.log 2>&1 && make -C build -j4 >> $O/build.log 2>&1 || { echo BUILDFAIL >> $ST; grep -n " error" $O/build.log | head -20 >> $ST; exit 1; }
cp build/forge $O/forge_tab; echo "BUILT 36cadb3b" >> $ST; rm -rf build/cuda_forge/CMakeFiles
OLD=$H/integ/forge_integ; NEW=$O/forge_tab
export FORGE_CUDA_BLOCKSIZE=128 FORGE_ALLOW_UNVERIFIED_SPECIES=1
S=$H/forge-r8/case/46.sern_design/run_1033_r11_1; C46=$R/case/46.sern_design
mk46() { local D=$C46/$1; mkdir $D || { echo "EXISTS $D" >> $ST; return 1; }
  for f in solverConfig.yaml bcondConfig.yaml probe.yaml species_meta.yaml sern.h5; do cp $S/$f $D/; done
  cp $S/resolved_species_*.yaml $D/ 2>/dev/null
  python3 $R/solver_density_cuda/tools/restart_field.py $S/res_12000.h5 $D/sern.h5 > $D/restart.log 2>&1 || echo "restart fail $1" >> $ST
  printf "plan thermophysics-solver-owned-species-db #7c 検証: $2\nIC: forge-r8 case/46 run_1033_r11_1 res_12000\nbuild 36cadb3b\n" > $D/IC_FROM.txt; }
mk46 run_1905_lumpdiff_tab_probe "V4f probe (表引き)" && (cd $C46/run_1905_lumpdiff_tab_probe && FORGE_DMIX_PROBE=$PWD/dmix.bin $NEW > forge_run.log 2>&1; echo "probe tab rc=$?" >> $ST)
python3 $R/solver_density_cuda/tests/unit/check_dmix_probe.py $C46/run_1905_lumpdiff_tab_probe/dmix.bin $(ls -t $C46/run_1905_lumpdiff_tab_probe/resolved_species_*.yaml | head -1) > $O/V4f_tab.txt 2>&1; echo "V4f tab $(tail -1 $O/V4f_tab.txt)" >> $ST
grep -h "lump diffusion" $C46/run_1905_lumpdiff_tab_probe/forge_run.log >> $ST
setsteps() { python3 - $1 <<'PY'
import sys, yaml
p=sys.argv[1]; c=yaml.safe_load(open(p)); c["time"]["last"]["nStepOuter"]=10000; c["time"]["outStepInterval"]=10000
yaml.safe_dump(c, open(p,"w"), default_flow_style=None, sort_keys=False, width=200)
PY
}
n=1906
for v in old1 form1 tab1 old2 form2 tab2; do
  D=run_${n}_lumpdiff_tabperf_$v; mk46 $D "SERN 2D 10000 step 継続 ($v)" || { n=$((n+1)); continue; }; setsteps $C46/$D/solverConfig.yaml
  B=$NEW; E=""; [[ $v == old* ]] && B=$OLD; [[ $v == form* ]] && E="FORGE_LUMPDIFF_TABLE=0"
  env $E FORGE_BIN=$B $R/solver_density_cuda/tools/run_case.sh $C46/$D > $C46/$D/run_case.out 2>&1
  echo "$D $(grep -h "wall, 10000 steps" $C46/$D/forge_run.log | tail -1)" >> $ST
  rm -f $C46/$D/sern.h5 $C46/$D/res_0.h5
  n=$((n+1))
done
python3 $R/solver_density_cuda/tools/check_field_regress.py --repeat $C46/run_1907_lumpdiff_tabperf_form1 $C46/run_1910_lumpdiff_tabperf_form2 --candidate $C46/run_1908_lumpdiff_tabperf_tab1 $C46/run_1911_lumpdiff_tabperf_tab2 --quantities ro,roUx,roUy,roe,P,T,roK,roOmega,roY0,roY1,Y0,Y1 > $O/tab_vs_formula.txt 2>&1
echo "tab vs formula $(grep VERDICT $O/tab_vs_formula.txt | tail -1)" >> $ST
echo DONE >> $ST
