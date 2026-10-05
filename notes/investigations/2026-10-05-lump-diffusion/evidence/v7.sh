#!/bin/bash
# plan thermophysics-solver-owned-species-db #7b の AWS 検証 (V4b CFD 回帰・V4f probe 照合・性能・V4g 記録)。
set +e
H=$HOME; O=$H/integ/v7; ST=$O/status.txt; mkdir -p $O; : > $ST
R=$H/forge-integ; T=$H/forge-species-tp
cd $R && git fetch -q origin +refs/heads/feature/species-transport:refs/remotes/origin/feature/species-transport; git checkout -q --detach 91a79925
[ "$(git -C $R rev-parse --short=8 HEAD)" = "91a79925" ] || { echo "WRONG HEAD $(git -C $R rev-parse --short HEAD)" >> $ST; exit 1; }
cd $R/solver_density_cuda && rm -rf build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_CUDA_ARCHITECTURES=86 "-DCMAKE_CXX_FLAGS=-I/usr/local/cuda/include -I/usr/local/cuda/include/cccl" > $O/build.log 2>&1 && make -C build -j4 >> $O/build.log 2>&1 || { echo BUILDFAIL >> $ST; grep -n " error" $O/build.log | head -20 >> $ST; exit 1; }
cp build/forge $O/forge_new; echo "BUILT $(git -C $R rev-parse --short HEAD)" >> $ST
rm -rf build/cuda_forge/CMakeFiles   # 空き容量 (オブジェクトは不要)
cuobjdump --dump-resource-usage $O/forge_new 2>/dev/null | grep "REG:" | grep -iE "species_diffusion|Dmix|audit|tpfd|slau" > $O/regs_new.txt
cuobjdump --dump-resource-usage $H/integ/forge_integ 2>/dev/null | grep "REG:" | grep -iE "species_diffusion|Dmix|audit|tpfd|slau" > $O/regs_old.txt
g++ -O2 -std=c++17 -I. -o $O/tld tests/unit/test_lump_diffusion_host.cpp > $O/tld.build 2>&1 && $O/tld > $O/tld.out 2>&1; echo "host unit rc=$? $(tail -1 $O/tld.out)" >> $ST
OLD=$H/integ/forge_integ; NEW=$O/forge_new
export FORGE_CUDA_BLOCKSIZE=128 FORGE_ALLOW_UNVERIFIED_SPECIES=1
setsteps() { python3 - $1 $2 $3 <<'PY'
import sys, yaml
p=sys.argv[1]; n=int(sys.argv[2]); oi=int(sys.argv[3]); c=yaml.safe_load(open(p)); c["time"]["last"]["nStepOuter"]=n; c["time"]["outStepInterval"]=oi
yaml.safe_dump(c, open(p,"w"), default_flow_style=None, sort_keys=False, width=200)
PY
}
# ---- V4f / 性能 / V4g: SERN 2D (run_1033_r11_1 の入力を読み取りで複製; res_12000 から継続)
S=$H/forge-r8/case/46.sern_design/run_1033_r11_1; C46=$R/case/46.sern_design; mkdir -p $C46
mk46() { local D=$C46/$1; mkdir $D || { echo "EXISTS $D" >> $ST; return 1; }
  for f in solverConfig.yaml bcondConfig.yaml probe.yaml species_meta.yaml sern.h5; do cp $S/$f $D/; done
  cp $S/resolved_species_*.yaml $D/ 2>/dev/null
  python3 $R/solver_density_cuda/tools/restart_field.py $S/res_12000.h5 $D/sern.h5 > $D/restart.log 2>&1 || echo "restart fail $1" >> $ST
  printf "plan thermophysics-solver-owned-species-db #7b 検証: $2\nIC: forge-r8 case/46 run_1033_r11_1 res_12000 (restart_field)\n" > $D/IC_FROM.txt; }
mk46 run_1900_lumpdiff_probe_new "V4f probe (新)" && (cd $C46/run_1900_lumpdiff_probe_new && FORGE_DMIX_PROBE=$PWD/dmix.bin $NEW > forge_run.log 2>&1; echo "probe new rc=$?" >> $ST)
python3 $R/solver_density_cuda/tests/unit/check_dmix_probe.py $C46/run_1900_lumpdiff_probe_new/dmix.bin $(ls -t $C46/run_1900_lumpdiff_probe_new/resolved_species_*.yaml | head -1) > $O/V4f_sern.txt 2>&1; echo "V4f sern $(tail -1 $O/V4f_sern.txt)" >> $ST
n=1901
for v in old1 new1 old2 new2; do
  D=run_${n}_lumpdiff_cont_$v; mk46 $D "SERN 2D 継続 10000 step ($v; 性能と変化量の記録)" || continue; setsteps $C46/$D/solverConfig.yaml 10000 10000
  B=$OLD; [[ $v == new* ]] && B=$NEW
  t0=$(date +%s); FORGE_BIN=$B $R/solver_density_cuda/tools/run_case.sh $C46/$D > $C46/$D/run_case.out 2>&1; t1=$(date +%s)
  echo "$D wall $((t1-t0)) s rows=$(grep -c . $C46/$D/residual_history.csv 2>/dev/null) loop_time=$(grep -iE "elapsed|wall|total time" $C46/$D/forge_run.log | tail -1)" >> $ST
  n=$((n+1))
done
python3 $R/solver_density_cuda/tools/check_field_regress.py --repeat $C46/run_1901_lumpdiff_cont_old1 $C46/run_1903_lumpdiff_cont_old2 --candidate $C46/run_1902_lumpdiff_cont_new1 $C46/run_1904_lumpdiff_cont_new2 > $O/SERN_change.txt 2>&1
echo "sern change $(grep VERDICT $O/SERN_change.txt | tail -1) (記録; lump 拡散の変更で差が出るのが期待)" >> $ST
# ---- V4b: case/16 (lump 無し) 旧 ×3 vs 新 ×2、200 step
C16=$T/case/16.nozzle_wys; IC0482=$H/forge-species/src_run0482/res_48000.h5; RC=$R/case/16.nozzle_wys; mkdir -p $RC
QS=ro,roUx,roUy,roe,P,T,h0,roK,roOmega,roY0,roY1,rog_0,roQ0_0,roQ1_0,roQ2_0,g_0
arm() { local TAG=$1 n=$2 SRC=$3 IC=$4 Q=$5; local ro="" rn=""
  for v in old1 old2 old3 new1 new2; do
    D=$RC/run_0${n}_lumpdiff_${TAG}_$v; mkdir $D || { echo "EXISTS $D" >> $ST; return; }
    for f in solverConfig.yaml bcondConfig.yaml probe.yaml species_db.yaml nozzle_user_2d.h5; do cp $SRC/$f $D/; done
    python3 $T/solver_density_cuda/tools/restart_field.py $IC $D/nozzle_user_2d.h5 --force-species > $D/restart.log 2>&1
    setsteps $D/solverConfig.yaml 200 200
    printf "plan thermophysics-solver-owned-species-db #7b V4b ($TAG, $v): lump 無しの不変、200 step\n" > $D/IC_FROM.txt
    B=$OLD; [[ $v == new* ]] && B=$NEW
    FORGE_BIN=$B $R/solver_density_cuda/tools/run_case.sh $D > $D/run_case.out 2>&1
    echo "$(basename $D) rows=$(grep -c . $D/residual_history.csv 2>/dev/null)" >> $ST
    [[ $v == old* ]] && ro="$ro $D" || rn="$rn $D"; n=$((n+1))
  done
  python3 $R/solver_density_cuda/tools/check_field_regress.py --repeat $ro --candidate $rn ${Q:+--quantities $Q} > $O/REGRESS_$TAG.txt 2>&1
  echo "regress $TAG $(grep VERDICT $O/REGRESS_$TAG.txt | tail -1)" >> $ST
  rm -f $RC/run_0*_lumpdiff_${TAG}_*/nozzle_user_2d.h5 $RC/run_0*_lumpdiff_${TAG}_*/res_0.h5
}
arm dry 985 $C16/run_0524_floor_dry_L1 $C16/run_0524_floor_dry_L1/res_48000.h5 ""
arm wetoff 990 $C16/run_0567_g2_off_r1 $IC0482 $QS
arm weton 995 $C16/run_0561_g2_onprod_r1 $IC0482 $QS
echo DONE >> $ST
