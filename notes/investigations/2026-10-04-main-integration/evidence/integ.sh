#!/bin/bash
# main 統合の検証 (notes/investigations/2026-10-04-main-integration/README.md の事前登録)。
set +e
H=$HOME; ST=$H/integ/status.txt; mkdir -p $H/integ; : > $ST
CM="-DCMAKE_BUILD_TYPE=Release -DCMAKE_CUDA_ARCHITECTURES=86 -DCMAKE_CXX_FLAGS=-I/usr/local/cuda/include\ -I/usr/local/cuda/include/cccl"
build() { # dir ref tag
  local D=$1 REF=$2 TAG=$3
  if [ ! -d $D ]; then git -C $H/forge-species-tp worktree add -q --detach $D $REF || { echo "WT FAIL $TAG" >> $ST; return 1; }; fi
  git -C $D checkout -q --detach $REF
  HF=$D/solver_density_cuda/third_party/HighFive; if [ -z "$(ls -A $HF 2>/dev/null)" ]; then rmdir $HF 2>/dev/null; ln -s $H/forge-species-tp/solver_density_cuda/third_party/HighFive $HF; fi
  cd $D/solver_density_cuda && rm -rf build
  cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_CUDA_ARCHITECTURES=86 "-DCMAKE_CXX_FLAGS=-I/usr/local/cuda/include -I/usr/local/cuda/include/cccl" > $H/integ/build_$TAG.log 2>&1 && make -C build -j4 >> $H/integ/build_$TAG.log 2>&1 || { echo "BUILDFAIL $TAG" >> $ST; grep -n " error" $H/integ/build_$TAG.log | head -20 >> $ST; return 1; }
  cp build/forge $H/integ/forge_$TAG; echo "BUILT $TAG $(git -C $D rev-parse --short HEAD)" >> $ST
}
cd $H/forge-species-tp && git fetch -q origin +refs/heads/integrate/main-2026-10-04:refs/remotes/origin/integrate/main-2026-10-04 +refs/heads/feature/sern-design:refs/remotes/origin/feature/sern-design
build $H/forge-integ origin/integrate/main-2026-10-04 integ || exit 1
cuobjdump --dump-resource-usage $H/integ/forge_integ 2>/dev/null | grep -E "Function|REG" | grep -B1 -E "REG" | grep -oE "Function [^ :]+|REG:[0-9]+" | paste - - | grep -E "slau|SLAU|species_advection_faceY|twophase|condensation_source" > $H/integ/regs_integ.txt
build $H/forge-integ-sern origin/feature/sern-design sern || exit 1
# CUDA 単体試験 (統合版)
cd $H/forge-integ/solver_density_cuda
for t in test_passive_fct test_passive_scalar test_renorm_gate test_twophase_kernel test_farfield_flux; do
  f=tests/unit/$t.cu; cmd=$(grep -m1 -oE "nvcc .*" $f | sed "s/ -o $t / -o $H\/integ\/$t /")
  [ -z "$cmd" ] && cmd="nvcc -O2 -arch=sm_86 -std=c++17 --expt-relaxed-constexpr -I. -Icuda_forge -o $H/integ/$t $f"
  eval "$cmd" > $H/integ/$t.build.log 2>&1 && { (cd $H/forge-integ/solver_density_cuda && timeout 600 $H/integ/$t > $H/integ/$t.out 2>&1); echo "unit $t rc=$? $(tail -1 $H/integ/$t.out)" >> $ST; } || echo "unit $t BUILDFAIL" >> $ST
done
export FORGE_CUDA_BLOCKSIZE=256 FORGE_ALLOW_UNVERIFIED_SPECIES=1
T=$H/forge-species-tp; R=$H/forge-integ; SP=$H/forge-species/forge_8a9b673b; SE=$H/integ/forge_sern; NEW=$H/integ/forge_integ
QS=ro,roUx,roUy,roe,P,T,h0,roK,roOmega,roY0,roY1,rog_0,roQ0_0,roQ1_0,roQ2_0,g_0
setsteps() { python3 - $1 "$2" <<'PY'
import sys, yaml
p=sys.argv[1]; pin=sys.argv[2]=="pin"; c=yaml.safe_load(open(p)); c["time"]["last"]["nStepOuter"]=200; c["time"]["outStepInterval"]=200
if pin:
    c.setdefault("space",{})["slauWallNormalChi"]=0; c.setdefault("mesh",{})["scalarGradient"]="gg"
yaml.safe_dump(c, open(p,"w"), default_flow_style=None, sort_keys=False, width=200)
PY
}
arm() { # tag case n src ic refbin pin(pin|-) quantities
  local TAG=$1 C=$R/case/$2 n=$3 SRC=$4 IC=$5 REF=$6 PIN=$7 Q=$8; local ro="" rn=""
  for v in ref1 ref2 ref3 new1 new2; do
    D=$C/run_0${n}_integ_${TAG}_$v; mkdir -p $C; mkdir $D || { echo "EXISTS $D" >> $ST; return; }
    cp -r $SRC/. $D/ 2>/dev/null; rm -f $D/res_*.h5 $D/residual_history* $D/forge_run*.log $D/*VERDICT* $D/stage_manifest.json $D/RUN_PROVENANCE.txt
    if [ -n "$IC" ]; then python3 $T/solver_density_cuda/tools/restart_field.py $IC $D/$(basename $(ls $D/*.h5 | grep -v -E "ic.h5|res_" | head -1)) --force-species > $D/restart.log 2>&1; fi
    setsteps $D/solverConfig.yaml $PIN
    printf "main 統合の検証 腕 $TAG ($v): 200 step。notes/investigations/2026-10-04-main-integration/README.md\nsrc: $SRC\n" > $D/IC_FROM.txt
    B=$REF; [[ $v == new* ]] && B=$NEW
    FORGE_BIN=$B $R/solver_density_cuda/tools/run_case.sh $D > $D/run_case.out 2>&1
    echo "$(basename $D) rows=$(grep -c . $D/residual_history.csv 2>/dev/null)" >> $ST
    [[ $v == ref* ]] && ro="$ro $D" || rn="$rn $D"; n=$((n+1))
  done
  python3 $R/solver_density_cuda/tools/check_field_regress.py --repeat $ro --candidate $rn ${Q:+--quantities $Q} > $H/integ/REGRESS_$TAG.txt 2>&1
  echo "regress $TAG $(grep VERDICT $H/integ/REGRESS_$TAG.txt | tail -1)" >> $ST
}
mkin() { local D=$H/integ/in_$1; mkdir -p $D; for f in solverConfig.yaml bcondConfig.yaml probe.yaml species_db.yaml nozzle_user_2d.h5; do cp $2/$f $D/; done; echo $D; }
C16=$T/case/16.nozzle_wys; IC0482=$H/forge-species/src_run0482/res_48000.h5
arm A1 16.nozzle_wys 970 $(mkin A1 $C16/run_0561_g2_onprod_r1) $IC0482 $SP pin $QS
arm A2 16.nozzle_wys 975 $(mkin A2 $C16/run_0567_g2_off_r1) $IC0482 $SP pin $QS
arm A3 16.nozzle_wys 980 $(mkin A3 $C16/run_0524_floor_dry_L1) $C16/run_0524_floor_dry_L1/res_48000.h5 $SE - ""
arm A4 44.vitiated_air_wt 540 $H/forge-species/g0in/c44 "" $SE - ""
echo DONE >> $ST
