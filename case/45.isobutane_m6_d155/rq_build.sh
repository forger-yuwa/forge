#!/bin/bash
# plan architecture-float-state-double-geometry §6.40: B で #20 の commit を FP64 でビルドし (~/forge-fgeom11-fp64)、量子化の単体試験を回す。
# 使い方: bash rq_build.sh <commit>。記録は rq_build.log と _band_ab/cold_pair/rq_pre.json。
set -u
SHA=$1
TOKEN=$(curl -s -X PUT http://169.254.169.254/latest/api/token -H "X-aws-ec2-metadata-token-ttl-seconds: 60")
[ "$(curl -s -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-id)" = "i-0ba2b91ba659254c4" ] || { echo "B ではない。中止"; exit 1; }
C45=~/forge-wallfit/case/45.isobutane_m6_d155; LOG=$C45/rq_build.log; PRE=$C45/_band_ab/cold_pair/rq_pre.json
D=~/forge-fgeom11-fp64; S=$D/solver_density_cuda; W=~/rq_unit; mkdir -p $W
echo "== $(date -Is) 取得 $SHA" > $LOG
git -C ~/forge-species fetch -q https://github.com/forger-yuwa/forge feature/nozzle-wall-fit-and-pipeline >> $LOG 2>&1 || { echo "取得に失敗" >> $LOG; exit 1; }
git -C ~/forge-species cat-file -e $SHA^{commit} || { echo "$SHA が無い" >> $LOG; exit 1; }
[ -e $D ] && { echo "$D が既にある (中止)" >> $LOG; exit 1; }
git -C ~/forge-species worktree add -q --detach $D $SHA >> $LOG 2>&1 || { echo "$D の作成に失敗" >> $LOG; exit 1; }
HF=~/forge-wallfit-bin-fp64/solver_density_cuda/third_party/HighFive
[ -e $S/third_party/HighFive/CMakeLists.txt ] || { rm -rf $S/third_party/HighFive; ln -s $HF $S/third_party/HighFive; }
sed -i 's#^typedef float flow_float;#//typedef float flow_float;#; s#^typedef float geom_float;#//typedef float geom_float;#; s#^typedef float3 flow_float3;#//typedef float3 flow_float3;#; s#^typedef float3 geom_float3;#//typedef float3 geom_float3;#;
        s#^//typedef double flow_float;#typedef double flow_float;#; s#^//typedef double geom_float;#typedef double geom_float;#; s#^//typedef double3 flow_float3;#typedef double3 flow_float3;#; s#^//typedef double3 geom_float3;#typedef double3 geom_float3;#' $S/flowFormat.hpp
grep -n "^typedef.*_float" $S/flowFormat.hpp >> $LOG
grep -q "^typedef double flow_float;" $S/flowFormat.hpp || { echo "FP64 への切り替えに失敗" >> $LOG; exit 1; }
mkdir -p $S/build && cd $S/build
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_CUDA_ARCHITECTURES=86 -DCMAKE_CXX_FLAGS="-I/usr/local/cuda/include -I/usr/local/cuda/include/cccl" > cmake.log 2>&1; echo "cmake rc=$?" >> $LOG
make -j4 forge > make.log 2>&1; rc=$?; echo "make rc=$rc $(date -Is) $(sha256sum forge 2>/dev/null)" >> $LOG
[ $rc = 0 ] || exit 1
cd $S; ut=False
/usr/local/cuda/bin/nvcc --expt-relaxed-constexpr -arch=sm_86 -I. -o $W/test_quantize_state tests/unit/test_quantize_state.cu > $W/build.log 2>&1 \
  && $W/test_quantize_state > $W/test_quantize_state.out 2>&1; v=$(grep -h "VERDICT" $W/test_quantize_state.out 2>/dev/null | tail -1); echo "test_quantize_state: $v" >> $LOG
echo "$v" | grep -q "VERDICT: PASS" && ut=True
python3 -c "import json; json.dump({'commit': '$SHA', 'unit_test_pass': $ut, 'forge_sha256': '$(sha256sum $S/build/forge | cut -c1-64)'}, open('$PRE', 'w'), ensure_ascii=False, indent=1)"
echo "== $(date -Is) 終わり" >> $LOG; touch $C45/rq_build.done
