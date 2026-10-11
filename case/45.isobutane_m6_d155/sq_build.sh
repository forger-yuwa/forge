#!/bin/bash
# plan architecture-float-state-double-geometry §6.36: A で #19 の commit を FP64 でビルドし (~/forge-fgeom10-fp64)、量子化の単体試験を回し、
# §6.26 の FP64 (~/forge-fgeom7-fp64) と SASS を照合して記録し (判定には使わない。新しい A* が対照)、**同じバイナリを B の同じパスへ写す**。
# 使い方: bash sq_build.sh <commit>。記録は sq_build.log と _band_ab/cold_pair/sq_pre.json。
set -u
SHA=$1
TOKEN=$(curl -s -X PUT http://169.254.169.254/latest/api/token -H "X-aws-ec2-metadata-token-ttl-seconds: 60")
[ "$(curl -s -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-id)" = "i-0b1a5e0b8dc152f00" ] || { echo "A ではない。中止"; exit 1; }
C45=~/forge-wallfit/case/45.isobutane_m6_d155; LOG=$C45/sq_build.log; PRE=$C45/_band_ab/cold_pair/sq_pre.json
D=~/forge-fgeom10-fp64; S=$D/solver_density_cuda; OLD=~/forge-fgeom7-fp64/solver_density_cuda
W=~/sq_sass; mkdir -p $W
BK="-i $HOME/.ssh/id_a2b -o StrictHostKeyChecking=no"; BR="ubuntu@172.31.14.139"
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
# 単体試験 (量子化の 1 セルの処理、FP64)
cd $S; ut=False
if [ -f tests/unit/test_quantize_state.cu ]; then
  /usr/local/cuda/bin/nvcc --expt-relaxed-constexpr -arch=sm_86 -I. -o $W/test_quantize_state tests/unit/test_quantize_state.cu > $W/test_quantize_state.build.log 2>&1 \
    && $W/test_quantize_state > $W/test_quantize_state.out 2>&1; v=$(grep -h "VERDICT" $W/test_quantize_state.out 2>/dev/null | tail -1); echo "test_quantize_state: $v" >> $LOG
  echo "$v" | grep -q "VERDICT: PASS" && ut=True
else echo "tests/unit/test_quantize_state.cu が無い" >> $LOG; fi
# SASS (記録)
CUOBJ=$(find ~/cuobj_venv -name cuobjdump -type f | head -1); export PATH=$(dirname $CUOBJ):$PATH
$CUOBJ -sass $OLD/build/forge > $W/base.sass; $CUOBJ -sass $S/build/forge > $W/new.sass
python3 -c "import json; json.dump({'commit': '$SHA', 'unit_test_pass': $ut, 'forge_sha256': '$(sha256sum $S/build/forge | cut -c1-64)', 'base_fp64_forge_sha256': '$(sha256sum $OLD/build/forge | cut -c1-64)'}, open('$PRE', 'w'), ensure_ascii=False, indent=1)"
python3 $C45/qa_sass_cmp.py $W/base.sass $W/new.sass $PRE > $W/cmp.json
python3 -c "import json; d=json.load(open('$PRE')); print('単体試験', d['unit_test_pass'], '関数', d['n_base'], d['n_new'], '一致', d['n_identical'], '基準だけ', d['only_base'], '新しいだけ', d['only_new'], '違う', d['differ'])" >> $LOG
# B へ同じバイナリを写す (B で別にビルドしない: nvcc のコード生成の揺れを持ち込まない)
ssh $BK $BR "mkdir -p ~/forge-fgeom10-fp64/solver_density_cuda/build && [ ! -e ~/forge-fgeom10-fp64/solver_density_cuda/build/forge ]" || { echo "B の写し先が既にある/作れない (中止)" >> $LOG; exit 1; }
scp $BK -q $S/build/forge $BR:forge-fgeom10-fp64/solver_density_cuda/build/forge || { echo "B への写しに失敗" >> $LOG; exit 1; }
ssh $BK $BR "echo 'A の ~/forge-fgeom10-fp64 (commit $SHA、FP64) のビルドの forge をそのまま写した (plan §6.36、$(date -Is))' > ~/forge-fgeom10-fp64/COPIED_FROM_A.txt"
echo "B へ写した: $(ssh $BK $BR 'sha256sum ~/forge-fgeom10-fp64/solver_density_cuda/build/forge')" >> $LOG
echo "== $(date -Is) 終わり" >> $LOG; touch $C45/sq_build.done
