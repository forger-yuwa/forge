#!/bin/bash
# plan architecture-float-state-double-geometry §6.34 の前提 1: A で #18 の commit を float でビルドし (~/forge-fgeom9-f32)、
# 単体試験 (test_qacc_axis・test_qacc_commit) を回し、§6.26 の float (~/forge-fgeom7-f32) と SASS を照合する。
# 使い方: bash qa_build.sh <commit>。結果は _band_ab/cold_pair/qa_pre.json と qa_build.log。
set -u
SHA=$1
TOKEN=$(curl -s -X PUT http://169.254.169.254/latest/api/token -H "X-aws-ec2-metadata-token-ttl-seconds: 60")
[ "$(curl -s -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-id)" = "i-0b1a5e0b8dc152f00" ] || { echo "A ではない。中止"; exit 1; }
C45=~/forge-wallfit/case/45.isobutane_m6_d155; LOG=$C45/qa_build.log; PRE=$C45/_band_ab/cold_pair/qa_pre.json
D=~/forge-fgeom9-f32; S=$D/solver_density_cuda; OLD=~/forge-fgeom7-f32/solver_density_cuda
W=~/qa_sass; mkdir -p $W/variants
echo "== $(date -Is) 取得 $SHA" > $LOG
git -C ~/forge-species fetch -q https://github.com/forger-yuwa/forge feature/nozzle-wall-fit-and-pipeline >> $LOG 2>&1 || { echo "取得に失敗" >> $LOG; exit 1; }
git -C ~/forge-species cat-file -e $SHA^{commit} || { echo "$SHA が無い" >> $LOG; exit 1; }
[ -e $D ] && { echo "$D が既にある (中止)" >> $LOG; exit 1; }
git -C ~/forge-species worktree add -q --detach $D $SHA >> $LOG 2>&1 || { echo "$D の作成に失敗" >> $LOG; exit 1; }
HF=~/forge-wallfit-bin-fp64/solver_density_cuda/third_party/HighFive
[ -e $S/third_party/HighFive/CMakeLists.txt ] || { rm -rf $S/third_party/HighFive; ln -s $HF $S/third_party/HighFive; }
grep -n "^typedef.*flow_float;" $S/flowFormat.hpp >> $LOG
mkdir -p $S/build && cd $S/build
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_CUDA_ARCHITECTURES=86 -DCMAKE_CXX_FLAGS="-I/usr/local/cuda/include -I/usr/local/cuda/include/cccl" > cmake.log 2>&1; echo "cmake rc=$?" >> $LOG
make -j4 forge > make.log 2>&1; rc=$?; echo "make rc=$rc $(date -Is) $(sha256sum forge 2>/dev/null)" >> $LOG
[ $rc = 0 ] || exit 1
# 単体試験
cd $S; ut=true
for t in test_qacc_axis test_qacc_commit; do
  /usr/local/cuda/bin/nvcc --expt-relaxed-constexpr -arch=sm_86 -I. -o $W/$t tests/unit/$t.cu > $W/$t.build.log 2>&1 || { echo "$t のビルドに失敗" >> $LOG; ut=false; continue; }
  $W/$t > $W/$t.out 2>&1; v=$(grep -h "VERDICT" $W/$t.out | tail -1); echo "$t: $v" >> $LOG
  echo "$v" | grep -q "VERDICT: PASS" || ut=false
done
# SASS の照合は qa_sass.sh (cuobjdump が A に無いので venv の CUDA 13 の cuobjdump を使う)
bash $C45/qa_sass.sh $SHA $( $ut && echo True || echo False ) >> $LOG 2>&1
python3 -c "import json; d=json.load(open('$PRE')); print('単体試験', d['unit_test_pass'], '照合で残った不一致', d['sass_mismatch_except_axis'], '揺れとみなした', list(d.get('differ_resolved_as_nvcc_variant', {})))" >> $LOG
echo "== $(date -Is) 終わり" >> $LOG; touch $C45/qa_build.done
