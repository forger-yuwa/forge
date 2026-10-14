#!/bin/bash
# plan architecture-float-state-double-geometry §5.1 #26 の回帰の準備: 移植の前 (PRE) と後 (POST) のソルバを、同じソースのパス ~/forge-port で
# float と FP64 の両方にクリーンビルドする (パスが違うと命令の並びが変わりビット比較の前提が崩れる、implementer の注意)。インスタンス B。
# 使い方: bash port_build.sh <PRE の commit> <POST の commit>。記録は port_build.log。
set -u
PRE=$1; POST=$2
C45=~/forge-wallfit/case/45.isobutane_m6_d155; LOG=$C45/port_build.log
D=~/forge-port; S=$D/solver_density_cuda
echo "== $(date -Is) PRE $PRE POST $POST" > $LOG
git -C ~/forge-species fetch -q https://github.com/forger-yuwa/forge feature/nozzle-wall-fit-and-pipeline >> $LOG 2>&1 || { echo "取得に失敗" >> $LOG; exit 1; }
for c in $PRE $POST; do git -C ~/forge-species cat-file -e $c^{commit} || { echo "$c が無い" >> $LOG; exit 1; }; done
[ -e $D ] && { echo "$D が既にある (中止)" >> $LOG; exit 1; }
git -C ~/forge-species worktree add -q --detach $D $PRE >> $LOG 2>&1 || { echo "$D の作成に失敗" >> $LOG; exit 1; }
HF=~/forge-wallfit-bin-fp64/solver_density_cuda/third_party/HighFive
build() {   # build <tag>  (今の checkout を float と FP64 で build-<tag>-f32・build-<tag>-fp64 にビルドする)
  local tag=$1
  [ -e $S/third_party/HighFive/CMakeLists.txt ] || { rm -rf $S/third_party/HighFive; ln -s $HF $S/third_party/HighFive; }
  for prec in f32 fp64; do
    if [ $prec = fp64 ]; then
      sed -i 's#^typedef float flow_float;#//typedef float flow_float;#; s#^typedef float geom_float;#//typedef float geom_float;#; s#^typedef float3 flow_float3;#//typedef float3 flow_float3;#; s#^typedef float3 geom_float3;#//typedef float3 geom_float3;#;
              s#^//typedef double flow_float;#typedef double flow_float;#; s#^//typedef double geom_float;#typedef double geom_float;#; s#^//typedef double3 flow_float3;#typedef double3 flow_float3;#; s#^//typedef double3 geom_float3;#typedef double3 geom_float3;#' $S/flowFormat.hpp
      grep -q "^typedef double flow_float;" $S/flowFormat.hpp || { echo "FP64 への切り替えに失敗" >> $LOG; exit 1; }
    fi
    B=$S/build-$tag-$prec; mkdir -p $B
    ( cd $B && cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_CUDA_ARCHITECTURES=86 -DCMAKE_CXX_FLAGS="-I/usr/local/cuda/include -I/usr/local/cuda/include/cccl" > cmake.log 2>&1 && make -j4 forge > make.log 2>&1 ); rc=$?
    echo "$tag $prec rc=$rc $(date -Is) $(sha256sum $B/forge 2>/dev/null | cut -c1-16)" >> $LOG
    [ $prec = fp64 ] && git -C $D checkout -q -- solver_density_cuda/flowFormat.hpp
    [ $rc = 0 ] || exit 1
  done
}
build pre
git -C $D status --porcelain -- . ':!solver_density_cuda/third_party/HighFive' ':!solver_density_cuda/build-*' | head -3 >> $LOG
git -C $D checkout -q --detach $POST || { echo "POST への切り替えに失敗" >> $LOG; exit 1; }
build post
echo "== $(date -Is) 終わり" >> $LOG; touch $C45/port_build.done
