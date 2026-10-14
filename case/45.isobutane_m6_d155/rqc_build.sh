#!/bin/bash
# plan architecture-float-state-double-geometry §6.54 の準備: 量子化の診断のビルド (~/forge-fgeom11-fp64、commit b0303547) と同じソースの float のビルド
# ~/forge-fgeom11-f32 を B で作る (陽性対照 F 用、codex diagnose 2026-10-12)。flowFormat.hpp は既定 (float) のまま。
# 使い方: bash rqc_build.sh b0303547。記録は rqc_build.log。
set -u
SHA=$1
TOKEN=$(curl -s -X PUT http://169.254.169.254/latest/api/token -H "X-aws-ec2-metadata-token-ttl-seconds: 60")
[ "$(curl -s -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-id)" = "i-0ba2b91ba659254c4" ] || { echo "B ではない。中止"; exit 1; }
C45=~/forge-wallfit/case/45.isobutane_m6_d155; LOG=$C45/rqc_build.log
D=~/forge-fgeom11-f32; S=$D/solver_density_cuda
echo "== $(date -Is) $SHA" > $LOG
git -C ~/forge-species cat-file -e $SHA^{commit} || { echo "$SHA が無い" >> $LOG; exit 1; }
[ "$(git -C ~/forge-fgeom11-fp64 rev-parse HEAD | cut -c1-8)" = "$(git -C ~/forge-species rev-parse $SHA | cut -c1-8)" ] || { echo "~/forge-fgeom11-fp64 のコミットが $SHA と違う" >> $LOG; exit 1; }
[ -e $D ] && { echo "$D が既にある (中止)" >> $LOG; exit 1; }
git -C ~/forge-species worktree add -q --detach $D $SHA >> $LOG 2>&1 || { echo "$D の作成に失敗" >> $LOG; exit 1; }
HF=~/forge-wallfit-bin-fp64/solver_density_cuda/third_party/HighFive
[ -e $S/third_party/HighFive/CMakeLists.txt ] || { rm -rf $S/third_party/HighFive; ln -s $HF $S/third_party/HighFive; }
grep -n "^typedef.*_float" $S/flowFormat.hpp >> $LOG
grep -q "^typedef float flow_float;" $S/flowFormat.hpp || { echo "flowFormat.hpp が既定の float でない" >> $LOG; exit 1; }
# FP64 のビルドとの差は flowFormat.hpp の typedef 4 行だけであることを確かめる
[ -z "$(git -C $D status --porcelain -- . ':!solver_density_cuda/third_party/HighFive')" ] || { echo "作業ツリーに変更がある (HighFive のリンクの差し替えを除く)" >> $LOG; git -C $D status --porcelain >> $LOG; exit 1; }
diff $S/flowFormat.hpp ~/forge-fgeom11-fp64/solver_density_cuda/flowFormat.hpp >> $LOG
mkdir -p $S/build && cd $S/build
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_CUDA_ARCHITECTURES=86 -DCMAKE_CXX_FLAGS="-I/usr/local/cuda/include -I/usr/local/cuda/include/cccl" > cmake.log 2>&1; echo "cmake rc=$?" >> $LOG
make -j4 forge > make.log 2>&1; rc=$?; echo "make rc=$rc $(date -Is) $(sha256sum forge 2>/dev/null)" >> $LOG
[ $rc = 0 ] || exit 1
echo "== $(date -Is) 終わり" >> $LOG; touch $C45/rqc_build.done
