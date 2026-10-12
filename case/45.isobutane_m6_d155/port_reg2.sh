#!/bin/bash
# plan architecture-float-state-double-geometry §6.60 (事前登録、codex diagnose 2026-10-13): 移植の回帰 R1・R3 の回し直し。インスタンス A で GPU を占有して逐次。
# バイナリは B の ~/forge-port で同じパスでビルドした移植の前・後 (port_build.log の sha) を A へ写して使う。
# R1: case/45 の G0・_tr/q32_init.h5 から point (run_0713_pt_d1 の設定) と B0 (run_0717_c50_d1 の設定) の 10 step、float と FP64 → run_0780〜0803。
# R3: case/48 の _smoke の設定・入力の 10 step、float と FP64 → case/48 の run_0071〜0082。
# 各組で前・後を交互に 3 本ずつ (pre1 post1 pre2 post2 pre3 post3)。判定は port_reg2_an.py。
set -uo pipefail
TOKEN=$(curl -s -X PUT http://169.254.169.254/latest/api/token -H "X-aws-ec2-metadata-token-ttl-seconds: 60")
[ "$(curl -s -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-id)" = "i-0b1a5e0b8dc152f00" ] || { echo "A ではない。中止"; exit 1; }
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
export FORGE_OMEGA_BUDGET=1
C45=~/forge-wallfit/case/45.isobutane_m6_d155; C48=~/forge-wallfit/case/48.flat_plate_cooled_m4; LOG=$C45/port_reg2.log
RC=~/forge-wallfit/solver_density_cuda/tools/run_case.sh; BIN=~/forge-port-bins
die() { echo "中止: $* ($(date -Is))" >> $LOG; exit 1; }
B_IP=172.31.14.139; K="-i $HOME/.ssh/id_a2b -o StrictHostKeyChecking=accept-new"
echo "== 準備 $(date -Is)" >> $LOG
mkdir -p $BIN $C45/_reg_src/pt $C45/_reg_src/b0
for t in pre post; do for p in f32 fp64; do
  mkdir -p $BIN/$t-$p; [ -e $BIN/$t-$p/forge ] || scp -q $K ubuntu@$B_IP:forge-port/solver_density_cuda/build-$t-$p/forge $BIN/$t-$p/forge || die "バイナリ $t-$p の取得"
done; done
declare -A WANT=([pre-f32]=9e203e1f03c87f60 [pre-fp64]=00f95b7b88c682c9 [post-f32]=a20f395f2a7f110f [post-fp64]=75ee4675b0c90f41)
for k in "${!WANT[@]}"; do s=$(sha256sum $BIN/$k/forge | cut -c1-16); echo "bin $k $s" >> $LOG; [ "$s" = "${WANT[$k]}" ] || die "バイナリ $k の sha が port_build.log と違う"; done
scp -q $K "ubuntu@$B_IP:forge-wallfit/case/45.isobutane_m6_d155/run_0713_pt_d1/*.yaml" $C45/_reg_src/pt/ || die "point の設定の取得"
scp -q $K "ubuntu@$B_IP:forge-wallfit/case/45.isobutane_m6_d155/run_0717_c50_d1/*.yaml" $C45/_reg_src/b0/ || die "B0 の設定の取得"
mkdir -p $C48/_smoke
for f in solverConfig.yaml bcondConfig.yaml probe.yaml mesh_init.h5; do   # _smoke にはサブディレクトリもあるので必要なファイルだけ (2026-10-13: まとめて写して失敗した)
  [ -e $C48/_smoke/$f ] || scp -q $K "ubuntu@$B_IP:forge-wallfit/case/48.flat_plate_cooled_m4/_smoke/$f" $C48/_smoke/$f || die "case/48 の _smoke/$f の取得"
done
[ "$(sha256sum $C45/_tr/q32_init.h5 | cut -c1-8)" = "ba859e02" ] || die "case/45 の入力が §6.26 と違う"
echo "入力: case/45 $(sha256sum $C45/_tr/q32_init.h5 | cut -c1-16)、case/48 $(sha256sum $C48/_smoke/mesh_init.h5 | cut -c1-16)" >> $LOG
waitfree() { for i in $(seq 1 120); do [ -z "$(pgrep -x forge)" ] && return 0; sleep 30; done; die "ほかの forge が 1 時間空かない (占有できない)"; }
one() {   # one <case_dir> <run> <bin> <cfg_dir> <入力>
  local cd=$1 r=$2 bin=$3 src=$4 inp=$5
  waitfree
  cd $cd; mkdir $r || die "$r が既にある"
  cp $src/*.yaml $r/
  local mf; mf=$(grep -o 'meshFileName: "\?[A-Za-z0-9_.]*' $r/solverConfig.yaml | sed 's/.*: "\?//'); cp $inp $r/$mf || die "$r の入力"
  sed -i -E "s/nStepOuter: [0-9]+\}/nStepOuter: 10}/; s/outStepInterval: [0-9]+$/outStepInterval: 10/" $r/solverConfig.yaml
  ( export FORGE_BIN=$bin; bash $RC $cd/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC )
  echo "$r $(basename $(dirname $bin)) RUN_RC=$(cat $r/RUN_RC) 設定の差 $(diff $src/solverConfig.yaml $r/solverConfig.yaml | grep -c '^[<>]') 行 入力 $(sha256sum $r/$mf | cut -c1-16) ほかの forge $(pgrep -x forge | wc -l)" >> $LOG
  cd $C45
}
n=780
for tag in pt b0; do for p in f32 fp64; do for k in 1 2 3; do for t in pre post; do
  one $C45 run_0${n}_q1_${tag}_${p}_${t}${k} $BIN/$t-$p/forge $C45/_reg_src/$tag $C45/_tr/q32_init.h5; n=$((n + 1))
done; done; done; done
n=71
for p in f32 fp64; do for k in 1 2 3; do for t in pre post; do
  one $C48 run_00${n}_q3_${p}_${t}${k} $BIN/$t-$p/forge $C48/_smoke $C48/_smoke/mesh_init.h5; n=$((n + 1))
done; done; done
python3 port_reg2_an.py >> $LOG 2>&1
echo "== 終了 $(date -Is)" >> $LOG; touch port_reg2.done
