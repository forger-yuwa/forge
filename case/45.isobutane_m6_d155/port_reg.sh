#!/bin/bash
# plan architecture-float-state-double-geometry §5.1 #26 の回帰 (§6.59 事前登録)。インスタンス B。port_build.sh の後に回す。
# バイナリ: ~/forge-port/solver_density_cuda/build-{pre,post}-{f32,fp64}/forge (同じソースのパスでビルドした移植の前と後)。
# R1 (case/45、lowMachPrecond 0 の経路が変わらない): G0・Q32 の起点、point (run_0713_pt_d1 の設定) と B0 (run_0717_c50_d1 の設定) の各 10 step、
#    float と FP64、移植の前を 2 本・後を 1 本 (再実行の差と比べる)。run_0763〜0774。
# R2 (case/45、前処理の経路の境界の修復): point・FP64・lowMachPrecond 3・precondEps 1.0 (β = 1) の 1 step、FORGE_DIAG_STAGE_DUMP、
#    移植の前 2 本・後 2 本、移植の後の P0 (lowMachPrecond 0・precondEps 1.0) 1 本。run_0775〜0779。
# R3 (case/48、CPG・node・lowMachPrecond 0 が変わらない): _smoke の設定・入力の 10 step、float と FP64、前 2・後 1。case/48 の run_0065〜0070。
# 1 本ずつ逐次に回す (GPU は試験 G と共有)。判定は port_reg_an.py。
set -uo pipefail
TOKEN=$(curl -s -X PUT http://169.254.169.254/latest/api/token -H "X-aws-ec2-metadata-token-ttl-seconds: 60")
[ "$(curl -s -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-id)" = "i-0ba2b91ba659254c4" ] || { echo "B ではない。中止"; exit 1; }
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
export FORGE_OMEGA_BUDGET=1
C45=~/forge-wallfit/case/45.isobutane_m6_d155; C48=~/forge-wallfit/case/48.flat_plate_cooled_m4; LOG=$C45/port_reg.log
RC=~/forge-wallfit/solver_density_cuda/tools/run_case.sh; BD=~/forge-port/solver_density_cuda
die() { echo "中止: $* ($(date -Is))" >> $LOG; exit 1; }
cd $C45
[ -e port_build.done ] || die "ビルドが終わっていない"
echo "== 開始 $(date -Is)" >> $LOG
for t in pre post; do for p in f32 fp64; do echo "bin $t $p $(sha256sum $BD/build-$t-$p/forge | cut -c1-16)" >> $LOG; done; done
[ "$(sha256sum _tr/q32_init.h5 | cut -c1-8)" = "ba859e02" ] || die "case/45 の入力が §6.26 と違う"
one() {   # one <case_dir> <run> <bin> <cfg_src> <nStep> <sed 追加> <入力> [stage]
  local cd=$1 r=$2 bin=$3 cfg=$4 ns=$5 ex=$6 inp=$7 st=${8:-0}
  cd $cd; mkdir $r || die "$r が既にある"
  cp $(dirname $cfg)/*.yaml $r/ 2>/dev/null; cp $cfg $r/solverConfig.yaml
  local mf; mf=$(grep -o 'meshFileName: "\?[A-Za-z0-9_.]*' $r/solverConfig.yaml | sed 's/.*: "\?//'); cp $inp $r/$mf || die "$r の入力"
  sed -i -E "s/nStepOuter: [0-9]+\}/nStepOuter: $ns}/; s/outStepInterval: [0-9]+$/outStepInterval: $ns/" $r/solverConfig.yaml
  [ -n "$ex" ] && sed -i "$ex" $r/solverConfig.yaml
  if [ $st = 1 ]; then ( export FORGE_BIN=$bin FORGE_DIAG_STAGE_DUMP=$cd/$r/stage.h5 FORGE_DIAG_STAGE_STEP=1; bash $RC $cd/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC )
  else ( export FORGE_BIN=$bin; bash $RC $cd/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC ); fi
  echo "$r $(basename $(dirname $bin)) 終了 RUN_RC=$(cat $r/RUN_RC) $(diff $cfg $r/solverConfig.yaml | grep -c '^[<>]') 行の差 $(sha256sum $r/$mf | cut -c1-16)" >> $LOG
  cd $C45
}
n=763
for cfg in run_0713_pt_d1 run_0717_c50_d1; do for p in f32 fp64; do for t in pre pre post; do
  tag=$([ $cfg = run_0713_pt_d1 ] && echo pt || echo b0)
  one $C45 run_0${n}_r1_${tag}_${p}_${t} $BD/build-$t-$p/forge $C45/$cfg/solverConfig.yaml 10 "" $C45/_tr/q32_init.h5; n=$((n + 1))
done; done; done
for spec in "pre 3" "pre 3" "post 3" "post 3" "post 0"; do
  set -- $spec; t=$1; m=$2
  ex="s/lowMachPrecond: 0, detectNaN: 1}/lowMachPrecond: $m, precondEps: 1.0, detectNaN: 1}/"
  one $C45 run_0${n}_r2_p${m}_${t} $BD/build-$t-fp64/forge $C45/run_0713_pt_d1/solverConfig.yaml 1 "$ex" $C45/_tr/q32_init.h5 1; n=$((n + 1))
done
n=65
for p in f32 fp64; do for t in pre pre post; do
  one $C48 run_00${n}_r3_${p}_${t} $BD/build-$t-$p/forge $C48/_smoke/solverConfig.yaml 10 "" $C48/_smoke/mesh_init.h5; n=$((n + 1))
done; done
python3 port_reg_an.py >> $LOG 2>&1
echo "== 終了 $(date -Is)" >> $LOG; touch port_reg.done
