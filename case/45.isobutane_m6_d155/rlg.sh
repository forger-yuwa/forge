#!/bin/bash
# plan architecture-float-state-double-geometry §6.57 (事前登録、codex diagnose 2026-10-12): 前処理の経路への切り替えが node の境界の更新を壊さないかの 1 step の確認。
# G0・Q32 の起点 (_tr/q32_init.h5)・point (§6.48 の run_0713_pt_d1 の設定)・FP64 (~/forge-fgeom11-fp64、FORGE_DIAG_STAGE_DUMP を持つ)。
# 両腕とも precondEps 1.0 (β = 1)、変えるのは lowMachPrecond 0 / 3 だけ。各 2 本・1 step、1 本ずつ (GPU を占有して逐次)。判定は rlg_an.py。
set -uo pipefail
TOKEN=$(curl -s -X PUT http://169.254.169.254/latest/api/token -H "X-aws-ec2-metadata-token-ttl-seconds: 60")
[ "$(curl -s -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-id)" = "i-0ba2b91ba659254c4" ] || { echo "B ではない。中止"; exit 1; }
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
export FORGE_OMEGA_BUDGET=1
C45=~/forge-wallfit/case/45.isobutane_m6_d155; LOG=$C45/rlg.log; RC=~/forge-wallfit/solver_density_cuda/tools/run_case.sh
F64=~/forge-fgeom11-fp64/solver_density_cuda/build/forge
cd $C45
die() { echo "中止: $* ($(date -Is))" >> $LOG; exit 1; }
echo "== 開始 $(date -Is) F64 $(sha256sum $F64 | cut -c1-16) 入力 $(sha256sum _tr/q32_init.h5 | cut -c1-16)" >> $LOG
[ "$(sha256sum $F64 | cut -c1-16)" = "74954e7be5ef38a8" ] || die "FP64 のビルドが登録 (fgeom11) と違う"
[ "$(sha256sum _tr/q32_init.h5 | cut -c1-8)" = "ba859e02" ] || die "共通の入力が §6.26 と違う"
grep -q "lowMachPrecond: 0, detectNaN: 1}" run_0713_pt_d1/solverConfig.yaml || die "point の設定の形が想定と違う"
[ -z "$(pgrep -x forge)" ] || die "ほかの forge が走っている (占有して逐次に回す)"
for spec in "run_0756_lg_p0a 0" "run_0757_lg_p0b 0" "run_0758_lg_p3a 3" "run_0759_lg_p3b 3"; do
  set -- $spec; r=$1; m=$2
  mkdir $r || die "$r が既にある"
  cp _tr/*.yaml $r/ && cp run_0713_pt_d1/solverConfig.yaml $r/solverConfig.yaml && cp _tr/q32_init.h5 $r/nozzle.h5 || die "$r の複製"
  sed -i "s/nStepOuter: 120000}/nStepOuter: 1}/; s/lowMachPrecond: 0, detectNaN: 1}/lowMachPrecond: $m, precondEps: 1.0, detectNaN: 1}/" $r/solverConfig.yaml
  [ "$(diff run_0713_pt_d1/solverConfig.yaml $r/solverConfig.yaml | grep -c '^[<>]')" = 4 ] && grep -q "lowMachPrecond: $m, precondEps: 1.0," $r/solverConfig.yaml && grep -q "nStepOuter: 1}" $r/solverConfig.yaml || die "$r の設定 (2 行だけのはず)"
  ( export FORGE_BIN=$F64 FORGE_DIAG_STAGE_DUMP=$C45/$r/stage.h5 FORGE_DIAG_STAGE_STEP=1; bash $RC $C45/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC )
  echo "$r lowMachPrecond $m 終了 $(date -Is) RUN_RC=$(cat $r/RUN_RC) stage.h5 $([ -e $r/stage.h5 ] && echo あり || echo なし) $(grep -h 'stageDump\|未対応' $r/forge_run.log | head -2 | cut -c1-160)" >> $LOG
done
python3 rlg_an.py >> $LOG 2>&1
echo "== 終了 $(date -Is)" >> $LOG; touch rlg.done
