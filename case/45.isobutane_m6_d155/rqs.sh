#!/bin/bash
# plan architecture-float-state-double-geometry §6.54 (事前登録、codex diagnose 2026-10-12): G0c・B0・継続場 `_ic/cont_q32.h5` から、
# 同じ FP64 のビルド ~/forge-fgeom11-fp64 で step 末尾の 9 量の量子化を OFF (C) / ON (Q)、陽性対照は同じソースの float ~/forge-fgeom11-f32 (F)。
# 3 本ずつ 2 回 (各回に C・Q・F を 1 本ずつ、同じ起動順、30 秒ずつ)、各 3,000 step・250 ごと。設定は run_0717_c50_d1 から nStepOuter・outStepInterval だけを変える。判定は rqs_an.py。
set -uo pipefail
TOKEN=$(curl -s -X PUT http://169.254.169.254/latest/api/token -H "X-aws-ec2-metadata-token-ttl-seconds: 60")
[ "$(curl -s -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-id)" = "i-0ba2b91ba659254c4" ] || { echo "B ではない。中止"; exit 1; }
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
export FORGE_OMEGA_BUDGET=1 FORGE_DIAG_COMMIT_LOSS=500
C45=~/forge-wallfit/case/45.isobutane_m6_d155; LOG=$C45/rqs.log; TL=~/forge-wallfit/solver_density_cuda/tools
RC=$TL/run_case.sh
B=solver_density_cuda/build/forge; F32=~/forge-fgeom11-f32/$B; F64=~/forge-fgeom11-fp64/$B
cd $C45
die() { echo "中止: $* ($(date -Is))" >> $LOG; exit 1; }
[ -e rqc_build.done ] || die "float のビルド (rqc_build.sh) が終わっていない"
echo "== 開始 $(date -Is) F32 $(sha256sum $F32 | cut -c1-16) F64 $(sha256sum $F64 | cut -c1-16) 継続場 $(sha256sum _ic/cont_q32.h5 | cut -c1-16)" >> $LOG
[ "$(sha256sum $F64 | cut -c1-16)" = "74954e7be5ef38a8" ] || die "FP64 のビルドが登録 (fgeom11) と違う"
[ "$(sha256sum _ic/cont_q32.h5 | cut -c1-16)" = "3baafa3d23d31e5e" ] || die "継続場の入力が登録と違う"
[ "$(df --output=avail -B1 . | tail -1 | tr -dc 0-9)" -ge 21474836480 ] || die "ディスクの空きが 20 GiB 未満"
arm() {   # arm <run> <バイナリ> <量子化 0|1>
  local r=$1 b=$2 q=$3
  mkdir $r || { echo "$r: 既にある (中止)" >> $LOG; return 1; }
  cp _g0c/*.yaml $r/ && cp _ic/cont_q32.h5 $r/nozzle.h5
  sed -i 's/nStepOuter: 200000}/nStepOuter: 3000}/; s/outStepInterval: 2500/outStepInterval: 250/; s/^output: .*/output: {level: 1, extraFields: [res_ro, res_roUx, res_roUy, res_roe, res_roK, res_roOmega, omg_prod, omg_dest, omg_cross, omg_trans, omg_axisym, dt_local]}/' $r/solverConfig.yaml
  sed 's/nStepOuter: 30000}/nStepOuter: 3000}/; s/outStepInterval: 500/outStepInterval: 250/' run_0717_c50_d1/solverConfig.yaml | cmp -s - $r/solverConfig.yaml || { echo "$r: 設定が登録 (run_0717 から 2 行だけ変更) と同一でない" >> $LOG; return 1; }
  if [ $q = 1 ]; then ( export FORGE_BIN=$b FORGE_DIAG_QUANTIZE_STATE=1 FORGE_DIAG_QUANTIZE_STATE_CHECK=1,2750,3000; bash $RC $C45/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC ) &
  else ( export FORGE_BIN=$b; bash $RC $C45/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC ) &
  fi
  echo "$r 起動 $(date -Is) $(basename $(dirname $(dirname $(dirname $b)))) 量子化 $q 入力 $(sha256sum $r/nozzle.h5 | cut -c1-16)" >> $LOG
}
batch() {   # batch <C> <Q> <F>
  local RUNS=("$@")
  for i in $(seq 1 30); do av=$(free -m | awk '/^Mem:/{print $7}'); [ "$av" -ge 8000 ] && break; sleep 60; done
  echo "available ${av} MB" >> $LOG; [ "$av" -ge 8000 ] || return 1
  arm $1 $F64 0 && sleep 30 && arm $2 $F64 1 && sleep 30 && arm $3 $F32 0 || return 1
  while true; do left=0; for r in "${RUNS[@]}"; do [ -f $r/RUN_RC ] || left=1; done; [ $left = 0 ] && break; sleep 30; done
  for r in "${RUNS[@]}"; do echo "$r 終了 $(date -Is) RUN_RC=$(cat $r/RUN_RC) 最後の行: $(grep -h '^step' $r/forge_run.log | tail -1 | cut -c1-100)" >> $LOG; done
}
batch run_0737_qs_c1 run_0738_qs_q1 run_0739_qs_f1 || die "1 回目"
batch run_0740_qs_c2 run_0741_qs_q2 run_0742_qs_f2 || die "2 回目"
python3 rqs_an.py >> $LOG 2>&1
echo "== 終了 $(date -Is)" >> $LOG
touch rqs.done
