#!/bin/bash
# plan architecture-float-state-double-geometry §6.51 の 2 回目だけを回す (rcap.sh は 1 回目の後に diff の終了コードで止まった、2026-10-11)。インスタンス B。
# §6.26 の B0 の設定・旧バイナリ・入力で、lineDtDirectionalCap だけを 50 / 1。各 FP64 ×2・float ×2、計 8 本。4 本ずつ 2 回 (各回に上限 50・1 と両精度を 1 本ずつ、同じ起動順)。
# 各 30,000 step・500 ごと、出力に dt_local を足す (両方の上限で同じ)。step 1〜19,500 の場は消す。判定は rcap_an.py。
set -uo pipefail
TOKEN=$(curl -s -X PUT http://169.254.169.254/latest/api/token -H "X-aws-ec2-metadata-token-ttl-seconds: 60")
[ "$(curl -s -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-id)" = "i-0ba2b91ba659254c4" ] || { echo "B ではない。中止"; exit 1; }
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
export FORGE_OMEGA_BUDGET=1 FORGE_DIAG_COMMIT_LOSS=500
C45=~/forge-wallfit/case/45.isobutane_m6_d155; LOG=$C45/rcap.log
RC=~/forge-wallfit/solver_density_cuda/tools/run_case.sh
B=solver_density_cuda/build/forge; F32=~/forge-fgeom7-f32/$B; F64=~/forge-fgeom7-fp64/$B
cd $C45
echo "== 2 回目の開始 $(date -Is) F32 $(sha256sum $F32 | cut -c1-16) F64 $(sha256sum $F64 | cut -c1-16) 共通の入力 $(sha256sum _tr/q32_init.h5 | cut -c1-16)" >> $LOG
[ "$(df --output=avail -B1 . | tail -1 | tr -dc 0-9)" -ge 42949672960 ] || { echo "ディスクの空きが 40 GiB 未満 (中止)" >> $LOG; exit 1; }
[ "$(sha256sum _tr/q32_init.h5 | cut -c1-8)" = "ba859e02" ] || { echo "共通の入力が §6.26 と違う (中止)" >> $LOG; exit 1; }
[ "$(sha256sum $F32 | cut -c1-16)" = "793215917d156680" ] || { echo "float のバイナリが §6.26 と違う (中止)" >> $LOG; exit 1; }
[ "$(sha256sum $F64 | cut -c1-16)" = "7e95e5bbd8f50c9e" ] || { echo "FP64 のバイナリが §6.26 と違う (中止)" >> $LOG; exit 1; }
arm() {   # arm <run> <バイナリ> <上限 50|1>
  local r=$1 b=$2 cap=$3
  mkdir $r || { echo "$r: 既にある (中止)" >> $LOG; return 1; }
  cp _tr/*.yaml $r/ && cp _tr/q32_init.h5 $r/nozzle.h5
  sed -i 's/nStepOuter: 200000}/nStepOuter: 30000}/; s/outStepInterval: 2500/outStepInterval: 500/; s/^output: .*/output: {level: 1, extraFields: [res_ro, res_roUx, res_roUy, res_roe, res_roK, res_roOmega, omg_prod, omg_dest, omg_cross, omg_trans, omg_axisym, dt_local]}/' $r/solverConfig.yaml
  sed -i "s/lineDtDirectionalCap: 50.0, /lineDtDirectionalCap: ${cap}.0, /" $r/solverConfig.yaml
  grep -q "lineDtDirectionalCap: ${cap}.0, implicitThermalJacobian: 5, lineImplicit: 1, lineDtDirectional: 1," $r/solverConfig.yaml && grep -q "nStepOuter: 30000}" $r/solverConfig.yaml && grep -q "dt_local\]" $r/solverConfig.yaml || { echo "$r: 設定の書き換えに失敗" >> $LOG; return 1; }
  ( export FORGE_BIN=$b; bash $RC $C45/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC ) &
  echo "$r 起動 $(date -Is) $(basename $(dirname $(dirname $(dirname $b)))) 上限 $cap" >> $LOG
}
batch() {   # batch <run50_d> <run1_d> <run50_f> <run1_f>
  local RUNS=("$@")
  arm $1 $F64 50 && sleep 120 && arm $2 $F64 1 && sleep 120 && arm $3 $F32 50 && sleep 120 && arm $4 $F32 1 || return 1
  while true; do
    left=0; for r in "${RUNS[@]}"; do [ -f $r/RUN_RC ] || left=1; done
    for r in "${RUNS[@]}"; do for f in $r/res_*.h5; do s=$(basename $f .h5); s=${s#res_}; [[ $s =~ ^[0-9]+$ ]] || continue
      if [ $s -gt 0 ] && [ $s -lt 20000 ]; then rm -f $f $r/res_$s.xmf; fi; done; done
    [ $left = 0 ] && break; sleep 120
  done
  for r in "${RUNS[@]}"; do echo "$r 終了 $(date -Is) RUN_RC=$(cat $r/RUN_RC) 最後の行: $(grep -h '^step' $r/forge_run.log | tail -1 | cut -c1-100)" >> $LOG; done
  # 次の回の前の dt の照合用に、50 と 1 の設定の差を記録
  diff $2/solverConfig.yaml $1/solverConfig.yaml > $2/config_vs_cap50.diff; true   # diff は差があると 1 を返す (2026-10-11: これで 1 回目の後に台本が止まった)
}
av=$(free -m | awk '/^Mem:/{print $7}'); echo "回 2 (rcap2.sh): available ${av} MB" >> $LOG
[ "$av" -ge 11000 ] || { echo "回 2: メモリが 11 GB 未満 (中止)" >> $LOG; exit 1; }
batch run_0721_c50_d2 run_0722_c1_d2 run_0723_c50_f2 run_0724_c1_f2 || { echo "回 2: batch が失敗 (中止)" >> $LOG; exit 1; }
python3 rcap_an.py >> $LOG 2>&1
echo "== 終了 $(date -Is)" >> $LOG
touch rcap.done
