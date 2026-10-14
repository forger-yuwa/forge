#!/bin/bash
# plan architecture-float-state-double-geometry §6.55 (事前登録、codex diagnose 2026-10-12): G0c・B0・継続場から space.pRef 0 / 5,499,000 を両精度で比べる。インスタンス B。
# 1. 関門: 診断 run_0751_pc_hoop (FP64 + pRef、1 step、FORGE_DIAG_HOOP_CLOSURE) の hoop.h5 で G0c の閉性 g ≤ 1e-12 (rpc_gate.py)。通らなければ本番は回さない。
# 2. 4 本ずつ 2 回 (各回に C・Cp・F・Fp を 1 本ずつ、この順に 30 秒ずつ)、各 3,000 step・250 ごと。バイナリは §6.54 と同じ fgeom11。判定は rpc_an.py。
set -uo pipefail
TOKEN=$(curl -s -X PUT http://169.254.169.254/latest/api/token -H "X-aws-ec2-metadata-token-ttl-seconds: 60")
[ "$(curl -s -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-id)" = "i-0ba2b91ba659254c4" ] || { echo "B ではない。中止"; exit 1; }
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
export FORGE_OMEGA_BUDGET=1 FORGE_DIAG_COMMIT_LOSS=500
C45=~/forge-wallfit/case/45.isobutane_m6_d155; LOG=$C45/rpc.log; TL=~/forge-wallfit/solver_density_cuda/tools
RC=$TL/run_case.sh
B=solver_density_cuda/build/forge; F32=~/forge-fgeom11-f32/$B; F64=~/forge-fgeom11-fp64/$B
PREF=5499000.0
cd $C45
die() { echo "中止: $* ($(date -Is))" >> $LOG; exit 1; }
echo "== 開始 $(date -Is) F32 $(sha256sum $F32 | cut -c1-16) F64 $(sha256sum $F64 | cut -c1-16) 継続場 $(sha256sum _ic/cont_q32.h5 | cut -c1-16) pRef $PREF" >> $LOG
[ "$(sha256sum $F32 | cut -c1-16)" = "1f9ff32232da8823" ] || die "float のビルドが §6.54 と違う"
[ "$(sha256sum $F64 | cut -c1-16)" = "74954e7be5ef38a8" ] || die "FP64 のビルドが §6.54 と違う"
[ "$(sha256sum _ic/cont_q32.h5 | cut -c1-16)" = "3baafa3d23d31e5e" ] || die "継続場の入力が登録と違う"
[ "$(df --output=avail -B1 . | tail -1 | tr -dc 0-9)" -ge 32212254720 ] || die "ディスクの空きが 30 GiB 未満"
mkcfg() {   # mkcfg <run> <nStep> <pRef 0|1>
  local r=$1 n=$2 p=$3
  cp _g0c/*.yaml $r/ && cp _ic/cont_q32.h5 $r/nozzle.h5
  sed -i "s/nStepOuter: 200000}/nStepOuter: $n}/; s/outStepInterval: 2500/outStepInterval: 250/; s/^output: .*/output: {level: 1, extraFields: [res_ro, res_roUx, res_roUy, res_roe, res_roK, res_roOmega, omg_prod, omg_dest, omg_cross, omg_trans, omg_axisym, dt_local]}/" $r/solverConfig.yaml
  [ $p = 1 ] && sed -i "s/^space: {/space: {pRef: $PREF, /" $r/solverConfig.yaml
  local want; want=$(sed "s/nStepOuter: 30000}/nStepOuter: $n}/; s/outStepInterval: 500/outStepInterval: 250/" run_0717_c50_d1/solverConfig.yaml)
  [ $p = 1 ] && want=$(echo "$want" | sed "s/^space: {/space: {pRef: $PREF, /")
  [ "$want" = "$(cat $r/solverConfig.yaml)" ] || { echo "$r: 設定が登録と同一でない" >> $LOG; return 1; }
}
# 1. 関門
r=run_0751_pc_hoop; mkdir $r || die "$r が既にある"; mkcfg $r 1 1 || die "$r の設定"
( export FORGE_BIN=$F64 FORGE_DIAG_HOOP_CLOSURE=$C45/$r/hoop.h5; bash $RC $C45/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC )
echo "$r 終了 RUN_RC=$(cat $r/RUN_RC) $(grep -h 'hoopClosure' $r/forge_run.log | head -1 | cut -c1-200)" >> $LOG
python3 rpc_gate.py >> $LOG 2>&1 || die "閉性の関門を通らない (本番は回さない)"
arm() {   # arm <run> <バイナリ> <pRef 0|1>
  local r=$1 b=$2 p=$3
  mkdir $r || { echo "$r: 既にある (中止)" >> $LOG; return 1; }
  mkcfg $r 3000 $p || return 1
  ( export FORGE_BIN=$b; bash $RC $C45/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC ) &
  echo "$r 起動 $(date -Is) $(basename $(dirname $(dirname $(dirname $b)))) pRef $p 入力 $(sha256sum $r/nozzle.h5 | cut -c1-16)" >> $LOG
}
batch() {   # batch <C> <Cp> <F> <Fp>
  local RUNS=("$@")
  for i in $(seq 1 30); do av=$(free -m | awk '/^Mem:/{print $7}'); [ "$av" -ge 10000 ] && break; sleep 60; done
  echo "available ${av} MB" >> $LOG; [ "$av" -ge 10000 ] || return 1
  arm $1 $F64 0 && sleep 30 && arm $2 $F64 1 && sleep 30 && arm $3 $F32 0 && sleep 30 && arm $4 $F32 1 || return 1
  while true; do left=0; for r in "${RUNS[@]}"; do [ -f $r/RUN_RC ] || left=1; done; [ $left = 0 ] && break; sleep 30; done
  for r in "${RUNS[@]}"; do echo "$r 終了 $(date -Is) RUN_RC=$(cat $r/RUN_RC) 最後の行: $(grep -h '^step' $r/forge_run.log | tail -1 | cut -c1-100)" >> $LOG; done
  for r in $2 $4; do grep -h "pRef" $r/forge_run.log | head -1 | cut -c1-200 | sed "s/^/$r: /" >> $LOG; done
}
batch run_0743_pc_c1 run_0744_pc_cp1 run_0745_pc_f1 run_0746_pc_fp1 || die "1 回目"
batch run_0747_pc_c2 run_0748_pc_cp2 run_0749_pc_f2 run_0750_pc_fp2 || die "2 回目"
python3 rpc_an.py >> $LOG 2>&1
echo "== 終了 $(date -Is)" >> $LOG
touch rpc.done
