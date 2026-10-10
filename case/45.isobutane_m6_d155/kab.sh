#!/bin/bash
# plan axisymmetric-freestream-hoop-gauge §4.13 (事前登録、codex diagnose 2026-10-11): 同じ保存場・新格子・同じ FP64 バイナリから、
# mesh.axisSegmentRWeight だけを 1/0 に変えて 60,000 step 固定で回す (出力ごとの揺れのキーだけの A/B)。インスタンス B で回す。
# 入力 (A から ~/kab_in/ に送っておく): run_0483_hp7_k1 の nozzle.h5 (新格子)・res_140000.h5 (起点)・*.yaml (設定・化学種)、
#   _band_ab/cold_pair/theta_run_018{1,3}_*_100000.npz (系列の抽出の帯)。判定は kab_an.py。
set -uo pipefail
TOKEN=$(curl -s -X PUT http://169.254.169.254/latest/api/token -H "X-aws-ec2-metadata-token-ttl-seconds: 60")
[ "$(curl -s -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-id)" = "i-0ba2b91ba659254c4" ] || { echo "B ではない。中止"; exit 1; }
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
C45=~/forge-wallfit/case/45.isobutane_m6_d155; IN=~/kab_in; LOG=$C45/kab.log
RC=~/forge-wallfit/solver_density_cuda/tools/run_case.sh; TL=~/forge-wallfit/solver_density_cuda/tools
F64=~/forge-fgeom7-fp64/solver_density_cuda/build/forge
cd $C45
echo "== 開始 $(date -Is) F64 $(sha256sum $F64 | cut -c1-16) 起点 res_140000 $(sha256sum $IN/res_140000.h5 | cut -c1-16) 格子 $(sha256sum $IN/nozzle.h5 | cut -c1-16)" >> $LOG
mkdir -p _band_ab/cold_pair && cp -n $IN/theta_run_018*_100000.npz _band_ab/cold_pair/
arm() {   # arm <run> <キー>
  local r=$1 k=$2
  mkdir $r || { echo "$r: 既にある (中止)" >> $LOG; return 1; }
  cp $IN/*.yaml $r/ && cp $IN/nozzle.h5 $r/nozzle.h5
  python3 $TL/restart_field.py $IN/res_140000.h5 $r/nozzle.h5 --dst-run $r --forge $F64 > $r/restart.log 2>&1 || { echo "$r: restart_field に失敗" >> $LOG; return 1; }
  grep -q "VERDICT: OK (9 量を移した、SRC とビット一致)" $r/restart.log || { echo "$r: 移した量が想定と違う: $(grep VERDICT $r/restart.log)" >> $LOG; return 1; }
  sed -i "s/^mesh: {axisSegmentRWeight: [01],/mesh: {axisSegmentRWeight: $k,/; s/nStepOuter: 200000}/nStepOuter: 60000}/" $r/solverConfig.yaml
  grep -q "^mesh: {axisSegmentRWeight: $k," $r/solverConfig.yaml && grep -q "nStepOuter: 60000}" $r/solverConfig.yaml && grep -q "outStepInterval: 2500" $r/solverConfig.yaml || { echo "$r: 設定の書き換えに失敗" >> $LOG; return 1; }
  ( export FORGE_BIN=$F64; bash $RC $C45/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC ) &
  echo "$r 起動 $(date -Is) キー $k" >> $LOG
}
arm run_0485_kab_k1 1 || exit 1
arm run_0486_kab_k0 0 || exit 1
sleep 180
for r in run_0485_kab_k1 run_0486_kab_k0; do echo "$r 起動の確認: $(grep -h 'axisSegmentRWeight' $r/forge_run.log | head -1 | cut -c1-110)" >> $LOG; done
wait
for r in run_0485_kab_k1 run_0486_kab_k0; do echo "$r 終了 $(date -Is) RUN_RC=$(cat $r/RUN_RC) 最後の行: $(grep -h '^step' $r/forge_run.log | tail -1 | cut -c1-100)" >> $LOG; done
python3 kab_an.py >> $LOG 2>&1
echo "== 終了 $(date -Is)" >> $LOG
touch kab.done
