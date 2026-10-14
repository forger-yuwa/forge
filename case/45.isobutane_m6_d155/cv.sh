#!/bin/bash
# plan axisymmetric-freestream-hoop-gauge §4.14 (事前登録、forge-25 の codex diagnose の勧め): 旧い変換の格子の腕。インスタンス A で回す。
# run_0490_cv_old_k0: 旧い変換の格子 (run_0452_v4_fp64 の nozzle.h5、rSurfVect なし)、キー 0、~/forge-fgeom7-fp64、run_0483 の res_140000 から 60,000 step・2,500 ごと。
# 新しい変換の格子の腕は B の run_0486_kab_k0 (§4.13 と共用)。判定は cv_an.py。
set -uo pipefail
TOKEN=$(curl -s -X PUT http://169.254.169.254/latest/api/token -H "X-aws-ec2-metadata-token-ttl-seconds: 60")
[ "$(curl -s -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-id)" = "i-0b1a5e0b8dc152f00" ] || { echo "A ではない。中止"; exit 1; }
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
C45=~/forge-wallfit/case/45.isobutane_m6_d155; LOG=$C45/cv.log
RC=~/forge-wallfit/solver_density_cuda/tools/run_case.sh; TL=~/forge-wallfit/solver_density_cuda/tools
F64=~/forge-fgeom7-fp64/solver_density_cuda/build/forge
SRCRUN=$C45/run_0483_hp7_k1; OLDMESH=$C45/run_0452_v4_fp64/nozzle.h5; r=run_0490_cv_old_k0
cd $C45
echo "== 開始 $(date -Is) F64 $(sha256sum $F64 | cut -c1-16) 起点 res_140000 $(sha256sum $SRCRUN/res_140000.h5 | cut -c1-16) 旧格子 $(sha256sum $OLDMESH | cut -c1-16)" >> $LOG
mkdir $r || { echo "$r: 既にある (中止)" >> $LOG; exit 1; }
cp $SRCRUN/*.yaml $r/ && cp $OLDMESH $r/nozzle.h5
python3 $TL/restart_field.py $SRCRUN/res_140000.h5 $r/nozzle.h5 --dst-run $r --forge $F64 > $r/restart.log 2>&1 || { echo "$r: restart_field に失敗 (中止)" >> $LOG; exit 1; }
grep -q "VERDICT: OK (9 量を移した、SRC とビット一致)" $r/restart.log || { echo "$r: 移した量が想定と違う: $(grep VERDICT $r/restart.log) (中止)" >> $LOG; exit 1; }
sed -i 's/^mesh: {axisSegmentRWeight: 1,/mesh: {axisSegmentRWeight: 0,/; s/nStepOuter: 200000}/nStepOuter: 60000}/' $r/solverConfig.yaml
grep -q "^mesh: {axisSegmentRWeight: 0," $r/solverConfig.yaml && grep -q "nStepOuter: 60000}" $r/solverConfig.yaml && grep -q "outStepInterval: 2500" $r/solverConfig.yaml || { echo "$r: 設定の書き換えに失敗 (中止)" >> $LOG; exit 1; }
diff <(grep -v "^mesh:" $r/solverConfig.yaml) <(ssh -o BatchMode=yes -i ~/.ssh/id_a2b ubuntu@172.31.14.139 "grep -v '^mesh:' ~/forge-wallfit/case/45.isobutane_m6_d155/run_0486_kab_k0/solverConfig.yaml") > $r/config_vs_run_0486.diff 2>&1; echo "run_0486 との設定の差 (mesh 行を除く): $(wc -l < $r/config_vs_run_0486.diff) 行" >> $LOG
( export FORGE_BIN=$F64; bash $RC $C45/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC ) &
echo "$r 起動 $(date -Is)" >> $LOG
sleep 180; echo "$r 起動の確認: $(grep -h 'axisSegmentRWeight' $r/forge_run.log | head -1 | cut -c1-120)" >> $LOG
wait
echo "$r 終了 $(date -Is) RUN_RC=$(cat $r/RUN_RC) 最後の行: $(grep -h '^step' $r/forge_run.log | tail -1 | cut -c1-100)" >> $LOG
touch cv.done
