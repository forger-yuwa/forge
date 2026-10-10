#!/bin/bash
# plan architecture-float-state-double-geometry §6.21 (事前登録、codex diagnose 2026-10-10): case/48 の float と FP64 の比較 (精度だけを変える A/B)。
# 同じソース (commit 0b4dff4e: ~/forge-fgeom7-{f32,fp64})・同じ double の幾何の格子 (case/45/_conv614/c48/new64/m.h5)・同じ初期保存量
# (主ワークツリーの run_0025_B_tw300_y3_fx05/res_48000.h5 の 7 量を restart_field.py で 1 回だけ移した入力 HDF5 を 4 本で共有)・同じ実効設定。
# 設定は run_0048 (B 300 K、生産 SST、blockDPLUR 1、cfl_pseudo 2) から、scalarGradient を gg → lsq (node の既定) に明示し、48,000 step・2,000 step ごとに出力。
# 腕: float ×2 (run_0057・0058)、FP64 ×2 (run_0059・0060) を順に回す。各 run の後に c48_prec_series.py で判定量の系列を作る。判定は tools/c48_prec_judge.py。
set -uo pipefail
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
C48=~/forge-wallfit/case/48.flat_plate_cooled_m4; C45=~/forge-wallfit/case/45.isobutane_m6_d155
RC=~/forge-wallfit/solver_density_cuda/tools/run_case.sh; TL=~/forge-wallfit/solver_density_cuda/tools
B=solver_density_cuda/build/forge; F32=~/forge-fgeom7-f32/$B; F64=~/forge-fgeom7-fp64/$B
SEED=$C48/_seed/run_0025_res_48000.h5; TPL=$C48/run_0048_fg6_dump32; LOG=$C48/c48_prec.log
cd $C48
echo "== 開始 $(date -Is) F32 $(sha256sum $F32 | cut -c1-16) F64 $(sha256sum $F64 | cut -c1-16) seed $(sha256sum $SEED | cut -c1-16)" >> $LOG
[ "$(df --output=avail -B1 . | tail -1 | tr -dc 0-9)" -ge 5368709120 ] || { echo "ディスクの空きが 5 GiB 未満 (中止)" >> $LOG; exit 1; }
# 共通の入力 HDF5 を 1 回だけ作る
mkdir -p _prec && cp $C45/_conv614/c48/new64/m.h5 _prec/mesh_init.h5
python3 $TL/restart_field.py $SEED _prec/mesh_init.h5 --dst-run $TPL > _prec/restart.log 2>&1 || { echo "restart_field に失敗 (中止)" >> $LOG; exit 1; }
grep -q "VERDICT: OK (7 量を移した、SRC とビット一致)" _prec/restart.log || { echo "移した量が想定と違う: $(grep VERDICT _prec/restart.log) (中止)" >> $LOG; exit 1; }
echo "入力 HDF5: $(grep VERDICT _prec/restart.log) sha $(sha256sum _prec/mesh_init.h5 | cut -c1-16)" >> $LOG
arm() {   # arm <run> <バイナリ>
  local r=$1 b=$2
  mkdir $r || { echo "$r: 既にある (中止)" >> $LOG; return 1; }
  cp $TPL/{solverConfig.yaml,bcondConfig.yaml,probe.yaml} $r/ && cp _prec/mesh_init.h5 $r/mesh.h5
  sed -i 's/nStepOuter: 1}/nStepOuter: 48000}/; s/outStepInterval: 1$/outStepInterval: 2000/; s/scalarGradient: gg/scalarGradient: lsq/; s/^output: {level: 2}/output: {level: 1}/' $r/solverConfig.yaml
  grep -q 'nStepOuter: 48000}' $r/solverConfig.yaml && grep -q 'outStepInterval: 2000$' $r/solverConfig.yaml && grep -q 'scalarGradient: lsq' $r/solverConfig.yaml && grep -q '^output: {level: 1}' $r/solverConfig.yaml \
    || { echo "$r: 設定の書き換えに失敗 (中止)" >> $LOG; return 1; }
  ( export FORGE_BIN=$b; bash $RC $C48/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC )
  echo "$r 終了 $(date -Is) RUN_RC=$(cat $r/RUN_RC) 最後の行: $(grep -h '^step' $r/forge_run.log | tail -1 | cut -c1-110)" >> $LOG
  python3 tools/c48_prec_series.py $r >> $LOG 2>&1; echo "  系列 rc=$?" >> $LOG
  for n in $(seq 2000 2000 22000); do rm -f $r/res_$n.h5 $r/res_$n.xmf; done   # 窓の前は系列に取り込んだので消す (res_0 と窓は残す)
}
arm run_0057_pv_f32a $F32 || exit 1
arm run_0058_pv_f32b $F32 || exit 1
arm run_0059_pv_f64a $F64 || exit 1
arm run_0060_pv_f64b $F64 || exit 1
python3 tools/c48_prec_judge.py >> $LOG 2>&1
echo "== 終了 $(date -Is)" >> $LOG
touch c48_prec.done
