#!/bin/bash
# plan architecture-float-state-double-geometry §6.37 (事前登録、codex diagnose 2026-10-11): case/45 の float の偏りが case/48 で再現するか (line implicit 系の設定を足す)。
# §6.21 と同じバイナリ (~/forge-fgeom7-{f32,fp64})・格子・入力 (_prec/mesh_init.h5、保存量 7 つはもともと float32 で表せる値)・設定 (run_0057) に、
# case/45 の陰解法の設定 (lineImplicit 1・lineDtDirectional 1・lineDtDirectionalCap 50・implicitThermalJacobian 5・cfl 4・cfl_pseudo 4) だけを足す。
# 腕: float ×2 (run_0061・0062)、FP64 ×2 (run_0063・0064)。48,000 step・2,000 ごと、AWS A で 4 本同時。系列は c48_prec_series.py (step 0〜48,000、共通の double の座標、je を step 0 に固定した抽出も)、判定は tools/c48_li_judge.py。
set -uo pipefail
TOKEN=$(curl -s -X PUT http://169.254.169.254/latest/api/token -H "X-aws-ec2-metadata-token-ttl-seconds: 60")
[ "$(curl -s -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-id)" = "i-0b1a5e0b8dc152f00" ] || { echo "A ではない。中止"; exit 1; }
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
C48=~/forge-wallfit/case/48.flat_plate_cooled_m4
RC=~/forge-wallfit/solver_density_cuda/tools/run_case.sh
B=solver_density_cuda/build/forge; F32=~/forge-fgeom7-f32/$B; F64=~/forge-fgeom7-fp64/$B
TPL=$C48/run_0057_pv_f32a; LOG=$C48/c48_li.log
cd $C48
echo "== 開始 $(date -Is) F32 $(sha256sum $F32 | cut -c1-16) F64 $(sha256sum $F64 | cut -c1-16) 入力 $(sha256sum _prec/mesh_init.h5 | cut -c1-16)" >> $LOG
[ "$(df --output=avail -B1 . | tail -1 | tr -dc 0-9)" -ge 10737418240 ] || { echo "ディスクの空きが 10 GiB 未満 (中止)" >> $LOG; exit 1; }
arm() {   # arm <run> <バイナリ>
  local r=$1 b=$2
  mkdir $r || { echo "$r: 既にある (中止)" >> $LOG; return 1; }
  cp $TPL/{solverConfig.yaml,bcondConfig.yaml,probe.yaml} $r/ && cp _prec/mesh_init.h5 $r/mesh.h5
  sed -i 's/cfl: 2.0, cfl_pseudo: 2.0,/cfl: 4.0, cfl_pseudo: 4.0,/; s/^  deltaT: {/  deltaT: {lineImplicit: 1, lineDtDirectional: 1, lineDtDirectionalCap: 50.0, implicitThermalJacobian: 5, /' $r/solverConfig.yaml
  grep -q "lineImplicit: 1, lineDtDirectional: 1, lineDtDirectionalCap: 50.0, implicitThermalJacobian: 5," $r/solverConfig.yaml && grep -q "cfl: 4.0, cfl_pseudo: 4.0," $r/solverConfig.yaml || { echo "$r: 設定の書き換えに失敗 (中止)" >> $LOG; return 1; }
  diff $r/solverConfig.yaml $TPL/solverConfig.yaml > $r/config_vs_run_0057.diff; echo "$r: run_0057 との設定の差 $(grep -c "^[<>]" $r/config_vs_run_0057.diff) 行" >> $LOG
  ( export FORGE_BIN=$b; bash $RC $C48/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC
    python3 tools/c48_prec_series.py $r --steps 0:48000:2000 --coords-from _prec/mesh_init.h5 > $r/series.log 2>&1; rc1=$?
    FJ=$(python3 -c "import csv,sys; r=next(csv.DictReader(open(sys.argv[1]))); print(','.join(f'{x}:{int(float(r[f\"je_{x}\"]))}' for x in ('0.3','0.6','0.9')))" $r/prec_series.csv)
    python3 tools/c48_prec_series.py $r --steps 0:48000:2000 --coords-from _prec/mesh_init.h5 --fix-je $FJ --out prec_series_fixje.csv >> $r/series.log 2>&1; rc2=$?
    echo "$rc1 $rc2 je0=$FJ" > $r/SERIES_RC ) &
  echo "$r 起動 $(date -Is) $(basename $(dirname $(dirname $(dirname $b))))" >> $LOG
}
arm run_0061_li_f32a $F32 || exit 1
arm run_0062_li_f32b $F32 || exit 1
arm run_0063_li_f64a $F64 || exit 1
arm run_0064_li_f64b $F64 || exit 1
sleep 180
for r in run_0061_li_f32a run_0062_li_f32b run_0063_li_f64a run_0064_li_f64b; do echo "$r 起動の確認: $(grep -h '\[lineImplicit\] lines' $r/forge_run.log | head -1 | cut -c1-120)" >> $LOG; done
wait
for r in run_0061_li_f32a run_0062_li_f32b run_0063_li_f64a run_0064_li_f64b; do
  echo "$r 終了 $(date -Is) RUN_RC=$(cat $r/RUN_RC) SERIES_RC=$(cat $r/SERIES_RC) 最後の行: $(grep -h '^step' $r/forge_run.log | tail -1 | cut -c1-110)" >> $LOG
  for n in $(seq 2000 2000 22000); do rm -f $r/res_$n.h5 $r/res_$n.xmf; done   # 窓の前は系列に取り込んだので消す (res_0 と窓は残す)
done
python3 tools/c48_li_judge.py >> $LOG 2>&1   # §6.37 の再現判定 (旧 §6.21 の判定器は記録として中から呼ぶ)
echo "== 終了 $(date -Is)" >> $LOG
touch c48_li.done
