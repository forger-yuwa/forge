#!/bin/bash
# plan architecture-float-state-double-geometry §6.40 (事前登録、codex diagnose 2026-10-11): FP64 のビルドで ρE だけを step 末尾に float32 に丸める短い試験。インスタンス B。
# 対照 C = FP64 OFF (run_0595〜0597)、介入 E = FP64 + ρE だけ量子化 (run_0598・0599・0620)、陽性対照 F = float (run_0621〜0623)。各 2,000 step・500 ごと。
# C・E のバイナリ ~/forge-fgeom11-fp64 (#20 の commit、rq_build.sh)、F は ~/forge-fgeom7-f32 (§6.26 の B・B′ と同じ)。1 回に C_i・E_i・F_i の 3 本、3 回に分ける。判定は rq_an.py。
set -uo pipefail
TOKEN=$(curl -s -X PUT http://169.254.169.254/latest/api/token -H "X-aws-ec2-metadata-token-ttl-seconds: 60")
[ "$(curl -s -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-id)" = "i-0ba2b91ba659254c4" ] || { echo "B ではない。中止"; exit 1; }
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
export FORGE_OMEGA_BUDGET=1 FORGE_DIAG_COMMIT_LOSS=500
C45=~/forge-wallfit/case/45.isobutane_m6_d155; LOG=$C45/rq.log
RC=~/forge-wallfit/solver_density_cuda/tools/run_case.sh
B=solver_density_cuda/build/forge; F64=~/forge-fgeom11-fp64/$B; F32=~/forge-fgeom7-f32/$B
cd $C45
echo "== 開始 $(date -Is) F64 $(sha256sum $F64 | cut -c1-16) F32 $(sha256sum $F32 | cut -c1-16) 共通の入力 $(sha256sum _tr/q32_init.h5 | cut -c1-16)" >> $LOG
[ "$(df --output=avail -B1 . | tail -1 | tr -dc 0-9)" -ge 21474836480 ] || { echo "ディスクの空きが 20 GiB 未満 (中止)" >> $LOG; exit 1; }
[ "$(sha256sum _tr/q32_init.h5 | cut -c1-8)" = "ba859e02" ] || { echo "共通の入力が §6.26 と違う (中止)" >> $LOG; exit 1; }
[ "$(sha256sum $F32 | cut -c1-16)" = "793215917d156680" ] || { echo "float のバイナリが §6.26 と違う (中止)" >> $LOG; exit 1; }
arm() {   # arm <run> <バイナリ> <量子化 0|1>
  local r=$1 b=$2 q=$3
  mkdir $r || { echo "$r: 既にある (中止)" >> $LOG; return 1; }
  cp _tr/*.yaml $r/ && cp _tr/q32_init.h5 $r/nozzle.h5
  sed -i 's/nStepOuter: 200000}/nStepOuter: 2000}/; s/outStepInterval: 2500/outStepInterval: 500/; s/^output: .*/output: {level: 1, extraFields: [res_ro, res_roUx, res_roUy, res_roe, res_roK, res_roOmega, omg_prod, omg_dest, omg_cross, omg_trans, omg_axisym]}/' $r/solverConfig.yaml
  grep -q "nStepOuter: 2000}" $r/solverConfig.yaml && grep -q "outStepInterval: 500$" $r/solverConfig.yaml && grep -q "^mesh: {axisSegmentRWeight: 1," $r/solverConfig.yaml || { echo "$r: 設定の書き換えに失敗" >> $LOG; return 1; }
  diff $r/solverConfig.yaml run_0488_tr_f32/solverConfig.yaml > $r/config_vs_run_0488.diff; echo "$r: run_0488 との設定の差 $(grep -c "^[<>]" $r/config_vs_run_0488.diff) 行 (nStepOuter の 1 行のはず)" >> $LOG
  if [ $q = 1 ]; then ( export FORGE_BIN=$b FORGE_DIAG_QUANTIZE_STATE=1 FORGE_DIAG_QUANTIZE_STATE_ONLY=roe FORGE_DIAG_QUANTIZE_STATE_CHECK=1,2000; bash $RC $C45/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC ) &
  else ( export FORGE_BIN=$b; bash $RC $C45/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC ) &
  fi
  echo "$r 起動 $(date -Is) $(basename $(dirname $(dirname $(dirname $b)))) 量子化 $q" >> $LOG
}
for k in 1 2 3; do
  case $k in 1) c=run_0595_rq_c1; e=run_0598_rq_e1; f=run_0621_rq_f1;; 2) c=run_0596_rq_c2; e=run_0599_rq_e2; f=run_0622_rq_f2;; 3) c=run_0597_rq_c3; e=run_0620_rq_e3; f=run_0623_rq_f3;; esac
  # メモリ: g5 は 16 GB・スワップなし。3 本ぶん (FP64 約 2.5 GB × 2 + float 約 2 GB) + 2 GB = 9 GB の available を 30 分まで待つ
  for w in $(seq 1 30); do av=$(free -m | awk '/^Mem:/{print $7}'); [ "$av" -ge 9000 ] && break; echo "$(date -Is) 回 $k: available ${av} MB < 9000 MB、待つ" >> $LOG; sleep 60; done
  [ "$av" -ge 9000 ] || { echo "回 $k: メモリが空かない (中止)" >> $LOG; exit 1; }
  echo "回 $k: available ${av} MB で投入" >> $LOG
  arm $c $F64 0 || exit 1; arm $e $F64 1 || exit 1; arm $f $F32 0 || exit 1
  wait
  for r in $c $e $f; do echo "$r 終了 $(date -Is) RUN_RC=$(cat $r/RUN_RC) 最後の行: $(grep -h '^step' $r/forge_run.log | tail -1 | cut -c1-100)" >> $LOG; done
  grep -h "\[quantizeState\]" $e/forge_run.log | cut -c1-300 >> $LOG
done
python3 rq_an.py >> $LOG 2>&1
echo "== 終了 $(date -Is)" >> $LOG
touch rq.done
