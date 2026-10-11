#!/bin/bash
# plan architecture-float-state-double-geometry §6.48 (事前登録、codex diagnose 2026-10-11): point の陰解法に固定した float と FP64 の比較 (案 P)。インスタンス B、4 本同時。
# FP64 ×2 (run_0713・0714)、float ×2 (run_0715・0716)。§6.26 の設定から lineImplicit・lineDtDirectional・lineDtDirectionalCap・implicitThermalJacobian の 4 キーを除く。
# 最大 120,000 step。30,000 step で段 1 の判定 (rpt_an.py --stage 1)、結果 A なら自分の forge を止めて終える。そうでなければ 120,000 step で段 2 の判定。
# 判定に使わない場 (step 1〜19,500、30,500〜99,500) は途中で消す。
set -uo pipefail
TOKEN=$(curl -s -X PUT http://169.254.169.254/latest/api/token -H "X-aws-ec2-metadata-token-ttl-seconds: 60")
[ "$(curl -s -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-id)" = "i-0ba2b91ba659254c4" ] || { echo "B ではない。中止"; exit 1; }
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
export FORGE_OMEGA_BUDGET=1 FORGE_DIAG_COMMIT_LOSS=500
C45=~/forge-wallfit/case/45.isobutane_m6_d155; LOG=$C45/rpt.log
RC=~/forge-wallfit/solver_density_cuda/tools/run_case.sh
B=solver_density_cuda/build/forge; F32=~/forge-fgeom7-f32/$B; F64=~/forge-fgeom7-fp64/$B
cd $C45
echo "== 開始 $(date -Is) F32 $(sha256sum $F32 | cut -c1-16) F64 $(sha256sum $F64 | cut -c1-16) 共通の入力 $(sha256sum _tr/q32_init.h5 | cut -c1-16)" >> $LOG
[ "$(df --output=avail -B1 . | tail -1 | tr -dc 0-9)" -ge 53687091200 ] || { echo "ディスクの空きが 50 GiB 未満 (中止)" >> $LOG; exit 1; }
[ "$(sha256sum _tr/q32_init.h5 | cut -c1-8)" = "ba859e02" ] || { echo "共通の入力が §6.26 と違う (中止)" >> $LOG; exit 1; }
[ "$(sha256sum $F32 | cut -c1-16)" = "793215917d156680" ] || { echo "float のバイナリが §6.26 と違う (中止)" >> $LOG; exit 1; }
[ "$(sha256sum $F64 | cut -c1-16)" = "7e95e5bbd8f50c9e" ] || { echo "FP64 のバイナリが §6.26 と違う (中止)" >> $LOG; exit 1; }
RUNS=()
arm() {   # arm <run> <バイナリ> <比べる旧 run>
  local r=$1 b=$2 base=$3
  mkdir $r || { echo "$r: 既にある (中止)" >> $LOG; return 1; }
  cp _tr/*.yaml $r/ && cp _tr/q32_init.h5 $r/nozzle.h5
  sed -i 's/nStepOuter: 200000}/nStepOuter: 120000}/; s/outStepInterval: 2500/outStepInterval: 500/; s/^output: .*/output: {level: 1, extraFields: [res_ro, res_roUx, res_roUy, res_roe, res_roK, res_roOmega, omg_prod, omg_dest, omg_cross, omg_trans, omg_axisym]}/' $r/solverConfig.yaml
  sed -i 's/lineDtDirectionalCap: 50.0, //; s/implicitThermalJacobian: 5, //; s/lineImplicit: 1, //; s/lineDtDirectional: 1, //' $r/solverConfig.yaml
  if grep -q "lineImplicit\|lineDtDirectional\|implicitThermalJacobian" $r/solverConfig.yaml; then echo "$r: point への書き換えに失敗 (中止)" >> $LOG; return 1; fi
  grep -q "nStepOuter: 120000}" $r/solverConfig.yaml && grep -q "outStepInterval: 500$" $r/solverConfig.yaml && grep -q "^mesh: {axisSegmentRWeight: 1," $r/solverConfig.yaml || { echo "$r: 設定の書き換えに失敗" >> $LOG; return 1; }
  diff $r/solverConfig.yaml $base/solverConfig.yaml > $r/config_vs_base.diff; echo "$base" > $r/CONFIG_BASE; echo "$r: $base との設定の差 $(grep -c "^[<>]" $r/config_vs_base.diff) 行" >> $LOG
  ( export FORGE_BIN=$b; bash $RC $C45/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC ) &
  echo "$r 起動 $(date -Is) $(basename $(dirname $(dirname $(dirname $b))))" >> $LOG; RUNS+=($r)
}
for spec in "run_0713_pt_d1 $F64 run_0487_tr_fp64" "run_0715_pt_f1 $F32 run_0488_tr_f32" "run_0714_pt_d2 $F64 run_0487_tr_fp64" "run_0716_pt_f2 $F32 run_0488_tr_f32"; do
  set -- $spec
  for w in $(seq 1 30); do av=$(free -m | awk '/^Mem:/{print $7}'); [ "$av" -ge 4500 ] && break; echo "$(date -Is) $1: available ${av} MB < 4500 MB、待つ" >> $LOG; sleep 60; done
  [ "$av" -ge 4500 ] || { echo "$1: メモリが空かない (中止)" >> $LOG; exit 1; }
  echo "$1: available ${av} MB で起動" >> $LOG
  arm $1 $2 $3 || exit 1
  sleep 120
done
keep() { local s=$1; [ $s -eq 0 ] && return 0; [ $s -ge 20000 ] && [ $s -le 30000 ] && return 0; [ $s -ge 100000 ] && return 0; return 1; }
stage1=0
while true; do
  left=0; for r in "${RUNS[@]}"; do [ -f $r/RUN_RC ] || left=1; done
  for r in "${RUNS[@]}"; do for f in $r/res_*.h5; do s=$(basename $f .h5); s=${s#res_}; [[ $s =~ ^[0-9]+$ ]] || continue
    keep $s || rm -f $f $r/res_$s.xmf; done; done
  if [ $stage1 = 0 ]; then
    all=1; for r in "${RUNS[@]}"; do [ -f $r/res_30000.h5 ] || all=0; done
    if [ $all = 1 ]; then
      sleep 60   # res_30000 の書き込みの完了を待つ
      echo "== 段 1 の判定 $(date -Is)" >> $LOG
      python3 rpt_an.py --stage 1 >> $LOG 2>&1; stage1=1
      if grep -q "== VERDICT §6.48 段 1: 結果 A" $LOG; then
        echo "段 1 で結果 A: 自分の run を止める $(date -Is)" >> $LOG
        for p in $(pgrep -x forge); do c=$(readlink /proc/$p/cwd); for r in "${RUNS[@]}"; do [ "$c" = "$C45/$r" ] && kill $p && echo "  止めた $p $r" >> $LOG; done; done
        touch rpt_stopped_after_stage1
      fi
    fi
  fi
  [ $left = 0 ] && break; sleep 120
done
for r in "${RUNS[@]}"; do echo "$r 終了 $(date -Is) RUN_RC=$(cat $r/RUN_RC) 最後の行: $(grep -h '^step' $r/forge_run.log | tail -1 | cut -c1-100)" >> $LOG; done
[ -f rpt_stopped_after_stage1 ] || python3 rpt_an.py --stage 2 >> $LOG 2>&1
echo "== 終了 $(date -Is)" >> $LOG
touch rpt.done
