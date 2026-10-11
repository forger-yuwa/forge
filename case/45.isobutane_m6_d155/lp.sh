#!/bin/bash
# plan architecture-float-state-double-geometry §6.46 (事前登録、codex diagnose 2026-10-11): pRef の長い試験 (30,000 step、下流の θ_r・Q_w)。インスタンス A、4 本同時。
# F0 = float・pRef なし (run_0709)、Bp = float + pRef (run_0710・0711)、Ap = FP64 + pRef (run_0712)。旧 A・B・B' の実バイナリ ~/forge-fgeom7-{f32,fp64}。
# 起動の順 F0 → Bp1 → Bp2 → Ap、各起動の前に available ≥ 4.5 GB を確かめ 3 分ずつずらす。step 1〜19,500 の場は途中で消す (判定は step 0 と 20,000〜30,000)。判定は lp_an.py。
set -uo pipefail
TOKEN=$(curl -s -X PUT http://169.254.169.254/latest/api/token -H "X-aws-ec2-metadata-token-ttl-seconds: 60")
[ "$(curl -s -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-id)" = "i-0b1a5e0b8dc152f00" ] || { echo "A ではない。中止"; exit 1; }
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
export FORGE_OMEGA_BUDGET=1 FORGE_DIAG_COMMIT_LOSS=500
PREF=5499000.0
C45=~/forge-wallfit/case/45.isobutane_m6_d155; LOG=$C45/lp.log
RC=~/forge-wallfit/solver_density_cuda/tools/run_case.sh
B=solver_density_cuda/build/forge; F32=~/forge-fgeom7-f32/$B; F64=~/forge-fgeom7-fp64/$B
cd $C45
echo "== 開始 $(date -Is) F32 $(sha256sum $F32 | cut -c1-16) F64 $(sha256sum $F64 | cut -c1-16) 共通の入力 $(sha256sum _tr/q32_init.h5 | cut -c1-16) pRef $PREF" >> $LOG
[ "$(df --output=avail -B1 . | tail -1 | tr -dc 0-9)" -ge 16106127360 ] || { echo "ディスクの空きが 15 GiB 未満 (中止)" >> $LOG; exit 1; }
[ "$(sha256sum _tr/q32_init.h5 | cut -c1-8)" = "ba859e02" ] || { echo "共通の入力が §6.26 と違う (中止)" >> $LOG; exit 1; }
[ "$(sha256sum $F32 | cut -c1-16)" = "793215917d156680" ] || { echo "float のバイナリが §6.26 と違う (中止)" >> $LOG; exit 1; }
[ "$(sha256sum $F64 | cut -c1-16)" = "7e95e5bbd8f50c9e" ] || { echo "FP64 のバイナリが §6.26 と違う (中止)" >> $LOG; exit 1; }
RUNS=()
arm() {   # arm <run> <バイナリ> <pRef 0|1> <比べる旧 run>
  local r=$1 b=$2 p=$3 base=$4
  mkdir $r || { echo "$r: 既にある (中止)" >> $LOG; return 1; }
  cp _tr/*.yaml $r/ && cp _tr/q32_init.h5 $r/nozzle.h5
  sed -i 's/nStepOuter: 200000}/nStepOuter: 30000}/; s/outStepInterval: 2500/outStepInterval: 500/; s/^output: .*/output: {level: 1, extraFields: [res_ro, res_roUx, res_roUy, res_roe, res_roK, res_roOmega, omg_prod, omg_dest, omg_cross, omg_trans, omg_axisym]}/' $r/solverConfig.yaml
  [ $p = 1 ] && sed -i "s/^space: {/space: {pRef: $PREF, /" $r/solverConfig.yaml
  grep -q "nStepOuter: 30000}" $r/solverConfig.yaml && grep -q "outStepInterval: 500$" $r/solverConfig.yaml && grep -q "^mesh: {axisSegmentRWeight: 1," $r/solverConfig.yaml || { echo "$r: 設定の書き換えに失敗" >> $LOG; return 1; }
  diff $r/solverConfig.yaml $base/solverConfig.yaml > $r/config_vs_base.diff; echo "$base" > $r/CONFIG_BASE; echo "$r: $base との設定の差 $(grep -c "^[<>]" $r/config_vs_base.diff) 行" >> $LOG
  ( export FORGE_BIN=$b; bash $RC $C45/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC ) &
  echo "$r 起動 $(date -Is) $(basename $(dirname $(dirname $(dirname $b)))) pRef $p" >> $LOG; RUNS+=($r)
}
for spec in "run_0709_lp_f0 $F32 0 run_0488_tr_f32" "run_0710_lp_bp1 $F32 1 run_0488_tr_f32" "run_0711_lp_bp2 $F32 1 run_0488_tr_f32" "run_0712_lp_ap $F64 1 run_0487_tr_fp64"; do
  set -- $spec
  for w in $(seq 1 30); do av=$(free -m | awk '/^Mem:/{print $7}'); [ "$av" -ge 4500 ] && break; echo "$(date -Is) $1: available ${av} MB < 4500 MB、待つ" >> $LOG; sleep 60; done
  [ "$av" -ge 4500 ] || { echo "$1: メモリが空かない (中止、判定不能として記録)" >> $LOG; exit 1; }
  echo "$1: available ${av} MB で起動" >> $LOG
  arm $1 $2 $3 $4 || exit 1
  sleep 180
done
# 途中の res (step 1〜19,500) を消しながら待つ
while true; do
  left=0; for r in "${RUNS[@]}"; do [ -f $r/RUN_RC ] || left=1; done
  for r in "${RUNS[@]}"; do for f in $r/res_*.h5; do s=$(basename $f .h5); s=${s#res_}; [[ $s =~ ^[0-9]+$ ]] || continue
    if [ $s -gt 0 ] && [ $s -lt 20000 ]; then rm -f $f $r/res_$s.xmf; fi; done; done
  [ $left = 0 ] && break; sleep 120
done
for r in "${RUNS[@]}"; do echo "$r 終了 $(date -Is) RUN_RC=$(cat $r/RUN_RC) 最後の行: $(grep -h '^step' $r/forge_run.log | tail -1 | cut -c1-100)" >> $LOG; done
python3 lp_an.py >> $LOG 2>&1
echo "== 終了 $(date -Is)" >> $LOG
touch lp.done
