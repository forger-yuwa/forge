#!/bin/bash
# plan architecture-float-state-double-geometry §6.52 (事前登録、codex diagnose 2026-10-11): 縮流部を粗くした格子 G0c で、B0 のまま float と FP64 を比べる。インスタンス B。
# §6.51 (rcap2.sh) の終了 (rcap.done) を待ち、上限 50 の陽性対照が成立したら (rg0c_an.py --pc-only)、G0c の共通の入力 _g0c/q32_init.h5 から
# §6.51 の上限 50 と同じ設定 (B0、出力に dt_local) で FP64 ×2・float ×2 を同時に 4 本、各 30,000 step・500 ごと。step 1〜19,500 の場は消す。判定は rg0c_an.py。
set -uo pipefail
TOKEN=$(curl -s -X PUT http://169.254.169.254/latest/api/token -H "X-aws-ec2-metadata-token-ttl-seconds: 60")
[ "$(curl -s -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-id)" = "i-0ba2b91ba659254c4" ] || { echo "B ではない。中止"; exit 1; }
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
export FORGE_OMEGA_BUDGET=1 FORGE_DIAG_COMMIT_LOSS=500
C45=~/forge-wallfit/case/45.isobutane_m6_d155; LOG=$C45/rg0c.log
RC=~/forge-wallfit/solver_density_cuda/tools/run_case.sh
B=solver_density_cuda/build/forge; F32=~/forge-fgeom7-f32/$B; F64=~/forge-fgeom7-fp64/$B
cd $C45
die() { echo "中止: $* ($(date -Is))" >> $LOG; exit 1; }
echo "== 待機の開始 $(date -Is)" >> $LOG
# 1. §6.51 の終了を待つ (最大 4 時間)
for i in $(seq 1 120); do
  [ -e rcap.done ] && break
  grep -q "中止" rcap.log && die "§6.51 の台本が中止した"
  sleep 120
done
[ -e rcap.done ] || die "§6.51 が 4 時間で終わらない"
for r in run_0717_c50_d1 run_0718_c1_d1 run_0719_c50_f1 run_0720_c1_f1 run_0721_c50_d2 run_0722_c1_d2 run_0723_c50_f2 run_0724_c1_f2; do
  [ "$(cat $r/RUN_RC 2>/dev/null)" = 0 ] || die "$r が正常に終わっていない"
done
# 2. 陽性対照
echo "== 陽性対照 (rg0c_an.py --pc-only) $(date -Is)" >> $LOG
python3 rg0c_an.py --pc-only >> $LOG 2>&1 || die "§6.51 の上限 50 の陽性対照が成立しない (G0c の試験は投入しない)"
# 3. 入力・バイナリ・準備
echo "== 開始 $(date -Is) F32 $(sha256sum $F32 | cut -c1-16) F64 $(sha256sum $F64 | cut -c1-16) 共通の入力 $(sha256sum _g0c/q32_init.h5 | cut -c1-16)" >> $LOG
[ -e g0c_prep.done ] || die "G0c の準備が終わっていない"
[ "$(sha256sum _g0c/q32_init.h5 | cut -c1-16)" = "6f8fcda6d97482e3" ] || die "G0c の入力が登録と違う"
[ "$(sha256sum $F32 | cut -c1-16)" = "793215917d156680" ] || die "float のバイナリが §6.26 と違う"
[ "$(sha256sum $F64 | cut -c1-16)" = "7e95e5bbd8f50c9e" ] || die "FP64 のバイナリが §6.26 と違う"
[ "$(df --output=avail -B1 . | tail -1 | tr -dc 0-9)" -ge 42949672960 ] || die "ディスクの空きが 40 GiB 未満"
# 4. メモリ (共有の B。4 本で約 9 GB。最大 1 時間待つ)
for i in $(seq 1 30); do av=$(free -m | awk '/^Mem:/{print $7}'); [ "$av" -ge 11000 ] && break; sleep 120; done
echo "available ${av} MB" >> $LOG; [ "$av" -ge 11000 ] || die "メモリが 11 GB 未満のまま"
arm() {   # arm <run> <バイナリ>
  local r=$1 b=$2
  mkdir $r || { echo "$r: 既にある (中止)" >> $LOG; return 1; }
  cp _g0c/*.yaml $r/ && cp _g0c/q32_init.h5 $r/nozzle.h5
  sed -i 's/nStepOuter: 200000}/nStepOuter: 30000}/; s/outStepInterval: 2500/outStepInterval: 500/; s/^output: .*/output: {level: 1, extraFields: [res_ro, res_roUx, res_roUy, res_roe, res_roK, res_roOmega, omg_prod, omg_dest, omg_cross, omg_trans, omg_axisym, dt_local]}/' $r/solverConfig.yaml
  cmp -s $r/solverConfig.yaml run_0717_c50_d1/solverConfig.yaml || { echo "$r: 設定が run_0717_c50_d1 (§6.51 の上限 50) と同一でない" >> $LOG; return 1; }
  ( export FORGE_BIN=$b; bash $RC $C45/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC ) &
  echo "$r 起動 $(date -Is) $(basename $(dirname $(dirname $(dirname $b))))" >> $LOG
}
RUNS=(run_0725_g0c_d1 run_0727_g0c_f1 run_0726_g0c_d2 run_0728_g0c_f2)
arm run_0725_g0c_d1 $F64 && sleep 120 && arm run_0727_g0c_f1 $F32 && sleep 120 && arm run_0726_g0c_d2 $F64 && sleep 120 && arm run_0728_g0c_f2 $F32 || die "起動の失敗"
while true; do
  left=0; for r in "${RUNS[@]}"; do [ -f $r/RUN_RC ] || left=1; done
  for r in "${RUNS[@]}"; do for f in $r/res_*.h5; do s=$(basename $f .h5); s=${s#res_}; [[ $s =~ ^[0-9]+$ ]] || continue
    if [ $s -gt 0 ] && [ $s -lt 20000 ]; then rm -f $f $r/res_$s.xmf; fi; done; done
  [ $left = 0 ] && break; sleep 120
done
for r in "${RUNS[@]}"; do echo "$r 終了 $(date -Is) RUN_RC=$(cat $r/RUN_RC) 最後の行: $(grep -h '^step' $r/forge_run.log | tail -1 | cut -c1-100)" >> $LOG; done
python3 rg0c_an.py >> $LOG 2>&1
echo "== 終了 $(date -Is)" >> $LOG
touch rg0c.done
