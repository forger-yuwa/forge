#!/bin/bash
# plan tooling-nozzle-core-grid §4.76 (事前登録、codex 諮問 2026-10-13 after-small-amp、ユーザ決定「いいよ」): 実行間差と摂動への応答の分離。
# Grefc の run_0625 の res_60000 から、基準 R・ビット一致の複製 Z・摂動 P (a = 1e−7、§4.65 と同じ形) の 3 本を同じバイナリ・設定 (上限 50) で同時に各 3000 step・250 ごと。
# 出力ごとに (R, Z) と (R, P) の差を cg_twin_diff.py で取り (タグ rzp_noise・rzp_sig)、0・1000 の倍数・3000 以外の場を消す。
# usage (AWS の ~/forge-coregrid/case/45.isobutane_m6_d155): bash cg_twin_rzp.sh
# **実行中にこのファイルを編集しないこと**。
set -uo pipefail
cd "$(dirname "$0")"
C45=$PWD
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
TL=$C45/../../solver_density_cuda/tools; RC=$TL/run_case.sh
F7_64=$HOME/forge-fgeom7-fp64/solver_density_cuda/build/forge
SRC=run_0625_cg_grefc; SEED=$SRC/res_60000.h5; REFC=run_0619_cg_gref_b0_ctl/res_60000.h5; LOG=$C45/cg_twin_rzp.log; AMP=1e-7; STEPS=3000; OUT=250
EXTRA="[res_ro, res_roUx, res_roUy, res_roe, res_roK, res_roOmega, volume, dt_local]"
free_gb() { df --output=avail -B1G . | tail -1 | tr -dc 0-9; }
echo "== 開始 $(date -Is) F7_64 $(sha256sum $F7_64 | cut -c1-16) 起点 $(sha256sum $SEED | cut -c1-16)" >> $LOG
[ "$(free_gb)" -ge 8 ] || { echo "ディスクの空きが 8 GB 未満 (中止)" >> $LOG; exit 1; }
mk() {  # mk <run> <摂動 0/1>
  local r=$1 pert=$2
  [ -e $r ] && { echo "$r: 既にある (中止)" >> $LOG; return 1; }
  mkdir $r; for f in $SRC/*; do case $(basename $f) in res_*|*.log|RUN_*|forge_*|residual_*|CONVERGENCE_VERDICT.txt|.done_*|*.png) ;; *) cp -a $f $r/;; esac; done
  python3 $TL/restart_field.py $SEED $r/nozzle.h5 --dst-run $r --forge $F7_64 > $r/restart.log 2>&1 || { echo "$r: restart_field に失敗" >> $LOG; return 1; }
  grep -q "VERDICT: OK (9 量を移した、SRC とビット一致)" $r/restart.log || { echo "$r: 移した量が想定と違う" >> $LOG; return 1; }
  sed -i "s/nStepOuter: [0-9]*}/nStepOuter: $STEPS}/; s/outStepInterval: [0-9]*/outStepInterval: $OUT/; s/output: {level: 1, extraFields: \[[^]]*\]}/output: {level: 1, extraFields: $EXTRA}/" $r/solverConfig.yaml
  local d=$(diff $SRC/solverConfig.yaml $r/solverConfig.yaml | grep -c '^>')
  [ "$d" = 3 ] && grep -q "lineDtDirectionalCap: 50.0," $r/solverConfig.yaml && grep -q "nStepOuter: $STEPS}" $r/solverConfig.yaml || { echo "$r: 設定の差が想定と違う ($d 行) (中止)" >> $LOG; return 1; }
  if [ $pert = 1 ]; then python3 cg_twin_perturb.py $r $AMP $SEED --ref-coords $REFC >> $LOG 2>&1 || { echo "$r: 摂動に失敗" >> $LOG; return 1; }; fi
  echo "$r: Grefc・上限 50・$STEPS step・出力 $OUT ごと・摂動 $pert" >> $LOG
}
R=run_0683_rzp_r; Z=run_0684_rzp_z; P=run_0685_rzp_p
mk $R 0 && mk $Z 0 && mk $P 1 || exit 1
# Z の入力が R とビット一致、P と R の差は roe だけ
python3 - $R $Z $P >> $LOG 2>&1 <<'PY' || { echo "入力の一致の検査に失敗 (中止)" >> $LOG; exit 1; }
import sys, h5py, numpy as np
R, Z, P = sys.argv[1:4]
def vals(r):
    with h5py.File(f"{r}/nozzle.h5", "r") as h: return {k: np.asarray(h["VALUE/" + k][:]) for k in h["VALUE"].keys()}
vr, vz, vp = vals(R), vals(Z), vals(P)
assert vr.keys() == vz.keys() == vp.keys()
bad = [k for k in vr if not np.array_equal(vr[k], vz[k])]
if bad: raise SystemExit(f"Z が R と違う: {bad}")
diff = [k for k in vr if not np.array_equal(vr[k], vp[k])]
if diff != ["roe"]: raise SystemExit(f"P と R の差が roe だけでない: {diff}")
print("入力の検査: Z は R とビット一致、P と R の差は roe だけ")
PY
run() { ( export FORGE_BIN=$F7_64; bash $RC $C45/$1 > $1/run_case_stdout.log 2>&1; echo $? > $1/RUN_RC ); }
( run $R ) & PR=$!; sleep 10; ( run $Z ) & PZ=$!; sleep 10; ( run $P ) & PP=$!; echo "$R・$Z・$P 起動 $(date -Is)" >> $LOG
while true; do
  for f in $R/res_*.h5; do
    n=$(basename $f .h5); n=${n#res_}; [[ $n =~ ^[0-9]+$ ]] || continue
    [ -f $R/.rzp_$n ] && continue; [ -f $Z/res_$n.h5 ] && [ -f $P/res_$n.h5 ] || continue
    if { [ -f $R/res_$((n + OUT)).h5 ] && [ -f $Z/res_$((n + OUT)).h5 ] && [ -f $P/res_$((n + OUT)).h5 ]; } || { [ -f $R/RUN_RC ] && [ -f $Z/RUN_RC ] && [ -f $P/RUN_RC ]; }; then
      python3 cg_twin_diff.py $R $Z $n rzp_noise && python3 cg_twin_diff.py $R $P $n rzp_sig && touch $R/.rzp_$n
      if [ $n -gt 0 ] && [ $((n % 1000)) -ne 0 ] && [ $n -ne $STEPS ]; then rm -f $R/res_$n.* $Z/res_$n.* $P/res_$n.*; fi
    fi
  done
  if [ -f $R/RUN_RC ] && [ -f $Z/RUN_RC ] && [ -f $P/RUN_RC ]; then left=0; for f in $R/res_*.h5; do n=$(basename $f .h5); n=${n#res_}; [[ $n =~ ^[0-9]+$ ]] && [ ! -f $R/.rzp_$n ] && left=1; done; [ $left = 0 ] && break; fi
  [ "$(free_gb)" -lt 3 ] && { echo "ディスク 3 GB 未満 (中止)" >> $LOG; kill $PR $PZ $PP 2>/dev/null; exit 1; }
  sleep 20
done
wait $PR $PZ $PP
echo "$R RUN_RC=$(cat $R/RUN_RC)、$Z RUN_RC=$(cat $Z/RUN_RC)、$P RUN_RC=$(cat $P/RUN_RC) $(date -Is)" >> $LOG
echo "== 終了 $(date -Is)" >> $LOG
