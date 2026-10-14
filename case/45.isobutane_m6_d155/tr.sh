#!/bin/bash
# plan architecture-float-state-double-geometry §6.26 (事前登録、codex diagnose 2026-10-11): 共通の Q32 の起点からの float と FP64 の軌跡の A/B。
# 起点 run_0483_hp7_k1 の res_145000 を新格子の nozzle.h5 に移し (9 量、ビット一致)、保存量 9 つだけを float32 に丸めて double に戻した共通の入力を作る。
# A = FP64、B = float、B' = float の再実行。各 30,000 step 固定・500 step ごと。インスタンス B で 3 本同時。判定は tr_an.py。
# 入力 (A から送っておく): ~/tr_in/res_145000.h5、~/kab_in/ の nozzle.h5・*.yaml (kab.sh と同じ)。
set -uo pipefail
TOKEN=$(curl -s -X PUT http://169.254.169.254/latest/api/token -H "X-aws-ec2-metadata-token-ttl-seconds: 60")
[ "$(curl -s -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-id)" = "i-0ba2b91ba659254c4" ] || { echo "B ではない。中止"; exit 1; }
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
export FORGE_OMEGA_BUDGET=1 FORGE_DIAG_COMMIT_LOSS=500
C45=~/forge-wallfit/case/45.isobutane_m6_d155; LOG=$C45/tr.log
RC=~/forge-wallfit/solver_density_cuda/tools/run_case.sh; TL=~/forge-wallfit/solver_density_cuda/tools
B=solver_density_cuda/build/forge; F64=~/forge-fgeom7-fp64/$B; F32=~/forge-fgeom7-f32/$B
SRC=~/tr_in/res_145000.h5
cd $C45
echo "== 開始 $(date -Is) F64 $(sha256sum $F64 | cut -c1-16) F32 $(sha256sum $F32 | cut -c1-16) 起点 $(sha256sum $SRC | cut -c1-16)" >> $LOG
[ "$(df --output=avail -B1 . | tail -1 | tr -dc 0-9)" -ge 32212254720 ] || { echo "ディスクの空きが 30 GiB 未満 (中止)" >> $LOG; exit 1; }
# 共通の入力 (Q32)
mkdir -p _tr && cp ~/kab_in/*.yaml _tr/ && cp ~/kab_in/nozzle.h5 _tr/q32_init.h5
python3 $TL/restart_field.py $SRC _tr/q32_init.h5 --dst-run _tr --forge $F64 > _tr/restart.log 2>&1 || { echo "restart_field に失敗 (中止)" >> $LOG; exit 1; }
grep -q "VERDICT: OK (9 量を移した、SRC とビット一致)" _tr/restart.log || { echo "移した量が想定と違う: $(grep VERDICT _tr/restart.log) (中止)" >> $LOG; exit 1; }
python3 - _tr/q32_init.h5 >> $LOG 2>&1 <<'EOF' || { echo "Q32 の丸めに失敗 (中止)" >> $LOG; exit 1; }
import sys, h5py, numpy as np
Q = ("ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega", "roY0", "roY1")
with h5py.File(sys.argv[1], "r+") as h:
    V = h["VALUE"]; n = 0
    for k in Q:
        a = np.asarray(V[k][:])
        if a.dtype != np.float64: raise SystemExit(f"{k} が float64 でない ({a.dtype})")
        b = a.astype(np.float32).astype(np.float64); n += int(np.count_nonzero(a != b)); V[k][...] = b
    print(f"[Q32] 保存量 {len(Q)} つを float32 に丸めた (変わった値 {n} 個)。wall_dist などほかは変えない")
EOF
echo "共通の入力 _tr/q32_init.h5 sha $(sha256sum _tr/q32_init.h5 | cut -c1-16)" >> $LOG
arm() {   # arm <run> <バイナリ>
  local r=$1 b=$2
  mkdir $r || { echo "$r: 既にある (中止)" >> $LOG; return 1; }
  cp _tr/*.yaml $r/ && cp _tr/q32_init.h5 $r/nozzle.h5
  sed -i 's/nStepOuter: 200000}/nStepOuter: 30000}/; s/outStepInterval: 2500/outStepInterval: 500/; s/^output: .*/output: {level: 1, extraFields: [res_ro, res_roUx, res_roUy, res_roe, res_roK, res_roOmega, omg_prod, omg_dest, omg_cross, omg_trans, omg_axisym]}/' $r/solverConfig.yaml
  grep -q "nStepOuter: 30000}" $r/solverConfig.yaml && grep -q "outStepInterval: 500$" $r/solverConfig.yaml && grep -q "^mesh: {axisSegmentRWeight: 1," $r/solverConfig.yaml && grep -q "omg_axisym" $r/solverConfig.yaml || { echo "$r: 設定の書き換えに失敗" >> $LOG; return 1; }
  ( export FORGE_BIN=$b; bash $RC $C45/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC ) &
  echo "$r 起動 $(date -Is) $(basename $(dirname $(dirname $(dirname $b))))" >> $LOG
}
arm run_0487_tr_fp64 $F64 || exit 1
arm run_0488_tr_f32 $F32 || exit 1
arm run_0489_tr_f32b $F32 || exit 1
sleep 240
for r in run_0487_tr_fp64 run_0488_tr_f32 run_0489_tr_f32b; do echo "$r 起動の確認: $(grep -h 'axisSegmentRWeight' $r/forge_run.log | head -1 | cut -c1-90) / $(grep -h "extraFields" $r/forge_run.log | head -1 | cut -c1-80)" >> $LOG; done
wait
for r in run_0487_tr_fp64 run_0488_tr_f32 run_0489_tr_f32b; do echo "$r 終了 $(date -Is) RUN_RC=$(cat $r/RUN_RC) 最後の行: $(grep -h '^step' $r/forge_run.log | tail -1 | cut -c1-100)" >> $LOG; done
python3 tr_an.py >> $LOG 2>&1
echo "== 終了 $(date -Is)" >> $LOG
touch tr.done
