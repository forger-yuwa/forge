#!/bin/bash
# plan architecture-float-state-double-geometry §6.53 (事前登録、codex diagnose 2026-10-12): G0c・B0 (上限 50) で初期場だけを変える A/B。インスタンス B。
# 1. 継続場の入力を作る: run_0725_g0c_d1/res_30000.h5 (G0c の FP64、30,000 step) の 9 保存量を restart_field.py (同一メッシュ) で
#    _g0c/q32_init.h5 の複製に写し、float32 に丸めて double に戻す → _ic/cont_q32.h5。
# 2. 4 本ずつ 2 回 (各回に旧起点・継続場 × FP64・float を 1 本ずつ、同じ起動順、30 秒ずつずらす)、各 6,000 step・250 ごと。
#    設定は run_0717_c50_d1 から nStepOuter (6000)・outStepInterval (250) だけを変える。判定は ric_an.py。
set -uo pipefail
TOKEN=$(curl -s -X PUT http://169.254.169.254/latest/api/token -H "X-aws-ec2-metadata-token-ttl-seconds: 60")
[ "$(curl -s -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-id)" = "i-0ba2b91ba659254c4" ] || { echo "B ではない。中止"; exit 1; }
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
export FORGE_OMEGA_BUDGET=1 FORGE_DIAG_COMMIT_LOSS=500
C45=~/forge-wallfit/case/45.isobutane_m6_d155; LOG=$C45/ric.log; TL=~/forge-wallfit/solver_density_cuda/tools
RC=$TL/run_case.sh
B=solver_density_cuda/build/forge; F32=~/forge-fgeom7-f32/$B; F64=~/forge-fgeom7-fp64/$B
cd $C45
die() { echo "中止: $* ($(date -Is))" >> $LOG; exit 1; }
echo "== 準備の開始 $(date -Is)" >> $LOG
[ "$(sha256sum _g0c/q32_init.h5 | cut -c1-16)" = "6f8fcda6d97482e3" ] || die "旧起点の入力が登録と違う"
[ -e run_0725_g0c_d1/res_30000.h5 ] || die "継続元 run_0725_g0c_d1/res_30000.h5 が無い"
[ -e _ic ] && die "_ic が既にある"
mkdir _ic && cp _g0c/*.yaml _ic/ && cp _g0c/q32_init.h5 _ic/cont_q32.h5 || die "_ic の準備"
python3 $TL/restart_field.py run_0725_g0c_d1/res_30000.h5 _ic/cont_q32.h5 --keep-src-dtype --dst-run _ic --forge $F64 > _ic/restart.log 2>&1 || die "restart_field.py (ログ _ic/restart.log)"
tail -3 _ic/restart.log >> $LOG
python3 - _ic/cont_q32.h5 run_0725_g0c_d1/res_30000.h5 >> $LOG 2>&1 <<'PY' || die "継続場の Q32 の丸め"
import sys, h5py, numpy as np
Q = ("ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega", "roY0", "roY1")
with h5py.File(sys.argv[1], "r+") as h, h5py.File(sys.argv[2], "r") as s:
    V = h["VALUE"]; n = 0
    for k in Q:
        a = np.asarray(V[k][:])
        if a.dtype != np.float64: raise SystemExit(f"{k} が float64 でない ({a.dtype})")
        if not np.array_equal(a, np.asarray(s["VALUE/" + k][:], dtype=np.float64)): raise SystemExit(f"{k} が継続元とビット一致しない")
        if not np.all(np.isfinite(a)): raise SystemExit(f"{k} に非有限")
        b = a.astype(np.float32).astype(np.float64); n += int(np.count_nonzero(a != b)); V[k][...] = b
    print(f"[Q32] 継続場の保存量 9 つを float32 に丸めた (変わった値 {n} 個)")
PY
SO=$(sha256sum _g0c/q32_init.h5 | cut -c1-16); SC=$(sha256sum _ic/cont_q32.h5 | cut -c1-16)
echo "== 開始 $(date -Is) F32 $(sha256sum $F32 | cut -c1-16) F64 $(sha256sum $F64 | cut -c1-16) 旧起点 $SO 継続場 $SC" >> $LOG
[ "$(sha256sum $F32 | cut -c1-16)" = "793215917d156680" ] || die "float のバイナリが §6.26 と違う"
[ "$(sha256sum $F64 | cut -c1-16)" = "7e95e5bbd8f50c9e" ] || die "FP64 のバイナリが §6.26 と違う"
[ "$(df --output=avail -B1 . | tail -1 | tr -dc 0-9)" -ge 42949672960 ] || die "ディスクの空きが 40 GiB 未満"
arm() {   # arm <run> <バイナリ> <入力>
  local r=$1 b=$2 inp=$3
  mkdir $r || { echo "$r: 既にある (中止)" >> $LOG; return 1; }
  cp _g0c/*.yaml $r/ && cp $inp $r/nozzle.h5
  sed -i 's/nStepOuter: 200000}/nStepOuter: 6000}/; s/outStepInterval: 2500/outStepInterval: 250/; s/^output: .*/output: {level: 1, extraFields: [res_ro, res_roUx, res_roUy, res_roe, res_roK, res_roOmega, omg_prod, omg_dest, omg_cross, omg_trans, omg_axisym, dt_local]}/' $r/solverConfig.yaml
  sed 's/nStepOuter: 30000}/nStepOuter: 6000}/; s/outStepInterval: 500/outStepInterval: 250/' run_0717_c50_d1/solverConfig.yaml | cmp -s - $r/solverConfig.yaml || { echo "$r: 設定が登録 (run_0717 から 2 行だけ変更) と同一でない" >> $LOG; return 1; }
  ( export FORGE_BIN=$b; bash $RC $C45/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC ) &
  echo "$r 起動 $(date -Is) $(basename $(dirname $(dirname $(dirname $b)))) 入力 $(sha256sum $r/nozzle.h5 | cut -c1-16)" >> $LOG
}
batch() {   # batch <old_d> <cont_d> <old_f> <cont_f>
  local RUNS=("$@")
  for i in $(seq 1 30); do av=$(free -m | awk '/^Mem:/{print $7}'); [ "$av" -ge 11000 ] && break; sleep 60; done
  echo "available ${av} MB" >> $LOG; [ "$av" -ge 11000 ] || return 1
  arm $1 $F64 _g0c/q32_init.h5 && sleep 30 && arm $2 $F64 _ic/cont_q32.h5 && sleep 30 && arm $3 $F32 _g0c/q32_init.h5 && sleep 30 && arm $4 $F32 _ic/cont_q32.h5 || return 1
  while true; do left=0; for r in "${RUNS[@]}"; do [ -f $r/RUN_RC ] || left=1; done; [ $left = 0 ] && break; sleep 30; done
  for r in "${RUNS[@]}"; do echo "$r 終了 $(date -Is) RUN_RC=$(cat $r/RUN_RC) 最後の行: $(grep -h '^step' $r/forge_run.log | tail -1 | cut -c1-100)" >> $LOG; done
}
batch run_0729_ic_old_d1 run_0730_ic_cont_d1 run_0731_ic_old_f1 run_0732_ic_cont_f1 || die "1 回目"
batch run_0733_ic_old_d2 run_0734_ic_cont_d2 run_0735_ic_old_f2 run_0736_ic_cont_f2 || die "2 回目"
python3 ric_an.py >> $LOG 2>&1
echo "== 終了 $(date -Is)" >> $LOG
touch ric.done
