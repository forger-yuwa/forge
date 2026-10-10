#!/bin/bash
# plan architecture-float-state-double-geometry §6.18 (事前登録): 精度の切り替えの A/B。起点は run_0451_v4_f32 の res_140000 (float の到達の場)。
#   A = float で続ける、A2 = A の再実行、B = FP64 に切り替える。各 30,000 step 固定、2,500 step ごと (res_ro 入り)。3 本を同時に投入。
# V4 と同じバイナリ (~/forge-fgeom5-*、56b1a6fe)。A には commit の診断 (FORGE_DIAG_COMMIT_LOSS=2500) を付ける。
set -uo pipefail
cd ~/forge-wallfit/case/45.isobutane_m6_d155
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
RC=~/forge-wallfit/solver_density_cuda/tools/run_case.sh
B32=~/forge-fgeom5-f32/solver_density_cuda/build/forge; B64=~/forge-fgeom5-fp64/solver_density_cuda/build/forge
export REAL_CONVERTER=~/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$PWD/conv_tolerant.sh
SRC=run_0183_ns_coldmesh_tw300_ext; FROM=run_0451_v4_f32
PREP64=~/forge-linespeed-fp64/solver_density_cuda/build/forge   # cold_cfl.py が登録と照合する準備用 (計算には使わない)
LOG=$PWD/v6ab.log
echo "== 開始 $(date -Is)" >> $LOG
for b in B32 B64; do echo "$b $(sha256sum ${!b} | cut -c1-16)" >> $LOG; done
[ "$(ls $FROM/res_*.h5 | grep -E 'res_[0-9]+\.h5$' | sed 's/.*res_//; s/\.h5//' | sort -n | tail -1)" = "140000" ] || { echo "$FROM の最終の res が 140000 でない (中止)" >> $LOG; exit 1; }
[ "$(df --output=avail -B1 . | tail -1 | tr -dc 0-9)" -ge 10737418240 ] || { echo "ディスクの空きが 10 GiB 未満 (中止)" >> $LOG; exit 1; }
arm() {   # arm <run> <バイナリ> [環境変数...]
  local r=$1 b=$2; shift 2
  ( export FORGE_BIN=$PREP64 COLD_ALT_BINARY=lineM_fp64; python3 cold_cfl.py prep $SRC $r --steps 30000 --out 2500 --cfl 4 --limiter-ref-from $SRC --line dir --itj 5 --cap 50 --extra res_ro --field-from $FROM > /dev/null ) || { echo "$r: 準備に失敗" >> $LOG; return 1; }
  rm -f $r/nozzle.msh
  grep -h "restart_field\|ビット\|bit" $r/restart_field.log 2>/dev/null | tail -2 | sed "s#^#  $r: #" >> $LOG
  ( export FORGE_BIN=$b; for kv in "$@"; do export "$kv"; done; bash $RC $PWD/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC ) &
  echo "$r 起動 $(date -Is) $(basename $(dirname $(dirname $(dirname $b)))) $*" >> $LOG
}
arm run_0460_ab_f32cont $B32 FORGE_DIAG_COMMIT_LOSS=2500
arm run_0461_ab_f32cont2 $B32
arm run_0462_ab_fp64switch $B64
wait
for r in run_0460_ab_f32cont run_0461_ab_f32cont2 run_0462_ab_fp64switch; do echo "$r 終了 $(date -Is) RUN_RC=$(cat $r/RUN_RC 2>/dev/null) 最後の res $(ls $r | grep -E '^res_[0-9]+\.h5$' | sed 's/res_//; s/\.h5//' | sort -n | tail -1)" >> $LOG; done
echo "== 終了 $(date -Is)" >> $LOG
touch v6ab.done
