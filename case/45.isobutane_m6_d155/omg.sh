#!/bin/bash
# ω の残差の場 (2026-10-10、ユーザ「ωがなんで高いんだろうか」): M64・B0・P の最終の状態から元の構成で 1 step、残差と Δτ を書き出す。
# 残差は状態だけで決まる (左辺に依らない)。バイナリ lineM_fp64 (LAYOUT2 既定)。
set -uo pipefail
cd "$(dirname "$0")"
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
export REAL_CONVERTER=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$PWD/conv_tolerant.sh
export FORGE_BIN=$HOME/forge-linespeed-fp64/solver_density_cuda/build/forge COLD_ALT_BINARY=lineM_fp64
LOG=$PWD/omg.log
SRC=run_0183_ns_coldmesh_tw300_ext
EX="--extra res_ro,res_roUx,res_roe,res_roK,res_roOmega,omega,k,vis_turb,dt_local"
B="--cfl 4 --limiter-ref-from $SRC"
echo "== 開始 $(date -Is)" >> $LOG
one() {   # one <run> <field-from> <構成の引数...> [env]
  local r=$1 ff=$2; shift 2
  python3 cold_cfl.py prep $SRC $r --steps 1 --out 1 $B $EX --field-from $ff "$@" > /dev/null || { echo "prep $r 失敗" >> $LOG; return; }
  rm -f $r/nozzle.msh
}
one run_0384_omg_M64 run_0378_tt_M64 --line dir --itj 5 --cap 50
( export FORGE_LINE_MAXLEN=64; python3 cold_cfl.py run run_0384_omg_M64 > run_0384_omg_M64/cold_pair_run_stdout.log 2>&1 ); echo "run_0384 rc=$?" >> $LOG
one run_0385_omg_B0 run_0252_ns_coldmesh_tw300_linedir_tj5_cap50_ext2 --line dir --itj 5 --cap 50
python3 cold_cfl.py run run_0385_omg_B0 > run_0385_omg_B0/cold_pair_run_stdout.log 2>&1; echo "run_0385 rc=$?" >> $LOG
one run_0386_omg_P run_0263_ns_coldmesh_tw300_cfl4_ext4
python3 cold_cfl.py run run_0386_omg_P > run_0386_omg_P/cold_pair_run_stdout.log 2>&1; echo "run_0386 rc=$?" >> $LOG
for r in run_0384_omg_M64 run_0385_omg_B0 run_0386_omg_P; do echo "$r $(grep -m1 'lineImplicit\] lines' $r/forge_run.log) $(grep -m1 'Thomas の配列の並び' $r/forge_run.log)" >> $LOG; head -1 $r/residual_history.csv | cut -d, -f1-12 >> $LOG; grep ",outer_end," $r/residual_history.csv | head -2 | cut -d, -f1-12 >> $LOG; done
echo "== 終了 $(date -Is)" >> $LOG
touch omg.done
