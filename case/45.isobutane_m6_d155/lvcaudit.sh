#!/bin/bash
# plan time_integration-line-viscous-jacobian-faceh §6.7 (2026-10-10): 製品の経路の照合の採取。
# run_0183 の res_100000 から値 3・マスク 7・キー 5・方向別・上限なし・cfl 4・ISP 0 を 1 step、ライン行列の書き出し (5 本のライン) つきで
#   run_0550_lvcaudit       監査用のビルド (-DFORGE_LINE_AUDIT、faceh_audit_fp64)
#   run_0551_lvcaudit_ctrl  同じソース・同じコンパイル条件の通常のビルド (faceh_ctrl_fp64)
#   run_0552_lvcaudit_ctrl2 通常のビルドの再実行 (atomicAdd による再実行の揺れを測る)
# AWS の自分の作業ツリー (~/forge-faceh-audit/case/45.isobutane_m6_d155) で動かす。run_0183 は ~/forge-wallfit の run へのリンク (読むだけ)。判定は lvcaudit_judge.py。
set -uo pipefail
cd "$(dirname "$0")"
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
export REAL_CONVERTER=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$PWD/conv_tolerant.sh
SRC=run_0183_ns_coldmesh_tw300_ext
LOG=$PWD/lvcaudit.log
NODES=1572,4113,7985,198560,264263
AUDIT_BIN=$HOME/forge-faceh-audit/solver_density_cuda/build/forge
CTRL_BIN=$HOME/forge-faceh-ctrl/solver_density_cuda/build/forge

echo "== 開始 $(date -Is)" >> $LOG
[ -e $SRC ] || ln -s $HOME/forge-wallfit/case/45.isobutane_m6_d155/$SRC $SRC
last=$(ls $SRC/ | grep -E '^res_[0-9]+\.h5$' | sort -t_ -k2 -n | tail -1)
[ "$last" = res_100000.h5 ] && [ "$(sha256sum $SRC/res_100000.h5 | cut -c1-16)" = 207d39f0e7f4aa03 ] || { echo "出発の場が res_100000.h5 (207d39f0…) でない ($last) — 止める" >> $LOG; exit 1; }

cap() {  # cap <run> <bin> <COLD_ALT_BINARY のキー> <監査か 0/1>
  local r=$1 bin=$2 key=$3 audit=$4
  ( export FORGE_BIN=$bin COLD_ALT_BINARY=$key
    python3 cold_cfl.py prep $SRC $r --steps 1 --out 1 --cfl 4 --limiter-ref-from $SRC --line dir --extra res_ro,volume --itj 5 --lvc 3 >> $LOG 2>&1 ) \
    || { echo "prep $r 失敗 — 止める" >> $LOG; return 1; }
  rm -f $r/nozzle.msh; mkdir -p $r/linedump
  ( export FORGE_BIN=$bin COLD_ALT_BINARY=$key FORGE_LINE_DUMP_DIR=$PWD/$r/linedump FORGE_LINE_DUMP_CALL=1 FORGE_LINE_DUMP_NODES=$NODES
    python3 cold_cfl.py run $r > $r/cold_pair_run_stdout.log 2>&1 ); local rc=$?
  echo "$r rc=$rc $(date -Is) $(grep -c '\[lineDump\] factor の直前を書いた' $r/forge_run.log) $(grep -c '\[lineAudit\] 記録を書いた' $r/forge_run.log)" >> $LOG
  [ $rc = 0 ] || { echo "$r: 終了コード $rc — 止める" >> $LOG; return 1; }
  grep -q '\[lineDump\] factor の直前を書いた' $r/forge_run.log || { echo "$r: 書き出しが無い — 止める" >> $LOG; return 1; }
  if [ $audit = 1 ]; then grep -q '\[lineAudit\] 記録を書いた' $r/forge_run.log || { echo "$r: 監査の記録が無い — 止める" >> $LOG; return 1; }
  else ! grep -q '\[lineAudit\]' $r/forge_run.log || { echo "$r: 通常のビルドなのに監査の記録がある — 止める" >> $LOG; return 1; }; fi
}

cap run_0550_lvcaudit $AUDIT_BIN faceh_audit_fp64 1 || exit 1
cap run_0551_lvcaudit_ctrl $CTRL_BIN faceh_ctrl_fp64 0 || exit 1
cap run_0552_lvcaudit_ctrl2 $CTRL_BIN faceh_ctrl_fp64 0 || exit 1
echo "== 終了 $(date -Is)" >> $LOG
touch lvcaudit.done
