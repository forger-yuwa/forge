#!/bin/bash
# plan time_integration-line-implicit-speed §6.19 の「既定化の後の確認」(2026-10-10): 新しいバイナリ (LAYOUT2 を既定にした版) の既定の経路・従来へ戻る分岐。
# FP64 = ~/forge-linespeed-fp64 (f9be0c4f + 既定化 + typedef double)、float = ~/forge-linespeed-f32 (f9be0c4f + 既定化)。
set -uo pipefail
cd "$(dirname "$0")"
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
LOG=$PWD/lay4.log
B64=$HOME/forge-linespeed-fp64/solver_density_cuda/build/forge
B32=$HOME/forge-linespeed-f32/solver_density_cuda/build/forge
RC=$HOME/forge-wallfit/solver_density_cuda/tools/run_case.sh
SRC=run_0183_ns_coldmesh_tw300_ext
G56=../56.gap_tp1187
FAILS=0
chk() { if eval "$2"; then echo "OK   $1" >> $LOG; else echo "NG   $1" >> $LOG; FAILS=$((FAILS + 1)); fi; }
echo "== 開始 $(date -Is)" >> $LOG
echo "FP64 sha256 $(sha256sum $B64 | cut -d' ' -f1)、float sha256 $(sha256sum $B32 | cut -d' ' -f1)" >> $LOG
prep45() { ( export FORGE_BIN=$(cat prof_target_linevisc.txt) COLD_ALT_BINARY=lineL_fp64 REAL_CONVERTER=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$PWD/conv_tolerant.sh
             python3 cold_cfl.py prep $SRC $1 --steps $2 --out $2 --cfl 4 --limiter-ref-from $SRC --line dir --itj 5 --cap 50 > /dev/null ) && rm -f $1/nozzle.msh; }
prep56() { mkdir $G56/$1 && cp $G56/run_0073_lay2cmp_varlen/{mesh.h5,bcondConfig.yaml,probe.yaml,case_setup.json,solverConfig.yaml} $G56/$1/; }
run() { local r=$1; shift; ( for kv in "$@"; do export "$kv"; done; bash $RC $r > $r/run_case_stdout.log 2>&1 ); echo $?; }
# R1 case/45 FP64 既定の経路 (比較あり)
r=run_0381_lay2def_cmp; prep45 $r 20; t0=$(date +%s); rc=$(run $r FORGE_BIN=$B64 FORGE_LINE_COMPARE=1)
chk "R1 rc=0" "[ $rc = 0 ]"
chk "R1 既定で LAYOUT2" "grep -q 'Thomas の配列の並び: LAYOUT2 (既定' $r/forge_run.log"
python3 line_cmp_judge.py $r --arm LAYOUT2 --eta-limit 1e-11 --bitwise --factors 20 --solves 100 --lines 4719 > /dev/null 2>&1
chk "R1 ビット一致 PASS" "grep -q '\"verdict\": \"PASS\"' $r/cmp_judge_LAYOUT2.json && [ \$(stat -c %Y $r/cmp_judge_LAYOUT2.json) -ge $t0 ]"
# R2 case/45 FP64 FORGE_LINE_LAYOUT=0
r=run_0382_lay2def_off; prep45 $r 20; rc=$(run $r FORGE_BIN=$B64 FORGE_LINE_LAYOUT=0)
chk "R2 rc=0" "[ $rc = 0 ]"
chk "R2 従来の並び" "grep -q 'Thomas の配列の並び: 従来 (節点番号の並び) (FORGE_LINE_LAYOUT=0' $r/forge_run.log && grep -q '^\[line\] factor 1 回目: モード LU$' $r/forge_run.log"
# R3 診断のスイッチ (未指定) → 従来
r=$G56/run_0074_lay2def_diag; prep56 run_0074_lay2def_diag; rc=$(run $r FORGE_BIN=$B32 FORGE_LINE_INV=1)
chk "R3 rc=0" "[ $rc = 0 ]"
chk "R3 診断で従来・INV" "grep -q 'Thomas の配列の並び: 従来 (診断のスイッチ FORGE_LINE_INV' $r/forge_run.log && grep -q '^\[line\] factor 1 回目: モード INV$' $r/forge_run.log"
# R4 明示の 2 と診断のスイッチ → 止まる
r=$G56/run_0075_lay2def_refuse; prep56 run_0075_lay2def_refuse; rc=$(run $r FORGE_BIN=$B32 FORGE_LINE_LAYOUT=2 FORGE_LINE_INV=1)
chk "R4 止まる (rc≠0・メッセージ)" "[ $rc != 0 ] && grep -rq 'FORGE_LINE_LAYOUT=2 と FORGE_LINE_INV は組み合わせない' $r/"
# R5 区画比の上限を下げる → 従来
r=$G56/run_0076_lay2def_ratio; prep56 run_0076_lay2def_ratio; rc=$(run $r FORGE_BIN=$B32 FORGE_LINE_LAYOUT_SLOT_RATIO=1.0)
chk "R5 rc=0・区画比で従来" "[ $rc = 0 ] && grep -q 'Thomas の配列の並び: 従来 (区画 107 × 643' $r/forge_run.log && grep -q '^\[line\] factor 1 回目: モード LU$' $r/forge_run.log"
# R6 確保の失敗の模擬 → 従来
r=$G56/run_0077_lay2def_oom; prep56 run_0077_lay2def_oom; rc=$(run $r FORGE_BIN=$B32 FORGE_LINE_LAYOUT_FAKE_OOM=1)
chk "R6 rc=0・確保の失敗で従来" "[ $rc = 0 ] && grep -q 'Thomas の配列の並び: 従来 (並べ替えた配列を確保できない)' $r/forge_run.log && grep -q '^\[line\] factor 1 回目: モード LU$' $r/forge_run.log"
# R7 float 既定の経路 (比較あり)、長さがばらつくライン
r=$G56/run_0078_lay2def_cmp; prep56 run_0078_lay2def_cmp; t0=$(date +%s); rc=$(run $r FORGE_BIN=$B32 FORGE_LINE_COMPARE=1)
chk "R7 rc=0・既定で LAYOUT2" "[ $rc = 0 ] && grep -q 'Thomas の配列の並び: LAYOUT2 (既定' $r/forge_run.log"
python3 line_cmp_judge.py $r --arm LAYOUT2 --eta-limit 1e-11 --bitwise --factors 20 --solves 80 --lines 643 > /dev/null 2>&1
LAY3_JUDGE_OUT=$PWD/_band_ab/cold_pair/lay4_R7_judge.json python3 lay3_judge.py $r $t0 > $r/lay3_judge_out.txt 2>&1; c=$?
chk "R7 ビット一致 (η の例外の規則)" "[ $c = 0 ]"
rm -f run_0381_lay2def_cmp/{nozzle.h5,res_0.h5,res_20.h5} run_0382_lay2def_off/{nozzle.h5,res_0.h5,res_20.h5}
echo "NG の数 $FAILS" >> $LOG
echo "== 終了 $(date -Is)" >> $LOG
echo $FAILS > lay4.rc
touch lay4.done
