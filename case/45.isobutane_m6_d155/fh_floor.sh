#!/bin/bash
# plan time_integration-implicit-thermal-jacobian §6.2 (2026-10-10): point 仕上げでの面エンタルピーの精度の A/B (§5.1 #5)。
# 出発 S = run_0354_m9_L5cut/res_40000.h5 (point 仕上げの最終の場)。設定は run_0354 と同じ (point・cfl 4・キー 0・リミッタの基準値は run_0183)。
# A = 既定 (面エンタルピーは float)、B = FORGE_DIAG_FACE_H_DOUBLE=1。
#   軌道 6 本 (A1 B1 A2 B2 A3 B3 の順、各 2000 step・500 ごと)。
#   共通の評価 28 本: 保存した状態から同じ設定で 1 step。step 0 の outer_begin の行が更新前の残差 R(Q)。
#     評価器 f = 切替なし、d = 切替あり。11 状態 × 2 + 再評価のノイズ 3 状態 × 2。
# バイナリは lineM_fp64 (sha256 05ad8bdf…) に固定。判定は fh_floor_judge.py (この台本は判定しない)。
# 場 (h5) は判定の後に消す (この台本では消さない。codex plan 段 2026-10-10 M1)。
set -uo pipefail
cd "$(dirname "$0")"
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
export REAL_CONVERTER=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$PWD/conv_tolerant.sh
export FORGE_BIN=$HOME/forge-linespeed-fp64/solver_density_cuda/build/forge COLD_ALT_BINARY=lineM_fp64
TOOLS=$(cd ../../solver_density_cuda/tools && pwd)
SRC=run_0183_ns_coldmesh_tw300_ext
S0=run_0354_m9_L5cut
ST=_fh_states
LOG=$PWD/fh_floor.log
B="--cfl 4 --limiter-ref-from $SRC --extra res_ro"

echo "== 開始 $(date -Is)" >> $LOG
if pgrep -x forge > /dev/null; then echo "他の forge が走っている — 止める ($(pgrep -x forge | tr '\n' ' '))" >> $LOG; exit 1; fi
[ "$(sha256sum $FORGE_BIN | cut -c1-16)" = 05ad8bdf6100874f ] || { echo "FORGE_BIN の sha256 が lineM_fp64 でない — 止める" >> $LOG; exit 1; }

# 格子の同一性: nozzle.h5 の MESH の全データセットのハッシュ (場は含めない)
meshsha() {  # meshsha <nozzle.h5> → stdout
  python3 -c '
import h5py, hashlib, sys
h = hashlib.sha256()
def visit(name, obj):
    if isinstance(obj, h5py.Dataset):
        h.update(name.encode()); h.update(obj[...].tobytes())
with h5py.File(sys.argv[1], "r") as f:
    f["MESH"].visititems(visit)
print(h.hexdigest())' "$1"
}

# 切替の表示を確かめる (B と評価器 d は表示がある、A と評価器 f は無い)
switch_ok() {  # switch_ok <run> <fh>
  if [ "$2" = 1 ]; then grep -q "FORGE_DIAG_FACE_H_DOUBLE" $1/forge_run.log
  else ! grep -q "FORGE_DIAG_FACE_H_DOUBLE" $1/forge_run.log; fi
}

state() {  # state <名前> <run> <res>: 状態だけを置いたディレクトリ (prep の --field-from は最後の res を取る)
  mkdir -p $ST/$1
  [ -f $2/$3 ] || { echo "状態 $2/$3 が無い" >> $LOG; return 1; }
  ln -sfn ../../$2/$3 $ST/$1/$3
  sha256sum $2/$3 | cut -d' ' -f1 > $ST/$1/SHA256
}

prep_run() {  # prep_run <run> <steps> <out> <状態の名前>
  python3 cold_cfl.py prep $SRC $1 --steps $2 --out $3 $B --field-from $ST/$4 >> $LOG 2>&1 || { echo "prep $1 失敗" >> $LOG; return 1; }
  rm -f $1/nozzle.msh
  meshsha $1/nozzle.h5 > $1/MESH_SHA.txt || { echo "$1: 格子のハッシュに失敗" >> $LOG; return 1; }
}

traj() {  # traj <run> <fh>
  local r=$1 fh=$2
  prep_run $r 2000 500 s0 || return 1
  ( [ "$fh" = 1 ] && export FORGE_DIAG_FACE_H_DOUBLE=1; python3 cold_cfl.py run $r > $r/cold_pair_run_stdout.log 2>&1 ); local rc=$?
  echo "$r fh=$fh rc=$rc $(date -Is) last=$(awk -F, '$3=="outer_begin"{s=$1} END{print s}' $r/residual_history.csv)" >> $LOG
  [ $rc = 0 ] || { echo "$r: 終了コード $rc — 止める" >> $LOG; return 1; }
  switch_ok $r $fh || { echo "$r: 切替の表示が腕と合わない — 止める" >> $LOG; return 1; }
  [ -f $r/res_2000.h5 ] || { echo "$r: res_2000.h5 が無い — 止める" >> $LOG; return 1; }
  ( cd $r && python3 $TOOLS/check_convergence.py . --segment > CONVERGENCE_SEGMENT.txt 2>&1 ); echo "$r check_convergence --segment rc=$?" >> $LOG
}

ev() {  # ev <run> <状態の名前> <fh>
  local r=$1 sd=$2 fh=$3
  prep_run $r 1 1 $sd || return 1
  ( [ "$fh" = 1 ] && export FORGE_DIAG_FACE_H_DOUBLE=1; python3 cold_cfl.py run $r > $r/cold_pair_run_stdout.log 2>&1 ); local rc=$?
  echo "$r state=$sd fh=$fh rc=$rc $(date -Is) step0=$(awk -F, '$1==0&&$3=="outer_begin"{print $8}' $r/residual_history.csv)" >> $LOG
  [ $rc = 0 ] || { echo "$r: 終了コード $rc — 止める" >> $LOG; return 1; }
  switch_ok $r $fh || { echo "$r: 切替の表示が評価器と合わない — 止める" >> $LOG; return 1; }
  grep -q "^0,-1,outer_begin" $r/residual_history.csv || { echo "$r: step 0 の行が無い — 止める" >> $LOG; return 1; }
}

mkdir -p $ST && meshsha $S0/nozzle.h5 > $ST/MESH_SHA_REF.txt 2>> $LOG || { echo "run_0354 の格子のハッシュに失敗" >> $LOG; exit 1; }
state s0 $S0 res_40000.h5 || exit 1
i=0
for arm in a1:0 b1:1 a2:0 b2:1 a3:0 b3:1; do
  n=${arm%%:*}; fh=${arm##*:}
  traj run_$(printf %04d $((500 + i)))_fh_$n $fh || exit 1
  i=$((i + 1))
done
state a1k1000 run_0500_fh_a1 res_1000.h5 && state a1k1500 run_0500_fh_a1 res_1500.h5 && state a1k2000 run_0500_fh_a1 res_2000.h5 &&
state b1k1000 run_0501_fh_b1 res_1000.h5 && state b1k1500 run_0501_fh_b1 res_1500.h5 && state b1k2000 run_0501_fh_b1 res_2000.h5 &&
state a2k2000 run_0502_fh_a2 res_2000.h5 && state b2k2000 run_0503_fh_b2 res_2000.h5 &&
state a3k2000 run_0504_fh_a3 res_2000.h5 && state b3k2000 run_0505_fh_b3 res_2000.h5 || exit 1

k=510
for sd in s0 a1k1000 a1k1500 a1k2000 b1k1000 b1k1500 b1k2000 a2k2000 a3k2000 b2k2000 b3k2000; do
  ev run_$(printf %04d $k)_fhe_${sd}_f $sd 0 || exit 1; k=$((k + 1))
  ev run_$(printf %04d $k)_fhe_${sd}_d $sd 1 || exit 1; k=$((k + 1))
done
for sd in s0 a1k2000 b1k2000; do   # 再評価のノイズ (同じ状態・同じ評価器をもう 1 回)
  ev run_$(printf %04d $k)_fhe_${sd}_fr $sd 0 || exit 1; k=$((k + 1))
  ev run_$(printf %04d $k)_fhe_${sd}_dr $sd 1 || exit 1; k=$((k + 1))
done
echo "== 終了 $(date -Is)" >> $LOG
touch fh_floor.done
