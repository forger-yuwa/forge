#!/bin/bash
# plan time_integration-line-viscous-jacobian-faceh §6 (2026-10-10): 熱伝導の近傍 K を入れたライン粘性 Jacobian (マスク 7) の破綻に、
# 面エンタルピーの float の評価が要るかの A/B。run_0183 の res_100000 から、キー 5・方向別・上限なし・cfl 4 (親 plan の run_0311 と同じ)。
#   主の 4 本だけ (値 3: a1 (既定)・b1 (FORGE_DIAG_FACE_H_DOUBLE=1)・a2・b2)。最大 2000 step・100 ごと (既知の破綻 122〜150 step の前の場を残す)。
#   序盤 200 step は E1b (run_0311) と同じ 112 節点の帳簿。値 2 の副の対は主の判定の後に決める (codex 諮問 2026-10-10)。
# 発散 (detectNaN で止まる、RUN_RC 1) は結果なので台本は止めずに次へ進む。起動の失敗・切替の表示の不一致では止める。
# バイナリは lineM_fp64 (sha256 05ad8bdf…) に固定。判定は lvcfh_judge.py (この台本は判定しない)。場・帳簿・res_nan は result 段のレビューまで残す。
set -uo pipefail
cd "$(dirname "$0")"
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
export REAL_CONVERTER=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$PWD/conv_tolerant.sh
export FORGE_BIN=$HOME/forge-linespeed-fp64/solver_density_cuda/build/forge COLD_ALT_BINARY=lineM_fp64
TOOLS=$(cd ../../solver_density_cuda/tools && pwd)
SRC=run_0183_ns_coldmesh_tw300_ext
LOG=$PWD/lvcfh.log

echo "== 開始 $(date -Is)" >> $LOG
if pgrep -x forge > /dev/null; then echo "他の forge が走っている — 止める ($(pgrep -x forge | tr '\n' ' '))" >> $LOG; exit 1; fi
[ "$(sha256sum $FORGE_BIN | cut -c1-16)" = 05ad8bdf6100874f ] || { echo "FORGE_BIN の sha256 が lineM_fp64 でない — 止める" >> $LOG; exit 1; }
# 出発の場: prep は $SRC の最後の res を取るので、それが res_100000.h5 で事前に固定した sha256 であることを先に確かめる
last=$(ls $SRC | grep -E '^res_[0-9]+\.h5$' | sort -t_ -k2 -n | tail -1)
[ "$last" = res_100000.h5 ] && [ "$(sha256sum $SRC/res_100000.h5 | cut -c1-16)" = 207d39f0e7f4aa03 ] || { echo "出発の場が res_100000.h5 (207d39f0…) でない ($last) — 止める" >> $LOG; exit 1; }
LEDGER_NODES=1572,1571,1570,1569,1568,1567,1566,1565,1564,1563,1562,1561,1560,1559,1558,1557,4113,4112,4111,4110,4109,4108,4107,4106,4105,4104,4103,4102,4101,4100,4099,4098,4960,4959,4958,4957,4956,4955,4954,4953,4952,4951,4950,4949,4948,4947,4946,4945,6170,6169,6168,6167,6166,6165,6164,6163,6162,6161,6160,6159,6158,6157,6156,6155,7985,7984,7983,7982,7981,7980,7979,7978,7977,7976,7975,7974,7973,7972,7971,7970,198560,198559,198558,198557,198556,198555,198554,198553,198552,198551,198550,198549,198548,198547,198546,198545,264263,264262,264261,264260,264259,264258,264257,264256,264255,264254,264253,264252,264251,264250,264249,264248

meshsha() {  # meshsha <nozzle.h5> → stdout (MESH の全データセットのハッシュ、場は含めない)
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

arm() {  # arm <run> <値 2|3> <fh 0|1>
  local r=$1 v=$2 fh=$3
  python3 cold_cfl.py prep $SRC $r --steps 2000 --out 100 --cfl 4 --limiter-ref-from $SRC --line dir --extra res_ro,volume --itj 5 --lvc $v >> $LOG 2>&1 \
    || { echo "prep $r 失敗 — 止める" >> $LOG; return 1; }
  rm -f $r/nozzle.msh
  meshsha $r/nozzle.h5 > $r/MESH_SHA.txt || { echo "$r: 格子のハッシュに失敗 — 止める" >> $LOG; return 1; }
  ( export FORGE_DUMP_LEDGER=$PWD/$r/ledger.csv FORGE_DUMP_LEDGER_NODES=$LEDGER_NODES FORGE_DUMP_LEDGER_CALLS=200
    [ "$fh" = 1 ] && export FORGE_DIAG_FACE_H_DOUBLE=1; python3 cold_cfl.py run $r > $r/cold_pair_run_stdout.log 2>&1 ); local rc=$?
  local nan=$(grep -m1 -o "Non-finite value detected in '[A-Za-z0-9_]*' at step [0-9]*" $r/forge_run.log)
  echo "$r v=$v fh=$fh rc=$rc $(date -Is) last=$(awk -F, '$3=="outer_begin"{s=$1} END{print s}' $r/residual_history.csv 2>/dev/null) ${nan:-有限}" >> $LOG
  [ -s $r/residual_history.csv ] && grep -q "'lineViscCoupling' in 'time.deltaT': $v" $r/forge_run.log \
    || { echo "$r: 起動していないか値の表示が違う — 止める" >> $LOG; return 1; }
  if [ "$fh" = 1 ]; then grep -q "FORGE_DIAG_FACE_H_DOUBLE" $r/forge_run.log || { echo "$r: 切替の表示が無い — 止める" >> $LOG; return 1; }
  else ! grep -q "FORGE_DIAG_FACE_H_DOUBLE" $r/forge_run.log || { echo "$r: 切替が入っている — 止める" >> $LOG; return 1; }; fi
  if [ -z "$nan" ] && [ $rc = 0 ]; then
    ( cd $r && python3 $TOOLS/check_convergence.py . --segment > CONVERGENCE_SEGMENT.txt 2>&1 ); echo "$r check_convergence --segment rc=$?" >> $LOG
  fi
}

arm run_0540_lvcfh_v3_a1 3 0 || exit 1
arm run_0541_lvcfh_v3_b1 3 1 || exit 1
arm run_0542_lvcfh_v3_a2 3 0 || exit 1
arm run_0543_lvcfh_v3_b2 3 1 || exit 1
echo "== 終了 $(date -Is)" >> $LOG
touch lvcfh.done
