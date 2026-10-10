#!/bin/bash
# plan time_integration-line-viscous-jacobian-dt §6 (2026-10-11): 値 3・マスク 7 (全部入り) の方向別 dt の上限の A/B。
# バイナリは全部元のセッションの段 ③ の FP64 (~/forge-fgeom3-fp64、読むだけ)。腕の違いは方向別 dt の上限だけ:
#   A = --line dir (上限なし)、B = --line dir --cap 50、C (条件付き) = --line only (方向別なし = point の dt)。FORGE_LVC_TERMS=7 は全腕で同じ。
# run_0183 の res_100000 から、値 3・キー 5・cfl 4 (faceh §6.9・§6.11 と同じ)。
# 順序: ① 1 step の書き出しを A・B・A の再実行・C で 4 本 (残差の場つき) → ② 事前のゲート lvcdt_pregate.py → 合格 (終了コード 0) のときだけ
# ③ A1・B1・A2・B2 を最大 2000 step・100 ごと、序盤 200 step は 112 節点の帳簿。ゲートが INVALID (1) か判別不能 (2) なら腕を回さずに止める。
# 本判定は lvcdt_judge.py (この台本は判定しない)。引数: なし = ①〜③、dumps = ①② だけ、arms = ③ だけ (lvcdt_pregate.py --verify が通るときだけ)、
# armsC = C1・C2 だけ (A/B だけの再判定 lvcdt_judge.py --ab-only が分岐 2「回避せず」で、--verify が通るときだけ)。
# AWS の自分の作業ツリー (~/forge-faceh-audit/case/45.isobutane_m6_d155) で動かす。run_0183 は ~/forge-wallfit の run へのリンク。
# 発散 (detectNaN で止まる、RUN_RC 1) は結果なので台本は止めずに次へ進む。起動の失敗・表示の不一致では止める。
set -uo pipefail
cd "$(dirname "$0")"
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
export REAL_CONVERTER=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$PWD/conv_tolerant.sh
BIN=$HOME/forge-fgeom3-fp64/solver_density_cuda/build/forge
TOOLS=$(cd ../../solver_density_cuda/tools && pwd)
SRC=run_0183_ns_coldmesh_tw300_ext
LOG=$PWD/lvcdt.log
NODES=1572,4113,7985,198560,264263
EXTRA=res_ro,res_roUx,res_roUy,res_roUz,res_roe,res_roK,res_roOmega,res_roY0,res_roY1,volume   # 書き出しの run だけ (残差の不変を全節点で見る)

echo "== 開始 $(date -Is)" >> $LOG
[ "$(sha256sum $BIN | cut -c1-16)" = 129de3f4e7f67aa3 ] || { echo "バイナリの sha256 が登録と違う — 止める" >> $LOG; exit 1; }
[ -e $SRC ] || ln -s $HOME/forge-wallfit/case/45.isobutane_m6_d155/$SRC $SRC
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

disk_ok() {  # 共有ディスクの空きが 4 GB 未満なら止める (2026-10-10: 他のセッションの run で 1 時間に 17 GB 減った。FINITE の腕 1 本は約 3.3 GB)
  local free=$(df --output=avail -B1G ~ | tail -1 | tr -d ' ')
  [ "$free" -ge 4 ] || { echo "ディスクの空き ${free} GB < 4 GB — $1 の前で止める" >> $LOG; return 1; }
}

lineargs() {  # lineargs <腕 a|b|c>: cold_cfl.py prep に渡すライン・dt の引数
  case $1 in a) echo "--line dir";; b) echo "--line dir --cap 50";; c) echo "--line only";; *) return 1;; esac
}

dump1() {  # dump1 <run> <腕 a|b|c>: 1 step のライン行列の書き出し (介入の成立の確認)
  local r=$1 arm=$2 mk=7 bin=$BIN key=fgeom3_fp64
  disk_ok $r || return 1
  ( export FORGE_BIN=$bin COLD_ALT_BINARY=$key
    python3 cold_cfl.py prep $SRC $r --steps 1 --out 1 --cfl 4 --limiter-ref-from $SRC $(lineargs $arm) --extra $EXTRA --itj 5 --lvc 3 >> $LOG 2>&1 ) \
    || { echo "prep $r 失敗 — 止める" >> $LOG; return 1; }
  rm -f $r/nozzle.msh; mkdir -p $r/linedump
  meshsha $r/nozzle.h5 > $r/MESH_SHA.txt || { echo "$r: 格子のハッシュに失敗 — 止める" >> $LOG; return 1; }
  echo $mk > $r/LVC_TERMS.txt
  ( export FORGE_BIN=$bin COLD_ALT_BINARY=$key FORGE_LVC_TERMS=$mk FORGE_LINE_DUMP_DIR=$PWD/$r/linedump FORGE_LINE_DUMP_CALL=1 FORGE_LINE_DUMP_NODES=$NODES
    python3 cold_cfl.py run $r > $r/cold_pair_run_stdout.log 2>&1 ); local rc=$?
  echo "$r rc=$rc $(date -Is) $(grep -c '\[lineDump\] factor の直前を書いた' $r/forge_run.log)" >> $LOG
  [ $rc = 0 ] && grep -q '\[lineDump\] factor の直前を書いた' $r/forge_run.log || { echo "$r: 書き出しが無いか失敗 — 止める" >> $LOG; return 1; }
}

arm() {  # arm <run> <腕 a|b|c>
  local r=$1 armk=$2 mk=7 v=3 fh=0 bin=$BIN key=fgeom3_fp64
  disk_ok $r || return 1
  export FORGE_BIN=$bin COLD_ALT_BINARY=$key
  python3 cold_cfl.py prep $SRC $r --steps 2000 --out 100 --cfl 4 --limiter-ref-from $SRC $(lineargs $armk) --extra res_ro,volume --itj 5 --lvc $v >> $LOG 2>&1 \
    || { echo "prep $r 失敗 — 止める" >> $LOG; return 1; }
  rm -f $r/nozzle.msh
  meshsha $r/nozzle.h5 > $r/MESH_SHA.txt || { echo "$r: 格子のハッシュに失敗 — 止める" >> $LOG; return 1; }
  echo $mk > $r/LVC_TERMS.txt
  ( export FORGE_LVC_TERMS=$mk FORGE_DUMP_LEDGER=$PWD/$r/ledger.csv FORGE_DUMP_LEDGER_NODES=$LEDGER_NODES FORGE_DUMP_LEDGER_CALLS=200
    [ "$fh" = 1 ] && export FORGE_DIAG_FACE_H_DOUBLE=1; python3 cold_cfl.py run $r > $r/cold_pair_run_stdout.log 2>&1 ); local rc=$?
  grep -q "Write failed" $r/forge_run.log && { echo "$r: 書き込みに失敗 (ディスク) — 止める" >> $LOG; return 1; }
  local nan=$(grep -m1 -o "Non-finite value detected in '[A-Za-z0-9_]*' at step [0-9]*" $r/forge_run.log)
  echo "$r v=$v arm=$armk rc=$rc $(date -Is) last=$(awk -F, '$3=="outer_begin"{s=$1} END{print s}' $r/residual_history.csv 2>/dev/null) ${nan:-有限}" >> $LOG
  [ -s $r/residual_history.csv ] && grep -q "'lineViscCoupling' in 'time.deltaT': $v" $r/forge_run.log \
    || { echo "$r: 起動していないか値の表示が違う — 止める" >> $LOG; return 1; }
  if [ $mk = 7 ]; then ! grep -q "FORGE_LVC_TERMS=" $r/forge_run.log || { echo "$r: マスク 7 なのに表示がある — 止める" >> $LOG; return 1; }
  else grep -q "診断のマスク FORGE_LVC_TERMS=$mk " $r/forge_run.log || { echo "$r: マスク $mk の表示が無い — 止める" >> $LOG; return 1; }; fi
  if [ "$fh" = 1 ]; then grep -q "FORGE_DIAG_FACE_H_DOUBLE" $r/forge_run.log || { echo "$r: 切替の表示が無い — 止める" >> $LOG; return 1; }
  else ! grep -q "FORGE_DIAG_FACE_H_DOUBLE" $r/forge_run.log || { echo "$r: 切替が入っている — 止める" >> $LOG; return 1; }; fi
  if [ -z "$nan" ] && [ $rc = 0 ]; then
    ( cd $r && python3 $TOOLS/check_convergence.py . --segment > CONVERGENCE_SEGMENT.txt 2>&1 ); echo "$r check_convergence --segment rc=$?" >> $LOG
  fi
}

MODE=${1:-all}
if [ $MODE = all ] || [ $MODE = dumps ]; then
dump1 run_0584_lvcdt_a_dump a || exit 1
dump1 run_0585_lvcdt_b_dump b || exit 1
dump1 run_0586_lvcdt_a_dump2 a || exit 1
dump1 run_0587_lvcdt_c_dump c || exit 1
python3 lvcdt_pregate.py > lvcdt_pregate.stdout 2>&1; prc=$?
echo "事前のゲート rc=$prc $(date -Is) $(grep '^VERDICT' lvcdt_pregate.stdout)" >> $LOG
[ $prc = 0 ] || { echo "事前のゲートが不合格 (1 = INVALID / 2 = 判別不能) — 腕を回さずに止める" >> $LOG; touch lvcdt.done; exit 1; }
fi
[ $MODE = dumps ] && { echo "== 書き出しと事前のゲートまで $(date -Is)" >> $LOG; touch lvcdt.done; exit 0; }
python3 lvcdt_pregate.py --verify >> $LOG 2>&1 \
  || { echo "事前のゲートの記録が PASS でないか、証拠が記録と違う — 腕を回さずに止める" >> $LOG; touch lvcdt.done; exit 1; }
if [ $MODE = armsC ]; then
  # C の起動の直前に、今の証拠で A/B だけを判定し直す (ゲート合格・INVALID なし・主判定「回避せず」のときだけ終了コード 0。plan-dt レビュー M5)
  python3 lvcdt_judge.py --ab-only > lvcdt_judge_ab.stdout 2>&1 \
    || { echo "A/B の再判定が分岐 2 (回避せず) でないか、ゲートが外れた — C を回さずに止める" >> $LOG; touch lvcdt.done; exit 1; }
  arm run_0588_lvcdt_c1 c || exit 1
  arm run_0589_lvcdt_c2 c || exit 1
  echo "== C の終了 $(date -Is)" >> $LOG
  touch lvcdt.done; exit 0
fi
arm run_0580_lvcdt_a1 a || exit 1
arm run_0581_lvcdt_b1 b || exit 1
arm run_0582_lvcdt_a2 a || exit 1
arm run_0583_lvcdt_b2 b || exit 1
echo "== 終了 $(date -Is)" >> $LOG
touch lvcdt.done
