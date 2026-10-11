#!/bin/bash
# plan time_integration-line-viscous-jacobian-dt §6.5 (2026-10-11): 全部入り (L7) の再実行 1 本 (run_0567_lvcmp_l7r)。
# ユーザ「やり残したことだった」。L7 (run_0553) と同じ設定・同じバイナリ・同じ出発の場・同じ見張りと記録係。§6.3 の書き出しと事前のゲート (lvcmp_pregate.json) をそのまま使う。
# B に他のセッションの forge が無くなるまで 10 分ごとに待ってから回す (ユーザ判断 2026-10-11 の運用)。判定は lvcmp_rerun_judge.py。
set -uo pipefail
cd "$(dirname "$0")"
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
export FORGE_CUDA_BLOCKSIZE=128
export REAL_CONVERTER=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$PWD/conv_tolerant.sh
BIN=$HOME/forge-fgeom3-fp64/solver_density_cuda/build/forge
KEY=fgeom3_fp64
TOOLS=$(cd ../../solver_density_cuda/tools && pwd)
SRC=run_0183_ns_coldmesh_tw300_ext
LOG=$PWD/lvcmp_rerun.log
NODES=1572,4113,7985,198560,264263
EXTRA=res_ro,res_roUx,res_roUy,res_roUz,res_roe,res_roK,res_roOmega,res_roY0,res_roY1,volume   # 書き出しの run だけ
BUDGET=150000

echo "== 開始 $(date -Is)" >> $LOG
mkdir -p _band_ab/cold_pair
[ "$(sha256sum $BIN | cut -c1-16)" = 129de3f4e7f67aa3 ] || { echo "バイナリの sha256 が登録と違う — 止める" >> $LOG; exit 1; }
[ -e $SRC ] || ln -s $HOME/forge-wallfit/case/45.isobutane_m6_d155/$SRC $SRC
last=$(ls $SRC | grep -E '^res_[0-9]+\.h5$' | sort -t_ -k2 -n | tail -1)
[ "$last" = res_100000.h5 ] && [ "$(sha256sum $SRC/res_100000.h5 | cut -c1-16)" = 207d39f0e7f4aa03 ] || { echo "出発の場が res_100000.h5 (207d39f0…) でない ($last) — 止める" >> $LOG; exit 1; }

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

disk_ok() {  # 共有ディスクの空きが 6 GB 未満なら止める
  local free=$(df --output=avail -B1G ~ | tail -1 | tr -d ' ')
  [ "$free" -ge 6 ] || { echo "ディスクの空き ${free} GB < 6 GB — $1 の前で止める" >> $LOG; return 1; }
}
armargs() {  # armargs <l7|l5|l0>: cold_cfl.py prep に渡す値の引数
  case $1 in l7|l5) echo "--lvc 3";; l0) echo "";; *) return 1;; esac
}
masknum() {  # masknum <l7|l5|l0>: FORGE_LVC_TERMS ("-" = 渡さない)
  case $1 in l7) echo 7;; l5) echo 5;; l0) echo -;; esac
}
modes_ok() {  # modes_ok <run> <l7|l5|l0>: 起動ログの実効のモード
  local f=$1/forge_run.log
  grep -q "^\[line\] factor 1 回目: モード LAYOUT2$" $f && grep -q "^\[line\] solve 1 回目: モード LAYOUT2$" $f || return 1
  ! grep -q "FORGE_DIAG_FACE_H_DOUBLE" $f || return 1
  case $2 in
    l7) grep -q "'lineViscCoupling' in 'time.deltaT': 3" $f && ! grep -q "FORGE_LVC_TERMS=" $f ;;
    l5) grep -q "'lineViscCoupling' in 'time.deltaT': 3" $f && grep -q "診断のマスク FORGE_LVC_TERMS=5 " $f ;;
    l0) ! grep -q "'lineViscCoupling' in 'time.deltaT'" $f && ! grep -q "FORGE_LVC_TERMS=" $f ;;
  esac
}

dump1() {  # dump1 <run> <l7|l5|l0>: 1 step のライン行列の書き出し
  local r=$1 a=$2 mk=$(masknum $2)
  disk_ok $r || return 1
  ( export FORGE_BIN=$BIN COLD_ALT_BINARY=$KEY
    python3 cold_cfl.py prep $SRC $r --steps 1 --out 1 --cfl 4 --limiter-ref-from $SRC --line dir --cap 50 --itj 5 $(armargs $a) --extra $EXTRA >> $LOG 2>&1 ) \
    || { echo "prep $r 失敗 — 止める" >> $LOG; return 1; }
  rm -f $r/nozzle.msh; mkdir -p $r/linedump
  meshsha $r/nozzle.h5 > $r/MESH_SHA.txt || { echo "$r: 格子のハッシュに失敗 — 止める" >> $LOG; return 1; }
  echo $mk > $r/LVC_TERMS.txt
  ( export FORGE_BIN=$BIN COLD_ALT_BINARY=$KEY FORGE_LINE_DUMP_DIR=$PWD/$r/linedump FORGE_LINE_DUMP_CALL=1 FORGE_LINE_DUMP_NODES=$NODES
    [ "$mk" != - ] && export FORGE_LVC_TERMS=$mk
    python3 cold_cfl.py run $r > $r/cold_pair_run_stdout.log 2>&1 ); local rc=$?
  echo "$r rc=$rc $(date -Is) $(grep -c '\[lineDump\] factor の直前を書いた' $r/forge_run.log)" >> $LOG
  [ $rc = 0 ] && grep -q '\[lineDump\] factor の直前を書いた' $r/forge_run.log || { echo "$r: 書き出しが無いか失敗 — 止める" >> $LOG; return 1; }
}

arm() {  # arm <run> <l7|l5|l0>: 水準に入るまで (見張りが止める)
  local r=$1 a=$2 mk=$(masknum $2)
  if [ -f $r/m9_watch.json ] && grep -q '"status": "\(REACHED\|CENSORED\|DIVERGED\)"' $r/m9_watch.json; then echo "$r は既に終わっている" >> $LOG; return 0; fi
  if [ -f $r/m9_watch.json ] && grep -q '"status": "\(EXEC_ERROR\|DATA_ERROR\)"' $r/m9_watch.json; then echo "$r は見張りの失敗で終わっている — 止める" >> $LOG; return 1; fi
  [ -e $r ] && { echo "$r が途中まである — 止める (使い直さない)" >> $LOG; return 1; }
  disk_ok $r || return 1
  ( export FORGE_BIN=$BIN COLD_ALT_BINARY=$KEY
    python3 cold_cfl.py prep $SRC $r --steps $BUDGET --out 2500 --cfl 4 --limiter-ref-from $SRC --line dir --cap 50 --itj 5 $(armargs $a) --extra res_ro >> $LOG 2>&1 ) \
    || { echo "prep $r 失敗 — 止める" >> $LOG; return 1; }
  rm -f $r/nozzle.msh
  meshsha $r/nozzle.h5 > $r/MESH_SHA.txt || { echo "$r: 格子のハッシュに失敗 — 止める" >> $LOG; return 1; }
  echo $mk > $r/LVC_TERMS.txt
  ( export FORGE_BIN=$BIN COLD_ALT_BINARY=$KEY
    [ "$mk" != - ] && export FORGE_LVC_TERMS=$mk
    python3 cold_cfl.py run $r > $r/cold_pair_run_stdout.log 2>&1; echo $? > $r/forge_rc ) &
  local job=$!
  python3 lvcmp_record.py $r >> $r/lvcmp_record.log 2>&1 &    # 見張りが消す前に各出力の sha256・有限性・ρ/P/T の正値を記録する (codex plan-2 m5)
  local recpid=$!
  sleep 120
  if ! modes_ok $r $a; then
    echo "$r: 実効のモードが期待と違う — 止める" >> $LOG
    python3 m9_stop.py $r || echo "$r: forge を止められない" >> $LOG
    wait $job; wait $recpid; return 1
  fi
  echo "$r: モード OK $(date -Is)" >> $LOG
  python3 m9_watch.py $r --phase line --budget $BUDGET >> $LOG 2>&1; local wrc=$?
  local wst=$(python3 -c "import json;print(json.load(open('$r/m9_watch.json'))['status'])" 2>/dev/null)
  case "$wst" in REACHED|CENSORED|DIVERGED|EXEC_ERROR|DATA_ERROR) ;;
    *) echo "$r: 見張りが最終の状態を残さずに終わった (rc $wrc) — forge と記録係を止める" >> $LOG
       python3 m9_stop.py $r || echo "$r: forge を止められない" >> $LOG
       kill $recpid 2>/dev/null ;;
  esac
  wait $job
  wait $recpid
  echo "$r 見張り rc=$wrc、forge rc=$(cat $r/forge_rc 2>/dev/null)、状態 $wst、記録係 $(python3 -c "import json;d=json.load(open('$r/lvcmp_outputs.json'));print(len(d['rows']), d.get('finished'))" 2>/dev/null) $(date -Is)" >> $LOG
  ( cd $r && python3 $TOOLS/check_convergence.py . --segment > CONVERGENCE_SEGMENT.txt 2>&1 ); echo "$r check_convergence --segment rc=$?" >> $LOG
  case "$wst" in REACHED|CENSORED|DIVERGED) return 0;; *) echo "$r: 見張りの状態 $wst — 止める (codex plan-2 m6)" >> $LOG; return 1;; esac
}

mkdir -p _band_ab/cold_pair
python3 lvcmp_pregate.py --verify >> $LOG 2>&1 \
  || { echo "事前のゲートの記録が PASS でないか、証拠が記録と違う — 回さずに止める" >> $LOG; exit 1; }
python3 -c "import sys; sys.path.insert(0, '.'); import cold_xcheck as XC; XC.common_yb()" >> $LOG 2>&1 \
  || { echo "見張りの初期化ができない — 回さずに止める" >> $LOG; exit 1; }
others() {   # 起動の前なので、B の forge は全部他の run (自分の作業ツリーの別の case も数える。codex plan-3 m6)
  pgrep -x forge | wc -l
}
t0=$(date +%s)
while [ "$(others)" -gt 0 ]; do
  if [ $(( $(date +%s) - t0 )) -gt 43200 ]; then echo "他のセッションの forge が 12 時間減らない — 回さずに止める $(date -Is)" >> $LOG; exit 1; fi
  sleep 600
done
echo "== 他のセッションの forge なし $(date -Is): L7r を回す" >> $LOG
arm run_0567_lvcmp_l7r l7 || exit 1
echo "== 終了 $(date -Is)" >> $LOG
touch lvcmp_rerun.done
