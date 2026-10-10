#!/bin/bash
# plan axisymmetric-freestream-hoop-gauge §4.6 の 4 (一様な静止場) の取り直し (2026-10-10)。
# fg9.sh の run_0463〜0468 は、新しい格子 (変換直後) に化学種のデータセット roY0・roY1 が無いまま restart_field.py で場を移したので
# 7 量しか移らず (12 量のうち化学種が落ちた)、still_field.py も roY* を書けなかった。res_0 は 570,878/570,999 節点で Y0 = 1 (純粋な化学種 0)、
# P 92.5 kPa・T 293.6 K で一様でなく、入口 (混合気の Pt 1e5・Tt 300) との食い違いで 10 step 後に 28 m/s の流れが出た (キー 0・1 で同じ値)。台本の誤り。
# ここでは場を移す前に新しい格子へ roY* のデータセットを作り (値は restart_field が上書きする)、移した量が 9 であること、
# still_field の後で P・T・Y が一様であることを確かめてから回す。判定は fg9_an.py (§4 の前提の検査つき)。
set -uo pipefail
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
C45=~/forge-wallfit/case/45.isobutane_m6_d155; W=$C45/_conv9; LOG=$C45/fg9.log
RC=~/forge-wallfit/solver_density_cuda/tools/run_case.sh; TL=~/forge-wallfit/solver_density_cuda/tools
F7_64=~/forge-fgeom7-fp64/solver_density_cuda/build/forge; F7_32=~/forge-fgeom7-f32/solver_density_cuda/build/forge
echo "== fg9d 開始 $(date -Is)" >> $LOG
cd $C45
SRC=run_0183_ns_coldmesh_tw300_ext
PREP64=~/forge-linespeed-fp64/solver_density_cuda/build/forge
export REAL_CONVERTER=~/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$PWD/conv_tolerant.sh
prepnew() {   # prepnew <run> <steps> <out>: B0 の構成で準備し、格子を新しい格子に替えて場を移す (化学種のデータセットを先に作る)
  ( export FORGE_BIN=$PREP64 COLD_ALT_BINARY=lineM_fp64; python3 cold_cfl.py prep $SRC $1 --steps $2 --out $3 --cfl 4 --limiter-ref-from $SRC --line dir --itj 5 --cap 50 > /dev/null ) || { echo "$1: 準備に失敗 (中止)" >> $LOG; exit 1; }
  rm -f $1/nozzle.msh; mv $1/nozzle.h5 $1/nozzle_oldmesh.h5; cp $W/c45/new64/nozzle.h5 $1/nozzle.h5
  python3 - $1/nozzle_oldmesh.h5 $1/nozzle.h5 <<'EOF' || { echo "$1: roY* のデータセットの作成に失敗 (中止)" >> $LOG; exit 1; }
import sys, re, h5py
with h5py.File(sys.argv[1], "r") as s, h5py.File(sys.argv[2], "r+") as d:
    ks = sorted(k for k in s["VALUE"] if re.fullmatch(r"roY\d+", k))
    if not ks: raise SystemExit("SRC に roY* が無い")
    for k in ks:
        if k not in d["VALUE"]: d["VALUE"].create_dataset(k, shape=s["VALUE/" + k].shape, dtype=s["VALUE/" + k].dtype)
    print("roY* を作った:", ks)
EOF
  python3 $TL/restart_field.py $1/nozzle_oldmesh.h5 $1/nozzle.h5 --dst-run $1 --forge $F7_64 > $1/restart_newmesh.log 2>&1 || { echo "$1: 新しい格子への場の移しに失敗 (中止)" >> $LOG; exit 1; }
  grep -q "VERDICT: OK (9 量を移した" $1/restart_newmesh.log || { echo "$1: 移した量が 9 でない: $(grep VERDICT $1/restart_newmesh.log) (中止)" >> $LOG; exit 1; }
  rm -f $1/nozzle_oldmesh.h5
}
key() {   # key <run> <0|1>
  grep -c '^mesh: {' $1/solverConfig.yaml | grep -qx 1 || { echo "$1: mesh の形が想定外 (中止)" >> $LOG; exit 1; }
  sed -i "s/^mesh: {/mesh: {axisSegmentRWeight: $2, /" $1/solverConfig.yaml
}
out() { grep -q '^output' $1/solverConfig.yaml && { echo "$1: output が既にある (中止)" >> $LOG; exit 1; }; echo "output: $2" >> $1/solverConfig.yaml; }
uniform() {   # uniform <run>: still_field の後の保存量が一様 (ro・roe・roY* が全節点で同じ値、運動量 0)
  python3 - $1/nozzle.h5 <<'EOF' >> $LOG 2>&1 || { echo "$1: 一様でない (中止)" >> $LOG; exit 1; }
import sys, re, h5py, numpy as np
with h5py.File(sys.argv[1], "r") as h:
    V = h["VALUE"]; bad = []
    for k in V:
        if k in ("ro", "roe", "roK", "roOmega") or re.fullmatch(r"roY\d+", k):
            a = np.asarray(V[k][:])
            if a.min() != a.max(): bad.append(f"{k} {a.min()}..{a.max()}")
        if k in ("roUx", "roUy", "roUz") and np.any(np.asarray(V[k][:]) != 0): bad.append(k)
    print(f"[uniform] {sys.argv[1]}: 保存量 {sorted(V.keys())}、一様でないもの {bad}")
    if bad or not any(re.fullmatch(r"roY\d+", k) for k in V): sys.exit(1)
EOF
}
run() { local r=$1 b=$2; shift 2; ( export FORGE_BIN=$b; for kv in "$@"; do export "$kv"; done; bash $RC $r > $r/run_case_stdout.log 2>&1 ); echo "$r rc=$? $(basename $(dirname $(dirname $(dirname $b)))) $(grep -h 'axisSegmentRWeight' $r/forge_run.log $r/run_case_stdout.log 2>/dev/null | sort -u | head -1 | cut -c1-140)" >> $LOG; }
STILL="run_0476_hp_still64_k1 run_0477_hp_still64_k0 run_0478_hp_still32_k1 run_0479_hp_still32_k1b run_0480_hp_still32_k0 run_0481_hp_still32_k0b"
for r in $STILL; do
  [ -e $r ] && { echo "$r: 既にある (中止)" >> $LOG; exit 1; }
  prepnew $r 10 1; python3 still_field.py $r --state-from $SRC/res_100000.h5 >> $LOG 2>&1 || { echo "$r: still_field に失敗 (中止)" >> $LOG; exit 1; }
  uniform $r
  case $r in *_k1*) key $r 1;; *) key $r 0;; esac; out $r '{level: 1, extraFields: [res_roUy, res_roUx]}'
done
run $C45/run_0476_hp_still64_k1 $F7_64 FORGE_DIAG_HOOP_CLOSURE=$C45/run_0476_hp_still64_k1/hoop.h5
run $C45/run_0477_hp_still64_k0 $F7_64 FORGE_DIAG_HOOP_CLOSURE=$C45/run_0477_hp_still64_k0/hoop.h5
for r in run_0478_hp_still32_k1 run_0479_hp_still32_k1b run_0480_hp_still32_k0 run_0481_hp_still32_k0b; do run $C45/$r $F7_32 FORGE_DIAG_HOOP_CLOSURE=$C45/$r/hoop.h5; done
for r in $STILL; do for n in 2 3 4 5 6 7 8 9; do rm -f $r/res_$n.h5 $r/res_$n.xmf $r/res_outlet_*_$n.h5 $r/res_outlet_*_$n.xmf; done; rm -f $r/nozzle.h5; done
echo "== fg9d 終了 $(date -Is)" >> $LOG
touch fg9d.done
