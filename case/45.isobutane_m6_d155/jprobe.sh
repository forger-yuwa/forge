#!/bin/bash
# plan time_integration-line-viscous-jacobian §6.13 (2026-10-09): 作用素の切り分け (方向微分)。バイナリ lineG、設定は e1b.sh と同じ (値 3・キー 5・方向別・上限なし・cfl 4)。
# S0 = run_0183 の res_100000 (方向は run_0313 = マスク 7 / run_0315 = マスク 5 の 1 step 目の書き出し)。
# S1 = マスク 7 の新しい run (run_0316_jp_m7_s1) の 20 step 目。その状態で 1 step 目の書き出し (run_0318 = 7 / run_0319 = 5) を取り、同じ手順を繰り返す。
# 各 1 step の評価 run は終了コード・出力の時刻・非有限を確かめ、外れたら全体を止める。重い h5 は npz に抜いてから消す (q0 の res_1 は残す)。
set -uo pipefail
cd "$(dirname "$0")"
export FORGE_CUDA_BLOCKSIZE=128 REAL_CONVERTER=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$PWD/conv_tolerant.sh
export FORGE_BIN=$(cat prof_target_linevisc.txt) COLD_ALT_BINARY=lineG_fp64
unset FORGE_ALLOW_UNVERIFIED_SPECIES FORGE_PROFILE FORGE_LINE_COMPARE FORGE_LINE_DUMP_DIR FORGE_DUMP_LEDGER FORGE_LVC_TERMS FORGE_WI_FORCE_DIAG
SRC=run_0183_ns_coldmesh_tw300_ext
B="--cfl 4 --limiter-ref-from $SRC --line dir --itj 5 --lvc 3"
EXTRA=res_ro,res_roUx,res_roUy,res_roUz,res_roe,wi_eheat,wi_ework
EPS=1e-2
LOG=$PWD/jprobe.log
mkdir -p _jprobe/logs
fail() { echo "失敗: $*" | tee -a $LOG; touch jprobe.done; exit 1; }

# 1 step の評価: probe <run> <field_from|-> [keep]
probe() {
  local r=$1 ff=$2 keep=${3:-}
  local opt=""; [ "$ff" != "-" ] && opt="--field-from $ff"
  python3 cold_cfl.py prep $SRC $r --steps 1 --out 1 $B --extra $EXTRA $opt > _jprobe/logs/$r.prep.log 2>&1 || return 1
  rm -f $r/nozzle.msh
  local t0; t0=$(date +%s)
  ( export FORGE_WI_FORCE_DIAG=1; python3 cold_cfl.py run $r > $r/cold_pair_run_stdout.log 2>&1 ); local rc=$?
  echo "$r rc=$rc" >> $LOG
  [ $rc -eq 0 ] || return 1
  [ -f $r/res_1.h5 ] && [ $(stat -c %Y $r/res_1.h5) -ge $t0 ] || { echo "$r: res_1 が無いか古い" >> $LOG; return 1; }
  python3 jprobe.py extract $r >> $LOG 2>&1 || return 1
  rm -f $r/nozzle.h5 $r/res_0.h5
  [ -z "$keep" ] && rm -f $r/res_1.h5
  return 0
}
# 摂動した場を作って 5 本 (±ε、±ε/2、緩和後の実補正) を評価: family <run 接頭辞> <tag>
family() {
  local pre=$1 tag=$2
  for k in pe me ph mh pp; do
    probe ${pre}_${k} _jprobe/${tag}_${k} || fail "${pre}_${k}"
    rm -f _jprobe/${tag}_${k}/res_0.h5          # nozzle.h5 へ写した後は不要 (sha256 は fields の記録にある)
  done
}

echo "== 開始 $(date -Is)" >> $LOG
# ---- S0 ----
python3 jprobe.py fields $SRC/res_100000.h5 run_0313_e1bdump_m7_linedump s0p7 --eps $EPS >> $LOG 2>&1 || fail "fields s0p7"
python3 jprobe.py fields $SRC/res_100000.h5 run_0315_e1bdump_m5_linedump s0p5 --eps $EPS >> $LOG 2>&1 || fail "fields s0p5"
probe run_0317_jp_s0_q0 - keep || fail q0
probe run_0317_jp_s0_q0b - || fail q0b
family run_0317_jp_s0p7 s0p7
family run_0317_jp_s0p5 s0p5
R0="q0=run_0317_jp_s0_q0,q0b=run_0317_jp_s0_q0b"
OPS0="7=run_0313_e1bdump_m7_linedump,5=run_0315_e1bdump_m5_linedump"
for t in p7 p5; do
  d=$([ $t = p7 ] && echo run_0313_e1bdump_m7_linedump || echo run_0315_e1bdump_m5_linedump)
  python3 jprobe.py compare $d s0$t --eps $EPS --ops $OPS0 \
    --runs "$R0,pe=run_0317_jp_s0${t}_pe,me=run_0317_jp_s0${t}_me,ph=run_0317_jp_s0${t}_ph,mh=run_0317_jp_s0${t}_mh,pp=run_0317_jp_s0${t}_pp" > _jprobe/s0${t}_compare.txt 2>&1 || fail "compare s0$t"
done
echo "== S0 終了 $(date -Is)" >> $LOG

# ---- S1: マスク 7 を 20 step ----
r=run_0316_jp_m7_s1
python3 cold_cfl.py prep $SRC $r --steps 20 --out 10 $B --extra $EXTRA > _jprobe/logs/$r.prep.log 2>&1 || fail "prep $r"
rm -f $r/nozzle.msh
python3 cold_cfl.py run $r > $r/cold_pair_run_stdout.log 2>&1; rc=$?; echo "$r rc=$rc" >> $LOG
[ $rc -eq 0 ] && [ -f $r/res_20.h5 ] || fail "$r"
python3 jprobe.py locate $r >> $LOG 2>&1 || fail "locate"
NODE=$(python3 -c "import json;print(json.load(open('_jprobe/${r}_locate.json'))['node'])")
for m in 7 5; do
  d=$([ $m = 7 ] && echo run_0318_jpdump_s1_m7 || echo run_0319_jpdump_s1_m5)
  mkdir ${d}_linedump || fail "${d}_linedump が既にある"
  python3 cold_cfl.py prep $SRC $d --steps 1 --out 1 $B --extra res_ro --field-from $r > _jprobe/logs/$d.prep.log 2>&1 || fail "prep $d"
  rm -f $d/nozzle.msh
  ( export FORGE_LVC_TERMS=$m FORGE_LINE_DUMP_DIR=$PWD/${d}_linedump FORGE_LINE_DUMP_CALL=1 FORGE_LINE_DUMP_NODES=1572,4113,7985,198560,264263,$NODE
    python3 cold_cfl.py run $d > $d/cold_pair_run_stdout.log 2>&1 ); rc=$?; echo "$d rc=$rc" >> $LOG
  [ $rc -eq 0 ] && [ -f ${d}_linedump/dqnew_s4.f64 ] || fail "$d"
  rm -f $d/nozzle.h5
done
# 2 つの書き出しが同じ状態・同じ rhs・同じ D (行 0〜3) であることを確かめる
python3 - >> $LOG 2>&1 <<'EOF' || fail "S1 の書き出しの照合"
import numpy as np, sys
sys.path.insert(0, "."); import jprobe as J
from pathlib import Path
a = J.load_dump(Path("run_0318_jpdump_s1_m7_linedump")); b = J.load_dump(Path("run_0319_jpdump_s1_m5_linedump"))
n = len(a["_nodes"])
ok = np.array_equal(a["node_line"], b["node_line"]) and np.array_equal(a["state_ro_roU_roe_cp_gamma"], b["state_ro_roU_roe_cp_gamma"])
drhs = float(np.max(np.abs(a["rhs_s0"] - b["rhs_s0"])) / np.max(np.abs(a["rhs_s0"])))
Da = a["D"].reshape(n, 5, 5); Db = b["D"].reshape(n, 5, 5)
dD = float(np.max(np.abs(Da[:, :4] - Db[:, :4])))
print(f"[S1 照合] 節点・状態の一致 {ok}、rhs の相対差 {drhs:.2e}、D の行 0〜3 の差 {dD:.2e}、ライン {sorted(set(a['_lines'].tolist()))}")
sys.exit(0 if ok and drhs <= 1e-9 and dD == 0 else 1)
EOF
python3 jprobe.py fields $r/res_20.h5 run_0318_jpdump_s1_m7_linedump s1p7 --eps $EPS >> $LOG 2>&1 || fail "fields s1p7"
python3 jprobe.py fields $r/res_20.h5 run_0319_jpdump_s1_m5_linedump s1p5 --eps $EPS >> $LOG 2>&1 || fail "fields s1p5"
probe run_0320_jp_s1_q0 $r keep || fail "s1 q0"
probe run_0320_jp_s1_q0b $r || fail "s1 q0b"
family run_0320_jp_s1p7 s1p7
family run_0320_jp_s1p5 s1p5
R1="q0=run_0320_jp_s1_q0,q0b=run_0320_jp_s1_q0b"
OPS1="7=run_0318_jpdump_s1_m7_linedump,5=run_0319_jpdump_s1_m5_linedump"
for t in p7 p5; do
  d=$([ $t = p7 ] && echo run_0318_jpdump_s1_m7_linedump || echo run_0319_jpdump_s1_m5_linedump)
  python3 jprobe.py compare $d s1$t --eps $EPS --ops $OPS1 \
    --runs "$R1,pe=run_0320_jp_s1${t}_pe,me=run_0320_jp_s1${t}_me,ph=run_0320_jp_s1${t}_ph,mh=run_0320_jp_s1${t}_mh,pp=run_0320_jp_s1${t}_pp" > _jprobe/s1${t}_compare.txt 2>&1 || fail "compare s1$t"
done
echo "== 終了 $(date -Is)" >> $LOG
touch jprobe.done
