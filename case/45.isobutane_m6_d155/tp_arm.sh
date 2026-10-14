#!/bin/bash
# plan architecture-float-state-double-geometry §6.62 (事前登録、codex diagnose `2026-10-13-test-p-design-diagnose.md`、2026-10-13 の資格の訂正の後): 試験 P の段 1。
# **実行中にこのファイルを編集しないこと** (bash は台本を読みながら実行するので、編集すると続きの位置がずれて同じ run を回し直すことがある)。
# point の陰解法 (run_0713_pt_d1 の設定) + G0c の yaml 一式 + G0c の Q32 の起点 (_g0c/q32_init.h5)、移植後の FP64 (~/forge-port の build-post-fp64) を両腕に使い、
# 変えるのは lowMachPrecond (0 / 3) の 1 行だけ。GPU を占有して逐次に回す:
#   (c) 資格の 1 step: run_0808_tp_q0 (0)・run_0809_tp_q3 (3)、nStepOuter 1、FORGE_DIAG_STAGE_DUMP (rlg.sh と同じ形)
#   (d) 段 1: run_0806_tp_p0 (0) → run_0807_tp_p3 (3)、各 20,000 step・出力 500 ごと。走行中に出力ごとに tp_region.py で局所の残差と物理性を
#       tp_series.jsonl に取り出し、取り出せたら step 0・500・2500・20000 以外の場 (.h5・.xmf) を消す (取り出しに失敗したら消さずにログ)。
#       待ち合わせの形は forge-25 の cg_pp3.sh (origin/feature/nozzle-core-grid) と同じ (次の出力か RUN_RC ができてから取り出す)。
#   (e) tp_check.py (資格 Q1〜Q6) と tp_judge.py (打ち切りの規則) → tp.done
# 設定は run_0713_pt_d1/solverConfig.yaml から sed で 3 行を書き換える: nStepOuter 120000 → 20000 (資格の run は 1)、deltaT の 2 行目の
#   lowMachPrecond の後に precondEps 0.15・implicitSolvePrecision 0 を足す (腕 3 では lowMachPrecond 3)、output の extraFields を §6.62 の 8 量に。
#   cfl・cfl_pseudo 4.0・implicitRelax 0.7・outStepInterval 500 は元の設定にあることを確かめる (違えば中止)。
# 環境変数は sp_arm.sh と同じ (FORGE_CUDA_BLOCKSIZE=128・FORGE_OMEGA_BUDGET=1・FORGE_DIAG_COMMIT_LOSS=500、ほかの FORGE_* は unset)。
# 資格の run が失敗したとき (RUN_RC が 0 でない、stage.h5 が無い) は段 1 を回さずに判定器だけ回す (資格不成立か判定不能になる)。
# 使い方 (AWS B の ~/forge-wallfit/case/45.isobutane_m6_d155): nohup bash tp_arm.sh > tp_arm.out 2>&1 < /dev/null &
set -uo pipefail
TOKEN=$(curl -s -X PUT http://169.254.169.254/latest/api/token -H "X-aws-ec2-metadata-token-ttl-seconds: 60")
[ "$(curl -s -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-id)" = "i-0ba2b91ba659254c4" ] || { echo "B ではない。中止"; exit 1; }
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
export FORGE_OMEGA_BUDGET=1 FORGE_DIAG_COMMIT_LOSS=500
C45=$(cd ~/forge-wallfit/case/45.isobutane_m6_d155 && pwd -P) || { echo "case/45 が無い。中止"; exit 1; }
LOG=$C45/tp.log; RC=~/forge-wallfit/solver_density_cuda/tools/run_case.sh
BIN=~/forge-port/solver_density_cuda/build-post-fp64/forge
INP=_g0c/q32_init.h5; BASE=run_0713_pt_d1
P0=run_0806_tp_p0; P3=run_0807_tp_p3; Q0=run_0808_tp_q0; Q3=run_0809_tp_q3
STEPS=20000; OUT=500; KEEP=" 0 500 2500 20000 "
MIN_START=$((10 * 1024 * 1024 * 1024)); MIN_RUN=$((2 * 1024 * 1024 * 1024))
EXTRA="[res_ro, res_roUx, res_roUy, res_roe, res_roK, res_roOmega, volume, dt_local]"
cd $C45
die() { echo "中止: $* ($(date -Is))" >> $LOG; exit 1; }
free_b() { df --output=avail -B1 . | tail -1 | tr -dc 0-9; }
wait_free() {   # 2026-10-13 の変更 (§6.62「GPU の共有」): ほかのセッションの forge が続くので占有を待たない。ほかの forge を記録し、メモリの空きだけ確かめる
  local o m
  o=$(for p in $(pgrep -x forge); do readlink /proc/$p/cwd; done | tr '\n' ' ')
  [ -n "$o" ] && echo "ほかの forge と GPU を共有して回す ($(date -Is)): $o" >> $LOG
  m=$(awk '/MemAvailable/ {print int($2 / 1024)}' /proc/meminfo)
  [ "$m" -ge 4608 ] || die "使えるメモリが 4.5 GiB 未満 ($m MiB)"
}
echo "== 開始 $(date -Is) バイナリ $(sha256sum $BIN | cut -c1-16) 入力 $(sha256sum $INP | cut -c1-16) 空き $(( $(free_b) / 1024 / 1024 )) MiB" >> $LOG

# (a) 準備の検査
[ "$(sha256sum $BIN | cut -c1-16)" = "75ee4675b0c90f41" ] || die "FP64 のバイナリが登録 (移植後 75ee4675b0c90f41) と違う"
[ "$(sha256sum $INP | cut -c1-16)" = "6f8fcda6d97482e3" ] || die "入力 $INP が登録 (6f8fcda6d97482e3) と違う"
[ -f $BASE/solverConfig.yaml ] || die "$BASE/solverConfig.yaml が無い"
for r in $P0 $P3 $Q0 $Q3; do [ -e $r ] && die "$r が既にある"; done
for f in bcondConfig.yaml probe.yaml species_meta.yaml resolved_species_fc8ba1a8828147e3.yaml resolved_species_fc8ba1a8828147e3_70cc58aa292c9504.yaml; do
  cmp -s run_0717_c50_d1/$f _g0c/$f || die "$f が run_0717_c50_d1 と _g0c で違う"     # sp_arm.sh と同じ照合
  cmp -s $BASE/$f _g0c/$f && echo "  $f: $BASE と _g0c で同一" >> $LOG || echo "  $f: $BASE と _g0c で違う (記録。登録は G0c の yaml 一式を使う)" >> $LOG
done
B=$BASE/solverConfig.yaml
[ "$(grep -c 'nStepOuter: 120000}' $B)" = 1 ] && [ "$(grep -c 'cfl: 4.0, cfl_pseudo: 4.0,' $B)" = 1 ] \
  && [ "$(grep -c 'implicitRelax: 0.7, lowMachPrecond: 0, detectNaN: 1}' $B)" = 1 ] && [ "$(grep -c 'lowMachPrecond' $B)" = 1 ] \
  && [ "$(grep -c 'outStepInterval: 500$' $B)" = 1 ] && [ "$(grep -c '^output: {level: 1, extraFields: \[' $B)" = 1 ] \
  && [ "$(grep -c 'precondEps\|implicitSolvePrecision' $B)" = 0 ] \
  || die "$BASE の設定の形が想定 (cfl・cfl_pseudo 4.0、implicitRelax 0.7、lowMachPrecond 0、outStepInterval 500、precondEps・implicitSolvePrecision 無し) と違う"
MF=$(grep -o 'meshFileName: "\?[A-Za-z0-9_.]*' $B | sed 's/.*: "\?//')
[ "$MF" = nozzle.h5 ] || die "meshFileName が nozzle.h5 でない ($MF)"
wait_free
[ "$(free_b)" -ge $MIN_START ] || die "ディスクの空きが 10 GiB 未満"
python3 tp_region.py --fix $INP >> $LOG 2>&1 || die "区間の節点の集合を固定できない (tp_region.py --fix)"

# (b) 設定
mkcfg() {   # mkcfg <run> <lowMachPrecond> <nStepOuter>
  local r=$1 m=$2 ns=$3
  mkdir $r || die "$r が作れない"
  cp _g0c/*.yaml $r/ && cp $B $r/solverConfig.yaml && cp $INP $r/nozzle.h5 || die "$r の複製"
  sed -i "s/nStepOuter: 120000}/nStepOuter: $ns}/; s/implicitRelax: 0.7, lowMachPrecond: 0, detectNaN: 1}/implicitRelax: 0.7, lowMachPrecond: $m, precondEps: 0.15, implicitSolvePrecision: 0, detectNaN: 1}/; s/^output: .*/output: {level: 1, extraFields: $EXTRA}/" $r/solverConfig.yaml
  grep -q "nStepOuter: $ns}" $r/solverConfig.yaml && grep -q "implicitRelax: 0.7, lowMachPrecond: $m, precondEps: 0.15, implicitSolvePrecision: 0, detectNaN: 1}" $r/solverConfig.yaml \
    && grep -q "cfl: 4.0, cfl_pseudo: 4.0," $r/solverConfig.yaml && grep -q "outStepInterval: 500$" $r/solverConfig.yaml \
    && grep -q "^output: {level: 1, extraFields: \[res_ro, res_roUx, res_roUy, res_roe, res_roK, res_roOmega, volume, dt_local\]}$" $r/solverConfig.yaml \
    && [ "$(grep -c lowMachPrecond $r/solverConfig.yaml)" = 1 ] && [ "$(sha256sum $r/nozzle.h5 | cut -c1-16)" = "6f8fcda6d97482e3" ] || die "$r の設定の書き換え"
  # 元の設定との差は 3 行 (nStepOuter・deltaT の 2 行目・output) だけ
  local dl; dl=$(diff $B $r/solverConfig.yaml | grep '^>')
  [ "$(diff $B $r/solverConfig.yaml | grep -c '^[<>]')" = 6 ] && [ -z "$(echo "$dl" | grep -v -e "nStepOuter: $ns}" -e "lowMachPrecond: $m, precondEps" -e '^> output: ')" ] || die "$r: $BASE との差が 3 行だけでない"
}
mkcfg $P0 0 $STEPS; mkcfg $P3 3 $STEPS; mkcfg $Q0 0 1; mkcfg $Q3 3 1
# P0 と P3 の設定の差は lowMachPrecond の 1 行だけ (資格の run も同じ)、ほかの yaml は同一
onediff() {   # onediff <a> <b> <a の文字列> <b の文字列>: 差が 1 行だけで、その行が文字列の置き換えだけで違う
  local a=$1 b=$2 sa=$3 sb=$4 la lb
  [ "$(diff $a/solverConfig.yaml $b/solverConfig.yaml | grep -c '^[<>]')" = 2 ] || return 1
  la=$(diff $a/solverConfig.yaml $b/solverConfig.yaml | grep '^<' | cut -c3-); lb=$(diff $a/solverConfig.yaml $b/solverConfig.yaml | grep '^>' | cut -c3-)
  [[ "$la" == *"$sa"* ]] && [[ "$lb" == *"$sb"* ]] && [ "${la/"$sa"/X}" = "${lb/"$sb"/X}" ]
}
onediff $P0 $P3 "lowMachPrecond: 0, " "lowMachPrecond: 3, " || die "P0 と P3 の設定の差が lowMachPrecond の 1 行だけでない"
onediff $Q0 $Q3 "lowMachPrecond: 0, " "lowMachPrecond: 3, " || die "Q0 と Q3 の設定の差が lowMachPrecond の 1 行だけでない"
onediff $P0 $Q0 "nStepOuter: $STEPS}" "nStepOuter: 1}" || die "P0 と Q0 の設定の差が nStepOuter の 1 行だけでない"
onediff $P3 $Q3 "nStepOuter: $STEPS}" "nStepOuter: 1}" || die "P3 と Q3 の設定の差が nStepOuter の 1 行だけでない"
for f in $P0/*.yaml; do b=$(basename $f); [ $b = solverConfig.yaml ] && continue
  for r in $P3 $Q0 $Q3; do cmp -s $f $r/$b || die "$b が $P0 と $r で違う"; done; done
echo "設定: $BASE から 3 行 (nStepOuter・deltaT の 2 行目・output)。P0 と P3 の差は lowMachPrecond の 1 行。$(diff $P0/solverConfig.yaml $P3/solverConfig.yaml | grep '^[<>]' | cut -c1-200 | tr '\n' '|')" >> $LOG

# (c) 資格の 1 step (rlg.sh と同じ形)
QOK=1
for r in $Q0 $Q3; do
  wait_free
  [ "$(free_b)" -ge $MIN_RUN ] || die "ディスクの空きが 2 GiB 未満 ($r の前)"
  echo "$r 起動 $(date -Is) lowMachPrecond $(grep -o 'lowMachPrecond: [0-9]' $r/solverConfig.yaml) 入力 $(sha256sum $r/nozzle.h5 | cut -c1-16)" >> $LOG
  ( export FORGE_BIN=$BIN FORGE_DIAG_STAGE_DUMP=$C45/$r/stage.h5 FORGE_DIAG_STAGE_STEP=1; bash $RC $C45/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC )
  echo "$r 終了 $(date -Is) RUN_RC=$(cat $r/RUN_RC) stage.h5 $([ -e $r/stage.h5 ] && echo あり || echo なし) $(grep -h 'stageDump\|未対応' $r/forge_run.log | head -2 | cut -c1-160)" >> $LOG
  { [ "$(cat $r/RUN_RC)" = 0 ] && [ -e $r/stage.h5 ]; } || QOK=0
done

# (d) 段 1
arm() {   # arm <run>
  local r=$1 f b n left p c FPID TPID
  wait_free
  [ "$(free_b)" -ge $MIN_RUN ] || die "ディスクの空きが 2 GiB 未満 ($r の前)"
  echo "$r 起動 $(date -Is) lowMachPrecond $(grep -o 'lowMachPrecond: [0-9]' $r/solverConfig.yaml) バイナリ $(sha256sum $BIN | cut -c1-16) 入力 $(sha256sum $r/nozzle.h5 | cut -c1-16)" >> $LOG
  date +%s.%N > $r/START_EPOCH
  ( export FORGE_BIN=$BIN; bash $RC $C45/$r > $r/run_case_stdout.log 2>&1; echo $? > $r/RUN_RC ) & FPID=$!
  ( while [ ! -e $r/RUN_RC ]; do   # sp_arm.sh と同じ: 時刻・残差の履歴の最終 step・forge の数・ほかの forge の cwd
      s=$(tail -1 $r/residual_history.csv 2>/dev/null | cut -d, -f1); n=$(pgrep -x forge | wc -l)
      o=$(for p in $(pgrep -x forge); do readlink /proc/$p/cwd; done | grep -vx "$C45/$r" | tr '\n' ' ')
      echo -e "$(date +%s.%N)\t${s:-}\t$n\t$o" >> $r/timing.tsv; sleep 30
    done ) & TPID=$!
  while true; do
    for f in $r/res_*.h5; do
      b=$(basename $f .h5); n=${b#res_}; [[ $n =~ ^[0-9]+$ ]] || continue     # res_outlet_*・res_nan_* は対象外
      [ -f $r/.tp_done_$n ] && continue
      if [ -f $r/res_$((n + OUT)).h5 ] || [ -f $r/RUN_RC ]; then           # 次の出力か終了を待って、書き終わった場だけ読む
        if python3 tp_region.py $r $n >> $r/tp_region.log 2>&1; then
          case "$KEEP" in *" $n "*) ;; *) rm -f $f $r/res_$n.xmf ;; esac
        else
          echo "$r: step $n の取り出しに失敗 (場は残す) $(date -Is)" >> $LOG
        fi
        touch $r/.tp_done_$n
      fi
    done
    if [ -f $r/RUN_RC ]; then
      left=0; for f in $r/res_*.h5; do b=$(basename $f .h5); n=${b#res_}; [[ $n =~ ^[0-9]+$ ]] && [ ! -f $r/.tp_done_$n ] && left=1; done
      [ $left = 0 ] && break
      continue
    fi
    if [ "$(free_b)" -lt $MIN_RUN ]; then   # 自分の forge だけ止める (cwd で絞る。pkill -f は使わない)
      for p in $(pgrep -x forge); do c=$(readlink /proc/$p/cwd); [ "$c" = "$C45/$r" ] && kill $p && echo "ディスクの空きが 2 GiB 未満で止めた: $r (pid $p) $(date -Is)" >> $LOG; done
    fi
    sleep 30
  done
  wait $FPID $TPID
  echo "$r 終了 $(date -Is) RUN_RC=$(cat $r/RUN_RC) 最後の出力 $(ls $r/res_*.h5 2>/dev/null | sed 's/.*res_//; s/\.h5$//' | grep -E '^[0-9]+$' | sort -n | tail -1) 取り出した出力 $(wc -l < $r/tp_series.jsonl 2>/dev/null) 空き $(( $(free_b) / 1024 / 1024 )) MiB" >> $LOG
}
if [ $QOK = 1 ]; then
  arm $P0
  arm $P3
else
  echo "資格の run が失敗した (RUN_RC が 0 でないか stage.h5 が無い)。段 1 を回さない $(date -Is)" >> $LOG
fi

# (e) 判定
echo "== 資格 (tp_check.py) $(date -Is)" >> $LOG; python3 tp_check.py >> $LOG 2>&1; echo "  rc=$?" >> $LOG
echo "== 段 1 の判定 (tp_judge.py) $(date -Is)" >> $LOG; python3 tp_judge.py >> $LOG 2>&1; echo "  rc=$?" >> $LOG
echo "== 終了 $(date -Is)" >> $LOG
touch tp.done
