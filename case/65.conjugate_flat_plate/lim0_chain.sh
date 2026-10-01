#!/bin/bash
# case/65 C の事後改訂設定 (limiter 0) チェーン — plan boundary-cht-conjugate-flat-plate.md §4.6 (事前登録どおり)
#   1. cfl 判別 run_0020 (run_0019 最終場・壁温固定、cfl 2、6000 step、毎 step 抜き出し) → ab_judge_cfl2lim0.py で cfl を決める
#   2. 本番再走 run_0022〜0027 を直列 (各本番 run の step 600000 から同一格子 restart、150000 step、連成あり)
#      停止則: run_0022 が STALLED/NaN なら以降を起動しない。どの run でも NaN/発散ならチェーンを止める。
#      延長則: still converging は 1 回だけ 150000 step 延長 (<run>_ext、2 run を連結して収束判定)。
# AWS 上で nohup 実行。進捗は ~/plate_lim0_chain.log、終了で ~/plate_lim0_chain.done
set -u
BIN=/home/ubuntu/forge-graetz/solver_density_cuda/.build-native/relwithdebinfo/forge
MD5=159220e4f9a06ca94df76abc78c49f66
C=~/forge-graetz/case/65.conjugate_flat_plate
T=~/forge-graetz/solver_density_cuda/tools
LOG=~/plate_lim0_chain.log
log() { echo "$(date '+%F %T') $*" | tee -a $LOG; }
[ "$(md5sum $BIN | cut -d' ' -f1)" = $MD5 ] || { log "バイナリの md5 が違う。中止"; exit 1; }
cd ~/forge-graetz && GIT_SSH_COMMAND="ssh -i ~/.ssh/forge-deploy -o StrictHostKeyChecking=no" git fetch -q origin feature/cht-conjugate-benchmarks && git checkout -q -f -B feature/cht-conjugate-benchmarks FETCH_HEAD
TREE=$(git log --oneline -1); log "tree $TREE"

prov() { # $1 run dir, $2 説明
  { echo "forge_bin $BIN"; echo "md5 $MD5"; echo "solver commit 86115cb1 (FP64 build)"; echo "tree $TREE"; echo "$2"; echo "FORGE_CUDA_BLOCKSIZE=128"; date; } > $1/RUN_PROVENANCE.txt
}
nanchk() { # $1 h5 → 0 なら有限
  python3 -c "
import h5py,numpy as np,sys
h=h5py.File('$1'); bad=[k for k in h['VALUE'] if not np.isfinite(h['VALUE/'+k][:]).all()]
print('nonfinite',bad,'Pmin',h['VALUE/P'][:].min(),'romin',h['VALUE/ro'][:].min(),'Tmin',h['VALUE/T'][:].min()); sys.exit(1 if bad else 0)"
}

# ---------- 1. cfl 判別 ----------
R=$C/run_0020_lim0_c1_n64_cfl2; S=$C/run_0019_lim0ext_c1_n64
if [ ! -f $R/AB_DONE ]; then
  mkdir $R || { log "$R が既にある。中止"; exit 1; }
  cp $S/mesh.h5 $S/solid.h5 $S/solverConfig.yaml $S/bcondConfig.yaml $S/probe.yaml $S/conjugate_state_5.h5 $S/wall_profile_5.csv $R/
  python3 $T/restart_field.py $S/res_36000.h5 $R/mesh.h5 --keep-src-dtype | tail -1 | tee -a $LOG
  sed -i "s/nStepOuter: 36000/nStepOuter: 6000/; s/cfl_pseudo: 0.5/cfl_pseudo: 2.0/" $R/solverConfig.yaml
  diff $S/solverConfig.yaml $R/solverConfig.yaml >> $LOG
  prov $R "restart from run_0019_lim0ext_c1_n64 res_36000.h5 (累計 48000) + conjugate_state_5.h5、固定壁温、limiter 0、cfl_pseudo 2.0、6000 step、毎 step 抜き出し"
  F=P,T,Uy,Ux,ro,dt_local,limiter_ro,limiter_Ux,limiter_Uy,limiter_P
  cd $R
  FORGE_CUDA_BLOCKSIZE=128 $BIN > forge_run.log 2>&1 < /dev/null &
  FP=$!
  python3 $C/ab_extract.py $R --until 6000 --fields $F --region-stats > ab_extract.log 2>&1 < /dev/null &
  cd ~
  while kill -0 $FP 2>/dev/null; do sleep 20; done
  while pgrep -f "[a]b_extract.py $R " > /dev/null; do sleep 5; done
  python3 $T/check_convergence.py $R > $R/CONVERGENCE_CHECK.txt 2>&1
  nanchk $R/res_6000.h5 >> $LOG 2>&1 || { log "run_0020 最終場に NaN。チェーン中止"; touch ~/plate_lim0_chain.done; exit 1; }
  (cd $C; python3 ab_judge_cfl2lim0.py $R > $R/AB_JUDGE_CFL2LIM0.txt 2>&1)
  touch $R/AB_DONE
fi
cat $R/ab_extract.log >> $LOG
CFL=$(grep -o "DECISION: cfl_pseudo=[0-9.]*" $R/AB_JUDGE_CFL2LIM0.txt | cut -d= -f2)
[ -n "$CFL" ] || { log "cfl 判別の判定が出ていない (REFUSED?)。チェーン中止"; tail -3 $R/AB_JUDGE_CFL2LIM0.txt >> $LOG; touch ~/plate_lim0_chain.done; exit 1; }
log "cfl 判別: $(grep VERDICT $R/AB_JUDGE_CFL2LIM0.txt) → 本番 cfl_pseudo $CFL"

# ---------- 2. 本番再走 ----------
gates() { # $1 run dir, $2 tol_solid
  local r=$1 tol=$2
  python3 $T/plot_residual.py $r/residual_history.csv -o $r/residual_history.png > /dev/null 2>&1
  python3 ~/forge-graetz/case/64.conjugate_pipe_wall/series_conj.py C $r > $r/SERIES_CONJ.txt 2>&1
  # check_cht_interface.py は run ディレクトリに CHT_INTERFACE_VERDICT.txt を自分で書く。旧閾値を先に回して別名で保存し、
  # 登録式の閾値で正本を書く (2026-10-01 result レビュー m4。初回のチェーンは 2 回目が 1 回目を上書きしていた)
  python3 $T/check_cht_interface.py $r --phys-id 5 --eps-abs 1.0 --eps-rel 1e-3 --dt-k 5e-4 --tol-solid 1e-9 --n-consec 80 > /dev/null 2>&1
  cp $r/CHT_INTERFACE_VERDICT.txt $r/CHT_INTERFACE_VERDICT_old_tol1e-9.txt
  python3 $T/check_cht_interface.py $r --phys-id 5 --eps-abs 1.0 --eps-rel 1e-3 --dt-k 5e-4 --tol-solid $tol --n-consec 80 > /dev/null 2>&1
  sed -i "1i # 登録式の tol_solid $tol (plan §4.6、事後改訂)。旧転記値 1e-9 の結果は CHT_INTERFACE_VERDICT_old_tol1e-9.txt" $r/CHT_INTERFACE_VERDICT.txt
  local Qt=$(python3 -c "
import h5py,glob,re,numpy as np
f=sorted(glob.glob('$r/res_solid_5_*.h5'),key=lambda s:int(re.search(r'_(\d+)\.h5',s).group(1)))[-1]
print(abs(float(np.asarray(h5py.File(f,'r')['VALUE/q_hole'][:]).sum())))")
  python3 $T/check_cht_balance.py $r --solid-mode fem2d --phys-id 5 --phys-name plate --q-floor $Qt --tol-rel 1e-3 --tol-abs $(python3 -c "print($Qt*1e-3)") > $r/CHT_BALANCE_VERDICT.txt 2>&1
  (cd $C; python3 ../64.conjugate_pipe_wall/eval_conj.py C $r > $r/EVAL_CONJ.txt 2>&1; echo "rc=$?" >> $r/EVAL_CONJ.txt)
  log "  $(basename $r): conv [$(grep -o '\-> [A-Za-z ()/—-]*' $r/CONVERGENCE_CHECK.txt | head -1)] | series [$(grep VERDICT $r/SERIES_CONJ.txt | tail -1)] | G-if 登録 [$(grep VERDICT $r/CHT_INTERFACE_VERDICT.txt)] 旧 [$(grep VERDICT $r/CHT_INTERFACE_VERDICT_old_tol1e-9.txt)] | G-cons [$(grep -E 'VERDICT|REFUSED' $r/CHT_BALANCE_VERDICT.txt | tail -1)] | 主判定 [$(grep VERDICT $r/EVAL_CONJ.txt | tail -1)]"
}

runone() { # $1 新 run 名, $2 元 run, $3 元 step, $4 tol_solid, $5 nStepOuter の元の値, $6 説明
  local R=$C/$1 S=$C/$2 st=$3 tol=$4 n0=$5
  mkdir $R || { log "$R が既にある。チェーン中止"; touch ~/plate_lim0_chain.done; exit 1; }
  cp $S/mesh.h5 $S/solid.h5 $S/solverConfig.yaml $S/bcondConfig.yaml $S/probe.yaml $S/conjugate_state_5.h5 $S/RUN_INPUTS.txt $R/
  cp $S/conjugate_Tw_5.csv $R/wall_profile_5.csv
  sed -i "s/ints: {conjugate: 1}/ints: {conjugate: 1, wallProfile: 1}/" $R/bcondConfig.yaml
  python3 $T/restart_field.py $S/res_${st}.h5 $R/mesh.h5 --keep-src-dtype | tail -1 >> $LOG
  sed -i "s/nStepOuter: $n0}/nStepOuter: 150000}/; s/^  limiter: 2/  limiter: 0/; s/cfl_pseudo: [0-9.]*,/cfl_pseudo: $CFL,/; s/tol_solid: [0-9.e-]*,/tol_solid: $tol,/" $R/solverConfig.yaml
  diff $S/solverConfig.yaml $R/solverConfig.yaml >> $LOG; diff $S/bcondConfig.yaml $R/bcondConfig.yaml >> $LOG
  prov $R "$6; restart from $2 res_${st}.h5 + conjugate_state_5.h5 + conjugate_Tw_5.csv (同一格子)、limiter 0、cfl_pseudo $CFL、tol_solid $tol、150000 step"
  log "起動 $1"
  (cd $R; FORGE_CUDA_BLOCKSIZE=128 $BIN > forge_run.log 2>&1 < /dev/null)
  python3 $T/check_convergence.py $R > $R/CONVERGENCE_CHECK.txt 2>&1
  if ! nanchk $R/res_150000.h5 >> $LOG 2>&1 || grep -q "NaN detected" $R/forge_run.log; then log "$1 で NaN/発散。チェーン中止"; touch ~/plate_lim0_chain.done; exit 1; fi
  gates $R $tol
}

extend() { # $1 run 名, $2 tol → still converging なら 1 回だけ延長し、連結した収束判定を書く
  local R=$C/$1
  grep -q "still converging" $R/CONVERGENCE_CHECK.txt || return 0
  local E=$C/${1}_ext
  mkdir $E || { log "$E が既にある"; return 0; }
  cp $R/mesh.h5 $R/solid.h5 $R/solverConfig.yaml $R/bcondConfig.yaml $R/probe.yaml $R/conjugate_state_5.h5 $R/RUN_INPUTS.txt $R/wall_profile_5.csv $E/
  cp $R/conjugate_Tw_5.csv $E/wall_profile_5.csv
  python3 $T/restart_field.py $R/res_150000.h5 $E/mesh.h5 --keep-src-dtype | tail -1 >> $LOG
  prov $E "延長則 (1 回だけ): restart from $1 res_150000.h5 + conjugate_state_5.h5 + conjugate_Tw_5.csv、同設定、150000 step"
  log "延長 ${1}_ext"
  (cd $E; FORGE_CUDA_BLOCKSIZE=128 $BIN > forge_run.log 2>&1 < /dev/null)
  if ! nanchk $E/res_150000.h5 >> $LOG 2>&1; then log "${1}_ext で NaN。チェーン中止"; touch ~/plate_lim0_chain.done; exit 1; fi
  local J=$(mktemp -d)/joined_$1; mkdir -p $J
  { cat $R/residual_history.csv; tail -n +2 $E/residual_history.csv | awk -F, -v OFS=, '{$1=$1+150000; print}'; } > $J/residual_history.csv
  { echo "### 連結判定 ($1 + ${1}_ext、limiter 0 の同一設定区間)"; python3 $T/check_convergence.py $J; } > $E/CONVERGENCE_CHECK.txt 2>&1
  rm -rf $(dirname $J)
  gates $E $2
}

# 名前 元 run 元 step tol_solid 元 nStepOuter 説明
runone run_0022_c1_n64_lim0 run_0007_c1_n64 600000 5.322e-07 600000 "C1 n64"
if grep -q "STALLED" $C/run_0022_c1_n64_lim0/CONVERGENCE_CHECK.txt; then
  log "停止則: run_0022 が STALLED。以降を起動しない"; touch ~/plate_lim0_chain.done; exit 0
fi
extend run_0022_c1_n64_lim0 5.322e-07
runone run_0023_c2_n64_lim0 run_0010_c2_n64 600000 7.832e-07 600000 "C2 n64"; extend run_0023_c2_n64_lim0 7.832e-07
runone run_0024_c1_n16_lim0 run_0005_c1_n16 600000 2.311e-06 600000 "C1 n16"; extend run_0024_c1_n16_lim0 2.311e-06
runone run_0025_c2_n16_lim0 run_0008_c2_n16 600000 3.409e-06 600000 "C2 n16"; extend run_0025_c2_n16_lim0 3.409e-06
runone run_0026_c1_n32_lim0 run_0006_c1_n32 600000 1.094e-06 600000 "C1 n32"; extend run_0026_c1_n32_lim0 1.094e-06
runone run_0027_c2_n32_lim0 run_0009_c2_n32 600000 1.610e-06 600000 "C2 n32"; extend run_0027_c2_n32_lim0 1.610e-06
log "チェーン完了"
touch ~/plate_lim0_chain.done
