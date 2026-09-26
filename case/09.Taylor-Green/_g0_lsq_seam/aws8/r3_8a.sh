#!/usr/bin/env bash
# #8a の R3 再確認: 既定 (ダンプ off) の新 f67fe877 が HEAD 前 bd22376d (solver 差分はダンプのみ) と同じ扱いか。
# 入力は前回 R3 の r3_prep (r3_prepare.py で作ったもの、ローカルから転送) をそのまま 4 run へ cp。旧 = bd22376d、新 = f67fe877。
set -uo pipefail
SC=~/pgrad8
H=~/forge-pgrad-new/case/09.Taylor-Green/_g0_lsq_seam
declare -A BIN=([old]=$HOME/forge-pgrad-bd/solver_density_cuda/build/forge [new]=$HOME/forge-pgrad-new/solver_density_cuda/build/forge)
for c in case48 case16 axi; do
  for v in old_a old_b new_a new_b; do
    R=$SC/r3_${c}_$v
    [ -e $R/res_1.h5 ] && { echo "exists $R"; continue; }
    mkdir -p $R
    for f in $SC/r3_prep/$c/*; do case $f in *.msh|*.geo|*.xmf|*.log|*quality.txt) ;; *) cp $f $R/ ;; esac; done
    env -u FORGE_DUMP_SCALARGRAD FORGE_BIN=${BIN[${v%_*}]} FORGE_CUDA_BLOCKSIZE=128 FORGE_CUDA_BLOCKSIZE_SMALL=128 \
      ~/forge-pgrad-new/solver_density_cuda/tools/run_case.sh $R > $R/run_case_stdout.log 2>&1
    echo "$c $v rc=$? $(grep forge_sha256 $R/RUN_PROVENANCE.txt)"
  done
done
cd $H && python3 r3_compare.py $SC --out $H/R8a_R3_dumpoff.txt | grep -n "VERDICT\|NaN\|FAIL"
echo R3_8A_DONE
