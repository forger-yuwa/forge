#!/usr/bin/env bash
# #8a: 面寄与ダンプ。usage: r8a.sh prepare | r8a.sh run <case> <bin...> | r8a.sh compare
set -uo pipefail
H=~/forge-pgrad-new/case/09.Taylor-Green/_g0_lsq_seam
SC=~/pgrad6_scratch
declare -A BIN=([new]=$HOME/forge-pgrad-new/solver_density_cuda/build/forge
               [olddump]=$HOME/forge-pgrad-olddump/solver_density_cuda/build/forge
               [wrong]=$HOME/forge-pgrad-wrong/solver_density_cuda/build/forge)
cd $H
case $1 in
  prepare) python3 r8a_facedump.py prepare $SC ;;
  run) c=$2; shift 2; for b in "$@"; do FORGE_BIN=${BIN[$b]} python3 r8a_facedump.py run $SC $c $b; done ;;
  compare) python3 r8a_facedump.py compare $SC --out $H/R8a_facedump.txt ;;
esac
