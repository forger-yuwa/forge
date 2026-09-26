#!/usr/bin/env bash
# #8c: g_suite 7 変種を新バイナリ (f67fe877、forge-pgrad-new の既定バイナリ、ダンプ off) で再実行 (本 run + GPU 定数場 run、定数場は 3 区分判定)
set -uo pipefail
export PATH=$HOME/venv-mesh/bin:$PATH
cd ~/forge-pgrad-new/case/09.Taylor-Green/_g0_lsq_seam
for v in ${VARIANTS:-tgv tgv_mirror tgv_shift tgv_repeat tgv_bcswap jitter32 channel}; do
  python3 g_suite.py $v --scratch ~/pgrad6_scratch > ~/pgrad8/gs8_$v.log 2>&1
  echo "$v rc=$?"; grep "^VERDICT" G_$v.txt
done
echo G8C_DONE
