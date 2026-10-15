#!/bin/bash
# plan tooling-nozzle-mknot-fast §6 V4・V5 の最終の版 (cba2914f = final) での証拠 (2026-10-15 result 段 NO-GO の M1)
#   (1) N 2400 の 1 評価の計時: 引き継ぎ時点 9468e0d0 (base) と final を同じ条件で各 3 回 (BLAS 1 スレッド)
#   (2) 1 評価のメモリのピーク (新しいプロセス、N 600・2400)
#   (3) カード 1 枚の生成 (4 並列): M5 c2b (case/70) と M5.5 の最終のカード (M6ノズルの設計 94790833 の card_m55.yaml、mesh.ni 3000 の明示入り)
set -u
cd ~/mknot_time
OUT=~/mknot_time/final_out; mkdir -p $OUT
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 REPS=3
PY=~/mknot_venv/bin/python
{ date -u; lscpu | grep -E "Model name|^CPU\(s\)"; $PY -c "import numpy, scipy, numba, llvmlite; print('numpy', numpy.__version__, 'scipy', scipy.__version__, 'numba', numba.__version__, 'llvmlite', llvmlite.__version__)"; } > $OUT/env.txt 2>&1
for v in base final; do
  (cd $v && CARD71=$HOME/mknot_time/card_m55.yaml $PY -W ignore ../time_eval.py $HOME/mknot_time/$v/design) > $OUT/time_$v.log 2>&1
done
for v in base final; do for w in a b; do for N in 600 2400; do
  CARD71=$HOME/mknot_time/card_m55.yaml /usr/bin/time -f "$v $w $N maxrss_kB %M wall %e" $PY -W ignore mem_one.py ~/mknot_time/$v/design $w $N 2>&1 | grep maxrss >> $OUT/mem.txt
done; done; done
for tag in a b; do
  if [ $tag = a ]; then card=$HOME/mknot_time/final/case/70.m5_air_wt_demo/card_m5_air_demo_c2b.yaml; else card=$HOME/mknot_time/card_m55_final.yaml; fi
  d=$OUT/card_$tag; mkdir -p $d; t0=$(date +%s)
  (cd ~/mknot_time/final/design && /usr/bin/time -v $PY -m forge_design.intake.wt_card $card --out $d/problem.yaml) > $d/gen.log 2>&1 < /dev/null
  echo "final $tag rc=$? wall_s=$(( $(date +%s) - t0 ))" >> $OUT/card_times.txt
done
echo ALL_DONE >> $OUT/card_times.txt
