#!/bin/bash
# case/71 δ* の補正の 1 回目を関門から NS の準備まで間を空けずに流す (AWS A、plan tooling-nozzle-card-local-mesh §5.1 #3e)。
# 起動直後から CPU・GPU の仕事を切らさない (idle 自動停止が関門の抽出の途中で A を止めた 2026-10-15 の対策; 偽のプロセスは使わない)。
set -uo pipefail
C=$HOME/forge-cardmesh/case/71.m55_combustion_wt
R=$C/run_0003_ns_cont
W=$HOME/m55_p1; mkdir -p "$W"
echo "[$(date -Iseconds)] (1) 関門 (上書き 3 つ)" | tee -a "$W/chain.log"
python3 $HOME/m55_diag_0003/gate_0003.py > $HOME/m55_diag_0003/gate_0003.log 2>&1
tail -1 $HOME/m55_diag_0003/gate_0003.log | tee -a "$W/chain.log"
test -f "$R/delta_r_next.csv" || { echo "STOP: 関門が通らない ($HOME/m55_diag_0003/gate_0003.log)" | tee -a "$W/chain.log"; exit 1; }
# (7) 容量: 途中の全場を消す (res_0・140k・160k・180k・200k は残す: #3f の再生用)
n=0; for f in "$R"/res_*.h5; do b=$(basename "$f"); s=${b#res_}; s=${s%.h5}
  [[ "$s" =~ ^[0-9]+$ ]] || continue
  case "$s" in 0|140000|160000|180000|200000) continue ;; esac
  rm -f "$f" "$R/res_$s.xmf"; n=$((n+1)); done
echo "[$(date -Iseconds)] run_0003_ns_cont の途中の全場 $n 枚を消した; $(df -h ~ | tail -1)" | tee -a "$W/chain.log"
bash $HOME/forge-cardmesh/case/71.m55_combustion_wt/_runs/run_0005_ns/p1_chain.sh
