#!/bin/bash
# plan time_integration-line-viscous-jacobian-dt §6.3 (2026-10-11): ユーザ判断「他の run が減るまで待つ」の実行。
# L7 (run_0553) の見張り・forge・記録係の終わりを PID で待ち、L7 の区間付きの収束 VERDICT を出したあと、
# B に他のセッションの forge が無くなるまで 10 分ごとに待ってから、lvcmp.sh arms で L0・L5 を順に回す (L7 は終わった状態として飛ばされる)。
# 使い方: lvcmp_chain.sh <見張りの PID> <forge の subshell の PID> <記録係の PID>   (pgrep -f は使わない)
set -uo pipefail
cd "$(dirname "$0")"
LOG=$PWD/lvcmp.log
TOOLS=$(cd ../../solver_density_cuda/tools && pwd)
W=$1 F=$2 R=$3
r=run_0553_lvcmp_l7
echo "== 連鎖の開始 $(date -Is): 見張り $W・forge $F・記録係 $R の終わりを待つ" >> $LOG
for p in $W $F $R; do while kill -0 $p 2>/dev/null; do sleep 30; done; done
wst=$(python3 -c "import json;print(json.load(open('$r/m9_watch.json'))['status'])" 2>/dev/null)
echo "$r 見張り (台本の外で終わった)、forge rc=$(cat $r/forge_rc 2>/dev/null)、状態 $wst、記録係 $(python3 -c "import json;d=json.load(open('$r/lvcmp_outputs.json'));print(len(d['rows']), d.get('finished'))" 2>/dev/null) $(date -Is)" >> $LOG
( cd $r && python3 $TOOLS/check_convergence.py . --segment > CONVERGENCE_SEGMENT.txt 2>&1 ); echo "$r check_convergence --segment rc=$?" >> $LOG
case "$wst" in REACHED|CENSORED|DIVERGED) ;; *) echo "$r: 見張りの状態 $wst — 残りを回さずに止める" >> $LOG; exit 1;; esac
others() {   # 他のセッションの forge の本数 (cwd が自分の作業ツリーの外)
  local n=0
  for p in $(pgrep -x forge); do case "$(readlink /proc/$p/cwd 2>/dev/null)" in $HOME/forge-faceh-audit/*) ;; *) n=$((n+1));; esac; done
  echo $n
}
t0=$(date +%s)
while [ "$(others)" -gt 0 ]; do
  if [ $(( $(date +%s) - t0 )) -gt 43200 ]; then echo "他のセッションの forge が 12 時間減らない — 残りを回さずに止める $(date -Is)" >> $LOG; exit 1; fi
  sleep 600
done
echo "== 他のセッションの forge なし $(date -Is): lvcmp.sh arms で L0・L5 を回す" >> $LOG
exec bash lvcmp.sh arms
