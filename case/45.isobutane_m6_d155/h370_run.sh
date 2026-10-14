#!/bin/bash
# plan tooling-nozzle-isothermal-wall-chain §5.1 #36f・§6「V-h370 の事前登録」(§4.9 codex diagnose 3) の実行の入口 (B)。
# 環境変数は V-tw370 (tw370_run.sh) と同じ (FORGE_CUDA_BLOCKSIZE=128・FORGE_OMEGA_BUDGET=1・FORGE_DIAG_COMMIT_LOSS=500・
# FORGE_BIN=~/forge-fgeom7-fp64 の FP64 (7e95e5bbd8f50c9e…)、ほかの FORGE_* は消す)。変換器は同じビルドの convertGmshToForge (711f4438…)。
# 中身は h370_run.py ((i) Euler 参照 run_0950 → 資格 → 固定 → (ii) NS run_0951 → 判定)。先に h370_prep.sh <commit> を通しておく。
# nohup で起動する前提 (stdin は /dev/null)。ログは case dir の h370_run.log (追記)。
# 使い方 (B):
#   cd ~/forge-wallfit/case/45.isobutane_m6_d155
#   nohup bash h370_run.sh --phase euler --commit <commit> --dry-run r1 > /dev/null 2>&1 < /dev/null &   # 回す前の確かめ (forge なし)
#   nohup bash h370_run.sh --phase all --commit <commit> > /dev/null 2>&1 < /dev/null &                  # (i) → (ii) → 判定
# 走行中にこの台本・h370_run.py・h370_judge.py を編集しないこと (bash は逐次読み; 判定は終わった後に読み直される)。
set -uo pipefail
C45=~/forge-wallfit/case/45.isobutane_m6_d155
LOG=$C45/h370_run.log
TOKEN=$(curl -s -X PUT http://169.254.169.254/latest/api/token -H "X-aws-ec2-metadata-token-ttl-seconds: 60")
[ "$(curl -s -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-id)" = "i-0ba2b91ba659254c4" ] \
  || { echo "B ではない。中止 ($(date -Is))" | tee -a $LOG; exit 1; }
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
export FORGE_OMEGA_BUDGET=1 FORGE_DIAG_COMMIT_LOSS=500
export FORGE_BIN=$HOME/forge-fgeom7-fp64/solver_density_cuda/build/forge
cd $C45 || exit 1
[ -e h370_prep.done ] || { echo "h370_prep.done が無い (先に h370_prep.sh)。中止 ($(date -Is))" >> $LOG; exit 1; }
echo "== h370_run.sh 開始 $(date -Is) 引数 $* 環境 $(compgen -e | grep '^FORGE_' | while read v; do printf '%s=%s ' $v "${!v}"; done)" >> $LOG
python3 h370_run.py "$@" >> $LOG 2>&1
rc=$?
echo "== h370_run.sh 終了 $(date -Is) h370_run.py rc=$rc" >> $LOG
exit $rc
