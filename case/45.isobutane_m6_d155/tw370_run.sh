#!/bin/bash
# plan tooling-nozzle-isothermal-wall-chain §5.1 #32・§6「V-tw370 の事前登録」(b5df7416) と「登録の訂正」(b6aeb263) の実行の入口 (B)。腕を 1 つずつ別のプロセスで。
# 環境変数は float plan §6.61 の sp_arm.sh と同じ (FORGE_CUDA_BLOCKSIZE=128・FORGE_OMEGA_BUDGET=1・FORGE_DIAG_COMMIT_LOSS=500、ほかの FORGE_* は消す)。
# FORGE_BIN は ~/forge-fgeom7-fp64 の FP64 (7e95e5bbd8f50c9e…) で、準備の変換器 (同じビルドの convertGmshToForge) もこれで決まる。
# 中身は tw370_run.py (腕 R0・C1: 準備 → 前提の検査 → S1 → S2 → 本段 + m9_watch.py; C50: C1 の受け渡しを待って段なしの本段 + m9_watch.py)。
# 先に tw370_prep.sh <commit> を通しておく (作業ツリー ~/forge-g0cchain の commit・入力の写し・q₃₀₀ の固定; tw370_prep.done)。
# 使い方 (B、GPU は共有。3 本を同時に投入してよい。C50 は C1 が本段に入るまで待つ):
#   cd ~/forge-wallfit/case/45.isobutane_m6_d155
#   for A in R0 C1 C50; do nohup bash tw370_run.sh --arm $A --commit <commit> > /dev/null 2>&1 < /dev/null & sleep 5; done
# 回す前の確かめ (run_08xx を使わない、forge・見張り・判定は回さない; 確かめた後は _tw370_dry_1_* を消してよい):
#   bash tw370_run.sh --arm R0 --commit <commit> --dry-prep _tw370_dry_1_R0
#   bash tw370_run.sh --arm C1 --commit <commit> --dry-prep _tw370_dry_1_C1
#   bash tw370_run.sh --arm C50 --commit <commit> --dry-prep _tw370_dry_1_C50 --c1-dir _tw370_dry_1_C1
# 記録: tw370_run_<腕>.log、run_08xx/ (TW370_PRECHECK.json・timing.tsv・START_EPOCH・MAIN_STAGE_STARTED・TW370_MAIN_START.json・m9_watch.json・
#   RUN_RC・RUN_END.json、C1 は C50_HANDOFF/、C50 は C1_HISTORY/)、_band_ab/cold_pair/tw370_judge.json。
# 走行中にこの台本・tw370_run.py を編集しないこと (bash は逐次読み)。
set -uo pipefail
C45=~/forge-wallfit/case/45.isobutane_m6_d155
ARM=none; prev=""
for x in "$@"; do [ "$prev" = "--arm" ] && ARM=$x; prev=$x; done
case $ARM in R0|C1|C50) ;; *) echo "--arm R0|C1|C50 が要る"; exit 1;; esac
LOG=$C45/tw370_run_$ARM.log
TOKEN=$(curl -s -X PUT http://169.254.169.254/latest/api/token -H "X-aws-ec2-metadata-token-ttl-seconds: 60")
[ "$(curl -s -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-id)" = "i-0ba2b91ba659254c4" ] \
  || { echo "B ではない。中止 ($(date -Is))" | tee -a $LOG; exit 1; }
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
export FORGE_OMEGA_BUDGET=1 FORGE_DIAG_COMMIT_LOSS=500
export FORGE_BIN=$HOME/forge-fgeom7-fp64/solver_density_cuda/build/forge
cd $C45 || exit 1
[ -e tw370_prep.done ] || { echo "tw370_prep.done が無い (先に tw370_prep.sh)。中止 ($(date -Is))" >> $LOG; exit 1; }
echo "== tw370_run.sh 開始 $(date -Is) 引数 $* 環境 $(compgen -e | grep '^FORGE_' | while read v; do printf '%s=%s ' $v "${!v}"; done)" >> $LOG
python3 tw370_run.py "$@" >> $LOG 2>&1
rc=$?
echo "== tw370_run.sh 終了 $(date -Is) tw370_run.py rc=$rc" >> $LOG
exit $rc
