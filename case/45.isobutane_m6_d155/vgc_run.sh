#!/bin/bash
# plan tooling-nozzle-isothermal-wall-chain §5.1 #31・§6「V-g0c-chain の事前登録」(commit 68d56723) の実行の入口 (B)。
# 環境変数は float plan §6.61 の sp_arm.sh と同じ (FORGE_CUDA_BLOCKSIZE=128・FORGE_OMEGA_BUDGET=1・FORGE_DIAG_COMMIT_LOSS=500、ほかの FORGE_* は消す)。
# FORGE_BIN は ~/forge-fgeom7-fp64 の FP64 (7e95e5bbd8f50c9e…) で、準備の変換器 (同じビルドの convertGmshToForge) もこれで決まる。
# 中身は vgc_run.py (準備 prepare_ns → 前提の検査 (1)〜(6) → 段階起動 full → 本段だけに m9_watch.py → vgc_judge.py)。
# 先に vgc_prep.sh <commit> を通しておく (作業ツリー ~/forge-g0cchain・入力の写し・台本の配置; vgc_prep.done)。
# 使い方 (B、GPU を占有して; ほかの forge が無くなるまで最大 --wait-max 秒待つ):
#   cd ~/forge-wallfit/case/45.isobutane_m6_d155 && nohup bash vgc_run.sh --commit <commit> > /dev/null 2>&1 < /dev/null &
# 回す前の確かめ (準備と前提の検査 (1)〜(6) だけ、run_0810 を使わず _vgc_dry_<名前>/ に作る。forge・見張り・判定は回さない、占有も待たない):
#   bash vgc_run.sh --commit <commit> --dry-prep _vgc_dry_1     (確かめた後は _vgc_dry_1/ を消してよい)
# 記録: vgc_run.log (この台本と vgc_run.py の出力)、run_0810_vgc_chain/ (VGC_PRECHECK.json・timing.tsv・START_EPOCH・MAIN_STAGE_STARTED・
#   VGC_MAIN_START.json・m9_watch.json・RUN_RC・RUN_END.json)、_band_ab/cold_pair/vgc_judge.json。
# 走行中にこの台本・vgc_run.py を編集しないこと (bash は逐次読み)。
set -uo pipefail
C45=~/forge-wallfit/case/45.isobutane_m6_d155; LOG=$C45/vgc_run.log
TOKEN=$(curl -s -X PUT http://169.254.169.254/latest/api/token -H "X-aws-ec2-metadata-token-ttl-seconds: 60")
[ "$(curl -s -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-id)" = "i-0ba2b91ba659254c4" ] \
  || { echo "B ではない。中止 ($(date -Is))" | tee -a $LOG; exit 1; }
export FORGE_CUDA_BLOCKSIZE=128
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
export FORGE_OMEGA_BUDGET=1 FORGE_DIAG_COMMIT_LOSS=500
export FORGE_BIN=$HOME/forge-fgeom7-fp64/solver_density_cuda/build/forge
cd $C45 || exit 1
[ -e vgc_prep.done ] || { echo "vgc_prep.done が無い (先に vgc_prep.sh)。中止 ($(date -Is))" >> $LOG; exit 1; }
echo "== vgc_run.sh 開始 $(date -Is) 引数 $* 環境 $(compgen -e | grep '^FORGE_' | while read v; do printf '%s=%s ' $v "${!v}"; done)" >> $LOG
python3 vgc_run.py "$@" >> $LOG 2>&1
rc=$?
echo "== vgc_run.sh 終了 $(date -Is) vgc_run.py rc=$rc" >> $LOG
exit $rc
