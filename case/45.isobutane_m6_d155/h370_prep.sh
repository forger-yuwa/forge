#!/bin/bash
# plan tooling-nozzle-isothermal-wall-chain §5.1 #36f・§6「V-h370 の事前登録」(§4.9 codex diagnose 3) の準備 (B)。forge は回さない。
#  (a) 作業ツリー ~/forge-g0cchain: ~/forge-wallfit の git で `git fetch origin feature/nozzle-wall-fit-and-pipeline` し、引数の commit (この台本を
#      commit した後の commit。fetch したブランチの祖先) に合わせる。無ければ `git worktree add --detach`。あって HEAD が違えば、commit されて
#      いない変更が無く、作業ツリーを読む台本 (vgc_*・tw370_*・h370_* の run・judge) が走っていないときだけ `git checkout --detach <commit>` で
#      動かす (元の HEAD を記録)。台本 (h370_prep.sh・h370_run.sh・h370_run.py・h370_judge.py・h370_problems.py) と入力 (問題 YAML 2 本・
#      δ_r の CSV) を作業ツリーから case dir へ置く (既にあって中身が違えば止める。消さない)。
#  (b) バイナリ (forge 7e95e5bbd8f50c9e・変換器 711f4438)・python3 の numpy 2.5.2 と h5py・yaml・scipy・matplotlib (判定の nozzle_report が読む)。
#  (c) B にある入力: CFD ピンの凍結元 run_0062 (問題 YAML の initial_line_run)・腕 A の旧い Euler run_0174 (prepare_info.json・nozzle.h5・
#      最後の res)・V5d の起動の記録 _band_ab/moc_v5d_run.log (h370_problems.py check が読む)。足りなければ写す手順を出して止める (exit 3)。
#  (d) h370_problems.py check を作業ツリーの版で (h370_run.py --phase check-problems → _band_ab/half_inlet/h370_problems_B.json)。
#  (e) 既にある run・準備・固定の記録 (h370_run.py は既存のものを使わずに止まる)、資源の記録。最後まで通ったら h370_prep.done (commit)。
# 使い方 (B、作業ツリーの commit を引数に):
#   git -C ~/forge-wallfit fetch origin feature/nozzle-wall-fit-and-pipeline
#   git -C ~/forge-wallfit show <commit>:case/45.isobutane_m6_d155/h370_prep.sh > /tmp/h370_prep.sh && bash /tmp/h370_prep.sh <commit>
set -uo pipefail
COMMIT=${1:-}
REPO=~/forge-wallfit; WT=~/forge-g0cchain; CASE=case/45.isobutane_m6_d155
C45=$REPO/$CASE; WT45=$WT/$CASE; LOG=$C45/h370_prep.log
BRANCH=feature/nozzle-wall-fit-and-pipeline
A_IP=172.31.5.94; B_IP=172.31.14.139; B_ID=i-0ba2b91ba659254c4
BIN=~/forge-fgeom7-fp64/solver_density_cuda/build/forge; CONV=~/forge-fgeom7-fp64/solver_density_cuda/build/convertGmshToForge
R062=run_0062_euler_wallfit_fit_r1_ext6k; R174=run_0174_euler_v5d_M_r1
R062_FILES="prepare_info.json nozzle.h5 res_6000.h5 solverConfig.yaml bcondConfig.yaml"
SCRIPTS="h370_prep.sh h370_run.sh h370_run.py h370_judge.py h370_problems.py"
INPUTS="problem_d155_euler_h_ref.yaml problem_d155_ns_prod_coldmesh_tw370_cgref_h_b0.yaml _band_ab/half_inlet/delta_r_initial_C.csv"
CSV_SHA=24251cb647226d6523109655c08c05210fdb8ba4a669aeaecb4407202945b918
V5D_LINE="problem_d155_euler_v5d_M.yaml cfee7f057648c21ef1cfd51a167397c4490ab2ee6a6f8fee31df7fcf84fdaf92"
log() { echo "$*" | tee -a $LOG; }
die() { log "中止: $* ($(date -Is))"; exit 1; }
sha16() { sha256sum "$1" 2>/dev/null | cut -c1-16; }
# 作業ツリーを読む台本が走っているか (/proc の cmdline を見る; pgrep -f は自分に一致するので使わない)
readers() {
  local p cmd
  for p in /proc/[0-9]*; do
    [ "${p#/proc/}" = "$$" ] && continue
    cmd=$(tr '\0' ' ' < $p/cmdline 2>/dev/null) || continue
    case "$cmd" in
      *python*vgc_run.py*|*python*vgc_judge.py*|*python*tw370_run.py*|*python*tw370_judge.py*|*python*h370_run.py*|*python*h370_judge.py*)
        echo "${p#/proc/}: $cmd";;
    esac
  done
}

TOKEN=$(curl -s -X PUT http://169.254.169.254/latest/api/token -H "X-aws-ec2-metadata-token-ttl-seconds: 60")
[ "$(curl -s -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-id)" = "$B_ID" ] || { echo "B ではない。中止"; exit 1; }
[[ "$COMMIT" =~ ^[0-9a-f]{7,40}$ ]] || { echo "使い方: bash h370_prep.sh <commit (16 進)>"; exit 1; }
[ -d $C45 ] || { echo "$C45 が無い"; exit 1; }
rm -f $C45/h370_prep.done
log "== h370_prep.sh 開始 $(date -Is) commit $COMMIT (台本 $0 sha $(sha16 $0))"

# --- (a) 作業ツリー ---------------------------------------------------------------------------------------------------
git -C $REPO fetch origin $BRANCH >> $LOG 2>&1 || die "git fetch origin $BRANCH"
FULL=$(git -C $REPO rev-parse --verify --quiet "$COMMIT^{commit}") || die "commit $COMMIT が $REPO に無い"
git -C $REPO merge-base --is-ancestor $FULL FETCH_HEAD || die "commit $FULL が origin/$BRANCH ($(git -C $REPO rev-parse FETCH_HEAD)) の祖先でない"
if [ -e $WT ]; then
  OLD=$(git -C $WT rev-parse HEAD 2>/dev/null) || die "$WT が git の作業ツリーでない"
  [ -z "$(git -C $WT status --porcelain --untracked-files=no)" ] || die "既存の $WT に commit されていない変更がある (消さない)"
  if [ "$OLD" != "$FULL" ]; then
    RUNNING=$(readers)
    [ -z "$RUNNING" ] || die "作業ツリーを読む台本が走っている (HEAD を動かさない): $RUNNING"
    git -C $WT checkout --detach $FULL >> $LOG 2>&1 || die "git -C $WT checkout --detach $FULL"
    log "作業ツリー $WT の HEAD を $OLD → $FULL に動かした"
  else
    log "作業ツリー $WT は既に $FULL (変更なし)"
  fi
else
  git -C $REPO worktree add --detach $WT $FULL >> $LOG 2>&1 || die "git worktree add --detach $WT $FULL"
  log "作業ツリー $WT を作った (HEAD $FULL)"
fi
cmp -s "$0" $WT45/h370_prep.sh || die "回している $0 が commit の h370_prep.sh と違う"
for f in $SCRIPTS $INPUTS; do
  [ -f $WT45/$f ] || die "作業ツリーに $f が無い"
  if [ -e $C45/$f ]; then
    cmp -s $WT45/$f $C45/$f || die "$C45/$f が commit の版と違う (消さない。中身を確かめてから退ける)"
    log "台本・入力 $f: case dir に同じものがある ($(sha16 $C45/$f))"
  else
    mkdir -p "$(dirname $C45/$f)" && cp $WT45/$f $C45/$f || die "$f の配置"; log "台本・入力 $f: case dir に置いた ($(sha16 $C45/$f))"
  fi
done
[ "$(sha256sum $C45/_band_ab/half_inlet/delta_r_initial_C.csv | cut -d' ' -f1)" = "$CSV_SHA" ] || die "δ_r の CSV の sha256 が登録 (24251cb6…) と違う"

# --- (b) バイナリ・python ----------------------------------------------------------------------------------------------------
s=$(sha16 $BIN); log "forge $BIN sha $s"; [ "$s" = 7e95e5bbd8f50c9e ] || die "forge が登録 (7e95e5bbd8f50c9e) と違う"
s=$(sha16 $CONV); log "変換器 $CONV sha $s"; [ "${s:0:8}" = 711f4438 ] || die "変換器が登録 (711f4438) と違う"
npv=$(python3 -c "import numpy; print(numpy.__version__)" 2>&1); log "python3 $(python3 -V 2>&1) numpy $npv"
[ "$npv" = 2.5.2 ] || die "numpy が登録 (2.5.2) と違う"
mods=$(python3 -c "import h5py, yaml, scipy, matplotlib; print('h5py', h5py.__version__, 'scipy', scipy.__version__, 'matplotlib', matplotlib.__version__)" 2>&1) \
  || die "python3 で h5py・yaml・scipy・matplotlib を読めない ($mods)"
log "python3 の部品: $mods"

# --- (c) B にある入力 ------------------------------------------------------------------------------------------------------
cd $C45 || die "cd $C45"
NEED=""
for f in $R062_FILES; do [ -f $R062/$f ] || NEED="$NEED $R062/$f"; done
L174=$(ls $R174 2>/dev/null | grep -E '^res_[0-9]+\.h5$' | sort -t_ -k2 -n | tail -1)
{ [ -f $R174/prepare_info.json ] && [ -f $R174/nozzle.h5 ] && [ -n "$L174" ]; } || NEED="$NEED $R174"
grep -qF "$V5D_LINE" _band_ab/moc_v5d_run.log 2>/dev/null || NEED="$NEED _band_ab/moc_v5d_run.log"
if [ -n "$NEED" ]; then
  {
    echo "== 足りない入力:$NEED"
    echo "   run_0062・run_0174 は A にある (tw370_prep.sh が A → B で写した経路と同じ)。A で:"
    echo "     K=\"-i \$HOME/.ssh/id_a2b -o StrictHostKeyChecking=accept-new\"; B=ubuntu@$B_IP; C=forge-wallfit/$CASE; cd ~/\$C"
    echo "     L=\$(ls $R174 | grep -E '^res_[0-9]+\\.h5\$' | sort -t_ -k2 -n | tail -1); ssh \$K \$B \"mkdir -p \$C/$R174 \$C/$R062\""
    echo "     scp \$K $R174/prepare_info.json $R174/nozzle.h5 $R174/\$L \$B:\$C/$R174/ && (cd $R062 && scp \$K $R062_FILES \$B:\$C/$R062/)"
    echo "   _band_ab/moc_v5d_run.log は手元 (commit されていない) にある。手元で:"
    echo "     scp -i ~/.ssh/test.pem case/45.isobutane_m6_d155/_band_ab/moc_v5d_run.log ubuntu@<B の公開 IP>:forge-wallfit/$CASE/_band_ab/"
    echo "   写した後に同じ引数でこの台本を回し直す。"
  } | tee -a $LOG
  log "h370_prep: 入力の写しを待つ (exit 3) $(date -Is)"; exit 3
fi
for f in $R062_FILES; do log "CFD ピン $R062/$f sha $(sha16 $R062/$f)"; done
log "腕 A の旧い Euler $R174: prepare_info.json $(sha16 $R174/prepare_info.json)、nozzle.h5 $(sha16 $R174/nozzle.h5)、最後の res $L174 $(sha16 $R174/$L174)"

# --- (d) 問題 YAML の検査 (作業ツリーの h370_problems.py check) -----------------------------------------------------------------
python3 h370_run.py --phase check-problems --commit $FULL > /tmp/h370_problems_check.out 2>&1 \
  || { tail -30 /tmp/h370_problems_check.out >> $LOG; die "h370_problems.py check (作業ツリーの版) が不合格"; }
log "問題 YAML の検査: $(grep -E '^VERDICT' /tmp/h370_problems_check.out | tail -1) → _band_ab/half_inlet/h370_problems_B.json"

# --- (e) 既存のもの・資源 ----------------------------------------------------------------------------------------------------
for r in run_0950_h370_euler_ref run_0951_h370_ns_b0 _h370_prelaunch_ns _band_ab/half_inlet/h370_euler_ref_pin.json _band_ab/half_inlet/h370_euler_prequalify.json; do
  [ -e $C45/$r ] && log "注意: $r が既にある (h370_run.py は既存のものを使わずに止まる; 固定 pin は同じ中身なら使う)"
done
log "資源: ディスクの空き $(df -h --output=avail $C45 | tail -1) (起動の前に 20 GiB 以上が要る)、メモリ $(free -m | awk 'NR==2{print $7" MB available"}')、GPU $(nvidia-smi --query-gpu=name,memory.used,memory.total,utilization.gpu --format=csv,noheader 2>&1)"
nvidia-smi --query-compute-apps=pid,process_name,used_memory --format=csv,noheader >> $LOG 2>&1
for p in $(pgrep -x forge); do log "  forge $p cwd $(readlink /proc/$p/cwd)"; done
printf 'commit %s\ndate %s\n' "$FULL" "$(date -Is)" > $C45/h370_prep.done
log "== h370_prep.sh 終了 $(date -Is) → h370_prep.done (次: 回す前の確かめ nohup bash h370_run.sh --phase euler --commit $FULL --dry-run r1 > /dev/null 2>&1 < /dev/null &、"
log "   合格なら _h370_dry_r1_* を退けて nohup bash h370_run.sh --phase all --commit $FULL > /dev/null 2>&1 < /dev/null &)"
