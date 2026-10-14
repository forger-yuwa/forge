#!/bin/bash
# plan tooling-nozzle-isothermal-wall-chain §5.1 #31・§6「V-g0c-chain の事前登録」(commit 68d56723) の準備 (B)。forge は回さない。
#  (a) 作業ツリー ~/forge-g0cchain: ~/forge-wallfit の git で `git fetch origin feature/nozzle-wall-fit-and-pipeline` し、引数の commit
#      (この台本を commit した後の commit。fetch したブランチの祖先であること) を `git worktree add --detach`。既にあれば HEAD と未 commit の変更なしを照合して使う。
#      台本 (vgc_prep.sh・vgc_run.sh・vgc_run.py・vgc_judge.py) と問題 YAML を作業ツリーから case dir (~/forge-wallfit/case/45…) へ置く
#      (既にあって中身が違えば止める。消さない)。見張り・判定が使う m9_watch.py・cold_series.py・cold_xcheck.py・sp_judge.py と cold_xcheck の
#      npz が作業ツリーの版と同じことを確かめる (違えば止める; ALLOW_TOOL_DIFF=1 なら差を記録して続ける)。
#  (b) 入力を A (172.31.5.94) から写し、登録の sha256 を照合する:
#      run_0179_ns_n012_N2_ext: res_20000.h5 (fb65c9b37d9bb475…) と prepare_ns の ic_from (interp_field.py) が読む隣の設定
#        (*.yaml = solverConfig・bcond・species_meta・resolved_species の記録など) と prepare_info.json・nozzle.msh (記録)。
#      _band_ab/prod_confirm/prep/ の一式 (delta_r_initial.csv 0bc9da5b6e8f0197…、wall_physical.csv 31b4fd80818b4a0c…)。
#      run_0062_euler_wallfit_fit_r1_ext6k (B に無いときだけ): 問題 YAML の geometry.initial_line_run (CFD ピンの初期線) で design_chain が読む
#        prepare_info.json・nozzle.h5・res_6000.h5・solverConfig.yaml・bcondConfig.yaml (登録の入力の一覧には無いが、prepare_ns が設計チェーンを
#        通すのに要る。登録の sha は無いので sha256 を記録し、前提の検査 (1) の .msh のビット一致で効き方を確かめる)。
#      run_0174_euler_v5d_M_r1 (euler_ref) は写さない: prepare_ns はこのパスを prepare_info.json に書くだけで、中のファイルを読まない
#        (design/forge_design/evaluate/runner_axismach.py の prepare_ns、`"euler_ref": str(euler_ref)`)。
#      写し方: B → A の ssh が通る鍵 (~/.ssh/id_b2a・id_a2b・id_ed25519・id_rsa・既定) があれば B から引く。無ければ「A で実行する A → B の
#      写しの手順」をログに出して止める (exit 3)。A から _vgc_in/ に置いた後、同じ引数でこの台本を回し直すと照合して取り込む。
#  (c) 写した後の検査の記録 (全ファイルの sha256、res_20000.h5 の量の有限性、判定の参照 T2a・T2b・run_0804 の存在、ディスク・メモリ・GPU) を
#      vgc_prep.log に書き、最後まで通ったら vgc_prep.done (commit) を置く。
# 使い方 (B、作業ツリーが無いときは git から取り出して回す):
#   git -C ~/forge-wallfit fetch origin feature/nozzle-wall-fit-and-pipeline
#   git -C ~/forge-wallfit show <commit>:case/45.isobutane_m6_d155/vgc_prep.sh > /tmp/vgc_prep.sh && bash /tmp/vgc_prep.sh <commit>
set -uo pipefail
COMMIT=${1:-}
REPO=~/forge-wallfit; WT=~/forge-g0cchain; CASE=case/45.isobutane_m6_d155
C45=$REPO/$CASE; WT45=$WT/$CASE; LOG=$C45/vgc_prep.log; STG=$C45/_vgc_in
BRANCH=feature/nozzle-wall-fit-and-pipeline
A_IP=172.31.5.94; A_ID=i-0b1a5e0b8dc152f00; B_IP=172.31.14.139; B_ID=i-0ba2b91ba659254c4
BIN=~/forge-fgeom7-fp64/solver_density_cuda/build/forge; CONV=~/forge-fgeom7-fp64/solver_density_cuda/build/convertGmshToForge
R179=run_0179_ns_n012_N2_ext; R062=run_0062_euler_wallfit_fit_r1_ext6k; R174=run_0174_euler_v5d_M_r1; PREP=_band_ab/prod_confirm/prep
PROBLEM=problem_d155_ns_prod_coldmesh_tw300_g0c_chain.yaml
R062_FILES="prepare_info.json nozzle.h5 res_6000.h5 solverConfig.yaml bcondConfig.yaml"
log() { echo "$*" | tee -a $LOG; }
die() { log "中止: $* ($(date -Is))"; exit 1; }
sha16() { sha256sum "$1" 2>/dev/null | cut -c1-16; }

TOKEN=$(curl -s -X PUT http://169.254.169.254/latest/api/token -H "X-aws-ec2-metadata-token-ttl-seconds: 60")
[ "$(curl -s -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-id)" = "$B_ID" ] || { echo "B ではない。中止"; exit 1; }
[[ "$COMMIT" =~ ^[0-9a-f]{7,40}$ ]] || { echo "使い方: bash vgc_prep.sh <commit (16 進)>"; exit 1; }
[ -d $C45 ] || { echo "$C45 が無い"; exit 1; }
rm -f $C45/vgc_prep.done
log "== vgc_prep.sh 開始 $(date -Is) commit $COMMIT (台本 $0 sha $(sha16 $0))"

# --- (a) 作業ツリー ---------------------------------------------------------------------------------------------------
git -C $REPO fetch origin $BRANCH >> $LOG 2>&1 || die "git fetch origin $BRANCH"
FULL=$(git -C $REPO rev-parse --verify --quiet "$COMMIT^{commit}") || die "commit $COMMIT が $REPO に無い"
git -C $REPO merge-base --is-ancestor $FULL FETCH_HEAD || die "commit $FULL が origin/$BRANCH ($(git -C $REPO rev-parse FETCH_HEAD)) の祖先でない"
if [ -e $WT ]; then
  [ "$(git -C $WT rev-parse HEAD 2>/dev/null)" = "$FULL" ] || die "既存の $WT の HEAD $(git -C $WT rev-parse HEAD 2>/dev/null) が $FULL でない (消さない)"
  [ -z "$(git -C $WT status --porcelain --untracked-files=no)" ] || die "既存の $WT に commit されていない変更がある"
  log "作業ツリー $WT は既にある (HEAD $FULL、変更なし)"
else
  git -C $REPO worktree add --detach $WT $FULL >> $LOG 2>&1 || die "git worktree add --detach $WT $FULL"
  log "作業ツリー $WT を作った (HEAD $FULL)"
fi
cmp -s "$0" $WT45/vgc_prep.sh || die "回している $0 が commit の vgc_prep.sh と違う"
for f in vgc_prep.sh vgc_run.sh vgc_run.py vgc_judge.py $PROBLEM; do
  [ -f $WT45/$f ] || die "作業ツリーに $f が無い"
  if [ -e $C45/$f ]; then
    cmp -s $WT45/$f $C45/$f || die "$C45/$f が commit の版と違う (消さない。中身を確かめてから退ける)"
    log "台本 $f: case dir に同じものがある ($(sha16 $C45/$f))"
  else
    cp $WT45/$f $C45/$f || die "$f の配置"; log "台本 $f: case dir に置いた ($(sha16 $C45/$f))"
  fi
done
TOOLDIFF=0
for f in m9_watch.py cold_series.py cold_xcheck.py sp_judge.py _band_ab/cold_pair/theta_run_0181_ns_coldmesh_ad_100000.npz \
         _band_ab/cold_pair/theta_run_0183_ns_coldmesh_tw300_ext_100000.npz; do
  if cmp -s $C45/$f $WT45/$f; then log "見張り・判定の部品 $f: 作業ツリーと同じ ($(sha16 $C45/$f))"
  else TOOLDIFF=1; log "見張り・判定の部品 $f: 作業ツリーと違う (case dir $(sha16 $C45/$f) / 作業ツリー $(sha16 $WT45/$f))"; diff $WT45/$f $C45/$f 2>/dev/null | head -20 >> $LOG; fi
done
[ $TOOLDIFF = 0 ] || [ "${ALLOW_TOOL_DIFF:-0}" = 1 ] || die "見張り・判定の部品が commit の版と違う (同じにするか ALLOW_TOOL_DIFF=1 で差を記録して続ける)"
s=$(sha16 $BIN); log "forge $BIN sha $s"; [ "$s" = 7e95e5bbd8f50c9e ] || die "forge が登録 (7e95e5bbd8f50c9e) と違う"
s=$(sha16 $CONV); log "変換器 $CONV sha $s"; [ "${s:0:8}" = 711f4438 ] || die "変換器が登録 (711f4438) と違う"
npv=$(python3 -c "import numpy; print(numpy.__version__)" 2>&1); log "python3 $(python3 -V 2>&1) numpy $npv"
[ "$npv" = 2.5.2 ] || die "numpy が登録 (2.5.2) と違う"
python3 - $WT/design $C45/$PROBLEM >> $LOG 2>&1 <<'PY' || die "問題 YAML の読み込み・main_line の検査 (作業ツリーの forge_design)"
import sys
sys.path.insert(0, sys.argv[1])
import forge_design
from forge_design.evaluate import runner_axismach as RA
p = RA.load_problem(sys.argv[2]); ml = RA._main_line_spec(p)
assert forge_design.__file__.startswith(sys.argv[1]), forge_design.__file__
print("[vgc_prep] forge_design", forge_design.__file__, "問題", p.name, "main_line", ml, "nStepOuter", p.evaluate["nStepOuter"], "outStepInterval", p.evaluate["outStepInterval"])
PY

# --- 判定・前提の検査が読む B の既存の run (無ければ回しても判定不能になるので止める) -------------------------------------------
cd $C45 || die "cd $C45"
[ -f run_0717_c50_d1/solverConfig.yaml ] || die "run_0717_c50_d1/solverConfig.yaml が無い (前提の検査 (5) の参照)"
python3 - >> $LOG 2>&1 <<'PY' || die "判定の参照 (T2a・T2b) か対照 run_0804 の記録が欠けている"
import json, hashlib, re
from pathlib import Path
bad = []
def sha16(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 22), b""): h.update(c)
    return h.hexdigest()[:16]
for r, need in (("run_0752_t2_d1", True), ("run_0754_t2_d2", True), ("run_0804_sp_g0c", True)):
    d = Path(r); m = d / "m9_watch.json"
    st = json.loads(m.read_text()) if m.is_file() else {}
    n = st.get("reach_step"); res = d / f"res_{n}.h5"; noz = d / "nozzle.h5"
    print(f"[vgc_prep] {r}: 見張り {st.get('status')} {n}、停止の出力 {res.name} {'あり' if res.is_file() else '無い'}、nozzle.h5 "
          f"{'読める' if noz.is_file() else '無い'}{' → ' + str(noz.resolve()) if noz.is_symlink() else ''}")
    if st.get("status") != "REACHED" or not res.is_file() or not noz.is_file():
        if r != "run_0804_sp_g0c" or st.get("status") != "REACHED":
            bad.append(r)
ctl = Path("run_0804_sp_g0c")
for f in ("bcondConfig.yaml", "species_meta.yaml"):
    if not (ctl / f).is_file(): bad.append(f"run_0804_sp_g0c/{f}")
rs = sorted(f.name for f in ctl.iterdir() if re.match(r"resolved_species_[0-9a-f]{16}(_[0-9a-f]{16})?\.yaml$", f.name))
print(f"[vgc_prep] run_0804_sp_g0c の物性の記録 {rs}")
cand = [p for p in (ctl / "nozzle.h5", Path("_g0c/q32_init.h5")) if p.is_file()]
shas = {str(p): sha16(p) for p in cand}
print(f"[vgc_prep] 格子の照合 (2) の参照の候補 {shas} (6f8fcda6d97482e3 であること)")
if "6f8fcda6d97482e3" not in shas.values(): bad.append("run_0804 の入力 (6f8fcda6d97482e3)")
if bad: raise SystemExit(f"欠け: {bad}")
PY
[ -e run_0810_vgc_chain ] && log "注意: run_0810_vgc_chain が既にある (vgc_run.py は既存の run を使わずに止まる)"

# --- (b) 入力 ------------------------------------------------------------------------------------------------------------
merge_stage() {   # A から _vgc_in/ に置かれたものを照合して取り込む (既にあって違えば止める。消さない)
  [ -d $STG ] || return 0
  local f src dst
  while IFS= read -r f; do
    src=$STG/$f; dst=$C45/$f
    if [ -e "$dst" ]; then
      cmp -s "$src" "$dst" || die "_vgc_in/$f が既にある $dst と違う (消さない)"
      rm -f "$src"; log "取り込み $f: 既にある同じもの ($(sha16 $dst))"
    else
      mkdir -p "$(dirname "$dst")" && mv "$src" "$dst" || die "_vgc_in/$f の取り込み"
      log "取り込み $f ($(sha16 $dst))"
    fi
  done < <(cd $STG && find . -type f | sed 's|^\./||' | sort)
  find $STG -depth -type d -empty -delete 2>/dev/null
}
missing() {       # 足りない入力を 1 行ずつ出す (R179 / PREP / R062)
  local r
  if [ "$(sha16 $R179/res_20000.h5)" != fb65c9b37d9bb475 ] || [ ! -f $R179/solverConfig.yaml ] || [ ! -f $R179/bcondConfig.yaml ] \
     || [ ! -f $R179/prepare_info.json ] || [ ! -f $R179/nozzle.msh ]; then echo R179; fi
  if [ "$(sha16 $PREP/delta_r_initial.csv)" != 0bc9da5b6e8f0197 ] || [ "$(sha16 $PREP/wall_physical.csv)" != 31b4fd80818b4a0c ] \
     || [ ! -f $PREP/wall_repr.json ]; then echo PREP; fi
  for r in $R062_FILES; do [ -f $R062/$r ] || { echo R062; break; }; done
}
merge_stage
NEED=$(missing | tr '\n' ' ')
if [ -n "$NEED" ]; then
  log "足りない入力: $NEED"
  AK=none
  for k in ~/.ssh/id_b2a ~/.ssh/id_a2b ~/.ssh/id_ed25519 ~/.ssh/id_rsa ""; do
    [ -z "$k" ] || [ -f "$k" ] || continue
    if ssh ${k:+-i $k} -o BatchMode=yes -o ConnectTimeout=10 -o StrictHostKeyChecking=accept-new ubuntu@$A_IP true < /dev/null > /dev/null 2>&1; then AK=$k; break; fi
  done
  if [ "$AK" = none ]; then
    {
      echo "== B → A の ssh が通らない (鍵 ~/.ssh/id_b2a・id_a2b・id_ed25519・id_rsa・既定のどれも)。A で次を実行して A → B へ写し、B で同じ引数でこの台本を回し直す:"
      echo "K=\"-i \$HOME/.ssh/id_a2b -o StrictHostKeyChecking=accept-new\"; B=ubuntu@$B_IP; C=forge-wallfit/$CASE; S=\$C/_vgc_in; cd ~/\$C"
      echo "sha256sum $R179/res_20000.h5 $PREP/delta_r_initial.csv $PREP/wall_physical.csv   # fb65c9b37d9bb475… / 0bc9da5b6e8f0197… / 31b4fd80818b4a0c… であること"
      case " $NEED " in *" R179 "*)
        echo "ssh \$K \$B \"mkdir -p \$S/$R179\" && scp \$K $R179/res_20000.h5 $R179/*.yaml $R179/prepare_info.json $R179/nozzle.msh \$B:\$S/$R179/";; esac
      case " $NEED " in *" PREP "*)
        echo "ssh \$K \$B \"mkdir -p \$S/_band_ab/prod_confirm\" && scp -r \$K $PREP \$B:\$S/_band_ab/prod_confirm/";; esac
      case " $NEED " in *" R062 "*)
        echo "ssh \$K \$B \"mkdir -p \$S/$R062\" && (cd $R062 && scp \$K $R062_FILES \$B:\$S/$R062/)";; esac
      echo "== (ここまで A で実行)"
    } | tee -a $LOG
    log "vgc_prep: A からの写しを待つ (exit 3) $(date -Is)"; exit 3
  fi
  SSH="ssh ${AK:+-i $AK} -o BatchMode=yes -o ConnectTimeout=10 ubuntu@$A_IP"; SCP="scp -q ${AK:+-i $AK} -o BatchMode=yes"
  log "B → A の ssh: 鍵 ${AK:-既定}"
  aid=$($SSH 'T=$(curl -s -X PUT http://169.254.169.254/latest/api/token -H "X-aws-ec2-metadata-token-ttl-seconds: 60"); curl -s -H "X-aws-ec2-metadata-token: $T" http://169.254.169.254/latest/meta-data/instance-id' < /dev/null)
  [ "$aid" = "$A_ID" ] || die "$A_IP の instance id ($aid) が A ($A_ID) でない"
  AC=forge-wallfit/$CASE
  pull() {   # pull <case dir からの相対の dir> <ファイル名...>
    local d=$1; shift; mkdir -p $STG/$d
    for f in "$@"; do $SCP "ubuntu@$A_IP:$AC/$d/$f" $STG/$d/$f < /dev/null || die "A の $d/$f を写せない"; done
  }
  case " $NEED " in *" R179 "*)
    $SSH "cd $AC/$R179 && sha256sum res_20000.h5 && ls res_*.h5" < /dev/null >> $LOG 2>&1
    pull $R179 res_20000.h5 prepare_info.json nozzle.msh $($SSH "cd $AC/$R179 && ls *.yaml" < /dev/null);;
  esac
  case " $NEED " in *" PREP "*)
    $SSH "cd $AC/$PREP && ls -la && sha256sum *" < /dev/null >> $LOG 2>&1
    pull $PREP $($SSH "cd $AC/$PREP && find . -maxdepth 1 -type f -printf '%f\n'" < /dev/null);;
  esac
  case " $NEED " in *" R062 "*)
    $SSH "cd $AC/$R062 && sha256sum $R062_FILES" < /dev/null >> $LOG 2>&1
    pull $R062 $R062_FILES;;
  esac
  merge_stage
  NEED=$(missing | tr '\n' ' ')
  [ -z "$NEED" ] || die "写した後も足りない: $NEED"
  log "写し: B → A の ssh (鍵 ${AK:-既定}) で B から引いた"
else
  log "入力はそろっている (写しは不要、または _vgc_in から取り込み済み)"
fi

# --- (c) 写した後の検査の記録 -----------------------------------------------------------------------------------------------
latest=$(ls $R179 | grep -E '^res_[0-9]+\.h5$' | sort -t_ -k2 -n | tail -1)
[ "$latest" = res_20000.h5 ] || die "$R179 の最後の res が $latest (res_20000.h5 であること; prepare_ns は最新の res を移す)"
log "== 入力の sha256"
( cd $C45 && sha256sum $R179/res_20000.h5 $R179/*.yaml $R179/prepare_info.json $R179/nozzle.msh $PREP/* \
     $(for f in $R062_FILES; do echo $R062/$f; done) ) >> $LOG 2>&1
[ -d $R174 ] && log "euler_ref $R174: B にある (prepare_ns は読まない)" || log "euler_ref $R174: B に無い (prepare_ns はパスを prepare_info.json に書くだけで読まないので写さない)"
python3 - $R179/res_20000.h5 >> $LOG 2>&1 <<'PY' || die "res_20000.h5 の検査"
import sys, h5py, numpy as np
with h5py.File(sys.argv[1], "r") as h:
    V = h["VALUE"]; nf = {k: int(np.count_nonzero(~np.isfinite(np.asarray(V[k][()])))) for k in V if V[k].dtype.kind == "f"}
    ro = np.asarray(V["ro"][()]); n = ro.size
    print(f"[vgc_prep] {sys.argv[1]}: 節点 {n}、量 {sorted(V)}、非有限 {sum(nf.values())}、ro min {ro.min():.4g}、属性 {dict(h.attrs)}")
    if sum(nf.values()) or ro.min() <= 0: raise SystemExit("非有限か ro ≤ 0")
PY
log "ディスクの空き $(df -h --output=avail $C45 | tail -1)、メモリ $(free -m | awk 'NR==2{print $7" MB available"}')、GPU $(nvidia-smi --query-gpu=name,memory.used,utilization.gpu --format=csv,noheader 2>&1)、forge $(pgrep -x forge | wc -l) 本"
for p in $(pgrep -x forge); do log "  forge $p cwd $(readlink /proc/$p/cwd)"; done
printf 'commit %s\ndate %s\n' "$FULL" "$(date -Is)" > $C45/vgc_prep.done
log "== vgc_prep.sh 終了 $(date -Is) → vgc_prep.done (次: nohup bash vgc_run.sh --commit $FULL > /dev/null 2>&1 < /dev/null &)"
