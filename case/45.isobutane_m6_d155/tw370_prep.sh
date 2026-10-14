#!/bin/bash
# plan tooling-nozzle-isothermal-wall-chain §5.1 #32・§6「V-tw370 の事前登録」(b5df7416) と「登録の訂正」(b6aeb263) の準備 (B)。forge は回さない。
#  (a) 作業ツリー ~/forge-g0cchain (登録の「コードは B の作業ツリー ~/forge-g0cchain (この登録の後の commit)」): ~/forge-wallfit の git で
#      `git fetch origin feature/nozzle-wall-fit-and-pipeline` し、引数の commit (この台本を commit した後の commit。fetch したブランチの祖先) に合わせる。
#      無ければ `git worktree add --detach`。あって HEAD が違えば、commit されていない変更が無く、vgc_*/tw370_* の台本が走っていないときだけ
#      `git checkout --detach <commit>` で動かす (元の HEAD を記録; V-g0c-chain の 018d993d は記録に残っている)。
#      台本 (tw370_prep.sh・tw370_run.sh・tw370_run.py・tw370_judge.py・tw370_yaml_check.py) と問題 YAML 3 本を作業ツリーから case dir へ置く
#      (既にあって中身が違えば止める。消さない)。見張り・判定の部品 (m9_watch.py・cold_series.py・cold_xcheck.py・sp_judge.py と cold_xcheck の npz)
#      が作業ツリーの版と同じことを確かめる (違えば止める; tw370_run.py も同じ照合をする)。バイナリ・numpy・YAML の差の検査 (tw370_yaml_check.py)。
#  (b) B にある入力の照合: ドナー run_0754_t2_d2 (m9_watch.json REACHED 142500、最後の res が res_142500.h5、sha が到達の記録と同じ、
#      nozzle.h5 → _tr/q32_init.h5 ba859e02…)、run_0752_t2_d1/solverConfig.yaml、run_0804_sp_g0c (bcond・物性・入力 6f8fcda6d97482e3…)、
#      _band_ab/prod_confirm/prep (delta_r_initial.csv 0bc9da5b6e8f0197…、wall_physical.csv 31b4fd80818b4a0c…)、初期線の凍結源 run_0062 (V-g0c-chain で
#      B に置いたもの)。足りないものと、δ_E・δ_target の抽出 (判定の記録) が読む Euler の参照 run_0174_euler_v5d_M_r1 の 3 つのファイル
#      (prepare_info.json・nozzle.h5・最後の res) は A から写す: B → A の ssh が通る鍵があれば B から引き、無ければ「A で実行する A → B の写しの手順」を
#      出して止める (exit 3)。A から _tw370_in/ に置いた後、同じ引数でこの台本を回し直すと照合して取り込む。
#  (c) 起動の前の固定: tw370_judge.py --fix-q300 で 300 K の参照 (run_0754 res_142500) の各断面・各量の値と a_q を
#      _band_ab/cold_pair/tw370_q300.json に固定する (登録の訂正 2)。300 K の帯の感度が登録の値 (θ +3.4509/+1.6904/−0.0580 %、
#      δ_loc −12.108/−4.639/+0.134 %) と合わなければ止める。資源 (GPU のメモリ・ホストの RAM・ディスク・既存のジョブ) を記録し、
#      最後まで通ったら tw370_prep.done (commit) を置く。
# 使い方 (B、作業ツリーの commit を引数に):
#   git -C ~/forge-wallfit fetch origin feature/nozzle-wall-fit-and-pipeline
#   git -C ~/forge-wallfit show <commit>:case/45.isobutane_m6_d155/tw370_prep.sh > /tmp/tw370_prep.sh && bash /tmp/tw370_prep.sh <commit>
set -uo pipefail
COMMIT=${1:-}
REPO=~/forge-wallfit; WT=~/forge-g0cchain; CASE=case/45.isobutane_m6_d155
C45=$REPO/$CASE; WT45=$WT/$CASE; LOG=$C45/tw370_prep.log; STG=$C45/_tw370_in
BRANCH=feature/nozzle-wall-fit-and-pipeline
A_IP=172.31.5.94; A_ID=i-0b1a5e0b8dc152f00; B_IP=172.31.14.139; B_ID=i-0ba2b91ba659254c4
BIN=~/forge-fgeom7-fp64/solver_density_cuda/build/forge; CONV=~/forge-fgeom7-fp64/solver_density_cuda/build/convertGmshToForge
DON=run_0754_t2_d2; R752=run_0752_t2_d1; R804=run_0804_sp_g0c; R062=run_0062_euler_wallfit_fit_r1_ext6k; R174=run_0174_euler_v5d_M_r1
PREP=_band_ab/prod_confirm/prep
PROBLEMS="problem_d155_ns_prod_coldmesh_tw370_g0_c1.yaml problem_d155_ns_prod_coldmesh_tw370_g0c_c1.yaml problem_d155_ns_prod_coldmesh_tw370_g0c_c50.yaml"
SCRIPTS="tw370_prep.sh tw370_run.sh tw370_run.py tw370_judge.py tw370_yaml_check.py"
R062_FILES="prepare_info.json nozzle.h5 res_6000.h5 solverConfig.yaml bcondConfig.yaml"
log() { echo "$*" | tee -a $LOG; }
die() { log "中止: $* ($(date -Is))"; exit 1; }
sha16() { sha256sum "$1" 2>/dev/null | cut -c1-16; }

TOKEN=$(curl -s -X PUT http://169.254.169.254/latest/api/token -H "X-aws-ec2-metadata-token-ttl-seconds: 60")
[ "$(curl -s -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-id)" = "$B_ID" ] || { echo "B ではない。中止"; exit 1; }
[[ "$COMMIT" =~ ^[0-9a-f]{7,40}$ ]] || { echo "使い方: bash tw370_prep.sh <commit (16 進)>"; exit 1; }
[ -d $C45 ] || { echo "$C45 が無い"; exit 1; }
rm -f $C45/tw370_prep.done
log "== tw370_prep.sh 開始 $(date -Is) commit $COMMIT (台本 $0 sha $(sha16 $0))"

# --- (a) 作業ツリー ---------------------------------------------------------------------------------------------------
git -C $REPO fetch origin $BRANCH >> $LOG 2>&1 || die "git fetch origin $BRANCH"
FULL=$(git -C $REPO rev-parse --verify --quiet "$COMMIT^{commit}") || die "commit $COMMIT が $REPO に無い"
git -C $REPO merge-base --is-ancestor $FULL FETCH_HEAD || die "commit $FULL が origin/$BRANCH ($(git -C $REPO rev-parse FETCH_HEAD)) の祖先でない"
if [ -e $WT ]; then
  OLD=$(git -C $WT rev-parse HEAD 2>/dev/null) || die "$WT が git の作業ツリーでない"
  [ -z "$(git -C $WT status --porcelain --untracked-files=no)" ] || die "既存の $WT に commit されていない変更がある (消さない)"
  if [ "$OLD" != "$FULL" ]; then
    RUNNING=$(pgrep -af "vgc_run.py|vgc_judge.py|tw370_run.py|tw370_judge.py" | grep -v pgrep || true)
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
cmp -s "$0" $WT45/tw370_prep.sh || die "回している $0 が commit の tw370_prep.sh と違う"
for f in $SCRIPTS $PROBLEMS; do
  [ -f $WT45/$f ] || die "作業ツリーに $f が無い"
  if [ -e $C45/$f ]; then
    cmp -s $WT45/$f $C45/$f || die "$C45/$f が commit の版と違う (消さない。中身を確かめてから退ける)"
    log "台本 $f: case dir に同じものがある ($(sha16 $C45/$f))"
  else
    cp $WT45/$f $C45/$f || die "$f の配置"; log "台本 $f: case dir に置いた ($(sha16 $C45/$f))"
  fi
done
for f in m9_watch.py cold_series.py cold_xcheck.py sp_judge.py problem_d155_ns_prod_coldmesh_tw300.yaml problem_d155_ns_prod_coldmesh_tw300_g0c_chain.yaml \
         _band_ab/cold_pair/theta_run_0181_ns_coldmesh_ad_100000.npz _band_ab/cold_pair/theta_run_0183_ns_coldmesh_tw300_ext_100000.npz; do
  cmp -s $C45/$f $WT45/$f || { diff $WT45/$f $C45/$f 2>/dev/null | head -20 >> $LOG; die "見張り・判定の部品 $f が作業ツリーの版と違う (case dir $(sha16 $C45/$f) / 作業ツリー $(sha16 $WT45/$f))"; }
  log "見張り・判定の部品 $f: 作業ツリーと同じ ($(sha16 $C45/$f))"
done
s=$(sha16 $BIN); log "forge $BIN sha $s"; [ "$s" = 7e95e5bbd8f50c9e ] || die "forge が登録 (7e95e5bbd8f50c9e) と違う"
s=$(sha16 $CONV); log "変換器 $CONV sha $s"; [ "${s:0:8}" = 711f4438 ] || die "変換器が登録 (711f4438) と違う"
npv=$(python3 -c "import numpy; print(numpy.__version__)" 2>&1); log "python3 $(python3 -V 2>&1) numpy $npv"
[ "$npv" = 2.5.2 ] || die "numpy が登録 (2.5.2) と違う"
( cd $C45 && python3 tw370_yaml_check.py --design $WT/design --json $C45/tw370_yaml_check.json ) > /tmp/tw370_yaml_check.out 2>&1 \
  || { tail -30 /tmp/tw370_yaml_check.out >> $LOG; die "問題 YAML 3 本の設定の差・main_line・壁温の検査 (tw370_yaml_check.py)"; }
log "問題 YAML の検査: $(tail -1 /tmp/tw370_yaml_check.out)"
[ -e $C45/_band_ab/cold_pair/tw370_judge.started ] && die "_band_ab/cold_pair/tw370_judge.started が既にある (前の試みの印; 中身を確かめてから退ける)"
for r in run_0811_tw370_g0_c1 run_0812_tw370_g0c_c1 run_0813_tw370_g0c_c50; do
  [ -e $C45/$r ] && log "注意: $r が既にある (tw370_run.py は既存の run を使わずに止まる)"
done

# --- (b) 入力 ------------------------------------------------------------------------------------------------------------
cd $C45 || die "cd $C45"
python3 - >> $LOG 2>&1 <<'PY' || die "B にある入力 (ドナー・参照・壁) が欠けているか登録と違う"
import json, hashlib, re
from pathlib import Path
def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 22), b""): h.update(c)
    return h.hexdigest()
bad = []
d = Path("run_0754_t2_d2"); m = json.loads((d / "m9_watch.json").read_text()) if (d / "m9_watch.json").is_file() else {}
rs = sorted((f for f in d.glob("res_[0-9]*.h5") if re.fullmatch(r"res_\d+\.h5", f.name)), key=lambda f: int(f.stem.split("_")[1]))
last = rs[-1].name if rs else None
s = sha(d / "res_142500.h5") if (d / "res_142500.h5").is_file() else None
print(f"[tw370_prep] ドナー {d}: 見張り {m.get('status')} {m.get('reach_step')}、最後の res {last}、res_142500.h5 sha {s}、到達の記録の sha と同じ {s == m.get('reach_sha256')}")
if m.get("status") != "REACHED" or m.get("reach_step") != 142500 or last != "res_142500.h5" or s != m.get("reach_sha256"): bad.append("ドナー")
noz = d / "nozzle.h5"
print(f"[tw370_prep] ドナーの nozzle.h5 → {noz.resolve() if noz.exists() else None} sha {sha(noz)[:16] if noz.exists() else None} (G0 の格子の照合 (2) は ba859e02)")
for f in ("solverConfig.yaml", "bcondConfig.yaml", "species_meta.yaml"):
    if not (d / f).is_file(): bad.append(f"{d}/{f}")
print(f"[tw370_prep] ドナーの物性の記録 {sorted(f.name for f in d.iterdir() if f.name.startswith('resolved_species_'))} (restart_field・interp_field の照合に要る)")
for p, want in ((Path("_tr/q32_init.h5"), "ba859e02"), (Path("run_0804_sp_g0c/nozzle.h5"), "6f8fcda6d97482e3"), (Path("_g0c/q32_init.h5"), "6f8fcda6d97482e3")):
    print(f"[tw370_prep] {p}: {sha(p)[:16] if p.is_file() else '無い'} ({want} であること)")
if not ((Path("_tr/q32_init.h5").is_file() and sha(Path("_tr/q32_init.h5")).startswith("ba859e02"))
        or (noz.is_file() and sha(noz).startswith("ba859e02"))): bad.append("G0 の格子の参照 ba859e02")
if not any(p.is_file() and sha(p).startswith("6f8fcda6d97482e3") for p in (Path("run_0804_sp_g0c/nozzle.h5"), Path("_g0c/q32_init.h5"))): bad.append("G0c の格子の参照")
for f in ("run_0752_t2_d1/solverConfig.yaml", "run_0804_sp_g0c/bcondConfig.yaml", "run_0804_sp_g0c/species_meta.yaml"):
    print(f"[tw370_prep] {f}: {sha(f)[:16] if Path(f).is_file() else '無い'}")
    if not Path(f).is_file(): bad.append(f)
for f, want in (("_band_ab/prod_confirm/prep/delta_r_initial.csv", "0bc9da5b6e8f0197"), ("_band_ab/prod_confirm/prep/wall_physical.csv", "31b4fd80818b4a0c")):
    s = sha(f)[:16] if Path(f).is_file() else None
    print(f"[tw370_prep] {f}: {s} ({want} であること)")
    if s != want: bad.append(f)
if bad: raise SystemExit(f"欠け・不一致: {bad}")
PY
merge_stage() {   # A から _tw370_in/ に置かれたものを照合して取り込む (既にあって違えば止める。消さない)
  [ -d $STG ] || return 0
  local f src dst
  while IFS= read -r f; do
    src=$STG/$f; dst=$C45/$f
    if [ -e "$dst" ]; then
      cmp -s "$src" "$dst" || die "_tw370_in/$f が既にある $dst と違う (消さない)"
      rm -f "$src"; log "取り込み $f: 既にある同じもの ($(sha16 $dst))"
    else
      mkdir -p "$(dirname "$dst")" && mv "$src" "$dst" || die "_tw370_in/$f の取り込み"
      log "取り込み $f ($(sha16 $dst))"
    fi
  done < <(cd $STG && find . -type f | sed 's|^\./||' | sort)
  find $STG -depth -type d -empty -delete 2>/dev/null
}
missing() {       # 足りない入力を 1 行ずつ出す (R174 / R062)
  local r
  if [ ! -f $R174/prepare_info.json ] || [ ! -f $R174/nozzle.h5 ] || [ -z "$(ls $R174 2>/dev/null | grep -E '^res_[0-9]+\.h5$')" ]; then echo R174; fi
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
      echo "K=\"-i \$HOME/.ssh/id_a2b -o StrictHostKeyChecking=accept-new\"; B=ubuntu@$B_IP; C=forge-wallfit/$CASE; S=\$C/_tw370_in; cd ~/\$C"
      case " $NEED " in *" R174 "*)
        echo "L=\$(ls $R174 | grep -E '^res_[0-9]+\\.h5\$' | sort -t_ -k2 -n | tail -1); echo \"最後の res \$L\"; sha256sum $R174/prepare_info.json $R174/nozzle.h5 $R174/\$L"
        echo "ssh \$K \$B \"mkdir -p \$S/$R174\" && scp \$K $R174/prepare_info.json $R174/nozzle.h5 $R174/\$L \$B:\$S/$R174/";; esac
      case " $NEED " in *" R062 "*)
        echo "ssh \$K \$B \"mkdir -p \$S/$R062\" && (cd $R062 && scp \$K $R062_FILES \$B:\$S/$R062/)";; esac
      echo "== (ここまで A で実行)"
    } | tee -a $LOG
    log "tw370_prep: A からの写しを待つ (exit 3) $(date -Is)"; exit 3
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
  case " $NEED " in *" R174 "*)
    L=$($SSH "cd $AC/$R174 && ls | grep -E '^res_[0-9]+\.h5$' | sort -t_ -k2 -n | tail -1" < /dev/null)
    [ -n "$L" ] || die "A の $R174 に res が無い"
    $SSH "cd $AC/$R174 && sha256sum prepare_info.json nozzle.h5 $L && ls -la" < /dev/null >> $LOG 2>&1
    pull $R174 prepare_info.json nozzle.h5 $L;;
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
  log "入力はそろっている (写しは不要、または _tw370_in から取り込み済み)"
fi
L174=$(ls $R174 | grep -E '^res_[0-9]+\.h5$' | sort -t_ -k2 -n | tail -1)
log "== Euler の参照 $R174 (δ_E・δ_target の抽出が読む prepare_info.json・nozzle.h5・最後の res $L174)"
( sha256sum $R174/prepare_info.json $R174/nozzle.h5 $R174/$L174 ) >> $LOG 2>&1
python3 - $R174 >> $LOG 2>&1 <<'PY' || die "Euler の参照 $R174 を読めない (δ_E の抽出の _load_structured と同じ読み方)"
import sys, json, h5py, numpy as np
from pathlib import Path
d = Path(sys.argv[1]); info = json.loads((d / "prepare_info.json").read_text())
ni, nj = int(info["mesh"]["ni"]), int(info["mesh"]["nj"])
res = sorted(d.glob("res_[0-9]*.h5"), key=lambda p: int(p.stem.split("_")[1]))[-1]
with h5py.File(d / "nozzle.h5") as f: nc = f["/MESH/COORD"].shape
with h5py.File(res) as f: q = np.asarray(f["/VALUE/ro"][:]) * np.asarray(f["/VALUE/Ux"][:])
print(f"[tw370_prep] {d}: ni {ni} nj {nj} scale_m {info.get('scale_m')}、COORD {nc}、{res.name} の ρU 非有限 {int(np.count_nonzero(~np.isfinite(q)))}")
if q.size != ni * nj or not np.all(np.isfinite(q)): raise SystemExit("Euler の参照が構造格子でないか非有限")
PY

# --- (c) 起動の前の固定と記録 --------------------------------------------------------------------------------------------------
python3 tw370_judge.py --fix-q300 >> $LOG 2>&1; rc=$?
[ $rc = 0 ] || die "q₃₀₀ の固定 (tw370_judge.py --fix-q300 rc $rc: 2 = 入力・既存のファイルと違う、3 = 300 K の帯の感度が登録の値と合わない)"
log "q₃₀₀: _band_ab/cold_pair/tw370_q300.json ($(sha16 _band_ab/cold_pair/tw370_q300.json))"
log "資源: ディスクの空き $(df -h --output=avail $C45 | tail -1)、メモリ $(free -m | awk 'NR==2{print $7" MB available"}')、GPU $(nvidia-smi --query-gpu=name,memory.used,memory.total,utilization.gpu --format=csv,noheader 2>&1)"
nvidia-smi --query-compute-apps=pid,process_name,used_memory --format=csv,noheader >> $LOG 2>&1
for p in $(pgrep -x forge); do log "  forge $p cwd $(readlink /proc/$p/cwd)"; done
printf 'commit %s\ndate %s\n' "$FULL" "$(date -Is)" > $C45/tw370_prep.done
log "== tw370_prep.sh 終了 $(date -Is) → tw370_prep.done (次: --dry-prep で R0・C1・C50 を確かめ、for A in R0 C1 C50; do nohup bash tw370_run.sh --arm \$A --commit $FULL > /dev/null 2>&1 < /dev/null & sleep 5; done)"
