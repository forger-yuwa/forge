#!/bin/bash
# plan verification-case45-euler-total-enthalpy §5.1 #3b・§6 E2 (2026-10-07 登録、改訂 a2259768・f158bac2・再改訂 a59e61b5): 半径方向の配点プロファイル
# への感度。変えるのは mesh.wall_first_frac と mesh.wall_first_frac_throat の 2 鍵だけ (連動して 1 要因):
#   run_0161_euler_t0cluster_g1    (A) problem_d155_euler_t0cluster_g1.yaml    G1 のまま (1.3e-5・4.5e-6)
#   run_0162_euler_t0cluster_u5em3 (B) problem_d155_euler_t0cluster_u5em3.yaml 両方 0.005
#   他は E1 の run_0153 と同じ (壁・2000 × 97・軸方向の節点・BC・熱物性・MOC・出口較正・バイナリ・等エントロピー IC の生成手順)。
# 手順:
#   1. 前提: run dir が無い (あれば番号の衝突として止める。既存 run は消さない)、forge の sha256 が run_0143 の RUN_PROVENANCE と同じ
#      (上書きの抜け道は作らない)、変換器の sha256 を記録。起動の記録 _band_ab/euler_t0_e2_launch.txt (評価器・スクリプトの sha256・git の状態)。
#   2. 準備 (euler_t0_e2.py prep): 問題の検査 → A・B の RA.prepare → 各腕の検査・A と B の照合 (列の x と壁の r が同一 ほか)・A と run_0153 の
#      格子の照合・起動前の IC の検査 (全節点 |T0 − 1600| ≤ 1 K)・実効のメッシュの記録 (E2_MESH.json・E2_MESH_sections.csv・
#      _band_ab/euler_t0_e2_mesh.json)。1 つでも不成立なら forge を 1 本も起動しない。
#   3. prep を run dir に複製し、入力が prep のままかを確かめ (verify-prep)、2 本を並列に回す (euler_t0_e2.py run =
#      runner_axismach.run_staged stages="soft": soft 3000 → 本段 54000、1000 ごと出力。soft 段の出力は <run>/_soft_stage/ に退避)。
#      本段の step 数は評価器の定数 (euler_t0_e2_eval.MAIN_NSTEPS) を prep・run が読み、solverConfig の実効値と窓を照合する (食い違えば止める)。
#      終了コードは <run>/RUN_RC。早期停止: 残差 CSV に NaN・Inf (または res_nan_*.h5) が出たら、その run の forge (このスクリプトが起動した
#      プロセスの子孫で cwd がその run のもの) を止めて <run>/EARLY_STOP.txt に記録する。他の run は続ける。他セッションの forge には触らない。
#   4. 残差図 (無ければ plot_residual.py)、check_convergence --segment (→ <run>/CONVERGENCE_VERDICT_segment.txt)、
#      評価器 (euler_t0_e2_eval.py → _band_ab/euler_t0_e2_eval.json)。
# usage: bash run_euler_t0_e2.sh
#        DRY=1 bash run_euler_t0_e2.sh   (乾式確認: _dry_e2/ の下に prep を作り、検査と記録まで。forge は起動しない [--resolve-species も
#                                         呼ばない]。本番の run 名の dir は作らない。ローカルでは CASE_RUNS・REAL_CONVERTER を渡す)
#   バイナリ: 既定は AWS の ~/forge-wallfit-bin (FORGE_BIN / REAL_CONVERTER で上書き可)。AWS のケース dir は ~/forge-wallfit/case/45.isobutane_m6_d155。
#   WATCH_SEC: 早期停止の見張りの間隔 [s] (既定 60)。
#   評価器 (euler_t0_e2_eval.py) は投入前に commit し、その commit と sha256 を plan §9 に書く (登録)。
#   **実行中にこのファイルを編集しないこと** (bash は逐次読み; run_case.sh の注記と同じ)。評価器・euler_t0_e2.py も実行中に変えない。
set -euo pipefail
: "${FORGE_BIN:=$HOME/forge-wallfit-bin/solver_density_cuda/build/forge}"
: "${REAL_CONVERTER:=$HOME/forge-wallfit-bin/solver_density_cuda/build/convertGmshToForge}"
: "${FORGE_CUDA_BLOCKSIZE:=128}"
: "${WATCH_SEC:=60}"
export FORGE_BIN REAL_CONVERTER FORGE_CUDA_BLOCKSIZE
export FORGE_CONVERTER="$(cd "$(dirname "$0")" && pwd)/conv_tolerant.sh"   # 終了時 GPUassert の既知の罠を許容
cd "$(dirname "$0")"
case "$WATCH_SEC" in ''|*[!0-9]*) echo "WATCH_SEC は正の整数: $WATCH_SEC"; exit 2;; esac
[ "$WATCH_SEC" -ge 1 ] || { echo "WATCH_SEC は 1 以上"; exit 2; }
TOOLS=../../solver_density_cuda/tools
PY=euler_t0_e2.py
EV=euler_t0_e2_eval.py
SELF=$(basename "$0")
SEG=$(python3 -c 'from euler_t0_e2_eval import SEGMENT_VERDICT_FILE as s; print(s)')   # 評価器と同じファイル名
REF_BIN_RUN=run_0143_euler_wallfit_monoG1_r1    # forge の sha256 の基準 (評価器の REF_RUN_BIN と同じ)
declare -A RUNDIR=([A]=run_0161_euler_t0cluster_g1 [B]=run_0162_euler_t0cluster_u5em3)
declare -A PREPOF=([A]=_prep_e2_g1 [B]=_prep_e2_u5em3)
KEYS=(A B)
if [ -n "${DRY:-}" ]; then OUTROOT=_dry_e2; MODE="DRY (乾式確認)"; else OUTROOT=.; MODE=本番; fi

# 1. 前提 -------------------------------------------------------------------------------------------------------------
for k in "${KEYS[@]}"; do
  if [ -e "${RUNDIR[$k]}" ]; then echo "${RUNDIR[$k]} が既にある — 番号の衝突。止める (既存 run は消さない)"; exit 2; fi
done
python3 "$PY" check-problems || { echo "問題の検査が不成立 — 止める"; exit 2; }
WANT=$(awk -F': ' '/^forge_sha256/{print $2; exit}' "${CASE_RUNS:-.}/$REF_BIN_RUN/RUN_PROVENANCE.txt" 2>/dev/null || true)
HAVE=$( [ -f "$FORGE_BIN" ] && sha256sum "$FORGE_BIN" | awk '{print $1}' || true )
if [ -z "$WANT" ] || [ "$WANT" != "$HAVE" ]; then
  if [ -n "${DRY:-}" ]; then
    echo "DRY: forge のバイナリの照合を省略 (FORGE_BIN ${HAVE:-無い} / $REF_BIN_RUN ${WANT:-記録なし}; 本番では止まる)"
  else
    echo "forge のバイナリ ($FORGE_BIN ${HAVE:-読めない}) が $REF_BIN_RUN (${WANT:-記録なし}) と違う — 同じバイナリでない。止める"; exit 2
  fi
fi
CONV_SHA=$( [ -f "$REAL_CONVERTER" ] && sha256sum "$REAL_CONVERTER" | awk '{print $1}' || true )
[ -n "$CONV_SHA" ] || { echo "変換器 $REAL_CONVERTER が無い — 止める"; exit 2; }
mkdir -p _band_ab
LAUNCH=_band_ab/euler_t0_e2_launch${DRY:+_dry}.txt
{
  echo "date        : $(date -Is)"
  echo "mode        : $MODE"
  echo "plan        : plans/active/verification-case45-euler-total-enthalpy.md §6 E2 (登録 a59e61b5)"
  echo "steps       : $(python3 -c 'import euler_t0_e2_eval as e; print("soft", e.SOFT_STEPS, "main", e.MAIN_NSTEPS, "out", e.OUT_INTERVAL, "win13", e.WIN13[0], e.WIN13[-1], "tail5", e.TAIL5[0], e.TAIL5[-1])')"
  echo "evaluator   : $EV $(sha256sum "$EV" | awk '{print $1}')"
  echo "prep_run    : $PY $(sha256sum "$PY" | awk '{print $1}')"
  echo "launcher    : $SELF $(sha256sum "$SELF" | awk '{print $1}')"
  for f in problem_d155_euler_t0cluster_g1.yaml problem_d155_euler_t0cluster_u5em3.yaml problem_d155_euler_pin_G1_recal_mono_moc.yaml \
           euler_t0_stage_ab.py moc_v5c_thermo_ab.py throat_mono_ab.py; do
    echo "input       : $f $(sha256sum "$f" | awk '{print $1}')"
  done
  echo "git_head    : $(git rev-parse HEAD 2>/dev/null || echo unknown)"
  echo "git_status  : $(git status --porcelain -- "$EV" "$PY" "$SELF" problem_d155_euler_t0cluster_*.yaml ../../design 2>/dev/null | tr '\n' ' ' || echo unknown)"
  echo "forge_bin   : $FORGE_BIN ${HAVE:-(無い)}"
  echo "forge_ref   : $REF_BIN_RUN ${WANT:-(記録なし)}"
  echo "converter   : $REAL_CONVERTER $CONV_SHA (wrapper $FORGE_CONVERTER)"
  echo "CASE_RUNS   : ${CASE_RUNS:-.}"
  echo "WATCH_SEC   : $WATCH_SEC  FORGE_CUDA_BLOCKSIZE: $FORGE_CUDA_BLOCKSIZE"
  echo "disk        : $(df -h . | tail -1)"
} > "$LAUNCH"
cat "$LAUNCH"

# 2. 準備 (2 本とも検査が通るまで forge を起動しない) --------------------------------------------------------------------
if [ -n "${DRY:-}" ]; then rm -rf "$OUTROOT"; mkdir "$OUTROOT"; else rm -rf "${PREPOF[A]}" "${PREPOF[B]}"; fi
DRYFLAG=${DRY:+--dry}
if ! python3 "$PY" prep --out-root "$OUTROOT" $DRYFLAG > "$OUTROOT/_prep_e2.log" 2>&1; then
  tail -25 "$OUTROOT/_prep_e2.log"; echo "準備が失敗 (検査の不成立を含む) — $OUTROOT/_prep_e2.log・_band_ab/euler_t0_e2_mesh${DRY:+_dry}.json。forge は起動していない"; exit 1
fi
grep -v "WARNING: cannot stamp" "$OUTROOT/_prep_e2.log" | tail -6
if [ -n "${DRY:-}" ]; then
  echo "DRY: 準備・検査・記録まで (forge は起動していない) — $OUTROOT/、_band_ab/euler_t0_e2_mesh_dry.json"
  exit 0
fi
for k in "${KEYS[@]}"; do
  cp -r "${PREPOF[$k]}" "${RUNDIR[$k]}"
  python3 "$PY" verify-prep "${PREPOF[$k]}" "${RUNDIR[$k]}" || { echo "${RUNDIR[$k]} の入力が prep と違う — 止める (forge は起動していない)"; exit 2; }
done

# 3. 実行 (2 本並列、残差の NaN・Inf で早期停止) ---------------------------------------------------------------------------
nan_in() {   # $1 = run dir。残差 CSV の rms_* 列 (rms_dq_* を除く; check_convergence の DIVERGED と同じ列) に nan・inf、または res_nan_*.h5 があれば真
  if compgen -G "$1/res_nan_*.h5" > /dev/null; then return 0; fi
  [ -f "$1/residual_history.csv" ] || return 1
  awk -F, 'NR == 1 { for (i = 1; i <= NF; i++) { h = $i; gsub(/[ \t\r]/, "", h); if (h ~ /^rms_/ && h !~ /^rms_dq_/) col[i] = 1 }; next }
           { for (i in col) { v = tolower($i); gsub(/[ \t\r]/, "", v); if (v ~ /^[-+]?(nan|inf)/) { bad = 1; exit } } }
           END { exit (bad ? 0 : 1) }' "$1/residual_history.csv"
}
forge_desc() {   # $1 の子孫のうち実行ファイル名が forge のもの (自分が起動したものだけを辿る。pgrep -x forge で全体を探さない)
  local todo="$1" out="" p k
  while [ -n "${todo// /}" ]; do
    set -- $todo; p=$1; shift; todo="$*"
    for k in $(pgrep -P "$p" 2>/dev/null || true); do
      todo="$todo $k"
      if [ "$(cat "/proc/$k/comm" 2>/dev/null)" = forge ]; then out="$out $k"; fi
    done
  done
  echo $out
}
early_stop_note() {   # $1 = run dir, $2 = いつ検出したか
  { echo "date  : $(date -Is)"; echo "reason: 残差に NaN・Inf (または res_nan_*.h5) — 早期停止 (plan §6 E2、$2)"
    echo "residual_tail:"; tail -n 3 "$1/residual_history.csv" 2>/dev/null || true
    ls "$1"/res_nan_*.h5 2>/dev/null || true; } > "$1/EARLY_STOP.txt"
}
watch_nan() {   # $1 = run dir, $2 = その run の python の pid
  local d=$1 py=$2 cwd fp
  cwd=$(cd "$d" && pwd -P)
  while kill -0 "$py" 2>/dev/null; do
    sleep "$WATCH_SEC"
    if nan_in "$d"; then
      early_stop_note "$d" "実行中に検出"
      for fp in $(forge_desc "$py"); do
        if [ "$(readlink "/proc/$fp/cwd" 2>/dev/null)" = "$cwd" ]; then
          if kill -TERM "$fp" 2>/dev/null; then echo "killed: forge pid $fp (cwd $cwd)" >> "$d/EARLY_STOP.txt"; fi
        fi
      done
      return 0
    fi
  done
}
run_one() {   # $1 = run dir。終了コードは RUN_RC に書く (wait -n で回収済みの pid を後から wait すると状態が取れないため)
  local d=$1 rc=0 py w
  python3 "$PY" run "$d" > "$d/run_stdout.log" 2>&1 &
  py=$!
  watch_nan "$d" "$py" &
  w=$!
  wait "$py" || rc=$?
  kill "$w" 2>/dev/null || true
  wait "$w" 2>/dev/null || true
  if nan_in "$d" && [ ! -f "$d/EARLY_STOP.txt" ]; then early_stop_note "$d" "終了後に検出"; fi
  echo "$rc" > "$d/RUN_RC"
}
for k in "${KEYS[@]}"; do
  run_one "${RUNDIR[$k]}" &
  echo "started ${RUNDIR[$k]} (pid $!)"
done
wait || true
FAILED=()
for k in "${KEYS[@]}"; do
  d=${RUNDIR[$k]}
  rc=$(cat "$d/RUN_RC" 2>/dev/null || echo missing)
  if [ "$rc" = 0 ] && [ ! -f "$d/EARLY_STOP.txt" ]; then
    echo "done $d $(tail -1 "$d/run_stdout.log")"
  else
    echo "FAIL $d (rc $rc)$( [ -f "$d/EARLY_STOP.txt" ] && echo ' 早期停止 (EARLY_STOP.txt)' || true)"; FAILED+=("$d")
  fi
done

# 4. 後処理: 残差図・収束判定 (判定区間 = stage_manifest の最後の区間 = 本段)、評価 ---------------------------------------------
for k in "${KEYS[@]}"; do
  d=${RUNDIR[$k]}
  [ -f "$d/residual_history.csv" ] || { echo "$d: residual_history.csv が無い"; continue; }
  [ -f "$d/residual_history.png" ] || python3 $TOOLS/plot_residual.py "$d/residual_history.csv" -o "$d/residual_history.png" > /dev/null 2>&1 \
    || echo "$d: 残差図の生成に失敗"
  python3 $TOOLS/check_convergence.py "$d" --segment > "$d/$SEG" 2>&1 || true
  echo "$d: $(grep -m1 -E '^=== .*-> ' "$d/$SEG" || tail -1 "$d/$SEG")"
done
python3 "$EV" . > _band_ab/euler_t0_e2_eval.log 2>&1 || echo "評価器が失敗 (rc $?) — _band_ab/euler_t0_e2_eval.log"
tail -30 _band_ab/euler_t0_e2_eval.log
[ ${#FAILED[@]} -eq 0 ] || { echo "失敗・早期停止した run: ${FAILED[*]}"; exit 1; }
echo ALLDONE
