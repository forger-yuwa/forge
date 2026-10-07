#!/bin/bash
# plan verification-case45-euler-total-enthalpy §5.1 #5b・§6 E4 (2026-10-07 登録 99431498): 出口較正のやり直し。1 回の起動で 1 段 (1 本) を回す。
#   段 1: bash run_e4_recal.sh d0          → run_0163_euler_e4_recal_d0 (problem_d155_euler_e4_recal_d0.yaml、Md_moc_offset δ₀ = +3.770e-4)
#   段 2: bash run_e4_recal.sh d1 <δ₁>     → run_0164_euler_e4_recal_d1 (段 1 の評価が「更新」で、<δ₁> が評価器の出した δ₁ の repr と同じときだけ。
#                                            e4_recal.py make-d1 が d0 から name と Md_moc_offset だけを変えた problem_d155_euler_e4_recal_d1.yaml を作る)
#   段 2 は主セッションが段 1 の判定 (_band_ab/e4_recal_eval.json) を見てから起動する。補正は 1 回だけ (登録)。
# 固定: 単調壁 [0, 1.5]・legacy + fixed2 (MOC のキー無し)・凍結の初期線 run_0062・r_t・BC・熱物性・バイナリ。格子は mesh_euler (2000 × 97・全域 0.005)。
# 手順 (E2 の run_euler_t0_e2.sh と同じ流儀):
#   1. 前提: run dir が無い (あれば番号の衝突として止める。既存 run は消さない)、問題の検査、forge の sha256 が run_0143 の RUN_PROVENANCE と同じ
#      (上書きの抜け道は作らない)、変換器の sha256 を記録。起動の記録 _band_ab/e4_recal_launch_<段>.txt (評価器・スクリプトの sha256・git の状態)。
#   2. 準備 (e4_recal.py prep): RA.prepare (等エントロピー IC) → 検査 (step 数の同期・MOC・mono_r2・格子の採用元と実効値・Md_moc_offset・
#      凍結の初期線・r_t・壁の証拠・メッシュ品質・起動前の IC |T0 − 1600| ≤ 1 K・段 2 は段 1 の prep との同一性) → 記録。不成立なら forge を起動しない。
#   3. prep を run dir に複製し、入力が prep のままかを確かめ (verify-prep)、回す (e4_recal.py run = runner_axismach.run_staged stages="soft":
#      soft 3000 → 本段 54000、1000 ごと出力。soft 段の出力は <run>/_soft_stage/ に退避)。終了コードは <run>/RUN_RC。
#      早期停止: 残差 CSV に NaN・Inf (または res_nan_*.h5) が出たら、その run の forge (このスクリプトが起動したプロセスの子孫で cwd がその run の
#      もの) を止めて <run>/EARLY_STOP.txt に記録する。他セッションの forge には触らない。
#   4. 残差図 (無ければ plot_residual.py)、check_convergence --segment (→ <run>/CONVERGENCE_VERDICT_segment.txt)、
#      評価器 (e4_recal_eval.py → _band_ab/e4_recal_eval.json; 段 1 と、あれば段 2 を評価する)。
# usage: bash run_e4_recal.sh d0 | d1 <δ₁>
#        DRY=1 bash run_e4_recal.sh d0     (乾式確認: _dry_e4/ の下に prep を作り、検査と記録まで。forge は起動しない [--resolve-species も
#                                           呼ばない]。本番の run 名の dir は作らない。ローカルでは CASE_RUNS・REAL_CONVERTER を渡す)
#   バイナリ: 既定は AWS の ~/forge-wallfit-bin (FORGE_BIN / REAL_CONVERTER で上書き可)。AWS のケース dir は ~/forge-wallfit/case/45.isobutane_m6_d155。
#   WATCH_SEC: 早期停止の見張りの間隔 [s] (既定 60)。
#   評価器 (e4_recal_eval.py) は投入前に commit し、その commit と sha256 を plan §9 に書く (登録)。
#   **実行中にこのファイルを編集しないこと** (bash は逐次読み)。評価器・e4_recal.py も実行中に変えない。
set -euo pipefail
: "${FORGE_BIN:=$HOME/forge-wallfit-bin/solver_density_cuda/build/forge}"
: "${REAL_CONVERTER:=$HOME/forge-wallfit-bin/solver_density_cuda/build/convertGmshToForge}"
: "${FORGE_CUDA_BLOCKSIZE:=128}"
: "${WATCH_SEC:=60}"
export FORGE_BIN REAL_CONVERTER FORGE_CUDA_BLOCKSIZE
export FORGE_CONVERTER="$(cd "$(dirname "$0")" && pwd)/conv_tolerant.sh"   # 終了時 GPUassert の既知の罠を許容
cd "$(dirname "$0")"
STAGE=${1:-}
DELTA=${2:-}
case "$STAGE" in
  d0) [ -z "$DELTA" ] || { echo "d0 は δ を取らない (δ₀ は登録の +3.770e-4)"; exit 2; } ;;
  d1) [ -n "$DELTA" ] || { echo "d1 は δ₁ (段 1 の評価器が出した repr) が必要: bash run_e4_recal.sh d1 <δ₁>"; exit 2; } ;;
  *) echo "usage: bash run_e4_recal.sh d0 | d1 <δ₁>"; exit 2 ;;
esac
case "$WATCH_SEC" in ''|*[!0-9]*) echo "WATCH_SEC は正の整数: $WATCH_SEC"; exit 2;; esac
[ "$WATCH_SEC" -ge 1 ] || { echo "WATCH_SEC は 1 以上"; exit 2; }
TOOLS=../../solver_density_cuda/tools
PY=e4_recal.py
EV=e4_recal_eval.py
SELF=$(basename "$0")
read -r RUN PREP PROB SEG < <(python3 -c "import e4_recal_eval as e, e4_recal as r; s='$STAGE'; print(e.RUNS[s], r.PREPS[s], e.PROBLEMS[s], e.SEGMENT_VERDICT_FILE)")
REF_BIN_RUN=$(python3 -c 'import e4_recal_eval as e; print(e.REF_RUN_BIN)')   # forge の sha256 の基準 (評価器と同じ)
if [ -n "${DRY:-}" ]; then OUTROOT=_dry_e4; MODE="DRY (乾式確認)"; else OUTROOT=.; MODE=本番; fi

# 1. 前提 -------------------------------------------------------------------------------------------------------------
if [ -e "$RUN" ]; then echo "$RUN が既にある — 番号の衝突。止める (既存 run は消さない)"; exit 2; fi
if [ "$STAGE" = d1 ]; then
  python3 "$PY" make-d1 --delta "$DELTA" || { echo "段 2 の問題を作れない (段 1 の判定・δ₁ を確かめる) — 止める"; exit 2; }
  python3 "$PY" check-problem d1 --delta "$DELTA" || { echo "問題の検査が不成立 — 止める"; exit 2; }
else
  python3 "$PY" check-problem d0 || { echo "問題の検査が不成立 — 止める"; exit 2; }
fi
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
LAUNCH=_band_ab/e4_recal_launch_${STAGE}${DRY:+_dry}.txt
{
  echo "date        : $(date -Is)"
  echo "mode        : $MODE"
  echo "plan        : plans/active/verification-case45-euler-total-enthalpy.md §6 E4 (登録 $(python3 -c 'import e4_recal_eval as e; print(e.PLAN_REG_COMMIT)'))"
  echo "stage       : $STAGE  run $RUN  problem $PROB  delta ${DELTA:-δ₀ (登録 +3.770e-4)}"
  echo "steps       : $(python3 -c 'import e4_recal_eval as e; print("soft", e.SOFT_STEPS, "main", e.MAIN_NSTEPS, "out", e.OUT_INTERVAL, "win13", e.WIN13[0], e.WIN13[-1], "tail5", e.TAIL5[0], e.TAIL5[-1], "prev5", e.PREV5[0], e.PREV5[-1])')"
  echo "evaluator   : $EV $(sha256sum "$EV" | awk '{print $1}')"
  echo "prep_run    : $PY $(sha256sum "$PY" | awk '{print $1}')"
  echo "launcher    : $SELF $(sha256sum "$SELF" | awk '{print $1}')"
  for f in "$PROB" problem_d155_euler_pin_G1_recal_mono.yaml problem_d155_ns_finemesh_recal_final_mono.yaml \
           euler_t0_e2.py euler_t0_e2_eval.py euler_t0_stage_ab.py moc_v5c_thermo_ab.py throat_mono_ab.py; do
    echo "input       : $f $(sha256sum "$f" | awk '{print $1}')"
  done
  echo "git_head    : $(git rev-parse HEAD 2>/dev/null || echo unknown)"
  echo "git_status  : $(git status --porcelain -- "$EV" "$PY" "$SELF" problem_d155_euler_e4_recal_*.yaml ../../design 2>/dev/null | tr '\n' ' ' || echo unknown)"
  echo "forge_bin   : $FORGE_BIN ${HAVE:-(無い)}"
  echo "forge_ref   : $REF_BIN_RUN ${WANT:-(記録なし)}"
  echo "converter   : $REAL_CONVERTER $CONV_SHA (wrapper $FORGE_CONVERTER)"
  echo "CASE_RUNS   : ${CASE_RUNS:-.}"
  echo "WATCH_SEC   : $WATCH_SEC  FORGE_CUDA_BLOCKSIZE: $FORGE_CUDA_BLOCKSIZE"
  echo "disk        : $(df -h . | tail -1)"
} > "$LAUNCH"
cat "$LAUNCH"

# 2. 準備 (検査が通るまで forge を起動しない) ---------------------------------------------------------------------------------
if [ -n "${DRY:-}" ]; then mkdir -p "$OUTROOT"; rm -rf "$OUTROOT/$PREP"; else rm -rf "$PREP"; fi
DRYFLAG=${DRY:+--dry}
if ! python3 "$PY" prep "$STAGE" --out-root "$OUTROOT" $DRYFLAG > "$OUTROOT/_prep_e4_${STAGE}.log" 2>&1; then
  tail -25 "$OUTROOT/_prep_e4_${STAGE}.log"; echo "準備が失敗 (検査の不成立を含む) — $OUTROOT/_prep_e4_${STAGE}.log・_band_ab/e4_recal_prep_${STAGE}${DRY:+_dry}.json。forge は起動していない"; exit 1
fi
grep -v "WARNING: cannot stamp" "$OUTROOT/_prep_e4_${STAGE}.log" | tail -4
if [ -n "${DRY:-}" ]; then
  echo "DRY: 準備・検査・記録まで (forge は起動していない) — $OUTROOT/$PREP、_band_ab/e4_recal_prep_${STAGE}_dry.json"
  exit 0
fi
cp -r "$PREP" "$RUN"
python3 "$PY" verify-prep "$PREP" "$RUN" || { echo "$RUN の入力が prep と違う — 止める (forge は起動していない)"; exit 2; }

# 3. 実行 (残差の NaN・Inf で早期停止) -----------------------------------------------------------------------------------------
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
  { echo "date  : $(date -Is)"; echo "reason: 残差に NaN・Inf (または res_nan_*.h5) — 早期停止 (plan §6 E4、$2)"
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
rc=0
python3 "$PY" run "$RUN" > "$RUN/run_stdout.log" 2>&1 &
PYPID=$!
echo "started $RUN (pid $PYPID)"
watch_nan "$RUN" "$PYPID" &
W=$!
wait "$PYPID" || rc=$?
kill "$W" 2>/dev/null || true
wait "$W" 2>/dev/null || true
if nan_in "$RUN" && [ ! -f "$RUN/EARLY_STOP.txt" ]; then early_stop_note "$RUN" "終了後に検出"; fi
echo "$rc" > "$RUN/RUN_RC"
FAILED=0
if [ "$rc" = 0 ] && [ ! -f "$RUN/EARLY_STOP.txt" ]; then
  echo "done $RUN $(tail -1 "$RUN/run_stdout.log")"
else
  echo "FAIL $RUN (rc $rc)$( [ -f "$RUN/EARLY_STOP.txt" ] && echo ' 早期停止 (EARLY_STOP.txt)' || true)"; FAILED=1
fi

# 4. 後処理: 残差図・収束判定 (判定区間 = stage_manifest の最後の区間 = 本段)、評価 ---------------------------------------------
if [ -f "$RUN/residual_history.csv" ]; then
  [ -f "$RUN/residual_history.png" ] || python3 $TOOLS/plot_residual.py "$RUN/residual_history.csv" -o "$RUN/residual_history.png" > /dev/null 2>&1 \
    || echo "$RUN: 残差図の生成に失敗"
  python3 $TOOLS/check_convergence.py "$RUN" --segment > "$RUN/$SEG" 2>&1 || true
  echo "$RUN: $(grep -m1 -E '^=== .*-> ' "$RUN/$SEG" || tail -1 "$RUN/$SEG")"
else
  echo "$RUN: residual_history.csv が無い"
fi
python3 "$EV" . > _band_ab/e4_recal_eval_${STAGE}.log 2>&1 || echo "評価器が失敗 (rc $?) — _band_ab/e4_recal_eval_${STAGE}.log"
tail -30 _band_ab/e4_recal_eval_${STAGE}.log
[ "$FAILED" = 0 ] || { echo "失敗・早期停止: $RUN"; exit 1; }
echo ALLDONE
