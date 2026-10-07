#!/bin/bash
# plan discretization-moc-axis-limit-and-corrector §5.1 #5 の 5a′・§6 V5b (別登録 2026-10-07、commit 4dd94af2): V5 の 7 本の延長の診断。
# 根拠 notes/reviews/2026-10-07-moc-v5-hold-next-step-diagnose.md。V5 の結果 (保留、_band_ab/moc_v5_euler_eval.json) は書き換えない。
#  親 (本段 18000 で終わった V5 の run) の res_18000.h5 から、同一メッシュの restart (restart_field.py、SRC とビット一致の検査つき) で
#  新しい run に継ぎ、本段と同じ設定 (2 次・cfl 2・implicitRelax 0.7)・同じバイナリ・soft 段なしで追加 36000 step (通算 54000) 回す。
#  出力は 1000 ごと (中間出力も残す)。再開は全 run 1 回。再延長はしない。
#    run_0154〜0156 ← 腕 B run_0143〜0145 / run_0157〜0159 ← 腕 M run_0150〜0152 / run_0160 ← ISEN run_0153
#    (対応の正本は moc_v5b_ext_eval.py の CHILDREN。このスクリプトは `moc_v5b_ext_eval.py pairs` で読む)
#  手順: 1. 前提: 子が無い、親の res_18000.h5 と V5 の時系列 (wallfit_series_v5.csv・_band_ab/wallfit_series_v5.json) がある、
#           forge の sha256 が run_0143 の RUN_PROVENANCE と同じ、親 7 本の RUN_PROVENANCE も同じ sha256。
#        2. 子を 7 本とも準備・検査してから起動する (1 本でも不成立なら forge を 1 本も起動しない):
#           prep-child (親の入力の複製。solverConfig は time.last.nStepOuter 18000 → 36000 だけを書き換え、親との差がそれだけで
#           出力間隔が 1000 であることを検査)
#           → restart_field.py <親>/res_18000.h5 <子>/nozzle.h5 (→ <子>/restart_field.log)
#           → verify-child (SRC とビット一致・restart 元の保存量 [roY* を含む] がすべて移ったか・移した量の数・メッシュが親と同じ・
#             設定 → <子>/V5B_EXTENSION.json)
#        3. NPAR 本ずつ並列に run (moc_v5b_ext_eval.py run = runner_axismach.run_staged stages="none")。終了コードは <子>/RUN_RC。
#           早期停止: 残差 CSV に NaN・Inf (または res_nan_*.h5) が出たら、その run の forge (このスクリプトが起動したプロセスの
#           子孫で、cwd がその run のもの) を止めて <子>/EARLY_STOP.txt に記録する。他の run は続ける。他セッションの forge には触らない。
#        4. 残差図 (無ければ plot_residual.py)、check_convergence --segment (→ <子>/CONVERGENCE_VERDICT_segment.txt)、
#           時系列 (eval_wallfit_euler.py --series: 子 7 本を同じ呼び出しで、評価座標は腕 B の r1 [run_0143]; tag v5b)、
#           判定 (moc_v5b_ext_eval.py → _band_ab/moc_v5b_ext_eval.json)。
# usage: bash run_moc_v5b_ext.sh [NPAR (同時実行本数, 既定 4)]
#        DRY=1 bash run_moc_v5b_ext.sh   (乾式確認: _dry_moc_v5b/ の下に子を作り、restart_field と検査までで止める。forge は起動しない
#                                         [--resolve-species も呼ばない]。本番の run 名の dir は作らない。乾式の子は DRY の印で run が拒否する)
#   バイナリ: 既定は AWS の ~/forge-wallfit-bin (FORGE_BIN で上書き可)。照合に上書きの抜け道は作らない (登録「同じバイナリ」)。
#   WATCH_SEC: 早期停止の見張りの間隔 [s] (既定 60)。
#   既存の子の run dir があれば番号の衝突として止める (既存 run は消さない。準備の途中で止まったときの作りかけの子も消さない)。
#   評価器 (moc_v5b_ext_eval.py) は投入前に commit し、その commit と sha256 を plan §9 に書く (登録)。起動時の sha256 と git の状態を
#   _band_ab/moc_v5b_launch.txt に残す。
#   **実行中にこのファイルを編集しないこと** (bash は逐次読み; run_case.sh の注記と同じ)。評価器・時系列の評価器も実行中に変えない。
set -euo pipefail
: "${FORGE_BIN:=$HOME/forge-wallfit-bin/solver_density_cuda/build/forge}"
: "${FORGE_CUDA_BLOCKSIZE:=128}"
: "${WATCH_SEC:=60}"
export FORGE_BIN FORGE_CUDA_BLOCKSIZE
cd "$(dirname "$0")"
NPAR=${1:-4}
case "$NPAR" in ''|*[!0-9]*) echo "NPAR は正の整数: $NPAR"; exit 2;; esac
[ "$NPAR" -ge 1 ] || { echo "NPAR は 1 以上"; exit 2; }
case "$WATCH_SEC" in ''|*[!0-9]*) echo "WATCH_SEC は正の整数: $WATCH_SEC"; exit 2;; esac
[ "$WATCH_SEC" -ge 1 ] || { echo "WATCH_SEC は 1 以上"; exit 2; }
TOOLS=../../solver_density_cuda/tools
EV=moc_v5b_ext_eval.py
SEG=$(python3 -c 'from throat_mono_judge import SEGMENT_VERDICT_FILE as s; print(s)')   # 評価器と同じファイル名
GEOM_REF=run_0143_euler_wallfit_monoG1_r1     # 評価座標と forge の sha256 の基準 (腕 B の r1)
PARENT_END=18000
CHILDREN=()
declare -A PARENT=()
while read -r c p role; do CHILDREN+=("$c"); PARENT[$c]=$p; done < <(python3 "$EV" pairs)
[ ${#CHILDREN[@]} -eq 7 ] || { echo "子の対応が 7 本でない (${#CHILDREN[@]}) — 止める"; exit 2; }
if [ -n "${DRY:-}" ]; then
  OUTDIR=_dry_moc_v5b; DRYFLAG=--dry; MODE="DRY (乾式確認)"
else
  OUTDIR=.; DRYFLAG=; MODE=本番
fi
SELF=$(basename "$0")

# 1. 前提 -------------------------------------------------------------------------------------------------------------
[ -f _band_ab/wallfit_series_v5.json ] || { echo "V5 の時系列の記録 _band_ab/wallfit_series_v5.json が無い — 止める (接続に使う)"; exit 2; }
for c in "${CHILDREN[@]}"; do
  p=${PARENT[$c]}
  if [ -e "$OUTDIR/$c" ] && [ -z "${DRY:-}" ]; then echo "$c が既にある — 番号の衝突。止める (既存 run は消さない)"; exit 2; fi
  [ -f "$p/res_$PARENT_END.h5" ] || { echo "親 $p/res_$PARENT_END.h5 が無い — 止める"; exit 2; }
  [ -f "$p/wallfit_series_v5.csv" ] || { echo "親 $p/wallfit_series_v5.csv (V5 の時系列) が無い — 止める (接続に使う)"; exit 2; }
done
psha() { awk -F': ' '/^forge_sha256/{print $2; exit}' "$1/RUN_PROVENANCE.txt" 2>/dev/null || true; }
WANT=$(psha "$GEOM_REF")
[ -n "$WANT" ] || { echo "$GEOM_REF/RUN_PROVENANCE.txt の forge_sha256 を読めない — 止める"; exit 2; }
for c in "${CHILDREN[@]}"; do
  got=$(psha "${PARENT[$c]}")
  [ "$got" = "$WANT" ] || { echo "親 ${PARENT[$c]} の forge_sha256 (${got:-読めない}) が $GEOM_REF ($WANT) と違う — 止める"; exit 2; }
done
HAVE=$( [ -f "$FORGE_BIN" ] && sha256sum "$FORGE_BIN" | awk '{print $1}' || true )
if [ "$HAVE" != "$WANT" ]; then
  if [ -n "${DRY:-}" ] && [ -z "$HAVE" ]; then
    echo "DRY: forge のバイナリ $FORGE_BIN が無い (ローカルの乾式確認) — バイナリの照合を省略 (本番では止まる)"
  else
    echo "forge のバイナリ ($FORGE_BIN ${HAVE:-読めない}) が $GEOM_REF ($WANT) と違う — 同じバイナリでない。止める"; exit 2
  fi
fi
mkdir -p _band_ab
LAUNCH=_band_ab/moc_v5b_launch${DRY:+_dry}.txt
{
  echo "date        : $(date -Is)"
  echo "mode        : $MODE"
  echo "evaluator   : $EV $(sha256sum "$EV" | awk '{print $1}')"
  echo "series_eval : eval_wallfit_euler.py $(sha256sum eval_wallfit_euler.py | awk '{print $1}')"
  echo "git_head    : $(git rev-parse HEAD 2>/dev/null || echo unknown)"
  echo "git_status  : $(git status --porcelain -- "$EV" "$SELF" eval_wallfit_euler.py moc_v5_euler_eval.py 2>/dev/null | tr '\n' ' ' || echo unknown)"
  echo "forge_bin   : $FORGE_BIN ${HAVE:-(無い)}"
  echo "forge_ref   : $GEOM_REF $WANT"
  echo "NPAR        : $NPAR  WATCH_SEC: $WATCH_SEC  FORGE_CUDA_BLOCKSIZE: $FORGE_CUDA_BLOCKSIZE"
  echo "disk        : $(df -h . | tail -1)"
} > "$LAUNCH"
cat "$LAUNCH"
if [ -n "${DRY:-}" ]; then
  rm -rf "$OUTDIR"; mkdir "$OUTDIR"
  # 乾式: forge を起動しない (restart_field の --resolve-species も)。化学種の属性は付かない (乾式の子は DRY の印で run が拒否する)
  export FORGE_BIN=/nonexistent/forge-dry-run-guard FORGE_ALLOW_UNVERIFIED_SPECIES=1
fi

# 2. 子を 7 本とも準備・検査する (forge はまだ起動しない) ------------------------------------------------------------------
for c in "${CHILDREN[@]}"; do
  p=${PARENT[$c]}; d="$OUTDIR/$c"
  python3 "$EV" prep-child "$p" "$d" $DRYFLAG || { echo "$c: 準備に失敗 — 止める (forge は起動していない)"; exit 1; }
  if ! python3 $TOOLS/restart_field.py "$p/res_$PARENT_END.h5" "$d/nozzle.h5" --dst-run "$d" > "$d/restart_field.log" 2>&1; then
    tail -5 "$d/restart_field.log"; echo "$c: restart_field が失敗 — 止める (forge は起動していない。作りかけの $d は消さない)"; exit 1
  fi
  tail -1 "$d/restart_field.log"
  python3 "$EV" verify-child "$p" "$d" $DRYFLAG || { echo "$c: 延長の検査が不成立 — 止める (forge は起動していない。$d/V5B_EXTENSION.json)"; exit 1; }
done
if [ -n "${DRY:-}" ]; then
  echo "DRY: 準備・restart・検査まで (forge は起動していない) — $OUTDIR/"
  exit 0
fi

# 3. 実行 (NPAR 本ずつ並列、残差の NaN・Inf で早期停止) ---------------------------------------------------------------
nan_in() {   # $1 = 子の dir。残差 CSV の rms_* 列 (rms_dq_* を除く; check_convergence の DIVERGED と同じ列) に nan・inf、または res_nan_*.h5 があれば真
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
early_stop_note() {   # $1 = 子の dir, $2 = いつ検出したか
  { echo "date  : $(date -Is)"; echo "reason: 残差に NaN・Inf (または res_nan_*.h5) — 早期停止 (plan §6 V5b、$2)"
    echo "residual_tail:"; tail -n 3 "$1/residual_history.csv" 2>/dev/null || true
    ls "$1"/res_nan_*.h5 2>/dev/null || true; } > "$1/EARLY_STOP.txt"
}
watch_nan() {   # $1 = 子の dir, $2 = その run の python の pid
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
run_one() {   # $1 = 子。終了コードは RUN_RC に書く (wait -n で回収済みの pid を後から wait すると状態が取れないため)
  local d=$1 rc=0 py w
  python3 "$EV" run "$d" > "$d/run_stdout.log" 2>&1 &
  py=$!
  watch_nan "$d" "$py" &
  w=$!
  wait "$py" || rc=$?
  kill "$w" 2>/dev/null || true
  wait "$w" 2>/dev/null || true
  if nan_in "$d" && [ ! -f "$d/EARLY_STOP.txt" ]; then early_stop_note "$d" "終了後に検出"; fi
  echo "$rc" > "$d/RUN_RC"
}
for c in "${CHILDREN[@]}"; do
  while [ "$(jobs -rp | wc -l)" -ge "$NPAR" ]; do wait -n || true; done
  run_one "$c" &
  echo "started $c (pid $!)"
done
wait || true
FAILED=()
for c in "${CHILDREN[@]}"; do
  rc=$(cat "$c/RUN_RC" 2>/dev/null || echo missing)
  if [ "$rc" = 0 ] && [ ! -f "$c/EARLY_STOP.txt" ]; then
    echo "done $c $(tail -1 "$c/run_stdout.log")"
  else
    echo "FAIL $c (rc $rc)$( [ -f "$c/EARLY_STOP.txt" ] && echo ' 早期停止 (EARLY_STOP.txt)' || true)"; FAILED+=("$c")
  fi
done

# 4. 後処理: 残差図・収束判定 (判定区間 = stage_manifest の最後の区間 = 子の main)、時系列、判定 ---------------------------
for c in "${CHILDREN[@]}"; do
  [ -f "$c/residual_history.csv" ] || { echo "$c: residual_history.csv が無い"; continue; }
  [ -f "$c/residual_history.png" ] || python3 $TOOLS/plot_residual.py "$c/residual_history.csv" -o "$c/residual_history.png" > /dev/null 2>&1 \
    || echo "$c: 残差図の生成に失敗"
  python3 $TOOLS/check_convergence.py "$c" --segment > "$c/$SEG" 2>&1 || true
  echo "$c: $(grep -m1 -E '^=== .*-> ' "$c/$SEG" || tail -1 "$c/$SEG")"
done
ALL7=$(IFS=,; echo "${CHILDREN[*]}")
python3 eval_wallfit_euler.py . --fixed-coef --series="$ALL7" --geom-ref="$GEOM_REF" --tag=v5b > _band_ab/moc_v5b_series.log 2>&1 \
  || echo "時系列の評価器が失敗 (rc $?) — _band_ab/moc_v5b_series.log"
tail -9 _band_ab/moc_v5b_series.log
python3 "$EV" . > _band_ab/moc_v5b_ext_eval.log 2>&1 || echo "判定の評価器が失敗 (rc $?) — _band_ab/moc_v5b_ext_eval.log"
tail -40 _band_ab/moc_v5b_ext_eval.log
[ ${#FAILED[@]} -eq 0 ] || { echo "失敗・早期停止した run: ${FAILED[*]}"; exit 1; }
echo ALLDONE
