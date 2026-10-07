#!/bin/bash
# plan discretization-moc-axis-limit-and-corrector §5.1 #5・§6 V5d (2026-10-07 登録 99431498): 新しい Euler の格子 (mesh_euler、2000 × 97・全域 0.005) での
# MOC の比較。腕 B = 単調壁・legacy + fixed2、腕 M = 単調壁・analytic + converge。両腕とも同じ出口較正の値 (E4/E4V の採用値) を Md_moc_offset に書き、
# 等エントロピー IC から各 3 本。
#   既定 (共用しない): 腕 B run_0171〜0173_euler_v5d_B_r{1,2,3}、腕 M run_0174〜0176_euler_v5d_M_r{1,2,3} の 6 本を新しく回す。
#   --share-e4v      : E4V の run_0164_euler_e4_recal_d1 を腕 B の 1 本として共用する (腕 B = run_0171・run_0172・run_0164; run_0173 は作らない)。
#                      共用の前提 (moc_v5d.py check-share・prep): run_0164 が完走・早期停止なし・判定窓の res と残差判定がある・Md_moc_offset が同じ・
#                      E4 の評価の参照 run が run_0164、腕 B の prep と問題 (name 以外)・壁・格子・設定ファイルが同じ。
# 手順 (E4 の run_e4_recal.sh と同じ流儀):
#   1. 前提: 新しく作る run dir が無い (あれば番号の衝突として止める。既存 run は消さない)、問題の生成と検査 (make-problems: 本番は E4 の評価
#      _band_ab/e4_recal_eval.json の採用値が --md-offset と同じこと)、共用なら check-share、forge の sha256 が run_0143 の RUN_PROVENANCE と同じ
#      (共用なら run_0164 とも同じ; 上書きの抜け道は作らない)、変換器の sha256 を記録、ディスクの空き (MIN_FREE_GB)。
#      起動の記録 _band_ab/moc_v5d_launch.txt・.json (評価器・スクリプトの sha256・腕の構成・git の状態)。
#   2. 準備 (moc_v5d.py prep B / prep M → cross-check): RA.prepare (等エントロピー IC) → 検査 → 腕 B と腕 M の照合。不成立なら forge を起動しない。
#   3. prep を run dir に複製し (腕ごとに同じ prep)、入力が prep のままかを確かめ (verify-prep)、NPAR 本ずつ並列に回す (moc_v5d.py run =
#      runner_axismach.run_staged stages="soft": soft 3000 → 本段 54000、1000 ごと出力。soft 段の出力は <run>/_soft_stage/ に退避)。
#      起動の順は腕 B と腕 M を交互に (B1 M1 B2 M2 B3 M3)。終了コードは <run>/RUN_RC。
#      早期停止: 残差 CSV に NaN・Inf (または res_nan_*.h5) が出たら、その run の forge (このスクリプトが起動したプロセスの子孫で cwd がその run の
#      もの) を止めて <run>/EARLY_STOP.txt に記録する。他セッションの forge には触らない。
#   4. 後処理: 新しい run ごとに残差図 (無ければ plot_residual.py)、check_convergence --segment (→ <run>/CONVERGENCE_VERDICT_segment.txt。
#      共用の run_0164 は E4 の投入スクリプトが書いたものを読むだけで書き直さない)、評価量の時系列 (eval_wallfit_euler.py --series: 6 本を
#      同じ呼び出し・同じ X_E・X_F [腕 B の r1 = run_0171] で; tag v5d)、判定 (moc_v5d_eval.py → _band_ab/moc_v5d_eval.json)。
# usage: bash run_moc_v5d.sh --md-offset <repr> [--share-e4v] [NPAR (同時実行本数, 既定 3)]
#        DRY=1 bash run_moc_v5d.sh --md-offset <repr> [--share-e4v]
#          (乾式確認: _dry_v5d/ の下に prep を作り、検査・照合・記録まで。forge は起動しない [--resolve-species も呼ばない]。本番の run 名の dir は
#           作らない。問題 YAML は case dir に書く [V5d の run dir が無い間は、値が変われば書き直す]。ローカルでは CASE_RUNS・REAL_CONVERTER を渡す)
#   バイナリ: 既定は AWS の ~/forge-wallfit-bin (FORGE_BIN / REAL_CONVERTER で上書き可)。AWS のケース dir は ~/forge-wallfit/case/45.isobutane_m6_d155。
#   WATCH_SEC: 早期停止の見張りの間隔 [s] (既定 60)。MIN_FREE_GB: 投入前に要るディスクの空き [GB] (既定 = 新しい run の本数 × 1.6 + 1)。
#   評価器 (moc_v5d_eval.py) とこのスクリプト・moc_v5d.py は投入前に commit し、その commit と sha256 を plan §9 に書く (登録)。
#   **実行中にこのファイルを編集しないこと** (bash は逐次読み)。評価器・moc_v5d.py も実行中に変えない。
set -euo pipefail
: "${FORGE_BIN:=$HOME/forge-wallfit-bin/solver_density_cuda/build/forge}"
: "${REAL_CONVERTER:=$HOME/forge-wallfit-bin/solver_density_cuda/build/convertGmshToForge}"
: "${FORGE_CUDA_BLOCKSIZE:=128}"
: "${WATCH_SEC:=60}"
export FORGE_BIN REAL_CONVERTER FORGE_CUDA_BLOCKSIZE
export FORGE_CONVERTER="$(cd "$(dirname "$0")" && pwd)/conv_tolerant.sh"   # 終了時 GPUassert の既知の罠を許容
cd "$(dirname "$0")"
MD=""; SHARE=""; NPAR=3
while [ $# -gt 0 ]; do
  case "$1" in
    --md-offset) [ $# -ge 2 ] || { echo "--md-offset に値が要る"; exit 2; }; MD=$2; shift 2 ;;
    --share-e4v) SHARE=--share-e4v; shift ;;
    -h|--help) sed -n '2,32p' "$0"; exit 0 ;;
    *) case "$1" in ''|*[!0-9]*) echo "不明な引数: $1 (usage: bash run_moc_v5d.sh --md-offset <repr> [--share-e4v] [NPAR])"; exit 2 ;; esac
       NPAR=$1; shift ;;
  esac
done
[ -n "$MD" ] || { echo "--md-offset <repr> が要る (E4/E4V の採用値。両腕に同じ値を書く)"; exit 2; }
[ "$NPAR" -ge 1 ] || { echo "NPAR は 1 以上"; exit 2; }
case "$WATCH_SEC" in ''|*[!0-9]*) echo "WATCH_SEC は正の整数: $WATCH_SEC"; exit 2;; esac
[ "$WATCH_SEC" -ge 1 ] || { echo "WATCH_SEC は 1 以上"; exit 2; }
TOOLS=../../solver_density_cuda/tools
PY=moc_v5d.py
EV=moc_v5d_eval.py
SELF=$(basename "$0")
if [ -n "${DRY:-}" ]; then OUTROOT=_dry_v5d; MODE="DRY (乾式確認)"; DRYF=--dry; else OUTROOT=.; MODE=本番; DRYF=""; fi
# 腕の構成・起動の順 (評価器と同じ定義を moc_v5d.py から読む)
PLAN_OUT=$(python3 "$PY" plan-runs $SHARE)
read -r -a ORDER <<< "$(awk '$1 == "ORDER" { $1 = ""; print }' <<< "$PLAN_OUT")"
ALL=$(awk '$1 == "ALL" { print $2 }' <<< "$PLAN_OUT")
GEOM_REF=$(awk '$1 == "GEOM_REF" { print $2 }' <<< "$PLAN_OUT")
SHARED=$(awk '$1 == "SHARED" { print $2 }' <<< "$PLAN_OUT")
declare -A PREPOF=()
while read -r tag rd arm prep; do [ "$tag" = RUN ] && PREPOF[$rd]=$prep; done <<< "$PLAN_OUT"
[ "${#ORDER[@]}" -ge 1 ] && [ -n "$ALL" ] && [ -n "$GEOM_REF" ] || { echo "腕の構成を読めない: $PLAN_OUT"; exit 2; }
SEG=$(python3 -c 'import e4_recal_eval as e; print(e.SEGMENT_VERDICT_FILE)')
REF_BIN_RUN=$(python3 -c 'import e4_recal_eval as e; print(e.REF_RUN_BIN)')   # forge の sha256 の基準 (E2・E4 と同じ run_0143)

# 1. 前提 -------------------------------------------------------------------------------------------------------------
for rd in "${ORDER[@]}"; do
  if [ -e "$rd" ]; then echo "$rd が既にある — 番号の衝突。止める (既存 run は消さない)"; exit 2; fi
done
python3 "$PY" make-problems --md-offset "$MD" $SHARE $DRYF || { echo "問題を作れない・検査が不成立 (E4 の採用値と --md-offset を確かめる) — 止める"; exit 2; }
python3 "$PY" check-problems --md-offset "$MD" $SHARE $DRYF || { echo "問題の検査が不成立 — 止める"; exit 2; }
if [ -n "$SHARE" ]; then
  python3 "$PY" check-share --md-offset "$MD" $DRYF || { echo "共用の前提 ($SHARED) が不成立 — 止める"; exit 2; }
fi
WANT=$(awk -F': ' '/^forge_sha256/{print $2; exit}' "${CASE_RUNS:-.}/$REF_BIN_RUN/RUN_PROVENANCE.txt" 2>/dev/null || true)
HAVE=$( [ -f "$FORGE_BIN" ] && sha256sum "$FORGE_BIN" | awk '{print $1}' || true )
WANT_S=""
if [ -n "$SHARE" ]; then WANT_S=$(awk -F': ' '/^forge_sha256/{print $2; exit}' "${CASE_RUNS:-.}/$SHARED/RUN_PROVENANCE.txt" 2>/dev/null || true); fi
if [ -z "$WANT" ] || [ "$WANT" != "$HAVE" ] || { [ -n "$SHARE" ] && [ "$WANT_S" != "$HAVE" ]; }; then
  if [ -n "${DRY:-}" ]; then
    echo "DRY: forge のバイナリの照合を省略 (FORGE_BIN ${HAVE:-無い} / $REF_BIN_RUN ${WANT:-記録なし}${SHARE:+ / $SHARED ${WANT_S:-記録なし}}; 本番では止まる)"
  else
    echo "forge のバイナリ ($FORGE_BIN ${HAVE:-読めない}) が $REF_BIN_RUN (${WANT:-記録なし})${SHARE:+・$SHARED (${WANT_S:-記録なし})} と違う — 同じバイナリでない。止める"; exit 2
  fi
fi
CONV_SHA=$( [ -f "$REAL_CONVERTER" ] && sha256sum "$REAL_CONVERTER" | awk '{print $1}' || true )
[ -n "$CONV_SHA" ] || { echo "変換器 $REAL_CONVERTER が無い — 止める"; exit 2; }
NNEW=${#ORDER[@]}
: "${MIN_FREE_GB:=$(awk -v n="$NNEW" 'BEGIN { printf "%.1f", n * 1.6 + 1.0 }')}"
FREE_GB=$(df -Pk . | awk 'NR == 2 { printf "%.1f", $4 / 1048576 }')
if awk -v f="$FREE_GB" -v m="$MIN_FREE_GB" 'BEGIN { exit !(f < m) }'; then
  if [ -n "${DRY:-}" ]; then echo "DRY: ディスクの空き ${FREE_GB} GB < ${MIN_FREE_GB} GB (本番では止まる)"
  else echo "ディスクの空き ${FREE_GB} GB < ${MIN_FREE_GB} GB (新しい run ${NNEW} 本 × 約 1.5 GB + 余裕) — 止める (Write failed で run が死ぬのを避ける)"; exit 2; fi
fi
mkdir -p _band_ab
LAUNCH=_band_ab/moc_v5d_launch${DRY:+_dry}.txt
{
  echo "date        : $(date -Is)"
  echo "mode        : $MODE"
  echo "plan        : plans/active/discretization-moc-axis-limit-and-corrector.md §6 V5d (登録 $(python3 -c 'import moc_v5d_eval as e; print(e.PLAN_REG_COMMIT)'))"
  echo "md_offset   : $MD  share_e4v: ${SHARE:-なし}  NPAR: $NPAR"
  echo "arms        : $ALL  (geom_ref $GEOM_REF)"
  echo "order       : ${ORDER[*]}"
  echo "steps       : $(python3 -c 'import moc_v5d_eval as e; print("soft", e.SOFT_STEPS, "main", e.MAIN_NSTEPS, "out", e.OUT_INTERVAL, "win13", e.WIN13[0], e.WIN13[-1], "tail5", e.TAIL5[0], e.TAIL5[-1])')"
  for f in "$EV" "$PY" "$SELF" eval_wallfit_euler.py e4_recal_eval.py moc_v5_euler_eval.py euler_t0_e2_eval.py euler_t0_e2.py \
           problem_d155_euler_v5d_B.yaml problem_d155_euler_v5d_M.yaml problem_d155_euler_e4_recal_d0.yaml problem_d155_euler_pin_G1_recal_mono_moc.yaml; do
    echo "input       : $f $(sha256sum "$f" | awk '{print $1}')"
  done
  echo "git_head    : $(git rev-parse HEAD 2>/dev/null || echo unknown)"
  echo "git_status  : $(git status --porcelain -- "$EV" "$PY" "$SELF" eval_wallfit_euler.py e4_recal_eval.py ../../design 2>/dev/null | tr '\n' ' ' || echo unknown)"
  echo "forge_bin   : $FORGE_BIN ${HAVE:-(無い)}"
  echo "forge_ref   : $REF_BIN_RUN ${WANT:-(記録なし)}${SHARE:+  $SHARED ${WANT_S:-(記録なし)}}"
  echo "converter   : $REAL_CONVERTER $CONV_SHA (wrapper $FORGE_CONVERTER)"
  echo "CASE_RUNS   : ${CASE_RUNS:-.}"
  echo "WATCH_SEC   : $WATCH_SEC  FORGE_CUDA_BLOCKSIZE: $FORGE_CUDA_BLOCKSIZE  MIN_FREE_GB: $MIN_FREE_GB"
  echo "disk        : $(df -h . | tail -1)"
} > "$LAUNCH"
cat "$LAUNCH"
python3 "$PY" launch-record --md-offset "$MD" $SHARE $DRYF \
  --extra "{\"forge_bin\": \"$FORGE_BIN\", \"forge_sha256\": \"${HAVE}\", \"forge_ref_sha256\": \"${WANT}\", \"converter_sha256\": \"$CONV_SHA\", \"npar\": $NPAR}"

# 2. 準備 (検査が通るまで forge を起動しない) ---------------------------------------------------------------------------------
if [ -n "${DRY:-}" ]; then mkdir -p "$OUTROOT"; rm -rf "$OUTROOT/_prep_v5d_B" "$OUTROOT/_prep_v5d_M"; else rm -rf _prep_v5d_B _prep_v5d_M; fi
for arm in B M; do
  if ! python3 "$PY" prep "$arm" --md-offset "$MD" $SHARE --out-root "$OUTROOT" $DRYF > "$OUTROOT/_prep_v5d_${arm}.log" 2>&1; then
    tail -25 "$OUTROOT/_prep_v5d_${arm}.log"
    echo "準備 ($arm) が失敗 (検査の不成立を含む) — $OUTROOT/_prep_v5d_${arm}.log・_band_ab/moc_v5d_prep_${arm}${DRY:+_dry}.json。forge は起動していない"; exit 1
  fi
  grep -v "WARNING: cannot stamp" "$OUTROOT/_prep_v5d_${arm}.log" | tail -3
done
python3 "$PY" cross-check "$OUTROOT/_prep_v5d_B" "$OUTROOT/_prep_v5d_M" --md-offset "$MD" $DRYF \
  || { echo "腕 B と腕 M の prep の照合が不成立 — 止める (forge は起動していない)"; exit 1; }
if [ -n "${DRY:-}" ]; then
  echo "DRY: 準備・検査・照合・記録まで (forge は起動していない) — $OUTROOT/_prep_v5d_{B,M}、_band_ab/moc_v5d_prep_{B,M,cross}_dry.json"
  exit 0
fi
for rd in "${ORDER[@]}"; do cp -r "${PREPOF[$rd]}" "$rd"; done
for rd in "${ORDER[@]}"; do
  python3 "$PY" verify-prep "${PREPOF[$rd]}" "$rd" || { echo "$rd の入力が prep と違う — 止める (forge は起動していない)"; exit 2; }
done

# 3. 実行 (NPAR 本ずつ並列、残差の NaN・Inf で早期停止) -------------------------------------------------------------------------
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
  { echo "date  : $(date -Is)"; echo "reason: 残差に NaN・Inf (または res_nan_*.h5) — 早期停止 (plan §6 V5d、$2)"
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
run_one() {   # $1 = run dir。python (forge の親) と見張りを起動して待ち、終了コードを <run>/RUN_RC に書く
  local rd=$1 rc=0 py w
  python3 "$PY" run "$rd" > "$rd/run_stdout.log" 2>&1 &
  py=$!
  watch_nan "$rd" "$py" &
  w=$!
  wait "$py" || rc=$?
  kill "$w" 2>/dev/null || true
  wait "$w" 2>/dev/null || true
  if nan_in "$rd" && [ ! -f "$rd/EARLY_STOP.txt" ]; then early_stop_note "$rd" "終了後に検出"; fi
  echo "$rc" > "$rd/RUN_RC"
}
for rd in "${ORDER[@]}"; do
  while [ "$(jobs -rp | wc -l)" -ge "$NPAR" ]; do wait -n || true; done
  run_one "$rd" &
  echo "started $rd (job $!)  $(date -Is)"
done
wait || true
FAILED=()
for rd in "${ORDER[@]}"; do
  rc=$(cat "$rd/RUN_RC" 2>/dev/null || echo missing)
  if [ "$rc" = 0 ] && [ ! -f "$rd/EARLY_STOP.txt" ]; then echo "done $rd $(tail -1 "$rd/run_stdout.log")"
  else echo "FAIL $rd (rc $rc)$( [ -f "$rd/EARLY_STOP.txt" ] && echo ' 早期停止 (EARLY_STOP.txt)' || true)"; FAILED+=("$rd"); fi
done

# 4. 後処理: 残差図・収束判定 (判定区間 = stage_manifest の最後の区間 = 本段)、時系列、評価 -------------------------------------------
for rd in "${ORDER[@]}"; do
  [ -f "$rd/residual_history.csv" ] || { echo "$rd: residual_history.csv が無い"; continue; }
  [ -f "$rd/residual_history.png" ] || python3 $TOOLS/plot_residual.py "$rd/residual_history.csv" -o "$rd/residual_history.png" > /dev/null 2>&1 \
    || echo "$rd: 残差図の生成に失敗"
  python3 $TOOLS/check_convergence.py "$rd" --segment > "$rd/$SEG" 2>&1 || true
  echo "$rd: $(grep -m1 -E '^=== .*-> ' "$rd/$SEG" || tail -1 "$rd/$SEG")"
done
if [ -n "$SHARE" ]; then echo "$SHARED: $(grep -m1 -E '^=== .*-> ' "$SHARED/$SEG" 2>/dev/null || echo "$SEG が無い") (E4 の投入スクリプトが書いたもの。書き直さない)"; fi
python3 eval_wallfit_euler.py . --fixed-coef --series="$ALL" --geom-ref="$GEOM_REF" --tag=v5d > _band_ab/moc_v5d_series.log 2>&1 \
  || echo "時系列の評価器が失敗 (rc $?) — _band_ab/moc_v5d_series.log"
tail -8 _band_ab/moc_v5d_series.log
python3 "$EV" . --md-offset "$MD" $SHARE > _band_ab/moc_v5d_eval.log 2>&1 || echo "判定の評価器が失敗 (rc $?) — _band_ab/moc_v5d_eval.log"
tail -30 _band_ab/moc_v5d_eval.log
[ ${#FAILED[@]} -eq 0 ] || { echo "失敗・早期停止した run: ${FAILED[*]}"; exit 1; }
echo ALLDONE
