#!/bin/bash
# plan tooling-nozzle-upstream-poly-and-throat-sizing §6 U4「NS の 3 条件」(N0・N1・N2) の投入 (AWS、case dir で実行)。
# 順序と根拠: plan discretization-moc-axis-limit-and-corrector §6 V5d (E3 → E4 → V5d → 較正値の確定 → N0 → N1 → N2)、
#   諮問 notes/reviews/2026-10-07-euler-grid-switch-plan-diagnose.md。合否は plan tooling-nozzle-throat-monotone-r2 §6 N・K の転記。
#   N0 = 今の生産 (単調壁・ramp・legacy MOC) で Md_moc_offset だけを新しい較正値に / N1 = N0 + pw_upstream poly / N2 = N1 + MOC analytic・converge。
#   較正値・k_f・r_t・NS の格子・実効の設定は 3 条件で共通 (make_ns_n012_problems.py が作った問題と ns_n012.py verify-set で照合)。
#
# usage:
#   bash run_ns_n012.sh main <MD_OFFSET> <N0_dry> <N1_dry> <N2_dry> <N0_cond> <N1_cond> <N2_cond> [NPAR (既定 3)]
#        dry 3 本: 準備 (prepare_ns + interp_field の cross-mesh 移送、IC = run_0149 の最終場) → IC の検査 → 3 条件の固定の検査
#        → 並列で段階起動 (soft 3000 → mid 3000 → 本段 2 次 cfl 1・relax 0.7・80000 step・5000 ごと) → 後処理。
#        凝縮 (COND の既定 after): dry の全ゲートに合格した条件だけ、dry の res_80000 から convert_species_field で 18000 step (1000 ごと)。
#        最後に評価器 (ns_n012_eval.py eval → _band_ab/ns_n012_eval.json)。
#   bash run_ns_n012.sh ext  <MD_OFFSET> <N0|N1|N2> <src_dry_run> <ext_run>     延長 1 回 (res_80000 から restart_field、20000 step)
#   bash run_ns_n012.sh cond <MD_OFFSET> <N0|N1|N2> <src_dry_or_ext_run> <cond_run>   凝縮だけ (src の最終 res から)
#   bash run_ns_n012.sh eval <MD_OFFSET> N0=<dry>[+<ext>] N1=... N2=... [K0=<cond> K1=... K2=...]   評価だけ
#   DRY=1 bash run_ns_n012.sh main ...   乾式確認: _dry_ns_n012/ の下に準備し、IC の検査・固定の検査・凝縮の準備 (dry の IC から、
#        convert_species_field の手前) まで。
#        forge は起動しない (FORGE_BIN を存在しない道に、FORGE_ALLOW_UNVERIFIED_SPECIES=1)。変換器は起動する (REAL_CONVERTER)。
#        DRY_IC_SRC_RUN (case dir からの相対) で IC のドナーを差し替えられる (乾式のみ。既定は run_0149)。
# 環境変数:
#   EULER_REF_LEGACY   N0・N1 の Euler 参照 (δ_E の抽出・報告; E4 で合格した run)。本番では必須
#   EULER_REF_ANALYTIC N2 の Euler 参照 (V5d の腕 M、同じ較正値の run)。本番では必須
#   COND=after|always|0  凝縮を続けて回す条件 (既定 after = dry の全ゲート合格の条件だけ。always は run_0147/0148 の投入と同じく判定を待たない)
#   WATCH_SEC (既定 5) 残差の NaN・Inf の見張りの間隔 [s] (序盤の確認: 残差が 300 行を超えた最初の見張りで <run>/EARLY_NAN_CHECK.txt)、MIN_FREE_GB (既定 6) 起動前に要る空き容量 [GB]
#   FORGE_BIN・REAL_CONVERTER (既定 ~/forge-wallfit-bin のビルド)、FORGE_CUDA_BLOCKSIZE (既定 128)
# 照合: forge の sha256 が run_0147 の RUN_PROVENANCE と同じ (上書きの抜け道は作らない)。既存の run dir があれば止める (消さない)。
# 早期停止: 残差 CSV の rms_* (rms_dq_* を除く) に NaN・Inf、または res_nan_*.h5 が出たら、その run の forge (このスクリプトの
#   子孫で cwd がその run) を止めて <run>/EARLY_STOP.txt に書く。他の run・他セッションの forge には触らない。
# **実行中にこのファイルを編集しないこと** (bash は逐次読み)。ns_n012.py・ns_n012_eval.py・make_ns_n012_problems.py も実行中に変えない。
set -euo pipefail
: "${FORGE_BIN:=$HOME/forge-wallfit-bin/solver_density_cuda/build/forge}"
: "${REAL_CONVERTER:=$HOME/forge-wallfit-bin/solver_density_cuda/build/convertGmshToForge}"
: "${FORGE_CUDA_BLOCKSIZE:=128}"; : "${WATCH_SEC:=5}"; : "${COND:=after}"; : "${MIN_FREE_GB:=6}"
export FORGE_BIN REAL_CONVERTER FORGE_CUDA_BLOCKSIZE
cd "$(dirname "$0")"
export FORGE_CONVERTER="$(pwd)/conv_tolerant.sh"
TOOLS=../../solver_density_cuda/tools
H=ns_n012.py; EV=ns_n012_eval.py; MK=make_ns_n012_problems.py
SEG=CONVERGENCE_VERDICT_segment.txt
REF_PROV=run_0147_ns_mono_final
IC_SRC=run_0149_ns_mono_final_ext; IC_RES=res_20000.h5
REF_DRY_A=run_0147_ns_mono_final; REF_DRY_B=run_0149_ns_mono_final_ext; REF_EU=run_0143_euler_wallfit_monoG1_r1
SELF=$(basename "$0")
case "$WATCH_SEC" in ''|*[!0-9]*) echo "WATCH_SEC は正の整数: $WATCH_SEC"; exit 2;; esac
case "$COND" in after|always|0) ;; *) echo "COND は after / always / 0: $COND"; exit 2;; esac
MODE=${1:-}; [ $# -ge 1 ] && shift
if [ -n "${DRY:-}" ]; then
  OUT=_dry_ns_n012; DRYF=--dry
  export FORGE_BIN=/nonexistent/forge-dry-run-guard FORGE_ALLOW_UNVERIFIED_SPECIES=1
  ICSRC=${DRY_IC_SRC_RUN:-$IC_SRC}
else
  OUT=.; DRYF=
  unset FORGE_ALLOW_UNVERIFIED_SPECIES
  ICSRC=$IC_SRC
fi

# --- 小道具 -----------------------------------------------------------------------------------------------------------
die() { echo "止める: $*"; exit 2; }
check_name() { [[ "$1" =~ ^run_[0-9]{4}_[A-Za-z0-9_.-]+$ ]] || die "run 名 '$1' が run_NNNN_<slug> でない"; }
check_new() {   # $1 = run 名, $2 = strict (乾式でも既存を拒否)。乾式の main は $OUT を作り直すので既存を問わない
  check_name "$1"
  if [ -z "${DRY:-}" ] || [ "${2:-}" = strict ]; then
    if [ -e "$OUT/$1" ]; then die "$OUT/$1 が既にある — 番号の衝突 (既存 run は消さない)"; fi
  fi
}
cond_of() { case "$1" in N0|N1|N2) ;; *) die "条件 '$1' は N0 / N1 / N2";; esac; }
euler_of() {   # 条件 → Euler 参照 (N0・N1 = legacy MOC の設計壁、N2 = analytic MOC の設計壁)
  if [ "$1" = N2 ]; then echo "${EULER_REF_ANALYTIC:-}"; else echo "${EULER_REF_LEGACY:-}"; fi
}
psha() { awk -F': ' '/^forge_sha256/{print $2; exit}' "$1/RUN_PROVENANCE.txt" 2>/dev/null || true; }
check_bin() {   # forge の sha256 が run_0147 と同じ (乾式でバイナリの無いときだけ省略)
  local want have
  want=$(psha "$REF_PROV")
  have=$( [ -f "$FORGE_BIN" ] && sha256sum "$FORGE_BIN" | awk '{print $1}' || true )
  if [ -n "${DRY:-}" ]; then echo "DRY: forge は起動しない (バイナリの照合を省略; 基準 ${want:-読めない})"; return 0; fi
  [ -n "$want" ] || die "$REF_PROV/RUN_PROVENANCE.txt の forge_sha256 を読めない"
  [ "$have" = "$want" ] || die "forge のバイナリ ($FORGE_BIN ${have:-読めない}) が $REF_PROV ($want) と違う"
  echo "forge sha256 $have (= $REF_PROV)"
}
check_disk() {   # $1 = 要る GB
  local avail
  avail=$(df --output=avail -BG . | tail -1 | tr -dc '0-9')
  [ "${avail:-0}" -ge "$1" ] || die "空き容量 ${avail} GB < $1 GB (共有ディスク。skill forge-aws-run §5)"
  echo "空き容量 ${avail} GB (要 $1 GB)"
}
check_euler() {   # $1 = 条件
  local eu; eu=$(euler_of "$1")
  if [ -n "${DRY:-}" ]; then return 0; fi
  [ -n "$eu" ] || die "条件 $1 の Euler 参照 (EULER_REF_$([ "$1" = N2 ] && echo ANALYTIC || echo LEGACY)) が未設定"
  [ -d "$eu" ] || die "Euler 参照 $eu が無い"
}
launch_record() {   # $1 = mode, 残り = 引数
  mkdir -p _band_ab
  local f="_band_ab/ns_n012_launch_$(date -u +%Y%m%dT%H%M%SZ)_$1${DRY:+_dry}.txt"
  {
    echo "date        : $(date -Is)"
    echo "mode        : $1 ${DRY:+(DRY 乾式確認)}"
    echo "args        : ${*:2}"
    for s in "$SELF" "$H" "$EV" "$MK"; do echo "script      : $s $(sha256sum "$s" | awk '{print $1}')"; done
    echo "git_head    : $(git rev-parse HEAD 2>/dev/null || echo unknown)"
    echo "git_status  : $(git status --porcelain -- "$SELF" "$H" "$EV" "$MK" ../../design ../../solver_density_cuda/tools 2>/dev/null | tr '\n' ' ' || echo unknown)"
    echo "forge_bin   : $FORGE_BIN $( [ -f "$FORGE_BIN" ] && sha256sum "$FORGE_BIN" | awk '{print $1}' || echo '(無い)')"
    echo "converter   : $REAL_CONVERTER (FORGE_CONVERTER=$FORGE_CONVERTER)"
    echo "euler_ref   : legacy ${EULER_REF_LEGACY:-未設定} / analytic ${EULER_REF_ANALYTIC:-未設定}"
    echo "env         : COND=$COND WATCH_SEC=$WATCH_SEC FORGE_CUDA_BLOCKSIZE=$FORGE_CUDA_BLOCKSIZE"
    echo "disk        : $(df -h . | tail -1)"
  } > "$f"
  cat "$f"
}
nan_csv() {   # $1 = 残差 CSV。rms_* (rms_dq_* を除く) の末尾 3000 行に nan・inf があれば真
  [ -f "$1" ] || return 1
  { head -1 "$1"; tail -n 3000 "$1"; } |
    awk -F, 'NR == 1 { for (i = 1; i <= NF; i++) { h = $i; gsub(/[ \t\r]/, "", h); if (h ~ /^rms_/ && h !~ /^rms_dq_/) col[i] = 1 }; next }
             { for (i in col) { v = tolower($i); gsub(/[ \t\r]/, "", v); if (v ~ /^[-+]?(nan|inf)/) { bad = 1; exit } } }
             END { exit (bad ? 0 : 1) }'
}
nan_in() {   # $1 = run dir。今の段の残差・済んだ段の残差 (residual_history_S*.csv) に nan・inf、または res_nan_*.h5 があれば真
  local f
  if compgen -G "$1/res_nan_*.h5" > /dev/null; then return 0; fi
  for f in "$1/residual_history.csv" "$1"/residual_history_S*.csv; do nan_csv "$f" && return 0; done
  return 1
}
early_check() {   # $1 = run dir。残差が 300 行を超えた最初の見張りで、序盤の NaN の確認を EARLY_NAN_CHECK.txt に残す (U4「序盤の NaN 確認」)
  local f n
  [ -f "$1/EARLY_NAN_CHECK.txt" ] && return 0
  for f in "$1"/residual_history_S1_soft.csv "$1/residual_history.csv"; do
    [ -f "$f" ] || continue
    n=$(( $(wc -l < "$f") - 1 ))
    if [ "$n" -ge 300 ]; then
      { echo "date : $(date -Is)"; echo "file : $(basename "$f") ($n 行)"
        if nan_in "$1"; then echo "VERDICT: NaN・Inf あり"; else echo "VERDICT: NaN・Inf なし"; fi; } > "$1/EARLY_NAN_CHECK.txt"
      return 0
    fi
  done
}
forge_desc() {   # $1 の子孫のうち実行ファイル名が forge のもの (自分が起動したものだけを辿る)
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
early_stop_note() {
  { echo "date  : $(date -Is)"; echo "reason: 残差に NaN・Inf (または res_nan_*.h5) — 早期停止 ($2)"
    echo "residual_tail:"; tail -n 3 "$1/residual_history.csv" 2>/dev/null || true
    ls "$1"/res_nan_*.h5 2>/dev/null || true; } > "$1/EARLY_STOP.txt"
}
watch_nan() {   # $1 = run dir, $2 = その run の python の pid
  local d=$1 py=$2 cwd fp
  cwd=$(cd "$d" && pwd -P)
  while kill -0 "$py" 2>/dev/null; do
    sleep "$WATCH_SEC"
    early_check "$d" || true
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
run_one() {   # $1 = run dir。終了コードは RUN_RC に書く
  local d=$1 rc=0 py w
  python3 "$H" run "$d" > "$d/run_stdout.log" 2>&1 &
  py=$!
  watch_nan "$d" "$py" &
  w=$!
  wait "$py" || rc=$?
  kill "$w" 2>/dev/null || true
  wait "$w" 2>/dev/null || true
  if nan_in "$d" && [ ! -f "$d/EARLY_STOP.txt" ]; then early_stop_note "$d" "終了後に検出"; fi
  echo "$rc" > "$d/RUN_RC"
}
run_many() {   # $1 = NPAR, 残り = run dir。並列に回して全部の終わりを待つ
  local np=$1 d; shift
  for d in "$@"; do
    while [ "$(jobs -rp | wc -l)" -ge "$np" ]; do wait -n || true; done
    run_one "$d" &
    echo "started $d (pid $!) $(date -Is)"
  done
  wait || true
}
ok_run() { [ "$(cat "$1/RUN_RC" 2>/dev/null)" = 0 ] && [ ! -f "$1/EARLY_STOP.txt" ]; }
report() {   # $1 = run, $2 = Euler 参照
  (cd ../../design && python3 -m forge_design.report.nozzle_report "../case/45.isobutane_m6_d155/$1" --euler "../case/45.isobutane_m6_d155/$2" \
     --no-pptx --wall-over-frac 5 > "../case/45.isobutane_m6_d155/$1/report_stdout.log" 2>&1) && echo "report $1 rc=0" || echo "report $1 rc=$?"
}
post_common() {   # $1 = run。NaN の全走査・残差図・判定区間の収束判定
  python3 "$H" nan-scan "$1" || echo "$1: NAN_SCAN が CLEAN でない"
  [ -f "$1/residual_history.png" ] || python3 $TOOLS/plot_residual.py "$1/residual_history.csv" -o "$1/residual_history.png" > /dev/null 2>&1 \
    || echo "$1: 残差図の生成に失敗"
  python3 $TOOLS/check_convergence.py "$1" --segment > "$1/$SEG" 2>&1 || true
  echo "$1: $(grep -m1 -E '^=== .*-> ' "$1/$SEG" || tail -1 "$1/$SEG")"
}
post_dry() {   # $1 = run (本段か延長), $2 = 条件。exitM_sampling_ab (quantities_series.csv)・報告・記録の時系列
  local r=$1 c=$2 eu prob
  eu=$(euler_of "$c"); prob=$(python3 -c "import make_ns_n012_problems as M; print(M.out_name('$c', 'dry'))")
  post_common "$r"
  report "$r" "$eu"
  EULER_REF=$eu SOLVE_JSON=c2pin_solve_recal.json FINAL_PROBLEM=$prob python3 exitM_sampling_ab.py "$r" "$IC_SRC" > "$r/exitM_sampling_stdout.log" 2>&1 \
    && echo "series $r ok" || echo "series $r rc=$?"
  python3 "$H" record-series "$r" --euler "$eu" > "$r/record_series_stdout.log" 2>&1 && echo "record $r ok" || echo "record $r rc=$?"
}
post_cond() {   # $1 = run, $2 = 条件
  post_common "$1"
  report "$1" "$(euler_of "$2")"
  python3 cond_series.py "$1" > "$1/cond_series_stdout.log" 2>&1 && echo "cond_series $1 ok" || echo "cond_series $1 rc=$?"
}
ref_records() {   # 生産 (run_0147 + 0149) の記録の時系列 (判定なし、無ければ作る。生産の Euler 参照 run_0143)
  mkdir -p _band_ab/ns_n012_ref
  for r in $REF_DRY_A $REF_DRY_B; do
    [ -f "_band_ab/ns_n012_ref/${r}_record.csv" ] && continue
    python3 "$H" record-series "$r" --euler "$REF_EU" --out "_band_ab/ns_n012_ref/${r}_record.csv" > "_band_ab/ns_n012_ref/${r}_record.log" 2>&1 \
      && echo "ref record $r ok" || echo "ref record $r rc=$?"
  done
}

# --- main ---------------------------------------------------------------------------------------------------------------
mode_main() {
  [ $# -ge 7 ] || die "main <MD_OFFSET> <N0_dry> <N1_dry> <N2_dry> <N0_cond> <N1_cond> <N2_cond> [NPAR]"
  local MD=$1 NPAR=${8:-3} c i
  local -a D=("$2" "$3" "$4") K=("$5" "$6" "$7") C=(N0 N1 N2)
  case "$NPAR" in ''|*[!0-9]*) die "NPAR は正の整数: $NPAR";; esac
  [ "$NPAR" -ge 1 ] || die "NPAR は 1 以上"
  for r in "${D[@]}" "${K[@]}"; do check_new "$r"; done
  [ "$(printf '%s\n' "${D[@]}" "${K[@]}" | sort -u | wc -l)" -eq 6 ] || die "run 名が重複している"
  python3 "$MK" --md-offset "$MD" --check || die "問題 YAML が較正値 $MD と一致しない (make_ns_n012_problems.py --md-offset $MD を先に)"
  [ -f "$ICSRC/$IC_RES" ] || [ -n "${DRY:-}" ] || die "IC のドナー $ICSRC/$IC_RES が無い"
  [ -d "$ICSRC" ] || die "IC のドナー $ICSRC が無い"
  for c in "${C[@]}"; do check_euler "$c"; done
  check_bin
  [ -n "${DRY:-}" ] || check_disk "$MIN_FREE_GB"
  launch_record main "$@"
  if [ -n "${DRY:-}" ]; then rm -rf "$OUT"; mkdir "$OUT"; fi
  # 1. dry 3 本の準備と検査 (1 本でも不成立なら forge を 1 本も起動しない)
  for i in 0 1 2; do
    python3 "$H" prep-dry "${C[$i]}" "$OUT/${D[$i]}" --md-offset "$MD" $DRYF ${DRY:+--ic-src-run "$ICSRC"} \
      || die "${D[$i]}: 準備に失敗 (forge は起動していない。作りかけの dir は消さない)"
    python3 "$H" ic-check "$OUT/${D[$i]}" $DRYF || die "${D[$i]}: IC の検査が不成立 ($OUT/${D[$i]}/IC_CHECK.json)"
  done
  python3 "$H" verify-set "$OUT/${D[0]}" "$OUT/${D[1]}" "$OUT/${D[2]}" --md-offset "$MD" $DRYF || die "3 条件の固定の検査が不成立"
  if [ -n "${DRY:-}" ]; then
    # 乾式: 凝縮の準備 (dry の移送後の IC を src に、convert_species_field まで) を通す
    for i in 0 1 2; do
      python3 "$H" prep-cond "${C[$i]}" "$OUT/${D[$i]}" "$OUT/${K[$i]}" --md-offset "$MD" --dry || die "${K[$i]}: 凝縮の準備に失敗"
    done
    echo "DRY: 準備・IC の検査・固定の検査・凝縮の準備まで (forge は起動していない) — $OUT/"
    return 0
  fi
  # 2. dry 3 本を並列に
  run_many "$NPAR" "${D[@]}"
  ref_records
  local -a OKD=()
  for i in 0 1 2; do
    if ok_run "${D[$i]}"; then post_dry "${D[$i]}" "${C[$i]}"; OKD+=("$i")
    else echo "FAIL ${D[$i]} (rc $(cat "${D[$i]}/RUN_RC" 2>/dev/null))$( [ -f "${D[$i]}/EARLY_STOP.txt" ] && echo ' 早期停止' || true)"
      post_common "${D[$i]}" || true; fi
  done
  # 3. 凝縮 (COND=after: dry の全ゲート合格の条件だけ)
  local -a KR=()
  if [ "$COND" != 0 ]; then
    for i in "${OKD[@]}"; do
      if [ "$COND" = after ] && ! python3 "$EV" gate-dry --md-offset "$MD" --cond-name "${C[$i]}" --run "${D[$i]}"; then
        echo "${C[$i]}: dry が全ゲート合格でない — 凝縮は保留 (延長: bash $SELF ext $MD ${C[$i]} ${D[$i]} <ext_run>、その後 bash $SELF cond ...)"
        continue
      fi
      check_disk 2
      python3 "$H" prep-cond "${C[$i]}" "${D[$i]}" "${K[$i]}" --md-offset "$MD" || { echo "${K[$i]}: 凝縮の準備に失敗 — この条件の凝縮は回さない"; continue; }
      KR+=("$i")
    done
    local -a KD=()
    for i in "${KR[@]}"; do KD+=("${K[$i]}"); done
    [ ${#KD[@]} -eq 0 ] || run_many "$NPAR" "${KD[@]}"
    for i in "${KR[@]}"; do
      if ok_run "${K[$i]}"; then post_cond "${K[$i]}" "${C[$i]}"; else echo "FAIL ${K[$i]}"; post_common "${K[$i]}" || true; fi
    done
  fi
  # 4. 評価
  local -a ARGS=(--set "N0=${D[0]}" --set "N1=${D[1]}" --set "N2=${D[2]}")
  for i in "${KR[@]}"; do ARGS+=(--cond "${C[$i]}=${K[$i]}"); done
  python3 "$EV" eval --md-offset "$MD" "${ARGS[@]}" > _band_ab/ns_n012_eval.log 2>&1 || echo "評価器が失敗 (rc $?) — _band_ab/ns_n012_eval.log"
  tail -60 _band_ab/ns_n012_eval.log
  echo ALLDONE
}

mode_ext() {
  [ $# -eq 4 ] || die "ext <MD_OFFSET> <N0|N1|N2> <src_dry_run> <ext_run>"
  local MD=$1 c=$2 src=$3 e=$4
  cond_of "$c"; check_name "$src"; check_new "$e" strict; check_euler "$c"; check_bin
  [ -n "${DRY:-}" ] || check_disk 1
  python3 "$MK" --md-offset "$MD" --check || die "問題 YAML が較正値 $MD と一致しない"
  launch_record ext "$@"
  python3 "$H" prep-ext "$OUT/$src" "$OUT/$e" $DRYF || die "$e: 延長の準備に失敗"
  if [ -n "${DRY:-}" ]; then echo "DRY: 延長の準備まで — $OUT/$e"; return 0; fi
  run_many 1 "$e"
  ok_run "$e" || { post_common "$e" || true; die "$e が失敗・早期停止 (rc $(cat "$e/RUN_RC" 2>/dev/null))"; }
  post_dry "$e" "$c"
  python3 "$EV" gate-dry --md-offset "$MD" --cond-name "$c" --run "$src" --ext "$e" || true
  echo "評価: bash $SELF eval $MD N0=...[+...] N1=... N2=... [K0=...]"
  echo ALLDONE
}

mode_cond() {
  [ $# -eq 4 ] || die "cond <MD_OFFSET> <N0|N1|N2> <src_dry_or_ext_run> <cond_run>"
  local MD=$1 c=$2 src=$3 k=$4
  cond_of "$c"; check_name "$src"; check_new "$k" strict; check_euler "$c"; check_bin
  [ -n "${DRY:-}" ] || check_disk 2
  python3 "$MK" --md-offset "$MD" --check || die "問題 YAML が較正値 $MD と一致しない"
  launch_record cond "$@"
  python3 "$H" prep-cond "$c" "$OUT/$src" "$OUT/$k" --md-offset "$MD" $DRYF || die "$k: 凝縮の準備に失敗"
  if [ -n "${DRY:-}" ]; then echo "DRY: 凝縮の準備まで — $OUT/$k"; return 0; fi
  run_many 1 "$k"
  ok_run "$k" || { post_common "$k" || true; die "$k が失敗・早期停止"; }
  post_cond "$k" "$c"
  echo ALLDONE
}

mode_eval() {
  [ $# -ge 4 ] || die "eval <MD_OFFSET> N0=<dry>[+<ext>] N1=... N2=... [K0=<cond> ...]"
  local MD=$1 a; shift
  local -a ARGS=()
  for a in "$@"; do
    case "$a" in N[012]=*) ARGS+=(--set "$a");; K[012]=*) ARGS+=(--cond "N${a:1}");; *) die "引数 '$a' は Nx=... か Kx=...";; esac
  done
  python3 "$EV" eval --md-offset "$MD" "${ARGS[@]}" $DRYF
}

case "$MODE" in
  main) mode_main "$@";;
  ext)  mode_ext "$@";;
  cond) mode_cond "$@";;
  eval) mode_eval "$@";;
  *) sed -n '2,30p' "$0"; exit 2;;
esac
