#!/bin/bash
# plan discretization-moc-axis-limit-and-corrector §5.1 #5 の 5a・§6 V5 (2026-10-07): MOC の軸処理を変えた壁 (腕 M) の Euler の確認。
# 生産 Euler 格子 G1。腕 B (現行の単調壁、run_0143〜0145) は既存の run を使う (このスクリプトでは回さない)。
#  M    : run_0150〜0152 euler_wallfit_mocG1_r{1,2,3}  (problem_d155_euler_pin_G1_recal_mono_moc.yaml。IC = run_0114 最終場 (res_6000) を、
#         IC の検査 (moc_v5_euler.py の C1〜C5、記録 IC_INSPECTION.json) が成立したときだけ検証付き番号写像で写す。3 本は同じ prep の複製)
#  ISEN : run_0153 euler_icdep_mocG1_isen               (同じ問題・同じ格子 [腕 M の prep を IC を入れる前に複製]、等エントロピー IC。IC 依存の確認)
#  段: soft (1 次 cfl 0.5、3000 step) → 本段 2 次 cfl 2・implicitRelax 0.7・18000 step・1000 ごと出力 (腕 B と同じ; forge は run_case.sh 経由)。
#  終了後: residual_history.png の確認 (無ければ plot_residual.py)、check_convergence --segment (→ 各 run の CONVERGENCE_VERDICT_segment.txt)、
#  評価量の時系列 (eval_wallfit_euler.py --series: 腕 B・M・ISEN の 7 本を同じ呼び出し・同じ X_E・X_F [腕 B の r1] で; tag v5)、
#  判定 (moc_v5_euler_eval.py → _band_ab/moc_v5_euler_eval.json)。
# usage: bash run_moc_v5_euler.sh [NPAR (同時実行本数, 既定 4)]
#        DRY=1 bash run_moc_v5_euler.sh   (乾式確認: prep --dry だけで止める。forge を起動しない。run dir は作らない。
#                                          ローカルでは CASE_RUNS・REAL_CONVERTER・ARMB_REF を渡す)
#   バイナリ: 既定は AWS の ~/forge-wallfit-bin (FORGE_BIN / REAL_CONVERTER で上書き可)。AWS のケース dir は ~/forge-wallfit/case/45.isobutane_m6_d155。
#   ARMB_REF: IC の検査 C4 で照合する腕 B の run (既定 run_0143_euler_wallfit_monoG1_r1)。
#   既存の run dir があれば番号の衝突として止める (既存 run は消さない)。
#   **実行中にこのファイルを編集しないこと** (bash は逐次読み; run_case.sh の注記と同じ)。
set -euo pipefail
: "${FORGE_BIN:=$HOME/forge-wallfit-bin/solver_density_cuda/build/forge}"
: "${REAL_CONVERTER:=$HOME/forge-wallfit-bin/solver_density_cuda/build/convertGmshToForge}"
: "${FORGE_CUDA_BLOCKSIZE:=128}"
export FORGE_BIN REAL_CONVERTER FORGE_CUDA_BLOCKSIZE
export FORGE_CONVERTER="$(cd "$(dirname "$0")" && pwd)/conv_tolerant.sh"   # 終了時 GPUassert の既知の罠を許容
cd "$(dirname "$0")"
NPAR=${1:-4}
case "$NPAR" in ''|*[!0-9]*) echo "NPAR は正の整数: $NPAR"; exit 2;; esac
[ "$NPAR" -ge 1 ] || { echo "NPAR は 1 以上"; exit 2; }
TOOLS=../../solver_density_cuda/tools
IC=run_0114_euler_pin_G1_recal_ext6k
ARMB_REF=${ARMB_REF:-run_0143_euler_wallfit_monoG1_r1}
[ -f "${CASE_RUNS:-.}/$IC/res_6000.h5" ] || { echo "IC の run ${CASE_RUNS:-.}/$IC/res_6000.h5 が無い — 止める"; exit 2; }
[ -f "$ARMB_REF/prepare_info.json" ] && [ -f "$ARMB_REF/nozzle.h5" ] && [ -f "$ARMB_REF/IC_MAP.json" ] \
  || { echo "腕 B の参照 $ARMB_REF (prepare_info.json・nozzle.h5・IC_MAP.json) が無い — 止める"; exit 2; }
SEG=$(python3 -c 'from throat_mono_judge import SEGMENT_VERDICT_FILE as s; print(s)')   # 評価器と同じファイル名
declare -A RUNDIR=([M1]=run_0150_euler_wallfit_mocG1_r1 [M2]=run_0151_euler_wallfit_mocG1_r2 [M3]=run_0152_euler_wallfit_mocG1_r3
                   [ISEN]=run_0153_euler_icdep_mocG1_isen)
declare -A PREPOF=([M1]=_prep_moc_v5_mocG1 [M2]=_prep_moc_v5_mocG1 [M3]=_prep_moc_v5_mocG1 [ISEN]=_prep_moc_v5_mocG1_isen)
ARMB=(run_0143_euler_wallfit_monoG1_r1 run_0144_euler_wallfit_monoG1_r2 run_0145_euler_wallfit_monoG1_r3)
KEYS=(M1 M2 M3 ISEN)
RUNS=()
for k in "${KEYS[@]}"; do
  rd=${RUNDIR[$k]}
  if [ -e "$rd" ]; then echo "$rd が既にある — 番号の衝突。止める (既存 run は消さない)"; exit 2; fi
  RUNS+=("$rd")
done
if [ -z "${DRY:-}" ]; then
  for rd in "${ARMB[@]}"; do [ -d "$rd" ] || { echo "腕 B の run $rd が無い — 止める"; exit 2; }; done
  # 同じ実効設定 (§6 V5): forge のバイナリが腕 B と同じか (評価器も同じ照合をして、違えば保留にする)。違えば投入前に止める
  WANT=$(awk -F': ' '/^forge_sha256/{print $2; exit}' "${ARMB[0]}/RUN_PROVENANCE.txt" 2>/dev/null || true)
  HAVE=$(sha256sum "$FORGE_BIN" 2>/dev/null | awk '{print $1}' || true)
  if [ -z "$WANT" ] || [ "$WANT" != "$HAVE" ]; then
    if [ -n "${BIN_OVERRIDE:-}" ]; then
      mkdir -p _band_ab
      echo "$(date -Is) BIN_OVERRIDE=${BIN_OVERRIDE} (FORGE_BIN $FORGE_BIN $HAVE / 腕 B ${WANT:-記録なし})" >> _band_ab/moc_v5_bin_override.log
      echo "forge のバイナリが腕 B と違うが BIN_OVERRIDE で進める (評価器の前提「同じ実効設定」は不成立になる): ${BIN_OVERRIDE}"
    else
      echo "forge のバイナリ ($FORGE_BIN ${HAVE:-読めない}) が腕 B ${ARMB[0]} (${WANT:-記録なし}) と違う — 同じ実効設定でない。止める"; exit 2
    fi
  fi
fi
# 1. 準備 (1 回)。_prep_moc_v5_* はこのスクリプト専用の作業 dir。IC の検査が不成立なら prep が例外で止まる (番号写像を使わない)
rm -rf _prep_moc_v5_mocG1 _prep_moc_v5_mocG1_isen
if [ -n "${DRY:-}" ]; then
  if ! python3 moc_v5_euler.py prep _prep_moc_v5_mocG1 --with-isen _prep_moc_v5_mocG1_isen --armB-ref "$ARMB_REF" --dry \
       > _prep_moc_v5_mocG1.log 2>&1; then
    tail -20 _prep_moc_v5_mocG1.log; echo "DRY: prep が失敗 — _prep_moc_v5_mocG1.log"; exit 1
  fi
  grep -v "WARNING: cannot stamp" _prep_moc_v5_mocG1.log | tail -4
  echo "DRY: prep だけで止める (記録 _band_ab/moc_v5_ic_inspection_dry.json)"; exit 0
fi
if ! python3 moc_v5_euler.py prep _prep_moc_v5_mocG1 --with-isen _prep_moc_v5_mocG1_isen --armB-ref "$ARMB_REF" \
     > _prep_moc_v5_mocG1.log 2>&1; then
  tail -20 _prep_moc_v5_mocG1.log; echo "prep が失敗 (IC の検査の不成立を含む) — _prep_moc_v5_mocG1.log・_band_ab/moc_v5_ic_inspection.json"; exit 1
fi
tail -4 _prep_moc_v5_mocG1.log
# 2. run dir を作り (同じ prep の複製)、入力が写像・IC 生成の直後のままかを確かめ、NPAR 本ずつ並列に回す
for k in "${KEYS[@]}"; do
  cp -r "${PREPOF[$k]}" "${RUNDIR[$k]}"
done
python3 moc_v5_euler.py verify-prep _prep_moc_v5_mocG1 "${RUNDIR[M1]}" "${RUNDIR[M2]}" "${RUNDIR[M3]}" || { echo "腕 M の run の入力が prep と違う — 止める"; exit 2; }
python3 moc_v5_euler.py verify-prep _prep_moc_v5_mocG1_isen "${RUNDIR[ISEN]}" || { echo "ISEN の run の入力が prep と違う — 止める"; exit 2; }
# 終了コードは run dir の RUN_RC に書く (wait -n で回収済みの pid を後から wait すると状態が取れないため)
run_one() { local rc=0; python3 moc_v5_euler.py run "$1" > "$1/run_stdout.log" 2>&1 || rc=$?; echo "$rc" > "$1/RUN_RC"; }
for rd in "${RUNS[@]}"; do
  while [ "$(jobs -rp | wc -l)" -ge "$NPAR" ]; do wait -n || true; done
  run_one "$rd" &
  echo "started $rd (pid $!)"
done
wait || true
FAILED=()
for rd in "${RUNS[@]}"; do
  rc=$(cat "$rd/RUN_RC" 2>/dev/null || echo missing)
  if [ "$rc" = 0 ]; then echo "done $rd $(tail -1 "$rd/run_stdout.log")"; else echo "FAIL $rd (rc $rc)"; FAILED+=("$rd"); fi
done
# 3. 後処理: 残差図・収束判定 (本段の区間 = stage_manifest の最後の区間; 評価器は $SEG を読む)
for rd in "${RUNS[@]}"; do
  [ -f "$rd/residual_history.csv" ] || { echo "$rd: residual_history.csv が無い"; continue; }
  [ -f "$rd/residual_history.png" ] || python3 $TOOLS/plot_residual.py "$rd/residual_history.csv" -o "$rd/residual_history.png" > /dev/null 2>&1 || echo "$rd: 残差図の生成に失敗"
  python3 $TOOLS/check_convergence.py "$rd" --segment > "$rd/$SEG" 2>&1 || true
  echo "$rd: $(grep -m1 -E '^=== .*-> ' "$rd/$SEG" || tail -1 "$rd/$SEG")"
done
# 4. 評価: 7 本の時系列を同じ呼び出しで (欠けた run は評価器が保留にする) → 判定
mkdir -p _band_ab
ALL7=$(IFS=,; echo "${ARMB[*]},${RUNDIR[M1]},${RUNDIR[M2]},${RUNDIR[M3]},${RUNDIR[ISEN]}")
python3 eval_wallfit_euler.py . --fixed-coef --series="$ALL7" --geom-ref="${ARMB[0]}" --tag=v5 > _band_ab/moc_v5_euler_series.log 2>&1 \
  || echo "時系列の評価器が失敗 (rc $?) — _band_ab/moc_v5_euler_series.log"
tail -9 _band_ab/moc_v5_euler_series.log
python3 moc_v5_euler_eval.py . > _band_ab/moc_v5_euler_eval.log 2>&1 || echo "判定の評価器が失敗 (rc $?) — _band_ab/moc_v5_euler_eval.log"
tail -22 _band_ab/moc_v5_euler_eval.log
[ ${#FAILED[@]} -eq 0 ] || { echo "失敗した run: ${FAILED[*]}"; exit 1; }
echo ALLDONE
