#!/bin/bash
# plan tooling-nozzle-throat-monotone-r2 §6 E1〜E4 (§5.1 #5・#5b・#6): 現行壁 (腕 A) と r″ 単調壁 (腕 B) の Euler A/B (生産 Euler 格子 G1、各腕 3 回)
# と、腕 B の IC 写像の予備 A/B。腕ごとに投入する (ARMS):
#  A   : run_0140〜0142 euler_wallfit_pinG1_r{1,2,3}  (problem_d155_euler_pin_G1_recal.yaml、IC = run_0114 最終場を restart_field)。prep を作り直す
#  B1  : run_0143 euler_wallfit_monoG1_r1             (problem_d155_euler_pin_G1_recal_mono.yaml、IC = 検証付き番号写像 ic_index_map --mode index)
#        run_0146 euler_icab_monoG1_nn                 (予備 A/B の α: 同じ B 格子 [B の prep を IC 前に複製]、最近傍対応で保存量を直接転送)
#        run_0143 が予備 A/B の β を兼ねる。B の prep を作り直す。終了後に予備 A/B の評価 (eval_wallfit_euler.py --icab) を書く
#  B23 : run_0144・0145 euler_wallfit_monoG1_r{2,3}   (B1 で作った B の prep をそのまま使う。verify-prep で run_0143 と同じ入力かを確認)
#        予備 A/B の VERDICT が「残留 IC 依存説を棄却」でなければ止める (§6 E1)。諮問のうえで進めるなら B23_OVERRIDE=<理由> を付ける
#        (理由は _band_ab/throat_mono_b23_override.log に残る)
#  B1 と B23 は同時に指定できない (予備 A/B の判定が先)。A は B1・B23 のどちらとも同時に指定できる。
#  段: soft (1 次 cfl 0.5、3000 step) → 本段 2 次 cfl 2・implicitRelax 0.7・18000 step・1000 ごと出力 (forge は run_case.sh 経由)。
#  終了後: residual_history.png の確認 (無ければ plot_residual.py)、check_convergence --segment (→ 各 run の CONVERGENCE_VERDICT_segment.txt;
#  評価器はこのファイルを読む)、評価器 (B1 は予備 A/B、6 本そろえば E3)。
# usage: ARMS=<A|B1|B23>[,…] bash run_throat_mono_ab.sh [NPAR (同時実行本数, 既定 3)]
#   例: ARMS=A,B1 bash run_throat_mono_ab.sh 5   /   ARMS=B23 bash run_throat_mono_ab.sh
#   バイナリ: 既定は AWS の ~/forge-wallfit-bin (FORGE_BIN / REAL_CONVERTER で上書き可)。AWS のケース dir は ~/forge-wallfit/case/45.isobutane_m6_d155。
#   既存の run dir があれば番号の衝突として止める (既存 run は消さない)。
#   **実行中にこのファイルを編集しないこと** (bash は逐次読み; run_case.sh の注記と同じ)。
set -euo pipefail
: "${FORGE_BIN:=$HOME/forge-wallfit-bin/solver_density_cuda/build/forge}"
: "${REAL_CONVERTER:=$HOME/forge-wallfit-bin/solver_density_cuda/build/convertGmshToForge}"
: "${FORGE_CUDA_BLOCKSIZE:=128}"
export FORGE_BIN REAL_CONVERTER FORGE_CUDA_BLOCKSIZE
export FORGE_CONVERTER="$(cd "$(dirname "$0")" && pwd)/conv_tolerant.sh"   # 終了時 GPUassert の既知の罠を許容
cd "$(dirname "$0")"
NPAR=${1:-3}
case "$NPAR" in ''|*[!0-9]*) echo "NPAR は正の整数: $NPAR"; exit 2;; esac
[ "$NPAR" -ge 1 ] || { echo "NPAR は 1 以上"; exit 2; }
[ -n "${ARMS:-}" ] || { echo "ARMS を指定する (A / B1 / B23 をカンマ区切り; 例 ARMS=A,B1)"; exit 2; }
PHASES=()
for p in ${ARMS//,/ }; do
  case "$p" in A|B1|B23) PHASES+=("$p");; *) echo "ARMS の値が不正: $p (A / B1 / B23)"; exit 2;; esac
done
has() { local x; for x in "${PHASES[@]}"; do [ "$x" = "$1" ] && return 0; done; return 1; }
if has B1 && has B23; then echo "B1 と B23 は同時に指定できない (予備 A/B の判定が先; plan §6 E1)"; exit 2; fi
TOOLS=../../solver_density_cuda/tools
IC=run_0114_euler_pin_G1_recal_ext6k
[ -f "$IC/res_6000.h5" ] || { echo "IC の run $IC/res_6000.h5 が無い — 止める"; exit 2; }
SEG=$(python3 -c 'from throat_mono_judge import SEGMENT_VERDICT_FILE as s; print(s)')   # 評価器と同じファイル名
declare -A RUNDIR=([A1]=run_0140_euler_wallfit_pinG1_r1 [A2]=run_0141_euler_wallfit_pinG1_r2 [A3]=run_0142_euler_wallfit_pinG1_r3
                   [B1]=run_0143_euler_wallfit_monoG1_r1 [B2]=run_0144_euler_wallfit_monoG1_r2 [B3]=run_0145_euler_wallfit_monoG1_r3
                   [NN]=run_0146_euler_icab_monoG1_nn)
declare -A PREPOF=([A1]=pinG1 [A2]=pinG1 [A3]=pinG1 [B1]=monoG1 [B2]=monoG1 [B3]=monoG1 [NN]=monoG1_nn)
KEYS=()
if has A; then KEYS+=(A1 A2 A3); fi
if has B1; then KEYS+=(B1 NN); fi
if has B23; then KEYS+=(B2 B3); fi
RUNS=()
for k in "${KEYS[@]}"; do
  rd=${RUNDIR[$k]}
  if [ -e "$rd" ]; then echo "$rd が既にある — 番号の衝突。止める (既存 run は消さない)"; exit 2; fi
  RUNS+=("$rd")
done
# 0. B23 の前提 (準備より先に確かめる): B の r1 と同じ prep が残っていること、予備 A/B の判定
if has B23; then
  [ -d "${RUNDIR[B1]}" ] || { echo "B23: ${RUNDIR[B1]} (腕 B の r1) が無い — 先に ARMS=B1"; exit 2; }
  python3 throat_mono_ab.py verify-prep _prep_throat_mono_monoG1 "${RUNDIR[B1]}" || { echo "B23: B の prep が r1 と同じ入力でない — 止める"; exit 2; }
  V=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["VERDICT"])' _band_ab/throat_mono_icab_fixedcoef.json 2>/dev/null || echo "missing (_band_ab/throat_mono_icab_fixedcoef.json)")
  case "$V" in
    残留\ IC\ 依存説を棄却*) echo "予備 A/B: $V";;
    *) if [ -n "${B23_OVERRIDE:-}" ]; then
         mkdir -p _band_ab
         echo "$(date -Is) B23_OVERRIDE=${B23_OVERRIDE} (予備 A/B VERDICT: $V)" >> _band_ab/throat_mono_b23_override.log
         echo "予備 A/B が棄却でない ($V) が B23_OVERRIDE で進める: ${B23_OVERRIDE}"
       else
         echo "B23: 予備 A/B の VERDICT が「残留 IC 依存説を棄却」でない ($V) — 諮問してから B23_OVERRIDE=<理由> で進める"; exit 2
       fi;;
  esac
fi
# 1. 準備 (腕ごとに 1 回)。_prep_throat_mono_* はこのスクリプト専用の作業 dir。B23 は B1 の prep を使う (作り直さない)
if has A; then
  rm -rf _prep_throat_mono_pinG1
  python3 throat_mono_ab.py prep _prep_throat_mono_pinG1 pinG1 2>&1 | tee _prep_throat_mono_pinG1.log | tail -3
fi
if has B1; then
  rm -rf _prep_throat_mono_monoG1 _prep_throat_mono_monoG1_nn
  python3 throat_mono_ab.py prep _prep_throat_mono_monoG1 monoG1 --with-nn _prep_throat_mono_monoG1_nn 2>&1 | tee _prep_throat_mono_monoG1.log | tail -6
fi
# 2. run dir を作り、NPAR 本ずつ並列に回す
for k in "${KEYS[@]}"; do
  cp -r "_prep_throat_mono_${PREPOF[$k]}" "${RUNDIR[$k]}"
done
if has B23; then
  python3 throat_mono_ab.py verify-prep _prep_throat_mono_monoG1 "${RUNDIR[B1]}" "${RUNDIR[B2]}" "${RUNDIR[B3]}" || { echo "B23: 複製した run の入力が r1 と違う — 止める"; exit 2; }
fi
# 終了コードは run dir の RUN_RC に書く (wait -n で回収済みの pid を後から wait すると状態が取れないため)
run_one() { local rc=0; python3 throat_mono_ab.py run "$1" > "$1/run_stdout.log" 2>&1 || rc=$?; echo "$rc" > "$1/RUN_RC"; }
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
# 4. 評価 (固定係数の刻み感度)
mkdir -p _band_ab
if has B1; then   # 予備 A/B (§6 E1): α 最近傍 run_0146 / β 番号写像 run_0143、本段の連続 10 枚
  python3 eval_wallfit_euler.py . --fixed-coef --icab="${RUNDIR[NN]},${RUNDIR[B1]}" > _band_ab/throat_mono_icab.log 2>&1 || echo "予備 A/B の評価器が失敗 (rc $?) — _band_ab/throat_mono_icab.log"
  tail -16 _band_ab/throat_mono_icab.log
fi
MISSING=()
for k in A1 A2 A3 B1 B2 B3; do
  if [ "$(cat "${RUNDIR[$k]}/RUN_RC" 2>/dev/null || echo x)" != 0 ]; then MISSING+=("${RUNDIR[$k]}"); fi
done
if [ ${#MISSING[@]} -eq 0 ]; then   # E3 (§6 E2〜E4): 腕 A・B の 6 本がそろったとき
  python3 eval_wallfit_euler.py . --fixed-coef --e3 --pair=pinG1,monoG1 > _band_ab/throat_mono_euler_ab.log 2>&1 || echo "E3 の評価器が失敗 (rc $?) — _band_ab/throat_mono_euler_ab.log"
  tail -16 _band_ab/throat_mono_euler_ab.log
else
  echo "E3 の評価は 6 本がそろってから (未完了: ${MISSING[*]})"
fi
[ ${#FAILED[@]} -eq 0 ] || { echo "失敗した run: ${FAILED[*]}"; exit 1; }
echo ALLDONE
