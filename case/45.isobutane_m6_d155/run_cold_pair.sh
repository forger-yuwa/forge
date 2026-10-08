#!/bin/bash
# 冷却壁の NS の対 (plan tooling-nozzle-isothermal-wall-chain §5.1 #14、§6 V-c45) の投入 (AWS の case dir で実行)。
# usage: bash run_cold_pair.sh <ad_run> <tw300_run>
#   2 本を準備 (cold_pair.py prep) してから並列に回す (cold_pair.py run → runner_axismach.run_staged_ns → run_case.sh)。
#   forge の sha256 は生産の run_0167 の RUN_PROVENANCE と同じであること (別のバイナリで回さない)。
#   残差の rms_* (rms_dq_* を除く) に NaN・Inf が出たら、その run の forge (cwd がその run の forge のうち、このスクリプトが起動したもの) を止める。
# **実行中にこのファイル・cold_pair.py を編集しないこと**。
set -euo pipefail
: "${FORGE_BIN:=$HOME/forge-wallfit-bin/solver_density_cuda/build/forge}"
: "${REAL_CONVERTER:=$HOME/forge-wallfit-bin/solver_density_cuda/build/convertGmshToForge}"
: "${FORGE_CUDA_BLOCKSIZE:=128}"; : "${MIN_FREE_GB:=8}"; : "${WATCH_SEC:=10}"
export FORGE_BIN REAL_CONVERTER FORGE_CUDA_BLOCKSIZE
unset FORGE_ALLOW_UNVERIFIED_SPECIES
cd "$(dirname "$0")"
export FORGE_CONVERTER="$(pwd)/conv_tolerant.sh"
die() { echo "止める: $*"; exit 2; }
[ $# -eq 2 ] || die "usage: bash run_cold_pair.sh <ad_run> <tw300_run>"
AD=$1; TW=$2
for r in "$AD" "$TW"; do [[ "$r" =~ ^run_[0-9]{4}_[A-Za-z0-9_.-]+$ ]] || die "run 名 '$r'"; [ ! -e "$r" ] || die "$r が既にある"; done
want=$(awk -F': ' '/^forge_sha256/{print $2; exit}' run_0167_ns_n012_N2/RUN_PROVENANCE.txt)
have=$(sha256sum "$FORGE_BIN" | awk '{print $1}')
[ -n "$want" ] && [ "$want" = "$have" ] || die "forge の sha256 が run_0167 と違う ($have vs $want)"
free=$(df -BG --output=avail . | tail -1 | tr -dc 0-9)
[ "$free" -ge "$MIN_FREE_GB" ] || die "空き容量 ${free} GB < ${MIN_FREE_GB} GB"
echo "forge $FORGE_BIN ($have)  空き ${free} GB"
python3 cold_pair.py prep ad "$AD"
python3 cold_pair.py prep tw300 "$TW"
pids=()
for r in "$AD" "$TW"; do
  ( python3 cold_pair.py run "$r" > "$r/cold_pair_run_stdout.log" 2>&1 ) &
  pids+=($!)
  echo "started $r (pid ${pids[-1]})"
done
# 見張り: 残差の NaN・Inf (段ごとの CSV も含む) で、その run の forge だけを止める
watch_nan() {
  local r=$1 f bad
  for f in "$r"/residual_history*.csv; do
    [ -f "$f" ] || continue
    bad=$(python3 - "$f" <<'PY'
import csv, math, sys
with open(sys.argv[1]) as fh:
    rd = csv.reader(fh); head = [h.strip() for h in next(rd, [])]
    cols = [i for i, h in enumerate(head) if h.startswith("rms_") and not h.startswith("rms_dq_")]
    for n, row in enumerate(rd, 1):
        for i in cols:
            try: v = float(row[i])
            except (ValueError, IndexError): continue
            if not math.isfinite(v): print(f"row {n} {head[i]}"); sys.exit()
PY
)
    if [ -n "$bad" ]; then
      echo "$(date -Is) NaN/Inf in $f: $bad" | tee -a "$r/EARLY_STOP.txt"
      for p in $(pgrep -x forge); do
        [ "$(readlink -f /proc/$p/cwd 2>/dev/null)" = "$(readlink -f "$r")" ] && kill "$p" && echo "killed forge $p ($r)" | tee -a "$r/EARLY_STOP.txt"
      done
      return 1
    fi
  done
  return 0
}
alive=1
while [ $alive -eq 1 ]; do
  alive=0
  for i in 0 1; do
    if kill -0 "${pids[$i]}" 2>/dev/null; then alive=1; fi
  done
  for r in "$AD" "$TW"; do [ -f "$r/EARLY_STOP.txt" ] || watch_nan "$r" || true; done
  sleep "$WATCH_SEC"
done
rc=0
for i in 0 1; do wait "${pids[$i]}" || rc=1; done
for r in "$AD" "$TW"; do python3 cold_pair.py nan-scan "$r" || true; done
echo "DONE rc=$rc"
