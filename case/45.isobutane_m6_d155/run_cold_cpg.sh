#!/bin/bash
# SU2 との照合の forge 側 (plan tooling-nozzle-isothermal-wall-chain §5.1 #26): 壁・格子は TP の対と同じ、CFD だけ CPG の 4 本を並列に回す (AWS の case dir)。
# usage: bash run_cold_cpg.sh
#   run_0184 断熱・素の SST / run_0185 300 K・素の SST / run_0186 断熱・生産の SST / run_0187 300 K・生産の SST
#   初期値は収束した TP の解 (断熱 run_0181 の res_100000、300 K run_0183 の res_100000) を CPG に組み直したもの。
# **実行中にこのファイル・cold_pair.py を編集しないこと**。
set -euo pipefail
: "${FORGE_BIN:=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/forge}"
: "${REAL_CONVERTER:=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge}"
: "${FORGE_CUDA_BLOCKSIZE:=128}"; : "${MIN_FREE_GB:=12}"; : "${WATCH_SEC:=15}"
export FORGE_BIN REAL_CONVERTER FORGE_CUDA_BLOCKSIZE
unset FORGE_ALLOW_UNVERIFIED_SPECIES
cd "$(dirname "$0")"
export FORGE_CONVERTER="$(pwd)/conv_tolerant.sh"
die() { echo "止める: $*"; exit 2; }
free=$(df -BG --output=avail . | tail -1 | tr -dc 0-9)
[ "$free" -ge "$MIN_FREE_GB" ] || die "空き容量 ${free} GB < ${MIN_FREE_GB} GB"
RUNS=(run_0184_ns_coldmesh_cpg_ad_plain run_0185_ns_coldmesh_cpg_tw300_plain run_0186_ns_coldmesh_cpg_ad_dilat2 run_0187_ns_coldmesh_cpg_tw300_dilat2)
for r in "${RUNS[@]}"; do [ ! -e "$r" ] || die "$r が既にある"; done
python3 cold_pair.py prep-cpg run_0181_ns_coldmesh_ad run_0181_ns_coldmesh_ad/res_100000.h5 "${RUNS[0]}" plain
python3 cold_pair.py prep-cpg run_0182_ns_coldmesh_tw300 run_0183_ns_coldmesh_tw300_ext/res_100000.h5 "${RUNS[1]}" plain
python3 cold_pair.py prep-cpg run_0181_ns_coldmesh_ad run_0181_ns_coldmesh_ad/res_100000.h5 "${RUNS[2]}" dilat2
python3 cold_pair.py prep-cpg run_0182_ns_coldmesh_tw300 run_0183_ns_coldmesh_tw300_ext/res_100000.h5 "${RUNS[3]}" dilat2
pids=()
for r in "${RUNS[@]}"; do
  ( python3 cold_pair.py run "$r" > "$r/cold_pair_run_stdout.log" 2>&1 ) &
  pids+=($!); echo "started $r (pid ${pids[-1]})"
done
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
  for p in "${pids[@]}"; do kill -0 "$p" 2>/dev/null && alive=1; done
  for r in "${RUNS[@]}"; do [ -f "$r/EARLY_STOP.txt" ] || watch_nan "$r" || true; done
  sleep "$WATCH_SEC"
done
rc=0
for p in "${pids[@]}"; do wait "$p" || rc=1; done
for r in "${RUNS[@]}"; do python3 cold_pair.py nan-scan "$r" || true; done
echo "DONE rc=$rc"
