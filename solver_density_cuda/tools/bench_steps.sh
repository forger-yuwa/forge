#!/usr/bin/env bash
# 速度計測用の短時間 run: run ディレクトリの solverConfig.yaml を nStepOuter=N・出力なしに書き換えて
# FORGE_PROFILE=1 で forge を回し、Time / ms/step / Runtime Profile Summary を bench_<label>.log に残す。
#
# 使い方: bench_steps.sh <run_dir> <nsteps> <label> [FORGE_BIN=...] [FORGE_ENV="K=V K=V"]
#   - run_dir は複製済みの run_* (入力 h5 + solverConfig/bcondConfig) であること。res_*.h5 は書かない。
#   - 同じ run_dir で label を変えて A/B する (バイナリは FORGE_BIN で切替)。
#   - 起動 (メッシュ読込・初期化) を相殺するには nsteps を 2 種類 (例 100/300) で回し差分を取る。
set -uo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
RUNDIR="$1"; NSTEP="$2"; LABEL="$3"
BIN="${FORGE_BIN:-$ROOT/solver_density_cuda/build/forge}"
export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/hdf5/serial:${LD_LIBRARY_PATH:-}
cd "$RUNDIR" || exit 1
[ -f solverConfig.yaml.orig ] || cp solverConfig.yaml solverConfig.yaml.orig
python3 - "$NSTEP" <<'PY'
import re, sys
n = sys.argv[1]
t = open("solverConfig.yaml.orig").read()
t, c1 = re.subn(r'(nStepOuter:\s*)\d+', r'\g<1>' + n, t, count=1)
t, c2 = re.subn(r'(outStepInterval:\s*)\d+', r'\g<1>1000000', t, count=1)
assert c1 == 1 and c2 == 1, (c1, c2)
open("solverConfig.yaml", "w").write(t)
PY
LOG="bench_${LABEL}_n${NSTEP}.log"
echo "[bench] $BIN  nstep=$NSTEP  label=$LABEL  env=${FORGE_ENV:-}" | tee "$LOG"
( export FORGE_PROFILE=1; for kv in ${FORGE_ENV:-}; do export "$kv"; done; /usr/bin/time -f "wall=%e s maxrss=%M KB" "$BIN" ) >> "$LOG" 2>&1
rc=$?
echo "[bench] exit=$rc"
grep -E "^Time = |wall=|Profiled steps|total_ms" "$LOG" | sed 's/^/  /'
python3 - "$LOG" "$NSTEP" <<'PY'
import re, sys
t = open(sys.argv[1]).read(); n = int(sys.argv[2])
m = re.search(r'^Time = ([0-9.]+) s', t, re.M)
if m: print(f"[bench] Time={float(m.group(1)):.2f}s  ->  {1000*float(m.group(1))/n:.2f} ms/step (incl. startup)")
PY
rm -f res_*.h5 res_*.xmf
exit $rc
