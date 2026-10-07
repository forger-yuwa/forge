#!/bin/bash
# 使い方: run_mem.sh <forge バイナリ> <格子 s070|s085> <tag>
# forgemem/prep_<grid> の入力を新しいディレクトリへ複製し、FORGE_MEMLOG=1・10 step で回す (計測専用)。
set -u
SP=/tmp/claude-1000/-home-sano-work-forge/60a080db-fad5-4e19-bdaa-3590a1035176/scratchpad
bin=$1; grid=$2; tag=$3
d=$SP/hostmem/mem_${grid}_${tag}
mkdir "$d" || exit 1
for f in solverConfig.yaml bcondConfig.yaml probe.yaml species_meta.yaml sern.h5; do cp "$SP/forgemem/prep_${grid}/$f" "$d/" || exit 1; done
sed -i 's/outStepInterval: 500/outStepInterval: 5/' "$d/solverConfig.yaml"
[ -n "${EXTRA_SED:-}" ] && sed -i "$EXTRA_SED" "$d/solverConfig.yaml"
grep -n "nStepOuter\|outStepInterval" "$d/solverConfig.yaml"
cd "$d"
sha256sum "$bin" > bin_sha256.txt
free -m > free_before.txt
export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/hdf5/serial:${LD_LIBRARY_PATH:-}
FORGE_MEMLOG=1 FORGE_ALLOW_UNVERIFIED_SPECIES=1 /usr/bin/time -v -o time_v.txt "$bin" > forge_run.log 2>&1
rc=$?
echo "rc=$rc" > rc.txt
free -m > free_after.txt
echo "$d rc=$rc maxRSS=$(grep 'Maximum resident' time_v.txt)"
