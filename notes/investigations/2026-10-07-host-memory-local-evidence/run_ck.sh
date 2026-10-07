#!/bin/bash
# 使い方: run_ck.sh <forge バイナリ> <tag> [ic_ckpt.h5 (再開時)]
# case/09 run_0160_passiveG_fct_ckpt100 の入力を scratch へ複製して 100 step 回す (元の case には書かない。計測専用)。
set -u
SP=/tmp/claude-1000/-home-sano-work-forge/60a080db-fad5-4e19-bdaa-3590a1035176/scratchpad
SRC=/home/sano/work/forge/case/09.Taylor-Green/run_0160_passiveG_fct_ckpt100
bin=$1; tag=$2; ic=${3:-}
d=$SP/hostmem/ck_${tag}
mkdir "$d" || exit 1
for f in solverConfig.yaml bcondConfig.yaml probe.yaml species_db.yaml Taylor-Green.h5 TG_stepxi_gaussY_seam.h5; do cp "$SRC/$f" "$d/" || exit 1; done
if [ -n "$ic" ]; then
  cp "$ic" "$d/ic_ckpt.h5" || exit 1
  sed -i 's/valueFileName: TG_stepxi_gaussY_seam.h5/valueFileName: ic_ckpt.h5/' "$d/solverConfig.yaml"
  echo "restart from $ic" > "$d/IC_FROM.txt"
fi
cd "$d"
sha256sum "$bin" > bin_sha256.txt
export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/hdf5/serial:${LD_LIBRARY_PATH:-}
FORGE_ALLOW_UNVERIFIED_SPECIES=1 /usr/bin/time -v -o time_v.txt "$bin" > forge_run.log 2>&1
rc=$?
echo "rc=$rc" > rc.txt
echo "$d rc=$rc $(grep -c . residual_history.csv 2>/dev/null) csv lines"
