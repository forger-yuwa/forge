#!/usr/bin/env bash
# 速度計測用の短時間 run。<template_run> の入力 (config/bcond/species/probe + メッシュ h5 のハードリンク) から
# 専用ディレクトリ <template_run>_bench/<label>_n<N>/ を作り、nStepOuter=N・出力なしで forge を回して
# "Time = ... ms/step" (時間ループ内の壁時計; 初期化・I/O を含まない) をログに残す。
#
# 使い方: bench_steps.sh <template_run> <nsteps> <label>
#   環境変数: FORGE_BIN=<binary>   (既定 build/forge)
#             FORGE_PROFILE=1      (セクション別計時。cudaEventSynchronize が入るので合否判定には使わない。既定 0)
#             FORGE_ENV="K=V ..."  (追加の環境変数)
#   - template_run の入力は変更しない (専用 dir に複製)。res_*.h5 は書かない (outStepInterval を巨大化)。
#   - A/B は基準/変更版を交互に 2 回以上回し中央値で判定する (plan performance-3d-node-sst-speedup §4.3)。
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
TPL="$1"; NSTEP="$2"; LABEL="$3"
BIN="${FORGE_BIN:-$ROOT/solver_density_cuda/build/forge}"
export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/hdf5/serial:${LD_LIBRARY_PATH:-}
[ -d "$TPL" ] || { echo "[bench] template run dir not found: $TPL"; exit 1; }
TPL="$(cd "$TPL" && pwd)"
SRC_CFG="$TPL/solverConfig.yaml.orig"; [ -f "$SRC_CFG" ] || SRC_CFG="$TPL/solverConfig.yaml"
D="${TPL}_bench/${LABEL}_n${NSTEP}"
if [ -e "$D" ] && [ "${FORGE_BENCH_OVERWRITE:-0}" != "1" ]; then echo "[bench] $D exists; use a new label or FORGE_BENCH_OVERWRITE=1"; exit 2; fi
mkdir -p "$D"
for f in bcondConfig.yaml species_db.yaml probe.yaml IC_FROM.txt; do [ -f "$TPL/$f" ] && cp "$TPL/$f" "$D/"; done
for h in "$TPL"/*.h5; do case "$(basename "$h")" in res_*) ;; *) [ -e "$D/$(basename "$h")" ] || ln "$h" "$D/$(basename "$h")" 2>/dev/null || cp "$h" "$D/";; esac; done
python3 - "$SRC_CFG" "$D/solverConfig.yaml" "$NSTEP" "${FORGE_CFG_SUB:-}" <<'PY'
import re, sys
src, dst, n, sub = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
t = open(src).read()
t, c1 = re.subn(r'(nStepOuter:\s*)\d+', r'\g<1>' + n, t, count=1)
t, c2 = re.subn(r'(outStepInterval:\s*)\d+', r'\g<1>1000000', t, count=1)
if c1 != 1 or c2 != 1:
    raise SystemExit(f"[bench] config rewrite failed (nStepOuter matches={c1}, outStepInterval matches={c2})")
if sub:   # FORGE_CFG_SUB="old=new" (A/B 用の config 差し替え, \n 可)
    old, new = sub.split("=", 1)
    old = old.encode().decode("unicode_escape"); new = new.encode().decode("unicode_escape")
    if old not in t: raise SystemExit("[bench] FORGE_CFG_SUB old text not found in solverConfig")
    t = t.replace(old, new); print("[bench] cfgsub applied:", repr(new))
open(dst, "w").write(t)
PY
cd "$D"
LOG="bench_$(date +%Y%m%d_%H%M%S).log"
{
  echo "[bench] bin=$BIN sha256=$(sha256sum "$BIN" | cut -c1-16) template=$TPL nstep=$NSTEP label=$LABEL"
  echo "[bench] FORGE_PROFILE=${FORGE_PROFILE:-0} FORGE_ENV=${FORGE_ENV:-} gpu=$(nvidia-smi --query-gpu=name,clocks.sm,temperature.gpu --format=csv,noheader 2>/dev/null | head -1)"
  echo "[bench] inputs: $(sha256sum solverConfig.yaml bcondConfig.yaml 2>/dev/null | awk '{print substr($1,1,12), $2}' | tr '\n' ' ')"
} | tee "$LOG"
set +e
( export FORGE_PROFILE="${FORGE_PROFILE:-0}"; for kv in ${FORGE_ENV:-}; do export "$kv"; done; "$BIN" ) >> "$LOG" 2>&1
rc=$?
set -e
echo "[bench] exit=$rc"
grep -E "^Time = |Profiled steps|total_ms" "$LOG" | sed 's/^/  /'
steps_done=$(grep -c "^Step : " "$LOG" || true)
echo "[bench] steps logged=$steps_done (expected $NSTEP) log=$D/$LOG"
exit $rc
