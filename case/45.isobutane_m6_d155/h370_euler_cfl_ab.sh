#!/usr/bin/env bash
# V-h370 の Euler 参照の擬似 CFL の A/B (plans/active/tooling-nozzle-isothermal-wall-chain.md §5.1 #36j・§6、
# 2026-10-14 事前登録 b1470a4a、codex diagnose notes/reviews/2026-10-14-h370-euler-ref-plateau-diagnose.md)。
# AWS B で nohup で回す: cd ~/forge-wallfit/case/45.isobutane_m6_d155 && nohup bash h370_euler_cfl_ab.sh > /dev/null 2>&1 < /dev/null &
# 両腕とも run_0950_h370_euler_ref/res_60000.h5 から restart_field.py (同一メッシュ、--keep-src-dtype) で始め、
# run_0950 の本段の設定から nStepOuter (60000 → 40000) と cfl・cfl_pseudo だけを変える。腕は 1 本ずつ回す。
# **実行中にこのファイルを編集しないこと** (bash は逐次読み)。
set -uo pipefail
C45=$HOME/forge-wallfit/case/45.isobutane_m6_d155
SRC=$C45/run_0950_h370_euler_ref
T=$HOME/forge-g0cchain/solver_density_cuda/tools
LOG=$C45/h370_euler_cfl_ab.log
export FORGE_CUDA_BLOCKSIZE=128 FORGE_OMEGA_BUDGET=1 FORGE_DIAG_COMMIT_LOSS=500
export FORGE_BIN=$HOME/forge-fgeom7-fp64/solver_density_cuda/build/forge
for v in $(compgen -e | grep '^FORGE_'); do
  case $v in FORGE_CUDA_BLOCKSIZE|FORGE_OMEGA_BUDGET|FORGE_DIAG_COMMIT_LOSS|FORGE_BIN) ;; *) unset "$v";; esac
done
log() { echo "[$(date -Is)] $*" | tee -a "$LOG"; }
die() { log "中止: $*"; exit 1; }

log "== h370_euler_cfl_ab.sh 開始 (台本 sha $(sha256sum "$0" | cut -c1-16)、環境 $(compgen -e | grep '^FORGE_' | while read -r v; do printf '%s=%s ' "$v" "${!v}"; done))"
[ "$(sha256sum "$FORGE_BIN" | cut -c1-16)" = 7e95e5bbd8f50c9e ] || die "forge の sha が登録 (7e95e5bbd8f50c9e) と違う"
[ "$(sha256sum "$SRC/solverConfig.yaml" | cut -c1-16)" = 0889d4d0704a69ae ] || die "run_0950 の本段の設定の sha が 0889d4d0704a69ae でない"
[ -f "$SRC/res_60000.h5" ] || die "run_0950 の res_60000.h5 が無い"
log "継続元 res_60000.h5 sha $(sha256sum "$SRC/res_60000.h5" | cut -c1-16)"
avail=$(df --output=avail -BG "$C45" | tail -1 | tr -dc 0-9); [ "${avail:-0}" -ge 10 ] || die "ディスクの空き ${avail} GiB < 10 GiB"

for arm in "run_0952_h370_euler_cfl2 2.0" "run_0953_h370_euler_cfl05 0.5"; do
  set -- $arm; R=$C45/$1; V=$2
  [ -e "$R" ] && die "$R が既にある (消さない)"
  mkdir "$R" || die "mkdir $R"
  for f in solverConfig.yaml bcondConfig.yaml nozzle.h5 species_meta.yaml probe.yaml prepare_info.json wall_repr.json; do
    cp "$SRC/$f" "$R/" || die "$f の複製"
  done
  cp "$SRC"/resolved_species_*.yaml "$R/" || die "resolved_species の複製"
  [ -f "$SRC/nozzle.xmf" ] && cp "$SRC/nozzle.xmf" "$R/"
  python3 "$T/restart_field.py" "$SRC/res_60000.h5" "$R/nozzle.h5" --keep-src-dtype --dst-run "$R" > "$R/restart_field.log" 2>&1 \
    || { tail -20 "$R/restart_field.log" >> "$LOG"; die "restart_field.py $1"; }
  grep -i "verdict" "$R/restart_field.log" | tail -2 | sed "s/^/  $1 restart: /" | tee -a "$LOG"
  sed -i "s/nStepOuter: 60000/nStepOuter: 40000/; s/cfl: 2.0, cfl_pseudo: 2.0/cfl: $V, cfl_pseudo: $V/" "$R/solverConfig.yaml"
  diff "$SRC/solverConfig.yaml" "$R/solverConfig.yaml" > "$R/config_diff.txt"
  nd=$(grep -c '^[<>]' "$R/config_diff.txt")
  grep -q "nStepOuter: 40000" "$R/solverConfig.yaml" && grep -q "cfl: $V, cfl_pseudo: $V" "$R/solverConfig.yaml" && [ "$nd" -le 4 ] \
    || { cat "$R/config_diff.txt" >> "$LOG"; die "$1 の設定の書き換えが想定と違う (差の行 $nd)"; }
  log "$1: 設定の差 (run_0950 の本段から) $(grep '^>' "$R/config_diff.txt" | tr '\n' ' ')"
  log "起動 $1 (cfl_pseudo $V、40000 step)"
  bash "$T/run_case.sh" "$R" > "$R/run_case_stdout.log" 2>&1
  log "終了 $1 run_case rc=$? ($(grep -h 'forge exit=' "$R/run_case_stdout.log" | tail -1))"
done
python3 "$C45/h370_euler_cfl_ab_judge.py" >> "$LOG" 2>&1
log "== 判定 rc=$? → _band_ab/half_inlet/h370_euler_cfl_ab.json"
