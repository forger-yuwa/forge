#!/usr/bin/env bash
# forge を run ディレクトリで実行し、終了後に必ず check_convergence.py を走らせて
# VERDICT を CONVERGENCE_VERDICT.txt に残す + 残差プロットを生成する。
#
# 使い方:  run_case.sh [run_dir]        (省略時は $PWD)
#   バックグラウンド長時間 run も:  nohup run_case.sh <run_dir> &
#
# これを唯一の forge 実行経路にすることで「forge を回したら収束チェック」を強制する
# (PreToolUse フックが直接 build/forge を弾く)。
set -uo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/hdf5/serial:${LD_LIBRARY_PATH:-}
RUNDIR="${1:-$PWD}"
cd "$RUNDIR" || { echo "run_case: bad run dir $RUNDIR"; exit 1; }

echo "[run_case] forge in $RUNDIR"
# FORGE_BIN で別ビルドの forge を指定できる (A/B 回帰で旧バイナリを回す用途。既定は build/forge)
"${FORGE_BIN:-$ROOT/solver_density_cuda/build/forge}" > forge_run.log 2>&1
rc=$?
echo "[run_case] forge exit=$rc"

echo "[run_case] convergence check ->"
# 判定の終了コードを **握りつぶさず記録する**。段階起動では段ごとに未収束でも先へ進むのが正しいので
# 既定では run_case 自体の終了コードには反映しないが、`FORGE_STRICT_CONVERGENCE=1` を立てると
# 未収束・判定不能を run_case の失敗として返す (2026-09-19 codex: `|| true` で飲み込んでいた)。
python3 "$ROOT/solver_density_cuda/tools/check_convergence.py" . > CONVERGENCE_VERDICT.txt 2>&1
conv_rc=$?
echo "CONVERGENCE_EXIT: $conv_rc" >> CONVERGENCE_VERDICT.txt
cat CONVERGENCE_VERDICT.txt
python3 "$ROOT/solver_density_cuda/tools/plot_implicit_residuals.py" \
    --input residual_history.csv --output residual_history.png > /dev/null 2>&1 || true
if [ "${FORGE_STRICT_CONVERGENCE:-0}" != "0" ] && [ "$rc" -eq 0 ] && [ "$conv_rc" -ne 0 ]; then
    echo "[run_case] FORGE_STRICT_CONVERGENCE=1 かつ収束判定が非ゼロ ($conv_rc) -> 失敗として返す"
    exit "$conv_rc"
fi
exit $rc
