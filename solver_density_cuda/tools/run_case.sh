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
python3 "$ROOT/solver_density_cuda/tools/check_convergence.py" . > CONVERGENCE_VERDICT.txt 2>&1 || true
cat CONVERGENCE_VERDICT.txt
python3 "$ROOT/solver_density_cuda/tools/plot_implicit_residuals.py" \
    --input residual_history.csv --output residual_history.png > /dev/null 2>&1 || true
exit $rc
