#!/usr/bin/env bash
# forge を run ディレクトリで実行し、終了後に必ず収束の判定 (check_convergence.py、settle_criteria.yaml のある run の本段は
# check_settled.py; 選択の規則は下の後処理) を走らせて VERDICT を CONVERGENCE_VERDICT.txt に残す + 残差プロットを生成する。
#
# 使い方:  run_case.sh [run_dir]        (省略時は $PWD)
#   バックグラウンド長時間 run も:  nohup run_case.sh <run_dir> &
#
# これを唯一の forge 実行経路にすることで「forge を回したら収束チェック」を強制する
# (PreToolUse フックが直接 build/forge を弾く)。
#
# **実行中にこのファイルを編集しないこと** (2026-09-20 に事故): bash はスクリプトを
# バイト位置で逐次読むので、長時間 run の最中に行を挿入すると forge 終了後の再開位置が
# ずれ、コメントの途中から実行されて **forge 起動行をもう一度実行する**。
# 実例: `case/56/run_0009_gap_t8_uniform` が 800k step 完了直後に 2 周目を始め、
# `residual_history.csv` が上書きされた (`launch.log` に `command not found` が残る)。
# 編集が要るときは、走っている run が無いことを確かめてからにする。
set -uo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/hdf5/serial:${LD_LIBRARY_PATH:-}
RUNDIR="${1:-$PWD}"
cd "$RUNDIR" || { echo "run_case: bad run dir $RUNDIR"; exit 1; }

echo "[run_case] forge in $RUNDIR"
# --- 来歴を残す (codex result-3 M6) ---
# 既定値はソルバ側で変わることがある (例: 2026-09-20 に space.limiterScaled の既定が 0→1、
# venkatK が 1.0→0.05)。config に書いていないキーは run のログに現れないので、
# **どのバイナリで回したか**が分からないと過去の run がどの作用素の結果か決まらない。
FORGE_EXE="${FORGE_BIN:-$ROOT/solver_density_cuda/build/forge}"
{
  echo "date        : $(date -Is)"
  echo "forge_bin   : $FORGE_EXE"
  echo "forge_mtime : $(date -Is -r "$FORGE_EXE" 2>/dev/null || echo unknown)"
  echo "forge_sha256: $(sha256sum "$FORGE_EXE" 2>/dev/null | cut -d' ' -f1)"
  echo "git_head    : $(git -C "$ROOT" rev-parse HEAD 2>/dev/null || echo unknown)"
  echo "git_dirty   : $(git -C "$ROOT" status --porcelain 2>/dev/null | wc -l) 件"
  echo "host        : $(hostname)"
  # forge の挙動を変える環境変数を**値ごと**残す (2026-09-23)。
  # FORGE_CUDA_BLOCKSIZE は結果に効く (case/56 で 100k あたりの減衰が 1.40 % vs 1.88 %) のに
  # どこにも記録されず、後から run_0020/run_0021 を突き合わせたときに揃っているか確認できなかった。
  # FORGE_TRANSPORT_TABLE (0 = 輸送物性を表でなく double 評価)・FORGE_DIAG_FACE_H_DOUBLE (1 = SLAU の面エンタルピーを double) も
  # 解を変えるので残す (2026-10-15、plan tooling-nozzle-card-local-mesh §5.1 #3d の全部 double の延長)。
  for v in FORGE_CUDA_BLOCKSIZE FORGE_KERNEL_SYNC FORGE_NODE_FX_HALF FORGE_AXIS_DIAG_ALPHA FORGE_RESID_SNAP FORGE_TRANSPORT_TABLE FORGE_DIAG_FACE_H_DOUBLE FORGE_BIN; do
    if [ -n "${!v-}" ]; then echo "env         : $v=${!v}"; fi
  done
  # 解を変える opt-in の**要求値**を config から拾う (実効値は forge のログが出す)。
  # plan time_integration-fp64-accumulator §4.4 が要求。要求どおり効いたかは
  # forge_run.log の "[qAccumulatorFP64] 有効" 行と突き合わせる。
  for key in qAccumulatorFP64 lineImplicit lowMachPrecond wallTreatmentSST sstEnergyIncludesK; do
    v=$(grep -o "${key}: *[0-9][0-9]*" solverConfig.yaml 2>/dev/null | head -1 | grep -o '[0-9][0-9]*$')
    if [ -n "$v" ]; then echo "cfg         : $key=$v (要求値)"; fi
  done
} > RUN_PROVENANCE.txt

# FORGE_BIN で別ビルドの forge を指定できる (A/B 回帰で旧バイナリを回す用途。既定は build/forge)
"${FORGE_BIN:-$ROOT/solver_density_cuda/build/forge}" > forge_run.log 2>&1
rc=$?
echo "[run_case] forge exit=$rc"
# 実効設定 (起動エコー) を RUN_PROVENANCE にも写す (plan convection-slau-wall-normal-chi-default §4.2)
grep -h "^'slauWallNormalChi' effective" forge_run.log >> RUN_PROVENANCE.txt 2>/dev/null || true
grep -h "^'scalarGradient' effective" forge_run.log >> RUN_PROVENANCE.txt 2>/dev/null || true

echo "[run_case] convergence check ->"
# 判定の終了コードを **握りつぶさず記録する**。段階起動では段ごとに未収束でも先へ進むのが正しいので
# 既定では run_case 自体の終了コードには反映しないが、`FORGE_STRICT_CONVERGENCE=1` を立てると
# 未収束・判定不能を run_case の失敗として返す (2026-09-19 codex: `|| true` で飲み込んでいた)。
# 判定の道具の選択 (2026-10-16、plan tooling-settle-convergence-and-monitor §4.9 の 2): settle_criteria.yaml の無い run は今どおり
# check_convergence.py。ある run は `check_settled.py --select` が決める (本段の記録のある run・手で組んだ run → check_settled.py
# (既定で本段の起動と照合; 照合の不一致・起動の記録の欠け・前の段と確かめられない run は INDETERMINATE)、前の段と確かめられた run
# (MAIN_STAGE_STARTED が無く stage_manifest.json に本段の段が無い) → check_convergence.py; plan §5.1 #7c)。基準のある run は選択に
# 失敗しても旧の道具へ落とさない。1 行目に使った道具を書く
TOOLS="$ROOT/solver_density_cuda/tools"
judge_tool="check_convergence.py"; judge_why="settle_criteria.yaml が無い"
if [ -e settle_criteria.yaml ] || [ -L settle_criteria.yaml ]; then
    sel=$(python3 "$TOOLS/check_settled.py" --select . 2>/dev/null)
    judge_tool=$(printf '%s\n' "$sel" | sed -n 1p)
    judge_why=$(printf '%s\n' "$sel" | sed -n 2p)
    if [ "$judge_tool" != "check_convergence.py" ] && [ "$judge_tool" != "check_settled.py" ]; then
        judge_why="道具の選択に失敗 (check_settled.py --select の出力 '${judge_tool}'); 基準のある run なので check_settled.py"
        judge_tool="check_settled.py"
    fi
fi
if [ "$judge_tool" = "check_settled.py" ]; then
    { echo "TOOL: check_settled.py ($judge_why)"; python3 "$TOOLS/check_settled.py" . --bind; } > CONVERGENCE_VERDICT.txt 2>&1
    conv_rc=$?
    if ! grep -q '^VERDICT: ' CONVERGENCE_VERDICT.txt; then
        echo "VERDICT: INDETERMINATE (check_settled.py が VERDICT を出さずに終わった; 終了コード $conv_rc)" >> CONVERGENCE_VERDICT.txt
        conv_rc=3
    fi
else
    { echo "TOOL: check_convergence.py ($judge_why)"; python3 "$TOOLS/check_convergence.py" .; } > CONVERGENCE_VERDICT.txt 2>&1
    conv_rc=$?
fi
echo "CONVERGENCE_EXIT: $conv_rc" >> CONVERGENCE_VERDICT.txt
cat CONVERGENCE_VERDICT.txt
python3 "$ROOT/solver_density_cuda/tools/plot_implicit_residuals.py" \
    --input residual_history.csv --output residual_history.png > /dev/null 2>&1 || true
if [ "${FORGE_STRICT_CONVERGENCE:-0}" != "0" ] && [ "$rc" -eq 0 ] && [ "$conv_rc" -ne 0 ]; then
    echo "[run_case] FORGE_STRICT_CONVERGENCE=1 かつ収束判定が非ゼロ ($conv_rc) -> 失敗として返す"
    exit "$conv_rc"
fi
exit $rc
