#!/bin/bash
# plan tooling-sern-te-wake-grid §5.1 R2c・§6 #11 (事前登録 2026-10-12、codex diagnose 2026-10-12 M3):
# 実 GPU で標準経路の判定区間の段が「未定常だけ → 自動延長 → PASS」を通るかを 1 件。起点 = run_1094_tewake_lsw10_m6/res_4000.h5
# (run_1097 の起点と同じ)、生産 YAML problem_3d_prod_3op_wallres_lswx10.yaml の m6_on、N 20000・W 10000・延長 1 回 (YAML から)。
#
# 使い方 (AWS-B の ~/forge-r8 = feature/sern-design の checkout、solver_density_cuda/build/forge あり。長いので nohup で):
#   nohup case/46.sern_design/diag/judgment_r2c.sh <新しい run 名> > case/46.sern_design/<新しい run 名>.r2c.log 2>&1 &
#   (追加のオプションは judgment_r2c.py に渡る: --src-run / --src-res / --problem / --op / --poll-s / --min-free-gb / --case-dir)
#
# 段: (1)〜(3) と (4) の前半 = prepare (格子生成・変換・品質、格子署名と双対幾何のハッシュの照合、本段の config の差の記録、
#     起点の場の配置。どれかが合わなければ時間積分の前に止める) → (4) judge (runner_sern.run_judgment、te_monitor 付き・強制の分岐なし)
#     → (5) post (check_convergence・check_quasisteady・check_floor_events を run に残し、R2C_SUMMARY.txt)
# 環境: FORGE_CUDA_BLOCKSIZE=128 を明示して渡す。FORGE_ALLOW_UNVERIFIED_SPECIES は既定で外し、起点の場に化学種の属性が無いとき
#     (prepare の判断、R2C_PREPARE.json の allow_unverified_species) だけ judge の実行に付ける (procedures/solver-settings.md の
#     「その実行だけ」の許可)。実際の値は JUDGMENT.json・runner_launches.jsonl・R2C_JUDGE.json に残る
# 終了コード: prepare で止めたらその番号 (2 既存・3 格子・4 config・5 化学種・6 入力/環境)、それ以外は judge の終了コード
set -u
ROOT=$(cd "$(dirname "$0")/../../.." && pwd)
cd "$ROOT" || exit 1
NEW=${1:?usage: judgment_r2c.sh <新しい run 名> [judgment_r2c.py のオプション]}
shift
PY=${PYTHON:-python3}
T=case/46.sern_design/diag/judgment_r2c.py
export FORGE_CUDA_BLOCKSIZE=128          # node SLAU のレジスタ上限 (既定 512 では起動できない)。検証の手順と同じ値
unset FORGE_ALLOW_UNVERIFIED_SPECIES     # 既定で付けない (要るときだけ judge の実行に付ける)
echo "R2C_START $NEW $(date) git $(git rev-parse HEAD 2>/dev/null)"
$PY $T prepare "$NEW" "$@"
rc=$?
if [ $rc -ne 0 ]; then echo "R2C: prepare で停止 (rc=$rc) $(date)"; exit $rc; fi
RD=$($PY $T rundir "$NEW" "$@") || exit 1
allow=$($PY -c "import json,sys; print(1 if json.load(open(sys.argv[1]))['allow_unverified_species'] else 0)" "$RD/R2C_PREPARE.json") || exit 1
if [ "$allow" = 1 ]; then
  echo "R2C: 起点の場に化学種の属性が無い -> judge の実行に限って FORGE_ALLOW_UNVERIFIED_SPECIES=1"
  FORGE_ALLOW_UNVERIFIED_SPECIES=1 $PY $T judge "$NEW" "$@"
else
  $PY $T judge "$NEW" "$@"
fi
jrc=$?
$PY $T post "$NEW" "$@"
prc=$?
echo "R2C_DONE $NEW judge rc=$jrc post rc=$prc $(date)"
exit $jrc
