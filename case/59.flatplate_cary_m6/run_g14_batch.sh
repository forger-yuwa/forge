#!/usr/bin/env bash
# G14 本体 (acceptance.json G14-E) の残り 8 系列を AWS で回す。2 レーン並列、1 本ごとに後片付けする。
#   bash run_g14_batch.sh <init_res.h5> <forge_tools_dir>
# 各 run: 細分化メッシュ fp_y1_0.35um_nx2800・全域 FP64・本段設定 40,000 step・出力 5000 毎。
# 完走後に check_convergence を回し、中間の全場スナップショット (res_<n>.h5) を消す
# (壁出力 res_wall_4_* と最終場・res_0 は残す。AWS のディスクが逼迫しているため)。
set -u
INIT=$1; TOOLS_DBL=$2
cd "$(dirname "$0")"
ROOT=$(git rev-parse --show-toplevel)

one() {  # run 名, 系列
  local run=$1 ser=$2
  if [ -e "$run" ]; then echo "[$run] 既にある。飛ばす"; return; fi
  python3 gen_runs.py --run "$run" --series "$ser" --mesh fp_y1_0.35um_nx2800 --init-from "$INIT" --init-same-mesh \
      --forge-tools "$TOOLS_DBL" --main-steps 40000 --out-int 5000 > "_gen_${run}.log" 2>&1
  python3 "$ROOT/solver_density_cuda/tools/check_convergence.py" "$run" > "$run/CONVERGENCE_VERDICT_check.txt" 2>&1
  local last; last=$(ls "$run" | grep -E '^res_[0-9]+\.h5$' | sed 's/res_//;s/\.h5//' | sort -n | tail -1)
  for f in "$run"/res_[0-9]*.h5; do
    b=$(basename "$f" .h5)
    [ "$b" = "res_$last" ] || [ "$b" = "res_0" ] || rm -f "$f" "$run/$b.xmf"
  done
  echo "[$run] 完了 $(head -1 "$run/CONVERGENCE_VERDICT_check.txt" | cut -c1-140)"
}

lane1() { one run_0006_tw07_re031 Re0.31_Tw0.7; one run_0008_tw05_re026 Re0.26_Tw0.5; one run_0010_tw03_re027 Re0.27_Tw0.3; one run_0012_tw04_re014 Re0.14_Tw0.4; }
lane2() { one run_0007_tw06_re026 Re0.26_Tw0.6; one run_0009_tw04_re027 Re0.27_Tw0.4; one run_0011_tw06_re014 Re0.14_Tw0.6; one run_0013_tw02_re014 Re0.14_Tw0.2; }
lane1 > _g14_lane1.log 2>&1 &
lane2 > _g14_lane2.log 2>&1 &
wait
echo "G14 BATCH DONE" | tee -a _g14_lane1.log
