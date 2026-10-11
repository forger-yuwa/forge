#!/bin/bash
# 計時用の包み (ncu、2026-10-09)。prof_target.txt の forge を、ライン Thomas の factor/solve を 1 回ずつ詳しく測る ncu の下で起動する。
T=$(cat "$(dirname "$(readlink -f "$0")")/prof_target.txt")
exec ncu --kernel-name regex:"lineThomas(Factor|Solve)_d|implicit_defect_correction_block_d" --launch-skip 30 --launch-count 6 \
  --section SpeedOfLight --section Occupancy --section LaunchStats --section ComputeWorkloadAnalysis --section MemoryWorkloadAnalysis \
  --section WarpStateStats --section InstructionStats -f -o "$PWD/ncu_prof" "$T" "$@"
