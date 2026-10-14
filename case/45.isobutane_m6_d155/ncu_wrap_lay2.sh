#!/bin/bash
# ncu の包み (plan time_integration-line-implicit-speed §5.1 #11、2026-10-10 ユーザ許可「sudo ncu」)。prof_target_linevisc.txt の forge を ncu の下で起動する。
# 性能カウンタは管理者限定 (RmProfilingAdminOnly 1) なので sudo -E で環境を引き継ぐ。並べ替え + 連鎖の短縮 (FORGE_LINE_LAYOUT=2) の 3 step 目の 1 step ぶん (分解 1・代入 5・block 5 = 11 回) を --set full で取る。
T=$(cat "$(dirname "$(readlink -f "$0")")/prof_target_linevisc.txt")
exec sudo -E env "PATH=$PATH" "LD_LIBRARY_PATH=${LD_LIBRARY_PATH:-}" /usr/local/cuda/bin/ncu --target-processes all \
  --kernel-name regex:"lineThomas(Factor|Solve)LP_d|implicit_defect_correction_block_d" --launch-skip 22 --launch-count 11 \
  --set full --import-source no -f -o "$PWD/ncu_prof" "$T" "$@"
