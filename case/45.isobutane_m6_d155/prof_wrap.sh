#!/bin/bash
# 計時用の包み (plan time_integration-line-implicit の速度の調査、2026-10-09)。run_case.sh の FORGE_BIN に渡すと、
# 同じディレクトリの prof_target.txt に書いた forge を nsys の下で起動する (CUDA のカーネル・memcpy を記録)。
T=$(cat "$(dirname "$(readlink -f "$0")")/prof_target.txt")
exec nsys profile -t cuda --stats=false -f true -o "$PWD/nsys_prof" "$T" "$@"
