#!/bin/bash
# 計時用の包み (plan time_integration-line-implicit-speed §6.2 の事後の内訳、2026-10-09)。prof_target_linevisc.txt の forge (lineI) を nsys の下で起動する。
T=$(cat "$(dirname "$(readlink -f "$0")")/prof_target_linevisc.txt")
exec nsys profile -t cuda --stats=false -f true -o "$PWD/nsys_prof" "$T" "$@"
