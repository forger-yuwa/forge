#!/bin/bash
# case/71 run_0002_ns の判定ツールを今の場で回す (AWS A、forge は起動しない)。plan tooling-nozzle-card-local-mesh §5.1 #3b。
# run dir には何も書かない (--segment は scratch の写しで回す: check_convergence --segment は residual_history_segment.csv を run dir に書くため)。
# 出力は $O (原文)。手元の case/71.m55_combustion_wt/_runs/run_0002_ns/judge/ に写して commit する。
set -uo pipefail
C=$HOME/forge-cardmesh/case/71.m55_combustion_wt
R=$C/run_0002_ns
E=$C/run_0001_euler
T=$HOME/forge-cardmesh/solver_density_cuda/tools
DS=$HOME/forge-cardmesh/design
W=$HOME/m55_diag_0002
O=$W/out
mkdir -p "$O"
S=$(dirname "$(readlink -f "$0")")
echo "# $(date -Iseconds)  git $(git -C $HOME/forge-cardmesh rev-parse HEAD)  python3 $(python3 -V 2>&1)" > "$O/_env.txt"

# 1a. --segment (記録された stage_manifest.json のまま; 本段はドライバが OSError で落ちて manifest に入っていない)
A1=$W/seg_asis/run_0002_ns; rm -rf "$W/seg_asis"; mkdir -p "$A1"
cp -p "$R/stage_manifest.json" "$R/forge_launches.jsonl" "$A1/"
for f in residual_history_S1_soft.csv residual_history_S2_mid.csv residual_history.csv; do ln -s "$R/$f" "$A1/$f"; done
{ python3 "$T/stage_manifest.py" "$A1" --segments; echo "STAGE_MANIFEST_EXIT: $?"; } > "$O/stage_manifest_segments_asis.txt" 2>&1
{ python3 "$T/check_convergence.py" --segment "$A1"; echo "CHECK_CONVERGENCE_EXIT: $?"; } > "$O/check_convergence_segment_asis.txt" 2>&1

# 1b. --segment (本段を manifest に足した写し: run_staged_ns の _record(cfg_main, "main") と同じ StageManifest.add を run の
#     solverConfig.yaml (本段) と bcondConfig.yaml で行う。履歴は residual_history.csv (本段だけ) を residual_history_main.csv として)
A2=$W/seg_with_main/run_0002_ns; rm -rf "$W/seg_with_main"; mkdir -p "$A2"
cp -p "$R/forge_launches.jsonl" "$A2/"
for f in residual_history_S1_soft.csv residual_history_S2_mid.csv; do ln -s "$R/$f" "$A2/$f"; done
ln -s "$R/residual_history.csv" "$A2/residual_history_main.csv"
python3 - "$T" "$R" "$A2" > "$O/stage_manifest_add_main.txt" 2>&1 <<'PY'
import json, sys
sys.path.insert(0, sys.argv[1])
from stage_manifest import StageManifest
R, A2 = sys.argv[2], sys.argv[3]
sm = StageManifest(A2)
sm.stages = json.load(open(R + "/stage_manifest.json"))["stages"]
sm.add("main", open(R + "/solverConfig.yaml").read(), open(R + "/bcondConfig.yaml").read(), history="residual_history_main.csv")
sm.write()
print(json.dumps(sm.stages[-1], indent=1, ensure_ascii=False))
PY
{ python3 "$T/stage_manifest.py" "$A2" --segments; echo "STAGE_MANIFEST_EXIT: $?"; } > "$O/stage_manifest_segments_with_main.txt" 2>&1
{ python3 "$T/check_convergence.py" --segment "$A2"; echo "CHECK_CONVERGENCE_EXIT: $?"; } > "$O/check_convergence_segment_with_main.txt" 2>&1
cp -p "$A2/stage_manifest.json" "$O/stage_manifest_with_main.json"

# 1c. 本段の履歴だけ (residual_history.csv) の通常判定 (run_case.sh の CONVERGENCE_VERDICT.txt と同じ入力)
{ python3 "$T/check_convergence.py" "$R"; echo "CHECK_CONVERGENCE_EXIT: $?"; } > "$O/check_convergence_main_only.txt" 2>&1

# 1d. Euler の参照 run_0001_euler の --segment (manifest は本段まである)
A3=$W/seg_euler/run_0001_euler; rm -rf "$W/seg_euler"; mkdir -p "$A3"
cp -p "$E/stage_manifest.json" "$E/forge_launches.jsonl" "$A3/"
for f in "$E"/residual_history*.csv; do case "$f" in *_segment.csv) ;; *) ln -s "$f" "$A3/$(basename "$f")" ;; esac; done
{ python3 "$T/stage_manifest.py" "$A3" --segments; echo "STAGE_MANIFEST_EXIT: $?"; } > "$O/stage_manifest_segments_euler.txt" 2>&1
{ python3 "$T/check_convergence.py" --segment "$A3"; echo "CHECK_CONVERGENCE_EXIT: $?"; } > "$O/check_convergence_segment_euler.txt" 2>&1

# 3. θ_r の測定帯の A/B (停止規則の系列の CSV も作る)
{ python3 "$S/theta_band_ab.py" --design "$DS" --run "$R" --euler-run "$E" --band-b-dir "$W/band_B_183000" \
    --out-json "$O/theta_band_ab.json" --series-csv "$O/stop_rule_series.csv"; echo "EXIT: $?"; } > "$O/theta_band_ab.stdout.txt" 2>&1

# 2a. 停止規則の系列の準定常 (tail は既定)
{ python3 "$T/check_quasisteady.py" --series-csv "$O/stop_rule_series.csv" \
    --series-cols r,mdot_in,mdot_out,theta_r_0.4,theta_r_0.7,theta_r_0.95,exit_core_M; echo "CHECK_QUASISTEADY_EXIT: $?"; } \
    > "$O/check_quasisteady_series.txt" 2>&1

# 2b. 最終の場の非有限・物理性
{ python3 "$S/final_field_check.py" --tools "$T" --run "$R" --step 183000 --out-json "$O/final_field_check.json"; echo "EXIT: $?"; } \
    > "$O/final_field_check.stdout.txt" 2>&1

# 2c. 壁解像 (no-slip 壁は bcond の wall (physID 3, wall_isothermal) だけ)
ls "$R"/res_wall_3_183000.h5 > "$O/check_wall_resolution.txt" 2>&1
{ python3 "$T/check_wall_resolution.py" "$R" --step 183000 --groups wall --weight area; echo "CHECK_WALL_RESOLUTION_EXIT: $?"; } \
    >> "$O/check_wall_resolution.txt" 2>&1

# 5. 実効の設定の照合 (問題 YAML → 生成器)
{ python3 "$S/config_regen_check.py" --design "$DS" --tools "$T" --problem "$C/problem_m55_final.yaml" --run "$R" \
    --out-json "$O/config_regen_check.json"; echo "EXIT: $?"; } > "$O/config_regen_check.stdout.txt" 2>&1

# 4. 容量
{ date -Iseconds; df -B1 "$HOME"; df -h "$HOME"; echo; du -sb "$R" "$E"; echo
  stat -c '%s %n' "$R/res_183000.h5" "$R/res_0.h5" "$R/res_wall_3_183000.h5" "$R/res_outlet_2_183000.h5" "$R/nozzle.h5" \
     "$R/residual_history.csv" "$R/forge_run.log" "$E/res_12000.h5"
  echo "res_N.h5 の数 $(ls "$R" | grep -cE '^res_[0-9]+\.h5$')  res_wall_3 の数 $(ls "$R" | grep -cE '^res_wall_3_[0-9]+\.h5$')  res_outlet_2 の数 $(ls "$R" | grep -cE '^res_outlet_2_[0-9]+\.h5$')"
} > "$O/disk.txt" 2>&1
echo DONE > "$O/_done.txt"
