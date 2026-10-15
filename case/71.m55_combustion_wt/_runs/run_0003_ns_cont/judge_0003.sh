#!/bin/bash
# case/71 run_0003_ns_cont の判定ツールを最後の場で回す (AWS A、forge は起動しない)。plan tooling-nozzle-card-local-mesh §5.1 #3d。
# run dir には何も書かない (--segment は scratch の写しで回す)。出力は $O (原文)。手元の _runs/run_0003_ns_cont/judge/ に写して commit する。
# 判定区間は延長の本段だけ (stages=none、設定が違うので run_0002_ns とは連結しない; plan #3d)。
set -uo pipefail
C=$HOME/forge-cardmesh/case/71.m55_combustion_wt
R=$C/run_0003_ns_cont
R2=$C/run_0002_ns
T=$HOME/forge-cardmesh/solver_density_cuda/tools
DS=$HOME/forge-cardmesh/design
W=$HOME/m55_diag_0003
O=$W/out
STEP=${STEP:-200000}
mkdir -p "$O"
S=$(dirname "$(readlink -f "$0")")
echo "# $(date -Iseconds)  git $(git -C $HOME/forge-cardmesh rev-parse HEAD)  python3 $(python3 -V 2>&1)  step $STEP" > "$O/_env.txt"

# 1a. --segment (記録された stage_manifest.json のまま、写しで)
A1=$W/seg/run_0003_ns_cont; rm -rf "$W/seg"; mkdir -p "$A1"
cp -p "$R/stage_manifest.json" "$R/forge_launches.jsonl" "$A1/" 2>/dev/null
for f in "$R"/residual_history*.csv; do case "$f" in *_segment.csv) ;; *) ln -s "$f" "$A1/$(basename "$f")" ;; esac; done
{ python3 "$T/stage_manifest.py" "$A1" --segments; echo "STAGE_MANIFEST_EXIT: $?"; } > "$O/stage_manifest_segments.txt" 2>&1
{ python3 "$T/check_convergence.py" --segment "$A1"; echo "CHECK_CONVERGENCE_EXIT: $?"; } > "$O/check_convergence_segment.txt" 2>&1
# 1b. 履歴そのまま (run_case.sh の CONVERGENCE_VERDICT.txt と同じ入力)
{ python3 "$T/check_convergence.py" "$R"; echo "CHECK_CONVERGENCE_EXIT: $?"; } > "$O/check_convergence_main_only.txt" 2>&1

# 2a. 停止規則の系列 (stop_rule.json の rows) の CSV と準定常 (tail は既定)
python3 - "$R/stop_rule.json" "$O/stop_rule_series.csv" <<'PY'
import json, sys
d = json.load(open(sys.argv[1]))
ks = ["step", "r", "mdot_in", "mdot_out", "theta_r_0.4", "theta_r_0.7", "theta_r_0.95", "exit_core_M", "mdot_rt"]
with open(sys.argv[2], "w") as fh:
    fh.write(",".join(ks) + "\n")
    for r in d["rows"]:
        fh.write(",".join(repr(r.get(k)) for k in ks) + "\n")
PY
{ python3 "$T/check_quasisteady.py" --series-csv "$O/stop_rule_series.csv" \
    --series-cols r,mdot_in,mdot_out,theta_r_0.4,theta_r_0.7,theta_r_0.95,exit_core_M; echo "CHECK_QUASISTEADY_EXIT: $?"; } \
    > "$O/check_quasisteady_series.txt" 2>&1

# 2b. 最後の場の非有限・物理性
{ python3 "$S/../run_0002_ns/final_field_check.py" --tools "$T" --run "$R" --step "$STEP" --out-json "$O/final_field_check.json"; echo "EXIT: $?"; } \
    > "$O/final_field_check.stdout.txt" 2>&1

# 2c. 壁解像 (no-slip 壁は wall (physID 3, wall_isothermal) だけ)
ls "$R"/res_wall_3_$STEP.h5 > "$O/check_wall_resolution.txt" 2>&1
{ python3 "$T/check_wall_resolution.py" "$R" --step "$STEP" --groups wall --weight area; echo "CHECK_WALL_RESOLUTION_EXIT: $?"; } \
    >> "$O/check_wall_resolution.txt" 2>&1

# 3. 監視の量 (S の h0・P0、壁の積分熱流、入口・出口のエネルギー流束、入口の列の静温) の時系列と準定常
STEPS=2000$(for s in $(seq 20000 20000 $STEP); do printf ",%d" $s; done)
{ python3 "$S/monitor_0003.py" --tools "$T" --design "$DS" --run "$R" --ref-run "$R2" --steps "$STEPS" \
    --out-csv "$O/monitor_series.csv" --out-json "$O/monitor_series.json"; echo "EXIT: $?"; } > "$O/monitor_series.stdout.txt" 2>&1
{ python3 "$T/check_quasisteady.py" --series-csv "$O/monitor_series.csv" \
    --series-cols n_P0_exceed,max_dP0_over_Pt,S_mean_dP0_over_Pt,S_mean_h0_minus_h0in,E_in,E_out,Q_wall,bal_plus,bal_minus,inlet_max_T_minus_Tt; \
  echo "CHECK_QUASISTEADY_EXIT: $?"; } > "$O/check_quasisteady_monitor.txt" 2>&1

# 4. 容量
{ date -Iseconds; df -h "$HOME"; du -sb "$R"; stat -c '%s %n' "$R/res_$STEP.h5" "$R/residual_history.csv" "$R/forge_run.log"
  echo "res_N.h5 の数 $(ls "$R" | grep -cE '^res_[0-9]+\.h5$')"; } > "$O/disk.txt" 2>&1
echo DONE > "$O/_done.txt"
