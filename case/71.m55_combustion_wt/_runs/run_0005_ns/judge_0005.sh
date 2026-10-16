#!/bin/bash
# case/71 run_0005_ns (δ* の補正の 1 回目の NS) の判定ツールを最後の場で回す (AWS A、forge は起動しない)。plan tooling-nozzle-card-local-mesh §5.1 #3e (6)。
# 判定区間は stage_manifest.json の最後の区間 (= 本段; `check_convergence --segment`)。run dir には何も書かない。出力は $O。
set -uo pipefail
C=$HOME/forge-cardmesh/case/71.m55_combustion_wt
R=$C/run_0005_ns
T=$HOME/forge-cardmesh/solver_density_cuda/tools
W=$HOME/m55_diag_0005
O=$W/out
STEP=${STEP:?STEP (最後の出力の step) を渡す}
mkdir -p "$O"
echo "# $(date -Iseconds)  git $(git -C $HOME/forge-cardmesh rev-parse HEAD)  step $STEP" > "$O/_env.txt"
# 1. 収束 (本段の区間)
{ python3 "$T/stage_manifest.py" "$R" --segments; echo "STAGE_MANIFEST_EXIT: $?"; } > "$O/stage_manifest_segments.txt" 2>&1
{ python3 "$T/check_convergence.py" --segment "$R"; echo "CHECK_CONVERGENCE_EXIT: $?"; } > "$O/check_convergence_segment.txt" 2>&1
# 2. 停止規則の系列の準定常
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
# 3. 最後の場の非有限・物理性 (run_0002 の台本)
{ python3 "$C/_runs/run_0002_ns/final_field_check.py" --tools "$T" --run "$R" --step "$STEP" --out-json "$O/final_field_check.json"; echo "EXIT: $?"; } \
    > "$O/final_field_check.stdout.txt" 2>&1
# 4. 壁解像 (no-slip 壁は wall だけ; 面積の重み)
{ python3 "$T/check_wall_resolution.py" "$R" --step "$STEP" --groups wall --weight area; echo "CHECK_WALL_RESOLUTION_EXIT: $?"; } \
    > "$O/check_wall_resolution.txt" 2>&1
# 5. 関門 (_after_ns が書いた HANDOFF_GATE.json の要点)
python3 - "$R/HANDOFF_GATE.json" > "$O/handoff_gate_summary.txt" 2>&1 <<'PY'
import json, sys
g = json.load(open(sys.argv[1]))
print(json.dumps({k: g.get(k) for k in ("proceed", "blocked_by", "settled")}, ensure_ascii=False, indent=1))
for k, v in (g.get("checks") or {}).items():
    print(k, json.dumps({kk: vv for kk, vv in (v or {}).items() if kk in ("pass", "kind", "measured", "r_mean", "r_amp", "worst", "why")},
                        ensure_ascii=False, default=str)[:1500])
PY
{ date -Iseconds; df -h "$HOME"; du -sb "$R"; } > "$O/disk.txt" 2>&1
echo DONE > "$O/_done.txt"
