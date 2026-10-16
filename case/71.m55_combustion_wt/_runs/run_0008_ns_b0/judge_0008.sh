#!/bin/bash
# case/71 run_0008_ns_b0 (2 係数の壁の NS、本段 B0; judge_0007.sh の写しで run と出力先だけ違う) の判定ツールを回す (AWS A、forge は起動しない)。plan tooling-nozzle-card-local-mesh §5.1 #3o・§6 V5。
# 判定区間は stage_manifest.json の最後の区間 (= 本段; `check_convergence --segment`)。run dir には何も書かない。出力は $O。
set -uo pipefail
C=$HOME/forge-cardmesh/case/71.m55_combustion_wt
R=$C/run_0008_ns_b0
T=$HOME/forge-cardmesh/solver_density_cuda/tools
W=$HOME/m55_diag_0008
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
    if not isinstance(v, dict):
        print(k, v); continue
    print(k, json.dumps({kk: vv for kk, vv in (v or {}).items() if kk in ("pass", "kind", "measured", "r_mean", "r_amp", "worst", "why")},
                        ensure_ascii=False, default=str)[:1500])
PY
{ date -Iseconds; df -h "$HOME"; du -sb "$R"; } > "$O/disk.txt" 2>&1
# 6. V5b 出口コア M (停止規則の量、最後の出力) と V5c 試験部の一様性 (全場の出力ごと → 両量の時系列の準定常)
PY=$HOME/forge-deliv/design/.venv-opt/bin/python
python3 - "$R/stop_rule.json" > "$O/v5b_exit_core_M.txt" 2>&1 <<'PYB'
import json, sys
r = json.load(open(sys.argv[1]))["rows"][-1]; M = r["exit_core_M"]; lo, hi = 5.5 * (1 - 5e-4), 5.5 * (1 + 5e-4)
print(f"V5b: step {r['step']} 出口コア M {M:.6f} (許容 {lo:.5f}〜{hi:.5f}) → {'PASS' if lo <= M <= hi else 'FAIL'}")
PYB
( cd $HOME/forge-cardmesh/design && PYTHONPATH=$PWD $PY "$C/_runs/run_0007_ns/v5c_series.py" "$R" "$O/v5c_series.csv" ) > "$O/v5c_series.txt" 2>&1
{ python3 "$T/check_quasisteady.py" --series-csv "$O/v5c_series.csv" --series-cols maxabs_eta0.1,halfamp_eta0.1,maxabs_axis,halfamp_axis; echo "CHECK_QUASISTEADY_EXIT: $?"; } \
    > "$O/check_quasisteady_v5c.txt" 2>&1
echo DONE > "$O/_done.txt"
