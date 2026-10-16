#!/bin/bash
# plan tooling-settle-convergence-and-monitor §6 V3(c) の判定 (事前登録 2026-10-16): run_0017 を v3 の基準で再生する。
# 道具は ~/forge-cardmesh-tools (e92f6b1e = #9〜#11 の修正の後)。基準は run_0015 と run_0017 の両方から同じ作り方で作り直して比べる。
set -uo pipefail
C=$HOME/forge-cardmesh/case/71.m55_combustion_wt
T=$HOME/forge-cardmesh-tools
PY=$HOME/venv-cardmesh/bin/python
W=$HOME/m55_p17/replay; mkdir -p $W
[ "$(git -C $T rev-parse --short=8 HEAD)" = e92f6b1e ] || { echo "STOP: 道具が e92f6b1e でない"; exit 1; }
for R in run_0015_ns_float_a run_0017_ns_float_free130; do
  ( cd $T/design && PYTHONPATH=$PWD $PY - "$C/$R" "$W/criteria_v3_$R.yaml" <<'PYC'
import sys, json, yaml, copy
from pathlib import Path
from forge_design.evaluate import settle_presets as SP
run, out = Path(sys.argv[1]), Path(sys.argv[2])
cur = yaml.safe_load((run / "settle_criteria.yaml").read_text())
ev = copy.deepcopy(cur["source"]["evaluate_settle"]); ev["preset"] = "wt_nozzle_v3"
info = json.loads((run / "prepare_info.json").read_text())
cfg = yaml.safe_load((run / "solverConfig.yaml").read_text()); bc = yaml.safe_load((run / "bcondConfig.yaml").read_text())
crit = SP.resolve(ev, info, cfg, bc, segment=SP.main_stage_segment())
out.write_text(yaml.safe_dump(crit, allow_unicode=True, sort_keys=False))
print("wrote", out, crit["source"].get("preset_version"))
PYC
  ) > $W/make_$R.log 2>&1 || { echo "STOP: 基準を作れない $R"; cat $W/make_$R.log; exit 1; }
done
$PY - "$W/criteria_v3_run_0015_ns_float_a.yaml" "$W/criteria_v3_run_0017_ns_float_free130.yaml" "$HOME/m55_p5/v3b_replay/criteria_v2_run_0015_ns_float_a.yaml" <<'PYC' | tee $W/criteria_compare.txt
import sys, yaml, hashlib
a, b, old = (yaml.safe_load(open(f)) for f in sys.argv[1:4])
def strip(c):
    c = dict(c); c.pop("source", None); return c
ok1 = strip(a) == strip(b)
print("run_0015 と run_0017 の基準 (source 以外):", "一致" if ok1 else "違う")
o = strip(old)
for q in o["quantities"]:
    pass
a2 = strip(a)
for q in a2["quantities"]:
    if q.get("kind") == "residual_finite":
        q = q
diffq = []
qa = {q["name"]: q for q in a2["quantities"]}; qo = {q["name"]: q for q in o["quantities"]}
for n in sorted(set(qa) | set(qo)):
    x, y = dict(qa.get(n, {})), dict(qo.get(n, {}))
    if n == "residual_finite":
        x.pop("required", None); x.pop("notes", None); y.pop("notes", None)
    if x != y:
        diffq.append(n)
rest = {k for k in set(a2) | set(o) if k != "quantities" and a2.get(k) != o.get(k)}
ok2 = not diffq and not rest and "required" in qa.get("residual_finite", {})
print("旧の run_0015 の v3b の基準との差が residual_finite の required (と notes) だけ:", "はい" if ok2 else f"いいえ (量 {diffq}、その他 {sorted(rest)})")
for f in sys.argv[1:4]:
    print(hashlib.sha256(open(f, "rb").read()).hexdigest(), f)
sys.exit(0 if (ok1 and ok2) else 1)
PYC
[ "${PIPESTATUS[0]}" = 0 ] || { echo "STOP: 基準の照合が事前登録と違う"; exit 1; }
( time $PY $T/solver_density_cuda/tools/settle_virtual_stop.py $C/run_0017_ns_float_free130 --criteria $W/criteria_v3_run_0017_ns_float_free130.yaml \
    --out $W/run_0017_ns_float_free130 --span 80000 --ts-max 50000 ) > $W/replay_run_0017.log 2>&1
echo "replay rc $?"
cut -c1-400 $W/run_0017_ns_float_free130/VERDICT.txt
echo REPLAY_DONE
