#!/usr/bin/env python3
"""正式な入口の回し直し (`wt_qualify.approval_rechecks`) の試験 (plan tooling-design-orchestrator §4.7.1 (6)、§5.1 #2c)。
run の側で作られた判定 (停止の判定・計算の格子の品質) を、判定の側のコードで回し直す。手元の case/71 run_0014_chain_v2 を読む (書かない)。

  (P) 停止の step 46171・格子の AR ≤ 5000 → 回し直しに問題なし (停止の判定 SETTLED・格子の品質 PASS)
  (N1) 停止の step を 20000 にした流れの記録 → 停止の判定を回し直すと SETTLED でない
  (N2) 承認の格子の AR の上限 10 → 格子の品質が PASS でない
  (N3) 停止の step が無い → 回し直せない

実行: /home/sano/work/forge/design/.venv-opt/bin/python design/tests/run_approval_rechecks_tests.py
所要: 1〜2 分 (check_settled を 2 回)。run が無ければ飛ばす (rc 0 と「飛ばした」)。
"""
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from forge_design.feedback import wt_qualify as Q  # noqa: E402

R14 = Path(os.environ.get("FORGE_CASE71_RUN_0014", Path(__file__).resolve().parents[2] / "case" / "71.m55_combustion_wt" / "run_0014_chain_v2"))
FAIL = 0


def check(name, cond, info=""):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name + (f"  [{info}]" if info else ""))
    if not cond:
        FAIL += 1


if not (R14 / "SETTLE_VERDICT.json").is_file() or not (R14 / "nozzle.h5").is_file():
    print(f"飛ばした: {R14} に run が無い")
    sys.exit(0)
APR = {"protected": {"evaluation": {"mesh_quality": {"ar_max": 5000.0, "skew_max": 0.9}}}}
r = Q.approval_rechecks(R14, {"stop_step": 46171}, APR)
check("(P) 停止の step 46171・AR ≤ 5000 → 問題なし", not r["problems"], f"{r.get('settle', {}).get('verdict_line')} / {r.get('mesh_quality', {}).get('verdict_line')}")
r = Q.approval_rechecks(R14, {"stop_step": 20000}, APR)
check("(N1) 停止の step 20000 → 停止の判定が SETTLED でない", any("SETTLED でない" in p for p in r["problems"]), f"{r.get('settle', {}).get('verdict_line')}")
r = Q.approval_rechecks(R14, {"stop_step": 46171}, {"protected": {"evaluation": {"mesh_quality": {"ar_max": 10.0, "skew_max": 0.9}}}})
check("(N2) AR の上限 10 → 格子の品質が PASS でない", any("品質" in p for p in r["problems"]), f"{r.get('mesh_quality', {}).get('verdict_line')}")
r = Q.approval_rechecks(R14, {}, APR)
check("(N3) 停止の step が無い → 回し直せない", any("停止の step が無い" in p for p in r["problems"]))
print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAIL'}")
sys.exit(1 if FAIL else 0)
