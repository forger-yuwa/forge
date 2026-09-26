#!/usr/bin/env python3
"""設計 DB の学習採否ゲート (plan gradient-scalar-lsq-unification #6、codex diagnose 2026-09-27 の A/B)。

forge は回さない。node・新 policy・数値ゲート PASS を固定し、作動点の実効 scalarGradient だけを変える。
合格 (測る前に固定): gg → 0 件、lsq → 1 件、不明 (None) → 0 件、2 作動点のうち 1 つだけ gg → 0 件、旧 policy の lsq → 0 件。
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import types  # noqa: E402
try:
    import pymoo  # noqa: F401
except ImportError:   # 選別関数だけを試すので、最適化モジュールは空のスタブで足りる
    _m = types.ModuleType("forge_design.opt.moo"); _m.propose_infill = None
    sys.modules["forge_design.opt.moo"] = _m
try:
    import smt  # noqa: F401
except ImportError:
    _s = types.ModuleType("forge_design.opt.surrogate"); _s.KrigingSet = object
    sys.modules["forge_design.opt.surrogate"] = _s
from forge_design.opt.driver_sern import FLAG_POLICY, SernCampaign  # noqa: E402


class _Fake:
    def __init__(self, rows):
        self.rows = rows


def row(sg_list, policy=FLAG_POLICY):
    return {"status": "PASS", "flag_policy": policy, "x": [0.0], "C_T_w": 0.9, "L_ramp": 1.0,
            "ops": {f"op{i}": {"scalar_gradient_effective": sg} for i, sg in enumerate(sg_list)}}


fails = 0
def check(name, rows, want):
    global fails
    X, _ = SernCampaign._XF(_Fake(rows))
    got = len(X)
    ok = got == want
    fails += not ok
    print(("PASS " if ok else "FAIL ") + f"{name}: 採用 {got} 件 (期待 {want})")


check("gg", [row(["gg"])], 0)
check("lsq", [row(["lsq"])], 1)
check("不明 (None)", [row([None])], 0)
check("2 作動点の片方だけ gg", [row(["lsq", "gg"])], 0)
check("旧 policy の lsq", [row(["lsq"], policy="2026-09-26")], 0)
print("VERDICT:", "PASS" if fails == 0 else f"FAIL ({fails})")
sys.exit(1 if fails else 0)
