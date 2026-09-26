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
from forge_design.evaluate import runner_sern as RS  # noqa: E402
import json, tempfile  # noqa: E402,E401


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

# 起動記録 → 実効値 (codex result 2026-09-27 M1): 最後の起動が読めなければ不明。前の lsq を引き継がない
def eff(lines):
    d = tempfile.mkdtemp()
    open(os.path.join(d, "forge_launches.jsonl"), "w").write("\n".join(lines) + "\n")
    return RS._last_launch_value(d, "scalarGradient", allowed=("gg", "lsq"))
L = lambda sg: json.dumps({"cfg_fnv": "x", "scalarGradient": sg})
for name, lines, want in (("lsq", [L("lsq")], "lsq"), ("gg→lsq", [L("gg"), L("lsq")], "lsq"),
                          ("lsq→キー欠落", [L("lsq"), json.dumps({"cfg_fnv": "x"})], None),
                          ("lsq→JSON 破損", [L("lsq"), "{broken"], None), ("未知値", [L("foo")], None)):
    got = eff(lines); ok = got == want; fails += not ok
    print(("PASS " if ok else "FAIL ") + f"起動記録 {name}: 実効値 {got!r} (期待 {want!r})")
    ok2 = SernCampaign._XF(_Fake([row([got])]))[0].shape[0] == (1 if want == "lsq" else 0); fails += not ok2
    print(("PASS " if ok2 else "FAIL ") + f"起動記録 {name}: 学習採否")

# Pareto の母集団 = 学習の母集団 (M2)
class _FakeS(_Fake):
    ref = (-0.90, 20.0); ops = []
rows = [dict(row(["lsq"]), tag="a", C_T_w=0.95, C_M_w=0.0), dict(row(["gg"]), tag="b", C_T_w=0.99, C_M_w=0.0),
        dict(row(["lsq"], policy="2026-09-26"), tag="c", C_T_w=1.00, C_M_w=0.0)]
summ = SernCampaign.summary(_FakeS(rows))
ok = summ["n_pass"] == 1 and [p["tag"] for p in summ["pareto"]] == ["a"] and summ["n_pass_excluded_by_policy"] == 2
fails += not ok
print(("PASS " if ok else "FAIL ") + f"Pareto: n_pass {summ['n_pass']} / 除外 {summ['n_pass_excluded_by_policy']} / Pareto {[p['tag'] for p in summ['pareto']]} (期待 1 / 2 / ['a'])")
print("VERDICT:", "PASS" if fails == 0 else f"FAIL ({fails})")
sys.exit(1 if fails else 0)
