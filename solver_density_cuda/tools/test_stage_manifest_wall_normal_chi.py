#!/usr/bin/env python3
"""stage_manifest の `space.slauWallNormalChi` 区間識別 (plan convection-slau-wall-normal-chi §5.1 #12, codex result-3 M5)。

合格 (測る前に固定): 省略 vs 0 → 同一 key・1 区間 / 0 vs 1 → 2 区間 / flow 形式と block 形式で同一 key。
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stage_manifest import stage_key, segments

BASE = "time: {{nStepInner: 5}}\n{space}\n"
FLOW = lambda extra: BASE.format(space="space: {convMethod: 1, limiter: 2%s}" % extra)
BLOCK = lambda extra: BASE.format(space="space:\n  convMethod: 1\n  limiter: 2\n" + (("  " + extra.strip(", ") + "\n") if extra else ""))
BC = "wall: {kind: wall}\n"

def segs(*cfgs):
    man = {"stages": [{"tag": str(i), "key": stage_key(c, BC)} for i, c in enumerate(cfgs)]}
    return len(segments(man))

fails = 0
def check(name, cond):
    global fails
    print(("PASS " if cond else "FAIL ") + name); fails += (not cond)

check("省略 vs 0 (flow) → 同一 key", stage_key(FLOW(""), BC) == stage_key(FLOW(", slauWallNormalChi: 0"), BC))
check("省略 vs 0 → 1 区間", segs(FLOW(""), FLOW(", slauWallNormalChi: 0")) == 1)
check("0 vs 1 → 2 区間", segs(FLOW(", slauWallNormalChi: 0"), FLOW(", slauWallNormalChi: 1")) == 2)
check("省略 vs 1 → 2 区間", segs(FLOW(""), FLOW(", slauWallNormalChi: 1")) == 2)
check("flow と block (1) → 同一 key", stage_key(FLOW(", slauWallNormalChi: 1"), BC) == stage_key(BLOCK(", slauWallNormalChi: 1"), BC))
check("flow と block (省略) → 同一 key", stage_key(FLOW(""), BC) == stage_key(BLOCK(""), BC))
check("キー無しの旧 manifest と省略が同一 (後方互換)", "space.slauWallNormalChi" not in stage_key(FLOW(""), BC))
# solver の区間分離 (plan convection-slau-wall-normal-chi-usage-rule §5.1 #5)
S = lambda solv: 'solver: "%s"\n' % solv + FLOW("")
check("SLAU vs ROE → 2 区間", segs(S("SLAU"), S("ROE")) == 2)
check("SLAU vs slau → 同一 key", stage_key(S("SLAU"), BC) == stage_key(S("slau"), BC))
check("SLAU vs SLAU2 → 2 区間", segs(S("SLAU"), S("SLAU2")) == 2)
print("VERDICT:", "PASS" if fails == 0 else "FAIL (%d)" % fails)
sys.exit(1 if fails else 0)
