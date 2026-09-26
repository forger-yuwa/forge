#!/usr/bin/env python3
"""stage_manifest の `mesh.scalarGradient` 区間識別と YAML 解析不能時の扱い
(plan gradient-scalar-lsq-unification §5.1 #5j、codex result-1 M5)。

合格 (測る前に固定): 明示 gg vs 明示 lsq → 2 区間 / 明示 lsq 同士 → 1 区間 / 省略 vs 明示 gg → 2 区間 (暫定、分けすぎる側) /
解析不能な YAML 2 段 (本文が違う) → 2 区間。
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stage_manifest import segments, stage_key  # noqa: E402

BC = "wall: {kind: wall}\n"
CFG = lambda grad: ('solver: "SLAU"\nmesh: {discretization: "node", nodeWallDirichlet: 1%s}\n'
                    'space: {convMethod: 1, limiter: 2}\n') % (f", scalarGradient: {grad}" if grad else "")


def segs(*cfgs):
    man = {"stages": [{"tag": str(i), "key": stage_key(c, BC)} for i, c in enumerate(cfgs)]}
    return len(segments(man))


fails = 0


def check(name, cond):
    global fails
    print(("PASS " if cond else "FAIL ") + name)
    fails += (not cond)


check("明示 gg vs 明示 lsq → 2 区間", segs(CFG("gg"), CFG("lsq")) == 2)
check("明示 lsq 同士 → 1 区間", segs(CFG("lsq"), CFG("lsq")) == 1)
check("省略 vs 明示 gg → 2 区間 (暫定: 分けすぎる側)", segs(CFG(None), CFG("gg")) == 2)
bad1 = CFG("gg") + "space: [unclosed\n"
bad2 = CFG("lsq") + "space: [unclosed\n"
check("解析不能な YAML (本文が違う) 2 段 → 2 区間", segs(bad1, bad2) == 2)
print("VERDICT:", "PASS" if fails == 0 else "FAIL (%d)" % fails)
sys.exit(1 if fails else 0)
