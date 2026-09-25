#!/usr/bin/env python3
"""stage_manifest の `space.slauWallNormalChi` 区間識別 (plan convection-slau-wall-normal-chi §5.1 #12, codex result-3 M5)。

合格 (測る前に固定): 省略 vs 0 → 同一 key・1 区間 / 0 vs 1 → 2 区間 / flow 形式と block 形式で同一 key。
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stage_manifest import stage_key, segments, fnv1a64, CHI_KEY as CHI

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

# ---- 既定化後 (plan convection-slau-wall-normal-chi-default §4.3): 実効値で区間を切る ----
K = lambda c: stage_key(c, BC)[CHI]
NODE = lambda extra, solv="SLAU", nwd=1: 'solver: "%s"\nmesh: {discretization: "node", nodeWallDirichlet: %d}\n' % (solv, nwd) + FLOW(extra)
CELL = lambda extra: 'solver: "SLAU"\nmesh: {discretization: "cell"}\n' + FLOW(extra)
check("node+SLAU で省略 → 実効 1", K(NODE("")) == "1")
check("node+SLAU で省略 と 明示 0 → 別区間", segs(NODE(""), NODE(", slauWallNormalChi: 0")) == 2)
check("node+SLAU で省略 と 明示 1 → 実効値は同じ (由来は inferred/explicit で別区間)", K(NODE("")) == K(NODE(", slauWallNormalChi: 1")))
check("cell で省略 → 実効 0", K(CELL("")) == "0")
check("非 SLAU (ROE) で省略 → 実効 0", K(NODE("", "ROE")) == "0")
check("nodeWallDirichlet 0 で省略 → 実効 0", K(NODE("", nwd=0)) == "0")
check("flow と block (1) → 同一 key", stage_key(FLOW(", slauWallNormalChi: 1"), BC) == stage_key(BLOCK(", slauWallNormalChi: 1"), BC))
# 旧形式 manifest の移行
old = {"stages": [{"tag": "a", "key": {"space.convMethod": "1"}}, {"tag": "b", "key": {"space.convMethod": "1"}},
                  {"tag": "c", "key": {"space.convMethod": "1", "space.slauWallNormalChi": "1"}}]}
sg = segments(old)
check("旧形式: キー無し 2 段は 0 で 1 区間、明示 1 の段は別区間", len(sg) == 2 and sg[0][0]["key"][CHI] == "0" and sg[1][0]["key"][CHI] == "1")
# 起動記録での確定: 同一 YAML でもバイナリ更新をまたいで実効値が違えば別区間
cfg = NODE("")
man = {"stages": [{"tag": "s1", "key": stage_key(cfg, BC), "cfg_fnv": "aa", "chi_source": "inferred"},
                  {"tag": "s2", "key": stage_key(cfg, BC), "cfg_fnv": "bb", "chi_source": "inferred"}]}
sg = segments(man, {"aa": "0", "bb": "1"})
check("同一 YAML・起動記録で実効 0→1 → 2 区間", len(sg) == 2)
sg = segments(man, {"aa": "1"})
check("推定の段と確定の段は値が同じでも連結しない", len(sg) == 2)
sg = segments(man, {"aa": "1", "bb": "1"})
check("両段とも確定で同値 → 1 区間", len(sg) == 1)
check("fnv1a64 が FNV-1a 64 の標準値 (空・\"a\")", fnv1a64(b"") == "cbf29ce484222325" and fnv1a64(b"a") == "af63dc4c8601ec8c")

# solver の区間分離 (plan convection-slau-wall-normal-chi-usage-rule §5.1 #5)
S = lambda solv: 'solver: "%s"\n' % solv + FLOW("")
check("SLAU vs ROE → 2 区間", segs(S("SLAU"), S("ROE")) == 2)
check("SLAU vs slau → 同一 key", stage_key(S("SLAU"), BC) == stage_key(S("slau"), BC))
check("SLAU vs SLAU2 → 2 区間", segs(S("SLAU"), S("SLAU2")) == 2)
print("VERDICT:", "PASS" if fails == 0 else "FAIL (%d)" % fails)
sys.exit(1 if fails else 0)
