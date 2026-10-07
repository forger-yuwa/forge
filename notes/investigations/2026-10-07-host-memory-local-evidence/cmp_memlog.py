# base と R1+R2 の FORGE_MEMLOG を工程別に並べ、2 格子 (s070/s085) の傾き (B/節点) を出す (計測専用)
import re, sys
SP = "/tmp/claude-1000/-home-sano-work-forge/60a080db-fad5-4e19-bdaa-3590a1035176/scratchpad/hostmem"
import os
grids = ["s070", "s085"]; tags = ["base", os.environ.get("NEWTAG", "new")]
def load(d):
    txt = open(f"{d}/forge_run.log").read()
    n = int(re.search(r"Number of Cells: (\d+)", txt).group(1))
    rows = {}; order = []; seen = {}
    for m in re.finditer(r"\[memlog\] (.+?) @(\S+) VmRSS=(\d+)MB VmHWM=(\d+)MB(?: \| (.*))?", txt):
        lab = m.group(1); seen[lab] = seen.get(lab, 0) + 1
        key = lab if seen[lab] == 1 else f"{lab}#{seen[lab]}"
        ex = m.group(5) or ""
        g = re.search(r"GPU used=(\d+)MB \(since first memlog \+(-?\d+)MB\)", ex)
        items = dict((a, float(b)) for a, b in re.findall(r"(\S+?)\[\d+\]=([\d.e+-]+)MB", ex))
        rows[key] = (int(m.group(3)), int(m.group(4)), int(g.group(2)) if g else None, items); order.append(key)
    return n, rows, order
D = {(g, t): load(f"{SP}/mem_{g}_{t}") for g in grids for t in tags}
N = {g: D[(g, "base")][0] for g in grids}
assert all(D[(g, tags[1])][0] == N[g] for g in grids)
print("nodes:", N)
def slope(y):  # 2 点の傾き B/節点 と切片 MB
    a = (y[1] - y[0]) * 1048576 / (N[grids[1]] - N[grids[0]]); b = y[0] - a * N[grids[0]] / 1048576; return a, b
order = D[("s085", "base")][2]
want = [k for k in order if "#" not in k or k.startswith("after step 1")]
print(f"{'stage':52s} | HWM MB base s070/s085 | new s070/s085 | HWM slope base -> new (B/node) | dslope | RSS slope base -> new | GPUd MB base/new s085")
for k in want:
    if not all(k in D[(g, t)][1] for g in grids for t in tags): print(f"{k:52s} | (missing)"); continue
    hb = [D[(g, "base")][1][k][1] for g in grids]; hn = [D[(g, tags[1])][1][k][1] for g in grids]
    rb = [D[(g, "base")][1][k][0] for g in grids]; rn = [D[(g, tags[1])][1][k][0] for g in grids]
    gb = D[("s085", "base")][1][k][2]; gn = D[("s085", tags[1])][1][k][2]
    sb, _ = slope(hb); sn, _ = slope(hn); qb, _ = slope(rb); qn, _ = slope(rn)
    print(f"{k:52s} | {hb[0]:5d} {hb[1]:5d} | {hn[0]:5d} {hn[1]:5d} | {sb:6.0f} -> {sn:6.0f} | {sn-sb:+6.0f} | {qb:6.0f} -> {qn:6.0f} | {gb} / {gn}")
print()
k = "end of initializeSimulation"
print("containers at", k)
for name in D[("s085", "base")][1][k][3]:
    yb = [D[(g, "base")][1][k][3].get(name, 0.0) for g in grids]; yn = [D[(g, tags[1])][1][k][3].get(name, 0.0) for g in grids]
    print(f"  {name:28s} base {yb[0]:8.1f} {yb[1]:8.1f} MB ({slope(yb)[0]:7.1f} B/node)  new {yn[0]:8.1f} {yn[1]:8.1f} MB ({slope(yn)[0]:7.1f} B/node)")
for g in grids:
    for t in tags:
        tv = open(f"{SP}/mem_{g}_{t}/time_v.txt").read()
        m = re.search(r"Maximum resident set size \(kbytes\): (\d+)", tv)
        print(f"  /usr/bin/time maxRSS {g} {t}: {int(m.group(1))/1024:.0f} MB")
