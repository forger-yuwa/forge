# forge の FORGE_MEMLOG 出力を工程別の表にし、3 格子から B/節点 (傾き) と定数を最小二乗で出す (計測専用)
import re,sys,numpy as np
runs={"s050":"mem_s050","s070":"mem_s070","s085":"mem_s085"}
data={}; N={}
for t,d in runs.items():
    txt=open(f"{d}/forge_run.log").read()
    N[t]=int(re.search(r"Number of Cells: (\d+)",txt).group(1))
    rows=[]; seen={}
    for m in re.finditer(r"\[memlog\] (.+?) @(\S+) VmRSS=(\d+)MB VmHWM=(\d+)MB(?: \| (.*))?",txt):
        lab=m.group(1); seen[lab]=seen.get(lab,0)+1
        key=lab if seen[lab]==1 else f"{lab}#{seen[lab]}"
        ex=m.group(5) or ""
        g=re.search(r"GPU used=(\d+)MB \(since first memlog \+(-?\d+)MB\)",ex)
        items=dict((a,float(b)) for a,b in re.findall(r"(\S+?)\[\d+\]=([\d.e+-]+)MB",ex))
        rows.append((key,int(m.group(3)),int(m.group(4)),int(g.group(2)) if g else None,int(g.group(1)) if g else None,items))
    data[t]=rows
keys=[r[0] for r in data["s085"]]
ns=np.array([N[t] for t in runs])
print("nodes:",N)
print(f"{'stage':58s} | RSS MB s050/s070/s085 | HWM MB s050/s070/s085 | GPUΔ MB | RSS fit B/node + MB | HWM fit | GPUΔ fit")
def fit(y):
    A=np.vstack([ns,np.ones_like(ns)]).T; c,res,_,_=np.linalg.lstsq(A,np.array(y,float)*1048576,rcond=None); return c[0],c[1]/1048576
for k in keys:
    vals=[]
    for t in runs:
        r=[x for x in data[t] if x[0]==k]
        vals.append(r[0] if r else None)
    if any(v is None for v in vals): 
        print(f"{k:58s} | (missing in some run)"); continue
    rss=[v[1] for v in vals]; hwm=[v[2] for v in vals]; g=[v[3] for v in vals]
    a,b=fit(rss); c,d=fit(hwm)
    gs="    -     -     -" if g[0] is None else f"{g[0]:5d} {g[1]:5d} {g[2]:5d}"
    gf="" if g[0] is None else "%6.0f +%5.0f"%fit(g)
    print(f"{k:58s} | {rss[0]:5d} {rss[1]:5d} {rss[2]:5d} | {hwm[0]:5d} {hwm[1]:5d} {hwm[2]:5d} | {gs} | {a:6.0f} +{b:5.0f} | {c:6.0f} +{d:5.0f} | {gf}")
print()
# コンテナ別 (end 時点、無ければ最後)
last={t:data[t][-1][5] for t in runs}
for name in last["s085"]:
    y=[last[t].get(name,0.0) for t in runs]; a,b=fit(y)
    print(f"{name:40s} {y[0]:9.1f} {y[1]:9.1f} {y[2]:9.1f} MB  -> {a:7.1f} B/node + {b:6.1f} MB")
print("\n== containers at 'end of initializeSimulation'")
tot=np.zeros(3)
for name in [x for x in data["s085"] if x[0]=="end of initializeSimulation"][0][5]:
    y=[[x for x in data[t] if x[0]=="end of initializeSimulation"][0][5].get(name,0.0) for t in runs]; a,b=fit(y)
    if "_d(" not in name: tot+=np.array(y)
    print(f"{name:40s} {y[0]:9.1f} {y[1]:9.1f} {y[2]:9.1f} MB  -> {a:7.1f} B/node + {b:6.1f} MB")
a,b=fit(tot); print(f"{'host containers total':40s} {tot[0]:9.1f} {tot[1]:9.1f} {tot[2]:9.1f} MB  -> {a:7.1f} B/node + {b:6.1f} MB")
