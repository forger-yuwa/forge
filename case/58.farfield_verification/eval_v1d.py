#!/usr/bin/env python3
"""V1(d) 判定 (plan boundary-node-farfield-characteristic §6 V1(d)): 一様状態 (SST 込み) での輸送残差の打消し。
帳簿ダンプ (FORGE_DUMP_LEDGER、全節点に印、初回 assembleResidual) と farfield の面ダンプ (FORGE_DUMP_FARFIELD) から、
  対流の寄与     = res_after_conv − res_before_conv            (ρ, ρu, ρE)
  k・ω 輸送の寄与 = res_after_rans_transport − res_after_conv  (ρk, ρω)
を全節点で出し、保存量ごとの規模 (その節点に接する面流束の絶対和) で割った値が 1e-5 以下を合格とする。
規模: ρ・ρE は各面流束の絶対和、運動量 3 成分は共通で Σ(|F_ρux|+|F_ρuy|+|F_ρuz|) (pRef = P∞ のゲージでは接線成分の
面流束が 0 になり成分別では割れないため)、k は Σ|ṁ| k∞、ω は Σ|ṁ| ω∞。
  python3 eval_v1d.py RUN_DIR
"""
import csv, glob, sys
from collections import defaultdict
import yaml

run = sys.argv[1]
fl = yaml.safe_load(open(run + "/bcondConfig.yaml"))["xmin"]["floats"]
kinf, oinf = float(fl["k"]), float(fl["omega"])

res = defaultdict(dict)   # (tag, field) -> {node: value}
for r in csv.DictReader(open(run + "/ledger.csv")):
    if r["call"] == "1" and r["tag"].startswith("res_"):
        res[(r["tag"], r["field"])][int(r["node"])] = float(r["value"])

MOM = ("F_roUx", "F_roUy", "F_roUz")
sc = defaultdict(lambda: defaultdict(float))  # node -> quantity -> scale
nffaces = 0
def add(node, F):
    s = sc[node]
    s["ro"] += abs(F["F_ro"]); s["roe"] += abs(F["F_roe"])
    s["mom"] += sum(abs(F[m]) for m in MOM)
    s["roK"] += abs(F["F_ro"]) * kinf; s["roOmega"] += abs(F["F_ro"]) * oinf
seen = set()
for r in csv.DictReader(open(run + "/ledger.csv.faces")):
    if r["call"] != "1" or r["ip"] in seen:
        continue
    seen.add(r["ip"])
    F = {k: float(r[k]) for k in ("F_ro", "F_roe") + MOM}
    add(int(r["ic0"]), F); add(int(r["ic1"]), F)
for p in glob.glob(run + "/ffdump.*.csv"):
    for r in csv.DictReader(open(p)):
        add(int(float(r["ic"])), {k: float(r[k]) for k in ("F_ro", "F_roe") + MOM})
        nffaces += 1
        if int(float(r["vacuum"])) or int(float(r["hll"])):
            print(f"  置換/退避あり: {p} ip {r['ip']}")

def contrib(field, t0, t1):
    a, b = res[(t1, field)], res[(t0, field)]
    return {n: a[n] - b[n] for n in a}

rows = [("ro", "ro", "res_before_conv", "res_after_conv"), ("roUx", "mom", "res_before_conv", "res_after_conv"),
        ("roUy", "mom", "res_before_conv", "res_after_conv"), ("roUz", "mom", "res_before_conv", "res_after_conv"),
        ("roe", "roe", "res_before_conv", "res_after_conv"),
        ("roK", "roK", "res_after_conv", "res_after_rans_transport"), ("roOmega", "roOmega", "res_after_conv", "res_after_rans_transport")]
worst_all = 0.0
print(f"{run}: 内部面 {len(seen)}、farfield 面 {nffaces}、節点 {len(sc)}")
for fld, skey, t0, t1 in rows:
    d = contrib("res_" + fld, t0, t1)
    w, wn = 0.0, -1
    for n, v in d.items():
        r = abs(v) / sc[n][skey] if sc[n][skey] > 0 else float("inf")
        if r > w:
            w, wn = r, n
    worst_all = max(worst_all, w)
    print(f"  {fld:8s}: 最大 |寄与|/規模 {w:.3e} (節点 {wn})")
print(f"VERDICT: {'PASS' if worst_all <= 1e-5 else 'FAIL'} (最大 {worst_all:.3e}、許容 1e-5)")
