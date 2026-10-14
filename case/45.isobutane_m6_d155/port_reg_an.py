"""plan architecture-float-state-double-geometry §6.59 (§5.1 #26 の回帰、事前登録) の判定。port_reg.sh の後に B で回す。
R1 (case/45、lowMachPrecond 0) と R3 (case/48 CPG、lowMachPrecond 0): 各組 (設定・精度) の移植の前 2 本 (a・b) と後 1 本について、res_<N>.h5 の
  VALUE の全浮動小数点の量と残差の履歴の rms 列の、max|後 − 前a| ≤ max|前b − 前a| (前の 2 本がビット一致なら後もビット一致) を全量で満たせば合格。
R2 (case/45、point・FP64・lowMachPrecond 3・β = 1、1 step、stage.h5 の d): 前処理の経路の境界の修復。
  (i) 移植の後の P3 で拘束の成分 (軸の ρv、壁の ρu・ρv、強いピンの等温壁の ρE) の max|d| が 0。
  (ii) 入口・出口の境界の節点の各成分で、移植の後の P0 に対する相対差 e が、e_前 − e_後 > 10 × 再実行の幅 (P3 の前・後の 2 本ずつの差の大きい方)。
  結果 A (修復の効果を支持): (i) と (ii) をすべて満たす。結果 B: 拘束の違反が残る、または凍結が改善しない。欠損は判定不能。
  移植の後も P0 と P3 に残る差 (軸の TP のソースの Jacobian の χ_hoop、解く精度、Jacobian の作り方) は領域ごとに記録だけ。
→ _band_ab/cold_pair/port_reg_judge.json"""
import csv, json, math, re
from pathlib import Path
import numpy as np, h5py
HERE = Path(__file__).resolve().parent; C48 = HERE.parent / "48.flat_plate_cooled_m4"
OUT = {"R1": {}, "R3": {}, "R2": {}, "undecidable": []}
def und(m): OUT["undecidable"].append(m); print("  [判定不能]", m)
def vals(d, n):
    p = d / f"res_{n}.h5"
    if not p.exists(): return None
    with h5py.File(p, "r") as h: return {k: np.asarray(h["VALUE/" + k][:], np.float64) for k in h["VALUE"] if h["VALUE"][k].dtype.kind == "f"}
def rms(d):
    p = d / "residual_history.csv"
    if not p.exists(): return None
    rows = list(csv.DictReader(open(p))); cols = [c for c in rows[0] if c.startswith("rms_")]
    return {c: np.array([float(x[c]) for x in rows]) for c in cols}
def triple(base, runs, n):
    a, b, c = (base / r for r in runs)
    if any(not (x / "RUN_RC").exists() or (x / "RUN_RC").read_text().strip() != "0" for x in (a, b, c)): return None, "RUN_RC が 0 でない run がある"
    va, vb, vc = vals(a, n), vals(b, n), vals(c, n); ra, rb, rc = rms(a), rms(b), rms(c)
    if None in (va, vb, vc, ra, rb, rc): return None, "出力が欠けている"
    res = {"bitident_pre": True, "over": []}
    for src, (A, B, C) in (("場", (va, vb, vc)), ("残差", (ra, rb, rc))):
        for k in A:
            if k not in B or k not in C or A[k].shape != C[k].shape: res["over"].append(f"{src} {k}: 形が違う"); continue
            dr = float(np.max(np.abs(B[k] - A[k]))) if A[k].size else 0.0; dp = float(np.max(np.abs(C[k] - A[k]))) if A[k].size else 0.0
            if dr != 0: res["bitident_pre"] = False
            if dp > dr: res["over"].append(f"{src} {k}: 後−前a {dp:.3e} > 前b−前a {dr:.3e}")
    res["pass"] = not res["over"]; return res, None
# R1
names = sorted(p.name for p in HERE.glob("run_07[67]*_r1_*") if p.is_dir())
for tag in ("pt", "b0"):
    for p in ("f32", "fp64"):
        g = [n for n in names if f"_r1_{tag}_{p}_" in n]
        pre = [n for n in g if n.endswith("_pre")]; post = [n for n in g if n.endswith("_post")]
        if len(pre) != 2 or len(post) != 1: und(f"R1 {tag} {p}: run がそろわない {g}"); continue
        r, e = triple(HERE, pre + post, 10)
        if e: und(f"R1 {tag} {p}: {e}"); continue
        OUT["R1"][f"{tag}|{p}"] = r; print("R1", tag, p, "合格" if r["pass"] else "不合格", "前の 2 本がビット一致" if r["bitident_pre"] else "前の 2 本に再実行の差", r["over"][:3])
# R3
names3 = sorted(p.name for p in C48.glob("run_00*_r3_*") if p.is_dir())
for p in ("f32", "fp64"):
    g = [n for n in names3 if f"_r3_{p}_" in n]; pre = [n for n in g if n.endswith("_pre")]; post = [n for n in g if n.endswith("_post")]
    if len(pre) != 2 or len(post) != 1: und(f"R3 {p}: run がそろわない {g}"); continue
    r, e = triple(C48, pre + post, 10)
    if e: und(f"R3 {p}: {e}"); continue
    OUT["R3"][p] = r; print("R3", p, "合格" if r["pass"] else "不合格", "前の 2 本がビット一致" if r["bitident_pre"] else "前の 2 本に再実行の差", r["over"][:3])
# R2
Q = ("ro", "roUx", "roUy", "roe")
def st(r):
    p = HERE / r / "stage.h5"
    if not p.exists(): return None
    with h5py.File(p, "r") as h:
        d = {q: np.asarray(h[f"/d/{q}"][:], np.float64) for q in Q}
        fl = {k: np.asarray(h[f"/cells/{k}_flag"][:]).astype(bool) for k in ("bnode", "axis", "wall", "iso_wall") if f"/cells/{k}_flag" in h}
    return d, fl
r2 = sorted(p.name for p in HERE.glob("run_07*_r2_*") if p.is_dir())
P3pre = [n for n in r2 if "_p3_pre" in n]; P3post = [n for n in r2 if "_p3_post" in n]; P0post = [n for n in r2 if "_p0_post" in n]
S = {n: st(n) for n in P3pre + P3post + P0post}
if len(P3pre) != 2 or len(P3post) != 2 or len(P0post) != 1 or any(v is None for v in S.values()): und(f"R2: run か stage.h5 がそろわない {r2}")
else:
    fl = S[P0post[0]][1]; n = len(S[P0post[0]][0]["ro"])
    ax, wa, iw, bn = (fl.get(k, np.zeros(n, bool)) for k in ("axis", "wall", "iso_wall", "bnode"))
    oth = bn & ~ax & ~wa
    nrm = lambda a: float(np.sqrt(np.sum(a * a)))  # noqa: E731
    cons = {"軸|roUy": (ax, "roUy"), "壁|roUx": (wa & ~ax, "roUx"), "壁|roUy": (wa & ~ax, "roUy"), "等温壁|roe": (iw & ~ax, "roe")}
    ci = {}
    for k, (m, q) in cons.items():
        ci[k] = {"前 P3 max|d|": float(max(np.max(np.abs(S[r][0][q][m])) for r in P3pre)) if m.any() else None,
                 "後 P3 max|d|": float(max(np.max(np.abs(S[r][0][q][m])) for r in P3post)) if m.any() else None,
                 "後 P0 max|d|": float(np.max(np.abs(S[P0post[0]][0][q][m]))) if m.any() else None}
    OUT["R2"]["拘束"] = ci
    ok_i = all(v["後 P3 max|d|"] == 0.0 for v in ci.values() if v["後 P3 max|d|"] is not None)
    REG = {"軸": ax, "壁": wa & ~ax, "入口・出口の境界": oth, "内部": ~bn & ~ax & ~wa}
    tab = {}; ok_ii = True
    for nm, m in REG.items():
        for q in Q:
            ref = S[P0post[0]][0][q][m]; n0 = nrm(ref)
            if n0 == 0: continue
            e_pre = nrm(S[P3pre[0]][0][q][m] - ref) / n0; e_post = nrm(S[P3post[0]][0][q][m] - ref) / n0
            s = max(nrm(S[P3pre[0]][0][q][m] - S[P3pre[1]][0][q][m]), nrm(S[P3post[0]][0][q][m] - S[P3post[1]][0][q][m])) / n0
            g_post = nrm(S[P3post[0]][0][q][m]) / n0
            tab[f"{nm}|{q}"] = dict(e_pre=e_pre, e_post=e_post, s=s, g_post=g_post)
            if nm == "入口・出口の境界" and not (e_pre - e_post > 10 * s): ok_ii = False
            print(f"R2 {nm:8s} {q:5s} e 前 {e_pre:.3e} → 後 {e_post:.3e}  幅 {s:.2e}  g 後 {g_post:.4f}")
    OUT["R2"]["領域"] = tab
    for k, v in ci.items(): print("R2 拘束", k, v)
    OUT["R2"]["verdict"] = ("結果 A (拘束の成分が 0 になり、入口・出口の凍結が直った)" if ok_i and ok_ii else
                            "結果 B (" + ("拘束の違反が残る" if not ok_i else "") + (" 入口・出口の凍結が直らない" if not ok_ii else "") + ")")
r1ok = OUT["R1"] and all(v["pass"] for v in OUT["R1"].values()) and len(OUT["R1"]) == 4
r3ok = OUT["R3"] and all(v["pass"] for v in OUT["R3"].values()) and len(OUT["R3"]) == 2
OUT["verdict"] = ("判定不能: " + "; ".join(OUT["undecidable"][:3])) if OUT["undecidable"] else \
    f"R1 {'合格' if r1ok else '不合格'}、R2 {OUT['R2'].get('verdict')}、R3 {'合格' if r3ok else '不合格'}"
print("== VERDICT §6.59:", OUT["verdict"])
(HERE / "_band_ab/cold_pair/port_reg_judge.json").write_text(json.dumps(OUT, indent=1, ensure_ascii=False, default=float)); print("→ _band_ab/cold_pair/port_reg_judge.json")
