#!/usr/bin/env python3
"""V2b 保存収支・V2f 局所逆流の判定 (plan boundary-node-farfield-characteristic §6 V2b/V2f)。
帳簿ダンプ (FORGE_DUMP_LEDGER、全節点) と farfield 面ダンプ (FORGE_DUMP_FARFIELD) から:

  identity RUN [--call N]   離散恒等式 (評価 N 回目、既定 1):
      Σ_節点 対流残差 (res_after_conv) = −Σ_farfield 面流束       (ρ, ρu, ρE)
      Σ_節点 化学種残差 (res_after_species) = −Σ ṁ Y_面  (Y_面 = ṁ≥0 で節点値、ṁ<0 で外側状態 R_Y0)
      Σ_節点 k・ω 輸送残差 (res_after_rans_transport) = −Σ ṁ φ_面
      判定 |差| ≤ 1e-5 × その保存量の全面流束の絶対和 (内部面 = 帳簿の面記録、スカラーは Σ|ṁ| × max|φ|)。
      寄与は「担当段の直後の値」そのもの (ρ・ρu・ρE = res_after_conv、k・ω = res_after_rans_transport、種 = res_after_species)。
      各残差配列は担当段の中でゼロに戻されてから書かれ、2 回目以降の評価では段の前の値に前回の最終残差が残っている
      (2026-09-29 run_0036_v2f で確認: 評価 2 の Σres_before_conv(ρ) = 評価 1 の対流残差和)。段の前後の差は 1 回目でしか成り立たない。
      Σ種残差 = 質量の対流残差 (相対 1e-6、Σ種流束 = 質量流束)。
  balance RUN_RESTART META_RUN   定常後の全体収支: 最終場から restart した 1 評価の res_final の全節点和 R を
      代表量 (ρ∞|u∞|A、(ρ∞|u∞|²+P∞)A、ρ∞|u∞|H∞A、ρ∞|u∞|k∞A、ρ∞|u∞|ω∞A、A = farfield 総面積) で割り ≤ 1e-4。
  replaced RUN              最終 res の |Y_EXH − 0| ≤ 1e-4 (内部が自由流に置き換わる)。
  v2f RUN                   面ダンプ全回で xmax の ṁ<0 の面数と、その面の外側組成 R_Y0 = 外気 (0)、ṁ>0 の面は節点値を運ぶことを
                            帳簿の化学種恒等式 (identity の化学種行) で確認する。
"""
import csv, glob, json, os, re, sys
from collections import defaultdict
import numpy as np

MOM = ("F_roUx", "F_roUy", "F_roUz")


def read_ledger(run, call):
    res = defaultdict(dict); st = defaultdict(dict)
    with open(run + "/ledger.csv") as fh:
        r = csv.reader(fh); next(r)
        c = str(call)
        for row in r:
            if row[0] != c:
                continue
            tag, node, fld, v = row[1], int(row[2]), row[3], float(row[4])
            if tag.startswith("res_"):
                res[(tag, fld)][node] = v
            elif tag == "after_eos_bc":
                st[fld][node] = v
    return res, st


def read_ff(run, call):
    rows = []
    for p in sorted(glob.glob(run + "/ffdump.*.csv")):
        m = re.search(r"ffdump\.(\d+)(?:\.(\d+))?\.csv$", p)
        pid, cl = int(m.group(1)), int(m.group(2) or 1)
        if cl != call:
            continue
        for r in csv.DictReader(open(p)):
            d = {k: float(v) for k, v in r.items()}; d["physID"] = pid
            rows.append(d)
    return rows


def faces_abs(run, call):
    s = defaultdict(float); seen = set(); mabs = 0.0
    with open(run + "/ledger.csv.faces") as fh:
        for r in csv.DictReader(fh):
            if r["call"] != str(call) or r["ip"] in seen:
                continue
            seen.add(r["ip"])
            s["ro"] += abs(float(r["F_ro"])); s["roe"] += abs(float(r["F_roe"]))
            s["mom"] += sum(abs(float(r[k])) for k in MOM)
            mabs += abs(float(r["F_ro"]))
    return s, mabs


def total(d):
    return sum(d.values())


def identity(run, call=1):
    res, st = read_ledger(run, call)
    ff = read_ff(run, call)
    if not ff:
        raise SystemExit(f"{run}: 面ダンプ (評価 {call}) が無い")
    fa, mabs_int = faces_abs(run, call)
    ok = True
    print(f"{run} 評価 {call}: farfield 面 {len(ff)}、置換 {sum(int(f['vacuum']) != 0 for f in ff)}、退避 {sum(int(f['hll']) != 0 for f in ff)}")
    conv = lambda q: dict(res[("res_after_conv", "res_" + q)])
    for q, fk, sk in (("ro", "F_ro", "ro"), ("roUx", "F_roUx", "mom"), ("roUy", "F_roUy", "mom"), ("roUz", "F_roUz", "mom"), ("roe", "F_roe", "roe")):
        a = total(conv(q)); b = -sum(f[fk] for f in ff)
        scale = fa[sk] + sum(abs(f[k]) for f in ff for k in ((fk,) if sk != "mom" else MOM))
        r = abs(a - b) / scale
        ok &= r <= 1e-5
        print(f"  {q:8s}: Σ残差 {a:+.6e}  −Σ面流束 {b:+.6e}  |差|/規模 {r:.2e} {'OK' if r <= 1e-5 else 'NG'}")
    mabs = mabs_int + sum(abs(f["F_ro"]) for f in ff)
    ic = lambda f: int(f["ic"])
    for q, t0, t1, face_val in (
            ("roY0", "res_after_rans_transport", "res_after_species", lambda f: st["Y0"][ic(f)] if f["F_ro"] >= 0 else f["R_Y0"]),
            ("roK", "res_after_conv", "res_after_rans_transport", lambda f: st_k(st, ic(f)) if f["F_ro"] >= 0 else f["R_k"]),
            ("roOmega", "res_after_conv", "res_after_rans_transport", lambda f: st_om(st, ic(f)) if f["F_ro"] >= 0 else f["R_om"])):
        if ("res_after_species", "res_" + q) not in res:
            continue
        d = dict(res[(t1, "res_" + q)])
        a = total(d); b = -sum(f["F_ro"] * face_val(f) for f in ff)
        phimax = max(abs(face_val(f)) for f in ff)
        scale = mabs * max(phimax, 1e-30)
        r = abs(a - b) / scale
        ok &= r <= 1e-5
        print(f"  {q:8s}: Σ残差 {a:+.6e}  −Σṁφ_面 {b:+.6e}  |差|/規模 {r:.2e} {'OK' if r <= 1e-5 else 'NG'}")
    if ("res_after_species", "res_roY1") in res:
        sp = total({n: res[("res_after_species", "res_roY0")][n] + res[("res_after_species", "res_roY1")][n]
                    for n in res[("res_after_species", "res_roY0")]})
        m = total(conv("ro"))
        r = abs(sp - m) / max(abs(m), mabs * 1e-12)
        ok &= r <= 1e-6 or abs(sp - m) <= 1e-6 * mabs
        print(f"  Σ種 vs 質量: {sp:+.6e} vs {m:+.6e} (相対 {r:.2e}、|差|/Σ|ṁ| {abs(sp - m) / mabs:.2e})")
    print(f"VERDICT identity: {'PASS' if ok else 'FAIL'}")
    return ok


def st_k(st, n):
    return st["roK"][n] / st["ro"][n]


def st_om(st, n):
    return st["roOmega"][n] / st["ro"][n]


def balance(run, meta_run):
    meta = json.load(open(meta_run + "/v2b_meta.json"))
    res, _ = read_ledger(run, 1)
    ff = read_ff(run, 1)
    A = sum(f["S"] for f in ff)
    rU = meta["ro"] * meta["U"]
    rep = {"ro": rU * A, "roUx": (rU * meta["U"] + meta["P"]) * A, "roUy": (rU * meta["U"] + meta["P"]) * A, "roUz": (rU * meta["U"] + meta["P"]) * A,
           "roe": rU * meta["H"] * A, "roY0": rU * A, "roY1": rU * A, "roK": rU * meta["k"] * A, "roOmega": rU * meta["omega"] * A}
    ok = True
    print(f"{run}: 定常後の全体収支 (res_final の全節点和 / 代表量、A = {A:.4g} m²)")
    for q, R in rep.items():
        key = ("res_final", "res_" + q)
        if key not in res:
            continue
        v = total(res[key]) / R
        ok &= abs(v) <= 1e-4
        print(f"  {q:8s}: {v:+.3e} {'OK' if abs(v) <= 1e-4 else 'NG'}")
    print(f"VERDICT balance: {'PASS' if ok else 'FAIL'} (≤ 1e-4)")
    return ok


def replaced(run):
    import h5py
    fs = sorted((int(re.search(r"res_(\d+)\.h5$", p).group(1)), p) for p in glob.glob(run + "/res_*.h5") if re.search(r"/res_\d+\.h5$", p))
    st, last = fs[-1]
    with h5py.File(last) as f:
        V = f["VALUE"]
        Y = V["roY0"][:] / V["ro"][:] if "roY0" in V else V["Y0"][:]
        nonfin = sum(int(np.sum(~np.isfinite(V[k][:]))) for k in V.keys())
    m = float(np.max(np.abs(Y)))
    ok = m <= 1e-4 and nonfin == 0
    print(f"{run} step {st}: max|Y_EXH| {m:.3e}、非有限 {nonfin} → VERDICT replaced: {'PASS' if ok else 'FAIL'} (≤ 1e-4)")
    return ok


def v2f(run):
    calls = sorted({int(re.search(r"ffdump\.2(?:\.(\d+))?\.csv$", p).group(1) or 1) for p in glob.glob(run + "/ffdump.2*.csv")})
    nin_tot = 0; bad = 0
    for c in calls:
        ff = [f for f in read_ff(run, c) if f["physID"] == 2]
        neg = [f for f in ff if f["F_ro"] < 0]
        nin_tot += len(neg)
        bad += sum(abs(f["R_Y0"]) > 1e-7 for f in neg)
        vac = sum(int(f["vacuum"]) != 0 or int(f["hll"]) != 0 for f in ff)
        print(f"  評価 {c:2d}: xmax {len(ff)} 面のうち ṁ<0 {len(neg)} 面 (外側組成 ≠ 外気 {sum(abs(f['R_Y0']) > 1e-7 for f in neg)})、置換/退避 {vac}")
    print(f"VERDICT v2f-faces: {'PASS' if nin_tot > 0 and bad == 0 else 'FAIL'} (ṁ<0 の面 延べ {nin_tot}、外側組成が外気でない面 {bad})")
    return nin_tot > 0 and bad == 0


if __name__ == "__main__":
    mode = sys.argv[1]
    call = int(sys.argv[sys.argv.index("--call") + 1]) if "--call" in sys.argv else 1
    ok = {"identity": lambda: identity(sys.argv[2], call), "balance": lambda: balance(sys.argv[2], sys.argv[3]),
          "replaced": lambda: replaced(sys.argv[2]), "v2f": lambda: v2f(sys.argv[2])}[mode]()
    sys.exit(0 if ok else 1)
