"""plan architecture-float-state-double-geometry §6.42 (事前登録、codex diagnose 2026-10-11) の判定: FP64 のビルドで保存量 9 つを全部 step 末尾に float32 に丸める短い試験 (G9 対 C)。A で回す。
腕: G9 = run_0706〜0708 (新)。C = run_0595〜0597、陽性対照 F = run_0621〜0623 (§6.40 の既存)。参照の型・指標 a は §6.40 (rq_an.py) と同じ定義。
C の再利用の条件: G9 と C で FP64 のバイナリの sha256 (rp.log・rq.log の見出し)・入力・実効設定 (run_0595 との差 0 行)・res_0 の 9 保存量がビット一致。
G9 の自己点検 (量の名前まで検査): (a) step 1・2000 の Q・Q_N・Q_M の全量で表せない値 0。(b) step 1 の roY{0,1}N・M は全 CV で 0、
  または「入口のピンを除く内点」で 0。
判定: e_G9 = mean(a_G9) − mean(a_C)、e_F = mean(a_F) − mean(a_C)、Welch、各比較の t の分位 1 − 0.05/4。
  結果 A: L_F ≥ 0.5 かつ L_G9 ≥ 0.5。結果 B: L_F ≥ 0.5 かつ [L_G9, U_G9] ⊂ [−0.1, 0.1]。それ以外は判別不能。→ _band_ab/cold_pair/rg_judge.json"""
import csv, json, math, re, sys
from pathlib import Path
import numpy as np, h5py
from scipy import stats
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cold_xcheck as XC  # noqa: E402
C = ("run_0595_rq_c1", "run_0596_rq_c2", "run_0597_rq_c3")
G9 = ("run_0706_rp_g9a", "run_0707_rp_g9b", "run_0708_rp_g9c")
F = ("run_0621_rq_f1", "run_0622_rq_f2", "run_0623_rq_f3")
REF = {"A": "run_0487_tr_fp64", "B": "run_0488_tr_f32", "Bp": "run_0489_tr_f32b"}
S = (1000, 1500, 2000)
Q9 = ("ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega", "roY0", "roY1")
QA = {"Q": list(Q9), "Q_N": [q + "N" for q in Q9[:7]], "Q_M": [q + "M" for q in Q9[:7]]}
QB = ["roY0N", "roY1N", "roY0M", "roY1M"]
OUT = {"undecidable": [], "pre": {}, "records": {}}
def und(m): OUT["undecidable"].append(m); print("  [判定不能]", m)
def sha_in(log):
    t = (HERE / log).read_text(errors="replace"); m = re.findall(r"F64 ([0-9a-f]{16})", t); return m[-1] if m else None
s_rq, s_rp = sha_in("rq.log"), sha_in("rp.log")
OUT["pre"]["FP64 のバイナリ (rq.log / rp.log)"] = [s_rq, s_rp]
if not (s_rq and s_rq == s_rp): und("G9 と C の FP64 のバイナリが同じと確かめられない")
for r in G9 + C + F:
    d = HERE / r
    rc = (d / "RUN_RC").read_text().strip() if (d / "RUN_RC").exists() else None
    if rc != "0": und(f"{r}: RUN_RC {rc}"); continue
    rows = [x for x in csv.DictReader(open(d / "residual_history.csv")) if x.get("phase", "outer_end") == "outer_end"]
    cols = [c for c in rows[0] if c.startswith("rms_")]
    if any(not math.isfinite(float(x[c])) for x in rows for c in cols): und(f"{r}: 残差に非有限")
    if int(rows[-1]["step"]) + 1 < 2000: und(f"{r}: 2,000 step に届かない")
    for s in (0, 500, 1000, 1500, 2000):
        with h5py.File(d / f"res_{s}.h5", "r") as h:
            if any(h["VALUE"][k].dtype.kind == "f" and not np.all(np.isfinite(h["VALUE"][k][:])) for k in h["VALUE"]): und(f"{r}: res_{s} に非有限")
for r in G9:
    df = HERE / r / "config_vs_run_0595.diff"
    n = sum(1 for l in df.read_text().splitlines() if l[:1] in "<>") if df.exists() else None
    if n != 0: und(f"{r}: run_0595 との設定の差が {n} 行")
def q9(run):
    with h5py.File(HERE / run / "res_0.h5", "r") as h: return {k: np.asarray(h["VALUE/" + k][:]) for k in Q9 if k in h["VALUE"]}
if not OUT["undecidable"]:
    r0 = q9(C[0])
    for r in C[1:] + G9:
        bad = [k for k, v in q9(r).items() if not np.array_equal(v, r0[k])]
        if bad or len(q9(r)) != 9: und(f"{r}: res_0 の 9 保存量が {C[0]} と違う {bad}")
def pairs(seg): return dict((k, int(v)) for k, v in re.findall(r"(ro[A-Za-z0-9]*) (\d+)", seg))
for r in G9:
    lines = [l for l in (HERE / r / "forge_run.log").read_text(errors="replace").splitlines() if "[quantizeState]" in l]
    OUT["pre"][f"{r} 自己点検"] = [l[:300] for l in lines if re.search(r"step \d+:", l)][:8]
    if not any("保存量 9 本" in l for l in lines): und(f"{r}: 9 本を丸める起動の記録が無い")
    for st in (1, 2000):
        la = [l for l in lines if f"step {st}:" in l and "表せない値" in l and "(Q)" in l]
        if not la: und(f"{r}: step {st} の (a) の行が無い"); continue
        segs = dict(zip(("Q", "Q_N", "Q_M"), re.split(r"\((?:Q|Q_N|Q_M)\)", la[0].split("表せない値", 1)[1])[:3]))
        for g, names in QA.items():
            pv = pairs(segs.get(g, ""))
            miss = [n for n in names if n not in pv]; nz = {n: pv[n] for n in names if pv.get(n)}
            if miss or nz: und(f"{r}: step {st} の (a) {g} で 名前が無い {miss}・表せない値 {nz}")
    lb = [l for l in lines if "step 1:" in l and "(roY_N)" in l]
    if not lb: und(f"{r}: step 1 の (b) の行が無い")
    else:
        pv = pairs(lb[0].split("表せない値", 1)[1].split("合計")[0])
        miss = [n for n in QB if n not in pv]
        if miss: und(f"{r}: (b) で名前が無い {miss}")
        elif any(pv[n] for n in QB):
            note = [l for l in lines if "step 1:" in l and "入口のピン" in l and "を除く内点では" in l]
            ok = bool(note) and re.search(r"を除く内点では.*合計 0(?!\d)", note[-1])
            OUT["pre"][f"{r} (b) 0 でない (全 CV)"] = {n: pv[n] for n in QB}
            if not ok: und(f"{r}: (b) で入口のピンを除く内点にも表せない値がある: {note[-1][:200] if note else '注意の行が無い'}")
verdict = None
if not OUT["undecidable"]:
    with h5py.File(HERE / REF["A"] / "nozzle.h5", "r") as h: xy0 = np.asarray(h["MESH/COORD"][:], float)
    ni, nj, Sc = XC.mesh_info(xy0.size // 3); X = (xy0.reshape(-1, 3)[:, 0] / Sc).reshape(ni, nj)
    I, J = np.meshgrid(np.arange(ni), np.arange(nj), indexing="ij")
    DOM = ((X < -2) & (I >= 1) & (J >= 1) & (J <= nj - 2)).ravel()
    def T(run, s):
        with h5py.File(HERE / run / f"res_{s}.h5", "r") as h: return np.asarray(h["VALUE/T"][:], float)
    TA0 = T(REF["A"], 0); z = lambda run, s: (T(run, s) - T(run, 0)) / TA0  # noqa: E731
    V = {s: 0.5 * (z(REF["B"], s) + z(REF["Bp"], s)) - z(REF["A"], s) for s in S}
    a = lambda run: float(np.mean([np.dot(V[s][DOM], z(run, s)[DOM]) / np.dot(V[s][DOM], V[s][DOM]) for s in S]))  # noqa: E731
    A_ = {g: [a(r) for r in runs] for g, runs in (("C", C), ("G9", G9), ("F", F))}
    OUT["a"] = A_; tq = 1 - 0.05 / 4; res = {}
    for g in ("G9", "F"):
        x, c = np.array(A_[g]), np.array(A_["C"]); e = x.mean() - c.mean(); vx, vc = x.var(ddof=1) / len(x), c.var(ddof=1) / len(c); se = math.sqrt(vx + vc)
        if se == 0: df_, t = float("inf"), 0.0
        else: df_ = (vx + vc) ** 2 / ((vx ** 2 / (len(x) - 1) if vx else 0) + (vc ** 2 / (len(c) - 1) if vc else 0)); t = float(stats.t.ppf(tq, df_))
        res[g] = dict(e=e, se=se, df=df_, t=t, L=e - t * se, U=e + t * se)
        print(f"e_{g} = {e:+.4f}  SE {se:.4f}  df {df_:.1f}  t {t:.2f}  区間 [{e - t * se:+.4f}, {e + t * se:+.4f}]  (a_{g} = {[round(v, 4) for v in A_[g]]}、a_C = {[round(v, 4) for v in A_['C']]})")
    OUT["compare"] = res
    pos = res["F"]["L"] >= 0.5
    if pos and res["G9"]["L"] >= 0.5: verdict = "結果 A: 9 保存量の末尾の量子化だけで、登録した成分の少なくとも半分を生じさせられる (十分性)"
    elif pos and -0.1 <= res["G9"]["L"] and res["G9"]["U"] <= 0.1: verdict = "結果 B: 9 保存量の末尾の量子化だけでは不十分 (第 2 仮説を退ける。残差の演算への局在化はしない)"
    else: verdict = "判別不能" + ("" if pos else " (陽性対照 F が成立しない)")
if verdict is None: verdict = "判定不能: " + "; ".join(OUT["undecidable"][:4])
OUT["verdict"] = verdict
print(f"== VERDICT §6.42: {verdict}")
(HERE / "_band_ab/cold_pair/rg_judge.json").write_text(json.dumps(OUT, indent=1, ensure_ascii=False, default=float)); print("→ _band_ab/cold_pair/rg_judge.json")
