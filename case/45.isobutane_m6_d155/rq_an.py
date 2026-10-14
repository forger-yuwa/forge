"""plan architecture-float-state-double-geometry §6.40 (事前登録、codex diagnose 2026-10-11) の判定: FP64 のビルドで ρE だけを step 末尾に float32 に丸める短い試験。A で回す (B が満杯のため run の前に変更)。
腕: C = run_0595〜0597 (FP64 OFF)、E = run_0598・0599・0620 (FP64 + ρE の量子化)、F = run_0621〜0623 (float、陽性対照)。参照: A = run_0487・B = run_0488・B' = run_0489。
参照の型 (固定): 評価域は x/r_t < −2 の内点 (軸 j = 0・壁 j = nj−1・入口 i = 0 を除く)。z_{r,s} = (T_{r,s} − T_{r,0}) / T_{A,0}、v_s = mean(z_B, z_B') − z_A (s = 1000・1500・2000)。
指標: a_r = mean_s ⟨v_s, z_{r,s}⟩ / ⟨v_s, v_s⟩。e_E = mean(a_E) − mean(a_C)、e_F = mean(a_F) − mean(a_C)。Welch の区間、2 比較で全体 95 % (各 t の分位 1 − 0.05/4)。
結果 A: e_F の下限 ≥ 0.5 かつ e_E の下限 ≥ 0.5。結果 B: e_F の下限 ≥ 0.5 かつ e_E の区間 ⊂ [−0.1, 0.1]。それ以外 (陽性対照の不成立を含む) は判別不能。
前提: C・E の res_0 の 9 保存量のビット一致、E の自己点検 (step 1・2000) の表せない値 0 と選ばなかった量の変化 0、全残差・保存した場が有限、座標の一致。
記録: 参照の A・B・B' の a、F の初期化の差、境界の節点の別集計、check_convergence、絶対値の分布。結果は _band_ab/cold_pair/rq_judge.json。"""
import csv, json, math, re, subprocess, sys
from pathlib import Path
import numpy as np, h5py
from scipy import stats
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cold_xcheck as XC  # noqa: E402
C = ("run_0595_rq_c1", "run_0596_rq_c2", "run_0597_rq_c3")
E = ("run_0598_rq_e1", "run_0599_rq_e2", "run_0620_rq_e3")
F = ("run_0621_rq_f1", "run_0622_rq_f2", "run_0623_rq_f3")
REF = {"A": "run_0487_tr_fp64", "B": "run_0488_tr_f32", "Bp": "run_0489_tr_f32b"}
S = (1000, 1500, 2000)
Q9 = ("ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega", "roY0", "roY1")
TOOLS = Path.home() / "forge-wallfit/solver_density_cuda/tools"
OUT = {"undecidable": [], "pre": {}, "records": {}}
def und(m): OUT["undecidable"].append(m); print("  [判定不能]", m)
with h5py.File(HERE / REF["A"] / "nozzle.h5", "r") as h: xy0 = np.asarray(h["MESH/COORD"][:], float)
def T(run, s):
    with h5py.File(HERE / run / f"res_{s}.h5", "r") as h: return np.asarray(h["VALUE/T"][:], float)
nn = xy0.size // 3; ni, nj, Sc = XC.mesh_info(nn)
X = (xy0.reshape(-1, 3)[:, 0] / Sc).reshape(ni, nj)
I, J = np.meshgrid(np.arange(ni), np.arange(nj), indexing="ij")
up = X < -2
DOM = (up & (I >= 1) & (J >= 1) & (J <= nj - 2)).ravel()
BND = (up & ~((I >= 1) & (J >= 1) & (J <= nj - 2))).ravel()
OUT["pre"]["評価域の節点数 (内点 / 境界)"] = [int(DOM.sum()), int(BND.sum())]
# 前提
for r in C + E + F:
    d = HERE / r
    rc = (d / "RUN_RC").read_text().strip() if (d / "RUN_RC").exists() else None
    if rc != "0": und(f"{r}: RUN_RC {rc}"); continue
    rows = [x for x in csv.DictReader(open(d / "residual_history.csv")) if x.get("phase", "outer_end") == "outer_end"]
    cols = [c for c in rows[0] if c.startswith("rms_")]
    if any(not math.isfinite(float(x[c])) for x in rows for c in cols): und(f"{r}: 残差に非有限")
    if int(rows[-1]["step"]) + 1 < 2000: und(f"{r}: 2,000 step に届かない")
    for s in (0, 500, 1000, 1500, 2000):
        p = d / f"res_{s}.h5"
        if not p.exists(): und(f"{r}: res_{s} が無い"); continue
        with h5py.File(p, "r") as h:
            if any(h["VALUE"][k].dtype.kind == "f" and not np.all(np.isfinite(h["VALUE"][k][:])) for k in h["VALUE"]): und(f"{r}: res_{s} に非有限")
    with h5py.File(d / "nozzle.h5", "r") as h:
        if not np.array_equal(np.asarray(h["MESH/COORD"][:], float), xy0): und(f"{r}: 座標が run_0487 と違う")
    p = subprocess.run([sys.executable, str(TOOLS / "check_convergence.py"), str(d)], capture_output=True, text=True)
    ov = [l for l in (p.stdout + p.stderr).splitlines() if "OVERALL" in l or "VERDICT" in l]
    OUT["records"][f"{r} check_convergence"] = ov[-1] if ov else "実行できない"
def q9(run):
    with h5py.File(HERE / run / "res_0.h5", "r") as h: return {k: np.asarray(h["VALUE/" + k][:]) for k in Q9 if k in h["VALUE"]}
if not OUT["undecidable"]:
    r0 = q9(C[0])
    if len(r0) != 9: und(f"res_0 の保存量が 9 でない: {list(r0)}")
    for r in C[1:] + E:
        bad = [k for k, v in q9(r).items() if not np.array_equal(v, r0[k])]
        if bad: und(f"{r}: res_0 の保存量が {C[0]} と違う {bad}")
    rf = q9(F[0])
    OUT["records"]["F の初期化の差 (res_0、C との最大相対差)"] = {k: float(np.max(np.abs(rf[k].astype(float) - r0[k].astype(float)) / np.maximum(np.abs(r0[k].astype(float)), 1e-30))) for k in r0}
    for r in E:
        lines = [l for l in (HERE / r / "forge_run.log").read_text(errors="replace").splitlines() if "[quantizeState]" in l]
        OUT["pre"][f"{r} 自己点検"] = [l[:300] for l in lines if re.search(r"step \d+:", l)][:8]
        for st in (1, 2000):
            ls = [l for l in lines if f"step {st}:" in l]
            rep = [l for l in ls if "表せない" in l]; chg = [l for l in ls if "選ばなかった量の丸め前後の変化" in l]
            if not rep or not chg: und(f"{r}: step {st} の自己点検の行が無い (表せない値 {len(rep)}、変化 {len(chg)})"); continue
            if any(not re.search(r"合計 0(?!\d)", l) for l in rep): und(f"{r}: step {st} で表せない値が 0 でない: {rep[0][:200]}")
            # 行は「選ばなかった量の丸め前後の変化 ro 0 … / 合計 0 (0 が期待。全 CV <節点数>、…)」。節点数を拾わないよう「合計 0」で見る
            if any(not re.search(r"合計 0(?!\d)", l) for l in chg): und(f"{r}: step {st} で選ばなかった量が丸めで変わった: {chg[0][:200]}")
        if not any("roe" in l and ("選んだ" in l or "ONLY" in l or "だけ" in l) for l in lines): und(f"{r}: ρE だけを選んだ起動の記録が無い")
# 参照の型と指標
verdict = None
if not OUT["undecidable"]:
    TA0 = T(REF["A"], 0)
    z = lambda run, s: (T(run, s) - T(run, 0)) / TA0  # noqa: E731
    V = {s: 0.5 * (z(REF["B"], s) + z(REF["Bp"], s)) - z(REF["A"], s) for s in S}
    def a(run, m=DOM): return float(np.mean([np.dot(V[s][m], z(run, s)[m]) / np.dot(V[s][m], V[s][m]) for s in S]))
    A = {g: [a(r) for r in runs] for g, runs in (("C", C), ("E", E), ("F", F))}
    OUT["a"] = A
    OUT["records"]["参照の a (A・B・B')"] = {k: a(v) for k, v in REF.items()}
    OUT["records"]["境界の節点の a (記録)"] = {g: [a(r, BND) for r in runs] for g, runs in (("C", C), ("E", E), ("F", F))}
    OUT["records"]["|v_s| の評価域の RMS"] = {s: float(np.sqrt(np.mean(V[s][DOM] ** 2))) for s in S}
    tq = 1 - 0.05 / 4; res = {}
    for g in ("E", "F"):
        x, c = np.array(A[g]), np.array(A["C"]); e = x.mean() - c.mean(); vx, vc = x.var(ddof=1) / len(x), c.var(ddof=1) / len(c); se = math.sqrt(vx + vc)
        if se == 0: df, t = float("inf"), 0.0
        else: df = (vx + vc) ** 2 / ((vx ** 2 / (len(x) - 1) if vx else 0) + (vc ** 2 / (len(c) - 1) if vc else 0)); t = float(stats.t.ppf(tq, df))
        res[g] = dict(e=e, se=se, df=df, t=t, L=e - t * se, U=e + t * se)
        print(f"e_{g} = {e:+.4f}  SE {se:.4f}  df {df:.1f}  t {t:.2f}  区間 [{e - t * se:+.4f}, {e + t * se:+.4f}]  (a_{g} = {[round(v, 4) for v in A[g]]}、a_C = {[round(v, 4) for v in A['C']]})")
    OUT["compare"] = res
    pos = res["F"]["L"] >= 0.5
    if pos and res["E"]["L"] >= 0.5: verdict = "結果 A: ρE の末尾の量子化だけで、登録した初期の温度差の成分を少なくとも半分の規模で生じさせられる (十分性。float の原因・必要性とは呼ばない)"
    elif pos and -0.1 <= res["E"]["L"] and res["E"]["U"] <= 0.1: verdict = "結果 B: ρE 単独の末尾の量子化では、その成分は出ない (第 1 仮説を退ける)"
    else: verdict = "判別不能" + ("" if pos else " (陽性対照 F が成立しない)")
if verdict is None: verdict = "判定不能: " + "; ".join(OUT["undecidable"][:4])
OUT["verdict"] = verdict
print(f"== VERDICT §6.40: {verdict}")
(HERE / "_band_ab/cold_pair/rq_judge.json").write_text(json.dumps(OUT, indent=1, ensure_ascii=False, default=float)); print("→ _band_ab/cold_pair/rq_judge.json")
