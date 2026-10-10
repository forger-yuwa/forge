"""plan architecture-float-state-double-geometry §6.6 (段 ④) の集計。fg5.sh の後に AWS の case/45 で回す (fg4_an.py を段 ④ の run に読み替えたもの)。
2 (起動時の照合) は fg5.log の [geomStage4] の行を表示する (判定は表を読んで行う: FP64 は全項目 0)。
欠損・形の不一致・非有限は判定不能 (= PASS でない)。結果を標準出力と fg5_verdict.json に書く。"""
import csv, json, h5py, numpy as np
from pathlib import Path
C45 = Path("/home/ubuntu/forge-wallfit/case/45.isobutane_m6_d155"); C56 = Path("/home/ubuntu/forge-wallfit/case/56.gap_tp1187")
V = {}
class Undecidable(Exception): pass
def verdict(key, ok, detail):
    V[key] = {"verdict": "PASS" if ok else "FAIL", "detail": detail}; print(f"  [{'PASS' if ok else 'FAIL'}] {key}: {detail}")
def judge(key, fn):
    try: fn()
    except (Undecidable, OSError, KeyError) as e:
        V[key] = {"verdict": "UNDECIDABLE", "detail": str(e)}; print(f"  [判定不能] {key}: {e}")
def vals(r, n, keys):
    with h5py.File(r / f"res_{n}.h5", "r") as h:
        out = {k: np.asarray(h["VALUE/" + k][:], np.float64) for k in keys if "VALUE/" + k in h}
    miss = [k for k in keys if k not in out]
    if miss: raise Undecidable(f"{r.name} res_{n} に {miss} が無い")
    for k, v in out.items():
        if not np.all(np.isfinite(v)): raise Undecidable(f"{r.name} res_{n} の {k} に非有限")
    return out
def hist(r):
    rows = [x for x in csv.DictReader(open(r / "residual_history.csv")) if x["phase"] == "outer_end"]
    cols = [c for c in rows[0] if c.startswith("rms_") and not c.startswith("rms_dq")]
    a = np.array([[float(x[c]) for c in cols] for x in rows])
    if not np.all(np.isfinite(a)): raise Undecidable(f"{r.name} の残差履歴に非有限")
    return a
KS = ("ro", "roUx", "roUy", "roe", "roK", "roOmega", "P", "T")
def ident(key, base, o, nw, o2, step=20, record=False):
    def fn():
        A, B, A2 = (vals(base / x, step, KS) for x in (o, nw, o2))
        worst = 0.0; parts = []
        for k in KS:
            s = np.max(np.abs(A[k])); dn = np.max(np.abs(B[k] - A[k])) / s; dr = np.max(np.abs(A2[k] - A[k])) / s
            worst = max(worst, dn / dr if dr > 0 else (0.0 if dn == 0 else np.inf)); parts.append(f"{k} {dn:.2e}/{dr:.2e}")
        ha, hb, ha2 = hist(base / o), hist(base / nw), hist(base / o2)
        rn = np.max(np.abs(hb - ha) / np.maximum(np.abs(ha), 1e-300)); rr = np.max(np.abs(ha2 - ha) / np.maximum(np.abs(ha), 1e-300))
        det = f"新−旧 / 再実行の差の最大比 {worst:.2f}; " + "; ".join(parts) + f"; 残差履歴 新−旧 {rn:.2e} / 再実行 {rr:.2e}"
        if record: print(f"  (記録) {key}: {det}"); V[key] = {"verdict": "RECORD", "detail": det}
        else: verdict(key, worst <= 10, det)
    judge(key, fn)
print("== 2 起動時の照合 (fg5.log の [geomStage4] の行)")
for l in open(C45 / "fg5.log"):
    if "geomStage4" in l: print(l.rstrip())
print("== 3 FP64 isp 1 の B0 20 step (旧 = 段 ③、新 = 段 ④)")
ident("3.fp64_isp1.case/45", C45, "run_0429_fg4_isp1_old", "run_0430_fg4_isp1_new", "run_0431_fg4_isp1_old2")
print("== 4 float case/56 20 step (旧 = 段 ③ の run_0085、新 = 段 ④)")
ident("4.float.case/56", C56, "run_0085_fg3_id32_new", "run_0088_fg4_id32_new", "run_0089_fg4_id32_old2")
print("== 4d FP64 isp 1 の 1 step: commit 前の残差と dq = Q₁ − Q₀")
def f4d():
    RK = ("res_ro", "res_roUx", "res_roUy", "res_roe", "res_roK", "res_roOmega"); QK = ("ro", "roUx", "roUy", "roe", "roK", "roOmega")
    o, nw, o2 = (C45 / x for x in ("run_0432_fg4_1s_old", "run_0433_fg4_1s_new", "run_0434_fg4_1s_old2"))
    R = {t: vals(r, 1, RK) for t, r in (("o", o), ("n", nw), ("o2", o2))}
    Q1 = {t: vals(r, 1, QK) for t, r in (("o", o), ("n", nw), ("o2", o2))}
    Q0 = {t: vals(r, 0, QK) for t, r in (("o", o), ("n", nw), ("o2", o2))}
    for k in QK:
        if not (np.array_equal(Q0["o"][k], Q0["n"][k]) and np.array_equal(Q0["o"][k], Q0["o2"][k])): raise Undecidable(f"Q₀ の {k} が腕で違う")
    worst = 0.0; parts = []
    for k in RK:
        s = np.max(np.abs(R["o"][k])); dn = np.max(np.abs(R["n"][k] - R["o"][k])) / s; dr = np.max(np.abs(R["o2"][k] - R["o"][k])) / s
        worst = max(worst, dn / dr if dr > 0 else (0.0 if dn == 0 else np.inf)); parts.append(f"{k} {dn:.2e}/{dr:.2e}")
    for k in QK:
        dqo, dqn, dqo2 = (Q1[t][k] - Q0[t][k] for t in ("o", "n", "o2"))
        s = np.max(np.abs(dqo)); dn = np.max(np.abs(dqn - dqo)) / s; dr = np.max(np.abs(dqo2 - dqo)) / s
        worst = max(worst, dn / dr if dr > 0 else (0.0 if dn == 0 else np.inf)); parts.append(f"dq_{k} {dn:.2e}/{dr:.2e}")
    verdict("4d.fp64_isp1_1step", worst <= 10, f"新−旧 / 再実行の差の最大比 {worst:.2f}; " + "; ".join(parts))
judge("4d.fp64_isp1_1step", f4d)
print("== 5 記録")
def f5b():
    A, B = vals(C45 / "run_0428_fg3_f32_new", 20, KS), vals(C45 / "run_0435_fg4_f32_new", 20, KS)
    parts = [f"{k} {np.max(np.abs(B[k] - A[k])) / np.max(np.abs(A[k])):.2e}" for k in KS]
    ha, hb = hist(C45 / "run_0428_fg3_f32_new"), hist(C45 / "run_0435_fg4_f32_new")
    det = "float case/45 B0 20 step の新−旧: " + "; ".join(parts) + f"; 最終 step の rms (旧 → 新): " + ", ".join(f"{a:.2e}→{b:.2e}" for a, b in zip(ha[-1], hb[-1]))
    print("  (記録) " + det); V["5.float_case45"] = {"verdict": "RECORD", "detail": det}
judge("5.float_case45", f5b)
json.dump(V, open(C45 / "fg5_verdict.json", "w"), ensure_ascii=False, indent=1)
bad = [k for k, v in V.items() if v["verdict"] not in ("PASS", "RECORD")]
print("== 総合:", "ALL PASS" if not bad else f"PASS でないもの {bad}")
