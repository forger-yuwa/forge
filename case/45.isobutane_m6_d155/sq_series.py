"""plan architecture-float-state-double-geometry §6.36: 機械ごとに、各 run の系列 (step 0・20,000〜30,000 の θ_r・Q_w) と前提の記録を作る。
使い方: python3 sq_series.py <A|B> <run> ... → _band_ab/cold_pair/sq_series_<A|B>.json。判定は sq_an.py (2 機械の記録を合わせる)。
記録するもの (判定は sq_an.py):
  - 残差の全行と使う出力 (step 0・20,000〜30,000) の有限性、check_convergence の VERDICT
  - res_0 の 9 保存量のハッシュ (10 本のビット一致と旧 A の res_0 との一致の照合用)
  - 使う出力ごとの、9 保存量で float32 で表せない値の個数 (G は 0 のはず)
  - forge_run.log の [quantizeState] の行 (起動・自己点検)
  - θ_r の系列 (20,000〜30,000) の check_quasisteady (--drift 0.0001 --osc 0.0001 --tail 1.0、記録だけ)"""
import csv, gzip, hashlib, json, math, re, subprocess, sys
from pathlib import Path
import numpy as np, h5py
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cold_series as CS  # noqa: E402
import cold_xcheck as XC  # noqa: E402
M = sys.argv[1]; RUNS = sys.argv[2:]
STEPS = [0] + list(range(20000, 30001, 500))
KEYS = ("theta_r_40", "theta_r_70", "theta_r_94", "Q_w")
TOOLS = Path.home() / "forge-wallfit/solver_density_cuda/tools"
Q9 = ("ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega", "roY0", "roY1")
yb_x, yb = XC.common_yb()
out = {"machine": M, "runs": {}}
def opn(p):
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p)
for name in RUNS:
    run = HERE / name; rec = {"problems": []}
    rh = run / "residual_history.csv"
    rows = [x for x in csv.DictReader(opn(rh if rh.exists() else rh.with_suffix(".csv.gz")))]
    rows = [x for x in rows if x.get("phase", "outer_end") == "outer_end"]
    cols = [c for c in rows[0] if c.startswith("rms_")]
    nonfin = [(x["step"], c) for x in rows for c in cols if not math.isfinite(float(x[c]))]
    rec["residual_rows"] = len(rows); rec["residual_last_step"] = int(rows[-1]["step"]); rec["residual_nonfinite"] = nonfin[:5]
    if nonfin: rec["problems"].append("残差に非有限")
    if int(rows[-1]["step"]) + 1 < 30000: rec["problems"].append("30,000 step に届かない")
    rec["RUN_RC"] = (run / "RUN_RC").read_text().strip() if (run / "RUN_RC").exists() else None
    if rec["RUN_RC"] != "0": rec["problems"].append(f"RUN_RC {rec['RUN_RC']}")
    p = subprocess.run([sys.executable, str(TOOLS / "check_convergence.py"), str(run)], capture_output=True, text=True)
    ov = [l for l in (p.stdout + p.stderr).splitlines() if "OVERALL" in l or "VERDICT" in l]
    rec["check_convergence"] = ov[-1] if ov else None
    if not ov: rec["problems"].append("check_convergence が実行できない")
    rec["config_diff_lines"] = sum(1 for l in (run / "config_vs_run_0487.diff").read_text().splitlines() if l[:1] in "<>")
    # 出力: 有限性・表せない値・res_0 のハッシュ
    rep = {}; ser = {}
    for s in STEPS:
        f = run / f"res_{s}.h5"
        if not f.exists(): rec["problems"].append(f"res_{s} が無い"); continue
        with h5py.File(f, "r") as h:
            V = h["VALUE"]
            if any(V[k].dtype.kind == "f" and not np.all(np.isfinite(V[k][:])) for k in V): rec["problems"].append(f"res_{s} に非有限")
            miss = [q for q in Q9 if q not in V]
            if miss: rec["problems"].append(f"res_{s} に {miss} が無い")
            cnt = {}
            for q in Q9:
                if q not in V: continue
                a = np.asarray(V[q][:]).astype(np.float64)
                cnt[q] = int(np.count_nonzero(a.astype(np.float32).astype(np.float64) != a))
            rep[s] = dict(total=sum(cnt.values()), by_qty=cnt, n=int(np.asarray(V["ro"][:]).size))
            if s == 0:
                hsh = hashlib.sha256()
                for q in Q9:
                    if q in V: hsh.update(q.encode()); hsh.update(np.ascontiguousarray(np.asarray(V[q][:])).tobytes())
                rec["res0_q9_sha256"] = hsh.hexdigest()
        r = CS.snapshot(run, s, yb_x, yb)
        if r.get("nonfinite") or any(not np.isfinite(r[k]) for k in KEYS): rec["problems"].append(f"step {s} の系列に非有限")
        ser[s] = {k: float(r[k]) for k in KEYS}
    rec["not_float32_representable"] = rep; rec["series"] = ser
    log = run / "forge_run.log"
    txt = opn(log if log.exists() else log.with_suffix(".log.gz")).read()
    rec["quantize_log"] = [l[:400] for l in txt.splitlines() if "[quantizeState]" in l][:20]
    # check_quasisteady (記録だけ)
    sc = run / "sq_theta_series.csv"
    with open(sc, "w", newline="") as f:
        w = csv.writer(f); w.writerow(["step"] + list(KEYS))
        for s in sorted(k for k in ser if k >= 20000): w.writerow([s] + [ser[s][k] for k in KEYS])
    q = subprocess.run([sys.executable, str(TOOLS / "check_quasisteady.py"), str(run), "--series-csv", str(sc), "--series-cols", ",".join(KEYS),
                        "--drift", "0.0001", "--osc", "0.0001", "--tail", "1.0"], capture_output=True, text=True)
    rec["check_quasisteady"] = [l for l in q.stdout.splitlines() if any(k in l for k in KEYS) or "VERDICT" in l][:12]
    out["runs"][name] = rec
    print(name, "問題:", rec["problems"] or "なし", "| 表せない値 step 30000:", rep.get(30000, {}).get("total"), "/", rep.get(30000, {}).get("n"), "| check_convergence:", rec["check_convergence"])
ref = HERE / "run_0487_tr_fp64/res_0.h5"
if ref.exists():
    with h5py.File(ref, "r") as h:
        hsh = hashlib.sha256()
        for q in Q9:
            if q in h["VALUE"]: hsh.update(q.encode()); hsh.update(np.ascontiguousarray(np.asarray(h["VALUE"][q][:])).tobytes())
    out["ref_run_0487_res0_sha256"] = hsh.hexdigest()
dst = HERE / f"_band_ab/cold_pair/sq_series_{M}.json"
dst.write_text(json.dumps(out, indent=1, ensure_ascii=False, default=float)); print("→", dst)
