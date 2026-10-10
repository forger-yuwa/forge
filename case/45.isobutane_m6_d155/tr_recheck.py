"""float plan §6.27 の確かめ直し (codex diagnose 2026-10-11 の Major 3): 判定の前提として、3 本の全期間・全出力の有限性と、全残差列の有限性を見る。
- residual_history.csv: outer_end の全行・全 rms_* 列が有限であること。
- 全出力 res_<n>.h5 (0〜30,000、500 ごと): VALUE の全データセットが有限であること。
- 正式ツール (check_convergence・check_quasisteady) は tr_judge.json に OVERALL の行が記録されていること (実行できたこと)。
結果は標準出力と _band_ab/cold_pair/tr_recheck.json。"""
import csv, json, math
from pathlib import Path
import numpy as np, h5py
HERE = Path(__file__).resolve().parent
RUNS = ("run_0487_tr_fp64", "run_0488_tr_f32", "run_0489_tr_f32b")
OUT = {"problems": []}
for r in RUNS:
    run = HERE / r
    rows = [x for x in csv.DictReader(open(run / "residual_history.csv")) if x.get("phase", "outer_end") == "outer_end"]
    cols = [c for c in rows[0] if c.startswith("rms_")]
    bad = next(((x["step"], c) for x in rows for c in cols if not math.isfinite(float(x[c]))), None)
    if bad: OUT["problems"].append(f"{r}: 残差に非有限 {bad}")
    if int(rows[-1]["step"]) + 1 < 30000: OUT["problems"].append(f"{r}: 残差が 30,000 step に届かない")
    nbad = 0; nfiles = 0
    for s in range(0, 30001, 500):
        p = run / f"res_{s}.h5"
        if not p.exists():
            if s == 0: continue
            OUT["problems"].append(f"{r}: res_{s} が無い"); continue
        nfiles += 1
        with h5py.File(p, "r") as h:
            V = h["VALUE"]
            for k in V:
                a = np.asarray(V[k][:])
                if a.dtype.kind == "f" and not np.all(np.isfinite(a)): nbad += 1; OUT["problems"].append(f"{r}: res_{s} の {k} に非有限")
    print(f"{r}: 残差 {len(rows)} 行・{len(cols)} 列 有限 {bad is None}、出力 {nfiles} 個の全データセットの非有限 {nbad}")
tj = json.loads((HERE / "_band_ab/cold_pair/tr_judge.json").read_text())
for k, v in tj.get("tools", {}).items():
    if "check_" in k and not (isinstance(v, str) and v.startswith("OVERALL")):
        OUT["problems"].append(f"ツールが実行できていない: {k} = {v}")
print("== tr_judge の VERDICT:", tj.get("verdict"))
OUT["verdict_recheck"] = "前提を満たす (判定は確定)" if not OUT["problems"] else "前提を満たさない: " + "; ".join(OUT["problems"][:5])
print("== 確かめ直し:", OUT["verdict_recheck"])
(HERE / "_band_ab/cold_pair/tr_recheck.json").write_text(json.dumps(OUT, indent=1, ensure_ascii=False))
