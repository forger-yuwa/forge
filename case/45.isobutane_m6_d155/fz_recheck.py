"""float plan §6.29 の確かめ直し (codex diagnose 2026-10-11 の Major 2): fz_an.py が見ていなかった前提を見る。判定の閾値は変えない。
- 抽出した系列 (θ_r(40/70/94)・Q_w) が 3 本・全 61 出力で有限で、点数がそろうこと (fz_judge.json の series)。
- 3 本の格子の座標が同じ (共通の入力) であること。
- 20,000〜30,000 step の系列の準定常の VERDICT (check_quasisteady、θ_r は --drift 0.0001 --osc 0.0001、Q_w は 0.0002、--tail 1.0)。記録であって判定の前提ではない。
結果は _band_ab/cold_pair/fz_recheck.json。"""
import json, subprocess, sys
from pathlib import Path
import numpy as np, h5py
HERE = Path(__file__).resolve().parent
ARMS = {"A": "run_0491_fz_fp64", "B": "run_0492_fz_f32", "Bp": "run_0493_fz_f32b"}
TOOLS = Path.home() / "forge-wallfit/solver_density_cuda/tools"
KEYS = ("theta_r_40", "theta_r_70", "theta_r_94"); QW = "Q_w"
fz = json.loads((HERE / "_band_ab/cold_pair/fz_judge.json").read_text())
OUT = {"problems": [], "quasisteady": {}}
coords = {}
for a, r in ARMS.items():
    ser = {int(k): v for k, v in fz["series"][a].items()}
    steps = sorted(ser)
    if steps != list(range(0, 30001, 500)): OUT["problems"].append(f"{a}: 系列の step がそろわない ({len(steps)} 点)")
    nf = [(s, k) for s in steps for k in KEYS + (QW,) if not np.isfinite(ser[s][k])]
    if nf: OUT["problems"].append(f"{a}: 系列に非有限 {nf[:3]}")
    with h5py.File(HERE / r / "nozzle.h5", "r") as h: coords[a] = np.asarray(h["MESH/COORD"][:])
    sc = HERE / r / "fz_series_20k_30k.csv"
    with open(sc, "w") as fh:
        fh.write("step," + ",".join(KEYS + (QW,)) + "\n")
        for s in steps:
            if 20000 <= s <= 30000: fh.write(f"{s}," + ",".join(f"{ser[s][k]:.10g}" for k in KEYS + (QW,)) + "\n")
    for cols, thr in ((",".join(KEYS), "0.0001"), (QW, "0.0002")):
        p = subprocess.run([sys.executable, str(TOOLS / "check_quasisteady.py"), "--series-csv", str(sc), "--series-cols", cols, "--tail", "1.0", "--drift", thr, "--osc", thr], capture_output=True, text=True)
        ov = [l for l in (p.stdout + p.stderr).splitlines() if l.startswith("OVERALL")]
        if not ov: OUT["problems"].append(f"{a}: check_quasisteady ({cols}) が実行できない")
        OUT["quasisteady"][f"{a} {cols} (閾値 {thr})"] = ov[-1] if ov else None
same = all(np.array_equal(coords["A"], coords[a]) for a in ARMS)
if not same: OUT["problems"].append("3 本の格子の座標が同じでない")
print("座標が同じ:", same)
for k, v in OUT["quasisteady"].items(): print(f"  {k}: {v}")
OUT["recheck"] = "前提を満たす (判定は判別不能のまま)" if not OUT["problems"] else "前提を満たさない: " + "; ".join(OUT["problems"])
print("== 確かめ直し:", OUT["recheck"])
(HERE / "_band_ab/cold_pair/fz_recheck.json").write_text(json.dumps(OUT, indent=1, ensure_ascii=False))
