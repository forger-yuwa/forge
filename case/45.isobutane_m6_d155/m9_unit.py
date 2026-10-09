"""§6.9 の 1 step の専有時間と起動時間 (lineK、1000 step × 3 本ずつ)。各 run の 101〜999 step の 1 step ごとの ms の平均 (step がちょうど 1 回ずつ、有限・正)、
起動時間 = (cold_cfl.py run の壁時計、run の wall.txt) − (全 step の ms の和)。投入前の GPU の他の計算プロセス (gpu_procs.txt) が 0 本でない run があれば判定不能 (終了コード 2)。
usage: python3 m9_unit.py → _band_ab/cold_pair/m9_unit.json"""
import json, math, re, statistics, sys
from collections import Counter
from pathlib import Path
HERE = Path(__file__).resolve().parent
OUT = HERE / "_band_ab" / "cold_pair" / "m9_unit.json"
OUT.unlink(missing_ok=True)
res = {}
for mode in ("P", "L0", "L5"):
    ms, st = [], []
    for i in (1, 2, 3):
        r = HERE / f"run_035{4 + i}_m9unit_{mode}_{i}"
        if (r / "gpu_procs.txt").read_text().strip() != "0": print(f"{r.name}: GPU が専有でない — 判定不能"); sys.exit(2)
        rows = [(int(k), float(v)) for k, v in re.findall(r"^step\s+(\d+) \| ([0-9.eE+-]+|nan|inf) ms/step", (r / "forge_run.log").read_text(errors="replace"), re.M)]
        sel = [v for k, v in rows if 101 <= k <= 999]
        if Counter(k for k, _ in rows if 101 <= k <= 999) != Counter(range(101, 1000)) or not all(math.isfinite(v) and v > 0 for v in sel):
            print(f"{r.name}: step 101..999 がそろわない — 判定不能"); sys.exit(2)
        ms.append(sum(sel) / len(sel))
        st.append(float((r / "wall.txt").read_text()) - sum(v for _, v in rows) / 1000.0)
    res[mode] = {"ms_per_step": ms, "median_ms": statistics.median(ms), "startup_s": st, "median_startup_s": statistics.median(st)}
OUT.write_text(json.dumps(res, indent=1, ensure_ascii=False)); print(json.dumps(res, indent=1, ensure_ascii=False))
