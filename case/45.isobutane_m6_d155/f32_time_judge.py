"""§6.4 (2) の判定: 3 組 (従来・F32c・F32cs) の 101〜999 step の 1 step ごとの ms の平均と、腕ごとの組の差 (腕 − 従来)。
usage: python3 f32_time_judge.py → _band_ab/cold_pair/f32_time_judge.json (run 名は run_034{1,2,3}_time3_{LU,F32c,F32cs}_{1,2,3})"""
import json, re, statistics
from pathlib import Path
HERE = Path(__file__).resolve().parent
def mean_ms(run):
    v = [float(m.group(2)) for m in re.finditer(r"^step\s+(\d+) \| ([0-9.]+) ms/step", (HERE / run / "forge_run.log").read_text(errors="replace"), re.M) if 101 <= int(m.group(1)) <= 999]
    if len(v) != 899: raise SystemExit(f"{run}: 101〜999 step の行が {len(v)} 本 — 判定不能")
    return sum(v) / len(v)
ms = {arm: [mean_ms(f"run_034{i}_time3_{arm}_{i}") for i in (1, 2, 3)] for arm in ("LU", "F32c", "F32cs")}
spread = max(ms["LU"]) - min(ms["LU"])
out = {"ms": ms, "LU_spread_ms": spread}
for arm in ("F32c", "F32cs"):
    d = [x - y for x, y in zip(ms[arm], ms["LU"])]; med = statistics.median(d)
    out[arm] = {"pair_diff_ms": d, "median_diff_ms": med,
                "verdict": "速い" if all(x < 0 for x in d) and abs(med) > spread else "遅い" if all(x > 0 for x in d) and abs(med) > spread else "判別不能"}
(HERE / "_band_ab" / "cold_pair" / "f32_time_judge.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
print(json.dumps(out, indent=1, ensure_ascii=False))
