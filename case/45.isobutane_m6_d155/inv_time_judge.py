"""§6.2 (3) の判定: 各 run の 101〜999 step の 1 step ごとの ms の平均 (forge_run.log の "step N | X ms/step" の行)、交互の 3 組の差 (INV − LU)。
usage: python3 inv_time_judge.py LU1 INV1 LU2 INV2 LU3 INV3 → 標準出力と _band_ab/cold_pair/inv_time_judge.json"""
import json, re, statistics, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
def mean_ms(run):
    v = [float(m.group(2)) for m in re.finditer(r"^step\s+(\d+) \| ([0-9.]+) ms/step", (HERE / run / "forge_run.log").read_text(errors="replace"), re.M) if 101 <= int(m.group(1)) <= 999]
    if len(v) != 899: raise SystemExit(f"{run}: 101〜999 step の行が {len(v)} 本 (899 本のはず) — 判定不能")
    return sum(v) / len(v)
runs = sys.argv[1:]; assert len(runs) == 6
lu = [mean_ms(r) for r in runs[0::2]]; inv = [mean_ms(r) for r in runs[1::2]]
d = [i - l for i, l in zip(inv, lu)]; spread = max(lu) - min(lu); med = statistics.median(d)
verdict = ("速い" if all(x < 0 for x in d) and abs(med) > spread else "遅い" if all(x > 0 for x in d) and abs(med) > spread else "判別不能")
out = {"runs": runs, "LU_ms": lu, "INV_ms": inv, "pair_diff_ms": d, "median_diff_ms": med, "LU_spread_ms": spread, "verdict": verdict}
(HERE / "_band_ab" / "cold_pair" / "inv_time_judge.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
print(json.dumps(out, indent=1, ensure_ascii=False))
