"""§6.4 (2) の判定: 3 組 (従来・F32c・F32cs) の 101〜999 step の 1 step ごとの ms の平均と、腕ごとの組の差 (腕 − 従来)。
各 run のログで step 101..999 がちょうど 1 回ずつ現れることを確かめ、外れたら判定不能 (終了コード 2)。
usage: python3 f32_time_judge.py → _band_ab/cold_pair/f32_time_judge.json (run 名は run_034{1,2,3}_time3_{LU,F32c,F32cs}_{1,2,3})"""
import json, re, statistics, sys
from collections import Counter
from pathlib import Path
HERE = Path(__file__).resolve().parent
OUT = HERE / "_band_ab" / "cold_pair" / "f32_time_judge.json"
def mean_ms(run):
    rows = [(int(m.group(1)), float(m.group(2))) for m in re.finditer(r"^step\s+(\d+) \| ([0-9.]+) ms/step", (HERE / run / "forge_run.log").read_text(errors="replace"), re.M)]
    sel = [(k, v) for k, v in rows if 101 <= k <= 999]
    if Counter(k for k, _ in sel) != Counter(range(101, 1000)):
        print(f"{run}: step 101..999 がちょうど 1 回ずつでない — 判定不能"); sys.exit(2)
    return sum(v for _, v in sel) / len(sel)
OUT.unlink(missing_ok=True)
ms = {arm: [mean_ms(f"run_034{i}_time3_{arm}_{i}") for i in (1, 2, 3)] for arm in ("LU", "F32c", "F32cs")}
spread = max(ms["LU"]) - min(ms["LU"])
out = {"ms": ms, "LU_spread_ms": spread}
for arm in ("F32c", "F32cs"):
    d = [x - y for x, y in zip(ms[arm], ms["LU"])]; med = statistics.median(d)
    out[arm] = {"pair_diff_ms": d, "median_diff_ms": med,
                "verdict": "速い" if all(x < 0 for x in d) and abs(med) > spread else "遅い" if all(x > 0 for x in d) and abs(med) > spread else "判別不能"}
out["note"] = "F32c と F32cs の短縮は各変更の実効の短縮量。律速 (演算かメモリか) の解釈は仮説で、#10 は配列の変更そのものの A/B で決める (codex 2026-10-09 m6)"
OUT.write_text(json.dumps(out, indent=1, ensure_ascii=False))
print(json.dumps(out, indent=1, ensure_ascii=False))
