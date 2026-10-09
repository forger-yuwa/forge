"""性能の組の判定 (plan time_integration-line-implicit-speed §6.2 (3) と同じ規則): 各 run の 101〜999 step の 1 step ごとの ms の平均、
組ごとの差 (B − A) が 3 組とも負でその中央値の大きさが A の 3 本の幅を超える → 速い、3 組とも正で同じ条件 → 遅い、ほかは判別不能。
各 run のログで step 101..999 がちょうど 1 回ずつ現れなければ判定不能 (終了コード 2)。
usage: python3 time_pairs_judge.py <out.json> A1 B1 A2 B2 A3 B3"""
import json, re, statistics, sys
from collections import Counter
from pathlib import Path
HERE = Path(__file__).resolve().parent
def mean_ms(run):
    rows = [(int(m.group(1)), float(m.group(2))) for m in re.finditer(r"^step\s+(\d+) \| ([0-9.]+) ms/step", (HERE / run / "forge_run.log").read_text(errors="replace"), re.M)]
    sel = [(k, v) for k, v in rows if 101 <= k <= 999]
    if Counter(k for k, _ in sel) != Counter(range(101, 1000)):
        print(f"{run}: step 101..999 がちょうど 1 回ずつでない — 判定不能"); sys.exit(2)
    return sum(v for _, v in sel) / len(sel)
out_path = HERE / sys.argv[1]; runs = sys.argv[2:]; assert len(runs) == 6
out_path.unlink(missing_ok=True)
A = [mean_ms(r) for r in runs[0::2]]; Bv = [mean_ms(r) for r in runs[1::2]]
d = [b - a for a, b in zip(A, Bv)]; spread = max(A) - min(A); med = statistics.median(d)
verdict = "速い" if all(x < 0 for x in d) and abs(med) > spread else "遅い" if all(x > 0 for x in d) and abs(med) > spread else "判別不能"
out = {"runs": runs, "A_ms": A, "B_ms": Bv, "pair_diff_ms": d, "median_diff_ms": med, "A_spread_ms": spread, "verdict": verdict}
out_path.write_text(json.dumps(out, indent=1, ensure_ascii=False)); print(json.dumps(out, indent=1, ensure_ascii=False))
