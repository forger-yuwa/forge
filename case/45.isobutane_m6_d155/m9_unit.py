"""§6.11 の 1 step の専有時間・出力 1 回の費用・起動時間 (lineL、1000 step × 3 本ずつ、codex plan-2 M2・M3 の反映)。
各 run: 1 step の時間 = ログの表示の 101〜999 step の 1 step ごとの ms の平均、出力 1 回の費用 = (表示の step 1000 の ms) − (1 step の時間)
(実測: run_0349_timeLAY_LU_1 で表示の step 999 が 34.65 ms、1000 が 274.75 ms。出力は表示の 1000 に入る — codex plan-2 M2 の「999 に入る」は実測と合わない)、
起動と終了の時間 = (wall.txt の壁時計) − (全 step の ms の和)。全部が有限で、1 step・出力・起動は負でないこと、表示の step 1..1000 がちょうど 1 回ずつ、
投入前の GPU の他の計算プロセスが 0 本 (gpu_procs.txt が "0")、実効のモードが期待どおり (m9_modes.txt) であることを確かめ、外れたら判定不能 (終了コード 2)。
usage: python3 m9_unit.py → _band_ab/cold_pair/m9_unit.json"""
import json, math, re, statistics, sys
from collections import Counter
from pathlib import Path
HERE = Path(__file__).resolve().parent
OUT = HERE / "_band_ab" / "cold_pair" / "m9_unit.json"
OUT.unlink(missing_ok=True)
def bad(msg): print(msg + " — 判定不能"); sys.exit(2)
res = {}
for mode in ("P", "L0", "L5"):
    ms, oc, su = [], [], []
    for i in (1, 2, 3):
        r = HERE / f"run_035{4 + i}_m9unit_{mode}_{i}"
        if (r / "gpu_procs.txt").read_text().strip() != "0": bad(f"{r.name}: GPU が専有でない")
        if (r / "m9_modes.txt").read_text().strip() != "OK": bad(f"{r.name}: 実効のモードが期待と違う")
        rows = [(int(k), float(v)) for k, v in re.findall(r"^step\s+(\d+) \| ([0-9.eE+-]+|nan|inf) ms/step", (r / "forge_run.log").read_text(errors="replace"), re.M)]
        if Counter(k for k, _ in rows) != Counter(range(1, 1001)): bad(f"{r.name}: 表示の step 1..1000 がそろわない")
        d = dict(rows)
        if not all(math.isfinite(v) and v > 0 for v in d.values()): bad(f"{r.name}: ms が有限・正でない")
        unit = sum(d[k] for k in range(101, 1000)) / 899
        out_c = d[1000] - unit
        wall = float((r / "wall.txt").read_text())
        start = wall - sum(d.values()) / 1000.0
        if not (math.isfinite(wall) and out_c >= 0 and start >= 0): bad(f"{r.name}: 出力の費用 {out_c} か起動の時間 {start} が負・非有限")
        ms.append(unit); oc.append(out_c); su.append(start)
    res[mode] = {"ms_per_step": ms, "median_ms": statistics.median(ms), "spread_ms": max(ms) - min(ms),
                 "output_ms": oc, "median_output_ms": statistics.median(oc), "startup_s": su, "median_startup_s": statistics.median(su)}
OUT.write_text(json.dumps(res, indent=1, ensure_ascii=False)); print(json.dumps(res, indent=1, ensure_ascii=False))
