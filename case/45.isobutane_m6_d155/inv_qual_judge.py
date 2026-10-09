"""§6.2 (4) の判定 (予備の選別): cold_series の時系列から |2πΣres_ro| ≤ 1.0 kg/s が step 0 を除く出力で 2 回続けて成り立つ最初の区間を探し、
手前の出力との線形補間の step と、そこまでの累積壁時計 (forge_run.log の各 step の ms の和、端数は補間) を出す。
usage: python3 inv_qual_judge.py <LU run> <INV run> → 標準出力と _band_ab/cold_pair/inv_qual_judge.json"""
import json, re, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent; OUTD = HERE / "_band_ab" / "cold_pair"
def reach(run):
    rows = [r for r in json.loads((OUTD / f"series_{run}.json").read_text())["rows"] if r["step"] > 0 and r.get("deficit") is not None]
    rows.sort(key=lambda r: r["step"])
    ms = {int(m.group(1)): float(m.group(2)) for m in re.finditer(r"^step\s+(\d+) \| ([0-9.]+) ms/step", (HERE / run / "forge_run.log").read_text(errors="replace"), re.M)}
    rec = {"run": run, "series": [(r["step"], r["deficit"]) for r in rows], "nonfinite": sum(r.get("nonfinite", 0) for r in rows)}
    for i in range(len(rows) - 1):
        a, b = abs(rows[i]["deficit"]), abs(rows[i + 1]["deficit"])
        if a <= 1.0 and b <= 1.0:
            if i == 0:
                st = float(rows[0]["step"]); rec["note"] = "最初の出力ですでに 1.0 以下 (補間の手前が無い)"
            else:
                s0, d0, s1, d1 = rows[i - 1]["step"], abs(rows[i - 1]["deficit"]), rows[i]["step"], a
                st = s0 + (d0 - 1.0) / (d0 - d1) * (s1 - s0) if d0 != d1 else float(s1)
            n = int(st); frac = st - n
            wall = sum(ms.get(k, 0.0) for k in range(1, n + 1)) + frac * ms.get(n + 1, 0.0)
            rec.update(reach_step=st, wall_ms=wall); return rec
    rec.update(reach_step=None, wall_ms=None, note="10000 step で到達しない (打ち切り)"); return rec
lu, inv = reach(sys.argv[1]), reach(sys.argv[2])
out = {"LU": lu, "INV": inv}
if lu["reach_step"] and inv["reach_step"]:
    out["reach_step_rel_diff"] = (inv["reach_step"] - lu["reach_step"]) / lu["reach_step"]
    out["wall_ratio_INV_over_LU"] = inv["wall_ms"] / lu["wall_ms"]
    out["verdict"] = ("判別不能 (到達の差が 10 % を超える)" if abs(out["reach_step_rel_diff"]) > 0.10
                      else ("逆行列の累積壁時計が短い" if inv["wall_ms"] < lu["wall_ms"] else "逆行列の累積壁時計が短くない"))
else:
    out["verdict"] = "打ち切り (片方または両方が未到達) — 採用の候補の判定は保留"
(OUTD / "inv_qual_judge.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
print(json.dumps({k: v for k, v in out.items() if k not in ("LU", "INV")} | {"LU_reach": lu["reach_step"], "INV_reach": inv["reach_step"], "LU_wall_s": (lu["wall_ms"] or 0) / 1e3, "INV_wall_s": (inv["wall_ms"] or 0) / 1e3}, indent=1, ensure_ascii=False))
