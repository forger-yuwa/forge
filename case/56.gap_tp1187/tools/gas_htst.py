#!/usr/bin/env python3
"""case/56 の気体 — Langley 8-ft HTST のメタン-空気燃焼ガス。

case/50 (TN D-5908, **同じ 8-ft HTST**) 向けに作った `CombustionProducts` をそのまま使う。
当量比は LHV 収支から Tt,c で逆算される。TP-1187 の Tt,c は 1707-1935 K。
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "50.deep_cavity_wieting_m7" / "tools"))
from gas_model import CombustionProducts                        # noqa: E402


def htst(Tt_c: float) -> CombustionProducts:
    return CombustionProducts(Tt_c)


if __name__ == "__main__":
    import json
    cond = json.loads((HERE.parent / "conditions.json").read_text(encoding="utf-8"))
    print(f"{'run':>4} {'Tt_c [K]':>9} {'phi':>7} {'R':>8} {'cp(Tt)':>9} {'gamma(Tt)':>10}")
    seen = set()
    for s in cond["series"]:
        Tt = s["Tt_c_K"]
        if Tt in seen:
            continue
        seen.add(Tt)
        g = htst(Tt)
        print(f"{s['run']:4d} {Tt:9.0f} {g.phi:7.4f} {g.R:8.2f} {g.cp(Tt):9.2f} {g.gamma(Tt):10.4f}")
