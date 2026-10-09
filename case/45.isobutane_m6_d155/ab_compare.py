"""旧・新バイナリの場の差を、同じバイナリの再実行の差と同じ尺度で並べる (v0_compare.py の一般化、2026-10-09)。
run 名が run_NNNN_<tag>_<old|new>_s<step><接尾辞> の形の run を集め、step ごとに組を作る。
usage (AWS の case dir): python3 ab_compare.py <tag> → _band_ab/cold_pair/AB_<tag>.json と標準出力
"""
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import v0_compare as V  # noqa: E402

if __name__ == "__main__":
    tag = sys.argv[1]
    steps = sorted({int(m.group(1)) for p in HERE.glob(f"run_*_{tag}_*_s*")
                    if (m := re.search(r"_s(\d+)[a-z]*$", p.name))})
    out = {f"s{st}": V.group(f"{tag} {st} step", rf"run_\d{{4}}_{tag}_(?P<b>old|new)_s{st}[a-z]*", st) for st in steps}
    (HERE / "_band_ab" / "cold_pair" / f"AB_{tag}.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
    for k, g in out.items():
        print(k, json.dumps(g["summary"], ensure_ascii=False))
