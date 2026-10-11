"""V0 (plan time_integration-implicit-thermal-jacobian §6) の比較: 同じバイナリの再実行の差と、旧・新バイナリの差を同じ尺度で並べる。
各組の res_<N>.h5 の VALUE/* について、場の差の相対 RMS (‖a − b‖₂ / ‖a‖₂) と最大絶対差を出す。
  - 同じバイナリの組 (旧 × 旧、新 × 新) = 再実行の差 (atomicAdd の加算順)
  - 旧 × 新 の組 = バイナリの差 + 再実行の差
usage (AWS の case dir): python3 v0_compare.py → _band_ab/cold_pair/V0_repeat.json
"""
import itertools
import json
import re
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
FIELDS = ("ro", "roUx", "roUy", "roe", "roK", "roOmega", "P", "T")


def runs(pattern):
    out = []
    for p in sorted(HERE.glob("run_*")):
        m = re.fullmatch(pattern, p.name)
        if m:
            out.append((p, m.group("b")))
    return out


def load(run: Path, step: int):
    f = run / f"res_{step}.h5"
    if not f.is_file():
        return None
    with h5py.File(f, "r") as h:
        return {k: np.asarray(h["VALUE"][k][:], float) for k in FIELDS if k in h["VALUE"]}


def diff(a, b):
    return {k: {"rel_rms": float(np.linalg.norm(a[k] - b[k]) / max(np.linalg.norm(a[k]), 1e-300)),
                "max_abs": float(np.max(np.abs(a[k] - b[k]))), "bit_identical": bool(np.array_equal(a[k], b[k]))}
            for k in FIELDS if k in a and k in b}


def group(tag, pattern, step):
    rs = [(p, b, load(p, step)) for p, b in runs(pattern)]
    rs = [(p, b, d) for p, b, d in rs if d is not None]
    pairs = {"old_old": [], "new_new": [], "old_new": []}
    for (p1, b1, d1), (p2, b2, d2) in itertools.combinations(rs, 2):
        key = "old_old" if b1 == b2 == "old" else "new_new" if b1 == b2 == "new" else "old_new"
        pairs[key].append({"a": p1.name, "b": p2.name, "diff": diff(d1, d2)})
    summ = {}
    for key, lst in pairs.items():
        if lst:
            summ[key] = {f: {"rel_rms_max": max(x["diff"][f]["rel_rms"] for x in lst),
                             "rel_rms_min": min(x["diff"][f]["rel_rms"] for x in lst),
                             "n_bit_identical": sum(x["diff"][f]["bit_identical"] for x in lst), "n_pairs": len(lst)}
                         for f in FIELDS if f in lst[0]["diff"]}
    print(f"== {tag} (res_{step}, {len(rs)} 本: {', '.join(p.name for p, _, _ in rs)})")
    for key, s in summ.items():
        print(f"  {key:8s} " + "  ".join(f"{f} {v['rel_rms_min']:.1e}〜{v['rel_rms_max']:.1e} (一致 {v['n_bit_identical']}/{v['n_pairs']})"
                                         for f, v in s.items() if f in ("ro", "roe", "roOmega", "T")))
    return {"step": step, "runs": [p.name for p, _, _ in rs], "pairs": pairs, "summary": summ}


if __name__ == "__main__":
    out = {"pt0_20": group("point cfl 4・ISP 0・20 step", r"run_\d{4}_v0_pt0_(?P<b>old|new)(_[bc])?", 20),
           "dir0_20": group("directional cfl 4・ISP 0・20 step", r"run_\d{4}_v0_dir0_(?P<b>old|new)(_[bc])?", 20),
           "pt0_1": group("point cfl 4・ISP 0・1 step", r"run_\d{4}_v0s1_pt0_(?P<b>old|new)_[ab]", 1),
           "dir0_1": group("directional cfl 4・ISP 0・1 step", r"run_\d{4}_v0s1_dir0_(?P<b>old|new)_[ab]", 1)}
    (HERE / "_band_ab" / "cold_pair" / "V0_repeat.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
