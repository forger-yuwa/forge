#!/usr/bin/env python3
r"""③ベル (M) L100: 壁圧差の分解 (`corrected_wall_pressure_split.py`) を NS の本段のスナップショット (res_0 〜 res_60000、21 枚) で順に出し、
指標の時系列を `check_quasisteady.classify_series` (tail 0.4・drift 0.005・osc 0.01・min 20 — δ_r の M3 と同じ線) で判定する
(codex (diagnose) 2026-10-13 の「圧力の指標の時系列に check_quasisteady を当てる、未定常なら原因の判定を保留」)。
Euler 側は最終場のまま。B で回す想定 (途中の場は B にだけある)。出力 corrected_wall_pressure_split_series.json。"""
import json
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "solver_density_cuda" / "tools"))
import check_quasisteady as cq  # noqa: E402
import corrected_wall_pressure_split as sp  # noqa: E402
import judge_corrected_wall as j  # noqa: E402

KEYS = ("rms_e_w", "rms_e_s", "rms_e_near", "mean_e_w", "mean_e_s", "mean_e_near", "ratio_rms_s_over_w")


def main():
    pairs = sys.argv[1:] or ["Vext", "Vzero"]
    ed = np.load(j.EDGE, allow_pickle=True)
    out = {"lines": {"tail": 0.4, "drift": 0.005, "osc": 0.01, "min": 20}, "pairs": {}}
    for n in pairs:
        nsr, eur = sp.PAIRS[n]
        run = j.CASE / nsr
        snaps = sorted(run.glob("res_[0-9]*.h5"), key=lambda p: int(p.stem.split("_")[1]))
        eu = sp.load(j.CASE / eur)
        steps, rows = [], []
        for s in snaps:
            r = sp.split(sp.load(run, s.name), eu, ed["x"], ed["band_y_b"])
            steps.append(int(s.stem.split("_")[1])); rows.append({k: r[k] for k in KEYS})
            print(n, s.name, " ".join(f"{k} {r[k]:+.5f}" for k in KEYS), flush=True)
        st = np.array(steps, float)
        ver = {}
        for k in KEYS:
            v, det, _ = cq.classify_series(st, np.array([r_[k] for r_ in rows]), 0.4, 0.005, 0.01, 20)
            m = re.search(r"最終比 ([+-]?[0-9.eE+-]+) %", det or "")
            ver[k] = {"verdict": v, "asym_final_ratio_pct": float(m.group(1)) if m else None, "detail": det}
        out["pairs"][n] = {"steps": steps, "series": rows, "verdicts": ver}
        print(n, {k: (ver[k]["verdict"], ver[k]["asym_final_ratio_pct"]) for k in KEYS}, flush=True)
    (HERE / "corrected_wall_pressure_split_series.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
