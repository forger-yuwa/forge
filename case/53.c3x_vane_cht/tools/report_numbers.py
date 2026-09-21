#!/usr/bin/env python3
r"""報告に載せる数字を run から一括で出す (領域別 $h$ 偏差・帯内率・準定常判定・収束判定・NaN)。

各 run について `compare_h.py` (領域別統計)、`h_series.py` + `check_quasisteady.py --series-csv`
(量ごとの準定常)、`check_convergence.py` (VERDICT) を回し、1 行ずつ表にする。
報告の数字はここから転記し、手で計算し直さない。

usage: python3 case/53.c3x_vane_cht/tools/report_numbers.py [--flux iface_q_eff] c3x:run_XXXX markii:run_YYYY ...
"""
import argparse, re, subprocess, sys, glob
from pathlib import Path
import numpy as np, h5py

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
T = ROOT / "solver_density_cuda/tools"
CASE = {"c3x": ("case/53.c3x_vane_cht", "run108"), "markii": ("case/54.markii_vane_cht", "run42")}


def sh(*a):
    return subprocess.run([sys.executable, *map(str, a)], capture_output=True, text=True).stdout


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("runs", nargs="+"); ap.add_argument("--flux", default="iface_q_eff")
    a = ap.parse_args()
    for spec in a.runs:
        vane, run = spec.split(":"); base, key = CASE[vane]; rd = ROOT / base / run
        out = sh(HERE / "compare_h.py", rd, "--run", key, "--flux", a.flux)
        reg = re.findall(r"^\s+(PS|SS laminar[^n]*|SS post-transition|all)\s+n=\s*(\d+)\s+bias\s+([+-][\d.]+)%\s+rms\s+([\d.]+)%\s+within report unc\.\s+([\d.]+)%", out, re.M)
        sh(HERE / "h_series.py", rd, "--run", key, "--flux", a.flux)
        qs = sh(T / "check_quasisteady.py", rd, "--series-csv", rd / f"h_series_{a.flux}.csv", "--series-cols", "PS,SS_lam,SS_post,all")
        ver = dict(re.findall(r"^\s+(PS|SS_lam|SS_post|all)\s+:.*?(STEADY|DRIFTING|OSCILLATING|TRANSIENT-UNSETTLED)", qs, re.M))
        cv = sh(T / "check_convergence.py", rd); m = re.search(r"->\s+(.*?) ===", cv); conv = m.group(1) if m else "?"
        fn = sorted(glob.glob(str(rd / "res_[0-9]*.h5")), key=lambda p: int(re.findall(r"(\d+)\.h5", p)[0]))[-1]
        with h5py.File(fn) as f:
            bad = [k for k in f["VALUE"] if not np.isfinite(np.array(f["VALUE"][k])).all()]
        print(f"\n## {base}/{run}   [{a.flux}]   convergence: {conv.split(' — ')[0]}   NaN fields: {bad or 'none'}")
        names = {"PS": "PS", "SS laminar": "SS_lam", "SS post-transition": "SS_post", "all": "all"}
        for r in reg:
            nm = r[0].strip().split(" (")[0]
            print(f"   {nm:<20} n={r[1]:>3}  bias {r[2]:>6} %  rms {r[3]:>5} %  in-band {r[4]:>5} %   quasi-steady: {ver.get(names.get(nm, nm), '?')}")


if __name__ == "__main__":
    main()
