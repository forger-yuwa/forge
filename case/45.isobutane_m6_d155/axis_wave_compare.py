"""試験部の軸 M の山の格子感度 (plan verification-m6-axis-wave-mesh-su2 §4.3)。

M は hypot(Ux, Uy) / sonic (VALUE/sonic が無ければ停止)。
- 副指標 b(η) = max_{x∈[60,80]} ΔM − median_{x∈[45,90]} ΔM [%pt]、ΔM = M/M_E − 1 (Euler 参照 run_0037)。
- 主指標 D̃(x, η) = D − median_{x∈[45,90]} D、D = 100·(M_B/M_A − 1) (NS どうしの直接差)。
- 時系列: 各 run の全 res_<n>.h5 について b(0)・ピーク位置を CSV に書く。

usage:
  python3 axis_wave_compare.py b    <run> [<run> ...]          # 各 run の最終場の b(η) とピーク位置
  python3 axis_wave_compare.py series <run> [...]              # res ごとの b(0) 時系列 → <run>/axis_wave_series.csv
  python3 axis_wave_compare.py diff <run_A> <run_B>            # 主指標 D̃ (B 対 A)
"""
import json
import os
import re
import sys
from pathlib import Path

import h5py
import numpy as np

C = Path(__file__).resolve().parent
EULER = "run_0037_euler_rt77p02"
ETAS = (0.0, 0.05, 0.10, 0.20)
XQ = np.linspace(10.0, 94.0, 841)


def res_files(rd):
    fs = [f for f in os.listdir(rd) if re.fullmatch(r"res_\d+\.h5", f) and f != "res_0.h5"]
    return sorted(fs, key=lambda f: int(re.findall(r"\d+", f)[0]))


def load(run, res=None):
    """(X, R, M) を (ni, nj) 形で返す。座標は r_t 単位。"""
    rd = C / run
    info = json.loads((rd / "prepare_info.json").read_text())
    S = info["scale_m"]
    ni = info["mesh"]["ni"]
    with h5py.File(rd / "nozzle.h5") as f:
        nc = f["/MESH/COORD"][:].reshape(-1, 3)
    nj = nc.shape[0] // ni
    res = res or res_files(rd)[-1]
    with h5py.File(rd / res) as f:
        if "sonic" not in f["/VALUE"]:
            raise SystemExit(f"{run}/{res}: VALUE/sonic が無い (M を作れない)")
        M = np.hypot(f["/VALUE/Ux"][:], f["/VALUE/Uy"][:]) / f["/VALUE/sonic"][:]
    X = (nc[:, 0] / S).reshape(ni, nj)
    R = (nc[:, 1] / S).reshape(ni, nj)
    Mm = M.reshape(ni, nj)
    if not (np.allclose(R[:, 0], 0.0) and np.all(np.diff(X[:, 0]) > 0)):
        raise SystemExit(f"{run}: 構造格子 (i*nj+j, j=0 が軸) を仮定できない")
    return X, R, Mm


def eta_profile(X, R, M, xq, eta):
    """各 station で η=r/r_w の行内線形補間 → x 方向線形補間。"""
    xa = X[:, 0]
    i = np.clip(np.searchsorted(xa, xq) - 1, 0, len(xa) - 2)
    w = (xq - xa[i]) / (xa[i + 1] - xa[i])
    out = np.empty(len(xq))
    for k in range(len(xq)):
        v = [np.interp(eta, R[ii] / R[ii, -1], M[ii]) for ii in (i[k], i[k] + 1)]
        out[k] = (1 - w[k]) * v[0] + w[k] * v[1]
    return out


def bump(d, lo=60, hi=80, blo=45, bhi=90):
    s = (XQ >= lo) & (XQ <= hi)
    b = (XQ >= blo) & (XQ <= bhi)
    k = np.argmax(d[s])
    return float(d[s][k] - np.median(d[b])), float(XQ[s][k])


_EULER = {}


def euler_profile(eta):
    if eta not in _EULER:
        _EULER[eta] = eta_profile(*load(EULER), XQ, eta)
    return _EULER[eta]


def b_of(run, res=None, etas=ETAS):
    F = load(run, res)
    out = {}
    for eta in etas:
        d = 100.0 * (eta_profile(*F, XQ, eta) / euler_profile(eta) - 1.0)
        out[eta] = bump(d) + bump(d, 25, 45, 20, 50)   # (b70, x70, b35, x35)
    return out


def main():
    mode, runs = sys.argv[1], sys.argv[2:]
    if mode == "b":
        for run in runs:
            for eta, (b70, x70, b35, x35) in b_of(run).items():
                print(f"{run} η={eta:.2f}: b70 {b70:.4f} %pt @x={x70:.1f}   b35 {b35:.4f} %pt @x={x35:.1f}")
    elif mode == "series":
        for run in runs:
            rows = []
            for r in res_files(C / run):
                step = int(re.findall(r"\d+", r)[0])
                v = b_of(run, r, etas=(0.0, 0.05))
                rows.append((step, *v[0.0][:2], *v[0.05][:2]))
            a = np.array(rows)
            np.savetxt(C / run / "axis_wave_series.csv", a, delimiter=",", comments="",
                       header="step,b0_pct,x0_peak,b005_pct,x005_peak", fmt="%.6g")
            tail = a[-5:, 1]
            print(f"{run}: n={len(a)}  b0 last={a[-1,1]:.4f} %pt @x={a[-1,2]:.1f}  "
                  f"tail5 range={tail.max()-tail.min():.4f} %pt")
    elif mode == "diff":
        A, B = load(runs[0]), load(runs[1])
        base = (XQ >= 45) & (XQ <= 90)
        win = (XQ >= 60) & (XQ <= 80)
        for eta in (0.0, 0.05, 0.10, 0.20):
            D = 100.0 * (eta_profile(*B, XQ, eta) / eta_profile(*A, XQ, eta) - 1.0)
            Dt = D - np.median(D[base])
            k = np.argmax(np.abs(Dt[win]))
            print(f"η={eta:.2f}: median D[45,90] {np.median(D[base]):+.4f} %   "
                  f"max|D̃|[60,80] {abs(Dt[win][k]):.4f} % @x={XQ[win][k]:.1f}")
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
