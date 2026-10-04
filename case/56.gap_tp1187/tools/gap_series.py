#!/usr/bin/env python3
"""すきま深部の量を全スナップショットで出し、準定常かを判定する。

`gap_eval.py` は 1 枚しか見ないので「深部が落ち着いたか」を判定できない。
AGENTS.md の準定常確認 (check_quasisteady.py) に食わせる CSV を作る。
"""
import subprocess, sys
from pathlib import Path
import numpy as np
import h5py

CASE = Path(__file__).resolve().parents[1]
TOOLS = CASE.parents[1] / "solver_density_cuda" / "tools"
W = 0.18e-2
DEPTHS = [(1.39, "zW139"), (2.83, "zW283"), (4.22, "zW422"), (8.44, "zW844"), (14.11, "zW1411")]


def wall_q(res):
    """前向き壁の q_w(深さ)。**forge の壁出力 `res_gap_6_*.h5` の `qwall` を使う**。

    場から幾何的に片側差分を取る方式は、襟が床付近で傾いていて同一 y の内側ノードが
    無く壁点を拾えない (`gap_eval.py` の注記どおり。2026-09-20 に再度踏んだ)。
    """
    step = int(Path(res).stem.split("_")[1])
    wf = Path(res).parent / f"res_gap_6_{step}.h5"
    if not wf.exists():
        raise SystemExit(f"REFUSED: {wf.name} が無い (bcond の outputHDFflg を確認)")
    with h5py.File(wf) as f:
        wc = f["/MESH/COORD"][:].reshape(-1, 3)
        qw = -f["/VALUE/qwall"][:].astype(float)      # 流体→壁を正
    m = wc[:, 0] > 1e-9                                # 下流側 (前向き壁 = 計測面)
    d, q = -wc[m, 1], qw[m]
    o = np.argsort(d)
    return d[o], q[o]


def deep_stats(res):
    with h5py.File(res) as h:
        c = h["/MESH/COORD"][:].reshape(-1, 3)
        T = h["/VALUE/T"][:].astype(np.float32)
        U = np.hypot(h["/VALUE/Ux"][:].astype(float), h["/VALUE/Uy"][:].astype(float))
    deep = (c[:, 1] < 0) & (np.abs(c[:, 1]) > 4.0 * W)
    return float(np.ptp(T[deep])), float(np.sqrt(np.mean(U[deep] ** 2)))


def main():
    run = sys.argv[1]
    rd = CASE / run
    files = sorted(rd.glob("res_[0-9]*.h5"), key=lambda p: int(p.stem.split("_")[1]))
    if len(files) < 4:
        raise SystemExit(f"{run}: スナップショットが {len(files)} 枚しかない")
    rows, hdr = [], ["step"] + [n for _, n in DEPTHS] + ["T_span_deep_K", "U_rms_deep"]
    # 壁出力が無い step (初期場など) は飛ばす
    files = [f for f in files
             if (f.parent / f"res_gap_6_{int(f.stem.split('_')[1])}.h5").exists()]
    if len(files) < 4:
        raise SystemExit(f"{run}: 壁出力つきスナップショットが {len(files)} 枚しかない")
    for f in files:
        d, q = wall_q(f)
        ts, us = deep_stats(f)
        rows.append([int(f.stem.split("_")[1])]
                    + [float(np.interp(z * W, d, q)) for z, _ in DEPTHS] + [ts, us])
    a = np.array(rows)
    out = rd / "gap_series.csv"
    np.savetxt(out, a, delimiter=",", header=",".join(hdr), comments="")
    print(f"{run}: {len(files)} 枚  step {a[0,0]:.0f}..{a[-1,0]:.0f}  → {out}")
    r = subprocess.run([sys.executable, str(TOOLS / "check_quasisteady.py"),
                        "--series-csv", str(out), "--series-cols", ",".join(hdr[1:])],
                       capture_output=True, text=True)
    print(r.stdout.rstrip())
    if r.stderr.strip():
        print(r.stderr.rstrip())


if __name__ == "__main__":
    main()
