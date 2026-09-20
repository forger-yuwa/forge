#!/usr/bin/env python3
r"""CHT 壁の界面量を**スナップショット時系列 CSV** にする (G-if / 準定常判定の入力)。

壁面ダンプ `res_<physName>_<physID>_<step>.h5` を全ステップ読み、次を 1 行 / 1 スナップショットで書く。
`check_quasisteady.py --series-csv <out> --series-cols Tw_mean,Tw_max,q_total` にそのまま渡せる。

| 列 | 意味 |
| --- | --- |
| `step` | ステップ |
| `Tw_mean` / `Tw_min` / `Tw_max` | 壁温 [K] (面積重み平均 / 最小 / 最大) |
| `Tw_at_x` | 指定 x に最も近い節点の壁温 [K] (`--x`) |
| `q_total` | 流体側の界面熱量 $\sum q\,A$ [W] (平面 2D は W/m) |
| `q_solid_total` | 固体側 $\sum (T_w-T_b)/R_{tot}\,A$ (`--Rtot`, `--Tb` を与えたときのみ) |
| `imbalance_max` | **局所の界面不釣合い** $\max|q_f-q_s|$ [W/m²] (G-if の絶対条件) |
| `imbalance_rel` | 同 / $\max|q_f|$ (G-if の相対条件) |

usage:
  python3 cht_wall_series.py <run_dir> --phys-id 4 [--phys-name wall] [--flux q_compact]
      [--Rtot 4.60829e-3 --Tb 300] [--x 0.5] [-o wall_series.csv]
"""
import argparse
import glob
import os
import re
import sys

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from cht_loop import read_wall_dump   # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("--phys-id", type=int, required=True)
    ap.add_argument("--phys-name", default="wall")
    ap.add_argument("--flux", default="q_compact", choices=["q_compact", "q_recon", "q_2nd"])
    ap.add_argument("--Rtot", type=float, default=None, help="固体の全抵抗 [m2K/W] (与えると固体側と比較する)")
    ap.add_argument("--Tb", type=float, default=300.0)
    ap.add_argument("--x", type=float, default=None, help="この x に最も近い節点の壁温を列に出す")
    ap.add_argument("-o", "--out", default=None)
    a = ap.parse_args()

    pat = re.compile(rf"^res_{re.escape(a.phys_name)}_{a.phys_id}_(\d+)\.h5$")
    files = []
    for p in glob.glob(os.path.join(a.run_dir, f"res_{a.phys_name}_{a.phys_id}_*.h5")):
        m = pat.match(os.path.basename(p))
        if m:
            files.append((int(m.group(1)), p))
    if not files:
        sys.exit(f"no wall dumps res_{a.phys_name}_{a.phys_id}_*.h5 in {a.run_dir}")
    files.sort()

    out = a.out or os.path.join(a.run_dir, "wall_series.csv")
    cols = ["step", "Tw_mean", "Tw_min", "Tw_max", "q_total"]
    if a.x is not None:
        cols.insert(4, "Tw_at_x")
    if a.Rtot is not None:
        cols += ["q_solid_total", "imbalance_max", "imbalance_rel"]

    rows = []
    for step, path in files:
        coords, faces, vals = read_wall_dump(path)
        # 面積重み: 線要素 (2D) / 三角・四角 (3D) を長さ・面積で分配
        w = np.zeros(len(coords))
        for f in faces:
            p = coords[f]
            if len(f) == 2:
                L = np.linalg.norm(p[1] - p[0]); w[f] += 0.5 * L
            else:
                for i in range(1, len(f) - 1):
                    tri = [f[0], f[i], f[i + 1]]
                    ar = 0.5 * np.linalg.norm(np.cross(coords[tri[1]] - coords[tri[0]],
                                                       coords[tri[2]] - coords[tri[0]]))
                    w[tri] += ar / 3.0
        Tw = np.array(vals["Ts"], float)
        q = np.array(vals["iface_" + a.flux], float)
        r = {"step": step,
             "Tw_mean": float(np.sum(Tw * w) / np.sum(w)),
             "Tw_min": float(Tw.min()), "Tw_max": float(Tw.max()),
             "q_total": float(np.sum(q * w))}
        if a.x is not None:
            i = int(np.argmin(np.abs(coords[:, 0] - a.x)))
            r["Tw_at_x"] = float(Tw[i])
        if a.Rtot is not None:
            qs = (Tw - a.Tb) / a.Rtot
            r["q_solid_total"] = float(np.sum(qs * w))
            r["imbalance_max"] = float(np.max(np.abs(q - qs)))
            r["imbalance_rel"] = float(np.max(np.abs(q - qs)) / max(np.max(np.abs(q)), 1e-30))
        rows.append(r)

    with open(out, "w") as f:
        f.write(",".join(cols) + "\n")
        for r in rows:
            f.write(",".join(f"{r[c]:.8g}" for c in cols) + "\n")
    print(f"[cht_wall_series] {len(rows)} snapshots -> {out}  (flux = iface_{a.flux})")
    for r in rows[-3:]:
        print("   " + "  ".join(f"{c}={r[c]:.6g}" for c in cols))
    print("\n判定は: python3 check_quasisteady.py --series-csv %s --series-cols %s"
          % (out, ",".join(c for c in cols if c != "step")))


if __name__ == "__main__":
    sys.exit(main())
