#!/usr/bin/env python3
"""2D 平板 run の 1 断面を `gen_inlet_profile.py gen --table` 用のテーブルに落とす。

3D すきま run (`gen_mesh_gap3d.py`) の入口は領域入口 ($x=x_{in}$) にあり、照合点の
横すきま ($x=0$) はそこから `x_in` ぶん下流にある。**照合点で $\\delta^*$ が実測の
位置 II (1.88 m) と揃うように、平板の $1.88 - |x_{in}|$ の断面**を抜く。

使い方:
    python3 tools/extract_inlet_table.py --run run_0002_fp_t8_long --x 1.7076 \\
        --out run_XXXX/_inlet_table.txt [--steady]

`--steady` を付けると、最後の 3 スナップショットで $\\delta^*$ と $\\theta$ が
動いていないかも出す (残差プラトーは量の定常化を意味しないため)。
"""
import argparse
from pathlib import Path
import numpy as np
import h5py

CASE = Path(__file__).resolve().parents[1]
COLS = ["y", "ro", "Ux", "Uy", "Uz", "Ps", "k", "omega"]


def latest_fields(run):
    fs = [p for p in (CASE / run).glob("res_[0-9]*.h5")]
    return sorted(fs, key=lambda p: int(p.stem.split("_")[1]))


def column(run, res, x_target):
    with h5py.File(CASE / run / "mesh.h5") as m:
        c = np.asarray(m["MESH/COORD"], dtype=float).reshape(-1, 3)
    xs = np.unique(c[:, 0])
    x0 = xs[np.argmin(np.abs(xs - x_target))]
    sel = np.where(np.abs(c[:, 0] - x0) < 1e-12)[0]
    order = np.argsort(c[sel, 1])
    sel = sel[order]
    with h5py.File(res) as f:
        v = {k: np.asarray(f["VALUE/" + k], dtype=float)[sel] for k in
             ["ro", "Ux", "Uy", "Uz", "P", "k", "omega"]}
        ny = [k for k in f["VALUE"] if k.startswith("Y") and k[1:].isdigit()]
        for k in sorted(ny, key=lambda s: int(s[1:])):
            v[k] = np.asarray(f["VALUE/" + k], dtype=float)[sel]
    v["y"] = c[sel, 1]
    v["Ps"] = v.pop("P")
    return v, x0, sorted(ny, key=lambda s: int(s[1:]))


def bl_thickness(v, edge_frac=0.995):
    """圧縮性の δ* と θ。

    **積分の上端を境界層端で切る**。自由流まで積分すると、平板上の僅かな圧縮で
    $\\rho u$ が端の値をわずかに上回り、0.35 m ぶんの符号違いが欠損を打ち消して
    $\\delta^*\\approx0$ になる (2026-09-25 に実際そうなった)。
    """
    y, ro, u = v["y"], v["ro"], v["Ux"]
    ru = ro * u
    ie = int(np.argmax(ru >= edge_frac * ru.max()))
    ie = max(ie, 2)
    re, ue = ro[ie], u[ie]
    f = ru[:ie + 1] / (re * ue)
    ds = np.trapz(1.0 - f, y[:ie + 1])
    th = np.trapz(f * (1.0 - u[:ie + 1] / ue), y[:ie + 1])
    i99 = int(np.argmax(u >= 0.99 * ue))
    return ds, th, y[i99]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    ap.add_argument("--x", type=float, required=True, help="抽出断面 [m]")
    ap.add_argument("--out", required=True)
    ap.add_argument("--steady", action="store_true", help="末尾 3 枚で δ*/θ の定常性を見る")
    a = ap.parse_args()

    fs = latest_fields(a.run)
    if not fs:
        raise SystemExit(f"{a.run}: res_*.h5 が無い")
    v, x0, ycols = column(a.run, fs[-1], a.x)
    ds, th, d99 = bl_thickness(v)
    print(f"[抽出] {a.run}/{fs[-1].name}  要求 x={a.x:.4f} m -> 実際 x={x0:.6f} m  節点 {len(v['y'])}")
    print(f"        δ*={ds*100:.4f} cm  θ={th*100:.4f} cm  H={ds/th:.3f}  δ99={d99*100:.3f} cm")
    print(f"        端: ro={v['ro'][-1]:.6g} Ux={v['Ux'][-1]:.6g} Ps={v['Ps'][-1]:.6g} "
          f"k={v['k'][-1]:.6g} omega={v['omega'][-1]:.6g}")
    if a.steady:
        for p in fs[-3:]:
            vv, _, _ = column(a.run, p, a.x)
            d2, t2, _ = bl_thickness(vv)
            print(f"        [定常性] {p.name:20s} δ*={d2*100:.4f} cm  θ={t2*100:.4f} cm")

    cols = COLS + ycols
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    arr = np.column_stack([v[c] for c in cols])
    np.savetxt(out, arr, header=" ".join(cols), comments="", fmt="%.8e")
    print(f"        -> {out}  ({arr.shape[0]} 行 x {arr.shape[1]} 列)")


if __name__ == "__main__":
    main()
