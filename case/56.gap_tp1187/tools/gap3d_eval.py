#!/usr/bin/env python3
"""3D すきま run の**前向き壁**の深さ分布を取り出す (A/B と D1 の照合に使う)。

前向き壁 = 横すきまの下流側の鉛直壁 ($x=+W/2$)。TP-1187 の熱電対 TC 87-92 は
縦すきま中心線から高々 6 mm 以内にあるので、**$z$ の帯で切って**深さ方向に並べる。

    python3 tools/gap3d_eval.py --run run_0049_ab3d_wall [--zband 0.006] [--csv out.csv]

出力: 各スナップショットについて深さ別の $q_w$ [W/m²] と前向き壁の最大静圧。
`--csv` を付けると `check_quasisteady.py --series-csv ... --series-cols ...` に
そのまま渡せる時系列 CSV を書く (派生量の定常性はツールで判定すること)。

**深さの定義**: ここでは上面 ($y=0$) からの鉛直距離を使う。TP-1187 の $s$ は
「タイル半径中心弧からの表面距離」なので、**実測との対応づけ (D1) は別途固定が要る**。
A/B は同じ定義で両腕を比べるだけなので、この定義のままで足りる。
"""
import argparse, glob, re
from pathlib import Path
import numpy as np
import h5py

CASE = Path(__file__).resolve().parents[1]
DEPTH_CM = [0.25, 0.51, 0.76, 1.52, 2.54, 3.81]      # TC 92,91,90,89,88,87
TC = [92, 91, 90, 89, 88, 87]


def snaps(run):
    fs = glob.glob(str(CASE / run / "res_gap_6_*.h5"))
    out = []
    for f in fs:
        m = re.search(r"res_gap_6_(\d+)\.h5$", f)
        if m:
            out.append((int(m.group(1)), f))
    return sorted(out)


def forward_wall(path, xd, r, D, zband, tol):
    with h5py.File(path) as h:
        c = np.asarray(h["MESH/COORD"], dtype=float).reshape(-1, 3)
        q = np.asarray(h["VALUE/qwall"], dtype=float)
        ps = np.asarray(h["VALUE/Ps"], dtype=float)
    sel = ((np.abs(c[:, 0] - xd) < tol) & (c[:, 1] < -r + tol) & (c[:, 1] > -D - tol)
           & (c[:, 2] <= zband))
    return c[sel], q[sel], ps[sel]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    ap.add_argument("--w", type=float, default=0.18e-2)
    ap.add_argument("--r", type=float, default=0.25e-2)
    ap.add_argument("--depth", type=float, default=6.35e-2)
    ap.add_argument("--zband", type=float, default=6.0e-3, help="z の帯 [m] (熱電対の横方向範囲)")
    ap.add_argument("--tol", type=float, default=1.0e-6)
    ap.add_argument("--csv", default=None)
    a = ap.parse_args()
    xd = 0.5 * a.w

    ss = snaps(a.run)
    if not ss:
        raise SystemExit(f"{a.run}: res_gap_6_*.h5 が無い (bcond gap の outputHDFflg を確認)")
    rows = []
    for step, f in ss:
        c, q, ps = forward_wall(f, xd, a.r, a.depth, a.zband, a.tol)
        if len(c) < 4:
            print(f"  step {step}: 前向き壁の節点が {len(c)} 個しか取れない"); continue
        o = np.argsort(-c[:, 1])                       # 浅い -> 深い
        y, qq = c[o, 1], q[o]
        vals = []
        for d in DEPTH_CM:
            vals.append(float(np.interp(-d * 1e-2, y[::-1], qq[::-1])))
        rows.append([step] + vals + [float(ps.max())])
        if step == ss[-1][0]:
            print(f"[{a.run}] step {step}  前向き壁の節点 {len(c)} (z<={a.zband*1e3:.1f} mm)")
            for tc, d, v in zip(TC, DEPTH_CM, vals):
                print(f"   TC {tc}  深さ {d:.2f} cm   q_w = {v:11.2f} W/m^2")
            print(f"   前向き壁の最大静圧 = {ps.max():.2f} Pa")
    cols = ["step"] + [f"q_tc{t}" for t in TC] + ["p_fwd_max"]
    if a.csv:
        out = Path(a.csv)
        np.savetxt(out, np.array(rows), delimiter=",", header=",".join(cols),
                   comments="", fmt="%.8g")
        print(f"   -> {out} ({len(rows)} スナップショット)")
        print(f"   定常性の判定: python3 solver_density_cuda/tools/check_quasisteady.py "
              f"--series-csv {out} --series-cols {','.join(cols[1:])}")


if __name__ == "__main__":
    main()
