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


def crossing_mach(run, step, W, r, zband):
    """交差部 (縦すきまが横すきまに入る所) の流れ状態。

    `diagnostician` (2026-09-26) の指摘: 腕 A の q_w スパイクは**交差部が遷音速化した
    スナップショット**と一致していた。生産 run では main 段の最初の 2-3 枚でこれを
    基準メッシュと比べ、早期に打ち切り判断ができるようにする。

    返り値: (口の下 y∈[-6,-2.5] mm の M_max, y<-2.5 mm で M>1 の節点数)
    """
    fld = CASE / run / f"res_{step}.h5"
    msh = CASE / run / "mesh.h5"
    if not fld.exists() or not msh.exists():
        return float("nan"), -1
    with h5py.File(msh) as m:
        c = np.asarray(m["MESH/COORD"], dtype=float).reshape(-1, 3)
    with h5py.File(fld) as h:
        if "VALUE/sonic" not in h:
            return float("nan"), -1
        a = np.asarray(h["VALUE/sonic"], dtype=float)
        u = np.sqrt(sum(np.asarray(h["VALUE/U" + k], dtype=float) ** 2 for k in "xyz"))
    M = u / np.maximum(a, 1e-30)
    col = (np.abs(c[:, 0]) < 0.5 * W + r) & (c[:, 2] <= zband)
    band = col & (c[:, 1] <= -r) & (c[:, 1] >= -6.0e-3)
    deep = col & (c[:, 1] < -r)
    return (float(M[band].max()) if band.any() else float("nan"),
            int((M[deep] > 1.0).sum()))


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
    ap.add_argument("--min-step", type=int, default=0,
                    help="この step 未満のスナップショットを捨てる。**段階起動では必須** "
                         "— 暖機段の出力が同じ step 番号で混ざり、定常性判定が別方程式の過渡を含む")
    ap.add_argument("--csv", default=None)
    a = ap.parse_args()
    xd = 0.5 * a.w

    ss = [(k, f) for k, f in snaps(a.run) if k >= a.min_step]
    if not ss:
        raise SystemExit(f"{a.run}: res_gap_6_*.h5 が無い (bcond gap の outputHDFflg を確認)")
    rows = []
    for step, f in ss:
        c, q, ps = forward_wall(f, xd, a.r, a.depth, a.zband, a.tol)
        if len(c) < 4:
            print(f"  step {step}: 前向き壁の節点が {len(c)} 個しか取れない"); continue
        # **z 線ごとに深さ補間してから帯で平均する**。帯内の全節点を y でソートして
        # np.interp に渡すと、同じ y に複数の z の値が重なって深さごとに 1 節点を
        # 拾うだけになり、時系列が桁で振れる (2026-09-25 に踏んだ。抽出アーチファクト)。
        zs = np.unique(np.round(c[:, 2], 9))
        prof = np.full((len(zs), len(DEPTH_CM)), np.nan)
        for iz, z0 in enumerate(zs):
            m = np.abs(c[:, 2] - z0) < 1e-9
            if m.sum() < 3:
                continue
            yy, qq = c[m, 1], q[m]
            o = np.argsort(yy)                          # 深い -> 浅い (y 昇順)
            prof[iz] = np.interp([-d * 1e-2 for d in DEPTH_CM], yy[o], qq[o])
        vals = np.nanmean(prof, axis=0)
        spread = np.nanstd(prof, axis=0)
        iz0 = int(np.argmin(np.abs(zs)))                # 縦すきま中心線に最も近い z 線
        vals0 = prof[iz0]
        mmax, nsup = crossing_mach(a.run, step, a.w, a.r, a.zband)
        rows.append([step] + list(vals) + list(vals0) + [float(ps.max()), mmax, nsup])
        if step == ss[-1][0]:
            print(f"[{a.run}] step {step}  前向き壁: z 線 {len(zs)} 本 x 節点 {len(c)} "
                  f"(z<={a.zband*1e3:.1f} mm)")
            for tc, d, v, sd, v0 in zip(TC, DEPTH_CM, vals, spread, vals0):
                print(f"   TC {tc}  深さ {d:.2f} cm   帯平均 q_w = {v:11.2f} +- {sd:9.2f}"
                      f"   z=0 線 {v0:11.2f} W/m^2")
            print(f"   前向き壁の最大静圧 = {ps.max():.2f} Pa")
            print(f"   交差部 (|x|<W/2+r, z<={a.zband*1e3:.1f} mm): 口下 y[-6,-2.5] mm の M_max "
                  f"= {mmax:.3f}、y<-2.5 mm で M>1 の節点 {nsup}")
    cols = (["step"] + [f"q_tc{t}" for t in TC]
            + [f"q0_tc{t}" for t in TC] + ["p_fwd_max", "M_cross_max", "n_supersonic"])
    if a.csv:
        out = Path(a.csv)
        np.savetxt(out, np.array(rows), delimiter=",", header=",".join(cols),
                   comments="", fmt="%.8g")
        print(f"   -> {out} ({len(rows)} スナップショット)")
        print(f"   定常性の判定: python3 solver_density_cuda/tools/check_quasisteady.py "
              f"--series-csv {out} --series-cols {','.join(cols[1:])}")


if __name__ == "__main__":
    main()
