#!/usr/bin/env python3
"""case/48 の精度の A/B (plan architecture-float-state-double-geometry §6.21・§6.23) の判定量の時系列を、各 res_<N>.h5 から作る。

    python3 c48_prec_series.py <run_dir> [--steps 2000:48000:2000] [--coords-from <double の座標の h5>]
                               [--fix-je 0.6:76] [--out prec_series.csv]

判定量は `cooled_plate_eval.py` の定義をそのまま使う:
- x = 0.3 / 0.6 / 0.9 の θ・δ*・Cf・q_w (`station()`: 壁法線の列、ρu が 99.5 % に達する最初の点 je を縁とする台形積分、壁の量は 3 点の片側差分)
- CD = ∫Cf dx・HF = ∫q_w dx (`plate_integrals()`: x ∈ [0.002, 0.998] の 200 点)
記録として各断面の縁の点 je と縁の高さ y_edge も出す。
- `--coords-from` (§6.23、codex diagnose の Major 1): 各 res の `/MESH/COORD` (float のビルドでは float32) の代わりに、共通の double の座標を使う。
  res の座標が、共通の座標を res の型に丸めたものと全節点で一致することを確かめてから置き換える (一致しなければ終了コード 2)。
- `--fix-je X:J` (§6.23 の B): 断面 X の縁の点を J に固定する (縁の値 ρe・ue も同じ節点から取る)。ほかの断面は動的のまま。
  固定しない断面は `cooled_plate_eval.station()` をそのまま呼ぶ。固定する断面は同じ式を je だけ差し替えて計算する (`station_je`。je を固定しないときに
  `station()` と全量が一致することを `--self-test` で確かめる)。
指定した step の res が 1 つでも欠ける・判定量に非有限があると、CSV を書かずに終了コード 2 で止める (判定不能にするため)。"""
import argparse, csv, math, sys
from pathlib import Path
import numpy as np
import h5py
sys.path.insert(0, str(Path(__file__).resolve().parent))
import cooled_plate_eval as E  # noqa: E402

XS = (0.3, 0.6, 0.9)
QTY = ("theta", "dstar", "Cf", "qw")


def station_je(D, xs, je_fixed=None, edge_frac=0.995):
    """cooled_plate_eval.station() と同じ式で、縁の点 je だけを固定できるもの (θ・δ*・Cf・q_w・je・y_edge を返す)。"""
    x, y = D["x"], D["y"]
    wallx = np.unique(np.round(x[(np.abs(y) < 1e-9) & (x >= -1e-9)], 7))
    xv = wallx[np.argmin(np.abs(wallx - xs))]
    col = np.where(np.abs(x - xv) < 1e-7)[0]
    col = col[np.argsort(y[col])]
    yy = y[col]; u = D["u"][col]; ro = D["ro"][col]; T = D["T"][col]; mu = D["mu"][col]
    ywin = 0.6 * (xv + 0.1) * math.tan(math.asin(1.0 / E.M_INF))
    m = yy <= ywin
    rou = ro * u
    imax = int(np.argmax(rou[m]))
    je = int(np.argmax(rou[m] >= edge_frac * rou[m][imax])) if je_fixed is None else int(je_fixed)
    ue, roe = u[je], ro[je]
    sl = slice(0, je + 1)
    dstar = float(np.trapezoid(1.0 - rou[sl] / (roe * ue), yy[sl]))
    theta = float(np.trapezoid(rou[sl] / (roe * ue) * (1.0 - u[sl] / ue), yy[sl]))
    muw = mu[0]
    dudy = E.deriv_wall(yy[:3], u[:3]); dTdy = E.deriv_wall(yy[:3], T[:3])
    tau = muw * dudy
    qw = muw * E.CP / E.PR * dTdy
    return dict(theta=theta, dstar=dstar, Cf=float(tau / (0.5 * E.RO_INF * E.U_INF ** 2)), qw=float(qw), je=je, y_edge=float(yy[je]))


def load(path, coords64):
    D = E.load_forge(path)
    if coords64 is not None:
        with h5py.File(path, "r") as h: c = np.asarray(h["/MESH/COORD"][:]).reshape(-1, 3)
        if c.shape != coords64.shape or not np.array_equal(c, coords64.astype(c.dtype)):
            sys.exit(print(f"[c48_prec_series] {path}: 座標が共通の double の座標と対応しない") or 2)
        D["x"], D["y"] = coords64[:, 0].copy(), coords64[:, 1].copy()
    return D


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("run"); ap.add_argument("--steps", default="2000:48000:2000")
    ap.add_argument("--coords-from", default=None); ap.add_argument("--fix-je", default=None); ap.add_argument("--out", default="prec_series.csv")
    ap.add_argument("--self-test", action="store_true", help="station_je(je 固定なし) が station() と全量で一致するかを、最後の step の場で確かめて終わる")
    a = ap.parse_args()
    s0, s1, ds = (int(v) for v in a.steps.split(":"))
    run = Path(a.run); steps = list(range(s0, s1 + 1, ds))
    coords64 = None
    if a.coords_from:
        with h5py.File(a.coords_from, "r") as h: coords64 = np.asarray(h["/MESH/COORD"][:], np.float64).reshape(-1, 3)
    fix = {}
    if a.fix_je:
        for kv in a.fix_je.split(","): xk, jk = kv.split(":"); fix[float(xk)] = int(jk)
    if a.self_test:
        D = load(run / f"res_{s1}.h5", coords64); bad = []
        for x in XS:
            s_ref = E.station(D, x); s_new = station_je(D, x)
            for k in QTY + ("je", "y_edge"):
                if s_ref[k] != s_new[k]: bad.append((x, k, s_ref[k], s_new[k]))
        print(f"[c48_prec_series] self-test {run.name} step {s1}: {'一致' if not bad else '不一致 ' + str(bad)}")
        sys.exit(0 if not bad else 2)
    miss = [n for n in steps if not (run / f"res_{n}.h5").exists()]
    if miss: sys.exit(print(f"[c48_prec_series] {run.name}: res が無い step {miss[:5]}… ({len(miss)} 個)") or 2)
    cols = ["step"] + [f"{q}_{x}" for x in XS for q in QTY] + ["CD", "HF"] + [f"je_{x}" for x in XS] + [f"yedge_{x}" for x in XS]
    rows = []
    for n in steps:
        D = load(run / f"res_{n}.h5", coords64)
        r = {"step": n}
        for x in XS:
            st = station_je(D, x, fix[x]) if x in fix else E.station(D, x)
            for q in QTY: r[f"{q}_{x}"] = st[q]
            r[f"je_{x}"] = st["je"]; r[f"yedge_{x}"] = st["y_edge"]
        pi = E.plate_integrals(D); r["CD"] = pi["CD"]; r["HF"] = pi["HF"]
        bad = [k for k, v in r.items() if not math.isfinite(float(v))]
        if bad: sys.exit(print(f"[c48_prec_series] {run.name}: step {n} の {bad} が非有限") or 2)
        rows.append(r)
    with open(run / a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols); w.writeheader(); [w.writerow(r) for r in rows]
    print(f"[c48_prec_series] {run.name}: {len(rows)} 点 (座標 {'共通の double' if coords64 is not None else '各 res'}、je 固定 {fix or 'なし'}) → {run / a.out}")


if __name__ == "__main__":
    main()
