#!/usr/bin/env python3
"""前駆 2D 平板の解から、境界層の積分量と 3D 入口プロファイル CSV を作る (plan §4.3)。

usage:
  # station 掃引 (delta_99 が目標になる x を探す)
  python3 tools/extract_profile.py scan --run precursor/run_0001_precursor_m5
  # 入口 CSV を書く (station は --x か、目標 delta から自動)
  python3 tools/extract_profile.py csv --run precursor/run_0001_precursor_m5 --out RUN/inlet_profile_1.csv
  # 準定常判定用の時系列 (delta*, theta を全スナップショットで)
  python3 tools/extract_profile.py series --run precursor/run_0001_precursor_m5 --x 0.12

**CSV は自分で直書きする** (plan §4.3, codex 2026-09-19 plan-2 M3):
`gen_inlet_profile.py gen --table` は `Ps` を出力列から落とすため (Tt/M 換算経路でしか
再生成されない)、そのまま渡すと圧力分布が BC の一様値に戻る。列は
`z ro Ux Uy Uz Ps k omega` で、**2D の壁法線速度 V を 3D の Uz** に写す (発達 BL では V≠0)。
"""
import argparse
import json
import sys
from pathlib import Path

import h5py
import numpy as np


def _trapz(y, x):
    """台形積分。**numpy 2 で `np.trapz` が削除された**ので互換に包む (AWS の numpy で落ちた)。"""
    f = getattr(np, "trapezoid", None) or getattr(np, "trapz")
    return f(y, x)

HERE = Path(__file__).resolve().parent
CASE = HERE.parent
sys.path.insert(0, str(CASE))
from setup import load as load_conditions  # noqa: E402


def snapshots(run):
    return sorted(Path(run).glob("res_[0-9]*.h5"), key=lambda p: int(p.stem.split("_")[1]))


def read(res):
    with h5py.File(res, "r") as f:
        c = np.array(f["MESH/COORD"]).reshape(-1, 3)
        v = {k: np.array(f["VALUE/" + k]) for k in
             ("ro", "Ux", "Uy", "Uz", "P", "T", "k", "omega", "vis_lam", "wall_dist")}
    return c, v


def column(c, v, x_target):
    """x=x_target に最も近いノード列を壁から順に返す (2D 平面メッシュ: y が壁法線)。"""
    xs = np.unique(c[:, 0])
    x0 = xs[np.argmin(np.abs(xs - x_target))]
    idx = np.where(np.abs(c[:, 0] - x0) <= 1e-9 * max(1.0, abs(x0)))[0]
    if idx.size == 0:
        raise SystemExit("station x=%g のノード列が空" % x0)
    idx = idx[np.argsort(c[idx, 1])]
    return x0, idx


def bl_quantities(c, v, idx):
    """delta_99 / delta* / theta / Re_theta / tau_w / edge 値。
    delta_99 は u/u_e、積分厚さは ro*u で定義 (case/48 tools と同じ)。"""
    y = c[idx, 1]
    u, ro, T = v["Ux"][idx], v["ro"][idx], v["T"][idx]
    ru = ro * u
    # 外縁: ru の**外側 30% の中央値**を基準にし、壁から見て最初に 0.995 を超える点。
    # (ru/ru[-1] に searchsorted を使うと、前縁波の反射などで ru が単調でないとき破綻する)
    ru_ref = float(np.median(ru[int(0.7 * len(ru)):]))
    above = np.where(ru >= 0.995 * ru_ref)[0]
    i_e = int(above[0]) if above.size else len(y) - 1
    i_e = min(max(i_e, 3), len(y) - 1)
    ue, roe_ = u[i_e], ro[i_e]
    d99 = np.interp(0.99 * ue, u[:i_e + 1], y[:i_e + 1])
    f1 = 1.0 - ro * u / (roe_ * ue)
    f2 = (ro * u) / (roe_ * ue) * (1.0 - u / ue)
    ds = _trapz(f1[:i_e + 1], y[:i_e + 1])
    th = _trapz(f2[:i_e + 1], y[:i_e + 1])
    mu_e = v["vis_lam"][idx][i_e]
    tau_w = v["vis_lam"][idx][0] * (u[1] - u[0]) / (y[1] - y[0])
    return dict(delta99=float(d99), delta_star=float(ds), theta=float(th),
                Re_theta=float(roe_ * ue * th / mu_e), tau_w=float(tau_w),
                Tw=float(T[0]), ue=float(ue), roe=float(roe_), i_edge=int(i_e),
                n_in_bl=int(np.sum(y <= d99)))


def cmd_scan(a):
    res = snapshots(a.run)[-1]
    c, v = read(res)
    d = load_conditions()
    tgt = d["delta_target_m"]
    xs = np.unique(c[:, 0])
    xs = xs[xs > 0.005]
    probe = xs[:: max(1, len(xs) // 40)]
    print("res=%s  target delta99=%.3f mm" % (res.name, tgt * 1e3))
    print("%10s %10s %10s %10s %10s %10s %8s" %
          ("x[mm]", "d99[mm]", "d*[mm]", "th[mm]", "Re_th", "tau_w", "pts<d99"))
    best, rows = None, []
    for x in probe:
        x0, idx = column(c, v, x)
        q = bl_quantities(c, v, idx)
        rows.append((x0, q))
        print("%10.2f %10.4f %10.4f %10.4f %10.0f %10.1f %8d"
              % (x0 * 1e3, q["delta99"] * 1e3, q["delta_star"] * 1e3, q["theta"] * 1e3,
                 q["Re_theta"], q["tau_w"], q["n_in_bl"]))
        if best is None or abs(q["delta99"] - tgt) < abs(best[1]["delta99"] - tgt):
            best = (x0, q)
    print("\n-> delta99=%.3f mm に最も近い station: x = %.2f mm (delta99 %.3f mm, Re_theta %.0f, Tw %.1f K)"
          % (tgt * 1e3, best[0] * 1e3, best[1]["delta99"] * 1e3, best[1]["Re_theta"], best[1]["Tw"]))
    print("   (断熱平板の回復温度: CPG %.1f K / TP %.1f K)" % (d["Taw_cpg"], d.get("Taw_tp", float("nan"))))
    return best


def cmd_csv(a):
    res = snapshots(a.run)[-1]
    c, v = read(res)
    d = load_conditions()
    x_target = a.x if a.x is not None else cmd_scan(a)[0]
    x0, idx = column(c, v, x_target)
    q = bl_quantities(c, v, idx)
    y = c[idx, 1]
    # 3D 入口は z が壁法線。2D の (x,y) -> 3D の (x,z)、V(壁法線) -> Uz、スパン Uy=0。
    cols = {"z": y, "ro": v["ro"][idx], "Ux": v["Ux"][idx],
            "Uy": np.zeros_like(y), "Uz": v["Uy"][idx],
            "Ps": v["P"][idx], "k": v["k"][idx], "omega": v["omega"][idx]}
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    hdr = ["z", "ro", "Ux", "Uy", "Uz", "Ps", "k", "omega"]
    with open(out, "w") as f:
        f.write(" ".join(hdr) + "\n")
        for i in range(len(y)):
            f.write(" ".join("%.9g" % cols[h][i] for h in hdr) + "\n")
    meta = dict(source=str(res), station_x_m=float(x0), n_rows=int(len(y)), **q)
    Path(str(out) + ".json").write_text(json.dumps(meta, indent=2))
    print("wrote %s (%d rows) at x=%.2f mm" % (out, len(y), x0 * 1e3))
    for h in hdr:
        print("  %-6s %12.6g .. %12.6g" % (h, cols[h].min(), cols[h].max()))
    print("  delta99 %.4f mm / delta* %.4f mm / theta %.4f mm / Re_theta %.0f / Tw %.1f K"
          % (q["delta99"] * 1e3, q["delta_star"] * 1e3, q["theta"] * 1e3, q["Re_theta"], q["Tw"]))


def cmd_series(a):
    rows = []
    for res in snapshots(a.run):
        c, v = read(res)
        _, idx = column(c, v, a.x)
        q = bl_quantities(c, v, idx)
        rows.append((int(res.stem.split("_")[1]), q))
    out = Path(a.out or (Path(a.run) / "bl_series.csv"))
    keys = ["delta99", "delta_star", "theta", "Re_theta", "tau_w", "Tw"]
    with open(out, "w") as f:
        f.write("step," + ",".join(keys) + "\n")
        for st, q in rows:
            f.write("%d," % st + ",".join("%.9g" % q[k] for k in keys) + "\n")
    print("wrote", out, "(%d snapshots at x=%.1f mm)" % (len(rows), a.x * 1e3))
    print("  use: check_quasisteady.py --series-csv %s --series-cols %s" % (out, ",".join(keys)))


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("scan", "csv", "series"):
        p = sub.add_parser(name)
        p.add_argument("--run", required=True)
        p.add_argument("--x", type=float, default=None, help="station [m]")
        p.add_argument("--out", default=None)
    a = ap.parse_args()
    if a.cmd == "scan":
        cmd_scan(a)
    elif a.cmd == "csv":
        if not a.out:
            raise SystemExit("--out が要る")
        cmd_csv(a)
    else:
        if a.x is None:
            raise SystemExit("--x が要る")
        cmd_series(a)


if __name__ == "__main__":
    main()
