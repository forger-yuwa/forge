#!/usr/bin/env python3
"""G14 本体 (acceptance.json の G14-E): 系列ごとに、事前登録した Re_H 点で R = St_CFD / St_exp を出す。

    python3 tools/cary_g14.py <run_dir> [<run_dir> ...] [--series-csv-dir DIR]

Re_H の作り方は tools/cary_energy_compare.py と同じ (前縁から積分、実験は遷移前の点の最小係数で前縁入熱を外挿)。
St は両者とも Re_H に対して線形補間する。x_b・遷移前の点・比較 Re_H は acceptance.json G14-E から読む
(結果を見てから変えない)。複数 run を渡すと群ごとの R(Tw/Tt) 表を出す。
"""
import argparse, glob, json, re
from pathlib import Path
import numpy as np
import h5py

CASE = Path(__file__).resolve().parents[1]
GAM, CP = 1.4, 1004.5


def mu_suth(T):
    return 1.716e-5 * (T / 273.0) ** 1.5 * (273.0 + 111.0) / (T + 111.0)


def cumtrapz(x, y):
    return np.concatenate([[0.0], np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(x))])


def wall_q(run, step):
    with h5py.File(Path(run) / f"res_wall_4_{step}.h5") as h:
        c = np.asarray(h["MESH/COORD"], float).reshape(-1, 3)
        q = -np.asarray(h["VALUE/qwall"], float)
    m = c[:, 0] >= 0.0
    x = np.round(c[m, 0], 9); ux = np.unique(x)
    return ux, np.array([q[m][x == u].mean() for u in ux])


def exp_curve(ser, su, pe, xs):
    T_inf, U, ro, M, Tw, Tt = su["T_inf"], su["U_inf"], su["ro_inf"], su["M"], su["Tw"], su["Tt"]
    rhoucp = ro * U * CP
    Taw = lambda r: T_inf * (1 + r * 0.5 * (GAM - 1) * M ** 2)
    st = np.array([np.nan if v is None else v for v in ser["St_inf"]], float)
    ok = np.isfinite(st)
    x_b = pe["xb_cm"] / 100.0
    r = np.where(xs < x_b, 0.845, 0.89)
    q = st * rhoucp * (Taw(r) - Tw)
    x_min = xs[int(np.nanargmin(np.where(xs < 0.30, st, np.nan)))]
    pre = [i for i in range(len(xs)) if ok[i] and xs[i] <= x_min + 1e-12][:5]
    C = min(q[i] * np.sqrt(xs[i]) for i in pre)
    xo, qo, so = xs[ok], q[ok], st[ok]
    den_H = mu_suth(T_inf) * CP * (Tt - Tw)
    ReH = (2 * C * np.sqrt(xo[0]) + cumtrapz(xo, qo)) / den_H
    # 比較は Re_v ≥ 1e6 の乱流点の範囲で (事前登録の実験 Re_H 範囲)
    sel = np.array([v in pe["pts_cm"] for v in np.round(xo * 100, 2)])
    return ReH[sel], so[sel], den_H, rhoucp * (Taw(0.89) - Tw), len(pre)


def evaluate(run, step, gate, cond):
    su = json.loads((Path(run) / "case_setup.json").read_text())
    sid = su["series"]
    grp = next(g for g, v in gate["groups"].items() if sid in v["series"])
    ser = next(s for s in cond["series"] if s["id"] == sid)
    xs = np.array(cond["x_cm"]) / 100.0
    ReH_e, St_e, den_H, den_St, npre = exp_curve(ser, su, gate["per_series_exp"][sid], xs)
    ux, q = wall_q(run, step)
    ReH_c = cumtrapz(ux, q) / den_H
    St_c = q / den_St
    pts = gate["groups"][grp]["ReH_points"]
    R = [np.interp(p, ReH_c, St_c) / np.interp(p, ReH_e, St_e) for p in pts]
    return dict(series=sid, group=grp, TwTt=su["Tw_over_Tt"], pts=pts, R=R, npre=npre,
                St_c=[float(np.interp(p, ReH_c, St_c)) for p in pts])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--series-csv-dir", default=None, help="各 run の St_CFD 時系列 CSV (check_quasisteady 用) の出力先")
    a = ap.parse_args()
    gate = next(g for g in json.loads((CASE / "acceptance.json").read_text())["gates"] if g["id"] == "G14-E")
    cond = json.loads((CASE / "conditions.json").read_text())
    rows = []
    for run in a.runs:
        steps = sorted(int(re.search(r"res_wall_4_(\d+)\.h5$", s).group(1)) for s in glob.glob(str(Path(run) / "res_wall_4_*.h5")))
        e = evaluate(run, steps[-1], gate, cond)
        rows.append(e)
        note = "" if e["npre"] > 1 else "  (前縁入熱の外挿感度は未評価: 遷移前 1 点)"
        print(f"[{Path(run).name}] step {steps[-1]}  {e['series']}  群 {e['group']}  Tw/Tt {e['TwTt']}{note}")
        print("   Re_H: " + "  ".join(f"{p:>6d}" for p in e["pts"]))
        print("   R   : " + "  ".join(f"{r:6.3f}" for r in e["R"]))
        if a.series_csv_dir:
            out = []
            for s in steps:
                out.append([s] + evaluate(run, s, gate, cond)["St_c"])
            f = Path(a.series_csv_dir) / f"_g14_{Path(run).name}.csv"
            np.savetxt(f, np.array(out), delimiter=",", header=",".join(["step"] + [f"St_ReH{p}" for p in e["pts"]]),
                       comments="", fmt="%.10g")
            print(f"   -> {f}")
    for g in ("H", "L"):
        rs = sorted([r for r in rows if r["group"] == g], key=lambda r: r["TwTt"])
        if len(rs) < 2:
            continue
        print(f"\n群 {g}: R(Tw/Tt) (行 = Tw/Tt、列 = Re_H)")
        for r in rs:
            print(f"   {r['TwTt']:4.2f}  " + "  ".join(f"{v:6.3f}" for v in r["R"]))
        lo, hi = rs[0], rs[-1]
        eff = [h / l for h, l in zip(hi["R"], lo["R"])]
        print(f"   効果量 R(Tw/Tt={hi['TwTt']}) / R(Tw/Tt={lo['TwTt']}) = " + "  ".join(f"{v:6.3f}" for v in eff))


if __name__ == "__main__":
    main()
