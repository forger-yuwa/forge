#!/usr/bin/env python3
"""D-7275 の中心線 St*_l (試験 26 = Fig 20、試験 28 = Fig 19(c) の digitize) と forge を比べる (case/60)。

    python3 tools/d7275_compare.py <run_dir> [--series-csv out.csv] [--x-map tc|R]

試験は run の case_setup.json の d7275_test で選ぶ (無ければ 26)。
還元は実験と同じにする: St*_l = q_w / [(ρVcp)*_l (Taw − Tw)]、(ρVcp)*_l と Taw は Table III の該当行。CFD の Tw は run の値 (300 K)。
x の対応 (2026-09-28 から既定 tc): Table I の中心線熱電対位置 − 0.0096 m (鋭い前縁基準)。R*/(係数) で出す x (--x-map R) は
試験 26 で熱電対位置の 1.03 倍、試験 28 で 0.98 倍と一様にずれる (係数側の丸め・読み取りの系統差)。
乱流の起点は 2 通り: A = 前縁 (x_CFD = x)、B = トリップ位置 12.7 cm (x_CFD = x − 0.127 m)。有効長の感度であって上下界ではない。
"""
import argparse, glob, json, re
from pathlib import Path
import numpy as np
import h5py

CASE = Path(__file__).resolve().parents[1]


def wall_q(run, step):
    with h5py.File(Path(run) / f"res_plate_4_{step}.h5") as h:
        c = np.asarray(h["MESH/COORD"], float).reshape(-1, 3)
        q = -np.asarray(h["VALUE/qwall"], float)
    m = (c[:, 0] >= 0.0) & (np.abs(c[:, 1]) < 1e-9)
    x = np.round(c[m, 0], 9); ux = np.unique(x)
    return ux, np.array([q[m][x == u].mean() for u in ux])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--series-csv", default=None)
    ap.add_argument("--x-map", choices=("tc", "R"), default="tc")
    a = ap.parse_args()
    run = Path(a.run)
    su = json.loads((run / "case_setup.json").read_text())
    test = int(su.get("d7275_test", 26))
    fn = {26: "digitize_d7275_fig20_test26.json", 28: "digitize_d7275_fig19c_test28.json"}[test]
    d = json.loads((CASE / fn).read_text())
    t3 = d[f"test{test}_table"]["TableIII"]; t2 = d[f"test{test}_table"]["TableII"]
    rvc = t3["rhoVcp_star_l_over_inf"] * t2["rhoVcp_inf_kW_m2K"] * 1e3
    den = rvc * (t3["Taw_K"] - su["Tw"])
    xfac = t3["RstarLl_over_RLinf"] * t2["Re_inf_per_m"]
    if a.x_map == "tc":
        pts = [(p["x_sharp_m"], p["St_star_l"]) for p in d["points"]]
    else:
        pts = [(p["R_star_l"] / xfac, p["St_star_l"]) for p in d["points"]]
    steps = sorted(int(re.search(r"res_plate_4_(\d+)\.h5$", s).group(1)) for s in glob.glob(str(run / "res_plate_4_*.h5")))
    ux, q = wall_q(run, steps[-1])
    St = q / den
    print(f"[{run.name}] 試験 {test}  step {steps[-1]}  x 対応 {a.x_map}  (ρVcp)*_l {rvc:.0f} W/m²K  Taw {t3['Taw_K']} K  Tw {su['Tw']} K")
    print("   x [m]   St*_exp    St*_CFD(A 前縁)  R_A     St*_CFD(B トリップ)  R_B")
    RA, RB = [], []
    def at(xx):   # 範囲外は NaN (np.interp は端点値を黙って返す、codex 2026-09-27)
        return float(np.interp(xx, ux, St)) if ux.min() <= xx <= ux.max() else float("nan")
    for x, se in pts:
        sa = at(x); sb = at(x - 0.127)
        RA.append(sa / se); RB.append(sb / se)
        print(f"  {x:6.3f}  {se:.3e}   {sa:.3e}      {sa/se:6.3f}   {sb:.3e}          {sb/se:6.3f}")
    RA, RB = np.array(RA), np.array(RB)
    print(f"  R_A 平均 {np.nanmean(RA):.3f} [{np.nanmin(RA):.3f}, {np.nanmax(RA):.3f}]   R_B 平均 {np.nanmean(RB):.3f} [{np.nanmin(RB):.3f}, {np.nanmax(RB):.3f}]")
    # x 対応 ±3 % の感度 (A)
    for f in (0.97, 1.03):
        r = np.array([at(x * f) / se for x, se in pts])
        print(f"  x × {f}: R_A 平均 {np.nanmean(r):.3f} (範囲外 {int(np.isnan(r).sum())} 点)")
    # パネル平均 (x 1.07–2.46 m) と carpet plot 比較用の q
    m = (ux >= 1.07) & (ux <= 2.46)
    integ = getattr(np, "trapezoid", None) or np.trapz
    q_avg = integ(q[m], ux[m]) / (ux[m][-1] - ux[m][0])
    print(f"  パネル平均 q_w (1.07–2.46 m) = {q_avg/1e3:.2f} kW/m²、位置 II (1.88 m) = {np.interp(1.88, ux, q)/1e3:.2f} kW/m²")
    if a.series_csv:
        rows = []
        for s in steps:
            u2, q2 = wall_q(run, s)
            rows.append([s] + [np.interp(x, u2, q2) / den for x, _ in pts])
        np.savetxt(a.series_csv, np.array(rows), delimiter=",",
                   header=",".join(["step"] + [f"St_x{x:.3f}" for x, _ in pts]), comments="", fmt="%.10g")
        print(f"  -> {a.series_csv}")


if __name__ == "__main__":
    main()
