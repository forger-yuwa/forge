#!/usr/bin/env python3
"""D-7275 試験 26 の中心線 St*_l (Fig 20 の digitize) と forge を比べる (case/60)。

    python3 tools/d7275_compare.py <run_dir> [--series-csv out.csv]

還元は実験と同じにする: St*_l = q_w / [(ρVcp)*_l (Taw − Tw)]、(ρVcp)*_l = 0.413 × 68.89 kW/m²K、Taw = 1728 K
(Table III の試験 26 行)。CFD の Tw は run の値 (300 K)。
x の対応: x = R*_l / (0.180 × 4.757e6 /m)。この換算は最下流点をパネル端より 3 % 下流に置く (±3 % の不確かさ)。
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
    a = ap.parse_args()
    run = Path(a.run)
    su = json.loads((run / "case_setup.json").read_text())
    d = json.loads((CASE / "digitize_d7275_fig20_test26.json").read_text())
    t3 = d["test26_table"]["TableIII"]; t2 = d["test26_table"]["TableII"]
    rvc = t3["rhoVcp_star_l_over_inf"] * t2["rhoVcp_inf_kW_m2K"] * 1e3
    den = rvc * (t3["Taw_K"] - su["Tw"])
    xfac = t3["RstarLl_over_RLinf"] * t2["Re_inf_per_m"]
    pts = [(p["R_star_l"] / xfac, p["St_star_l"]) for p in d["points"]]
    steps = sorted(int(re.search(r"res_plate_4_(\d+)\.h5$", s).group(1)) for s in glob.glob(str(run / "res_plate_4_*.h5")))
    ux, q = wall_q(run, steps[-1])
    St = q / den
    print(f"[{run.name}] step {steps[-1]}  (ρVcp)*_l {rvc:.0f} W/m²K  Taw {t3['Taw_K']} K  Tw {su['Tw']} K")
    print("   x [m]   St*_exp    St*_CFD(A 前縁)  R_A     St*_CFD(B トリップ)  R_B")
    RA, RB = [], []
    for x, se in pts:
        sa = np.interp(x, ux, St); sb = np.interp(x - 0.127, ux, St)
        RA.append(sa / se); RB.append(sb / se)
        print(f"  {x:6.3f}  {se:.3e}   {sa:.3e}      {sa/se:6.3f}   {sb:.3e}          {sb/se:6.3f}")
    RA, RB = np.array(RA), np.array(RB)
    print(f"  R_A 平均 {RA.mean():.3f} [{RA.min():.3f}, {RA.max():.3f}]   R_B 平均 {RB.mean():.3f} [{RB.min():.3f}, {RB.max():.3f}]")
    # x 対応 ±3 % の感度 (A)
    for f in (0.97, 1.03):
        r = np.array([np.interp(x * f, ux, St) / se for x, se in pts])
        print(f"  x × {f}: R_A 平均 {r.mean():.3f}")
    # パネル平均 (x 1.07–2.46 m) と carpet plot 比較用の q
    m = (ux >= 1.07) & (ux <= 2.46)
    print(f"  パネル平均 q_w (1.07–2.46 m) = {np.trapz(q[m], ux[m])/(ux[m][-1]-ux[m][0])/1e3:.2f} kW/m²、位置 II (1.88 m) = {np.interp(1.88, ux, q)/1e3:.2f} kW/m²")
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
