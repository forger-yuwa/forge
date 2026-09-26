#!/usr/bin/env python3
"""Cary TN D-5863 の St と forge の St を比べる (case/59、acceptance.json の T4-0a-0 / G14)。

    python3 tools/cary_compare.py <run_dir> [--step N] [--series-csv out.csv]

St = q_w / (ρ∞ u∞ cp (Taw − Tw))、Taw = T∞(1 + 0.89·(γ−1)/2·M∞²) で Cary と同じに還元する。
比較原点: A = x_CFD = x_exp、B = x_CFD = x_exp − x_peak (acceptance.json の値)。
--series-csv を付けると、全スナップショットについて比較点の St_forge (原点 A/B) を書く (check_quasisteady 用)。
"""
import argparse, glob, json, re
from pathlib import Path
import numpy as np
import h5py

CASE = Path(__file__).resolve().parents[1]
GAM, CP = 1.4, 1004.5


def wall_profile(run, step):
    with h5py.File(Path(run) / f"res_wall_4_{step}.h5") as h:
        c = np.asarray(h["MESH/COORD"], float).reshape(-1, 3)
        q = -np.asarray(h["VALUE/qwall"], float)
    m = c[:, 0] >= 0.0
    x = np.round(c[m, 0], 9); ux = np.unique(x)
    return ux, np.array([q[m][x == u].mean() for u in ux])


def steps(run):
    return sorted(int(re.search(r"res_wall_4_(\d+)\.h5$", p).group(1))
                  for p in glob.glob(str(Path(run) / "res_wall_4_*.h5")))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--step", type=int, default=None)
    ap.add_argument("--series-csv", default=None)
    a = ap.parse_args()
    run = Path(a.run)
    su = json.loads((run / "case_setup.json").read_text())
    cond = json.loads((CASE / "conditions.json").read_text())
    acc = next(g for g in json.loads((CASE / "acceptance.json").read_text())["gates"] if g["id"] == "T4-0a-0")
    ser = next(s for s in cond["series"] if s["id"] == su["series"])
    xs = np.array(cond["x_cm"]); st = {x: v for x, v in zip(cond["x_cm"], ser["St_inf"])}
    pts = [p for p in acc["compare_points_cm"] if st.get(p)]
    xpk = 22.22
    Taw = su["T_inf"] * (1 + 0.89 * 0.5 * (GAM - 1) * su["M"] ** 2)
    den = su["ro_inf"] * su["U_inf"] * CP * (Taw - su["Tw"])
    sts = steps(run)
    step = a.step or sts[-1]

    def St_at(x_cm, ux, q):
        return np.interp(x_cm / 100.0, ux, q) / den

    ux, q = wall_profile(run, step)
    print(f"[{run.name}] step {step}  series {su['series']}  Taw(r0.89) {Taw:.1f} K  Tw {su['Tw']:.1f} K  "
          f"ρ∞u∞cp(Taw−Tw) {den:.4g} W/m²K")
    print("  x_exp[cm]  St_exp     St_A(x)    St_B(x−22.22)  R_A     R_B     D=R_B/R_A−1")
    RA, RB = [], []
    for p in pts:
        sa, sb = St_at(p, ux, q), St_at(p - xpk, ux, q)
        ra, rb = sa / st[p], sb / st[p]; RA.append(ra); RB.append(rb)
        print(f"  {p:8.2f}  {st[p]:.3e}  {sa:.3e}  {sb:.3e}     {ra:6.3f}  {rb:6.3f}  {rb/ra-1:+7.3f}")
    RA, RB = np.array(RA), np.array(RB); D = RB / RA - 1
    print(f"  平均 R_A {RA.mean():.3f} [{RA.min():.3f}, {RA.max():.3f}]  R_B {RB.mean():.3f} [{RB.min():.3f}, {RB.max():.3f}]  "
          f"D {D.mean():+.3f} [{D.min():+.3f}, {D.max():+.3f}]")
    if np.all(np.abs(D) > 0.10):
        v = "全点 |D| > 10 % → 『原点感度は小さい』を棄却 (比較座標を先に固定する)"
    elif np.all(np.abs(D) < 0.10):
        v = "全点 |D| < 10 % → 原点の影響は小さい"
    else:
        v = "点によって 10 % をまたぐ → 判定不能"
    print("  T4-0a-0 の読み (acceptance.json の事前登録どおり): " + v)
    # 参考: CFD の加熱の立ち上がり (前縁からの St の形)
    Sx = q / den
    print(f"  参考: CFD St は x=1 cm で {np.interp(0.01, ux, Sx):.2e}、5 cm で {np.interp(0.05, ux, Sx):.2e}、"
          f"最大は x={ux[np.argmax(Sx)]*100:.2f} cm")
    if a.series_csv:
        rows = []
        for s in sts:
            u2, q2 = wall_profile(run, s)
            rows.append([s] + [St_at(p, u2, q2) for p in pts] + [St_at(p - xpk, u2, q2) for p in pts])
        cols = ["step"] + [f"stA_{p:g}" for p in pts] + [f"stB_{p:g}" for p in pts]
        np.savetxt(a.series_csv, np.array(rows), delimiter=",", header=",".join(cols), comments="", fmt="%.8g")
        print(f"  -> {a.series_csv} ({len(rows)} 枚)  check_quasisteady --series-cols {','.join(cols[1:])}")


if __name__ == "__main__":
    main()
