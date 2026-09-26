#!/usr/bin/env python3
"""Cary TN D-5863 と forge を**エネルギー厚さ Reynolds 数 Re_H** の座標で比べる
(case/59、acceptance.json の T4-0a-E。codex diagnose 2026-09-27 の提案)。

    python3 tools/cary_energy_compare.py <run_dir> [--series Re0.27_Tw0.2] [--series-csv out.csv]

    Re_H(x) = [∫₀ˣ q_w(s) ds] / [μ_e (h₀,e − h_w)],   h₀,e − h_w = cp (Tt − Tw)   (熱量的完全の空気)

- 実験側: q_w,exp = St_exp · ρ∞u∞cp (Taw(r) − Tw)。r は原報の還元規約 (層流 0.845 / 乱流 0.89)。
  遷移開始 x_b より上流を層流、それ以降を乱流として戻す (遷移部の規約は原報に明記が無いので、
  全区間 0.89 とした場合の差を別に表示する — 判定には使わない)。
- **未計測の前縁入熱 Q₀ = ∫₀^{7.30 cm} q_w dx** は、遷移前の 5 点 (7.30–13.02 cm) の q_w√x の
  最小係数 (腕 A) / 最大係数 (腕 B) で q_w = C/√x と外挿して 2C√x₁ とする。7.30 cm 以降は台形積分。
- CFD 側: 前縁 (x=0) から台形積分。**加熱ピークで積分をゼロに戻さない**。
- 比較: 実験点 x_exp の Re_H,exp と等しい Re_H,CFD の位置 x_CFD で St_CFD を取り、R = St_CFD / St_exp。
  比較点は acceptance.json T4-0a-0 の 17 点、重み 1/17。
"""
import argparse, glob, json, re
from pathlib import Path
import numpy as np
import h5py

CASE = Path(__file__).resolve().parents[1]
GAM, CP = 1.4, 1004.5
R_LAM, R_TURB = 0.845, 0.89


def mu_suth(T):
    return 1.716e-5 * (T / 273.0) ** 1.5 * (273.0 + 111.0) / (T + 111.0)


def wall_q(run, step):
    with h5py.File(Path(run) / f"res_wall_4_{step}.h5") as h:
        c = np.asarray(h["MESH/COORD"], float).reshape(-1, 3)
        q = -np.asarray(h["VALUE/qwall"], float)
    m = c[:, 0] >= 0.0
    x = np.round(c[m, 0], 9); ux = np.unique(x)
    return ux, np.array([q[m][x == u].mean() for u in ux])


def cumtrapz(x, y):
    return np.concatenate([[0.0], np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(x))])


def exp_ReH(x_m, q, x_b, C, den_H):
    """実験の Re_H。x_m [m] は測定点 (昇順)、q は q_w,exp。前縁 0–x_m[0] は C/√x。"""
    Q0 = 2.0 * C * np.sqrt(x_m[0])
    return (Q0 + cumtrapz(x_m, q)) / den_H


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--series", default=None, help="省略時は run の case_setup.json の系列")
    ap.add_argument("--xb-cm", type=float, default=None, help="遷移開始 x_b [cm] (省略時は acceptance.json)")
    ap.add_argument("--step", type=int, default=None)
    ap.add_argument("--series-csv", default=None)
    a = ap.parse_args()
    run = Path(a.run)
    su = json.loads((run / "case_setup.json").read_text())
    cond = json.loads((CASE / "conditions.json").read_text())
    acc = {g["id"]: g for g in json.loads((CASE / "acceptance.json").read_text())["gates"]}
    gE = acc["T4-0a-E"]
    sid = a.series or su["series"]
    ser = next(s for s in cond["series"] if s["id"] == sid)
    x_b = (a.xb_cm or gE["x_b_cm"][sid]) / 100.0
    pre = [p / 100.0 for p in gE["pre_transition_points_cm"]]
    pts = acc["T4-0a-0"]["compare_points_cm"]

    T_inf, U, ro, M, Tw, Tt = su["T_inf"], su["U_inf"], su["ro_inf"], su["M"], su["Tw"], su["Tt"]
    rhoucp = ro * U * CP
    Taw = lambda r: T_inf * (1 + r * 0.5 * (GAM - 1) * M ** 2)
    den_H = mu_suth(T_inf) * CP * (Tt - Tw)                  # μ_e (h0,e − h_w)、M_e = M_∞ (α=0)
    den_St = rhoucp * (Taw(R_TURB) - Tw)                     # 比較点 (乱流) の St の分母

    xs = np.array(cond["x_cm"]) / 100.0
    st = np.array([np.nan if v is None else v for v in ser["St_inf"]], float)
    ok = np.isfinite(st)
    xs_ok, st_ok = xs[ok], st[ok]

    def q_exp(r_rule):
        r = np.where(xs_ok < x_b, R_LAM, R_TURB) if r_rule == "split" else np.full_like(xs_ok, R_TURB)
        return st_ok * rhoucp * (Taw(r) - Tw)

    def coef(q):
        k = np.array([np.interp(p, xs_ok, q) * np.sqrt(p) for p in pre])
        return k.min(), k.max()

    q_split = q_exp("split")
    Cmin, Cmax = coef(q_split)
    arms = {"A": Cmin, "B": Cmax}
    ReH_exp = {k: exp_ReH(xs_ok, q_split, x_b, C, den_H) for k, C in arms.items()}
    # 参考: 遷移部も乱流 r で戻した場合 (判定に使わない)
    q_turb = q_exp("turb")
    ReH_exp_rT = exp_ReH(xs_ok, q_turb, x_b, coef(q_turb)[0], den_H)

    def evaluate(step):
        ux, q = wall_q(run, step)
        ReH_cfd = cumtrapz(ux, q) / den_H
        out = {}
        for k, R_e in list(ReH_exp.items()) + [("A_rT", ReH_exp_rT)]:
            rows = []
            for p in pts:
                xe = p / 100.0
                if not np.isfinite(np.interp(xe, xs, st)):
                    continue
                target = np.interp(xe, xs_ok, R_e)
                xc = np.interp(target, ReH_cfd, ux)                # Re_H,CFD は x に単調増加
                stc = np.interp(xc, ux, q) / den_St
                rows.append((p, target, xc * 100, stc, stc / np.interp(xe, xs_ok, st_ok)))
            out[k] = rows
        return out

    steps = sorted(int(re.search(r"res_wall_4_(\d+)\.h5$", s).group(1)) for s in glob.glob(str(run / "res_wall_4_*.h5")))
    step = a.step or steps[-1]
    res = evaluate(step)
    print(f"[{run.name}] step {step}  系列 {sid}  x_b {x_b*100:.2f} cm  μ_e(h0e−hw) {den_H:.4g}  "
          f"C_A {Cmin:.4g}  C_B {Cmax:.4g} (W/m^1.5、max/min {Cmax/Cmin:.3f})")
    print("  x_exp  Re_H,exp(A)  x_CFD(A) cm  St_CFD(A)   R_A     | Re_H,exp(B)  x_CFD(B)   R_B    R_B/R_A−1 | R (A, 遷移部も r0.89)")
    RA = np.array([r[4] for r in res["A"]]); RB = np.array([r[4] for r in res["B"]]); RT = np.array([r[4] for r in res["A_rT"]])
    for ra, rb, rt in zip(res["A"], res["B"], res["A_rT"]):
        print(f"  {ra[0]:5.2f}  {ra[1]:10.4g}  {ra[2]:8.2f}   {ra[3]:.3e}  {ra[4]:6.3f}  | {rb[1]:10.4g}  {rb[2]:8.2f}  {rb[4]:6.3f}  {rb[4]/ra[4]-1:+7.4f}  | {rt[4]:6.3f}")
    d = RB / RA - 1
    print(f"  R_A 平均 {RA.mean():.3f} [{RA.min():.3f}, {RA.max():.3f}]  R_B 平均 {RB.mean():.3f}  "
          f"|R_B/R_A−1| 最大 {np.abs(d).max():.4f}  (参考: 遷移部 r 規約の差 最大 {np.abs(RT/RA-1).max():.4f})")
    v = ("全 17 点で |R_B/R_A−1| ≤ 0.05 → 『前縁欠測が比較を支配する』を退ける"
         if np.all(np.abs(d) <= 0.05) else "1 点以上で 0.05 超 → この比較方法の本採用を保留")
    print("  T4-0a-E の読み (事前登録どおり): " + v)
    if a.series_csv:
        rows = []
        for s in steps:
            e = evaluate(s)
            rows.append([s] + [r[3] for r in e["A"]] + [r[3] for r in e["B"]])
        cols = ["step"] + [f"stA_{r[0]:g}" for r in res["A"]] + [f"stB_{r[0]:g}" for r in res["B"]]
        np.savetxt(a.series_csv, np.array(rows), delimiter=",", header=",".join(cols), comments="", fmt="%.10g")
        print(f"  -> {a.series_csv} ({len(rows)} 枚)")


if __name__ == "__main__":
    main()
