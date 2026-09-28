#!/usr/bin/env python3
r"""case/63 Graetz の運動量収支の後処理 A/B (plan `boundary-cht-axisymmetric-graetz.md` §5.1 #6c、事前登録 2026-09-28)。

V-g1 ① (−dp/dx と 8μŪ/R² の差 2.25 %) の由来を、既存の保存場だけで調べる (追加計算 0 step)。
区間 [x1, x2] = [0.25, 0.75] L_heat (最も近い節点列)、全項 per rad:

    A: 断面積分した圧力力 ΔF_p = [∫p r dr]_{x2}^{x1} と、Poiseuille 抵抗の区間積分 4μ∫Ū dx の差
    B: 定常 x 運動量の積分形
         [∫(ρu² + p − τ_xx) r dr]_{x1}^{x2} − R ∫ τ_rx,w dx = 0
       の不一致 (τ_xx = 2μ u_x − (2/3)μ ∇·u、∇·u = u_x + v_r + v/r、τ_rx,w = μ ∂u/∂r|_{r=R})

事前基準: B の不一致 ≤ 0.1 % of ΔF_p、かつ抽出の不確かさ (壁せん断の片側差分を 2 次 → 3 次に替えた B の差。**2026-09-28 結果前に修正**: 初版の 2 次 → 1 次は厳密な Poiseuille の合成場でも 3.1 % を出し、1 次差分自体の誤差を測るだけだった) ≤ その 1/3。
判定: A で差が残り B が閉じて追加項が差を説明 → 第 1 仮説 (参照式の欠落項) を支持 / B が閉じない → 棄却 /
抽出精度・定常性が不足 → 判定不能。**B が閉じても Nu 精度や V-g1 ② は合格にならない。**

    python3 case/63.graetz_cht/momentum_balance.py <run> [--last 25]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

import graetz_common as gc
from eval_graetz import Run, rint


def ddx(Rn: Run, f, key, xv, order):
    """節点列 x = xv での ∂(key)/∂x (列ごと、非等間隔)。order 2: 中心 3 点、1: 前進 2 点。"""
    k = int(np.searchsorted(Rn.xcols, xv))
    xm, x0, xp = Rn.xcols[k - 1], Rn.xcols[k], Rn.xcols[k + 1]
    fm, f0, fp = (f[key][Rn.col[x]] for x in (xm, x0, xp))
    if order == 3:
        order = 2           # ∂u/∂x は項が小さい (十分発達流) ので 3 点中心のまま。比較は壁せん断の次数だけ
    a, b = x0 - xm, xp - x0
    return (-b / (a * (a + b))) * fm + ((b - a) / (a * b)) * f0 + (a / (b * (a + b))) * fp


def section_terms(Rn: Run, f, xv, order):
    i = Rn.col[np.round(xv, 10)]
    r = Rn.xyz[i, 1]
    ro, u, v, p = f["ro"][i], f["Ux"][i], f["Uy"][i], f["P"][i]
    ux = ddx(Rn, f, "Ux", xv, order)
    vr = np.gradient(v, r, edge_order=2)
    v_over_r = np.where(r > 0, v / np.where(r > 0, r, 1.0), vr)      # 軸上は極限 ∂v/∂r
    div = ux + vr + v_over_r
    txx = 2 * gc.MU * ux - (2.0 / 3.0) * gc.MU * div
    return dict(mom=rint(ro * u * u * r, r), p=rint(p * r, r), txx=rint(txx * r, r),
                ubar=rint(u * r, r) / (gc.R ** 2 / 2))


def wall_shear(Rn: Run, f, xv, order):
    i = Rn.col[np.round(xv, 10)]
    r = Rn.xyz[i, 1]
    u = f["Ux"][i]
    h = r[-1] - r[-2]
    if order == 2:
        dudr = (3 * u[-1] - 4 * u[-2] + u[-3]) / (2 * h)
    else:   # order 3: 4 点片側 (3 次式まで厳密)
        dudr = (11 * u[-1] - 18 * u[-2] + 9 * u[-3] - 2 * u[-4]) / (6 * h)
    return gc.MU * dudr


def balance(Rn: Run, st, order):
    f = Rn.field(st)
    x1 = Rn.xcols[np.argmin(np.abs(Rn.xcols - 0.25 * gc.L_HEAT))]
    x2 = Rn.xcols[np.argmin(np.abs(Rn.xcols - 0.75 * gc.L_HEAT))]
    s1, s2 = section_terms(Rn, f, x1, order), section_terms(Rn, f, x2, order)
    xs = Rn.xcols[(Rn.xcols >= x1) & (Rn.xcols <= x2)]
    tw = np.array([wall_shear(Rn, f, x, order) for x in xs])
    ub = np.array([section_terms(Rn, f, x, order)["ubar"] for x in xs]) if order == 2 else None
    trap = getattr(np, "trapezoid", None) or np.trapz
    Fw = gc.R * trap(tw, xs)                                    # 壁から受ける x 方向の力 (負)
    dP = s1["p"] - s2["p"]                                       # 圧力力 (正)
    dM = s2["mom"] - s1["mom"]                                   # 運動量流束の増分
    dT = s2["txx"] - s1["txx"]                                   # 軸方向粘性応力の増分
    resid = (s2["mom"] + s2["p"] - s2["txx"]) - (s1["mom"] + s1["p"] - s1["txx"]) - Fw
    out = dict(step=st, x1=x1, x2=x2, dP=dP, dM=dM, dTxx=dT, Fw=Fw, resid=resid)
    if ub is not None:
        out["A_pois"] = 4 * gc.MU * trap(ub, xs)                # 4μ ∫Ū dx (= ∫8μŪ/R² dx · R²/2)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run")
    ap.add_argument("--last", type=int, default=25)
    a = ap.parse_args()
    Rn = Run(a.run)
    steps = sorted(int(p.name[4:-3]) for p in Path(a.run).glob("res_[0-9]*.h5"))[-a.last:]
    rows = []
    for st in steps:
        b2, b1 = balance(Rn, st, 2), balance(Rn, st, 3)
        rows.append([st, b2["dP"], b2["A_pois"], b2["dM"], b2["dTxx"], b2["Fw"], b2["resid"], b1["resid"]])
    R = np.array(rows)
    cols = "step,dP,A_pois,dM,dTxx,Fw,resid_o2,resid_o1"
    out = Path(a.run) / "momentum_balance_series.csv"
    np.savetxt(out, R, delimiter=",", header=cols, comments="", fmt=["%d"] + ["%.12e"] * 7)
    last = R[-1]
    dP, A, dM, dT, Fw, r2, r1 = last[1:]
    b_rel = abs(r2) / dP
    ext = abs(r2 - r1) / dP
    print(f"=== 運動量収支 {a.run}  (step {int(last[0])}、区間 x {balance(Rn, int(last[0]), 2)['x1']*1e3:.2f}–"
          f"{balance(Rn, int(last[0]), 2)['x2']*1e3:.2f} mm、per rad)")
    print(f"  圧力力 ΔF_p           {dP:.6e}")
    print(f"  A: Poiseuille 4μ∫Ū dx {A:.6e}   → (ΔF_p − A)/ΔF_p = {(dP - A)/dP:+.4e}")
    print(f"  B の項: 運動量流束の増分 {dM:+.4e} ({dM/dP:+.3e} of ΔF_p)、軸方向粘性応力の増分 {dT:+.4e} ({dT/dP:+.3e})、"
          f"壁せん断の力 {Fw:+.4e} ({Fw/dP:+.3e})")
    print(f"  B の不一致 (壁せん断 2 次) {r2:+.4e} = {r2/dP:+.4e} of ΔF_p   (3 次 {r1/dP:+.4e})")
    print(f"  抽出の不確かさ |B(2 次) − B(3 次)|/ΔF_p = {ext:.3e}")
    # 系列の変動 (末尾 25 枚の max−min を ΔF_p で割る)
    for j, nm in ((1, "ΔF_p"), (2, "A"), (3, "dM"), (5, "Fw"), (6, "B")):
        print(f"  参考: 末尾 {len(R)} 枚の {nm} の変動 (max−min)/ΔF_p = {np.ptp(R[:, j])/dP:.3e}")
    ok_b = b_rel <= 1e-3
    ok_e = ext <= 1e-3 / 3
    print(f"\n  {'PASS' if ok_b else 'FAIL'}  B の不一致 ≤ 0.1 % of ΔF_p   ({b_rel:.3e})")
    print(f"  {'PASS' if ok_e else 'FAIL'}  抽出の不確かさ ≤ 1/3 × 0.1 %  ({ext:.3e})")
    verdict = "PASS" if (ok_b and ok_e) else ("判定不能" if not ok_e else "FAIL")
    print(f"VERDICT (B の閉じ): {verdict}   系列 → {out}  (準定常は check_quasisteady.py --series-csv で別に)")
    return 0 if verdict == "PASS" else (2 if verdict == "判定不能" else 1)


if __name__ == "__main__":
    sys.exit(main())
