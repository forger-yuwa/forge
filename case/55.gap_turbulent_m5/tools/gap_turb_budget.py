#!/usr/bin/env python3
"""すきま内で乱流モデルがどれだけ仕事をしているかを、収束場から直接測る。

深部の熱流束が文献より 1-2 桁低い件で、格子・非定常・前処理は既に潰れている。
残る候補のうち「近淀み域の乱流モデル」は**新たに run を回さなくても測れる**:
壁熱流束の熱伝導率は $\\lambda + \\lambda_t$ ($\\lambda_t = c_p\\mu_t/Pr_t$) なので、
深さごとの $\\mu_t/\\mu$ から**乱流が担う割合**が出る。さらに、文献帯に届かせるには
$\\mu_t/\\mu$ がいくつ必要かを逆算できる。

出力: 深さ z/W ごとに
  q_w/q_fp          (後壁・前壁)
  mu_t/mu           壁隣接ノードの値
  lambda_t/lambda   = (mu_t/mu)(Pr/Pr_t)  ← 乱流が熱流束に効く比率
  必要 mu_t/mu      文献帯 (既定 0.10 q_fp @ z/W<=6) に届くのに要する値
"""
import argparse, json, sys
from pathlib import Path
import numpy as np
import h5py

HERE = Path(__file__).resolve().parent
CASE = HERE.parent
PR_T = 0.9


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--res", default=None)
    ap.add_argument("--qfp-from", default="run_0005_T0_iso500_main")
    ap.add_argument("--band", type=float, default=0.10,
                    help="文献帯の目標 q/q_fp (既定 0.10 = TP-1187 の『乱流は深さ 60 % で 0.1 q_FP』)")
    a = ap.parse_args()

    rd = CASE / a.run
    if sorted(rd.glob("res_nan_*.h5")):
        sys.exit("REFUSED: res_nan_* がある (発散 run)")
    res = Path(a.res) if a.res else sorted(
        rd.glob("res_[0-9]*.h5"), key=lambda f: int(f.stem.split("_")[1]))[-1]

    geom = json.loads((CASE / "geometry.json").read_text(encoding="utf-8"))
    setup = json.loads((rd / "case_setup.json").read_text(encoding="utf-8"))
    w = geom["cavity"]["widths"][str(setup["series"]["w_over_d"])] * 1e-3
    d = geom["cavity"]["depth"] * 1e-3
    xr = geom["cavity"]["x_rear_wall_from_le"] * 1e-3
    xf = xr - w

    arr = np.loadtxt(CASE / a.qfp_from / "wall_q.csv", delimiter=",", skiprows=1)
    q_fp = float(np.interp(xr - 0.5 * w, arr[:, 0], arr[:, 1]))

    with h5py.File(res, "r") as f:
        c = np.asarray(f["MESH"]["COORD"]).reshape(-1, 3)
        V = f["VALUE"]
        mut = np.asarray(V["vis_turb"]); mul = np.asarray(V["vis_lam"])
        lam = np.asarray(V["thermCond"]); T = np.asarray(V["T"])
        cp_over = None
    x, y = c[:, 0], c[:, 1]
    # Pr = cp mu / lambda を場から復元 (cp は出力していないので lambda, mu, Pr の関係で逆算しない。
    # かわりに lambda_t/lambda = (mu_t/mu)(Pr/Pr_t) の Pr を層流値から作る)
    # Pr_lam = cp*mu/lam。cp は T から出せないので、代表値を case_setup から取る。
    cp = 1005.0
    Pr_lam = cp * mul / np.maximum(lam, 1e-30)

    # すきま内の x 列を並べ、**第一内部列**を取る。
    # 壁ノードそのものは SST が k をピン留めするので mu_t = 0 で、乱流の寄与を測れない。
    tol = 1e-9
    ingap = (y < tol) & (y > -d - tol) & (x >= xf - tol) & (x <= xr + tol)
    xs = np.unique(np.round(x[ingap], 9))
    cols = {"後壁 第一内部": xs[-2], "前壁 第一内部": xs[1], "すきま中央": xs[len(xs) // 2]}

    def column(xw):
        sel = (np.abs(x - xw) < tol) & (y < tol) & (y > -d - tol)
        idx = np.where(sel)[0]
        return idx[np.argsort(-y[idx])]

    out = {}
    for name, xw in cols.items():
        idx = column(xw)
        if len(idx) == 0:
            continue
        zw = -y[idx] / w
        out[name] = dict(z_over_w=zw, mut_mu=mut[idx] / np.maximum(mul[idx], 1e-30),
                         Pr=Pr_lam[idx], T=T[idx])

    print(f"run {a.run}  res {res.name}   q_fp = {q_fp*1e-3:.2f} kW/m²  (W = {w*1e3:.2f} mm)")
    names = list(cols.keys())
    print("第一内部列の値 (壁ノードは k ピン留めで µt=0 のため)")
    print(f"{'z/W':>6} |" + "".join(f" {n+' µt/µ':>16} {'λt/λ':>9} |" for n in names))
    targets = [0.5, 1, 2, 3, 4, 6, 8, 12, 20]
    for zt in targets:
        row = f"{zt:6.1f} |"
        for name in names:
            o = out.get(name)
            if o is None:
                row += f" {'-':>16} {'-':>9} |"; continue
            i = int(np.argmin(np.abs(o["z_over_w"] - zt)))
            r = o["mut_mu"][i]; lt = r * o["Pr"][i] / PR_T
            row += f" {r:16.3e} {lt:9.3e} |"
        print(row)

    # --- 深さごとの q/q_fp と、文献帯に届くのに必要な µt/µ ---
    prof = np.loadtxt(rd / "cavity_rear.csv", delimiter=",", skiprows=1)
    zq = prof[:, 0] * d / w            # x_over_d -> z/W
    qq = prof[:, 1] / q_fp
    o = out[names[0]]
    print(f"\n必要 µt/µ の逆算 (目標 q/q_fp = {a.band:.2f}; λ_eff/λ を倍率ぶん上げる近似)")
    print(f"{'z/W':>6} {'q/q_fp (後壁)':>14} {'不足倍率':>9} {'必要 µt/µ':>11} {'現状 µt/µ':>11} {'桁差':>6}")
    for zt in (1, 2, 3, 4, 6):
        j = int(np.argmin(np.abs(zq - zt)))
        i = int(np.argmin(np.abs(o["z_over_w"] - zt)))
        q = qq[j]
        if q <= 0:
            continue
        need_ratio = a.band / q                       # λ_eff/λ をこの倍率に
        need_mut = (need_ratio - 1.0) * PR_T / o["Pr"][i]
        cur = o["mut_mu"][i]
        dec = np.log10(need_mut / cur) if cur > 0 else float("inf")
        print(f"{zt:6.1f} {q:14.3e} {need_ratio:9.1f} {need_mut:11.1f} {cur:11.3e} "
              f"{dec:6.1f}" if np.isfinite(dec) else
              f"{zt:6.1f} {q:14.3e} {need_ratio:9.1f} {need_mut:11.1f} {cur:11.3e}   inf")
    print("  ※ 温度勾配が変わらないとした粗い近似。桁の議論にのみ使う。")

    # 外部 BL の代表値
    sel = (y > 0) & (y < 0.006) & (np.abs(x - (xr - 0.05)) < 2e-4)
    if sel.sum():
        print(f"\n参考: 開口 50 mm 上流の外部 BL 内 µt/µ 最大 = "
              f"{np.max(mut[sel]/np.maximum(mul[sel],1e-30)):.1f}")
    print(f"\nPr_lam (代表) = {np.median(out[names[0]]['Pr']):.3f},  Pr_t = {PR_T}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
