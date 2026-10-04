#!/usr/bin/env python3
"""T2' の入口 BL プロファイルを `inletProfile` CSV にする。

**`inlet_bl.json` は使わない**。あれは訂正前の台帳 (θ = 1.42 cm の誤読値) で δ* と θ を
同時に拘束した結果で、δ = 54.3 cm・残差 0.554 (未収束)・Π = −0.54 と、README が
「不可能」と判定した側の解になっている。訂正後の正本は `bl_derived.json` で、
そこでは **δ* を錨にして θ は出力**として扱う (README「訂正後の台帳」)。

ここでは同じプロファイル族 (van Driest + Spalding 内層 + Coles wake、Walz 温度) を使い、
  - δ* = 12.80 cm (Fig 4 の実測。TM X-71945 本文の独立記述と 5 % 一致)
  - u_tau = 66.31 m/s (`bl_derived.json`)
を与えて δ を解く。出力 θ が台帳の 0.511 cm に戻ることを**自己整合の確認**として印字する。

usage: python3 tools/gen_inlet_csv.py --out run_XXXX/inlet_profile_1.csv
"""
import argparse, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import importlib.util as ilu
_sp = ilu.spec_from_file_location("inlet_bl", HERE / "inlet_bl.py")
_ib = ilu.module_from_spec(_sp); _sp.loader.exec_module(_ib)
from gas_air import dry_air                                        # noqa: E402
_sc = ilu.spec_from_file_location("t2_conditions", HERE / "conditions.py")
_cm = ilu.module_from_spec(_sc); _sc.loader.exec_module(_cm)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--re", type=float, default=1.47e6, help="Re'/m の系列")
    ap.add_argument("--Tw", type=float, default=300.0)
    ap.add_argument("--H", type=float, default=0.3937, help="領域高さ [m]")
    ap.add_argument("--ny", type=int, default=600)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    root = HERE.parent
    cond = json.loads((root / "conditions.json").read_text(encoding="utf-8"))
    der = json.loads((root / "derived.json").read_text(encoding="utf-8"))
    bld = json.loads((root / "bl_derived.json").read_text(encoding="utf-8"))
    s = next(x for x in der["series"] if abs(x["Re_m"] - a.re) / a.re < 1e-6)
    b = next(x for x in bld["rows"] if abs(x["Re_m"] - a.re) / a.re < 1e-6)

    air = _cm.Air(dry_air(1000.0))
    fs = dict(ro=s["ro_inf"], U=s["U_inf"], mu=s["mu_inf"], T=s["T_inf"], p=s["p_inf"])
    prof = _ib.Profile(fs, a.Tw, s["T_aw"], air, cond["M"])
    u_tau = b["u_tau"]
    delta = prof.solve_one(b["dstar_cm"] * 1e-2, "dstar", u_tau)
    pr = prof.build(delta, u_tau, ny=a.ny)
    ds, th = prof.thicknesses(pr)
    k, om = prof.turbulence(pr)

    print(f"=== T2' 入口 BL (Re' = {a.re:.2e}/m, T_w = {a.Tw:.0f} K) ===")
    print(f"  入力  δ* = {b['dstar_cm']:.2f} cm (錨),  u_tau = {u_tau:.2f} m/s (bl_derived)")
    print(f"  解    δ  = {delta*100:.3f} cm   [台帳 {b['delta_cm']:.3f} cm, "
          f"差 {abs(delta*100-b['delta_cm'])/b['delta_cm']*100:.2f} %]")
    print(f"  出力  δ* = {ds*100:.3f} cm  θ = {th*100:.4f} cm   "
          f"[台帳 θ {b['theta_cm_derived']:.4f} cm, 差 {abs(th*100-b['theta_cm_derived'])/b['theta_cm_derived']*100:.2f} %]")
    print(f"  H = {ds/th:.2f} [台帳 {b['H'] if 'H' in b else float('nan'):.2f}],  "
          f"Re_theta = {fs['ro']*fs['U']*th/fs['mu']:.3e} [台帳 {b['Re_theta']:.3e}]")

    # 領域高さまで自由流で埋める (δ より上は一様)
    y = pr["y"]
    keep = y <= a.H
    y = y[keep]; k = k[keep]; om = om[keep]
    ro = pr["ro"][keep]; ux = pr["u"][keep]
    if y[-1] < a.H:
        y = np.append(y, a.H); k = np.append(k, k[-1]); om = np.append(om, om[-1])
        ro = np.append(ro, fs["ro"]); ux = np.append(ux, fs["U"])
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as f:
        f.write("y k omega ro Ux Uy Uz\n")
        for i in range(len(y)):
            f.write(f"{y[i]:.9e} {k[i]:.9e} {om[i]:.9e} {ro[i]:.9e} {ux[i]:.9e} 0 0\n")
    print(f"  -> {out}  ({len(y)} 点, y {y[0]:.3e}..{y[-1]:.4f} m)")
    print(f"     ro {ro.min():.4e}..{ro.max():.4e},  Ux {ux.min():.1f}..{ux.max():.1f},  "
          f"k {k.min():.3e}..{k.max():.3e},  omega {om.min():.3e}..{om.max():.3e}")


if __name__ == "__main__":
    main()
