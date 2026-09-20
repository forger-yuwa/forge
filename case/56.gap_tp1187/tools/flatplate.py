#!/usr/bin/env python3
"""T4 Gate B — 較正パネルの $q_{FP}$ を独立に組んで台帳と突き合わせる。

TP-1187 の分母は**同一パネルホルダ上の 2D 平板較正パネルの冷壁実測**で、Table III に
run ごとの絶対値が kW/m² で載っている (層流 6.55 / 9.08、乱流 66.85 / 70.25 / 136.98)。

ここでは Eckert 参照温度法 (層流 Blasius / 乱流 1/5 乗則) で独立に予測し、比を出す。
**本文は理論の偏りを明記している**:
    層流 — 理論が実測より約 10 % 高い
    乱流 — 理論が実測より約 30 % 高い
したがって本ゲートの合格条件は「1.0 に合うこと」ではなく、**この偏りを再現すること**。
再現できれば、Gate A の自由流・燃焼ガス物性・$q_{FP}$ の読み取りが一貫していることになる。

幾何 (Fig 6): 前縁から 位置 I = 117 cm、位置 II = 188 cm。乱流はトリップが前縁から 13 cm。
壁は冷壁 (模型は室温から挿入) なので $T_w$ = 300 K を既定とする。

**迎角を落とすと乱流が合わない**。Table II の乱流 run は $\alpha$=7.2-7.6° のものが多く、
自由流条件のまま計算すると理論/実測が 0.50-0.65 に落ちる ($\alpha\approx0.2°$ の run 8/14 だけ
1.21 で本文の 1.30 に近い)。パネルは平板なので**斜め衝撃波で圧縮された背後**が縁条件になる。
局所 $\gamma(T_\infty)$ の完全気体斜め衝撃波で縁条件を作り直す。
"""
import argparse, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from gas_htst import htst                                       # noqa: E402
import importlib.util as _ilu                                   # noqa: E402
_sp = _ilu.spec_from_file_location("t4_cond", HERE / "conditions.py")
_m = _ilu.module_from_spec(_sp); _sp.loader.exec_module(_m)
Gas = _m.Gas

X_LOC = {"I": 1.17, "II": 1.88}     # 前縁からの距離 [m]
X_TRIP = 0.13                        # 乱流トリップ位置 [m]


def oblique(M1, gam, theta_deg):
    """完全気体の斜め衝撃波。弱解の衝撃波角 beta を解き、背後の比を返す。"""
    th = np.radians(theta_deg)
    if th < 1e-6:
        return dict(p=1.0, T=1.0, ro=1.0, M2=M1, beta=np.nan)
    from scipy.optimize import brentq
    def f(b):
        return (2.0 / np.tan(b) * (M1 ** 2 * np.sin(b) ** 2 - 1.0)
                / (M1 ** 2 * (gam + np.cos(2.0 * b)) + 2.0) - np.tan(th))
    b = brentq(f, np.arcsin(1.0 / M1) + 1e-9, np.radians(60.0))
    Mn1 = M1 * np.sin(b)
    p = 1.0 + 2.0 * gam / (gam + 1.0) * (Mn1 ** 2 - 1.0)
    ro = (gam + 1.0) * Mn1 ** 2 / ((gam - 1.0) * Mn1 ** 2 + 2.0)
    T = p / ro
    Mn2 = np.sqrt((1.0 + 0.5 * (gam - 1.0) * Mn1 ** 2) / (gam * Mn1 ** 2 - 0.5 * (gam - 1.0)))
    return dict(p=p, T=T, ro=ro, M2=Mn2 / np.sin(b - th), beta=np.degrees(b))


def qdot(fs, gas, Tt, T_w, x, regime, alpha=0.0):
    """Eckert 参照温度法。冷壁熱流束 [W/m2]。alpha>0 なら斜め衝撃波で縁条件を作る。"""
    gam_inf = gas.cp(fs["T_inf"]) / (gas.cp(fs["T_inf"]) - gas.R)
    M1 = fs["U_inf"] / gas.a(fs["T_inf"])
    sh = oblique(M1, gam_inf, alpha)
    T_inf = fs["T_inf"] * sh["T"]
    p_e = fs["p_inf"] * sh["p"]
    U = sh["M2"] * gas.a(T_inf)
    r = np.sqrt(gas.Pr(T_inf)) if regime == "laminar" else 0.89
    T_aw = T_inf + r * (Tt - T_inf)
    T_ref = T_inf + 0.5 * (T_w - T_inf) + 0.22 * (T_aw - T_inf)
    ro_ref = p_e / (gas.R * T_ref)
    mu_ref, cp_ref, Pr_ref = gas.mu(T_ref), gas.cp(T_ref), gas.Pr(T_ref)
    Re_ref = ro_ref * U * x / mu_ref
    if regime == "laminar":
        cf2 = 0.332 / np.sqrt(Re_ref)
    else:
        cf2 = 0.0296 * Re_ref ** -0.2
    St = cf2 * Pr_ref ** (-2.0 / 3.0)
    return St * ro_ref * U * cp_ref * (T_aw - T_w), T_aw, T_ref, Re_ref, sh["beta"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tw", type=float, default=300.0, help="冷壁温度 [K]")
    a = ap.parse_args()
    cond = json.loads((HERE.parent / "conditions.json").read_text(encoding="utf-8"))
    der = json.loads((HERE.parent / "derived.json").read_text(encoding="utf-8"))
    by_run = {s["run"]: s for s in cond["series"]}

    print(f"=== T4 Gate B : 較正パネル q_FP の独立検算 (T_w = {a.tw:.0f} K) ===")
    print("本文の既知の偏り: 層流 理論が +10 %, 乱流 理論が +30 %\n")
    print(f"{'run':>4} {'流れ':>6} {'位置':>4} {'a[deg]':>7} {'beta':>6} {'q 理論':>8} {'q 実測':>8} "
          f"{'理論/実測':>9} | {'T_aw':>7} {'Re_ref':>10}")
    rows = []
    for o in der["series"]:
        s = by_run[o["run"]]
        for loc in ("I", "II"):
            q_meas = s.get(f"qfp_loc_{loc}_kW_m2")
            bl = s[f"bl_loc_{loc}"]
            if q_meas is None or bl not in ("LA", "T"):
                continue
            regime = "laminar" if bl == "LA" else "turbulent"
            gas = Gas(htst(s["Tt_c_K"]))
            x = X_LOC[loc] - (X_TRIP if regime == "turbulent" else 0.0)
            q, T_aw, T_ref, Re_ref, beta = qdot(o, gas, s["Tt_c_K"], a.tw, x, regime,
                                                alpha=s["alpha_deg"])
            ratio = q * 1e-3 / q_meas
            print(f"{o['run']:4d} {regime[:4]:>6} {loc:>4} {s['alpha_deg']:7.1f} {beta:6.1f} "
                  f"{q*1e-3:8.2f} {q_meas:8.2f} {ratio:9.3f} | {T_aw:7.1f} {Re_ref:10.3e}")
            rows.append(dict(run=o["run"], regime=regime, loc=loc, q_theory_kW=q * 1e-3,
                             q_meas_kW=q_meas, ratio=ratio, T_aw=T_aw, T_ref=T_ref,
                             Re_ref=Re_ref, alpha_deg=s["alpha_deg"], beta_deg=float(beta)))
    lam = np.array([r["ratio"] for r in rows if r["regime"] == "laminar"])
    tur = np.array([r["ratio"] for r in rows if r["regime"] == "turbulent"])
    print(f"\n層流 理論/実測 : {lam.min():.3f} - {lam.max():.3f} (中央 {np.median(lam):.3f})"
          f"   ← 本文 1.10")
    print(f"乱流 理論/実測 : {tur.min():.3f} - {tur.max():.3f} (中央 {np.median(tur):.3f})"
          f"   ← 本文 1.30")
    # 表示文と実装を一致させる (codex result M4: 旧実装は 0.85-1.45 / 0.95-1.75 で ±0.25 でなかった)
    ok_l = 1.10 - 0.25 < np.median(lam) < 1.10 + 0.25
    ok_t = 1.30 - 0.25 < np.median(tur) < 1.30 + 0.25
    verdict = "PASS" if (ok_l and ok_t) else "FAIL"
    print(f"GATE B VERDICT: {verdict}  (本文の偏り 1.10 / 1.30 を ±0.25 で挟めるか)")
    (HERE.parent / "gateB.json").write_text(json.dumps(
        dict(_gate="B", T_w=a.tw, verdict=verdict,
             laminar_median=float(np.median(lam)), turbulent_median=float(np.median(tur)),
             rows=rows), indent=2), encoding="utf-8")
    return 0 if verdict == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
