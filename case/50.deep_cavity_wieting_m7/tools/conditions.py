#!/usr/bin/env python3
"""case/50 の条件再構成 — ゲート A (再現確認) と ゲート B (照合)。

plan `case-hypersonic-gap-heating-validation.md` §4.2:
  ゲート A: 組成と輸送則を固定したうえで、M・Re'・Tt から (T∞,p∞,ρ∞,U∞) を**一意に**決める。
           これは入力の再現であって検証ではない。
  ゲート B: 決めた状態から Pohlhausen + Eckert 基準温度で q_fp を計算し、Table IV と照合する
           (文献の q_fp も解析値なので「基準式と物性の整合確認」)。
  診断:    再構成した全圧と実測 pt の差を残す (組成誤差を全圧損失に吸収させない)。

使い方: python3 tools/conditions.py [--pr 0.75] [--json out.json]
"""
import argparse, json, math, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from gas_model import CombustionProducts                                   # noqa: E402

CASE = HERE.parent


def solve_state(gas, M, Re_m, Tt):
    """h(Tt) = h(T) + U²/2, U = M a(T) を T について解く → (T,U,rho,p)。"""
    def f(T):
        U = M * gas.a(T)
        return gas.h(Tt) - gas.h(T) - 0.5 * U * U
    lo, hi = 60.0, Tt - 1.0
    flo, fhi = f(lo), f(hi)
    if flo * fhi > 0:
        raise ValueError("解が範囲外")
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(lo) * f(mid) <= 0:
            hi = mid
        else:
            lo = mid
    T = 0.5 * (lo + hi)
    U = M * gas.a(T)
    mu = gas.mu(T)
    rho = Re_m * mu / U
    p = rho * gas.R * T
    return T, U, rho, p, mu


def pt_isentropic(gas, T, p, Tt, n=4000):
    """ln(pt/p) = ∫_T^Tt cp/(R T') dT' (熱的完全気体・等エントロピー)。"""
    Ts = np.linspace(T, Tt, n)
    cps = np.array([gas.cp(t) for t in Ts])
    integ = np.trapz(cps / (gas.R * Ts), Ts)
    return p * math.exp(integ)


def q_flat_plate(gas, T, U, p, x, Tw, pr_override=None, iters=6):
    """Pohlhausen 層流平板 + Eckert 基準温度。q_w [W/m²] と中間量を返す。"""
    Tstar = 0.5 * (T + Tw)
    Taw = T
    for _ in range(iters):
        Pr_s = pr_override if pr_override else gas.Pr(Tstar)
        r = math.sqrt(Pr_s)
        cp_s = gas.cp(Tstar)
        Taw = T + r * U * U / (2.0 * cp_s)
        Tstar = T + 0.5 * (Tw - T) + 0.22 * (Taw - T)
    Pr_s = pr_override if pr_override else gas.Pr(Tstar)
    mu_s = gas.mu(Tstar)
    cp_s = gas.cp(Tstar)
    rho_s = p / (gas.R * Tstar)
    Re_s = rho_s * U * x / mu_s
    St = 0.332 / math.sqrt(Re_s) * Pr_s ** (-2.0 / 3.0)
    q = St * rho_s * U * cp_s * (Taw - Tw)
    return dict(q=q, Taw=Taw, Tstar=Tstar, Pr_star=Pr_s, Re_star_x=Re_s, St=St, rho_star=rho_s)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pr", type=float, default=None,
                    help="Pr を固定する (例 0.75 = W70 本文の値)。既定はモデルの Pr(T*)")
    ap.add_argument("--json", default=None, help="結果を JSON で書き出す")
    a = ap.parse_args()

    cond = json.loads((CASE / "conditions.json").read_text(encoding="utf-8"))
    geom = json.loads((CASE / "geometry.json").read_text(encoding="utf-8"))
    Tw = cond["wall"]["Tw"]
    x_rear = geom["cavity"]["x_rear_wall_from_le"] * 1e-3
    widths = {k: v * 1e-3 for k, v in geom["cavity"]["widths"].items()}

    print(f"{'id':>3} {'w/d':>6} {'Tt[K]':>6} {'M':>5} {'Re/m':>9} | "
          f"{'T∞[K]':>7} {'p∞[Pa]':>8} {'ρ∞':>9} {'U∞[m/s]':>8} | "
          f"{'pt_rec/pt':>9} | {'Taw[K]':>7} {'T*[K]':>6} {'Pr*':>6} | "
          f"{'q calc':>7} {'q rep':>6} {'差%':>7}")
    out, gate_a_ok, gate_b_ok = [], True, True
    for s in cond["series"]:
        gas = CombustionProducts(s["Tt_K"], T_react=cond["gas"]["T_react"])
        T, U, rho, p, mu = solve_state(gas, s["M"], s["Re_m"], s["Tt_K"])
        pt_rec = pt_isentropic(gas, T, p, s["Tt_K"])
        w = widths[f"{s['w_over_d']:.3f}".rstrip('0').rstrip('.')] if False else widths[str(s["w_over_d"])]
        x_mid = x_rear - 0.5 * w
        fp = q_flat_plate(gas, T, U, p, x_mid, Tw, pr_override=a.pr)
        d = 100.0 * (fp["q"] / (s["qfp_kW"] * 1e3) - 1.0)
        gate_b_ok &= abs(d) <= 5.0
        print(f"{s['id']:>3} {s['w_over_d']:6.3f} {s['Tt_K']:6.0f} {s['M']:5.2f} {s['Re_m']:9.2e} | "
              f"{T:7.2f} {p:8.1f} {rho:9.5f} {U:8.1f} | "
              f"{pt_rec/(s['pt_MPa']*1e6):9.3f} | {fp['Taw']:7.1f} {fp['Tstar']:6.1f} {fp['Pr_star']:6.4f} | "
              f"{fp['q']*1e-3:7.2f} {s['qfp_kW']:6.1f} {d:+7.2f}")
        out.append(dict(id=s["id"], w_over_d=s["w_over_d"], Tt=s["Tt_K"], M=s["M"], Re_m=s["Re_m"],
                        phi=gas.phi, Y=gas.Y, R=gas.R, T_inf=T, p_inf=p, rho_inf=rho, U_inf=U,
                        mu_inf=mu, pt_reconstructed=pt_rec, pt_measured=s["pt_MPa"] * 1e6,
                        x_midpoint=x_mid, q_fp_calc=fp["q"], q_fp_reported=s["qfp_kW"] * 1e3,
                        diff_percent=d, **{k: fp[k] for k in ("Taw", "Tstar", "Pr_star", "Re_star_x")}))
    print()
    pt_ratios = [o["pt_reconstructed"] / o["pt_measured"] for o in out]
    print(f"ゲート A (再現確認): T∞/p∞/ρ∞/U∞ は M・Re'・Tt から一意に解けた ({len(out)} 系列) → PASS")
    print(f"  診断: 再構成 pt / 実測 pt = {min(pt_ratios):.3f} – {max(pt_ratios):.3f} "
          f"(1 から外れる分は全圧損失 + 組成・物性の誤差。**合わせに行かない**)")
    dmax = max(abs(o["diff_percent"]) for o in out)
    print(f"ゲート B (照合): q_fp の差 最大 {dmax:.2f} % → {'PASS' if gate_b_ok else 'FAIL'} (判定ライン ±5 %)")
    if a.pr:
        print(f"  (Pr を {a.pr} に固定した場合)")
    if a.json:
        Path(a.json).write_text(json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
        print(f"  → {a.json}")


if __name__ == "__main__":
    main()
