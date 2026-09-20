#!/usr/bin/env python3
"""T2-G0 — 試験の流入乱流 BL を再構成する。

すきまに近づくのは**トンネル壁の厚い乱流 BL** ($\\delta^*$ = 10.5-12.8 cm = すきま幅の 46-56 倍)
なので、前縁から発達させず入口にプロファイルを与える (§4.2b)。

**$C_f$ を仮定しない**。Gate B で報告の平板相関と Van Driest II が 1.4 倍・冪も違うと
分かったので、どちらかを入力にすると流入場がその選択に汚染される。代わりに
**Fig 4 の実測 $\\delta^*$ と $\\theta$ (2 条件) から $(\\delta, u_\\tau)$ を解く**。
得られた $C_f$ は出力であり、VD-II と相関のどちらに近いかを**判定材料**にできる。

構成:
  速度  van Driest 変換した $u_{vd}^+$ に Spalding 内層 + Coles wake
  温度  Walz / Crocco-Busemann (回復係数 r=0.89)
  密度  等圧 $\\rho = p/(RT)$
  k, ω  対数層の平衡値を壁で漸近させる
"""
import json, sys
from pathlib import Path
import numpy as np
from scipy.optimize import brentq, fsolve

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from gas_air import dry_air                                      # noqa: E402
import importlib.util as _ilu                                    # noqa: E402
_sp = _ilu.spec_from_file_location("t2_conditions", HERE / "conditions.py")
_m = _ilu.module_from_spec(_sp); _sp.loader.exec_module(_m)
Air = _m.Air

KAPPA, B_LOG, BETA_STAR = 0.41, 5.2, 0.09
R_FACTOR = 0.89


def spalding_table(u_max=40.0, n=4000):
    """Spalding の内層則 y+(u+) を表にする (u+ 単調なので後で反転できる)。"""
    up = np.linspace(1e-6, u_max, n)
    ku = KAPPA * up
    yp = up + np.exp(-KAPPA * B_LOG) * (
        np.exp(ku) - 1.0 - ku - ku ** 2 / 2.0 - ku ** 3 / 6.0)
    return yp, up


_YP, _UP = spalding_table()


def u_inner(yp):
    return np.interp(yp, _YP, _UP)


def wake(eta, Pi):
    return (2.0 * Pi / KAPPA) * np.sin(0.5 * np.pi * np.clip(eta, 0, 1)) ** 2


class Profile:
    """(delta, u_tau) を与えると BL 全体を返す。"""

    def __init__(self, fs, T_w, T_aw, air, M):
        self.fs, self.T_w, self.T_aw, self.air, self.M = fs, T_w, T_aw, air, M
        self.ro_w = fs["p"] / (air.R * T_w)
        self.mu_w = air.mu(T_w)
        # Crocco:  T/T_w = 1 + B xi - A^2 xi^2
        # Crocco: T/T_w = 1 + Bc xi - A2 xi^2。**端点整合で書く** (xi=1 で T=T_inf)。
        # gamma 一定の r(g-1)/2 M^2 T_inf/T_w を使うと半完全気体の T_aw と矛盾し、
        # 端で T が 0 を切って arcsin が定義域を外れる。
        self.Bc = (T_aw - T_w) / T_w
        self.A2 = (T_aw - fs["T"]) / T_w
        self.A = np.sqrt(self.A2)
        self.D = np.sqrt(self.Bc ** 2 + 4.0 * self.A2)
        # van Driest 変換の端値
        self.uvd_e = (fs["U"] / self.A) * (
            np.arcsin((2.0 * self.A2 - self.Bc) / self.D) + np.arcsin(self.Bc / self.D))

    def xi_of_uvd(self, uvd):
        """van Driest 速度 -> xi = u/Ue (Crocco の逆変換)。"""
        s = np.sin(self.A * uvd / self.fs["U"] - np.arcsin(self.Bc / self.D))
        return (self.D * s + self.Bc) / (2.0 * self.A2)

    def build(self, delta, u_tau, ny=600):
        dp = delta * u_tau * self.ro_w / self.mu_w              # delta+
        Pi = KAPPA * 0.5 * (self.uvd_e / u_tau - u_inner(dp))   # 端で閉じるよう Pi を決める
        y = np.concatenate([[0.0], np.geomspace(1e-9, delta, ny - 1)])
        yp = y * u_tau * self.ro_w / self.mu_w
        uvd = u_tau * (u_inner(yp) + wake(y / delta, Pi))
        uvd = np.minimum(uvd, self.uvd_e)
        xi = np.clip(self.xi_of_uvd(uvd), 0.0, 1.0)
        u = xi * self.fs["U"]
        T = self.T_w * (1.0 + self.Bc * xi - self.A2 * xi ** 2)
        ro = self.fs["p"] / (self.air.R * T)
        return dict(y=y, u=u, T=T, ro=ro, xi=xi, Pi=Pi, delta_plus=dp,
                    delta=delta, u_tau=u_tau)

    def thicknesses(self, pr):
        roe, Ue = self.fs["ro"], self.fs["U"]
        f = pr["ro"] * pr["u"] / (roe * Ue)
        dstar = np.trapz(1.0 - f, pr["y"])
        theta = np.trapz(f * (1.0 - pr["u"] / Ue), pr["y"])
        return dstar, theta

    def solve(self, dstar_t, theta_t, guess=(0.30, 50.0)):
        def res(p):
            d, ut = abs(p[0]), abs(p[1])
            ds, th = self.thicknesses(self.build(d, ut))
            return [ds / dstar_t - 1.0, th / theta_t - 1.0]
        sol = fsolve(res, guess, full_output=True)
        p, info, ier, msg = sol
        d, ut = abs(p[0]), abs(p[1])
        return d, ut, ier, np.max(np.abs(res(p)))

    def solve_one(self, target, which, u_tau):
        """u_tau を外から固定し、delta だけを 1 つの積分厚さに合わせる。"""
        def f(d):
            ds, th = self.thicknesses(self.build(d, u_tau))
            return (ds if which == "dstar" else th) - target
        return brentq(f, 1e-3, 3.0, xtol=1e-9)

    def turbulence(self, pr):
        """平衡対数層の k, ω。壁では omega の粘性極限、外縁では小さく落とす。"""
        ut, y = pr["u_tau"], pr["y"]
        yp = np.maximum(y * ut * self.ro_w / self.mu_w, 1e-12)
        damp = (1.0 - np.exp(-yp / 25.0)) ** 2                 # 粘性底層で k を落とす
        k = (ut ** 2 / np.sqrt(BETA_STAR)) * damp * np.clip(1.0 - y / pr["delta"], 0, 1) ** 2
        nu_w = self.mu_w / self.ro_w
        om_vis = 6.0 * nu_w / (0.075 * np.maximum(y, 1e-9) ** 2)
        om_log = ut / (np.sqrt(BETA_STAR) * KAPPA * np.maximum(y, 1e-9))
        om = np.minimum(om_vis, 1e30)
        om = np.sqrt(om_log ** 2 + np.minimum(om_vis, 1e12) ** 2)
        k = np.maximum(k, 1e-6 * self.fs["U"] ** 2 * 1e-6)
        om = np.maximum(om, 1e-3)
        return k, om


def main():
    cond = json.loads((HERE.parent / "conditions.json").read_text(encoding="utf-8"))
    der = json.loads((HERE.parent / "derived.json").read_text(encoding="utf-8"))
    air = Air(dry_air(1000.0))
    M = cond["M"]
    T_w = float(sys.argv[1]) if len(sys.argv) > 1 else 300.0

    print(f"=== T2-G0 : 実測 δ*, θ から流入 BL を再構成 (T_w = {T_w:.0f} K) ===\n")
    print(f"{'Re_inf':>10} {'δ* 目標':>9} {'θ 目標':>8} | {'δ [cm]':>8} {'u_tau':>7} "
          f"{'Π':>6} {'δ+':>9} {'Cf':>9} {'Cf/VDII':>8} {'Cf/相関':>8} {'残差':>8}")
    rows = []
    for s, bl in zip(der["series"], [cond["boundary_layer"]["series"][i] for i in (0, 1, 3)]):
        fs = dict(ro=s["ro_inf"], U=s["U_inf"], mu=s["mu_inf"], T=s["T_inf"], p=s["p_inf"])
        prof = Profile(fs, T_w, s["T_aw"], air, M)
        ds_t, th_t = bl["dstar_cm"] * 1e-2, bl["theta_cm"] * 1e-2
        d, ut, ier, resid = prof.solve(ds_t, th_t)
        pr = prof.build(d, ut)
        tau_w = prof.ro_w * ut ** 2
        Cf = 2.0 * tau_w / (fs["ro"] * fs["U"] ** 2)
        # 比較対象
        import importlib.util as il
        sp = il.spec_from_file_location("fp", HERE / "flatplate.py")
        fpm = il.module_from_spec(sp); sp.loader.exec_module(fpm)
        Cf_vd, _, _, _ = fpm.van_driest_ii(M, fs["T"], T_w, air,
                                           fs["ro"] * fs["U"] * th_t / fs["mu"],
                                           T_aw=s["T_aw"])
        # 相関の h から逆算した等価 Cf (Reynolds 相似)
        C = cond["reference"]["C_digitized"]
        h_corr = C * s["Re_m"] ** 0.69
        T_ref = fs["T"] + 0.5 * (T_w - fs["T"]) + 0.22 * (s["T_aw"] - fs["T"])
        Cf_corr = 2.0 * h_corr / (fs["ro"] * fs["U"] * air.cp(T_ref)) * air.Pr(T_ref) ** (2 / 3)
        print(f"{s['Re_m']:10.2e} {ds_t*100:9.2f} {th_t*100:8.2f} | {d*100:8.2f} {ut:7.2f} "
              f"{pr['Pi']:6.3f} {pr['delta_plus']:9.2e} {Cf:9.3e} {Cf/Cf_vd:8.3f} "
              f"{Cf/Cf_corr:8.3f} {resid:8.1e}")
        rows.append(dict(Re_m=s["Re_m"], delta=d, u_tau=ut, Pi=float(pr["Pi"]),
                         delta_plus=float(pr["delta_plus"]), Cf=float(Cf),
                         Cf_over_VDII=float(Cf / Cf_vd), Cf_over_corr=float(Cf / Cf_corr),
                         dstar_target=ds_t, theta_target=th_t, residual=float(resid),
                         T_w=T_w, converged=int(ier)))

    print(f"\n判定材料:")
    print(f"  Π (Coles) = {[round(r['Pi'],3) for r in rows]}  "
          f"— 平衡平板の標準値は 0.4-0.7")
    print(f"  Cf / VD-II  = {[round(r['Cf_over_VDII'],3) for r in rows]}")
    print(f"  Cf / 相関等価 = {[round(r['Cf_over_corr'],3) for r in rows]}")
    ok = all(r["residual"] < 1e-6 for r in rows)
    print(f"\nG0-a (δ*, θ の同時一致): {'PASS' if ok else 'FAIL'}  "
          f"(最大残差 {max(r['residual'] for r in rows):.1e})")

    (HERE.parent / "inlet_bl.json").write_text(json.dumps(
        dict(_gate="T2-G0a", T_w=T_w, rows=rows), indent=2), encoding="utf-8")
    print("書き出し: case/51.gap_turbulent_thd8233/inlet_bl.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
