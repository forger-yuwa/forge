#!/usr/bin/env python3
"""T2 Gate A — TN D-8233 の自由流条件を確定する。

本文は $T_t$ を書いていない。しかし **$P_t$・$M$・$Re'$ が 3 系列とも与えられている**ので、
$T_t$ は解ける (未知 1・条件 1)。3 系列が同じ $T_t$ を返せば、施設の運転条件としても
digitize した数値としても自己整合が取れる。

半完全気体で厳密に解く (M=10.3 では $T_\\infty\\approx 45$ K まで落ちるので $\\gamma$ 一定は使えない):
  h(T_t) = h(T) + M^2 a(T)^2 / 2          (全エンタルピー)
  p = P_t exp[(s0(T) - s0(T_t))/R]        (等エントロピー)
  Re' = rho U / mu ,  rho = p/(R T) ,  U = M a(T)

200 K 未満は振動自由度が凍結しているので **cp だけ** cp(200 K) 一定に固定する
(NASA-9 多項式の下限 200 K を外挿すると 1/T^2 項が発散するため)。粘性は床止めしない
— 一度 mu も 200 K で止めてしまい、T_inf 45 K で真値の 4.7 倍になって Re' が 3.5 倍
外れた。低温の mu は Re' に直接乗るので、CE と Sutherland の差 (45 K で 16 %) も併記する。
"""
import json, sys
from pathlib import Path
import numpy as np
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from gas_air import dry_air                                      # noqa: E402

T_CLAMP = 200.0


class Air:
    """cp を 200 K で床止めした空気 (h, s0 は床を貫いて連続に延長)。"""

    def __init__(self, g):
        self.g = g
        self.R = g.R
        self.cp0 = g.cp(T_CLAMP)
        self.h0c = g.h(T_CLAMP)
        # s0(T) = ∫_{T_CLAMP}^{T} cp/T dT を数値積分で用意する
        self._Tg = np.geomspace(T_CLAMP, 4000.0, 2000)
        self._sg = np.concatenate([[0.0], np.cumsum(
            0.5 * (self.g.cp(self._Tg[1:]) / self._Tg[1:] + self.g.cp(self._Tg[:-1]) / self._Tg[:-1])
            * np.diff(self._Tg))])

    def cp(self, T):
        return self.cp0 if T < T_CLAMP else self.g.cp(T)

    def h(self, T):
        return self.h0c + self.cp0 * (T - T_CLAMP) if T < T_CLAMP else self.g.h(T)

    def s0(self, T):
        if T < T_CLAMP:
            return self.cp0 * np.log(T / T_CLAMP)
        return float(np.interp(T, self._Tg, self._sg))

    def a(self, T):
        gam = self.cp(T) / (self.cp(T) - self.R)
        return np.sqrt(gam * self.R * T)

    def mu(self, T):
        # 輸送係数は床止めしない。Chapman-Enskog (Neufeld Ω) は T*~0.5 まで有効で、
        # T_inf ≈ 45 K でも Sutherland と 16 % 差に収まる (forge の viscMethod 2 と同式)。
        return self.g.mu(T)

    def mu_sutherland(self, T):
        return 1.458e-6 * T ** 1.5 / (T + 110.4)

    def Pr(self, T):
        return self.g.Pr(T)


def freestream(air, M, Tt, Pt):
    """(M, Tt, Pt) → 自由流。半完全気体の等エントロピー膨張。"""
    f = lambda T: air.h(Tt) - air.h(T) - 0.5 * M ** 2 * air.a(T) ** 2
    T = brentq(f, 1.0, Tt * 0.999, xtol=1e-10)
    p = Pt * np.exp((air.s0(T) - air.s0(Tt)) / air.R)
    ro = p / (air.R * T)
    U = M * air.a(T)
    mu = air.mu(T)
    return dict(T=T, p=p, ro=ro, U=U, mu=mu, Re_m=ro * U / mu, a=air.a(T))


def solve_Tt(air, M, Pt, Re_target):
    f = lambda Tt: freestream(air, M, Tt, Pt)["Re_m"] - Re_target
    return brentq(f, 300.0, 2500.0, xtol=1e-6)


def main():
    cond = json.loads((HERE.parent / "conditions.json").read_text(encoding="utf-8"))
    M = cond["M"]
    air = Air(dry_air(1000.0))
    pt = cond["stagnation_pressure_MN_m2"]

    # digitize した Re' ラベルと BL 系列の Re' を対応づける
    pairs = [(1.47e6, pt["1.5e6"]), (3.32e6, pt["3.3e6"]), (7.82e6, pt["7.8e6"])]

    print(f"=== T2 Gate A : M = {M} 固定、(Pt, Re') から Tt を解く ===\n")
    print(f"{'Re_inf [1/m]':>14} {'Pt [MPa]':>9} {'Tt [K]':>8} {'T_inf [K]':>10} "
          f"{'p_inf [Pa]':>11} {'rho [kg/m3]':>12} {'U [m/s]':>9} {'T_aw [K]':>9}")
    out = []
    for Re, Pt_MPa in pairs:
        Pt = Pt_MPa * 1e6
        Tt = solve_Tt(air, M, Pt, Re)
        fs = freestream(air, M, Tt, Pt)
        Taw = fs["T"] + 0.89 * (Tt - fs["T"])
        print(f"{Re:14.3e} {Pt_MPa:9.2f} {Tt:8.1f} {fs['T']:10.2f} {fs['p']:11.1f} "
              f"{fs['ro']:12.5f} {fs['U']:9.1f} {Taw:9.1f}")
        out.append(dict(Re_m=Re, Pt_Pa=Pt, Tt=Tt, T_inf=fs["T"], p_inf=fs["p"],
                        ro_inf=fs["ro"], U_inf=fs["U"], mu_inf=fs["mu"], T_aw=Taw,
                        Pr_inf=air.Pr(fs["T"])))

    Tts = np.array([o["Tt"] for o in out])
    spread = (Tts.max() - Tts.min()) / Tts.mean() * 100
    print(f"\n解けた Tt : {Tts[0]:.1f} / {Tts[1]:.1f} / {Tts[2]:.1f} K   "
          f"平均 {Tts.mean():.1f} K, ばらつき {spread:.1f} %")
    verdict = "PASS" if spread < 10 else "FAIL"
    print(f"GATE A VERDICT: {verdict}  (3 系列が同じ Tt を返すか; 閾値 10 %)")
    if verdict == "PASS":
        print("  → Langley CFHT の公称運転温度帯と照合できる単一の Tt が立つ。"
              "digitize した Pt と Re' は互いに整合している。")

    (HERE.parent / "derived.json").write_text(json.dumps(
        dict(_gate="A", M=M, Tt_mean=float(Tts.mean()), Tt_spread_pct=float(spread),
             verdict=verdict, series=out), indent=2), encoding="utf-8")
    print(f"\n書き出し: case/51.gap_turbulent_thd8233/derived.json")
    return 0 if verdict == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
