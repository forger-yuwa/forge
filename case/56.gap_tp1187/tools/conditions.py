#!/usr/bin/env python3
"""T4 Gate A — 8-ft HTST の自由流を、台帳の過決定性で検算しながら確定する。

Table II は run ごとに $T_{t,c}$, $P_{t,c}$, 動圧 $q$, 単位 Reynolds 数 $R$ を**両方**載せている。

**燃焼器全圧 $P_{t,c}$ をそのまま等エントロピー膨張させると合わない**。$M$ を自由にして
$q$ と $Re'$ から別々に解くと 9.2 と 10.5 になり (公称 7.0)、互いにも 16 % 食い違う。
これは 8-ft HTST が燃焼駆動の blowdown で、**燃焼器から試験部までに大きな全圧損失がある**ため。

そこで **$M$=7 を固定し、$q$ から試験部全圧 $P_{t,t}$ を解いて、$Re'$ を予測して実測と比べる**
(未知 1・条件 2 の過決定)。$P_{t,t}/P_{t,c}$ が 17 run で一定なら、固定ノズルの回復率として
筋が通る — これが本ゲートの判定。

気体は case/50 と同じ `CombustionProducts` (メタン-空気燃焼生成物、当量比は $T_{t,c}$ から逆算)。
半完全気体で等エントロピー膨張を厳密に解く。
"""
import json, sys
from pathlib import Path
import numpy as np
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from gas_htst import htst                                       # noqa: E402

T_CLAMP = 200.0


class Gas:
    """cp を 200 K で床止め (NASA-9 下限)。輸送は床止めしない。"""

    def __init__(self, g):
        self.g = g; self.R = g.R
        self.cp0 = g.cp(T_CLAMP); self.h0c = g.h(T_CLAMP)
        self._Tg = np.geomspace(T_CLAMP, 4000.0, 2000)
        self._sg = np.concatenate([[0.0], np.cumsum(
            0.5 * (g.cp(self._Tg[1:]) / self._Tg[1:] + g.cp(self._Tg[:-1]) / self._Tg[:-1])
            * np.diff(self._Tg))])

    def cp(self, T): return self.cp0 if T < T_CLAMP else self.g.cp(T)
    def h(self, T):  return self.h0c + self.cp0 * (T - T_CLAMP) if T < T_CLAMP else self.g.h(T)
    def s0(self, T): return self.cp0 * np.log(T / T_CLAMP) if T < T_CLAMP else float(
        np.interp(T, self._Tg, self._sg))
    def a(self, T):
        gam = self.cp(T) / (self.cp(T) - self.R); return np.sqrt(gam * self.R * T)
    def mu(self, T): return self.g.mu(T)
    def Pr(self, T): return self.g.Pr(T)


def freestream(gas, M, Tt, Pt):
    f = lambda T: gas.h(Tt) - gas.h(T) - 0.5 * M ** 2 * gas.a(T) ** 2
    T = brentq(f, 1.0, Tt * 0.999, xtol=1e-10)
    p = Pt * np.exp((gas.s0(T) - gas.s0(Tt)) / gas.R)
    ro = p / (gas.R * T); U = M * gas.a(T); mu = gas.mu(T)
    return dict(T=T, p=p, ro=ro, U=U, mu=mu,
                q_dyn=0.5 * ro * U ** 2, Re_m=ro * U / mu)


def solve_Ptt(gas, M, Tt, q_dyn):
    """M 固定で、動圧に合う試験部全圧を解く。"""
    f = lambda P: freestream(gas, M, Tt, P)["q_dyn"] - q_dyn
    return brentq(f, 1e4, 5e7, xtol=1.0)


def main():
    cond = json.loads((HERE.parent / "conditions.json").read_text(encoding="utf-8"))
    M = cond["facility"]["M_nominal"]
    print(f"=== T4 Gate A : M = {M} 固定、q から試験部全圧 Pt_t を解き Re' で検算 ===\n")
    print(f"{'run':>4} {'Pt_c':>6} {'Pt_t':>7} {'Pt_t/Pt_c':>10} | {'Re 予測':>10} {'Re 実測':>10} "
          f"{'比':>6} | {'T_inf':>7} {'p_inf':>8} {'U_inf':>7} {'T_aw':>7}")
    out = []
    for s_ in cond["series"]:
        Tt, Ptc = s_["Tt_c_K"], s_["Pt_c_MPa"] * 1e6
        gas = Gas(htst(Tt))
        Ptt = solve_Ptt(gas, M, Tt, s_["q_dyn_kPa"] * 1e3)
        fs = freestream(gas, M, Tt, Ptt)
        Taw = fs["T"] + 0.89 * (Tt - fs["T"])
        rec, rr = Ptt / Ptc, fs["Re_m"] / s_["Re_m"]
        print(f"{s_['run']:4d} {Ptc*1e-6:6.1f} {Ptt*1e-6:7.3f} {rec:10.3f} | "
              f"{fs['Re_m']:10.3e} {s_['Re_m']:10.3e} {rr:6.3f} | "
              f"{fs['T']:7.1f} {fs['p']:8.1f} {fs['U']:7.1f} {Taw:7.1f}")
        out.append(dict(run=s_["run"], M=M, Pt_test_Pa=Ptt, recovery=rec,
                        Re_pred=fs["Re_m"], Re_ledger=s_["Re_m"], Re_ratio=rr,
                        T_inf=fs["T"], p_inf=fs["p"], ro_inf=fs["ro"], U_inf=fs["U"],
                        mu_inf=fs["mu"], T_aw=Taw, Tt_c=Tt, Pt_c_Pa=Ptc,
                        Pr_inf=gas.Pr(fs["T"])))
    rec = np.array([o["recovery"] for o in out])
    rr = np.array([o["Re_ratio"] for o in out])
    spread = 100 * (rec.max() - rec.min()) / np.median(rec)
    print(f"\nPt_t/Pt_c : {rec.min():.3f} - {rec.max():.3f} (中央 {np.median(rec):.3f}), "
          f"ばらつき {spread:.0f} %  ← Pt_c は 4.2-17.3 MPa と 4 倍振れている")
    print(f"Re 予測/実測 : {rr.min():.3f} - {rr.max():.3f} (中央 {np.median(rr):.3f})")
    worst = out[int(np.argmax(np.abs(rr - 1.0)))]
    print(f"最悪 run {worst['run']}: Re 比 {worst['Re_ratio']:.3f}")
    verdict = "PASS" if spread < 30 and abs(np.median(rr) - 1.0) < 0.15 else "FAIL"
    print(f"GATE A VERDICT: {verdict}  (回復率のばらつき < 30 % かつ Re 比の中央が 1 の 15 % 以内)")
    if verdict == "PASS":
        print("  → 固定ノズルの全圧回復率として一定で、独立な Re' も再現する。"
              "台帳・気体モデル・M=7 が互いに整合している。")
    (HERE.parent / "derived.json").write_text(json.dumps(
        dict(_gate="A", M=M, verdict=verdict, recovery_median=float(np.median(rec)),
             Re_ratio_median=float(np.median(rr)), series=out), indent=2), encoding="utf-8")
    print("\n書き出し: case/56.gap_tp1187/derived.json")
    return 0 if verdict == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
