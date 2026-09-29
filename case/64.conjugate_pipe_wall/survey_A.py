#!/usr/bin/env python3
r"""A (厚肉管) の条件選定 (plan §5.1 #3b)。解析流速 (Poiseuille)・源項なしの参照解で、候補条件ごとに
上流へ回り込む熱の割合 Q_up/Q_total、固体の厚さ方向温度差 ΔT_s、固体断面の軸方向熱量 Q_ax を出す。forge の結果は見ない。"""
import itertools
import numpy as np
import conjugate_ref as cr

R, k_f, cp, rho, Um = 1e-3, 5.700284e-2, 1004.5, 1.176829, 17.35944
r_o = 2 * R
def run(ks_ratio, Lh_R, h, Tc, L_u, nx_scale=1.0, nf=32, ns=16):
    """外面の h・T_c を直接与える (A1/A2 で固定し k_s だけ変える: codex plan レビュー M3)。"""
    k_s = ks_ratio * k_f
    P, xs, ys = cr.pipe_problem(nf, ns, nx_scale, R, r_o, k_f, k_s, rho, cp, Um, L_u, Lh_R * R, 10e-3, (h, Tc))
    P.T_in = 0.0; T = P.solve()
    jw = nf
    # 界面熱流束 (固体 → 流体): 固体側の 2 次片側差分 q = k_s dT/dr|_R (外向き r 正、流体へは −) → 流体へ = k_s*(−dT/dr)… 符号: 固体は外から温められ内壁で流体へ放熱
    dr = ys[jw + 1] - ys[jw]
    q_i = -k_s * (-3 * T[:, jw] + 4 * T[:, jw + 1] - T[:, jw + 2]) / (2 * dr) * (-1)   # 流体へ向かう (−r 方向) 熱流束
    up = xs <= 0.0                                              # x=0 の節点を含めて積分端を閉じる (初版は最後の負の節点で切れていた)
    trap = getattr(np, "trapezoid", None) or np.trapz
    Q_up = R * trap(q_i[up], xs[up]); Q_tot = P.robin_heat()
    heat = (xs >= 0) & (xs <= Lh_R * R)
    dTs = (T[:, -1] - T[:, jw])                                  # 外面 − 内壁
    # 固体断面の軸方向熱量 (per rad): −∫ k_s dT/dx r dr
    Tx = np.gradient(T, xs, axis=0)
    Qax = np.array([-trap(k_s * Tx[i, jw:] * ys[jw:], ys[jw:]) for i in range(len(xs))])
    Tw_rise = trap(T[heat, jw], xs[heat]) / (xs[heat][-1] - xs[heat][0])     # 長さ重みの平均 (初版は節点の単純平均)
    return dict(Q_up_frac=Q_up / Q_tot, Tw_rise=Tw_rise, dTs_max=dTs[heat].max(), Qax_max=np.abs(Qax).max() / Q_tot,
                Tw_up_at_minus5R=np.interp(-5 * R, xs, T[:, jw]) / Tw_rise, xs=xs, T=T)

if __name__ == "__main__":
  print(" ks/kf  Lh/R  h    Q_up/Q_tot  Tw上昇[K]  ΔT_s max[K]  |Q_ax|max/Q_tot  x=−5R の予熱/上昇")
  for ksr, LhR, h in itertools.product((10, 100), (10,), (570.0,)):
    d = run(ksr, LhR, h, 10.0, 80e-3)
    print(f" {ksr:5d}  {LhR:4d}  {h:5.0f}   {d['Q_up_frac']:9.4f}   {d['Tw_rise']:8.3f}   {d['dTs_max']:9.4f}    {d['Qax_max']:9.4f}        {d['Tw_up_at_minus5R']:.4f}")
