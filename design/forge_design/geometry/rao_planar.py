r"""平面 (2D) 最大推力理論の断片: 縁条件・理想推力・Rao 軌跡。

plan: plans/active/tooling-nozzle-sern-chain.md §4.3。
Guderley–Hantsch (1955) / Rao (1958) の長さ拘束最大推力の平面版。制御面 = ランプ後縁 e から
出る最終 C⁺ (c–e)。平面では C⁺ 上の適合関係 θ−ν = const と乗数関係が位置に依らないので
制御面上の状態は一様 (直線 Mach 線) になり、縁 e では

    ½ ρ_e w_e² sin 2θ_e = (p_e − p_a) cot μ_e      (Cain 2010, 式 4.3)

すなわち sin 2θ_e = 2 (1 − p_a/p_e) cot μ_e / (γ M_e²)。p_e = p_a で θ_e = 0 (一様平行流、
非粘性の絶対最大)、p_e > p_a (短いノズル) で θ_e > 0。
理論本文 (Rao) は未入手のため、これらは調査ノート (sern-design-method-survey.md §5) の通り
自前導出扱い — テスト (run_sern_moc_tests.py) の掃引で最適性を数値検算する。
"""
from __future__ import annotations

import numpy as np

from .moc_kernel import pm_nu, mu_vec


def p_over_pt(M, g=1.4):
    M = np.asarray(M, dtype=float)
    return (1.0 + 0.5 * (g - 1.0) * M * M) ** (-g / (g - 1.0))


def mach_from_p_over_pt(r, g=1.4):
    r = np.asarray(r, dtype=float)
    return np.sqrt(2.0 / (g - 1.0) * (r ** (-(g - 1.0) / g) - 1.0))


def massflux_ratio(M, M_ref, g=1.4):
    """ρw / (ρw)_ref (同一全状態)。"""
    def f(m):
        return m * (1.0 + 0.5 * (g - 1.0) * m * m) ** (-(g + 1.0) / (2.0 * (g - 1.0)))
    return f(np.asarray(M, dtype=float)) / f(M_ref)


def theta_lip(M_e, pe_over_pa, g=1.4):
    """縁条件から θ_e [rad]。pe/pa < 1 (過膨張) は 0 を返す (理論の適用外)。"""
    M_e = float(M_e)
    r = 1.0 - 1.0 / float(pe_over_pa)
    if r <= 0.0:
        return 0.0
    s = 2.0 * r / np.tan(float(mu_vec(M_e))) / (g * M_e * M_e)
    if s >= 1.0:
        return 0.25 * np.pi
    return 0.5 * float(np.arcsin(s))


def ideal_gross_thrust(M_in, pa_over_pin, g=1.4):
    """入口一様流 (p_in=1, H=1) を p_a まで等エントロピー膨張した理想総推力 F/(p_in H) = ṁ u_e。"""
    pe_pt = pa_over_pin * float(p_over_pt(M_in, g))
    M_e = float(mach_from_p_over_pt(pe_pt, g))
    T_ratio = (1.0 + 0.5 * (g - 1.0) * M_in ** 2) / (1.0 + 0.5 * (g - 1.0) * M_e ** 2)
    u_ratio = M_e / M_in * np.sqrt(T_ratio)
    return g * M_in ** 2 * u_ratio, M_e


def inlet_stream_thrust(M_in, pa_over_pin, g=1.4):
    """入口面の流れ推力 (ṁ u + (p − p_a) H)/(p_in H)。"""
    return g * M_in ** 2 + (1.0 - pa_over_pin)


def lip_residual(M_c, theta_c, p_c_over_pa, g=1.4):
    """key point 状態の縁条件残差 θ_c − θ_lip [rad] (0 が Rao 最適)。"""
    return float(theta_c) - theta_lip(M_c, p_c_over_pa, g)
