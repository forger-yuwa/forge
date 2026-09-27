#!/usr/bin/env python3
r"""V-g1 の参考比較用: 有限 Pe (720) の ellip 解の加熱開始断面 x⁺ = 0 の θ(ξ) を CSV に書く (合否に使わない)。

    python3 case/63.graetz_cht/ref_x0_profile.py     # → ref_x0_profile.csv (xi, theta)。T = T_c − ΔT θ
"""
import numpy as np

import graetz_ref

x, nu, tb, xi, th = graetz_ref.solve_ellip(120, 720.0, 0.01, 0.12, 0.01, 2400, 400, 400, return_field=True)
i0 = int(np.argmin(np.abs(x)))
assert abs(x[i0]) < 1e-15
np.savetxt("ref_x0_profile.csv", np.c_[xi, th[i0]], delimiter=",", header="xi,theta", comments="", fmt="%.10e")
print(f"x⁺=0 断面: θ 軸 {th[i0,0]:.6f} / ξ=0.9 {np.interp(0.9, xi, th[i0]):.6f} / 壁 {th[i0,-1]:.6f}")
