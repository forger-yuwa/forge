#!/usr/bin/env python3
r"""case/63 Graetz (軸対称 CHT・流れあり) の共通定数と部品。

plan [`boundary-cht-axisymmetric-graetz.md`](../../plans/accepted/boundary-cht-axisymmetric-graetz.md) §4。
case/61 の `axcht.py` (変換・固体帯・壁ダンプ・固体 Q_sol) を import して流用する。
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "61.conjugate_annulus"))
import axcht  # noqa: E402  (流用: gmsh_and_convert / solid_strip / write_solid / mesh_quality / 壁ダンプ)

# ---- 登録値 (plan §4.1) ----
R = 1.0e-3                 # 管半径 [m]
D = 2.0 * R
L_UP, L_HEAT, L_DOWN = 10.0e-3, 0.1728, 10.0e-3   # 上流断熱 / 加熱 (x⁺ 0.12) / 下流断熱 [m]
CP, GAMMA = 1004.5, 1.4
RGAS = CP * (GAMMA - 1.0) / GAMMA
T_IN, P_OUT = 300.0, 101325.0
MACH, RE, PR = 0.05, 1000.0, 0.72
RHO = P_OUT / (RGAS * T_IN)
C_IN = (GAMMA * RGAS * T_IN) ** 0.5
U_M = MACH * C_IN
MU = RHO * U_M * D / RE
K_F = MU * CP / PR
PE = RE * PR

# 物理 ID (gen_mesh.py の Physical Curve)
PID = {"inlet": 1, "outlet": 2, "wall_up": 3, "wall_heat": 4, "wall_down": 5, "axis": 6}


def summary() -> str:
    return (f"R {R*1e3:g} mm, M {MACH}, Re {RE:g}, Pr {PR}, Pe {PE:g}: rho {RHO:.6f}, c {C_IN:.4f}, "
            f"U_m {U_M:.5f} m/s, mu {MU:.6e}, k_f {K_F:.6e}")


if __name__ == "__main__":
    print(summary())
