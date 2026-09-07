"""⑤ SERN の幾何コンテナ (無次元: 燃焼器出口高さ H = 1、入口面 x = 0)。

座標: x = 流れ方向、y = 上向き。ランプ = 上壁 ((0,1) の角部から θ_r0 で膨張)、
カウル = 下壁 ((0,0) から −θ_c0、長さ L_cowl)。角度は radian。
plan: plans/active/tooling-nozzle-sern-chain.md §4.1。
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class SernEnvelope:
    """幾何包絡 (spec)。None は無制約。"""
    L_ramp_max: float | None = None
    L_cowl_max: float | None = None
    y_ramp_max: float | None = None   # 後胴線 (定数近似)。ランプがこれを越えたら違反


@dataclass
class SernDesign:
    """逆設計の結果 (moc_sern.design_ramp の戻り値)。全て無次元。"""
    ramp_xy: np.ndarray           # (n,2) 角部 → a → … → e
    ramp_M: np.ndarray            # 壁上の M (同じ長さ)
    ramp_theta: np.ndarray        # 壁角 [rad]
    cowl_xy: np.ndarray           # (m,2) (0,0) → TE
    cowl_M: np.ndarray
    cowl_theta: np.ndarray
    key_point: tuple              # (x_c, y_c, M_c, theta_c)
    foot_a: tuple                 # (x_a, y_a) C⁻ の壁足
    lip_e: tuple                  # (x_e, y_e)
    mass_fraction: float          # 指定 f
    mass_fraction_check: float    # 終端特性線 c–e を横切る質量流量比 (検算)
    p_ext_over_p_in: float
    info: dict = field(default_factory=dict)

    @property
    def L_ramp(self) -> float:
        return float(self.ramp_xy[-1, 0])

    @property
    def H_exit(self) -> float:
        """ランプ後縁とカウル後縁 (または自由境界) の高さ差 (参考値)。"""
        return float(self.ramp_xy[-1, 1] - self.cowl_xy[-1, 1])

    def check_envelope(self, env: SernEnvelope) -> list[str]:
        v = []
        if env.L_ramp_max is not None and self.L_ramp > env.L_ramp_max:
            v.append(f"L_ramp {self.L_ramp:.4f} > L_ramp_max {env.L_ramp_max:.4f}")
        if env.L_cowl_max is not None and self.cowl_xy[-1, 0] > env.L_cowl_max + 1e-12:
            v.append(f"L_cowl {self.cowl_xy[-1, 0]:.4f} > L_cowl_max {env.L_cowl_max:.4f}")
        if env.y_ramp_max is not None and float(self.ramp_xy[:, 1].max()) > env.y_ramp_max:
            v.append(f"ramp y_max {self.ramp_xy[:, 1].max():.4f} > y_ramp_max {env.y_ramp_max:.4f}")
        return v
