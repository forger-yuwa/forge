#!/usr/bin/env python3
r"""case/64 厚肉管の共役 Graetz (A1/A2) の登録値と共通部品。plan
[`boundary-cht-conjugate-benchmarks.md`](../../plans/active/boundary-cht-conjugate-benchmarks.md) §4.2・§4.5。
流体の物性・流れは case/63 と同じ (graetz_common を流用)。"""
from __future__ import annotations
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "63.graetz_cht"))
import graetz_common as gc          # noqa: E402  (R, D, MU, K_F, CP, RHO, U_M, T_IN, P_OUT, axcht)
axcht = gc.axcht

R = gc.R
R_O = 2.0 * R                        # 固体殻の外半径
L_UP, L_HEAT, L_DOWN = 80.0 * R, 10.0 * R, 10.0 * R
H_O, DT_C = 570.0, 10.0              # 外面 Robin (加熱区間だけ): h_o [W/m²K]、T_c − T_in [K]
KS_RATIO = {"A1": 10.0, "A2": 100.0}
PID = {"inlet": 1, "outlet": 2, "wall": 3, "axis": 4}
