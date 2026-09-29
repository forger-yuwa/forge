#!/usr/bin/env python3
r"""case/65 共役平板 (C1/C2) の登録値。plan [`boundary-cht-conjugate-benchmarks.md`](../../plans/active/boundary-cht-conjugate-benchmarks.md) §4.3・§4.5。
流体の物性は case/63 と同じ定数物性 (Pr 0.72)、M 0.1、Re_L 1e4。"""
from __future__ import annotations
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "63.graetz_cht"))
import graetz_common as gc          # noqa: E402
axcht = gc.axcht
U_INF = 0.1 * gc.C_IN                # M 0.1
L = 1.0e4 * gc.MU / (gc.RHO * U_INF) # Re_L = 1e4 → 10 mm
B = 0.2 * L                          # 板厚 b/L 0.2
L_UP, L_DOWN = 0.5 * L, 0.5 * L
H_TOP = 6 * 5 * L / 1.0e4 ** 0.5     # 6 δ(L)
H_BACK, DT_H = 1.0e8, 10.0           # 下面 Robin
KS_RATIO = {"C1": 10.0, "C2": 100.0}
PID = {"inlet": 1, "outlet": 2, "top": 3, "slip_up": 4, "plate": 5, "slip_down": 6}
