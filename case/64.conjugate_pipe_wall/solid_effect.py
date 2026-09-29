#!/usr/bin/env python3
r"""case/64 (A): forge 側の固体の効果の指標 (plan §4.4・§6「固体が効いていることの確認」)。
固体ダンプ `res_solid_3_<step>.h5` の温度から、(i) 加熱区間の厚さ方向の温度差 ΔT_s = T(r_o) − T(R) の最大、
(ii) 固体断面を通る軸方向熱量 Q_ax(x) = −∫ k_s ∂T/∂x r dr (per rad) の最大を Q_tot (外面の Robin 入熱) で割った値、
(iii) 上流へ回り込む熱 Q_up/Q_tot を出し、参照解 (eval_conj の主参照) の同じ量と並べる。

    python3 solid_effect.py <run>
"""
from __future__ import annotations
import re, sys
from pathlib import Path
import h5py
import numpy as np
import pipe_common as pc
import eval_conj as ev

trap = getattr(np, "trapezoid", None) or np.trapz


def solid_indicators(xs, rs, T, k_s):
    """テンソル格子 (x, r) の固体温度から ΔT_s と Q_ax を作る。"""
    heat = (xs >= -1e-12) & (xs <= pc.L_HEAT + 1e-9)
    dTs = (T[:, -1] - T[:, 0])[heat].max()
    Tx = np.gradient(T, xs, axis=0)
    Qax = np.array([-trap(k_s * Tx[i] * rs, rs) for i in range(len(xs))])
    return dTs, Qax


def main():
    run = Path(sys.argv[1])
    case = re.search(r"kind cht (A1|A2)", (run / "RUN_INPUTS.txt").read_text()).group(1); pc.CASE = case
    k_s = pc.KS_RATIO[case] * pc.gc.K_F
    st = ev.last_step(run)
    with h5py.File(run / "solid.h5", "r") as f:
        C = np.asarray(f["MESH/COORD"][:], float); perm = np.asarray(f["MESH/PERM"][:], int)
    with h5py.File(run / f"res_solid_3_{st}.h5", "r") as f:
        Tn = np.asarray(f["VALUE/T"][:], float); qh = float(np.asarray(f["VALUE/q_hole"][:], float).sum())
    xr, rr = np.round(C[:, 0], 10), np.round(C[:, 1], 10)
    xs, rs = np.unique(xr), np.unique(rr)
    if len(xs) * len(rs) != len(C): print("REFUSED: 固体がテンソル格子でない"); return 2
    Tg = np.full((len(xs), len(rs)), np.nan); Tg[np.searchsorted(xs, xr), np.searchsorted(rs, rr)] = Tn
    if np.isnan(Tg).any():                                   # 固体ダンプが RCM 順なら PERM で戻す
        Tg[np.searchsorted(xs, xr), np.searchsorted(rs, rr)] = Tn[perm]
    dTs_f, Qax_f = solid_indicators(xs, rs, Tg, k_s)
    Qtot = abs(qh)
    # 参照 (中間水準) の固体温度
    xs_f, ys_f, G, W = ev.load(run, st, "wall_3")
    P, xsr, ysr, jw, ks = ev.build("A", xs_f, ys_f, G, 1, pc=pc); Tr = ev._solve_rowwise(P, P.T_in_profile)
    dTs_r, Qax_r = solid_indicators(xsr, ysr[jw:], Tr[:, jw:], ks)
    print(f"=== 固体の効果 {run} ({case}) step {st}")
    print(f"  厚さ方向の温度差 max ΔT_s [K]: forge {dTs_f:.4f} / 参照 {dTs_r:.4f}")
    print(f"  固体の軸方向熱量 max |Q_ax|/Q_tot: forge {np.abs(Qax_f).max()/Qtot:.4f} / 参照 {np.abs(Qax_r).max()/abs(P.robin_heat()):.4f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
