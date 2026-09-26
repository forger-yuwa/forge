#!/usr/bin/env python3
r"""V6′ の $K$/$D_f$ 感度を**共通帯・節点ごと**に取る (plan §5.1 #100 ⑤ / #101 ④)。

各 run の解析解からの誤差を並べても run 間感度にはならない (帯も run ごとに違う)。ここでは
基準 run (`Df_scale` 5 / `interval` 50) と比較 run の**同じ $y$ の界面節点**で最終壁ダンプを突き合わせ、

    温度差  max |T_w1(run) - T_w1(ref)| / (固体の温度上昇 94.18 K)   <= 0.5 %
    熱流束差 max |q_w1(run) - q_w1(ref)| / q*                         <= 0.5 %

を出す。**共通帯** = 各 run の `v6p_band.json` (eval_v6p.py が書く) の交わり。

使い方:
  python3 case/58.conjugate_slot/sens_v6p.py <ref_run> <run> [<run> ...] [--step 100000]
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import h5py
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from eval_v6p import analytic, last_step  # noqa: E402


def wall(run, step):
    with h5py.File(f"{run}/res_slot_front_5_{step}.h5", "r") as f:
        y = np.asarray(f["MESH/COORD"][:], float).reshape(-1, 3)[:, 1]
        return y, np.asarray(f["VALUE/Ts"][:], float), np.abs(np.asarray(f["VALUE/iface_q_eff"][:], float))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("ref")
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--step", type=int, default=None)
    ap.add_argument("--tol", type=float, default=0.5, help="[%%] (登録値)")
    ap.add_argument("--fix-band", type=float, nargs=2, metavar=("YTOP", "YBOT"), default=None,
                    help="共通帯を外から固定する [m] (省略時は各 run の v6p_band.json の交わり)")
    a = ap.parse_args()
    qs, _, drop = analytic()
    allr = [a.ref] + a.runs
    bands = []
    for r in allr:
        fn = f"{r.rstrip('/')}/v6p_band.json"
        if not os.path.exists(fn):
            raise SystemExit(f"{fn} が無い (先に eval_v6p.py を回すこと)")
        bands.append(json.load(open(fn)))
    ytop = min(b["ytop"] for b in bands); ybot = max(b["ybot"] for b in bands)
    if a.fix_band is not None:
        ytop, ybot = max(a.fix_band), min(a.fix_band)
    # 各 run は**それぞれの最終 step** で比べる (延長 run は基準より長い。--step は全 run 共通に固定するとき)
    step = a.step if a.step else last_step(a.ref.rstrip("/"))
    y0, T0, q0 = wall(a.ref.rstrip("/"), step)
    sel = np.where((y0 <= ytop + 1e-9) & (y0 >= ybot - 1e-9))[0]
    print(f"共通帯 y {ytop*1e3:.3f} .. {ybot*1e3:.3f} mm ({len(sel)} 節点), 基準 {a.ref}, step {step}")
    print(f"{'run':<40}{'温度差 max [%]':>16}{'平均':>10}{'q 差 max [%]':>15}{'平均':>10}  判定")
    worst = 0
    for r in a.runs:
        y, T, q = wall(r.rstrip("/"), a.step if a.step else last_step(r.rstrip("/")))
        j = np.array([np.argmin(np.abs(y - y0[i])) for i in sel])
        if np.abs(y[j] - y0[sel]).max() > 1e-9:
            raise SystemExit(f"{r}: 壁節点が基準と座標で対応しない")
        dT = np.abs(T[j] - T0[sel]) / drop * 100
        dq = np.abs(q[j] - q0[sel]) / qs * 100
        ok = dT.max() <= a.tol and dq.max() <= a.tol
        worst |= not ok
        print(f"{os.path.basename(r.rstrip('/')):<40}{dT.max():16.4f}{dT.mean():10.4f}"
              f"{dq.max():15.4f}{dq.mean():10.4f}  {'PASS' if ok else 'FAIL'}")
    print(f"\nVERDICT: {'PASS' if not worst else 'FAIL'} (最終スナップショットの差。局所の静定は別に判定する)")
    return int(worst)


if __name__ == "__main__":
    raise SystemExit(main())
