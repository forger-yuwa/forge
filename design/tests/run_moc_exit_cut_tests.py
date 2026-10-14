#!/usr/bin/env python3
"""cplus の閉包の物理出口での壁の表の切り方 `moc_inverse._cut_wall_at_exit` の試験
(plan tooling-nozzle-deltastar-band-root-zone §5.1 #4c、deliverables plan §5.1 #15 の諮問の方向)。

切る規則 `wall[:, 0] <= x_F` は、丸めで x_F が終端の壁点の 1e-13 下になるとその点を落とし、壁の終点が 1 点間隔短くなっていた。
今は切った後、x_F が残った最後の点より 1e-9 r_t 以上下流で、元の壁に x_F を挟む区間があれば、その区間で全列を線形補間した点
(x_F, ...) を末尾に足す。

  (a) 合成の壁: x_F が壁点の 1e-13 下 → 補間した点が足され、終点 = x_F・値はその壁点と丸めの差だけ
  (b) x_F がちょうど壁点の上・1e-13 上 → 足さない (今の切断とビット同一)
  (c) x_F が点の間 → その区間の線形補間の点が足される (r は r_F の定義 np.interp と同じ値)
  (d) しきい値 1e-9 r_t の両側、x_F が壁の範囲の外 (挟む区間なし)
  (e) 実際の逆 MOC (cplus・characteristic・M4/R2・n500): 壁の表の終点が x_F から 1e-9 r_t 以内、exit に記録

実行 (forge-deliv の直下で): /home/sano/work/forge/design/.venv-opt/bin/python design/tests/run_moc_exit_cut_tests.py
所要: 約 10 秒 ((e) の design_chain 1 回)。
"""
import sys
import tempfile
from pathlib import Path

import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from forge_design.evaluate.runner_axismach import design_chain  # noqa: E402
from forge_design.geometry import moc_inverse as MI  # noqa: E402
from forge_design.probdef import load_problem  # noqa: E402

FAIL = 0


def check(name, cond, info=""):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name + (f"  [{info}]" if info else ""))
    if not cond:
        FAIL += 1


def old_cut(w, x_F):
    return w[w[:, 0] <= x_F]


# 合成の壁 [x, r, θ, M]: 単調な曲線 (列ごとに違う非線形の関数)
x = np.linspace(0.0, 10.0, 101)
W = np.c_[x, 1.0 + 0.3 * np.tanh(0.4 * x), 0.2 * np.exp(-0.3 * x), 1.0 + 0.4 * x + 0.01 * x ** 2]
k = 50
xk = float(W[k, 0])
check("合成の壁: 50 番目の点の x がちょうど 5.0", xk == 5.0, f"{xk!r}")

print("== (a) x_F が壁点の 1e-13 下")
x_F = xk - 1e-13
w, info = MI._cut_wall_at_exit(W, x_F)
check("今の切断はその点を落とす (終点が 1 点間隔短い)", old_cut(W, x_F)[-1, 0] == W[k - 1, 0])
check("点が 1 つ足され、終点 = x_F", len(w) == k + 1 and w[-1, 0] == x_F, f"n {len(w)}, 終点 {w[-1, 0]!r}")
check("残った点は元の表とビット同一", w[:k].tobytes() == W[:k].tobytes())
check("足した点の r = np.interp(x_F) (r_F の定義と同じ値)", w[-1, 1] == float(np.interp(x_F, W[:, 0], W[:, 1])))
check("足した点の全列が落ちた壁点と 1e-12 以内", np.max(np.abs(w[-1] - W[k])) <= 1e-12, f"{np.max(np.abs(w[-1] - W[k])):.2e}")
check("exit の記録: 足した・間隔 ≈ 1 点間隔 (0.1)", info["x_F_appended"] is True and abs(info["x_F_minus_last_node_rt"] - 0.1) < 1e-12,
      str(info))

print("== (b) x_F がちょうど壁点の上・1e-13 上")
for d, lab in ((0.0, "ちょうど"), (1e-13, "1e-13 上")):
    x_F = xk + d
    w, info = MI._cut_wall_at_exit(W, x_F)
    check(f"{lab}: 足さない (今の切断とビット同一、終点 = その壁点)",
          w.shape == old_cut(W, x_F).shape and w.tobytes() == old_cut(W, x_F).tobytes() and w[-1, 0] == xk
          and info["x_F_appended"] is False and info["x_F_minus_last_node_rt"] == x_F - xk, str(info))

print("== (c) x_F が点の間")
x_F = 5.03
w, info = MI._cut_wall_at_exit(W, x_F)
t = (x_F - W[k, 0]) / (W[k + 1, 0] - W[k, 0])
lin = W[k] + t * (W[k + 1] - W[k])
check("点が 1 つ足され、終点 = x_F", len(w) == k + 2 and w[-1, 0] == x_F and w[:k + 1].tobytes() == W[:k + 1].tobytes())
check("足した点の全列が区間 [5.0, 5.1] の線形補間 (相対 1e-14)", np.max(np.abs(w[-1] - lin) / np.maximum(1.0, np.abs(lin))) <= 1e-14,
      f"{np.max(np.abs(w[-1] - lin)):.2e}")
check("足した点の r = np.interp(x_F)", w[-1, 1] == float(np.interp(x_F, W[:, 0], W[:, 1])))
check("exit の記録: 足した・間隔 0.03", info["x_F_appended"] is True and abs(info["x_F_minus_last_node_rt"] - 0.03) < 1e-12, str(info))

print("== (d) しきい値 1e-9 r_t の両側、壁の範囲の外")
w, info = MI._cut_wall_at_exit(W, xk + 0.5e-9)
check("間隔 0.5e-9 < 1e-9: 足さない (落ちていない扱い)", len(w) == k + 1 and w[-1, 0] == xk and info["x_F_appended"] is False)
w, info = MI._cut_wall_at_exit(W, xk + 2e-9)
check("間隔 2e-9 ≥ 1e-9: 足す", len(w) == k + 2 and w[-1, 0] == xk + 2e-9 and info["x_F_appended"] is True)
w, info = MI._cut_wall_at_exit(W, 10.05)
check("x_F が壁の終点より下流 (挟む区間なし): 足さない", w.tobytes() == W.tobytes() and info["x_F_appended"] is False, str(info))

print("== (e) 実際の逆 MOC (cplus・characteristic)")
BASE = {"name": "exit_cut", "type": "wind_tunnel_axisym_axismach", "gas": {"gamma": 1.4, "cp": 1004.5},
        "spec": {"Pt": 1.0e6, "Tt": 300.0, "M_design": 4.0, "r_throat": 0.05},
        "dv": {"L_total": {"value": 20.0}},
        "geometry": {"R": 2.0, "start_line": "throat_char", "wall_mode": "cplus", "n_axis_inv": 500, "n_start": 41,
                     "Lc_mode": "from_length"},
        "mesh": {}, "evaluate": {}}
with tempfile.TemporaryDirectory() as td:
    path = Path(td) / "exit_cut.yaml"
    path.write_text(yaml.safe_dump(BASE, sort_keys=False, allow_unicode=True))
    d = design_chain(load_problem(path))
x_F, x_e = float(d["exit"]["x_F"]), float(d["wall_inv"][-1, 0])
ex = d["exit"]
check("壁の表の終点が x_F から 1e-9 r_t 以内 (下流に出ない)", 0.0 <= x_F - x_e < 1e-9, f"x_F {x_F!r}, x_e {x_e!r}")
check("exit に足したか・間隔を記録", isinstance(ex.get("x_F_appended"), bool) and np.isfinite(ex.get("x_F_minus_last_node_rt")), str(
    {k_: ex.get(k_) for k_ in ("x_F_appended", "x_F_minus_last_node_rt")}))
check("設計壁の終点 = 壁の表の終点", float(d["wall"].x_e) == x_e)

print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAIL'}")
sys.exit(1 if FAIL else 0)
