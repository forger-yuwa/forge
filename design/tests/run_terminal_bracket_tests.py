#!/usr/bin/env python3
"""cplus × characteristic の終端 C⁺ の補間 (`moc_inverse._terminal_bracket`・`_cplus_wall_cross`) の試験
(plan tooling-nozzle-card-local-mesh §4.5・§6 V9 (a)・(a′))。

物理出口 F は、x_E を挟む 2 節点から出る網の C⁺ と壁の交点を x_E で線形補間した点。x_E が節点に乗るときはその 1 本の値
(= 以前の「最寄りの節点」の値とビット同一)。

  (a) 降順の軸 [3, 2, 1] で、節点に乗る x_E (最上流の 1 を含む) は 1 本、間は 2 本と重み
  (a′) 範囲の外・非有限・降順でない軸は ValueError
  (b) 実際の逆 MOC (M4/R2・n500、`run_moc_exit_cut_tests.py` (e) と同じ問題): `_design_cplus` の引数を捕まえ、
      x_E を節点に置くと x_F = その節点の C⁺ と壁の交点 (以前の最寄りの値) とビット同一、
      節点の中点に置くと x_F = 2 本の交点の平均 (1e-12 r_t)、範囲の外は ValueError

実行: /home/sano/work/forge/design/.venv-opt/bin/python design/tests/run_terminal_bracket_tests.py
所要: 約 10 秒 (design_chain 1 回と `_design_cplus` 4 回)。
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


def raises(fn, exc=ValueError):
    try:
        fn()
    except exc:
        return True
    except Exception:
        return False
    return False


print("== (a) 降順の軸 [3, 2, 1]")
ax = np.array([3.0, 2.0, 1.0])
check("x_E = 1 (最上流の節点): 1 本", MI._terminal_bracket(ax, 1.0) == [(2, 1.0)], str(MI._terminal_bracket(ax, 1.0)))
check("x_E = 3 (最下流の節点): 1 本", MI._terminal_bracket(ax, 3.0) == [(0, 1.0)])
check("x_E = 2 (中の節点): 1 本", MI._terminal_bracket(ax, 2.0) == [(1, 1.0)])
br = MI._terminal_bracket(ax, 2.25)
check("x_E = 2.25: 下流 0 (重み 0.25)・上流 1 (重み 0.75)", br == [(0, 0.25), (1, 0.75)], str(br))
br = MI._terminal_bracket(ax, 1.0 + 1e-12)
check("x_E が最上流の節点の 1e-12 下流: 下流 1・上流 2", [m for m, _ in br] == [1, 2] and 0.0 < br[0][1] < 1e-11, str(br))

print("== (a′) 異常な入力")
check("範囲の外 (下流): ValueError", raises(lambda: MI._terminal_bracket(ax, 3.0 + 1e-9)))
check("範囲の外 (上流): ValueError", raises(lambda: MI._terminal_bracket(ax, 1.0 - 1e-9)))
check("NaN: ValueError", raises(lambda: MI._terminal_bracket(ax, float("nan"))))
check("inf: ValueError", raises(lambda: MI._terminal_bracket(ax, float("inf"))))
check("昇順の軸: ValueError", raises(lambda: MI._terminal_bracket(ax[::-1].copy(), 2.0)))
check("同じ値の節点を含む軸: ValueError", raises(lambda: MI._terminal_bracket(np.array([3.0, 2.0, 2.0, 1.0]), 2.5)))

print("== (b) 実際の逆 MOC (cplus・characteristic)")
BASE = {"name": "terminal_bracket", "type": "wind_tunnel_axisym_axismach", "gas": {"gamma": 1.4, "cp": 1004.5},
        "spec": {"Pt": 1.0e6, "Tt": 300.0, "M_design": 4.0, "r_throat": 0.05},
        "dv": {"L_total": {"value": 20.0}},
        "geometry": {"R": 2.0, "start_line": "throat_char", "wall_mode": "cplus", "n_axis_inv": 500, "n_start": 41,
                     "Lc_mode": "from_length"},
        "mesh": {}, "evaluate": {}}
captured = {}
orig = MI._design_cplus


def capture(*a, **k):
    captured["a"], captured["k"] = a, k
    return orig(*a, **k)


MI._design_cplus = capture
try:
    with tempfile.TemporaryDirectory() as td:
        path = Path(td) / "terminal_bracket.yaml"
        path.write_text(yaml.safe_dump(BASE, sort_keys=False, allow_unicode=True))
        d = design_chain(load_problem(path))
finally:
    MI._design_cplus = orig
ex = d["exit"]
check("設計の exit に挟む節点・重み・2 本の交点を記録", len(ex.get("x_E_bracket", [])) in (1, 2) and "x_E_weight" in ex and "x_F_nodes" in ex,
      str({k_: ex.get(k_) for k_ in ("x_E_bracket", "x_E_weight", "x_F_nodes")}))
a, k = list(captured["a"]), dict(captured["k"])
# 引数の並び: (inv, init, n_ax, ax, mdot_star, g, x0, x_wall0, exit_mode, x_E, M_d, start_line)
axv = np.asarray(a[3], dtype=float)
m_hi = int(ex["x_E_bracket"][0])


def run_at(x_E):
    b = list(a)
    b[9] = x_E
    return orig(*b, **k)


r = run_at(float(axv[m_hi]))
lev, wall_full = r["levels"], r["wall_full"]
x_old = MI._cplus_wall_cross(lev, wall_full, m_hi)
check("x_E を節点に置く: 1 本・x_F = その節点の C⁺ の交点 (以前の最寄りの値) とビット同一",
      r["exit"]["x_E_bracket"] == [m_hi] and r["exit"]["x_F"] == x_old and r["exit"]["x_E_index"] == m_hi,
      f"x_F {r['exit']['x_F']!r}, 最寄り {x_old!r}")
xm = 0.5 * (float(axv[m_hi]) + float(axv[m_hi + 1]))
r = run_at(xm)
x_hi, x_lo = MI._cplus_wall_cross(lev, wall_full, m_hi), MI._cplus_wall_cross(lev, wall_full, m_hi + 1)
check("x_E を節点の中点に置く: x_F = 2 本の交点の平均 (1e-12 r_t)",
      abs(r["exit"]["x_F"] - 0.5 * (x_hi + x_lo)) <= 1e-12 and abs(r["exit"]["x_E_weight"] - 0.5) <= 1e-9,
      f"x_F {r['exit']['x_F']!r}, 平均 {0.5 * (x_hi + x_lo)!r}, 重み {r['exit']['x_E_weight']!r}")
check("2 本の交点の差は網の間隔の程度 (0 < 差 < 2 × 間隔)",
      0.0 < x_hi - x_lo < 2.0 * (float(axv[m_hi]) - float(axv[m_hi + 1])), f"{x_hi - x_lo:.3e}")
check("x_E が最下流の節点より下流: ValueError", raises(lambda: run_at(float(axv[0]) + 1.0)))

print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAIL'}")
sys.exit(1 if FAIL else 0)
