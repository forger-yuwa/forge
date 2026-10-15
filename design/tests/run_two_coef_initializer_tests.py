"""CONTUR × 倍率 c + 一定のずらし δ₀ (plan tooling-nozzle-card-local-mesh §4.3・§5.1 #3n) と 3 係数 (c + c₁·x)·δ_C + δ₀ (#3t) の試験。

(1) 既定 (キー無し・(1, 0) の明示) は `integral_delta_r` の δ_r・導関数・記録が今までとビット同一。
(2) 2 係数: δ_r = c·δ_C + δ₀/scale (r_t 単位) が表・関数の値で成り立ち、導関数は c·δ_C′ (5 次補間の表の誤差の範囲)。
    記録に係数と係数前のスロート・出口の値が残る。寸法の逆算の記録にも係数が出る。
(3) 不正な値 (c ≤ 0・数でない) は ValueError。
(4) case/71 の pass 0 の物理壁の δ (run_0003_ns_cont の delta_in_prev = 格子の NS の壁 − Euler の壁を測った値、手元に記録があるとき) と
    (1) の既定が格子の測定の差 (≤ 0.002 mm) で一致する。
使い方: design/.venv-opt/bin/python tests/run_two_coef_initializer_tests.py
"""
import gzip
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "design"))
from forge_design.evaluate import runner_axismach as RA  # noqa: E402
from forge_design.feedback import deltastar_loop as DL  # noqa: E402

FAIL = []


def check(name, ok, detail=""):
    print(("ok   " if ok else "FAIL ") + name + (f"  ({detail})" if detail else ""))
    if not ok:
        FAIL.append(name)


YAML = ROOT / "case/71.m55_combustion_wt/problem_m55_final.yaml"
p = RA.load_problem(YAML)
d = RA.design_chain(p)
S = float(p.spec["r_throat"])

r0, f0, i0 = RA.integral_delta_r(p, d, {"model": "contur"})
r1, f1, i1 = RA.integral_delta_r(p, d, {"model": "contur", "delta_scale": 1.0, "delta_offset_m": 0.0})
xq = np.linspace(-0.5, float(r0["x"][-1]), 3001)
same = (np.array_equal(r0["delta_r"], r1["delta_r"]) and np.array_equal(f0(xq), f1(xq))
        and all(np.array_equal(f0(xq, k), f1(xq, k)) for k in (1, 2, 3)) and i0 == i1 and "delta_r_contur" not in r1)
check("(1) 既定の明示 (1, 0) はキー無しとビット同一 (表・値・導関数 1〜3・記録)", same)

c, d0 = 0.852073, -0.32201e-3
r2, f2, i2 = RA.integral_delta_r(p, d, {"model": "contur", "delta_scale": c, "delta_offset_m": d0})
x = r2["x"]
check("(2) 係数前の CONTUR は既定の δ_r とビット同一", np.array_equal(r2["delta_r_contur"], r0["delta_r"]))
check("(2) 表: δ_r = c·δ_C + δ₀/scale", np.max(np.abs(r2["delta_r"] - (c * r0["delta_r"] + d0 / S))) <= 1e-15,
      f"最大差 {np.max(np.abs(r2['delta_r'] - (c * r0['delta_r'] + d0 / S))):.2e} r_t")
ev = np.max(np.abs(f2(xq) - (c * f0(xq) + d0 / S)))
check("(2) 関数の値: c·δ_C + δ₀/scale (5 次補間の表の誤差 ≤ 1e−9 r_t)", ev <= 1e-9, f"最大差 {ev:.2e} r_t")
ed = max(np.max(np.abs(f2(xq, k) - c * f0(xq, k))) / max(np.max(np.abs(f0(xq, k))), 1e-30) for k in (1, 2))
check("(2) 導関数 1・2 階は c·δ_C′ (相対 ≤ 1e−6、ずらしは導関数に入らない)", ed <= 1e-6, f"相対の最大差 {ed:.2e}")
tc = i2.get("two_coef") or {}
check("(2) 記録に係数・係数前のスロートと出口の値・scale が残る",
      tc.get("delta_scale") == c and tc.get("delta_offset_m") == d0 and tc.get("scale_m") == S
      and abs(tc.get("delta_r_throat_contur") - i0["delta_r_throat"]) == 0.0 and tc.get("delta_r_exit_contur") == i0["delta_r_exit"]
      and abs(i2["delta_r_throat"] - (c * i0["delta_r_throat"] + d0 / S)) <= 1e-12,
      f"スロート δ_C {i0['delta_r_throat'] * S * 1e3:.5f} mm → {i2['delta_r_throat'] * S * 1e3:.5f} mm、出口 {i0['delta_r_exit'] * S * 1e3:.4f} → {i2['delta_r_exit'] * S * 1e3:.4f} mm")
sup = DL._sizing_delta_supplier(p, {"model": "contur", "delta_scale": c, "delta_offset_m": d0})
_drx, rec = sup(d, S)
sup0 = DL._sizing_delta_supplier(p, {"model": "contur"})
_drx0, rec0 = sup0(d, S)
check("(2) 寸法の逆算の記録に係数が出る (既定では出ない)",
      rec["settings"].get("delta_scale") == c and rec["settings"].get("delta_offset_m") == d0
      and "delta_scale" not in rec0["settings"] and "delta_offset_m" not in rec0["settings"])
check("(2) 寸法の逆算の δ_r 関数 = prepare_ns の経路の δ_r 関数 (同じ r_t)", np.array_equal(_drx(xq), f2(xq)))

for bad in ({"delta_scale": 0.0}, {"delta_scale": -1.0}, {"delta_scale": "x"}, {"delta_offset_m": float("nan")}):
    try:
        RA.integral_delta_r(p, d, {"model": "contur", **bad})
        check(f"(3) 不正な値 {bad} は ValueError", False, "例外が出ない")
    except ValueError:
        check(f"(3) 不正な値 {bad} は ValueError", True)

ref = ROOT / "case/71.m55_combustion_wt/_runs/run_0005_ns/judge/diag_wave/pass0_run_0003_delta_r_next.csv.gz"
if ref.is_file():
    with gzip.open(ref, "rt") as fh:
        a = np.loadtxt(fh, delimiter=",", skiprows=1)
    m = (a[:, 0] >= 0.0)
    e = np.max(np.abs(f0(a[m, 0]) - a[m, 2])) * S * 1e3
    check("(4) 既定の CONTUR と pass 0 の物理壁の δ (run_0003_ns_cont の delta_in_prev、格子から測った値、x ≥ 0) が ≤ 0.002 mm で一致", e <= 2e-3, f"最大差 {e:.2e} mm")
else:
    print("[-- ] (4) pass 0 の記録が無い — 飛ばす")

# (5) 3 係数 (plan §4.3 の末尾・#3t): δ = (c + c₁·x_m)·δ_C + δ₀/scale (x_m = x·scale)
c0_, c1_, d0_ = 0.6743229906561874, 0.02148862611164012, -3.0538768354168355e-04
r3, f3, i3 = RA.integral_delta_r(p, d, {"model": "contur", "delta_scale": c0_, "delta_scale_slope_per_m": c1_, "delta_offset_m": d0_})
want = (c0_ + c1_ * r0["x"] * S) * r0["delta_r"] + d0_ / S
check("(5) 3 係数の表: δ = (c + c₁·x_m)·δ_C + δ₀/scale", np.max(np.abs(r3["delta_r"] - want)) <= 1e-15,
      f"最大差 {np.max(np.abs(r3['delta_r'] - want)):.2e} r_t")
ev3 = np.max(np.abs(f3(xq) - ((c0_ + c1_ * np.asarray(xq) * S) * f0(xq) + d0_ / S)))
check("(5) 3 係数の関数の値 (5 次補間の表の誤差 ≤ 1e−9 r_t)", ev3 <= 1e-9, f"最大差 {ev3:.2e} r_t")
xin = xq[(xq > 0.5) & (xq < float(r0["x"][-1]) - 0.5)]
d1w = c1_ * S * f0(xin) + (c0_ + c1_ * S * xin) * f0(xin, 1)
ed3 = np.max(np.abs(f3(xin, 1) - d1w)) / max(np.max(np.abs(d1w)), 1e-30)
check("(5) 3 係数の導関数 = c₁·scale·δ_C + (c + c₁·x_m)·δ_C′ (相対 ≤ 1e−6)", ed3 <= 1e-6, f"相対の最大差 {ed3:.2e}")
check("(5) 記録に c₁ が残る (2 係数の記録には無い)", i3.get("delta_scale_slope_per_m") == c1_ and (i3.get("two_coef") or {}).get("delta_scale_slope_per_m") == c1_
      and "delta_scale_slope_per_m" not in i2 and "delta_scale_slope_per_m" not in (i2.get("two_coef") or {}))
r4, f4, i4 = RA.integral_delta_r(p, d, {"model": "contur", "delta_scale": c, "delta_offset_m": d0, "delta_scale_slope_per_m": 0.0})
check("(5) c₁ = 0 の明示は 2 係数とビット同一", np.array_equal(r4["delta_r"], r2["delta_r"]) and np.array_equal(f4(xq), f2(xq)) and i4 == i2)
_dx, rec3 = DL._sizing_delta_supplier(p, {"model": "contur", "delta_scale": c0_, "delta_scale_slope_per_m": c1_, "delta_offset_m": d0_})(d, S)
check("(5) 寸法の逆算の記録に c₁ が出て、δ_r 関数は prepare_ns の経路と同じ", rec3["settings"].get("delta_scale_slope_per_m") == c1_ and np.array_equal(_dx(xq), f3(xq)))
for bad in ({"delta_scale_slope_per_m": "x"}, {"delta_scale_slope_per_m": float("inf")}):
    try:
        RA.integral_delta_r(p, d, {"model": "contur", **bad})
        check(f"(5) 不正な c₁ {bad} は ValueError", False, "例外が出ない")
    except ValueError:
        check(f"(5) 不正な c₁ {bad} は ValueError", True)

print("\nALL PASS" if not FAIL else f"\nFAIL {len(FAIL)}: {FAIL}")
sys.exit(1 if FAIL else 0)
