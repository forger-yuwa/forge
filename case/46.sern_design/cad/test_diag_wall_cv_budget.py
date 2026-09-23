#!/usr/bin/env python3
"""`diag_wall_cv_budget.py` の `slau_mdot` が実装と一致することの単体試験。

**なぜ要るか** (codex result 段 M1, 2026-09-23): ツールが変更前の chi (速度の大きさ基準) のまま
`--equilibrium` を回していたため、**検証したい補充項そのものが診断から抜けていた**。
「予測と一致して PASS」という報告が無効になった。同じ事故を防ぐため、
**カーネルの式 (`convectiveFlux_slau_d.inc.cuh:538-579`) と項ごとに合うこと**を試験に固定する。

  python3 test_diag_wall_cv_budget.py
"""
import sys, math
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import numpy as np
from diag_wall_cv_budget import slau_mdot

FAIL = []
def check(name, got, want, tol=1e-9):
    ok = abs(got - want) <= tol * max(1.0, abs(want))
    print(f"{'ok  ' if ok else 'FAIL'} {name}: got {got:.6g}, want {want:.6g}")
    if not ok: FAIL.append(name)

def kernel_mdot(A, n, L, R, wall_normal_chi):
    """`convectiveFlux_slau_d.inc.cuh` の mdot を独立に書き下したもの (ツールとは別実装)"""
    VnL = L["Ux"]*n[0] + L["Uy"]*n[1] + L["Uz"]*n[2]
    VnR = R["Ux"]*n[0] + R["Uy"]*n[1] + R["Uz"]*n[2]
    c = 0.5*(L["sonic"] + R["sonic"])
    Mp, Mm = VnL/c, VnR/c
    g = -max(min(Mp, 0.0), -1.0) * min(max(Mm, 0.0), 1.0)
    vh = (L["ro"]*abs(VnL) + R["ro"]*abs(VnR)) / (L["ro"] + R["ro"])
    vhp, vhm = (1-g)*vh + g*abs(VnL), (1-g)*vh + g*abs(VnR)
    u2L = L["Ux"]**2 + L["Uy"]**2 + L["Uz"]**2
    u2R = R["Ux"]**2 + R["Uy"]**2 + R["Uz"]**2
    if wall_normal_chi:
        Mh = min(1.0, math.sqrt(0.5*(VnL**2 + VnR**2))/c)
    else:
        Mh = min(1.0, math.sqrt(0.5*(u2L + u2R))/c)
    chi = (1-Mh)**2
    return A*0.5*((L["ro"]*(VnL+vhp) + R["ro"]*(VnR-vhm)) - chi/c*(R["P"]-L["P"]))

def st(ro, Ux, Uy, Uz, P, son): return {"ro":ro,"Ux":Ux,"Uy":Uy,"Uz":Uz,"P":P,"sonic":son}

print("=== codex result M1 の反例 (壁 0 / 内点は接線のみ 1000 m/s / c 700 / Δp 900 Pa / A 1 m²) ===")
# 面法線 x、内点の速度は y 方向 (接線) のみ -> 面法線速度は両側 0
L = st(1.0e-3, 0.0, 0.0, 0.0, 100.0, 700.0)        # 壁ノード
R = st(1.0e-2, 0.0, 1000.0, 0.0, 1000.0, 700.0)    # 内点 (接線のみ)
n = (1.0, 0.0, 0.0); A = 1.0
m_off, *_ = slau_mdot(A, *n, L, R, False)
m_on,  *_ = slau_mdot(A, *n, L, R, True)
check("flag off: chi=0 で補充項が消える (mdot=0)", m_off, 0.0)
check("flag on : 面法線 Mach=0 -> chi=1 で壁へ流入", m_on, -0.642857142857, 1e-6)
check("flag off がカーネル (旧) と一致", m_off, kernel_mdot(A, n, L, R, False))
check("flag on  がカーネル (新) と一致", m_on,  kernel_mdot(A, n, L, R, True))

print("\n=== 面向き反転で mdot が符号反転のみ (V0) ===")
for wf in (False, True):
    a_, *_ = slau_mdot(A, *n, L, R, wf)
    b_, *_ = slau_mdot(A, -n[0], -n[1], -n[2], R, L, wf)
    check(f"反転 (wall_face={wf})", b_, -a_, 1e-9)

print("\n=== 等状態では移流のみ・圧力差項ゼロ (V0) ===")
S = st(2.0e-2, 300.0, 50.0, 0.0, 5000.0, 600.0)
for wf in (False, True):
    m, adv, pre, *_ = slau_mdot(A, *n, S, S, wf)
    check(f"等状態の圧力差項 (wall_face={wf})", pre, 0.0)
    check(f"等状態の mdot = A*ro*Vn (wall_face={wf})", m, A*S["ro"]*S["Ux"], 1e-9)

print("\n=== flag の差は Δp に比例する (M5: 剥離縁だけに効くのではない) ===")
# 壁隣接面で接線速度差があれば、付着境界層でも chi_n != chi になる。
# 差の大きさを決めるのは **Δp** であって剥離の有無ではない。
prev = None
for dp in (1.0, 10.0, 100.0):
    L2 = st(5.0e-2, 0.0, 0.0, 0.0, 10000.0, 500.0)
    R2 = st(5.1e-2, 0.0, 400.0, 0.0, 10000.0 + dp, 500.0)   # 接線流のみ
    d_off, *_ = slau_mdot(A, *n, L2, R2, False)
    d_on,  *_ = slau_mdot(A, *n, L2, R2, True)
    d = abs(d_on - d_off)
    print(f"     Δp={dp:6.1f} Pa: mdot off {d_off:.6e} / on {d_on:.6e} / 差 {d:.4e}")
    if prev is not None:
        check(f"差が Δp に比例 (Δp {prev[0]:g} -> {dp:g})", d / prev[1], dp / prev[0], 1e-6)
    prev = (dp, d)
# 同じ接線速度でも Δp=0 なら差はゼロ
L3 = st(5.0e-2, 0.0, 0.0, 0.0, 10000.0, 500.0)
R3 = st(5.1e-2, 0.0, 400.0, 0.0, 10000.0, 500.0)
e_off, *_ = slau_mdot(A, *n, L3, R3, False)
e_on,  *_ = slau_mdot(A, *n, L3, R3, True)
check("Δp=0 なら flag の差はゼロ (接線速度があっても)", e_on - e_off, 0.0, 1e-12)

print()
if FAIL:
    print("FAILED:", ", ".join(FAIL)); raise SystemExit(1)
print("ALL PASS")
