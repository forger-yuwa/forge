#!/usr/bin/env python3
"""③ベル (M) モードの TP の Rao の生成器 `geometry/rao_axisym.py` と等エントロピーの accessor `gas/isentrope.py` の単体試験
(plan tooling-nozzle-bell-two-mode §5.1 #5・§6 V1 の T1–T6、導出は notes/reviews/2026-10-11-bell-s1-rao-tp-derivation.md)。

**この版で回すもの** (式と線のレベル、前進 MOC を使わないか短い):

- accessor: TP の経路を定比熱の擬似種で回すと CPG の閉形式と 1e-12 で一致 (状態・ν・M(ν)・質量流束密度)、基準ケース EXH で
  d ln V/dM = 1/[M(1−k)]・dν/dM = √(M²−1)/[M(1−k)]・d ln ρ/dM = −M/(1−k) を差分で 1e-8、1D の C_F が排気ファイルの
  1D 理想値 (F_1D,vac 1159.1 kN @ p_c 10 MPa・r_t 0.14 m → C_F 1.8824) と 1e-4。
- **T1b** (平面縮約): 縁の関数 `edge_function` の根が `rao_planar.theta_lip` と 1e-8 (γ 1.2・1.4、M_e 2〜4.5、p_e/p_a 1.5〜10)。
  加えて DE (Φ₂・Φ₃ 一定で引いた線) が軸対称の C⁺ の適合式 d(θ−ν) = −sinμ sinθ/r dσ を CPG と EXH で満たす
  (残差 ≤ 1e-12 × 尺度、源項の符号を反転した対照は ≥ 1e-2 × 尺度)。
- **T2** (包絡線の恒等式、前進 MOC なし): 適合式を満たす合成の C⁻ データ (θ を与え ν を C⁻ の適合式で積分) の上で D を動かし、
  dF/dη_D = 2π[Φ₃ dL/dη_D + F(E) dr_E/dη_D] (F = BD の下端 → D + DE の推力、F(E) = `edge_function`) が相対 1e-5。
  DE の M 格子 8e-3・4e-3・2e-3 のどれでも成り立ち、格子を倍にして悪化しない。CPG と EXH、p_a 0 と p_a > 0。
- **T1d** (定比熱の擬似種で CPG 生成器と 1e-8): 同じ D からの DE の E と推力、同じ壁・初期線の kernel 場 (前進 MOC の網)。

**この版で回さないもの (方針の穴、主セッションに返した)**: T1a・T1c・T3・T4・T5・T6 は `design_rao` を基準ケースの L_k で回す
試験だが、記録 (構成手順 4、plan §4.4) の「DE で |c| > 1 は解なし」の規則の下では θ_B > μ_B の kernel で D が θ_D = μ_D の点より
B 寄りになる解 (L ≲ 0.9 L_cone) がすべて INFEASIBLE になる。|c| = 1 は DE の M のパラメータ化の折り返し (θ = μ で M が極値) で、
r をパラメータにすれば DE は滑らかに続く (Jacobian は正則)。規則を変えるかは上位の判断。
`--blocked` を付けると、その事実 (規則の下の L_k 5 本の可否) を記録として出力する (数分)。

実行: design/.venv-opt の Python で `python design/tests/run_rao_axisym_tests.py`。
"""
import sys
import time
import warnings
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "design"))
from forge_design.gas.isentrope import Isentrope, constant_cp_gas  # noqa: E402
from forge_design.gas.rocket_exhaust import RocketExhaust  # noqa: E402
from forge_design.geometry import rao_axisym as ra  # noqa: E402
from forge_design.geometry.rao_planar import theta_lip  # noqa: E402

EXHAUST = ROOT / "case" / "68.bell_two_mode" / "exhaust_lox_lh2_pc100_of6.yaml"
FAIL = 0


def check(name, cond, detail=""):
    global FAIL
    print(("ok   " if cond else "FAIL ") + name + (f"  [{detail}]" if detail else ""), flush=True)
    if not cond:
        FAIL += 1


EXH = RocketExhaust.from_file(EXHAUST)
IG_EXH = ra.as_isentrope(EXH.moc_gas())

# ================================================================ accessor
for g in (1.2, 1.3):
    ic = Isentrope(g)
    it = Isentrope(constant_cp_gas(g, Tt=3000.0))
    M = np.linspace(1.0, 7.5, 1301)
    sc, stt = ic.state(M), it.state(M)
    worst = max(float(np.max(np.abs(sc[k] - stt[k]) / np.maximum(np.abs(sc[k]), 1e-300)))
                for k in ("T", "V", "a", "rho", "p", "k", "mu"))
    check(f"accessor: 定比熱の擬似種の TP = CPG 閉形式 (γ {g}、状態 7 量)", worst <= 1e-12, f"最大相対差 {worst:.2e}")
    nu = ic.nu(M)
    d_nu = float(np.max(np.abs(it.nu(M) - nu)))
    d_inv = float(np.max(np.abs(it.mach_of_nu(nu) - M)))
    d_G = float(np.max(np.abs(it.mass_flux_density(M) / ic.mass_flux_density(M) - 1.0)))
    check(f"accessor: ν・M(ν)・ρV/(ρ*a*) の TP = CPG (γ {g})", max(d_nu, d_inv, d_G) <= 1e-12,
          f"ν {d_nu:.1e}・M(ν) {d_inv:.1e}・質量流束 {d_G:.1e}")
M = np.array([1.3, 2.0, 3.0, 4.5, 6.0])
h = 1e-5
s0, sp, sm = IG_EXH.state(M), IG_EXH.state(M + h), IG_EXH.state(M - h)
e1 = np.max(np.abs((np.log(sp["V"]) - np.log(sm["V"])) / (2 * h) / s0["dlnV_dM"] - 1))
e2 = np.max(np.abs((IG_EXH.nu(M + h) - IG_EXH.nu(M - h)) / (2 * h) / (np.sqrt(M * M - 1) * s0["dlnV_dM"]) - 1))
e3 = np.max(np.abs((np.log(sp["rho"]) - np.log(sm["rho"])) / (2 * h) / (-M / (1 - s0["k"])) - 1))
e4 = np.max(np.abs((np.log(sp["a"]) - np.log(sm["a"])) / (np.log(sp["V"]) - np.log(sm["V"])) / s0["k"] - 1))
check("accessor: EXH の関係式 (d ln V/dM・dν/dM・d ln ρ/dM・k = d ln a/d ln V) を差分で 1e-8",
      max(e1, e2, e3, e4) <= 1e-8, f"{e1:.1e} / {e2:.1e} / {e3:.1e} / {e4:.1e}")
F1 = IG_EXH.ideal_1d(50.0)
CF_file = 1159.1e3 / (EXH.Pc * np.pi * 0.14 ** 2)
check("accessor: 1D の C_F (ε 50) が F_1D,vac 1159.1 kN と 1e-4 (plan §5.1 #2a)", abs(F1["CF"] / CF_file - 1) <= 1e-4,
      f"C_F {F1['CF']:.6f} / {CF_file:.6f}、M_e {F1['M_e']:.4f}")

# ================================================================ T1b
worst = 0.0
for g in (1.2, 1.4):
    ic = Isentrope(g)
    for Me in (2.0, 3.0, 4.5):
        st = ic.state(np.array([Me]))
        for pe_pa in (1.5, 3.0, 10.0):
            pa = float(st["p"][0]) / pe_pa
            th = brentq(lambda t: ra.edge_function(ic, Me, t, 1.0, pa), 1e-9, 0.25 * np.pi - 1e-9, xtol=1e-15)
            worst = max(worst, abs(th - theta_lip(Me, pe_pa, g)))
check("T1b: 縁条件 edge_function = 0 の根が rao_planar.theta_lip と 1e-8 (γ 1.2・1.4、M_e 2・3・4.5、p_e/p_a 1.5・3・10)",
      worst <= 1e-8, f"最大差 {worst:.2e} rad")


def synth_bd(ig, xB=0.7, rB=1.25, thB=np.deg2rad(30.0), MB=1.9, sig_max=3.0, N=1500):
    """合成の C⁻ データ: θ(σ) を与え、ν を C⁻ の適合式 d(θ+ν)/dσ = sinμ sinθ/r で積分 (DOP853、rtol 1e-13)。
    η = 0 が下端、η = N が B。"""
    nuB = float(ig.nu(np.array([MB]))[0])

    def thf(s):
        return thB * (1 - 0.8 * s / sig_max) - 0.05 * np.sin(2 * s)

    def dthf(s):
        return -thB * 0.8 / sig_max - 0.1 * np.cos(2 * s)

    def rhs(s, y):
        x, r, nu = y
        Mm = float(ig.mach_of_nu(np.array([nu]))[0])
        mu = np.arcsin(1 / Mm)
        th = thf(s)
        return [np.cos(th - mu), np.sin(th - mu), np.sin(mu) * np.sin(th) / r - dthf(s)]
    sol = solve_ivp(rhs, (0, sig_max), [xB, rB, nuB], method="DOP853", rtol=1e-13, atol=1e-14, dense_output=True)
    ss = np.linspace(sig_max, 0, N + 1)
    y = sol.sol(ss)
    return ra.BDLine(np.column_stack([y[0], y[1], thf(ss), y[2]]), ig)


for name, ig in (("CPG γ 1.2", ra.as_isentrope(1.2)), ("EXH", IG_EXH)):
    bd = synth_bd(ig)
    d = ra._de_at(bd, 0.4 * bd.N, 0.0, 4e-3)
    rl = d["rl"]
    Ms = np.linspace(d["D"]["M"], d["E"]["M"], 200)
    f = rl.eval(Ms)
    st = ig.state(Ms)
    lhs = f["thp"] - np.sqrt(Ms * Ms - 1) * st["dlnV_dM"]
    src = -np.sin(f["mu"]) * np.sin(f["th"]) / f["r"] * f["rp"] / np.sin(f["th"] + f["mu"])
    scale = float(np.max(np.abs(lhs)))
    res = float(np.max(np.abs(lhs - src)))
    flip = float(np.max(np.abs(lhs + src)))
    check(f"T1b: DE (Φ₂・Φ₃ 一定) が軸対称の C⁺ の適合式を満たす ({name})", res <= 1e-12 * scale and flip >= 1e-2 * scale,
          f"残差 {res:.1e}、符号反転の対照 {flip:.2e}、尺度 {scale:.2f}")

# ================================================================ T2
for name, ig in (("CPG γ 1.2", ra.as_isentrope(1.2)), ("EXH", IG_EXH)):
    bd = synth_bd(ig)
    for p_a in (0.0, 2e-3):
        rels = {}
        for dM in (8e-3, 4e-3, 2e-3):
            rr = []
            for eta in (0.3 * bd.N, 0.5 * bd.N, 0.7 * bd.N):
                hh = 0.5
                v = {k: ra._de_at(bd, eta + k * hh, p_a, dM) for k in (-2, -1, 0, 1, 2)}
                D5 = lambda q: (-q[2] + 8 * q[1] - 8 * q[-1] + q[-2]) / (12 * hh)
                lhs = D5({k: w["F_part"] for k, w in v.items()})
                rhs = 2 * np.pi * (v[0]["rl"].Phi3 * D5({k: w["E"]["x"] for k, w in v.items()})
                                   + v[0]["edge"] * D5({k: w["E"]["r"] for k, w in v.items()}))
                rr.append(abs(lhs - rhs) / abs(lhs))
            rels[dM] = max(rr)
        ok = all(v <= 1e-5 for v in rels.values()) and rels[2e-3] <= max(rels[8e-3], 1e-9) * 10
        check(f"T2: 包絡線の恒等式 dF/dη = 2π[Φ₃ dL/dη + F(E) dr_E/dη] ({name}、p_a {p_a:g})", ok,
              "相対差 " + "・".join(f"dM {k:g}: {v:.1e}" for k, v in rels.items()))

# ================================================================ T1d (線と場のレベル)
g = 1.2
icpg = ra.as_isentrope(g)
itp = ra.as_isentrope(constant_cp_gas(g, Tt=3000.0))
bd_c, bd_t = synth_bd(icpg), synth_bd(itp)
dc = ra._de_at(bd_c, 0.45 * bd_c.N, 0.0, 4e-3)
dt = ra._de_at(bd_t, 0.45 * bd_t.N, 0.0, 4e-3)
worst = max(abs(dc["E"][k] - dt["E"][k]) / max(abs(dc["E"][k]), 1.0) for k in ("x", "r", "th", "M", "T", "mdot"))
check("T1d: 定比熱の擬似種の TP と CPG で、同じ合成 C⁻ からの DE の E と推力が 1e-8", worst <= 1e-8, f"最大相対差 {worst:.1e}")
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    line_c = ra._hall_line(1.0, icpg, 81)
    line_t = ra._hall_line(1.0, itp, 81)
kc = ra.kernel_field(icpg, 1.0, 0.0, line_c, 0.55, L_march=40.0)
kt = ra.kernel_field(itp, 1.0, 0.0, line_t, 0.55, L_march=40.0)
K = min(kc.K, kt.K)
m = np.isfinite(kc.A[:K + 1]) & np.isfinite(kt.A[:K + 1])
dnet = float(np.max(np.abs(kc.A[:K + 1][m] - kt.A[:K + 1][m])))
check("T1d: 定比熱の擬似種の TP と CPG で、同じ円弧 kernel・Hall 初期線 (n 81) の前進 MOC の網が 1e-8", dnet <= 1e-8,
      f"網の (x, r, θ, ν) の最大差 {dnet:.1e} (前線 {K + 1} 本)")

# ================================================================ 回さないもの
BLOCKED = {
    "T1a": "CPG・円弧 kernel・R_t 0.382 のチャート照合 (80 % ベル) — 解が DE の |c| = 1 を越える側にあり、規則の下では根が無い",
    "T1c": "自由 ε と固定 ε の等価 — 同上 (L_k の範囲)",
    "T3": "保存 (閉合と推力の 2 会計) — design_rao が基準ケースの L_k で成立しない",
    "T4": "停留性 — 同上",
    "T5": "C2 (|g̃| ≤ 1e-8) — 同上",
    "T6": "分解能 (n 161/321/641) — 同上。加えて Hall 初期線は R_t ≤ 1 で n ≤ 321 (moc_forward.HALL_N_MAX_RT_LE1)",
}
for k, v in BLOCKED.items():
    print(f"未実施 {k}: {v}", flush=True)

if "--blocked" in sys.argv:
    from forge_design.geometry.wall_bell_c2 import cone15_length
    Lc = cone15_length(50.0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        line = ra._hall_line(1.0, IG_EXH, 161)
    t0 = time.time()
    net = ra.kernel_field(IG_EXH, 1.0, 0.0, line, 0.75, L_march=2.5 * Lc)
    for fr in (0.6, 0.7, 0.8, 0.9, 1.0):
        mid = ra.solve_thetaB(net, fr * Lc, eps=50.0, s_range=(0.40, min(0.75, net.s_need)), n_scan=13, n_scan_D=30)
        sol = [(round(np.rad2deg(r["theta_B"]), 3), round(np.rad2deg(r["D_root"]["E"]["th"]), 3)) for r in mid["roots"]]
        print(f"記録 (円弧 kernel・EXH・n 161・規則 |c| ≤ 1): L {fr:.0%} ({fr * Lc:.3f} r_t) → (θ_B, θ_E) {sol or '根なし'}"
              + ("" if sol else f" — {mid['why'][:120]}"), flush=True)
    print(f"  ({time.time() - t0:.0f} s)")

print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAIL'} (未実施 {len(BLOCKED)} 件: 方針の判断待ち)")
sys.exit(1 if FAIL else 0)
