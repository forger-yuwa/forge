#!/usr/bin/env python3
"""③ベル (M) モードの TP の Rao の生成器 `geometry/rao_axisym.py` と等エントロピーの accessor `gas/isentrope.py` の単体試験
(plan tooling-nozzle-bell-two-mode §5.1 #5・§6 V1 の T1–T6、導出は notes/reviews/2026-10-11-bell-s1-rao-tp-derivation.md)。

**この版で回すもの** (式と線のレベル、前進 MOC を使わないか短い):

- accessor: TP の経路を定比熱の擬似種で回すと CPG の閉形式と 1e-12 で一致 (状態・ν・M(ν)・質量流束密度)、基準ケース EXH で
  d ln V/dM = 1/[M(1−k)]・dν/dM = √(M²−1)/[M(1−k)]・d ln ρ/dM = −M/(1−k) を差分で 1e-8、1D の C_F が排気ファイルの
  1D 理想値 (F_1D,vac 1159.1 kN @ p_c 10 MPa・r_t 0.14 m → C_F 1.8824) と 1e-4。
- **T1b** (平面縮約): 縁の関数 `edge_function` の根が `rao_planar.theta_lip` と 1e-8 (γ 1.2・1.4、M_e 2〜4.5、p_e/p_a 1.5〜10)。
  加えて DE (Φ₂・Φ₃ 一定で r をパラメータに引いた線) が軸対称の C⁺ の適合式 d(θ−ν) = −sinμ sinθ/r dσ を CPG と EXH で満たす
  (残差 ≤ 1e-12 × 尺度、源項の符号を反転した対照は ≥ 1e-2 × 尺度)。**折り返し (θ = μ) を跨ぐ D** (`synth_bd` の M_B 2.2・
  θ_B 32°、D を B 寄り) を必須にし (plan §6 V1、diagnostician 2026-10-11)、線が実際に θ − μ の符号を変えること、折り返し点で
  dθ/dr = −tanθ/(2r) (C⁺ の適合式と Φ₃ 一定が同じ値を与える点検値) と dM/dr = 0 を相対 1e-9 で確かめる。
- **T2** (包絡線の恒等式、前進 MOC なし): 適合式を満たす合成の C⁻ データ (θ を与え ν を C⁻ の適合式で積分) の上で D を動かし、
  dF/dη_D = 2π[Φ₃ dL/dη_D + F(E) dr_E/dη_D] (F = BD の下端 → D + DE の推力、F(E) = `edge_function`) が相対 1e-5。
  DE の r 格子 `dr` 0.04・0.02・0.01 のどれでも成り立ち、格子を倍にして悪化しない。CPG と EXH、p_a 0 と p_a > 0。
  **跨ぐ D の組** (上と同じ合成、η_D 0.72〜0.78 N — 5 点の差分の全点で θ_D > μ_D・θ_E < μ_E) も同じ条件で必須。
- **(b) 枝の追跡** (`solve_thetaB`、plan §4.4 の訂正): `_link_branches` が η_D の連続性で結ぶこと (添字がずれる場合) と、
  det_D の符号 (DE の族) の違う根を結ばないこと。
- **(e) 受理規則と kernel の既定** (plan §4.4「主枝の受理規則」「(M) は円弧 kernel・B で C1」、diagnostician 2026-10-11):
  `_accept_candidate` は det_D > 0 だけで受理し、θ_E ≤ θ_B と de < 0 は記録のフラグ、det_D > 0 かつ θ_E > θ_B は受理して
  `policy_hole` を立てる (4 象限)。`design_rao` の既定は円弧 (`kernel="arc"`)、`linear` は c の明示が要る、`solve_kappaB` は
  c を明示パラメータとして受け取るだけ (外側の走査・brentq の引数が無い)。
- **T1d** (定比熱の擬似種で CPG 生成器と 1e-8): 同じ D からの DE の E と推力、同じ壁・初期線の kernel 場 (前進 MOC の網)。

**この版で回さないもの**: T3・T5・T6 は `design_rao` を基準ケースの L_k で回す run で、ここでは回さない —
`case/68.bell_two_mode/rao_m/run_rao_t6.py` (円弧 kernel・L100・L70 × n 161/241/321 + n 321 の dr 半分・m・q 倍) が
表 `rao_t6.json`・`rao_t6.csv` を出し、合否は親と上位が判定する (plan §6 V1 の T3・T5・T6、2026-10-11 改)。
T5 は C2 を課さない記録 (g̃ と Δκ_B)。T1c・T4・T1a′ は T6 の後 (§5.1 #5 (g))。T1a は記録 (合否なし) に格下げ (#6 で CFD ピンの後)。

実行: design/.venv-opt の Python で `python design/tests/run_rao_axisym_tests.py`。
"""
import sys
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


def de_compat(d, n_pts=400):
    """DE 上の r 等間隔の点で C⁺ の適合式の残差 (dθ/dr − ν' dM/dr と −sinμ sinθ/(r sin(θ+μ)))。戻り (残差, 対照, 尺度, 状態)。"""
    rl = d["rl"]
    rs = np.linspace(d["D"]["r"], d["E"]["r"], n_pts)
    f = rl.state_at(rs)
    lhs = f["dth_dr"] - np.sqrt(f["M"] ** 2 - 1) * f["dlnV"] * f["dM_dr"]
    src = -np.sin(f["mu"]) * np.sin(f["th"]) / (f["r"] * np.sin(f["th"] + f["mu"]))
    scale = float(np.max(np.abs(lhs)))
    return float(np.max(np.abs(lhs - src))), float(np.max(np.abs(lhs + src))), scale, f


CROSS = dict(MB=2.2, thB=np.deg2rad(32.0))          # θ_B > μ_B (μ(2.2) = 27.0°) — D を B 寄りに取ると DE が θ = μ を跨ぐ
for name, ig in (("CPG γ 1.2", ra.as_isentrope(1.2)), ("EXH", IG_EXH)):
    bd = synth_bd(ig)
    d = ra._de_at(bd, 0.4 * bd.N, 0.0, 0.02)
    res, flip, scale, _ = de_compat(d)
    check(f"T1b: DE (Φ₂・Φ₃ 一定、r パラメータ) が軸対称の C⁺ の適合式を満たす ({name})",
          res <= 1e-12 * scale and flip >= 1e-2 * scale, f"残差 {res:.1e}、符号反転の対照 {flip:.2e}、尺度 {scale:.2f}")
    # 跨ぐ D
    bdx = synth_bd(ig, **CROSS)
    dx_ = ra._de_at(bdx, 0.75 * bdx.N, 0.0, 0.02)
    ok_line = bool(dx_["E"].get("ok"))
    res, flip, scale, f = de_compat(dx_) if ok_line else (np.inf, 0.0, 1.0, None)
    psi = (f["th"] - f["mu"]) if ok_line else np.array([0.0])
    crosses = ok_line and psi[0] > 0.0 and psi[-1] < 0.0
    check(f"T1b: 跨ぐ D (M_B 2.2・θ_B 32°、η_D 0.75 N) の DE が θ = μ を通って C⁺ の適合式を満たす ({name})",
          crosses and res <= 1e-12 * scale and flip >= 1e-2 * scale,
          f"θ−μ {np.rad2deg(psi[0]):+.3f}° → {np.rad2deg(psi[-1]):+.3f}°、残差 {res:.1e}、対照 {flip:.2e}、尺度 {scale:.2f}")
    if crosses:
        rl = dx_["rl"]
        i = int(np.flatnonzero(np.diff(np.sign(psi)) != 0)[0])
        rr = np.linspace(dx_["D"]["r"], dx_["E"]["r"], len(psi))
        g = lambda r: float((lambda q: q["th"][0] - q["mu"][0])(rl.state_at([r])))
        r_f = brentq(g, rr[i], rr[i + 1], xtol=1e-15, rtol=4 * np.finfo(float).eps)
        q = rl.state_at([r_f])
        ref = -np.tan(q["th"][0]) / (2.0 * r_f)
        e_th = abs(q["dth_dr"][0] / ref - 1.0)
        e_M = abs(q["dM_dr"][0]) * abs(1.0 / ref)
        # 陰関数の微分の式は θ = μ で代数的に −tanθ/(2r) になるので、解いた線そのものの傾き (中心差分) でも確かめる
        hh = 1e-4 * r_f
        qq = rl.state_at([r_f - hh, r_f + hh])
        fd_th = (qq["th"][1] - qq["th"][0]) / (2 * hh)
        fd_M = (qq["M"][1] - qq["M"][0]) / (2 * hh)
        e_fd = abs(fd_th / ref - 1.0)
        e_fdM = abs(fd_M / ref)
        check(f"T1b: 折り返し点 (θ = μ) で dθ/dr = −tanθ/(2r)・dM/dr = 0 ({name})",
              e_th <= 1e-9 and e_M <= 1e-9 and e_fd <= 1e-6 and e_fdM <= 1e-6,
              f"r_f {r_f:.6f}・M_f {q['M'][0]:.6f}・θ_f {np.rad2deg(q['th'][0]):.4f}°、点検値 {ref:.10f}: 式 {q['dth_dr'][0]:.10f} "
              f"(相対 {e_th:.1e})・解いた線の中心差分 {fd_th:.10f} (相対 {e_fd:.1e})、|dM/dr|/|点検値|: 式 {e_M:.1e}・差分 {e_fdM:.1e}")

# ================================================================ T2
DR_T2 = (0.04, 0.02, 0.01)


def t2_rel(bd, etas, p_a, dr):
    """包絡線の恒等式の相対差の最大 (η_D の 5 点差分、hh 0.5) と、差分の全点の (θ_D − μ_D, θ_E − μ_E) の符号。"""
    rr, sgn = [], []
    for eta in etas:
        hh = 0.5
        v = {k: ra._de_at(bd, eta + k * hh, p_a, dr) for k in (-2, -1, 0, 1, 2)}
        D5 = lambda q: (-q[2] + 8 * q[1] - 8 * q[-1] + q[-2]) / (12 * hh)
        lhs = D5({k: w["F_part"] for k, w in v.items()})
        rhs = 2 * np.pi * (v[0]["rl"].Phi3 * D5({k: w["E"]["x"] for k, w in v.items()})
                           + v[0]["edge"] * D5({k: w["E"]["r"] for k, w in v.items()}))
        rr.append(abs(lhs - rhs) / abs(lhs))
        for w in v.values():
            sgn.append((np.sign(w["D"]["th"] - np.arcsin(1 / w["D"]["M"])), np.sign(w["E"]["th"] - np.arcsin(1 / w["E"]["M"]))))
    return max(rr), sgn


for name, ig in (("CPG γ 1.2", ra.as_isentrope(1.2)), ("EXH", IG_EXH)):
    for label, bd, etas in (("", synth_bd(ig), None), ("跨ぐ D、", synth_bd(ig, **CROSS), (0.72, 0.75, 0.78))):
        etas = [e * bd.N for e in (etas or (0.3, 0.5, 0.7))]
        for p_a in (0.0, 2e-3):
            rels, sg = {}, []
            for dr in DR_T2:
                rels[dr], sgn = t2_rel(bd, etas, p_a, dr)
                sg += sgn
            ok = all(v <= 1e-5 for v in rels.values()) and rels[DR_T2[-1]] <= max(rels[DR_T2[0]], 1e-9) * 10
            extra = ""
            if label:
                cross_all = all(a > 0 and b < 0 for a, b in sg)
                ok = ok and cross_all
                extra = f"、差分の全点で θ_D > μ_D・θ_E < μ_E: {cross_all}"
            check(f"T2: 包絡線の恒等式 dF/dη = 2π[Φ₃ dL/dη + F(E) dr_E/dη] ({label}{name}、p_a {p_a:g})", ok,
                  "相対差 " + "・".join(f"dr {k:g}: {v:.1e}" for k, v in rels.items()) + extra)

# ================================================================ (b) 中の根の枝の追跡
# 添字でなく η_D の連続性で結ぶ: 低い η_D に根が 1 つ増えても、元の根は近い η_D の根へつながる (添字なら 10 → 5 へ跳ぶ)。
# det_D の符号 (DE の族) が違う根は η_D が近くても結ばない (det_D = 0 を跨ぐと E が跳ぶので連続な枝ではない)。
lab = ra._link_branches([[(10.0, 1)], [(5.0, 1), (12.0, 1)], [(4.0, 1), (14.0, 1)]], tol=5.0)
check("(b) 枝: 根が増えて添字がずれても η_D の近い根を同じ枝に結ぶ", lab[0][0] == lab[1][1] == lab[2][1] and lab[1][0] == lab[2][0]
      and lab[1][0] != lab[0][0], f"枝番号 {lab}")
lab = ra._link_branches([[(34.4, 1)], [(52.9, -1)], [(53.6, -1)]], tol=24.0)
check("(b) 枝: det_D の符号 (族) が違う根は |Δη_D| ≤ tol でも結ばない", lab[0][0] != lab[1][0] and lab[1][0] == lab[2][0],
      f"枝番号 {lab} (θ_B 37°→38° の円弧 kernel の主枝の終わり η 34.4 と偽枝 η 52.9 の型)")

# ================================================================ (e) 受理規則と kernel の既定
# 受理は det_D > 0 だけ。θ_E ≤ θ_B・de < 0 はフラグ。det_D > 0 かつ θ_E > θ_B は受理 + policy_hole (方針の穴として諮る)。
import inspect  # noqa: E402

d2r = np.deg2rad
quad = {}
for name, (det, thB, thE, de) in {"det>0・θE≤θB": (0.1, d2r(36.0), d2r(11.0), -0.12),
                                  "det<0・θE≤θB (θ_B ≥ 40.8° の偽枝の型)": (-0.005, d2r(41.0), d2r(40.5), 8.8),
                                  "det>0・θE>θB": (0.1, d2r(36.0), d2r(37.0), -0.12),
                                  "det<0・θE>θB": (-0.3, d2r(38.0), d2r(41.0), 2.0)}.items():
    cd = {}
    quad[name] = (ra._accept_candidate(cd, thB, thE, det, de), cd)
ok_e = (quad["det>0・θE≤θB"][0] and not quad["det>0・θE≤θB"][1].get("policy_hole")
        and not quad["det<0・θE≤θB (θ_B ≥ 40.8° の偽枝の型)"][0] and "det_D" in quad["det<0・θE≤θB (θ_B ≥ 40.8° の偽枝の型)"][1]["why"]
        and quad["det>0・θE>θB"][0] and bool(quad["det>0・θE>θB"][1].get("policy_hole"))
        and not quad["det<0・θE>θB"][0]
        and all(set(("det_D", "de", "theta_E_le_theta_B", "de_lt_0")) <= set(c) for _, c in quad.values())
        and [c["theta_E_le_theta_B"] for _, c in quad.values()] == [True, True, False, False]
        and [c["de_lt_0"] for _, c in quad.values()] == [True, False, True, False])
check("(e) 受理規則: det_D > 0 だけで受理、θ_E ≤ θ_B・de < 0 はフラグ、det_D > 0 かつ θ_E > θ_B は受理 + policy_hole (4 象限)",
      ok_e, "; ".join(f"{k}: 受理 {v[0]}・hole {bool(v[1].get('policy_hole'))}" for k, v in quad.items()))
sig_d = inspect.signature(ra.design_rao).parameters
sig_k = inspect.signature(ra.solve_kappaB).parameters
api_ok = (sig_d["kernel"].default == "arc" and sig_d["c"].default is None and "c_scan" not in sig_d
          and sig_k["c"].default is inspect.Parameter.empty and sig_k["c"].kind is inspect.Parameter.KEYWORD_ONLY
          and not any(k in sig_k for k in ("c_scan", "xtol_c")) and not hasattr(ra, "C_SCAN"))
raised = {}
for label, kw in (("linear・c 無し", {"kernel": "linear"}), ("arc・c ≠ 0", {"kernel": "arc", "c": -0.35})):
    try:
        ra.design_rao(IG_EXH, R_t=1.0, L=22.66, eps=50.0, **kw)
        raised[label] = False
    except ValueError:
        raised[label] = True
check("(e) kernel の既定: design_rao は円弧が既定、linear は c の明示が要る、solve_kappaB は c を受け取るだけ (外側の走査・brentq なし)",
      api_ok and all(raised.values()), f"既定 kernel {sig_d['kernel'].default!r}、ValueError {raised}")

# ================================================================ T1d (線と場のレベル)
g = 1.2
icpg = ra.as_isentrope(g)
itp = ra.as_isentrope(constant_cp_gas(g, Tt=3000.0))
bd_c, bd_t = synth_bd(icpg), synth_bd(itp)
dc = ra._de_at(bd_c, 0.45 * bd_c.N, 0.0, 0.02)
dt = ra._de_at(bd_t, 0.45 * bd_t.N, 0.0, 0.02)
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
NOT_RUN = {
    "T1a": "記録 (合否なし) に格下げ — #6 (CFD ピン) の中で case/29 の R_t 0.382 の Euler を凍結源に CPG γ 1.2・1.3・1.4 (diagnostician 2026-10-11)",
    "T1a′": "(D) 壁の MOC 推力 ≤ (M) の MOC 推力 — T6 の後 (§5.1 #5 (g))",
    "T1c": "自由 ε と固定 ε の等価 — T6 の後",
    "T3": "保存 (推力の 2 会計の差・閉合 |W_m − E|) — T6 の系列で判定 (run_rao_t6.py の表、合否は親と上位)",
    "T4": "停留性 — T6 の後 (§5.1 #5 (g))",
    "T5": "C2 は課さない (2026-10-11 改) — g̃ と Δκ_B を T6 の系列 (n 161/241/321) で記録 (run_rao_t6.py)",
    "T6": "分解能 n 161/241/321 + n 321 の dr 半分・m・q 倍 (円弧・L100・L70) — run_rao_t6.py の表、合否は親と上位",
}
for k, v in NOT_RUN.items():
    print(f"未実施 {k}: {v}", flush=True)

print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAIL'} (未実施 {len(NOT_RUN)} 件: 理由は各行)")
sys.exit(1 if FAIL else 0)
