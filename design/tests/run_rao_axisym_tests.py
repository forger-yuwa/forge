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
- **(g1) Goursat の既定** (plan §5.1 #5 (g1)、2026-10-12): `design_rao` の m・q の既定が初期線の点数 n (旧 n/2)、明示した値は
  そのまま。根の探索を差し替えて (根なしを返す) 戻りの記録で確かめる。
- **(W3) 初期線を外から渡す口** (plan §5.1 #6 W3、2026-10-12): `design_rao` の既定 (`line=None`) が solve_kappaB に渡す線と
  引数が `_hall_line` の線・従来の既定と厳密に一致し外の線の検査を呼ばない (既定の経路はビット同一)、渡した線は同じ dict の
  まま検査 `check_initial_line` を 1 回通る、検査は Hall の線を通し壁から外れた線・θ_r 無し・亜音速の点・r が単調でない線を
  solve_kappaB の前に ValueError で止める。検証 run (case/68 の run_0015) が手元にあれば CFD ピンの線 (n 321) が
  `moc_forward._check_line` を通ること (無ければ未実施に記録)。
- **(g) 与えた壁の推力 `forward_thrust`** (T4・T1a′ の評価器、前進 MOC で場を解き直す): 放射源流 (γ 1.4・半角 10°・厳密な
  C⁻ 初期線) の円錐壁で、推力 (初期線 + 壁の圧力積分) が厳密値に n 49/97/193 で単調に近づき最細区間の次数 ≥ 1.8・
  n 193 で相対 1e-6 以下、ゲート通過・壁の依存域に反転なし。衝撃波を要する壁 (曲げ戻し θ₁ 0°) は依存域に反転ありと返す。
  `ExtendedWall` は x ≤ L で元の壁とビット同一・先は接線。
- **(g) (M) の壁の表現 `RaoWall`**: 円弧 kernel の B が厳密、端の拘束 (B の r・r'、終点の r・r') が 1e-12。滑らかな合成の壁
  (B で C1) の点列に位置の雑音 1e-5 を足したとき、罰則付きの当てはめ (`"pspline"`) は真の壁から 1e-5 以内で曲率の振れ
  (dκ の符号変化) が少なく、点列の補間 (`"interp"`) は曲率が大きく振れる (前進 MOC が BDE の中で反転を出した原因の再現)。
- **(g) 摂動 `PerturbedWall`**: 族 `sin`・`sin_c1` が両端で 0 (B・終点の位置が不変 → L・ε 不変)、`sin_c1` は端で勾配も不変、
  r' が差分と 1e-7。
- **(h) 整合な前線の再標本化 `resampled_front_line`** (T4 の評価器の初期線、plan §5.1 #5 (h)、diagnostician 2026-10-12):
  円弧 kernel (R_t 1) + 接線の壁で Hall n 321 の前進 MOC を回し、壁点 x ≈ 0.3 の前線 k₀ が軸まで ST_OK で、前線そのものを
  初期線にしたマーチが元の網の前線を 1e-12 で再現すること (整合な初期線)。N 321・641・1281 の線が `_check_line` を通り
  (壁点は壁に厳密に乗る)、r 等間隔・端が前線のまま。N 321 の線が前線を再現する往復の誤差を記録し N で単調に減ること
  (前線には Hall のマーチの出だしの軸点の θ_r の振動から来る 2Δ の鋸歯の帯 [k₀ − j ∈ −6〜12] があり、往復の誤差はそこで
  最大)。滑らかな合成の前線で補間誤差が O(h⁴) で減ること。軸まで OK でない前線は ValueError。

**この版で回さないもの**: T3・T5・T6 は `design_rao` を基準ケースの L_k で回す run で、ここでは回さない —
`case/68.bell_two_mode/rao_m/run_rao_t6.py` (円弧 kernel・L100・L70 × n 161/241/321 + n 321 の dr 半分・m・q 倍) が
表 `rao_t6.json`・`rao_t6.csv` を出し、合否は親と上位が判定する (plan §6 V1 の T3・T5・T6、2026-10-11 改)。
T5 は C2 を課さない記録 (g̃ と Δκ_B)。T4・T1a′ は `run_rao_t4_gen.py` (壁の生成) + `run_rao_t4.py` (表 `rao_t4.csv`・
`rao_t1ap.csv`)、T1c は `run_rao_t4_gen.py` (表 `rao_t1c.csv`) で、合否は親と上位 (§5.1 #5 (g))。T1a は記録 (合否なし) に
格下げ (#6 で CFD ピンの後)。

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

# ================================================================ (g1) Goursat の既定 m = q = n
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    line81 = ra._hall_line(1.0, IG_EXH, 81)
_orig_skb = ra.solve_kappaB
ra.solve_kappaB = lambda *a, **k: {"c": 0.0, "ok": False, "why": "試験の差し替え (根なし)", "roots": []}
try:
    r_def = ra.design_rao(IG_EXH, R_t=1.0, L=22.66, eps=50.0, line=line81)
    r_exp = ra.design_rao(IG_EXH, R_t=1.0, L=22.66, eps=50.0, line=line81, m=40, q=30)
finally:
    ra.solve_kappaB = _orig_skb
check("(g1) design_rao の Goursat の既定 m = q = n (初期線の点数)、明示した m・q はそのまま",
      r_def["status"] == "infeasible" and (r_def["n"], r_def["m"], r_def["q"]) == (81, 81, 81) and (r_exp["m"], r_exp["q"]) == (40, 30),
      f"既定 (n, m, q) = ({r_def['n']}, {r_def['m']}, {r_def['q']})、明示 ({r_exp['m']}, {r_exp['q']})")

# ================================================================ (W3) 初期線を外から渡す口 (plan §5.1 #6 W3、2026-10-12)
# 既定 (line=None) は Hall の線を作って solve_kappaB に渡すだけで、外から渡す線の検査 `check_initial_line` を通らない。
# 既定の経路がビット同一であることは、solve_kappaB に渡る線と引数が `_hall_line` の線・従来の既定 (dr 3.2/n・m = q = n) と
# 厳密に一致することで確かめる (solve_kappaB 以降は変えていない)。
_cap = {}


def _capture(*a, **k):
    _cap.clear()
    _cap.update(k)
    _cap["args"] = a
    return {"c": 0.0, "ok": False, "why": "試験の差し替え (根なし)", "roots": []}


_orig_cil = ra.check_initial_line
_n_cil = [0]


def _count_cil(*a, **k):
    _n_cil[0] += 1
    return _orig_cil(*a, **k)


ra.solve_kappaB, ra.check_initial_line = _capture, _count_cil
try:
    r_w3 = ra.design_rao(IG_EXH, R_t=1.0, L=22.66, eps=50.0, n=81)
    cap_def = dict(_cap)
    n_cil_def = _n_cil[0]
    ref81 = ra._hall_line(1.0, IG_EXH, 81)
    r_w3g = ra.design_rao(IG_EXH, R_t=1.0, L=22.66, eps=50.0, line=ref81)
    cap_giv = dict(_cap)
    n_cil_giv = _n_cil[0] - n_cil_def
finally:
    ra.solve_kappaB, ra.check_initial_line = _orig_skb, _orig_cil
same_line = all(np.array_equal(cap_def["line"][k], ref81[k]) for k in ("x", "r", "M", "theta")) and \
    cap_def["line"]["thr_axis"] == ref81["thr_axis"]
check("(W3) design_rao の既定 (line=None) は Hall の線 (`_hall_line` と全配列・θ_r が厳密一致) を solve_kappaB に渡し、"
      "dr 3.2/n・m = q = n・line_source 'hall'、外の線の検査は呼ばない (既定の経路はビット同一)",
      same_line and cap_def["args"][1:] == (1.0, 22.66) and cap_def["dr"] == 3.2 / 81 and cap_def["c"] == 0.0
      and cap_def["eps"] == 50.0 and (r_w3["m"], r_w3["q"]) == (81, 81) and r_w3["line_source"] == "hall" and n_cil_def == 0,
      f"線の一致 {same_line}・dr {cap_def['dr']}・(m, q) ({r_w3['m']}, {r_w3['q']})・検査の呼び出し {n_cil_def}")
check("(W3) 外から渡した線はそのまま (同じ dict) solve_kappaB に渡り、検査 `check_initial_line` を 1 回通り、"
      "line_source は線の source",
      cap_giv["line"] is ref81 and n_cil_giv == 1 and r_w3g["line_source"] == "hall_throat_characteristic",
      f"line_source {r_w3g['line_source']}・検査 {n_cil_giv} 回")


def _cil_ok(line_):
    try:
        ra.check_initial_line(line_, 1.0, IG_EXH)
        return True
    except ValueError:
        return False


# 検査: Hall の線は通る (n 81/161/321)、壁から外れた線・θ_r の無い線・亜音速の点を含む線・r が単調でない線は ValueError
# (solve_kappaB の前に止まる)
ok_cil = all(_cil_ok(ra._hall_line(1.0, IG_EXH, n_)) for n_ in (81, 161, 321))
bad_cases = {}
for lab, edit in (("壁から外れた (r × 1.001)", lambda d_: d_.update({"r": d_["r"] * 1.001})),
                  ("θ_r 無し", lambda d_: d_.pop("thr_axis")),
                  ("亜音速の点", lambda d_: d_.update({"M": np.where(np.arange(len(d_["M"])) == 40, 0.98, d_["M"])})),
                  ("r が単調でない", lambda d_: d_.update({"r": np.r_[d_["r"][:40], d_["r"][41], d_["r"][40], d_["r"][42:]]}))):
    bl = {k: (np.array(v, dtype=float).copy() if isinstance(v, np.ndarray) else v) for k, v in ref81.items()}
    edit(bl)
    ra.solve_kappaB = _capture
    _cap.clear()
    try:
        ra.design_rao(IG_EXH, R_t=1.0, L=22.66, eps=50.0, line=bl)
        bad_cases[lab] = "通った"
    except ValueError as e:
        bad_cases[lab] = ("ValueError" if not _cap else "ValueError (solve_kappaB の後)") \
            if "初期線の検査に不合格" in str(e) else f"別の ValueError: {e}"
    finally:
        ra.solve_kappaB = _orig_skb
check("(W3) check_initial_line: Hall の線 (n 81/161/321) は通り、壁から外れた線・θ_r の無い線・亜音速の点・r が単調でない線は "
      "design_rao が solve_kappaB の前に ValueError で止める",
      ok_cil and all(v == "ValueError" for v in bad_cases.values()), str(bad_cases))
# CFD ピンの線 (検証 run が手元にあるときだけ): run_0015 (Hall 版 L100 の TP Euler) から n 321・r 等間隔の線が検査を通る
RUN_PIN = ROOT / "case" / "68.bell_two_mode" / "run_0015_tp_m_euler_l100"
PIN_RAN = (RUN_PIN / "res_30000.h5").exists()
if PIN_RAN:
    from forge_design.feedback.cfd_initial_line import expected_bell_source, moc_initial_line, pinned_factory  # noqa: E402
    from forge_design.geometry import moc_forward as mf_w3  # noqa: E402
    from forge_design.probdef import load_problem  # noqa: E402
    eb_w3 = expected_bell_source(load_problem(ROOT / "case" / "68.bell_two_mode" / "problem_bell_tp_m_euler_l100.yaml"))
    th_w3 = pinned_factory(RUN_PIN, "res_30000.h5", eb_w3["expect"])(eb_w3["R"], float(IG_EXH.moc.gamma_throat()))
    pin_w3 = moc_initial_line(th_w3, IG_EXH.moc, 321)
    info_w3 = ra.check_initial_line(pin_w3, 1.0, IG_EXH)
    W_w3 = mf_w3._Wall(ra.KernelWall(1.0, c=0.0, s_lin=0.5, theta_cap=np.deg2rad(40.0)), 10.0)
    codes_w3 = mf_w3._check_line(pin_w3, W_w3, IG_EXH.moc)[0]
    m_w3 = ra.line_integrals(IG_EXH, pin_w3["x"], pin_w3["r"], pin_w3["theta"], pin_w3["M"], 0.0)[0] / np.pi
    check("(W3) CFD ピンの線 (run_0015、n 321・r 等間隔) が `moc_forward._check_line` と `check_initial_line` を通る",
          not codes_w3 and info_w3["wall_dr"] == 0.0 and pin_w3["source"] == "cfd_pinned" and len(pin_w3["x"]) == 321,
          f"M 最小 {info_w3['M_min']:.5f}・θ_r {info_w3['thr_axis']:.6f}・軸端 x {pin_w3['x'][0]:.6f}・m*/ṁ_1D {m_w3:.7f} (記録)")

# ================================================================ (g) 与えた壁の推力 forward_thrust (放射源流の円錐壁)
G_RAD, T10, X1_RAD, L_RAD = 1.4, np.tan(np.deg2rad(10.0)), 1.2, 3.0
IG_RAD = ra.as_isentrope(G_RAD)


def _mach_of_ar(AR):
    AR = np.asarray(AR, dtype=float)
    e = (G_RAD + 1.0) / (2.0 * (G_RAD - 1.0))
    M = 2.0 + 1.5 * np.sqrt(np.maximum(AR - 1.0, 0.0))
    for _ in range(80):
        t = (2.0 / (G_RAD + 1.0)) * (1.0 + 0.5 * (G_RAD - 1.0) * M * M)
        f = t ** e / M - AR
        df = -t ** e / M ** 2 + e * t ** (e - 1.0) * (2.0 / (G_RAD + 1.0)) * (G_RAD - 1.0)
        dM = f / df
        M = M - dM
        if np.all(np.abs(dM) < 1e-15 * M):
            break
    return M


def _src_M(x, r):
    return _mach_of_ar(np.hypot(x, r) ** 2)


class _Ray:
    def r(self, x, deriv=0):
        x = np.asarray(x, dtype=float)
        return x * T10 if deriv == 0 else np.full_like(x, T10)


class _TurnWall:
    """レイ (10°) → [x1, x1+dx] で勾配を smoothstep で tan θ₁ へ → 直線 (run_moc_forward_tests.py と同じ)。"""

    def __init__(self, x1, dx, th1_deg):
        self.x1, self.dx, self.s1 = float(x1), float(dx), float(np.tan(np.deg2rad(th1_deg)))

    def r(self, x, deriv=0):
        x = np.asarray(x, dtype=float)
        t = np.clip((x - self.x1) / self.dx, 0.0, 1.0)
        S = 3.0 * t * t - 2.0 * t ** 3
        IS = t ** 3 - 0.5 * t ** 4
        if deriv == 1:
            return T10 + (self.s1 - T10) * S
        r1 = self.x1 * T10
        r_in = r1 + T10 * (x - self.x1) + (self.s1 - T10) * self.dx * IS
        r_after = r1 + T10 * self.dx + 0.5 * (self.s1 - T10) * self.dx + self.s1 * (x - self.x1 - self.dx)
        return np.where(x < self.x1, x * T10, np.where(x <= self.x1 + self.dx, r_in, r_after))


def _radial_line(n):
    rw = X1_RAD * T10

    def rhs(r, y):
        M = float(_src_M(y[0], r))
        return [1.0 / np.tan(np.arctan2(r, y[0]) - np.arcsin(1.0 / M))]
    sol = solve_ivp(rhs, (rw, 0.0), [X1_RAD], method="DOP853", rtol=1e-13, atol=1e-14, dense_output=True)
    rq = np.linspace(0.0, rw, n)
    xq = sol.sol(rq)[0]
    xq[-1] = X1_RAD
    return {"x": xq, "r": rq, "M": _src_M(xq, rq), "theta": np.arctan2(rq, xq), "thr_axis": 1.0 / xq[0]}


from scipy.integrate import quad  # noqa: E402

lf = _radial_line(4001)
_, T_line_ex = ra.line_integrals(IG_RAD, lf["x"], lf["r"], lf["theta"], lf["M"], 0.0, spline=True)
T_wall_ex = quad(lambda x: float(IG_RAD.state(np.array([float(_src_M(x, x * T10))]))["p"][0]) * 2.0 * np.pi * x * T10 * T10,
                 X1_RAD, L_RAD, epsabs=1e-14, epsrel=1e-13, limit=200)[0]          # ∫ p 2π r r' dx (r = x tan10°)
F_ex = T_line_ex + T_wall_ex
errs, ok_all = {}, True
for n in (49, 97, 193):
    rr_ = ra.forward_thrust(_Ray(), L_RAD, G_RAD, _radial_line(n))
    ok_all = ok_all and rr_["ok"] and rr_["gate"]["gate_pass"] and rr_["dod_clean"]
    errs[n] = (rr_["F"] - F_ex) / F_ex
od = [np.log2(abs(errs[49] / errs[97])), np.log2(abs(errs[97] / errs[193]))]
check("(g) forward_thrust: 放射源流の円錐壁で推力が厳密値に単調に近づき、最細区間の次数 ≥ 1.8・n 193 で相対 1e-6 以下、"
      "ゲート通過・壁の依存域に反転なし",
      ok_all and abs(errs[49]) > abs(errs[97]) > abs(errs[193]) and od[-1] >= 1.8 and abs(errs[193]) <= 1e-6,
      f"相対誤差 n 49 {errs[49]:+.2e}・97 {errs[97]:+.2e}・193 {errs[193]:+.2e}、次数 {od[0]:.2f}・{od[1]:.2f} (F 厳密 {F_ex:.10f})")
rr_ = ra.forward_thrust(_TurnWall(1.6, 0.3, 0.0), L_RAD, G_RAD, _radial_line(97))
check("(g) forward_thrust: 衝撃波を要する壁 (曲げ戻し θ₁ 0°) は壁の依存域に反転ありと返す (dod_clean False)",
      rr_["ok"] and not rr_["dod_clean"] and rr_["dod"]["n_bad_in_dod"] > 0 and not rr_["gate"]["gate_pass"],
      f"依存域の反転 {rr_['dod']['n_bad_in_dod']} 点、ゲート {rr_['gate']['verdict']} {rr_['gate']['codes']}")
ew = ra.ExtendedWall(_Ray(), 2.0)
xq = np.array([0.5, 1.3, 2.0, 2.5, 3.0])
check("(g) ExtendedWall: x ≤ L は元の壁とビット同一、L の先は接線",
      np.array_equal(ew.r(xq[:3]), _Ray().r(xq[:3])) and np.allclose(ew.r(xq[3:]), 2.0 * T10 + T10 * (xq[3:] - 2.0), rtol=0, atol=1e-15)
      and np.all(ew.r(xq, 1) == T10))

# ================================================================ (g) (M) の壁の表現 RaoWall と摂動
thB_s, thE_s, x_e = np.deg2rad(30.0), np.deg2rad(6.0), 20.0     # 合成の壁: θ が B の θ_B から θ_E へ滑らかに (B で C1)
xB_s, rB_s = np.sin(thB_s), 2.0 - np.cos(thB_s)


def _th_true(x):
    u = (x - xB_s) / 3.0
    return thE_s + (thB_s - thE_s) * np.exp(-u) * (1.0 + u)


xx = np.linspace(xB_s, x_e, 40001)
gx_, gw_ = np.polynomial.legendre.leggauss(6)
dxs = np.diff(xx)
mid = 0.5 * (xx[:-1] + xx[1:])[:, None] + 0.5 * dxs[:, None] * gx_
rr_true = rB_s + np.concatenate([[0.0], np.cumsum(np.sum(np.tan(_th_true(mid)) * 0.5 * dxs[:, None] * gw_, axis=1))])
from scipy.interpolate import CubicHermiteSpline  # noqa: E402

r_true = CubicHermiteSpline(xx, rr_true, np.tan(_th_true(xx)))
xs_pts = xB_s + (x_e - xB_s) * (np.arange(161) / 160.0) ** 2          # B の近くが密 (生成器の壁点と同じ偏り)
rng = np.random.default_rng(0)
noise = 1e-5 * rng.standard_normal(len(xs_pts))
noise[0] = noise[-1] = 0.0
pts = np.column_stack([xs_pts, r_true(xs_pts) + noise])
fits = {}
xd = np.linspace(xB_s + 1e-6, x_e - 1e-6, 20001)
kap_true = np.gradient(_th_true(xd), xd) * np.cos(_th_true(xd))       # κ = dθ/ds = cosθ dθ/dx
for fit in ("pspline", "interp"):
    w_ = ra.RaoWall(1.0, thB_s, pts, float(_th_true(x_e)), fit=fit)
    rp_ = w_.r(xd, 1)
    kap_ = w_.spl.derivative(2)(xd) / (1.0 + rp_ * rp_) ** 1.5
    fits[fit] = (w_, float(np.max(np.abs(w_.r(xd) - r_true(xd)))), float(np.max(np.abs(kap_ - kap_true))),
                 w_.fit_info["n_dkappa_sign_changes"])
wp_ = fits["pspline"][0]
end_ok = (abs(float(wp_.r(np.array([xB_s]))[0]) - rB_s) <= 1e-12 and abs(float(wp_.r(np.array([xB_s]), 1)[0]) - np.tan(thB_s)) <= 1e-12
          and abs(float(wp_.spl(np.array([xB_s + 1e-12]))[0]) - rB_s) <= 1e-11 and abs(float(wp_.dspl(np.array([xB_s]))[0]) - np.tan(thB_s)) <= 1e-12
          and abs(wp_.r_end - pts[-1, 1]) <= 1e-12 and abs(wp_.rp_end - np.tan(_th_true(x_e))) <= 1e-12
          and abs(wp_.r_B - rB_s) <= 1e-15 and abs(wp_.x_B - xB_s) <= 1e-15)
check("(g) RaoWall: B (円弧の端) が厳密、端の拘束 (B の r・r'、終点の r・r') が 1e-12、kernel と BE が B で C1", end_ok)
check("(g) RaoWall: 雑音 σ 1e-5 の点列 (B の近くが密) から、罰則付きの当てはめは真の壁と位置 3σ・曲率 3e-3 以内、"
      "点列の補間は曲率の誤差がその 100 倍以上",
      fits["pspline"][1] <= 3e-5 and fits["pspline"][2] <= 3e-3 and fits["interp"][2] >= 100 * fits["pspline"][2],
      f"pspline: 位置 {fits['pspline'][1]:.1e}・曲率 {fits['pspline'][2]:.1e}・dκ の符号変化 {fits['pspline'][3]} "
      f"(λ {wp_.fit_info['lam']:.1e}); interp: 位置 {fits['interp'][1]:.1e}・曲率 {fits['interp'][2]:.2f}・符号変化 {fits['interp'][3]}; "
      f"真の κ の幅 {kap_true.max() - kap_true.min():.3f}")
pw_ok, pw_det = True, []
for fam in ("sin", "sin_c1"):
    for k in (1, 2, 3):
        eta, deta = ra.perturbation_family(fam, k)
        pw = ra.PerturbedWall(wp_, fam, k, 4e-3)
        xe_ = np.array([wp_.x_B, wp_.x_end])
        same_end = np.array_equal(pw.r(xe_), wp_.r(xe_))
        xs_ = np.linspace(wp_.x_B + 0.01, wp_.x_end - 0.01, 9)
        hh = 1e-6
        fd = (pw.r(xs_ + hh) - pw.r(xs_ - hh)) / (2 * hh)
        e_fd = float(np.max(np.abs(fd - pw.r(xs_, 1))))
        end0 = abs(eta(0.0)) <= 1e-15 and abs(eta(1.0)) <= 1e-12
        slope_end = (abs(deta(0.0)) <= 1e-12 and abs(deta(1.0)) <= 1e-12) if fam == "sin_c1" else True
        pw_ok = pw_ok and same_end and e_fd <= 1e-7 and end0 and slope_end
        pw_det.append(f"{fam} k{k}: 端 {same_end}・r' 差分 {e_fd:.0e}")
check("(g) PerturbedWall: 両端で 0 (B・終点の位置が不変 → L・ε 不変)、sin_c1 は端で勾配も不変、r' が差分と 1e-7", pw_ok,
      "; ".join(pw_det))
t_ = wp_.t_of_x(np.array([wp_.x_B, wp_.x_end]))
check("(g) RaoWall の弧長の正規化 t が B で 0・終点で 1", abs(t_[0]) <= 1e-15 and abs(t_[1] - 1.0) <= 1e-12, f"{t_}")

# ================================================================ (h) 整合な前線の再標本化 resampled_front_line
from types import SimpleNamespace  # noqa: E402

from scipy.interpolate import CubicSpline  # noqa: E402

from forge_design.geometry import moc_forward as mf  # noqa: E402


class _ArcTangent:
    """円弧 kernel (R_t 1、閉じた式) を θ_B 30° まで → 接線。前線 k₀ (壁点 x ≈ 0.3) は円弧の上なので B の先の形に依らない。"""

    def __init__(self, thB):
        self.xB, self.rB, self.t = np.sin(thB), 2.0 - np.cos(thB), np.tan(thB)

    def r(self, x, deriv=0):
        x = np.asarray(x, dtype=float)
        xa = np.atleast_1d(x)
        out = np.empty_like(xa)
        k = xa <= self.xB
        root = np.sqrt(np.maximum(1.0 - xa[k] ** 2, 0.0))
        out[k] = (2.0 - root) if deriv == 0 else xa[k] / root
        out[~k] = (self.rB + self.t * (xa[~k] - self.xB)) if deriv == 0 else self.t
        return out.reshape(np.shape(x))


W_AT, L_AT = _ArcTangent(np.deg2rad(30.0)), 8.0
res_h = mf.forward_moc(W_AT, L_AT, IG_EXH.moc, ra._hall_line(1.0, IG_EXH, 321))
net_h, n_h = res_h["net"], res_h["n"]
k0_h = int(np.argmin(np.abs(net_h.a["x"][:, n_h] - 0.3)))
fr_h = {f: net_h.a[f][k0_h].copy() for f in net_h.FIELDS}
front_ok = bool(np.all(net_h.status[k0_h] == mf.ST_OK))
lf_h = {"x": fr_h["x"], "r": fr_h["r"], "M": fr_h["M"], "theta": fr_h["th"], "thr_axis": float(net_h.thr_axis[k0_h])}
net2_h = mf.forward_moc(W_AT, L_AT, IG_EXH.moc, lf_h)["net"]
K_h = 30
d_march = max(float(np.nanmax(np.abs(net_h.a[f][k0_h:k0_h + K_h] - net2_h.a[f][:K_h]))) for f in net_h.FIELDS)
st_same = bool(np.array_equal(net_h.status[k0_h:k0_h + K_h], net2_h.status[:K_h]))
check("(h) 前線 k₀ (Hall n 321 のマーチ、壁点 x ≈ 0.3 の円弧の上) は軸まで ST_OK で、前線そのものを初期線にしたマーチが元の網の"
      "前線 k₀+1〜k₀+30 を 1e-12 で再現する (整合な初期線 — Hall の初期線と違い前線 1 で不整合を出さない)",
      front_ok and st_same and d_march <= 1e-12,
      f"k₀ {k0_h}・壁点 x {fr_h['x'][-1]:.6f}・軸点 x {fr_h['x'][0]:.6f}、マーチの差 {d_march:.1e}、状態一致 {st_same}")
Wc_h = mf._Wall(W_AT, L_AT)
lines_h, cl_det, cl_ok = {}, [], True
for N in (321, 641, 1281):
    ln = ra.resampled_front_line(net_h, k0_h, N)
    lines_h[N] = ln
    codes, _, info = mf._check_line(ln, Wc_h, IG_EXH.moc)
    ends = (ln["r"][0] == 0.0 and ln["theta"][0] == 0.0 and ln["x"][0] == fr_h["x"][0] and ln["M"][0] == fr_h["M"][0]
            and ln["x"][-1] == fr_h["x"][-1] and ln["r"][-1] == fr_h["r"][-1] and ln["theta"][-1] == fr_h["th"][-1]
            and ln["M"][-1] == fr_h["M"][-1] and ln["thr_axis"] == net_h.thr_axis[k0_h])
    uni = float(np.max(np.abs(np.diff(ln["r"]) - fr_h["r"][-1] / N)))
    ok_n = not codes and len(ln["r"]) == N + 1 and info["wall_dr"] == 0.0 and info["wall_dtheta"] == 0.0 and ends and uni <= 1e-15
    cl_ok = cl_ok and ok_n
    cl_det.append(f"N {N}: 検査の符号 {codes or 'なし'}・壁 Δr {info['wall_dr']:.0e}・端 {ends}")
check("(h) resampled_front_line: N 321・641・1281 の線が `_check_line` を通り (亜音速・壁から外れる・軸・順序なし、壁点は壁に"
      "厳密に乗る)、r 等間隔の N+1 点、両端と軸端の θ_r が前線のまま", cl_ok, "; ".join(cl_det))
# N = 321 で元の前線を再現するか: 再標本化した線を通る 3 次スプラインを前線の点で評価した往復の誤差 (記録) と、N で減ること
rt_h = {}
for N, ln in lines_h.items():
    rt_h[N] = {nm: np.abs(CubicSpline(ln["r"], ln[key])(fr_h["r"]) - fr_h[nm]) for nm, key in (("x", "x"), ("th", "theta"), ("M", "M"))}
band = (k0_h - np.arange(n_h + 1) <= 12) & (k0_h - np.arange(n_h + 1) >= -6)
jmax = int(np.argmax(rt_h[321]["th"]))
mono = all(float(rt_h[a][nm].max()) > float(rt_h[b][nm].max()) for nm in ("x", "th", "M") for a, b in ((321, 641), (641, 1281)))
check("(h) resampled_front_line: N = 321 の線が前線の点を再現する往復の誤差 (記録) が N 321 → 641 → 1281 で単調に減る",
      mono,
      "θ の最大 " + "・".join(f"N {N} {float(v['th'].max()):.1e}" for N, v in rt_h.items())
      + "、M の最大 " + "・".join(f"N {N} {float(v['M'].max()):.1e}" for N, v in rt_h.items())
      + f"; N 321 の θ の最大は j {jmax} (k₀ − j = {k0_h - jmax}・r {fr_h['r'][jmax]:.3f}) — 前線の k₀ − j ∈ [−6, 12] の帯 "
      f"(Hall のマーチの出だしの軸点の θ_r の振動から来る C⁺ の 2Δ の鋸歯) の外では θ {float(rt_h[321]['th'][~band].max()):.1e}・"
      f"M {float(rt_h[321]['M'][~band].max()):.1e}")
# 実装の検査: 同じ非等間隔の r の点に滑らかな関数を置いた合成の前線で、N 321 の補間誤差が O(h⁴) で減る (前線の点の間隔を
# すべて半分にすると誤差が 1/12 以下 = 次数 3.6 以上)。誤差の大きさは参考に 3 次スプライン (端を固定) の上界 (5/384) h⁴ max|f⁗| と並べる
fx = lambda r: 4.27 - 3.6 * r + 0.3 * np.sin(2.0 * r)  # noqa: E731
fth = lambda r: 0.24 * r + 0.05 * r ** 3  # noqa: E731
fM = lambda r: 3.26 - 1.2 * r ** 2 + 0.2 * np.cos(3.0 * r)  # noqa: E731
M4_s = max(0.3 * 16.0, 0.0, 0.2 * 81.0)                              # max|f⁗| (x: 0.3·2⁴、θ: 0、M: 0.2·3⁴)


def _synth_err(rr):
    net_s = SimpleNamespace(status=np.zeros((1, len(rr)), dtype=np.int8), thr_axis=np.array([0.24]),
                            a={"x": fx(rr)[None, :], "r": rr[None, :], "th": fth(rr)[None, :], "M": fM(rr)[None, :]})
    ln_s = ra.resampled_front_line(net_s, 0, 321)
    return max(float(np.max(np.abs(ln_s["x"] - fx(ln_s["r"])))), float(np.max(np.abs(ln_s["theta"] - fth(ln_s["r"])))),
               float(np.max(np.abs(ln_s["M"] - fM(ln_s["r"])))))


rr_s = fr_h["r"]
rr_s2 = np.sort(np.concatenate([rr_s, 0.5 * (rr_s[1:] + rr_s[:-1])]))
e_s1, e_s2 = _synth_err(rr_s), _synth_err(rr_s2)
hmax_s = float(np.max(np.diff(rr_s)))
check("(h) resampled_front_line: 滑らかな合成の前線 (前線 k₀ と同じ非等間隔の r) で N 321 の補間誤差が前線の点の間隔の半分で "
      "1/12 以下 (O(h⁴))", e_s1 / e_s2 >= 12.0,
      f"最大誤差 {e_s1:.1e} → {e_s2:.1e} (比 {e_s1 / e_s2:.1f})、参考: 端を固定した 3 次スプラインの上界 (5/384) h⁴ max|f⁗| = "
      f"{5.0 / 384.0 * hmax_s ** 4 * M4_s:.1e} (h {hmax_s:.4f})")
k_bad = int(np.flatnonzero([not np.all(net_h.status[k] == mf.ST_OK) for k in range(net_h.status.shape[0])])[0])
try:
    ra.resampled_front_line(net_h, k_bad, 321)
    raised = False
except ValueError:
    raised = True
check("(h) resampled_front_line: 軸まで ST_OK でない前線 (出口の外の軸点) は ValueError", raised, f"前線 {k_bad}")

# ================================================================ 回さないもの
NOT_RUN = {
    "T1a": "記録 (合否なし) に格下げ — #6 (CFD ピン) の中で case/29 の R_t 0.382 の Euler を凍結源に CPG γ 1.2・1.3・1.4 (diagnostician 2026-10-11)",
    "T1a′": "(D) 壁の MOC 推力 ≤ (M) の MOC 推力 — run の表 (case/68.bell_two_mode/rao_m/run_rao_t4.py → rao_t1ap.csv)、合否は親と上位",
    "T1c": "自由 ε → 固定 ε の解き直し — run の表 (run_rao_t4_gen.py → rao_t1c.csv)、合否は親と上位",
    "T3": "保存 (推力の 2 会計の差・閉合 |W_m − E|) — T6 の系列で判定 (run_rao_t6.py の表、合否は親と上位)",
    "T4": "停留性 — run の表 (run_rao_t4_gen.py の壁 + run_rao_t4.py → rao_t4.csv、(h) の A/B は run_rao_t4h.py → rao_t4h.csv)、"
          "合否は親と上位。評価器と壁の表現の単体はこの試験の (g)、再標本化した初期線は (h)",
    "T5": "C2 は課さない (2026-10-11 改) — g̃ と Δκ_B を T6 の系列 (n 161/241/321) で記録 (run_rao_t6.py)",
    "T6": "分解能 n 161/241/321 + n 321 の dr 半分・m・q 倍 (円弧・L100・L70) — run_rao_t6.py の表、合否は親と上位",
    "W3 ピン版の (M)": "CFD ピンの初期線で L100 を生成し Hall 版との差を記録 — run の表 (rao_m/run_pin_line.py → rao_pin.json・"
                       "rao_pin_compare.csv)、合否は親と上位",
}
if not PIN_RAN:
    NOT_RUN["W3 CFD ピンの線"] = f"検証 run {RUN_PIN} が無い (git 管理外) — 線の検査を省略"
for k, v in NOT_RUN.items():
    print(f"未実施 {k}: {v}", flush=True)

print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAIL'} (未実施 {len(NOT_RUN)} 件: 理由は各行)")
sys.exit(1 if FAIL else 0)
