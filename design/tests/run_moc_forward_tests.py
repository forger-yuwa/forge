#!/usr/bin/env python3
"""壁を与えた前進 MOC と前進 MOC ゲート `geometry/moc_forward.py` の試験 (plan tooling-nozzle-bell-two-mode §5.1 #3b)。

合格条件 (§5.1 #3b):
1. **厳密解の再現**: 風洞の放射源流の試験 (`moc_axis_limit_radial.py`) と同じ放射源流 (γ 1.4・半角 10°・A/A* = R²) を
   前進で解く。壁 = 原点からのレイ、初期線 = x 1.2 の壁点から厳密場の中で下ろした C⁻ (高次の ODE 積分)。
   ゲートを通り、誤差 (壁の ν・軸の ν・場の ν と θ・前線を通る流量) が解像度とともに減り、最細区間の次数が
   壁 ≥ 1.8・場 ≥ 1.5・流量 ≥ 1.8。
2. **衝撃波を要する壁の負例**: レイの途中で壁を曲げ戻す (10° → θ₁ を Δx 0.3 で)。弱い戻し (9°・7°) は解像度 3 水準で合格、
   強い戻し (5°・0°) は 3 水準とも「網の反転 (C⁺ 沿い)」と「C⁻ 族の交差」で不受理、交差の始まりの x が水準間で 0.03 r_t 以内。
   内部衝撃波の指標 (C⁻ 族の間隔の比の最小) が戻しの強さに単調 (9° > 7° > 0 > 5°)。
3. **亜音速入力の負例**: 初期線に M ≤ 1 の点 → `init_subsonic` (クリップ前に検査、マーチしない)。マーチ中に ν < 0
   (強い戻し) → `unit_subsonic`。
4. 不受理の理由の対応: C⁺ が壁に届かない (急拡大) → `unit_wall_miss`、修正子の上限 → `unit_iter_maxiter`
   (`interior_vec` の対の分類)、初期線の非有限・壁から外れ・順序。
5. `moc_kernel.py`・`moc_inverse.py`・`moc_diagnostics.py` を変更していない (既定はビット不変。風洞の既存試験は別に回す:
   `run_inverse_tests.py`・`run_moc_axis_limit_tests.py`)。
6. 決定性 (同じ入力で網とゲートが同一)。
7. **ガード 3 点** (§5.1 #3c、diagnostician 2026-10-11): (a) Hall の初期線は R_t ≤ 1 で n > 321 を拒否 (321 可・641 不可)、
   (b) 前線 1 (k = 1) の網の反転は壁の不受理でなく「判定不能 (初期線不整合)」(Hall R 1・n 641 の実例と、網に反転を
   入れた合成例)、(c) 不受理の根拠の点に軸からの添字 j を記録し、根拠が j ≤ 2 だけなら「判定不能 (軸近傍)」(合成例)、
   j > 2 の根拠があれば不受理のまま。既存の負例 (曲げ戻し 5°・0°) は verdict が reject。
info: ③ベルの基準ガス (EXH) と Hall の初期線 (R 1) での動作、Hall の初期線の C⁻ 適合式の不整合の大きさ。

実行: design/.venv-opt の Python で `python design/tests/run_moc_forward_tests.py`。
"""
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "design"))
from forge_design.geometry import moc_kernel as MK  # noqa: E402
from forge_design.geometry import moc_forward as MF  # noqa: E402
from forge_design.geometry.moc_forward import (HALL_N_MAX_RT_LE1, NEAR_AXIS_J, _hall_line_raw,  # noqa: E402
                                               axis_points, forward_moc, gate_record, hall_initial_line,
                                               wall_points)
from forge_design.geometry.moc_kernel import pm_mach_vec, pm_nu  # noqa: E402

FAIL = 0


def check(name, cond, detail=""):
    global FAIL
    print(("ok   " if cond else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
    if not cond:
        FAIL += 1


# ---------------------------------------------------------------- 放射源流 (厳密解)
G = 1.4
THC = np.deg2rad(10.0)
T10 = np.tan(THC)
X1 = 1.2                       # 初期線の壁点 (moc_axis_limit_radial.py の初期線の x と同じ)
L_RAD = 3.0
X_EVAL = 2.4                   # 誤差を測る範囲 x ≤ 2.4 (moc_axis_limit_radial.py の壁の標本 1.4〜2.4 と同じ上端)


def mach_of_ar(AR):
    """A/A* → 超音速の M (Newton、丸めまで)。"""
    AR = np.asarray(AR, dtype=float)
    e = (G + 1.0) / (2.0 * (G - 1.0))
    M = 2.0 + 1.5 * np.sqrt(np.maximum(AR - 1.0, 0.0))
    for _ in range(80):
        t = (2.0 / (G + 1.0)) * (1.0 + 0.5 * (G - 1.0) * M * M)
        f = t ** e / M - AR
        df = -t ** e / M ** 2 + e * t ** (e - 1.0) * (2.0 / (G + 1.0)) * (G - 1.0)
        dM = f / df
        M = M - dM
        if np.all(np.abs(dM) < 1e-15 * M):
            break
    return M


def src_M(x, r):
    return mach_of_ar(np.hypot(x, r) ** 2)


class Ray:
    def r(self, x, deriv=0):
        x = np.asarray(x, dtype=float)
        return x * T10 if deriv == 0 else np.full_like(x, T10)


class TurnWall:
    """レイ (10°) → [x1, x1+dx] で勾配を smoothstep で tan θ₁ へ → 直線。C1 (勾配連続)。"""

    def __init__(self, x1, dx, th1_deg):
        self.x1, self.dx, self.s1 = float(x1), float(dx), float(np.tan(np.deg2rad(th1_deg)))

    def r(self, x, deriv=0):
        x = np.asarray(x, dtype=float)
        t = np.clip((x - self.x1) / self.dx, 0.0, 1.0)
        S = 3.0 * t * t - 2.0 * t ** 3
        IS = t ** 3 - 0.5 * t ** 4
        slope = T10 + (self.s1 - T10) * S
        if deriv == 1:
            return slope
        r1 = self.x1 * T10
        r_in = r1 + T10 * (x - self.x1) + (self.s1 - T10) * self.dx * IS
        r_after = r1 + T10 * self.dx + 0.5 * (self.s1 - T10) * self.dx + self.s1 * (x - self.x1 - self.dx)
        return np.where(x < self.x1, x * T10, np.where(x <= self.x1 + self.dx, r_in, r_after))


def radial_line(n):
    """x = X1 の壁点から厳密場の C⁻ を軸まで (DOP853、r を独立変数)。n 点 r 等間隔、軸→壁。軸端の θ_r = 1/x (厳密)。"""
    rw = X1 * T10

    def rhs(r, y):
        M = float(src_M(y[0], r))
        return [1.0 / np.tan(np.arctan2(r, y[0]) - np.arcsin(1.0 / M))]
    sol = solve_ivp(rhs, (rw, 0.0), [X1], method="DOP853", rtol=1e-13, atol=1e-14, dense_output=True)
    rq = np.linspace(0.0, rw, n)
    xq = sol.sol(rq)[0]
    xq[-1] = X1
    return {"x": xq, "r": rq, "M": src_M(xq, rq), "theta": np.arctan2(rq, xq), "thr_axis": 1.0 / xq[0]}


def radial_errors(res):
    net = res["net"]
    X, R, TH, NU = (net.a[f] for f in ("x", "r", "th", "nu"))
    m = (net.status == 0) & (X <= X_EVAL)
    w = wall_points(res)
    a = axis_points(res)
    mw, ma = w[:, 0] <= X_EVAL, a[:, 0] <= X_EVAL
    return {"wall_nu": float(np.max(np.abs(w[mw, 3] - pm_nu(src_M(w[mw, 0], w[mw, 1]), G)))),
            "axis_nu": float(np.max(np.abs(a[ma, 3] - pm_nu(src_M(a[ma, 0], 0.0 * a[ma, 0]), G)))),
            "field_nu": float(np.max(np.abs(NU[m] - pm_nu(src_M(X[m], R[m]), G)))),
            "field_th": float(np.max(np.abs(TH[m] - np.arctan2(R[m], X[m])))),
            "mdot": res["gate"]["indicators"]["mdot_err_max"]}


print("== 1. 放射源流 (厳密解) の再現")
LEVELS = (25, 49, 97, 193, 385)
err = {"analytic": [], "legacy": []}
for n in LEVELS:
    line = radial_line(n)
    for al in ("analytic", "legacy"):
        t0 = time.time()
        res = forward_moc(Ray(), L_RAD, G, line, axis_limit=al)
        e = radial_errors(res)
        err[al].append(e)
        if al == "analytic":
            check(f"放射源流 n {n}: ゲート合格・出口まで届く", res["gate"]["pass"] and res["end"]["reason"] == "exit_reached",
                  f"{res['gate']['codes']} 前線 {res['n_fronts']}")
        print(f"     n {n:4d} {al:8s} 壁 ν {e['wall_nu']:.3e}  軸 ν {e['axis_nu']:.3e}  場 ν {e['field_nu']:.3e}  "
              f"場 θ {e['field_th']:.3e}  流量 {e['mdot']:.3e}  ({time.time() - t0:.2f} s)")
for al in ("analytic", "legacy"):
    od = {k: [np.log2(err[al][i][k] / err[al][i + 1][k]) for i in range(len(LEVELS) - 1)]
          for k in ("wall_nu", "axis_nu", "field_nu", "field_th", "mdot")}
    print(f"info 次数 ({al}): " + "  ".join(f"{k} {['%.2f' % v for v in vs]}" for k, vs in od.items()))
    if al == "analytic":
        mono = all(err[al][i][k] > err[al][i + 1][k] for i in range(len(LEVELS) - 1) for k in od)
        check("放射源流: 全ての誤差が解像度とともに単調に減る", mono)
        check("放射源流: 最細区間の次数 壁 ν ≥ 1.8", od["wall_nu"][-1] >= 1.8, f"{od['wall_nu'][-1]:.2f}")
        check("放射源流: 最細区間の次数 場 ν ≥ 1.5 (最大誤差は軸上)", od["field_nu"][-1] >= 1.5, f"{od['field_nu'][-1]:.2f}")
        check("放射源流: 最細区間の次数 流量 ≥ 1.8 (前進 MOC は流量を構成的に保存しない — 非循環の指標)",
              od["mdot"][-1] >= 1.8, f"{od['mdot'][-1]:.2f}")

# ---------------------------------------------------------------- 衝撃波を要する壁
print("== 2. 衝撃波を要する壁 (レイを x 1.6 から Δx 0.3 で曲げ戻す)")
scan = {}
for th1 in (9.0, 7.0, 5.0, 0.0):
    scan[th1] = []
    for n in (49, 97, 193):
        rec = gate_record(forward_moc(TurnWall(1.6, 0.3, th1), L_RAD, G, radial_line(n)))
        scan[th1].append(rec)
        print(f"     θ₁ {th1:4.1f}° n {n:3d}: 合否 {rec['gate_pass']}  {rec['codes']}  C⁻ 間隔比の最小 "
              f"{rec['cminus_sep_ratio_min']:.3f} @ ({rec['cminus_sep_x']:.3f}, {rec['cminus_sep_r']:.3f})  "
              f"交差の始まり x {rec['reversal_onset_x']}")
for th1 in (9.0, 7.0):
    check(f"弱い戻し θ₁ {th1}°: 3 水準とも合格", all(r["gate_pass"] for r in scan[th1]))
for th1 in (5.0, 0.0):
    recs = scan[th1]
    check(f"強い戻し θ₁ {th1}°: 3 水準とも不受理、理由に網の反転 (C⁺ 沿い) と C⁻ 族の交差",
          all((not r["gate_pass"]) and "reversal_cplus_step" in r["codes"] and "crossing_cminus" in r["codes"] for r in recs))
    check(f"強い戻し θ₁ {th1}°: verdict は reject (判定不能にしない)・根拠の j の最大 > {NEAR_AXIS_J}",
          all(r["verdict"] == "reject" and r["evidence_j_max"] > NEAR_AXIS_J for r in recs),
          f"{[(r['verdict'], r['evidence_j_min'], r['evidence_j_max']) for r in recs]}")
    xs = [r["reversal_onset_x"] for r in recs]
    check(f"強い戻し θ₁ {th1}°: 交差の始まりの x が 3 水準で 0.03 r_t 以内 (物理の包絡線)",
          max(xs) - min(xs) <= 0.03, f"{['%.4f' % v for v in xs]}")
s = {th1: scan[th1][-1]["cminus_sep_ratio_min"] for th1 in scan}
check("内部衝撃波の指標 (C⁻ 間隔比の最小、n 193) が戻しの強さに単調: 9° > 7° > 0 > 5°",
      s[9.0] > s[7.0] > 0.0 > s[5.0], f"9° {s[9.0]:.3f}・7° {s[7.0]:.3f}・5° {s[5.0]:.3f}・0° {s[0.0]:.3f}")
check("θ₁ 7° は詰まっているが交差していない (0 < 比 < 0.5、3 水準)",
      all(0.0 < r["cminus_sep_ratio_min"] < 0.5 for r in scan[7.0]), f"{[round(r['cminus_sep_ratio_min'], 3) for r in scan[7.0]]}")

# ---------------------------------------------------------------- 亜音速
print("== 3. 亜音速")
line97 = radial_line(97)
for Mbad in (0.98, 1.0):
    l2 = {k: (np.array(v, dtype=float).copy() if isinstance(v, np.ndarray) else v) for k, v in line97.items()}
    l2["M"][40] = Mbad
    res = forward_moc(Ray(), L_RAD, G, l2)
    g_ = res["gate"]
    check(f"亜音速入力 (初期線の 1 点が M {Mbad}): 理由 init_subsonic だけ・マーチしない",
          (not g_["pass"]) and g_["codes"] == ["init_subsonic"] and g_["stage"] == "initial_line" and res["net"] is None,
          g_["reasons"][0] if g_["reasons"] else "")
check("(クリップの確認) pm_nu は M < 1 を ν = 0 にする — 前進ゲートはその前に検査する",
      float(pm_nu(0.98, G)) == 0.0 and float(pm_mach_vec(np.array([-0.01]), G)[0]) <= 1.0 + 1e-9)
rec = gate_record(forward_moc(TurnWall(1.6, 0.05, -45.0), L_RAD, G, line97))
check("マーチ中の亜音速 (壁を −45° へ急に戻す): 理由に unit_subsonic、不受理", (not rec["gate_pass"]) and "unit_subsonic" in rec["codes"],
      rec["codes"])

# ---------------------------------------------------------------- 理由の対応
print("== 4. 不受理の理由の対応")
rec = gate_record(forward_moc(TurnWall(1.6, 0.02, 80.0), L_RAD, G, line97))
check("急拡大 (壁が 80° へ): C⁺ が壁に届かない unit_wall_miss・wall_not_reached",
      (not rec["gate_pass"]) and "unit_wall_miss" in rec["codes"] and "wall_not_reached" in rec["codes"], rec["codes"])
rec = gate_record(forward_moc(Ray(), L_RAD, G, line97, max_corr=1))
check("修正子の上限 1 回 (収束できない): interior_vec の分類 → unit_iter_maxiter", "unit_iter_maxiter" in rec["codes"], rec["codes"])
for name, mod, code in [("非有限", lambda l: l["theta"].__setitem__(10, np.nan), "init_nonfinite"),
                        ("壁点が壁から外れる (Δr 1e-6)", lambda l: l["r"].__setitem__(-1, l["r"][-1] + 1e-6), "init_off_wall"),
                        ("r が軸→壁の順でない", lambda l: l["r"].__setitem__(5, l["r"][3]), "init_order")]:
    l2 = {k: (np.array(v, dtype=float).copy() if isinstance(v, np.ndarray) else v) for k, v in line97.items()}
    mod(l2)
    g_ = forward_moc(Ray(), L_RAD, G, l2)["gate"]
    check(f"初期線の検査: {name} → {code}", (not g_["pass"]) and code in g_["codes"], g_["codes"])
try:
    forward_moc(Ray(), L_RAD, G, line97, corrector="fixed2")
    check("corrector='fixed2' は受けない (収束の分類を理由に使うため)", False)
except ValueError:
    check("corrector='fixed2' は受けない (収束の分類を理由に使うため)", True)

# ---------------------------------------------------------------- 既定のビット不変・決定性
print("== 5. 既定のビット不変・決定性")
r = subprocess.run(["git", "-C", str(ROOT), "diff", "--quiet", "HEAD", "--", "design/forge_design/geometry/moc_kernel.py",
                    "design/forge_design/geometry/moc_inverse.py", "design/forge_design/geometry/moc_diagnostics.py"])
check("moc_kernel.py・moc_inverse.py・moc_diagnostics.py は HEAD から無変更 (風洞の既定はビット不変)", r.returncode == 0)
check("moc_kernel の定数 (AXIS_LIMIT_FRAC 0.05・CORR_TOL 1e-12・CORR_MAX 50・RESID_TOL 1e-10) が既定のまま",
      (MK.AXIS_LIMIT_FRAC, MK.CORR_TOL, MK.CORR_MAX, MK.RESID_TOL) == (0.05, 1e-12, 50, 1e-10))
ra = forward_moc(TurnWall(1.6, 0.3, 5.0), L_RAD, G, line97)
rb = forward_moc(TurnWall(1.6, 0.3, 5.0), L_RAD, G, line97)
same = all(np.array_equal(ra["net"].a[f], rb["net"].a[f], equal_nan=True) for f in ("x", "r", "th", "nu", "M")) \
    and np.array_equal(ra["net"].status, rb["net"].status) and ra["gate"]["codes"] == rb["gate"]["codes"] \
    and ra["gate"]["reasons"] == rb["gate"]["reasons"]
check("同じ入力で網 (全点・状態) とゲートが同一", same)

# ---------------------------------------------------------------- 7. ガード 3 点 (§5.1 #3c)
print("== 7. ガード 3 点")
# (a) Hall の初期線の点数
for R_, n_, ok_ in ((1.0, 161, True), (1.0, 321, True), (1.0, 322, False), (1.0, 641, False), (0.5, 641, False),
                    (2.0, 641, True)):
    try:
        hl_ = hall_initial_line(R_, 1.4, n_)
        got = len(hl_["x"]) == n_
    except ValueError:
        got = False
    check(f"(a) Hall の初期線 R {R_:g}・n {n_}: {'受ける' if ok_ else '拒否 (ValueError)'}", got == ok_)
check(f"(a) 上限の定数 HALL_N_MAX_RT_LE1 = 321", HALL_N_MAX_RT_LE1 == 321)
check("(a) 制限のない本体 _hall_line_raw は hall_initial_line (n 321) と同じ線",
      all(np.array_equal(_hall_line_raw(1.0, 1.4, 321)[k], hall_initial_line(1.0, 1.4, 321)[k]) for k in ("x", "r", "M", "theta")))
# (b) 前線 1 の反転 → 判定不能 (初期線不整合): 実例 (CPG γ 1.4・Hall R 1・n 641 を (D) の壁で)
from forge_design.geometry.wall_bell_c2 import BellC2Spec, cone15_length, design_bell_c2  # noqa: E402
dsg = design_bell_c2(0.7 * cone15_length(50.0), 10.0, 32.5, BellC2Spec())
for n_, exp in ((641, "undecidable"), (321, None)):
    res_ = forward_moc(dsg.wall, dsg.L, 1.4, _hall_line_raw(1.0, 1.4, n_))
    g_ = res_["gate"]
    if exp == "undecidable":
        check(f"(b) 実例: γ 1.4・Hall R 1・n {n_}・(D) L 70 %・θe 10°・θmax 32.5° → 判定不能 (初期線不整合)、pass は False",
              g_["verdict"] == "undecidable" and (g_["undecidable"] or {}).get("kind") == "initial_line" and not g_["pass"],
              f"{g_['verdict']} {(g_['undecidable'] or {}).get('kind')} {g_['codes']}")
    else:
        check(f"(b) 同じ壁で n {n_} は初期線不整合にならない", (g_["undecidable"] or {}).get("kind") != "initial_line",
              f"{g_['verdict']} {g_['codes']}")


def _regate(res):
    return MF._gate(res, MF._Wall(Ray(), L_RAD), G, MF.N_CROSS_SAMPLE)


def _perturbed(n_, k_, j_, dx):
    """放射源流 (合格) の網の点 (k_, j_) を x 方向に dx ずらしてゲートだけ評価し直す (合成の負例)。"""
    res_ = forward_moc(Ray(), L_RAD, G, radial_line(n_))
    res_["net"].a["x"][k_, j_] += dx
    return res_, _regate(res_)


base_ = forward_moc(Ray(), L_RAD, G, radial_line(49))
check("(合成の前提) 放射源流 n 49 は通過、verdict pass・根拠 0 点",
      base_["gate"]["verdict"] == "pass" and base_["gate"]["evidence_j"]["n"] == 0)
dx_back = -0.05
res_b, g_b = _perturbed(49, 1, 30, dx_back)
check("(b) 合成: 前線 1 の点 (k 1・j 30) を C⁺ の担体より手前へ → 判定不能 (初期線不整合)、他の理由があっても優先",
      g_b["verdict"] == "undecidable" and g_b["undecidable"]["kind"] == "initial_line" and not g_b["pass"],
      f"{g_b['verdict']} {g_b['codes']}")
res_c, g_c = _perturbed(49, 1, 30, +0.0)
check("(b) 合成の対照: ずらさなければ通過", g_c["verdict"] == "pass")
# (c) 根拠の j の記録と軸近傍
kmid = 20
res_n, g_n = _perturbed(49, kmid, 1, dx_back)
ej = g_n["evidence_j"]
check(f"(c) 合成: 前線 {kmid} の j 1 の点を手前へ → 根拠はすべて j ≤ {NEAR_AXIS_J} → 判定不能 (軸近傍)",
      g_n["verdict"] == "undecidable" and g_n["undecidable"]["kind"] == "near_axis" and ej["j_max"] <= NEAR_AXIS_J
      and ej["n"] > 0 and ej["n_unlocated"] == 0, f"{g_n['verdict']} {g_n['codes']} j {ej}")
check("(c) 根拠の点に種類・k・j・x・r が入る", all({"kind", "k", "j", "x", "r"} <= set(e) for e in g_n["evidence"])
      and any(e["kind"] == "reversal_cplus_step" and e["k"] == kmid and e["j"] == 1 for e in g_n["evidence"]),
      f"{g_n['evidence'][:3]}")
res_f, g_f = _perturbed(49, kmid, 24, dx_back)
check("(c) 合成: 同じずらしを j 24 に → 不受理 (reject)、根拠の j の最大 > 2",
      g_f["verdict"] == "reject" and g_f["evidence_j"]["j_max"] > NEAR_AXIS_J, f"{g_f['verdict']} {g_f['evidence_j']}")
res_g, g_g = _perturbed(49, kmid, NEAR_AXIS_J + 1, dx_back)
check(f"(c) 合成: j {NEAR_AXIS_J + 1} (境界の 1 つ外) は不受理 (reject)", g_g["verdict"] == "reject",
      f"{g_g['verdict']} {g_g['evidence_j']}")
rec_ = gate_record({**res_n, "gate": g_n})
check("(c) gate_record に verdict・判定不能の種類・根拠の j・最小間隔比の j が入る",
      rec_["verdict"] == "undecidable" and rec_["undecidable_kind"] == "near_axis" and rec_["evidence_j_max"] <= NEAR_AXIS_J
      and rec_["cminus_sep_j"] is not None and "reversal_onset_j" in rec_, f"{ {k: rec_[k] for k in ('verdict', 'undecidable_kind', 'evidence_j_max', 'cminus_sep_j', 'reversal_onset_j')} }")
# 位置の分からない根拠 (前線数の上限) が 1 つでもあれば軸近傍の判定不能にしない
res_cap = forward_moc(Ray(), L_RAD, G, radial_line(49), max_fronts=5)
check("(c) 位置の分からない根拠 (前線数の上限で出口に届かない) は不受理のまま",
      res_cap["gate"]["verdict"] == "reject" and res_cap["gate"]["evidence_j"]["n_unlocated"] >= 1,
      f"{res_cap['gate']['verdict']} {res_cap['gate']['codes']} {res_cap['gate']['evidence_j']}")
# 初期線の検査で止まったときも判定不能 (初期線の検査)
l2 = {k: (np.array(v, dtype=float).copy() if isinstance(v, np.ndarray) else v) for k, v in radial_line(49).items()}
l2["M"][10] = 0.98
g_s = forward_moc(Ray(), L_RAD, G, l2)["gate"]
check("初期線の検査で止まった (亜音速入力) → verdict undecidable・種類 initial_line_check・pass False",
      g_s["verdict"] == "undecidable" and g_s["undecidable"]["kind"] == "initial_line_check" and not g_s["pass"])

# ---------------------------------------------------------------- info: ③ベルの基準ガス・Hall の初期線
print("== info: ③ベルの基準ガス (EXH) と Hall の初期線")
try:
    from forge_design.gas.rocket_exhaust import RocketExhaust
    from forge_design.geometry.wall_bell_c2 import BellC2Wall, cone15_length
    gas = RocketExhaust.from_file(ROOT / "case/68.bell_two_mode/exhaust_lox_lh2_pc100_of6.yaml").moc_gas()
    for R_ in (1.0, 2.0):
        hl = _hall_line_raw(R_, gas, 2001)        # 不整合の大きさを測る診断 (点数の制限なし)
        nu = np.asarray(pm_nu(hl["M"], gas), dtype=float)
        mu = np.arcsin(1.0 / hl["M"])
        s_ = np.where(hl["r"] > 1e-9, np.sin(hl["theta"]) / np.maximum(hl["r"], 1e-300), hl["thr_axis"])
        f_ = np.sin(mu) * s_ / np.cos(hl["theta"] - mu)
        dS = 0.5 * (f_[1:] + f_[:-1]) * (hl["x"][:-1] - hl["x"][1:])
        J = hl["theta"] + nu - np.concatenate([np.cumsum(dS[::-1])[::-1], [0.0]])
        print(f"info Hall の初期線 (R {R_}, γ* {hl['gamma_hall']:.4f}) の C⁻ 適合式の不整合 "
              f"(θ+ν−∫源項 の壁と軸端の差) {J[-1] - J[0]:+.3e} rad")
    hl = hall_initial_line(1.0, gas, 161)
    Lc = cone15_length(50.0)
    w = BellC2Wall(0.7 * Lc, 10.0, 32.5, 0.1)
    t0 = time.time()
    rec = gate_record(forward_moc(w, w.L, gas, hl))
    print(f"info EXH・Hall R 1・(D) L 70 %・θ_e 10°・θ_max 32.5°・x_m/L 0.1 (n 161): 合否 {rec['gate_pass']} {rec['codes']}・"
          f"前線 {rec['n_fronts']}・C⁻ 間隔比 {rec['cminus_sep_ratio_min']:.3f}・流量 {rec['mdot_err_max']:.2e} "
          f"({time.time() - t0:.1f} s)")
    rec = gate_record(forward_moc(w, w.L, gas, _hall_line_raw(1.0, gas, 641)))
    check("(b) 実例 (EXH): Hall R 1・n 641・同じ壁 → 判定不能 (初期線不整合) (V0-a で全点不受理 → 分類 B と誤判定された件)",
          rec["verdict"] == "undecidable" and rec["undecidable_kind"] == "initial_line", f"{rec['verdict']} {rec['codes']}")
except FileNotFoundError as e:
    print(f"info 基準ガスの排気ファイルが無いので飛ばす ({e})")

print("\nALL PASS" if FAIL == 0 else f"\n{FAIL} FAIL")
sys.exit(0 if FAIL == 0 else 1)
