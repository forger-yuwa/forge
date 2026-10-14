#!/usr/bin/env python3
"""H の固定の Euler の参照 (補正 δ_r ≡ 0 の物理壁) の準備の経路と、抽出器の δ_in と壁の入力の照合の試験
(plans/active/tooling-nozzle-isothermal-wall-chain.md §4.9 codex diagnose 2・§5.1 #36b)。

壁は case/45 の単調壁の問題 (`problem_d155_ns_finemesh_recal_final_mono.yaml`) の初期線だけ Hall に差し替えて作る (CFD ピンの凍結源 run に
依存しない; run_physical_upstream_tests.py と同じ)。H は配管 r 3.2425・L_U 9・L_pipe 3.5。NS の δ_r は physical_upstream の無い問題の
積分法 (`integral_delta_r`) で作る。CFD・変換器は使わない (格子は生成器の座標)。

1. `delta_r_zero`: 値と導関数が厳密に 0、`delta_r_from_table` と同じ属性
2. `euler_reference_wall`: キー無しは設計壁そのもの (今と同じ)。H (δ_r = 0 の物理壁) は、設計壁 S の出口端の微小な減少 (設計壁の
   ゲートは合格) が物理壁の厳密なゲートで落ち、§5.1 #36b の判断待ちの理由つきで止まる (ゲートを緩めない; 判断が出たら書き換える)。
   参照の壁の定義 (直管・S(0) への 5 次 Hermite・S; 試験の中の `HRefWall`) と NS の H の物理壁の差は x ≥ 0 で δ_r、x < 0 で
   x = 0 の δ_r の (値, 1 階, 2 階) への 5 次 Hermite
3. Euler の prepare の入口の検査: H + ramp・H で mesh_euler 無し・キー無しで single_bspline / poly の明示は run dir を作る前に例外
4. 照合の核 `compare_delta_in` の δ_in が抽出器 (`deltastar_from_core_matched_euler`) の δ_in とビット同一 (小さな合成の run)、
   `check_runs` が run の nozzle.h5 から同じ値を出す
5. H の NS の格子 (壁法線の近壁層つき) と参照の壁の Euler の格子 (mesh_euler 2000 × 97・全域 0.005・throat_refine 4・width 3) で照合:
   (i) NS の列 = Euler の列なら x ≥ 0 の |δ_in − δ_r| ≤ 1e-9 mm、(ii) 別の列なら δ_in − δ_r が Euler の壁の列間の線形補間の誤差と一致
   (独立の np.interp で)。負例: Euler の範囲が出口に届かない・Euler の壁を 2 µm ずらす・設計壁の Euler (旧い参照の形)

usage: design/.venv-opt/bin/python design/tests/run_euler_h_ref_tests.py   (所要 1〜2 分)
"""
import json
import sys
import tempfile
from pathlib import Path

import h5py
import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from forge_design.evaluate import runner_axismach as RA  # noqa: E402
from forge_design.evaluate.runner_axismach import (build_physical_wall, delta_r_from_table, delta_r_zero, design_chain,  # noqa: E402
                                                   euler_reference_wall, integral_delta_r, load_problem)
from forge_design.geometry.wall_axismach import monotone_gate  # noqa: E402
from forge_design.geometry.wall_walldriven import _hermite_quintic, _poly_eval  # noqa: E402
from forge_design.meshing.mesh2d import generate_axisym_mesh  # noqa: E402
from forge_design.metrics.delta_in_check import DELTA_IN_TOL_MM, check_runs, compare_delta_in  # noqa: E402
from forge_design.metrics.deltastar import deltastar_from_core_matched_euler  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
PROB = ROOT / "case/45.isobutane_m6_d155/problem_d155_ns_finemesh_recal_final_mono.yaml"
PU = {"r_inlet": 3.2425, "L_U": 9.0, "L_pipe": 3.5}
MESH_EULER = {"ni": 2000, "nj": 97, "throat_refine": 4.0, "throat_width": 3.0, "wall_first_frac": 0.005}   # V5d の Euler (run_0174) と同じ
FAIL = 0


def check(name, cond):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name)
    if not cond:
        FAIL += 1


def raises(fn, exc=Exception):
    try:
        fn()
    except exc as e:  # noqa: PERF203
        return f"{type(e).__name__}: {e}"
    return None


if not PROB.exists():
    print(f"問題 YAML {PROB} が無い — 全体を飛ばす (合格扱いにしない)")
    sys.exit(2)

TMP = Path(tempfile.mkdtemp(prefix="eulerhref_"))


def hall_yaml(name: str, mut=None) -> Path:
    """生産問題の写し (初期線を Hall に、pw_ramp・pw_upstream を外して既定 poly、mesh_euler を足す)。mut(raw) で書き換える。"""
    raw = yaml.safe_load(PROB.read_text())
    g = raw["geometry"]
    g["initial_line"] = "hall"
    for k in ("initial_line_run", "initial_line_res", "pw_ramp", "pw_upstream"):
        g.pop(k, None)
    raw["mesh_euler"] = dict(MESH_EULER)
    if mut is not None:
        mut(raw)
    out = TMP / f"{name}.yaml"
    out.write_text(yaml.safe_dump(raw, allow_unicode=True, sort_keys=False))
    return out


# --- 1. delta_r_zero ---------------------------------------------------------------------------------------------------
z = delta_r_zero(-12.5, 95.3)
xz = np.linspace(-13.0, 96.0, 20001)
check("delta_r_zero: 値と 1〜3 階の導関数が厳密に 0 (範囲外も)", all(np.all(z(xz, n) == 0.0) for n in range(4)))
check("delta_r_zero: delta_r_from_table と同じ属性 (supports_deriv・spline・x_range)・is_zero",
      z.supports_deriv and hasattr(z, "spline") and z.x_range == (-12.5, 95.3) and z.is_zero)

# --- 2. euler_reference_wall --------------------------------------------------------------------------------------------
yC = hall_yaml("c45_C")
yH = hall_yaml("c45_H", lambda r: r["geometry"].__setitem__("physical_upstream", dict(PU)))
yHL = hall_yaml("c45_H_legacy", lambda r: (r["geometry"].__setitem__("physical_upstream", dict(PU)),
                                           r["geometry"].__setitem__("physical_wall_repr", "legacy")))
pC, pH, pHL = load_problem(yC), load_problem(yH), load_problem(yHL)
dC = design_chain(pC)
dH = design_chain(pH)
scale = float(pC.spec["r_throat"])
S = dC["wall"]
w0, r0 = euler_reference_wall(pC, dC, scale)
check("キー無し: 設計壁そのもの (同じ物体)・記録なし (今と同じ)", w0 is dC["wall"] and r0 is None)
# H: δ_r ≡ 0 の物理壁は x ≥ 0 が設計壁 S そのもので、S の出口端の微小な減少 (設計壁のゲート monotone_gate は合格) が物理壁の厳密な
# ゲート (スロートの後で r′ ≥ −1e-12) で落ちる。どちらで判定するかは plan §5.1 #36b の判断待ちなので、ゲートを緩めずに理由つきで
# 止まることを確かめる (判断が出たらこの検査を書き換える)
mgS = monotone_gate(S._spl, S.x0, S.x_e)
e = raises(lambda: euler_reference_wall(pH, dH, scale), ValueError)
check(f"H: δ_r = 0 の物理壁は設計壁 S の出口端 (r′ 最小 {mgS['min_r1']:.2e} at x {mgS['x_min_theta']:.2f}・D {mgS['D']:.1e} r_t、"
      f"設計壁のゲート {mgS['pass']}) で物理壁の厳密なゲートに落ち、判断待ちの理由つきで止まる",
      e is not None and "#36b の判断待ち" in e and "設計壁のゲート" in e and mgS["pass"] and mgS["min_r1"] < -1e-12)


class HRefWall:
    """Euler の参照の壁の定義 (試験の中の参照の経路): 直管 r_inlet ([x_in, −L_U))・(−L_U, r_inlet, 0, 0) → x = 0 の S の (r, r′, r″) の
    5 次 Hermite ([−L_U, 0))・S (x ≥ 0)。#34 の A と同じ定義 (half_inlet_ab.py の HalfInletWall を S に当てたもの)。"""

    def __init__(self, design, r_U, L_U):
        self.d, self.r_U, self.L_U = design, float(r_U), float(L_U)
        self.x_in, self.x_e = float(design.x_in), float(design.x_e)
        self.c = _hermite_quintic(-self.L_U, 0.0, self.r_U, 0.0, 0.0, *[float(design.r(np.r_[0.0], n)[0]) for n in range(3)])

    def r(self, x, deriv: int = 0):
        x = np.asarray(x, dtype=float)
        out = np.array(self.d.r(x, deriv), dtype=float, copy=True)
        mp, mq = x < -self.L_U, (x >= -self.L_U) & (x < 0.0)
        out[mp] = self.r_U if deriv == 0 else 0.0
        if mq.any():
            out[mq] = _poly_eval(self.c, (x[mq] + self.L_U) / self.L_U, deriv) / self.L_U ** deriv
        return out


WEL = HRefWall(S, PU["r_inlet"], PU["L_U"])
xd = np.r_[0.0, np.linspace(0.0, float(S.x_e), 200001)[1:]]
xu = np.r_[np.linspace(-9.0, 0.0, 90001)[:-1], -1e-13]
# NS の H の物理壁 (build_physical_wall、区分表現) − 参照の壁: x ≥ 0 は δ_r、x < 0 は x = 0 の δ_r の (値, 1 階, 2 階) への 5 次 Hermite
# (C² でつながる)。抽出器の δ_in は x < 0 でもこの差になる (旧い参照の −248 mm の形の差は入らない)
_, drx, _ = integral_delta_r(pC, dC, pC.raw["deltastar_initializer"])
WNL = build_physical_wall(pHL, dH, scale, delta_r_x=drx, offset="radial")
gap_dn = float(np.max(np.abs((WNL.r(xd) - WEL.r(xd)) - drx(xd))))
cq = _hermite_quintic(-9.0, 0.0, 0.0, 0.0, 0.0, *[float(drx(np.r_[0.0], n)[0]) for n in range(3)])
gap_up = max(float(np.max(np.abs((WNL.r(xu, n) - WEL.r(xu, n)) - _poly_eval(cq, (xu + 9.0) / 9.0, n) / 9.0 ** n))) for n in range(3))
check(f"NS の H の物理壁 − 参照の壁: x ≥ 0 は δ_r (差 {gap_dn:.1e} ≤ 1e-14)、x < 0 は δ_r(0) への 5 次 Hermite (差 {gap_up:.1e} ≤ 1e-12)",
      gap_dn <= 1e-14 and gap_up <= 1e-12)

# --- 3. Euler の prepare の入口の検査 (run dir を作る前) -----------------------------------------------------------------------
rd = TMP / "never_euler"
cases = {
    "H + pw_upstream ramp": (hall_yaml("bad_ramp", lambda r: r["geometry"].update({"physical_upstream": dict(PU), "pw_upstream": "ramp",
                                                                                     "pw_ramp": [-11.0, -6.0]})), "poly"),
    "H で mesh_euler 無し": (hall_yaml("bad_noeu", lambda r: (r["geometry"].__setitem__("physical_upstream", dict(PU)), r.pop("mesh_euler"))),
                          "mesh_euler"),
    "H で入口の位置が違う": (hall_yaml("bad_xin", lambda r: r["geometry"].__setitem__("physical_upstream", {**PU, "L_pipe": 3.0})), "入口の位置"),
    "キー無しで physical_wall_repr: single_bspline (回帰)": (hall_yaml("bad_sb", lambda r: r["geometry"].__setitem__("physical_wall_repr",
                                                                                                                   "single_bspline")), "single_bspline"),
    "キー無しで pw_upstream: poly の明示 (回帰)": (hall_yaml("bad_poly", lambda r: r["geometry"].__setitem__("pw_upstream", "poly")), "poly"),
}
for lab, (y, key) in cases.items():
    e = raises(lambda y=y: RA.prepare(y, rd))
    check(f"Euler の prepare: {lab} → run dir を作る前に例外 ({(e or '')[:60]}…)", e is not None and key in e and not rd.exists())


# --- 4. 照合の核と抽出器の δ_in -------------------------------------------------------------------------------------------
def write_run(d: Path, X, R, q, sc: float, info_extra: dict) -> None:
    d.mkdir()
    ni, nj = X.shape
    with h5py.File(d / "nozzle.h5", "w") as f:
        f["/MESH/COORD"] = np.c_[X.ravel() * sc, R.ravel() * sc, np.zeros(ni * nj)].ravel()
    with h5py.File(d / "res_1.h5", "w") as f:
        f["/VALUE/ro"] = np.ones(ni * nj)
        f["/VALUE/Ux"] = q.ravel()
    (d / "prepare_info.json").write_text(json.dumps({"mesh": {"ni": ni, "nj": nj}, "scale_m": sc, "x_E": 10.0, **info_extra}))


sc = 0.0766
rE_f = lambda x: 1.0 + 0.03 * x * x  # noqa: E731
dl_f = lambda x: 0.01 + 0.002 * x    # noqa: E731
xN = np.sort(np.r_[np.linspace(0.0, 12.0, 37), [0.37, 5.55, 11.9]])
xEu = np.r_[np.linspace(-0.5, 6.0, 9), np.linspace(6.6, 12.5, 6)]
nj = 160
s = 1.0 - (1.0 - np.linspace(0.0, 1.0, nj)) ** 2          # 壁へ寄せた半径方向の配置 (帯の中に点が要る)
RN = (rE_f(xN) + dl_f(xN))[:, None] * s[None, :]
XN = np.repeat(xN[:, None], nj, axis=1)
qN = 1.0 - np.exp(-((rE_f(xN) + dl_f(xN))[:, None] - RN) / 0.004)   # 境界層の厚さ ~0.02 (δ_in 0.01〜0.034 と同じ桁)
RE = rE_f(xEu)[:, None] * s[None, :]
XE = np.repeat(xEu[:, None], nj, axis=1)
write_run(TMP / "syn_ns", XN, RN, qN, sc, {"throat_physical": {"r": 1.01}})
write_run(TMP / "syn_eu", XE, RE, np.ones_like(RE), sc, {})
ex = deltastar_from_core_matched_euler(TMP / "syn_ns", TMP / "syn_eu")
cmp = compare_delta_in(XN[:, 0] * sc / sc, (RN[:, -1] * sc) / sc, XE[:, 0] * sc / sc, (RE[:, -1] * sc) / sc, dl_f, sc, return_arrays=True)
A = cmp["arrays"]
idx = [int(np.where(A["x"] == xv)[0][0]) for xv in ex["x"]]
check(f"核の δ_in が抽出器の δ_in とビット同一 (抽出器が返した {len(ex['x'])} / {len(xN)} 列)",
      len(ex["x"]) >= 0.8 * len(xN) and np.array_equal(A["delta_in"][idx], ex["delta_in"]))
csv = TMP / "syn_dr.csv"
np.savetxt(csv, np.c_[np.linspace(-1.0, 13.0, 50), dl_f(np.linspace(-1.0, 13.0, 50))], delimiter=",", header="x_rt,delta_r_rt", comments="")
cr = check_runs(TMP / "syn_ns", TMP / "syn_eu", csv)
# check_runs は δ_r を CSV から delta_r_from_table で読む (直線の表の 5 次補間なので丸めの差だけ)
check(f"check_runs: run の nozzle.h5 から核と同じ最大差 ({cr['judged']['max_abs_rt']:.6e} ≈ {cmp['judged']['max_abs_rt']:.6e})・csv の sha256",
      abs(cr["judged"]["max_abs_rt"] - cmp["judged"]["max_abs_rt"]) <= 1e-13 and len(cr["delta_r_csv_sha256"]) == 64 and cr["scale_same"]
      and cr["covered"])

# --- 5. H の NS と H の Euler の格子で照合 -------------------------------------------------------------------------------------
# 抽出器は Euler の壁を列の間で線形に補間するので、NS の列が Euler の列と別の x にあると δ_in に S″·h²/8 の補間の誤差が入る。
# (i) NS の列 = Euler の列 (同じ配置): 補間を通らず δ_in = δ_r (丸めだけ)。(ii) 別の配置 (ni 3001・throat_refine 3・width 1.5):
# δ_in − δ_r が Euler の壁の列間の線形補間の誤差と一致することを、独立に計算した np.interp(x, xE, r_E(xE)) − r_E(x) で確かめる
# (合否は tol との比較のまま。この配置での値は記録)。どちらも NS は壁法線の近壁層つき
mpE = RA.mesh_params_euler(pH, scale)
cE, _, _ = generate_axisym_mesh(WEL, mpE)
XEs, REs = cE[:, 0].reshape(mpE.ni, -1) / scale, cE[:, 1].reshape(mpE.ni, -1) / scale


def ns_cols(mesh_upd):
    y = hall_yaml("c45_H_ns_" + "_".join(map(str, mesh_upd.values())).replace(".", "p"),
                  lambda r: (r["geometry"].__setitem__("physical_upstream", dict(PU)), r["mesh"].update(mesh_upd)))
    pn = load_problem(y)
    W = build_physical_wall(pn, dH, scale, delta_r_x=drx, offset="radial")
    mp = RA.mesh_params(pn, scale, ni=561, nj=97, wall_first_frac=4.5e-5)
    c, _, _ = generate_axisym_mesh(W, mp)
    X, R = c[:, 0].reshape(mp.ni, -1) / scale, c[:, 1].reshape(mp.ni, -1) / scale
    return X, R, W


XNa, RNa, WNa = ns_cols({"wall_normal_layer": [0.02, 0.8]})                                   # ni 2000・throat_refine 4・width 3 (= Euler)
XNs, RNs, WN = ns_cols({"ni": 3001, "throat_refine": 3.0, "throat_width": 1.5, "wall_normal_layer": [0.02, 0.8]})
check("NS の格子 (壁法線の近壁層): 壁の節点の x = 軸の節点の x (列の x)。(i) は Euler と同じ列、(ii) は別の列",
      np.array_equal(XNs[:, -1], XNs[:, 0]) and np.array_equal(XNa[:, -1], XNa[:, 0]) and np.array_equal(XNa[:, 0], XEs[:, 0])
      and len(np.intersect1d(XNs[:, 0], XEs[:, 0])) <= 2)
ra = compare_delta_in(XNa[:, 0], RNa[:, -1], XEs[:, 0], REs[:, -1], drx, scale, r_ns_exact=WNa.r, r_eu_exact=WEL.r)
check(f"(i) 同じ列: x ≥ 0 の max|δ_in − δ_r| = {ra['judged']['max_abs_mm']:.1e} mm ≤ 1e-9 mm・全列の被覆・合格", ra["pass"]
      and ra["judged"]["max_abs_mm"] <= 1e-9)
res = compare_delta_in(XNs[:, 0], RNs[:, -1], XEs[:, 0], REs[:, -1], drx, scale, r_ns_exact=WN.r, r_eu_exact=WEL.r, return_arrays=True)
j, dec = res["judged"], res["decomposition"]
xs_, sel = res["arrays"]["x"], res["arrays"]["x"] >= 0.0
lin = np.interp(xs_, XEs[:, 0], WEL.r(XEs[:, 0])) - WEL.r(xs_)     # 独立に計算した線形補間の誤差 (x = 0 の直前の Euler の列は Q の上)
dev = float(np.max(np.abs((res["arrays"]["delta_in"] - res["arrays"]["delta_r"])[sel] + lin[sel])))
check(f"(ii) 別の列: δ_in − δ_r = −(Euler の壁の列間の線形補間の誤差) (差 {dev:.1e} r_t ≤ 1e-12)。この配置の max {j['max_abs_mm']:.2e} mm "
      f"(x {j['x_at_max']:.3f}) → {'合格' if res['pass'] else '不合格'} (tol {DELTA_IN_TOL_MM} mm; 内訳: Euler の補間 "
      f"{dec['euler_interp_minus_exact']['max_abs_mm']:.1e} mm・NS の節点 {dec['ns_nodes_minus_exact']['max_abs_mm']:.1e} mm・真の間隔 − δ_r "
      f"{dec['exact_gap_minus_delta_r_x_ge_min']['max_abs_mm']:.1e} mm)",
      dev <= 1e-12 and res["covered"] and res["pass"] == (j["max_abs_mm"] <= DELTA_IN_TOL_MM))
# 負例 ((i) の同じ列の配置で)
cut = XEs[:, 0] <= float(S.x_e) - 0.5
r1 = compare_delta_in(XNa[:, 0], RNa[:, -1], XEs[cut, 0], REs[cut, -1], drx, scale)
check(f"負例: Euler の範囲が出口に 0.5 r_t 届かない → 不合格 (外れ {r1['n_uncovered']} 列)", not r1["pass"] and r1["n_uncovered"] > 0)
shift = 2e-3 / 1e3 / scale
r2 = compare_delta_in(XNa[:, 0], RNa[:, -1], XEs[:, 0], REs[:, -1] + shift * (XEs[:, 0] >= 40.0), drx, scale)
check(f"負例: Euler の壁を試験部 (x ≥ 40) だけ 2 µm 外へ → 不合格 (max {r2['judged']['max_abs_mm']:.2e} mm)", not r2["pass"])
cD, _, _ = generate_axisym_mesh(S, mpE)
XDs, RDs = cD[:, 0].reshape(mpE.ni, -1) / scale, cD[:, 1].reshape(mpE.ni, -1) / scale
r3 = compare_delta_in(XNa[:, 0], RNa[:, -1], XDs[:, 0], RDs[:, -1], drx, scale)
up = r3["upstream_x_lt_min"]
check(f"設計壁の Euler (旧い参照の形): x ≥ 0 は同じ判定 ({r3['judged']['max_abs_mm']:.1e} mm)、x < 0 の δ_in は H と設計の縮流部の差 "
      f"(max {up['max_abs_mm']:.1f} mm at x {up['x_at_max']:.2f}) — 判定しない区間に出る",
      r3["pass"] and up["max_abs_mm"] > 100.0)

print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAIL'}")
sys.exit(1 if FAIL else 0)
