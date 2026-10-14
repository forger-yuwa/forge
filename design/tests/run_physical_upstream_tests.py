#!/usr/bin/env python3
"""物理壁だけの上流の形 (`geometry.physical_upstream`) の試験 (plans/active/tooling-nozzle-isothermal-wall-chain.md §4.9・§5.1 #34)。

壁は case/45 の単調壁の生産問題 (`problem_d155_ns_finemesh_recal_final_mono.yaml`) の初期線だけ Hall に差し替えて作る (CFD ピンの
凍結源 run に依存しない; run_pw_upstream_poly_tests.py と同じ)。δ_r は physical_upstream の無い問題の積分法 (`integral_delta_r`) で作り、
両方の壁に同じものを渡す。H は core-grid plan §4.98 の形 (配管 r 3.2425・L_U 9・L_pipe 3.5、入口の位置 −12.5 は同じ)。

1. キー: `_physical_upstream` (キー無し → None・正しい値・不正値・欠け・未知のキー・真偽値・非有限・r_inlet ≤ 1・joint でない壁・
   pw_upstream ramp・入口の位置の違い・μ > 20)、`PhysicalNozzleWall` の直接の引数の検査
2. 回帰: キー無しの壁は `physical_upstream` None・`up_ref` が設計壁の Hermite・壁ファイルに physical_upstream が無い・ゲートに H_ref が無い
3. H の壁 (区分表現と 1 本の B-spline): x ≥ 0 が C (キー無し) と同じ、直管 = r_inlet、Q の端条件 (−L_U で (r_U, 0, 0)、0 で S + δ_r)、
   既往の定義 (forge-25 の cg_hi_prep.py: (−9, 3.2425, 0, 0) → x = 0 の C の (r, r′, r″) の 5 次 Hermite) と 1e-12 以内、
   ゲート (|Q″ − H″| の H は物理壁の上流の名目の Hermite)、μ、凍結した設計の記録、validate
4. 壁ファイル: legacy と single_bspline の往復 (legacy は physical_upstream で直管・Q の区間を復元)、壊れた記録は例外
5. チェーンの関数: 積分法 (`integral_delta_r`) と prepare_ns (delta_r_csv 無し・Tw_profile) は run dir を作る前に例外、
   Euler の prepare (δ_r = 0 の物理壁の経路、§5.1 #36b; 詳しくは run_euler_h_ref_tests.py) は mesh_euler の無い問題で run dir を作る前に止まる、
   design_chain は設計壁の値のまま (physical_upstream で設計が変わらない)

usage: design/.venv-opt/bin/python design/tests/run_physical_upstream_tests.py   (所要 1 分程度)
"""
import json
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from forge_design.evaluate import runner_axismach as RA  # noqa: E402
from forge_design.evaluate.runner_axismach import (_gam_or_gas, _physical_upstream, build_physical_wall, design_chain,  # noqa: E402
                                                   integral_delta_r, load_problem)
from forge_design.geometry.wall_axismach import (PhysicalNozzleWall, SingleBSplinePhysicalWall, check_required_attrs,  # noqa: E402
                                                 load_wall_file, save_wall_file, wall_file_record)
from forge_design.geometry.wall_walldriven import _hermite_quintic, _poly_eval  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
PROB = ROOT / "case/45.isobutane_m6_d155/problem_d155_ns_finemesh_recal_final_mono.yaml"
FAIL = 0
PU = {"r_inlet": 3.2425, "L_U": 9.0, "L_pipe": 3.5}


def check(name, cond):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name)
    if not cond:
        FAIL += 1


def raises(fn, exc=ValueError):
    try:
        fn()
    except exc as e:  # noqa: PERF203
        return str(e)
    return None


if not PROB.exists():
    print(f"問題 YAML {PROB} が無い — 全体を飛ばす (合格扱いにしない)")
    sys.exit(2)

TMP = Path(tempfile.mkdtemp(prefix="physup_"))


def hall_yaml(name: str, mut=None) -> Path:
    """生産問題の写し (初期線を Hall に、pw_ramp・pw_upstream を外して既定 poly)。mut(g, spec) で geometry・spec を書き換える。"""
    import yaml
    raw = yaml.safe_load(PROB.read_text())
    g = raw["geometry"]
    g["initial_line"] = "hall"
    for k in ("initial_line_run", "initial_line_res", "pw_ramp", "pw_upstream"):
        g.pop(k, None)
    if mut is not None:
        mut(g, raw["spec"])
    out = TMP / f"{name}.yaml"
    out.write_text(yaml.safe_dump(raw, allow_unicode=True, sort_keys=False))
    return out


# --- 1. キー ---------------------------------------------------------------------------------------------------------
J = {"wall_repr": "joint", "r_inlet": 6.485, "L_U": 12.0, "L_pipe": 0.5, "R": 2.0}
check("キー無し → None (今と同じ)", _physical_upstream(J) is None)
r_ = _physical_upstream({**J, "physical_upstream": dict(PU)})
check(f"正しい値 → 解決済みの値と凍結した設計の値 (μ {r_['mu']:.4f}、x_in {r_['x_in']})",
      r_["r_inlet"] == 3.2425 and r_["L_U"] == 9.0 and r_["L_pipe"] == 3.5 and r_["x_in"] == -12.5
      and r_["frozen_design"] == {"r_inlet": 6.485, "L_U": 12.0, "L_pipe": 0.5} and abs(r_["mu"] - 81.0 / (2.0 * 2.2425)) <= 1e-12)
check("整数の値も受ける (YAML の 9 など)", _physical_upstream({**J, "physical_upstream": {"r_inlet": 3.2425, "L_U": 9, "L_pipe": 3.5}})["L_U"] == 9.0)
bad = {"mapping でない": [3.2425, 9.0, 3.5], "null": None, "欠け": {"r_inlet": 3.2425, "L_U": 9.0},
       "未知のキー": {**PU, "r_U": 3.2425}, "真偽値": {**PU, "L_pipe": True}, "文字列": {**PU, "L_U": "9.0"}, "非有限": {**PU, "r_inlet": float("nan")},
       "r_inlet ≤ 1": {**PU, "r_inlet": 1.0}, "L_U ≤ 0": {"r_inlet": 3.2425, "L_U": 0.0, "L_pipe": 12.5},
       "L_pipe ≤ 0": {"r_inlet": 3.2425, "L_U": 12.5, "L_pipe": 0.0}, "入口の位置が違う": {**PU, "L_pipe": 3.0},
       "μ > 20 (L_U² / (R Δr))": {"r_inlet": 2.0, "L_U": 9.0, "L_pipe": 3.5}}
for lab, v in bad.items():
    e = raises(lambda v=v: _physical_upstream({**J, "physical_upstream": v}))
    check(f"不正値 ({lab}) → 例外 ({(e or '')[:60]}…)", e is not None)
e = raises(lambda: _physical_upstream({**J, "wall_repr": "interp", "physical_upstream": dict(PU)}))
check(f"joint でない壁 → 例外 ({(e or '')[:40]}…)", e is not None and "joint" in e)
e = raises(lambda: _physical_upstream({**J, "pw_upstream": "ramp", "pw_ramp": [-11.0, -6.0], "physical_upstream": dict(PU)}))
check(f"pw_upstream ramp → 例外 ({(e or '')[:40]}…)", e is not None and "poly" in e)

# --- 2・3. 壁 ------------------------------------------------------------------------------------------------------
yC = hall_yaml("c45_C")
yH = hall_yaml("c45_H", lambda g, s: g.__setitem__("physical_upstream", dict(PU)))
pC, pH = load_problem(yC), load_problem(yH)
dC = design_chain(pC)
dH = design_chain(pH)
xq = np.linspace(-12.5, float(dC["wall"].x_e), 40001)
check("design_chain は physical_upstream で変わらない (設計壁 r、r′、r″ がビット同一・縮流部は凍結した設計の値)",
      all(np.array_equal(dC["wall"].r(xq, n), dH["wall"].r(xq, n)) for n in range(3)) and float(dH["wall"].up.r_U) == 6.485
      and float(dH["wall"].up.L_U) == 12.0)
_, drx, _ = integral_delta_r(pC, dC, pC.raw["deltastar_initializer"])
scale = float(pC.spec["r_throat"])
args = (dC["wall"], dC["wall_inv"], scale, float(pC.spec["Pt"]), float(pC.spec["Tt"]), _gam_or_gas(pC), pC.cp)
PC = PhysicalNozzleWall(*args, offset="radial", delta_r_x=drx)
PH = PhysicalNozzleWall(*args, offset="radial", delta_r_x=drx, physical_upstream={"r_U": 3.2425, "L_U": 9.0, "L_pipe": 3.5})
check("回帰: キー無しの壁は physical_upstream None・up_ref は設計壁の Hermite・ゲートに H_ref が無い・既定ランプとの差の診断あり",
      PC.physical_upstream is None and PC.up_ref is dC["wall"].up and "H_ref" not in PC.upstream_gate
      and "max_abs_diff_rt" in PC.upstream_gate["vs_ramp_default"])
check("回帰: build_physical_wall (キー無し) は直接の PhysicalNozzleWall と同じ壁 (1 本の B-spline の元の区分表現と r がビット同一)",
      np.array_equal(SingleBSplinePhysicalWall(PC).r(xq), build_physical_wall(pC, dC, scale, delta_r_x=drx, offset="radial").r(xq)))
xd = np.linspace(0.0, PC.x_e, 200001)
check("H: x ≥ 0 が C とビット同一 (r, r′, r″, r‴)", all(np.array_equal(PH.r(xd, n), PC.r(xd, n)) for n in range(4)))
xp = np.linspace(-12.5, -9.0, 2001)[:-1]
check("H: 直管 [−12.5, −9) は r = 3.2425、r′ = r″ = 0", np.all(PH.r(xp) == 3.2425) and np.all(PH.r(xp, 1) == 0.0) and np.all(PH.r(xp, 2) == 0.0))
check("H: x_in・r_U・L_U・_herm_x0 は物理壁の上流の値", PH.x_in == -12.5 and PH.r_U == 3.2425 and PH.L_U == 9.0 and PH._herm_x0 == -9.0)
e0 = [float(PC.r(np.r_[0.0], n)[0]) for n in range(3)]
c_ref = _hermite_quintic(-9.0, 0.0, 3.2425, 0.0, 0.0, *e0)        # 既往の定義 (forge-25 cg_hi_prep.py の HalfInletWall)
xu = np.r_[np.linspace(-9.0, 0.0, 90001)[:-1], -1e-13]
ref = [_poly_eval(c_ref, (xu + 9.0) / 9.0, n) / 9.0 ** n for n in range(3)]
dref = [float(np.max(np.abs(PH.r(xu, n) - ref[n]))) for n in range(3)]
check(f"H: [−9, 0) が既往の H の定義 (x = 0 の C の (r, r′, r″) への 5 次 Hermite) と一致 (差 {dref[0]:.1e}・{dref[1]:.1e}・{dref[2]:.1e} ≤ 1e-12)",
      max(dref) <= 1e-12)
g = PH.upstream_gate
check(f"H: ゲート合格・継ぎ目の跳び {g['max_seam_jump']:.1e} ≤ 1e-8・H_ref は物理壁の上流の名目の Hermite",
      g["pass"] and g["max_seam_jump"] <= 1e-8 and "physical_upstream" in g["H_ref"] and "status" in g["vs_ramp_default"])
from forge_design.geometry.wall_walldriven import UpstreamThroatPoly  # noqa: E402
Hn = UpstreamThroatPoly(r_U=3.2425, R_t=2.0, L_U=9.0)
d2n = float(np.max(np.abs(PH.r(xu, 2) - Hn.r(xu, 2))))
check(f"H: |Q″ − H″| は名目の Hermite (3.2425, L_U 9) に対して ({g['max_abs_d2_change_vs_H']:.6e} ≥ 密な点 {d2n:.6e}、差 ≤ 1e-9)",
      g["max_abs_d2_change_vs_H"] >= d2n - 1e-15 and g["max_abs_d2_change_vs_H"] - d2n <= 1e-9)
pu = PH.physical_upstream
check(f"H: 記録 (μ {pu['mu']:.4f}、凍結した設計 {pu['frozen_design']})",
      abs(pu["mu"] - Hn.mu) <= 1e-15 and pu["frozen_design"] == {"r_U": 6.485, "L_U": 12.0, "L_pipe": 0.5, "x_in": -12.5}
      and pu["x_contraction"] == -9.0)
check(f"H: validate() が空 ({PH.validate()})・必須属性", PH.validate() == [] and check_required_attrs(PH) is None)
xs_m = np.linspace(-12.5, 0.0, 30001)
check("H: 縮流部で単調に縮む (r′ ≤ 0、x < −0.01)", float(np.max(PH.r(xs_m[xs_m < -0.01], 1))) <= 0.0)
check(f"H: 物理スロートは大域最小 (x_t {PH.x_throat:.6f}、r_t {PH.r_throat:.7f}; C は x_t {PC.x_throat:.6f})",
      abs(PH.r_throat - float(np.min(PH.r(np.linspace(-1, 1, 200001))))) <= 1e-9)
SH = SingleBSplinePhysicalWall(PH)
dS = [float(np.max(np.abs(SH.r(xq, n) - PH.r(xq, n)))) for n in range(3)]
check(f"H (1 本の B-spline): 区分表現との差 {dS[0]:.1e}・{dS[1]:.1e}・{dS[2]:.1e} (≤ 1e-9・1e-8・1e-6)、ゲート合格・H_ref・記録を写す",
      dS[0] <= 1e-9 and dS[1] <= 1e-8 and dS[2] <= 1e-6 and SH.upstream_gate["pass"] and "H_ref" in SH.upstream_gate
      and SH.physical_upstream is pu and SH.up_ref is PH.up_ref and SH.validate() == [] and SH.x_in == -12.5)
WB = build_physical_wall(pH, dC, scale, delta_r_x=drx, offset="radial")
check("build_physical_wall (H の問題) は 1 本の B-spline の H の壁 (直接の構築とビット同一)",
      isinstance(WB, SingleBSplinePhysicalWall) and np.array_equal(WB.r(xq), SH.r(xq)) and WB.physical_upstream["r_U"] == 3.2425)
e = raises(lambda: PhysicalNozzleWall(*args, offset="radial", delta_r_x=drx, physical_upstream={"r_U": 3.2425, "L_U": 9.0, "L_pipe": 3.5},
                                      ramp=(-11.0, -6.0), upstream="ramp"))
check(f"PhysicalNozzleWall: ramp + physical_upstream → 例外 ({(e or '')[:40]}…)", e is not None and "poly" in e)
e = raises(lambda: PhysicalNozzleWall(*args, offset="radial", delta_r_x=drx, physical_upstream={"r_U": 3.2425, "L_U": 9.0}))
check("PhysicalNozzleWall: キーの欠け → 例外", e is not None)
e = raises(lambda: PhysicalNozzleWall(*args, offset="radial", delta_r_x=drx, physical_upstream={"r_U": 3.2425, "L_U": 9.0, "L_pipe": 3.0}))
check("PhysicalNozzleWall: 入口の位置の違い → 例外", e is not None and "入口の位置" in e)

def _mu_bad():
    W = PhysicalNozzleWall(*args, offset="radial", delta_r_x=drx, physical_upstream={"r_U": 1.5, "L_U": 5.0, "L_pipe": 7.5})
    return W.validate()


e = raises(_mu_bad)
if e is None:
    msgs = _mu_bad()
    check(f"PhysicalNozzleWall: μ > 20 の上流は validate が不合格 ({[m[:50] for m in msgs][:1]})", any("μ" in m for m in msgs))
else:
    check(f"PhysicalNozzleWall: μ > 20 の上流は構築時の上流のゲートで不合格 ({e[:60]}…)", "前で r′ ≤ 0 False" in e)

# --- 4. 壁ファイル -------------------------------------------------------------------------------------------------
for kind, W in (("legacy", PH), ("single_bspline", SH)):
    d_ = TMP / f"wf_{kind}"
    d_.mkdir()
    path, rec, _ = save_wall_file(d_, W, scale, kind)
    back = load_wall_file(path)
    check(f"壁ファイル ({kind}): 往復で r がビット同一・physical_upstream を記録・設計壁は凍結した値 (r_U {rec['design_wall']['pipe']['r_U']})",
          np.array_equal(back["physical"].r(xq), W.r(xq)) and rec["physical_wall"]["physical_upstream"]["r_U"] == 3.2425
          and rec["design_wall"]["pipe"]["r_U"] == 6.485 and rec["domain"] == [-12.5, float(W.x_e)])
recC = wall_file_record(SingleBSplinePhysicalWall(PC), scale, "single_bspline")
check("回帰: キー無しの壁の壁ファイルに physical_upstream が無い", "physical_upstream" not in recC["physical_wall"])
rec_bad = json.loads((TMP / "wf_legacy" / "wall_repr.json").read_text())
rec_bad["physical_wall"]["physical_upstream"].pop("L_U")
(TMP / "wf_bad").mkdir()
(TMP / "wf_bad" / "wall_repr.json").write_text(json.dumps(rec_bad))
e = raises(lambda: load_wall_file(TMP / "wf_bad"))
check(f"壁ファイル (legacy): physical_upstream の欠けは例外 ({(e or '')[:40]}…)", e is not None)
rec_nopu = json.loads((TMP / "wf_legacy" / "wall_repr.json").read_text())
rec_nopu["physical_wall"].pop("physical_upstream")
(TMP / "wf_nopu").mkdir()
(TMP / "wf_nopu" / "wall_repr.json").write_text(json.dumps(rec_nopu))
e = raises(lambda: load_wall_file(TMP / "wf_nopu"))
check(f"壁ファイル (legacy): physical_upstream を消すと Q の区間が設計壁と合わず例外 (黙って別の壁にしない) ({(e or '')[:40]}…)", e is not None)

# --- 5. チェーンの関数 -----------------------------------------------------------------------------------------------
e = raises(lambda: integral_delta_r(pH, dH, pH.raw["deltastar_initializer"]))
check(f"積分法の初期壁 + physical_upstream → 例外 ({(e or '')[:40]}…)", e is not None and "physical_upstream" in e)
rd = TMP / "never_ns"
e = raises(lambda: RA.prepare_ns(yH, rd))
check(f"prepare_ns: physical_upstream で delta_r_csv が無い → run dir を作る前に例外 ({(e or '')[:40]}…)",
      e is not None and "physical_upstream" in e and not rd.exists())
csv = TMP / "dr.csv"
xt = np.linspace(-12.5 + 1e-6, float(dC["wall"].x_e), 1500)
np.savetxt(csv, np.c_[xt, drx(xt)], delimiter=",", header="x_rt,delta_r_rt", comments="")
yT = hall_yaml("c45_H_tw", lambda g, s: (g.__setitem__("physical_upstream", dict(PU)),
                                          s.__setitem__("wall_thermal", {"mode": "isothermal", "Tw_profile": {"s": [-2.0, 0.0, 1.0],
                                                                                                             "Tw_K": [300.0, 300.0, 300.0]}})))
e = raises(lambda: RA.prepare_ns(yT, rd, delta_r_csv=str(csv)))
check(f"prepare_ns: physical_upstream + Tw_profile → run dir を作る前に例外 ({(e or '')[:50]}…)",
      e is not None and "physical_upstream" in e and "Tw_profile" in e and not rd.exists())
yB = hall_yaml("c45_H_bad", lambda g, s: g.__setitem__("physical_upstream", {**PU, "L_pipe": 3.0}))
e = raises(lambda: RA.prepare_ns(yB, rd, delta_r_csv=str(csv)))
check(f"prepare_ns: 不正な physical_upstream → run dir を作る前に例外 ({(e or '')[:40]}…)",
      e is not None and "physical_upstream" in e and not rd.exists())
e = raises(lambda: RA.prepare(yH, TMP / "never_euler"))
# 2026-10-14 (§5.1 #36b): Euler の prepare は physical_upstream の問題を δ_r = 0 の物理壁で準備する経路に変えた (試験は
# run_euler_h_ref_tests.py)。この問題には mesh_euler が無いので、physical_upstream の検査を通った後、run dir を作る前に mesh_euler で止まる
check(f"Euler の prepare: physical_upstream の問題 (mesh_euler 無し) は run dir を作る前に mesh_euler で止まる ({(e or '')[:40]}…)",
      e is not None and "mesh_euler" in e and not (TMP / "never_euler").exists())

print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAIL'}")
sys.exit(1 if FAIL else 0)
