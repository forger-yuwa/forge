#!/usr/bin/env python3
"""CFD ピン抽出器 (`feedback/cfd_initial_line.py`) の単体試験 = plan tooling-nozzle-cfd-pinned-initial-line §6 V0 (事前登録)。

case/45 の Euler 格子に Hall 場を載せた**合成 Hall 場**から、Hall の線とアンカーを再現できるかを見る:
  線 max|ΔM| ≤ 2e-4・max|Δθ| ≤ 0.004°・|Δx₀| ≤ 2e-4、アンカー |ΔM| ≤ 2e-5・|ΔM′| ≤ 1e-4・|ΔM″| ≤ 1e-3、
  m* の再現 ≤ 1e-4 (相対)、ν↔M 往復 ≤ 2e-4。
加えて入力の拒否 (亜音速の場・アンカー位置の取り違え) と design_chain のキー検査。
(試作 case/45.isobutane_m6_d155/test_cfd_initial_line.py の移植。)

usage: design/.venv-opt/bin/python design/tests/run_cfd_initial_line_tests.py [格子の run]
  格子の run の既定: $FORGE_CFDPIN_GRID_RUN か /home/sano/work/forge/case/45.isobutane_m6_d155/run_0062_euler_wallfit_fit_r1_ext6k
  (検証 run は git 管理外。無ければ exit 2 で止める — 合格扱いにしない)
"""
import json
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "design"))

from forge_design.evaluate.runner_axismach import _gam_or_gas, design_chain, load_problem  # noqa: E402
from forge_design.feedback.cfd_initial_line import CFDPinnedThroat  # noqa: E402
from forge_design.geometry.moc_inverse import _flux_along, _Pt  # noqa: E402
from forge_design.geometry.moc_kernel import pm_mach, pm_nu  # noqa: E402
from forge_design.geometry.transonic import HallThroat  # noqa: E402

FAIL = 0


def check(name, cond):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name)
    if not cond:
        FAIL += 1


PROB = ROOT / "case/45.isobutane_m6_d155/problem_d155_euler_c2final_n2400.yaml"
run = Path(sys.argv[1] if len(sys.argv) > 1 else os.environ.get(
    "FORGE_CFDPIN_GRID_RUN", "/home/sano/work/forge/case/45.isobutane_m6_d155/run_0062_euler_wallfit_fit_r1_ext6k"))
if not (run / "nozzle.h5").exists():
    print(f"格子の run が無い: {run} (引数か FORGE_CFDPIN_GRID_RUN で与える)")
    sys.exit(2)

import h5py  # noqa: E402

info = json.loads((run / "prepare_info.json").read_text())
S = float(info["scale_m"])
ni = int(info["mesh"]["ni"])
with h5py.File(run / "nozzle.h5") as f:
    nc = f["/MESH/COORD"][:].reshape(-1, 3)
X = (nc[:, 0] / S).reshape(ni, -1)
R = (nc[:, 1] / S).reshape(ni, -1)

# --- V0: 合成 Hall 場 ------------------------------------------------------------
g = 1.2735422978978068
ht = HallThroat(R=2.0, gamma=g)
m = (X[:, 0] > -2.2) & (X[:, 0] < 2.7)
Mh = np.ones(X.shape)
Th = np.zeros(X.shape)
Mh[m] = ht.mach(X[m], R[m])
Th[m] = ht.theta(X[m], R[m])
cfd = CFDPinnedThroat(2.0, g, X, R, Mh, Th)
xc, rc, Mc, tc = cfd.throat_characteristic(n=41)
xh, rh, Mhl, thl = ht.throat_characteristic(n=41)
aC = cfd.axis_anchor(cfd.x0_cfd)
aH = ht.axis_anchor(cfd.x0_cfd)
p = load_problem(PROB)
gas = _gam_or_gas(p)
gg = gas if hasattr(gas, "nu") else g


def flux(x, r, M, t):
    return _flux_along([_Pt(float(a), float(b), float(d), float(pm_nu(float(c), gg)), gg)
                        for a, b, c, d in zip(x, r, M, t)], gg)[-1]


mC, mH = flux(xc, rc, Mc, tc), flux(xh, rh, Mhl, thl)
Ms = np.linspace(1.0005, 1.3, 200)
rt = np.array([float(pm_mach(float(pm_nu(float(v), gg)), gg)) for v in Ms])
out = dict(line_dM=float(np.abs(Mc - Mhl).max()), line_dtheta_deg=float(np.degrees(np.abs(tc - thl)).max()),
           dx0=float(xc[0] - xh[0]), anchor_dM=aC[0] - aH[0], anchor_dMp=aC[1] - aH[1], anchor_dMpp=aC[2] - aH[2],
           mstar_rel=float(mC / mH - 1), nu_roundtrip=float(np.abs(rt - Ms).max()))
gates = dict(line_dM=2e-4, line_dtheta_deg=0.004, dx0=2e-4, anchor_dM=2e-5, anchor_dMp=1e-4, anchor_dMpp=1e-3,
             mstar_rel=1e-4, nu_roundtrip=2e-4)
for k, v in gates.items():
    check(f"V0 {k}: {out[k]:.3e} (≤ {v:g})", abs(out[k]) <= v)

# --- 入力の拒否 --------------------------------------------------------------------
try:
    CFDPinnedThroat(2.0, g, X, R, np.full(X.shape, 0.9), np.zeros(X.shape))
    check("一様 M=0.9 (亜音速) の場を拒否", False)
except ValueError as e:
    check(f"一様 M=0.9 (亜音速) の場を拒否 ({e})", "亜音速" in str(e))
try:
    cfd.axis_anchor(cfd.x0_cfd + 0.01)
    check("軸着地 x0 以外でのアンカー要求を拒否", False)
except ValueError:
    check("軸着地 x0 以外でのアンカー要求を拒否", True)

# --- design_chain のキー -----------------------------------------------------------
for key, val in (("initial_line", "spline"), ("wall_repr", "smooth")):
    pp = load_problem(PROB)
    pp.geometry[key] = val
    try:
        design_chain(pp)
        check(f"geometry.{key}: {val} を拒否", False)
    except ValueError:
        check(f"geometry.{key}: {val} を拒否", True)
pp = load_problem(PROB)
pp.geometry["initial_line"] = "cfd"
try:
    design_chain(pp)
    check("initial_line: cfd で initial_line_run 無しを拒否", False)
except ValueError:
    check("initial_line: cfd で initial_line_run 無しを拒否", True)
pp = load_problem(PROB)
pp.geometry.update({"initial_line": "cfd", "initial_line_res": "res_6000.h5", "Md_moc_offset": -4.16e-4,
                    "initial_line_run": os.path.relpath(run, PROB.parent)})      # problem ファイルからの相対パス
if (run / "res_6000.h5").exists():
    d = design_chain(pp)
    il = d["initial_line"]
    check(f"相対パスの initial_line_run を解決 ({il['run']})", Path(il["run"]).resolve() == run.resolve())
    check(f"報告の Md は spec のまま ({d['Md']}), MOC は較正値 ({d['Md_moc']})",
          d["Md"] == float(pp.spec["M_design"]) and d["Md_moc"] == float(pp.spec["M_design"]) - 4.16e-4)
    check(f"x0 = CFD 線の軸着地 ({il['x0']:.6f})", abs(il["x0"] - d["x_A"]) == 0.0 and il["source"] == "cfd")
    check(f"アンカーの出所を成分別に記録 ({d['anchor_source']})",
          d["anchor_source"] == {"M": "cfd", "Mp": "cfd", "Mpp": "hall@x0_cfd"})
    hs = il.get("sha256_16") or {}
    check(f"snapshot・抽出線・熱力学条件のハッシュを記録 ({hs})",
          all(isinstance(hs.get(k), str) and len(hs[k]) == 16 for k in ("snapshot", "line", "thermo")))
    check(f"凍結源の形・ガスを照合した ({il.get('match')})",
          set(("R", "gamma_hall", "gas", "L_U", "r_U", "L_pipe", "contraction_shape")) <= set((il.get("match") or {}).get("checked", [])))
    d_hall = design_chain(load_problem(PROB))
    check(f"Hall 経路のアンカー出所 ({d_hall['anchor_source']})", d_hall["anchor_source"] == {"M": "hall", "Mp": "hall", "Mpp": "hall"})

    # --- 凍結入力契約 (codex result M2 の再現を拒否する) -----------------------------
    from forge_design.feedback.cfd_initial_line import pinned_factory  # noqa: E402
    try:
        pinned_factory(run, "res_6000.h5")(3.0, 1.4)
        check("R=3・γ=1.4 (凍結源と別の形・ガス) を拒否", False)
    except ValueError as e:
        check(f"R=3・γ=1.4 (凍結源と別の形・ガス) を拒否 ({str(e)[:90]}…)", True)
    try:
        pinned_factory(run, None)
        check("snapshot (res) 省略を拒否", False)
    except ValueError:
        check("snapshot (res) 省略を拒否", True)
    for key, val in (("L_U", 11.5), ("r_inlet", 6.4)):
        pp = load_problem(PROB)
        pp.geometry.update({"initial_line": "cfd", "initial_line_res": "res_6000.h5", "initial_line_run": str(run), key: val})
        try:
            design_chain(pp)
            check(f"使う側の {key}={val} (凍結源の縮流部と違う) を拒否", False)
        except ValueError as e:
            check(f"使う側の {key}={val} (凍結源の縮流部と違う) を拒否 ({str(e)[:80]}…)", "凍結源" in str(e))
    pp = load_problem(PROB)
    pp.geometry.update({"initial_line": "cfd", "initial_line_res": "res_6000.h5", "initial_line_run": str(run)})
    pp.spec["Tt"] = 1500.0
    try:
        design_chain(pp)
        check("使う側の Tt (熱力学条件) 違いを拒否", False)
    except ValueError as e:
        check(f"使う側の Tt (熱力学条件) 違いを拒否 ({str(e)[:80]}…)", True)
    # 実効設定: 読み取り内容だけを差し替えた凍結源 (h5 は symlink、書き込みはしない)
    import tempfile  # noqa: E402
    import yaml  # noqa: E402
    for label, edit_cfg, edit_bc in (
            ("viscMethod 2・SST", lambda c: (c["physProp"].update({"viscMethod": 2, "visc": 1.8e-5}),
                                             c.update({"turbulence": {"model": "SST"}})), None),
            ("viscMethod 0 だが visc ≠ 0", lambda c: c["physProp"].update({"visc": 1.8e-5}), None),
            ("乱流 SST のみ", lambda c: c.update({"turbulence": {"model": "SST"}}), None),
            ("壁が no-slip", None, lambda b: b["wall"].update({"kind": "wall"})),
            ("cell 離散化", lambda c: c["mesh"].update({"discretization": "cell"}), None)):
        td = Path(tempfile.mkdtemp(prefix="cfdpin_bad_"))
        for fn in ("nozzle.h5", "res_6000.h5"):
            (td / fn).symlink_to(run / fn)
        (td / "prepare_info.json").write_text((run / "prepare_info.json").read_text())
        cfg = yaml.safe_load((run / "solverConfig.yaml").read_text())
        bc = yaml.safe_load((run / "bcondConfig.yaml").read_text())
        if edit_cfg:
            edit_cfg(cfg)
        if edit_bc:
            edit_bc(bc)
        (td / "solverConfig.yaml").write_text(yaml.safe_dump(cfg))
        (td / "bcondConfig.yaml").write_text(yaml.safe_dump(bc))
        try:
            pinned_factory(td, "res_6000.h5")
            check(f"凍結源の実効設定 {label} を拒否", False)
        except ValueError as e:
            check(f"凍結源の実効設定 {label} を拒否 ({str(e)[:100]}…)", True)
    # 対照: 同じ手順で書き直しただけの (変更なし) 設定は受理する
    td = Path(tempfile.mkdtemp(prefix="cfdpin_ok_"))
    for fn in ("nozzle.h5", "res_6000.h5"):
        (td / fn).symlink_to(run / fn)
    (td / "prepare_info.json").write_text((run / "prepare_info.json").read_text())
    (td / "solverConfig.yaml").write_text(yaml.safe_dump(yaml.safe_load((run / "solverConfig.yaml").read_text())))
    (td / "bcondConfig.yaml").write_text(yaml.safe_dump(yaml.safe_load((run / "bcondConfig.yaml").read_text())))
    ok_t = pinned_factory(td, "res_6000.h5")(2.0, d["gamma_hall"])
    check("対照: 設定を変えずに書き直した凍結源は受理", ok_t.x0_cfd == il["x0"])

    # --- 実効入口・熱力学条件 (codex diagnose 2026-10-06 result-interpretation の判別 A/B を回帰化) ----------
    # 凍結源の入力を一時ディレクトリへ写し (h5 は symlink、元 run は書き換えない)、入口 BC・physProp だけを変える。
    # A = 入口 Tt 1600 (そのまま) は受理、B = 1500 は拒否。使う側は同じ問題 (spec.Tt 1600)。
    def _src_copy(edit_bc=None, edit_cfg=None, prefix="cfdpin_eff_"):
        td = Path(tempfile.mkdtemp(prefix=prefix))
        for fn in ("nozzle.h5", "res_6000.h5"):
            (td / fn).symlink_to(run / fn)
        (td / "prepare_info.json").write_text((run / "prepare_info.json").read_text())
        cfg = yaml.safe_load((run / "solverConfig.yaml").read_text())
        bc = yaml.safe_load((run / "bcondConfig.yaml").read_text())
        if edit_cfg:
            edit_cfg(cfg)
        if edit_bc:
            edit_bc(bc)
        (td / "solverConfig.yaml").write_text(yaml.safe_dump(cfg))
        (td / "bcondConfig.yaml").write_text(yaml.safe_dump(bc))
        return td

    def _chain_with(td):
        pp = load_problem(PROB)
        pp.geometry.update({"initial_line": "cfd", "initial_line_res": "res_6000.h5", "Md_moc_offset": -4.16e-4,
                            "initial_line_run": str(td)})
        return design_chain(pp)
    td = _src_copy(prefix="cfdpin_tt1600_")
    d_a = _chain_with(td)
    ma = d_a["initial_line"].get("match") or {}
    want = {"inlet.Pt", "inlet.Tt", "inlet.Y", "physProp.thermalMethod", "physProp.species", "physProp.thermoHrefTemp",
            "meta.gas.Tt", "meta.species.Y_transport", "meta.species.transported"}
    check(f"A: 入口 Tt 1600 (凍結源そのまま) は受理し実効入口・熱力学条件を照合済みと記録 (欠け {sorted(want - set(ma.get('checked', [])))})",
          want <= set(ma.get("checked", [])) and d_a["initial_line"]["x0"] == il["x0"])
    check("A: 熱力学ハッシュは元の凍結源と同じ (書き直しだけでは変わらない)",
          d_a["initial_line"]["sha256_16"]["thermo"] == il["sha256_16"]["thermo"])
    for label, edit_bc, edit_cfg, key in (
            ("B: 入口 Tt 1500", lambda b: b["inlet"]["floats"].update({"Tt": 1500.0}), None, "入口 Tt"),
            ("入口 Pt 5.0 MPa", lambda b: b["inlet"]["floats"].update({"Pt": 5.0e6}), None, "入口 Pt"),
            ("入口組成 Y0/Y1 の入れ替え", lambda b: b["inlet"]["floats"].update({"Y0": 0.0858, "Y1": 0.9142}), None, "入口 Y"),
            ("thermalMethod 0 (CPG)", None, lambda c: c["physProp"].update({"thermalMethod": 0}), "thermalMethod"),
            ("thermoHrefTemp 0", None, lambda c: c["physProp"].update({"thermoHrefTemp": 0.0}), "thermoHrefTemp"),
            ("species の lump 分率違い", None,
             lambda c: c["physProp"]["species"][0]["lump"].update({"CO2": 0.13}), "physProp.species"),
            ("入口が inlet_Pressure でない", lambda b: b["inlet"].update({"kind": "inlet_uniformVelocity"}), None,
             "inlet_Pressure")):
        td = _src_copy(edit_bc, edit_cfg)
        try:
            _chain_with(td)
            check(f"凍結源の実効入口・熱力学条件 {label} を拒否", False)
        except ValueError as e:
            check(f"凍結源の実効入口・熱力学条件 {label} を拒否 ({str(e)[:110]}…)", key in str(e))
    # 使う側の問題を通さない呼び方 (pinned_factory 単独) でも、凍結源の prepare_info と実効 BC の食い違いは拒否
    td = _src_copy(lambda b: b["inlet"]["floats"].update({"Tt": 1500.0}))
    try:
        pinned_factory(td, "res_6000.h5")(2.0, d["gamma_hall"])
        check("pinned_factory 単独でも入口 Tt 1500 (prepare_info.gas.Tt 1600 と不整合) を拒否", False)
    except ValueError as e:
        check(f"pinned_factory 単独でも入口 Tt 1500 (prepare_info.gas.Tt 1600 と不整合) を拒否 ({str(e)[:90]}…)",
              "prepare_info" in str(e))
else:
    print("skip: res_6000.h5 が無いので design_chain (cfd) と凍結入力契約の検査を省略")

print(json.dumps(out, indent=1))
print("FAIL 件数:", FAIL)
sys.exit(1 if FAIL else 0)
