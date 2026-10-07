#!/usr/bin/env python3
"""物理壁の全域 1 本の 5 次 B-spline と保存した壁の復元・STEP の試験 (plans/active/tooling-nozzle-wall-single-bspline.md §4・§6 W4・W5)。

壁は case/45 の単調壁の生産問題 (`problem_d155_ns_finemesh_recal_final_mono.yaml`) の初期線だけ Hall に差し替えて作る
(CFD ピンの凍結源 run に依存しない)。δ_r は prepare_ns と同じ積分法の経路 (`integral_delta_r`、YAML の deltastar_initializer)。

1. 構築: 継ぎ目・定義域が壁の属性から決まる / ノットの重複度 (端 6・継ぎ目 3・他 1) / 許容誤差内 / スロート / ランプのゲート / 必須属性
2. 一般性 (plan §6 W4): 違う scale_m・L_U・L_pipe と既定ランプ (pw_ramp 無し) で、継ぎ目・定義域・mm 換算が壁の属性から決まる
3. 拒否: δ_r の表がランプ開始〜出口を覆わない (出口側・入口側) / 解析経路でない壁 / キーの不正値 / 定義域の外 (float32 丸めを超える)
4. ic.py: 物理壁のスロート属性が欠けたら例外 (黙って (0, 1) に戻らない)、設計壁は従来どおり (0, 1)
5. 壁ファイル: 書いて読み直し (single_bspline・legacy) が元の壁とビット一致 / 復元した評価関数は有効域の外で例外 /
   要素の欠損・版・種類・形式で例外
6. 報告の単体検査 (plan §6 W4 の 4 種): 旧 run (壁ファイル無し → 旧経路、図にその旨) / 新しい形式 (保存した係数から評価、CSV を読まない) /
   係数の欠損 (例外) / 上流を含む差分図 (x = −6 で設計壁が上流 Hermite の値、S を有効域の外で外挿しない)
7. STEP: グレビル点で x(u) = u・ノットの mm 換算 (1000·scale_m)・自前の de Boor と scipy の一致 /
   freecadcmd があれば実寸・全制御点で書き出し → 読み直しの構造と転送誤差・回転面 (無ければ飛ばして exit 2 — 合格扱いにしない)

usage: design/.venv-opt/bin/python design/tests/run_wall_single_bspline_tests.py   (所要 1 分程度)
"""
import copy
import json
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from forge_design.evaluate import ic as ic_mod  # noqa: E402
from forge_design.evaluate.runner_axismach import (_gam_or_gas, _physical_wall_repr, delta_r_from_table,  # noqa: E402
                                                   design_chain, integral_delta_r, load_problem)
from forge_design.geometry.wall_axismach import (WALL_FILE, PhysicalNozzleWall, SingleBSplinePhysicalWall,  # noqa: E402
                                                 check_required_attrs, default_pw_ramp, load_wall_file, save_wall_file)

ROOT = Path(__file__).resolve().parents[2]
PROB = ROOT / "case/45.isobutane_m6_d155/problem_d155_ns_finemesh_recal_final_mono.yaml"
FAIL = 0
SKIPPED = []


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


def build(cfg: str):
    p = load_problem(PROB)
    g = p.geometry
    g["initial_line"] = "hall"
    g.pop("initial_line_run", None)
    g.pop("initial_line_res", None)
    if cfg == "gen":
        g.pop("pw_ramp")
        g["L_U"], g["L_pipe"], g["n_axis_inv"] = 10.0, 1.0, 800
        p.spec["r_throat"] = 0.05
    d = design_chain(p)
    res_init, drx, _ = integral_delta_r(p, d, p.raw["deltastar_initializer"])
    pw = g.get("pw_ramp")
    args = (d["wall"], d["wall_inv"], float(p.spec["r_throat"]), float(p.spec["Pt"]), float(p.spec["Tt"]), _gam_or_gas(p), p.cp)
    PW = PhysicalNozzleWall(*args, offset="radial", delta_r_x=drx, ramp=(None if pw is None else tuple(float(v) for v in pw)))
    return p, d, res_init, drx, args, PW, SingleBSplinePhysicalWall(PW)


if not PROB.exists():
    print(f"問題 YAML {PROB} が無い — 全体を飛ばす (合格扱いにしない)")
    sys.exit(2)

# --- 1. 構築 (case/45 の形、pw_ramp [−11, −6]) -------------------------------------------------------------------
p, d, res_init, drx, args, PW, B = build("c45")
w = d["wall"]
S = float(p.spec["r_throat"])
t = np.asarray(B.spline.t)
dist, mult = np.unique(t, return_counts=True)
jt = [-float(w.up.L_U), -11.0, -6.0, 0.0]
check(f"継ぎ目 = [−L_U, ランプ両端, 0] (壁の属性から): {B.joints}", B.joints == jt)
check(f"定義域 = 壁の属性 [x_in, x_e] = [{B.x_in}, {B.x_e:.6f}]", (B.x_in, B.x_e) == (float(PW.x_in), float(PW.x_e)) and B.x_in == -12.5)
isj = np.isin(dist, jt)
check(f"重複度: 端 6・継ぎ目 3・他 1 (異なるノット {len(dist)}, 係数 {len(B.spline.c)})",
      mult[0] == 6 and mult[-1] == 6 and np.all(mult[isj] == 3) and np.all(mult[~isj][1:-1] == 1))
fe = B.fit_diag["max_err"]
check(f"元の壁との差 (密な点 + 区間多項式の極値): r {fe['d0']:.2e} ≤ 1.3e-7, r′ {fe['d1']:.2e} ≤ 1e-7, r″ {fe['d2']:.2e} ≤ 1e-5",
      fe["d0"] <= 1.3e-7 and fe["d1"] <= 1e-7 and fe["d2"] <= 1e-5)
xs = np.linspace(-12.5, B.x_e, 20001)
check(f"r(x) が元の壁と一致 (均等 20001 点の max |Δr| {np.abs(B.r(xs) - PW.r(xs)).max():.2e} ≤ 1.3e-7)",
      np.abs(B.r(xs) - PW.r(xs)).max() <= 1.3e-7)
jj = B.fit_diag["joint_jumps"]
check(f"継ぎ目の r・r′・r″ の跳び ≤ 1e-8 (区間多項式の左右の極限): max {max(max(j['jump_d0'], j['jump_d1'], j['jump_d2']) for j in jj):.1e}",
      all(max(j["jump_d0"], j["jump_d1"], j["jump_d2"]) <= 1e-8 for j in jj))
check(f"スロートを 1 本から求め直す: Δx {B.x_throat - PW.x_throat:.1e}, Δr {B.r_throat - PW.r_throat:.1e}, Δκ {B.kappa_throat - PW.kappa_throat:.1e}",
      abs(B.x_throat - PW.x_throat) <= 1e-6 and abs(B.r_throat - PW.r_throat) <= 1.3e-7 and abs(B.kappa_throat - PW.kappa_throat) <= 1e-5)
check(f"ランプのゲートを 1 本で評価し直して合格 ({B.ramp_gate['max_abs_d2_change']:.4e} ≤ {B.ramp_gate['limit_d2']}, max r′ {B.ramp_gate['max_r1']:.3f})",
      B.ramp_gate["pass"] and B.ramp_gate["repr"] == "single_bspline" and B.ramp_gate["source_wall"] == PW.ramp_gate)
check(f"validate() が空 ({B.validate()})", B.validate() == [])
check("必須属性 (REQUIRED_ATTRS) が全部ある", raises(lambda: check_required_attrs(B)) is None)
Bm = copy.copy(B)
del Bm.ramp_gate
check("必須属性が欠けたら check_required_attrs が例外 (ramp_gate を消す)", raises(lambda: check_required_attrs(Bm)) is not None)
check("validate() も必須属性の欠落を報告する", any("必須属性" in m for m in Bm.validate()))
check("_dstar_hist は元の壁のもの (prepare_info の dstar_throat_correlation が同じ)", float(B._dstar_hist(0.0)) == float(PW._dstar_hist(0.0)))
check(f"r の定義域外: float32 丸め分 ({B._dom_tol:.1e}) の内側は端で評価", float(B.r(np.array([B.x_in - 0.5 * B._dom_tol]))[0]) == float(B.r(np.array([B.x_in]))[0]))
check("r の定義域外: それより外は例外 (外挿しない)", raises(lambda: B.r(np.array([B.x_e + 1e-3]))) is not None)

# --- 2. 一般性: 違う scale_m (0.05 m)・L_U 10・L_pipe 1・既定ランプ ----------------------------------------------
pg, dg, _, drg, argsg, PWg, Bg = build("gen")
lo_d, hi_d = default_pw_ramp(dg["wall"])
check(f"既定ランプ (pw_ramp 無し) の継ぎ目 = [−L_U, 既定 lo {lo_d:.5f}, −0.5·L_U, 0]: {Bg.joints}",
      Bg.joints == [-10.0, lo_d, -5.0, 0.0] and PWg._ramp_source == "default")
check(f"定義域 [x_in, x_e] = [−L_U − L_pipe, 設計壁の x_e] = [{Bg.x_in}, {Bg.x_e:.5f}]", Bg.x_in == -11.0 and Bg.x_e == float(dg["wall"].x_e))
check(f"一般性の壁も許容誤差内 ({Bg.fit_diag['max_err']['d0']:.1e}, {Bg.fit_diag['max_err']['d1']:.1e}, {Bg.fit_diag['max_err']['d2']:.1e})",
      Bg.fit_diag["max_err"]["d0"] <= 1.3e-7 and Bg.fit_diag["max_err"]["d1"] <= 1e-7 and Bg.fit_diag["max_err"]["d2"] <= 1e-5)

# --- 3. 拒否 ------------------------------------------------------------------------------------------------------
xt, dt_ = res_init["x"], drx(res_init["x"])
short_hi = delta_r_from_table(xt[xt < B.x_e - 1.0], dt_[xt < B.x_e - 1.0])
PWs = PhysicalNozzleWall(*args, offset="radial", delta_r_x=short_hi, ramp=(-11.0, -6.0))
e = raises(lambda: SingleBSplinePhysicalWall(PWs))
check(f"δ_r の表が出口まで覆わない → 例外 ({(e or '')[:60]}…)", e is not None and "覆わない" in e)
short_lo = delta_r_from_table(xt[xt > -10.0], dt_[xt > -10.0])
PWl = PhysicalNozzleWall(*args, offset="radial", delta_r_x=short_lo, ramp=(-11.0, -6.0))
check("δ_r の表がランプ開始 (−11) より後から始まる → 例外", raises(lambda: SingleBSplinePhysicalWall(PWl)) is not None)
PWo = PhysicalNozzleWall(*args, offset="radial", delta_r_x=drx, analytic=False)
check("解析経路でない PhysicalNozzleWall → 例外", raises(lambda: SingleBSplinePhysicalWall(PWo)) is not None)
nodr = lambda x, deriv=0: drx(x, deriv)  # noqa: E731
nodr.supports_deriv = True
PWn = PhysicalNozzleWall(*args, offset="radial", delta_r_x=nodr, ramp=(-11.0, -6.0))
check("δ_r の関数に spline・x_range が無い → 例外 (黙って既定にしない)", raises(lambda: SingleBSplinePhysicalWall(PWn)) is not None)
check("キー無し → None (既定 = legacy、壁ファイルなし)", _physical_wall_repr({}) is None)
check("キー legacy / single_bspline を受ける", _physical_wall_repr({"physical_wall_repr": "legacy"}) == "legacy"
      and _physical_wall_repr({"physical_wall_repr": "single_bspline"}) == "single_bspline")
check("キーの不正値 (null・大文字・空白・数値・真偽値) は例外",
      all(raises(lambda v=v: _physical_wall_repr({"physical_wall_repr": v})) is not None
          for v in (None, "Legacy", " legacy", "single-bspline", 1, True)))

# --- 4. ic.py のスロート ------------------------------------------------------------------------------------------
check("ic: SingleBSplinePhysicalWall のスロート = (x_throat, r_throat)", ic_mod._throat_of(B) == (B.x_throat, B.r_throat))
check("ic: PhysicalNozzleWall のスロート = (x_throat, r_throat) (従来と同じ値)", ic_mod._throat_of(PW) == (PW.x_throat, PW.r_throat))
check("ic: 設計壁 (物理スロートを持たない) は従来どおり (0, 1)", ic_mod._throat_of(w) == (0.0, 1.0))
Bx = copy.copy(B)
del Bx.x_throat
check("ic: 物理壁のスロート属性が欠けたら例外 (黙って (0, 1) に戻らない)", raises(lambda: ic_mod._throat_of(Bx)) is not None)
Bn = copy.copy(B)
Bn.r_throat = float("nan")
check("ic: 物理壁のスロートが非有限なら例外", raises(lambda: ic_mod._throat_of(Bn)) is not None)

# --- 5. 壁ファイル ------------------------------------------------------------------------------------------------
with tempfile.TemporaryDirectory() as td:
    td = Path(td)
    (td / "sb").mkdir(); (td / "lg").mkdir()
    path_sb, rec_sb, sha_sb = save_wall_file(td / "sb", B, S, "single_bspline")
    path_lg, rec_lg, sha_lg = save_wall_file(td / "lg", PW, S, "legacy")
    Wsb, Wlg = load_wall_file(td / "sb"), load_wall_file(path_lg)
    xq = np.linspace(B.x_in, B.x_e, 50001)
    check("single_bspline: 復元した物理壁が元とビット一致 (r, r′, r″)",
          all(np.array_equal(Wsb["physical"].r(xq, n), B.r(xq, n)) for n in range(3)))
    check("legacy: 復元した物理壁が PhysicalNozzleWall とビット一致 (r, r′, r″, r‴)",
          all(np.array_equal(Wlg["physical"].r(xq, n), PW.r(xq, n)) for n in range(4)))
    check("復元した設計壁 (直管 + 上流 Hermite + S) が設計壁とビット一致 (r, r′, r″)",
          all(np.array_equal(Wsb["design"].r(xq, n), w.r(xq, n)) for n in range(3)))
    check("記録: 形式・版・表現の種類・単位 (r_t と scale_m)・有効域",
          rec_sb["format"] == "forge_design.nozzle_wall" and rec_sb["version"] == 1 and rec_sb["physical_wall_repr"] == "single_bspline"
          and rec_sb["units"]["length"] == "r_t" and rec_sb["units"]["scale_m"] == S
          and rec_sb["domain"] == [B.x_in, B.x_e] and rec_sb["domain_m"] == [B.x_in * S, B.x_e * S])
    check("復元した物理壁は有効域の外で例外 (外挿しない)", raises(lambda: Wsb["physical"].r(np.array([B.x_e + 1e-9]))) is not None
          and raises(lambda: Wlg["physical"].r(np.array([B.x_in - 1e-9]))) is not None)
    check("復元した設計壁 S は自分の有効域 [0, x_e] の外で例外", raises(lambda: Wsb["design"].S(np.array([-6.0]))) is not None)

    def broken(mut, name):
        r = json.loads(path_sb.read_text())
        mut(r)
        q = td / name
        q.write_text(json.dumps(r))
        return raises(lambda: load_wall_file(q))
    check("要素の欠損: physical_wall.c → 例外", broken(lambda r: r["physical_wall"].pop("c"), "m1.json") is not None)
    check("要素の欠損: 係数を 1 個減らす → 例外", broken(lambda r: r["physical_wall"]["c"].pop(), "m2.json") is not None)
    check("要素の欠損: design_wall.upstream_hermite → 例外", broken(lambda r: r["design_wall"].pop("upstream_hermite"), "m3.json") is not None)
    check("要素の欠損: design_wall.S.t → 例外", broken(lambda r: r["design_wall"]["S"].pop("t"), "m4.json") is not None)
    check("上流 Hermite の係数が (r_U, R_t, L_U) と合わない → 例外",
          broken(lambda r: r["design_wall"]["upstream_hermite"]["coef"].__setitem__(1, 1e-3), "m5.json") is not None)
    check("版が違う → 例外", broken(lambda r: r.__setitem__("version", 2), "m6.json") is not None)
    check("表現の種類が違う → 例外", broken(lambda r: r.__setitem__("physical_wall_repr", "bspline"), "m7.json") is not None)
    check("形式が違う → 例外", broken(lambda r: r.__setitem__("format", "x"), "m8.json") is not None)
    check("physical_wall.kind と physical_wall_repr の食い違い → 例外",
          broken(lambda r: r["physical_wall"].__setitem__("kind", "legacy"), "m9.json") is not None)
    check("壁ファイルの書き出しは legacy に解析経路でない壁を受けない", raises(lambda: save_wall_file(td, PWo, S, "legacy")) is not None)

    # --- 6. 報告の単体検査 (W4) ---------------------------------------------------------------------------------
    from forge_design.report.nozzle_report import fig_wall_shape
    F = {"S": S}
    old = td / "old_run"; old.mkdir()
    xs_p = np.linspace(PW.x_throat, PW.x_e, 1500)
    np.savetxt(old / "wall_physical.csv", np.c_[xs_p * S, PW.r(xs_p) * S], delimiter=",", header="x_m,r_m", comments="")
    np.savetxt(old / "wall_design.csv", np.c_[d["wall_inv"] * [S, S, 1.0, 1.0]], delimiter=",", header="x_m,r_m,theta_rad,M_wall", comments="")
    o1 = fig_wall_shape(old, F, old / "fig.png")
    check(f"報告 (旧 run): 壁ファイル無し → 旧経路 (CSV)、図にその旨 ({o1.get('wall_source')})",
          o1.get("wall_source") == "legacy_csv" and "旧経路" in o1.get("wall_source_note", "") and (old / "fig.png").exists())
    new = td / "new_run"; new.mkdir()
    save_wall_file(new, B, S, "single_bspline")
    # 新しい形式: わざと違う CSV を置いても、保存した係数から評価する (CSV を読まない)
    np.savetxt(new / "wall_physical.csv", np.c_[xs_p * S, 1.01 * PW.r(xs_p) * S], delimiter=",", header="x_m,r_m", comments="")
    o2 = fig_wall_shape(new, F, new / "fig.png")
    check(f"報告 (新しい形式): 保存した係数から評価 ({o2.get('wall_source')}; 出口半径 {o2.get('exit_radius_m'):.9f} m = 係数の値、CSV の 1.01 倍ではない)",
          o2.get("wall_source") == "saved_coefficients" and o2["exit_radius_m"] == float(B.r(np.array([B.x_e]))[0] * S))
    check(f"報告 (新しい形式): 評価量が旧経路と同じ定義で近い (r″ 高周波 {o2['r2_highfreq_max_x_gt2']:.3e} vs CSV 経路 {o1['r2_highfreq_max_x_gt2']:.3e})",
          abs(o2["r2_highfreq_max_x_gt2"] - o1["r2_highfreq_max_x_gt2"]) <= 1e-3 * max(o1["r2_highfreq_max_x_gt2"], 1e-12) + 1e-9)
    bad = td / "bad_run"; bad.mkdir()
    r_ = json.loads((new / WALL_FILE).read_text()); r_["physical_wall"].pop("c")
    (bad / WALL_FILE).write_text(json.dumps(r_))
    check("報告 (係数の欠損): 例外 (旧経路に黙って落ちない)", raises(lambda: fig_wall_shape(bad, F, bad / "fig.png")) is not None)
    Wn = load_wall_file(new)
    h6 = float(Wn["design"].r(np.array([-6.0]))[0])
    from scipy.interpolate import BSpline
    ext = float(BSpline(w._spl.t, w._spl.c, 5)(-6.0))
    check(f"報告 (上流を含む差分図): x = −6 の設計壁 = 上流 Hermite {h6:.6f} (= {float(w.up.r(np.array([-6.0]))[0]):.6f})、"
          f"S を外挿した {ext:.4f} ではない", h6 == float(w.up.r(np.array([-6.0]))[0]) and abs(h6 - ext) > 1.0)
    check("報告 (上流を含む差分図): 図は入口から (x_in の差分 = 物理壁 − 設計壁 が評価できる)",
          np.isfinite(Wn["physical"].r(np.array([B.x_in]))[0] - Wn["design"].r(np.array([B.x_in]))[0]))

    # --- 7. STEP ------------------------------------------------------------------------------------------------
    from forge_design.export.wall_step import (FREECADCMD, DeBoor, read_step, sample_params, step_curve_data,
                                               transfer_check, write_step)
    for lab, W_, Bx_ in (("case/45 の形", Wsb, B), ("一般性 (scale 0.05 m)", None, Bg)):
        if W_ is None:
            (td / "gen").mkdir()
            save_wall_file(td / "gen", Bx_, float(pg.spec["r_throat"]), "single_bspline")
            W_ = load_wall_file(td / "gen")
        dat = step_curve_data(W_["record"])
        s_ = 1000.0 * W_["scale_m"]
        check(f"STEP ({lab}): mm 換算 = 1000·scale_m = {s_:g} (壁の属性から)、ノット = t·{s_:g}",
              dat["scale_mm_per_rt"] == s_ and np.array_equal(dat["t_mm"], np.asarray(Bx_.spline.t) * s_))
        db = DeBoor(dat["t_mm"], dat["poles_mm"], 5, nder=2)
        uu = np.linspace(dat["domain_mm"][0], dat["domain_mm"][1], 997)
        ev = np.array([db(u)[0] for u in uu])
        check(f"STEP ({lab}): グレビル点の制御点で x(u) = u (max |x − u| {np.abs(ev[:, 0] - uu).max():.1e} mm ≤ 1e-9)",
              np.abs(ev[:, 0] - uu).max() <= 1e-9)
        check(f"STEP ({lab}): 自前の de Boor と保存した B-spline (scipy) の r が一致 (max {np.abs(ev[:, 1] - Bx_.spline(uu / s_) * s_).max():.1e} mm)",
              np.abs(ev[:, 1] - Bx_.spline(uu / s_) * s_).max() <= 1e-9)
    if Path(FREECADCMD).is_file():
        dat = step_curve_data(Wsb["record"])
        stp = td / "wall.step"
        wr = write_step(dat, stp)
        pts = sample_params(dat, 2)
        rd = read_step(stp, [u for u, _, _ in pts], revolve=True)
        tr = transfer_check(dat, rd, pts)
        check(f"STEP 書き出し: 次数 {wr['degree']}・制御点 {wr['n_poles']}・非有理", wr["degree"] == 5 and wr["n_poles"] == dat["n_poles"] and not wr["is_rational"])
        check(f"STEP 読み直し: 構造一致 {tr['structure']}", all(tr["structure"].values()))
        check(f"STEP 転送誤差: 位置 {tr['pos_max_mm']:.1e} mm ≤ 1e-6、角度 {tr['angle_max_rad']:.1e} ≤ 1e-9、曲率 {tr['kappa_abs_max_per_mm']:.1e} /mm",
              tr["pass"])
        rv = rd.get("revolve") or {}
        check(f"回転面 (内面) を作れる: {rv}", rv.get("is_valid") is True and rv.get("area_mm2", 0) > 0)
    else:
        SKIPPED.append(f"STEP の FreeCAD 往復 ({FREECADCMD} が無い)")

print(f"FAIL 件数: {FAIL}")
if SKIPPED:
    print("飛ばした項目 (合格扱いにしない):", SKIPPED)
sys.exit(1 if FAIL else (2 if SKIPPED else 0))
