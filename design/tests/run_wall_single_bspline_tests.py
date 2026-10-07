#!/usr/bin/env python3
"""物理壁の全域 1 本の 5 次 B-spline と保存した壁の復元・STEP の試験 (plans/accepted/tooling-nozzle-wall-single-bspline.md §4・§6 W4・W5、
係数の求め方は plans/accepted/tooling-nozzle-upstream-poly-and-throat-sizing.md §4.1b のノット挿入 — 2026-10-07 に最小二乗の版から直した)。

壁は case/45 の単調壁の生産問題 (`problem_d155_ns_finemesh_recal_final_mono.yaml`) の初期線だけ Hall に差し替え、`pw_ramp`・
`pw_upstream` を外して上流を既定の `poly` にして作る (CFD ピンの凍結源 run に依存しない)。δ_r は prepare_ns と同じ積分法の経路
(`integral_delta_r`、YAML の deltastar_initializer)。

1. 構築: 継ぎ目 (−L_U・0)・定義域が壁の属性から決まる / ノットの重複度 (端 6・継ぎ目 3・他 1) / 元の区分表現との差が丸めの範囲
   (半径 ≤ 1e-12・r′ ≤ 1e-10・r″ ≤ 1e-8) / ノット除去の誤差 ≤ 1e-12 / スロート (大域最小) / 上流のゲート / 必須属性 (方式別)
2. 一般性 (plan §6 W4): 違う scale_m・L_U・L_pipe で、継ぎ目・定義域・mm 換算が壁の属性から決まる
3. 拒否: δ_r の表が [0, x_e] を覆わない (出口側・入口側) / 解析経路でない壁 / ランプの壁 (ramp + single_bspline) / キーの不正値 /
   定義域の外 (float32 丸めを超える)
4. ic.py: 物理壁のスロート属性が欠けたら例外 (黙って (0, 1) に戻らない)、設計壁は従来どおり (0, 1)
5. 壁ファイル (版 2): 書いて読み直し (single_bspline・legacy poly・legacy ramp) が元の壁とビット一致 / 版 1 はランプとして読む /
   復元した評価関数は有効域の外で例外 / 要素の欠損・版・種類・形式・方式で例外
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
                                                 check_required_attrs, load_wall_file, required_attrs, save_wall_file)

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
    g.pop("pw_ramp", None)                 # 上流は既定の poly (ノット挿入の版は poly の壁だけ)
    g.pop("pw_upstream", None)
    if cfg == "gen":
        g["L_U"], g["L_pipe"], g["n_axis_inv"] = 10.0, 1.0, 800
        p.spec["r_throat"] = 0.05
    d = design_chain(p)
    res_init, drx, _ = integral_delta_r(p, d, p.raw["deltastar_initializer"])
    args = (d["wall"], d["wall_inv"], float(p.spec["r_throat"]), float(p.spec["Pt"]), float(p.spec["Tt"]), _gam_or_gas(p), p.cp)
    PW = PhysicalNozzleWall(*args, offset="radial", delta_r_x=drx)
    return p, d, res_init, drx, args, PW, SingleBSplinePhysicalWall(PW)


if not PROB.exists():
    print(f"問題 YAML {PROB} が無い — 全体を飛ばす (合格扱いにしない)")
    sys.exit(2)

# --- 1. 構築 (case/45 の形、上流 poly) -----------------------------------------------------------------------------
p, d, res_init, drx, args, PW, B = build("c45")
w = d["wall"]
S = float(p.spec["r_throat"])
t = np.asarray(B.spline.t)
dist, mult = np.unique(t, return_counts=True)
jt = [-float(w.up.L_U), 0.0]
check(f"元の壁は poly (既定): {PW.pw_upstream} / {PW.pw_upstream_source}", PW.pw_upstream == "poly" and PW.pw_upstream_source == "default")
check(f"継ぎ目 = [−L_U, 0] (壁の属性から): {B.joints}", B.joints == jt)
check(f"定義域 = 壁の属性 [x_in, x_e] = [{B.x_in}, {B.x_e:.6f}]", (B.x_in, B.x_e) == (float(PW.x_in), float(PW.x_e)) and B.x_in == -12.5)
isj = np.isin(dist, jt)
check(f"重複度: 端 6・継ぎ目 3・他 1 (異なるノット {len(dist)}, 係数 {len(B.spline.c)})",
      mult[0] == 6 and mult[-1] == 6 and np.all(mult[isj] == 3) and np.all(mult[~isj][1:-1] == 1))
fe = B.fit_diag["max_err"]
check(f"元の区分表現との差 (密な点 + 区間多項式の極値) が丸めの範囲: r {fe['d0']:.2e} ≤ 1e-12, r′ {fe['d1']:.2e} ≤ 1e-10, r″ {fe['d2']:.2e} ≤ 1e-8",
      fe["d0"] <= 1e-12 and fe["d1"] <= 1e-10 and fe["d2"] <= 1e-8)
kr = B.fit_diag["knot_removal"]
check(f"継ぎ目のノット除去の誤差 ≤ 1e-12 r_t: {[(k_['x'], k_['max_rt']) for k_ in kr]}", all(k_["max_rt"] <= 1e-12 for k_ in kr) and len(kr) == 2)
check(f"構成はノット挿入 (最小二乗なし): {B.fit_diag['construction']}, 挿入 {B.fit_diag['insertion']}",
      B.fit_diag["construction"] == "knot_insertion" and "lstsq_rank" not in B.fit_diag)
xs = np.linspace(-12.5, B.x_e, 20001)
check(f"r(x) が元の壁と一致 (均等 20001 点の max |Δr| {np.abs(B.r(xs) - PW.r(xs)).max():.2e} ≤ 1e-12)",
      np.abs(B.r(xs) - PW.r(xs)).max() <= 1e-12)
jj = B.fit_diag["joint_jumps"]
check(f"継ぎ目の r・r′・r″ の跳び ≤ 1e-8 (区間多項式の左右の極限): max {max(max(j['jump_d0'], j['jump_d1'], j['jump_d2']) for j in jj):.1e}",
      all(max(j["jump_d0"], j["jump_d1"], j["jump_d2"]) <= 1e-8 for j in jj))
check(f"継ぎ目の左右の極限が元の壁の左右の極限と一致 (≤ 1e-10): max {max(v for j in jj for k_, v in j.items() if '_vs_source_' in k_):.1e}",
      all(v <= 1e-10 for j in jj for k_, v in j.items() if "_vs_source_" in k_))
check(f"スロートを 1 本から大域最小で求め直す: Δx {B.x_throat - PW.x_throat:.1e}, Δr {B.r_throat - PW.r_throat:.1e}, Δκ {B.kappa_throat - PW.kappa_throat:.1e}",
      abs(B.x_throat - PW.x_throat) <= 1e-9 and abs(B.r_throat - PW.r_throat) <= 1e-12 and abs(B.kappa_throat - PW.kappa_throat) <= 1e-8)
ug = B.upstream_gate
check(f"上流のゲートを 1 本で評価し直して合格 (|r″ − H″| {ug['max_abs_d2_change_vs_H']:.4e} ≤ {ug['limit_d2']}, 継ぎ目 {ug['max_seam_jump']:.1e}, "
      f"一意 {ug['unique_min']}, 単調 {ug['monotone_before']}/{ug['monotone_after']})",
      ug["pass"] and ug["repr"] == "single_bspline" and ug["source_wall"]["pass"]
      and abs(ug["max_abs_d2_change_vs_H"] - PW.upstream_gate["max_abs_d2_change_vs_H"]) <= 1e-9)
check(f"validate() が空 ({B.validate()})", B.validate() == [])
check("必須属性 (poly + 1 本の B-spline) が全部ある", raises(lambda: check_required_attrs(B)) is None
      and set(required_attrs(B)) >= {"upstream_gate", "upstream_poly", "spline", "joints", "fit_diag"} and "ramp_gate" not in required_attrs(B))
Bm = copy.copy(B)
del Bm.upstream_gate
check("必須属性が欠けたら check_required_attrs が例外 (upstream_gate を消す)", raises(lambda: check_required_attrs(Bm)) is not None)
check("validate() も必須属性の欠落を報告する", any("必須属性" in m for m in Bm.validate()))
Bp = copy.copy(B)
del Bp.pw_upstream
check("pw_upstream が無ければ必須属性の検査が例外 (方式が分からないまま既定に落とさない)", raises(lambda: check_required_attrs(Bp)) is not None)
check("_dstar_hist は元の壁のもの (prepare_info の dstar_throat_correlation が同じ)", float(B._dstar_hist(0.0)) == float(PW._dstar_hist(0.0)))
check(f"r の定義域外: float32 丸め分 ({B._dom_tol:.1e}) の内側は端で評価", float(B.r(np.array([B.x_in - 0.5 * B._dom_tol]))[0]) == float(B.r(np.array([B.x_in]))[0]))
check("r の定義域外: それより外は例外 (外挿しない)", raises(lambda: B.r(np.array([B.x_e + 1e-3]))) is not None)

# --- 2. 一般性: 違う scale_m (0.05 m)・L_U 10・L_pipe 1 ------------------------------------------------------------
pg, dg, _, drg, argsg, PWg, Bg = build("gen")
check(f"継ぎ目 = [−L_U, 0] = {Bg.joints}", Bg.joints == [-10.0, 0.0])
check(f"定義域 [x_in, x_e] = [−L_U − L_pipe, 設計壁の x_e] = [{Bg.x_in}, {Bg.x_e:.5f}]", Bg.x_in == -11.0 and Bg.x_e == float(dg["wall"].x_e))
check(f"一般性の壁も丸めの範囲 ({Bg.fit_diag['max_err']['d0']:.1e}, {Bg.fit_diag['max_err']['d1']:.1e}, {Bg.fit_diag['max_err']['d2']:.1e})",
      Bg.fit_diag["max_err"]["d0"] <= 1e-12 and Bg.fit_diag["max_err"]["d1"] <= 1e-10 and Bg.fit_diag["max_err"]["d2"] <= 1e-8)

# --- 3. 拒否 ------------------------------------------------------------------------------------------------------
xt, dt_ = res_init["x"], drx(res_init["x"])
short_hi = delta_r_from_table(xt[xt < B.x_e - 1.0], dt_[xt < B.x_e - 1.0])
# 表が出口の手前で終わると、その先は δ_r 一定で物理壁 = 設計壁の r′ が出口直前でわずかに負になり、poly の壁そのものがゲートで止まる
e = raises(lambda: PhysicalNozzleWall(*args, offset="radial", delta_r_x=short_hi))
check(f"δ_r の表が出口まで覆わない poly の壁: ゲート (後で r′ ≥ 0) で止まる ({(e or '')[:40]}…)", e is not None and "後で r′ ≥ 0 False" in e)
# 1 本の B-spline 側の範囲の検査は、元の壁の δ_r だけを差し替えて試す (元の壁のゲートとは独立に)
PWs = copy.copy(PW); PWs._dr = short_hi
e = raises(lambda: SingleBSplinePhysicalWall(PWs))
check(f"δ_r の表が出口まで覆わない → 1 本の B-spline は例外 ({(e or '')[:60]}…)", e is not None and "覆わない" in e)
short_lo = delta_r_from_table(xt[xt > 0.5], dt_[xt > 0.5])
PWl = copy.copy(PW); PWl._dr = short_lo
e = raises(lambda: SingleBSplinePhysicalWall(PWl))
check(f"δ_r の表が設計スロート 0 より後 (0.5) から始まる → 例外 ({(e or '')[:50]}…)", e is not None and "覆わない" in e)
PWo = PhysicalNozzleWall(*args, offset="radial", delta_r_x=drx, analytic=False)
check("解析経路でない PhysicalNozzleWall → 例外", raises(lambda: SingleBSplinePhysicalWall(PWo)) is not None)
PWr = PhysicalNozzleWall(*args, offset="radial", delta_r_x=drx, ramp=(-11.0, -6.0), upstream="ramp")
e = raises(lambda: SingleBSplinePhysicalWall(PWr))
check(f"ランプの壁 (pw_upstream ramp) → 例外 ({(e or '')[:70]}…)", e is not None and "ramp" in e)
nodr = lambda x, deriv=0: drx(x, deriv)  # noqa: E731
nodr.supports_deriv = True
check("δ_r の関数に spline・x_range が無い → poly の壁を作らない (区切りが決まらない、黙って既定にしない)",
      raises(lambda: PhysicalNozzleWall(*args, offset="radial", delta_r_x=nodr)) is not None)
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

# --- 5. 壁ファイル (版 2) -----------------------------------------------------------------------------------------
with tempfile.TemporaryDirectory() as td:
    td = Path(td)
    (td / "sb").mkdir(); (td / "lg").mkdir(); (td / "lr").mkdir()
    path_sb, rec_sb, sha_sb = save_wall_file(td / "sb", B, S, "single_bspline")
    path_lg, rec_lg, sha_lg = save_wall_file(td / "lg", PW, S, "legacy")
    path_lr, rec_lr, sha_lr = save_wall_file(td / "lr", PWr, S, "legacy")
    Wsb, Wlg, Wlr = load_wall_file(td / "sb"), load_wall_file(path_lg), load_wall_file(path_lr)
    xq = np.linspace(B.x_in, B.x_e, 50001)
    check("single_bspline: 復元した物理壁が元とビット一致 (r, r′, r″)",
          all(np.array_equal(Wsb["physical"].r(xq, n), B.r(xq, n)) for n in range(3)))
    check("legacy (poly): 復元した物理壁が PhysicalNozzleWall とビット一致 (r, r′, r″, r‴)",
          all(np.array_equal(Wlg["physical"].r(xq, n), PW.r(xq, n)) for n in range(4)))
    check("legacy (ramp): 復元した物理壁が PhysicalNozzleWall とビット一致 (r, r′, r″, r‴)",
          all(np.array_equal(Wlr["physical"].r(xq, n), PWr.r(xq, n)) for n in range(4)))
    check("復元した設計壁 (直管 + 上流 Hermite + S) が設計壁とビット一致 (r, r′, r″)",
          all(np.array_equal(Wsb["design"].r(xq, n), w.r(xq, n)) for n in range(3)))
    check("記録: 形式・版 2・表現の種類・上流の方式・単位 (r_t と scale_m)・有効域",
          rec_sb["format"] == "forge_design.nozzle_wall" and rec_sb["version"] == 2 and rec_sb["physical_wall_repr"] == "single_bspline"
          and rec_sb["pw_upstream"] == "poly" and rec_sb["physical_wall"]["pw_upstream"] == "poly"
          and rec_sb["units"]["length"] == "r_t" and rec_sb["units"]["scale_m"] == S
          and rec_sb["domain"] == [B.x_in, B.x_e] and rec_sb["domain_m"] == [B.x_in * S, B.x_e * S])
    up = rec_lg["physical_wall"]["upstream_poly"]
    check(f"legacy (poly): Q の係数・基底・区間を保存 ({up['basis']}, {up['domain']})",
          rec_lg["pw_upstream"] == "poly" and len(up["coef"]) == 6 and up["domain"] == [-float(w.up.L_U), 0.0]
          and up["coef"] == [float(v) for v in PW._q_c] and "ramp" not in rec_lg["physical_wall"])
    check("読み込みの戻りに方式と版", (Wsb["pw_upstream"], Wlg["pw_upstream"], Wlr["pw_upstream"]) == ("poly", "poly", "ramp")
          and Wsb["version"] == 2)
    # 版 1 (方式の記録が無い) はランプとして読む
    r1 = json.loads(path_lr.read_text())
    r1["version"] = 1
    r1.pop("pw_upstream"); r1["physical_wall"].pop("pw_upstream")
    (td / "v1").mkdir(); (td / "v1" / WALL_FILE).write_text(json.dumps(r1))
    W1 = load_wall_file(td / "v1")
    check("版 1 (方式の記録なし) は ramp として読み、復元した壁がランプの壁とビット一致",
          W1["pw_upstream"] == "ramp" and W1["version"] == 1 and np.array_equal(W1["physical"].r(xq), PWr.r(xq)))
    check("復元した物理壁は有効域の外で例外 (外挿しない)", raises(lambda: Wsb["physical"].r(np.array([B.x_e + 1e-9]))) is not None
          and raises(lambda: Wlg["physical"].r(np.array([B.x_in - 1e-9]))) is not None)
    check("復元した設計壁 S は自分の有効域 [0, x_e] の外で例外", raises(lambda: Wsb["design"].S(np.array([-6.0]))) is not None)

    def broken(src, mut, name):
        r = json.loads(src.read_text())
        mut(r)
        q = td / name
        q.write_text(json.dumps(r))
        return raises(lambda: load_wall_file(q))
    check("要素の欠損: physical_wall.c → 例外", broken(path_sb, lambda r: r["physical_wall"].pop("c"), "m1.json") is not None)
    check("要素の欠損: 係数を 1 個減らす → 例外", broken(path_sb, lambda r: r["physical_wall"]["c"].pop(), "m2.json") is not None)
    check("要素の欠損: design_wall.upstream_hermite → 例外", broken(path_sb, lambda r: r["design_wall"].pop("upstream_hermite"), "m3.json") is not None)
    check("要素の欠損: design_wall.S.t → 例外", broken(path_sb, lambda r: r["design_wall"]["S"].pop("t"), "m4.json") is not None)
    check("上流 Hermite の係数が (r_U, R_t, L_U) と合わない → 例外",
          broken(path_sb, lambda r: r["design_wall"]["upstream_hermite"]["coef"].__setitem__(1, 1e-3), "m5.json") is not None)
    check("版が未対応 (3) → 例外", broken(path_sb, lambda r: r.__setitem__("version", 3), "m6.json") is not None)
    check("表現の種類が違う → 例外", broken(path_sb, lambda r: r.__setitem__("physical_wall_repr", "bspline"), "m7.json") is not None)
    check("形式が違う → 例外", broken(path_sb, lambda r: r.__setitem__("format", "x"), "m8.json") is not None)
    check("physical_wall.kind と physical_wall_repr の食い違い → 例外",
          broken(path_sb, lambda r: r["physical_wall"].__setitem__("kind", "legacy"), "m9.json") is not None)
    check("版 2 で pw_upstream が無い → 例外", broken(path_sb, lambda r: r.pop("pw_upstream"), "m10.json") is not None)
    check("版 2 の single_bspline で pw_upstream ramp → 例外 (ramp + single_bspline は作れない)",
          broken(path_sb, lambda r: (r.__setitem__("pw_upstream", "ramp"), r["physical_wall"].__setitem__("pw_upstream", "ramp")), "m11.json") is not None)
    check("pw_upstream の不正値 → 例外", broken(path_lg, lambda r: r.__setitem__("pw_upstream", "Poly"), "m12.json") is not None)
    check("最上位と物理壁の pw_upstream の食い違い → 例外",
          broken(path_lg, lambda r: r["physical_wall"].__setitem__("pw_upstream", "ramp"), "m13.json") is not None)
    check("legacy (poly) で upstream_poly が無い → 例外", broken(path_lg, lambda r: r["physical_wall"].pop("upstream_poly"), "m14.json") is not None)
    check("legacy (poly) の Q の係数が端条件と合わない (継ぎ目の跳び) → 例外",
          broken(path_lg, lambda r: r["physical_wall"]["upstream_poly"]["coef"].__setitem__(3, r["physical_wall"]["upstream_poly"]["coef"][3] + 1e-6),
                 "m15.json") is not None)
    check("legacy (ramp) で ramp が無い → 例外", broken(path_lr, lambda r: r["physical_wall"].pop("ramp"), "m16.json") is not None)
    check("版 1 なのに pw_upstream がある → 例外",
          broken(path_lr, lambda r: r.__setitem__("version", 1), "m17.json") is not None)
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
    # 判定器の負例 (FreeCAD 不要; 2026-10-07 result 段レビュー M1): 保存した曲線から作った「完全な読み直し」を基準に、
    # z の移動・途中で切れた辺・頂点の欠落・回転面の不正を検出すること
    import copy
    from forge_design.export.wall_step import revolve_ok
    dat = step_curve_data(Wsb["record"])
    pts = sample_params(dat, 1)
    dbp = DeBoor(dat["t_mm"], dat["poles_mm"], dat["degree"], nder=2)
    k0, k1 = float(dat["knots_mm"][0]), float(dat["knots_mm"][-1])
    e0, e1 = dbp(k0, "right")[0], dbp(k1, "left")[0]
    rd0 = {"n_edges": 1, "curve_type": "BSplineCurve", "degree": dat["degree"], "n_poles": dat["n_poles"], "is_rational": False,
           "knots": [float(v) for v in dat["knots_mm"]], "mults": [int(m) for m in dat["mults"]],
           "poles": [[float(x), float(y), 0.0] for x, y in dat["poles_mm"]], "edge_first": k0, "edge_last": k1,
           "vertices": [[float(e0[0]), float(e0[1]), 0.0], [float(e1[0]), float(e1[1]), 0.0]],
           "eval": [[[*map(float, r[0]), 0.0], [*map(float, r[1]), 0.0], [*map(float, r[2]), 0.0]] for r in (dbp(u, sd) for u, sd, _ in pts)]}
    check("STEP 判定器: 保存した曲線そのものの読み直しは合格", transfer_check(dat, rd0, pts)["pass"])
    rz = copy.deepcopy(rd0)
    for q in rz["poles"] + rz["vertices"]:
        q[2] += 1.0
    for ev_ in rz["eval"]:
        ev_[0][2] += 1.0
    trz = transfer_check(dat, rz, pts)
    check(f"STEP 判定器 (負例): 全体を z に 1 mm 動かすと不合格 (平面性 {trz['structure']['planar']}・位置 {trz['pos_max_mm']:.2g} mm)",
          (not trz["pass"]) and not trz["structure"]["planar"] and trz["pos_max_mm"] >= 0.999)
    rm = copy.deepcopy(rd0)
    rm["edge_last"] = 0.5 * (k0 + k1)
    trm = transfer_check(dat, rm, pts)
    check("STEP 判定器 (負例): 辺の終わりを定義域の中央にすると不合格", (not trm["pass"]) and not trm["structure"]["edge_full_range"])
    rv_ = copy.deepcopy(rd0)
    rv_["vertices"] = [rd0["vertices"][0], [0.5 * (k0 + k1), float(dbp(0.5 * (k0 + k1))[0][1]), 0.0]]
    trv = transfer_check(dat, rv_, pts)
    check("STEP 判定器 (負例): 実際の端点が曲線の端でないと不合格", (not trv["pass"]) and not trv["structure"]["edge_endpoints"])
    rn = copy.deepcopy(rd0)
    rn.pop("vertices")
    check("STEP 判定器 (負例): 頂点が無ければ不合格", not transfer_check(dat, rn, pts)["pass"])
    check("回転面の判定: 妥当・面 1 枚以上・面積 > 0 だけ合格",
          revolve_ok({"is_valid": True, "n_faces": 1, "area_mm2": 1.0})
          and not revolve_ok({"is_valid": False, "n_faces": 1, "area_mm2": 1.0})
          and not revolve_ok({"is_valid": True, "n_faces": 0, "area_mm2": 1.0})
          and not revolve_ok({"is_valid": True, "n_faces": 1, "area_mm2": 0.0}) and not revolve_ok(None))
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
        check(f"回転面 (内面) を作れる: {rv}", revolve_ok(rv))
    else:
        SKIPPED.append(f"STEP の FreeCAD 往復 ({FREECADCMD} が無い)")

print(f"FAIL 件数: {FAIL}")
if SKIPPED:
    print("飛ばした項目 (合格扱いにしない):", SKIPPED)
sys.exit(1 if FAIL else (2 if SKIPPED else 0))
