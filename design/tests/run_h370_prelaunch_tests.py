#!/usr/bin/env python3
"""V-h370 の投入前の実装 (plans/active/tooling-nozzle-isothermal-wall-chain.md §5.1 #36e、§4.9 codex diagnose 3・§6 V-h370) の試験。
CFD・変換器は使わない (格子は生成器の座標、場は合成)。

  A. Euler の参照の壁 (案 A、`UpstreamDesignWall`): 凍結した生産の設計 (CFD ピン run_0062) で H の上流を付けた設計壁が設計壁のゲート
     (`monotone_gate`) で通る。物理壁のゲートは今と同じ (`THROAT_R1_TOL` = 1e-12、物理壁のクラスに δ_r = 0 を入れた壁は今も落ちる)。
     負例: S に 1e-5 r_t の凹み → 設計壁のゲートで不合格、上流が単調でない H → 不合格。物理壁の表現の明示は Euler の prepare で例外。
     壁ファイル (design_upstream) の保存 → 復元で Euler の格子 (mesh_euler、throat_refine 5) の座標・接続がビット一致、壊れた壁ファイルは例外。
  B. 生成の座標の δ_in の照合: NS (Gref H・370 K の問題、δ_r = delta_r_initial_C.csv) と Euler (throat_refine 5) で x ≥ 0 の
     max|δ_in − δ_r| ≤ 0.001 mm・全列の被覆 (throat_refine 4 では超える = 旧い FAIL の再現)。NS の格子で一定 x の断面への補間の端の不足が 0 列。
  C. 停止規則のキーの解決 (370 K の NS の問題 → `_stop_rule_settings`): 上限 200,000・保持 40,000・閾値・断面 (40・70・94 r_t)・抽出の条件。
     キーの無い問題は今と同じ (None)。
  D. 抽出: 一定 x の断面への補間 (remap_constant_x=True) が帯の生成 (`fix_band`)・関門の抽出・出力ごとの δ の抽出に渡る。端の補間範囲の
     不足は remap_edge="missing" で欠測 (not_extracted; 既定の clamp は今と同じく端値で全列)。Euler の参照の固定: 指定の res を読む・sha 不一致は
     例外 (抽出器・prepare_ns の停止規則の設定・本段の起動前・CLI)。ディレクトリだけは最新の res (今と同じ)。
  E. 見張りの記録: `delta_snapshot` が使用域の全列の δ_E・δ_target を返し、`DeltaRecorder` が出力ごとに記録する (古い出力は飛ばす)。
  F. 受け渡しの関門の各条件 (合成の run・本物の抽出、収束・y₁⁺・次の壁は偽物で切り替え): 全部満たせば進む・表を run の直下に写す。
     負例で止まる: 保存則・y₁⁺ (FAIL・判定不能)・被覆 (端の欠測)・使用域の held・2 窓 (窓の中の振れ・窓の平均の差・記録の欠け・覆わない)・
     次の壁・Euler の参照の違い・ω の違い・抽出の例外 (sha の不一致)。上書きの 2 つでは足した条件を通さない。
  G. 次の壁の幾何 (`next_wall_check`、本物の壁の構築): 積分法の δ_r で合格、出口側に凹みを入れた δ で不合格。
  H. 1 出力あたりの抽出の時間 (`delta_snapshot`) を実寸の格子 (NS 4719 × 160・Euler 2000 × 97) で測る (--quick で省く)。

usage: /home/sano/work/forge/design/.venv-opt/bin/python design/tests/run_h370_prelaunch_tests.py [--quick] [--pin <run_0062 のディレクトリ>]
  --pin が無ければ case/45 の run_0062_euler_wallfit_fit_r1_ext6k、無ければ /home/sano/work/forge の同じ run (環境変数 H370_PIN も可)。
  ピンが無ければ A・B・H を飛ばし、ALL PASS にしない (終了コード 2)。所要: --quick 約 1 分、全部で約 3〜4 分。
"""
from __future__ import annotations

import copy
import json
import os
import shutil
import sys
import tempfile
import time
from pathlib import Path

import h5py
import numpy as np
import yaml

DESIGN = Path(__file__).resolve().parents[1]
ROOT = DESIGN.parent
CASE = ROOT / "case" / "45.isobutane_m6_d155"
sys.path.insert(0, str(DESIGN))
sys.path.insert(0, str(CASE))
os.environ.setdefault("FORGE_NO_DESIGN_REPORT", "1")

from forge_design.evaluate import runner_axismach as RA  # noqa: E402
from forge_design.evaluate import stop_rule as SR  # noqa: E402
from forge_design.feedback import deltastar_loop as DL  # noqa: E402
from forge_design.geometry import wall_axismach as WA  # noqa: E402
from forge_design.meshing.mesh2d import generate_axisym_mesh  # noqa: E402
from forge_design.metrics import deltastar as DS  # noqa: E402
from forge_design.metrics.delta_in_check import DELTA_IN_TOL_MM, compare_delta_in  # noqa: E402

QUICK = "--quick" in sys.argv[1:]
TMP = Path(tempfile.mkdtemp(prefix="h370_prelaunch_"))
FAIL = 0
SKIPPED = []


def check(name, cond, detail=""):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name + (f"  [{detail}]" if detail else ""), flush=True)
    if not cond:
        FAIL += 1


def skip(name, why):
    SKIPPED.append(name)
    print(f"skip {name}  [{why}]", flush=True)


def raises(fn, exc=Exception):
    try:
        fn()
    except exc as e:  # noqa: BLE001
        return f"{type(e).__name__}: {e}"
    return None


def find_pin():
    a = sys.argv[1:]
    cand = [Path(a[a.index("--pin") + 1])] if "--pin" in a else []
    if os.environ.get("H370_PIN"):
        cand.append(Path(os.environ["H370_PIN"]))
    cand += [CASE / "run_0062_euler_wallfit_fit_r1_ext6k",
             Path("/home/sano/work/forge/case/45.isobutane_m6_d155/run_0062_euler_wallfit_fit_r1_ext6k")]
    return next((c.resolve() for c in cand if (c / "prepare_info.json").exists()), None)


PIN = find_pin()


# ====================================================================================== A・B. 生産の設計 (CFD ピン)
def production():
    """凍結した生産の設計 (h370_problems の NS・Euler の問題をピンで読む)。"""
    import h370_problems as HP
    td = TMP / "pinned"
    td.mkdir(exist_ok=True)
    pn, pe = HP.load_pinned(HP.NS, PIN, td), HP.load_pinned(HP.EU, PIN, td)
    d, de = RA.design_chain(pn), RA.design_chain(pe)
    return HP, pn, pe, d, de


def test_wall_and_delta_in():
    print("== A. Euler の参照の壁 (案 A) と壁ファイル / B. 生成の座標の δ_in の照合", flush=True)
    HP, pn, pe, d, de = production()
    S = de["wall"]
    scale = float(pe.spec["r_throat"])
    wE, er = RA.euler_reference_wall(pe, de, scale)
    g = er["gate"]
    mg = WA.monotone_gate(S._spl, S.x0, S.x_e)
    check(f"案 A: H の上流を付けた設計壁が設計壁のゲートで合格 (S の出口端 r′ 最小 {mg['min_r1']:.4e}・累積減少 {mg['D']:.4e} r_t; "
          f"codex の再現値 −4.0164e−8・6.1293e−9)、上流は単調な収縮・正の半径、一意な大域最小 (x_t {g['global_min']['x_t']:.1e})、C²",
          isinstance(wE, WA.UpstreamDesignWall) and g["pass"] and g["downstream_monotone_gate"]["pass"]
          and abs(mg["min_r1"] - (-4.0164e-8)) < 5e-12 and abs(mg["D"] - 6.1293e-9) < 5e-13
          and g["upstream"]["monotone_contraction"] and g["upstream"]["positive_radius"] and g["global_min"]["pass"]
          and g["max_seam_jump"] <= WA.POLY_SEAM_TOL and er["repr"] == WA.DESIGN_UPSTREAM_KIND and not wE.validate())
    xs = np.r_[0.0, np.linspace(0.0, float(S.x_e), 100001)[1:]]
    check("案 A: x ≥ 0 は設計壁 S とビット同一 (値・1〜3 階)", all(np.array_equal(wE.r(xs, n), S.r(xs, n)) for n in range(4)))
    # 物理壁のゲートは今と同じ (緩めていない)
    check("物理壁のゲートの定数は今と同じ (THROAT_R1_TOL 1e-12・設計壁のゲート r′ ≥ −1.74533e-6・D ≤ 1e-7)",
          WA.THROAT_R1_TOL == 1e-12 and WA.MONO_R1_MIN == -1.74533e-6 and WA.MONO_D_MAX == 1e-7)
    e = raises(lambda: RA.build_physical_wall(pe, de, scale, delta_r_x=RA.delta_r_zero(S.x_in, S.x_e), offset="radial"), ValueError)
    check(f"物理壁のクラスに δ_r = 0 を入れた H の壁は今も物理壁の厳密なゲート (スロートの後で r′ ≥ −1e-12) で落ちる ({(e or '')[:60]}…)",
          e is not None and "後で r′ ≥ 0 False" in e)
    # 負例: S に凹み → 設計壁のゲートで不合格 / 上流が単調でない H → 不合格
    from scipy.interpolate import BSpline
    Sb = copy.copy(S)
    c = np.array(S._spl.c, dtype=float, copy=True)
    c[-4] -= 1e-5                                          # 出口の近く (r′ ≈ 0) の係数を下げて凹みを作る
    Sb._spl = BSpline(S._spl.t, c, S._spl.k)
    pu = {"r_U": 3.2425, "L_U": 9.0, "L_pipe": 3.5}
    e = raises(lambda: WA.UpstreamDesignWall(Sb, pu), ValueError)
    check(f"負例: S の出口の近くに 1e-5 r_t の凹み → 設計壁のゲート (monotone_gate) で不合格 ({(e or '')[:70]}…)", e is not None and "monotone_gate" in e)
    e = raises(lambda: WA.UpstreamDesignWall(S, {"r_U": 1.002, "L_U": 9.0, "L_pipe": 3.5}), ValueError)
    check(f"負例: 上流の直管の半径 1.002 (Q が単調な収縮にならない) → 不合格 ({(e or '')[:70]}…)", e is not None and ("単調な収縮" in e or "大域最小" in e))
    # Euler の prepare: 物理壁の表現の明示は run dir を作る前に例外 (案 A の壁は物理壁ではない)
    raw = yaml.safe_load(HP.EU.read_text())
    raw["geometry"]["initial_line_run"] = str(PIN)
    raw["geometry"]["physical_wall_repr"] = "legacy"
    yb = TMP / "eu_bad_repr.yaml"
    yb.write_text(yaml.safe_dump(raw, allow_unicode=True, sort_keys=False))
    rd = TMP / "never_euler"
    e = raises(lambda: RA.prepare(yb, rd), ValueError)
    check(f"Euler の prepare: H + physical_wall_repr の明示 → run dir を作る前に例外 ({(e or '')[:60]}…)",
          e is not None and "physical_wall_repr" in e and not rd.exists())
    # 壁ファイルの保存 → 復元 → 格子の座標がビット一致
    mpE = RA.mesh_params_euler(pe, scale)
    cE, qE, bE = generate_axisym_mesh(wE, mpE)
    wd = TMP / "wallfile"
    wd.mkdir()
    wpath, wrec, wsha = WA.save_wall_file(wd, wE, scale, WA.DESIGN_UPSTREAM_KIND)
    back = WA.load_wall_file(wpath)
    cB, qB, bB = generate_axisym_mesh(back["physical"], mpE)
    xa = np.linspace(wE.x_in, wE.x_e, 200001)
    check(f"壁ファイル (design_upstream) の保存 → 復元: 壁の値・1〜3 階がビット一致、Euler の格子 (throat_refine {mpE.throat_refine:g}、"
          f"{mpE.ni} × {cE.shape[0] // mpE.ni}) の座標・接続・境界がビット一致",
          back["kind"] == WA.DESIGN_UPSTREAM_KIND and all(np.array_equal(back["physical"].r(xa, n), wE.r(xa, n)) for n in range(4))
          and np.array_equal(cB, cE) and np.array_equal(qB, qE) and all(np.array_equal(np.asarray(bB[k]), np.asarray(bE[k])) for k in bE)
          and mpE.throat_refine == 5.0 and mpE.ni == 2000)
    bad = json.loads(wpath.read_text())
    bad["physical_wall"]["upstream_poly"]["coef"][2] += 1e-6
    (wd / "bad").mkdir()
    (wd / "bad" / WA.WALL_FILE).write_text(json.dumps(bad))
    e = raises(lambda: WA.load_wall_file(wd / "bad"), ValueError)
    check(f"壊れた壁ファイル (Q の係数を 1e-6 ずらす) は復元で例外 ({(e or '')[:50]}…)", e is not None and "継ぎ目" in e)
    e = raises(lambda: WA.wall_file_record(wE, scale, "legacy"), ValueError)
    check("UpstreamDesignWall を物理壁の種類 (legacy) で書こうとすると例外", e is not None)
    # B. 生成の座標の δ_in (NS = Gref H 370 K の問題 + delta_r_initial_C.csv、Euler = 案 A の壁 + mesh_euler)
    tbl = np.loadtxt(HP.CSV, delimiter=",", skiprows=1)
    drx = RA.delta_r_from_table(tbl[:, 0], tbl[:, 1])
    W = RA.build_physical_wall(pn, d, scale, delta_r_x=drx, offset="radial")
    mp = RA.mesh_params(pn, scale, ni=561, nj=97, wall_first_frac=4.5e-5)
    c, q, b = generate_axisym_mesh(W, mp)
    ni, nj = mp.ni, c.shape[0] // mp.ni
    XN, RN = c[:, 0].reshape(ni, nj) / scale, c[:, 1].reshape(ni, nj) / scale
    XE, RE = cE[:, 0].reshape(mpE.ni, -1) / scale, cE[:, 1].reshape(mpE.ni, -1) / scale
    res5 = compare_delta_in(XN[:, 0], RN[:, -1], XE[:, 0], RE[:, -1], drx, scale, r_ns_exact=W.r, r_eu_exact=wE.r)
    mp4 = copy.copy(mpE)
    mp4.throat_refine = 4.0
    c4, _, _ = generate_axisym_mesh(wE, mp4)
    res4 = compare_delta_in(XN[:, 0], RN[:, -1], c4[:, 0].reshape(mp4.ni, -1)[:, 0] / scale, c4[:, 1].reshape(mp4.ni, -1)[:, -1] / scale,
                            drx, scale)
    j5, j4 = res5["judged"], res4["judged"]
    check(f"δ_in の照合 (生成の座標): throat_refine 5 で x ≥ 0 の max|δ_in − δ_r| = {j5['max_abs_mm']:.6e} mm ≤ {DELTA_IN_TOL_MM} mm "
          f"(x {j5['x_at_max']:.4f}; codex の再計算 7.96404e−4)・全列の被覆", res5["pass"] and res5["covered"]
          and j5["max_abs_mm"] <= DELTA_IN_TOL_MM and abs(j5["max_abs_mm"] - 7.96404e-4) < 1e-9)
    check(f"  throat_refine 4 (旧い FAIL の再現) は {j4['max_abs_mm']:.6e} mm > {DELTA_IN_TOL_MM} mm (記録 1.136651e−3)",
          not res4["pass"] and abs(j4["max_abs_mm"] - 1.136651e-3) < 1e-9)
    miss = np.zeros(ni, dtype=bool)
    for jj in range(nj):
        miss |= (XN[:, 0] < XN[0, jj]) | (XN[:, 0] > XN[-1, jj])
    check(f"NS の格子 (Gref H、壁法線の近壁層) で一定 x の断面への補間の端の不足は 0 列 (列内の x の幅の最大 {np.ptp(XN, axis=1).max():.4f} r_t)",
          not miss.any())
    return {"pn": pn, "pe": pe, "d": d, "de": de, "wE": wE, "W": W, "scale": scale, "ns_coords": c, "ni": ni, "nj": nj,
            "eu_coords": cE, "niE": mpE.ni}


# ====================================================================================== C. 停止規則のキーの解決
def test_stop_rule_keys():
    print("== C. 370 K の NS の問題の停止規則のキーの解決", flush=True)
    import h370_problems as HP
    pn = RA.load_problem(HP.NS)
    eu = TMP / "euler_dummy_c"
    eu.mkdir(exist_ok=True)
    sr = RA._stop_rule_settings(pn, None, str(eu), None)
    want = {"requested": True, "use": True, "cap": 200000, "hold_steps": 40000, "out_step_interval": 2500, "fallback_steps": 200000,
            "window": 20000, "min_points": 9, "r_mean_max": 5e-4, "r_amp_max": 5e-3, "theta_drift_max": 5e-4, "theta_amp_max": 5e-3,
            "x_fracs": [0.419738616669543, 0.734542579171700, 0.986385749173426]}
    bad = {k: sr.get(k) for k, v in want.items() if sr.get(k) != v}
    check("NS の問題 (stop_rule・stop_rule_cap 200000・stop_rule_hold_steps 40000・stop_rule_thresholds) → 登録どおりに解決 (既定と同じ値も明示)",
          not bad and sr["cap_source"] == "evaluate.stop_rule_cap"
          and sorted(sr["thresholds_source"]["evaluate.stop_rule_thresholds"]) == sorted(SR.thresholds().keys()), str(bad))
    xw = 95.29740274407894
    check("x_fracs × x_w = 40・70・94 r_t (x_w = 修正後の設計壁の出口 95.2974)",
          all(abs(f * xw - x) < 1e-6 for f, x in zip(sr["x_fracs"], (40.0, 70.0, 94.0))))
    ex = {k: v for k, v in sr["extract"].items() if k != "source"}
    check("抽出の条件 (prepare_info の stop_rule.extract) = 受け渡しの固定の条件 (edge・remap_constant_x=True・端は欠測・ω 0.5・λ 1.0・ノット 2.0・昇格の上限 ×1000)",
          ex == DL.handoff_extract_settings() == {"band_select": "edge", "remap_constant_x": True, "remap_edge": "missing", "omega": 0.5,
                                                   "smooth_lam": 1.0, "knot_spacing": 2.0, "max_lam_factor": 1000.0}, str(ex))
    p3 = RA.load_problem(HP.H300)
    check("停止規則のキーの無い問題 (300 K の H) は今と同じ (None = prepare_info に書かない)", RA._stop_rule_settings(p3, None, str(eu), None) is None)


# ====================================================================================== 合成の構造格子 (D〜G)
S_SC = 0.0766
NI, NJ = 56, 48
X0, X1 = -2.0, 20.0


def rwE_f(x):
    return 1.0 + 0.002 * x * x


def dl_f(x):
    return 0.01 + 0.002 * (x + 2.0)


def core_q(x, r, rw):
    return 100.0 * (1.0 + 0.02 * np.tanh(x / 10.0)) * (1.0 - 0.01 * (r / rw) ** 2)


def radial(n):
    return 1.0 - (1.0 - np.linspace(0.0, 1.0, n)) ** 2


def write_run(rd: Path, X, R, Q, info: dict, steps=(2500,), scale=S_SC, perturb=None):
    """node 構造格子 run (nozzle.h5 の COORD [m]・res_<N>.h5 の ro・Ux・Uy・sonic・prepare_info.json)。perturb(step) → Q の倍率。"""
    rd.mkdir(parents=True)
    ni, nj = X.shape
    with h5py.File(rd / "nozzle.h5", "w") as f:
        f["/MESH/COORD"] = np.c_[X.ravel() * scale, R.ravel() * scale, np.zeros(ni * nj)].ravel()
    for st in steps:
        k = 1.0 if perturb is None else perturb(st)
        with h5py.File(rd / f"res_{st}.h5", "w") as f:
            f["/VALUE/ro"] = np.ones(ni * nj)
            f["/VALUE/Ux"] = (Q * k).ravel()
            f["/VALUE/Uy"] = np.zeros(ni * nj)
            f["/VALUE/sonic"] = np.full(ni * nj, 300.0)
    (rd / "prepare_info.json").write_text(json.dumps({"mesh": {"ni": ni, "nj": nj}, "scale_m": scale, "x_E": 15.0, **info}))


def synth_ns(sign: float = 1.0):
    """NS: 壁 r = rwE + δ_r、壁法線に沿わせた近壁層 (sign = +1 は mesher と同じ向き、−1 は逆向き = 両端の列で行の x が範囲を覆わない)。"""
    x = np.linspace(X0, X1, NI)
    rw = rwE_f(x) + dl_f(x)
    R = rw[:, None] * radial(NJ)[None, :]
    d = rw[:, None] - R
    slope = (0.004 * x + 0.002) / np.sqrt(1.0 + (0.004 * x + 0.002) ** 2)
    beta = np.clip(1.0 - d / (0.1 * rw[:, None]), 0.0, 1.0)
    X = x[:, None] + sign * beta * d * slope[:, None]
    y = rw[:, None] - R
    Q = core_q(X, R, rw[:, None]) * np.minimum(1.0, np.maximum(y / (5.0 * dl_f(X)), 0.0) ** (1.0 / 7.0))
    return X, R, Q


def synth_euler():
    x = np.linspace(X0 - 0.5, X1 + 0.5, 70)
    rw = rwE_f(x)
    R = rw[:, None] * radial(NJ)[None, :]
    X = np.repeat(x[:, None], NJ, axis=1)
    return X, R, core_q(X, R, rw[:, None])


def ns_info(euler_ref):
    sr = SR.resolve_settings({"stop_rule": True}, nsteps_arg=None, n_main=200000, n_main_source="試験", out_int=2500, euler_ref=euler_ref)
    return {"wall_thermal": {"mode": "isothermal", "Tw": 370.0}, "physical_wall": {"domain": [X0, X1]}, "throat_physical": {"r": 1.012},
            "euler_ref": SR.euler_ref_spec(euler_ref) if isinstance(euler_ref, dict) else str(euler_ref), "stop_rule": sr}


def make_euler(name="eu"):
    """Euler の参照 run (res_1000 と res_2000 を別の場で; 固定の試験用)。戻り (run dir, 固定の指定 res_1000)。"""
    X, R, Q = synth_euler()
    rd = TMP / name
    write_run(rd, X, R, Q, {}, steps=(1000, 2000), perturb=lambda st: 1.0 if st == 1000 else 1.003)
    return rd, {"run": str(rd), "res": "res_1000.h5", "sha256": DS._sha256_file(rd / "res_1000.h5")}


def test_extract_and_pin():
    print("== D. 抽出の一定 x の断面への補間・端の欠測・Euler の参照の固定", flush=True)
    eu, pin = make_euler("eu_d")
    Xo, Ro, Qo = synth_ns(+1.0)
    Xb, Rb, Qb = synth_ns(-1.0)
    write_run(TMP / "ns_ok_d", Xo, Ro, Qo, ns_info(pin))
    write_run(TMP / "ns_bad_d", Xb, Rb, Qb, ns_info(pin))
    kw = dict(band_select="edge", remap_constant_x=True)
    a = DS.deltastar_from_core_matched_euler(TMP / "ns_bad_d", pin, **kw)
    b = DS.deltastar_from_core_matched_euler(TMP / "ns_bad_d", pin, remap_edge="missing", **kw)
    ne = b["not_extracted"]
    check("既定 (remap_edge='clamp') は今と同じ: 端の補間範囲の不足は端値で埋めて全列を抽出 (not_extracted 0)",
          a["not_extracted"]["n"] == 0 and len(a["x"]) == NI and a["settings"].get("remap_clamped_points", 0) > 0)
    check(f"remap_edge='missing': 端の補間範囲の足りない列 (入口・出口) を欠測にして抽出しない (not_extracted {ne['n']}、理由 {sorted(set(ne['reason']))})",
          ne["n"] == 2 and set(ne["reason"]) == {"remap_missing"} and np.allclose(ne["x"], [X0, X1], rtol=0.0, atol=1e-12) and len(b["x"]) == NI - 2
          and np.array_equal(b["x"], a["x"][1:-1]))
    c_ = DS.deltastar_from_core_matched_euler(TMP / "ns_ok_d", pin, remap_edge="missing", **kw)
    check("壁法線の向きが mesher と同じ格子は欠測 0 列 (全列を抽出)", c_["not_extracted"]["n"] == 0 and len(c_["x"]) == NI)
    # Euler の参照の固定
    check("固定の指定 {run, res: res_1000.h5, sha256} → その res を読む (最新の res_2000 ではない)・記録に残す",
          c_["settings"]["euler_res"] == "res_1000.h5" and c_["settings"]["euler_ref_pinned"]["res"] == "res_1000.h5")
    d_ = DS.deltastar_from_core_matched_euler(TMP / "ns_ok_d", str(eu), remap_edge="missing", **kw)
    check("ディレクトリだけの指定は今と同じ (最新の res_2000)", d_["settings"]["euler_res"] == "res_2000.h5" and "euler_ref_pinned" not in d_["settings"]
          and d_["settings"]["euler_run"] == str(eu))
    wrong = {**pin, "sha256": "0" * 64}
    e = raises(lambda: DS.deltastar_from_core_matched_euler(TMP / "ns_ok_d", wrong, **kw), ValueError)
    check(f"sha の不一致 → 抽出器が例外 ({(e or '')[:50]}…)", e is not None and "sha256" in e)
    e = raises(lambda: DS.massflow_ratio(TMP / "ns_ok_d", wrong), ValueError)
    check("sha の不一致 → 流量の帳簿 (massflow_ratio) も例外", e is not None)
    e = raises(lambda: DS.resolve_euler_ref({**pin, "res": "../res_1000.h5"}), ValueError)
    check("res にディレクトリを含む指定は例外 (ファイル名だけ)", e is not None)
    P = lambda ev: type("P", (), {"evaluate": ev})()  # noqa: E731
    ev = {"stop_rule": True, "nStepOuter": 200000, "outStepInterval": 2500}
    e = raises(lambda: RA._stop_rule_settings(P(ev), None, wrong, None), ValueError)
    check(f"prepare_ns の停止規則の設定: 固定の sha の不一致 → run dir を作る前に例外 ({(e or '')[:50]}…)", e is not None and "sha256" in e)
    sr_ok = RA._stop_rule_settings(P(ev), None, pin, None)
    check("prepare_ns の停止規則の設定: 正しい固定 → euler_ref に {run, res, sha256} を残す", sr_ok["euler_ref"] == pin and sr_ok["use"])
    rdw = TMP / "ns_wrong_pin"
    shutil.copytree(TMP / "ns_ok_d", rdw)
    info = json.loads((rdw / "prepare_info.json").read_text())
    info["stop_rule"]["euler_ref"] = wrong
    (rdw / "prepare_info.json").write_text(json.dumps(info))
    e = raises(lambda: RA._stop_rule_for_main(rdw, None), RuntimeError)
    check(f"本段の起動前: 固定の sha の不一致 → forge を起動する前に RuntimeError ({(e or '')[:50]}…)", e is not None and "固定" in e)
    e = raises(lambda: DL.main(["--problem", "x.yaml", "--euler-ref", str(eu), "--euler-ref-res", "res_1000.h5", "--euler-ref-sha256", "0" * 64,
                                "--run-dir", str(TMP / "never_cli"), "--delta-r-csv", str(TMP / "never.csv")]), ValueError)
    check("CLI: --euler-ref-res・--euler-ref-sha256 の不一致 → 例外 (run dir を作らない)", e is not None and not (TMP / "never_cli").exists())
    e = raises(lambda: DL.main(["--problem", "x.yaml", "--euler-ref", str(eu), "--euler-ref-res", "res_1000.h5", "--run-dir", str(TMP / "never_cli")]),
               SystemExit)
    check("CLI: --euler-ref-res だけ (sha256 なし) → 使い方の誤り", e is not None)
    # remap_constant_x=True が帯の生成 (fix_band) に渡る
    calls = []
    real = DL.extract_and_merge

    def spy(*a, **k):
        calls.append(k)
        return real(*a, **k)
    DL.extract_and_merge = spy
    try:
        band = SR.fix_band(TMP / "ns_ok_d", 2500, pin, TMP / "band_d", **SR.band_extract_kw(ns_info(pin)["stop_rule"]))
    finally:
        DL.extract_and_merge = real
    check("帯の生成 (fix_band): band_select='edge'・remap_constant_x=True・remap_edge='missing'・固定の Euler の参照で抽出",
          len(calls) == 1 and calls[0]["band_select"] == "edge" and calls[0]["remap_constant_x"] is True and calls[0]["remap_edge"] == "missing"
          and band["euler_ref"]["res"] == "res_1000.h5" and len(band["x"]) == NI, str(calls))
    e = raises(lambda: SR.fix_band(TMP / "ns_ok_d", 2500, wrong, TMP / "band_d2", **SR.band_extract_kw(ns_info(pin)["stop_rule"])), SR.Unavailable)
    check("帯の生成: 固定の sha の不一致 → Unavailable (停止規則を使えない)", e is not None and "sha256" in e)
    return eu, pin


def test_recorder(pin):
    print("== E. 見張りの中の出力ごとの δ_E・δ_target の記録", flush=True)
    Xo, Ro, Qo = synth_ns(+1.0)
    steps = (2500, 5000, 7500)
    rd = TMP / "ns_rec"
    write_run(rd, Xo, Ro, Qo, ns_info(pin), steps=steps, perturb=lambda st: 1.0 + 1e-4 * st / 2500)
    uc = SR.use_columns(rd)
    ex = {**DL.handoff_extract_settings(), "source": "試験"}
    calls = []
    real = DL.extract_and_merge

    def spy(*a, **k):
        calls.append(k)
        return real(*a, **k)
    DL.extract_and_merge = spy
    try:
        t0 = time.time()
        one = SR.delta_snapshot(rd, 2500, pin, ex, uc["x"])
        t1 = time.time() - t0
    finally:
        DL.extract_and_merge = real
    xin = np.linspace(X0, X1, NI)
    check(f"delta_snapshot: 使用域 [0, x_w] の全列 ({len(uc['x'])} 列) の δ_E・δ_target が有限、欠測 0 (合成の格子で {t1:.1f} s)",
          len(one["delta_E"]) == len(uc["x"]) == int(np.sum(xin >= 0.0)) and np.all(np.isfinite(one["delta_E"]))
          and np.all(np.isfinite(one["delta_target"])) and one["n_missing"] == 0 and one["res_sha256"] == DS._sha256_file(rd / "res_2500.h5"))
    check("delta_snapshot: 受け渡しと同じ条件 (edge・remap_constant_x=True・端は欠測・ω 0.5・λ 1.0・昇格の上限 ×1000) で抽出",
          len(calls) == 1 and {k: calls[0][k] for k in DL.handoff_extract_settings()} == DL.handoff_extract_settings(), str(calls))
    rec = SR.DeltaRecorder(rd, pin, ex, uc["x"], span=40000)
    for st in steps:
        rec.submit(st)
    got = rec.finish(300.0)
    check("DeltaRecorder: 出力ごとに抽出して結果を返す (3 出力、どれも δ_E・δ_target・held)",
          not got.pop("_unfinished") and sorted(got) == list(steps) and all("delta_E" in got[s] and "delta_target" in got[s] for s in steps)
          and got[2500]["delta_E"] == one["delta_E"])
    rec = SR.DeltaRecorder(rd, pin, ex, uc["x"], span=2000)
    real_snap = SR.delta_snapshot
    SR.delta_snapshot = lambda *a, **k: (time.sleep(0.5), {"delta_E": [0.0]})[1]
    try:
        for st in steps:
            rec.submit(st)
        got = rec.finish(30.0)
    finally:
        SR.delta_snapshot = real_snap
    got.pop("_unfinished")
    check("DeltaRecorder: 最新の出力から 2 窓 (span) より古い出力は抽出せず skipped (追いつかないとき)",
          "skipped" in got.get(5000, {}) and "delta_E" in got.get(7500, {}), str({k: list(v) for k, v in got.items()}))
    rows = [{"step": s} for s in steps]
    SR._attach_deltas(rows, {2500: {"delta_E": [1.0]}, 9999: {}})
    check("_attach_deltas: 同じ step の行の delta に付ける", rows[0].get("delta") == {"delta_E": [1.0]} and "delta" not in rows[1])


# ====================================================================================== F. 受け渡しの関門
def stub_conv(kind="PASS"):
    return lambda run_dir: {"cmd": ["check_convergence.py", "--segment", "--drop", "3", str(run_dir)], "rc": 0 if kind == "PASS" else 1,
                            "kind": kind, "verdict": [f"OVERALL: {kind} (試験の偽物)"], "output": f"OVERALL: {kind} (試験の偽物)\n"}


def stub_wr(kind="PASS"):
    return lambda run_dir, step=None: {"cmd": ["check_wall_resolution.py"], "rc": {"PASS": 0, "FAIL": 1}.get(kind, 2), "kind": kind,
                                       "verdict": [f"VERDICT: {kind}"], "output": f"VERDICT: {kind} (試験の偽物)\n"}


def stub_nw(ok=True):
    return lambda problem, delta_next: ({"pass": True, "wall": "試験の偽物"} if ok else {"pass": False, "why": "試験の偽物: 次の壁が非単調"})


def gate_run(name, pin, ns=(+1.0), r_fn=None, delta_fn=None, s_stop=62500, steps=None, status="REACHED", euler_ref=None):
    """合成の受け渡しの run: 場 (res_<s>.h5、本物の抽出に使う)・prepare_info (停止規則の設定)・STOP_RULE_STATE.json・stop_rule.json
    (行の r と出力ごとの δ の記録)。r_fn(step) → r、delta_fn(step, x) → (δ_E, δ_target) [r_t] (None の行は記録なし)。"""
    X, R, Q = synth_ns(ns)
    rd = TMP / name
    info = ns_info(pin if euler_ref is None else euler_ref)
    write_run(rd, X, R, Q, info, steps=(s_stop,))
    sr = info["stop_rule"]
    uc = SR.use_columns(rd)
    steps = steps or list(range(2500, s_stop + 1, 2500))
    rows = []
    for st in steps:
        row = {"step": st, "nonfinite_all": 0, "r": float(r_fn(st)) if r_fn else 0.0}
        dd = (lambda s_, x_: (np.full(len(x_), 1e-3), np.full(len(x_), 2e-3))) if delta_fn is None else delta_fn
        v = dd(st, uc["x"])
        if v is not None:
            row["delta"] = {"delta_E": list(map(float, v[0])), "delta_target": list(map(float, v[1])), "held": [0] * len(uc["x"]),
                            "missing": [], "n_missing": 0, "n_held": 0, "seconds": 0.0}
        rows.append(row)
    rec = {"version": 1, "settings": sr, "status": status, "rows": rows, "stop_step": s_stop,
           "delta_series": {"x": uc["x"].tolist(), "use_domain": [0.0, uc["x_w"]], "scale_m": uc["scale_m"], "extract": sr["extract"],
                            "euler_ref": sr["euler_ref"], "span_steps": 40000}}
    SR.save_record(rd, rec)
    SR.write_state(rd, {"status": status, "stop_step": s_stop, "forge_rc": 143})
    return rd


def test_gate(eu, pin):
    print("== F. 受け渡しの関門 (§4.9 の更新に進む条件、上書きでは通さない)", flush=True)
    real = {k: getattr(DL, k) for k in ("convergence_verdict", "wall_resolution_verdict", "next_wall_check", "extract_and_merge")}
    calls = []

    def spy(*a, **k):
        calls.append(k)
        return real["extract_and_merge"](*a, **k)

    def run(rd, conv="PASS", wr="PASS", nw=True, **kw):
        DL.convergence_verdict, DL.wall_resolution_verdict, DL.next_wall_check = stub_conv(conv), stub_wr(wr), stub_nw(nw)
        calls.clear()
        return DL.handoff_gate(rd, problem="試験.yaml", **kw)
    DL.extract_and_merge = spy
    try:
        rd = gate_run("g_ok", pin)
        g = run(rd, euler_ref=pin)
        saved = json.loads((rd / DL.HANDOFF_GATE).read_text())
        check("全部の条件を満たす → 進む (抽出は edge・remap_constant_x=True・端は欠測・ω 0.5・λ 1.0・昇格の上限 ×1000、固定の Euler の参照)",
              g["proceed"] and not g["blocked_by"] and len(calls) == 1 and calls[0]["band_select"] == "edge" and calls[0]["remap_constant_x"] is True
              and calls[0]["remap_edge"] == "missing" and calls[0]["omega"] == 0.5 and calls[0]["smooth_lam"] == 1.0
              and calls[0]["max_lam_factor"] == 1000.0 and g["checks"]["euler_ref"]["pinned"], f"{g['blocked_by']}")
        check("  HANDOFF_GATE.json に各条件・抽出の条件 (λ・昇格の上限・ω) を記録、表を run の直下に写す (handoff_extract にも残る)",
              saved["proceed"] and saved["extract"]["settings"]["smooth_lam"] == 1.0 and saved["extract"]["settings"]["max_lam_factor"] == 1000.0
              and saved["extract"]["settings"]["omega"] == 0.5 and all(saved["checks"][c]["pass"] for c in
                                                                    ("conservation", "wall_resolution", "two_windows", "coverage", "next_wall", "euler_ref", "omega"))
              and (rd / "delta_r_next.csv").exists() and (rd / DL.HANDOFF_EXTRACT_DIR / "delta_r_next.csv").exists()
              and (rd / DL.HANDOFF_WALL_RESOLUTION).exists() and "delta_r_next.csv" in saved["extract"]["made"])
        cases = [
            ("保存則: r の窓の平均 1e-3 (> 5e-4)", dict(r_fn=lambda s: 1e-3), {}, "保存則"),
            ("保存則: r の振幅 1 % (> 0.5 %)", dict(r_fn=lambda s: 1e-2 * (-1) ** (s // 2500)), {}, "保存則"),
            ("局所 y₁⁺: check_wall_resolution が FAIL", {}, dict(wr="FAIL"), "局所 y₁⁺"),
            ("局所 y₁⁺: 判定不能 (壁のダンプが無い)", {}, dict(wr="UNDETERMINED"), "局所 y₁⁺"),
            ("被覆: 端の補間範囲の不足 (入口・出口の列が欠測)", dict(ns=-1.0), {}, "被覆"),
            ("2 窓: 後の窓の中で δ_E が 0.013 mm 振れる",
             dict(delta_fn=lambda s, x: (np.full(len(x), 1e-3 + (1.7e-4 if s == 52500 else 0.0)), np.full(len(x), 2e-3))), {}, "2 窓"),
            ("2 窓: 後の窓の中で δ_target が 0.0115 mm 振れる",
             dict(delta_fn=lambda s, x: (np.full(len(x), 1e-3), np.full(len(x), 2e-3 + (1.5e-4 if s > 42500 else 0.0)))), {}, "2 窓"),
            ("2 窓: 窓の中の出力の δ の記録が無い (抽出の失敗)", dict(delta_fn=lambda s, x: None if s == 30000 else (np.full(len(x), 1e-3),) * 2), {}, "2 窓"),
            ("2 窓: 履歴が 2 窓を覆わない (停止 30,000)", dict(s_stop=30000), {}, "2 窓"),
            ("次の壁の幾何の健全性が不合格", {}, dict(nw=False), "次の壁"),
        ]
        for lab, mk, rkw, key in cases:
            rd = gate_run("g_" + str(abs(hash(lab)) % 10 ** 8), pin, **mk)
            g = run(rd, euler_ref=pin, **rkw)
            check(f"負例 {lab} → 止める ({'; '.join(g['blocked_by'])[:110]}…)", not g["proceed"] and any(key in b for b in g["blocked_by"])
                  and not (rd / "delta_r_next.csv").exists())
        # 2 窓: 各窓の中の振れは 0.0095 mm (≤ 0.01) だが窓の平均どうしの差 0.0169 mm (共有の点 s − W = 42,500 を挟んで段差)
        r_t_mm = S_SC * 1e3
        a1, a2 = 0.0095 / r_t_mm, 0.019 / r_t_mm
        rd = gate_run("g_meandiff", pin, delta_fn=lambda s, x: (np.full(len(x), 1e-3),
                                                                np.full(len(x), 2e-3 + (0.0 if s < 42500 else (a1 if s == 42500 else a2)))))
        g = run(rd, euler_ref=pin)
        w = g["checks"]["two_windows"].get("worst", {})
        check(f"負例 2 窓: 窓の中の振れ {w.get('delta_target.range_W1', {}).get('max_mm', float('nan')):.4f}・"
              f"{w.get('delta_target.range_W2', {}).get('max_mm', float('nan')):.4f} mm (≤ 0.01) でも窓の平均どうしの差 "
              f"{w.get('delta_target.mean_diff', {}).get('max_mm', float('nan')):.4f} mm > 0.01 → 止める",
              not g["proceed"] and w["delta_target.range_W1"]["max_mm"] <= 0.01 and w["delta_target.range_W2"]["max_mm"] <= 0.01
              and w["delta_target.mean_diff"]["max_mm"] > 0.01 and any("2 窓" in b for b in g["blocked_by"]))
        # 使用域の held: 抽出の表の held を立てる (偽の抽出)
        rd = gate_run("g_held", pin)

        def held_extract(run_dir, euler_run, out_dir=None, **k):
            s_ = real["extract_and_merge"](run_dir, euler_run, out_dir=out_dir, **k)
            nx = DL.read_delta_r_next(Path(out_dir) / "delta_r_next.csv")
            h = np.asarray(nx["held"]).copy()
            h[np.argmin(np.abs(nx["x_rt"] - 10.0))] = 1
            np.savetxt(Path(out_dir) / "delta_r_next.csv", np.c_[nx["x_rt"], nx["delta_target"], nx["delta_in_prev"], nx["delta_E"], h],
                       delimiter=",", header=",".join(DL.DELTA_R_NEXT_COLUMNS), comments="")
            return s_
        DL.extract_and_merge = held_extract
        g = run(rd, euler_ref=pin)
        DL.extract_and_merge = spy
        check(f"負例 使用域の held 1 点 (x ≈ 10) → 止める ({'; '.join(g['blocked_by'])[:80]}…)",
              not g["proceed"] and g["checks"]["coverage"]["n_held_use"] == 1 and not (rd / "delta_r_next.csv").exists())
        rd = gate_run("g_eref", pin)
        g = run(rd, euler_ref=str(eu))
        check("負例 呼び出し側の Euler の参照が本段の記録 (固定の res_1000) と違う (ディレクトリだけ) → 止める",
              not g["proceed"] and not g["checks"]["euler_ref"]["pass"] and not calls)
        g = run(rd, euler_ref=pin, omega=0.4)
        check("負例 ω 0.4 (固定の 0.5 と違う) → 止める", not g["proceed"] and not g["checks"]["omega"]["pass"] and not calls)
        wrong = {**pin, "sha256": "f" * 64}
        rd = gate_run("g_sha", pin, euler_ref=wrong)
        g = run(rd)
        check(f"負例 本段の記録の Euler の参照の sha が合わない → 抽出の例外で止める ({'; '.join(g['blocked_by'])[:70]}…)",
              not g["proceed"] and any("抽出が例外" in b for b in g["blocked_by"]))
        # 上書きの 2 つでは足した条件を通さない
        rd = gate_run("g_override", pin, r_fn=lambda s: 1e-3)
        g = run(rd, conv="NOT_CONVERGED", euler_ref=pin, allow_not_reached=True, allow_unconverged=True)
        check("上書き (allow_not_reached・allow_unconverged) は収束の判定だけを通し、保存則の不合格は通さない",
              not g["proceed"] and any("保存則" in b for b in g["blocked_by"]) and not any("収束の判定" in b for b in g["blocked_by"]))
        rd = gate_run("g_override2", pin, status="CENSORED")
        g = run(rd, wr="FAIL", euler_ref=pin, allow_not_reached=True, allow_unconverged=True)
        check("上書きは CENSORED を通すが、局所 y₁⁺ の FAIL は通さない",
              not g["proceed"] and any("局所 y₁⁺" in b for b in g["blocked_by"]) and not any("停止の状態" in b for b in g["blocked_by"]))
    finally:
        for k, v in real.items():
            setattr(DL, k, v)


# ====================================================================================== G. 次の壁の幾何 (本物の構築)
def test_next_wall():
    print("== G. 次の壁の幾何の健全性 (next_wall_check、本物の壁の構築)", flush=True)
    prob = CASE / "problem_d155_ns_finemesh_recal_final_mono.yaml"
    if not prob.exists():
        skip("G. 次の壁の幾何", f"{prob.name} が無い")
        return
    raw = yaml.safe_load(prob.read_text())
    g = raw["geometry"]
    g["initial_line"] = "hall"
    for k in ("initial_line_run", "initial_line_res", "pw_ramp", "pw_upstream"):
        g.pop(k, None)
    y = TMP / "g_hall.yaml"
    y.write_text(yaml.safe_dump(raw, allow_unicode=True, sort_keys=False))
    p = RA.load_problem(y)
    d = RA.design_chain(p)
    _, drx, _ = RA.integral_delta_r(p, d, p.raw["deltastar_initializer"])
    xt = np.linspace(float(d["wall"].x_in) + 1e-6, float(d["wall"].x_e), 1500)
    ok = TMP / "next_ok.csv"
    np.savetxt(ok, np.c_[xt, drx(xt), drx(xt), drx(xt), np.zeros_like(xt)], delimiter=",", header=",".join(DL.DELTA_R_NEXT_COLUMNS), comments="")
    r_ok = DL.next_wall_check(y, ok)
    check(f"積分法の δ_r の次の壁 → 合格 ({r_ok.get('wall')}、物理スロート x {r_ok.get('throat', {}).get('x', float('nan')):.4f})", r_ok["pass"], str(r_ok.get("why")))
    dip = drx(xt) - 5e-3 * np.exp(-((xt - 90.0) / 0.5) ** 2)          # 出口の近く (S′ が小さい) に深さ 5e-3 r_t の凹み
    bad = TMP / "next_bad.csv"
    np.savetxt(bad, np.c_[xt, dip, dip, dip, np.zeros_like(xt)], delimiter=",", header=",".join(DL.DELTA_R_NEXT_COLUMNS), comments="")
    r_bad = DL.next_wall_check(y, bad)
    check(f"x = 90 に 5e-3 r_t の凹みを入れた δ の次の壁 → 不合格 ({(r_bad.get('why') or '')[:80]}…)", not r_bad["pass"])
    check("問題 YAML が無い → 不合格 (次の壁を作れない)", not DL.next_wall_check(None, ok)["pass"])


# ====================================================================================== H. 1 出力あたりの抽出の時間 (実寸)
def test_timing(prod):
    print("== H. 1 出力あたりの抽出の時間 (実寸の格子: NS 4719 × 160・Euler 2000 × 97、合成の場)", flush=True)
    sc, ni, nj, niE = prod["scale"], prod["ni"], prod["nj"], prod["niE"]
    cN, cE = prod["ns_coords"], prod["eu_coords"]
    XN, RN = cN[:, 0].reshape(ni, nj) / sc, cN[:, 1].reshape(ni, nj) / sc
    XE, RE = cE[:, 0].reshape(niE, -1) / sc, cE[:, 1].reshape(niE, -1) / sc
    rwN = RN[:, -1:]
    yN = rwN - RN
    dl = 0.02 + 0.003 * np.maximum(XN, 0.0) ** 0.8
    QN = core_q(XN, RN, rwN) * np.minimum(1.0, np.maximum(yN / dl, 0.0) ** (1.0 / 7.0))
    QE = core_q(XE, RE, RE[:, -1:])
    eu = TMP / "eu_full"
    write_run(eu, XE, RE, QE, {}, steps=(1000,), scale=sc)
    pin = {"run": str(eu), "res": "res_1000.h5", "sha256": DS._sha256_file(eu / "res_1000.h5")}
    xw = float(prod["W"].x_e)
    ns = TMP / "ns_full"
    info = {"wall_thermal": {"mode": "isothermal", "Tw": 370.0}, "physical_wall": {"domain": [float(prod["W"].x_in), xw]},
            "throat_physical": {"r": float(prod["W"].r_throat)}}
    write_run(ns, XN, RN, QN, {**info, "x_E": 39.0}, steps=(2500,), scale=sc)
    uc = SR.use_columns(ns)
    t0 = time.time()
    rec = SR.delta_snapshot(ns, 2500, pin, {**DL.handoff_extract_settings(), "source": "試験"}, uc["x"])
    dt = time.time() - t0
    interval = 2500 * 0.073
    print(f"    1 出力あたりの抽出の時間 = {dt:.1f} s (delta_snapshot の記録 {rec['seconds']} s)、使用域の列 {len(uc['x'])}、欠測 {rec['n_missing']}、"
          f"held {rec['n_held']}、λ {rec['lam']}; 出力の間隔 2,500 step × 73〜129 ms/step (B0 の実績) = {interval:.0f}〜{2500 * 0.129:.0f} s",
          flush=True)
    check(f"実寸の格子で 1 出力の抽出 {dt:.1f} s < 出力の間隔の下の見積もり {interval:.0f} s (別スレッドで回し、停止の判定は待たせない)・"
          f"使用域の欠測 0 (生産の NS の格子で一定 x の断面への補間の端の不足が無い)",
          dt < interval and rec["n_missing"] == 0 and len(rec["delta_E"]) == len(uc["x"]) and np.all(np.isfinite(rec["delta_E"])))
    return dt


def main() -> int:
    print(f"出力: {TMP}  ピン: {PIN}", flush=True)
    prod = None
    if PIN is None:
        skip("A・B・H (生産の設計)", "CFD ピンの凍結元 run_0062 が見つからない (--pin・H370_PIN で渡す)")
    else:
        prod = test_wall_and_delta_in()
    test_stop_rule_keys()
    eu, pin = test_extract_and_pin()
    test_recorder(pin)
    test_gate(eu, pin)
    test_next_wall()
    if prod is not None and not QUICK:
        test_timing(prod)
    elif QUICK:
        skip("H. 実寸の抽出の時間", "--quick")
    print(f"\n{'ALL PASS' if FAIL == 0 and not SKIPPED else (f'{FAIL} FAIL' if FAIL else 'PASS (飛ばしあり: ' + ', '.join(SKIPPED) + ')')}")
    if FAIL == 0:
        shutil.rmtree(TMP, ignore_errors=True)
    return 1 if FAIL else (2 if SKIPPED else 0)


if __name__ == "__main__":
    sys.exit(main())
