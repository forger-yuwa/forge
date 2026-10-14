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
     落ち着きの条件の上書き (§5.1 #36h): `allow_unsettled` (理由の文字列) は保存則・2 窓の**実測の超過**だけを通し、settled: false・
     unsettled_overrides (理由・r の平均と振幅と窓・2 窓の変動の最大とその量・窓・x)・δ_E の作り方 (最新の出力 1 枚) を HANDOFF_GATE.json と
     抽出の要約 (run の直下・handoff_extract) に残す。上書き 3 つを全部渡しても通さない: 2 窓の判定不能 (記録なし・行の欠け・覆わない・
     δ の記録の held・欠測)・保存則の判定不能 (非有限・点数の不足)・停止の状態 DIVERGED・y₁⁺ (FAIL・判定不能)・被覆 (端の欠測)・使用域の
     held・次の壁・Euler の参照・ω・sha の不一致。上書きを渡しても要らなかった (全部満たす) run は settled: true。呼び出し側の作り方の口
     (`extractor` + `delta_e_source`): 作り方を記録する、片方だけ・不正な記録は ValueError、その抽出にも被覆・held の条件が掛かる。
  G. 次の壁の幾何 (`next_wall_check`、本物の壁の構築): 積分法の δ_r で合格、出口側に凹みを入れた δ で不合格。
  H. 1 出力あたりの抽出の時間 (`delta_snapshot`) を実寸の格子 (NS 4719 × 160・Euler 2000 × 97) で測る (--quick で省く)。
  I. 出口の列の被覆 (§5.1 #36g、B の dry-run の発見): Euler の補間できる範囲 (`deltastar.euler_cover_range`) の幾何の許容 (端の外側だけ、
     端の列の間隔 × 1e−4)。Euler の端が NS の最後の列より 1.68e−9 r_t 短い → 覆う・端の列の値で評価 (抽出器・δ_in の照合・関門の被覆・
     δ の記録・帯の生成)、許容のちょうど外 → 欠測、0.07 r_t の不足 → 欠測、内側の列には効かない (ビット同一)、既定 (clamp) は ±1e−9 のまま。
     2 窓に入るべき出力の列挙 (`stop_rule.required_window_steps`)。

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
            ("2 窓: 記録の行そのものが欠けた出力 47,500 (入るべき出力を出力の間隔から列挙)",
             dict(steps=[k for k in range(2500, 62501, 2500) if k != 47500]), {}, "記録の行が無い出力"),
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
        # 上書きの 2 つでは足した条件を通さない (#36h から上書きは理由の文字列)
        rd = gate_run("g_override", pin, r_fn=lambda s: 1e-3)
        g = run(rd, conv="NOT_CONVERGED", euler_ref=pin, allow_not_reached="試験: 上書き", allow_unconverged="試験: 上書き")
        check("上書き (allow_not_reached・allow_unconverged) は収束の判定だけを通し、保存則の不合格は通さない",
              not g["proceed"] and any("保存則" in b for b in g["blocked_by"]) and not any("収束の判定" in b for b in g["blocked_by"]))
        rd = gate_run("g_override2", pin, status="CENSORED")
        g = run(rd, wr="FAIL", euler_ref=pin, allow_not_reached="試験: 上書き", allow_unconverged="試験: 上書き")
        check("上書きは CENSORED を通すが、局所 y₁⁺ の FAIL は通さない",
              not g["proceed"] and any("局所 y₁⁺" in b for b in g["blocked_by"]) and not any("停止の状態" in b for b in g["blocked_by"]))
        test_gate_unsettled(eu, pin, run, held_extract, real, spy)
    finally:
        for k, v in real.items():
            setattr(DL, k, v)


def test_gate_unsettled(eu, pin, run, held_extract, real, spy):
    """F の続き: 落ち着きの条件の上書き allow_unsettled (§5.1 #36h、forge-f5 と合意)。"""
    print("== F2. 落ち着きの条件の上書き (allow_unsettled、§5.1 #36h)", flush=True)
    US = "試験: 実用として 5 枚平均で区切る (case/70 相当)"
    ALL = dict(allow_not_reached="試験: 未到達", allow_unconverged="試験: 未収束", allow_unsettled=US)
    mm = S_SC * 1e3
    jitter = lambda s, x: (np.full(len(x), 1e-3 + (1.7e-4 if s == 52500 else 0.0)), np.full(len(x), 2e-3))  # noqa: E731

    def marks(rd):
        """run の直下と handoff_extract の抽出の要約の印 (無ければ None)。"""
        out = []
        for f in (rd / "delta_r_extract_summary.json", rd / DL.HANDOFF_EXTRACT_DIR / "delta_r_extract_summary.json"):
            out.append(json.loads(f.read_text()) if f.exists() else None)
        return out

    # --- 正例: 保存則の実測の超過 (r の窓の平均 1e-3 > 5e-4) ---
    rd = gate_run("u_cons", pin, r_fn=lambda s: 1e-3)
    g = run(rd, euler_ref=pin, allow_unsettled=US)
    saved = json.loads((rd / DL.HANDOFF_GATE).read_text())
    uo = g["unsettled_overrides"]
    m_top, m_ext = marks(rd)
    check("allow_unsettled (理由つき): 保存則の実測の超過 (r の平均 1e-3) を通す → 進む・表を run の直下に写す",
          g["proceed"] and not g["blocked_by"] and (rd / "delta_r_next.csv").exists() and g["override"]["applied"] == ["allow_unsettled"],
          f"{g['blocked_by']}")
    check("  記録: settled: false・unsettled_overrides に保存則の実測 (r の平均 1e-3・振幅 0・窓 [42,500, 62,500]・閾値・点数 9・理由)",
          saved["settled"] is False and len(uo) == 1 and uo[0]["override"] == "allow_unsettled" and uo[0]["condition"] == "conservation"
          and abs(uo[0]["r_mean"] - 1e-3) < 1e-12 and uo[0]["r_amp"] == 0.0 and uo[0]["window"] == [42500, 62500]
          and uo[0]["r_mean_max"] == 5e-4 and uo[0]["r_amp_max"] == 5e-3 and uo[0]["n"] == 9 and uo[0]["reason"] == US
          and saved["override"]["allow_unsettled_reason"] == US and saved["unsettled_overrides"] == uo, f"{uo}")
    des = saved["delta_e_source"]
    check("  δ_E の作り方: 最新の出力 1 枚 (kind single_snapshot・res_62500.h5・停止の出力と同じ・関門の抽出)",
          des["kind"] == "single_snapshot" and des["n_snapshots"] == 1 and des["res"] == "res_62500.h5" and des["steps"] == [62500]
          and des["is_stop_output"] is True and "extract_and_merge" in des["provider"], f"{des}")
    check("  印を抽出の要約に写す (run の直下と handoff_extract の delta_r_extract_summary.json; grep '\"settled\": false')",
          m_top is not None and m_ext is not None and m_top["settled"] is False and m_top["unsettled_overrides"] == uo
          and m_top["delta_e_source"] == des and m_ext["settled"] is False
          and all('"settled": false' in f.read_text() for f in (rd / DL.HANDOFF_GATE, rd / "delta_r_extract_summary.json")))
    # --- 正例: 2 窓の実測の超過 (後の窓の中で δ_E が 0.013 mm 振れる) ---
    rd = gate_run("u_2w", pin, delta_fn=jitter)
    g = run(rd, euler_ref=pin, allow_unsettled=US)
    uo = g["unsettled_overrides"]
    check("allow_unsettled: 2 窓の実測の超過を通す → 進む、記録は変動の最大 (0.013 mm)・どの量か (delta_E.range_W2)・窓・x・超えた量の全部",
          g["proceed"] and g["settled"] is False and len(uo) == 1 and uo[0]["condition"] == "two_windows"
          and uo[0]["quantity"] == "delta_E.range_W2" and abs(uo[0]["max_mm"] - 1.7e-4 * mm) < 1e-9
          and uo[0]["windows"] == [[22500, 42500], [42500, 62500]] and uo[0]["tol_mm"] == 0.01
          and set(uo[0]["exceeded"]) == {"delta_E.range_W2"} and uo[0]["x_at_max"] is not None and uo[0]["reason"] == US,
          f"{g['blocked_by']} {uo}")
    # --- 正例: 3 つ全部 (CENSORED + 保存則 + 2 窓) ---
    rd = gate_run("u_all", pin, status="CENSORED", r_fn=lambda s: 1e-3, delta_fn=jitter)
    g = run(rd, euler_ref=pin, allow_not_reached="試験: 未到達", allow_unsettled=US)
    check("allow_not_reached + allow_unsettled: CENSORED・保存則・2 窓をまとめて通す → 進む、unsettled_overrides は停止の状態・保存則・2 窓の順",
          g["proceed"] and g["settled"] is False and [u["condition"] for u in g["unsettled_overrides"]] == ["stop_state", "conservation", "two_windows"]
          and g["unsettled_overrides"][0]["status"] == "CENSORED" and g["override"]["applied"] == ["allow_not_reached", "allow_unsettled"],
          f"{g['blocked_by']} {g['unsettled_overrides']}")
    # --- 上書きを渡しても要らなかった run は settled: true ---
    rd = gate_run("u_notneeded", pin)
    g = run(rd, euler_ref=pin, **ALL)
    check("上書き 3 つを渡しても全部の条件を満たす run → settled: true・unsettled_overrides: []・applied: [] (使わずに通った)",
          g["proceed"] and g["settled"] is True and g["unsettled_overrides"] == [] and g["override"]["applied"] == []
          and marks(rd)[0]["settled"] is True, f"{g['blocked_by']}")
    # --- 既定 (上書きなし) は今までどおり: 保存則の超過で止める ---
    rd = gate_run("u_default", pin, r_fn=lambda s: 1e-3)
    g = run(rd, euler_ref=pin)
    check("既定 (allow_unsettled なし) は今までどおり保存則の超過で止める (settled: false の記録、抽出しない)",
          not g["proceed"] and any("保存則" in b for b in g["blocked_by"]) and g["settled"] is False and g["unsettled_overrides"] == []
          and not (rd / "delta_r_next.csv").exists())

    # --- 負例: 上書き 3 つを全部渡しても通さない ---
    def held_rec(step, key):
        def mod(rd):
            rec = json.loads((rd / SR.RECORD).read_text())
            for r in rec["rows"]:
                if r["step"] == step:
                    r["delta"][key] = 1
            SR.save_record(rd, rec)
        return mod
    neg = [
        ("2 窓の判定不能: 窓の中の出力の δ の記録が無い (抽出の失敗)",
         dict(delta_fn=lambda s, x: None if s == 30000 else (np.full(len(x), 1e-3),) * 2), {}, None, "2 窓"),
        ("2 窓の判定不能: 記録の行が欠けた出力 47,500", dict(steps=[k for k in range(2500, 62501, 2500) if k != 47500]), {}, None,
         "記録の行が無い出力"),
        ("2 窓の判定不能: 履歴が 2 窓を覆わない (停止 30,000)", dict(s_stop=30000), {}, None, "2 窓"),
        ("2 窓: δ の記録の held (出力 50,000)", {}, {}, held_rec(50000, "n_held"), "2 窓"),
        ("2 窓: δ の記録の欠測 (出力 50,000)", {}, {}, held_rec(50000, "n_missing"), "2 窓"),
        ("保存則の判定不能: r の非有限 (出力 50,000)", dict(r_fn=lambda s: float("nan") if s == 50000 else 0.0), {}, None, "保存則"),
        ("保存則の判定不能: 窓の点数の不足 (出力の間隔 5,000)", dict(steps=list(range(5000, 62501, 5000))), {}, None, "保存則"),
        ("停止の状態 DIVERGED", dict(status="DIVERGED"), {}, None, "停止の状態"),
        ("局所 y₁⁺: FAIL", {}, dict(wr="FAIL"), None, "局所 y₁⁺"),
        ("局所 y₁⁺: 判定不能", {}, dict(wr="UNDETERMINED"), None, "局所 y₁⁺"),
        ("被覆: 端の補間範囲の不足 (欠測)", dict(ns=-1.0), {}, None, "被覆"),
        ("次の壁の幾何の健全性", {}, dict(nw=False), None, "次の壁"),
    ]
    for i, (lab, mk, rkw, mod, key) in enumerate(neg):
        rd = gate_run(f"u_neg{i}", pin, **mk)
        if mod is not None:
            mod(rd)
        g = run(rd, euler_ref=pin, **{**ALL, **rkw})
        check(f"負例 (上書き 3 つ) {lab} → 止める ({'; '.join(g['blocked_by'])[:100]}…)",
              not g["proceed"] and any(key in b for b in g["blocked_by"]) and not (rd / "delta_r_next.csv").exists())
    rd = gate_run("u_neg_held", pin)
    DL.extract_and_merge = held_extract
    g = run(rd, euler_ref=pin, **ALL)
    DL.extract_and_merge = spy
    check("負例 (上書き 3 つ) 使用域の held 1 点 → 止める", not g["proceed"] and g["checks"]["coverage"]["n_held_use"] == 1
          and not (rd / "delta_r_next.csv").exists())
    rd = gate_run("u_neg_eref", pin)
    g = run(rd, euler_ref=str(eu), **ALL)
    check("負例 (上書き 3 つ) Euler の参照が本段の記録と違う → 止める", not g["proceed"] and not g["checks"]["euler_ref"]["pass"])
    g = run(rd, euler_ref=pin, omega=0.4, **ALL)
    check("負例 (上書き 3 つ) ω 0.4 → 止める", not g["proceed"] and not g["checks"]["omega"]["pass"])
    rd = gate_run("u_neg_sha", pin, euler_ref={**pin, "sha256": "f" * 64})
    g = run(rd, **ALL)
    check("負例 (上書き 3 つ) Euler の参照の sha の不一致 → 抽出の例外で止める", not g["proceed"] and any("抽出が例外" in b for b in g["blocked_by"]))

    # --- 呼び出し側の δ_E の作り方の口 (extractor + delta_e_source) ---
    win = {"kind": "window_mean", "n_snapshots": 5, "steps": [42500, 47500, 52500, 57500, 62500], "tool": "試験の偽物"}
    rd = gate_run("u_caller", pin, r_fn=lambda s: 1e-3)
    g = run(rd, euler_ref=pin, allow_unsettled=US, extractor=lambda *a, **k: real["extract_and_merge"](*a, **k), delta_e_source=win)
    m_top = marks(rd)[0]
    check("口: extractor + delta_e_source (5 枚の窓の平均の記録) → その抽出で進み、作り方を HANDOFF_GATE.json と抽出の要約に残す (provider caller)",
          g["proceed"] and g["delta_e_source"]["kind"] == "window_mean" and g["delta_e_source"]["n_snapshots"] == 5
          and g["delta_e_source"]["steps"] == win["steps"] and "caller" in g["delta_e_source"]["provider"]
          and m_top["delta_e_source"] == g["delta_e_source"] and m_top["settled"] is False, f"{g['blocked_by']} {g.get('delta_e_source')}")
    for lab, kw in (("extractor だけ", dict(extractor=spy)), ("delta_e_source だけ", dict(delta_e_source=win)),
                    ("kind が空", dict(extractor=spy, delta_e_source={**win, "kind": " "})),
                    ("n_snapshots 0", dict(extractor=spy, delta_e_source={**win, "n_snapshots": 0})),
                    ("n_snapshots が bool", dict(extractor=spy, delta_e_source={**win, "n_snapshots": True})),
                    ("delta_e_source が mapping でない", dict(extractor=spy, delta_e_source="window_mean"))):
        e = raises(lambda kw=kw: DL.handoff_gate(rd, problem="試験.yaml", euler_ref=pin, write=False, **kw), ValueError)
        check(f"口の負例: {lab} → ValueError (作り方の記録と抽出を食い違わせない)", e is not None, str(e)[:100])
    rd = gate_run("u_caller_held", pin)
    g = run(rd, euler_ref=pin, extractor=held_extract, delta_e_source=win, **ALL)
    check("口の負例: 呼び出し側の抽出にも被覆・held の条件が掛かる (使用域の held 1 点 → 上書き 3 つでも止める)",
          not g["proceed"] and g["checks"]["coverage"]["n_held_use"] == 1 and not (rd / "delta_r_next.csv").exists())
    rd = gate_run("u_caller_bad", pin)
    g = run(rd, euler_ref=pin, extractor=lambda *a, **k: None, delta_e_source=win)
    check("口の負例: 呼び出し側の抽出が要約を返さない → 抽出の例外で止める", not g["proceed"] and any("抽出が例外" in b for b in g["blocked_by"]))


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


# ====================================================================================== I. 出口の列の被覆 (幾何の許容)
def euler_ending(name: str, short: float):
    """Euler の参照 run (列の x = linspace(X0 − 0.5, X1, 69) の最後の列だけを X1 − short に)。戻り (run dir, 固定の指定, 列の x)。"""
    x = np.linspace(X0 - 0.5, X1, 69)
    x[-1] = X1 - short
    rw = rwE_f(x)
    R = rw[:, None] * radial(NJ)[None, :]
    X = np.repeat(x[:, None], NJ, axis=1)
    rd = TMP / name
    write_run(rd, X, R, core_q(X, R, rw[:, None]), {}, steps=(1000,))
    return rd, {"run": str(rd), "res": "res_1000.h5", "sha256": DS._sha256_file(rd / "res_1000.h5")}, x


def test_exit_cover(pin0):
    print("== I. 出口の列の被覆: Euler の補間できる範囲の幾何の許容 (§5.1 #36g)", flush=True)
    from forge_design.metrics.delta_in_check import compare_delta_in, wall_columns_from_run
    xE = np.linspace(X0 - 0.5, X1 + 0.5, 70)
    r = DS.euler_cover_range(xE)
    r0 = DS.euler_cover_range(xE, edge_tol=False)
    check("euler_cover_range: 許容は端の外側だけ (tol = 1e−4 × 端の列の間隔)、既定 (edge_tol=False) は ±1e−9",
          r["tol_hi"] == 1e-4 * (xE[-1] - xE[-2]) and r["tol_lo"] == 1e-4 * (xE[1] - xE[0]) and r["hi"] == xE[-1] + r["tol_hi"]
          and r["lo"] == xE[0] - r["tol_lo"] and r0["hi"] == xE[-1] + 1e-9 and r0["lo"] == xE[0] - 1e-9)
    xx = np.array([xE[0] - 2 * r["tol_lo"], xE[0] - 0.5 * r["tol_lo"], xE[0], 3.0, xE[-1], xE[-1] + 0.5 * r["tol_hi"], xE[-1] + 2 * r["tol_hi"]])
    check("euler_covers: 端の外側は許容の内だけを覆う", list(DS.euler_covers(xx, r)) == [False, True, True, True, True, True, False]
          and DS.euler_covers(3.0, r) is True)
    # NS (最後の列 x = X1) と、端が X1 − d の Euler
    Xo, Ro, Qo = synth_ns(+1.0)
    h0 = X1 - np.linspace(X0 - 0.5, X1, 69)[-2]
    d_out = 1.001e-4 * h0 / (1.0 + 1.001e-4)            # 許容 1e−4 × (h0 − d) のちょうど外 (0.1 %)
    cases = {"over": 1.68e-9, "just_out": d_out, "short": 0.07}
    runs = {}
    for k, d in cases.items():
        eu, pin, x = euler_ending(f"eu_i_{k}", d)
        rg = DS.euler_cover_range(wall_columns_from_run(eu)["x"])
        rd = TMP / f"ns_i_{k}"
        write_run(rd, Xo, Ro, Qo, ns_info(pin), steps=(2500,))
        runs[k] = (eu, pin, rd, rg)
    eu, pin, rd, rg = runs["just_out"]
    check(f"前提: 「ちょうど外」の Euler の端は NS の最後の列より許容の 1.001 倍だけ短い (hi {rg['hi']:.15g} < X1 {X1}、差 / 許容 "
          f"{(X1 - rg['x1']) / rg['tol_hi']:.4f})", rg["hi"] < X1 and 1.0005 < (X1 - rg["x1"]) / rg["tol_hi"] < 1.0015)
    kw = dict(band_select="edge", remap_constant_x=True, remap_edge="missing")
    eu, pin, rd, rg = runs["over"]
    a = DS.deltastar_from_core_matched_euler(rd, pin, **kw)
    ne = a["not_extracted"]
    E = DS._load_structured(eu, res="res_1000.h5")
    check(f"1.68e−9 r_t の越え: 全列を抽出 (not_extracted 0)、許容の内側で端を越えた列 {ne['edge_tol_cols']['x']}・越え "
          f"{ne['edge_tol_cols']['beyond'][0]:.3g} r_t を記録、Euler の壁は端の列の値",
          ne["n"] == 0 and len(a["x"]) == NI and ne["edge_tol_cols"]["n"] == 1 and ne["edge_tol_cols"]["x"] == [X1]
          and abs(ne["edge_tol_cols"]["beyond"][0] - 1.68e-9) < 1e-12 and a["r_wall_euler"][-1] == E["r"][-1, -1]
          and a["settings"]["euler_cover"]["edge_tol"], str(ne)[:300])
    b = DS.deltastar_from_core_matched_euler(rd, pin, band_select="edge", remap_constant_x=True)
    check("既定 (remap_edge='clamp') の経路は今と同じ ±1e−9: 1.68e−9 r_t の越えは Euler の範囲の外 (outside_euler)",
          b["not_extracted"]["n"] == 1 and b["not_extracted"]["reason"] == ["outside_euler"] and "euler_cover" not in b["settings"])
    for k, lab in (("just_out", "許容のちょうど外"), ("short", "0.07 r_t の不足")):
        _, pk, rdk, _ = runs[k]
        c_ = DS.deltastar_from_core_matched_euler(rdk, pk, **kw)
        check(f"{lab}: 最後の列は欠測 (outside_euler、外挿しない)", c_["not_extracted"]["n"] == 1 and c_["not_extracted"]["x"] == [X1]
              and c_["not_extracted"]["reason"] == ["outside_euler"] and c_["not_extracted"]["edge_tol_cols"]["n"] == 0)
    # 内側の列には効かない: 許容をほぼ 0 にした抽出と、列ごとの値 (adaptive の帯 = 列ごとに独立) がビット同一
    kw_a = dict(remap_constant_x=True, remap_edge="missing")
    full = DS.deltastar_from_core_matched_euler(rd, pin, **kw_a)
    old_frac = DS.EULER_EDGE_TOL_FRAC
    DS.EULER_EDGE_TOL_FRAC = 1e-12
    try:
        tiny = DS.deltastar_from_core_matched_euler(rd, pin, **kw_a)
    finally:
        DS.EULER_EDGE_TOL_FRAC = old_frac
    same = all(np.array_equal(full[q][:-1], tiny[q], equal_nan=True) for q in ("x", "delta_r_raw", "delta_in", "r_wall_euler", "alpha", "band_y_b"))
    check("内側の列には効かない: 許容あり (最後の列を端の値で抽出) と許容ほぼ 0 (最後の列は欠測) の、内側の全列の列ごとの値がビット同一",
          len(full["x"]) == NI and len(tiny["x"]) == NI - 1 and same)
    # δ_in の照合 (delta_in_check.compare_delta_in): 同じ共通の範囲
    xn = Xo[:, 0]
    rwn = Ro[:, -1]
    for k, want in (("over", True), ("just_out", False), ("short", False)):
        eu, _, _, _ = runs[k]
        W = wall_columns_from_run(eu)
        rr = compare_delta_in(xn, rwn, W["x"], W["rw"], lambda x: dl_f(x), S_SC)
        check(f"δ_in の照合 ({k}): 被覆 {'あり (許容の内側の列 ' + str(rr['n_edge_tol']) + ')' if want else 'なし (最後の列が外)'}",
              rr["covered"] == want and rr["n_uncovered"] == (0 if want else 1) and rr["n_edge_tol"] == (1 if want else 0)
              and rr["euler_cover"]["edge_tol"], str({q: rr[q] for q in ("covered", "n_uncovered", "n_edge_tol", "edge_tol_beyond")}))
    eu, _, _, _ = runs["over"]
    W = wall_columns_from_run(eu)
    r1 = compare_delta_in(xn, rwn, W["x"], W["rw"], lambda x: dl_f(x), S_SC, return_arrays=True)
    check("δ_in の照合: 許容の内側で越えた列の Euler の壁は端の列の値 (δ_in = NS の壁 − Euler の最後の列の壁)",
          r1["arrays"]["r_wall_euler"][-1] == W["rw"][-1] and r1["arrays"]["delta_in"][-1] == rwn[-1] - W["rw"][-1])
    # 関門の抽出と被覆・δ の記録・帯の生成 (受け渡しの条件 handoff_extract_settings)
    eu, pin, rd, _ = runs["over"]
    od = TMP / "i_handoff"
    summ = DL.extract_and_merge(rd, pin, out_dir=od, **DL.handoff_extract_settings())
    cov = DL.coverage_check(summ, od / "delta_r_next.csv", X1)
    check("関門の抽出と被覆 (1.68e−9 r_t の越え): 全域の被覆の合格・許容の内側の列を記録", cov["pass"] and cov["not_extracted"]["n"] == 0
          and cov["not_extracted"]["edge_tol_cols"]["n"] == 1, str(cov.get("why")))
    _, pk, rdk, _ = runs["short"]
    summ = DL.extract_and_merge(rdk, pk, out_dir=TMP / "i_handoff_short", **DL.handoff_extract_settings())
    cov = DL.coverage_check(summ, TMP / "i_handoff_short" / "delta_r_next.csv", X1)
    check("関門の被覆 (0.07 r_t の不足): 不合格 (outside_euler)", not cov["pass"] and cov["not_extracted"]["reason"] == ["outside_euler"])
    uc = SR.use_columns(rd)
    one = SR.delta_snapshot(rd, 2500, pin, {**DL.handoff_extract_settings(), "source": "試験"}, uc["x"])
    band = SR.fix_band(rd, 2500, pin, TMP / "i_band", **SR.band_extract_kw(ns_info(pin)["stop_rule"]))
    check("δ の記録 (delta_snapshot) は使用域の最後の列 (x_w) まで欠測 0、帯の生成 (fix_band) も全列",
          uc["x"][-1] == X1 and one["n_missing"] == 0 and np.all(np.isfinite(one["delta_E"])) and len(band["x"]) == NI and band["x"][-1] == X1)
    # 2 窓に入るべき出力の列挙
    check("required_window_steps: 停止 62,500・窓 20,000・間隔 2,500 → 17 出力 (22,500〜62,500)、不正な引数は例外",
          SR.required_window_steps(62500, 20000, 2500) == list(range(22500, 62501, 2500))
          and raises(lambda: SR.required_window_steps(62500, 20000, None), ValueError) is not None
          and raises(lambda: SR.required_window_steps(62500, 0, 2500), ValueError) is not None)
    _ = pin0


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
    test_exit_cover(pin)
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
