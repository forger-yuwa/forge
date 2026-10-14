#!/usr/bin/env python3
r"""③ベル (M) モード: CFD ピンの初期線 (W3) — 凍結源の照合・ピン線の抽出・ピン版の (M) 壁と Hall 版との差 (記録、合否は親と上位)。

計画: plans/active/tooling-nozzle-bell-two-mode.md §5.1 #6 (W3 の前段・記録項目)・§6 V2 (ピン版の事前登録)。判断の記録:
notes/reviews/2026-10-12-bell-m-v2-euler-diagnostician.md (凍結源 = run_0015、前段 (a) 照合・(b) ṁ のドリフト・(c) 格子 ×1.5)。

段 (サブコマンド):

- `lines` (数十秒): 凍結源の候補 run (`SOURCE_RUNS`; 無い run は飛ばす) ごとに
  (a) `cfd_initial_line.check_source_matches` の照合 (node・軸対称・Euler [viscMethod 0・visc 0・乱流なし・物理壁 slip]・
  ガス [排気ファイルの記録・熱力学の参照値・入口 Pt・Tt・組成・physProp]・寸法 [R_t・L_U・L_pipe・r_U]・収縮部の壁節点と
  設計の `UpstreamThroatPoly`)。照合できない項目は理由つきで `unchecked_reasons` に、
  (b) ṁ の末尾 (本段 15000〜30000) の相対ドリフト (thrust_history.csv)、
  ピン線 (`moc_initial_line`、n 321・r 等間隔) の検査 (`rao_axisym.check_initial_line` = `moc_forward._check_line`)・
  線を横切る流量 m*/ṁ_1D (等エントロピー写像、`rao_axisym.line_integrals`)・C⁻ の適合式の不整合・Hall n 321 との差、
  鋸歯: (i) 線そのものの 2Δ 指標、(ii) Hall 版 L100 の壁 (`rao_t6_parts/L100_n321_mq2`) の前進 MOC の前線 k₀ (壁点 x ≈ 0.3) の
  帯 k₀ − j ∈ [−6, 12] (`run_rao_t4h.py` の `reproduction` と同じ方法。Hall の線で記録 1.45e-3 rad を再現する) と出だしの軸点の θ_r。
  run の組の差 (0015 vs 0016 = 照合、0015 vs 0017 = 格子確認 (c): 登録の線 m*/ṁ_1D の差 ≤ 1e-4・θ の最大差 ≤ 1e-4 rad)。
  出力 `rao_pin_line.json`。
- `design` (2 本並列、各 25〜30 分): (1) ピン版 L100 (凍結源 `FREEZE` のピン線 n 321、円弧 kernel、Goursat m = q = n = 321 [既定]) と
  同じ根で m = q = 320 に作り直した壁 (Hall 版 `L100_n321_mq2` と同じ m・q での比較用)、(2) Hall の回し直し (既定の初期線
  `line=None`・n 321、m = q = 160 で `rao_t6_parts/L100_n321` と、作り直しの 320 で `L100_n321_mq2` と値・壁の点列を比べる —
  記録の Hall 版を今のコードが再現するか。321 の壁は `rao_pin_hallref/` に置く)。出力 `rao_pin_parts/<tag>.json`・
  `<tag>_wall.npz` (`run_rao_t6.py` と同じ形式 — runner の `geometry.rao_wall_file` で読める) と `rao_pin.json` (キャンペーンの要約:
  `setup` は `load_rao_wall` が WALL_DESIGN.json へ写す)。
- `compare` (数十秒): Hall 版とピン版の (M) の差の表 (θ_B・θ_E・θ_max・x_E・x_D・F_cs・2 会計の差・流量・壁の r の差 [全域と
  x 3.4〜5.3]・BE の 2Δ 指標)、目安 Δθ_B・Δθ_E ≤ 0.1°・壁 Δr ≤ 1e-3 r_t (超えたら記録、止めない)、ピン版の壁ファイルが
  `load_rao_wall` で読めること、V2 の問題 YAML (`../problem_bell_tp_m_euler_l100_pin.yaml`) で runner が壁を作り
  WALL_DESIGN.json に moc_reference (F_cs・η_cs・ṁ/ṁ_1D) を書くこと (変換器・forge は回さない)。`rao_pin.json` に追記し
  `rao_pin_compare.csv` を書く。

実行 (手元 PC はユーザの作業機 — 並列 2 本まで、nice):
`nice -n 10 /home/sano/work/forge/design/.venv-opt/bin/python case/68.bell_two_mode/rao_m/run_pin_line.py lines|design|compare`
"""
from __future__ import annotations

import os

for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_k, "1")

import argparse  # noqa: E402
import csv  # noqa: E402
import hashlib  # noqa: E402
import json  # noqa: E402
import subprocess  # noqa: E402
import sys  # noqa: E402
import tempfile  # noqa: E402
import time  # noqa: E402
import traceback  # noqa: E402
import warnings  # noqa: E402
from concurrent.futures import ProcessPoolExecutor, as_completed  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
from scipy.interpolate import CubicSpline  # noqa: E402

HERE = Path(__file__).resolve().parent
CASE = HERE.parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "design"))
sys.path.insert(0, str(HERE))
import run_rao_t6 as t6  # noqa: E402  (記録の形式 `_record`・m・q の作り直し `_rebuild_mq`・`_jsonable` を共用)
from forge_design.feedback.cfd_initial_line import (expected_bell_source, moc_initial_line,  # noqa: E402
                                                    pinned_factory)
from forge_design.gas.rocket_exhaust import RocketExhaust  # noqa: E402
from forge_design.geometry import moc_forward as mf  # noqa: E402
from forge_design.geometry import rao_axisym as ra  # noqa: E402
from forge_design.geometry.moc_kernel import pm_nu  # noqa: E402
from forge_design.geometry.wall_bell_c2 import cone15_length  # noqa: E402
from forge_design.probdef import load_problem  # noqa: E402

EXHAUST = CASE / "exhaust_lox_lh2_pc100_of6.yaml"
PROBLEM = CASE / "problem_bell_tp_m_euler_l100.yaml"           # 使う側 (寸法・ガス・入口) = 照合の期待値
PROBLEM_PIN = CASE / "problem_bell_tp_m_euler_l100_pin.yaml"    # ピン版 V2 の問題 (壁ファイルだけ差し替え)
SOURCE_RUNS = ("run_0015_tp_m_euler_l100", "run_0016_tp_m_euler_l70", "run_0017_tp_m_euler_l100_fine")
FREEZE = "run_0015_tp_m_euler_l100"                             # 凍結源 (diagnostician 2026-10-12)
RES = "res_30000.h5"
PAIRS = (("run_0015_tp_m_euler_l100", "run_0016_tp_m_euler_l70", "照合 (依存域は x_B より上流 — 一致するはず)"),
         ("run_0015_tp_m_euler_l100", "run_0017_tp_m_euler_l100_fine", "格子確認 (c): ni 1801・nj 181 (×1.5)"))
GRID_LINE = {"mstar_ratio": 1e-4, "theta_rad": 1e-4}           # (c) の登録の線 (§5.1 #6)
MSTAR_EXPECT = (0.9926, 2e-4)                                   # 期待 m*/ṁ_1D (§5.1 #6 W3 の記録項目)
MDOT_DRIFT_LINE, MDOT_STEP0 = 1e-5, 15000                        # (b) の線と区間の始まり
N_LINE = 321
R_T, EPS, P_A = 1.0, 50.0, 0.0
L100 = cone15_length(EPS)
X_FRONT = 0.3                                                    # 鋸歯を測る前線の壁点 (`run_rao_t4h.py` と同じ)
BAND = (-6, 12)                                                  # 前線の鋸歯の帯 k₀ − j (同)
HALL_WALL = HERE / "rao_t6_parts" / "L100_n321_mq2_wall.npz"     # Hall 版 (run_0015 が回した壁)
HALL_REC = {"mq160": HERE / "rao_t6_parts" / "L100_n321", "mq320": HERE / "rao_t6_parts" / "L100_n321_mq2"}
PARTS = HERE / "rao_pin_parts"
HALLREF = HERE / "rao_pin_hallref"
OUT_LINE = HERE / "rao_pin_line.json"
OUT = HERE / "rao_pin.json"
OUT_CSV = HERE / "rao_pin_compare.csv"
TAG_PIN = "L100_pin_n321"                                       # 生産の既定 (m = q = n = 321) — V2 の壁
TAG_PIN_320 = "L100_pin_n321_mq320"                              # 同じ根で m = q = 320 (Hall 版 mq2 と同じ m・q)
TAG_HALL_RERUN = "L100_hall_n321_rerun"
TAG_HALL_321 = "L100_hall_n321_mq321"
X_SAW = (3.4, 5.3)                                               # Hall 版の壁で鋸歯が乗る区間 (plan §5.1 #5 の欠陥票)
TARGET = {"dtheta_deg": 0.1, "dr": 1e-3}                         # Hall 版とピン版の差の目安 (diagnostician 2026-10-12)
SOURCES = ("design/forge_design/geometry/rao_axisym.py", "design/forge_design/feedback/cfd_initial_line.py",
           "design/forge_design/gas/isentrope.py", "design/forge_design/geometry/moc_forward.py",
           "design/forge_design/geometry/moc_kernel.py", "design/forge_design/evaluate/runner.py")
CMP_KEYS = ("theta_B_deg", "theta_E_deg", "theta_max_deg", "x_E", "x_E_minus_L", "x_B", "r_B", "s_B", "x_D", "r_D", "M_D",
            "theta_D_deg", "eta_D", "det_D", "de", "M_E", "F_cs", "F_wall", "thrust_rel_diff", "mdot_rel_diff", "T_line",
            "T_kernel", "T_BE", "T_closure", "T_DE", "CF", "isp_ratio", "closure_norm", "kappaW_B", "dkappa_B", "g_tilde",
            "goursat_nfail", "n_wall_missing", "streamline_resid_max")


def git_info() -> dict:
    head = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "-C", str(ROOT), "status", "--porcelain", "design/forge_design", "design/tests"],
                           capture_output=True, text=True).stdout
    sha = {s: hashlib.sha256((ROOT / s).read_bytes()).hexdigest()[:16] for s in SOURCES}
    return {"commit": head, "design_tree_dirty": [ln for ln in dirty.splitlines() if ln.strip()], "source_sha256_16": sha}


def gas():
    """生成器の気体 (`run_rao_t6.py` と同じ: 排気ファイルの `RocketExhaust.moc_gas`)。"""
    return ra.as_isentrope(RocketExhaust.from_file(EXHAUST).moc_gas())


def _dump(path: Path, obj) -> None:
    path.write_text(json.dumps(t6._jsonable(obj), ensure_ascii=False, indent=1))


# ================================================================================================ 線
def pinned_throat(run: str, ig):
    """凍結源 run の CFD ピン (照合つき)。戻り (throat, 照合の記録)。"""
    eb = expected_bell_source(load_problem(PROBLEM))
    g_hall = float(ig.moc.gamma_throat())
    th = pinned_factory(CASE / run, RES, eb["expect"])(eb["R"], g_hall)
    match = dict(th.source.get("match") or {})
    reasons = {}
    for k in match.get("unchecked", []):
        reasons[k] = {"gamma_hall": "③の凍結源の prepare_info は Hall の γ を持たない (evaluate/runner.py)。γ_Hall を使うのは "
                                    "M″_A (Hall の式) だけで、③の初期線 (`moc_initial_line`) は M_A・M′_A だけを使う — 照合の対象外"}.get(
            k, "照合の記録が凍結源に無い")
    match["unchecked_reasons"] = reasons
    match["R_used"], match["gamma_hall_used"] = eb["R"], g_hall
    return th, match


def mdot_drift(run: str) -> dict:
    """(b) thrust_history.csv の本段の末尾 (step ≥ 15000) の ṁ の相対ドリフト。"""
    rows = list(csv.DictReader(open(CASE / run / "thrust_history.csv")))
    rows = [r for r in rows if int(r["step"]) >= MDOT_STEP0]
    out = {"steps": [int(rows[0]["step"]), int(rows[-1]["step"])], "n": len(rows), "line": MDOT_DRIFT_LINE}
    for k in ("mdot_exit_kg_s", "mdot_in_kg_s", "mdot_ratio_1d"):
        v = np.array([float(r[k]) for r in rows])
        out[k] = {"first": float(v[0]), "last": float(v[-1]), "mean": float(v.mean()),
                  "rel_last_minus_first": float((v[-1] - v[0]) / v.mean()), "rel_max_minus_min": float((v.max() - v.min()) / v.mean())}
    out["within_line"] = bool(all(out[k]["rel_max_minus_min"] <= MDOT_DRIFT_LINE for k in ("mdot_exit_kg_s", "mdot_in_kg_s")))
    out["final_mdot_ratio_1d"] = float(rows[-1]["mdot_ratio_1d"])
    return out


def compat_mismatch(x, r, th, M, thr_axis, g) -> float:
    """C⁻ の適合式 d(θ+ν) = sinμ (sinθ/r) dσ の線に沿った不整合 (θ+ν−∫源項 の壁と軸端の差) [rad]
    (`design/tests/run_moc_forward_tests.py` の info と同じ式。Hall R 1 で −1.09e-2)。"""
    nu = np.asarray(pm_nu(M, g), dtype=float)
    mu = np.arcsin(1.0 / M)
    s_ = np.where(r > 1e-9, np.sin(th) / np.maximum(r, 1e-300), thr_axis)
    f_ = np.sin(mu) * s_ / np.cos(th - mu)
    dS = 0.5 * (f_[1:] + f_[:-1]) * (x[:-1] - x[1:])
    J = th + nu - np.concatenate([np.cumsum(dS[::-1])[::-1], [0.0]])
    return float(J[-1] - J[0])


def saw2(f, idx=None) -> dict:
    """2Δ の鋸歯の指標: a_j = (f_{j−1} − 2f_j + f_{j+1})/4 (純粋な 2Δ の振幅 A なら |a| = A、滑らかなら h² f''/4)。
    最大・中央値・符号の交代の割合 (2Δ なら 1 に近い)。`idx` は a の位置 (j) の部分集合のマスク。"""
    f = np.asarray(f, dtype=float)
    a = (f[:-2] - 2.0 * f[1:-1] + f[2:]) / 4.0
    j = np.arange(1, len(f) - 1)
    if idx is not None:
        a, j = a[idx[1:-1]], j[idx[1:-1]]
    if len(a) < 2:
        return {"n": int(len(a))}
    alt = np.sign(a[1:]) != np.sign(a[:-1])
    k = int(np.argmax(np.abs(a)))
    return {"n": int(len(a)), "max_abs": float(np.abs(a).max()), "j_at_max": int(j[k]), "median_abs": float(np.median(np.abs(a))),
            "alt_frac": float(np.mean(alt))}


def front_sawtooth(line: dict, ig) -> dict:
    """Hall 版 L100 の壁 (B より上流は円弧 kernel) の前進 MOC を前線 120 本まで回し、壁点 x ≈ 0.3 の前線 k₀ の鋸歯の帯
    (`run_rao_t4h.py` の `reproduction` と同じ: 帯 k₀ − j ∈ [−6, 12] の外の点を通る 3 次スプラインからのずれ) と、前線全体の
    2Δ 指標、出だしの前線の軸点の θ_r (Hall では符号が振動する — 鋸歯の起点)。前線 k₀ は B (x 0.496) より上流なので壁の
    B より先には依らない (打ち切りの前線数にも依らない — 前線は順に作る)。"""
    w, _ = ra.load_rao_wall(HALL_WALL)
    L = float(w.x_end)
    res = mf.forward_moc(ra.ExtendedWall(w, L), L, ig.moc, line, max_fronts=120)
    net, n = res["net"], res["n"]
    xw = net.a["x"][:, n]
    k0 = int(np.nanargmin(np.abs(xw - X_FRONT)))
    fr = {f: net.a[f][k0] for f in net.FIELDS}
    h = hashlib.sha256(b"".join(np.ascontiguousarray(net.a[f][k0]).tobytes() for f in net.FIELDS)
                       + np.float64(net.thr_axis[k0]).tobytes()).hexdigest()[:16]
    band = (k0 - np.arange(n + 1) <= BAND[1]) & (k0 - np.arange(n + 1) >= BAND[0])
    dev = {}
    for nm in ("th", "M"):
        d = fr[nm] - CubicSpline(fr["r"][~band], fr[nm][~band])(fr["r"])
        jj = np.flatnonzero(band)[int(np.argmax(np.abs(d[band])))]
        dd = np.diff(d[band])
        # 2Δ の鋸歯は「片側のずれ + 点ごとに大小が交代する成分」(Hall: −3.1e-4, −6.9e-4, −5.6e-4, −1.0e-3, …) なので、ずれの符号でなく
        # ずれの増分の符号の交代率で見る (2Δ なら → 1、滑らかなら → 0)。前線の点の r の間隔は C⁺ の境 (k₀ − j = 0、初期線の軸端からの
        # C⁺) で跳ぶので、添字の 2 階差分 (`saw2`) は前線には使わない
        dev[nm] = {"max_abs_dev_in_band": float(np.abs(d[band]).max()), "k0_minus_j_at_max": int(k0 - jj),
                   "max_abs_dev_outside_band": float(np.abs(d[~band]).max()),
                   "ddev_sign_alt_frac_in_band": float(np.mean(np.sign(dd[1:]) != np.sign(dd[:-1])))}
    # 網の点の r の間隔の 2Δ (Hall では帯で 1 点おきに 0.0084 ↔ 0.0090 と交代した): 間隔の増分の符号の交代率
    ddr = np.diff(np.diff(fr["r"]))
    bm = band[1:-1]
    return {"k0": k0, "x_wall": float(xw[k0]), "all_ok_to_axis": bool(np.all(net.status[k0] == mf.ST_OK)),
            "front_sha256_16": h, "march_end": res["end"], "band_k0_minus_j": list(BAND),
            "band_r_range": [float(fr["r"][band].min()), float(fr["r"][band].max())], "dev_from_smooth": dev,
            "dr_increment_sign_alt_frac_in_band": float(np.mean(np.sign(ddr[bm][1:]) != np.sign(ddr[bm][:-1]))),
            "dr_increment_sign_alt_frac_outside_band": float(np.mean(np.sign(ddr[~bm][1:]) != np.sign(ddr[~bm][:-1]))),
            "thr_axis_first_fronts": [float(v) for v in net.thr_axis[1:13]]}


def line_record(line: dict, ig, hall: dict | None) -> dict:
    x, r, th, M = (np.asarray(line[k], dtype=float) for k in ("x", "r", "theta", "M"))
    m_s, T_s = ra.line_integrals(ig, x, r, th, M, P_A, spline=True)
    m_t = mf._line_flux(x, r, th, M, ig.moc)
    m1d = float(ig.ideal_1d(EPS, P_A)["mdot"]) * np.pi                 # ṁ_1D (ρ*a*·πr_t²、2π 込みの単位)
    out = {"source": line.get("source"), "n": int(len(x)), "x_axis": float(x[0]), "M_axis": float(M[0]),
           "M_wall": float(M[-1]), "M_min": float(M.min()), "thr_axis": float(line["thr_axis"]),
           "Mp_axis_anchor": line.get("Mp_axis_anchor"), "check_initial_line": ra.check_initial_line(line, R_T, ig),
           "mstar_nd": m_s, "mstar_ratio_1d": m_s / m1d, "mstar_ratio_1d_trapezoid": m_t / m1d, "T_line_nd": T_s,
           "saw2_theta": saw2(th), "saw2_M": saw2(M),
           "sha256_16": hashlib.sha256(np.concatenate([x, r, M, th]).tobytes()).hexdigest()[:16]}
    if hall is not None:
        out["minus_hall"] = line_diff(line, hall, ig)
    return out


def line_diff(a: dict, b: dict, ig) -> dict:
    """同じ r の点 (どちらも r 等間隔 n 点) で比べた差 a − b (θ・M・x の最大、m*/ṁ_1D、軸端)。"""
    if not np.array_equal(np.asarray(a["r"]), np.asarray(b["r"])):
        raise ValueError("line_diff: r の点が違う")
    m1d = float(ig.ideal_1d(EPS, P_A)["mdot"]) * np.pi
    ma = ra.line_integrals(ig, a["x"], a["r"], a["theta"], a["M"], P_A, spline=True)[0]
    mb = ra.line_integrals(ig, b["x"], b["r"], b["theta"], b["M"], P_A, spline=True)[0]
    out = {}
    for k in ("theta", "M", "x"):
        d = np.asarray(a[k], dtype=float) - np.asarray(b[k], dtype=float)
        j = int(np.argmax(np.abs(d)))
        out[f"max_abs_d{k}"] = float(np.abs(d[j]))
        out[f"d{k}_at_max"] = float(d[j])
        out[f"r_at_max_d{k}"] = float(a["r"][j])
    out.update({"dmstar_ratio_1d": (ma - mb) / m1d, "mstar_rel": ma / mb - 1.0, "dx_axis": float(a["x"][0] - b["x"][0]),
                "dthr_axis": float(a["thr_axis"] - b["thr_axis"]), "dM_axis": float(a["M"][0] - b["M"][0])})
    return out


def cmd_lines() -> dict:
    t0 = time.time()
    ig = gas()
    hall = ra._hall_line(R_T, ig, N_LINE)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        hall_2001 = mf._hall_line_raw(R_T, ig.moc, 2001)
    rec = {"plan": "plans/active/tooling-nozzle-bell-two-mode.md §5.1 #6 (W3)・§6 V2",
           "decision_record": "notes/reviews/2026-10-12-bell-m-v2-euler-diagnostician.md", **git_info(),
           "setup": {"problem_expect": str(PROBLEM.relative_to(ROOT)), "res": RES, "n_line": N_LINE, "freeze": FREEZE,
                     "mstar_expect": list(MSTAR_EXPECT), "grid_line": GRID_LINE, "mdot_drift_line": MDOT_DRIFT_LINE,
                     "mstar_definition": "線を横切る等エントロピー写像の流量 (rao_axisym.line_integrals、3 次スプライン + GL) / "
                                         "ṁ_1D (ρ*a*πr_t²)。trapezoid は moc_forward._line_flux (中点則)"},
           "hall": {**line_record(hall, ig, None),
                    "compat_mismatch_2001": compat_mismatch(hall_2001["x"], hall_2001["r"], hall_2001["theta"], hall_2001["M"],
                                                            hall_2001["thr_axis"], ig.moc),
                    "front": front_sawtooth(hall, ig)},
           "runs": {}, "pairs": {}}
    lines = {}
    for run in SOURCE_RUNS:
        if not (CASE / run / RES).exists():
            rec["runs"][run] = {"skipped": f"{run}/{RES} が無い"}
            continue
        try:
            th, match = pinned_throat(run, ig)
            line = moc_initial_line(th, ig.moc, N_LINE)
            l2001 = moc_initial_line(th, ig.moc, 2001)
            lines[run] = line
            lr = line_record(line, ig, hall)
            lr["compat_mismatch_2001"] = compat_mismatch(l2001["x"], l2001["r"], l2001["theta"], l2001["M"], l2001["thr_axis"], ig.moc)
            lr["front"] = front_sawtooth(line, ig)
            dr_ = mdot_drift(run)
            m_eu = dr_["final_mdot_ratio_1d"]
            lr["mstar_vs_euler_rel"] = lr["mstar_ratio_1d"] / m_eu - 1.0
            lr["mstar_minus_expect"] = lr["mstar_ratio_1d"] - MSTAR_EXPECT[0]
            lr["mstar_within_expect"] = bool(abs(lr["mstar_ratio_1d"] - MSTAR_EXPECT[0]) <= MSTAR_EXPECT[1])
            rec["runs"][run] = {"match": match, "cfd": line["cfd"], "mdot_drift": dr_, "line": lr,
                                "anchor": {"M_A": line["M_axis_anchor"], "Mp_A": line["Mp_axis_anchor"], "x_A": line["x_axis"]}}
            np.savez(HERE / f"rao_pin_line_{run.split('_')[1]}_n{N_LINE}.npz", x=line["x"], r=line["r"], M=line["M"],
                     theta=line["theta"], thr_axis=line["thr_axis"])
            print(f"{run}: 照合 {len(match['checked'])} 項目 (照合できない {match['unchecked']})、m*/ṁ_1D {lr['mstar_ratio_1d']:.7f} "
                  f"(Euler {m_eu:.7f}、相対 {lr['mstar_vs_euler_rel']:+.2e})、Hall との差 {lr['minus_hall']['dmstar_ratio_1d']:+.3e}・"
                  f"θ 最大 {lr['minus_hall']['max_abs_dtheta']:.2e} rad・M 最大 {lr['minus_hall']['max_abs_dM']:.2e}、"
                  f"ṁ ドリフト {dr_['mdot_exit_kg_s']['rel_max_minus_min']:.1e}、前線の帯 θ {lr['front']['dev_from_smooth']['th']['max_abs_dev_in_band']:.2e}",
                  flush=True)
        except Exception as e:  # noqa: BLE001
            rec["runs"][run] = {"error": f"{type(e).__name__}: {e}", "traceback": traceback.format_exc()}
            print(f"{run}: ERROR {type(e).__name__}: {e}", flush=True)
    for a, b, what in PAIRS:
        if a in lines and b in lines:
            d = line_diff(lines[b], lines[a], ig)
            d["what"] = what
            if "run_0017" in b:
                d["registered_line"] = GRID_LINE
                d["within_registered_line"] = bool(abs(d["dmstar_ratio_1d"]) <= GRID_LINE["mstar_ratio"]
                                                   and d["max_abs_dtheta"] <= GRID_LINE["theta_rad"])
            rec["pairs"][f"{b}_minus_{a}"] = d
            print(f"{b} − {a} ({what}): m*/ṁ_1D {d['dmstar_ratio_1d']:+.2e}・θ 最大 {d['max_abs_dtheta']:.2e} rad "
                  f"(r {d['r_at_max_dtheta']:.3f})・M 最大 {d['max_abs_dM']:.2e}・x 最大 {d['max_abs_dx']:.2e}", flush=True)
    rec["time_s"] = time.time() - t0
    _dump(OUT_LINE, rec)
    print(f"{OUT_LINE.name} ({rec['time_s']:.0f} s)", flush=True)
    return rec


# ================================================================================================ 設計
def _setup(line: dict | None = None) -> dict:
    s = {"gas": "EXH (RocketExhaust.moc_gas)", "exhaust": str(EXHAUST.relative_to(ROOT)), "R_t": R_T, "eps": EPS, "p_a": P_A,
         "kernel": "arc", "accept_rule": "det_D > 0", "L_rt": L100, "L_frac": 1.0, "n_line": N_LINE, "dr_default": "3.2/n",
         "mq_default": "n (2026-10-12)",
         "initial_line": f"CFD ピン ({FREEZE}/{RES}、moc_initial_line n {N_LINE}・r 等間隔)",
         "note": "合否は判定しない (親と上位)。目安 Δθ_B・Δθ_E ≤ 0.1°・壁 Δr ≤ 1e-3 r_t (diagnostician 2026-10-12)"}
    if line is not None:
        s["initial_line_detail"] = {"run": FREEZE, "res": RES, "cfd": line.get("cfd"), "x_axis": line["x_axis"],
                                    "M_axis_anchor": line["M_axis_anchor"], "Mp_axis_anchor": line["Mp_axis_anchor"],
                                    "thr_axis": line["thr_axis"],
                                    "sha256_16": hashlib.sha256(np.concatenate([line[k] for k in ("x", "r", "M", "theta")]).tobytes()).hexdigest()[:16]}
    return s


def _save(dirp: Path, rec: dict, res: dict | None) -> None:
    dirp.mkdir(exist_ok=True)
    _dump(dirp / f"{rec['tag']}.json", rec)
    if res is not None and res.get("status") == "ok":
        np.savez(dirp / f"{rec['tag']}_wall.npz", BE=res["BE"], kernel_wall=res["kernel_wall"], wall=res["wall"])


def job_pin() -> list:
    t0 = time.time()
    ig = gas()
    th, match = pinned_throat(FREEZE, ig)
    line = moc_initial_line(th, ig.moc, N_LINE)
    spec = {"tag": TAG_PIN, "L_frac": 1.0, "variant": "pin"}
    res = ra.design_rao(ig, R_t=R_T, L=L100, eps=EPS, kernel="arc", line=line, p_a=P_A)
    rec = t6._record(res, spec, L100, time.time() - t0)
    lr = {"source": "cfd_pinned", "run": FREEZE, "res": RES, "n": N_LINE, "match": match, **_setup(line)["initial_line_detail"]}
    rec.update({"line_source": res.get("line_source"), "initial_line": lr,
                "note": f"CFD ピンの初期線 ({FREEZE}/{RES}、n {N_LINE}・r 等間隔) の (M)。Goursat m = q = n (既定)"})
    _save(PARTS, rec, res)
    out = [rec]
    print(t6._line(rec), flush=True)
    if res["status"] == "ok":
        t1 = time.time()
        r2 = t6._rebuild_mq(res, ig, 320, 320)
        rec2 = t6._record(r2, {**spec, "tag": TAG_PIN_320, "variant": "pin_mq320"}, L100, time.time() - t1)
        rec2.update({"line_source": res.get("line_source"), "initial_line": lr,
                     "note": f"{TAG_PIN} と同じ根 (kernel 場・BD・DE) で Goursat の m・q を 320 (Hall 版 L100_n321_mq2 と同じ) に。"
                             f"time_s は作り直しの分だけ (設計 1 本は {rec['time_s']:.0f} s)"})
        _save(PARTS, rec2, r2)
        out.append(rec2)
        print(t6._line(rec2), flush=True)
    return out


def _same_rec(new: dict, old: dict) -> dict:
    keys = [k for k in CMP_KEYS if k in old and k in new]
    diff = {k: (new[k] - old[k]) if isinstance(new[k], (int, float)) and isinstance(old[k], (int, float)) else
            (None if new[k] == old[k] else "differs") for k in keys}
    nz = {k: v for k, v in diff.items() if v not in (0, 0.0, None)}
    return {"n_keys": len(keys), "identical": not nz, "nonzero": nz}


def _same_npz(arrs: dict, path: Path) -> dict:
    with np.load(path) as z:
        out = {k: bool(np.array_equal(arrs[k], z[k])) for k in ("BE", "kernel_wall", "wall")}
        out["BE_max_abs_diff"] = float(np.max(np.abs(arrs["BE"] - z["BE"]))) if arrs["BE"].shape == z["BE"].shape else None
    return out


def job_hall_rerun() -> list:
    """既定の Hall の初期線 (`line=None`) で L100 n 321 を回し直し、記録 (`rao_t6_parts/L100_n321` [m = q = 160] と
    `L100_n321_mq2` [320]) と値・壁の点列を比べる。321 (今の既定) の壁は `rao_pin_hallref/` に置く。"""
    t0 = time.time()
    ig = gas()
    spec = {"tag": TAG_HALL_RERUN, "L_frac": 1.0, "variant": "hall_rerun_mq160"}
    res = ra.design_rao(ig, R_t=R_T, L=L100, eps=EPS, kernel="arc", n=N_LINE, m=160, q=160, p_a=P_A)
    rec = t6._record(res, spec, L100, time.time() - t0)
    rec["line_source"] = res.get("line_source")
    out = [rec]
    if res["status"] == "ok":
        old = json.loads((HALL_REC["mq160"].with_suffix(".json")).read_text())
        rec["vs_record_mq160"] = {**_same_rec(rec, old), "npz": _same_npz(res, Path(str(HALL_REC["mq160"]) + "_wall.npz"))}
        r2 = t6._rebuild_mq(res, ig, 320, 320)
        rec2 = t6._record(r2, {**spec, "tag": TAG_HALL_RERUN + "_mq320", "variant": "hall_rerun_mq320"}, L100, 0.0)
        old2 = json.loads((HALL_REC["mq320"].with_suffix(".json")).read_text())
        rec["vs_record_mq320"] = {**_same_rec(rec2, old2), "npz": _same_npz(r2, Path(str(HALL_REC["mq320"]) + "_wall.npz"))}
        r3 = t6._rebuild_mq(res, ig, 321, 321)
        rec3 = t6._record(r3, {**spec, "tag": TAG_HALL_321, "variant": "hall_mq321"}, L100, 0.0)
        rec3.update({"line_source": "hall", "note": "Hall n 321 の (M) を今の既定 m = q = n = 321 で (ピン版 L100_pin_n321 と同じ m・q の比較用)"})
        _save(HALLREF, rec3, r3)
        out.append(rec3)
    _save(HALLREF, rec, None)
    print(t6._line(rec) + f" | 記録との一致 mq160 {rec.get('vs_record_mq160', {}).get('identical')}・"
          f"mq320 {rec.get('vs_record_mq320', {}).get('identical')}", flush=True)
    return out


def write_campaign(extra: dict | None = None) -> dict:
    prev = json.loads(OUT.read_text()) if OUT.exists() else {}
    line = None
    if (PARTS / f"{TAG_PIN}.json").exists():
        line = json.loads((PARTS / f"{TAG_PIN}.json").read_text()).get("initial_line")
    s = _setup()
    if line is not None:
        s["initial_line_detail"] = {k: line.get(k) for k in ("run", "res", "cfd", "x_axis", "M_axis_anchor", "Mp_axis_anchor",
                                                              "thr_axis", "sha256_16")}
    rows = [json.loads(p.read_text()) for p in sorted(PARTS.glob("*.json"))] if PARTS.exists() else []
    out = {**prev, "plan": "plans/active/tooling-nozzle-bell-two-mode.md §5.1 #6 (W3)・§6 V2",
           "decision_record": "notes/reviews/2026-10-12-bell-m-v2-euler-diagnostician.md", "setup": s, **git_info(),
           "rows": rows, **(extra or {})}
    _dump(OUT, out)
    return out


def cmd_design(workers: int) -> None:
    if workers > 2:
        raise SystemExit("--workers は 2 まで (手元 PC はユーザの作業機)")
    t0 = time.time()
    jobs = {"pin": job_pin, "hall_rerun": job_hall_rerun}
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(f): k for k, f in jobs.items()}
        for fu in as_completed(futs):
            try:
                fu.result()
            except Exception as e:  # noqa: BLE001
                print(f"{futs[fu]}: ERROR {type(e).__name__}: {e}\n{traceback.format_exc()}", flush=True)
            write_campaign({"design_total_s": time.time() - t0})
    print(f"design 完了 ({time.time() - t0:.0f} s)", flush=True)


# ================================================================================================ 比較
def _wall_r(path: Path):
    w, info = ra.load_rao_wall(path)
    return w, info


def be_saw(BE: np.ndarray) -> dict:
    x, th = BE[:, 0], BE[:, 2]
    out = {}
    for lo, hi, nm in ((BE[0, 0], X_SAW[0], "B_to_3.4"), (X_SAW[0], X_SAW[1], "3.4_to_5.3"), (X_SAW[1], 10.0, "5.3_to_10"),
                       (10.0, BE[-1, 0] + 1.0, "10_to_E")):
        out[nm] = saw2(th, (x >= lo) & (x < hi))
    return out


def wall_diff(pa: Path, pb: Path) -> dict:
    """壁ファイル a − b の r の差: runner が使う表現 (`load_rao_wall` の RaoWall、罰則付き 5 次 B-spline) と点列の 3 次補間。"""
    wa, ia = _wall_r(pa)
    wb, ib = _wall_r(pb)
    L = min(float(wa.x_end), float(wb.x_end))
    x = np.linspace(0.0, L, 40001)
    out = {"L_a": float(wa.x_end), "L_b": float(wb.x_end)}
    for nm, (ra_, rb_) in (("fit", (wa.r(x), wb.r(x))),):
        d = np.asarray(ra_) - np.asarray(rb_)
        for lab, msk in (("all", np.ones_like(x, bool)), ("x3.4_5.3", (x >= X_SAW[0]) & (x <= X_SAW[1]))):
            j = int(np.argmax(np.abs(d[msk])))
            out[f"{nm}_{lab}"] = {"max_abs_dr": float(np.abs(d[msk][j])), "dr_at_max": float(d[msk][j]),
                                  "x_at_max": float(x[msk][j])}
    with np.load(pa) as za, np.load(pb) as zb:
        Wa, Wb = za["wall"], zb["wall"]
    xa_ok = np.all(np.diff(Wa[:, 0]) > 0) and np.all(np.diff(Wb[:, 0]) > 0)
    if xa_ok:
        xx = np.linspace(max(Wa[0, 0], Wb[0, 0]), min(Wa[-1, 0], Wb[-1, 0]), 40001)
        d = CubicSpline(Wa[:, 0], Wa[:, 1])(xx) - CubicSpline(Wb[:, 0], Wb[:, 1])(xx)
        for lab, msk in (("all", np.ones_like(xx, bool)), ("x3.4_5.3", (xx >= X_SAW[0]) & (xx <= X_SAW[1]))):
            j = int(np.argmax(np.abs(d[msk])))
            out[f"points_{lab}"] = {"max_abs_dr": float(np.abs(d[msk][j])), "dr_at_max": float(d[msk][j]),
                                    "x_at_max": float(xx[msk][j])}
    out["fit_info_a"], out["fit_info_b"] = ia["fit_info"], ib["fit_info"]
    return out


def derived(rec: dict, ig) -> dict:
    """記録から η_cs・η_wall (F/(π F_1D,nd))・流量 m* (F_cs/(isp_ratio·F_1D,nd)、runner.rao_moc_reference と同じ逆算)。"""
    i1 = ig.ideal_1d(EPS, P_A)
    F1 = float(i1["F"])
    m_line = float(rec["F_cs"]) / (float(rec["isp_ratio"]) * F1)
    return {"eta_cs": float(rec["F_cs"]) / (np.pi * F1), "eta_wall": float(rec["F_wall"]) / (np.pi * F1),
            "mstar_nd": m_line, "mstar_ratio_1d": m_line / np.pi, "F_per_mdot_ratio": float(rec["isp_ratio"])}


TABLE = (("theta_B_deg", "°", "dtheta_deg"), ("theta_E_deg", "°", "dtheta_deg"), ("theta_max_deg", "°", None),
         ("x_B", "r_t", None), ("x_D", "r_t", None), ("r_D", "r_t", None), ("M_D", "", None), ("x_E", "r_t", None),
         ("M_E", "", None), ("F_cs", "ρ*a*²r_t²", None), ("F_wall", "ρ*a*²r_t²", None), ("thrust_rel_diff", "", None),
         ("mdot_rel_diff", "", None), ("eta_cs", "", None), ("eta_wall", "", None), ("mstar_ratio_1d", "", None),
         ("F_per_mdot_ratio", "", None), ("dkappa_B", "1/r_t", None), ("g_tilde", "rad/r_t", None),
         ("closure_norm", "r_t", None))


def cmd_compare() -> dict:
    ig = gas()
    pairs = {"pin_mq320_vs_hall_mq2 (Euler run_0015 の壁と同じ m・q)": (PARTS / TAG_PIN_320, HALL_REC["mq320"]),
             "pin_mq321_vs_hall_mq321 (今の既定 m = q = n)": (PARTS / TAG_PIN, HALLREF / TAG_HALL_321)}
    cmp_, csv_rows = {}, []
    for name, (pa, pb) in pairs.items():
        ja, jb = Path(str(pa) + ".json"), Path(str(pb) + ".json")
        if not (ja.exists() and jb.exists()):
            cmp_[name] = {"skipped": f"{ja.name if not ja.exists() else jb.name} が無い"}
            continue
        a, b = json.loads(ja.read_text()), json.loads(jb.read_text())
        if a.get("status") != "ok" or b.get("status") != "ok":
            cmp_[name] = {"skipped": f"status {a.get('status')} / {b.get('status')} ({a.get('why')} / {b.get('why')})"}
            continue
        a.update(derived(a, ig))
        b.update(derived(b, ig))
        tab = {}
        for k, unit, tgt in TABLE:
            d = float(a[k]) - float(b[k])
            row = {"hall": float(b[k]), "pin": float(a[k]), "pin_minus_hall": d, "unit": unit}
            if tgt is not None:
                row["target_abs"] = TARGET[tgt]
                row["within_target"] = bool(abs(d) <= TARGET[tgt])
            tab[k] = row
            csv_rows.append({"comparison": name, "quantity": k, "unit": unit, "hall": b[k], "pin": a[k], "pin_minus_hall": d,
                             "target_abs": row.get("target_abs", ""), "within_target": row.get("within_target", "")})
        wd = wall_diff(Path(str(pa) + "_wall.npz"), Path(str(pb) + "_wall.npz"))
        for lab in ("all", "x3.4_5.3"):
            v = wd[f"fit_{lab}"]
            csv_rows.append({"comparison": name, "quantity": f"wall_r_fit_{lab}", "unit": "r_t", "hall": "", "pin": "",
                             "pin_minus_hall": v["dr_at_max"], "target_abs": TARGET["dr"],
                             "within_target": bool(v["max_abs_dr"] <= TARGET["dr"])})
        wd["within_target_fit_all"] = bool(wd["fit_all"]["max_abs_dr"] <= TARGET["dr"])
        with np.load(Path(str(pa) + "_wall.npz")) as za, np.load(Path(str(pb) + "_wall.npz")) as zb:
            saw = {"pin": be_saw(za["BE"]), "hall": be_saw(zb["BE"])}
        cmp_[name] = {"tags": {"pin": a["tag"], "hall": b["tag"]}, "table": tab, "wall_r": wd, "BE_saw2_theta": saw}
        print(f"== {name}", flush=True)
        for k, row in tab.items():
            print(f"  {k:18s} Hall {row['hall']:+.8g}  pin {row['pin']:+.8g}  差 {row['pin_minus_hall']:+.3e}"
                  + (f"  (目安 {row['target_abs']:g}: {'内' if row['within_target'] else '外'})" if "target_abs" in row else ""),
                  flush=True)
        print(f"  壁 r の差 (当てはめ) 全域 最大 {wd['fit_all']['max_abs_dr']:.3e} @ x {wd['fit_all']['x_at_max']:.3f}・"
              f"x 3.4〜5.3 最大 {wd['fit_x3.4_5.3']['max_abs_dr']:.3e} @ x {wd['fit_x3.4_5.3']['x_at_max']:.3f}", flush=True)
        print(f"  BE の θ の 2Δ 指標 x 3.4〜5.3: Hall 最大 {saw['hall']['3.4_to_5.3'].get('max_abs', float('nan')):.2e}・交代 "
              f"{saw['hall']['3.4_to_5.3'].get('alt_frac', float('nan')):.2f} / ピン 最大 "
              f"{saw['pin']['3.4_to_5.3'].get('max_abs', float('nan')):.2e}・交代 {saw['pin']['3.4_to_5.3'].get('alt_frac', float('nan')):.2f}",
              flush=True)
    # ピン版の壁ファイルが runner の読み口で読めること・V2 の問題 YAML の参照値
    v2 = {}
    try:
        w, info = ra.load_rao_wall(PARTS / f"{TAG_PIN}_wall.npz")
        v2["load_rao_wall"] = {"ok": True, "R_t": info["R_t"], "B_replaced": info["B_replaced"],
                               "last_point_minus_E": info["last_point_minus_E"], "fit_info": info["fit_info"],
                               "campaign_setup_initial_line": ((info.get("campaign") or {}).get("setup") or {}).get("initial_line")}
        from forge_design.evaluate import runner as rn
        p = load_problem(PROBLEM_PIN)
        wall, wrec = rn._build_wall_rao(p)
        with tempfile.TemporaryDirectory(prefix="rao_pin_v2_") as td:
            summ = rn._rao_wall_records(p, Path(td), wall, wrec)
            wd_json = json.loads((Path(td) / "WALL_DESIGN.json").read_text())
        ref = wd_json["moc_reference"]
        v2["problem"] = str(PROBLEM_PIN.relative_to(ROOT))
        v2["wall_design_checks"] = wd_json["checks"]
        v2["joints"] = wd_json["joints"]
        v2["moc_reference"] = ref
        v2["prepare_info_summary"] = summ
        v2["campaign_in_wall_design"] = (wd_json.get("campaign") or {}).get("setup")
        eu = json.loads((CASE / FREEZE / "WALL_DESIGN.json").read_text())["moc_reference"]
        m_eu = mdot_drift(FREEZE)["final_mdot_ratio_1d"]
        v2["preview"] = {
            "note": "ピン版 V2 の内訳 1 の見込み: ピン版の壁の Euler はまだ無いので、同じスロート・収縮部の run_0015 の Euler ṁ で代用 "
                    "(判定は親が V2 の run の後に)",
            "mdot_moc_over_mdot_euler_run0015_minus_1": ref["mdot_ratio_1d"] / m_eu - 1.0, "mdot_euler_run0015": m_eu,
            "hall_moc_reference": {k: eu.get(k) for k in ("F_cs_N", "eta_F_cs", "mdot_ratio_1d")}}
        print(f"V2 の YAML: moc_reference F_cs {ref['F_cs_N']:.1f} N・η_cs {ref['eta_F_cs']:.6f}・ṁ/ṁ_1D {ref['mdot_ratio_1d']:.7f} "
              f"(Hall 版 η_cs {eu['eta_F_cs']:.6f}・ṁ/ṁ_1D {eu['mdot_ratio_1d']:.7f})、ṁ_MOC/ṁ_Euler(run_0015) − 1 "
              f"{v2['preview']['mdot_moc_over_mdot_euler_run0015_minus_1']:+.2e}", flush=True)
    except Exception as e:  # noqa: BLE001
        v2["error"] = f"{type(e).__name__}: {e}"
        v2["traceback"] = traceback.format_exc()
        print(f"V2 の YAML の確認: ERROR {v2['error']}", flush=True)
    hall_rr = HALLREF / f"{TAG_HALL_RERUN}.json"
    rerun = None
    if hall_rr.exists():
        r = json.loads(hall_rr.read_text())
        rerun = {k: r.get(k) for k in ("vs_record_mq160", "vs_record_mq320", "time_s", "status")}
    out = write_campaign({"comparison": cmp_, "v2_problem": v2, "hall_rerun_vs_record": rerun, "target": TARGET})
    with open(OUT_CSV, "w", newline="") as fh:
        cols = ("comparison", "quantity", "unit", "hall", "pin", "pin_minus_hall", "target_abs", "within_target")
        wr = csv.DictWriter(fh, fieldnames=list(cols))
        wr.writeheader()
        for row in csv_rows:
            wr.writerow({k: (f"{row[k]:.10g}" if isinstance(row[k], float) else row[k]) for k in cols})
    print(f"{OUT.name}・{OUT_CSV.name}", flush=True)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=("lines", "design", "compare"))
    ap.add_argument("--workers", type=int, default=2)
    a = ap.parse_args()
    if a.step == "lines":
        cmd_lines()
    elif a.step == "design":
        cmd_design(a.workers)
    else:
        cmd_compare()


if __name__ == "__main__":
    main()
