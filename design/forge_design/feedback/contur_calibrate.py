r"""CONTUR の積分式の境界層の厚みを NS の δ_E に合わせる 3 係数の較正 (plan tooling-nozzle-card-local-mesh §4.3・§5.1 #6a)。

壁の境界層の厚み (半径方向の δ_r) を
    δ(x) = (c₀ + c₁·x_m)·δ_C(x) + δ₀
とする (δ_C = CONTUR の積分法 `runner_axismach.integral_delta_r` の平滑化後の値 (係数前)、x_m = 設計スロートからの実寸 [m]、
c₀ 無次元・c₁ [1/m]・δ₀ [m] で全長一定)。係数は 1 本の NS run の δ_E (停止規則が止めた場から抽出した平滑化後の値、未緩和) に
  - スロート (x = 0、列の直線補間) と出口 (最後の列) の 2 点で**等式**に合わせ、
  - 残る 1 自由度 c₁ を x ≥ 0 の列の **Δx 重み (台形) の最小二乗**で決める (モデルは c₁ に線形なので閉形式;
    codex diagnose 2026-10-15 m55-three-coef M1)。
式 (実寸、C = δ_C、E = δ_E、L = 出口の x_m、添字 t・e はスロート・出口):
    a = (E_e − E_t)/(C_e − C_t)、b(x) = x·C(x) − L·C_e·(C(x) − C_t)/(C_e − C_t)、r(x) = E(x) − E_t − a·(C(x) − C_t)
    c₁ = Σ wᵢ bᵢ rᵢ / Σ wᵢ bᵢ²、c₀ = a − c₁·L·C_e/(C_e − C_t)、δ₀ = E_t − c₀·C_t
(wᵢ は x ≥ 0 の列の台形の重み)。

入力 (run の中のファイル):
  - δ_E の表 (`delta_r_next.csv` の形式: x_rt, delta_target, delta_in_prev, delta_E, held)。探す順:
      1. `--delta-next PATH` (明示; `.gz` も可。抽出を run の外でした古い run 用)
      2. 停止規則 v2 の抽出 `STOP_RULE_STATE.json` の delta_E (状態 ok・`stop_rule_v2/delta_E_<step>/delta_r_next.csv`・sha256 を照合)
      3. run 直下の `delta_r_next.csv` (受け渡しの関門が写した表)
    状態に delta_E の記録があって ok でなければ止める (別の表に黙って落ちない)。x ≥ 0 に held・非有限があれば止める。
  - 入れた δ: 既定は表の `delta_in_prev` (= 格子の NS の壁 − Euler の壁を列ごとに測った値; δ_E と同じ格子・同じ列で測った値。
    §4.3 の 3 係数の値はこれで求めた)。`--delta-in function` で `delta_r_initial.csv` (積分法の関数の値) の 5 次補間を使う
    (両者の差 (x ≥ 0) は既定でも `delta_in_check` に記録する; case/71 run_0008 で 1.2e−4 mm、run_0010 で 1.1e−4 mm)。
  - 係数前の δ_C: `prepare_info.json` の `initializer` (= `delta_r_initial.json`; 両方あれば係数の一致を確かめる) の
    delta_scale・delta_scale_slope_per_m・delta_offset_m で δ_C = (δ_in − δ₀)/(c₀ + c₁·x_m) と戻す (`integral_delta_r` の式の逆)。
    係数が無い (pass 0 = CONTUR の式のまま) なら (1, 0, 0) で δ_C = δ_in。壁が積分法の初期壁でない run (`delta_r_csv` 等) は止める。
  - scale (= その run の r_t [m]) は `prepare_info.json` の `scale_m`。

出力 (`--out-dir`、既定 `<run>/contur_calibrate/`; `--no-write` で書かない):
  - `CONTUR_CALIBRATION.json`: 係数 (c₀・c₁ [1/m]・δ₀ [m])、今の係数、入力のパスと sha256、誤差の要約
    (今の壁: 最大 |δ_in − δ_E| と位置・x ≥ 3 m の最大・スロートと出口の差; 較正後の射影: 同じ量)、確かめたこと
  - `contur_overrides.yaml`: カードの `overrides` にそのまま貼る断片 (`deltastar_initializer.model: contur` と 3 係数;
    model を書かないと係数のキーだけでは `resolve_integral_initializer` が ValueError — codex 2026-10-15 ①)
係数は**この run の r_t の CONTUR に対する値**。新しい r_t の設計では生成器が CONTUR を新しい r_t で計算し直してから係数を当てる
(§4.3「CONTUR の r_t の依存は式が持つ」)。「較正後の射影の誤差」はこの NS の δ_E に対する誤差で、新しい壁の NS の誤差ではない。

使い方 (リポジトリの design/ で、または PYTHONPATH=design):
  design/.venv-opt/bin/python -m forge_design.feedback.contur_calibrate RUN_DIR [--delta-next PATH] [--delta-in grid|function]
      [--x-far-m 3.0] [--out-dir DIR | --no-write] [--json]
終了コード: 0 = 係数を出した、1 = 入力の不備で止めた (理由を stderr)。
関数: `calibrate_run(run_dir, ...)` (記録の dict)、`calibrate_coefficients(x_m, C_m, E_m)`・`error_summary(...)`・`locate_delta_e(...)`。
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

FORMAT = "forge.contur_calibration/v1"
PLAN = "plans/accepted/tooling-nozzle-card-local-mesh.md §4.3・§5.1 #6a"
OUT_DIRNAME = "contur_calibrate"
OUT_JSON = "CONTUR_CALIBRATION.json"
OUT_YAML = "contur_overrides.yaml"
X_FAR_M = 3.0                  # 誤差の要約の「下流側」の始まり [m] (case/71 の報告の慣例; 判定には使わない)
DELTA_IN_WARN_MM = 0.005       # 格子で測った δ_in と関数の δ_in の差 (x ≥ 0) がこれを超えたら警告 (記録だけ)
OVERRIDE_KEYS = ("deltastar_initializer.model", "deltastar_initializer.delta_scale",
                 "deltastar_initializer.delta_scale_slope_per_m", "deltastar_initializer.delta_offset_m")


def _sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def _git_head() -> str | None:
    try:
        r = subprocess.run(["git", "-C", str(Path(__file__).resolve().parent), "rev-parse", "HEAD"], capture_output=True, text=True,
                           timeout=10)
        return r.stdout.strip() or None
    except (OSError, subprocess.SubprocessError):
        return None


# ------------------------------------------------------------------ δ_E の表
def read_delta_next(path) -> dict:
    """`delta_r_next.csv` (`.gz` も可) を `deltastar_loop.read_delta_r_next` で列名で読む (列の並びの検査も同じ)。"""
    from .deltastar_loop import read_delta_r_next
    path = Path(path)
    if path.suffix != ".gz":
        return read_delta_r_next(path)
    with tempfile.TemporaryDirectory(prefix="contur_calibrate_") as td:
        tmp = Path(td) / path.name[:-3]
        with gzip.open(path, "rb") as fi, open(tmp, "wb") as fo:
            fo.write(fi.read())
        return read_delta_r_next(tmp)


def locate_delta_e(run_dir, delta_next=None) -> dict:
    """δ_E の表を探す (順はモジュールの docstring)。戻り: {path, sha256, source, state (停止規則 v2 の記録の要点 or None)}。
    停止規則の状態に delta_E の記録があって ok でない・表の sha256 が記録と違う → ValueError (別の表に落ちない)。"""
    from ..evaluate import stop_rule as SR
    run_dir = Path(run_dir)
    if delta_next is not None:
        p = Path(delta_next)
        if not p.is_file():
            raise ValueError(f"--delta-next {p} が無い")
        return {"path": str(p.resolve()), "sha256": _sha256(p), "source": "explicit (--delta-next)", "state": None}
    st_path = run_dir / SR.STATE
    if st_path.is_file():
        st = json.loads(st_path.read_text())
        de = st.get("delta_E")
        if de is not None:
            if not isinstance(de, dict) or de.get("status") != "ok":
                raise ValueError(f"{SR.STATE} の δ_E の抽出が ok でない ({(de or {}).get('status') if isinstance(de, dict) else de!r}: "
                                 f"{(de or {}).get('error') or (de or {}).get('reason') if isinstance(de, dict) else ''}) — "
                                 "抽出の失敗を別の表で埋めない (--delta-next で明示する)")
            k = de.get("step")
            dname = f"delta_E_{k}"
            cand = [Path(str(de.get("dir") or "")), run_dir / SR.V2_DIR / dname]
            od = next((c for c in cand if c.name == dname and (c / "delta_r_next.csv").is_file()), None)
            if od is None:
                raise ValueError(f"{SR.STATE} の δ_E の置き場 ({de.get('dir')!r}、{run_dir / SR.V2_DIR / dname}) に delta_r_next.csv が無い")
            p = od / "delta_r_next.csv"
            sha = _sha256(p)
            if sha != de.get("delta_r_next_sha256"):
                raise ValueError(f"{p} の sha256 {sha[:16]}… が {SR.STATE} の記録 {str(de.get('delta_r_next_sha256'))[:16]}… と違う")
            return {"path": str(p.resolve()), "sha256": sha, "source": f"停止規則 v2 の抽出 ({SR.STATE} の delta_E、step {k})",
                    "state": {"status": st.get("status"), "method": st.get("method"), "stop_step": st.get("stop_step"),
                              "delta_E_step": k, "delta_E_res_sha256": de.get("res_sha256")}}
    p = run_dir / "delta_r_next.csv"
    if p.is_file():
        return {"path": str(p.resolve()), "sha256": _sha256(p), "source": "run 直下の delta_r_next.csv", "state": None}
    raise ValueError(f"{run_dir} に δ_E の表が無い ({SR.STATE} の delta_E・run 直下の delta_r_next.csv のどちらも無い; --delta-next で与える)")


# ------------------------------------------------------------------ 今の係数 (入れた δ の式)
def coefficient_record(run_dir) -> dict:
    """run の壁の δ の式: prepare_info.json の scale_m・initializer (= delta_r_initial.json) の係数。
    戻り: {scale_m, c0, c1_per_m, d0_m, pass0 (係数なし), initializer (設定の要点), files (パスと sha256), contur_ref (係数前の
    スロート・出口の記録 [r_t])}。積分法の初期壁でない run は ValueError。"""
    run_dir = Path(run_dir)
    pi_path = run_dir / "prepare_info.json"
    if not pi_path.is_file():
        raise ValueError(f"{pi_path} が無い")
    info = json.loads(pi_path.read_text())
    S = float(info["scale_m"])
    files = {"prepare_info.json": _sha256(pi_path)}
    ini_pi = info.get("initializer")
    js_path = run_dir / "delta_r_initial.json"
    ini_js = json.loads(js_path.read_text()) if js_path.is_file() else None
    if ini_js is not None:
        files["delta_r_initial.json"] = _sha256(js_path)
    ini = ini_pi if isinstance(ini_pi, dict) else ini_js
    if not isinstance(ini, dict):
        raise ValueError(f"{run_dir} の壁は積分法の初期壁 (CONTUR の式) でない (prepare_info.json の initializer・delta_r_initial.json が無い; "
                         f"dstar_source {info.get('dstar_source')!r}) — 係数前の δ_C を作れない")
    if info.get("delta_r_csv"):
        raise ValueError(f"{run_dir} の壁は delta_r_csv ({info.get('delta_r_csv')}) で作った — CONTUR の式の壁でない")

    def coefs(d):
        return (float(d.get("delta_scale", 1.0)), float(d.get("delta_scale_slope_per_m", 0.0)), float(d.get("delta_offset_m", 0.0)))

    c = coefs(ini)
    if isinstance(ini_pi, dict) and isinstance(ini_js, dict) and coefs(ini_pi) != coefs(ini_js):
        raise ValueError(f"prepare_info.json の initializer の係数 {coefs(ini_pi)} と delta_r_initial.json の係数 {coefs(ini_js)} が違う")
    tc = ini.get("two_coef") if isinstance(ini.get("two_coef"), dict) else None
    if tc is not None:
        if (float(tc.get("delta_scale", 1.0)), float(tc.get("delta_scale_slope_per_m", 0.0)), float(tc.get("delta_offset_m", 0.0))) != c:
            raise ValueError(f"initializer.two_coef の係数 {tc} が initializer の係数 {c} と違う")
        if tc.get("scale_m") is not None and float(tc["scale_m"]) != S:
            raise ValueError(f"initializer.two_coef.scale_m {tc['scale_m']!r} が prepare_info.json の scale_m {S!r} と違う")
    if ini.get("rt_m") is not None and float(ini["rt_m"]) != S:
        raise ValueError(f"initializer.rt_m {ini['rt_m']!r} が prepare_info.json の scale_m {S!r} と違う")
    pass0 = c == (1.0, 0.0, 0.0)
    if pass0:
        ref = {"throat_rt": ini.get("delta_r_throat"), "exit_rt": ini.get("delta_r_exit"), "from": "initializer.delta_r_throat・delta_r_exit"}
    else:
        ref = {"throat_rt": (tc or {}).get("delta_r_throat_contur"), "exit_rt": (tc or {}).get("delta_r_exit_contur"),
               "from": "initializer.two_coef.delta_r_throat_contur・delta_r_exit_contur"}
    keep = ("model", "closure", "closure_version", "cf_scale", "n_scale", "a_crocco", "thermal_bc", "x_virtual_m", "theta0_m")
    return {"scale_m": S, "c0": c[0], "c1_per_m": c[1], "d0_m": c[2], "pass0": pass0,
            "initializer": {k: ini.get(k) for k in keep if k in ini}, "files": files, "contur_ref": ref,
            "x_E": info.get("x_E"), "x_wall_end_rt": ((info.get("physical_wall") or {}).get("domain") or [None, None])[1]}


def function_delta_in(run_dir, x_rt) -> np.ndarray | None:
    """`delta_r_initial.csv` (`.gz` も可; x_rt, delta_r_rt, …) の δ (係数を当てた後の関数の値) を 5 次補間で x_rt に評価 [r_t]。無ければ None。"""
    from scipy.interpolate import make_interp_spline
    run_dir = Path(run_dir)
    for name in ("delta_r_initial.csv", "delta_r_initial.csv.gz"):
        p = run_dir / name
        if p.is_file():
            op = gzip.open if name.endswith(".gz") else open
            with op(p, "rt") as fh:
                a = np.loadtxt(fh, delimiter=",", skiprows=1)
            spl = make_interp_spline(a[:, 0], a[:, 1], k=5)
            return spl(np.clip(np.asarray(x_rt, dtype=float), a[0, 0], a[-1, 0]))
    return None


# ------------------------------------------------------------------ 較正と誤差
def trapezoid_weights(x) -> np.ndarray:
    """列の x に対する台形の積分の重み (端は半分)。"""
    x = np.asarray(x, dtype=float)
    w = np.zeros_like(x)
    if len(x) >= 2:
        dx = np.diff(x)
        w[:-1] += 0.5 * dx
        w[1:] += 0.5 * dx
    return w


def calibrate_coefficients(x_m, C_m, E_m) -> dict:
    """δ = (c₀ + c₁·x)·C + δ₀ をスロート (x = 0、列の直線補間) と出口 (最後の列) の E に等式で合わせ、c₁ を x ≥ 0 の列の台形の重みの
    最小二乗で決める (閉形式; モジュールの docstring の式)。x_m [m] は昇順で 0 を挟むこと。戻り: {c0, c1_per_m, d0_m, 端点の値, 重みの説明}。"""
    x = np.asarray(x_m, dtype=float); C = np.asarray(C_m, dtype=float); E = np.asarray(E_m, dtype=float)
    if not (len(x) >= 3 and np.all(np.diff(x) > 0)):
        raise ValueError("列の x が昇順でない・点が足りない")
    if not (x[0] <= 0.0 < x[-1]):
        raise ValueError(f"列の x [{x[0]:.6g}, {x[-1]:.6g}] m が設計スロート x = 0 を挟まない")
    Ct = float(np.interp(0.0, x, C)); Et = float(np.interp(0.0, x, E))
    Ce = float(C[-1]); Ee = float(E[-1]); L = float(x[-1])
    if not (L > 0.0 and np.isfinite([Ct, Et, Ce, Ee]).all()) or Ce == Ct:
        raise ValueError(f"端点の値が不正 (L {L!r}、C_t {Ct!r}、C_e {Ce!r}、E_t {Et!r}、E_e {Ee!r}; C_e = C_t では解けない)")
    m = x >= 0.0
    xs, Cs, Es = x[m], C[m], E[m]
    if not (np.isfinite(Cs).all() and np.isfinite(Es).all()):
        raise ValueError("x ≥ 0 に非有限の δ_C・δ_E がある")
    w = trapezoid_weights(xs)
    a = (Ee - Et) / (Ce - Ct)
    b = xs * Cs - L * Ce * (Cs - Ct) / (Ce - Ct)
    r = Es - Et - a * (Cs - Ct)
    den = float(np.sum(w * b * b))
    if not den > 0.0:
        raise ValueError("最小二乗の分母 Σ w b² が 0 (c₁ を決められない)")
    c1 = float(np.sum(w * b * r) / den)
    c0 = float(a - c1 * L * Ce / (Ce - Ct))
    d0 = float(Et - c0 * Ct)
    return {"c0": c0, "c1_per_m": c1, "d0_m": d0, "throat": {"C_m": Ct, "E_m": Et}, "exit": {"x_m": L, "C_m": Ce, "E_m": Ee},
            "n_lsq": int(m.sum()), "weights": "x ≥ 0 の列の台形 (Δx) の重み", "constraints": "x = 0 (列の直線補間) と最後の列で等式"}


def model_delta(x_m, C_m, c0, c1, d0) -> np.ndarray:
    """δ = (c₀ + c₁·x)·C + δ₀ [m]。"""
    x = np.asarray(x_m, dtype=float)
    return (c0 + c1 * x) * np.asarray(C_m, dtype=float) + d0


def error_summary(x_m, delta_m, E_m, x_far_m: float = X_FAR_M) -> dict:
    """x ≥ 0 の |δ − δ_E| の最大と位置 (と符号つきの値 δ − δ_E)、x ≥ x_far の最大、スロート (x = 0) と出口 (最後の列) の差、
    Δx 重みの rms。単位 mm (位置は m)。"""
    x = np.asarray(x_m, dtype=float)
    e = np.asarray(delta_m, dtype=float) - np.asarray(E_m, dtype=float)
    m = x >= 0.0
    xs, es = x[m], e[m]
    i = int(np.argmax(np.abs(es)))
    far = x >= float(x_far_m)
    w = trapezoid_weights(xs)
    out = {"max_abs_mm": float(abs(es[i]) * 1e3), "x_m_at_max": float(xs[i]), "signed_at_max_mm": float(es[i] * 1e3),
           f"max_abs_x_ge_{x_far_m:g}m_mm": float(np.max(np.abs(e[far])) * 1e3) if far.any() else None,
           "throat_mm": float(np.interp(0.0, x, e) * 1e3), "exit_mm": float(e[-1] * 1e3),
           "rms_dx_mm": float(np.sqrt(np.sum(w * es * es) / np.sum(w)) * 1e3) if np.sum(w) > 0 else None,
           "sign": "δ − δ_E (正 = 入れた δ が NS より厚い)"}
    out["x_far_m"] = float(x_far_m)
    return out


def calibrate_run(run_dir, delta_next=None, delta_in: str = "grid", x_far_m: float = X_FAR_M) -> dict:
    """run の δ_E と今の壁の δ から 3 係数を求め、記録の dict を返す (ファイルは書かない; `write_outputs`)。
    delta_in: "grid" (既定、表の delta_in_prev) / "function" (delta_r_initial.csv の 5 次補間)。"""
    run_dir = Path(run_dir).resolve()
    if delta_in not in ("grid", "function"):
        raise ValueError(f"delta_in は grid / function ({delta_in!r})")
    loc = locate_delta_e(run_dir, delta_next)
    tbl = read_delta_next(loc["path"])
    co = coefficient_record(run_dir)
    S = co["scale_m"]
    x_rt = np.asarray(tbl["x_rt"], dtype=float)
    E_rt = np.asarray(tbl["delta_E"], dtype=float)
    din_g = np.asarray(tbl["delta_in_prev"], dtype=float)
    held = np.asarray(tbl["held"], dtype=float)
    use = x_rt >= 0.0
    if np.any(held[use] != 0):
        raise ValueError(f"δ_E の表の x ≥ 0 に held の列が {int(np.sum(held[use] != 0))} 本ある — 抽出が欠けた δ_E で較正しない")
    if not (np.isfinite(E_rt[use]).all() and np.isfinite(din_g[use]).all()):
        raise ValueError("δ_E の表の x ≥ 0 に非有限の δ_E・delta_in_prev がある")
    din_f = function_delta_in(run_dir, x_rt)
    check = {"function_source": None}
    if din_f is not None:
        dd = float(np.max(np.abs(din_f[use] - din_g[use])) * S * 1e3)
        check = {"max_abs_grid_minus_function_x_ge_0_mm": dd, "warn_mm": DELTA_IN_WARN_MM, "warn": dd > DELTA_IN_WARN_MM,
                 "function_source": "delta_r_initial.csv の 5 次補間"}
    if delta_in == "function":
        if din_f is None:
            raise ValueError(f"--delta-in function: {run_dir} に delta_r_initial.csv が無い")
        din_rt = din_f
    else:
        din_rt = din_g
    x_m = x_rt * S
    din_m = din_rt * S
    E_m = E_rt * S
    c0, c1, d0 = co["c0"], co["c1_per_m"], co["d0_m"]
    # 係数前の δ_C (integral_delta_r の δ = (c₀ + c₁·x_m)·δ_C + δ₀ の逆)。使う列 (x ≥ 0 とスロートを挟む 1 本) で倍率が正であること
    fac = c0 + c1 * x_m
    i0 = int(np.searchsorted(x_m, 0.0))
    need = np.zeros(len(x_m), dtype=bool)
    need[max(i0 - 1, 0):] = True
    if np.any(fac[need] <= 0.0):
        raise ValueError(f"今の係数の倍率 c₀ + c₁·x が x ≥ 0 で正でない (最小 {float(np.min(fac[need])):.6g}) — δ_C に戻せない")
    C_m = np.full_like(x_m, np.nan)
    C_m[need] = (din_m[need] - d0) / fac[need]
    cal = calibrate_coefficients(x_m[need], C_m[need], E_m[need])
    n0, n1, nd = cal["c0"], cal["c1_per_m"], cal["d0_m"]
    new_fac = n0 + n1 * x_m[use]
    warns = []
    if check.get("warn"):
        warns.append(f"格子で測った δ_in と関数の δ_in の差 (x ≥ 0) が {check['max_abs_grid_minus_function_x_ge_0_mm']:.3g} mm "
                     f"(> {DELTA_IN_WARN_MM} mm)")
    if not n0 > 0.0:
        warns.append(f"新しい c₀ = {n0!r} ≤ 0 — 生成器 (runner_axismach._two_coef) は delta_scale ≤ 0 を拒む")
    if np.any(new_fac <= 0.0):
        warns.append(f"新しい倍率 c₀ + c₁·x が x ≥ 0 で正でない (最小 {float(np.min(new_fac)):.6g})")
    if float(co["initializer"].get("cf_scale") or 1.0) != 1.0 or (co["initializer"].get("closure_version") or "contur_v1") != "contur_v1":
        warns.append(f"基底の CONTUR の設定が既定でない ({ {k: co['initializer'].get(k) for k in ('cf_scale', 'closure_version')} }) — "
                     "新しい係数は同じ設定の CONTUR に当てる前提 (カード・問題 YAML に残すこと)")
    ref = co["contur_ref"]
    Cchk = {}
    for nm, xv, rv in (("throat", 0.0, ref.get("throat_rt")), ("exit", x_m[-1], ref.get("exit_rt"))):
        cv = float(np.interp(xv, x_m[need], C_m[need]))
        Cchk[nm] = {"delta_C_from_delta_in_mm": cv * 1e3, "recorded_contur_mm": None if rv is None else float(rv) * S * 1e3,
                    "diff_mm": None if rv is None else (cv - float(rv) * S) * 1e3}
    Cchk["note"] = f"δ_in から戻した δ_C と記録の CONTUR の値 ({ref.get('from')}) の差 (格子の測定の差の大きさ)"
    cur = error_summary(x_m, din_m, E_m, x_far_m)
    proj = error_summary(x_m[need], model_delta(x_m[need], C_m[need], n0, n1, nd), E_m[need], x_far_m)
    overrides = {OVERRIDE_KEYS[0]: "contur", OVERRIDE_KEYS[1]: n0, OVERRIDE_KEYS[2]: n1, OVERRIDE_KEYS[3]: nd}
    x_end = co.get("x_wall_end_rt")
    return {
        "format": FORMAT, "plan": PLAN, "created": datetime.now(timezone.utc).isoformat(timespec="seconds"), "git_head": _git_head(),
        "run": str(run_dir),
        "inputs": {"delta_E": loc, "delta_E_columns": f"delta_E (未緩和の平滑化抽出; {tbl['header_kind']} header)",
                   "delta_in": {"source": delta_in, "column": "delta_in_prev (格子の NS の壁 − Euler の壁、列ごと)" if delta_in == "grid"
                                else "delta_r_initial.csv の 5 次補間"},
                   "files": co["files"], "scale_m": S, "n_cols": int(len(x_m)), "x_range_m": [float(x_m[0]), float(x_m[-1])],
                   "x_wall_end_m": None if x_end is None else float(x_end) * S, "initializer": co["initializer"]},
        "current": {"c0": c0, "c1_per_m": c1, "d0_m": d0, "pass0": co["pass0"],
                    "form": "δ = (c₀ + c₁·x_m)·δ_C + δ₀" + (" (係数なし = CONTUR の式のまま)" if co["pass0"] else "")},
        "coefficients": {"c0": n0, "c1_per_m": n1, "d0_m": nd},
        "calibration": {k: v for k, v in cal.items() if k not in ("c0", "c1_per_m", "d0_m")},
        "overrides": overrides,
        "overrides_yaml": overrides_yaml(overrides, run_dir, S),
        "errors": {"current_wall": cur, "projected_after_calibration": proj,
                   "note": "「較正後の射影」はこの NS の δ_E に対する誤差 (新しい壁の NS の誤差ではない)"},
        "checks": {"delta_in": check, "delta_C_vs_recorded_contur": Cchk, "n_held_x_ge_0": 0},
        "warnings": warns,
    }


def yaml_float(v) -> str:
    """float を PyYAML (YAML 1.1) が**数として**読む文字列にする (桁は repr のまま、往復で同じ float)。repr の `-5e-05` のように
    仮数に小数点の無い指数の形は、PyYAML の float の規則 (小数点が要る) に合わず文字列と読まれるので `-5.0e-05` にする
    (plan tooling-nozzle-card-local-mesh §5.1 #6d)。有限でない値と、読み直して同じ float にならない値は ValueError。"""
    import yaml
    x = float(v)
    if not np.isfinite(x):
        raise ValueError(f"yaml_float: 有限でない値 {x!r}")
    r = repr(x)
    m, sep, e = r.partition("e")
    if sep and "." not in m:
        r = f"{m}.0e{e}"
    back = yaml.safe_load(r)
    if not isinstance(back, float) or back != x:
        raise ValueError(f"yaml_float: {r!r} を YAML で読み直すと {back!r} ({type(back).__name__}) で {x!r} に戻らない")
    return r


def overrides_yaml(overrides: dict, run_dir, scale_m: float) -> str:
    """カードの `overrides` に貼る断片 (YAML)。数は `yaml_float` (PyYAML が数として読み、往復で同じ float)。書いた断片を読み直して
    overrides と同じ値・型 (数は float) に戻ることを確かめる (戻らなければ RuntimeError)。"""
    import yaml
    lines = [f"# contur_calibrate ({PLAN}): {Path(run_dir).name} (r_t {scale_m!r} m) の δ_E に合わせた 3 係数",
             "# δ = (c₀ + c₁·x_m)·δ_CONTUR + δ₀ (x_m は設計スロートからの実寸 [m])。新しい r_t でも同じ係数を CONTUR に当てる",
             "overrides:"]
    for k in OVERRIDE_KEYS:
        v = overrides[k]
        lines.append(f"  {k}: {v if isinstance(v, str) else yaml_float(v)}")
    text = "\n".join(lines) + "\n"
    back = (yaml.safe_load(text) or {}).get("overrides") or {}
    def same(k):
        v, w = overrides[k], back.get(k)
        return w == v if isinstance(v, str) else (isinstance(w, float) and w == float(v))
    bad = [k for k in OVERRIDE_KEYS if not same(k)]
    if bad:
        raise RuntimeError(f"overrides の断片を読み直すと値・型が違う: {[(k, back.get(k), overrides[k]) for k in bad]}")
    return text


def write_outputs(rec: dict, out_dir) -> dict:
    """記録 (JSON) と overrides の断片 (YAML) を書く。戻り: {json, yaml} のパス。"""
    od = Path(out_dir)
    od.mkdir(parents=True, exist_ok=True)
    pj, py = od / OUT_JSON, od / OUT_YAML
    pj.write_text(json.dumps(rec, indent=1, ensure_ascii=False))
    py.write_text(rec["overrides_yaml"])
    return {"json": str(pj), "yaml": str(py)}


def summary_lines(rec: dict) -> list[str]:
    co, cu, er = rec["coefficients"], rec["current"], rec["errors"]
    a, b = er["current_wall"], er["projected_after_calibration"]
    xf = a["x_far_m"]
    k = f"max_abs_x_ge_{xf:g}m_mm"
    L = [f"run: {rec['run']}",
         f"δ_E: {rec['inputs']['delta_E']['source']} — {rec['inputs']['delta_E']['path']}",
         f"今の係数: c₀ {cu['c0']!r}・c₁ {cu['c1_per_m']!r} /m・δ₀ {cu['d0_m']!r} m" + (" (pass 0)" if cu["pass0"] else ""),
         f"新しい係数: c₀ {co['c0']!r}・c₁ {co['c1_per_m']!r} /m・δ₀ {co['d0_m']!r} m",
         f"今の壁の誤差 |δ_in − δ_E|: 最大 {a['max_abs_mm']:.4f} mm @ x {a['x_m_at_max']:.3f} m、x ≥ {xf:g} m {a[k]:.4f} mm、"
         f"スロート {a['throat_mm']:+.4f}・出口 {a['exit_mm']:+.4f} mm",
         f"較正後の射影 (この NS の δ_E に対して): 最大 {b['max_abs_mm']:.4f} mm @ x {b['x_m_at_max']:.3f} m、x ≥ {xf:g} m {b[k]:.4f} mm、"
         f"スロート {b['throat_mm']:+.2e}・出口 {b['exit_mm']:+.2e} mm"]
    chk = rec["checks"]["delta_in"]
    if chk.get("function_source"):
        L.append(f"δ_in の照合 (格子 − 関数, x ≥ 0): 最大 {chk['max_abs_grid_minus_function_x_ge_0_mm']:.2e} mm")
    L += [f"警告: {w}" for w in rec["warnings"]]
    L += ["カードの overrides の断片:", rec["overrides_yaml"].rstrip()]
    return L


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        prog="python -m forge_design.feedback.contur_calibrate",
        description="CONTUR の積分式の δ を NS の δ_E に合わせる 3 係数 δ = (c₀ + c₁·x_m)·δ_C + δ₀ を求める "
                    "(スロートと出口で等式、c₁ は x ≥ 0 の Δx 重みの最小二乗; plan tooling-nozzle-card-local-mesh §4.3・§5.1 #6a)。",
        epilog="出力: <out-dir>/CONTUR_CALIBRATION.json (係数・誤差の要約・入力の sha256) と contur_overrides.yaml (カードの overrides の断片)。"
               "δ_E は --delta-next → 停止規則 v2 の抽出 (STOP_RULE_STATE.json の delta_E) → run 直下の delta_r_next.csv の順に探す。"
               "終了コード 0 = 係数を出した、1 = 入力の不備で止めた。")
    ap.add_argument("run", help="NS の run のディレクトリ (prepare_info.json・delta_r_initial.json と δ_E の表)")
    ap.add_argument("--delta-next", default=None, metavar="CSV", help="δ_E の表 (delta_r_next.csv の形式、.gz 可) を明示する")
    ap.add_argument("--delta-in", default="grid", choices=("grid", "function"),
                    help="入れた δ: grid (既定、表の delta_in_prev = 格子で測った値) / function (delta_r_initial.csv の関数の値)")
    ap.add_argument("--x-far-m", type=float, default=X_FAR_M, help=f"誤差の要約の下流側の始まり [m] (既定 {X_FAR_M}; 記録だけ)")
    ap.add_argument("--out-dir", default=None, help=f"出力の置き場 (既定 <run>/{OUT_DIRNAME})")
    ap.add_argument("--no-write", action="store_true", help="ファイルを書かない (標準出力だけ)")
    ap.add_argument("--json", action="store_true", help="記録の JSON を標準出力に出す")
    a = ap.parse_args(argv)
    try:
        rec = calibrate_run(a.run, delta_next=a.delta_next, delta_in=a.delta_in, x_far_m=a.x_far_m)
    except (ValueError, OSError, KeyError) as e:
        print(f"contur_calibrate: 止めた — {e}", file=sys.stderr)
        return 1
    if a.json:
        print(json.dumps(rec, indent=1, ensure_ascii=False))
    else:
        print("\n".join(summary_lines(rec)))
    if not a.no_write:
        paths = write_outputs(rec, a.out_dir if a.out_dir else Path(a.run) / OUT_DIRNAME)
        print(f"書いた: {paths['json']}、{paths['yaml']}", file=sys.stderr if a.json else sys.stdout)
    return 0


if __name__ == "__main__":
    sys.exit(main())
