#!/usr/bin/env python3
r"""③ベル (M) モード §5.1 #5 (h) — T4 (停留性) の評価器の A/B (事前登録): Hall の初期線をやめ、n 321 の Hall マーチの
整合な C⁻ 前線 k₀ を r 等間隔に N ∈ {321, 641, 1281} へ 3 次補間した初期線で、(M) の壁の推力の摂動応答を測る
(数値を出すだけ、判定 A/B/C は親と上位)。

plan: plans/active/tooling-nozzle-bell-two-mode.md §5.1 #5 (h)・§6 V1 の T4。
判断の記録: notes/reviews/2026-10-12-bell-rao-t4-stationarity-diagnostician.md (「次の一手」と「やらない」)。

**壁**: `rao_t6_parts/L100_n321_mq2_wall.npz` の BE (n 321・Goursat の m = q 倍の生成器の壁点列) を、前任 (`run_rao_t4.py`
の `base_wall`) と同じ `RaoWall` の罰則付き 5 次 B-spline (48 区間・λ GCV、B で C1) にしたもの。点列の末尾は生成器の E
(`L100_n321_mq2.json` の x_E・r_E・θ_E) に置き換える。

**初期線** (`rao_axisym.resampled_front_line`): その壁と Hall n 321 の初期線で前進 MOC (`forward_moc`、analytic・converge) を
L まで回し、壁点が x ≈ 0.3 (`X_FRONT`、円弧 kernel の上 — B は x 0.496) の前線 k₀ を取る (軸まで ST_OK を確かめる)。
前線の x・θ・M を r の 3 次スプライン (not-a-knot) にして r 等間隔の N+1 点で評価、軸端の θ_r は `net.thr_axis[k₀]`。
B の上流なので BE の摂動は初期線に効かない。前線は円弧の上だけに依るので L70 の壁でも同じ (follow 段で確かめる)。

**推力**: `rao_axisym.forward_thrust(wall, L, gas, line, splits=(x_B,))` — 前進 MOC で場を解き直し、F = 初期線の運動量 +
圧力 + 壁の圧力の軸力 (x₀ = 前線の壁点 → L)。会計は 2 つ: **スプライン** (`F`: 壁点の p の区間ごとの 3 次スプライン ×
壁の式の r・r'、Gauss–Legendre) と **台形** (`F_trap`: 壁点の台形 [r で])。

**摂動**: 族 `sin_c1` = sin(πt) sin(kπt) (t = B → E の弧長の正規化、B・E の位置と接線を保つ)、ε ∈ {0, ±1, ±2, ±4}×10⁻³ r_t
の 7 点。F(ε) を u = ε/ε₀ (ε₀ = 4e-3) の 3 次で最小二乗 → c₁ε₀・c₂ε₀²、ε* = |c₁/(2c₂)| [r_t]、懸かる推力
c₁²/(4|c₂|) / F_1D (F_1D = 同じ凍結 TP の 1D 理想推力、A* = π r_t²)。両会計で別々に当てはめる。

**事前登録の判定** (diagnostician 2026-10-12、plan §5.1 #5 (h)。ここでは数値と条件の真偽を出すだけ):
- **A**: N 1281 で |c₁ε₀| ≤ 1e-7 かつ両会計の c₁ε₀ の差 ≤ 1e-7 かつ 641 → 1281 で |c₁ε₀| が減る → `--follow` で k2・k3・L70・
  `ctrl_plus` を N 1281 で回す。
- **B**: 両会計の差 ≤ 1e-7 だが |c₁ε₀| が 641 と 1281 で 30 % 以内に一致し両方 > 1e-7 (壁は ε* だけ非停留)。
- **C**: N 1281 でも 2 会計の差 > 1e-7 (評価器はこの線を解像できない)。
登録文は |c₁ε₀| の会計を指定していないので、条件の真偽は会計ごとに出す。「30 % 以内」は |c₁ε₀(1281)/c₁ε₀(641) − 1| ≤ 0.3 で
数える (符号も一致する読み)。

**参照の行** (A/B の判定には使わない、解釈の材料): `front` = 前線 k₀ そのもの (321 点、r 非等間隔) を初期線にしたもの
(再標本化なしの整合な線)、`hall` = Hall n 321 の初期線 (前任の評価器、同じ壁)。

**再現の検査** (最初に確かめる項目、`reproduction`): 前線 k₀ の状態、各 N の線の `_check_line` の結果、再標本化した線を通る
3 次スプラインを前線の点で評価した往復の誤差 (x・θ・M の最大と位置)、前線そのものを初期線にしたマーチが元の網の前線を
再現するか。

**follow 段** (`--follow`、判定 A のときだけ回す): N 1281 で L100 の k2・k3、L70 (`rao_t6_parts/L70_n321_mq2_wall.npz`) の
k1・k2・k3、`ctrl_plus` (L100・L70) の k1・k2・k3。`ctrl_plus` の壁は前任の `run_rao_t4_gen.py` と同じ作り方 (`_control`:
根の BD の上で σ_D を 1.02 倍した D⁺ から Rao 線を引き直し → Goursat → 壁、x = L で切る) で、生成器を n 321・m = q = 320
(t6 の mq2 と同じ設定) で回し直して作る (1 本 約 26 分)。生成し直した基準の壁が t6 の mq2 の BE と一致するかを記録する。

出力 (このディレクトリ): `rao_t4h_parts/` (ジョブごと、再開用)、`rao_t4h.csv` (壁 × 初期線 × N × 族 × k の 1 行、両会計の列)、
`rao_t4h.json` (全 F(ε)・当てはめ・判定の数値・再現の検査・来歴)。既存の `rao_t4*.csv/json` は上書きしない。
実行: `nice -n 10 /home/sano/work/forge/design/.venv-opt/bin/python case/68.bell_two_mode/rao_m/run_rao_t4h.py`
(2 並列まで)。`--follow` で follow 段を足す、`--only` (ジョブ名の部分一致)、`--assemble-only`。
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
import time  # noqa: E402
import traceback  # noqa: E402
import warnings  # noqa: E402
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
from scipy.interpolate import CubicSpline  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "design"))
sys.path.insert(0, str(HERE))
from forge_design.gas.rocket_exhaust import RocketExhaust  # noqa: E402
from forge_design.geometry import moc_forward as mf  # noqa: E402
from forge_design.geometry import rao_axisym as ra  # noqa: E402

EXHAUST = HERE.parent / "exhaust_lox_lh2_pc100_of6.yaml"
T6 = HERE / "rao_t6_parts"
PARTS = HERE / "rao_t4h_parts"
EPS = 50.0
R_T = 1.0
P_A = 0.0
EPS0 = 4e-3
AMPS = (1e-3, 2e-3, 4e-3)
EPSS = tuple(sorted([-a for a in AMPS] + [0.0] + list(AMPS)))
FAMILY = "sin_c1"
N_LIST = (321, 641, 1281)
N_HALL = 321
X_FRONT = 0.3                         # 前線 k₀ の壁点の x (円弧 kernel の上、B より上流)
LINE_ABS = 1e-7                       # 登録の線 (|c₁ε₀|・両会計の差)
B_REL = 0.30                          # 判定 B の「641 と 1281 で 30 % 以内」
CTRL_SHIFT = 0.02                     # ctrl_plus: σ_D を 1 + 0.02 倍 (run_rao_t4_gen.py と同じ)
WALLS = {"L100": "L100_n321_mq2", "L70": "L70_n321_mq2"}
SOURCES = ("design/forge_design/geometry/rao_axisym.py", "design/forge_design/geometry/moc_forward.py",
           "design/forge_design/geometry/moc_kernel.py", "design/forge_design/gas/isentrope.py",
           "case/68.bell_two_mode/rao_m/run_rao_t4_gen.py")


def git_info() -> dict:
    head = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "-C", str(ROOT), "status", "--porcelain", "design/forge_design", "design/tests"],
                           capture_output=True, text=True).stdout
    sha = {s: hashlib.sha256((ROOT / s).read_bytes()).hexdigest()[:16] for s in SOURCES}
    return {"commit": head, "design_tree_dirty": [ln for ln in dirty.splitlines() if ln.strip()], "source_sha256_16": sha}


def _jsonable(v):
    if isinstance(v, dict):
        return {str(k): _jsonable(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [_jsonable(x) for x in v]
    if isinstance(v, (bool, np.bool_)):
        return bool(v)
    if isinstance(v, (np.floating, float)):
        return None if not np.isfinite(v) else float(v)
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, np.ndarray):
        return _jsonable(v.tolist())
    return v


_CACHE: dict = {}


def gas():
    if "ig" not in _CACHE:
        _CACHE["ig"] = ra.as_isentrope(RocketExhaust.from_file(EXHAUST).moc_gas())
    return _CACHE["ig"]


def hall(n: int) -> dict:
    key = ("hall", n)
    if key not in _CACHE:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            _CACHE[key] = ra._hall_line(R_T, gas(), n)
    return _CACHE[key]


def F1D_abs() -> float:
    """1D 理想推力 (同じ凍結 TP で ε まで、A* = π r_t²)。相対量の分母 (run_rao_t4.py と同じ)。"""
    return float(np.pi * gas().ideal_1d(EPS, P_A)["F"])


# ================================================================================================ 壁
def _rao_wall(theta_B: float, BE: np.ndarray, E_xrth, x_end_if_beyond: float | None = None):
    """run_rao_t4.py の `base_wall` と同じ: BE の点列の末尾を E に置き換え、RaoWall の pspline (既定)。
    `x_end_if_beyond` を与え E の x がそれより先なら、そこで切る (ctrl_plus)。"""
    P = np.asarray(BE, dtype=float)[:, :2].copy()
    last_gap = (float(P[-1, 0] - E_xrth[0]), float(P[-1, 1] - E_xrth[1]))
    P[-1] = (float(E_xrth[0]), float(E_xrth[1]))
    x_end = x_end_if_beyond if (x_end_if_beyond is not None and E_xrth[0] > x_end_if_beyond) else None
    w = ra.RaoWall(R_T, float(theta_B), P, float(E_xrth[2]), x_end=x_end)
    return w, last_gap


def load_wall(wall: str):
    """`L100` / `L70` (t6 の mq2 の壁)、`L100_ctrl_plus` / `L70_ctrl_plus` (follow 段の生成の部品)。戻り (RaoWall, 情報)。"""
    if wall in WALLS:
        tag = WALLS[wall]
        rec = json.loads((T6 / f"{tag}.json").read_text())
        z = np.load(T6 / f"{tag}_wall.npz")
        E = (rec["x_E"], rec["r_E"], np.deg2rad(rec["theta_E_deg"]))
        w, gap = _rao_wall(np.deg2rad(rec["theta_B_deg"]), z["BE"], E)
        info = {"source": f"rao_t6_parts/{tag}", "L_wall": float(w.x_end), "n_gen": rec["n"], "m_gen": rec["m"], "q_gen": rec["q"],
                "F_cs_gen": rec["F_cs"], "F_wall_gen": rec["F_wall"]}
    else:
        base, which = wall.split("_", 1)
        p = PARTS / f"gen_{base}_rec.json"
        z = PARTS / f"gen_{base}_walls.npz"
        if not (p.exists() and z.exists()):
            raise FileNotFoundError(f"{wall}: 生成の部品 gen_{base} が無い (follow 段の gen ジョブを先に)")
        rec = json.loads(p.read_text())
        arr = np.load(z)
        if f"{which}_BE" not in arr:
            raise FileNotFoundError(f"{wall}: {which} が無い ({(rec.get('controls') or {}).get(which, {}).get('why')})")
        L = float(rec["L_rt"])
        w, gap = _rao_wall(rec["theta_B"], arr[f"{which}_BE"], arr[f"{which}_E_xrthM"][:3], x_end_if_beyond=L)
        info = {"source": f"rao_t4h_parts/gen_{base} ({which})", "L_wall": float(w.x_end), "n_gen": rec["n"], "m_gen": rec["m"],
                "q_gen": rec["q"], "F_cs_gen": None, "F_wall_gen": None}
    info.update({"x_B": w.x_B, "theta_end_deg": float(np.rad2deg(w.theta_end)), "eps_wall": float(w.r_end) ** 2,
                 "B_replaced": w.B_replaced, "last_point_minus_E": gap, "fit_info": w.fit_info})
    return w, info


# ================================================================================================ 初期線
def front(wall: str):
    """Hall n 321 でその壁の前進 MOC を L まで回し、壁点 x ≈ X_FRONT の前線 k₀ を返す (網・k₀・情報)。"""
    key = ("front", wall)
    if key not in _CACHE:
        w, info = load_wall(wall)
        L = info["L_wall"]
        res = mf.forward_moc(ra.ExtendedWall(w, L), L, gas().moc, hall(N_HALL))
        net, n = res["net"], res["n"]
        xw = net.a["x"][:, n]
        k0 = int(np.nanargmin(np.abs(xw - X_FRONT)))
        st = net.status[k0]
        h = hashlib.sha256(b"".join(np.ascontiguousarray(net.a[f][k0]).tobytes() for f in net.FIELDS)
                           + np.float64(net.thr_axis[k0]).tobytes()).hexdigest()[:16]
        finfo = {"k0": k0, "x_wall": float(xw[k0]), "r_wall": float(net.a["r"][k0, n]), "x_axis": float(net.a["x"][k0, 0]),
                 "n_front": int(n), "all_ok_to_axis": bool(np.all(st == mf.ST_OK)),
                 "n_not_ok": int(np.sum(st != mf.ST_OK)), "thr_axis": float(net.thr_axis[k0]),
                 "front_dr_min": float(np.min(np.diff(net.a["r"][k0]))), "front_dr_max": float(np.max(np.diff(net.a["r"][k0]))),
                 "march_end": res["end"], "march_gate_verdict": res["gate"]["verdict"], "front_sha256_16": h,
                 "x_B": info["x_B"], "front_wall_upstream_of_B": bool(xw[k0] < info["x_B"])}
        _CACHE[key] = (net, k0, finfo)
    return _CACHE[key]


def make_line(wall: str, kind: str, N: int | None):
    if kind == "hall":
        return hall(N_HALL)
    net, k0, _ = front(wall)
    if kind == "front":
        return {"x": net.a["x"][k0].copy(), "r": net.a["r"][k0].copy(), "M": net.a["M"][k0].copy(),
                "theta": net.a["th"][k0].copy(), "thr_axis": float(net.thr_axis[k0]), "source": "front_k0_itself"}
    if kind == "resampled":
        return ra.resampled_front_line(net, k0, int(N))
    raise ValueError(kind)


# ================================================================================================ 当てはめ
def fit_cubic(eps, F) -> dict:
    """F(ε) を u = ε/ε₀ の 3 次で最小二乗 (7 点)。c_i = a_i/ε₀^i (run_rao_t4.py の `fit_cubic` と同じ式)。中心差分も振幅ごとに。"""
    if any(f is None for f in F):
        return {"ok": False, "why": "推力が取れない点がある"}
    e = np.asarray(eps, dtype=float)
    f = np.asarray(F, dtype=float)
    u = e / EPS0
    A = np.vander(u, 4, increasing=True)
    a, *_ = np.linalg.lstsq(A, f, rcond=None)
    res = A @ a - f
    c = [a[0], a[1] / EPS0, a[2] / EPS0 ** 2, a[3] / EPS0 ** 3]
    F0 = float(f[np.flatnonzero(e == 0.0)[0]])
    cd = {}
    for amp in AMPS:
        fp, fm = float(f[np.flatnonzero(np.isclose(e, amp))[0]]), float(f[np.flatnonzero(np.isclose(e, -amp))[0]])
        cd[f"{amp:g}"] = {"c1": (fp - fm) / (2 * amp), "c2": (fp + fm - 2 * F0) / (2 * amp * amp)}
    t1, t2 = abs(c[1] * EPS0), abs(c[2] * EPS0 ** 2)
    eps_star = abs(c[1] / (2.0 * c[2])) if c[2] != 0 else None
    hung = c[1] ** 2 / (4.0 * abs(c[2])) if c[2] != 0 else None
    return {"ok": True, "c0": float(c[0]), "c1": float(c[1]), "c2": float(c[2]), "c3": float(c[3]),
            "c1_eps0": float(c[1] * EPS0), "c2_eps0sq": float(c[2] * EPS0 ** 2), "c3_eps0cu": float(c[3] * EPS0 ** 3),
            "ratio_c1_c2": (t1 / t2) if t2 > 0 else None, "c2_sign": int(np.sign(c[2])),
            "eps_star": eps_star, "hung_thrust": hung, "hung_over_F1D": (hung / F1D_abs()) if hung is not None else None,
            "fit_rms": float(np.sqrt(np.mean(res ** 2))), "fit_max": float(np.max(np.abs(res))), "central": cd, "F0": F0}


# ================================================================================================ ジョブ
def t4_job(spec: dict) -> dict:
    t0 = time.time()
    ig = gas()
    w, info = load_wall(spec["wall"])
    L = info["L_wall"]
    line = make_line(spec["wall"], spec["line"], spec.get("N"))
    lw = mf._Wall(ra.ExtendedWall(w, L), L)
    codes, reasons, linfo = mf._check_line(line, lw, ig.moc)
    if codes:
        return {**spec, "status": "error", "why": f"初期線の検査: {codes} {reasons}", "info": info, "line_check": linfo}
    base = ra.forward_thrust(w, L, ig, line, splits=(w.x_B,), gate=True)
    if not base.get("ok"):
        return {**spec, "status": "error", "why": f"基準の壁の推力: {base.get('why')}", "info": info}
    fam = {}
    for k in spec["ks"]:
        Fs, Ft, recs = [], [], []
        for e in EPSS:
            r = base if e == 0.0 else ra.forward_thrust(ra.PerturbedWall(w, FAMILY, k, e), L, ig, line, splits=(w.x_B,), gate=False)
            Fs.append(r.get("F") if r.get("ok") else None)
            Ft.append(r.get("F_trap") if r.get("ok") else None)
            recs.append({"eps": e, "ok": bool(r.get("ok")), "why": r.get("why"), "F": r.get("F"), "F_trap": r.get("F_trap"),
                         "T_wall": r.get("T_wall"), "dod_clean": r.get("dod_clean"),
                         "dod_n_bad": (r.get("dod") or {}).get("n_bad_in_dod"), "end_reason": r.get("end_reason"),
                         "n_nodes_le_L": r.get("n_nodes_le_L"), "L_eval": r.get("L_eval")})
        fam[f"{FAMILY}_k{k}"] = {"family": FAMILY, "k": k, "eps": list(EPSS), "F": Fs, "F_trap": Ft, "runs": recs,
                                 "fit": fit_cubic(EPSS, Fs), "fit_trap": fit_cubic(EPSS, Ft)}
    return {**spec, "status": "ok", "info": info, "time_s": time.time() - t0,
            "line_info": {k: line.get(k) for k in ("source", "k0", "N", "n_front", "front_dr_min", "front_dr_max", "x_wall",
                                                    "r_wall", "x_axis", "thr_axis") if k in line} | {"n_points": int(len(line["r"])),
                                                                                                      "x_wall_line": float(line["x"][-1])},
            "line_check": linfo, "front_info": (front(spec["wall"])[2] if spec["line"] in ("front", "resampled") else None),
            "base": {"F": base["F"], "F_trap": base["F_trap"], "T_line": base["T_line"], "T_wall": base["T_wall"],
                     "T_wall_trap": base["T_wall_trap"], "T_wall_seg": base["T_wall_seg"], "mdot_line": base["mdot_line"],
                     "gate": base["gate"], "dod": base["dod"], "dod_clean": base["dod_clean"],
                     "nodes_max_abs_diff": base.get("nodes_max_abs_diff"), "mdot_err_max": base.get("mdot_err_max"),
                     "L_eval": base["L_eval"], "n_nodes_le_L": base["n_nodes_le_L"], "x_last_node_le_L": base["x_last_node_le_L"],
                     "n_fronts_eval": base.get("n_fronts_eval")},
            "families": fam}


def gen_job(spec: dict) -> dict:
    """follow 段: 生成器を t6 の mq2 と同じ設定 (Hall n 321・m = q = 320・dr 3.2/321) で回し直し、ctrl_plus・ctrl_minus の壁を
    run_rao_t4_gen.py の `_control` で作る。"""
    import run_rao_t4_gen as g4
    t0 = time.time()
    ig = gas()
    tag = WALLS[spec["base"]]
    t6 = json.loads((T6 / f"{tag}.json").read_text())
    z6 = np.load(T6 / f"{tag}_wall.npz")
    L = float(t6["L_rt"])
    res = ra.design_rao(ig, R_t=R_T, L=L, eps=EPS, kernel="arc", line=hall(N_HALL), p_a=P_A, m=int(t6["m"]), q=int(t6["q"]),
                        dr=float(t6["dr"]))
    rec = {"tag": f"gen_{spec['base']}", "base": spec["base"], "L_rt": L, "n": res["n"], "m": res["m"], "q": res["q"],
           "dr": res["dr"], "status": res["status"], "why": res.get("why")}
    arrays = {}
    if res["status"] == "ok":
        same = res["BE"].shape == z6["BE"].shape
        rec.update({"theta_B": float(res["theta_B"]), "theta_B_deg": float(np.rad2deg(res["theta_B"])),
                    "theta_E_deg": float(np.rad2deg(res["theta_E"])), "x_E": float(res["E"]["x"]), "r_E": float(res["E"]["r"]),
                    "F_cs": res["thrust"]["F_cs"],
                    "vs_t6": {"theta_B_deg": float(np.rad2deg(res["theta_B"])) - t6["theta_B_deg"],
                              "theta_E_deg": float(np.rad2deg(res["theta_E"])) - t6["theta_E_deg"],
                              "x_E": float(res["E"]["x"]) - t6["x_E"], "r_E": float(res["E"]["r"]) - t6["r_E"],
                              "F_cs": res["thrust"]["F_cs"] - t6["F_cs"], "BE_same_shape": bool(same),
                              "BE_max_abs_diff": float(np.max(np.abs(res["BE"] - z6["BE"]))) if same else None}})
        arrays["BE"] = res["BE"]
        arrays["E_xrthM"] = np.array([res["E"]["x"], res["E"]["r"], res["E"]["th"], res["E"]["M"]])
        ctrl = {}
        for name, fac in (("ctrl_plus", 1.0 + CTRL_SHIFT), ("ctrl_minus", 1.0 - CTRL_SHIFT)):
            t1 = time.time()
            try:
                c = g4._control(ig, res, fac)
            except Exception as e:  # noqa: BLE001
                c = {"ok": False, "why": f"{type(e).__name__}: {e}", "traceback": traceback.format_exc()}
            c["time_s"] = time.time() - t1
            if c.get("ok"):
                arrays[f"{name}_BE"] = c.pop("_BE")
                arrays[f"{name}_E_xrthM"] = np.array([c["x_E"], c["r_E"], c["th_E"], c["M_E"]])
            ctrl[name] = c
        rec["controls"] = ctrl
    rec["time_s"] = time.time() - t0
    PARTS.mkdir(exist_ok=True)
    if arrays:
        np.savez(PARTS / f"gen_{spec['base']}_walls.npz", **arrays)
    (PARTS / f"gen_{spec['base']}_rec.json").write_text(json.dumps(_jsonable(rec), ensure_ascii=False, indent=1))
    return {**spec, "status": rec["status"], "why": rec.get("why"), "gen": rec, "time_s": rec["time_s"]}


def jobs(follow: bool) -> list:
    J = []
    for N in N_LIST:
        J.append({"name": f"ab_L100_resampled_N{N}", "kind": "t4", "phase": "ab", "wall": "L100", "line": "resampled", "N": N,
                  "ks": [1]})
    J.append({"name": "ref_L100_front", "kind": "t4", "phase": "ref", "wall": "L100", "line": "front", "N": None, "ks": [1]})
    J.append({"name": "ref_L100_hall321", "kind": "t4", "phase": "ref", "wall": "L100", "line": "hall", "N": None, "ks": [1]})
    if follow:
        J.append({"name": "follow_L100_resampled_N1281_k23", "kind": "t4", "phase": "follow", "wall": "L100", "line": "resampled",
                  "N": 1281, "ks": [2, 3]})
        J.append({"name": "follow_L70_resampled_N1281", "kind": "t4", "phase": "follow", "wall": "L70", "line": "resampled",
                  "N": 1281, "ks": [1, 2, 3]})
        for b in ("L100", "L70"):
            J.append({"name": f"gen_{b}", "kind": "gen", "phase": "follow", "base": b})
            J.append({"name": f"follow_{b}_ctrl_plus_resampled_N1281", "kind": "t4", "phase": "follow", "wall": f"{b}_ctrl_plus",
                      "line": "resampled", "N": 1281, "ks": [1, 2, 3], "needs": f"gen_{b}"})
    return J


def run(spec: dict) -> dict:
    try:
        out = t4_job(spec) if spec["kind"] == "t4" else gen_job(spec)
    except Exception as e:  # noqa: BLE001
        out = {**spec, "status": "error", "why": f"{type(e).__name__}: {e}", "traceback": traceback.format_exc()}
    PARTS.mkdir(exist_ok=True)
    (PARTS / f"{spec['name']}.json").write_text(json.dumps(_jsonable(out), ensure_ascii=False, indent=1))
    print(summary_line(out), flush=True)
    return out


def summary_line(o: dict) -> str:
    if o.get("status") != "ok":
        return f"{o['name']}: {o.get('status')} — {o.get('why')}"
    if o["kind"] == "gen":
        g = o["gen"]
        return f"{o['name']}: ok vs t6 {g.get('vs_t6')} ctrl {[(k, v.get('ok')) for k, v in (g.get('controls') or {}).items()]} ({o['time_s']:.0f} s)"
    s = (f"{o['name']}: F0 {o['base']['F']:.10f} (台形 {o['base']['F_trap']:.10f}) 壁点 {o['base']['n_nodes_le_L']} "
         f"依存域 clean {o['base']['dod_clean']} ゲート {o['base']['gate']['verdict']} ({o['time_s']:.0f} s)")
    for key, v in o["families"].items():
        for acc, f in (("spline", v["fit"]), ("trap", v["fit_trap"])):
            s += (f"\n   {key} {acc}: " + (f"c1ε0 {f['c1_eps0']:+.3e} c2ε0² {f['c2_eps0sq']:+.3e} ε* {f['eps_star']:.2e} "
                                           f"当てはめ rms {f['fit_rms']:.1e}" if f.get("ok") else f"不可 {f.get('why')}"))
    return s


# ================================================================================================ 再現の検査
def reproduction() -> dict:
    """最初に確かめる項目: 前線 k₀ の状態、各 N の線の `_check_line`、往復の誤差、前線そのものからのマーチの再現。"""
    ig = gas()
    out = {}
    for wall in ("L100",):
        net, k0, finfo = front(wall)
        w, info = load_wall(wall)
        L = info["L_wall"]
        lw = mf._Wall(ra.ExtendedWall(w, L), L)
        fr = {f: net.a[f][k0] for f in net.FIELDS}
        n = finfo["n_front"]
        band = (k0 - np.arange(n + 1) <= 12) & (k0 - np.arange(n + 1) >= -6)
        per_N = {}
        for N in N_LIST:
            ln = ra.resampled_front_line(net, k0, N)
            codes, reasons, linfo = mf._check_line(ln, lw, ig.moc)
            rt = {}
            for nm, key in (("x", "x"), ("th", "theta"), ("M", "M")):
                e = np.abs(CubicSpline(ln["r"], ln[key])(fr["r"]) - fr[nm])
                j = int(np.argmax(e))
                rt[nm] = {"max": float(e.max()), "j": j, "k0_minus_j": int(k0 - j), "r": float(fr["r"][j]),
                          "max_outside_band": float(e[~band].max())}
            per_N[N] = {"check_line_codes": codes, "check_line_reasons": reasons, "check_line_info": linfo, "round_trip": rt}
        # 前線そのものを初期線にしたマーチが元の網の前線 k₀+1.. を再現するか
        lf = make_line(wall, "front", None)
        r2 = mf.forward_moc(ra.ExtendedWall(w, L), L, ig.moc, lf)
        K = 30
        d = max(float(np.nanmax(np.abs(net.a[f][k0:k0 + K] - r2["net"].a[f][:K]))) for f in net.FIELDS)
        # 前線の鋸歯の帯 (k₀ − j ∈ [−6, 12]) の振幅: 帯の外の点を通る 3 次スプラインからのずれ
        saw = {}
        for nm in ("th", "M"):
            dev = fr[nm] - CubicSpline(fr["r"][~band], fr[nm][~band])(fr["r"])
            saw[nm] = {"max_abs_dev_in_band": float(np.abs(dev[band]).max()),
                       "k0_minus_j_at_max": int(k0 - np.flatnonzero(band)[int(np.argmax(np.abs(dev[band])))])}
        out[wall] = {"front": finfo, "per_N": per_N, "front_as_line_march_max_diff_30_fronts": d,
                     "front_as_line_status_same": bool(np.array_equal(net.status[k0:k0 + K], r2["net"].status[:K])),
                     "sawtooth_band": {"k0_minus_j": [-6, 12], "r_range": [float(fr["r"][band].min()), float(fr["r"][band].max())],
                                       "dev_from_smooth": saw},
                     "thr_axis_first_fronts": [float(v) for v in net.thr_axis[1:13]]}
    return out


# ================================================================================================ 集計
def _get(parts, name):
    o = parts.get(name)
    return o if (o is not None and o.get("status") == "ok") else None


def judgement(parts: dict) -> dict:
    F1 = F1D_abs()
    c = {"spline": {}, "trap": {}}
    c2 = {"spline": {}, "trap": {}}
    F0 = {"spline": {}, "trap": {}}
    es, hg = {"spline": {}, "trap": {}}, {"spline": {}, "trap": {}}
    for N in N_LIST:
        o = _get(parts, f"ab_L100_resampled_N{N}")
        if o is None:
            return {"why": f"ab_L100_resampled_N{N} が無い / ok でない"}
        v = o["families"][f"{FAMILY}_k1"]
        for acc, f in (("spline", v["fit"]), ("trap", v["fit_trap"])):
            if not f.get("ok"):
                return {"why": f"N {N} の {acc} の当てはめが不可"}
            c[acc][N], c2[acc][N], F0[acc][N] = f["c1_eps0"], f["c2_eps0sq"], f["F0"]
            es[acc][N], hg[acc][N] = f["eps_star"], f["hung_over_F1D"]
    diff = {N: abs(c["spline"][N] - c["trap"][N]) for N in N_LIST}
    conv = {}
    for acc in ("spline", "trap"):
        d1, d2 = F0[acc][321] - F0[acc][641], F0[acc][641] - F0[acc][1281]      # h の比は 641/321・1281/641 (≈ 2)
        conv[acc] = {"F0": F0[acc], "F0_over_F1D": {N: F0[acc][N] / F1 for N in N_LIST},
                     "d_321_641": d1, "d_641_1281": d2, "d_over_F1D": [d1 / F1, d2 / F1],
                     "observed_order": (float(np.log(abs(d1 / d2)) / np.log(1281.0 / 641.0)) if d2 != 0 else None)}
    cond = {}
    for acc in ("spline", "trap"):
        a, b = c[acc][641], c[acc][1281]
        A1 = abs(b) <= LINE_ABS
        A2 = diff[1281] <= LINE_ABS
        A3 = abs(b) < abs(a)
        rel = (b / a - 1.0) if a != 0 else None
        Bc = bool(A2 and abs(a) > LINE_ABS and abs(b) > LINE_ABS and rel is not None and abs(rel) <= B_REL)
        cond[acc] = {"A1_abs_c1eps0_1281_le_1e-7": bool(A1), "A2_acc_diff_1281_le_1e-7": bool(A2),
                     "A3_abs_c1eps0_decreases_641_to_1281": bool(A3), "A_all": bool(A1 and A2 and A3),
                     "B_rel_641_1281": rel, "B_all": Bc, "C_acc_diff_1281_gt_1e-7": bool(diff[1281] > LINE_ABS)}
    return {"F1D_abs": F1, "c1_eps0": c, "c2_eps0sq": c2, "acc_diff_c1_eps0": diff, "eps_star_rt": es, "hung_over_F1D": hg,
            "F0_convergence": conv, "conditions": cond,
            "note": "判定 (A/B/C) は親と上位。登録文は |c₁ε₀| の会計を指定していないので条件は会計ごと。"
                    "B の 30 % は |c₁ε₀(1281)/c₁ε₀(641) − 1| ≤ 0.3。"}


ROW_KEYS = ("c1_eps0", "c2_eps0sq", "c3_eps0cu", "ratio_c1_c2", "eps_star", "hung_over_F1D", "fit_rms", "fit_max")


def assemble(total_s=None) -> dict:
    parts = {}
    if PARTS.exists():
        for p in sorted(PARTS.glob("*.json")):
            parts[p.stem] = json.loads(p.read_text())
    rows = []
    for name, o in sorted(parts.items()):
        if o.get("kind") != "t4":
            continue
        if o.get("status") != "ok":
            rows.append({"job": name, "status": o.get("status"), "why": o.get("why")})
            continue
        b, info, li = o["base"], o["info"], o["line_info"]
        for key, v in o["families"].items():
            row = {"job": name, "status": "ok", "phase": o["phase"], "wall": o["wall"], "line": o["line"], "N": o.get("N"),
                   "n_line_points": li["n_points"], "x_line_wall": li["x_wall_line"], "family": v["family"], "k": v["k"],
                   "L_wall": info["L_wall"], "F0": b["F"], "F0_trap": b["F_trap"], "T_line": b["T_line"],
                   "n_nodes_le_L": b["n_nodes_le_L"], "base_gate": b["gate"]["verdict"], "base_dod_clean": b["dod_clean"],
                   "n_pert_dod_clean": sum(1 for r in v["runs"] if r.get("dod_clean")), "n_pert": len(v["runs"]),
                   "time_s": o["time_s"]}
            for acc, f in (("", v["fit"]), ("_trap", v["fit_trap"])):
                if f.get("ok"):
                    row.update({f"{k}{acc}": f[k] for k in ROW_KEYS})
                    row[f"c2_sign{acc}"] = f["c2_sign"]
            if v["fit"].get("ok") and v["fit_trap"].get("ok"):
                row["acc_diff_c1_eps0"] = abs(v["fit"]["c1_eps0"] - v["fit_trap"]["c1_eps0"])
            rows.append(row)
    rep = parts.get("_reproduction")
    out = {"plan": "plans/active/tooling-nozzle-bell-two-mode.md §5.1 #5 (h)・§6 V1 T4",
           "record": "notes/reviews/2026-10-12-bell-rao-t4-stationarity-diagnostician.md",
           "method": {"wall": "rao_t6_parts/L100_n321_mq2 の BE → RaoWall pspline (5 次・48 区間・λ GCV、B で C1)、末尾は E",
                      "line": f"Hall n {N_HALL} の前進 MOC の前線 k₀ (壁点 x ≈ {X_FRONT}) を r 等間隔 N+1 点へ 3 次補間 "
                              "(rao_axisym.resampled_front_line)",
                      "thrust": "rao_axisym.forward_thrust (splits = (x_B,))、会計 = スプライン (F) と台形 (F_trap)",
                      "family": FAMILY, "eps": EPSS, "eps0": EPS0, "N": N_LIST, "fit": "u = ε/ε₀ の 3 次の最小二乗 (7 点)",
                      "F1D_abs": F1D_abs(),
                      "registered": {"A": "N 1281 で |c₁ε₀| ≤ 1e-7 かつ両会計の c₁ε₀ の差 ≤ 1e-7 かつ 641 → 1281 で |c₁ε₀| が減る",
                                     "B": "両会計の差 ≤ 1e-7 だが |c₁ε₀| が 641 と 1281 で 30 % 以内に一致し両方 > 1e-7",
                                     "C": "N 1281 でも 2 会計の差 > 1e-7"},
                      "references": "ref_L100_front (前線 k₀ そのもの、再標本化なし)・ref_L100_hall321 (Hall n 321、前任の評価器) — 判定に使わない"},
           **git_info(), "reproduction": rep, "judgement": judgement(parts), "rows": rows,
           "parts": {k: v for k, v in parts.items() if k != "_reproduction"}, "total_s": total_s}
    (HERE / "rao_t4h.json").write_text(json.dumps(_jsonable(out), ensure_ascii=False, indent=1))
    write_csv(HERE / "rao_t4h.csv", rows)
    return out


def write_csv(path: Path, rows: list) -> None:
    cols = []
    for r in rows:
        for k in r:
            if k not in cols and not isinstance(r[k], (dict, list)):
                cols.append(k)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: ("" if r.get(k) is None else (f"{r[k]:.10g}" if isinstance(r[k], float) else r[k])) for k in cols})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default=None)
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--follow", action="store_true", help="判定 A のときだけ: k2・k3・L70・ctrl_plus を N 1281 で")
    ap.add_argument("--assemble-only", action="store_true")
    a = ap.parse_args()
    if a.workers > 2:
        raise SystemExit("--workers は 2 まで (手元 PC はユーザの作業機)")
    t0 = time.time()
    PARTS.mkdir(exist_ok=True)
    if not a.assemble_only:
        rp = PARTS / "_reproduction.json"
        if not rp.exists():
            rep = reproduction()
            rp.write_text(json.dumps(_jsonable(rep), ensure_ascii=False, indent=1))
            f = rep["L100"]["front"]
            print(f"再現の検査: 前線 k₀ {f['k0']} (壁点 x {f['x_wall']:.6f}・軸点 x {f['x_axis']:.6f}) 軸まで OK {f['all_ok_to_axis']}、"
                  f"_check_line {[(N, v['check_line_codes']) for N, v in rep['L100']['per_N'].items()]}、"
                  f"往復 θ {[(N, round(v['round_trip']['th']['max'], 12)) for N, v in rep['L100']['per_N'].items()]}、"
                  f"前線そのもののマーチの差 {rep['L100']['front_as_line_march_max_diff_30_fronts']:.1e}", flush=True)
            if not f["all_ok_to_axis"] or any(v["check_line_codes"] for v in rep["L100"]["per_N"].values()):
                print("止める: 前線 k₀ が軸まで OK でない / 再標本化した線が _check_line を通らない", flush=True)
                assemble(time.time() - t0)
                return
        todo = [j for j in jobs(a.follow) if (a.only is None or any(s in j["name"] for s in a.only.split(",")))
                and not (PARTS / f"{j['name']}.json").exists()]
        print(f"回すジョブ {[j['name'] for j in todo]} (並列 {a.workers})", flush=True)
        with ProcessPoolExecutor(max_workers=a.workers) as ex:
            pend = {}
            waiting = list(todo)

            def submit_ready():
                for j in list(waiting):
                    need = j.get("needs")
                    if need is None or (PARTS / f"{need}.json").exists():
                        pend[ex.submit(run, j)] = j
                        waiting.remove(j)
            submit_ready()
            while pend:
                done, _ = wait(list(pend), return_when=FIRST_COMPLETED)
                for fu in done:
                    pend.pop(fu)
                    fu.result()
                submit_ready()
                assemble(time.time() - t0)
            for j in waiting:
                print(f"{j['name']}: 回さない (前提 {j.get('needs')} が無い)", flush=True)
    out = assemble(time.time() - t0)
    print(json.dumps(_jsonable(out["judgement"]), ensure_ascii=False, indent=1), flush=True)
    print(f"rao_t4h.csv {len(out['rows'])} 行 ({time.time() - t0:.0f} s)", flush=True)


if __name__ == "__main__":
    main()
