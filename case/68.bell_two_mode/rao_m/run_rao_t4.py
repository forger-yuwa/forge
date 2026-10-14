#!/usr/bin/env python3
r"""③ベル (M) モード §5.1 #5 (g2) T4 停留性・(g3) T1a′ — 与えた壁で前進 MOC を解き直した推力の表 (数値を出すだけ、
合否の判定は親と上位)。

plan: plans/active/tooling-nozzle-bell-two-mode.md §5.1 #5 (g)・§6 V1 の T1a′・T4・「TP の停留性」。
(M) の壁と負対照は `run_rao_t4_gen.py` の部品 (`rao_t4_parts/gen_*`) から読む (生成器の既定 m = q = n)。

**推力の求め方 (M・D 共通)**: `rao_axisym.forward_thrust` — 壁を与えた前進 MOC (`moc_forward.forward_moc`、analytic・converge、
初期線 = Hall の n_eval 点、ガス = 基準ケースの凍結 TP の `Isentrope`) で場を解き、F = 初期線の運動量 + 圧力 (`line_integrals`) +
壁の圧力の軸力 ∫(p − p_a) 2π r r' dx (壁点の p の区間ごとの 3 次スプライン × 壁の式の r・r'、Gauss–Legendre)。(M) は B で区間を
分ける (B で曲率が跳ぶ。摂動 `sin` では勾配も跳ぶ)。推力に使った壁点の依存域に網の反転・失敗が無いことを `dod_clean` に記録。
(M) の BE は生成器の壁点列の罰則付き 5 次 B-spline (`RaoWall` の `"pspline"`、λ は GCV) — 点列をそのまま補間すると B の近くで
曲率が ±0.1 振動し、前進 MOC が BDE の中で C⁻ 前線の交差を出す (2026-10-12 実測、Goursat 網には反転なし)。

**T4** (§6 V1): 基準の壁 = (M) の壁 (円弧 kernel・EXH・R_t 1・ε 50、L100・L70 × n 161、L100 × n 241)。壁 BE に
η_k(t) (t = B → E の弧長の正規化) を半径方向に ε ∈ {±1, ±2, ±4}×10⁻³ r_t で足す (B・E は動かない → L・ε・kernel・初期線 [流量]
は同じ)。族は `sin` = sin(kπt) (B・E で勾配が変わる) と `sin_c1` = sin(πt) sin(kπt) (B・E の接線連続を保つ)、k = 1, 2, 3。
F(ε) の 7 点を u = ε/ε₀ (ε₀ = 4e-3) の 3 次で最小二乗し c₁ = a₁/ε₀・c₂ = a₂/ε₀²。登録の線: |c₁ε₀| ≤ 0.05 |c₂ε₀²|、
|c₁(M)| ≤ 0.05 |c₁(負対照)|、c₂ の符号 (負なら局所最大)。n_eval = 161・241・321 (評価器の分解能) × 生成器の n (161・241) で
c₁ の n 依存を見る。
**負対照** (§6 V1 T4「D を 2 % ずらした非 Rao 線から作った壁」): 根の BD 上で σ_D (B → D の弧長) を 1.02 倍した D⁺ から
Φ₂・Φ₃ 一定の線を引き直した壁 (x_E' > L) を x = L で切ったもの (`ctrl_plus`)。切り口の制御面 (C⁺) は Φ₂・Φ₃ 一定でないので
端を固定した摂動に停留でない。参考に 0.98 倍の D⁻ の壁を切らずに使ったもの (`ctrl_minus`、その (L', ε') の Rao 輪郭 = 停留の
見込み) も n_eval = n_gen で回す。
**表現の感度** (L100 n 161・n_eval 161): BE の当てはめの λ を GCV の 1/10・10 倍にした基準の壁で同じ T4。

**T1a′** (§6 V1): 同じ L_k (100 %・70 %)・ε で、前進 MOC ゲートを通った (D) 壁の推力 ≤ (M) の推力。(D) は `design_bell_c2`
(θ_max {20, 25, 32.5, 40}° × θ_e {5, 10, 15}° の V0 の格子のうち成立する点 + (M) の (θ_max, θ_e) の近くの 4 点)。ゲートは同じ
前進 MOC (n_eval、L まで) の `gate_record`。V0-c の n 321 の判定 (`../v0_moc_gate/v0c_gate.csv`、ガスは `moc_gas` そのもの) も
並べる。(M) の推力は T4 の ε = 0 と同じ評価。

出力 (このディレクトリ): `rao_t4_eval_parts/` (ジョブごと、再開用)、`rao_t4.csv` (T4: 基準の壁 × n_eval × 族 × k の 1 行)、
`rao_t4.json` (全 F(ε)・当てはめ・ゲート・依存域・基準の壁の推力と生成器の F_cs の差・来歴)、`rao_t1ap.csv`・`rao_t1ap.json`。
実行: `nice -n 10 /home/sano/work/forge/design/.venv-opt/bin/python case/68.bell_two_mode/rao_m/run_rao_t4.py` (2 並列)。
`--only` (ジョブ名の部分一致、カンマ区切り)、`--assemble-only`。
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
from concurrent.futures import ProcessPoolExecutor, as_completed  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "design"))
from forge_design.gas.rocket_exhaust import RocketExhaust  # noqa: E402
from forge_design.geometry import rao_axisym as ra  # noqa: E402
from forge_design.geometry.wall_bell_c2 import BellC2Spec, cone15_length, design_bell_c2  # noqa: E402

EXHAUST = HERE.parent / "exhaust_lox_lh2_pc100_of6.yaml"
GEN = HERE / "rao_t4_parts"
PARTS = HERE / "rao_t4_eval_parts"
V0C = HERE.parent / "v0_moc_gate" / "v0c_gate.csv"
EPS = 50.0
R_T = 1.0
P_A = 0.0
EPS0 = 4e-3
AMPS = (1e-3, 2e-3, 4e-3)
EPSS = tuple(sorted([-a for a in AMPS] + [0.0] + list(AMPS)))
FAMILIES = ("sin", "sin_c1")
KS = (1, 2, 3)
N_EVALS = (161, 241, 321)
RATIO_LINE = 0.05
SOURCES = ("design/forge_design/geometry/rao_axisym.py", "design/forge_design/geometry/moc_forward.py",
           "design/forge_design/geometry/moc_kernel.py", "design/forge_design/gas/isentrope.py",
           "design/forge_design/geometry/wall_bell_c2.py")
D_GRID = {"theta_max": (20.0, 25.0, 32.5, 40.0), "theta_e": (5.0, 10.0, 15.0)}


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
    key = ("line", n)
    if key not in _CACHE:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            _CACHE[key] = ra._hall_line(R_T, gas(), n)
    return _CACHE[key]


def F1D_abs() -> float:
    """1D 理想推力 (同じ凍結 TP で ε まで、A* = π r_t²)。相対量の分母。"""
    return float(np.pi * gas().ideal_1d(EPS, P_A)["F"])


def load_gen(tag: str):
    j = GEN / f"{tag}.json"
    z = GEN / f"{tag}_walls.npz"
    if not (j.exists() and z.exists()):
        return None
    rec = json.loads(j.read_text())
    if rec.get("status") != "ok":
        return None
    return rec, dict(np.load(z))


def base_wall(gen_tag: str, which: str = "M", lam_scale: float | None = None):
    """(M) の壁 (`which` = "M")、負対照 `ctrl_plus` (x = L で切る)・`ctrl_minus` (切らない)。戻り (RaoWall, 情報 dict)。"""
    g = load_gen(gen_tag)
    if g is None:
        raise FileNotFoundError(f"{gen_tag} の部品が無い / ok でない")
    rec, arr = g
    L = float(rec["L_rt"])
    if which == "M":
        BE, E = arr["BE"], arr["E_xrthM"]
    else:
        if f"{which}_BE" not in arr:
            raise FileNotFoundError(f"{gen_tag} の {which} が無い ({(rec.get('controls') or {}).get(which, {}).get('why')})")
        BE, E = arr[f"{which}_BE"], arr[f"{which}_E_xrthM"]
    P = BE[:, :2].copy()
    last_gap = (float(P[-1, 0] - E[0]), float(P[-1, 1] - E[1]))
    P[-1] = E[:2]
    x_end = L if (which == "ctrl_plus" and E[0] > L) else None
    lam = None
    if lam_scale is not None:
        w0 = ra.RaoWall(R_T, rec["theta_B"], P, float(E[2]), x_end=x_end)
        lam = w0.fit_info["lam"] * lam_scale
    w = ra.RaoWall(R_T, rec["theta_B"], P, float(E[2]), x_end=x_end, lam=lam)
    L_w = float(w.x_end)
    info = {"gen_tag": gen_tag, "which": which, "L_wall": L_w, "eps_wall": float(w.r_end) ** 2, "theta_end_deg": float(np.rad2deg(w.theta_end)),
            "x_B": w.x_B, "theta_B_deg": float(np.rad2deg(rec["theta_B"])), "B_replaced": w.B_replaced,
            "last_point_minus_E": last_gap, "fit_info": w.fit_info, "n_gen": rec["n"], "m_gen": rec["m"],
            "F_cs_gen": rec.get("F_cs") if which == "M" else None, "F_wall_gen": rec.get("F_wall") if which == "M" else None,
            "lam_scale": lam_scale}
    return w, info


# ================================================================================================ T4
def t4_job(spec: dict) -> dict:
    t0 = time.time()
    ig = gas()
    line = hall(spec["n_eval"])
    w, info = base_wall(spec["gen"], spec["which"], spec.get("lam_scale"))
    L = info["L_wall"]
    base = ra.forward_thrust(w, L, ig, line, splits=(w.x_B,), gate=True)
    if not base.get("ok"):
        return {**spec, "status": "error", "why": f"基準の壁の推力: {base.get('why')}", "info": info}
    vals = {}
    for fam in FAMILIES:
        for k in KS:
            Fs = []
            recs = []
            for e in EPSS:
                if e == 0.0:
                    r = base
                else:
                    r = ra.forward_thrust(ra.PerturbedWall(w, fam, k, e), L, ig, line, splits=(w.x_B,), gate=False)
                Fs.append(r.get("F") if r.get("ok") else None)
                recs.append({"eps": e, "ok": bool(r.get("ok")), "why": r.get("why"), "F": r.get("F"),
                             "dod_clean": r.get("dod_clean"), "dod_n_bad": (r.get("dod") or {}).get("n_bad_in_dod"),
                             "end_reason": r.get("end_reason"), "F_trap": r.get("F_trap")})
            vals[f"{fam}_k{k}"] = {"family": fam, "k": k, "eps": list(EPSS), "F": Fs, "runs": recs, "fit": fit_cubic(EPSS, Fs)}
    out = {**spec, "status": "ok", "info": info, "time_s": time.time() - t0,
           "base": {"F": base["F"], "F_trap": base["F_trap"], "T_line": base["T_line"], "T_wall": base["T_wall"],
                    "T_wall_seg": base["T_wall_seg"], "mdot_line": base["mdot_line"], "gate": base["gate"],
                    "dod": base["dod"], "dod_clean": base["dod_clean"], "nodes_max_abs_diff": base.get("nodes_max_abs_diff"),
                    "mdot_err_max": base.get("mdot_err_max"), "L_eval": base["L_eval"], "n_nodes_le_L": base["n_nodes_le_L"]},
           "families": vals}
    return out


def fit_cubic(eps, F) -> dict:
    """F(ε) を u = ε/ε₀ の 3 次で最小二乗 (7 点)。c_i = a_i/ε₀^i。中心差分の c₁・c₂ も振幅ごとに。"""
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
    t1 = abs(c[1] * EPS0)
    t2 = abs(c[2] * EPS0 ** 2)
    return {"ok": True, "c0": float(c[0]), "c1": float(c[1]), "c2": float(c[2]), "c3": float(c[3]),
            "c1_eps0": float(c[1] * EPS0), "c2_eps0sq": float(c[2] * EPS0 ** 2), "c3_eps0cu": float(c[3] * EPS0 ** 3),
            "ratio_c1_c2": (t1 / t2) if t2 > 0 else None, "c2_sign": int(np.sign(c[2])),
            "fit_rms": float(np.sqrt(np.mean(res ** 2))), "fit_max": float(np.max(np.abs(res))), "central": cd, "F0": F0}


# ================================================================================================ T1a′
def d_points(L_frac: float, gm: dict | None) -> list:
    pts = [(tm, te, "grid") for tm in D_GRID["theta_max"] for te in D_GRID["theta_e"]]
    if gm is not None:
        tmM, teM = gm["theta_max_deg"], gm["theta_E_deg"]
        for tm, te in ((tmM, teM), (tmM, teM + 2.0), (tmM - 2.5, teM), (tmM - 2.5, teM + 2.0)):
            pts.append((round(tm, 4), round(te, 4), "near_M"))
    return pts


def t1ap_job(spec: dict) -> dict:
    t0 = time.time()
    ig = gas()
    n = spec["n_eval"]
    line = hall(n)
    L = spec["L_frac"] * cone15_length(EPS)
    rows = []
    gm = load_gen(spec["gen"])
    for tm, te, src in d_points(spec["L_frac"], gm[0] if gm else None):
        row = {"L_frac": spec["L_frac"], "L_rt": L, "n_eval": n, "family": "D", "source": src, "theta_max_deg": tm, "theta_e_deg": te}
        try:
            d = design_bell_c2(L, te, tm, BellC2Spec(eps=EPS, R_t=R_T))
        except ValueError as e:
            rows.append({**row, "design_status": "error", "why": str(e)})
            continue
        row.update({"design_status": d.status, "x_m": d.x_m})
        if not d.ok:
            rows.append({**row, "why": d.reason})
            continue
        r = ra.forward_thrust(d.wall, L, ig, line, splits=(), gate=True)
        if not r.get("ok"):
            rows.append({**row, "why": r.get("why")})
            continue
        g = r["gate"]
        rows.append({**row, "F": r["F"], "F_trap": r["F_trap"], "T_wall": r["T_wall"], "gate_verdict": g["verdict"],
                     "gate_pass": g["gate_pass"], "gate_codes": g["codes"], "cminus_sep_ratio_min": g["cminus_sep_ratio_min"],
                     "dod_clean": r["dod_clean"], "dod_n_bad": r["dod"]["n_bad_in_dod"], "mdot_err_max": r.get("mdot_err_max")})
    return {**spec, "status": "ok", "rows": rows, "time_s": time.time() - t0}


# ================================================================================================ ジョブ
def jobs() -> list:
    J = []
    for gen, n_gen in (("gen_L100_n161", 161), ("gen_L70_n161", 161), ("gen_L100_n241", 241)):
        for ne in N_EVALS:
            J.append({"name": f"t4_{gen}_M_ne{ne}", "kind": "t4", "gen": gen, "which": "M", "n_eval": ne})
            J.append({"name": f"t4_{gen}_ctrl_plus_ne{ne}", "kind": "t4", "gen": gen, "which": "ctrl_plus", "n_eval": ne})
        J.append({"name": f"t4_{gen}_ctrl_minus_ne{n_gen}", "kind": "t4", "gen": gen, "which": "ctrl_minus", "n_eval": n_gen})
    for ls in (0.1, 10.0):
        J.append({"name": f"t4_gen_L100_n161_M_ne161_lam{ls:g}", "kind": "t4", "gen": "gen_L100_n161", "which": "M", "n_eval": 161,
                  "lam_scale": ls})
    for fr, gen in ((1.0, "gen_L100_n161"), (0.7, "gen_L70_n161")):
        for ne in N_EVALS:
            J.append({"name": f"t1ap_L{int(round(fr * 100))}_ne{ne}", "kind": "t1ap", "L_frac": fr, "gen": gen, "n_eval": ne})
    return J


def run(spec: dict) -> dict:
    try:
        out = t4_job(spec) if spec["kind"] == "t4" else t1ap_job(spec)
    except Exception as e:  # noqa: BLE001
        out = {**spec, "status": "error", "why": f"{type(e).__name__}: {e}", "traceback": traceback.format_exc()}
    PARTS.mkdir(exist_ok=True)
    (PARTS / f"{spec['name']}.json").write_text(json.dumps(_jsonable(out), ensure_ascii=False, indent=1))
    print(summary_line(out), flush=True)
    return out


def summary_line(o: dict) -> str:
    if o.get("status") != "ok":
        return f"{o['name']}: {o.get('status')} — {o.get('why')}"
    if o["kind"] == "t1ap":
        ps = [r for r in o["rows"] if r.get("gate_pass")]
        return f"{o['name']}: D {len(o['rows'])} 点 (ゲート通過 {len(ps)}) ({o['time_s']:.0f} s)"
    s = f"{o['name']}: F0 {o['base']['F']:.9f} 依存域 clean {o['base']['dod_clean']} ({o['time_s']:.0f} s)"
    for key, v in o["families"].items():
        f = v["fit"]
        s += (f"\n   {key}: " + (f"c1ε0 {f['c1_eps0']:+.3e} c2ε0² {f['c2_eps0sq']:+.3e} 比 {f['ratio_c1_c2']:.3g} 当てはめ rms {f['fit_rms']:.1e}"
                                 if f.get("ok") else f"不可 {f.get('why')}"))
    return s


# ================================================================================================ 集計
def assemble() -> dict:
    parts = {}
    if PARTS.exists():
        for p in sorted(PARTS.glob("*.json")):
            parts[p.stem] = json.loads(p.read_text())
    F1 = F1D_abs()
    # --- T4 ---
    t4rows = []
    ctrl = {}
    for name, o in parts.items():
        if o.get("kind") == "t4" and o.get("status") == "ok" and o["which"] == "ctrl_plus":
            ctrl[(o["gen"], o["n_eval"])] = o
    bases = {}
    for name, o in sorted(parts.items()):
        if o.get("kind") != "t4":
            continue
        if o.get("status") != "ok":
            t4rows.append({"job": name, "status": o.get("status"), "why": o.get("why")})
            continue
        info, b = o["info"], o["base"]
        if o["which"] == "M" and o.get("lam_scale") is None:
            bases[(o["gen"], o["n_eval"])] = b["F"]
        for key, v in o["families"].items():
            f = v["fit"]
            row = {"job": name, "status": "ok", "gen": o["gen"], "which": o["which"], "lam_scale": o.get("lam_scale"),
                   "n_gen": info["n_gen"], "m_gen": info["m_gen"], "n_eval": o["n_eval"], "L_wall": info["L_wall"],
                   "eps_wall": info["eps_wall"], "family": v["family"], "k": v["k"], "F0": b["F"],
                   "F0_minus_Fcs_gen": (b["F"] - info["F_cs_gen"]) if info.get("F_cs_gen") is not None else None,
                   "F0_minus_Fwall_gen": (b["F"] - info["F_wall_gen"]) if info.get("F_wall_gen") is not None else None,
                   "base_gate": b["gate"]["verdict"], "base_dod_clean": b["dod_clean"],
                   "n_pert_dod_clean": sum(1 for r in v["runs"] if r.get("dod_clean")), "n_pert": len(v["runs"]),
                   "fit_ok": f.get("ok")}
            if f.get("ok"):
                row.update({k: f[k] for k in ("c1", "c2", "c3", "c1_eps0", "c2_eps0sq", "c3_eps0cu", "ratio_c1_c2", "c2_sign",
                                              "fit_rms", "fit_max")})
                row["c1_cd_1e-3"] = f["central"]["0.001"]["c1"]
                row["c1_cd_4e-3"] = f["central"]["0.004"]["c1"]
                row["c2_cd_1e-3"] = f["central"]["0.001"]["c2"]
                row["c2_cd_4e-3"] = f["central"]["0.004"]["c2"]
                row["ratio_le_0.05"] = bool(f["ratio_c1_c2"] is not None and f["ratio_c1_c2"] <= RATIO_LINE)
                if info.get("F_cs_gen") is not None:
                    row["absdiff_F0_Fcs_over_abs_c2eps0sq"] = abs(b["F"] - info["F_cs_gen"]) / abs(f["c2_eps0sq"])
                c = ctrl.get((o["gen"], o["n_eval"]))
                if o["which"] == "M" and c is not None:
                    cf = c["families"][key]["fit"]
                    if cf.get("ok"):
                        row["ctrl_c1"] = cf["c1"]
                        row["c1_over_ctrl_c1"] = abs(f["c1"]) / abs(cf["c1"]) if cf["c1"] != 0 else None
                        row["ctrl_ratio_le_0.05"] = bool(row["c1_over_ctrl_c1"] is not None and row["c1_over_ctrl_c1"] <= RATIO_LINE)
            t4rows.append(row)
    # --- T1a′ ---
    v0c = {}
    if V0C.exists():
        with open(V0C, newline="") as fh:
            for r in csv.DictReader(fh):
                if r["config"] == "hall_n321" and r["family"] == "D":
                    v0c[(round(float(r["L_frac"]), 3), float(r["theta_max_or_n_deg"]), float(r["theta_e_deg"]))] = r["verdict"]
    t1rows = []
    t1sum = {}
    for name, o in sorted(parts.items()):
        if o.get("kind") != "t1ap" or o.get("status") != "ok":
            continue
        FM = bases.get((o["gen"], o["n_eval"]))
        gm = load_gen(o["gen"])
        Mrow = {"L_frac": o["L_frac"], "n_eval": o["n_eval"], "family": "M", "source": o["gen"], "F": FM,
                "theta_max_deg": gm[0]["theta_max_deg"] if gm else None, "theta_e_deg": gm[0]["theta_E_deg"] if gm else None}
        t4M = parts.get(f"t4_{o['gen']}_M_ne{o['n_eval']}")
        if t4M and t4M.get("status") == "ok":
            Mrow.update({"gate_verdict": t4M["base"]["gate"]["verdict"], "gate_codes": t4M["base"]["gate"]["codes"],
                         "dod_clean": t4M["base"]["dod_clean"], "cminus_sep_ratio_min": t4M["base"]["gate"]["cminus_sep_ratio_min"]})
        t1rows.append(Mrow)
        worst = None
        for r in o["rows"]:
            rr = dict(r)
            rr["v0c_n321_verdict"] = v0c.get((round(o["L_frac"], 3), float(r["theta_max_deg"]), float(r["theta_e_deg"])))
            if FM is not None and r.get("F") is not None:
                rr["F_minus_FM"] = r["F"] - FM
                rr["F_minus_FM_over_F1D"] = (r["F"] - FM) / F1
                if r.get("gate_pass"):
                    if worst is None or rr["F_minus_FM"] > worst["F_minus_FM"]:
                        worst = rr
            t1rows.append(rr)
        npass = sum(1 for r in o["rows"] if r.get("gate_pass"))
        t1sum[name] = {"L_frac": o["L_frac"], "n_eval": o["n_eval"], "F_M": FM, "n_D": len(o["rows"]), "n_D_gate_pass": npass,
                       "max_F_D_minus_F_M_over_F1D (gate pass)": None if worst is None else worst["F_minus_FM_over_F1D"],
                       "argmax": None if worst is None else (worst["theta_max_deg"], worst["theta_e_deg"]),
                       "all_D_le_M (gate pass)": None if worst is None else bool(worst["F_minus_FM"] <= 0.0)}
    out = {"plan": "plans/active/tooling-nozzle-bell-two-mode.md §5.1 #5 (g2)(g3)・§6 V1 T4・T1a′",
           "method": {"thrust": "rao_axisym.forward_thrust (前進 MOC + 初期線 + 壁の圧力積分、壁の依存域の検査)",
                      "M_wall": "RaoWall pspline (5 次・48 区間・λ GCV、B で C1)", "families": FAMILIES, "k": KS, "eps": EPSS,
                      "eps0": EPS0, "fit": "u = ε/ε₀ の 3 次の最小二乗 (7 点)", "n_eval": N_EVALS, "F1D_abs": F1,
                      "control": "ctrl_plus = σ_D ×1.02 の D⁺ の Rao 線の壁を x = L で切る / ctrl_minus = σ_D ×0.98 (切らない、参考)",
                      "line": f"|c1 ε0| ≤ {RATIO_LINE} |c2 ε0²|、|c1(M)| ≤ {RATIO_LINE} |c1(ctrl_plus)|、c2 の符号"},
           **git_info(), "t4_rows": t4rows, "t1ap_rows": t1rows, "t1ap_summary": t1sum, "parts": parts}
    (HERE / "rao_t4.json").write_text(json.dumps(_jsonable({k: v for k, v in out.items() if k not in ("t1ap_rows", "t1ap_summary")}),
                                                 ensure_ascii=False, indent=1))
    (HERE / "rao_t1ap.json").write_text(json.dumps(_jsonable({k: out[k] for k in ("plan", "method", "commit", "source_sha256_16",
                                                                                    "design_tree_dirty", "t1ap_rows", "t1ap_summary")}),
                                                   ensure_ascii=False, indent=1))
    write_csv(HERE / "rao_t4.csv", t4rows)
    write_csv(HERE / "rao_t1ap.csv", t1rows)
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
    ap.add_argument("--assemble-only", action="store_true")
    a = ap.parse_args()
    if a.workers > 2:
        raise SystemExit("--workers は 2 まで (手元 PC はユーザの作業機)")
    t0 = time.time()
    if not a.assemble_only:
        todo = [j for j in jobs() if (a.only is None or any(s in j["name"] for s in a.only.split(",")))
                and not (PARTS / f"{j['name']}.json").exists()]
        ready = [j for j in todo if load_gen(j["gen"]) is not None]
        print(f"回すジョブ {len(ready)} (生成の部品が無くて待つもの {len(todo) - len(ready)})", flush=True)
        with ProcessPoolExecutor(max_workers=a.workers) as ex:
            futs = [ex.submit(run, j) for j in ready]
            for fu in as_completed(futs):
                fu.result()
    out = assemble()
    print(f"rao_t4.csv {len(out['t4_rows'])} 行・rao_t1ap.csv {len(out['t1ap_rows'])} 行 ({time.time() - t0:.0f} s)", flush=True)


if __name__ == "__main__":
    main()
