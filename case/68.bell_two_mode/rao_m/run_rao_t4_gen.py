#!/usr/bin/env python3
r"""③ベル (M) モード §5.1 #5 (g): T4 停留性・T1a′ の (M) 壁と負対照の生成、T1c (自由 ε ↔ 固定 ε) の run。

plan: plans/active/tooling-nozzle-bell-two-mode.md §5.1 #5 (g)・§6 V1 の T1a′・T1c・T4・「TP の停留性」。
判断の記録: notes/reviews/2026-10-11-bell-rao-mainbranch-c2-diagnostician.md ((M) = 円弧 kernel・B で C1、受理 det_D > 0)。

生成器 `geometry/rao_axisym.py` の `design_rao` を (g1) の後の既定 (**Goursat の m = q = n**、円弧 kernel、受理 det_D > 0) で回す。
設定は T6 (`run_rao_t6.py`) と同じ: ガス = 基準ケースの凍結 TP (`RocketExhaust.moc_gas` → `Isentrope`)、R_t 1、ε 50、
初期線 = Hall の n 点、L_k = 等面積比 15° 円錐の長さ (√ε − 1)/tan15° の 100 %・70 %、p_a 0。

run (1 行 = 1 設計、2 本ずつ並列):
- `gen_L100_n161`・`gen_L70_n161`: T4 の主設定と T1a′ の (M) 壁。(g1) の記録として T6 の n 161 (m = q = 80) との差を出す。
- `gen_L100_n241`: T4 の「n 241 でも 1 組」。
- `t1c_free_L70_n161`: 自由 ε (`eps=None`、縁条件 `edge_function` = 0、p_a 0)。L70 は固定 ε 50 の解で縁の関数が
  −9.5e-5 (rao_lk_n161.json) と 0 に近いので、自由 ε の根が ε 50 の近くにある。
- `t1c_fixed_L70_n161`: 自由 ε の解の ε_free = r_E² で固定 ε として解き直す (free の後に投入)。

**負対照** (plan §6 V1 T4「D を 2 % ずらした非 Rao 線から作った壁」、gen_* の各 run で同じプロセス内で作る): 根の BD (kernel 場の
C⁻) の上で、B から D までの BD の弧長 σ_D を 1.02 倍にした点 D⁺ (軸寄り) と 0.98 倍の点 D⁻ (B 寄り) から DE を Φ₂・Φ₃ 一定で
引き直し (流量は ṁ_BD(D'→B))、領域 BD'E' を同じ m・q で Goursat、壁を流量で取る。D'E' はその (L', ε') の Rao 線で
(L, ε) の Rao 線ではない。評価 (`run_rao_t4.py`) は x_E' > L の側 (D⁺ の見込み) を x = L で切った壁を負対照に使う
(切った壁の制御面 = 切り口を通る C⁺ は Φ₂・Φ₃ 一定でないので、端を固定した摂動に対して停留でない)。他方の側は
「その (L', ε') の停留輪郭」として記録だけ残す。

出力 (このディレクトリ):
- `rao_t4_parts/<tag>.json`・`<tag>_walls.npz` (BE の点列 [x, r, θ, M]、B・E、負対照 D⁺・D⁻ の BE' と E')。既にあれば回さない。
- `rao_t4_gen.json`: 全 run の要約、(g1) の T6 n 161 との差、T1c の比較 (自由 ε → 固定 ε のビット一致の検査)、来歴。
- `rao_t1c.csv`: T1c の 2 行と差の行。
実行: `nice -n 10 /home/sano/work/forge/design/.venv-opt/bin/python case/68.bell_two_mode/rao_m/run_rao_t4_gen.py`
(n 161・241 は 1 本 25〜30 分)。`--only` で一部、`--assemble-only` で部品から表だけ作り直す。
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
from scipy.optimize import brentq  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "design"))
from forge_design.gas.rocket_exhaust import RocketExhaust  # noqa: E402
from forge_design.geometry import rao_axisym as ra  # noqa: E402
from forge_design.geometry.wall_bell_c2 import cone15_length  # noqa: E402

EXHAUST = HERE.parent / "exhaust_lox_lh2_pc100_of6.yaml"
PARTS = HERE / "rao_t4_parts"
T6_PARTS = HERE / "rao_t6_parts"
EPS = 50.0
R_T = 1.0
P_A = 0.0
CTRL_SHIFT = 0.02                    # 負対照: BD の弧長 σ_D (B → D) を 1 ± 0.02 倍
SOURCES = ("design/forge_design/geometry/rao_axisym.py", "design/forge_design/gas/isentrope.py",
           "design/forge_design/geometry/moc_forward.py", "design/forge_design/geometry/moc_kernel.py")

JOBS = [{"tag": "gen_L100_n161", "kind": "gen", "L_frac": 1.0, "n": 161},
        {"tag": "gen_L70_n161", "kind": "gen", "L_frac": 0.7, "n": 161},
        {"tag": "t1c_free_L70_n161", "kind": "t1c_free", "L_frac": 0.7, "n": 161},
        {"tag": "gen_L100_n241", "kind": "gen", "L_frac": 1.0, "n": 241}]
T1C_FIXED = {"tag": "t1c_fixed_L70_n161", "kind": "t1c_fixed", "L_frac": 0.7, "n": 161}


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


def _gas_and_line(n: int):
    ig = ra.as_isentrope(RocketExhaust.from_file(EXHAUST).moc_gas())
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        line = ra._hall_line(R_T, ig, n)
    return ig, line


def _record(res: dict, spec: dict, L: float, t_s: float) -> dict:
    """`design_rao` の戻りから表に載せる量 (run_rao_t6.py の `_record` と同じ列 + 縁の関数)。"""
    rec = {"tag": spec["tag"], "kind": spec["kind"], "L_frac": spec["L_frac"], "L_rt": L, "n": res["n"], "eps_in": res["eps"],
           "dr": res["dr"], "m": res["m"], "q": res["q"], "kernel": res["kernel"], "c": res["c"],
           "accept_rule": res.get("accept_rule"), "status": res["status"], "why": res.get("why"),
           "rejected": res.get("rejected", []), "lost": res.get("lost", []), "policy_holes": res.get("policy_holes", []),
           "time_s": t_s}
    if res["status"] != "ok":
        return rec
    E, th, cl, gt = res["E"], res["thrust"], res["closure"], res["g_tilde"]
    rec.update({
        "theta_B_deg": float(np.rad2deg(res["theta_B"])), "theta_B": float(res["theta_B"]), "s_B": res["s_B"],
        "eta_D": res["eta_D"], "eta0_BD": res["eta0_BD"], "N_BD": res["N_BD"], "det_D": res["det_D"], "de": res["de"],
        "theta_E_le_theta_B": res["theta_E_le_theta_B"], "de_lt_0": res["de_lt_0"],
        "M_D": float(res["D"]["M"]), "theta_D_deg": float(np.rad2deg(res["D"]["th"])), "x_D": float(res["D"]["x"]),
        "r_D": float(res["D"]["r"]),
        "theta_E_deg": float(np.rad2deg(res["theta_E"])), "theta_E": float(res["theta_E"]),
        "theta_max_deg": float(np.rad2deg(res["theta_max"])),
        "x_E": float(E["x"]), "r_E": float(E["r"]), "eps_E": float(E["r"]) ** 2, "M_E": float(E["M"]),
        "x_E_minus_L": float(E["x"]) - L, "x_B": res["x_B"], "r_B": res["r_B"], "edge": res.get("edge"),
        "F_cs": th["F_cs"], "F_wall": th["F_wall"], "thrust_rel_diff": th["rel_diff"], "mdot_rel_diff": th["mdot_rel_diff"],
        "T_line": th["T_line"], "T_kernel": th["T_kernel"], "T_BE": th["T_BE"], "mdot_line": th["mdot_line"],
        "CF": res["CF"], "F_1D": res["F_1D"], "isp_ratio": res["isp_ratio"],
        "closure_dx": cl["dx"], "closure_dr": cl["dr"], "closure_dth": cl["dth"], "closure_norm": float(np.hypot(cl["dx"], cl["dr"])),
        "kappa_B": res["kappa_B"], "kappaW_B": res["kappaW_B"], "dkappa_B": res["dkappa_B"], "g_tilde": gt["g"],
        "goursat_nfail": res["goursat_nfail"], "n_wall_missing": res["n_wall_missing"],
        "streamline_resid_max": res["streamline_resid_max"], "n_candidates": len(res["candidates"]),
    })
    return rec


def _bd_sigma(bd, eta: float, n_pan: int = 400) -> float:
    """BD の弧長 σ (η → B = N)。"""
    eg, wg = ra._gl(np.linspace(eta, bd.N, n_pan + 1)[:-1], np.linspace(eta, bd.N, n_pan + 1)[1:])
    d = bd.dspl(eg)
    return float(np.sum(np.hypot(d[..., 0], d[..., 1]) * wg))


def _control(ig, res: dict, factor: float) -> dict:
    """負対照: BD の弧長 σ_D を factor 倍にした D' から DE (Φ₂・Φ₃ 一定) → Goursat → 壁。戻りは記録と壁の点列。"""
    mr = res["_mid"]
    bd = mr["bd"]
    eta_D = float(mr["D_root"]["eta"])
    sig_D = _bd_sigma(bd, eta_D)
    target = factor * sig_D
    lo, hi = bd.eta0 + 1e-9, float(bd.N) - 1e-9
    f = lambda e: _bd_sigma(bd, e) - target  # noqa: E731
    if not f(lo) * f(hi) < 0.0:
        return {"ok": False, "why": f"σ = {target:.6g} の点が BD [{bd.eta0}, {bd.N}] に無い"}
    eta_c = float(brentq(f, lo, hi, xtol=1e-13))
    dd = ra._de_at(bd, eta_c, P_A, res["dr"])
    E = dd["E"]
    out = {"factor": factor, "eta_D": eta_c, "sigma_D": sig_D, "sigma_Dp": _bd_sigma(bd, eta_c), "eta_D_root": eta_D,
           "D": {k: float(v) for k, v in dd["D"].items()}, "det_D": float(dd["rl"].det_D) if dd["rl"].ok_D else None}
    if not E.get("ok"):
        return {**out, "ok": False, "why": f"DE' が引けない: {E.get('why')}"}
    D_root = {"eta": eta_c, "D": dd["D"], "E": E, "rl": dd["rl"]}
    P = ra.goursat_fill(ig, bd, D_root, res["m"], res["q"])
    nfail = ra.goursat_fill.last_nfail
    wl = ra.wall_from_massflow(ig, bd, D_root, P, P_A)
    out.update({"ok": True, "x_E": float(E["x"]), "r_E": float(E["r"]), "eps_E": float(E["r"]) ** 2, "th_E": float(E["th"]),
                "M_E": float(E["M"]), "mdot": float(E["mdot"]), "x_E_minus_L": float(E["x"]) - res["L"],
                "goursat_nfail": nfail, "n_wall_missing": wl["n_missing"], "closure": wl["closure"],
                "streamline_resid_max": wl["streamline_resid_max"], "_BE": wl["wall"]})
    return out


def run_job(spec: dict) -> dict:
    t0 = time.time()
    L = spec["L_frac"] * cone15_length(EPS)
    try:
        ig, line = _gas_and_line(spec["n"])
        eps = None if spec["kind"] == "t1c_free" else float(spec.get("eps", EPS))
        res = ra.design_rao(ig, R_t=R_T, L=L, eps=eps, kernel="arc", line=line, p_a=P_A)
        rec = _record(res, spec, L, time.time() - t0)
        arrays = {}
        if res["status"] == "ok":
            arrays["BE"] = res["BE"]
            arrays["kernel_wall"] = res["kernel_wall"]
            arrays["B_xr"] = np.array([res["x_B"], res["r_B"]])
            arrays["E_xrthM"] = np.array([res["E"]["x"], res["E"]["r"], res["E"]["th"], res["E"]["M"]])
            if spec["kind"] == "gen":
                ctrl = {}
                for name, fac in (("ctrl_plus", 1.0 + CTRL_SHIFT), ("ctrl_minus", 1.0 - CTRL_SHIFT)):
                    t1 = time.time()
                    try:
                        c = _control(ig, res, fac)
                    except Exception as e:  # noqa: BLE001
                        c = {"ok": False, "why": f"{type(e).__name__}: {e}", "traceback": traceback.format_exc()}
                    c["time_s"] = time.time() - t1
                    if c.get("ok"):
                        arrays[f"{name}_BE"] = c.pop("_BE")
                        arrays[f"{name}_E_xrthM"] = np.array([c["x_E"], c["r_E"], c["th_E"], c["M_E"]])
                    ctrl[name] = c
                rec["controls"] = ctrl
        PARTS.mkdir(exist_ok=True)
        (PARTS / f"{spec['tag']}.json").write_text(json.dumps(_jsonable(rec), ensure_ascii=False, indent=1))
        if arrays:
            np.savez(PARTS / f"{spec['tag']}_walls.npz", **arrays)
    except Exception as e:  # noqa: BLE001
        rec = {"tag": spec["tag"], "kind": spec["kind"], "L_frac": spec["L_frac"], "n": spec["n"], "status": "error",
               "why": f"{type(e).__name__}: {e}", "traceback": traceback.format_exc(), "time_s": time.time() - t0}
        PARTS.mkdir(exist_ok=True)
        (PARTS / f"{spec['tag']}.json").write_text(json.dumps(_jsonable(rec), ensure_ascii=False, indent=1))
    print(_line(rec), flush=True)
    return rec


def _line(rec: dict) -> str:
    if rec.get("status") != "ok":
        return f"{rec['tag']}: {rec.get('status')} — {rec.get('why')} ({rec.get('time_s', 0):.0f} s)"
    s = (f"{rec['tag']}: ok m {rec['m']} θ_B {rec['theta_B_deg']:.6f}° θ_E {rec['theta_E_deg']:.6f}° θ_max {rec['theta_max_deg']:.6f}° "
         f"x_E {rec['x_E']:.9f} ε_E {rec['eps_E']:.9f} 縁 {rec['edge']:+.4e} F_cs {rec['F_cs']:.9f} 2 会計の差 "
         f"{rec['thrust_rel_diff']:+.3e} 閉合 {rec['closure_norm']:.3e} Δκ_B {rec['dkappa_B']:+.5f} ({rec['time_s']:.0f} s)")
    for name, c in (rec.get("controls") or {}).items():
        s += (f"\n    {name}: " + (f"ok η_D' {c['eta_D']:.5f} (根 {c['eta_D_root']:.5f}) x_E' {c['x_E']:.5f} (−L {c['x_E_minus_L']:+.4f}) "
                                  f"ε_E' {c['eps_E']:.4f} θ_E' {np.rad2deg(c['th_E']):.4f}°" if c.get("ok") else f"不可 — {c.get('why')}"))
    return s


# ================================================================================================ 集計
G1_KEYS = ("theta_B_deg", "theta_E_deg", "eta_D", "x_E", "r_E", "F_cs", "theta_max_deg", "F_wall", "thrust_rel_diff",
           "closure_norm", "kappaW_B", "dkappa_B", "g_tilde", "streamline_resid_max")
T1C_KEYS = ("theta_B_deg", "theta_B", "s_B", "eta_D", "det_D", "theta_E_deg", "theta_E", "theta_max_deg", "x_E", "r_E", "eps_E",
            "M_E", "edge", "F_cs", "F_wall", "thrust_rel_diff", "closure_dx", "closure_dr", "kappaW_B", "dkappa_B")


def g1_vs_t6(rows: dict) -> dict:
    """(g1) m = q = n の n 161 と T6 の n 161 (m = q = 80) の差 (新 − 旧)。"""
    out = {}
    for lab in ("L100", "L70"):
        new = rows.get(f"gen_{lab}_n161")
        p = T6_PARTS / f"{lab}_n161.json"
        if new is None or new.get("status") != "ok" or not p.exists():
            out[lab] = {"why": "どちらかが無い"}
            continue
        old = json.loads(p.read_text())
        out[lab] = {"m_new": new["m"], "m_old": old["m"],
                    "diff_new_minus_old": {k: new[k] - old[k] for k in G1_KEYS if new.get(k) is not None and old.get(k) is not None},
                    "new": {k: new.get(k) for k in G1_KEYS}, "old": {k: old.get(k) for k in G1_KEYS}}
    return out


def t1c_compare(rows: dict) -> dict:
    fr, fx = rows.get("t1c_free_L70_n161"), rows.get("t1c_fixed_L70_n161")
    if fr is None or fx is None or fr.get("status") != "ok" or fx.get("status") != "ok":
        return {"why": "自由 ε か固定 ε の run が無い / ok でない",
                "free_status": None if fr is None else fr.get("status"), "fixed_status": None if fx is None else fx.get("status"),
                "free_why": None if fr is None else fr.get("why"), "fixed_why": None if fx is None else fx.get("why")}
    d = {k: fx[k] - fr[k] for k in T1C_KEYS if fr.get(k) is not None and fx.get(k) is not None}
    bit = {k: bool(fx[k] == fr[k]) for k in T1C_KEYS if fr.get(k) is not None and fx.get(k) is not None}
    return {"eps_free": fr["eps_E"], "eps_fixed_in": fx["eps_in"], "diff_fixed_minus_free": d, "bit_identical": bit,
            "all_bit_identical": all(bit.values()), "max_abs_diff_angles_deg": max(abs(d[k]) for k in ("theta_B_deg", "theta_E_deg")),
            "free": {k: fr.get(k) for k in T1C_KEYS}, "fixed": {k: fx.get(k) for k in T1C_KEYS}}


def assemble(total_s=None) -> dict:
    rows = {}
    for j in JOBS + [T1C_FIXED]:
        p = PARTS / f"{j['tag']}.json"
        if p.exists():
            rows[j["tag"]] = json.loads(p.read_text())
    out = {"plan": "plans/active/tooling-nozzle-bell-two-mode.md §5.1 #5 (g)・§6 V1 T1a′・T1c・T4",
           "setup": {"gas": "EXH (RocketExhaust.moc_gas → Isentrope)", "exhaust": str(EXHAUST.relative_to(ROOT)), "R_t": R_T,
                     "eps": EPS, "p_a": P_A, "initial_line": "Hall n 点", "kernel": "arc", "accept_rule": "det_D > 0",
                     "mq_default": "n (g1)", "dr_default": "3.2/n", "L_cone15_rt": cone15_length(EPS),
                     "control": f"BD の弧長 σ_D (B → D) を 1 ± {CTRL_SHIFT} 倍した D' から Rao 線 (その L', ε' の) → Goursat → 壁"},
           **git_info(), "rows": rows, "g1_vs_t6_n161": g1_vs_t6(rows), "t1c": t1c_compare(rows), "total_s": total_s}
    (HERE / "rao_t4_gen.json").write_text(json.dumps(_jsonable(out), ensure_ascii=False, indent=1))
    # T1c の表
    cols = ("tag", "eps_in") + T1C_KEYS + ("time_s", "why")
    with open(HERE / "rao_t1c.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for t in ("t1c_free_L70_n161", "t1c_fixed_L70_n161"):
            r = rows.get(t)
            if r is not None:
                w.writerow([("" if r.get(c) is None else (f"{r[c]:.17g}" if isinstance(r[c], float) else r[c])) for c in cols])
        tc = out["t1c"]
        if "diff_fixed_minus_free" in tc:
            w.writerow(["fixed_minus_free", ""] + [("" if tc["diff_fixed_minus_free"].get(c) is None
                                                    else f"{tc['diff_fixed_minus_free'][c]:.3e}") for c in T1C_KEYS] + ["", ""])
    return out


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
        todo = [j for j in JOBS if (a.only is None or j["tag"] in a.only.split(",")) and not (PARTS / f"{j['tag']}.json").exists()]
        want_fixed = (a.only is None or T1C_FIXED["tag"] in a.only.split(",")) and not (PARTS / f"{T1C_FIXED['tag']}.json").exists()
        print(f"回す run: {[j['tag'] for j in todo]} + {'t1c_fixed (free の後)' if want_fixed else ''} (並列 {a.workers})", flush=True)
        with ProcessPoolExecutor(max_workers=a.workers) as ex:
            pend = {ex.submit(run_job, j): j for j in todo}
            fixed_sent = False

            def try_fixed():
                nonlocal fixed_sent
                p = PARTS / "t1c_free_L70_n161.json"
                if want_fixed and not fixed_sent and p.exists():
                    fr = json.loads(p.read_text())
                    if fr.get("status") == "ok":
                        spec = {**T1C_FIXED, "eps": fr["eps_E"]}
                        pend[ex.submit(run_job, spec)] = spec
                    else:
                        print(f"t1c_fixed は回さない: 自由 ε の run が {fr.get('status')} — {fr.get('why')}", flush=True)
                    fixed_sent = True
            try_fixed()
            while pend:
                done, _ = wait(list(pend), return_when=FIRST_COMPLETED)
                for fu in done:
                    pend.pop(fu)
                    fu.result()
                try_fixed()
                assemble(time.time() - t0)
    out = assemble(time.time() - t0)
    print(f"rao_t4_gen.json ({len(out['rows'])} run) total {time.time() - t0:.0f} s", flush=True)
    print(json.dumps(_jsonable({"g1_vs_t6_n161": {k: v.get("diff_new_minus_old", v) for k, v in out["g1_vs_t6_n161"].items()},
                                "t1c": {k: v for k, v in out["t1c"].items() if k not in ("free", "fixed")}}),
                     ensure_ascii=False, indent=1), flush=True)


if __name__ == "__main__":
    main()
