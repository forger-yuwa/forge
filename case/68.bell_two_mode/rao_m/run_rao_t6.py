#!/usr/bin/env python3
r"""③ベル (M) モード: 円弧 kernel の TP の Rao 輪郭の分解能 (T6) と保存 (T3) の系列 — 記録、合否の判定はしない
(plan tooling-nozzle-bell-two-mode §5.1 #5 (f)・§6 V1 の T3・T5・T6 [2026-10-11 改]。判断の記録
notes/reviews/2026-10-11-bell-rao-mainbranch-c2-diagnostician.md)。

生成器 `geometry/rao_axisym.py` の `design_rao` を §5.1 #5 (e) の後 (受理条件 det_D > 0、θ_E ≤ θ_B・de < 0 はフラグ、
既定 = 円弧 kernel、C2 の外の根なし) で回す。

設定: ガス = 基準ケースの凍結 TP (`RocketExhaust.moc_gas`、../exhaust_lox_lh2_pc100_of6.yaml)、R_t 1、ε 50 (固定)、
初期線 = Hall の n 点、kernel = 円弧 (曲率半径 R_t)、L_k = 等面積比 15° 円錐の長さ (√ε − 1)/tan15° の 100 % と 70 %
(成立域の両端)。DE の r 格子 dr・Goursat の点数 m・q・走査点数は `design_rao` の既定 (dr 3.2/n、m = q = n//2)。

run (1 行 = 1 設計):
- `base`: n {161, 241, 321} × L {100, 70} %。n 161 は `rao_lk_n161.*` (前の規則 θ_E ≤ θ_B) と同じ設定を今のコードで回し直し、
  値の一致を `n161_vs_rao_lk` に記録する (規則の変更で n 161 の解が変わっていないことの確認)。
- `dr_half` (n 321): DE の r 格子だけ半分 (dr = 3.2/321/2)。中の根から解き直す (dr は E の位置、したがって根に効く)。
- `mq_double` (n 321): Goursat の m・q だけ倍 (320)。m・q は根の探索に入らない (`design_rao` で Goursat → 壁 → 推力の段にだけ
  渡る) ので、`base` n 321 と同じ根 (同じ kernel 場・BD・DE) で Goursat → 壁 → 推力 → κ_W(B⁺) を作り直す
  (`design_rao` を m・q を倍にして回したのと同じ計算)。kernel 場 (n) と DE/Goursat (dr, m, q) の分解能を分けて見る。

出力 (このディレクトリ):
- `rao_t6_parts/<tag>.json`・`<tag>_wall.npz`: run ごとの記録 (途中で止まっても再開できる。既にあれば回さない) と壁の点列。
- `rao_t6.json`: 全行 + 解析 (T3: 推力の 2 会計の差と閉合 |W_m − E| の n 161/241/321 の値・単調性・観測次数・n 321 で
  ≤ 1e-4 か。T6: θ_B・θ_E・θ_max・F の 3 水準の Richardson 外挿 [Celik 2008 の非一様比]、n 321 の dr 半分・m・q 倍の差と
  n 241 → 321 の差。T5: g̃・Δκ_B の系列)、n 161 の旧結果との一致、来歴 (commit・ソースの sha256)。
- `rao_t6.csv`: 1 行 1 run の表。

並列は 2 本まで (手元 PC はユーザの作業機)。`nice -n 10` で起動する:
`nice -n 10 design/.venv-opt/bin/python case/68.bell_two_mode/rao_m/run_rao_t6.py` (n 321 は 1 本 1〜1.5 時間、dr 半分は約倍)。
`--only L100_n161,...` で一部、`--assemble-only` で部品から表だけ作り直す。
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
from forge_design.geometry.wall_bell_c2 import cone15_length  # noqa: E402

EXHAUST = HERE.parent / "exhaust_lox_lh2_pc100_of6.yaml"
PARTS = HERE / "rao_t6_parts"
EPS = 50.0
R_T = 1.0
P_A = 0.0
L_FRACS = (1.0, 0.7)
NS = (161, 241, 321)
N_REF = 321
OLD_N161 = HERE / "rao_lk_n161.json"
SOURCES = ("design/forge_design/geometry/rao_axisym.py", "design/forge_design/gas/isentrope.py",
           "design/forge_design/geometry/moc_forward.py", "design/forge_design/geometry/moc_kernel.py")


def lab(frac: float) -> str:
    return f"L{int(round(100 * frac))}"


# 回す順 (2 本ずつ同じ長さの組): n 161 → 241 → 321 (+ m・q 倍) → 321 の dr 半分
JOBS = ([{"tag": f"{lab(f)}_n{n}", "L_frac": f, "n": n, "variant": "base"} for n in NS for f in L_FRACS]
        + [{"tag": f"{lab(f)}_n{N_REF}_drhalf", "L_frac": f, "n": N_REF, "variant": "dr_half"} for f in L_FRACS])


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
    """`design_rao` の戻りから表に載せる量 (スカラー) と中の根の記録を抜く。"""
    rec = {"tag": spec["tag"], "L_frac": spec["L_frac"], "L_rt": L, "n": res["n"], "variant": spec["variant"],
           "dr": res["dr"], "m": res["m"], "q": res["q"], "kernel": res["kernel"], "c": res["c"],
           "accept_rule": res.get("accept_rule"), "status": res["status"], "why": res.get("why"),
           "eta_tol": res.get("eta_tol"), "rejected": res.get("rejected", []), "lost": res.get("lost", []),
           "policy_holes": res.get("policy_holes", []), "branches": res.get("branches", {}),
           "scan_mid": res.get("scan_mid", []), "time_s": t_s}
    if res["status"] != "ok":
        return rec
    E, th, cl, gt = res["E"], res["thrust"], res["closure"], res["g_tilde"]
    rec.update({
        "theta_B_deg": float(np.rad2deg(res["theta_B"])), "s_B": res["s_B"], "eta_D": res["eta_D"],
        "eta_D_over_nm1": res["eta_D"] / (res["n"] - 1), "eta0_BD": res["eta0_BD"], "N_BD": res["N_BD"],
        "det_D": res["det_D"], "de": res["de"], "theta_E_le_theta_B": res["theta_E_le_theta_B"], "de_lt_0": res["de_lt_0"],
        "branch": res["branch"], "M_D": float(res["D"]["M"]), "theta_D_deg": float(np.rad2deg(res["D"]["th"])),
        "x_D": float(res["D"]["x"]), "r_D": float(res["D"]["r"]),
        "theta_E_deg": float(np.rad2deg(res["theta_E"])), "theta_max_deg": float(np.rad2deg(res["theta_max"])),
        "x_E": float(E["x"]), "r_E": float(E["r"]), "M_E": float(E["M"]), "x_E_minus_L": float(E["x"]) - L,
        "x_B": res["x_B"], "r_B": res["r_B"],
        "F_cs": th["F_cs"], "F_wall": th["F_wall"], "thrust_rel_diff": th["rel_diff"], "mdot_rel_diff": th["mdot_rel_diff"],
        "T_line": th["T_line"], "T_kernel": th["T_kernel"], "T_BE": th["T_BE"], "T_closure": th["T_closure"], "T_DE": th["T_DE"],
        "CF": res["CF"], "isp_ratio": res["isp_ratio"],
        "closure_dx": cl["dx"], "closure_dr": cl["dr"], "closure_dth": cl["dth"], "closure_norm": float(np.hypot(cl["dx"], cl["dr"])),
        "kappa_B": res["kappa_B"], "kappaW_B": res["kappaW_B"], "dkappa_B": res["dkappa_B"],
        "g_tilde": gt["g"], "g_de": gt["de"], "g_kernel": gt["kernel"], "g_cplus_angle_err": gt["cplus_angle_err"],
        "goursat_nfail": res["goursat_nfail"], "n_wall_missing": res["n_wall_missing"],
        "streamline_resid_max": res["streamline_resid_max"], "candidates": res["candidates"], "net_K": res["net_K"],
    })
    return rec


def _rebuild_mq(res: dict, ig, m2: int, q2: int) -> dict:
    """同じ根 (`res` の kernel 場・BD・DE) で Goursat を m2・q2 にして壁 → 推力 → κ_W(B⁺) を作り直す (`design_rao` の後段と同じ式)。"""
    net, mr = res["_net"], res["_mid"]
    P = ra.goursat_fill(ig, mr["bd"], mr["D_root"], m2, q2)
    nfail = ra.goursat_fill.last_nfail
    wl = ra.wall_from_massflow(ig, mr["bd"], mr["D_root"], P, P_A)
    th = ra.thrust_accounts(net, mr, wl, P_A)
    BE = wl["wall"]
    kW = ra.wall_curvature_B(BE)
    F1 = ig.ideal_1d(EPS, P_A)
    kern = res["kernel_wall"]
    full = np.vstack([kern, BE[1:]]) if abs(BE[0, 0] - kern[-1, 0]) < 1e-6 else np.vstack([kern, BE])
    out = dict(res)
    out.update({"m": m2, "q": q2, "BE": BE, "wall": full, "thrust": th, "closure": wl["closure"], "n_wall_missing": wl["n_missing"],
                "streamline_resid_max": wl["streamline_resid_max"], "goursat_nfail": nfail, "kappaW_B": kW,
                "dkappa_B": kW - mr["kappa_B"],
                "theta_max": max(float(np.max(res["kernel_wall"][:, 2])), float(np.max(BE[:, 2]))),
                "CF": th["F_cs"] / (np.pi * ig.p0),
                "isp_ratio": (th["F_cs"] / th["mdot_line"]) / (F1["F"] / F1["mdot"])})
    return out


def _save(rec: dict, res: dict | None) -> None:
    PARTS.mkdir(exist_ok=True)
    (PARTS / f"{rec['tag']}.json").write_text(json.dumps(_jsonable(rec), ensure_ascii=False, indent=1))
    if res is not None and res.get("status") == "ok":
        np.savez(PARTS / f"{rec['tag']}_wall.npz", BE=res["BE"], kernel_wall=res["kernel_wall"], wall=res["wall"])


def _line(rec: dict) -> str:
    if rec["status"] != "ok":
        return f"{rec['tag']}: {rec['status']} — {rec.get('why')} ({rec['time_s']:.0f} s)"
    hole = f" POLICY_HOLE {len(rec['policy_holes'])}" if rec.get("policy_holes") else ""
    return (f"{rec['tag']}: ok θ_B {rec['theta_B_deg']:.4f}° θ_E {rec['theta_E_deg']:.4f}° θ_max {rec['theta_max_deg']:.4f}° "
            f"det_D {rec['det_D']:+.4f} de {rec['de']:+.4g} x_E {rec['x_E']:.6f} 推力差 {rec['thrust_rel_diff']:+.3e} "
            f"閉合 {rec['closure_norm']:.3e} Δκ_B {rec['dkappa_B']:+.4f} g̃ {rec['g_tilde']:+.4e} "
            f"候補 {len(rec['candidates'])}・捨て {len(rec['rejected'])}・lost {len(rec['lost'])}{hole} ({rec['time_s']:.0f} s)")


def run_job(spec: dict) -> list:
    """1 run (`base` n 321 は `mq_double` の行も)。部品を書いて記録の list を返す。例外は status `error` の行にする。"""
    out = []
    t0 = time.time()
    L = spec["L_frac"] * cone15_length(EPS)
    try:
        ig, line = _gas_and_line(spec["n"])
        dr = 3.2 / spec["n"] / 2.0 if spec["variant"] == "dr_half" else None
        res = ra.design_rao(ig, R_t=R_T, L=L, eps=EPS, kernel="arc", line=line, p_a=P_A, dr=dr)
        rec = _record(res, spec, L, time.time() - t0)
        _save(rec, res)
        print(_line(rec), flush=True)
        out.append(rec)
        if spec["variant"] == "base" and spec["n"] == N_REF and res["status"] == "ok":
            t1 = time.time()
            r2 = _rebuild_mq(res, ig, 2 * res["m"], 2 * res["q"])
            spec2 = {**spec, "tag": spec["tag"] + "_mq2", "variant": "mq_double"}
            rec2 = _record(r2, spec2, L, time.time() - t1)
            rec2["note"] = (f"{spec['tag']} と同じ根 (kernel 場・BD・DE) で Goursat の m・q だけ倍。time_s は作り直しの分だけ "
                            f"(根の探索を含む設計 1 本は {rec['time_s']:.0f} s)")
            _save(rec2, r2)
            print(_line(rec2), flush=True)
            out.append(rec2)
    except Exception as e:  # noqa: BLE001
        rec = {"tag": spec["tag"], "L_frac": spec["L_frac"], "L_rt": L, "n": spec["n"], "variant": spec["variant"],
               "status": "error", "why": f"{type(e).__name__}: {e}", "traceback": traceback.format_exc(),
               "time_s": time.time() - t0}
        _save(rec, None)
        print(_line(rec), flush=True)
        out.append(rec)
    return out


# ================================================================================================ 解析
def _h(n: int) -> float:
    """初期線の区間幅 (r_t 単位、Hall n 点 → n−1 区間)。dr = 3.2/n と m = q = n//2 も n に比例して細かくなる。"""
    return 1.0 / (n - 1)


def richardson(f1: float, f2: float, f3: float, h1: float, h2: float, h3: float) -> dict:
    """3 水準 (f1 が最も細かい) の観測次数と外挿値 (Celik et al. 2008、非一様の細分比を固定点反復で)。"""
    r21, r32 = h2 / h1, h3 / h2
    e21, e32 = f2 - f1, f3 - f2
    out = {"f_fine": f1, "e21": e21, "e32": e32, "ratio_e32_e21": (e32 / e21) if e21 != 0 else None,
           "monotone": bool(e21 * e32 > 0.0), "r21": r21, "r32": r32}
    if e21 == 0.0 or e32 == 0.0:
        out.update({"p": None, "f_ext": None, "why": "差が 0"})
        return out
    s = np.sign(e32 / e21)
    p = abs(np.log(abs(e32 / e21))) / np.log(r21)
    conv = False
    for _ in range(200):
        try:
            qv = np.log((r21 ** p - s) / (r32 ** p - s))
        except FloatingPointError:
            break
        p_new = abs(np.log(abs(e32 / e21)) + qv) / np.log(r21)
        if not np.isfinite(p_new):
            break
        if abs(p_new - p) <= 1e-12 * max(1.0, p):
            p, conv = p_new, True
            break
        p = p_new
    if not (conv and np.isfinite(p) and p > 0):
        out.update({"p": float(p) if np.isfinite(p) else None, "f_ext": None, "why": "次数の反復が収束しない"})
        return out
    f_ext = (r21 ** p * f1 - f2) / (r21 ** p - 1.0)
    out.update({"p": float(p), "f_ext": float(f_ext), "f_fine_minus_ext": float(f1 - f_ext),
                "rel_fine_minus_ext": float((f1 - f_ext) / f_ext) if f_ext != 0 else None})
    return out


def _series_order(vals: dict) -> dict:
    """|e(n)| の単調減少と観測次数 ln(|e_a|/|e_b|)/ln(h_a/h_b) (隣どうしと 161 → 321)。"""
    ns = sorted(vals)
    a = {n: abs(vals[n]) for n in ns}
    orders = {}
    for na, nb in list(zip(ns[:-1], ns[1:])) + [(ns[0], ns[-1])]:
        if a[na] > 0 and a[nb] > 0:
            orders[f"{na}->{nb}"] = float(np.log(a[na] / a[nb]) / np.log(_h(na) / _h(nb)))
    return {"values": {str(n): vals[n] for n in ns}, "abs": {str(n): a[n] for n in ns},
            "monotone_decreasing": bool(all(a[x] > a[y] for x, y in zip(ns[:-1], ns[1:]))), "orders": orders,
            "abs_at_n321_le_1e-4": bool(a[ns[-1]] <= 1e-4), "min_order": min(orders.values()) if orders else None}


T6_QTY = ("theta_B_deg", "theta_E_deg", "theta_max_deg", "F_cs", "F_wall", "CF", "eta_D_over_nm1", "dkappa_B", "g_tilde",
          "kappaW_B")
REF_QTY = T6_QTY + ("thrust_rel_diff", "closure_norm", "mdot_rel_diff")


def analyse(rows: list) -> dict:
    by = {r["tag"]: r for r in rows}
    out = {}
    for f in L_FRACS:
        L_lab = lab(f)
        base = {n: by.get(f"{L_lab}_n{n}") for n in NS}
        ok = all(b is not None and b.get("status") == "ok" for b in base.values())
        a = {"complete": ok}
        if not ok:
            a["why"] = {str(n): (None if b is None else b.get("why")) for n, b in base.items() if b is None or b.get("status") != "ok"}
            out[L_lab] = a
            continue
        # 枝の一貫性 (n で解の枝が変わっていないか): det_D の符号・候補数・θ_B と η_D/(n−1) の変化
        a["branch_check"] = {
            "det_D": {str(n): base[n]["det_D"] for n in NS}, "all_det_D_pos": all(base[n]["det_D"] > 0 for n in NS),
            "n_candidates": {str(n): len(base[n]["candidates"]) for n in NS},
            "n_rejected": {str(n): len(base[n]["rejected"]) for n in NS},
            "n_lost": {str(n): len(base[n]["lost"]) for n in NS},
            "n_policy_holes": {str(n): len(base[n]["policy_holes"]) for n in NS},
            "theta_E_le_theta_B": {str(n): base[n]["theta_E_le_theta_B"] for n in NS},
            "de_lt_0": {str(n): base[n]["de_lt_0"] for n in NS},
            "eta_D_over_nm1": {str(n): base[n]["eta_D_over_nm1"] for n in NS},
        }
        a["T3"] = {"thrust_rel_diff": _series_order({n: base[n]["thrust_rel_diff"] for n in NS}),
                   "closure_norm": _series_order({n: base[n]["closure_norm"] for n in NS}),
                   "mdot_rel_diff": _series_order({n: base[n]["mdot_rel_diff"] for n in NS})}
        a["T6_richardson"] = {k: richardson(base[321][k], base[241][k], base[161][k], _h(321), _h(241), _h(161))
                              for k in T6_QTY}
        a["T5_record"] = {k: {"values": {str(n): base[n][k] for n in NS},
                              "d_161_241": base[241][k] - base[161][k], "d_241_321": base[321][k] - base[241][k]}
                          for k in ("g_tilde", "dkappa_B", "kappaW_B")}
        ref = {}
        for v, tag in (("dr_half", f"{L_lab}_n{N_REF}_drhalf"), ("mq_double", f"{L_lab}_n{N_REF}_mq2")):
            r = by.get(tag)
            if r is None or r.get("status") != "ok":
                ref[v] = {"why": None if r is None else r.get("why"), "status": None if r is None else r.get("status")}
                continue
            ref[v] = {k: {"value": r[k], "minus_base321": r[k] - base[321][k], "base321_minus_base241": base[321][k] - base[241][k]}
                      for k in REF_QTY}
        a["n321_refinement"] = ref
        out[L_lab] = a
    return out


def compare_old_n161(rows: list) -> dict:
    """n 161 の回し直しと `rao_lk_n161.json` (前の規則) の値の差。"""
    if not OLD_N161.exists():
        return {"why": "rao_lk_n161.json が無い"}
    old = {round(r["L_frac"], 3): r for r in json.loads(OLD_N161.read_text())["L_k"]}
    keys = ("theta_B_deg", "eta_D", "det_D", "theta_E_deg", "theta_max_deg", "x_E", "r_E", "M_E", "thrust_rel_diff",
            "mdot_rel_diff", "F_cs", "F_wall", "kappaW_B", "n_wall_missing", "goursat_nfail")
    out = {}
    for r in rows:
        if r.get("variant") != "base" or r.get("n") != 161:
            continue
        o = old.get(round(r["L_frac"], 3))
        if o is None or r.get("status") != "ok" or o.get("status") != "ok":
            out[r["tag"]] = {"why": "旧結果に同じ L が無いか、どちらかが ok でない",
                             "new_status": r.get("status"), "old_status": None if o is None else o.get("status")}
            continue
        d = {}
        for k in keys:
            if r.get(k) is not None and o.get(k) is not None:
                d[k] = r[k] - o[k]
        if o.get("closure"):
            d["closure_dx"] = r["closure_dx"] - o["closure"]["dx"]
            d["closure_dr"] = r["closure_dr"] - o["closure"]["dr"]
        d["dkappa_B_vs_old_kappaW_minus_kappa"] = r["dkappa_B"] - o["kappaW_B_minus_kappa_B"]
        out[r["tag"]] = {"diff_new_minus_old": d, "max_abs": max(abs(v) for v in d.values()),
                         "identical": all(v == 0 for v in d.values()), "old_status": o["status"], "new_status": r["status"],
                         "old_n_candidates": len(o.get("candidates") or []), "old_n_rejected": len(o.get("rejected") or [])}
    return out


CSV_COLS = ("tag", "L_frac", "L_rt", "n", "variant", "dr", "m", "q", "status", "theta_B_deg", "theta_E_deg", "theta_max_deg",
            "x_E", "x_E_minus_L", "r_E", "M_E", "eta_D", "eta_D_over_nm1", "det_D", "de", "theta_E_le_theta_B", "de_lt_0",
            "F_cs", "F_wall", "thrust_rel_diff", "mdot_rel_diff", "CF", "isp_ratio", "closure_dx", "closure_dr",
            "closure_norm", "kappaW_B", "dkappa_B", "g_tilde", "g_de", "g_kernel", "n_candidates", "n_rejected", "n_lost",
            "n_policy_holes", "goursat_nfail", "n_wall_missing", "streamline_resid_max", "time_s", "why")


def csv_row(r: dict) -> dict:
    row = {k: r.get(k) for k in CSV_COLS}
    row["n_candidates"] = len(r.get("candidates") or [])
    row["n_rejected"] = len(r.get("rejected") or [])
    row["n_lost"] = len(r.get("lost") or [])
    row["n_policy_holes"] = len(r.get("policy_holes") or [])
    return row


def write_csv(path: Path, rows: list) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(CSV_COLS))
        w.writeheader()
        for r in rows:
            w.writerow({k: ("" if r.get(k) is None else (f"{r[k]:.10g}" if isinstance(r[k], float) else r[k])) for k in CSV_COLS})


def assemble(total_s: float | None = None) -> dict:
    order = [j["tag"] for j in JOBS]
    order = order[:4] + order[4:6] + [order[4] + "_mq2", order[5] + "_mq2"] + order[6:]
    # 並べ替え: 表は L ごと・n 昇順・細分の順
    order = sorted(order, key=lambda t: (-int(t.split("_")[0][1:]), int(t.split("_")[1][1:]), len(t)))
    rows = []
    for t in order:
        p = PARTS / f"{t}.json"
        if p.exists():
            rows.append(json.loads(p.read_text()))
    out = {"plan": "plans/active/tooling-nozzle-bell-two-mode.md §5.1 #5 (f)・§6 V1 T3・T5・T6 (2026-10-11 改)",
           "decision_record": "notes/reviews/2026-10-11-bell-rao-mainbranch-c2-diagnostician.md",
           "setup": {"gas": "EXH (RocketExhaust.moc_gas)", "exhaust": str(EXHAUST.relative_to(ROOT)), "R_t": R_T, "eps": EPS,
                     "p_a": P_A, "initial_line": "Hall n 点", "kernel": "arc", "accept_rule": "det_D > 0",
                     "L_cone15_rt": cone15_length(EPS), "L_fracs": list(L_FRACS), "ns": list(NS),
                     "dr_default": "3.2/n", "mq_default": "n//2", "h_for_order": "1/(n-1)",
                     "N_SCAN_SB": ra.N_SCAN_SB, "N_SCAN_D": ra.N_SCAN_D, "ETA_LINK_FRAC": ra.ETA_LINK_FRAC,
                     "THETA_B_RANGE_DEG": ra.THETA_B_RANGE_DEG,
                     "note": "合否は判定しない (親と上位)。T3 の線: n 321 で ≤ 1e-4・n 161/241/321 で単調減少・次数 ≥ 1"},
           **git_info(), "rows": rows, "analysis": analyse(rows), "n161_vs_rao_lk": compare_old_n161(rows),
           "policy_holes_any": any(bool(r.get("policy_holes")) for r in rows), "total_s": total_s}
    (HERE / "rao_t6.json").write_text(json.dumps(_jsonable(out), ensure_ascii=False, indent=1))
    write_csv(HERE / "rao_t6.csv", [csv_row(r) for r in rows])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default=None, help="回す tag (カンマ区切り)")
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--assemble-only", action="store_true")
    a = ap.parse_args()
    if a.workers > 2:
        raise SystemExit("--workers は 2 まで (手元 PC はユーザの作業機)")
    t0 = time.time()
    if not a.assemble_only:
        jobs = [j for j in JOBS if a.only is None or j["tag"] in a.only.split(",")]
        jobs = [j for j in jobs if not (PARTS / f"{j['tag']}.json").exists()]
        print(f"回す run: {[j['tag'] for j in jobs]} (並列 {a.workers})", flush=True)
        with ProcessPoolExecutor(max_workers=a.workers) as ex:
            futs = [ex.submit(run_job, j) for j in jobs]
            for fu in as_completed(futs):
                for rec in fu.result():
                    if rec.get("policy_holes"):
                        print(f"!!! POLICY_HOLE {rec['tag']}: det_D > 0 なのに θ_E > θ_B の候補 "
                              f"{[h.get('policy_hole') for h in rec['policy_holes']]}", flush=True)
                assemble(time.time() - t0)
    total = time.time() - t0
    if a.assemble_only and (HERE / "rao_t6.json").exists():
        total = json.loads((HERE / "rao_t6.json").read_text()).get("total_s", total)   # run の総時間を保つ
    out = assemble(total)
    print(f"rao_t6.json・rao_t6.csv ({len(out['rows'])} 行、policy_holes {out['policy_holes_any']}) total {time.time() - t0:.0f} s",
          flush=True)


if __name__ == "__main__":
    main()
