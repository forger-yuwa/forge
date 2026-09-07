"""⑤ SERN の力係数 (forge 壁面出力 res_wall_<physID>_<step>.h5 から)。

plan §4.7。壁面出力 (outputHDFflg: 1) は MESH/COORD (n,3) と VALUE/Ps (n) [node: ノード値 / cell: 面値]
を持つ。ランプ・カウルは x に単調なので x ソートした折れ線で (p − p_a) を台形積分する。
規約 (moc_sern.wall_forces と同じ): 推力 = 壁力の −x 成分、揚力 = +y、モーメント = 頭上げ正 (基準点指定)。
壁摩擦 (twall_x/y) は forge の「流体に働く traction」規約 (viscousFlux_d.cu) に従い、壁力 = −Σ twall·|dl| として別枠で加算する。
"""
from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np


def polyline_forces(h5file, fluid_below: bool, p_a: float, x_ref: float, y_ref: float) -> dict:
    with h5py.File(h5file, "r") as f:
        xy = f["MESH/COORD"][:].reshape(-1, 3)[:, :2].astype(float)
        conn = f["MESH/CONNE"][:].reshape(-1, 4)[:, 2:4].astype(int)   # [型, 節点数, n_a, n_b]
        p = f["VALUE/Ps"][:].astype(float)
        tw = None
        if "VALUE/twall_x" in f:
            tw = np.column_stack([f["VALUE/twall_x"][:], f["VALUE/twall_y"][:]]).astype(float)
    a, b = xy[conn[:, 0]], xy[conn[:, 1]]
    flip = a[:, 0] > b[:, 0]                # 接線を +x 向きに揃える
    a2 = np.where(flip[:, None], b, a); b2 = np.where(flip[:, None], a, b)
    d = b2 - a2
    ln = np.hypot(d[:, 0], d[:, 1])
    n_out = np.column_stack([-d[:, 1], d[:, 0]]) if fluid_below else np.column_stack([d[:, 1], -d[:, 0]])
    if len(p) == len(conn):                 # cell: 面値
        pf = p
        tf = tw
    elif len(p) == len(xy):                 # node: 節点値 → 面平均
        pf = 0.5 * (p[conn[:, 0]] + p[conn[:, 1]])
        tf = None if tw is None else 0.5 * (tw[conn[:, 0]] + tw[conn[:, 1]])
    else:
        raise ValueError(f"{h5file}: Ps 長 {len(p)} が面数 {len(conn)} とも節点数 {len(xy)} とも合わない")
    F = ((pf - p_a)[:, None] * n_out)       # |dl| 込み
    xm = 0.5 * (a2 + b2)
    Mz = (xm[:, 0] - x_ref) * F[:, 1] - (xm[:, 1] - y_ref) * F[:, 0]
    out = {"Fx_p": float(F[:, 0].sum()), "Fy_p": float(F[:, 1].sum()), "M_noseup_p": float(-Mz.sum()),
           "length": float(ln.sum()), "n_faces": int(len(conn)), "p_mean": float(np.sum(pf * ln) / ln.sum())}
    if tf is not None:
        tm = tf * ln[:, None]
        out.update({"Fx_tau": float(tm[:, 0].sum()), "Fy_tau": float(tm[:, 1].sum())})
        # 剥離指標: 壁接線 (+x 向き) 方向のせん断が負 (逆流) な面の長さ割合。符号規約は呼び出し側で補正 (tau_sign)
        tvec = np.column_stack([d[:, 0], d[:, 1]]) / ln[:, None]
        tt = (tf[:, 0] * tvec[:, 0] + tf[:, 1] * tvec[:, 1])
        out["_tau_t"] = tt; out["_ln"] = ln; out["_xm"] = xm[:, 0]
    return out


def sern_forces(run_dir, step: int, p_a: float, F_ideal: float, H: float, x_ref: float = 0.0, y_ref: float = 0.0,
                ids=None, mdot_u_in: float = 0.0, p_in: float = 0.0, twall_on_fluid: bool = True) -> dict:
    """1 出力ステップの力係数。F_ideal [N/m]、H [m]。総推力 = 入口流れ推力 + 壁推力。
    twall_on_fluid: forge の twall の規約。cell (viscousFlux_d.cu L527: tau = −τ_w ê_t を流体へ) は True、
    node (L679: W の CV に入る力 = 壁ノードに働く力) は False。case/26 平板 (node, 流れ +x) で twall_x>0、
    case/46 run_0011 (cell) で twall_x<0 を確認 (2026-09-04)。"""
    ids = ids or {"ramp": 4, "cowl_in": 5, "cowl_out": 6}
    run_dir = Path(run_dir)
    parts = {}
    for name, pid in ids.items():
        cands = list(run_dir.glob(f"res_*_{pid}_{step}.h5"))
        if not cands:
            continue
        fpath = cands[0]
        parts[name] = polyline_forces(fpath, fluid_below=(name != "cowl_in"), p_a=p_a, x_ref=x_ref, y_ref=y_ref)
    Fx = sum(v["Fx_p"] for v in parts.values()); Fy = sum(v["Fy_p"] for v in parts.values())
    Mn = sum(v["M_noseup_p"] for v in parts.values())
    T_wall = -Fx
    F_gross = mdot_u_in + (p_in - p_a) * H + T_wall
    out = {"step": step, "T_wall": T_wall, "L": Fy, "M_noseup": Mn, "F_gross": F_gross, "F_ideal": F_ideal,
           "C_T": F_gross / F_ideal, "C_T_wall": T_wall / F_ideal, "C_L": Fy / F_ideal, "C_M": Mn / (F_ideal * H),
           "parts": parts}
    if any("Fx_tau" in v for v in parts.values()):
        # forge の twall は「流体に働く」traction (cell: tau = −τ_w ê_t、node: W の CV に入る力 —
        # viscousFlux_d.cu L527/L679、2026-09-04 確認)。壁に働く摩擦力 = −Σ twall·|dl| で、流れ方向 (+x)
        # の抗力になる → 推力寄与は +Σ twall_x·|dl| (負)。揚力寄与は −Σ twall_y·|dl|。
        sgn = 1.0 if twall_on_fluid else -1.0     # 壁に働く摩擦力 = −sgn·Σ twall·|dl|
        Fx_t = sgn * sum(v.get("Fx_tau", 0.0) for v in parts.values()); Fy_t = sgn * sum(v.get("Fy_tau", 0.0) for v in parts.values())
        out["C_T_with_shear"] = (F_gross + Fx_t) / F_ideal
        out["C_L_with_shear"] = (Fy - Fy_t) / F_ideal
        out["C_T_friction"] = Fx_t / F_ideal
        # 剥離: 壁に働く接線力が −x (逆流) の長さ割合と最上流位置 (ランプ・カウル内面)
        for name in ("ramp", "cowl_in"):
            v = parts.get(name)
            if v and "_tau_t" in v:
                on_wall = -sgn * v["_tau_t"]          # 壁に働く接線力 (+x = 付着流)
                rev = on_wall < 0.0
                v["sep_frac"] = float(v["_ln"][rev].sum() / v["_ln"].sum())
                v["sep_x_min"] = float(v["_xm"][rev].min()) if rev.any() else None
                for k in ("_tau_t", "_ln", "_xm"):
                    v.pop(k, None)
        out["sep_frac_ramp"] = parts.get("ramp", {}).get("sep_frac")
        out["sep_x_min_ramp"] = parts.get("ramp", {}).get("sep_x_min")
    for v in parts.values():                      # 作業配列は JSON に出さない
        for k in ("_tau_t", "_ln", "_xm"):
            v.pop(k, None)
    return out


def force_history(run_dir, **kw) -> list:
    run_dir = Path(run_dir)
    steps = sorted({int(f.stem.split("_")[-1]) for f in run_dir.glob("res_*_[0-9]*_[0-9]*.h5") if not f.stem.startswith("res_nan")})
    return [sern_forces(run_dir, s, **kw) for s in steps]


def steadiness(values, tail: float = 0.4, drift: float = 0.02, osc: float = 0.05) -> dict:
    """末尾 tail 割合の時系列で頭打ち判定 (check_quasisteady と同じ語彙)。"""
    v = np.asarray(values, dtype=float)
    if len(v) < 4:
        return {"verdict": "TRANSIENT-UNSETTLED", "reason": "snapshots < 4"}
    n = max(3, int(np.ceil(len(v) * tail)))
    t = v[-n:]
    ref = max(abs(np.mean(t)), 1e-12)
    slope = np.polyfit(np.arange(n), t, 1)[0] * (n - 1) / ref
    fluct = (t.max() - t.min()) / ref
    if abs(slope) > drift:
        return {"verdict": "DRIFTING", "trend_frac": float(slope), "fluct_frac": float(fluct), "mean": float(t.mean())}
    if fluct > osc:
        return {"verdict": "OSCILLATING", "trend_frac": float(slope), "fluct_frac": float(fluct),
                "mean": float(t.mean()), "amp": float(0.5 * (t.max() - t.min()))}
    return {"verdict": "STEADY", "trend_frac": float(slope), "fluct_frac": float(fluct), "mean": float(t.mean())}
