"""case/71 M5.5 の手組みの NS 格子 (plan tooling-nozzle-card-local-mesh §5.1 #3a、ユーザ決定「a」2026-10-15)。生成器の既定は変えない。

build: カードから作った問題 YAML (problem_m55_card.yaml) の mesh・mesh_euler を、場所ごとの第一層の表と流れ方向の密度の表の格子に置き換えた
       問題 YAML (problem_m55_ns.yaml) を書く。
  1. 要求の第一層 y1,req(x): `wt_card.estimate_mesh` と同じ式 (設計壁の 1 次元等エントロピー・Eckert の参照温度・平板の Cf、y1+ 目標
     Y1PLUS_TARGET 0.8)。入口の端は s = S_C r_t の値で頭打ち (y1,req ∝ s^0.1 → 0 のため)。式の写しが estimate_mesh の断面 (8001 点) と
     一致することを確かめてから細かい刻みで使う。
  2. 係数 c(x) (初期格子の経験の係数; M6 の H の表 ÷ 見積もりの比の区間ごとの最小以下、`c45/h_over_est_ratio_summary.json`):
     配管 1.0・縮流部 0.6・スロート付近 0.5・拡大部 C_DIV。区間の境は smoothstep で、比の最小が大きい側 (上流側) の区間の中でつなぐ。
  3. 正規化した要求 f_req(x) = c(x)·y1,req(x) / max(r_設計壁(x), r_物理壁(x)) (物理壁 = prepare_ns の CFD 前の積分法の壁; 格子は物理壁に
     置くので実寸 = 表 × r_物理壁 ≤ c·y1,req になるように)。`wall_first_frac_table` = f_req の下包絡 (節点 0.125 r_t (x < 10)・1 r_t、
     各節点の値は前後の区間の f_req の最小を有効 6 桁で切り下げ; log 線形の補間が区間の最小を超えない)。表の値と補間の区間の中を
     細かい刻み (CHECK_DX) で検査し、超える区間は点を足す。
  4. `x_density_table`: 目標の間隔 h(x) = min(AR_MARGIN·ar_max·y1,表·cosθ_w, h_cap(x)) (物理壁; h_cap は雛形の分布
     (throat_refine 4・throat_width 3) で ni 2000 にしたときの各位置の間隔) の逆数の上包絡。AR (壁際のセルの max 辺 / min 辺) と
     dx ≤ h_cap を満たす最小の ni を二分法で求め、10 の倍数に切り上げる。throat_refine 1・wall_normal_layer [0.02, 0.8]・msh_digits 17。
  5. mesh_euler の ni は NS と同じ (雛形の規則)。max_nodes・max_nodes_reason を両方に書く: --approval (承認の記録) があればカードの
     承認済みの mesh_limits の値と承認元 (数えた節点数が超えれば問題 YAML を書かずに止める)、無ければ入力の問題 YAML の上限のまま
     (超えれば止める; 自動で広げない。plan tooling-design-orchestrator §4.6・§5.1 #8、2026-10-18 から)。
  5b. --axis-cap-frac c を渡したとき (#3b、ユーザ 2026-10-15「上限3%」): NS の径方向を mesh2d の axis_cap_frac (nj 固定版: 各断面で壁の
     第一セルから比 q_i の等比、間隔 / r_w が c に達したら軸まで一様) にし、nj は全断面で `mesh2d._radial_fracs_capfixed` の隣接比の
     最大 ≤ --radial-ratio-max (既定 1.2) になる最小の値。流れ方向・第一層の表・mesh_euler (nj を含む) は変えない。
  6. 最終の雛形 (axismach_final_v1: physical_wall_repr・pw_upstream を持つ; §5.1 #3b) の問題なら、Euler 用の写し <out の stem>_euler.yaml を
     生成器と同じ `wt_card.write_euler_companion` で書き、生成器が書いた <problem の stem>_euler.yaml と mesh・mesh_euler 以外が同じことを確かめる。
check: 置き換えた問題 YAML の NS の格子を prepare_ns と同じ経路 (設計チェーン → 積分法の δ_r → build_physical_wall → mesh_params →
       generate_axisym_mesh → msh 17 桁 → 変換) で作り、forge を起動せずに検査する: check_mesh_quality --ar-max・mesh_ar_relax_check・変換後の
       第一層厚と表の照合・見積もりの y1+・節点数・dx の比。Euler の格子 (prepare と同じ経路) も節点数と品質。
AWS A で ~/forge-cardmesh/design を import して走らせる (`--design` で変えられる)。"""
from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

S_C = 0.011            # 入口の端の頭打ち [r_t] (見積もりの断面の 1 点目)
C_PIPE, C_CONTR, C_THROAT = 1.0, 0.6, 0.5
C_DIV = 0.49           # 拡大部 (M6 の比の最小 0.49874 を有効 2 桁で切り下げ; 0.5 は比の最小を超える)
BLEND = 1.0            # 区間の境の smoothstep の幅 [r_t]
X_THROAT_HALF = 2.0    # スロート付近 |x| ≤ 2 (見積もり・c45_hcheck.py の区間と同じ)
NI_REF, REF_REFINE, REF_WIDTH = 2000, 4.0, 3.0      # h_cap の基準 (雛形の流れ方向の分布)
BUILD_DX, CHECK_DX = 2.0e-4, 2.0e-5                  # 下包絡を作る刻み・検査の刻み [r_t]
WALL_NORMAL_LAYER = [0.02, 0.8]
NORM = "max"           # f_req の半径: max(r_設計壁, r_物理壁)。--sanity (見積もりの再現) では設計壁 (design)
MSH_DIGITS = 17


def _imports(design: str):
    sys.path.insert(0, str(Path(design).expanduser()))
    global W, M2, RA, load_problem
    from forge_design.intake import wt_card as W               # noqa: F401
    from forge_design.meshing import mesh2d as M2              # noqa: F401
    from forge_design.evaluate import runner_axismach as RA    # noqa: F401
    from forge_design.probdef import load_problem              # noqa: F401


def smooth(t):
    t = np.clip(t, 0.0, 1.0)
    return t * t * (3.0 - 2.0 * t)


def coef(x, x_U: float):
    """c(x): 配管 C_PIPE → [x_U − BLEND, x_U] で C_CONTR → [−2 − BLEND, −2] で C_THROAT → [2 − BLEND, 2] で C_DIV。"""
    x = np.asarray(x, float)
    c = np.full(x.shape, C_PIPE)
    c = c + (C_CONTR - C_PIPE) * smooth((x - (x_U - BLEND)) / BLEND)
    c = c + (C_THROAT - C_CONTR) * smooth((x - (-X_THROAT_HALF - BLEND)) / BLEND)
    c = c + (C_DIV - C_THROAT) * smooth((x - (X_THROAT_HALF - BLEND)) / BLEND)
    return c


def y1req_design(x, wall, card, rt: float):
    """`wt_card.estimate_mesh` と同じ式の y1,req [m] (x は設計壁の入口端から昇順、s は x の格子で累積する)。等温一様の壁だけ。"""
    from forge_design.evaluate.ic import _invert_area_ratio_gas, _isentropic_pressure_ratio_gas
    if card["wall"]["thermal"] != "isothermal" or "Tw_profile" in card["wall"]:
        raise ValueError("等温一様の壁だけ (estimate_mesh の他の枝は写していない)")
    gas = W._gas(card)
    Pt, Tt = float(card["conditions"]["Pt"]), float(card["conditions"]["Tt"])
    rw = wall.r(x)
    s = np.concatenate([[0.0], np.cumsum(np.hypot(np.diff(x), np.diff(rw)))]) * rt
    AR = np.maximum(rw, 1.0) ** 2
    Mx = _invert_area_ratio_gas(AR, x >= 0.0, gas)
    T = np.where(Mx >= 1.0, gas.T_of_M(np.maximum(Mx, 1.0)), np.interp(np.minimum(Mx, 1.0), gas._Mu, gas._Tu))
    P = Pt * _isentropic_pressure_ratio_gas(T, gas)
    Rg = float(gas.R)
    u = Mx * np.sqrt(gas.gamma(T) * Rg * T)
    Tw = np.full_like(T, float(card["wall"]["Tw"]))
    Ts = T * (1.0 + 0.032 * Mx ** 2 + 0.58 * (Tw / T - 1.0))
    rho_s, mu_s = P / (Rg * Ts), W._mu_sutherland(Ts)
    rho_w, mu_w = P / (Rg * Tw), W._mu_sutherland(Tw)
    with np.errstate(divide="ignore", invalid="ignore"):
        Re = rho_s * u * s / mu_s
        Cf = 0.0592 * Re ** -0.2
        tau = Cf * 0.5 * rho_s * u * u
        utau = np.sqrt(tau / rho_w)
        y1 = W.Y1PLUS_TARGET * mu_w / (rho_w * utau)
    return y1, s / rt


def walls(problem_path):
    """(問題, 設計チェーン, 設計壁, 物理壁 (prepare_ns の CFD 前の積分法の壁), 積分法の記録)。prepare_ns・deltastar_loop --init-integral と同じ構築。"""
    from forge_design.feedback.deltastar_loop import resolve_integral_initializer
    p = load_problem(problem_path)
    d = RA.design_chain(p)
    RA.require_moc_gate(p, d)
    scale = float(p.spec["r_throat"])
    init = resolve_integral_initializer(problem_path)
    _res, delta_r_x, info = RA.integral_delta_r(p, d, init)
    pwu = RA._pw_upstream(p.geometry)
    wp = RA.build_physical_wall(p, d, scale, delta_r_x=delta_r_x, dstar_x=None, offset="radial", pwu=pwu)
    msgs = wp.validate()
    if msgs:
        raise ValueError("物理壁フィルタ不合格: " + "; ".join(msgs))
    return p, d, d["wall"], wp, {"initializer": init, "model": info.get("model") if isinstance(info, dict) else None}


class Req:
    """要求 f_req(x) = c(x)·y1,req(x)/max(r_d, r_p) を任意の刻みの格子で作る (設計壁の範囲の外は端の値で保持)。"""

    def __init__(self, card, wd, wp, rt, x_U):
        self.card, self.wd, self.wp, self.rt, self.x_U = card, wd, wp, rt, x_U
        self.xd0, self.xd1 = float(wd.x_in), float(wd.x_e)
        self.xp0, self.xp1 = float(wp.x_in), float(wp.x_e)
        self.X0, self.X1 = min(self.xd0, self.xp0), max(self.xd1, self.xp1)
        # 入口の端の頭打ちの値: 配管は平らなので s = x − x_in。[x_in, x_in + S_C] の 2 点の式で s = S_C の y1 を取る
        yc, sc = y1req_design(np.array([self.xd0, self.xd0 + S_C]), wd, card, rt)
        if abs(sc[1] - S_C) > 1e-12 or abs(float(wd.r(np.array([self.xd0 + S_C]))[0]) - float(wd.r(np.array([self.xd0]))[0])) > 1e-12:
            raise ValueError("入口の端が平らな配管でない (s = x − x_in の前提が崩れる)")
        self.y1_cap_raw = float(yc[1])

    def on(self, dx):
        n = int(math.ceil((self.xd1 - self.xd0) / dx)) + 1
        xd = np.linspace(self.xd0, self.xd1, n)
        y1_raw, s = y1req_design(xd, self.wd, self.card, self.rt)
        # 注: _isentropic_pressure_ratio_gas の積分の格子は評価する点の T の最小に依るので、端の値 (s = S_C) は別の 2 点で取った値
        #     ではなく、この格子の上で s = S_C を補間した値に揃える (頭打ちの値は格子に依らずほぼ同じ; 差は記録する)
        y1c = y1_raw.copy()
        k = int(np.searchsorted(s, S_C))
        y_sc = float(np.interp(S_C, s[k - 1:k + 1], y1_raw[k - 1:k + 1]))
        y1c[s < S_C] = y_sc
        ext_lo = np.arange(self.X0, self.xd0, dx) if self.X0 < self.xd0 else np.empty(0)
        ext_hi = (np.arange(self.xd1, self.X1, dx)[1:] if self.X1 > self.xd1 else np.empty(0))
        x = np.concatenate([ext_lo, xd, ext_hi, [self.X1] if self.X1 > xd[-1] + 1e-12 else []])
        y1 = np.concatenate([np.full(len(ext_lo), y1c[0]), y1c, np.full(len(x) - len(ext_lo) - len(xd), y1c[-1])])
        y1r = np.concatenate([np.full(len(ext_lo), np.nan), y1_raw, np.full(len(x) - len(ext_lo) - len(xd), np.nan)])
        rd = self.wd.r(np.clip(x, self.xd0, self.xd1))
        rp = self.wp.r(np.clip(x, self.xp0, self.xp1))
        c = coef(x, self.x_U)
        f = c * y1 / ((np.maximum(rd, rp) if NORM == "max" else rd) * self.rt)
        return {"x": x, "y1": y1, "y1_raw": y1r, "rd": rd, "rp": rp, "c": c, "f": f, "y_sc": y_sc}


def table_nodes(X0, X1):
    tx = np.concatenate([np.arange(X0, min(10.0, X1), 0.125), np.arange(10.0, X1, 1.0)]) if X1 > 10 else np.arange(X0, X1, 0.125)
    return np.unique(np.concatenate([tx[tx < X1 - 1e-6], [X1]]))


def window_env(tx, x, v, mode):
    """各節点の前後の区間 [tx_{k−1}, tx_{k+1}] の最小 (mode min) / 最大 (max)。"""
    idx = np.searchsorted(x, tx)
    out = np.empty(len(tx))
    for k in range(len(tx)):
        lo = idx[max(k - 1, 0)]
        hi = idx[min(k + 1, len(tx) - 1)]
        lo = max(lo - 1, 0)
        hi = min(hi + 1, len(x))
        seg = v[lo:hi]
        out[k] = seg.min() if mode == "min" else seg.max()
    return out


def floor_sig(v, n):
    e = np.floor(np.log10(v))
    q = 10.0 ** (e - (n - 1))
    return np.floor(v / q * (1 + 1e-15)) * q


def ceil_sig(v, n):
    e = np.floor(np.log10(v))
    q = 10.0 ** (e - (n - 1))
    return np.ceil(v / q * (1 - 1e-15)) * q


def fmt_table(t, sig):
    return "[" + ", ".join(f"[{float(a)!r}, {float(b):.{sig}g}]" for a, b in t) + "]"


def parse_table(txt):
    return np.asarray(json.loads(txt), dtype=float)


def interp_f(t, x):
    return np.exp(np.interp(x, t[:, 0], np.log(t[:, 1])))


def build(a):
    t0 = time.time()
    # 承認の照合は最初の全格子の生成 (下の estimate_mesh) の前 (plan tooling-design-orchestrator §4.6「照合の場所と順」の 1; 2026-10-18 codex
    # 諮問 mesh-limits の Major 1): カードが承認の原本か承認の記録を指す試行の入力か、カードだけで決まる保護の条件が承認と同じかを確かめて、
    # 節点数の上限 (mesh_limits) を取る。結び付かなければここで止まる (何も作らない)
    lim = None
    if a.approval:
        from forge_design.intake import approval as AP
        lim = AP.mesh_limits_for_card(a.approval, a.card)
    card = W.load_card(a.card)
    import yaml
    prob_dict = yaml.safe_load(open(a.problem))
    p, d, wd, wp, init_rec = walls(a.problem)
    rt = float(p.spec["r_throat"])
    x_U = -float(wd.up.L_U)
    ar_max = float(prob_dict["mesh"].get("ar_max", 5000))
    ARM = W.AR_MARGIN * ar_max
    nj = int(prob_dict["mesh"]["nj"])
    log = lambda *s: print(*s, flush=True)  # noqa: E731
    log(f"r_t {rt!r} m  設計壁 [{wd.x_in:.6f}, {wd.x_e:.6f}]  物理壁 [{wp.x_in:.6f}, {wp.x_e:.6f}]  x_U {x_U:.6f}  ar_max {ar_max:g} (AR 上限 {ARM:g})"
        f"  初期化 {init_rec}")
    # ---- 1. 式の写しの照合 (estimate_mesh の断面と同じ 8001 点で)
    # 断面の値 (sections) だけを使うので全格子は作らない (full_grid=False; plan tooling-design-orchestrator §5.1 #8、2026-10-18 codex 諮問 v2c-fail:
    # 上限を超える依頼でも、上限の照合の前に入力の雛形の全格子を作っていた)
    est = W.estimate_mesh(card, prob_dict, chain=d, full_grid=False)
    sec = est.pop("sections")
    y_chk, _ = y1req_design(sec["x_rt"], wd, card, rt)
    m = np.isfinite(sec["y1_m"]) & (sec["y1_m"] > 0)
    rel = float(np.max(np.abs(y_chk[m] / sec["y1_m"][m] - 1.0)))
    log(f"式の写し vs estimate_mesh (8001 点): 相対差の最大 {rel:.3e}")
    if rel > 1e-12:
        raise RuntimeError(f"y1,req の式の写しが estimate_mesh と一致しない ({rel:.3e})")
    # ---- 2. 係数の検査 (M5.5 の区間ごとの c の最大 ≤ M6 の比の区間の最小)
    rs = json.load(open(a.ratio_summary))["segments"]
    rmin = {s["seg"]: s["ratio_min"] for s in rs}
    m6 = {"pipe": rmin["入口の配管"], "contr": rmin["縮流部 x<-2"], "throat": rmin["スロート付近 |x|<=2"],
          "div": min(rmin["拡大部 2<x<=x_e/2"], rmin["拡大部 x>x_e/2"])}
    xx = np.linspace(min(wd.x_in, wp.x_in), max(wd.x_e, wp.x_e), 200001)
    cc = coef(xx, x_U)
    seg_m55 = {"pipe": (xx <= x_U), "contr": (xx >= x_U) & (xx <= -X_THROAT_HALF), "throat": np.abs(xx) <= X_THROAT_HALF,
               "div": xx >= X_THROAT_HALF}
    c_chk = {k: {"c_max": float(cc[mk].max()), "c_min": float(cc[mk].min()), "m6_ratio_min": m6[k],
                 "ok": bool(cc[mk].max() <= m6[k])} for k, mk in seg_m55.items()}
    for k, v in c_chk.items():
        log(f"係数 {k:<7s} c {v['c_min']:.4f}–{v['c_max']:.4f}  M6 の比の最小 {v['m6_ratio_min']:.5f}  {'OK' if v['ok'] else 'NG'}")
    if not all(v["ok"] for v in c_chk.values()) and not a.sanity:
        raise RuntimeError("補間後の c(x) が M6 の比の区間の最小を超える")
    # ---- 3. 第一層の表 (下包絡) と検査
    R = Req(card, wd, wp, rt, x_U)
    B = R.on(BUILD_DX)
    C = R.on(CHECK_DX)
    log(f"要求の格子: 作る {len(B['x'])} 点・検査 {len(C['x'])} 点; 頭打ちの値 y1(s = {S_C}) {B['y_sc']:.6e} / {C['y_sc']:.6e} m "
        f"(2 点の式 {R.y1_cap_raw:.6e} m)")
    tx = table_nodes(R.X0, R.X1)
    added = []
    for it in range(20):
        tv = floor_sig(window_env(tx, C["x"], C["f"], "min"), 6)
        ft = parse_table(fmt_table(np.stack([tx, tv], 1), 6))
        q = interp_f(ft, C["x"]) / C["f"]
        if q.max() <= 1.0:
            break
        # 超える区間に点を足す (その区間で比が最大の位置)
        bad_iv = np.unique(np.searchsorted(tx, C["x"][q > 1.0]) - 1)
        new = []
        for iv in bad_iv:
            mm = (C["x"] >= tx[iv]) & (C["x"] <= tx[min(iv + 1, len(tx) - 1)])
            new.append(float(C["x"][mm][np.argmax(q[mm])]))
        added += new
        tx = np.unique(np.concatenate([tx, new]))
    else:
        raise RuntimeError("第一層の表の検査が 20 回で収まらない")
    q_b = interp_f(ft, B["x"]) / B["f"]
    log(f"第一層の表 {len(ft)} 点 (足した点 {len(added)}): 表/要求 (検査の刻み {CHECK_DX}) の最大 {q.max():.6f}・最小 {q.min():.4f}; "
        f"(作る刻み) 最大 {q_b.max():.6f}")
    # ---- 4. 流れ方向の密度の表と ni (物理壁)
    xp0, xp1 = float(wp.x_in), float(wp.x_e)
    xg = np.linspace(xp0, xp1, 400001)
    rp = wp.r(xg)
    try:
        drp = np.asarray(wp.r(xg, 1), float)
    except TypeError:
        drp = np.gradient(rp, xg)
    cos = 1.0 / np.sqrt(1.0 + drp ** 2)
    dens_ref = 1.0 + (REF_REFINE - 1.0) * np.exp(-((xg / REF_WIDTH) ** 2))
    I_ref = float(np.trapezoid(dens_ref, xg))
    h_cap = I_ref / ((NI_REF - 1) * dens_ref)
    y1_tab = interp_f(ft, xg) * rp
    h_AR = ARM * y1_tab * cos
    h = np.minimum(h_AR, h_cap)
    dens = 1.0 / h
    tdx = table_nodes(R.X0, R.X1)
    dv = window_env(tdx, xg, dens, "max")
    dv = ceil_sig(dv / dv.min(), 5)
    dt = parse_table(fmt_table(np.stack([tdx, dv], 1), 5))

    def stations(ni):
        return M2._x_stations(xp0, xp1, int(ni), 1.0, 1.0, density_table=dt)

    def wall_ar(xs):
        r = wp.r(xs)
        y1 = interp_f(ft, xs) * r
        ch = np.hypot(np.diff(xs), np.diff(r))
        e = np.stack([ch, y1[:-1], y1[1:]])
        return e.max(0) / e.min(0)

    def ok(ni):
        xs = stations(ni)
        if wall_ar(xs).max() > ARM:
            return False
        xm = 0.5 * (xs[1:] + xs[:-1])
        return bool(np.all(np.diff(xs) <= np.interp(xm, xg, h_cap) * (1 + 1e-9)))
    lo, hi = 2, 64
    while not ok(hi):
        lo, hi = hi, hi * 2
    while hi - lo > 1:
        mid = (lo + hi) // 2
        lo, hi = (lo, mid) if ok(mid) else (mid, hi)
    ni_min = hi
    ni = int(math.ceil(ni_min / 10.0) * 10)
    if not ok(ni):
        raise RuntimeError(f"10 の倍数に切り上げた ni {ni} が条件を満たさない")
    xs = stations(ni)
    ar = wall_ar(xs)
    dxs = np.diff(xs)
    qd = dxs[1:] / dxs[:-1]
    kq = int(np.argmax(np.maximum(qd, 1 / qd)))
    ni_ideal = 1 + float(np.trapezoid(dens, xg))
    n_nodes = ni * nj
    log(f"ni 最小 {ni_min} → {ni} (理想 {ni_ideal:.0f})  節点 {n_nodes:,}  壁際 AR 最大 {ar.max():.1f} @ x {0.5 * (xs[np.argmax(ar)] + xs[np.argmax(ar) + 1]):.3f}"
        f"  dx {dxs.min():.5f}–{dxs.max():.4f} r_t  dx 比最大 {max(qd[kq], 1 / qd[kq]):.4f} @ x {xs[kq + 1]:.3f}")
    if a.sanity:
        log(f"--sanity: 係数 {a.sanity}・半径 {NORM} で ni {ni_min} → {ni}、節点 {n_nodes:,} (YAML は書かない)")
        return
    # ---- 4b. 径方向 (--axis-cap-frac): 全断面で隣接比の最大 ≤ radial_ratio_max の最小の nj
    radial_rec = None
    if a.axis_cap_frac is not None:
        cap, qlim = float(a.axis_cap_frac), float(a.radial_ratio_max)
        fr_st = interp_f(ft, xs)            # 各断面の第一セル / 局所半径 (mesh2d._first_frac_profile と同じ log 線形)

        def rad_ratio(nj_, frs):
            return max(M2._radial_fracs_capfixed(int(nj_), float(f), cap)[1] for f in frs)
        nj_lo = int(math.ceil(1.0 + (1.0 - float(fr_st.max())) / cap)) + 1      # 和が 1 になる最小 (fr + (n−1)c ≥ 1) の 1 つ上から
        fmin = [float(fr_st.min())]          # 比は第一セルが薄いほど大きい: 最小の第一セルで探し、全断面で確かめる
        lo_n, hi_n = nj_lo, nj_lo
        while rad_ratio(hi_n, fmin) > qlim:
            lo_n, hi_n = hi_n, hi_n * 2
        while hi_n - lo_n > 1:
            mid = (lo_n + hi_n) // 2
            lo_n, hi_n = (lo_n, mid) if rad_ratio(mid, fmin) <= qlim else (mid, hi_n)
        nj_cap = hi_n
        q_all = rad_ratio(nj_cap, fr_st)
        q_prev = rad_ratio(nj_cap - 1, fr_st)
        if not (q_all <= qlim < q_prev):
            raise RuntimeError(f"径方向の nj の決定が全断面で成り立たない (nj {nj_cap}: 比の最大 {q_all}、nj − 1: {q_prev})")
        nj = nj_cap
        n_nodes = ni * nj
        radial_rec = {"axis_cap_frac": cap, "radial_ratio_max": qlim, "nj": nj, "ratio_max_all_stations": q_all,
                      "ratio_max_nj_minus_1": q_prev, "first_frac_min": float(fr_st.min()), "first_frac_max": float(fr_st.max()),
                      "x_first_frac_min": float(xs[int(np.argmin(fr_st))]), "nj_template": int(prob_dict["mesh"]["nj"])}
        log(f"径方向: axis_cap_frac {cap}・比 ≤ {qlim} → nj {nj} (比の最大 {q_all:.5f}; nj {nj - 1} では {q_prev:.5f}; 最も薄い第一セル "
            f"{fr_st.min():.4e} @ x {radial_rec['x_first_frac_min']:.3f})  節点 {n_nodes:,} (雛形の nj {radial_rec['nj_template']})")
    # ---- 5. 置き換えた問題 YAML
    me = dict(prob_dict["mesh_euler"])
    # ---- 節点数の上限 (plan tooling-design-orchestrator §4.6・§5.1 #8): 書く上限は、承認 (--approval) があればカードの承認済みの mesh_limits と
    # 承認元、無ければ入力の問題 YAML の上限のまま (明示が無ければ既定 30 万)。実効の節点数 (axis_gap_frac・axis_cap_frac の導出の nj を含む) は
    # 下で、書く予定の YAML の節から共通の check_node_budget で数え、超えれば問題 YAML を書かない。以前は数えた点数から上限を自動で広げ、
    # 固定の文「ユーザ決定「a」2026-10-15」を書いていた (別の依頼でも過去の承認の文で上限を通せた)
    if lim is not None:
        src = (f"承認の記録 {Path(lim['approval_file']).name} (sha256 {lim['approval_sha256'][:16]}…、{lim['user_ref']}) の "
               "カードの mesh_limits.{key}")
        max_ns, max_eu = lim["max_nodes_ns"], lim["max_nodes_euler"]
        reason_ns = src.format(key="max_nodes_ns")
        reason_eu = src.format(key="max_nodes_euler")
    else:
        b_ns, b_eu = RA.node_budget(prob_dict.get("mesh") or {}, "mesh"), RA.node_budget(me, "mesh_euler")
        max_ns, max_eu = b_ns["max_nodes"], b_eu["max_nodes"]
        reason_ns = b_ns["max_nodes_reason"] or "既定 (runner_axismach.MAX_NODES_2D)"
        reason_eu = b_eu["max_nodes_reason"] or "既定 (runner_axismach.MAX_NODES_2D)"
    src = Path(a.problem).read_text()
    lines = src.splitlines(keepends=True)

    def block_span(key):
        i0 = next(i for i, ln in enumerate(lines) if ln.startswith(f"{key}:"))
        i1 = next((i for i in range(i0 + 1, len(lines)) if lines[i].strip() and not lines[i].startswith((" ", "#"))), len(lines))
        # 次のトップレベルの直前のコメント行・空行は次のブロックのもの
        while i1 - 1 > i0 and (not lines[i1 - 1].strip() or lines[i1 - 1].startswith("#")):
            i1 -= 1
        return i0, i1
    old_mesh = dict(prob_dict["mesh"])
    mesh_txt = (
        "mesh:\n"
        f"  # 手組みの格子 (plan tooling-nozzle-card-local-mesh {a.plan_item}、ユーザ決定「a」2026-10-15)。_mesh/hand_mesh.py build が作る (手で編集しない)。\n"
        "  # 生成器の見積もり (第一層 2 値) は使わない。第一層の表 = f_req = c(x)·y1,req(x)/max(r_設計壁, r_物理壁) の下包絡 (log 線形; 補間の区間の中も\n"
        f"  # 刻み {CHECK_DX} r_t で表 ≤ f_req を検査)。y1,req は estimate_mesh と同じ式 (y1+ 目標 0.8、入口の端は s = {S_C} r_t の値で頭打ち)、\n"
        f"  # c: 配管 {C_PIPE}・縮流部 {C_CONTR}・スロート付近 {C_THROAT}・拡大部 {C_DIV} (初期格子の経験の係数: M6 の H の表 ÷ 見積もりの比の区間の最小以下、\n"
        "  # _mesh/c45/h_over_est_ratio_summary.json)。物理壁は CFD 前の積分法の壁 (prepare_ns と同じ構築)。流れ方向の密度の表 = 目標の間隔\n"
        f"  # min({W.AR_MARGIN}·ar_max·y1,表·cosθ_w, 雛形の分布 (throat_refine {REF_REFINE}・width {REF_WIDTH}) で ni {NI_REF} の間隔) の逆数の上包絡、\n"
        f"  # ni は壁際の AR ≤ {ARM:g} と間隔の上限を満たす最小 ({ni_min}) を 10 の倍数に切り上げ。壁法線の近壁層・17 桁の座標 (第一層 / 半径 ≈ 4e-7)\n"
        "  discretization: node\n"
        f"  ni: {ni}\n"
        f"  nj: {nj}\n"
        f"  wall_first_frac: {float(ft[-1, 1]):.6g}   # 表があれば使われない (記録・既定値の読み手のため出口の表の値)\n"
        f"  wall_first_frac_table: {fmt_table(ft, 6)}\n"
        f"  x_density_table: {fmt_table(dt, 5)}\n"
        f"  wall_normal_layer: [{WALL_NORMAL_LAYER[0]}, {WALL_NORMAL_LAYER[1]}]   # 近壁を壁法線に (局所半径比 d_n, d_b; M6 の H と同じ)\n"
        + ("" if radial_rec is None else
           f"  axis_cap_frac: {radial_rec['axis_cap_frac']:g}   # 径方向 (ユーザ 2026-10-15「上限3%」): 壁の第一セルから等比、間隔 / r_w が {radial_rec['axis_cap_frac']:g} で"
           f"軸まで一様。nj {nj} は全断面で隣接比 ≤ {radial_rec['radial_ratio_max']:g} の最小 (mesh2d._radial_fracs_capfixed、比の最大 "
           f"{radial_rec['ratio_max_all_stations']:.5f})\n")
        + (f"  msh_digits: {MSH_DIGITS}\n"
        "  throat_refine: 1.0\n"
        "  throat_width: 1.0\n"
        f"  ar_max: {int(ar_max)}\n"
        f"  max_nodes: {max_ns}\n"
        f"  max_nodes_reason: {json.dumps(reason_ns, ensure_ascii=False)}\n"))
    me_new = {k: v for k, v in me.items() if k not in ("ni", "max_nodes", "max_nodes_reason")}
    eu_lines = ["mesh_euler:\n",
                "  # ni は NS (mesh) と同じ (雛形の規則; plan tooling-nozzle-mesh-size-guard §4.1)。_mesh/hand_mesh.py build が書く\n",
                f"  ni: {ni}\n"]
    for k, v in me_new.items():
        eu_lines.append(f"  {k}: {v}\n")
    eu_lines += [f"  max_nodes: {max_eu}\n", f"  max_nodes_reason: {json.dumps(reason_eu, ensure_ascii=False)}\n"]
    i0, i1 = block_span("mesh")
    lines[i0:i1] = [mesh_txt]
    j0, j1 = block_span("mesh_euler")
    lines[j0:j1] = ["".join(eu_lines)]
    final = any(k in prob_dict.get("geometry", {}) for k in ("physical_wall_repr", "pw_upstream"))   # 最終の雛形 (Euler は写しで回す)
    use = (f"NS の問題 (Euler は {Path(a.out_yaml).stem}_euler.yaml)" if final else "Euler・NS の問題")
    hdr = (f"# case/71 M5.5 の {use} (plan tooling-nozzle-card-local-mesh {a.plan_item})。{Path(a.problem).name} (要求カードの生成器の出力、\n"
           "# card_m55.yaml の overrides mesh.ni 3000 は生成のためだけ) の mesh・mesh_euler を _mesh/hand_mesh.py build で置き換えたもの (手で編集しない)。\n"
           "# 他のキー・_card の記録 (mesh_estimate を含む) は生成器の出力のまま (_card.mesh_estimate は旧い見積もりで、この格子ではない)。\n")
    out = Path(a.out_yaml)
    # 一時のファイルに書いて実効の節点数を数え (格子の配列は作らない; mesh2d.node_count)、上限の内なら置き換える (codex 諮問 mesh-limits の
    # Major 2: 指定の ni × nj では Euler の axis_gap_frac の導出の点数を取り逃がす)。超えれば問題 YAML を書かずに止める
    tmp_out = out.with_name(out.stem + ".tmp" + out.suffix)
    tmp_out.write_text(hdr + "".join(lines))
    try:
        p2 = load_problem(tmp_out)
        mp2 = RA.mesh_params(p2, 1.0, ni=561, nj=97, wall_first_frac=4.5e-5)
        hint = ("ユーザの了承を取ってカードの mesh_limits を直し、承認の記録を作り直して --approval で渡す" if lim is not None else
                "ユーザの了承を取ってカードの mesh_limits に書き、承認の記録を --approval で渡す (自動で広げない)")
        nodes_ns = RA.check_node_budget(mp2, p2.mesh, "mesh", allow_hint=hint)
        nodes_eu = RA.check_node_budget(RA.mesh_params_euler(p2, 1.0), RA.mesh_euler_block(p2), "mesh_euler", allow_hint=hint)
    except ValueError as ex:
        tmp_out.unlink(missing_ok=True)
        raise RuntimeError(f"節点数が上限を超える — 問題 YAML ({out}) を書かずに止めた: {ex}") from ex
    import os
    os.replace(tmp_out, out)
    limits_rec = None if lim is None else {**lim, "nodes_ns": nodes_ns["n_nodes"], "nodes_euler": nodes_eu["n_nodes"]}
    # 読み直した表が作った表と同じ
    t_ft = np.asarray(p2.mesh["wall_first_frac_table"], float)
    t_dt = np.asarray(p2.mesh["x_density_table"], float)
    if not (np.array_equal(t_ft, ft) and np.array_equal(t_dt, dt)):
        raise RuntimeError("書いた YAML の表が作った表と一致しない")
    changed = sorted(k for k in set(old_mesh) | set(p2.mesh) if old_mesh.get(k) != p2.mesh.get(k))
    # ---- 6. 最終の雛形なら Euler 用の写し (生成器と同じ関数) と、生成器の写しとの照合 (mesh・mesh_euler 以外が同じ)
    eu_rec = None
    if final:
        eo = W.write_euler_companion(out)
        ref = Path(a.problem).with_name(Path(a.problem).stem + "_euler" + Path(a.problem).suffix)
        if eo is None or not ref.exists():
            raise RuntimeError(f"最終の雛形の問題なのに Euler の写しが書けない・生成器の写し {ref} が無い")
        d_eo, d_ref, d_out = (yaml.safe_load(open(f)) for f in (eo, ref, out))
        diff = sorted(k for k in set(d_eo) | set(d_ref) if k not in ("mesh", "mesh_euler") and d_eo.get(k) != d_ref.get(k))
        if diff or d_eo["mesh_euler"] != d_out["mesh_euler"] or d_eo["mesh"] != d_out["mesh"]:
            raise RuntimeError(f"Euler の写しが生成器の写しと mesh・mesh_euler 以外で違う ({diff}) か、NS の格子の節と違う")
        if any(k in d_eo["geometry"] for k in ("physical_wall_repr", "pw_upstream")):
            raise RuntimeError("Euler の写しに物理壁のキーが残っている")
        eu_rec = {"file": str(eo), "generator_copy": str(ref), "same_except_mesh": True,
                  "geometry_keys_dropped": sorted(set(d_out["geometry"]) - set(d_eo["geometry"]))}
        log(f"Euler の写し: {eo} (生成器の写し {ref.name} と mesh・mesh_euler 以外が同じ; 除いたキー {eu_rec['geometry_keys_dropped']})")
    rec = {"plan": f"tooling-nozzle-card-local-mesh {a.plan_item}", "problem_in": str(a.problem), "problem_out": str(out),
           "euler_companion": eu_rec, "radial": radial_rec,
           "r_t_m": rt, "x_U": x_U, "design_wall": [float(wd.x_in), float(wd.x_e)], "physical_wall": [xp0, xp1],
           "initializer": init_rec, "ar_max": ar_max, "ar_limit": ARM, "y1plus_target": W.Y1PLUS_TARGET, "s_c_rt": S_C,
           "formula_check_rel_max": rel, "coef": {"pipe": C_PIPE, "contr": C_CONTR, "throat": C_THROAT, "div": C_DIV, "blend": BLEND,
                                                  "check": c_chk},
           "first_table": {"n": len(ft), "added_points": added, "table_over_req_max_check": float(q.max()),
                           "table_over_req_min_check": float(q.min()), "table_over_req_max_build": float(q_b.max()),
                           "check_dx": CHECK_DX, "build_dx": BUILD_DX},
           "density_table": {"n": len(dt), "h_cap_ref": {"ni": NI_REF, "throat_refine": REF_REFINE, "throat_width": REF_WIDTH,
                                                         "min": float(h_cap.min()), "max": float(h_cap.max())},
                             "h_AR_min": float(h_AR.min()), "x_h_AR_min": float(xg[np.argmin(h_AR)])},
           "ni_min": ni_min, "ni": ni, "ni_ideal": ni_ideal, "nj": nj, "nodes": n_nodes, "wall_ar_max_stations": float(ar.max()),
           "dx_min": float(dxs.min()), "dx_max": float(dxs.max()), "dx_ratio_max": float(max(qd[kq], 1 / qd[kq])),
           "x_dx_ratio_max": float(xs[kq + 1]), "nodes_check": {"mesh": nodes_ns, "mesh_euler": nodes_eu}, "node_limits": limits_rec,
           "mesh_keys_changed": changed, "y_sc_build_m": B["y_sc"], "y_sc_2pt_m": R.y1_cap_raw, "time_s": time.time() - t0}
    Path(a.record).write_text(json.dumps(W._json_safe(rec), ensure_ascii=False, indent=1))
    # 断面ごとの値 (作る刻みの 10 点おき)
    with open(Path(a.record).with_suffix(".csv"), "w") as fh:
        fh.write("x_rt,r_design_rt,r_physical_rt,c,y1_req_raw_m,y1_req_capped_m,f_req,f_table,table_over_req\n")
        fB = interp_f(ft, B["x"])
        for i in list(range(0, len(B["x"]), 100)) + [len(B["x"]) - 1]:
            fh.write(f"{B['x'][i]:.6f},{B['rd'][i]:.8f},{B['rp'][i]:.8f},{B['c'][i]:.5f},{B['y1_raw'][i]:.6e},{B['y1'][i]:.6e},"
                     f"{B['f'][i]:.6e},{fB[i]:.6e},{fB[i] / B['f'][i]:.6f}\n")
    log(f"書いた: {out}  記録 {a.record}  ({time.time() - t0:.0f} s)")


def _qc_and_relax(rd: Path, p, mp, wall, cfg_qc: str, bcond: str, ar_max: float, digits: int, nodes, block):
    """prepare / prepare_ns と同じ: 格子を作り msh を書き、cell 変換の写しで check_mesh_quality と mesh_ar_relax_check (rd は作成済み)。"""
    from forge_design.evaluate.runner import FORGE_TOOLS, _ENV, converter_path, mesh_ar_relax_check, run_converter
    coords, quads, bedges = M2.generate_axisym_mesh(wall, mp)
    RA._check_generated_nodes(coords, nodes, block)
    M2.write_msh41_2d(rd / "nozzle.msh", coords, quads, bedges, digits=digits)
    (rd / "bcondConfig.yaml").write_text(bcond)
    (rd / "probe.yaml").write_text(RA.PROBE_STUB)
    (rd / "solverConfig.yaml").write_text(cfg_qc)
    run_converter([str(converter_path()), "nozzle.msh", "nozzle_qc.h5"], cwd=rd, env=_ENV)
    q = subprocess.run([sys.executable, str(FORGE_TOOLS / "check_mesh_quality.py"), "nozzle_qc.h5", "--ar-max", str(int(ar_max))],
                       cwd=rd, env=_ENV, capture_output=True, text=True)
    (rd / "MESH_QUALITY.txt").write_text("# cell 変換コピーで検査 (品質は primal の性質)\n" + q.stdout + q.stderr)
    relax = mesh_ar_relax_check(rd / "nozzle_qc.h5", coords, bedges, float(mp.scale), ar_max)
    (rd / "MESH_AR_RELAX.json").write_text(json.dumps(relax, indent=1))
    verdict = next((ln for ln in reversed(q.stdout.splitlines()) if ln.startswith("VERDICT")), None)
    return coords, quads, bedges, q.returncode, verdict, relax


def check(a):
    t0 = time.time()
    from forge_design.evaluate.runner import converter_path
    card = W.load_card(a.card)
    p, d, wd, wp, init_rec = walls(a.problem)
    rt = float(p.spec["r_throat"])
    x_U = -float(wd.up.L_U)
    log = lambda *s: print(*s, flush=True)  # noqa: E731
    out = {"problem": str(a.problem), "converter": str(converter_path()), "initializer": init_rec, "r_t_m": rt}
    out["converter_sha256"] = subprocess.run(["sha256sum", str(converter_path())], capture_output=True, text=True).stdout.split()[0]
    # ---- NS (prepare_ns の格子の経路; forge は起動しない)
    nodes = RA.check_node_budget(RA.mesh_params(p, 1.0, ni=561, nj=97, wall_first_frac=4.5e-5), p.mesh, "mesh")
    mp = RA.mesh_params(p, rt, ni=561, nj=97, wall_first_frac=4.5e-5)
    n, out_int, _ = RA._ns_steps(p, None)
    rd = Path(a.run_dir) / "ns"
    ar_max = float(p.mesh.get("ar_max", 1000))
    rd.mkdir(parents=True, exist_ok=False)
    cfg_qc = RA._apply_gas_to_config(RA._config_euler(p, n, out_int, 4.0, 1), p, rd)     # 種のファイルを rd に書く (prepare_ns と同じ)
    bc = RA._bcond_with_species(RA._bcond(p, euler=False, profile_ok=True), RA._tp_species_Y(p))
    coords, quads, bedges, rc, verdict, relax = _qc_and_relax(rd, p, mp, wp, cfg_qc, bc, ar_max, int(p.mesh.get("msh_digits", 10)),
                                                              nodes, "mesh")
    log(f"NS: 節点 {coords.shape[0]:,}  品質 rc {rc}  {verdict}")
    log(f"NS: AR 緩和 違反 {relax['n_violations']} (AR > {relax['ar_threshold']:.0f} の {relax['n_high_ar']} セル、スキュー最大 "
        f"{relax['max_skew_high_ar']})  入口面 {'平面' if relax['inlet_planar'] else '非平面'}・出口面 {'平面' if relax['outlet_planar'] else '非平面'} "
        f"(出口の x の広がり {relax['outlet_x_spread_m']:.3e} m)")
    # node 変換 (prepare_ns と同じ config; wall_dist は no-slip 壁で作られる)
    from forge_design.evaluate.runner import _ENV, run_converter
    main_line = RA._main_line_spec(p, None, None)
    cfg_ns, _cfg_main = RA._ns_configs(p, rd, n, out_int, None, None, main_line)
    (rd / "solverConfig.yaml").write_text(cfg_ns)
    run_converter([str(converter_path()), "nozzle.msh", "nozzle.h5"], cwd=rd, env=_ENV)
    import h5py
    with h5py.File(rd / "nozzle.h5", "r") as f:
        C3 = np.asarray(f["MESH/COORD"][:], dtype=float).reshape(-1, 3)
        cdtype = str(f["MESH/COORD"].dtype)
    ni, nj = mp.ni, coords.shape[0] // mp.ni
    same_order = C3.shape == coords.shape and float(np.max(np.abs(C3[:, :2] - coords[:, :2]))) <= 1e-15 * float(np.abs(coords).max())
    diff_max = float(np.max(np.abs(C3[:, :2] - coords[:, :2]))) if C3.shape == coords.shape else float("nan")
    log(f"変換後の座標: dtype {cdtype}  msh と同じ順序・値 {same_order} (最大差 {diff_max:.3e} m)")
    if not same_order:
        raise RuntimeError("変換後の座標の順序・値が msh と違う (第一層の照合の前提が崩れる)")
    X = C3[:, 0].reshape(ni, nj) / rt
    Rr = C3[:, 1].reshape(ni, nj) / rt
    xw, rw_ = X[:, -1], Rr[:, -1]
    thick = np.hypot(X[:, -1] - X[:, -2], Rr[:, -1] - Rr[:, -2]) * rt       # 壁法線の層: 壁点と内側の点の距離 = 法線の第一層厚 [m]
    ft = np.asarray(p.mesh["wall_first_frac_table"], float)
    xs = xw
    expect = interp_f(ft, xs) * wp.r(xs) * rt
    # 要求 c·y1,req (設計壁の式、格子と同じ刻みの式で; 頭打ちあり) を壁点の x で
    R = Req(card, wd, wp, rt, x_U)
    G = R.on(CHECK_DX)
    req = np.interp(xs, G["x"], G["c"] * G["y1"])
    y1raw = np.interp(xs, G["x"], np.where(np.isfinite(G["y1_raw"]), G["y1_raw"], G["y1"]))
    ratio = thick / req
    y1p_est = W.Y1PLUS_TARGET * thick[1:] / y1raw[1:]          # 見積もりの y1+ (補正なし; 入口の節点を除く)
    y1p_cal = W.Y1PLUS_TARGET * thick / req                    # M6 の比で補正した要求に対する y1+ (≤ 0.8 なら要求以下)
    k1 = int(np.argmax(ratio))
    kp = int(np.argmax(y1p_est))
    dxs = np.diff(X[:, 0])
    qd = dxs[1:] / dxs[:-1]
    kq = int(np.argmax(np.maximum(qd, 1 / qd)))
    chord = np.hypot(np.diff(xw), np.diff(rw_))
    qc = chord[1:] / chord[:-1]
    segs = {"pipe": xs <= x_U, "contr": (xs >= x_U) & (xs <= -X_THROAT_HALF), "throat": np.abs(xs) <= X_THROAT_HALF, "div": xs >= X_THROAT_HALF}
    # 径方向の配点 (変換後の座標): 軸側のセルの高さ / 局所半径、j の線に沿った隣接間隔の比 (断面ごとの最大) を区間ごとに
    dsj = np.hypot(np.diff(X, axis=1), np.diff(Rr, axis=1))           # (ni, nj − 1)、軸 → 壁
    qj = dsj[:, 1:] / dsj[:, :-1]
    q_st = np.max(np.maximum(qj, 1.0 / qj), axis=1)
    h_ax = (Rr[:, 1] - Rr[:, 0]) / Rr[:, -1]
    x_half = 0.5 * float(wp.x_e)
    rsegs = {"pipe": xs <= x_U, "contr": (xs > x_U) & (xs < -X_THROAT_HALF), "throat": np.abs(xs) <= X_THROAT_HALF,
             "div_front": (xs > X_THROAT_HALF) & (xs <= x_half), "div_rear": xs > x_half}
    radial = {k: {"stations": int(mk.sum()), "axis_cell_over_rw_min": float(h_ax[mk].min()), "axis_cell_over_rw_max": float(h_ax[mk].max()),
                  "ratio_max_min": float(q_st[mk].min()), "ratio_max_max": float(q_st[mk].max())} for k, mk in rsegs.items()}
    radial["all"] = {"stations": len(xs), "axis_cell_over_rw_min": float(h_ax.min()), "axis_cell_over_rw_max": float(h_ax.max()),
                     "ratio_max_min": float(q_st.min()), "ratio_max_max": float(q_st.max()), "x_ratio_max": float(xs[int(np.argmax(q_st))]),
                     "axis_cap_frac": p.mesh.get("axis_cap_frac")}
    by_seg = {}
    for k, mk in segs.items():
        mk1 = mk.copy()
        mk1[0] = False
        by_seg[k] = {"stations": int(mk.sum()), "thick_over_req_max": float(ratio[mk].max()), "thick_um_min": float(thick[mk].min() * 1e6),
                     "thick_um_max": float(thick[mk].max() * 1e6), "y1plus_est_max": float((W.Y1PLUS_TARGET * thick[mk1] / y1raw[mk1]).max()),
                     "y1plus_cal_max": float(y1p_cal[mk].max())}
    out["ns"] = {"run_dir": str(rd), "n_nodes": int(coords.shape[0]), "ni": ni, "nj": nj, "nodes_budget": nodes,
                 "quality_rc": rc, "quality_verdict": verdict, "ar_relax": relax,
                 "coord_dtype": cdtype, "coord_same_as_msh": bool(same_order), "coord_max_diff_m": diff_max,
                 "first_layer": {"thick_over_req_max": float(ratio[k1]), "x_at_max": float(xs[k1]),
                                 "thick_over_table_minus1_absmax": float(np.max(np.abs(thick / expect - 1.0))),
                                 "thick_um_min": float(thick.min() * 1e6), "thick_um_max": float(thick.max() * 1e6),
                                 "n_over_req": int((ratio > 1.0).sum())},
                 "y1plus_est_max": float(y1p_est[kp]), "x_y1plus_est_max": float(xs[1:][kp]),
                 "y1plus_cal_max": float(y1p_cal.max()), "x_y1plus_cal_max": float(xs[int(np.argmax(y1p_cal))]),
                 "by_segment": by_seg, "radial": radial,
                 "dx_ratio_max": float(max(qd[kq], 1 / qd[kq])), "x_dx_ratio_max": float(X[kq + 1, 0]),
                 "wall_chord_ratio_max": float(np.max(np.maximum(qc, 1 / qc))),
                 "dx_min_rt": float(dxs.min()), "dx_max_rt": float(dxs.max())}
    log(f"NS: 第一層 実寸/(c·y1,req) の最大 {ratio[k1]:.6f} @ x {xs[k1]:.3f} (超える断面 {(ratio > 1).sum()})  実寸/(表×r_物理壁) − 1 の最大 "
        f"{out['ns']['first_layer']['thick_over_table_minus1_absmax']:.2e}  厚さ {thick.min() * 1e6:.4f}–{thick.max() * 1e6:.3f} µm")
    log(f"NS: 見積もりの y1+ の最大 {y1p_est[kp]:.4f} @ x {xs[1:][kp]:.3f}  補正した要求に対する y1+ の最大 {y1p_cal.max():.4f}  "
        f"dx 比の最大 {out['ns']['dx_ratio_max']:.4f} @ x {out['ns']['x_dx_ratio_max']:.3f}  壁の弦の比の最大 {out['ns']['wall_chord_ratio_max']:.4f}")
    for k, v in by_seg.items():
        log(f"  {k:<7s} 断面 {v['stations']:5d}  実寸/(c·y1,req) 最大 {v['thick_over_req_max']:.5f}  厚さ {v['thick_um_min']:.4f}–{v['thick_um_max']:.3f} µm  "
            f"見積もりの y1+ 最大 {v['y1plus_est_max']:.4f}  補正 y1+ 最大 {v['y1plus_cal_max']:.4f}")
    log(f"NS: 径方向 (axis_cap_frac {p.mesh.get('axis_cap_frac')}、nj {nj}): 断面ごとの隣接比の最大の範囲・軸側のセルの高さ / 局所半径")
    for k, v in radial.items():
        log(f"  {k:<9s} 断面 {v['stations']:5d}  隣接比の最大 {v['ratio_max_min']:.4f}–{v['ratio_max_max']:.4f}  "
            f"軸側のセル / r_w {v['axis_cell_over_rw_min']:.5f}–{v['axis_cell_over_rw_max']:.5f}")
    np.savetxt(rd / "first_layer_check.csv", np.c_[xs, wp.r(xs), thick, req, y1raw, ratio], delimiter=",",
               header="x_rt,r_physical_rt,thick_m,c_y1req_m,y1req_raw_m,thick_over_req", comments="")
    # ---- Euler (prepare の格子の経路)
    m_eu = RA.mesh_euler_block(p)
    nodes_e = RA.check_node_budget(RA.mesh_params_euler(p, 1.0), m_eu, "mesh_euler")
    we, _ew = RA.euler_reference_wall(p, d, rt, pwu=RA._pw_upstream(p.geometry))
    mpe = RA.mesh_params_euler(p, rt)
    n_e = int(p.evaluate.get("nStepOuter", 12000))
    oi_e = int(p.evaluate.get("outStepInterval", max(n_e // 3, 1)))
    rde = Path(a.run_dir) / "euler"
    rde.mkdir(parents=True, exist_ok=False)
    cfg_e = RA._apply_gas_to_config(RA._config_euler(p, n_e, oi_e, 4.0, 1), p, rde)
    ce, qe, be, rce, ve, relax_e = _qc_and_relax(rde, p, mpe, we, cfg_e, RA._bcond_with_species(RA._bcond(p, euler=True), RA._tp_species_Y(p)),
                                                 float(m_eu.get("ar_max", 1000)), 10, nodes_e, "mesh_euler")
    xe = ce[:, 0].reshape(mpe.ni, -1)[:, 0] / rt
    dxe = np.diff(xe)
    qe_ = dxe[1:] / dxe[:-1]
    out["euler"] = {"run_dir": str(rde), "n_nodes": int(ce.shape[0]), "ni": mpe.ni, "nj": int(ce.shape[0] // mpe.ni), "nodes_budget": nodes_e,
                    "quality_rc": rce, "quality_verdict": ve, "ar_relax": relax_e,
                    "dx_ratio_max": float(np.max(np.maximum(qe_, 1 / qe_))), "dx_min_rt": float(dxe.min()), "dx_max_rt": float(dxe.max())}
    log(f"Euler: 節点 {ce.shape[0]:,}  品質 rc {rce}  {ve}  dx {dxe.min():.5f}–{dxe.max():.4f} r_t  dx 比最大 {out['euler']['dx_ratio_max']:.4f}")
    out["time_s"] = time.time() - t0
    Path(a.record).write_text(json.dumps(W._json_safe(out), ensure_ascii=False, indent=1))
    log(f"書いた: {a.record} ({time.time() - t0:.0f} s)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=("build", "check"))
    ap.add_argument("--design", default="~/forge-cardmesh/design")
    ap.add_argument("--card", required=True)
    ap.add_argument("--problem", required=True)
    ap.add_argument("--ratio-summary", default=str(Path(__file__).resolve().parent / "c45" / "h_over_est_ratio_summary.json"))
    ap.add_argument("--out-yaml")
    ap.add_argument("--record", required=True)
    ap.add_argument("--run-dir")
    ap.add_argument("--axis-cap-frac", type=float, default=None,
                    help="NS の径方向を mesh2d の axis_cap_frac にする (#3b、ユーザ 2026-10-15「上限3%」は 0.03)。nj は --radial-ratio-max で決める")
    ap.add_argument("--radial-ratio-max", type=float, default=1.2, help="--axis-cap-frac の nj を決める全断面の隣接比の上限")
    ap.add_argument("--plan-item", default="§5.1 #3a", help="YAML の注記と記録に書く plan の項目 (#3a の再現は既定のまま)")
    ap.add_argument("--approval", default=None,
                    help="build: 承認の記録 card_approval.json。節点数の上限をカードの承認済みの mesh_limits と照合して書く (plan "
                         "tooling-design-orchestrator §4.6・§5.1 #8)。無ければ入力の問題 YAML の上限 (既定 30 万) を超えたら止める")
    ap.add_argument("--sanity", type=float, default=None,
                    help="見積もり (plan §3 の (C)) の再現の確認: 全区間の係数をこの値、半径は設計壁で ni だけを出す (YAML は書かない)")
    a = ap.parse_args()
    _imports(a.design)
    if a.sanity is not None:
        global C_PIPE, C_CONTR, C_THROAT, C_DIV, NORM
        C_PIPE = C_CONTR = C_THROAT = C_DIV = float(a.sanity)
        NORM = "design"
    if a.mode == "build":
        build(a)
    else:
        check(a)


if __name__ == "__main__":
    main()
