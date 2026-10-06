"""plan tooling-nozzle-throat-monotone-r2 §5.1 #10 (b) / §3 仮説 H: MOC の自己整合 (再出発) 試験 (CFD 0 step)。

目的: スロート始点直後の壁角差 Δθ = θ_wall − atan(x/R) の「x → 0 に残る約 +0.04° の成分」が、MOC の始点処理で
生じるのか、入力 (初期線・m*・アンカー) の不整合で生じるのかを分ける材料を取る。
MOC 解の内部の C⁻ (網の列 = 軸節点から壁へ上る C⁻) を初期線として MOC をやり直し、新しい始点の直後に元の壁との差が出るかを見る。

生産コード (design/forge_design/) は変えない。`design_chain` は生産のまま回し、`moc_inverse._design_cplus` と
`runner_axismach.inverse_design` を一時的に包んで MOC 網 (levels)・初期前線・m*・throat・target を捕まえるだけ。
再出発は生産の部品 (`InverseMOC.fill_levels`・`_flux_along`・`cplus_flux_wall`・`inverse_design`) をそのまま使う。
捕まえた前線を自前の手順で回し直した壁が生産の壁とビット一致することを各解で確かめる (selfcheck)。

MOC の始点処理 (コードを読んだ結果。cplus 経路):
  - 初期前線 = 軸節点 (下流 → 上流、x0 は含めない) + 初期線 (軸 → 壁、r 等間隔 n_start 点)。初期線の足は壁の始点そのもの。
  - fill_levels: L_k[i] = 単位過程 (B = L_{k−1}[i] が C⁻ 担体、A = L_{k−1}[i+1] が C⁺ 担体)。列 = C⁻、反対角線 = C⁺。
    初期線内の隣接対は A が B の C⁻ 上にあるので P ≈ A に縮退し (縮退カスケード)、初期線の写しが段ごとに作り直される。
  - 壁: 各 C⁺ に沿って累積流束 (起点までの初期線の流束から) が m* に達する点。足の C⁺ は起点で m* なので足自身が壁の第 0 点。
  - m* = 初期線を横切る流束 (_flux_along) そのもの。外から与えない。
  - 軸端: x0・M は初期線の軸端の値。アンカー M_A はその同じ値 (CFD ピン)、M′_A は軸 evenfit の窓フィット、M″_A は Hall の式。

単位過程の差し替え (この試験用。生産は prod):
  prod  生産のまま (予測子-修正子 n_corr 2、軸上の端点の sinθ/r は相手の端点の値で代用。相手も軸上なら 0)。
  conv  修正子だけを収束させる (n_corr 20)。軸の 1 段目 (軸節点 2 個から作る点) で修正子が振動 (比 ≈ −1/2) し、
        n_corr 2 では θ が収束値から ~0.02〜0.05° ずれるため。
  K2c   conv + 軸上の端点の sinθ/r を常に 0 (相手を借りない)。軸節点 2 個の対では生産と同じで、軸節点と軸外の点の対
        (再出発の初期線の軸端) でも同じ扱いになる → 網の列を初期線に戻したとき単位過程が恒等になる (自己整合)。
  prod で V1 が網を再現しないのは、列の 1 段目を「軸節点 2 個の対」で作ったのに、再出発では「軸節点と列の点の対」で
  作り直すため (代用が対の組み方で変わる + 修正子の非収束)。K2c はその 2 つを除いた試験。

再出発の変種 (列 i、壁との交点 W):
  V1n  列 i の網の節点 (軸 → 壁の直下の節点) をそのまま初期線にする。足 = 壁の直下の節点。
       比較相手 = 元の網で同じ節点を通る流線 (cplus_flux_wall の閾値をその節点の累積流束にしたもの)。
  V1w  V1n の線に、列 i と元の壁の交点 W (列の線分上で線形補間) を足す。比較相手 = 元の壁。
  V2   V1w の線を r 等間隔 n 点 (n = 元の n_start) に線形補間し直す (生産の throat_characteristic と同じ扱い)。軸節点は元と同じ。
  V3   V2 の線 + 軸格子を新しい足から作り直す (`inverse_design` をそのまま呼ぶ。生産の始点と同じ組み立て)。
軸節点が同じ V1n・V1w では、列 < i の網の節点が元と一致するか (net_diff、W から出る C⁺ 上の点は除く) も見る。

補助 (§5.1 #10 の残り): 初期線の軸端の (x, M) とアンカー (x_A, M_A)、m* と各種の流量積分の相対差、
初期線の縮退カスケード (L_j[n_ax] = 初期線の j 番目の点を j 回単位過程に通し直した写しと元の点の差)、
軸の 1 段目の点の修正子の非収束量と、列の対を単位過程に通し直したときの差。

usage: [CASE_RUNS=<run_0062 のある case dir>] python3 throat_moc_restart_test.py [OUT_JSON]
       既定の CASE_RUNS はローカル主ツリー /home/sano/work/forge/case/45.isobutane_m6_d155 (throat_start_limit.py 等と同じ参照)。
       出力: _band_ab/throat_moc_restart_test.json (所要 約 35 分 [2026-10-06 実測 2085 s]、メモリ 約 3 GB)
"""
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

C = Path(__file__).resolve().parent
ROOT = C.parents[1]
sys.path.insert(0, str(ROOT / "design"))
from forge_design.evaluate import runner_axismach as RA  # noqa: E402
from forge_design.geometry import moc_inverse as MI  # noqa: E402
from forge_design.geometry import moc_kernel as MK  # noqa: E402
from forge_design.geometry.moc_kernel import _Pt, pm_mach_vec, pm_nu  # noqa: E402

RUNS = Path(os.environ.get("CASE_RUNS", "/home/sano/work/forge/case/45.isobutane_m6_d155"))
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else C / "_band_ab/throat_moc_restart_test.json"
PROBS = {"pin": "problem_d155_ns_finemesh_recal_final.yaml",      # 最終問題 (CFD ピン初期線)
         "hall": "problem_d155_euler_c2final_n2400.yaml"}          # Hall 初期線 (同じ形・同じ軸格子設定)
# 分解能 (生産は n_start 41 / n_axis_inv 2400 / axis_dx0 0.03)
LV = {"L1": dict(n_start=161, n_axis=2400, dx0=0.03), "L2": dict(n_start=321, n_axis=2400, dx0=0.03),
      "L3": dict(n_start=321, n_axis=4800, dx0=0.03), "L4": dict(n_start=321, n_axis=4800, dx0=0.015),
      "L5": dict(n_start=321, n_axis=4800, dx0=0.0075)}
# (問題, 単位過程, 分解能, 再出発をするか)。L3・L5 は元の解の当てはめだけ (軸間隔の感度)
SPEC = ([("pin", "prod", L, L in ("L1", "L2", "L4")) for L in LV]
        + [("pin", "conv", "L1", True), ("pin", "conv", "L4", False)]
        + [("pin", "K2c", L, L in ("L1", "L2", "L4")) for L in LV]
        + [("hall", "prod", L, L in ("L1", "L2", "L4")) for L in LV]
        + [("hall", "K2c", L, L in ("L1", "L2")) for L in ("L1", "L2", "L4")])
TARGETS = (0.01, 0.05, 0.12, 0.3)    # 再出発の足 (列と元の壁の交点) の x の目標 [r_t]。0.01 は始点直近の補助
FIT = (0.01, 0.04)                   # 試験 (c) と同じ当てはめ窓 (始点からの距離)
N_CORR_CONV = 20                     # conv・K2c の修正子回数 (軸の 1 段目の振動が 1e-7° 以下に収まる)
deg = np.degrees

# ---- 単位過程の差し替え (この試験用。with kernel(...) を抜けると必ず生産に戻す) --------------------------------
_SOR, _IV = MK._sin_over_r_vec, MI.interior_vec


def _sor_axis_zero(r_p, th_p, r_o, th_o):
    """K2c: 軸上 (生産と同じ判定) の端点の sinθ/r は相手から借りずに 0。軸外の点は自身の値 (生産と同じ)。"""
    on_axis = (r_p <= 1e-9) | (r_p < MK.AXIS_LIMIT_FRAC * r_o)
    ok = ~on_axis & (r_p > 1e-9)
    return np.where(ok, np.sin(th_p) / np.where(ok, r_p, 1.0), 0.0)


def _iv_conv(*a):
    return MK.interior_vec(*a[:-1], N_CORR_CONV)          # 最後の引数 (n_corr) だけ差し替え


class kernel:
    def __init__(self, name):
        if name not in ("prod", "conv", "K2c"):
            raise ValueError(name)
        self.name = name

    def __enter__(self):
        if self.name in ("conv", "K2c"):
            MI.interior_vec = _iv_conv
        if self.name == "K2c":
            MK._sin_over_r_vec = _sor_axis_zero
        return self

    def __exit__(self, *exc):
        MK._sin_over_r_vec, MI.interior_vec = _SOR, _IV
        return False


# ---- 生産の設計チェーンから網を捕まえる -----------------------------------------------------------------
def solve(prob, n_start, n_axis, dx0):
    """生産の design_chain を回し、MOC 網・初期前線・m*・throat・target を捕まえる (単位過程は呼び出し側の kernel)。"""
    p = RA.load_problem(C / prob)
    p.geometry.update(n_start=int(n_start), n_axis_inv=int(n_axis), axis_dx0=float(dx0))
    if p.geometry.get("initial_line") == "cfd":
        p.geometry["initial_line_run"] = str(RUNS / p.geometry["initial_line_run"])
    cap = {"n_calls": 0}
    o_dc, o_inv = MI._design_cplus, RA.inverse_design

    def dc(inv, init, n_ax, ax, mdot_star, g, *a, **k):
        out = o_dc(inv, init, n_ax, ax, mdot_star, g, *a, **k)
        cap.update(inv=inv, init=init, n_ax=int(n_ax), ax=np.asarray(ax, float), mstar=float(mdot_star), g=g,
                   lev=out["levels"], wall=out["wall_full"])
        return out

    def invd(throat, target, **kw):
        cap["n_calls"] += 1
        cap.update(throat=throat, target=target, kw=dict(kw))
        return o_inv(throat, target, **kw)
    MI._design_cplus, RA.inverse_design = dc, invd
    try:
        d, err = RA.design_chain(p), None
    except Exception as e:          # 壁 QA・当てはめ等の下流ゲートで落ちても MOC 網は捕まえてある (記録して続ける)
        d, err = None, f"{type(e).__name__}: {str(e)[:300]}"
    finally:
        MI._design_cplus, RA.inverse_design = o_dc, o_inv
    if "lev" not in cap or cap["n_calls"] != 1:
        raise RuntimeError(f"MOC 網を 1 回だけ捕まえられなかった (呼び出し {cap['n_calls']} 回, err={err})")
    cap.update(d=d, err=err, R=float(p.geometry.get("R", 2.0)))
    return cap


def run_net(inv, init, n_ax, g):
    """`_design_cplus` と同じ手順で網と壁を作る (exit 処理は除く)。戻り: levels, cum0, m*, wall_full。"""
    lev = inv.fill_levels(init)
    cum0 = np.zeros(len(init))
    cum0[n_ax:] = MI._flux_along(init[n_ax:], g)
    ms = float(MI._flux_along(init[n_ax:], g)[-1])   # inverse_design の mdot_star と同じ式
    return lev, cum0, ms, MI.cplus_flux_wall(lev, cum0, ms, g)


def pts_exact(rows, g):
    """網の節点 [x, r, θ, ν, M] → _Pt (ν・M はそのまま)。"""
    return [_Pt(float(a[0]), float(a[1]), float(a[2]), float(a[3]), g, float(a[4])) for a in rows]


def pts_like_production(x, r, M, th, g):
    """inverse_design と同じ組み立て (M → ν = pm_nu(M)、M は ν から引き直す)。"""
    return [_Pt(float(x[q]), float(r[q]), float(th[q]), float(pm_nu(float(M[q]), g)), g) for q in range(len(x))]


def column(lev, i):
    """列 i (軸節点 L_0[i] から上る C⁻) の節点列。最初の欠損の手前まで。"""
    col = lev[:, i]
    ok = np.isfinite(col[:, 0]) & np.isfinite(col[:, 1])
    cut = int(np.argmin(ok)) if not ok.all() else len(ok)
    return np.array(col[:cut])


def crossing(col, wall, g):
    """列と元の壁 (x の折れ線) の交点。戻り: (k*, W[x,r,θ,ν,M]) — k* は壁以上に出た最初の節点。"""
    rw = np.interp(col[:, 0], wall[:, 0], wall[:, 1])
    above = col[:, 1] - rw >= 0.0
    above[0] = False
    if not above.any():
        return None, None
    k = int(np.argmax(above))
    a, b = col[k - 1], col[k]

    def f(s):
        return (a[1] + s * (b[1] - a[1])) - np.interp(a[0] + s * (b[0] - a[0]), wall[:, 0], wall[:, 1])
    s = brentq(f, 0.0, 1.0, xtol=1e-15, rtol=1e-15)
    W = a + s * (b - a)
    W[4] = float(pm_mach_vec(np.array([W[3]]), g)[0])        # M は補間した ν から
    return k, W


def flux_to_node(lev, cum0, m, k, g):
    """元の網で、C⁺ (起点 L_0[m]) の k 番目の節点までの累積流束 (cplus_flux_wall と同じ式・同じ演算順)。"""
    line = MI.cplus_lines(lev, m)
    if len(line) <= k:
        raise RuntimeError(f"C⁺ {m} が {k} 点に届かない")
    x, r, th, M = line[:, 0], line[:, 1], line[:, 2], line[:, 4]
    Mm, thm, rm = 0.5 * (M[1:] + M[:-1]), 0.5 * (th[1:] + th[:-1]), 0.5 * (r[1:] + r[:-1])
    d = (MI._mass_flux_density(Mm, g) * (np.cos(thm) * np.diff(r) - np.sin(thm) * np.diff(x)) * 2.0 * np.pi * rm)
    cum = cum0[m] + np.concatenate([[0.0], np.cumsum(d)])
    return float(cum[k]), line[k]


def cascade(lev, n_ax, n_line):
    """初期線内の隣接対の縮退カスケード: L_j[n_ax] (= 線の j 番目の点を j 回単位過程に通し直した写し) と元の点の差。
    網の列 n_ax−1 はこの写しを C⁺ 担体として使う。入力が網の離散 C⁻ 関係と整合なら差は 0。θ・ν は度。"""
    j = np.arange(n_line)
    dg, l0 = np.array(lev[j, n_ax]), np.array(lev[0, n_ax + j])
    ok = np.isfinite(dg[:, 0])
    dpos = np.hypot(dg[:, 0] - l0[:, 0], dg[:, 1] - l0[:, 1])
    dth, dnu = deg(dg[:, 2] - l0[:, 2]), deg(dg[:, 3] - l0[:, 3])
    rt = l0[:, 1] / l0[-1, 1]
    prof = {f"{f:g}": {"dth_deg": float(np.interp(f, rt[ok], dth[ok])), "dnu_deg": float(np.interp(f, rt[ok], dnu[ok]))}
            for f in (0.02, 0.1, 0.25, 0.5, 0.75, 0.9, 0.97)}
    top = [{"r": float(l0[q, 1]), "dpos": float(dpos[q]), "dth_deg": float(dth[q]), "dnu_deg": float(dnu[q])}
           for q in range(max(n_line - 3, 0), n_line) if ok[q]]
    q_max = int(np.nanargmax(np.where(ok, np.abs(dth), -1.0)))
    return {"n_line": int(n_line), "n_finite": int(ok.sum()), "max_dpos": float(np.nanmax(dpos[ok])),
            "max_abs_dth_deg": float(np.abs(dth[q_max])), "r_at_max_abs_dth": float(l0[q_max, 1]),
            "max_abs_dnu_deg": float(np.nanmax(np.abs(dnu[ok]))), "profile_r_over_rtop": prof, "top3": top}


def first_row(cap, cols):
    """軸の 1 段目 (列の最初の点 = 軸節点 2 個の対から作る点) の検査。θ は度。
    nc2_minus_nc40: 修正子 2 回と 40 回の差 (いまの sinθ/r の扱いで)。reproc_dth_k1: 列の対 (軸節点, 1 段目の点) を
    いまの単位過程に通し直したときの 1 段目の点の差 (自己整合なら 0)。reproc_dth_k5plus_max: 2 段目以降の同じ量の最大。"""
    lev, g = cap["lev"], cap["g"]
    out = {}
    for i in cols:
        col = column(lev, i)
        A, B = lev[0, i + 1], lev[0, i]
        r2 = MK.interior_vec(*A[[0, 1, 2, 3, 4]][:, None], *B[[0, 1, 2, 3, 4]][:, None], g, 1.0, 2)
        r40 = MK.interior_vec(*A[[0, 1, 2, 3, 4]][:, None], *B[[0, 1, 2, 3, 4]][:, None], g, 1.0, 40)
        Aq, Bq = col[1:], col[:-1]
        q = MI.interior_vec(Aq[:, 0], Aq[:, 1], Aq[:, 2], Aq[:, 3], Aq[:, 4], Bq[:, 0], Bq[:, 1], Bq[:, 2], Bq[:, 3], Bq[:, 4],
                            g, 1.0, 2)
        dq = deg(q[2] - col[1:, 2])
        out[str(i)] = {"axis_x": float(col[0, 0]), "axis_dx": float(B[0] - A[0]), "r1": float(col[1, 1]),
                       "theta1_deg": float(deg(col[1, 2])),
                       "nc2_minus_nc40_dth_deg": float(deg(r2[2][0] - r40[2][0])),
                       "nc2_minus_nc40_dnu_deg": float(deg(r2[3][0] - r40[3][0])),
                       "reproc_dth_k1_deg": float(dq[0]), "reproc_dth_k5plus_max_deg": float(np.abs(dq[4:]).max())}
    return out


def compare(new, ref, x_f, R):
    """新しい壁 (x_f より下流の点) と比較相手の壁 (x で線形補間) の差。θ は度、r は r_t。"""
    m = (new[:, 0] > x_f + 1e-12) & (new[:, 0] <= ref[-1, 0])
    x, r, th = new[m, 0], new[m, 1], new[m, 2]
    rr, tr = np.interp(x, ref[:, 0], ref[:, 1]), np.interp(x, ref[:, 0], ref[:, 2])
    dth, dr, xi = deg(th - tr), r - rr, x - x_f
    first = [{"x": float(x[q]), "x_minus_xf": float(xi[q]), "dth_deg": float(dth[q]), "dr": float(dr[q]),
              "Dth_new_deg": float(deg(th[q] - np.arctan(x[q] / R))), "Dth_ref_deg": float(deg(tr[q] - np.arctan(x[q] / R)))}
             for q in range(min(6, len(x)))]
    w = (xi >= FIT[0]) & (xi <= FIT[1])
    fit = None
    if w.sum() >= 3:
        b, a = np.polyfit(xi[w], dth[w], 1)
        fit = {"a_deg": float(a), "b_deg_per_rt": float(b), "n": int(w.sum())}
    win = {}
    for lab, lo, hi in (("near_0_0.05", 0.0, 0.05), ("mid_0.05_0.3", 0.05, 0.3)):
        s = (xi > lo) & (xi <= hi)
        win[lab] = None if not s.any() else {"max_abs_dth_deg": float(np.abs(dth[s]).max()),
                                             "max_abs_dr": float(np.abs(dr[s]).max()), "n": int(s.sum())}
    for lab, lo, hi in (("far_to_3", x_f + 0.3, 3.0), ("far_3_end", 3.0, np.inf)):
        s = (x > lo) & (x <= hi)
        win[lab] = None if not s.any() else {"max_abs_dth_deg": float(np.abs(dth[s]).max()),
                                             "max_abs_dr": float(np.abs(dr[s]).max()), "n": int(s.sum())}
    return {"n_compared": int(m.sum()), "first6": first, "fit_window_from_xf": list(FIT), "fit": fit, "windows": win,
            "end": {"x": float(x[-1]), "dth_deg": float(dth[-1]), "dr": float(dr[-1])} if len(x) else None}


def net_diff(lev_new, lev, i, m_max):
    """軸節点が同じ再出発で、列 < i の網の節点が元と一致するか (両方に有る節点の最大差)。θ・ν は度。
    節点 L_k[q] は起点 L_0[q+k] の C⁺ 上にあるので、q + k ≤ m_max (= 元の網と同じ起点を持つ C⁺) の節点だけ比べる。"""
    n2 = min(lev_new.shape[0], lev.shape[0])
    same_cplus = (np.arange(n2)[:, None] + np.arange(i)[None, :]) <= m_max
    out = {}
    for c, name in ((0, "x"), (1, "r"), (2, "th_deg"), (3, "nu_deg")):
        a, b = lev_new[:n2, :i, c], lev[:n2, :i, c]
        both = np.isfinite(a) & np.isfinite(b) & same_cplus
        dd = np.abs(a[both] - b[both])
        out[name] = float(deg(dd.max()) if c >= 2 else dd.max()) if dd.size else float("nan")
        if c == 0:
            out["n_both"] = int(both.sum())
            out["n_only_new"] = int((np.isfinite(a) & ~np.isfinite(b) & same_cplus).sum())
    return out


class _Line:
    """`inverse_design` に渡す throat の代役 (throat_characteristic だけ使われる)。"""

    def __init__(self, x, r, M, th):
        self.a = (np.asarray(x, float), np.asarray(r, float), np.asarray(M, float), np.asarray(th, float))

    def throat_characteristic(self, n=None, **_):
        return self.a


def resample(rows, n):
    """列の線 (軸 → 足) を r 等間隔 n 点に線形補間 (CFDPinnedThroat.throat_characteristic と同じく x・M・θ を r で補間)。"""
    r = rows[:, 1]
    if not np.all(np.diff(r) > 0):
        raise RuntimeError("列の r が単調増加でない")
    rq = np.linspace(0.0, r[-1], n)
    return np.interp(rq, r, rows[:, 0]), rq, np.interp(rq, r, rows[:, 4]), np.interp(rq, r, rows[:, 2])


def restart_one(cap, i, n_rs):
    inv, g, init, n_ax, lev, wall, R = (cap[k] for k in ("inv", "g", "init", "n_ax", "lev", "wall", "R"))
    col = column(lev, i)
    k, W = crossing(col, wall, g)
    if k is None:
        return None
    x_w = float(W[0])
    rec = {"column": int(i), "axis_x": float(col[0, 0]), "x_w": x_w, "r_w": float(W[1]), "k_first_above": int(k),
           "M_w": float(W[4]), "theta_w_deg": float(deg(W[2])),
           "wall_vs_column_dth_deg": float(deg(W[2] - np.interp(x_w, wall[:, 0], wall[:, 2]))), "variants": {}}
    axis_pts = init[:i]
    cum0_o = np.zeros(len(init))
    cum0_o[n_ax:] = MI._flux_along(init[n_ax:], g)
    # --- V1n: 足 = 壁の直下の網節点 (列 i の k*−1 番) -------------------------------------------
    line = col[:k]
    lev1, _, ms1, wall1 = run_net(inv, axis_pts + pts_exact(line, g), i, g)
    F, node = flux_to_node(lev, cum0_o, i + (k - 1), k - 1, g)
    if np.abs(node[:2] - line[-1, :2]).max() > 0:
        raise RuntimeError("C⁺ の添字と列の節点が食い違う")
    ref1 = MI.cplus_flux_wall(lev, cum0_o, F, g)
    rec["variants"]["V1n"] = {"x_f": float(line[-1, 0]), "r_f": float(line[-1, 1]), "mstar_new": ms1, "mstar_ref": F,
                              "mstar_rel": ms1 / F - 1.0, "cascade": cascade(lev1, i, len(line)),
                              "net_diff": net_diff(lev1, lev, i, i + k - 1), "cmp": compare(wall1, ref1, float(line[-1, 0]), R)}
    del lev1
    # --- V1w: 足 = 列と元の壁の交点 W -------------------------------------------------------------
    linew = np.vstack([col[:k], W])
    lev2, _, ms2, wall2 = run_net(inv, axis_pts + pts_exact(linew, g), i, g)
    rec["variants"]["V1w"] = {"x_f": x_w, "mstar_new": ms2, "mstar_ref": cap["mstar"], "mstar_rel": ms2 / cap["mstar"] - 1.0,
                              "cascade": cascade(lev2, i, len(linew)), "net_diff": net_diff(lev2, lev, i, i + k - 1),
                              "cmp": compare(wall2, wall, x_w, R)}
    del lev2
    # --- V2: W までの線を r 等間隔 n 点に補間し直す (軸節点は元のまま) ---------------------------------
    xq, rq, Mq, tq = resample(linew, n_rs)
    lev3, _, ms3, wall3 = run_net(inv, axis_pts + pts_like_production(xq, rq, Mq, tq, g), i, g)
    rec["variants"]["V2"] = {"x_f": x_w, "n_line": int(n_rs), "mstar_new": ms3, "mstar_rel": ms3 / cap["mstar"] - 1.0,
                             "cascade": cascade(lev3, i, n_rs), "cmp": compare(wall3, wall, x_w, R)}
    del lev3
    # --- V3: 同じ線 + 軸格子を新しい足から作り直す (inverse_design をそのまま) --------------------------
    kw = dict(cap["kw"])
    kw.update(n_axis=int(i) + 1, n_start=int(n_rs))      # 軸節点数は元の列より下流の数に合わせる
    try:
        r3 = MI.inverse_design(_Line(xq, rq, Mq, tq), cap["target"], **kw)
        ms4 = float(r3["mdot_start"])
        rec["variants"]["V3"] = {"x_f": x_w, "n_line": int(n_rs), "n_axis": int(i) + 1, "axis_dx0": kw.get("axis_dx0"),
                                 "mstar_new": ms4, "mstar_rel": ms4 / cap["mstar"] - 1.0,
                                 "cascade": cascade(r3["levels"], i, n_rs), "cmp": compare(r3["wall_full"], wall, x_w, R)}
        del r3
    except Exception as e:
        rec["variants"]["V3"] = {"error": f"{type(e).__name__}: {str(e)[:200]}"}
    return rec


def axis_anchor_check(cap):
    """3a: 初期線の軸端 (MOC が使う x・M) とアンカー (x_A・M_A)、law(x_A)、最初の軸節点。"""
    init, n_ax, d, ht = cap["init"], cap["n_ax"], cap["d"], cap["throat"]
    x_l, M_l = init[n_ax].x, init[n_ax].M
    out = {"x_line_axis": x_l, "M_line_axis": M_l, "nu_line_axis_deg": float(deg(init[n_ax].nu))}
    if d is not None:
        M_A, Mp_A, Mpp_A = d["anchor"]
        out.update(x_A=d["x_A"], M_A=M_A, Mp_A=Mp_A, Mpp_A=Mpp_A, dx=x_l - d["x_A"], dM=M_l - M_A,
                   law_xA_minus_M_A=float(d["law"](d["x_A"])) - M_A, anchor_source=d["anchor_source"])
    a = init[n_ax - 1]
    out["first_axis_node"] = {"x": a.x, "dx_from_line": a.x - x_l, "M": a.M,
                              "law_minus_M": (float(d["law"](a.x)) - a.M) if d is not None else None}
    if hasattr(ht, "_ax_x"):        # CFD ピン: M′_A を取る 4 次窓フィットの x0 での値 (M_A は線の軸端から)
        w = np.abs(ht._ax_x - x_l) <= ht._anchor_hw
        c = np.polyfit(ht._ax_x[w] - x_l, ht._ax_M[w], 4)
        out["cfd_axis_window_fit"] = {"M_at_x0": float(c[-1]), "Mp_at_x0": float(c[-2]), "M_at_x0_minus_M_A": float(c[-1]) - M_l,
                                      "evenfit_interp_at_x0": float(np.interp(x_l, ht._ax_x, ht._ax_M)),
                                      "field_M_axis_at_x0": float(ht._field(x_l, 0.0)[0])}
    else:                            # Hall: 場の軸上の M (同じ級数の u(x, 0))
        out["hall_mach_axis_at_x0"] = float(ht.mach(x_l, 0.0))
        out["hall_anchor_M_at_x0"] = float(ht.axis_anchor(x_l)[0])
    return out


def mstar_check(cap):
    """3b: m* (= 初期線を横切る流量。inverse_design はこれを定義に使い、外から与えない) と他の流量積分の相対差。"""
    ht, g, n_ax, init = cap["throat"], cap["g"], cap["n_ax"], cap["init"]
    ms = cap["mstar"]
    out = {"n_start": int(len(init) - n_ax), "mstar": ms, "mstar_over_pi": ms / np.pi,
           "note": "mass_flux_density は ρV/(ρ*a*) (semi-perfect 表)。1 次元の理想 ṁ* = π"}
    cum = MI._flux_along(init[n_ax:], g)
    out["top_segment_share"] = float((cum[-1] - cum[-2]) / ms)
    conv = {}
    for n in (41, 81, 161, 321, 641, 1281, 2561):
        xs, rr, MM, tt = ht.throat_characteristic(n=n)
        conv[str(n)] = float(MI._flux_along(pts_like_production(xs, rr, MM, tt, g), g)[-1]) / ms - 1.0
    out["line_flux_n_rel_to_mstar"] = conv
    if hasattr(ht, "_line_x"):     # CFD ピン: 追跡した細かい線そのもの (ds 2e-4) と、凍結源の場の縦断面 (全域超音速のもの)
        fine = pts_like_production(ht._line_x, ht._line_r, ht._line_M, ht._line_T, g)
        out["fine_traced_line_rel"] = float(MI._flux_along(fine, g)[-1]) / ms - 1.0
        out["fine_traced_line_npts"] = int(len(fine))
        x0 = init[n_ax].x
        sec = {}
        for xs_ in (x0, x0 + 0.25, x0 + 0.5):
            rw = float(ht._rw(xs_))
            rr = np.linspace(0.0, rw, 4001)
            Mv = np.ravel(ht._Ms.ev(np.full(rr.size, xs_), rr / rw))
            tv = np.ravel(ht._Ts.ev(np.full(rr.size, xs_), rr / rw))
            if Mv.min() < 1.0:
                sec[f"{xs_:.4f}"] = {"skip": f"亜音速を含む (M min {Mv.min():.4f})"}
                continue
            f = MI._mass_flux_density(Mv, g) * np.cos(tv) * 2.0 * np.pi * rr
            sec[f"{xs_:.4f}"] = {"rel": float(np.trapezoid(f, rr)) / ms - 1.0, "r_w": rw, "M_min": float(Mv.min())}
        out["cfd_vertical_sections_rel"] = sec
    else:
        out["hall_cd_series"] = float(ht.cd_series())
        out["mstar_over_pi_minus_cd_series"] = ms / np.pi - float(ht.cd_series())
    return out


def orig_summary(cap):
    wall, R = cap["wall"], cap["R"]
    x, th = wall[:, 0], wall[:, 2]
    D = deg(th - np.arctan(x / R))
    w = (x >= FIT[0]) & (x <= FIT[1])
    b, a = np.polyfit(x[w], D[w], 1)
    return {"first6": [{"x": float(x[q]), "Dth_deg": float(D[q])} for q in range(1, 7)],
            "fit": {"a_deg": float(a), "b_deg_per_rt": float(b), "n": int(w.sum())}}


def pick_columns(cap):
    """元の壁との交点 x が各目標に最も近い列 (重複は除く)。"""
    cand = []
    for i in range(cap["n_ax"] - 1, 0, -1):
        kk, W = crossing(column(cap["lev"], i), cap["wall"], cap["g"])
        if kk is None:
            continue
        cand.append((i, float(W[0])))
        if W[0] > max(TARGETS) * 1.5:
            break
    out = []
    for tx in TARGETS:
        i, _ = min(cand, key=lambda c: abs(c[1] - tx))
        if i not in [c[0] for c in out]:
            out.append((i, tx))
    return out


def main():
    t0 = time.time()
    commit = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "-C", str(ROOT), "status", "--porcelain", "design/forge_design"],
                           capture_output=True, text=True).stdout.strip()
    runs, mst = [], {}
    for lab, ker, L, do_rs in SPEC:
        lv = LV[L]
        t = time.time()
        with kernel(ker):
            cap = solve(PROBS[lab], **lv)
            n_ax = cap["n_ax"]
            # 自前の手順が生産 (同じ単位過程) の壁を再現するか (ビット一致)
            _, _, ms_chk, wall_chk = run_net(cap["inv"], cap["init"], n_ax, cap["g"])
            same = (wall_chk.shape == cap["wall"].shape and float(np.abs(wall_chk - cap["wall"]).max()) == 0.0
                    and ms_chk == cap["mstar"])
            if not same:
                raise RuntimeError(f"{lab} {ker} {L}: 自前の網が設計チェーンの壁を再現しない")
            del wall_chk
            cols = pick_columns(cap)
            s = {"label": lab, "problem": PROBS[lab], "kernel": ker, "level": L, **lv,
                 "design_chain_error": cap["err"], "selfcheck_bitwise": same,
                 "n_ax": n_ax, "x0": cap["init"][n_ax].x, "mstar": cap["mstar"],
                 "orig": orig_summary(cap), "cascade_orig": cascade(cap["lev"], n_ax, len(cap["init"]) - n_ax),
                 "first_row": first_row(cap, [n_ax - 1] + [c[0] for c in cols]),
                 "axis_anchor": axis_anchor_check(cap), "restarts": []}
            if ker == "prod" and L == "L1":
                mst[lab] = mstar_check(cap)
            if do_rs:
                for i, tx in cols:
                    rec = restart_one(cap, i, lv["n_start"])
                    rec["target_x"] = tx
                    s["restarts"].append(rec)
                    v = rec["variants"]
                    print(f"  {lab} {ker} {L}: 列 {i} x_w {rec['x_w']:.4f}  "
                          + "  ".join(f"{k}: a={v[k]['cmp']['fit']['a_deg']:+.5f}° d1={v[k]['cmp']['first6'][0]['dth_deg']:+.5f}°"
                                      if "cmp" in v[k] and v[k]["cmp"]["fit"] else
                                      (f"{k}: d1={v[k]['cmp']['first6'][0]['dth_deg']:+.5f}°" if "cmp" in v[k] else f"{k}: ERR")
                                      for k in v), flush=True)
        o = s["orig"]["fit"]
        print(f"{lab} {ker} {L} {lv}: 元の壁 Δθ [0.01,0.04] a={o['a_deg']:+.4f}° b={o['b_deg_per_rt']:+.3f}°/r_t  "
              f"第 1 点 {s['orig']['first6'][0]['Dth_deg']:+.4f}° (x {s['orig']['first6'][0]['x']:.4f})  "
              f"カスケード max|dθ| {s['cascade_orig']['max_abs_dth_deg']:.5f}°  err={cap['err']}  ({time.time() - t:.0f} s)",
              flush=True)
        runs.append(s)
        del cap
    out = {"plan": "plans/active/tooling-nozzle-throat-monotone-r2.md §5.1 #10 (b)", "commit": commit,
           "design_tree_dirty": dirty, "case_runs": str(RUNS), "R": 2.0, "targets_x": list(TARGETS), "fit_window": list(FIT),
           "levels": LV, "n_corr_conv": N_CORR_CONV, "runs": runs, "mstar": mst, "elapsed_s": time.time() - t0}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1, ensure_ascii=False, default=float))
    print(f"-> {OUT} ({time.time() - t0:.0f} s)")


if __name__ == "__main__":
    main()
