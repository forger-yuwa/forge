"""case/70: 設計壁の単調性の検査の A/B (plan tooling-nozzle-deltastar-band-root-zone §5.1 #4b、
諮問 notes/reviews/2026-10-13-design-wall-monotone-gate-diagnose.md、forge 0 step)。

壁は固定し、変えるのは判定方式だけ:
  A = 今の検査 (`AxisMachCFDWall.validate` の改修前の式: [x_in, x_e] の 4000 点の等間隔標本のうち x ≥ x0 の隣接差分 < −1e-9)
  B = 連続関数での複合ゲート `wall_axismach.monotone_gate` (設計壁 `_spl` の [x0, x_e] で 最小の r′ ≥ −1.74533e-6 = −tan(1e-4°)
      かつ 最大の累積減少 D ≤ 1e-7 r_t)
対象:
  - 走査した設計: c2b の要求カード (`card_m5_air_demo_c2b.yaml`) から r_t ごとに問題を作って design_chain
    (前の担当の走査と同じ r_t の列: c2b の r_t・0.02500〜0.02512 の 25 点・2 点)
  - 既存の通過設計: case/70 の problem_m5_air_c2b・final・c2・walltemp、case/45 の生産 YAML (problem_d155_ns_prod.yaml;
    凍結源 run_0062 はこの作業ツリーに無いので initial_line_run だけ /home/sano/work/forge/case/45 の絶対パスに置き換えて読む)
  - 事前固定の負例と対照 (c2b の設計壁のスプラインに局所の変形を足した区間多項式; 下の `deformed`)
壁は validate を呼ばずに組む (design_chain の間だけ `AxisMachCFDWall.validate` を空にする; 検査は副作用を持たないので壁は同じ)。
MOC の壁点列 (当てはめ前) の減少も別に記録する (当てはめ・補間が発生源かの切り分け)。

usage (case dir, design/.venv-opt):  python monotone_gate_ab.py
出力: _monotone_gate_ab/monotone_gate_ab.json と標準出力の表
"""
import copy
import inspect
import json
import sys
from pathlib import Path

import numpy as np
import yaml
from scipy.interpolate import PPoly

C = Path(__file__).resolve().parent
sys.path.insert(0, str(C.parents[1] / "design"))
from forge_design.evaluate.runner_axismach import design_chain  # noqa: E402
from forge_design.geometry import wall_axismach as WA  # noqa: E402
from forge_design.intake import wt_card as W  # noqa: E402
from forge_design.probdef import load_problem  # noqa: E402

OUT = C / "_monotone_gate_ab"
OUT.mkdir(exist_ok=True)
CARD = C / "card_m5_air_demo_c2b.yaml"
RT_SCAN = [0.024973561273227432] + [float(v) for v in np.linspace(0.02500, 0.02512, 25)] + [0.025077342164181705, 0.02507572544991335]
EXISTING = [("c70_c2b", C / "problem_m5_air_c2b.yaml"), ("c70_final", C / "problem_m5_air_final.yaml"),
            ("c70_c2", C / "problem_m5_air_c2.yaml"), ("c70_walltemp", C / "problem_m5_air_walltemp.yaml")]
C45_PROD = C.parents[1] / "case/45.isobutane_m6_d155/problem_d155_ns_prod.yaml"
C45_RUNS = Path("/home/sano/work/forge/case/45.isobutane_m6_d155")
A_N, A_TOL = 4000, 1e-9
DEG = np.pi / 180.0

_ORIG_VALIDATE = WA.AxisMachCFDWall.validate


def chain_without_validate(fn):
    """design_chain の間だけ設計壁の validate を空にして壁を取り出す (validate の例外の前の壁)。"""
    WA.AxisMachCFDWall.validate = lambda self, n=4000: []
    try:
        return fn()
    finally:
        WA.AxisMachCFDWall.validate = _ORIG_VALIDATE


def arm_a(wall) -> dict:
    """改修前の `AxisMachCFDWall.validate` の単調性の式そのまま (点の間隔に依存する)。"""
    xs = np.linspace(wall.x_in, wall.x_e, A_N)
    rv = wall.r(xs)
    m = xs >= wall.x0
    dr = np.diff(rv[m])
    xm = xs[m][1:]
    bad = dr < -A_TOL
    return {"pass": not bool(np.any(bad)), "n_bad": int(bad.sum()), "min_dr": float(dr.min()),
            "x_min_dr": float(xm[int(np.argmin(dr))]), "dx": float(xs[1] - xs[0])}


def arm_b(wall) -> dict:
    g = WA.monotone_gate(wall._spl, wall.x0, wall.x_e)
    if g["applicable"]:
        g["fail_slope"] = not g["min_r1"] >= g["r1_min"]
        g["fail_D"] = not g["D"] <= g["d_max"]
    return g


def moc_points(tbl) -> dict:
    """MOC の壁点列 (当てはめ前) の減少: 隣接差分の最小・負の数・点列の累積減少・点列の θ の最小。"""
    tb = np.asarray(tbl, dtype=float)
    x, r, th = tb[:, 0], tb[:, 1], tb[:, 2]
    dr = np.diff(r)
    gap = np.maximum.accumulate(r) - r
    i = int(np.argmin(dr))
    return {"n": int(len(r)), "min_dr": float(dr[i]), "x_min_dr": float(x[i + 1]), "n_neg_dr": int(np.sum(dr < 0.0)),
            "D_pts": float(gap.max()), "x_D_pts": float(x[int(np.argmax(gap))]),
            "min_theta_deg": float(np.degrees(th.min())), "x_min_theta": float(x[int(np.argmin(th))]),
            "n_neg_theta": int(np.sum(th < 0.0))}


def record(name, kind, d, src, rt=None) -> dict:
    w = d["wall"]
    rec = {"name": name, "kind": kind, "source": src, "r_t_m": rt, "wall_class": type(w).__name__,
           "spl_class": type(w._spl).__name__, "x0": float(w.x0), "x_e": float(w.x_e), "x_in": float(w.x_in),
           "L_c": float(d["L_c"]), "Lc_solve": d.get("Lc_solve"),
           "fit_max_dr_pts": (w.fit_diag or {}).get("max_dr_pts") if getattr(w, "fit_diag", None) else None,
           "A": arm_a(w), "B": arm_b(w), "moc": moc_points(d["wall_inv"]),
           "validate_now": _ORIG_VALIDATE(w)}
    return rec


# --- 負例と対照 (事前固定; plan §5.1 #4b) -------------------------------------------------------------------------------
def _shift(c, d):
    """係数 c (高次から, 変数 s) の多項式 f(s) → g(s') = f(s' + d) の係数。"""
    return np.poly1d(np.asarray(c, dtype=float))(np.poly1d([1.0, float(d)])).coeffs


def _smooth_up(x0_, T):
    """[x0_, x0_+T] で 0 → 1 の 5 次 smoothstep (x の多項式の係数を返す関数: 区間の始点 a で局所座標に写す)。"""
    h = np.poly1d([6.0, -15.0, 10.0, 0.0, 0.0, 0.0])
    return lambda a: h(np.poly1d([1.0 / T, (a - x0_) / T])).coeffs


def _smooth_down(x1_, T):
    h = np.poly1d([6.0, -15.0, 10.0, 0.0, 0.0, 0.0])
    return lambda a: h(np.poly1d([-1.0 / T, (x1_ - a) / T])).coeffs


def _bump(xc, half, amp):
    """amp·(1 − u²)², u = (x − xc)/half の局所座標の係数。"""
    Q = np.poly1d([1.0, 0.0, -2.0, 0.0, 1.0])
    return lambda a: (amp * Q(np.poly1d([1.0 / half, (a - xc) / half]))).coeffs


def deformed(S, x0, x_e, plateau, T, q_amp, comp_side, comp_len, comp_gap):
    r"""設計壁 S (BSpline) に局所の変形 φ を足した区間多項式 r = S + φ (PPoly) を返す。

    窓 [p0 − T, p1 + T] で r′ を S′ から q に置き換える: φ′ = χ (q − S′) + comp。χ は両端 5 次 smoothstep・平坦部 [p0, p1] で 1、
    q = q_amp (1 − u²)² (平坦部だけ; q_amp < 0 で収縮)。平坦部の r′ は q そのもの (最小 q_amp)、累積の変化は q_amp·(8/15)·L。
    遷移部の r′ = S′ (1 − χ) ≥ 0 なので減少は平坦部だけ。comp ≥ 0 は窓で失った高さを窓の前 (pre) か後 (post) で戻す
    5 次の山 (∫comp = ∫χ (S′ − q)) で、φ は窓と comp の外で 0 (壁の他の部分は元のまま)。"""
    p0, p1 = plateau
    L = p1 - p0
    win = (p0 - T, p1 + T)
    if comp_side == "post":
        cr = (win[1] + comp_gap, win[1] + comp_gap + comp_len)
    else:
        cr = (win[0] - comp_gap - comp_len, win[0] - comp_gap)
    assert x0 < min(win[0], cr[0]) and max(win[1], cr[1]) < x_e, (win, cr, x0, x_e)
    Spp = PPoly.from_spline(S)
    bks = np.unique(np.concatenate([Spp.x[(Spp.x >= x0) & (Spp.x <= x_e)], [x0, x_e, win[0], p0, p1, win[1], cr[0], cr[1]]]))
    bks = bks[(bks >= x0) & (bks <= x_e)]
    up, dn = _smooth_up(win[0], T), _smooth_down(win[1], T)
    q = _bump(0.5 * (p0 + p1), 0.5 * L, q_amp)

    def chi(a, mid):
        if mid < win[0] or mid > win[1]:
            return np.array([0.0])
        if mid < p0:
            return up(a)
        if mid > p1:
            return dn(a)
        return np.array([1.0])

    def qq(a, mid):
        return q(a) if p0 < mid < p1 else np.array([0.0])

    def S_loc(a, mid):
        j = int(np.searchsorted(Spp.x, mid, side="right") - 1)
        return _shift(Spp.c[:, j], a - Spp.x[j])

    # 窓で失う高さ ∫χ (S′ − q) → comp の振幅
    loss = 0.0
    for a, b in zip(bks[:-1], bks[1:]):
        mid = 0.5 * (a + b)
        f = np.polymul(chi(a, mid), np.polysub(np.polyder(S_loc(a, mid)), qq(a, mid)))
        loss += float(np.polyval(np.polyint(f), b - a))
    K = loss / (comp_len * 8.0 / 15.0)
    comp = _bump(0.5 * (cr[0] + cr[1]), 0.5 * comp_len, K)
    pieces, phi = [], 0.0
    for a, b in zip(bks[:-1], bks[1:]):
        mid = 0.5 * (a + b)
        Sl = S_loc(a, mid)
        dphi = np.polymul(chi(a, mid), np.polysub(qq(a, mid), np.polyder(Sl)))
        if cr[0] < mid < cr[1]:
            dphi = np.polyadd(dphi, comp(a))
        ph = np.polyint(dphi, k=phi)
        phi = float(np.polyval(ph, b - a))
        pieces.append(np.polyadd(Sl, ph))
    deg = max(len(p) for p in pieces) - 1
    c = np.zeros((deg + 1, len(pieces)))
    for j, p in enumerate(pieces):
        c[deg + 1 - len(p):, j] = p
    return PPoly(c, bks), {"plateau": [p0, p1], "T": T, "window": list(win), "comp": list(cr), "comp_amp": K,
                           "height_loss": loss, "phi_end": phi, "q_amp": q_amp,
                           "design_min_theta_deg": float(np.degrees(np.arctan(min(q_amp, 0.0)))),
                           "design_D": float(-q_amp * 8.0 / 15.0 * L) if q_amp < 0 else 0.0}


def dense_check(pp, x0, x_e, window, n=4_000_001, nw=2_000_001) -> dict:
    """負例の独立な確かめ (gate を使わない): 全域と窓の密な標本の r′ の最小と累積減少 (標本は真値の min r′ を上から、D を下から近似)。"""
    xs = np.unique(np.concatenate([np.linspace(x0, x_e, n), np.linspace(window[0], window[1], nw)]))
    r, r1 = pp(xs), pp(xs, 1)
    gap = np.maximum.accumulate(r) - r
    i = int(np.argmin(r1))
    return {"min_r1": float(r1[i]), "min_theta_deg": float(np.degrees(np.arctan(r1[i]))), "x_min_theta": float(xs[i]),
            "D": float(gap.max()), "x_D": float(xs[int(np.argmax(gap))]), "n": int(len(xs))}


def main():
    out = {"plan": "plans/active/tooling-nozzle-deltastar-band-root-zone.md §5.1 #4b",
           "consult": "notes/reviews/2026-10-13-design-wall-monotone-gate-diagnose.md",
           "A": f"[x_in, x_e] の {A_N} 点の等間隔標本、x ≥ x0 の隣接差分 < −{A_TOL:g} で不合格",
           "B": f"monotone_gate: [x0, x_e] で min r′ ≥ {WA.MONO_R1_MIN:g} (= −tan(1e-4°)) かつ D ≤ {WA.MONO_D_MAX:g} r_t",
           "validate_uses_gate": "monotone_gate(" in inspect.getsource(_ORIG_VALIDATE),
           "designs": [], "negatives": []}
    card = W.load_card(CARD)
    for rt in RT_SCAN:
        prob = W.card_to_problem(card, float(rt))
        try:
            d = chain_without_validate(lambda: W._design_chain_of(prob))
        except Exception as ex:      # noqa: BLE001 — 設計が組めないものは記録して続ける
            out["designs"].append({"name": f"scan_rt{rt:.9f}", "kind": "scan", "r_t_m": rt, "error": str(ex)[:300]})
            print(f"scan rt {rt:.9f}: design_chain error {str(ex)[:120]}", flush=True)
            continue
        out["designs"].append(record(f"scan_rt{rt:.9f}", "scan", d, str(CARD.name), rt))
        print(_row(out["designs"][-1]), flush=True)
    for name, path in EXISTING:
        p = load_problem(path)
        d = chain_without_validate(lambda: design_chain(p))
        out["designs"].append(record(name, "existing", d, str(path.relative_to(C.parents[1])), float(p.spec["r_throat"])))
        print(_row(out["designs"][-1]), flush=True)
    # case/45 の生産 YAML (凍結源の run だけ別の作業ツリー)
    raw = yaml.safe_load(C45_PROD.read_text())
    raw["geometry"]["initial_line_run"] = str(C45_RUNS / raw["geometry"]["initial_line_run"])
    tmp = OUT / "_c45_prod_abs_initial_line.yaml"
    tmp.write_text(yaml.safe_dump(raw, sort_keys=False, allow_unicode=True))
    p = load_problem(tmp)
    d = chain_without_validate(lambda: design_chain(p))
    out["designs"].append(record("c45_prod", "existing", d, str(C45_PROD.relative_to(C.parents[1])) +
                                 f" (initial_line_run → {raw['geometry']['initial_line_run']})", float(p.spec["r_throat"])))
    print(_row(out["designs"][-1]), flush=True)

    # 負例と対照: c2b の設計壁に局所の変形
    base_d = chain_without_validate(lambda: design_chain(load_problem(C / "problem_m5_air_c2b.yaml")))
    bw = base_d["wall"]
    S, x0, x_e = bw._spl, float(bw.x0), float(bw.x_e)
    a_slope = float(np.tan(2e-4 * DEG))          # min θ = −2e-4°
    a_cum = 1.5e-6                               # min θ = −8.594e-5° ≥ −1e-4°
    L_cum = 2e-7 / (a_cum * 8.0 / 15.0)          # D = 2e-7 r_t
    L_slope = 0.01                               # 格子点の間 (A の標本間隔 ≈ 0.017 r_t) の収縮、D = 1.86e-8 < 1e-7
    xm = 24.0                                    # 拡大部の途中 (θ ≈ 3.4°)
    pe = x_e - 0.06                              # 出口直前: 平坦部の終わり
    spec = [
        ("neg_slope_mid", "負例: 格子点の間の収縮 (傾き)", (xm - L_slope / 2, xm + L_slope / 2), 0.05, -a_slope, "post", "slope"),
        ("neg_slope_exit", "負例: 格子点の間の収縮 (傾き)", (pe - L_slope, pe), 0.03, -a_slope, "pre", "slope"),
        ("neg_cum_mid", "負例: 緩い収縮 (累積)", (xm - L_cum / 2, xm + L_cum / 2), 0.05, -a_cum, "post", "D"),
        ("neg_cum_exit", "負例: 緩い収縮 (累積)", (pe - L_cum, pe), 0.03, -a_cum, "pre", "D"),
        ("ctl_mono_mid", "対照: 厳密単調 (平坦部 q > 0)", (xm - L_cum / 2, xm + L_cum / 2), 0.05, +a_cum, "post", None),
        ("ctl_mono_exit", "対照: 厳密単調 (平坦部 q > 0)", (pe - L_cum, pe), 0.03, +a_cum, "pre", None),
        ("ctl_par_mid", "対照: 平行 (平坦部 r′ = 0)", (xm - L_cum / 2, xm + L_cum / 2), 0.05, 0.0, "post", None),
        ("ctl_par_exit", "対照: 平行 (平坦部 r′ = 0)", (pe - L_cum, pe), 0.03, 0.0, "pre", None),
    ]
    out["negatives"].append({"name": "ctl_base_c2b", "desc": "対照: 変形なしの c2b の設計壁", "expect": None,
                             "A": arm_a(bw), "B": arm_b(bw)})
    print(_nrow(out["negatives"][-1]), flush=True)
    for name, desc, plateau, T, qa, side, expect in spec:
        pp, info = deformed(S, x0, x_e, plateau, T, qa, side, comp_len=0.4, comp_gap=0.1)
        wn = copy.copy(bw)
        wn._spl = pp
        rec = {"name": name, "desc": desc, "expect": expect, "construction": info,
               "max_abs_phi_outside": float(np.max(np.abs(pp(np.r_[np.linspace(x0, min(info["window"][0], info["comp"][0]), 2001),
                                                                    np.linspace(max(info["window"][1], info["comp"][1]), x_e, 2001)])
                                                         - S(np.r_[np.linspace(x0, min(info["window"][0], info["comp"][0]), 2001),
                                                                   np.linspace(max(info["window"][1], info["comp"][1]), x_e, 2001)])))),
               "dense": dense_check(pp, x0, x_e, info["window"]), "A": arm_a(wn), "B": arm_b(wn)}
        b = rec["B"]
        if expect is None:
            rec["ok"] = bool(b["pass"])
        elif expect == "slope":
            rec["ok"] = (not b["pass"]) and b["fail_slope"] and not b["fail_D"]
        else:
            rec["ok"] = (not b["pass"]) and b["fail_D"] and not b["fail_slope"]
        out["negatives"].append(rec)
        print(_nrow(rec), flush=True)
    out["negatives"][0]["ok"] = bool(out["negatives"][0]["B"]["pass"])

    # 判定 (plan §5.1 #4b; 結果を見て閾値を広げない)
    real = [r for r in out["designs"] if "error" not in r]
    a_fail = [r["name"] for r in real if not r["A"]["pass"]]
    a_fail_b_fail = [r["name"] for r in real if not r["A"]["pass"] and not r["B"]["pass"]]
    a_pass_b_fail = [r["name"] for r in real if r["A"]["pass"] and not r["B"]["pass"]]
    neg_ok = all(r["ok"] for r in out["negatives"] if r["expect"] is not None)
    ctl_ok = all(r["ok"] for r in out["negatives"] if r["expect"] is None)
    adopt = bool(a_fail) and not a_fail_b_fail and not a_pass_b_fail and neg_ok and ctl_ok
    out["judgement"] = {"n_real": len(real), "n_error": len(out["designs"]) - len(real), "A_fail": a_fail,
                        "A_fail_B_fail": a_fail_b_fail, "A_pass_B_fail": a_pass_b_fail,
                        "negatives_all_rejected_for_expected_reason": neg_ok, "controls_all_pass": ctl_ok,
                        "adopt_B": adopt,
                        "rule": "A で落ちる実壁が B で通り、負例を全部 (予定の理由で) 拒否し、対照が通る → B を採用。"
                                "実壁が B でも落ちる → 「微小だから許容」を棄却。負例が通る → B の実装を不採用"}
    print("\n判定:", json.dumps(out["judgement"], ensure_ascii=False, indent=1))
    (OUT / "monotone_gate_ab.json").write_text(json.dumps(out, indent=1, ensure_ascii=False, default=float))


def _row(r):
    a, b, m = r["A"], r["B"], r["moc"]
    rt = f"{r['r_t_m']:.9f}" if r["r_t_m"] is not None else "-"
    return (f"{r['name']:<22} rt {rt} x_e {r['x_e']:.6f} | A {'pass' if a['pass'] else 'FAIL'} n_bad {a['n_bad']:2d} "
            f"min Δr {a['min_dr']:+.3e} @ {a['x_min_dr']:.4f} | B {'pass' if b['pass'] else 'FAIL'} "
            f"min θ {b['min_theta_deg']:+.3e}° @ {b['x_min_theta']:.4f} D {b['D']:.3e} @ {b['x_D']:.4f} | "
            f"MOC min Δr {m['min_dr']:+.3e} n_neg {m['n_neg_dr']} D_pts {m['D_pts']:.2e} min θ {m['min_theta_deg']:+.2e}°")


def _nrow(r):
    a, b = r["A"], r["B"]
    dn = r.get("dense")
    s = (f"{r['name']:<16} expect {str(r['expect']):<5} | A {'pass' if a['pass'] else 'FAIL'} n_bad {a['n_bad']:2d} "
         f"min Δr {a['min_dr']:+.3e} | B {'pass' if b['pass'] else 'FAIL'} min θ {b['min_theta_deg']:+.4e}° "
         f"@ {b['x_min_theta']:.5f} D {b['D']:.4e} @ {b['x_D']:.5f}")
    if dn:
        s += f" | dense min θ {dn['min_theta_deg']:+.4e}° D {dn['D']:.4e}"
    if "ok" in r:
        s += f" | {'ok' if r['ok'] else 'NG'}"
    return s


if __name__ == "__main__":
    main()
