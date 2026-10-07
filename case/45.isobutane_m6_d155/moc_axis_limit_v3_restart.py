"""plan discretization-moc-axis-limit-and-corrector §6 V3: 再出発の自己整合 (記録)。CFD 0 step。

`throat_moc_restart_test.py` (monotone plan §5.1 #10 (b)) と同じ再出発試験を、新しい手順 (geometry.moc_axis_limit: analytic +
moc_corrector: converge) で回す。単位過程の差し替え (monkeypatch) はせず、問題 YAML のキーで選ぶ (生産経路そのもの)。
**記録であって物理的な正しさの証明ではない** (K2c の誤った源項でも通るため。plan §6 V3)。

試験 (列 i = 軸節点 L_0[i] から上る C⁻、W = 列と元の壁の交点):
  V1n  列 i の網の節点 (軸 → 壁の直下) を**補間し直さず**初期線にして充填し直す。軸節点 init[:i] は元のまま。
       軸端点の θ_r は元の値 (analytic のとき)。列 < i の網の節点 (同じ C⁺ を持つもの) の x・r・θ・ν の再現を見る。
  V1w  V1n の線に W を足したもの (比較相手 = 元の壁)。
  軸端のメタデータ: 列の軸端 (= L_0[i]) の x・r・θ・ν・M と θ_r が元の軸節点と一致し、θ_r が軸則から作り直した値と一致するか。
  列の 1 段目の作り直し: 対 (L_1[i], L_0[i]) を同じ単位過程に通し直したときの L_1[i] との θ の差 (全列の最大と、選んだ列)。
条件: 単位過程 {current = legacy + fixed2 (生産)、new = analytic + converge、ref = legacy + converge (修正子だけ)} ×
      分解能 {P = 生産 (n_start 41・n_axis 2400・dx0 0.03)、L1 = n_start 161 (restart test の L1)}。
列の選び方・比較の式は throat_moc_restart_test.py の関数をそのまま使う (pick_columns・crossing・cascade・compare・net_diff)。

usage: [CASE_RUNS=<run_0062 のある case dir>] python3 moc_axis_limit_v3_restart.py [OUT_JSON]
       出力 _band_ab/moc_axis_limit_v3_restart.json
"""
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

C = Path(__file__).resolve().parent
ROOT = C.parents[1]
sys.path.insert(0, str(ROOT / "design"))
sys.path.insert(0, str(C))
from forge_design.evaluate import runner_axismach as RA  # noqa: E402
from forge_design.geometry import moc_inverse as MI  # noqa: E402
from forge_design.geometry.moc_kernel import axis_theta_r, interior_vec  # noqa: E402
import throat_moc_restart_test as T  # noqa: E402  (関数だけ使う。main は走らない)

RUNS = Path(os.environ.get("CASE_RUNS", "/home/sano/work/forge/case/45.isobutane_m6_d155"))
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else C / "_band_ab/moc_axis_limit_v3_restart.json"
PROB = "problem_d155_ns_finemesh_recal_final.yaml"     # restart test の "pin" と同じ (単調壁の問題と MOC 網は同一)
KERNELS = {"current": None, "new": ("analytic", "converge"), "ref": ("legacy", "converge")}
LV = {"P": dict(n_start=41, n_axis=2400, dx0=0.03), "L1": dict(n_start=161, n_axis=2400, dx0=0.03)}
deg = np.degrees


def solve(keys, n_start, n_axis, dx0):
    """生産の design_chain (キーで単位過程を選ぶ) を回し、MOC 網・初期前線・θ_r・target を捕まえる。"""
    p = RA.load_problem(C / PROB)
    p.geometry.update(n_start=int(n_start), n_axis_inv=int(n_axis), axis_dx0=float(dx0))
    p.geometry["initial_line_run"] = str((RUNS / p.geometry["initial_line_run"]).resolve())
    if keys is not None:
        p.geometry["moc_axis_limit"], p.geometry["moc_corrector"] = keys
    cap = {"n_calls": 0}
    o_dc, o_inv = MI._design_cplus, RA.inverse_design

    def dc(inv, init, n_ax, ax, mdot_star, g, *a, **k):
        out = o_dc(inv, init, n_ax, ax, mdot_star, g, *a, **k)
        cap.update(inv=inv, init=init, n_ax=int(n_ax), mstar=float(mdot_star), g=g, thr=k.get("axis_thr"),
                   lev=out["levels"], wall=out["wall_full"])
        return out

    def invd(throat, target, **kw):
        cap["n_calls"] += 1
        cap.update(throat=throat, target=target, kw=dict(kw))
        return o_inv(throat, target, **kw)
    MI._design_cplus, RA.inverse_design = dc, invd
    try:
        d, err = RA.design_chain(p), None
    except Exception as e:  # noqa: BLE001  下流のゲートで落ちても網は捕まえてある (記録して続ける)
        d, err = None, f"{type(e).__name__}: {str(e)[:300]}"
    finally:
        MI._design_cplus, RA.inverse_design = o_dc, o_inv
    if "lev" not in cap or cap["n_calls"] != 1:
        raise RuntimeError(f"MOC 網を 1 回だけ捕まえられなかった (呼び出し {cap['n_calls']} 回, err={err})")
    cap.update(d=d, err=err, R=float(p.geometry.get("R", 2.0)))
    return cap


def new_inv(cap):
    inv = cap["inv"]
    return MI.InverseMOC(gamma=inv.g, delta=inv.delta, n_corr=inv.n_corr, axis_limit=inv.axis_limit,
                         corrector=inv.corrector, tol=inv.tol, max_corr=inv.max_corr)


def run_net(cap, init, thr):
    """`_design_cplus` と同じ手順で網と壁を作る (exit 処理は除く)。"""
    inv = new_inv(cap)
    g = cap["g"]
    n_ax = cap["_n_ax_run"]
    lev = inv.fill_levels(init, axis_thr=thr)
    cum0 = np.zeros(len(init))
    cum0[n_ax:] = MI._flux_along(init[n_ax:], g)
    ms = float(MI._flux_along(init[n_ax:], g)[-1])
    return lev, ms, MI.cplus_flux_wall(lev, cum0, ms, g), inv.last_diag


def reproc_rows(cap, cols):
    """列の 1 段目の作り直し: 対 (A = L_1[i]、B = L_0[i]) を同じ単位過程に通し直した L_1[i] との差 (θ・ν は度)。"""
    lev, g, inv, thr = cap["lev"], cap["g"], cap["inv"], cap["thr"]
    ii = np.array([i for i in cols if np.isfinite(lev[1, i, 0])])
    A, B = lev[1, ii], lev[0, ii]
    kw = {}
    if inv.axis_limit == "analytic":
        kw.update(axis_limit="analytic", thrA=np.full(len(ii), np.nan), thrB=thr[ii])
    if inv.corrector == "converge":
        kw.update(corrector="converge", tol=inv.tol, max_corr=inv.max_corr)
    q = interior_vec(*(A[:, c] for c in range(5)), *(B[:, c] for c in range(5)), g, inv.delta, inv.n_corr, **kw)
    return ii, deg(q[2] - A[:, 2]), deg(q[3] - A[:, 3])


def axis_meta(cap, i, line):
    """列の軸端のメタデータ: 再出発の初期線の軸端 (= L_0[i]) が元の軸節点と一致し、θ_r が軸則から作り直した値と一致するか。"""
    lev, thr, kw = cap["lev"], cap["thr"], cap["kw"]
    a = lev[0, i]
    out = {"x": float(a[0]), "same_xrthnuM": bool(np.array_equal(line[0], a))}
    if thr is not None:
        M_t = float(cap["target"](float(a[0])))
        Mp_t = float(kw["target_dM"](float(a[0])))
        th_re = float(axis_theta_r(M_t, Mp_t, cap["g"]))
        out.update(theta_r=float(thr[i]), theta_r_rebuilt=th_re, dtheta_r=float(thr[i] - th_re),
                   M_target_minus_M_node=M_t - float(a[4]))
    return out


def restart_one(cap, i):
    lev, init, wall, g, R, thr = cap["lev"], cap["init"], cap["wall"], cap["g"], cap["R"], cap["thr"]
    n_ax = cap["n_ax"]
    col = T.column(lev, i)
    k, W = T.crossing(col, wall, g)
    if k is None:
        return None
    rec = {"column": int(i), "axis_x": float(col[0, 0]), "x_w": float(W[0]), "k_first_above": int(k), "variants": {}}
    cum0_o = np.zeros(len(init))
    cum0_o[n_ax:] = MI._flux_along(init[n_ax:], g)
    axis_pts = init[:i]
    for name, line in (("V1n", col[:k]), ("V1w", np.vstack([col[:k], W]))):
        new = axis_pts + T.pts_exact(line, g)
        thr2 = None
        if thr is not None:
            thr2 = np.full(len(new), np.nan)
            thr2[:i + 1] = thr[:i + 1]                    # 列の軸端 = 元の軸節点 L_0[i] (同じ θ_r)
        cap["_n_ax_run"] = i
        lev2, ms2, wall2, dg2 = run_net(cap, new, thr2)
        v = {"n_line": int(len(line)), "mstar_new": ms2, "axis_meta": axis_meta(cap, i, line),
             "cascade": T.cascade(lev2, i, len(line)), "net_diff": T.net_diff(lev2, lev, i, i + k - 1),
             "unconverged_pairs": int(dg2["pairs"]["iter_nonfinite"] + dg2["pairs"]["iter_maxiter"]),
             "iters_max": dg2["iters"]["max"]}
        if name == "V1n":
            F, node = T.flux_to_node(lev, cum0_o, i + (k - 1), k - 1, g)
            ref = MI.cplus_flux_wall(lev, cum0_o, F, g)
            v.update(mstar_ref=F, mstar_rel=ms2 / F - 1.0, cmp=T.compare(wall2, ref, float(line[-1, 0]), R))
        else:
            v.update(mstar_ref=cap["mstar"], mstar_rel=ms2 / cap["mstar"] - 1.0, cmp=T.compare(wall2, wall, float(W[0]), R))
        rec["variants"][name] = v
        del lev2
    return rec


def main():
    t0 = time.time()
    commit = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "-C", str(ROOT), "status", "--porcelain", "design/forge_design"],
                           capture_output=True, text=True).stdout.strip()
    runs = []
    for L, lv in LV.items():
        for kname, keys in KERNELS.items():
            t = time.time()
            cap = solve(keys, **lv)
            n_ax = cap["n_ax"]
            cap["_n_ax_run"] = n_ax
            _, ms_chk, wall_chk, _ = run_net(cap, cap["init"], cap["thr"])
            same = (wall_chk.shape == cap["wall"].shape and float(np.abs(wall_chk - cap["wall"]).max()) == 0.0
                    and ms_chk == cap["mstar"])
            if not same:
                raise RuntimeError(f"{kname} {L}: 自前の網が設計チェーンの壁を再現しない")
            cols = T.pick_columns(cap)
            ii, dth, dnu = reproc_rows(cap, range(1, n_ax))
            s = {"kernel": kname, "keys": keys, "level": L, **lv, "design_chain_error": cap["err"],
                 "selfcheck_bitwise": same, "n_ax": n_ax, "x0": cap["init"][n_ax].x,
                 "moc_gate": (cap["d"]["moc"]["gate"] if cap["d"] is not None else None),
                 "reproc_first_row_all": {"n": int(len(ii)), "max_abs_dth_deg": float(np.abs(dth).max()),
                                          "x_at_max": float(cap["lev"][0, ii[int(np.argmax(np.abs(dth)))], 0]),
                                          "max_abs_dnu_deg": float(np.abs(dnu).max())},
                 "orig": T.orig_summary(cap), "restarts": []}
            for i, tx in cols:
                rec = restart_one(cap, i)
                if rec is None:
                    continue
                rec["target_x"] = tx
                j = int(np.flatnonzero(ii == i)[0]) if np.any(ii == i) else None
                rec["reproc_first_row_dth_deg"] = None if j is None else float(dth[j])
                s["restarts"].append(rec)
                v = rec["variants"]
                print(f"  {kname} {L}: 列 {i} x_w {rec['x_w']:.4f}  1段目作り直し dθ {rec['reproc_first_row_dth_deg']:+.2e}°  "
                      + "  ".join(f"{nm}: 網 dθ {v[nm]['net_diff']['th_deg']:.2e}° dν {v[nm]['net_diff']['nu_deg']:.2e}° "
                                  f"dx {v[nm]['net_diff']['x']:.1e} 軸端一致 {v[nm]['axis_meta']['same_xrthnuM']} "
                                  f"d1 {v[nm]['cmp']['first6'][0]['dth_deg']:+.2e}°" for nm in v), flush=True)
            print(f"{kname} {L}: 1 段目の作り直し max|dθ| {s['reproc_first_row_all']['max_abs_dth_deg']:.3e}° "
                  f"(x {s['reproc_first_row_all']['x_at_max']:.3f})  元の壁 第 1 点 Δθ {s['orig']['first6'][0]['Dth_deg']:+.4f}°  "
                  f"({time.time() - t:.0f} s)", flush=True)
            runs.append(s)
            del cap
    out = {"plan": "plans/active/discretization-moc-axis-limit-and-corrector.md §6 V3 (記録)", "commit": commit,
           "design_tree_dirty": dirty, "case_runs": str(RUNS), "problem": PROB, "levels": LV,
           "kernels": {k: v for k, v in KERNELS.items()}, "targets_x": list(T.TARGETS), "fit_window": list(T.FIT),
           "runs": runs, "elapsed_s": time.time() - t0}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1, ensure_ascii=False, default=float))
    print(f"-> {OUT} ({time.time() - t0:.0f} s)")


if __name__ == "__main__":
    main()
