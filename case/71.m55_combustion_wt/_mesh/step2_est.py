"""M5.5 (case/71) 見積もり step 2: 3 つの作り方 (A) 2 値 + ブレンド・(B) 第 1 層の表 + 今の流れ方向の分布・(C) 第 1 層の表 + 密度の表
の最小の ni。使い捨て (コードは変えない)。step1_chain.py の出力 (wall.pkl・sections.npz・est.json・prob_parts.json) を読む。
引数: --tag <名前> --wall wall.pkl --sections sections.npz --est est.json --parts prob_parts.json --nj 97 --sc 0.01 [--build]"""
import sys, json, math, time, pickle, argparse
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path.home() / "forge-mknot/design"))
from forge_design.intake import wt_card as W
from forge_design.meshing import mesh2d as M2

ap = argparse.ArgumentParser()
ap.add_argument("--dir", default=str(Path.home() / "m55_mesh_est"))
ap.add_argument("--tag", default="m55")
ap.add_argument("--nj", type=int, default=97)
ap.add_argument("--sc", type=float, default=0.01, help="入口の前縁の切り落とし: 入口端からの弧長 s < sc [r_t] は sc の値")
ap.add_argument("--target_scale", type=float, default=1.0, help="第 1 層の表の倍率 (1.25 = y1+ 目標 1.0)")
ap.add_argument("--ni_ref", type=int, default=2000, help="dx の上限の基準の ni (雛形の throat_refine・throat_width で)")
ap.add_argument("--refine", type=float, default=4.0)
ap.add_argument("--width", type=float, default=3.0)
ap.add_argument("--skipA", action="store_true")
ap.add_argument("--frac_table_json", default=None, help="第 1 層の表を外から与える ([[x, frac], ...] の JSON)")
ap.add_argument("--LU", type=float, default=None)
ap.add_argument("--build", action="store_true", help="(C) の格子を mesh2d で実際に作って測る")
ap.add_argument("--buildA", action="store_true", help="(A) の格子も実際に作って測る")
a = ap.parse_args()
D = Path(a.dir)
ARM = W.AR_MARGIN * 5000.0
log = lambda *x: print(*x, flush=True)  # noqa: E731

wall = pickle.load(open(D / "wall.pkl", "rb"))
sec = np.load(D / "sections.npz")
est = json.load(open(D / "est.json"))
parts = json.load(open(D / "prob_parts.json"))
mesh = dict(parts["mesh"])
x = sec["x_rt"]; f = sec["y1_over_rw"].copy(); y1m = sec["y1_m"]; Tw = sec["Tw_K"]
rt = float(parts["spec"]["r_throat"])
x_in, x_e = float(wall.x_in), float(wall.x_e)
rw = np.asarray(wall.r(x), dtype=float)
drw = np.gradient(rw, x)
s_rt = np.concatenate([[0.0], np.cumsum(np.hypot(np.diff(x), np.diff(rw)))])
LU = a.LU if a.LU is not None else (float(wall.up.L_U) if hasattr(wall, "up") else None)
log(f"[{a.tag}] r_t {rt} m  x_in {x_in:.4f}  x_e {x_e:.4f}  L_U {LU}  nj {a.nj}  AR 上限 {ARM:g}")

# ---- 入口の前縁の切り落とし (平板の起点で要る第 1 層 → 0、y1 ∝ s^0.1)
ic = int(np.searchsorted(s_rt, a.sc))
f_cl = f.copy(); f_cl[:ic] = f[ic]
log(f"[{a.tag}] 前縁の切り落とし s < {a.sc} r_t ({ic} 点): f(s_c) {f[ic]:.3e}、f(第 2 点) {f[1]:.3e}")

# ---- 表の x 節点: x < 10 は 0.125 r_t、以降 1 r_t、両端は x_in・x_e ちょうど
tx = np.concatenate([np.arange(x_in, min(10.0, x_e), 0.125), np.arange(10.0, x_e, 1.0)]) if x_e > 10 else np.arange(x_in, x_e, 0.125)
tx = np.unique(np.concatenate([tx[tx < x_e - 1e-6], [x_e]]))
def env(vals, mode):
    out = np.empty(len(tx))
    for k in range(len(tx)):
        lo = tx[max(k - 1, 0)]; hi = tx[min(k + 1, len(tx) - 1)]
        m = (x >= lo - 1e-12) & (x <= hi + 1e-12)
        out[k] = vals[m].min() if mode == "min" else vals[m].max()
    return out
tv = env(f_cl, "min") * a.target_scale         # 第 1 層の表 (下包絡: log 線形の補間が f_cl を超えない)
if a.frac_table_json:
    _t = np.asarray(json.load(open(a.frac_table_json)), dtype=float)
    tv = np.exp(np.interp(tx, _t[:, 0], np.log(_t[:, 1])))
    log(f"[{a.tag}] 第 1 層の表を {a.frac_table_json} から ({len(_t)} 点 → {len(tx)} 点に写す)")
fr_tab = lambda xs: np.exp(np.interp(xs, tx, np.log(tv)))   # noqa: E731  (mesh2d._first_frac_profile と同じ log 線形)
ratio = fr_tab(x) / (f_cl * a.target_scale)
log(f"[{a.tag}] 第 1 層の表 {len(tx)} 点: 表/要求 の最大 {ratio.max():.4f}・最小 {ratio.min():.4f}・平均 {ratio.mean():.4f}")

def stations(ni, refine, width, dtab=None):
    return M2._x_stations(x_in, x_e, int(ni), refine, width, density_table=dtab)

def wall_ar(xs, frs):
    r = np.asarray(wall.r(xs), dtype=float)
    y1 = frs * r
    ch = np.hypot(np.diff(xs), np.diff(r))
    e = np.stack([ch, y1[:-1], y1[1:]])
    return e.max(axis=0) / e.min(axis=0)

def ni_min(ok):
    lo, hi = 2, 64
    while not ok(hi):
        lo, hi = hi, hi * 2
        if hi > 10 ** 6:
            raise RuntimeError("ni > 1e6")
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if ok(mid):
            hi = mid
        else:
            lo = mid
    return hi

def y1p_est(xs, frs):
    """生成器と同じ目安: 0.8 × (配点の第 1 層比) / (y1+ 0.8 に要る第 1 層比)。入口の節点を除く (切り落としの無い f で)。"""
    v = W.Y1PLUS_TARGET * frs[1:] / np.interp(xs[1:], x, f)
    k = int(np.argmax(v))
    return float(v[k]), float(xs[1:][k])

def dx_ratio_max(xs):
    d = np.diff(xs); q = d[1:] / d[:-1]
    k = int(np.argmax(np.maximum(q, 1 / q)))
    return float(max(q[k], 1 / q[k])), float(xs[k + 1])

res = {"tag": a.tag, "r_t_m": rt, "x_in": x_in, "x_e": x_e, "L_U": LU, "nj": a.nj, "ar_limit": ARM, "s_c_rt": a.sc,
       "target_scale": a.target_scale, "y1plus_target_generator": W.Y1PLUS_TARGET, "table_points": len(tx)}

# ---- 区間ごとの y1,req
segs = [("入口の配管", x_in, -LU if LU else -2.0), ("縮流部 (x < -2)", -LU if LU else x_in, -2.0), ("スロート付近 |x| ≤ 2", -2.0, 2.0),
        ("拡大部の前半 2 < x ≤ x_e/2", 2.0, 0.5 * x_e), ("拡大部の後半 x > x_e/2", 0.5 * x_e, x_e)]
rows = []
for nm, lo, hi in segs:
    m = (x >= lo) & (x <= hi)
    m_noLE = m & (s_rt >= a.sc)
    rows.append({"seg": nm, "x_lo": lo, "x_hi": hi, "rw_min": float(rw[m].min()), "rw_max": float(rw[m].max()),
                 "y1_um_min": float(y1m[m_noLE].min() * 1e6), "y1_um_max": float(y1m[m].max() * 1e6),
                 "y1_um_median": float(np.median(y1m[m]) * 1e6),
                 "frac_min": float(f[m_noLE].min()), "frac_max": float(f[m].max()),
                 "x_at_frac_min": float(x[m_noLE][np.argmin(f[m_noLE])]),
                 "Tw_K": float(np.median(Tw[m]))})
res["segments"] = rows
for r in rows:
    log(f"[{a.tag}] {r['seg']:<28s} x [{r['x_lo']:7.2f},{r['x_hi']:7.2f}] r_w {r['rw_min']:.3f}–{r['rw_max']:.3f}  "
        f"y1,req {r['y1_um_min']:.3f}–{r['y1_um_max']:.3f} µm (中央 {r['y1_um_median']:.3f})  y1/r_w {r['frac_min']:.3e}–{r['frac_max']:.3e} "
        f"(最小 @ x {r['x_at_frac_min']:.2f})")

# ---- (A) 今のやり方 (2 値 + ブレンド; estimate_mesh が選んだ上流のブレンド)
frac, frac_t = est["wall_first_frac"], est["wall_first_frac_throat"]
mesh_c = dict(mesh, wall_first_up_x0=est["wall_first_up_x0"], wall_first_up_x1=est["wall_first_up_x1"], throat_refine=a.refine,
              throat_width=a.width)
if not a.skipA and frac is not None:
    niA = W._ni_for_ar(wall, mesh_c, frac, frac_t, ARM)
    arA, xA = W._wall_ar_fast(wall, niA, mesh_c, frac, frac_t)
    xsA = stations(niA, a.refine, a.width)
    frA = W._wall_first_profile(xsA, mesh_c, frac, frac_t)
    ypA = y1p_est(xsA, frA)
    res["A"] = {"ni": niA, "nodes": niA * a.nj, "ar_wall_max": arA, "x_ar_max": xA, "frac": frac, "frac_t": frac_t,
                "blend": est.get("blend_choice"), "up_x0": est["wall_first_up_x0"], "up_x1": est["wall_first_up_x1"],
                "y1plus_est_max": ypA[0], "x_y1plus_est_max": ypA[1], "dx_ratio_max": dx_ratio_max(xsA)}
    log(f"[{a.tag}] (A) ni {niA}  節点 {niA * a.nj:,}  壁際 AR 最大 {arA:.1f} @ x {xA:.3f}  y1+ 目安最大 {ypA[0]:.3f} @ {ypA[1]:.2f}")
    if a.buildA:
        t0 = time.time()
        niA10 = int(math.ceil(niA / 10) * 10)
        full = W._mesh_ar_full(wall, niA10, mesh_c, frac, frac_t)
        res["A"]["full_mesh_ni10"] = {"ni": niA10, **full}
        log(f"[{a.tag}] (A) 実格子 ni {niA10}: AR 最大 {full['ar_max']:.1f} @ x {full['x_at_max']:.3f} ({time.time() - t0:.0f}s)")

# ---- (B) 第 1 層の表 + 今の流れ方向の分布 (throat_refine・throat_width)
okB = lambda ni: wall_ar(xs := stations(ni, a.refine, a.width), fr_tab(xs)).max() <= ARM   # noqa: E731
niB = ni_min(okB)
xsB = stations(niB, a.refine, a.width)
arB = wall_ar(xsB, fr_tab(xsB)); kB = int(np.argmax(arB))
ypB = y1p_est(xsB, fr_tab(xsB))
res["B"] = {"ni": niB, "nodes": niB * a.nj, "ar_wall_max": float(arB[kB]), "x_ar_max": float(0.5 * (xsB[kB] + xsB[kB + 1])),
            "y1plus_est_max": ypB[0], "x_y1plus_est_max": ypB[1], "dx_ratio_max": dx_ratio_max(xsB)}
log(f"[{a.tag}] (B) ni {niB}  節点 {niB * a.nj:,}  壁際 AR 最大 {arB[kB]:.1f} @ x {res['B']['x_ar_max']:.3f}  y1+ 目安最大 {ypB[0]:.3f} @ {ypB[1]:.2f}")
# (B) で AR を決めている場所の内訳 (区間ごとに必要な ni)
def ni_needed_region(lo, hi, refine, width, dtab=None):
    def ok(ni):
        xs = stations(ni, refine, width, dtab)
        ar = wall_ar(xs, fr_tab(xs)); xm = 0.5 * (xs[1:] + xs[:-1])
        m = (xm >= lo) & (xm <= hi)
        return (not m.any()) or ar[m].max() <= ARM
    return ni_min(ok)
res["B"]["ni_by_segment"] = {r["seg"]: ni_needed_region(r["x_lo"], r["x_hi"], a.refine, a.width) for r in rows}
log(f"[{a.tag}] (B) 区間ごとに要る ni: {res['B']['ni_by_segment']}")

# ---- (C) 第 1 層の表 + 密度の表 dx ≤ ARM·y1·cosθ、dx ≤ (ni_ref の雛形の配点の間隔)
dens_ref = 1.0 + (a.refine - 1.0) * np.exp(-((x / a.width) ** 2))
I_ref = float(np.trapezoid(dens_ref, x))
h_ref_local = I_ref / ((a.ni_ref - 1) * dens_ref)          # ni_ref 列の雛形の配点の局所の間隔 [r_t]
h_ref_glob = float(h_ref_local.max())
h_AR = ARM * fr_tab(x) * rw / np.sqrt(1.0 + drw ** 2)
res["C_inputs"] = {"h_ref_local_min": float(h_ref_local.min()), "h_ref_local_max": h_ref_glob, "h_AR_min": float(h_AR.min()),
                   "x_h_AR_min": float(x[np.argmin(h_AR)])}
log(f"[{a.tag}] (C) ni_ref {a.ni_ref} の間隔 {h_ref_local.min():.4f}–{h_ref_glob:.4f} r_t、AR からの間隔の最小 {h_AR.min():.5f} r_t @ x {x[np.argmin(h_AR)]:.2f}")

def run_C(label, h_cap):
    h = np.minimum(h_AR, h_cap)
    d = 1.0 / h
    ni_ideal = 1 + float(np.trapezoid(d, x))
    dv = env(d, "max")                                     # 密度の表 (上包絡: 線形補間が 1/h を下回らない)
    dtab = np.stack([tx, dv / dv.min()], axis=1)
    hcap_of = lambda xm: np.interp(xm, x, h_cap if np.ndim(h_cap) else np.full_like(x, h_cap))   # noqa: E731
    def ok(ni):
        xs = stations(ni, 1.0, 1.0, dtab)
        if wall_ar(xs, fr_tab(xs)).max() > ARM:
            return False
        xm = 0.5 * (xs[1:] + xs[:-1])
        return bool(np.all(np.diff(xs) <= hcap_of(xm) * (1 + 1e-9)))
    ni = ni_min(ok)
    xs = stations(ni, 1.0, 1.0, dtab)
    ar = wall_ar(xs, fr_tab(xs)); k = int(np.argmax(ar))
    yp = y1p_est(xs, fr_tab(xs))
    xm = 0.5 * (xs[1:] + xs[:-1]); dx = np.diff(xs)
    bind_AR = dx >= 0.98 * ARM * np.minimum(fr_tab(xs)[:-1] * wall.r(xs)[:-1], fr_tab(xs)[1:] * wall.r(xs)[1:]) / np.hypot(1.0, np.diff(wall.r(xs)) / dx)
    out = {"label": label, "ni": ni, "nodes": ni * a.nj, "ni_ideal": ni_ideal, "ar_wall_max": float(ar[k]),
           "x_ar_max": float(0.5 * (xs[k] + xs[k + 1])), "y1plus_est_max": yp[0], "x_y1plus_est_max": yp[1],
           "dx_ratio_max": dx_ratio_max(xs), "dx_min_rt": float(dx.min()), "dx_max_rt": float(dx.max()),
           "x_dx_min": float(xm[np.argmin(dx)]),
           "frac_AR_bound_stations": float(bind_AR.mean()),
           "x_density_table": dtab.tolist()}
    # 区間ごとの列の数
    out["stations_by_segment"] = {r["seg"]: int(((xs >= r["x_lo"]) & (xs < r["x_hi"])).sum()) for r in rows}
    log(f"[{a.tag}] (C{label}) ni {ni} (理想 {ni_ideal:.0f})  節点 {ni * a.nj:,}  壁際 AR 最大 {ar[k]:.1f} @ x {out['x_ar_max']:.3f}  "
        f"y1+ 目安最大 {yp[0]:.3f} @ {yp[1]:.2f}  dx 比最大 {out['dx_ratio_max'][0]:.5f} @ {out['dx_ratio_max'][1]:.2f}  "
        f"dx {dx.min():.5f}–{dx.max():.4f} r_t  AR 律速の列 {100 * out['frac_AR_bound_stations']:.0f} %")
    log(f"[{a.tag}] (C{label}) 区間ごとの列: {out['stations_by_segment']}")
    return out, dtab

res["C_local"], dtab_local = run_C("local", h_ref_local)
res["C_global"], _ = run_C("global", h_ref_glob)
res["first_frac_table"] = np.stack([tx, tv], axis=1).tolist()

if a.build:
    t0 = time.time()
    c = res["C_local"]
    prm = M2.Mesh2DParams(ni=c["ni"], nj=a.nj, wall_first_frac=float(tv.max()), throat_refine=1.0, throat_width=1.0, scale=1.0,
                          wall_first_frac_table=res["first_frac_table"], x_density_table=dtab_local.tolist())
    coords, quads, _ = M2.generate_axisym_mesh(wall, prm)
    P = coords[:, :2][quads]
    e = np.linalg.norm(np.roll(P, -1, axis=1) - P, axis=2)
    ar = e.max(axis=1) / e.min(axis=1)
    k = int(np.argmax(ar))
    X = coords[:, 0].reshape(c["ni"], a.nj); R = coords[:, 1].reshape(c["ni"], a.nj)
    xs = X[:, 0]
    g = np.diff(R, axis=1)                    # 軸 → 壁
    q = g[:, :-1] / g[:, 1:]                  # 隣の (軸側 / 壁側) の比
    ar_wall_row = ar.reshape(c["ni"] - 1, a.nj - 1)[:, -1]
    first = (R[:, -1] - R[:, -2]) / R[:, -1]
    res["C_local"]["built"] = {"n_nodes": int(coords.shape[0]), "n_cells": int(len(quads)), "ar_max": float(ar[k]),
                               "x_ar_max": float(P[k, :, 0].mean()), "r_ar_max_over_rw": float(P[k, :, 1].mean() / wall.r(P[k, :, 0].mean())),
                               "ar_wall_row_max": float(ar_wall_row.max()), "n_cells_ar_gt_1000": int((ar > 1000).sum()),
                               "dx_ratio_max": dx_ratio_max(xs), "radial_growth_max": float(np.max(np.maximum(q, 1 / q))),
                               "first_frac_built_min": float(first.min()), "first_frac_built_max": float(first.max()),
                               "time_s": time.time() - t0}
    b = res["C_local"]["built"]
    log(f"[{a.tag}] (C local) 実格子 {b['n_nodes']:,} 点: AR 最大 {b['ar_max']:.1f} @ x {b['x_ar_max']:.3f} (r/r_w {b['r_ar_max_over_rw']:.4f})  "
        f"壁の列の AR 最大 {b['ar_wall_row_max']:.1f}  AR>1000 のセル {b['n_cells_ar_gt_1000']:,}  dx 比最大 {b['dx_ratio_max'][0]:.5f} @ {b['dx_ratio_max'][1]:.2f}  "
        f"径方向の比の最大 {b['radial_growth_max']:.4f}  ({b['time_s']:.0f}s)")
    np.savez_compressed(D / f"{a.tag}_C_local_stations.npz", xs=xs, first_frac=first, rw=R[:, -1])

json.dump(res, open(D / f"{a.tag}_result.json", "w"), ensure_ascii=False, indent=1)
# y1,req の断面ごとの値 (間引き) を CSV に
with open(D / f"{a.tag}_y1req.csv", "w") as fh:
    fh.write("x_rt,x_m,s_rt,r_w_rt,Tw_K,y1_req_m,y1_req_over_rw,first_frac_table\n")
    for i in range(0, len(x), 10):
        fh.write(f"{x[i]:.6f},{x[i] * rt:.6f},{s_rt[i]:.6f},{rw[i]:.6f},{Tw[i]:.3f},{y1m[i]:.6e},{f[i]:.6e},{fr_tab(np.array([x[i]]))[0]:.6e}\n")
log(f"[{a.tag}] 書いた: {a.tag}_result.json・{a.tag}_y1req.csv")
