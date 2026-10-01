#!/usr/bin/env python3
r"""case/64 (A) と case/65 (C) の評価器 (plan `boundary-cht-conjugate-benchmarks.md` §4.1・§4.6・§4.7)。

forge の run の最終スナップショットの流れ場 (ρ・u・v・p) を固定し、`conjugate_ref.Problem` で流体と固体の温度を独立に解き直した
**主参照**と、forge の界面温度 T_i・界面熱流束 q_i (= −`iface_q_eff`、流体へ向かう熱を正) を壁節点で比べる。

主参照の不確かさ U (§4.6) = (a) 参照格子の細分化 (forge 格子を 1 回・2 回 2 等分した 2 水準の差) + (b) 写像 (双線形 ↔ 3 次) の差
+ (c) 散逸・圧力仕事の微分 (2 次 ↔ 3 次スプライン) の差 + (e) 界面熱流束の取り出し (固体側片側差分 2 次 ↔ 3 次) の差。
(d) 上流を 1.5 倍に延ばした参照の変化は **U に入れず別掲** する (2026-09-30 disposition M2 (ii)、plan §4.6 の事後改訂):
forge と参照は同一の有限領域・BC を解いており、延長は問題そのもの (入口位置・固体端面・源項の範囲) を変える。
A の温度は加熱区間・予熱域に加えて**全長の壁節点**でも判定する (M2)。熱流束は登録どおり加熱区間と予熱域 (−40R ≤ x < 0) の別判定。
各項目の判定を**すべて**表示し、総合は FAIL が 1 つでもあれば「FAIL (一部判定不能)」のように併記する (M3)。
合否は |観測差| + U ≤ 許容、U > 許容/3 なら判定不能 (終了コード 2)。

    python3 eval_conj.py A <run>          # case/64
    python3 eval_conj.py C <run>          # case/65
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent / "65.conjugate_flat_plate"))
import conjugate_ref as cr  # noqa: E402

trap = getattr(np, "trapezoid", None) or np.trapz


def refuse(msg):
    print(f"REFUSED: {msg}\nVERDICT: REFUSED"); sys.exit(2)


def last_step(run):
    s = [int(m.group(1)) for p in run.glob("res_[0-9]*.h5") for m in [re.match(r"res_(\d+)\.h5$", p.name)] if m]
    if not s: refuse("res_*.h5 が無い")
    return max(s)


def load(run, st, wall_name):
    with h5py.File(run / "mesh.h5", "r") as m:
        xyz = np.asarray(m["MESH/COORD"][:], float).reshape(-1, 3)
    with h5py.File(run / f"res_{st}.h5", "r") as h:
        F = {k: np.asarray(h["VALUE"][k][:], float) for k in ("ro", "Ux", "Uy", "P", "T")}
    for k, v in F.items():
        if not np.isfinite(v).all(): refuse(f"{k} に非有限値")
    xr, yr = np.round(xyz[:, 0], 10), np.round(xyz[:, 1], 10)
    xs, ys = np.unique(xr), np.unique(yr)
    if len(xs) * len(ys) != len(xyz): refuse("流体の節点がテンソル格子でない")
    ix = np.searchsorted(xs, xr); iy = np.searchsorted(ys, yr)
    G = {}
    for k, v in F.items():
        a = np.full((len(xs), len(ys)), np.nan); a[ix, iy] = v; G[k] = a
    with h5py.File(run / f"res_{wall_name}_{st}.h5", "r") as w:
        c = np.asarray(w["MESH/COORD"][:], float).reshape(-1, 3)
        q = -np.asarray(w["VALUE/iface_q_eff"][:], float); tw = np.asarray(w["VALUE/iface_Tw_bc"][:], float)
        ok = np.asarray(w["VALUE/iface_ok"][:], float)
    o = np.argsort(c[:, 0])
    return np.asarray(xs, float), np.asarray(ys, float), G, dict(x=c[o, 0], y=c[o, 1], q=q[o], Tw=tw[o], ok=ok[o])


def load_fields_only(run, st):
    """合成試験用: 壁ダンプ無しで流体場だけ読む。"""
    with h5py.File(run / "mesh.h5", "r") as m:
        xyz = np.asarray(m["MESH/COORD"][:], float).reshape(-1, 3)
    with h5py.File(run / f"res_{st}.h5", "r") as h:
        F = {k: np.asarray(h["VALUE"][k][:], float) for k in ("ro", "Ux", "Uy", "P", "T")}
    xr, yr = np.round(xyz[:, 0], 10), np.round(xyz[:, 1], 10)
    xs, ys = np.unique(xr), np.unique(yr); ix = np.searchsorted(xs, xr); iy = np.searchsorted(ys, yr)
    G = {}
    for k, v in F.items():
        a = np.full((len(xs), len(ys)), np.nan); a[ix, iy] = v; G[k] = a
    return xs, ys, G, None


def bisect(a, k):
    for _ in range(k):
        a = np.sort(np.concatenate([a, 0.5 * (a[1:] + a[:-1])]))
    return a


def build(case, xs_f, ys_f, G, k, method="linear", deriv="o2", pc=None, extend=False, solid_axial=1.0):
    """forge の流体格子 (xs_f, ys_f) を k 回 2 等分した参照格子に固体を足して Problem を作る。"""
    gc = pc.gc
    xs = bisect(xs_f, k); yfl = bisect(ys_f, k)
    n_ext = 0
    if extend and case == "A":                                        # 登録 (§4.6 (d)): 上流を 1.5 倍 (−80R → −120R) に延ばす
        dx0 = xs[1] - xs[0]
        xe = np.arange(xs[0] - dx0, -120 * pc.R - 1e-12, -dx0)[::-1]
        n_ext = len(xe); xs = np.concatenate([xe, xs])
    ymap, n_ext_y = yfl, 0
    if extend and case == "C":                                        # 登録 (§4.6 (d)): 上境界を 1.5 倍に (後継 plan §4.6.2 M2、A と同じ作法)
        dy = yfl[-1] - yfl[-2]
        ye = np.arange(yfl[-1] + dy, 1.5 * pc.H_TOP + 1e-12, dy)
        n_ext_y = len(ye); yfl = np.concatenate([yfl, ye])
    if case == "A":
        ns_f = len(ys_f) - 1                                           # 固体は流体と同じ半径間隔 (gen_solid の既定)
        ysol = np.linspace(pc.R, pc.R_O, ns_f * 2 ** k + 1)[1:]
        ys = np.concatenate([yfl, ysol]); jw = len(yfl) - 1
        mat = np.full((len(xs) - 1, len(ys) - 1), cr.FLUID); mat[:, jw:] = cr.SOLID
        axisym = True; k_s = pc.KS_RATIO[pc.CASE] * gc.K_F
        rob = lambda xm, y: (pc.H_O, gc.T_IN + pc.DT_C) if (-1e-12 <= xm <= pc.L_HEAT + 1e-9 and abs(y - ys[-1]) < 1e-12) else None
        fl_rows = slice(0, jw + 1)
    else:
        ns_f = len(ys_f) - 1
        ysol = np.linspace(-pc.B, 0.0, ns_f * 2 ** k + 1)[:-1]
        ys = np.concatenate([ysol, yfl]); jw = len(ysol)
        mat = np.full((len(xs) - 1, len(ys) - 1), cr.FLUID)
        for i in range(len(xs) - 1):
            xm = 0.5 * (xs[i] + xs[i + 1]); mat[i, :jw] = cr.SOLID if (0 <= xm <= pc.L) else cr.VOID
        axisym = False; k_s = pc.KS_RATIO[pc.CASE] * gc.K_F
        rob = lambda xm, y: (pc.H_BACK, gc.T_IN + pc.DT_H) if (0 <= xm <= pc.L and abs(y - ys[0]) < 1e-15) else None
        fl_rows = slice(jw, len(ys))
    nx, ny = len(xs), len(ys)
    Fm = {kk: cr.map_field(xs_f, ys_f, G[kk], xs[n_ext:], ymap, method) for kk in ("ro", "Ux", "Uy", "P", "T")}
    uy_top = None
    if n_ext_y:
        # 延長部の固定流れ場: forge の上端の行の ρ・u・v・P・T をそのまま y 方向に繰り返す (源項は延長後の場から同じ flow_source で作る)
        uy_top = float(np.abs(Fm["Uy"][:, -1]).max())
        for kk in Fm:
            Fm[kk] = np.concatenate([Fm[kk], np.repeat(Fm[kk][:, -1:], n_ext_y, axis=1)], axis=1)
    if n_ext:
        # 延長部の固定流れ場: forge の入口列の ρ・u・v・T をそのまま延ばし、圧力は入口の勾配で線形に延ばす (十分発達した流れ)
        dpdx = (Fm["P"][1] - Fm["P"][0]) / (xs[n_ext + 1] - xs[n_ext])
        for kk in ("ro", "Ux", "Uy", "T"):
            Fm[kk] = np.concatenate([np.repeat(Fm[kk][:1], n_ext, axis=0), Fm[kk]])
        Fm["P"] = np.concatenate([Fm["P"][0][None, :] + dpdx[None, :] * (xs[:n_ext, None] - xs[n_ext]), Fm["P"]])
    ro = np.ones((nx, ny)); u = np.zeros((nx, ny)); v = np.zeros((nx, ny)); S = np.zeros((nx, ny))
    ro[:, fl_rows] = Fm["ro"]; u[:, fl_rows] = Fm["Ux"]; v[:, fl_rows] = Fm["Uy"]
    phi, work = cr.flow_source(xs, yfl, Fm["ro"], Fm["Ux"], Fm["Uy"], Fm["P"], gc.MU, axisym, deriv=deriv)
    S[:, fl_rows] = phi + work
    P = cr.Problem(xs, ys, mat, gc.K_F, k_s, gc.CP, axisym, 0.0, ro=ro, u=u, v=v, source=S, robin=rob, solid_axial=solid_axial)
    P.uy_top = uy_top
    # 入口は forge の入口列の温度分布 (流体部分) を Dirichlet に: T_in を行ごとに
    Tin_col = np.full(ny, gc.T_IN); Tin_col[fl_rows] = Fm["T"][0]
    P.T_in_profile = Tin_col
    return P, xs, ys, jw, k_s


def _solve_rowwise(P, prof):
    """入口 Dirichlet を行ごとの値 (forge の入口列の温度) にして解く。"""
    P.inlet_profile = prof
    return P.solve()


def interface_q(P, T, jw, k_s, case, order):
    ys = P.ys
    if case == "A":                                                # 固体は r > R (上側)。流体へ = k_s ∂T/∂r|_{R+}
        h = ys[jw + 1] - ys[jw]
        d = (-3 * T[:, jw] + 4 * T[:, jw + 1] - T[:, jw + 2]) / (2 * h) if order == 2 else \
            (-11 * T[:, jw] + 18 * T[:, jw + 1] - 9 * T[:, jw + 2] + 2 * T[:, jw + 3]) / (6 * h)
        return k_s * d
    h = ys[jw] - ys[jw - 1]                                        # 固体は y < 0 (下側)。流体へ = −k_s ∂T/∂y|_{0−}
    d = (3 * T[:, jw] - 4 * T[:, jw - 1] + T[:, jw - 2]) / (2 * h) if order == 2 else \
        (11 * T[:, jw] - 18 * T[:, jw - 1] + 9 * T[:, jw - 2] - 2 * T[:, jw - 3]) / (6 * h)
    return -k_s * d


def metrics(case, pc, xs, T, q, jw, xw=None):
    gc = pc.gc
    Ti = T[:, jw]
    if case == "A":
        heat = (xs >= -1e-12) & (xs <= pc.L_HEAT + 1e-9)
        rise = trap(Ti[heat] - gc.T_IN, xs[heat]) / pc.L_HEAT
        return dict(rise=rise, Ti=Ti, q=q)
    win = (xs >= 0.2 * pc.L - 1e-12) & (xs <= 0.9 * pc.L + 1e-12)
    return dict(theta=(Ti - gc.T_IN) / pc.DT_H, q=q, qmean=trap(q[win], xs[win]) / (xs[win][-1] - xs[win][0]))


def expected_wall_x(case, xs_f, ys_f, pc):
    """メッシュから決まる評価対象の壁節点の x (A: r = R の全節点、C: y = 0 かつ 0 ≤ x ≤ L)。"""
    if case == "A":
        return xs_f.copy()
    return xs_f[(xs_f >= -1e-12) & (xs_f <= pc.L + 1e-9)]


def solid_from_dump(run, st, pid, k_s):
    """forge の固体ダンプ (テンソル格子の帯) から (x, y, T[nx, ny]) と外面 Robin 入熱 |Σ q_hole|。"""
    with h5py.File(run / "solid.h5", "r") as f:
        C = np.asarray(f["MESH/COORD"][:], float)
    with h5py.File(run / f"res_solid_{pid}_{st}.h5", "r") as f:
        Tn = np.asarray(f["VALUE/T"][:], float); qh = float(np.asarray(f["VALUE/q_hole"][:], float).sum())
    if not (np.isfinite(Tn).all() and np.isfinite(qh)): refuse("固体ダンプに非有限値")
    xr, yr = np.round(C[:, 0], 10), np.round(C[:, 1], 10)
    xs, ys = np.unique(xr), np.unique(yr)
    if len(xs) * len(ys) != len(C): refuse("固体がテンソル格子でない")
    Tg = np.full((len(xs), len(ys)), np.nan); Tg[np.searchsorted(xs, xr), np.searchsorted(ys, yr)] = Tn
    if np.isnan(Tg).any(): refuse("固体の温度を格子に並べられない")
    return xs, ys, Tg, abs(qh)


def solid_indicators(case, pc, xs, ys_s, Ts, k_s, Qtot):
    """固体の効果: 加熱区間 (A) / 評価窓 (C) の厚さ方向の温度差の最大と、軸方向熱量の最大 / Q_tot。
    ys_s・Ts は固体部分だけ (A: r = R → r_o、C: y = −b → 0)。"""
    Tx = np.gradient(Ts, xs, axis=0)
    w = ys_s if case == "A" else np.ones_like(ys_s)
    Qax = np.array([-trap(k_s * Tx[i] * w, ys_s) for i in range(len(xs))])
    if case == "A":
        sel = (xs >= -1e-12) & (xs <= pc.L_HEAT + 1e-9)
    else:
        sel = (xs >= 0.2 * pc.L - 1e-12) & (xs <= 0.9 * pc.L + 1e-12)
    dTs = np.abs(Ts[:, -1] - Ts[:, 0])[sel].max()
    return float(dTs), float(np.abs(Qax[sel]).max() / Qtot)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case", choices=["A", "C"]); ap.add_argument("run")
    ap.add_argument("--levels", type=int, default=2, help="参照格子の細分回数の最大 (登録 2 → 3 水準: 0,1,2)")
    a = ap.parse_args()
    run = Path(a.run)
    txt = (run / "RUN_INPUTS.txt").read_text()
    m = re.search(r"kind cht (A1|A2|C1|C2)", txt)
    if not m: refuse("RUN_INPUTS.txt に kind cht A1/A2/C1/C2 が無い")
    cname = m.group(1)
    if a.case == "A":
        import pipe_common as pc
    else:
        import plate_common as pc
    pc.CASE = cname
    gc = pc.gc
    pid = 3 if a.case == "A" else 5
    st = last_step(run)
    xs_f, ys_f, G, W = load(run, st, "wall_3" if a.case == "A" else "plate_5")
    xw = W["x"]
    exp = expected_wall_x(a.case, xs_f, ys_f, pc)
    if len(xw) != len(np.unique(np.round(xw, 10))): refuse("壁ダンプの節点に重複")
    if len(xw) != len(exp) or not np.allclose(xw, exp, rtol=0, atol=1e-9): refuse(f"壁ダンプの節点 ({len(xw)}) が期待集合 ({len(exp)}) と一致しない")
    if (W["ok"] != 1).any(): refuse("iface_ok = 0 の壁節点がある")
    for kk in ("q", "Tw"):
        if not np.isfinite(W[kk]).all(): refuse(f"壁ダンプの {kk} に非有限値")
    k_s_reg = pc.KS_RATIO[cname] * gc.K_F
    # ---- 参照の各変種 (§4.6)
    variants = {}
    def run_variant(key, **kw):
        P, xs, ys, jw, k_s = build(a.case, xs_f, ys_f, G, kw.pop("k"), pc=pc, **kw)
        T = _solve_rowwise(P, P.T_in_profile)
        variants[key] = dict(P=P, xs=xs, ys=ys, jw=jw, k_s=k_s, T=T, q2=interface_q(P, T, jw, k_s, a.case, 2),
                             q3=interface_q(P, T, jw, k_s, a.case, 3), Qtot=abs(P.robin_heat()))
    for k in range(a.levels + 1):
        run_variant(("grid", k), k=k)
    run_variant("cubic", k=1, method="cubic"); run_variant("spline", k=1, deriv="spline")
    run_variant("ext", k=1, extend=True)
    if a.case == "C":                                                  # 後継 plan §4.6.2 M1: 軸方向伝導なしの参照 (solid_axial=0)
        run_variant(("noax", a.levels), k=a.levels, solid_axial=0.0); run_variant(("noax", a.levels - 1), k=a.levels - 1, solid_axial=0.0)
        run_variant(("noax", 1), k=1, solid_axial=0.0)
        run_variant("noax_cubic", k=1, method="cubic", solid_axial=0.0); run_variant("noax_spline", k=1, deriv="spline", solid_axial=0.0)

    def wall_vals(key, qkey="q2"):
        V = variants[key]; xs = V["xs"]
        idx = [int(np.argmin(np.abs(xs - x))) for x in xw]
        if max(abs(xs[i] - x) for i, x in zip(idx, xw)) > 1e-9: refuse("参照格子に forge の壁節点の x が無い")
        return V["T"][idx, V["jw"]], V[qkey][idx]

    def ratio(q, Qtot):
        up = xw <= 0
        return pc.R * trap(q[up], xw[up]) / Qtot

    def solid_ind(key):
        V = variants[key]
        sl = slice(V["jw"], None) if a.case == "A" else slice(0, V["jw"] + 1)
        return solid_indicators(a.case, pc, V["xs"], V["ys"][sl], V["T"][:, sl], V["k_s"], V["Qtot"])

    L = a.levels
    Tr, qr = wall_vals(("grid", L)); TrP, qrP = wall_vals(("grid", L - 1)); T1, q1 = wall_vals(("grid", 1))
    Tc_, qc_ = wall_vals("cubic"); Ts_, qs_ = wall_vals("spline"); _, q3 = wall_vals(("grid", L), "q3")
    uT = np.abs(Tr - TrP) + np.abs(Tc_ - T1) + np.abs(Ts_ - T1)
    uq = np.abs(qr - qrP) + np.abs(qc_ - q1) + np.abs(qs_ - q1) + np.abs(q3 - qr)
    # (d) 上流延長 (−120R) の差は**比較の U に入れない** (2026-09-30 disposition M2 (ii) — 事後改訂): forge と参照は同一の有限領域・
    # 同一 BC の問題を解いており、延長は入口位置・固体端面・源項の範囲を変える「問題定義の感度」。別掲する。
    Te_, qe_ = wall_vals("ext")
    ext_sens = (np.abs(Te_ - T1), np.abs(qe_ - q1))
    Tf, qf = W["Tw"], W["q"]
    # forge の固体
    xs_s, ys_s, Ts_f, Qtot_f = solid_from_dump(run, st, pid, k_s_reg)
    dTs_f, Qax_f = solid_indicators(a.case, pc, xs_s, ys_s if a.case == "A" else ys_s, Ts_f, k_s_reg, Qtot_f)
    rows = []
    if a.case == "A":
        heat = (xw >= -1e-12) & (xw <= pc.L_HEAT + 1e-9); pre = (xw < 0) & (xw >= -40 * pc.R)
        base = variants[("grid", L)]
        rise = trap(Tr[heat] - gc.T_IN, xw[heat]) / pc.L_HEAT
        qo = base["Qtot"] / (pc.R_O * pc.L_HEAT)
        def region(mask, nm, with_q=True):
            rows.append((f"温度 (T_i−T_in) {nm}: max |Δ|/上昇", (np.abs(Tf - Tr)[mask] / rise).max(), (uT[mask] / rise).max(), 0.01))
            if with_q:
                rows.append((f"熱流束 q_i {nm}: max |Δ|/q_o", (np.abs(qf - qr)[mask] / qo).max(), (uq[mask] / qo).max(), 0.02))
        region(heat, "加熱区間"); region(pre, "予熱域 (−40R ≤ x < 0)"); region(np.ones_like(heat), "全長", with_q=False)
        # Q_up/Q_tot: forge は forge の総入熱、参照は各変種の総入熱で割る (M1)
        r_f = ratio(qf, Qtot_f)
        rv = {key: ratio(wall_vals(key)[1], variants[key]["Qtot"]) for key in variants}
        r_r = rv[("grid", L)]
        U_r = (abs(rv[("grid", L)] - rv[("grid", L - 1)]) + abs(rv["cubic"] - rv[("grid", 1)]) + abs(rv["spline"] - rv[("grid", 1)])
               + abs(ratio(q3, base["Qtot"]) - r_r))
        ext_ratio = abs(rv["ext"] - rv[("grid", 1)])
        rows.append(("上流へ回り込む熱 Q_up/Q_tot の差 (絶対)", abs(r_f - r_r), U_r, 0.005))
        info = (f"壁温上昇 (参照、長さ平均) {rise:.4f} K、q_o {qo:.2f} W/m²、Q_tot forge {Qtot_f:.6e} / 参照 {base['Qtot']:.6e} W/rad、"
                f"Q_up/Q_tot forge {r_f:.5f} / 参照 {r_r:.5f}")
        tolT_abs, tolI = 0.01 * rise, 0.005
    else:
        win = (xw >= 0.2 * pc.L - 1e-12) & (xw <= 0.9 * pc.L + 1e-12)
        th_f = (Tf - gc.T_IN) / pc.DT_H; th_r = (Tr - gc.T_IN) / pc.DT_H
        qm = trap(qr[win], xw[win]) / (xw[win][-1] - xw[win][0])
        rows.append(("界面温度 θ_i: 窓内 max |Δθ|", np.abs(th_f - th_r)[win].max(), (uT / pc.DT_H)[win].max(), 0.01))
        rows.append(("熱流束 q_i: 窓内 max |Δ|/窓内平均", (np.abs(qf - qr)[win] / qm).max(), (uq[win] / qm).max(), 0.03))
        info = f"窓内 θ_i 参照 {th_r[win].min():.4f}…{th_r[win].max():.4f}、q 平均 {qm:.2f} W/m²"
        tolT_abs, tolI = 0.01 * pc.DT_H, 0.03
    # 固体の効果 (§4.4・§6): 効果 ≥ 5U かつ効果 > 許容幅
    si = {key: solid_ind(key) for key in variants}
    dTs_r, Qax_r = si[("grid", L)]
    def U_of(j):
        u = abs(si[("grid", L)][j] - si[("grid", L - 1)][j]) + abs(si["cubic"][j] - si[("grid", 1)][j]) + abs(si["spline"][j] - si[("grid", 1)][j])
        return u
    eff = [("固体の厚さ方向の温度差 ΔT_s [K]", dTs_f, dTs_r, U_of(0), tolT_abs)]
    if a.case == "A":
        eff.append(("固体の軸方向熱量 max|Q_ax|/Q_tot", Qax_f, Qax_r, U_of(1), tolI))
    print(f"=== {a.case} ({cname}) {run}  step {st}、壁節点 {len(xw)} (期待集合と一致)、参照 {L + 1} 水準")
    print(f"  {info}")
    bad = und = False
    for nm, d, U, tol in rows:
        if not (np.isfinite(d) and np.isfinite(U)):
            v = "判定不能"; und = True
        elif U > tol / 3:
            v = "判定不能"; und = True
        elif d + U <= tol:
            v = "PASS"
        else:
            v = "FAIL"; bad = True
        print(f"  {v:5s} {nm:<42} 差 {d:.4e} + U {U:.4e}  (許容 {tol:g}、U 上限 {tol/3:.3g})")
    print("  --- 固体が効いていることの確認 (効果 ≥ 5U かつ 効果 > 許容幅)")
    for nm, vf, vr, U, tol in eff:
        good = np.isfinite(vf) and vf >= 5 * U and vf > tol
        # 発注元 plan §6: 「C1 の軸方向熱量は許容と同程度なので、軸方向伝導の判定は C2 だけで行う」(登録文どおり、C1 は参考表示)
        if cname == "C1" and "軸方向" in nm:
            print(f"  参考  {nm:<34} forge {vf:.4e} / 参照 {vr:.4e}、U {U:.3e} (効果/U {vf/max(U,1e-300):.1f})、許容幅 {tol:.3e} (効果/許容 {vf/tol:.1f}) — C1 は判定しない (§6)")
            continue
        bad |= not good
        print(f"  {'PASS' if good else 'FAIL'}  {nm:<34} forge {vf:.4e} / 参照 {vr:.4e}、U {U:.3e} (効果/U {vf/max(U,1e-300):.1f})、許容幅 {tol:.3e} (効果/許容 {vf/tol:.1f})")
    if a.case == "C":
        # 後継 plan §4.6.2 M1: 登録量 (発注元 §4.5 C2 行) = 評価窓内の軸方向伝導あり/なしの θ_i 差。C2 はゲート、C1 は参考
        def dth(ka, kn):
            Ta_, qa_ = wall_vals(ka); Tn_, qn_ = wall_vals(kn)
            return float((np.abs(Ta_ - Tn_)[win] / pc.DT_H).max()), float((np.abs(qa_ - qn_)[win]).max() / qm)
        dL, dqL = dth(("grid", L), ("noax", L)); dLm, _ = dth(("grid", L - 1), ("noax", L - 1)); d1, _ = dth(("grid", 1), ("noax", 1))
        dc, _ = dth("cubic", "noax_cubic"); ds, _ = dth("spline", "noax_spline")
        U_ax = abs(dL - dLm) + abs(dc - d1) + abs(ds - d1)
        Tn, _ = wall_vals(("noax", L)); th_n = (Tn - gc.T_IN) / pc.DT_H
        D_f = float(np.abs(th_f - th_n)[win].max()); U_th = float((uT / pc.DT_H)[win].max())
        good = np.isfinite(dL) and dL >= 5 * U_ax and dL > 0.01
        tag = ("PASS" if good else "FAIL") if cname == "C2" else "参考"
        if cname == "C2":
            bad |= not good
        print(f"  {tag:5s} 軸方向伝導あり/なしの θ_i 差 (参照)       Δθ_ax {dL:.4e}、U_ax {U_ax:.3e} (効果/U {dL/max(U_ax,1e-300):.1f})、許容幅 0.01 (効果/許容 {dL/0.01:.2f})"
              + ("" if cname == "C2" else " — C1 は判定しない (§6)"))
        print(f"  報告  forge と「なし」参照の距離 D_f {D_f:.4e} (主判定の U_θ {U_th:.3e}、D_f − U_θ {D_f - U_th:.4e})")
        print(f"  参考  あり/なしの q 差 max/q_mean {dqL:.4e}、固体の軸方向熱量 max|Q_ax|/Q_tot forge {Qax_f:.4e} / 参照 {Qax_r:.4e} (U {U_of(1):.3e})")
    if ext_sens is not None and a.case == "C":
        uyt = variants["ext"]["P"].uy_top
        print(f"  --- 別掲 (比較の U に含めない): 上境界を 1.5 倍に延ばしたときの参照の変化 (問題定義の感度、後継 plan §4.6.2 M2)")
        print(f"      θ_i 窓内 max {(ext_sens[0][win] / pc.DT_H).max():.4e} (許容の 1/3 = 3.3e-3)、q_i 窓内 max {(ext_sens[1][win]).max()/qm:.4e} of q_mean (許容の 1/3 = 1e-2)、"
              f"上端行の max|Uy| {uyt:.3e} m/s ({uyt/pc.U_INF:.2e} U∞)")
    if ext_sens is not None and a.case == "A":
        print(f"  --- 別掲 (比較の U に含めない): 上流を 1.5 倍に延ばしたときの参照の変化 (問題定義の感度)")
        print(f"      壁温 max {ext_sens[0].max():.4e} K ({ext_sens[0].max()/rise*100:.3f} % of 上昇)、q_i max {ext_sens[1].max()/qo*100:.3f} % of q_o、"
              f"Q_up/Q_tot {ext_ratio:.3e}")
    out = run / f"eval_conj_{st}_L{L}.csv"
    np.savetxt(out, np.c_[xw, Tf, Tr, qf, qr, uT, uq], delimiter=",", comments="",
               header=f"x,T_forge,T_ref,q_forge,q_ref,U_T,U_q  # levels {L}", fmt="%.10e")
    verdict = ("FAIL (一部判定不能)" if und else "FAIL") if bad else ("判定不能" if und else "PASS")
    print(f"  CSV {out.name}")
    print(f"VERDICT: {verdict}")
    return 1 if bad else (2 if und else 0)


if __name__ == "__main__":
    sys.exit(main())
