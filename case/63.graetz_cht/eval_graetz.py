#!/usr/bin/env python3
r"""case/63 Graetz の評価器 (plan `boundary-cht-axisymmetric-graetz.md` §4.5・§6)。

加熱 run と対照 run (ΔT = 0、同じ格子) から、各壁節点の差し引き局所 Nu を作り、古典 Graetz (march) と比べる。

    Nu_i = −(q_i − q0_i) D / ( k_f [ (Tw_i − Tb_i) − (Tw0_i − Tb0_i) ] )

- q = `iface_q_eff` (固体向き正。加熱では負)、Tw = `iface_Tw_bc`。`iface_ok` = 0 の節点が窓内にあれば FAIL。
- Tb は同じ x の節点列の混合平均 ∫ρ u_x T r dr / ∫ρ u_x r dr (Simpson。半径方向は一様格子)。節点の対応は `MESH/COORD` の x で取る。
- x⁺ = x / (D Pe_meas) (**原点は加熱開始点 x = 0**)、Pe_meas は実測の質量流束から Re·Pr で出す。
- 基準: march (nr 400) の q ∝ Nu θ_b を**各壁節点の双対面区間** [x_{i−1/2}, x_{i+1/2}] で平均し、節点の θ_b で割った Nu。

モード:
    snap   : 最終 (または --step) スナップショットで V-g1 (流れ・入口・軸) と V-g2 (Nu の誤差) を判定し、節点ごとの CSV を書く
    series : 両 run に共通する全スナップショットで、x⁺ = 3e-3, 1e-2, 3e-2, 0.1 の Nu・分子・分母と総入熱を系列 CSV に書く
             (check_quasisteady.py --series-csv に渡す)。対照 run 単独の q0・(Tw0−Tb0) の系列も別 CSV に書く

    python3 case/63.graetz_cht/eval_graetz.py snap   <heated_run> <control_run> [--step N] [--win-lo 3e-3 --win-hi 0.1]
    python3 case/63.graetz_cht/eval_graetz.py series <heated_run> <control_run>
"""
from __future__ import annotations

import argparse
import functools
import re
import sys
from pathlib import Path

import h5py
import numpy as np
from scipy.integrate import simpson

_trapz = getattr(np, "trapezoid", None) or np.trapz   # numpy 2 は trapz を削除 (AWS で落ちた 2026-09-28)

import graetz_common as gc
import graetz_ref

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "solver_density_cuda" / "tools"))
from check_quasisteady import _monotone_limit  # noqa: E402  (漸近値の推定を判定ツールと揃える)

PTS = (3e-3, 1e-2, 3e-2, 0.1)


# ------------------------------------------------------------------ 基準解
@functools.lru_cache(maxsize=1)
def reference():
    xs, nu, tb = graetz_ref.solve_march(400, 0.2, 3000, return_tb=True)
    return xs, nu, tb


def ref_nu_segment(xlo, xhi, xnode):
    """区間 [xlo, xhi] (x⁺、加熱区間の外は q=0) の q ∝ Nu θ_b の平均を、節点 x⁺ の θ_b で割った Nu。"""
    xs, nu, tb = reference()
    q = nu * tb                                     # q D/(k ΔT_in)
    lo, hi = max(xlo, 0.0), xhi
    if hi <= lo:
        return np.nan
    # 対数刻みの点列で台形積分 (x⁺ → 0 の q ∝ x^{-1/3} は可積分。区間に 0 を含むなら 1e-9 から)
    g = np.geomspace(max(lo, 1e-9), hi, 400)
    qg = np.interp(np.log(g), np.log(xs), q)
    qbar = _trapz(qg, g) / (xhi - xlo)            # 区間全体 (加熱外は 0) の平均
    tbn = np.interp(np.log(max(xnode, 1e-9)), np.log(xs), tb)
    return qbar / tbn


# ------------------------------------------------------------------ run の読み込み
def refuse(msg):
    print(f"REFUSED: {msg}")
    print("VERDICT: REFUSED")
    sys.exit(2)


class Run:
    def __init__(self, path: Path):
        self.p = Path(path)
        with h5py.File(self.p / "mesh.h5", "r") as m:
            self.xyz = np.asarray(m["MESH/COORD"][:], float).reshape(-1, 3)
        x = np.round(self.xyz[:, 0], 10)
        self.xcols = np.unique(x)
        self.col = {xv: np.where(x == xv)[0][np.argsort(self.xyz[x == xv, 1])] for xv in self.xcols}

    def steps(self):
        out = []
        for f in self.p.glob("res_[0-9]*.h5"):
            out.append(int(re.search(r"res_(\d+)\.h5$", f.name).group(1)))
        return sorted(s for s in out if (self.p / f"res_wall_heat_4_{s}.h5").exists())

    def field(self, st):
        with h5py.File(self.p / f"res_{st}.h5", "r") as h:
            V = h["VALUE"]
            d = {k: np.asarray(V[k][:], float) for k in ("ro", "Ux", "Uy", "T", "P")}
        for k, v in d.items():
            if not np.isfinite(v).all():
                raise SystemExit(f"{self.p}/res_{st}.h5: {k} に非有限値")
        return d

    def expected_wall_x(self):
        """メッシュから決まる加熱壁の節点の x (r = R、0 ≤ x ≤ L_heat)。壁ダンプはこれと 1 対 1 でなければならない。"""
        rmax = self.xyz[:, 1].max()
        on = np.abs(self.xyz[:, 1] - rmax) < 1e-12
        x = self.xyz[on, 0]
        return np.sort(x[(x >= -1e-7) & (x <= gc.L_HEAT + 1e-7)])      # 座標は float32 保存 (0.1728 → +4.5e-9)

    def wall(self, st):
        with h5py.File(self.p / f"res_wall_heat_4_{st}.h5", "r") as w:
            c = np.asarray(w["MESH/COORD"][:], float).reshape(-1, 3)
            q = np.asarray(w["VALUE/iface_q_eff"][:], float)
            tw = np.asarray(w["VALUE/iface_Tw_bc"][:], float)
            ok = np.asarray(w["VALUE/iface_ok"][:], float)
        o = np.argsort(c[:, 0])
        xw = c[o, 0]
        # 期待節点集合との照合 (2026-09-27 codex plan レビュー M3: 欠損・重複を NaN 検査では拾えない)
        exp = self.expected_wall_x()
        if len(xw) != len(np.unique(np.round(xw, 12))):
            refuse(f"{self.p} step {st}: 壁ダンプの節点に重複がある")
        if len(xw) != len(exp) or not np.allclose(xw, exp, rtol=0, atol=1e-9):
            refuse(f"{self.p} step {st}: 壁ダンプの節点 ({len(xw)}) がメッシュの加熱壁節点 ({len(exp)}) と一致しない (欠損?)")
        return xw, q[o], tw[o], ok[o]


def rint(y, r):
    """半径方向の積分。格子は一様・区間数が偶数なので Simpson 則 (3 次式まで厳密: 放物速度 × r の ṁ が誤差なく出る)。
    台形則だと N_r=16 で ṁ が (h/R)² = 0.39 % ずれ、x⁺ を通じて入口近くの Nu を動かす (2026-09-27)。"""
    h = np.diff(r)
    if len(r) % 2 == 0 or np.ptp(h) > 1e-4 * h.mean():          # 座標は float32 保存 (幅の相対ばらつき ~2e-6)
        raise SystemExit(f"半径方向の節点が一様・偶数区間でない (点 {len(r)}、幅 {h.min():.3e}..{h.max():.3e})")
    return simpson(y, x=r)


def col_integrals(run: Run, f, xv):
    """節点列 x = xv の ∫ρu r dr, ∫ρuT r dr, ∫u r dr (Simpson、r は 0 → R)。"""
    i = run.col[np.round(xv, 10)]
    r = run.xyz[i, 1]
    ru = f["ro"][i] * f["Ux"][i]
    m = rint(ru * r, r)
    mt = rint(ru * f["T"][i] * r, r)
    uu = rint(f["Ux"][i] * r, r)
    return m, mt, uu, r, i


def nu_profile(hr: Run, cr: Run, st: int):
    """加熱 run と対照 run の同一 step から、壁節点ごとの x, x⁺, Nu, 分子, 分母, ok を返す。"""
    fh, fc = hr.field(st), cr.field(st)
    xw, qh, twh, okh = hr.wall(st)
    xw0, qc, twc, okc = cr.wall(st)
    if not np.array_equal(xw, xw0):
        raise SystemExit("加熱 run と対照 run の壁節点が一致しない (同じ格子か?)")
    tbh, tbc, mdot = [], [], []
    for xv in xw:
        m, mt, _, _, _ = col_integrals(hr, fh, xv)
        m0, mt0, _, _, _ = col_integrals(cr, fc, xv)
        tbh.append(mt / m); tbc.append(mt0 / m0); mdot.append(m)
    tbh, tbc, mdot = map(np.asarray, (tbh, tbc, mdot))
    # 質量流束 (per rad) m' = ∫ρu r dr → 断面平均 ρU = 2 m'/R² → Re = ρU D/μ
    Re = (2.0 * mdot.mean() / gc.R ** 2) * gc.D / gc.MU
    Pe = Re * gc.PR
    num = -(qh - qc)
    den = (twh - tbh) - (twc - tbc)
    nu = num * gc.D / (gc.K_F * den)
    xp = xw / (gc.D * Pe)
    return dict(x=xw, xp=xp, nu=nu, num=num, den=den, ok=(okh == 1) & (okc == 1), Pe=Pe, Re=Re,
                mdot=mdot, q=qh, q0=qc, tw=twh, tb=tbh, tw0=twc, tb0=tbc)


def ref_for_nodes(xp):
    """各壁節点の双対面区間 (隣接節点との中点) で平均した基準 Nu。端の節点は片側だけ。"""
    mid = 0.5 * (xp[1:] + xp[:-1])
    lo = np.concatenate([[xp[0] - (mid[0] - xp[0])], mid])     # x=0 の節点は上流側の半区間を含む (加熱外 q=0)
    hi = np.concatenate([mid, [xp[-1]]])
    return np.array([ref_nu_segment(a, b, c) for a, b, c in zip(lo, hi, xp)])


# ------------------------------------------------------------------ V-g1 (流れ・入口・軸)
def flow_checks(run: Run, st: int, rows: list):
    f = run.field(st)
    # 加熱区間中央の列
    xs = run.xcols
    xm = xs[np.argmin(np.abs(xs - gc.L_HEAT / 2))]
    m, mt, uu, r, i = col_integrals(run, f, xm)
    ubar = uu / (gc.R ** 2 / 2)
    upar = 2 * ubar * (1 - (r / gc.R) ** 2)
    dev = np.abs(f["Ux"][i] - upar).max() / (2 * ubar)
    rows.append(("V-g1 u(r) と放物線 (加熱区間中央, 同流量) max/U_max", dev, 1e-2))
    # 質量流束の変動 (全列)
    ms = np.array([col_integrals(run, f, xv)[0] for xv in xs])
    rows.append(("V-g1 ṁ(x) の変動 (max−min)/mean", float(np.ptp(ms) / ms.mean()), 1e-3))
    # 圧力勾配 (加熱区間の中央 1/2、軸上の節点)
    ax = np.where(np.abs(run.xyz[:, 1]) < 1e-12)[0]
    sel = ax[(run.xyz[ax, 0] > gc.L_HEAT / 4) & (run.xyz[ax, 0] < 3 * gc.L_HEAT / 4)]
    slope = np.polyfit(run.xyz[sel, 0], f["P"][sel], 1)[0]
    rows.append(("V-g1 −dp/dx と 8μŪ/R² の差 (相対)", float(abs(-slope / (8 * gc.MU * ubar / gc.R ** 2) - 1)), 2e-2))
    # 入口の速度分布 (入口面の節点列、指定の放物線)
    xin = xs[0]
    _, _, _, rin, iin = col_integrals(run, f, xin)
    upin = 2 * gc.U_M * (1 - (rin / gc.R) ** 2)
    rows.append(("V-g1 入口 u(r) と指定放物線 max/(2U_m)", float(np.abs(f["Ux"][iin] - upin).max() / (2 * gc.U_M)), 1e-2))
    rows.append(("V-g1 入口面 x=−L_up の T(r) の max−min [K]", float(np.ptp(f["T"][iin])), 0.05))
    # 軸近傍の u_r (r < R/4)
    near = np.where(run.xyz[:, 1] < gc.R / 4)[0]
    rows.append(("V-g1 参考: 軸近傍 (r<R/4) max|u_r| / U_m", float(np.abs(f["Uy"][near]).max() / gc.U_M), None))
    return f


def inlet_T_check(cr: Run, st: int, rows: list):
    """対照 run (ΔT=0) の加熱開始断面 x=0 の内部節点の T の max − min (入口のエントロピー固定が作る非一様)。"""
    f = cr.field(st)
    x0 = cr.xcols[np.argmin(np.abs(cr.xcols))]
    i = cr.col[x0][:-1]                               # 壁節点を除く
    rows.append(("V-g1 対照 run の x=0 断面 T の max−min [K]", float(np.ptp(f["T"][i])), 0.05))


# ------------------------------------------------------------------ main
def gate_print(rows):
    bad = []
    for name, v, tol in rows:
        if tol is None:
            print(f"  ----  {name:<58} {v:.4e}")
            continue
        ok = np.isfinite(v) and v <= tol
        bad += [] if ok else [name]
        print(f"  {'PASS' if ok else 'FAIL'}  {name:<58} {v:.4e}  (許容 {tol:g})")
    return bad


def cmd_snap(a):
    hr, cr = Run(a.heated), Run(a.control)
    common = sorted(set(hr.steps()) & set(cr.steps()))
    if not common:
        raise SystemExit("両 run に共通するスナップショットが無い")
    st = a.step if a.step is not None else common[-1]
    if st not in common:
        raise SystemExit(f"step {st} が両 run に無い (共通: {common[-5:]})")
    d = nu_profile(hr, cr, st)
    ref = ref_for_nodes(d["xp"])
    err = d["nu"] / ref - 1.0
    W = (d["xp"] >= a.win_lo) & (d["xp"] <= a.win_hi)
    W2 = (d["xp"] >= 0.03) & (d["xp"] <= a.win_hi)
    out = Path(a.heated) / f"graetz_nu_{st}.csv"
    np.savetxt(out, np.c_[d["x"], d["xp"], d["nu"], ref, err, d["num"], d["den"], d["q"], d["q0"], d["tw"], d["tb"],
                          d["tw0"], d["tb0"], d["ok"]],
               delimiter=",", comments="", fmt="%.10e",
               header="x,xplus,Nu,Nu_ref,err,num,den,q,q0,Tw,Tb,Tw0,Tb0,ok")
    print(f"=== Graetz 評価: heated {a.heated} / control {a.control} / step {st}")
    print(f"  実測 Re {d['Re']:.4f} (設定 {gc.RE:g})、Pe {d['Pe']:.3f}、窓内 {W.sum()} 節点、CSV {out}")
    rows = []
    if (~d["ok"][W]).any():
        refuse(f"窓内に iface_ok = 0 の節点が {int((~d['ok'][W]).sum())} 個 (Tw が初期値のまま: conjugateWall.cpp:273)")
    # np.max (nanmax でない): 窓内に 1 つでも非有限があれば max が nan になり FAIL する
    rows.append(("V-g2 窓内の非有限 Nu の節点数", float((~np.isfinite(err[W])).sum()), 0.0))
    rows.append(("V-g2 max_W |Nu/Nu_Graetz − 1|", float(np.max(np.abs(err[W]))) if W.any() else np.nan, 2e-2))
    rows.append(("V-g2 max_{x⁺∈[0.03,0.1]} |Nu/Nu_Graetz − 1|", float(np.max(np.abs(err[W2]))) if W2.any() else np.nan, 1e-2))
    for p in PTS:
        e = np.interp(np.log(p), np.log(d["xp"][d["xp"] > 0]), err[d["xp"] > 0])
        rows.append((f"参考: x⁺={p:g} の誤差 (対数補間)", float(e), None))
    flow_checks(hr, st, rows)
    inlet_T_check(cr, st, rows)
    # 参考 (合否に使わない、plan §6 V-g1): 加熱 run の x=0 断面 T(r) と ellip (Pe 720) の場 T_w − ΔT θ(ξ) の差
    refp = Path(__file__).resolve().parent / "ref_x0_profile.csv"
    if refp.exists():
        rp = np.genfromtxt(refp, delimiter=",", names=True)
        fh = hr.field(st)
        x0 = hr.xcols[np.argmin(np.abs(hr.xcols))]
        i = hr.col[x0]
        xi = hr.xyz[i, 1] / hr.xyz[i, 1].max()
        dT = d["tw"][0] - gc.T_IN
        Tref = d["tw"][0] - dT * np.interp(xi, rp["xi"], rp["theta"])
        rows.append(("参考: 加熱 x=0 断面 T(r) と ellip の差 max [K]", float(np.abs(fh["T"][i] - Tref).max()), None))
        rows.append(("V-g1 参考: 入口面 T の 300 K からの平均偏差 [K]", float(fh["T"][hr.col[hr.xcols[0]]].mean() - gc.T_IN), None))
    bad = gate_print(rows)
    print(f"\nVERDICT (閾値のみ; 収束・準定常・G-if は別ツール): {'PASS' if not bad else 'FAIL'}"
          + (f"  ({'; '.join(bad)})" if bad else ""))
    return 0 if not bad else 1


def cmd_series(a):
    hr, cr = Run(a.heated), Run(a.control)
    common = sorted(set(hr.steps()) & set(cr.steps()))
    if len(common) < 4:
        refuse(f"共通スナップショットが {len(common)} 個しかない")
    rows, rows0, pern = [], [], []
    Wfix = None
    for st in common:
        d = nu_profile(hr, cr, st)
        xp = d["xp"]
        W = (xp >= a.win_lo) & (xp <= a.win_hi)
        if Wfix is None:
            Wfix = W.copy()                               # 窓内の節点集合は最初のスナップショットで固定 (以後同じ節点を追う)
        if (~d["ok"][Wfix]).any():
            refuse(f"step {st}: 窓内に iface_ok = 0 の節点")
        m = xp > 0
        vals = []
        for p in PTS:
            for key in ("nu", "num", "den"):
                vals.append(float(np.interp(np.log(p), np.log(xp[m]), d[key][m])))
        Q = float(_trapz(d["num"], d["x"]) * gc.R)
        rows.append([st] + vals + [Q])
        pern.append(np.concatenate([[st], d["nu"][Wfix], d["num"][Wfix], d["den"][Wfix]]))
        v0 = []
        for p in PTS:
            v0.append(float(np.interp(np.log(p), np.log(xp[m]), d["q0"][m])))
            v0.append(float(np.interp(np.log(p), np.log(xp[m]), (d["tw0"] - d["tb0"])[m])))
        rows0.append([st] + v0)
    cols = ["step"] + [f"{k}_{p:g}" for p in PTS for k in ("nu", "num", "den")] + ["Qnet"]
    out = Path(a.heated) / "graetz_series.csv"
    np.savetxt(out, np.array(rows), delimiter=",", comments="", header=",".join(cols),
               fmt=["%d"] + ["%.12e"] * (len(cols) - 1))
    cols0 = ["step"] + [f"{k}_{p:g}" for p in PTS for k in ("q0", "dTwb0")]
    out0 = Path(a.control) / "graetz_control_series.csv"
    np.savetxt(out0, np.array(rows0), delimiter=",", comments="", header=",".join(cols0),
               fmt=["%d"] + ["%.12e"] * (len(cols0) - 1))
    # 窓内の全節点 (codex plan レビュー M3: 4 観測点を避けた局所ドリフト・節点間振動を拾う)
    nW = int(Wfix.sum())
    colsN = ["step"] + [f"nu_n{i}" for i in range(nW)] + [f"num_n{i}" for i in range(nW)] + [f"den_n{i}" for i in range(nW)]
    P = np.array(pern)
    outN = Path(a.heated) / "graetz_series_nodes.csv"
    np.savetxt(outN, P, delimiter=",", comments="", header=",".join(colsN), fmt=["%d"] + ["%.12e"] * (len(colsN) - 1))
    # 末尾半分 (スナップショット数) の各節点の Nu の (max−min)/|mean|
    tail = P[len(P) // 2:, 1:1 + nW]
    var = np.ptp(tail, axis=0) / np.abs(tail.mean(axis=0))
    worst = int(np.argmax(var))
    tool = "python3 solver_density_cuda/tools/check_quasisteady.py"
    print(f"{len(common)} スナップショット (step {common[0]}–{common[-1]})、窓内 {nW} 節点")
    print(f"  → {out}\n  → {out0}\n  → {outN}")
    print(f"準定常 (加熱・4 点): {tool} --series-csv {out} --series-cols {','.join(cols[1:])} --tail 0.5 --drift 0.001 --osc 0.001")
    print(f"準定常 (加熱・全節点): {tool} --series-csv {outN} --series-cols {','.join(colsN[1:])} --tail 0.5 --drift 0.001 --osc 0.001")
    print(f"準定常 (対照): {tool} --series-csv {out0} --series-cols {','.join(cols0[1:])} --tail 0.5 --drift 0.001 --osc 0.001")
    # 反復の不確かさ u_it (plan §4.6): 単調な節点は check_quasisteady と同じ漸近値推定の |最終/漸近 − 1|、
    # 単調でない節点は末尾変動。末尾変動 (ゲート) とは別の量 (単調な残存過渡は末尾変動が小さくても漸近値から離れうる)
    steps_arr = P[:, 0]
    uit = np.empty(nW)
    kinds = []
    for j in range(nW):
        mono = _monotone_limit(steps_arr, P[:, 1 + j])
        if mono is not None:
            uit[j] = abs(mono[2]); kinds.append("単調")
        else:
            uit[j] = var[j]; kinds.append("非単調")
    ju = int(np.argmax(uit))
    print(f"  反復の不確かさ u_it = max_node {uit.max():.3e} (窓内 #{ju}, {kinds[ju]}; 単調 {kinds.count('単調')} / {nW} 節点)"
          f" — 比較の前提: V-g4 ≤ 1e-3、固体層 ≤ 1.5e-4、V-g3 は格子間差がこれ以下なら保留")
    with open(Path(a.heated) / "graetz_uit.txt", "w") as fo:
        fo.write(f"u_it {uit.max():.6e}\ntail_var_max {var.max():.6e}\nsnapshots {len(common)} steps {common[0]}-{common[-1]}\n")
    ok = bool(np.isfinite(var).all() and var.max() <= a.tail_var)
    print(f"  {'PASS' if ok else 'FAIL'}  末尾半分 ({len(tail)} スナップショット) の各節点 Nu の (max−min)/|mean| の最大 "
          f"{var.max():.3e} (節点 x⁺ 窓内 #{worst})  (許容 {a.tail_var:g})")
    print(f"VERDICT (末尾変動のみ; check_quasisteady は上のコマンドで別に回す): {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mode", choices=["snap", "series"])
    ap.add_argument("heated")
    ap.add_argument("control")
    ap.add_argument("--step", type=int, default=None)
    ap.add_argument("--win-lo", type=float, default=3e-3)
    ap.add_argument("--win-hi", type=float, default=0.1)
    ap.add_argument("--tail-var", type=float, default=3e-4, help="series: 末尾半分の各節点 Nu の変動の許容 (登録 0.03 %%)")
    a = ap.parse_args()
    return cmd_snap(a) if a.mode == "snap" else cmd_series(a)


if __name__ == "__main__":
    sys.exit(main())
