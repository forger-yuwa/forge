#!/usr/bin/env python3
r"""case/65 B1 (壁温の更新の有無 A/B) の観測系列と予備確認 (plan `boundary-cht-conjugate-flat-plate.md` §4.1.1・§4.1.2)。

    python3 ab_series.py pre  <run_A> <run_B>           # 予備確認: 受入条件と採取間隔の検定
    python3 ab_series.py series <run> [--every N]       # 系列 CSV (節点群ごと) を書く

節点集合は mesh.h5 の `BCONDS/<physID>/iCells` と `MESH/COORD` から作る (`CELLS/centCoords` は使わない)。
- 板の界面 (physID 5): 前縁帯 x/L∈[0,0.02]・後縁帯 [0.98,1] (各 10 点を期待)、全界面の熱流束も保持
- slip: 上流 (physID 4) x/L∈[−0.1,0)・下流 (physID 6) (1,1.1] (各 21 点を期待)
- 上記の直上の内部 1 列 (同じ x で y が最小の正の節点)
尺度: P は q∞ = ½ρ∞U∞² = 709.275 Pa、T は 10 K、Uy は U∞、界面熱流束は元 run (run_0007) の step 600000 の
評価窓 x/L∈[0.2,0.9] の台形平均の絶対値。欠落・重複・非有限・期待数の不一致は REFUSED。
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import plate_common as pc  # noqa: E402

trap = getattr(np, "trapezoid", None) or np.trapz
L = pc.L
Q_INF = 709.275
T_SCALE = 10.0
EXPECT = {"le": 10, "te": 10, "up": 21, "down": 21}


def refuse(msg):
    print(f"VERDICT: REFUSED ({msg})")
    sys.exit(2)


def steps_of(run, pat):
    return sorted(int(m.group(1)) for p in Path(run).glob("*.h5") for m in [re.fullmatch(pat, p.name)] if m)


def node_sets(run):
    with h5py.File(Path(run) / "mesh.h5", "r") as m:
        ic = {p: np.asarray(m[f"BCONDS/{p}/iCells"][:], int) for p in (4, 5, 6)}
    with h5py.File(Path(run) / "res_0.h5", "r") as r:
        c = np.asarray(r["MESH/COORD"][:], float).reshape(-1, 3)
    x, y = c[:, 0] / L, c[:, 1]
    tol = 1e-9
    def band(ids, lo, hi, lo_open=False, hi_open=False):
        xs = x[ids]
        ok = (xs > lo + tol if lo_open else xs >= lo - tol) & (xs < hi - tol if hi_open else xs <= hi + tol)
        s = ids[ok]
        return s[np.argsort(x[s])]
    sets = {"le": band(ic[5], 0.0, 0.02), "te": band(ic[5], 0.98, 1.0),
            "up": band(ic[4], -0.1, 0.0, hi_open=True), "down": band(ic[6], 1.0, 1.1, lo_open=True)}
    for k, s in sets.items():
        if len(s) != EXPECT[k]:
            refuse(f"節点群 {k}: {len(s)} 点 (期待 {EXPECT[k]})")
        if len(np.unique(s)) != len(s) or np.abs(y[s]).max() > tol:
            refuse(f"節点群 {k}: 重複または y≠0")
    # 直上の内部 1 列: 同じ x で y が最小の正の節点
    pos = np.where(y > tol)[0]
    for k in list(sets):
        above = []
        for n in sets[k]:
            cand = pos[np.abs(x[pos] - x[n]) < 1e-7]
            if len(cand) == 0:
                refuse(f"節点群 {k}: x/L={x[n]:.5f} の直上節点が無い")
            above.append(cand[np.argmin(y[cand])])
        sets[k + "_in"] = np.array(above)
    return sets, c, ic[5]


def fluid_series(run, sets, steps):
    out = {k: {q: [] for q in ("P", "T", "Uy")} for k in sets}
    for st in steps:
        with h5py.File(Path(run) / f"res_{st}.h5", "r") as r:
            v = {q: np.asarray(r[f"VALUE/{q}"][:], float) for q in ("P", "T", "Uy")}
        for k, s in sets.items():
            for q in v:
                out[k][q].append(v[q][s])
    for k in out:
        for q in out[k]:
            a = np.array(out[k][q])
            if not np.isfinite(a).all():
                refuse(f"{run}: 流体系列 {k}/{q} に非有限値")
            out[k][q] = a
    return out


def wall_series(run, steps):
    xs, Ts, qs = None, [], []
    for st in steps:
        with h5py.File(Path(run) / f"res_plate_5_{st}.h5", "r") as w:
            c = np.asarray(w["MESH/COORD"][:], float).reshape(-1, 3)
            T = np.asarray(w["VALUE/iface_Tw_bc"][:], float)
            q = -np.asarray(w["VALUE/iface_q_eff"][:], float)
        o = np.argsort(c[:, 0])
        if xs is None:
            xs = c[o, 0]
        elif not np.allclose(xs, c[o, 0], rtol=0, atol=1e-12):
            refuse(f"{run} step {st}: 壁ダンプの節点が変わった")
        Ts.append(T[o]); qs.append(q[o])
    Ts, qs = np.array(Ts), np.array(qs)
    if not (np.isfinite(Ts).all() and np.isfinite(qs).all()):
        refuse(f"{run}: 壁ダンプに非有限値")
    return xs, Ts, qs


def q_scale():
    src = HERE / "run_0007_c1_n64" / "res_plate_5_600000.h5"
    if not src.exists():
        refuse(f"熱流束の尺度の元 {src} が無い")
    with h5py.File(src, "r") as w:
        c = np.asarray(w["MESH/COORD"][:], float).reshape(-1, 3)
        q = -np.asarray(w["VALUE/iface_q_eff"][:], float)
    o = np.argsort(c[:, 0]); x = c[o, 0]; q = q[o]
    win = (x >= 0.2 * L - 1e-12) & (x <= 0.9 * L + 1e-12)
    s = abs(trap(q[win], x[win]) / (x[win][-1] - x[win][0]))
    if not (np.isfinite(s) and s > 0):
        refuse(f"熱流束の尺度 {s}")
    return s


def amp(a):
    """節点ごとの時間標準偏差の最大 (a: [時刻, 節点])。"""
    return float(np.std(a, axis=0).max())


def pre(a):
    runA, runB = Path(a.runA), Path(a.runB)
    sets, c, _ = node_sets(runA)
    stF = steps_of(runA, r"res_(\d+)\.h5"); stW = steps_of(runA, r"res_plate_5_(\d+)\.h5")
    if stF != steps_of(runB, r"res_(\d+)\.h5") or stW != steps_of(runB, r"res_plate_5_(\d+)\.h5"):
        refuse("A と B の出力 step が揃っていない")
    if np.diff(stW).max() != 1:
        refuse("毎 step 出力になっていない")
    qs = q_scale()
    print(f"=== 予備確認 A={runA.name} B={runB.name}  壁ダンプ {len(stW)} 枚 (step {stW[0]}–{stW[-1]})、流体 {len(stF)} 枚")
    print(f"  節点群: " + ", ".join(f"{k} {len(v)}" for k, v in sets.items()) + f"、熱流束の尺度 {qs:.6g} W/m²")
    xA, TA, qA = wall_series(runA, stW); xB, TB, qB = wall_series(runB, stW)
    ok = True
    # 受入条件
    prof = np.loadtxt(runB / "wall_profile_5.csv", skiprows=1)
    pT = prof[np.argsort(prof[:, 0]), 3]
    if len(pT) != len(xB) or not np.allclose(np.sort(prof[:, 0]), xB, atol=1e-9):
        refuse("wall_profile_5.csv と壁ダンプの節点が一致しない")
    d_init = float(np.abs(TA[0] - TB[0]).max()); d_prof = float(np.abs(TB[0] - pT).max())
    d_fix = float(np.abs(TB - TB[0]).max()); d_A = float(np.abs(TA - TA[0]).max())
    for nm, v, lim, good in (("初期壁温の A/B 差 (最初のダンプ) [K]", d_init, 1e-9, d_init <= 1e-9),
                             ("B の初期壁温 − wall_profile [K]", d_prof, 1e-9, d_prof <= 1e-9),
                             ("B の壁温の時間変化 max [K]", d_fix, 0.0, d_fix == 0.0)):
        ok &= good
        print(f"  {'PASS' if good else 'FAIL'}  {nm}: {v:.3e} (≤ {lim:g})")
    hist = np.genfromtxt(runA / "conjugate_history.csv", delimiter=",", names=True)
    nup = int(np.atleast_1d(hist["update"]).size)
    goodA = nup >= 3 and d_A > 0
    ok &= goodA
    print(f"  {'PASS' if goodA else 'FAIL'}  A の連成継続: 更新 {nup} 回 (step {', '.join(str(int(s)) for s in np.atleast_1d(hist['step']))})、壁温の時間変化 max {d_A:.3e} K")
    # 採取間隔の検定: 毎 step と N step 間引きの標準偏差 (量ごと・節点群ごと)
    N = a.every
    print(f"  採取間隔の検定 (毎 step vs {N} step 間引き、標準偏差の節点最大、10 % 以内):")
    fA = fluid_series(runA, sets, stF); fB = fluid_series(runB, sets, stF)
    le = (xA / L <= 0.02 + 1e-9); te = (xA / L >= 0.98 - 1e-9)
    rows = []
    for tag, F, q in (("A", fA, qA), ("B", fB, qB)):
        for k in sets:
            for qn, sc in (("P", Q_INF), ("T", T_SCALE), ("Uy", pc.U_INF)):
                rows.append((tag, k, qn, F[k][qn] / sc))
        for k, m in (("iface_le", le), ("iface_te", te), ("iface_all", np.ones_like(le))):
            rows.append((tag, k, "q", q[:, m] / qs))
    for tag, k, qn, s in rows:
        full = amp(s); sub = amp(s[::N])
        rel = abs(sub - full) / full if full > 0 else 0.0
        floor = full < 1e-12
        good = floor or rel <= 0.10
        ok &= good
        print(f"    {'PASS' if good else 'FAIL'}  {tag} {k:9s} {qn:2s}: 毎 step {full:.3e} / 間引き {sub:.3e}  (差 {rel*100:.1f} %{', 測定床' if floor else ''})")
    print(f"VERDICT: {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    p = sp.add_parser("pre"); p.add_argument("runA"); p.add_argument("runB"); p.add_argument("--every", type=int, default=11)
    a = ap.parse_args()
    return pre(a)


if __name__ == "__main__":
    sys.exit(main())
