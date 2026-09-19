#!/usr/bin/env python3
"""case/50 平板評価 — 壁熱流束 q_w(x) を場から出し、解析値 (Pohlhausen + Eckert) と Table IV に照合する。

plan §4.8: **抽出器は平板の解析解で先に検証してからキャビティに適用する**。本スクリプトがその検証。
  q_w = λ_w (dT/dn)_w   — λ は forge の出力 `thermCond` (viscMethod 2 の Mason–Saxena) をそのまま使う
  dT/dn は壁ノードから内側 2 点の**2 次片側差分** (非等間隔の 2 次多項式当てはめ)

usage: python3 tools/plate_eval.py RUN [--res res_8000.h5] [--series] [--out prefix]
"""
import argparse, json, sys
from pathlib import Path
import numpy as np
import h5py

HERE = Path(__file__).resolve().parent
CASE = HERE.parent
sys.path.insert(0, str(HERE))
from gas_model import CombustionProducts                     # noqa: E402
from conditions import solve_state, q_flat_plate             # noqa: E402


def wall_columns(coord, tol=1e-9):
    """構造格子の y 列を x でまとめる → {x: [node index (y 昇順)]}"""
    x = np.round(coord[:, 0], 9)
    order = np.lexsort((coord[:, 1], x))
    cols, cur, cx = {}, [], None
    for i in order:
        if cx is None or abs(x[i] - cx) > tol:
            if cur:
                cols[cx] = np.array(cur)
            cur, cx = [], x[i]
        cur.append(i)
    if cur:
        cols[cx] = np.array(cur)
    return cols


def dTdn_wall(y, T):
    """非等間隔 3 点の 2 次当てはめで壁 (y[0]) における勾配。"""
    h1, h2 = y[1] - y[0], y[2] - y[0]
    a = (T[1] - T[0]) / h1
    b = ((T[2] - T[0]) / h2 - a) / (h2 - h1)
    return a - b * h1


def evaluate(res, setup):
    with h5py.File(res, "r") as f:
        c = f["/MESH/COORD"][:].reshape(-1, 3)
        T = f["/VALUE/T"][:].astype(float)
        lam = f["/VALUE/thermCond"][:].astype(float)
        P = f["/VALUE/P"][:].astype(float)
        U = f["/VALUE/Ux"][:].astype(float)
    cols = wall_columns(c)
    xs, qs, Tws, lws, Pws = [], [], [], [], []
    for x, idx in sorted(cols.items()):
        if x < 1e-9 or x > setup["x_plate_end"] + 1e-9:
            continue
        col = idx[c[idx, 1] >= -1e-12]                    # 平板上 (y>=0) のみ
        if len(col) < 3:
            continue
        y = c[col, 1]
        o = np.argsort(y)
        col, y = col[o], y[o]
        if y[0] > 1e-9:
            continue
        xs.append(x); qs.append(lam[col[0]] * dTdn_wall(y[:3], T[col[:3]]))
        Tws.append(T[col[0]]); lws.append(lam[col[0]]); Pws.append(P[col[0]])
    return np.array(xs), np.array(qs), np.array(Tws), np.array(lws), np.array(Pws)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--res", default=None)
    ap.add_argument("--series", action="store_true", help="全 res_* で x_ref の q_w 時系列 CSV を書く")
    ap.add_argument("--x-ref", type=float, default=None, help="照合位置 [m] (既定 = キャビティ中点)")
    a = ap.parse_args()

    rd = CASE / a.run
    setup = json.loads((rd / "case_setup.json").read_text(encoding="utf-8"))
    geom = json.loads((CASE / "geometry.json").read_text(encoding="utf-8"))
    cond = json.loads((CASE / "conditions.json").read_text(encoding="utf-8"))
    s = setup["series"]
    w = geom["cavity"]["widths"][str(s["w_over_d"])] * 1e-3
    x_ref = a.x_ref or (geom["cavity"]["x_rear_wall_from_le"] * 1e-3 - 0.5 * w)
    setup["x_plate_end"] = 0.197

    gas = CombustionProducts(s["Tt_K"], T_react=cond["gas"]["T_react"])
    T0, U0, rho0, p0, mu0 = solve_state(gas, s["M"], s["Re_m"], s["Tt_K"])
    ana = q_flat_plate(gas, T0, U0, p0, x_ref, setup["Tw"])

    files = sorted(rd.glob("res_[0-9]*.h5"), key=lambda f: int(f.stem.split("_")[1]))
    res = rd / a.res if a.res else (files[-1] if files else rd / "mesh.h5")
    x, q, Tw, lam, Pw = evaluate(res, setup)
    i = int(np.argmin(np.abs(x - x_ref)))
    print(f"run {a.run}  res {res.name}   壁 {len(x)} 点  (x {x.min()*1e3:.1f}–{x.max()*1e3:.1f} mm)")
    print(f"照合位置 x = {x[i]*1e3:.2f} mm (目標 {x_ref*1e3:.2f})   T_w = {Tw[i]:.2f} K  λ_w = {lam[i]:.5f}  P_w = {Pw[i]:.1f} Pa")
    print()
    print(f"  forge  q_w = {q[i]*1e-3:7.2f} kW/m²")
    print(f"  解析   q_w = {ana['q']*1e-3:7.2f} kW/m²  (Pohlhausen + Eckert, Pr* {ana['Pr_star']:.3f}, Taw {ana['Taw']:.1f} K)")
    print(f"  文献   q_fp= {s['qfp_kW']:7.2f} kW/m²  (Table IV)")
    print()
    print(f"  forge / 解析 = {q[i]/ana['q']:.4f}   forge / 文献 = {q[i]/(s['qfp_kW']*1e3):.4f}   解析 / 文献 = {ana['q']/(s['qfp_kW']*1e3):.4f}")
    np.savetxt(rd / "wall_q.csv", np.c_[x, q, Tw, lam, Pw], delimiter=",",
               header="x_m,q_w_W_m2,Tw_K,lambda_w,P_w", comments="")
    print(f"  → {rd/'wall_q.csv'}")

    if a.series and files:
        rows = []
        for fpath in files:
            xx, qq, *_ = evaluate(fpath, setup)
            j = int(np.argmin(np.abs(xx - x_ref)))
            rows.append((int(fpath.stem.split("_")[1]), qq[j]))
        np.savetxt(rd / "q_series.csv", np.array(rows), delimiter=",",
                   header="step,q_w_ref_W_m2", comments="")
        print(f"  → {rd/'q_series.csv'}  ({len(rows)} スナップショット)")


if __name__ == "__main__":
    main()
