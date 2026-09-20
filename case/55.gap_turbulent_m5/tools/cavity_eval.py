#!/usr/bin/env python3
"""case/50 キャビティ評価 — W70 (TN D-5908) の座標規約で q_s/q_fp と平均比・圧力比を出す。

W70 の規約 (図 1・図 6):
  後壁: x を**上端から下向き**に取り x/d、床: y を**後壁から前向き**に取り y/d
  平均比: q̄_c/q_fp = Q_c/Q_fp で、**分母は開口面積** w·ℓ (2D なので単位スパン当たり [W/m] を w で割る)
  圧力: キャビティ床から 0.05 in の点 ÷ 前縁から 3.5 in の模型面圧

分母は 2 通り出す (plan §4.4 / §4.2c):
  (a) 文献分母  = Table IV の q_fp
  (b) forge 分母 = 同一 x の smooth run (T0) の q_w   ← 分母バイアスを打ち消した比

usage: python3 tools/cavity_eval.py RUN --t0 run_0003_T0_A1_main [--series]
"""
import argparse, json, sys
from pathlib import Path
import numpy as np
import h5py

HERE = Path(__file__).resolve().parent
CASE = HERE.parent
sys.path.insert(0, str(HERE))


def quad_grad(s, T):
    """非等間隔 3 点 (s[0] が壁) の 2 次当てはめによる壁での dT/ds。"""
    h1, h2 = s[1] - s[0], s[2] - s[0]
    a = (T[1] - T[0]) / h1
    b = ((T[2] - T[0]) / h2 - a) / (h2 - h1)
    return a - b * h1


def load(res):
    with h5py.File(res, "r") as f:
        return (f["/MESH/COORD"][:].reshape(-1, 3), f["/VALUE/T"][:].astype(float),
                f["/VALUE/thermCond"][:].astype(float), f["/VALUE/P"][:].astype(float))


def wall_flux_line(c, T, lam, fixed_axis, fixed_val, span_axis, span_lo, span_hi, inward_sign, tol=1e-9):
    """fixed_axis = fixed_val の壁線上の各点で、内側 3 点から q_w = λ dT/dn を出す。"""
    sel = (np.abs(c[:, fixed_axis] - fixed_val) < tol) & \
          (c[:, span_axis] > span_lo - tol) & (c[:, span_axis] < span_hi + tol)
    idx = np.where(sel)[0]
    idx = idx[np.argsort(c[idx, span_axis])]
    # 内側方向の座標候補 (構造格子なので fixed_axis の一意値を使う)
    axis_vals = np.unique(np.round(c[:, fixed_axis], 9))
    if inward_sign > 0:
        cand = axis_vals[axis_vals > fixed_val + tol][:3]
    else:
        cand = axis_vals[axis_vals < fixed_val - tol][-3:][::-1]
    out_s, out_q, out_T = [], [], []
    for i in idx:
        sv = c[i, span_axis]
        pts = [i]
        for av in cand:
            m = np.where((np.abs(c[:, fixed_axis] - av) < tol) & (np.abs(c[:, span_axis] - sv) < tol))[0]
            if len(m):
                pts.append(int(m[0]))
        if len(pts) < 3:
            continue
        n = np.abs(c[pts, fixed_axis] - fixed_val)
        out_s.append(sv)
        out_q.append(lam[pts[0]] * quad_grad(n[:3], T[np.array(pts[:3])]))
        out_T.append(T[pts[0]])
    return np.array(out_s), np.array(out_q), np.array(out_T)


def evaluate(res, xf, xr, d):
    c, T, lam, P = load(res)
    # 後壁 (x=xr, y<0): 内側は -x 方向
    yb, qb, _ = wall_flux_line(c, T, lam, 0, xr, 1, -d, 0.0, inward_sign=-1)
    rear = (np.abs(yb) / d, qb)                       # x/d (上端 0 → 床 1)
    # 前壁 (x=xf, y<0): 内側は +x
    yf, qf, _ = wall_flux_line(c, T, lam, 0, xf, 1, -d, 0.0, inward_sign=+1)
    front = (np.abs(yf) / d, qf)
    # 床 (y=-d): 内側は +y。W70 は後壁から前向きに y/d
    xb, qfl, _ = wall_flux_line(c, T, lam, 1, -d, 0, xf, xr, inward_sign=+1)
    floor = ((xr - xb) / d, qfl)
    return dict(rear=rear, front=front, floor=floor, c=c, P=P)


def integrate(s, q):
    o = np.argsort(s)
    return float(np.trapz(q[o], s[o]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--t0", default="run_0003_T0_A1_main", help="分母 (b) に使う smooth run")
    ap.add_argument("--res", default=None)
    ap.add_argument("--series", action="store_true")
    ap.add_argument("--qfp-from", default=None,
                    help="分母 q_fp を取る等温平板 run (wall_q.csv を開口中点で内挿)。"
                         "h 比は駆動温度差 (T_aw - T_w) が分子分母で共通なので q 比と一致する")
    a = ap.parse_args()

    rd = CASE / a.run
    setup = json.loads((rd / "case_setup.json").read_text(encoding="utf-8"))
    geom = json.loads((CASE / "geometry.json").read_text(encoding="utf-8"))
    cond = json.loads((CASE / "conditions.json").read_text(encoding="utf-8"))
    s = setup["series"]
    w = geom["cavity"]["widths"][str(s["w_over_d"])] * 1e-3
    d = geom["cavity"]["depth"] * 1e-3
    xr = geom["cavity"]["x_rear_wall_from_le"] * 1e-3
    xf = xr - w
    x_mid = xr - 0.5 * w

    # 分母 (a) 文献 / (b) forge smooth / (参考) 解析
    if a.qfp_from:
        arr0 = np.loadtxt(CASE / a.qfp_from / "wall_q.csv", delimiter=",", skiprows=1)
        q_lit = float(np.interp(xr - 0.5 * w, arr0[:, 0], arr0[:, 1]))
        print(f"分母: 等温平板 run {a.qfp_from} の q_w(x={1e3*(xr-0.5*w):.1f} mm) = {q_lit*1e-3:.3f} kW/m²")
        LBL_A = f"forge 等温平板 ({a.qfp_from})"
    else:
        q_lit = s["qfp_kW"] * 1e3
        LBL_A = "文献"
    q_ana = float("nan")   # case/55 は乱流なので層流平板の解析値は使わない
    t0csv = CASE / a.t0 / "wall_q.csv"
    if t0csv.exists():
        arr = np.loadtxt(t0csv, delimiter=",", skiprows=1)
        q_cfd = float(np.interp(x_mid, arr[:, 0], arr[:, 1]))
    else:
        q_cfd = float("nan")

    nan = sorted(rd.glob("res_nan_*.h5"))
    if nan:
        raise SystemExit(f"REFUSED: {rd.name} は発散している ({[f.name for f in nan]})。"
                         " 評価しない (初期場を読んで『結果』にしないため)")
    files = sorted(rd.glob("res_[0-9]*.h5"), key=lambda f: int(f.stem.split("_")[1]))
    res = rd / a.res if a.res else (files[-1] if files else rd / "mesh.h5")
    ev = evaluate(res, xf, xr, d)

    print(f"run {a.run}  res {res.name}   w = {w*1e3:.3f} mm (w/d = {s['w_over_d']}), d = {d*1e3:.2f} mm")
    print(f"分母: {LBL_A} {q_lit*1e-3:.2f} / forge smooth {q_cfd*1e-3:.2f} / 解析 {q_ana*1e-3:.2f} kW/m²")
    print()
    for name, (xs, q) in (("後壁", ev["rear"]), ("前壁", ev["front"]), ("床", ev["floor"])):
        print(f"[{name}] {len(xs)} 点  q_w: {np.min(q)*1e-3:+.3f} … {np.max(q)*1e-3:+.3f} kW/m²")
        for tgt in (0.0125, 0.05, 0.1, 0.25, 0.5, 0.75, 1.0):
            if xs.min() <= tgt <= xs.max():
                qq = float(np.interp(tgt, np.sort(xs), q[np.argsort(xs)]))
                print(f"    {'x/d' if name!='床' else 'y/d'} = {tgt:5.3f} : "
                      f"q/q_lit = {qq/q_lit:7.4f}   q/q_cfd = {qq/q_cfd:7.4f}   q = {qq*1e-3:8.4f} kW/m²")
    # 総入熱 (単位スパン) と開口面積平均
    Qc = (integrate(ev["rear"][0] * d, ev["rear"][1]) + integrate(ev["front"][0] * d, ev["front"][1])
          + integrate(ev["floor"][0] * d, ev["floor"][1]))
    print()
    print(f"総入熱 Q_c' = {Qc:.3f} W/m (後壁 {integrate(ev['rear'][0]*d, ev['rear'][1]):.3f} + "
          f"前壁 {integrate(ev['front'][0]*d, ev['front'][1]):.3f} + 床 {integrate(ev['floor'][0]*d, ev['floor'][1]):.3f})")
    for lbl, qd in ((LBL_A, q_lit), ("forge smooth", q_cfd)):
        print(f"  開口面積平均 q̄_c/q_fp [{lbl}] = {Qc/w/qd:.4f}   (= Q_c/Q_fp)")
    out = dict(run=a.run, w=w, d=d, q_lit=q_lit, q_cfd=q_cfd, q_ana=q_ana, Qc_per_span=Qc,
               qbar_over_lit=Qc / w / q_lit, qbar_over_cfd=Qc / w / q_cfd)
    (rd / "cavity_eval.json").write_text(json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
    np.savetxt(rd / "cavity_rear.csv", np.c_[ev["rear"][0], ev["rear"][1]], delimiter=",",
               header="x_over_d,q_w_W_m2", comments="")
    np.savetxt(rd / "cavity_floor.csv", np.c_[ev["floor"][0], ev["floor"][1]], delimiter=",",
               header="y_over_d,q_w_W_m2", comments="")
    print(f"  → {rd/'cavity_eval.json'}, cavity_rear.csv, cavity_floor.csv")

    if a.series and files:
        rows = []
        for fp in files:
            e = evaluate(fp, xf, xr, d)
            Q = (integrate(e["rear"][0] * d, e["rear"][1]) + integrate(e["front"][0] * d, e["front"][1])
                 + integrate(e["floor"][0] * d, e["floor"][1]))
            qdeep = float(np.interp(0.75, np.sort(e["rear"][0]), e["rear"][1][np.argsort(e["rear"][0])]))
            rows.append((int(fp.stem.split("_")[1]), Q, qdeep))
        np.savetxt(rd / "cavity_series.csv", np.array(rows), delimiter=",",
                   header="step,Qc_per_span_W_m,q_rear_xd075_W_m2", comments="")
        print(f"  → {rd/'cavity_series.csv'} ({len(rows)} スナップショット)")
        # 量が落ち着いていない run の平均比を「結果」として引用させない (2026-09-19 の失敗の実体化)
        import subprocess
        tool = CASE.parents[1] / "solver_density_cuda" / "tools" / "check_quasisteady.py"
        r = subprocess.run([sys.executable, str(tool), "--series-csv", str(rd / "cavity_series.csv"),
                            "--series-cols", "Qc_per_span_W_m"], capture_output=True, text=True)
        verdict = [l for l in r.stdout.splitlines() if "Qc_per_span" in l]
        print("  " + (verdict[0].strip() if verdict else "(準定常判定を取得できず)"))
        if "STEADY" not in (verdict[0] if verdict else ""):
            print("  " + "!" * 72)
            print("  !! この run は準定常でない。上の平均比・分布を『結果』として引用しないこと。")
            print("  !! 落ち着くまで継続してから比べる (case/50 実測: 一様 IC から ~160k step、")
            print("  !!   収束場からのスキーム変更でも 20k–95k step かかった)。")
            print("  " + "!" * 72)


if __name__ == "__main__":
    main()
