#!/usr/bin/env python3
"""case/56 横すきま評価 — 前向き壁の $q_w$ を深さ別に出し、TP-1187 の実測と比べる。

この 2D 計算は縦すきまの衝突を含まない。実測は衝突域で測られているので値は大きく違うが、
**2D と 3D に単調な大小関係の保証は無く、差を衝突の寄与として分離する根拠も無い**
(2026-09-20 codex result 1 巡目 M8 で撤回)。**探索的な参考値**として扱う。

$q_w=\\lambda_w(dT/dn)_w$。壁ノードから内側 2 点の 2 次片側差分 (case/50/55 と同式)。
下流側 (前向き壁 = 計測面) の輪郭: 円弧 arc_d → 鉛直壁 → 床。
深さ $z$ は上面 (y=0) からの深さで、TP-1187 Fig 5 の熱電対深さと同じ定義。
"""
import argparse, json, sys
from pathlib import Path
import numpy as np, h5py

CASE = Path(__file__).resolve().parents[1]
TC_DEPTH_CM = {"92": 0.25, "91": 0.51, "90": 0.76, "89": 1.52, "88": 2.54, "87": 3.81}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--res", default=None)
    ap.add_argument("--qfp", type=float, default=None, help="既定は conditions の run 位置 II")
    a = ap.parse_args()
    rd = CASE / a.run
    if list(rd.glob("res_nan_*.h5")):
        raise SystemExit(f"REFUSED: {a.run} は発散している")
    files = sorted(rd.glob("res_[0-9]*.h5"), key=lambda f: int(f.stem.split("_")[1]))
    res = rd / a.res if a.res else files[-1]
    setup = json.loads((rd / "case_setup.json").read_text(encoding="utf-8"))
    cond = json.loads((CASE / "conditions.json").read_text(encoding="utf-8"))
    geom = json.loads((CASE / "geometry.json").read_text(encoding="utf-8"))
    ser = [s for s in cond["series"] if s["run"] == setup["tp1187_run"]][0]
    qfp = (a.qfp if a.qfp else ser["qfp_loc_II_kW_m2"]) * 1e3
    W = 0.18e-2; r = geom["tile"]["edge_radius_cm"] * 1e-2
    xd = 0.5 * W

    # forge の壁出力 (res_gap_6_*.h5) の qwall を使う。
    # 場から幾何的に差分を取る方式は、襟が床付近で傾いているため同一 y の内側ノードが
    # 無く、壁点 2 個しか拾えなかった (2026-09-20 に踏んだ)。
    step = int(res.stem.split("_")[1])
    wf = rd / f"res_gap_6_{step}.h5"
    if not wf.exists():
        raise SystemExit(f"REFUSED: {wf.name} が無い (bcond の outputHDFflg を確認)")
    with h5py.File(wf, "r") as f:
        wc = f["/MESH/COORD"][:].reshape(-1, 3)
        qw = -f["/VALUE/qwall"][:].astype(float)     # 流体→壁を正にする
        tw = f["/VALUE/Ts"][:].astype(float)
    # 下流側 (前向き壁 = 計測面) は x > 0。深さ = -y
    m = wc[:, 0] > 1e-9
    dep, q = -wc[m, 1], qw[m]
    o = np.argsort(dep)
    A = np.column_stack([dep[o], q[o], tw[m][o]])
    np.savetxt(rd / "gap_wall_q.csv", A, delimiter=",",
               header="depth_m,q_w_W_m2,Tw_K", comments="")

    print(f"run {a.run}  res {res.name}  TP-1187 run {setup['tp1187_run']} 位置 II")
    print(f"q_FP = {qfp*1e-3:.2f} kW/m²   W = {W*1e3:.2f} mm   壁点 {len(A)}")
    print(f"\n{'TC':>4} {'深さ[cm]':>8} {'z/W':>6} {'q forge':>9} {'q/q_FP':>8} "
          f"{'実測 run8':>9} {'forge/実測':>10}")
    meas = {"87": 0.00, "88": 0.05, "90": 2.44, "91": 3.21, "92": 4.06}
    for tc, dcm in sorted(TC_DEPTH_CM.items(), key=lambda kv: kv[1]):
        dep = dcm * 1e-2
        if dep < A[0, 0] or dep > A[-1, 0]:
            print(f"{tc:>4} {dcm:8.2f} {dep/W:6.2f}   (壁データ範囲外)")
            continue
        q = float(np.interp(dep, A[:, 0], A[:, 1]))
        mv = meas.get(tc)
        rr = f"{q/qfp/mv:10.3f}" if mv else f"{'-':>10}"
        print(f"{tc:>4} {dcm:8.2f} {dep/W:6.2f} {q*1e-3:9.3f} {q/qfp:8.4f} "
              f"{(f'{mv:9.2f}' if mv is not None else '        -')} {rr}")
    print("\n注: 2D と 3D に単調な大小関係の保証はない (縦すきまを除くと流入・圧力・再循環・"
          "熱輸送がすべて変わる)。**この比は片側検証ではなく探索的な参考値**で、"
          "差を衝突加熱の寄与として分離する根拠も無い (2026-09-20 codex result 1 巡目 M8)。")
    print(f"  → {rd/'gap_wall_q.csv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
