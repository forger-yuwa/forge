#!/usr/bin/env python3
"""積分量での比較 — **観測モデルに鈍い比較量**を作る。

局所プロファイルは薄板内の横方向伝導に汚染され、リップ端の扱い (原典から決まらない) が
効果より大きく効くので判定に使えない (#30 の結論)。一方**横方向伝導は熱を再配分するが
保存する**ので、深さ方向の積分 $\\int (q/q_{fp})\\,d(x/d)$ は局所値より鈍いはず。

そこで:
  (a) 実測 (Fig 6 digitize) を積分
  (b) CFD 生を同じ区間で積分
  (c) CFD に観測モデルを掛けたものを同じ区間で積分 (リップ端 断熱 / Dirichlet の両極)
を並べ、**(c) の幅が (a)-(b) 差より小さいか**を見る。小さければ積分量は使える。
"""
import argparse, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
CASE = HERE.parent
sys.path.insert(0, str(HERE))
import importlib.util as _ilu
_sp = _ilu.spec_from_file_location("skin_model", HERE / "skin_model.py")
sm = _ilu.module_from_spec(_sp); _sp.loader.exec_module(sm)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--t-eval", type=float, default=0.8)
    ap.add_argument("--n", type=int, default=1201)
    a = ap.parse_args()
    geom = json.loads((CASE / "geometry.json").read_text(encoding="utf-8"))
    d = geom["cavity"]["depth"] * 1e-3
    cases = [("run_0006_T1_wd0063_long", 0.063, "ref/w70_fig6a_wd0063_rear.csv"),
             ("run_0008_T1_wd0211_long", 0.211, "ref/w70_fig6b_wd0211_rear.csv")]
    print(f"深さ積分 I = ∫ (q/q_fp) d(x/d)。積分区間は**実測が存在する範囲**に揃える。\n")
    for run, wd, ref in cases:
        setup = json.loads((CASE / run / "case_setup.json").read_text(encoding="utf-8"))
        W = geom["cavity"]["widths"][str(wd)] * 1e-3
        qfp = setup["series"]["qfp_kW"] * 1e3
        R = sm.read_ref(CASE / ref)
        o = np.argsort(R[:, 0]); R = R[o]
        lo, hi = R[0, 0], R[-1, 0]
        s, qc = sm.build_path(run, d, W, a.n)
        xd = s / d
        m = (xd >= lo) & (xd <= hi)
        T_lip = qfp * a.t_eval / (sm.RHO_C * sm.TAU)
        app_d = sm.solve(s, qc, a.t_eval, T_lip, lip="dirichlet")[1]
        app_a = sm.solve(s, qc, a.t_eval, 0.0, lip="adiabatic")[1]
        I_meas = np.trapz(np.interp(xd[m], R[:, 0], R[:, 1]), xd[m])
        I_raw = np.trapz(qc[m] / qfp, xd[m])
        I_ad = np.trapz(app_a[m] / qfp, xd[m])
        I_di = np.trapz(app_d[m] / qfp, xd[m])
        obs = max(I_ad, I_di) - min(I_ad, I_di)
        gap = abs(I_meas - I_raw)
        print(f"=== {run}  w/d={wd}  積分区間 x/d = {lo:.3f}–{hi:.3f} ===")
        print(f"  実測                    I = {I_meas:.5f}")
        print(f"  CFD 生                  I = {I_raw:.5f}   (実測比 {I_raw/I_meas:.3f})")
        print(f"  CFD + 観測モデル 断熱   I = {I_ad:.5f}   (実測比 {I_ad/I_meas:.3f})")
        print(f"  CFD + 観測モデル Dir    I = {I_di:.5f}   (実測比 {I_di/I_meas:.3f})")
        print(f"  → 観測モデルによる幅 {obs:.5f} vs CFD-実測差 {gap:.5f}  "
              f"= {obs/gap if gap>0 else float('inf'):.2f} 倍")
        verdict = "使える (モデル幅 < 差)" if obs < gap else "**使えない (モデル幅 >= 差)**"
        print(f"  積分量は比較量として: {verdict}\n")
    print("局所プロファイルとの対比: 局所値ではリップ端の扱いだけで z/W=2 が 0.112 ↔ 0.316 (3 倍) 動く。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
