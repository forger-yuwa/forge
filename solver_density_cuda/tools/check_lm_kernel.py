#!/usr/bin/env python3
"""LM2009 遷移モデルのカーネル (`cuda_forge/transition_d.cu`) を numpy 参照実装 (`lm2009_reference.py`) と節点ごとに突き合わせる。

`output: {level: 2}` で回した LM run の res_*.h5 から、カーネルが書いた診断場 (lmFonset, lmFlength, lmFtheta, lmRethCorr,
lmGammaSep, lmPgamma, lmEgamma, lmPtheta) と、同じファイルの状態量 (ro, U, 速度勾配, k, omega, vis_lam, wall_dist, gammaTr, reTheta)
から参照実装で再計算した値を比べる。

注意: 診断場は**その step の更新前**の状態で計算され、出力の状態量は**更新後**なので、厳密には 1 step ぶんずれる。
収束場か、cfl を極小にした 1 step run で使うこと (plan turbulence-transition-lm2009 §5.1 #3)。

usage: check_lm_kernel.py RUN_DIR [--step N] [--tol 1e-3]
VERDICT PASS / FAIL (exit 0/1)。判定は「参照値の最大絶対値で正規化した差の 99.9 パーセンタイル」が tol 以下。
"""
import argparse, glob, os, sys
import h5py, numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lm2009_reference as lm


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("run"); ap.add_argument("--step", type=int, default=None); ap.add_argument("--tol", type=float, default=1e-3)
    a = ap.parse_args()
    files = sorted(glob.glob(os.path.join(a.run, "res_[0-9]*.h5")), key=lambda f: int(os.path.basename(f)[4:-3]))
    f = os.path.join(a.run, f"res_{a.step}.h5") if a.step is not None else files[-1]
    with h5py.File(f, "r") as h:
        V = {k: np.array(h["VALUE"][k], dtype=float) for k in h["VALUE"]}
    need = ["lmFonset", "lmPgamma", "dUxdx", "gammaTr", "reTheta", "vis_lam", "wall_dist"]
    lack = [n for n in need if n not in V]
    if lack:
        print(f"[check_lm_kernel] {f}: 必要な場が無い {lack} (output.level 2 の LM run が要る)"); print("VERDICT: FAIL"); return 1
    n = len(V["ro"])
    U = np.stack([V["Ux"], V["Uy"], V["Uz"]], -1)
    G = np.zeros((n, 3, 3))
    for i, ui in enumerate("xyz"):
        for j, xj in enumerate("xyz"):
            G[:, i, j] = V[f"dU{ui}d{xj}"]
    # カーネルがソースを評価しない節点 (壁・停滞点・入口ピン) は診断が 0。参照側も同じ節点を外す。
    active = V["lmFlength"] > 0.0
    r = lm.sources(V["ro"], U, G, np.maximum(V["k"], 0.0), np.maximum(V["omega"], 1e-12), V["vis_lam"], np.maximum(V["wall_dist"], 1e-30), V["gammaTr"], V["reTheta"])
    pg = r["f_length"] * lm.C_A1 * V["ro"] * np.sqrt(2.0 * np.sum((0.5 * (G + np.swapaxes(G, -1, -2))) ** 2, axis=(-1, -2))) * np.sqrt(np.maximum(r["f_onset"] * V["gammaTr"], 0)) * (1 - lm.C_E1 * V["gammaTr"])
    pairs = [("F_onset", V["lmFonset"], r["f_onset"]), ("F_length", V["lmFlength"], r["f_length"]), ("F_theta", V["lmFtheta"], r["f_theta"]),
             ("Re_theta_corr", V["lmRethCorr"], r["re_theta_corr"]), ("gamma_sep", V["lmGammaSep"], r["gamma_sep"]),
             ("P_gamma", V["lmPgamma"], pg), ("P_gamma-E_gamma", V["lmPgamma"] - V["lmEgamma"], r["src_gamma"]), ("P_theta", V["lmPtheta"], r["src_reth"])]
    print(f"[check_lm_kernel] {f}: {int(active.sum())}/{n} 節点でソースを評価")
    print(f"{'量':<18}{'max|参照|':>14}{'p99.9 差/max':>16}{'最大差/max':>14}  判定")
    ok = True
    for name, ker, ref in pairs:
        k_, r_ = ker[active], ref[active]
        sc = max(np.abs(r_).max(), 1e-30); d = np.abs(k_ - r_) / sc
        p = float(np.percentile(d, 99.9)); good = p <= a.tol; ok &= good
        print(f"{name:<18}{sc:14.5e}{p:16.3e}{d.max():14.3e}  {'ok' if good else 'NG'}")
    ge = np.maximum(V["gammaTr"], np.where(active, V["lmGammaSep"], 0.0))
    dge = np.abs(V["gammaEff"] - ge).max(); print(f"gammaEff = max(gamma, gamma_sep): 最大差 {dge:.3e}")
    if "lmCorrIter" in V:
        it = V["lmCorrIter"][active]; print(f"相関の不動点反復: 平均 {it.mean():.1f} 回, 最大 {it.max():.0f} 回, 上限 100 到達 {int((it >= 100).sum())} 節点")
    print("VERDICT:", "PASS" if ok else "FAIL"); return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
