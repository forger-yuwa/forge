#!/usr/bin/env python3
"""LM2009 遷移モデルのカーネル (`cuda_forge/transition_d.cu`) を numpy 参照実装 (`lm2009_reference.py`、倍精度) と**節点ごと**に突き合わせる。

`output: {level: 2}` で回した LM run の res_*.h5 から、カーネルが書いた診断場 (lmFonset, lmFlength, lmFtheta, lmRethCorr, lmGammaSep,
lmPgamma, lmEgamma, lmPtheta, src_jac_gamma, src_jac_reth, gammaEff, lmCorrIter) と、同じファイルの状態量から参照実装で再計算した値を比べる。

判定 (全部そろって PASS。codex result M3 で強化):
  * 各量について**全節点**で |カーネル − 参照| <= rtol·|参照| + atol·max|参照|  (既定 rtol 1e-3, atol 1e-5。float32 の相殺誤差ぶんの絶対項)
  * gammaEff = max(gamma, gamma_sep) が全節点で一致 (同じ絶対・相対許容)
  * 相関の不動点反復が上限 100 回に達した節点が 0
上下限の作動 (gamma が 1e-4 / 1、Re_theta_t が 20 に張り付いた節点の割合) は判定せず表示する。

注意: 診断場は**その step の更新前**の状態で計算され、出力の状態量は**更新後**なので、厳密には 1 step ぶんずれる。
準定常に達した場か、cfl を極小にした 1 step run で使うこと。

usage: check_lm_kernel.py RUN_DIR [--step N] [--rtol 1e-3] [--atol 1e-5]
"""
import argparse, glob, os, sys
import h5py, numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lm2009_reference as lm


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("run"); ap.add_argument("--step", type=int, default=None)
    ap.add_argument("--rtol", type=float, default=1e-3); ap.add_argument("--atol", type=float, default=1e-5)
    a = ap.parse_args()
    files = sorted(glob.glob(os.path.join(a.run, "res_[0-9]*.h5")), key=lambda f: int(os.path.basename(f)[4:-3]))
    f = os.path.join(a.run, f"res_{a.step}.h5") if a.step is not None else files[-1]
    with h5py.File(f, "r") as h:
        V = {k: np.array(h["VALUE"][k], dtype=float) for k in h["VALUE"]}
    V = check_inputs(V, f)
    if V is None:
        print("VERDICT: FAIL"); return 1
    return judge(V, f, a.rtol, a.atol)


NEED = ["ro", "Ux", "Uy", "Uz", "k", "omega", "sonic", "vis_lam", "wall_dist", "gammaTr", "reTheta", "gammaEff",
        "lmFonset", "lmFlength", "lmFtheta", "lmRethCorr", "lmGammaSep", "lmPgamma", "lmEgamma", "lmPtheta", "lmCorrIter",
        "src_jac_gamma", "src_jac_reth"] + [f"dU{i}d{j}" for i in "xyz" for j in "xyz"]


def check_inputs(V, f="(memory)"):
    """必須場の有無・形状・有限性を先に見る (codex result-2 M1: NaN や診断場の欠落が PASS になっていた)。不備なら None。"""
    lack = [n for n in NEED if n not in V]
    if lack:
        print(f"[check_lm_kernel] {f}: 必要な場が無い {lack} (output.level 2 の LM run が要る)"); return None
    n = len(V["ro"]); bad = [k for k in NEED if np.shape(V[k]) != (n,)]
    if bad or n == 0:
        print(f"[check_lm_kernel] {f}: 形状が合わない場 {bad} (n={n})"); return None
    nonfin = [k for k in NEED if not np.all(np.isfinite(V[k]))]
    if nonfin:
        print(f"[check_lm_kernel] {f}: 非有限値 (NaN/Inf) を含む場 {nonfin}"); return None
    return V


def judge(V, f, rtol, atol):
    class A: pass
    a = A(); a.rtol, a.atol = rtol, atol
    n = len(V["ro"])
    U = np.stack([V["Ux"], V["Uy"], V["Uz"]], -1)
    G = np.zeros((n, 3, 3))
    for i, ui in enumerate("xyz"):
        for j, xj in enumerate("xyz"):
            G[:, i, j] = V[f"dU{ui}d{xj}"]
    # 評価対象は**検査される側の出力でなく入力から**決める (codex result-2 M1): 壁距離 > 1e-10 かつ |U| >= 1e-6 a。
    # カーネルはこのほか入口ピン節点 (gamma = 1 に固定) も評価しない。入力側で対象なのにカーネルが評価していない節点は、
    # 「gamma = 1 の入口ピン」として説明できる数 (全体の 2 % 以下) でなければ FAIL。
    umag = np.linalg.norm(U, axis=-1)
    expected = (V["wall_dist"] > 1e-10) & (umag >= 1e-6 * V["sonic"])
    kern = V["lmFlength"] > 0.0
    skipped = expected & ~kern
    ok_mask = (skipped.sum() <= 0.02 * n) and bool(np.all(np.abs(V["gammaTr"][skipped] - 1.0) < 1e-6)) and not bool((kern & ~expected).any())
    active = expected & kern
    if active.sum() == 0:
        print("[check_lm_kernel] 評価対象の節点が 0"); print("VERDICT: FAIL"); return 1
    r = lm.sources(V["ro"], U, G, np.maximum(V["k"], 0.0), np.maximum(V["omega"], 1e-12), V["vis_lam"], np.maximum(V["wall_dist"], 1e-30), V["gammaTr"], V["reTheta"])
    pairs = [("F_onset", V["lmFonset"], r["f_onset"]), ("F_length", V["lmFlength"], r["f_length"]), ("F_theta", V["lmFtheta"], r["f_theta"]),
             ("Re_theta_corr", V["lmRethCorr"], r["re_theta_corr"]), ("gamma_sep", V["lmGammaSep"], r["gamma_sep"]),
             ("P_gamma", V["lmPgamma"], r["P_gamma"]), ("E_gamma", V["lmEgamma"], r["E_gamma"]), ("P_theta", V["lmPtheta"], r["src_reth"]),
             ("gammaEff", V["gammaEff"], r["gamma_eff"])]
    pairs += [("diag_gamma (forge)", V["src_jac_gamma"], r["jac_gamma"]), ("diag_Re_theta", V["src_jac_reth"], r["jac_reth"])]
    print(f"[check_lm_kernel] {f}: {int(active.sum())}/{n} 節点でソースを評価  (許容 rtol {a.rtol:g}, atol {a.atol:g}·max)")
    print(f"{'量':<20}{'max|参照|':>14}{'最大 |差|/許容':>16}{'超過節点':>10}  判定")
    ok = ok_mask
    print(f"評価対象: 入力から {int(expected.sum())} 節点、うちカーネルが評価しなかった {int(skipped.sum())} 節点 (入口ピン gamma=1 のはず)  {'ok' if ok_mask else 'NG'}")
    for name, ker, ref in pairs:
        k_, r_ = ker[active], ref[active]
        sc = max(np.abs(r_).max(), 1e-30); ratio = np.abs(k_ - r_) / (a.rtol * np.abs(r_) + a.atol * sc)
        nbad = int((~(ratio <= 1.0)).sum()); good = nbad == 0; ok &= good      # NaN も違反に数える
        print(f"{name:<20}{sc:14.5e}{ratio.max():16.3e}{nbad:10d}  {'ok' if good else 'NG'}")
    # ソースを評価しない節点では gammaEff = gamma のはず
    dge = np.abs(V["gammaEff"][~active] - V["gammaTr"][~active]).max() if (~active).any() else 0.0
    good = dge <= 1e-6; ok &= good
    print(f"非評価節点の gammaEff - gamma: 最大 {dge:.3e}  {'ok' if good else 'NG'}")
    if True:
        it = V["lmCorrIter"][active]; ncap = int((it >= 100).sum()); ok &= (ncap == 0)
        print(f"相関の不動点反復: 平均 {it.mean():.1f} 回, 最大 {it.max():.0f} 回, 上限 100 到達 {ncap} 節点  {'ok' if ncap == 0 else 'NG'}")
    g, t = V["gammaTr"], V["reTheta"]
    print(f"上下限の作動 (表示のみ): gamma<=1e-4 が {100*np.mean(g <= 1.0001e-4):.2f} %, gamma>=1 が {100*np.mean(g >= 1.0):.2f} %, Re_theta_t<=20 が {100*np.mean(t <= 20.0):.2f} %")
    if True:
        s2 = r["su2_jac_gamma"][active]; fo = r["jac_gamma"][active]; m = np.isfinite(s2)
        print(f"参考: forge の対角 / SU2 のソース微分の負部 = 中央値 {np.median(fo[m][-s2[m] > 1e-12] / (-s2[m])[-s2[m] > 1e-12]):.2f} (項ごとに負部を取るぶん forge が大きい。定常解には影響しない)")
    print("VERDICT:", "PASS" if ok else "FAIL"); return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
