#!/usr/bin/env python3
r"""壁半 CV に入る**内部面の粘性エネルギー流束を項別に再構成**する (viscousFlux_d の内部面の式を後処理で再現)。

level 2 + `dTdx` 出力の `res_N.h5` と壁ダンプ `res_wall_<phys>_N.h5` (interfaceDiag: 1) を読み、
各壁節点について内部面の寄与を
  heat_main = k_f (T1-T0)/|d| * delta            (over-relaxed の法線項)
  heat_corr = k_f (grad T)_f . (S - delta_vec)   (非直交補正項)
  work_main = U_f . [mu (U1-U0)/|d| * delta]
  work_corr = U_f . [mu (grad U)_f.k + mu (grad U)_f^T.S - 2/3 mu div u S]
に分けて積算する。和が壁ダンプの (Rraw - Rpre - Fw) と一致することを検算してから内訳を読む。
面値・面勾配は solver と同じ幾何重み fx (= 面重心から両節点への法線距離の比) で補間する。

usage: python3 case/53.c3x_vane_cht/tools/wallcv_flux_terms.py <run_dir> --step N [--out fig.png]
"""
import argparse, sys
from pathlib import Path
import numpy as np, h5py
sys.path.insert(0, str(Path(__file__).parent))
from compare_h import arc_map
from resid_split import smooth

CP, PRT = 1004.5, 0.9


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("run"); ap.add_argument("--step", required=True)
    ap.add_argument("--phys", default="5"); ap.add_argument("--out"); ap.add_argument("--f-half", action="store_true",
                    help="fx を 0.5 に固定して再構成 (幾何重みの影響を見る)")
    a = ap.parse_args(); rd = Path(a.run)
    with h5py.File(rd / "mesh.h5") as f:
        X = np.array(f["MESH/COORD"], float).reshape(-1, 3)
        nN = int(f["MESH"].attrs["nNormalPlanes"])
        pl = np.array(f["PLANES/STRUCT"])[:6 * nN].reshape(-1, 6); c0, c1 = pl[:, 4], pl[:, 5]
        S = np.array(f["PLANES/surfVect"], float).reshape(-1, 3)[:nN]
        pc = np.array(f["PLANES/centCoords"], float).reshape(-1, 3)[:nN]
        iw = np.array(f[f"BCONDS/{a.phys}/iCells"])
    with h5py.File(rd / f"res_{a.step}.h5") as f:
        V = {k: np.array(f["VALUE"][k], float) for k in f["VALUE"]}
    with h5py.File(rd / f"res_wall_{a.phys}_{a.step}.h5") as f:
        W = {k: np.array(f["VALUE"][k], float) for k in f["VALUE"]}
    isw = np.zeros(len(X), bool); isw[iw] = True
    sel = np.where(isw[c0] | isw[c1])[0]
    i0, i1, Sv, P = c0[sel], c1[sel], S[sel], pc[sel]
    ss_ = np.linalg.norm(Sv, axis=1); n = Sv / ss_[:, None]
    d0 = np.abs(((P - X[i0]) * n).sum(1)); d1 = np.abs(((P - X[i1]) * n).sum(1))
    fx = np.where(d0 + d1 > 0, d1 / (d0 + d1), 0.5)
    if a.f_half: fx[:] = 0.5
    F = lambda q: fx * V[q][i0] + (1 - fx) * V[q][i1]
    d = X[i1] - X[i0]; dcc = np.linalg.norm(d, axis=1)
    D = np.maximum(np.abs((d * Sv).sum(1)), 1e-30)
    delta = dcc * ss_ ** 2 / D; dvec = d * (ss_ ** 2 / D)[:, None]; kv = Sv - dvec
    mu = F("vis_lam") + F("vis_turb"); tc = F("thermCond") + CP * F("vis_turb") / PRT
    gT = np.stack([F("dTdx"), F("dTdy"), F("dTdz")], 1)
    heat_main = tc * (V["T"][i1] - V["T"][i0]) / dcc * delta
    heat_corr = tc * (gT * kv).sum(1)
    U = "Ux Uy Uz".split(); X3 = "xyz"
    G = np.stack([np.stack([F(f"d{u}d{x}") for x in X3], 1) for u in U], 1)      # G[:,i,j] = dUi/dxj
    Uf = np.stack([F(u) for u in U], 1); dU = np.stack([V[u][i1] - V[u][i0] for u in U], 1)
    divu = G[:, 0, 0] + G[:, 1, 1] + G[:, 2, 2]
    tau_main = mu[:, None] * dU / dcc[:, None] * delta[:, None]
    tau_corr = mu[:, None] * (np.einsum("fij,fj->fi", G, kv) + np.einsum("fji,fj->fi", G, Sv) - (2 / 3) * divu[:, None] * Sv)
    work_main = (tau_main * Uf).sum(1); work_corr = (tau_corr * Uf).sum(1)
    terms = {"heat_main": heat_main, "heat_corr": heat_corr, "work_main": work_main, "work_corr": work_corr}
    acc = {k: np.zeros(len(X)) for k in terms}
    for k, v in terms.items():
        np.add.at(acc[k], i0, v); np.add.at(acc[k], i1, -v)
    A = (W["ifaceRraw"] - W["ifaceFw"]) / W["iface_q_eff"]
    tgt = (W["ifaceRraw"] - W["ifaceRpre"] - W["ifaceFw"]) / A
    T = {k: acc[k][iw] / A for k in terms}; tot = sum(T.values())
    s, ss = arc_map(X[iw][:, :2]); o = np.argsort(s[ss]); sS = s[ss][o]
    g = lambda y: y[ss][o]
    print(f"[wallcv_flux_terms] {rd}/res_{a.step}.h5  fx={'0.5 fixed' if a.f_half else 'geometric'}")
    print(f"  check: reconstructed - solver (Rraw-Rpre-Fw)/A : rms {np.sqrt(((tot-tgt)[ss]**2).mean()):.1f} W/m2 "
          f"(solver mean {tgt[ss].mean():.0f})")
    for lo, hi in ((0.30, 0.40), (0.45, 0.95)):
        m = (sS >= lo) & (sS <= hi)
        print(f"  s/S {lo}-{hi}:   mean [W/m2]   wave rms (y - smooth41)   corr with q_eff wave")
        qw = g(W["iface_q_eff"]) - smooth(g(W["iface_q_eff"]), 41)
        for k in list(T) + ["total", "solver"]:
            y = g(T[k]) if k in T else g(tot if k == "total" else tgt); w = y - smooth(y, 41)
            print(f"    {k:<10} {y[m].mean():12.1f} {np.sqrt((w[m]**2).mean()):12.1f} {np.corrcoef(w[m], qw[m])[0,1]:+12.2f}")
    if a.out:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
        for k in T:
            ax[0].plot(sS, (g(T[k]) - smooth(g(T[k]), 41)) / 1e3, lw=1, label=k)
        ax[0].plot(sS, (g(tgt) - smooth(g(tgt), 41)) / 1e3, "k", lw=1, label="solver interior viscous")
        ax[1].plot(sS, g(fx_wall(fx, isw, i0, i1, dcc, iw)), lw=1, label="wall-side weight fx on the wall-normal face")
        ax[1].axhline(0.5, color="gray", lw=.5)
        for x in ax: x.legend(loc="upper left", bbox_to_anchor=(1.0, 1.0), fontsize=8); x.grid(alpha=.3)
        ax[0].set_ylabel("detrended [kW/m$^2$]"); ax[0].set_ylim(-6, 6); ax[1].set_xlabel("s/S (suction side)"); ax[0].set_xlim(0.2, 1)
        fig.tight_layout(); fig.savefig(a.out, dpi=110)


def fx_wall(fx, isw, i0, i1, dcc, iw):
    out = {}
    for k in np.where(isw[i0] ^ isw[i1])[0]:
        w = i0[k] if isw[i0[k]] else i1[k]; fw = fx[k] if isw[i0[k]] else 1 - fx[k]
        if w not in out or dcc[k] < out[w][0]: out[w] = (dcc[k], fw)
    return np.array([out[w][1] for w in iw])


if __name__ == "__main__":
    main()
