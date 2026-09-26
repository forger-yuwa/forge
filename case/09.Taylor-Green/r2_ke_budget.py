#!/usr/bin/env python3
"""R2 (plan boundary-node-periodic-gradient-fix §5.1 #5d) の運動エネルギー収支。

    python3 r2_ke_budget.py OLD_RUN NEW_RUN [--png r2_ke_budget.png] [--csv r2_ke_budget.csv]

一意 DOF (33³ の node 格子から周期重複の最後の層を落とした 32³) で、各 snapshot について
  K    = Σ ½ρ|u|² h³
  ε    = Σ τ_ij ∂_j u_i h³        (τ_ij = μ(∂_j u_i + ∂_i u_j − ⅔δ_ij ∇·u)、μ 一定)
  Π    = Σ p ∇·u h³               (圧力の膨張仕事。M0=0.4 なので無視できない)
を出し、収支残差を K0/t_c で正規化する:
  r_ε  = (dK/dt + ε) / (K0/t_c)             (diagnostician の定義)
  r    = (dK/dt − Π + ε) / (K0/t_c)         (圧縮性の収支。判定はこちら)
微分は周期 2 次中心差分 (solver の離散化とは別物なので、新でも r は 0 にならない)。dK/dt は snapshot 間の中心差分。

判定 (測る前に固定、plan §5.1 #5d): 新は |r| ≤ 0.05。旧は t ≲ 7 で r > 0 が持続すれば「継ぎ目由来の注入」、
符号不定なら「収支が閉じない」に留める。
"""
import argparse
import glob
import os
import re

import h5py
import numpy as np

MU = 2.5e-4
DT = 0.007
L = 2.0 * np.pi
U0 = 0.4
TC = 1.0 / U0          # t_c = L_ref/U0 (L_ref = 1)


def step_of(p):
    m = re.search(r"res_(\d+)\.h5$", os.path.basename(p))
    return int(m.group(1)) if m else -1


def grid_order(coord):
    c = coord.reshape(-1, 3).astype(np.float64)
    n1 = round(len(c) ** (1 / 3))
    h = L / (n1 - 1)
    idx = np.rint(c / h).astype(int)
    assert idx.min() == 0 and idx.max() == n1 - 1, (idx.min(), idx.max())
    order = np.full((n1, n1, n1), -1)
    order[idx[:, 0], idx[:, 1], idx[:, 2]] = np.arange(len(c))
    assert (order >= 0).all()
    return order[:-1, :-1, :-1], h          # 周期重複の最後の層を落とす


def ddx(a, ax, h):
    return (np.roll(a, -1, ax) - np.roll(a, 1, ax)) / (2 * h)


def budget(run):
    files = sorted(glob.glob(os.path.join(run, "res_*.h5")), key=step_of)
    order, h = None, None
    rows = []
    for p in files:
        with h5py.File(p, "r") as f:
            if order is None:
                order, h = grid_order(f["MESH/COORD"][()])
            g = lambda k: f["VALUE/" + k][()].astype(np.float64)[order]
            ro, P = g("ro"), g("P")
            u = [g("roUx") / ro, g("roUy") / ro, g("roUz") / ro]
        du = [[ddx(u[i], j, h) for j in range(3)] for i in range(3)]      # du[i][j] = ∂_j u_i
        div = du[0][0] + du[1][1] + du[2][2]
        eps = 0.0
        for i in range(3):
            for j in range(3):
                tau = MU * (du[i][j] + du[j][i] - (2.0 / 3.0 if i == j else 0.0) * div)
                eps += np.sum(tau * du[i][j])
        K = np.sum(0.5 * ro * (u[0] ** 2 + u[1] ** 2 + u[2] ** 2))
        rows.append((step_of(p) * DT, K * h ** 3, eps * h ** 3, np.sum(P * div) * h ** 3))
    a = np.array(rows)
    t, K, eps, Pi = a.T
    dKdt = np.gradient(K, t)
    K0 = K[0]
    r_eps = (dKdt + eps) / (K0 / TC)
    r = (dKdt - Pi + eps) / (K0 / TC)
    return t, K / K0, eps / (K0 / TC), Pi / (K0 / TC), r_eps, r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("old")
    ap.add_argument("new")
    ap.add_argument("--png", default=None)
    ap.add_argument("--csv", default=None)
    a = ap.parse_args()
    res = {"old": budget(a.old), "new": budget(a.new)}
    out = []
    for tag in ("old", "new"):
        t, Kn, en, Pin, r_eps, r = res[tag]
        # 両端は片側差分なので除く
        s = slice(1, -1)
        early = (t > 0) & (t <= 7.0)
        early_i = early[s]
        rr = r[s]
        line = (f"{tag}: |r| max {np.abs(rr).max():.4f} (t={t[s][np.argmax(np.abs(rr))]:.2f})、"
                f"|r_eps| max {np.abs(r_eps[s]).max():.4f}、max K/K0 {Kn.max():.5f}、"
                f"t≤7 の r: 平均 {rr[early_i].mean():+.4f} / 正の割合 {np.mean(rr[early_i] > 0):.2f} / "
                f"min {rr[early_i].min():+.4f} max {rr[early_i].max():+.4f}")
        out.append(line)
    t_new, _, _, _, _, r_new = res["new"]
    ok = np.abs(r_new[1:-1]).max() <= 0.05
    out.append(f"VERDICT (新 |r| ≤ 0.05): {'PASS' if ok else 'FAIL'}")
    print("\n".join(out))
    if a.csv:
        with open(a.csv, "w") as f:
            f.write("run,t,K_over_K0,eps_n,Pi_n,r_eps,r\n")
            for tag in ("old", "new"):
                for row in zip(*res[tag]):
                    f.write(tag + "," + ",".join(f"{v:.8g}" for v in row) + "\n")
    if a.png:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
        for tag, c in (("old", "C3"), ("new", "C0")):
            t, Kn, en, Pin, r_eps, r = res[tag]
            ax[0].plot(t, Kn, c, label=f"{tag} K/K0")
            ax[1].plot(t[1:-1], r[1:-1], c, label=f"{tag} r (with p div u)")
            ax[1].plot(t[1:-1], r_eps[1:-1], c + "--", lw=0.8, label=f"{tag} r_eps (no p div u)")
        ax[1].axhline(0.05, color="k", lw=0.5, ls=":")
        ax[1].axhline(-0.05, color="k", lw=0.5, ls=":")
        ax[0].set_ylabel("K/K0")
        ax[1].set_ylabel("(dK/dt - Pi + eps)/(K0/t_c)")
        ax[1].set_xlabel("t")
        for x in ax:
            x.legend(fontsize=8)
            x.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(a.png, dpi=130)


if __name__ == "__main__":
    main()
