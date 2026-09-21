#!/usr/bin/env python3
r"""node (median-dual) の拡散演算子で、**面補間重み fx の選び方が解の収束次数をどう変えるか**を製造解で測る。

plan discretization-node-face-weight-midpoint §5.1 #7 (codex plan M2 / result M2)。ソルバ本体は製造解のソース項を
入れられないので、`viscousFlux_d.cu` の内部面の熱伝導の式

    F = k_f [ (T_B - T_A)/|d| * delta + g_f . (S - delta_vec) ],   delta = |d||S|^2/(d.S),
    k_f = f k_A + (1-f) k_B,   g_f = f g_A + (1-f) g_B

を**同じ形で** Python に写し、円環 (内壁 = 曲面壁) の O 型構造格子の双対メッシュ上で定常問題
$-\nabla\cdot(k\nabla T)=s$ を解く。格子は実機 (case/53 C3X) に寄せる: 第一層 2 µm、壁沿い 0.7 mm (AR 350)、
壁法線方向は等比 (成長率 1.15)、**壁沿い間隔は滑らかな変調 + 節点ごとの交番** (面重心の接線ずれを作る。
これが無いと回転非不変の欠陥は現れない)。円環なので壁の向きは全方位を含む。

fx の 3 通りを比べる:
  half : 0.5 (辺中点。現行 = 2026-09-22 以降)
  code : calcStructualVariables_d.cu の旧式  d = || n (.) Delta ||  (成分ごとの積のノルム。回転不変でない)
  proj : 正しい射影  d = | n . Delta |

細分化は両方向 2 倍、第一層 1/2、成長率は平方根 (層の物理的な厚さ分布を保つ)。勾配 g は解析値を使う
(重みの効果だけを切り出すため。O 型格子では非直交補正項自体が小さい)。ソースは各双対面の厳密流束を
Gauss 求積して作るので、ソースの求積誤差は入らない。

usage: python3 solver_density_cuda/tools/mms_face_weight.py [--levels 4] [--growth 1.15] [--out fig.png]
"""
import argparse
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

R_IN, R_OUT = 0.05, 0.06
TW, DT, DELTA = 566.0, 200.0, 5.0e-5          # 壁温、温度差、製造解の境界層厚 [m]


def T_exact(x, y):
    r = np.hypot(x, y); th = np.arctan2(y, x)
    return TW + DT * (1.0 - np.exp(-(r - R_IN) / DELTA)) * (1.0 + 0.2 * np.cos(3.0 * th))


def gradT_exact(x, y):
    r = np.hypot(x, y); th = np.arctan2(y, x)
    e = np.exp(-(r - R_IN) / DELTA); a = 1.0 + 0.2 * np.cos(3.0 * th)
    dTdr = DT * e / DELTA * a
    dTdth = DT * (1.0 - e) * (-0.6 * np.sin(3.0 * th))
    er = np.stack([np.cos(th), np.sin(th)], -1); et = np.stack([-np.sin(th), np.cos(th)], -1)
    return dTdr[..., None] * er + (dTdth / r)[..., None] * et


def k_of_T(T):                                 # 空気の k(T) に近い強さの温度依存 (T^0.75)
    return 0.04 * (T / TW) ** 0.75


def mesh(level, growth, d1_0=2.0e-6, nth0=225):
    nth = nth0 * 2 ** level; d1 = d1_0 / 2 ** level; g = growth ** (1.0 / 2 ** level)
    # 壁法線: 等比で R_OUT を超えるまで積み、最後を R_OUT に合わせて一様に縮める
    h = [d1]
    while sum(h) < R_OUT - R_IN:
        h.append(h[-1] * g)
    h = np.array(h); h *= (R_OUT - R_IN) / h.sum()
    r = R_IN + np.r_[0.0, np.cumsum(h)]
    # 壁沿い: 滑らかな変調 (±30 %) + 節点ごとの交番 (±0.3 %)。細分化しても同じ相対振幅
    i = np.arange(nth)
    w = (1.0 + 0.3 * np.sin(5.0 * 2 * np.pi * (i + 0.5) / nth)) * (1.0 + 0.003 * (-1.0) ** i)
    th = np.r_[0.0, np.cumsum(w)][:-1] / w.sum() * 2 * np.pi + 0.123
    X = r[None, :] * np.cos(th)[:, None]; Y = r[None, :] * np.sin(th)[:, None]
    return X, Y                                  # [nth, nr]


GP = np.array([-0.8611363116, -0.3399810436, 0.3399810436, 0.8611363116])
GW = np.array([0.3478548451, 0.6521451549, 0.6521451549, 0.3478548451])


def seg_exact_flux(m, c, nvec_sign_ref):
    """線分 m->c を通る厳密流束 ∫ k gradT . n ds (n は nvec_sign_ref 側を向く単位法線)。"""
    t = c - m; L = np.linalg.norm(t, axis=-1)
    n = np.stack([t[..., 1], -t[..., 0]], -1) / L[..., None]
    sgn = np.sign(np.sum(n * nvec_sign_ref, -1)); n = n * sgn[..., None]
    F = np.zeros(L.shape)
    for gp, gw in zip(GP, GW):
        p = m + (0.5 * (gp + 1.0))[..., None] * t if np.ndim(gp) else m + 0.5 * (gp + 1.0) * t
        T = T_exact(p[..., 0], p[..., 1]); g = gradT_exact(p[..., 0], p[..., 1])
        F += 0.5 * gw * L * k_of_T(T) * np.sum(g * n, -1)
    return F, n * L[..., None], 0.5 * (m + c), L


def solve(level, growth, mode):
    X, Y = mesh(level, growth); nth, nr = X.shape
    P = np.stack([X, Y], -1)
    Pn = np.roll(P, -1, axis=0)                                  # θ 方向の次の節点
    C = 0.25 * (P[:, :-1] + P[:, 1:] + Pn[:, :-1] + Pn[:, 1:])   # quad(i,j) の重心 [nth, nr-1]
    Te = T_exact(X, Y); ke = k_of_T(Te); ge = gradT_exact(X, Y)
    idx = -np.ones((nth, nr), int); inner = np.zeros((nth, nr), bool); inner[:, 1:-1] = True
    idx[inner] = np.arange(inner.sum()); N = inner.sum()
    rows, cols, vals = [], [], []; rhs = np.zeros(N)

    def add_edges(A_ij, B_ij, c1, c2):
        """辺 A->B (節点添字の組) と両隣の quad 重心 c1, c2 から面を作って系に足す。"""
        nonlocal rows, cols, vals, rhs
        xa, xb = P[A_ij], P[B_ij]; d = xb - xa; m = 0.5 * (xa + xb)
        F1, S1, pm1, L1 = seg_exact_flux(m, c1, d); F2, S2, pm2, L2 = seg_exact_flux(m, c2, d)
        Fex = F1 + F2; S = S1 + S2; pc = (pm1 * L1[..., None] + pm2 * L2[..., None]) / (L1 + L2)[..., None]
        ss = np.linalg.norm(S, axis=-1); n = S / ss[..., None]; dS = np.sum(d * S, -1)
        if mode == "half":
            f = np.full(ss.shape, 0.5)
        else:
            if mode == "code":
                d0 = np.linalg.norm(n * (pc - xa), axis=-1); d1_ = np.linalg.norm(n * (pc - xb), axis=-1)
            else:
                d0 = np.abs(np.sum(n * (pc - xa), -1)); d1_ = np.abs(np.sum(n * (pc - xb), -1))
            f = d1_ / (d0 + d1_)
        kf = f * ke[A_ij] + (1 - f) * ke[B_ij]
        gf = f[..., None] * ge[A_ij] + (1 - f)[..., None] * ge[B_ij]
        coef = kf * ss ** 2 / dS                                   # (T_B - T_A) の係数
        corr = kf * np.sum(gf * (S - d * (ss ** 2 / dS)[..., None]), -1)
        ia, ib = idx[A_ij], idx[B_ij]; TA, TB = Te[A_ij], Te[B_ij]
        # 節点 A の収支: sum_e F_e(discrete) = sum_e F_e(exact)。B は逆符号
        for (me, other, sgn, Tother_known) in ((ia, ib, +1.0, TB), (ib, ia, -1.0, TA)):
            sel = me >= 0
            r_ = me[sel]; o_ = other[sel]; c_ = coef[sel]
            # sgn*[c (T_o - T_me)*s + corr] = sgn*Fex, s = +1 (A 視点) / -1 (B 視点では T_A - T_B = -(T_B-T_A))
            # 離散流束 (A->B 向き) は c (T_B - T_A) + corr。A には +、B には - で入る。
            rows.append(r_); cols.append(r_); vals.append(-c_ if sgn > 0 else -c_)
            known = o_ < 0
            rows.append(r_[~known]); cols.append(o_[~known]); vals.append(c_[~known])
            b = (Fex - corr)[sel] * sgn
            # A 視点: c(T_B - T_A) = Fex - corr。B 視点: c(T_A - T_B) = -(Fex - corr)  → 同じ形 c(T_o - T_me) = sgn*(Fex-corr)
            b = b - np.where(known, c_ * Tother_known[sel], 0.0)
            np.add.at(rhs, r_, b)

    I, J = np.meshgrid(np.arange(nth), np.arange(nr), indexing="ij")
    Im = (I - 1) % nth; Ip = (I + 1) % nth
    # 半径方向の辺 (i,j)->(i,j+1): 隣接 quad は (i-1,j) と (i,j)
    a = (I[:, :-1], J[:, :-1]); b = (I[:, :-1], J[:, :-1] + 1)
    add_edges(a, b, C[Im[:, :-1], J[:, :-1]], C[I[:, :-1], J[:, :-1]])
    # 周方向の辺 (i,j)->(i+1,j), j=1..nr-2: 隣接 quad は (i,j-1) と (i,j)
    a = (I[:, 1:-1], J[:, 1:-1]); b = (Ip[:, 1:-1], J[:, 1:-1])
    add_edges(a, b, C[I[:, 1:-1], J[:, 1:-1] - 1], C[I[:, 1:-1], J[:, 1:-1]])
    A = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(N, N))
    Th = Te.copy(); Th[inner] = spla.spsolve(A.tocsc(), rhs)
    err = (Th - Te)[:, 1:-1]
    # 壁面の勾配流束 (ソルバの iface_q_compact と同じ作り) の誤差: 厳密な k dT/dr との相対差
    d1 = np.hypot(X[:, 1] - X[:, 0], Y[:, 1] - Y[:, 0])
    q_h = ke[:, 0] * (Th[:, 1] - Th[:, 0]) / d1; q_e0 = ke[:, 0] * (Te[:, 1] - Te[:, 0]) / d1
    return dict(n=nth, nr=nr, l2=float(np.sqrt(np.mean(err ** 2))), linf=float(np.abs(err).max()),
                t1=float(np.sqrt(np.mean((Th[:, 1] - Te[:, 1]) ** 2))),
                qrel=float(np.sqrt(np.mean(((q_h - q_e0) / q_e0) ** 2))))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--levels", type=int, default=4)
    ap.add_argument("--growth", type=float, nargs="+", default=[1.1, 1.2]); ap.add_argument("--out")
    a = ap.parse_args()
    res = {}
    for g in a.growth:
        print(f"\n== 成長率 {g} (最粗格子。細分化ごとに平方根)、第一層 2 µm / 2^L、壁沿い 225 x 2^L 節点 ==")
        print(f"{'fx':<6}{'L':>2}{'nth x nr':>13}{'T err L2 [K]':>15}{'p':>6}{'T err Linf':>13}{'p':>6}"
              f"{'T1 err [K]':>13}{'p':>6}{'q_w err rel':>13}{'p':>6}")
        for mode in ("half", "proj", "code"):
            prev = None
            for L in range(a.levels):
                r = solve(L, g, mode); res[(g, mode, L)] = r
                p = {k: (np.log2(prev[k] / r[k]) if prev else float("nan")) for k in ("l2", "linf", "t1", "qrel")}
                print(f"{mode:<6}{L:>2}{r['n']:>7} x{r['nr']:>4}{r['l2']:>15.3e}{p['l2']:>6.2f}{r['linf']:>13.3e}{p['linf']:>6.2f}"
                      f"{r['t1']:>13.3e}{p['t1']:>6.2f}{r['qrel']:>13.3e}{p['qrel']:>6.2f}")
                prev = r
    if a.out:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, len(a.growth), figsize=(5.2 * len(a.growth), 4.2), squeeze=False)
        for x, g in zip(ax[0], a.growth):
            for mode, c in (("half", "k"), ("proj", "tab:blue"), ("code", "tab:red")):
                hs = [1.0 / 2 ** L for L in range(a.levels)]
                x.loglog(hs, [res[(g, mode, L)]["l2"] for L in range(a.levels)], "o-", color=c, label=f"fx = {mode}")
            x.loglog(hs, [res[(g, "half", 0)]["l2"] * h ** 2 for h in hs], "k:", lw=0.8, label="slope 2")
            x.set_xlabel("relative mesh size"); x.set_ylabel("temperature error, L2 [K]"); x.set_title(f"growth {g}")
            x.grid(alpha=.3, which="both"); x.legend(fontsize=8)
        fig.tight_layout(); fig.savefig(a.out, dpi=110)


if __name__ == "__main__":
    main()
