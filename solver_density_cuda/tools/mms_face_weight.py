#!/usr/bin/env python3
r"""node (median-dual) の粘性流束で、**面補間重み fx の選び方が収束次数をどう変えるか**を製造解で測る。

plan discretization-node-face-weight-midpoint §5.1 #7 (codex plan M2 / result M2 / result-2 M1・M2)。
ソルバ本体は製造解のソース項を入れられないので、内部面の式は **`node_visc_face.visc_face()`** を使う。
これは `viscousFlux_d.cu` の既定経路を写したもので、実機 (case/53 C3X) の場に当てるとカーネルが積んだ
熱伝導・粘性仕事を 480 壁節点で 0.28 / 0.00 W/m² の差で再現する (`wallcv_flux_terms.py` と同じ関数)。

解く問題 (定常、**半径 0.2 m の円環の 15° 扇形**、内壁 = 曲面壁で $u=0$・$T=T_w$、外周と扇形の両端は厳密解の Dirichlet)。
半径は「弦のたるみ $\Delta s^2/8r$ / 第一層の半分」が最粗で 0.3 (実機 C3X は 0.1–1) になるように選んだ。
これを大きくしすぎる (たるみ > 第一層) と双対面が解析解の壁の内側に入り、厳密流束が意味を失う:
    熱      $-\nabla\cdot(k(T)\nabla T)=s_T$
    運動量  $-\nabla\cdot\tau(u)=s_u$,  $\tau=\mu(T)\,[\nabla u+\nabla u^\mathsf{T}-\tfrac23(\nabla\cdot u)I]$  (2 成分、接線 + 法線)
ソースは各双対面の厳密流束の Gauss 求積から作る。勾配 (非直交補正・転置項・発散項に入る) は解析値を使う
(重みの効果だけを切り出す)。

格子は実機に寄せる: 第一層 2 µm、壁沿い 0.69 mm (最粗で AR 約 350)、壁法線方向は等比、**壁沿い間隔に滑らかな変調
±30 % と節点ごとの交番 ±0.3 %** (面重心の接線ずれを作る。これが無いと旧式の回転非不変は現れない)。細分化は両方向 2 倍・
第一層 1/2・成長率は平方根。

測る量 (どれも誤差の rms と、隣り合う 2 水準からの観測次数 p):
    T, u           … 節点値の誤差
    wallQ, wallW   … **壁の検査体積に内部面から入る熱伝導 / 粘性仕事** (= ソルバの `iface_q_eff` を作る 2 項) の、
                     同じ面を通る厳密流束に対する誤差。面積あたり
    q2pt, tau2pt   … 2 点差分の壁面熱流束 / 壁せん断 (`iface_q_compact` と同じ作り) の、**解析的な壁面値**に対する誤差。
                     これは 2 点差分そのものが 1 次なので fx に依らず p≈1 になる (codex result-2 M1)
    q2pt_c         … 同じ 2 点差分を**厳密温度に当てた値**に対する誤差 (離散どうしの比較。旧版が「壁面熱流束」と呼んでいた量)

usage: python3 solver_density_cuda/tools/mms_face_weight.py [--levels 4] [--growth 1.1 1.2] [--out fig.png]
"""
import argparse
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from node_visc_face import visc_face, face_weight

R_IN, R_OUT = 0.20, 0.21
TH0, DTH = np.radians(15.0), np.radians(15.0)   # 扇形 15°–30° (旧式の回転非不変が最大になる 22.5° を中心に)
TW, DT, DELTA = 566.0, 200.0, 5.0e-5          # 壁温、温度差、製造解の境界層厚 [m]
U0, DELTA_U = 300.0, 5.0e-5


def polar(p):
    r = np.hypot(p[..., 0], p[..., 1]); th = np.arctan2(p[..., 1], p[..., 0])
    er = np.stack([np.cos(th), np.sin(th)], -1); et = np.stack([-np.sin(th), np.cos(th)], -1)
    return r, th, er, et


def fields(p):
    """厳密解: T, gradT, U (2), G[i,j]=dUi/dxj, k, mu。"""
    r, th, er, et = polar(p); s = r - R_IN
    e = np.exp(-s / DELTA); a = 1.0 + 0.2 * np.cos(24 * th); da = -4.8 * np.sin(24 * th)
    T = TW + DT * (1 - e) * a
    gT = (DT * e / DELTA * a)[..., None] * er + (DT * (1 - e) * da / r)[..., None] * et
    g = 1 - np.exp(-s / DELTA_U); dg = np.exp(-s / DELTA_U) / DELTA_U
    ut, ut_r, ut_t = U0 * g * a, U0 * dg * a, U0 * g * da
    b = np.sin(24 * th); db = 24 * np.cos(24 * th)
    ur, ur_r, ur_t = 0.02 * U0 * g * g * b, 0.04 * U0 * g * dg * b, 0.02 * U0 * g * g * db
    U = ur[..., None] * er + ut[..., None] * et
    dUdr = ur_r[..., None] * er + ut_r[..., None] * et
    dUdt = (ur_t - ut)[..., None] * er + (ut_t + ur)[..., None] * et
    G = dUdr[..., :, None] * er[..., None, :] + (dUdt / r[..., None])[..., :, None] * et[..., None, :]
    return T, gT, U, G, 0.04 * (T / TW) ** 0.75, 3.0e-5 * (T / TW) ** 0.7


GP = np.array([-0.8611363116, -0.3399810436, 0.3399810436, 0.8611363116])
GW = np.array([0.3478548451, 0.6521451549, 0.6521451549, 0.3478548451])


def seg_exact(m, c, ref):
    """線分 m->c を通る厳密流束 (熱 [W/m], 応力 [N/m] 2 成分, 仕事 [W/m])、面ベクトル、線分中点、長さ。"""
    t = c - m; L = np.linalg.norm(t, axis=-1)
    n = np.stack([t[..., 1], -t[..., 0]], -1) / L[..., None]
    n = n * np.sign(np.sum(n * ref, -1))[..., None]
    Fh = np.zeros(L.shape); Ft = np.zeros(L.shape + (2,)); Fw = np.zeros(L.shape)
    for gp, gw in zip(GP, GW):
        T, gT, U, G, k, mu = fields(m + 0.5 * (gp + 1.0) * t)
        tau = mu[..., None, None] * (G + np.swapaxes(G, -1, -2)
                                     - (2.0 / 3.0) * np.trace(G, axis1=-2, axis2=-1)[..., None, None] * np.eye(2))
        tn = np.einsum("...ij,...j->...i", tau, n); w = 0.5 * gw * L
        Fh += w * k * np.sum(gT * n, -1); Ft += w[..., None] * tn; Fw += w * np.sum(tn * U, -1)
    return Fh, Ft, Fw, n * L[..., None], 0.5 * (m + c), L


def mesh(level, growth, d1_0=2.0e-6, nth0=76):
    nth = nth0 * 2 ** level + 1; d1 = d1_0 / 2 ** level; g = growth ** (1.0 / 2 ** level)
    h = [d1]
    while sum(h) < R_OUT - R_IN:
        h.append(h[-1] * g)
    h = np.array(h); h *= (R_OUT - R_IN) / h.sum()
    r = R_IN + np.r_[0.0, np.cumsum(h)]
    i = np.arange(nth - 1)
    w = (1.0 + 0.3 * np.sin(2.0 * 2 * np.pi * (i + 0.5) / (nth - 1))) * (1.0 + 0.003 * (-1.0) ** i)
    th = TH0 + np.r_[0.0, np.cumsum(w)] / w.sum() * DTH
    return np.stack([r[None, :] * np.cos(th)[:, None], r[None, :] * np.sin(th)[:, None]], -1)   # [nth, nr, 2]


def build_edges(P):
    """辺 (A->B) ごとの幾何と厳密流束。扇形なので周期の折り返しは無い。壁の行 (j=0) の周方向の辺は片側の線分だけ。"""
    nth, nr, _ = P.shape
    C = 0.25 * (P[:-1, :-1] + P[:-1, 1:] + P[1:, :-1] + P[1:, 1:])      # quad(i,j) の重心 [nth-1, nr-1]
    I, J = np.meshgrid(np.arange(nth), np.arange(nr), indexing="ij")
    flat = lambda i, j: i * nr + j
    out = []
    # 半径方向 (i,j)->(i,j+1), i=1..nth-2: 隣接 quad (i-1,j), (i,j)
    a = (I[1:-1, :-1], J[1:-1, :-1]); b = (a[0], a[1] + 1)
    out.append((a, b, C[a[0] - 1, a[1]], C[a[0], a[1]]))
    # 周方向 (i,j)->(i+1,j), i=0..nth-2, j=1..nr-2: 隣接 quad (i,j-1), (i,j)
    a = (I[:-1, 1:-1], J[:-1, 1:-1]); b = (a[0] + 1, a[1])
    out.append((a, b, C[a[0], a[1] - 1], C[a[0], a[1]]))
    # 壁の行の周方向 (i,0)->(i+1,0): quad (i,0) だけ
    a = (I[:-1, :1], J[:-1, :1]); b = (a[0] + 1, a[1])
    out.append((a, b, None, C[a[0], a[1]]))
    E = dict(A=[], B=[], d=[], S=[], pc=[], Fh=[], Ft=[], Fw=[])
    for a, b, c1, c2 in out:
        xa, xb = P[a], P[b]; d = xb - xa; m = 0.5 * (xa + xb)
        Fh, Ft, Fw, S, pm, L = seg_exact(m, c2, d); pc = pm * L[..., None]; Ls = L.copy()
        if c1 is not None:
            Fh1, Ft1, Fw1, S1, pm1, L1 = seg_exact(m, c1, d)
            Fh, Ft, Fw, S = Fh + Fh1, Ft + Ft1, Fw + Fw1, S + S1; pc = pc + pm1 * L1[..., None]; Ls = Ls + L1
        pc = pc / Ls[..., None]
        for k_, v in (("A", flat(*a)), ("B", flat(*b)), ("d", d), ("S", S), ("pc", pc), ("Fh", Fh), ("Ft", Ft), ("Fw", Fw)):
            E[k_].append(v.reshape(-1, *v.shape[2:]))
    return {k_: np.concatenate(v) for k_, v in E.items()}, nth, nr


def solve_scalar(N, unk, A, B, coef_w, rhs_edge, known):
    """辺ごとの係数 c と既知右辺 b で  sum_e c (phi_B - phi_A) = sum_e b  を A に、逆符号を B に課す。"""
    rows, cols, vals = [], [], []; rhs = np.zeros(N)
    for me, ot, sg in ((A, B, 1.0), (B, A, -1.0)):
        sel = unk[me] >= 0
        r_, o_, c_ = unk[me][sel], ot[sel], coef_w[sel]
        rows.append(r_); cols.append(r_); vals.append(-c_)
        kn = unk[o_] < 0
        rows.append(r_[~kn]); cols.append(unk[o_][~kn]); vals.append(c_[~kn])
        np.add.at(rhs, r_, sg * rhs_edge[sel] - np.where(kn, c_ * known[o_], 0.0))
    M = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(N, N))
    return spla.spsolve(M.tocsc(), rhs)


def run(level, growth, mode):
    P = mesh(level, growth); E, nth, nr = build_edges(P)
    X = P.reshape(-1, 2); Te, gTe, Ue, Ge, ke, mue = fields(X)
    A, B = E["A"], E["B"]
    n = E["S"] / np.linalg.norm(E["S"], axis=-1)[:, None]
    f = face_weight(mode, n, E["pc"], X[A], X[B])
    jj = np.arange(nth * nr) % nr
    ii = np.arange(nth * nr) // nr
    unk = -np.ones(nth * nr, int); inner = (jj > 0) & (jj < nr - 1) & (ii > 0) & (ii < nth - 1)
    unk[inner] = np.arange(inner.sum()); N = inner.sum()
    args = lambda U, T: (X[B] - X[A], E["S"], f, U[A], U[B], T[A], T[B], Ge[A], Ge[B], gTe[A], gTe[B], mue[A], mue[B], ke[A], ke[B])
    # 既知項 = 厳密流束 - (勾配に依る補正項)。補正項は「節点差を 0 にした」流束として取り出す
    Z = np.zeros_like(Ue); ZT = np.zeros_like(Te)
    tau_c, heat_c, _, coef = visc_face(*args(Z, ZT))
    kf = f * ke[A] + (1 - f) * ke[B]; muf = f * mue[A] + (1 - f) * mue[B]
    Th = Te.copy(); Th[inner] = solve_scalar(N, unk, A, B, kf * coef, E["Fh"] - heat_c, Te)
    Uh = Ue.copy()
    for c in range(2):
        Uh[inner, c] = solve_scalar(N, unk, A, B, muf * coef, E["Ft"][:, c] - tau_c[:, c], Ue[:, c])
    # ---- 壁の検査体積に内部面から入る熱伝導・粘性仕事 (解いた場で評価) ----
    tau, heat, work, _ = visc_face(*args(Uh, Th))
    wall = (jj == 0) & (ii >= 2) & (ii <= nth - 3)      # 扇形の端から 2 節点は除く (端の Dirichlet の影響を見ない)
    def into_wall(F):
        acc = np.zeros(nth * nr); np.add.at(acc, A, F); np.add.at(acc, B, -F); return acc[wall]
    Pw = P[:, 0]; seg = np.linalg.norm(Pw[1:] - Pw[:-1], axis=-1); Aw = (0.5 * (seg[1:] + seg[:-1]))[1:-1]
    qh, qe = into_wall(heat) / Aw, into_wall(E["Fh"]) / Aw
    wh, we = into_wall(work) / Aw, into_wall(E["Fw"]) / Aw
    # ---- 2 点差分の壁面値 (iface_q_compact と同じ作り) ----
    sl = slice(2, nth - 2)
    d1 = np.linalg.norm(P[sl, 1] - P[sl, 0], axis=-1); r_, th_, er, et = polar(Pw[sl])
    q2 = ke.reshape(nth, nr)[sl, 0] * (Th.reshape(nth, nr)[sl, 1] - TW) / d1
    q2c = ke.reshape(nth, nr)[sl, 0] * (Te.reshape(nth, nr)[sl, 1] - TW) / d1
    qan = ke.reshape(nth, nr)[sl, 0] * np.sum(gTe.reshape(nth, nr, 2)[sl, 0] * er, -1)
    ut1 = np.sum(Uh.reshape(nth, nr, 2)[sl, 1] * et, -1)
    t2 = mue.reshape(nth, nr)[sl, 0] * ut1 / d1
    tan = mue.reshape(nth, nr)[sl, 0] * np.einsum("ni,nij,nj->n", et, Ge.reshape(nth, nr, 2, 2)[sl, 0], er)
    rms = lambda x: float(np.sqrt(np.mean(x ** 2)))
    return dict(n=nth, nr=nr, T=rms((Th - Te)[inner]), u=rms(np.linalg.norm(Uh - Ue, axis=-1)[inner]),
                wallQ=rms(qh - qe) / abs(qe).mean(), wallW=rms(wh - we) / abs(we).mean(),
                q2pt=rms((q2 - qan) / qan), q2pt_c=rms((q2 - q2c) / q2c), tau2pt=rms((t2 - tan) / tan),
                f_min=float(f[(jj[A] == 0) & (jj[B] == 1)].min()), f_max=float(f[(jj[A] == 0) & (jj[B] == 1)].max()),
                work_share=float(abs(we).mean() / abs(qe).mean()))


KEYS = ("T", "u", "wallQ", "wallW", "q2pt", "q2pt_c", "tau2pt")


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--levels", type=int, default=4)
    ap.add_argument("--growth", type=float, nargs="+", default=[1.1, 1.2]); ap.add_argument("--out")
    a = ap.parse_args(); res = {}
    for g in a.growth:
        print(f"\n== 成長率 {g} (最粗。細分化ごとに平方根)、第一層 2 µm / 2^L、壁沿い 76 x 2^L 区間 (0.69 mm / 2^L) ==")
        print(f"{'fx':<5}{'L':>2}{'nth x nr':>12}{'f at wall':>13} |" + "".join(f"{k:>11}{'p':>6}" for k in KEYS))
        for mode in ("half", "proj", "code"):
            prev = None
            for L in range(a.levels):
                r = run(L, g, mode); res[(g, mode, L)] = r
                cells = "".join(f"{r[k]:>11.3e}{(np.log2(prev[k] / r[k]) if prev else float('nan')):>6.2f}" for k in KEYS)
                print(f"{mode:<5}{L:>2}{r['n']:>6} x{r['nr']:>4}{r['f_min']:>7.3f}-{r['f_max']:.3f} |{cells}")
                prev = r
        print(f"   (壁の検査体積に入る粘性仕事 / 熱伝導 = {res[(g, 'half', 0)]['work_share']:.3f})")
    if a.out:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        g = a.growth[0]; hs = np.array([1.0 / 2 ** L for L in range(a.levels)])
        fig, ax = plt.subplots(1, 4, figsize=(15, 4.0))
        for x, k, lab in zip(ax, ("T", "u", "wallQ", "wallW"),
                             ("temperature error [K]", "velocity error [m/s]",
                              "conduction into the wall CV, rel. error", "viscous work into the wall CV, rel. error")):
            for mode, c in (("half", "k"), ("proj", "tab:blue"), ("code", "tab:red")):
                x.loglog(hs, [res[(g, mode, L)][k] for L in range(a.levels)], "o-", color=c, label=f"fx = {mode}")
            x.loglog(hs, res[(g, "half", 0)][k] * hs ** 2, "k:", lw=0.8, label="slope 2")
            x.loglog(hs, res[(g, "half", 0)][k] * hs, "k--", lw=0.6, label="slope 1")
            x.set_xlabel("relative mesh size"); x.set_title(lab, fontsize=9); x.grid(alpha=.3, which="both")
        ax[0].legend(fontsize=7); fig.suptitle(f"growth ratio {g}", fontsize=10); fig.tight_layout(); fig.savefig(a.out, dpi=110)


if __name__ == "__main__":
    main()
