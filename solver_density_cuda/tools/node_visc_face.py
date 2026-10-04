#!/usr/bin/env python3
r"""node (median-dual) の**内部面の粘性流束**を `viscousFlux_d.cu` と同じ式で評価する (numpy、ベクトル化)。

`heatCorrSU2 == 0` の既定経路 (over-relaxed 法線項 + 面平均勾配の補正):

    delta = |d||S|^2 / |d.S|,   k = S - d |S|^2/|d.S|
    tau   = mu_f [ (U1-U0)/|d| * delta + G_f k + G_f^T S - (2/3) tr(G_f) S ]
    heat  = k_f  [ (T1-T0)/|d| * delta + gT_f . k ]
    work  = tau . U_f
    ( )_f = f ( )_0 + (1-f) ( )_1

これを**2 か所から使う**ことで、製造解の試験が実カーネルと同じ式を測っていることを担保する:
  - `case/53.c3x_vane_cht/tools/wallcv_flux_terms.py` … 実機の場に当て、カーネルが積んだ `wi_eheat`/`wi_ework` と照合
  - `mms_face_weight.py` … 製造解の格子系列
面重み f の旧式・射影・中点は `face_weight()`。
"""
import numpy as np


def face_weight(mode, n, pc, x0, x1):
    """面補間重み (節点 0 側の重み)。mode: half | code (旧実装: 成分ごとの積のノルム) | proj (法線射影)。"""
    if mode == "half":
        return np.full(n.shape[:-1], 0.5)
    if mode == "code":
        d0 = np.linalg.norm(n * (pc - x0), axis=-1); d1 = np.linalg.norm(n * (pc - x1), axis=-1)
    elif mode == "proj":
        d0 = np.abs(np.sum(n * (pc - x0), -1)); d1 = np.abs(np.sum(n * (pc - x1), -1))
    else:
        raise ValueError(mode)
    return np.where(d0 + d1 > 0, d1 / np.maximum(d0 + d1, 1e-300), 0.5)


def visc_face(d, S, f, U0, U1, T0, T1, G0, G1, gT0, gT1, mu0, mu1, k0, k1):
    """内部面 1 枚ぶんの (tau [N], heat [W], work [W], coef)。ベクトルの最後の軸が空間成分。

    G[..., i, j] = dU_i/dx_j。戻り値の符号は「節点 0 の残差に + で入る向き」(カーネルの res[ic0] += と同じ)。
    coef は主項の係数 |S|^2/|d.S| (= delta/|d|)。陰的に解くときの行列係数に使う。
    """
    dcc = np.linalg.norm(d, axis=-1); ss = np.linalg.norm(S, axis=-1)
    D = np.maximum(np.abs(np.sum(d * S, -1)), 1e-300)
    coef = ss ** 2 / D
    kv = S - d * coef[..., None]
    fe = f[..., None]
    Uf = fe * U0 + (1 - fe) * U1
    Gf = fe[..., None] * G0 + (1 - fe[..., None]) * G1
    gTf = fe * gT0 + (1 - fe) * gT1
    mu = f * mu0 + (1 - f) * mu1; kk = f * k0 + (1 - f) * k1
    div = np.trace(Gf, axis1=-2, axis2=-1)
    tau = mu[..., None] * ((U1 - U0) * coef[..., None] + np.einsum("...ij,...j->...i", Gf, kv)
                           + np.einsum("...ji,...j->...i", Gf, S) - (2.0 / 3.0) * div[..., None] * S)
    heat = kk * ((T1 - T0) * coef + np.sum(gTf * kv, -1))
    work = np.sum(tau * Uf, -1)
    return tau, heat, work, coef
