#!/usr/bin/env python3
"""後縁の冷点の LSQ 勾配を forge の gradLSQ=2 と同じ式で再現し、近傍ごとの寄与と辺ごとの比を出す (読むだけ)。

forge: d_ij = x_j - x_i (node モードは centCoords を節点座標に置換 = mesh.cpp nodeValueAtNode), w = 1/|d|^2, M = Σ w d dᵀ,
       M⁺ = スペクトル打ち切り擬似逆 (λ < 1e-2 λmax を落とす), g = M⁺ Σ w d (φ_j - φ_i)。
       内部面の隣接のみ (境界半割面は入れない)。再構成 u_f = u_i + ψ g·(d/2)。
使い方: te_grad_diag.py mesh.h5 res.h5 ledger.csv.faces node1 [node2 ...]
"""
import sys
import numpy as np
import h5py

mesh_h5, res_h5, faces_csv = sys.argv[1:4]
targets = [int(a) for a in sys.argv[4:]]

with h5py.File(mesh_h5, "r") as f:
    cc = f["CELLS/centCoords"][...].astype(np.float64).reshape(-1, 3)
    xn = f["MESH/COORD"][...].astype(np.float64).reshape(-1, 3)
    vol = f["CELLS/volume"][...].astype(np.float64)
    st = f["PLANES/STRUCT"][...]
    sv = f["PLANES/surfVect"][...].astype(np.float64).reshape(-1, 3)
    sa = f["PLANES/surfArea"][...].astype(np.float64)
N = cc.shape[0]
with h5py.File(res_h5, "r") as f:
    V = f["VALUE"]
    ro = V["ro"][...].astype(np.float64)
    U = {c: V["ro" + c][...].astype(np.float64) / ro for c in ("Ux", "Uy", "Uz")}
    T = V["T"][...].astype(np.float64) if "T" in V else None

# PLANES/STRUCT: [nn, nodes..., nc, cells...] を順に
tset = set(targets)
planes = {}   # node -> list of (ip, other, sign, internal)
pos = 0
nP = sa.shape[0]
stl = st.tolist()
for ip in range(nP):
    nn = stl[pos]; pos += 1 + nn
    nc = stl[pos]; cells = stl[pos + 1: pos + 1 + nc]; pos += 1 + nc
    for k, c in enumerate(cells):
        if c in tset:
            other = cells[1 - k] if nc == 2 else -1
            internal = (nc == 2 and other < N and other >= 0)
            planes.setdefault(c, []).append((ip, other, 1.0 if k == 0 else -1.0, internal))

led = {}
try:
    import csv
    for r in csv.DictReader(open(faces_csv)):
        if r["call"] != "1":
            continue
        led[int(r["ip"])] = r
except FileNotFoundError:
    pass


def pinv(M, thresh=1e-2):
    lam, vec = np.linalg.eigh(M)
    keep = lam >= thresh * lam[-1]
    Mi = sum(np.outer(vec[:, k], vec[:, k]) / lam[k] for k in range(3) if keep[k])
    return Mi, lam, int((~keep).sum())


np.set_printoptions(linewidth=200, precision=4, suppress=False)
for i in targets:
    print("=" * 100)
    print(f"node {i}: x_node {xn[i]}  cc {cc[i]}  |cc-x| {np.linalg.norm(cc[i]-xn[i]):.3e}  vol {vol[i]:.3e}"
          + (f"  T {T[i]:.1f}" if T is not None else ""))
    print(f"  U_i = ({U['Ux'][i]:.2f}, {U['Uy'][i]:.2f}, {U['Uz'][i]:.2f})")
    nb = [(ip, j, s) for (ip, j, s, internal) in planes.get(i, []) if internal]
    bd = [(ip, j, s) for (ip, j, s, internal) in planes.get(i, []) if not internal]
    M = np.zeros((3, 3)); D = []
    for ip, j, s in nb:
        d = xn[j] - xn[i]; w = 1.0 / d.dot(d)
        M += w * np.outer(d, d); D.append((ip, j, d, w))
    Mi, lam, ndeg = pinv(M)
    print(f"  M eig {lam}  cond {lam[-1]/max(lam[0],1e-300):.3e}  dropped {ndeg}  boundary half-faces {len(bd)}")
    # GG: (1/V) Σ φ_f S_f (内部面 = 算術平均、境界半割面 = 節点値)
    for comp in ("Ux", "Uy", "Uz"):
        phi = U[comp]
        g = np.zeros(3); contrib = []
        for ip, j, d, w in D:
            c = Mi @ (w * d) * (phi[j] - phi[i]); g += c; contrib.append((ip, j, d, c))
        gg = np.zeros(3)
        for ip, j, s in nb:
            gg += 0.5 * (phi[i] + phi[j]) * s * sv[ip] * sa[ip]
        for ip, j, s in bd:
            gg += phi[i] * s * sv[ip] * sa[ip]
        gg /= vol[i]
        print(f"  [{comp}] g_LSQ {g}  g_GG {gg}")
        print(f"    {'nbr':>8} {'ip':>8} {'d (m)':>36} {'|d|':>9} {'phi_j':>9} {'dphi':>9} {'g.d':>10} {'ratio':>7} {'GG.d':>10} {'contrib(c_ij dphi).d_k for k=this edge':>12}  ledger: L-recon, psi_eff")
        for ip, j, d, c in contrib:
            dphi = phi[j] - phi[i]
            gd = g.dot(d); ggd = gg.dot(d)
            ratio = gd / dphi if abs(dphi) > 1e-9 else float("nan")
            extra = ""
            r = led.get(ip)
            if r is not None:
                side = "L" if int(r["ic0"]) == i else "R"
                uf = float(r[f"{comp}_{side}"])
                psi = (uf - phi[i]) / (0.5 * gd) if abs(gd) > 1e-9 else float("nan")
                extra = f"  {uf:10.2f}  psi_eff {psi:7.4f}"
            print(f"    {j:>8} {ip:>8} {np.array2string(d, precision=3, separator=','):>36} {np.linalg.norm(d):9.3e} {phi[j]:9.2f} {dphi:9.2f} {gd:10.2f} {ratio:7.2f} {ggd:10.2f}{extra}")
        # どの近傍が各辺方向の g.d を作っているか
        print("    寄与行列 (行 = 近傍 j の寄与 c_j dphi_j を、列 = 各辺 d_k へ射影):")
        hdr = " ".join(f"{k:>9}" for _, k, _, _ in contrib)
        print(f"    {'from\\to':>8} {hdr}")
        for ip, j, d, c in contrib:
            row = " ".join(f"{c.dot(dk):9.1f}" for _, _, dk, _ in contrib)
            print(f"    {j:>8} {row}")
