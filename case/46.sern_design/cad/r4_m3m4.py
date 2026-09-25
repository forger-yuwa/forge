#!/usr/bin/env python3
"""usage-rule plan §5.1 #8 (codex result M3/M4) の測定。判定条件は plan に先に固定した (commit 6b785918)。

M3: run_0440 (outflow + flag 0、convMethod 0 区間、step 815 で NaN) の NaN 前の連続 3 dump で、3 CV の
    (i) ρ_w が単調減少し最後の dump で ρ_w/ρ_i < 0.1 (ρ_i = 非壁隣接ノードの平均)、
    (ii) 全接続面の正味流出 Σṁ (流出正、1 次 SLAU・flag 0) が 3 dump とも > 0。
M4: 前 plan の「膨張角窓」を旧定義の窓のまま、式だけ全 dump 対 l_ij に差し替えて再判定。
    x_f = −0.000877683 m (旧抽出)、t = 2 mm (実厚 0.5 mm と異なる登録済み評価長さ)、座標・重みは旧定義と同じ
    CELLS/centCoords の x と台形則、末尾 29 dump。帯内: max l ≤ 0.01 / 差が残る: min l > 0.01 / 他: 判定不能。

usage: r4_m3m4.py m3 RUN0440_DIR --ids 153797,153880,189814 --wall-phys-ids 1,2,3,4,10,11,12,13,15
       r4_m3m4.py m4 RUN0_EXT RUN1_EXT
"""
import argparse, glob, re, sys
from pathlib import Path
import h5py, numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from diag_wall_cv_budget import slau_mdot, parse_struct


def dumps(d, nan_ok=False):
    out = {}
    for f in glob.glob(str(Path(d) / "res_[0-9]*.h5")):
        m = re.fullmatch(r"res_(\d+)\.h5", Path(f).name)
        if m and int(m.group(1)) > 0:
            out[int(m.group(1))] = f
    return dict(sorted(out.items()))


def m3(a):
    ids = [int(s) for s in a.ids.split(",")]
    mesh = str(Path(a.run) / "sern.h5")
    with h5py.File(mesh) as f:
        S = np.asarray(f["PLANES/surfVect"], np.float64).reshape(-1, 3)
        own, nei, _ = parse_struct(np.asarray(f["PLANES/STRUCT"]), len(S))
        wall = set()
        for pid in a.wall_phys_ids.split(","):
            if pid in f["BCONDS"]:
                c = np.asarray(f["BCONDS"][pid]["iCells"]).ravel(); wall.update(c[c >= 0].tolist())
    own = own.astype(int); nei = nei.astype(int)
    d = dumps(a.run)
    good = []
    for s, fn in d.items():
        with h5py.File(fn) as g:
            ro = np.asarray(g["VALUE/ro"])
        if np.all(np.isfinite(ro)):
            good.append(s)
    use = good[-3:]
    print(f"dumps (有限) {good}; 判定に使う連続 3 dump {use} (NaN step 815 の前)")
    res = {cv: [] for cv in ids}
    for s in use:
        with h5py.File(d[s]) as g:
            V = {k: np.asarray(g["VALUE"][k], np.float64) for k in ("ro", "Ux", "Uy", "Uz", "P", "sonic")}
        for cv in ids:
            faces = np.where(((own == cv) | (nei == cv)) & (nei >= 0))[0]
            tot = 0.0; nb = []
            for k in faces:
                A = float(np.linalg.norm(S[k])); n = S[k] / A
                L = {q: V[q][own[k]] for q in V}; R = {q: V[q][nei[k]] for q in V}
                m = slau_mdot(A, *n, L, R)[0]
                tot += m if own[k] == cv else -m
                j = nei[k] if own[k] == cv else own[k]
                if j not in wall:
                    nb.append(j)
            res[cv].append((s, V["ro"][cv], float(np.mean(V["ro"][nb])), tot))
    allok = True
    for cv, rows in res.items():
        rho = [r[1] for r in rows]; out = [r[3] for r in rows]
        mono = all(rho[i + 1] < rho[i] for i in range(len(rho) - 1))
        ratio = rows[-1][1] / rows[-1][2]
        c1 = mono and ratio < 0.1; c2 = all(o > 0 for o in out)
        allok &= c1 and c2 and len(rows) == 3
        print(f"CV {cv}:")
        for s, rw, ri, o in rows:
            print(f"   step {s:6d}  ρ_w {rw:.5e}  ρ_i {ri:.5e}  ρ_w/ρ_i {rw/ri:.4f}  Σṁ(流出) {o:+.4e} kg/s")
        print(f"   (i) 単調減少 {mono} かつ 最後の ρ_w/ρ_i {ratio:.4f} < 0.1 → {c1};  (ii) Σṁ>0 が 3 dump → {c2}")
    print("VERDICT (M3):", "3 条件のうち (i)(ii) 成立" if allok else "(i)(ii) 不成立 — 限定文のまま")


def m4(a):
    t, xf = 0.002, -0.000877683
    def ramp(rd):
        with h5py.File(str(Path(rd) / "sern.h5")) as f:
            ic = np.asarray(f["BCONDS"]["4"]["iCells"]).ravel(); ic = ic[ic >= 0]
            cc = np.asarray(f["CELLS/centCoords"], np.float64).reshape(-1, 3)
        o = np.argsort(cc[ic, 0]); return ic[o], cc[ic[o], 0]
    nodes, x = ramp(a.run0)
    win = (x >= xf - 2 * t) & (x <= xf + 5 * t); wx = x[win]
    w = np.zeros(wx.size); dx = np.diff(wx); w[:-1] += dx / 2; w[1:] += dx / 2
    P = {}
    for tag, rd in (("flag0", a.run0), ("flag1", a.run1)):
        d = dumps(rd); st = list(d)[-29:]
        P[tag] = [np.asarray(h5py.File(d[s])["/VALUE/P"], np.float64)[nodes][win] for s in st]
        print(f"{tag}: 末尾 29 dump step {st[0]}..{st[-1]}")
    pbar = np.mean(P["flag0"], axis=0); nrm = np.sqrt(np.sum(w * pbar ** 2))
    l = np.array([[np.sqrt(np.sum(w * (p1 - p0) ** 2)) / nrm for p1 in P["flag1"]] for p0 in P["flag0"]])
    v = "帯内" if l.max() <= 0.01 else ("差が残る" if l.min() > 0.01 else "判定不能")
    print(f"窓 [{xf-2*t:.6f}, {xf+5*t:.6f}] m、窓内ノード {int(win.sum())}")
    print(f"l_ij (全 {l.size} 対): min {100*l.min():.4f} %  max {100*l.max():.4f} %  → {v}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); sp = ap.add_subparsers(dest="cmd", required=True)
    p3 = sp.add_parser("m3"); p3.add_argument("run"); p3.add_argument("--ids", required=True); p3.add_argument("--wall-phys-ids", required=True)
    p4 = sp.add_parser("m4"); p4.add_argument("run0"); p4.add_argument("run1")
    a = ap.parse_args()
    m3(a) if a.cmd == "m3" else m4(a)
