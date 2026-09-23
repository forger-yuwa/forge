#!/usr/bin/env python3
"""V5 (#10d): 面流束 massflux のビット比較と tool 照合。plan §6 V5 の P0-P5 を判定する。

  python3 _v5_compare.py --mesh mesh/fp_y1_12um.h5 --state RUN/mesh.h5 \
      --m0 a.bin --m0b b.bin --m1 c.bin [--m1-256 d.bin] --wall-phys-ids 4
"""
import argparse, sys
from pathlib import Path
import h5py, numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "46.sern_design" / "cad"))
from diag_wall_cv_budget import slau_mdot, parse_struct

GAM, CP = 1.4, 1004.5
R = CP * (GAM - 1) / GAM

def load_bin(p, n):
    v = np.fromfile(p, dtype=np.float32)
    if v.size != n: raise SystemExit(f"{p}: {v.size} 面 != メッシュ {n} 面")
    return v

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mesh", required=True); ap.add_argument("--state", required=True)
    ap.add_argument("--m0", required=True); ap.add_argument("--m0b", required=True)
    ap.add_argument("--m1", required=True); ap.add_argument("--m1-256", default=None)
    ap.add_argument("--wall-phys-ids", required=True)
    a = ap.parse_args()

    with h5py.File(a.mesh) as f:
        S = np.asarray(f["PLANES/surfVect"], np.float64).reshape(-1, 3)
        st = np.asarray(f["PLANES/STRUCT"])
        ncv = len(np.asarray(f["CELLS/volume"]))
        wallf = np.zeros(ncv, bool)
        for k in a.wall_phys_ids.split(","):
            ic = np.asarray(f["BCONDS"][k.strip()]["iCells"]).ravel()
            wallf[ic[(ic >= 0) & (ic < ncv)]] = True
    own, nei, got = parse_struct(st, len(S))
    if got != len(S): raise SystemExit("PLANES/STRUCT を全面読めない")
    nP = len(S); A = np.linalg.norm(S, axis=1)
    print(f"# faces {nP} (内部 {(nei>=0).sum()}, 境界半割 {(nei<0).sum()})  壁ノード {int(wallf.sum())}/{ncv}")

    # **状態はカーネルが読む配列のダンプから読む** (`<massflux>.state`: ro,Ux,Uy,Uz,P,sonic x nCells)。
    # res_0.h5 の値を代理に使うと P が 1716/89440 節点で **1 ulp** ずれ、P2 の母集団を誤る
    # (§6 V5 初版が 3 面で FAIL した原因。2026-09-24)。
    if a.state.endswith(".state"):
        v = np.fromfile(a.state, np.float32)
        if v.size % 6: raise SystemExit(f"{a.state}: 6 配列に割れない")
        nc = v.size // 6
        v = v.reshape(6, nc)
        P32 = v[4].copy()
        ro, Ux, Uy, Uz, P, son = (v[i].astype(np.float64) for i in range(6))
    else:
        raise SystemExit("--state には <massflux>.state (カーネルのダンプ) を渡すこと")

    m0, m0b, m1 = (load_bin(getattr(a, k), nP) for k in ("m0", "m0b", "m1"))

    intf = nei >= 0
    o, ne = own.astype(int), nei.astype(int)
    tgt = np.zeros(nP, bool)
    tgt[intf] = wallf[o[intf]] | wallf[ne[intf]]

    def bitsame(x, y, m):
        return int((x[m] != y[m]).sum())

    print("\n=== P0 前提: flag 0 x2 が全面ビット同一 ===")
    nd = bitsame(m0, m0b, np.ones(nP, bool))
    print(f"  不一致 {nd}/{nP}  -> {'PASS' if nd==0 else 'FAIL (試験不能)'}")
    if nd: raise SystemExit("P0 不成立: 面流束が非決定的。合否を判定しない。")

    print("\n=== P1 非対象面 (両端とも非壁) のビット同一 ===")
    m_non = intf & ~tgt
    nd = bitsame(m0, m1, m_non)
    print(f"  対象 {int(m_non.sum())} 面  不一致 {nd}  -> {'PASS' if nd==0 else 'FAIL'}")

    # tool の float64 予測
    def tool(flag):
        out = np.zeros(nP)
        ii = np.where(intf)[0]
        for ip in ii:
            i0, i1 = o[ip], ne[ip]
            nx, ny, nz = S[ip]/A[ip]
            L = dict(ro=ro[i0], P=P[i0], Ux=Ux[i0], Uy=Uy[i0], Uz=Uz[i0], sonic=son[i0])
            Rr = dict(ro=ro[i1], P=P[i1], Ux=Ux[i1], Uy=Uy[i1], Uz=Uz[i1], sonic=son[i1])
            out[ip] = slau_mdot(A[ip], nx, ny, nz, L, Rr, wall_face=bool(flag and tgt[ip]))[0]
        return out
    t0, t1 = tool(0), tool(1)
    dt = t1 - t0
    robar = np.zeros(nP); robar[intf] = 0.5*(ro[o[intf]]+ro[ne[intf]])
    chat = np.zeros(nP); chat[intf] = 0.5*(son[o[intf]]+son[ne[intf]])
    tau = 1e-5 * A * robar * chat

    print("\n=== P2 対象面のうち ΔP=0 (float32 同値) の面 ===")
    dP0 = np.zeros(nP, bool)
    dP0[intf] = P32[o[intf]] == P32[ne[intf]]   # **カーネルが読む Ps** で母集団を定義する
    m_p2 = tgt & dP0
    nd = bitsame(m0, m1, m_p2)
    print(f"  対象 {int(m_p2.sum())} 面  不一致 {nd}  -> {'PASS' if nd==0 else 'FAIL'}")

    print("\n=== P3 対象面 S (壁隣接 ∧ ΔP≠0) の差を tool と照合 ===")
    S_ = tgt & ~dP0 & (np.abs(dt) > 0)
    dk = (m1 - m0).astype(np.float64)
    err = np.abs(dk - dt)
    ok_abs = err[S_] <= 2*tau[S_]
    big = S_ & (np.abs(dt) >= 100*tau)
    rel = np.abs(dk[big]-dt[big])/np.abs(dt[big]) if big.sum() else np.array([0.0])
    print(f"  S {int(S_.sum())} 面  max|Δk−Δt|/τ_b = {(err[S_]/tau[S_]).max():.3f} (許容 2)  -> {'PASS' if ok_abs.all() else 'FAIL'}")
    print(f"  |Δt|>=100τ_b の {int(big.sum())} 面  max 相対差 {rel.max():.3e} (許容 1e-2)  -> {'PASS' if rel.max()<=1e-2 else 'FAIL'}")

    print("\n=== P4 全内部面の絶対照合 |m_k − m_t| <= τ_b ===")
    for lbl, mk, tk in (("flag0", m0, t0), ("flag1", m1, t1)):
        e = np.abs(mk.astype(np.float64) - tk)[intf] / np.maximum(tau[intf], 1e-300)
        bad = int((e > 1.0).sum())
        wbad = int((e[tgt[intf]] > 1.0).sum()) if tgt[intf].any() else 0
        print(f"  {lbl}: max e/τ_b {e.max():.3f}  超過 {bad}/{int(intf.sum())} (うち対象面 {wbad})  -> {'PASS' if bad==0 else 'FAIL'}")

    if a.m1_256:
        print("\n=== P5 ブロックサイズ 128 vs 256 (flag 1) ===")
        m1b = load_bin(a.m1_256, nP)
        nd = bitsame(m1, m1b, np.ones(nP, bool))
        print(f"  不一致 {nd}/{nP}  -> {'PASS' if nd==0 else 'FAIL (面流束に blockDim 依存経路)'}")

if __name__ == "__main__":
    main()
