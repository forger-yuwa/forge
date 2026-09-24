#!/usr/bin/env python3
"""V7 (#10c): 固定点の判定。plan convection-slau-wall-normal-chi §6 V7 の定義どおりに距離を出す。

  D_intra(cfl) = 種 dump -> 延長終端     (その run が静止しているか)
  D_inter      = 両延長終端どうし        (固定点が CFL に依るか)

ノルム: 壁隣接集合が主・全域が副。分母は**同じ集合上のその列自身の RMS**、
ただし運動量 3 成分は共通分母 ||rho|u|||_2 (壁集合で roUz≈0 だと自分の RMS で発散する)。
節点数等重み。列ごとに個別判定 (束ねない)。
"""
import argparse, h5py, numpy as np

WALL_IDS = [1,2,3,4,10,11,12,13,15]          # cowl_base,cowl_in,cowl_out,cowl_side,ramp,sidewall_end,sidewall_in,sidewall_out,vehicle
COLS = ["ro","roUx","roUy","roUz","roe"]

def load(fn):
    with h5py.File(fn,"r") as f:
        d = {k: np.asarray(f["/VALUE/"+k], np.float64) for k in COLS}
        if "P" in f["/VALUE"]: d["P"] = np.asarray(f["/VALUE/P"], np.float64)
    return d

def adj_mask(mesh, ncv, wm):
    """壁ノードに面で隣接する**内部**ノード集合 (第一内部ノード)。
    壁ノード自身は nodeWallDirichlet で u=0 にピンされ roU が恒等 0 なので相対量が定義できない。"""
    import sys
    sys.path.insert(0, "cad")
    from diag_wall_cv_budget import parse_struct
    with h5py.File(mesh,"r") as f:
        S = np.asarray(f["PLANES/surfVect"]).reshape(-1,3)
        st = np.asarray(f["PLANES/STRUCT"])
    own, nei, _ = parse_struct(st, len(S))
    o, n = own.astype(int), nei.astype(int)
    m = np.zeros(ncv, bool)
    ok = n >= 0
    sel = ok & wm[o] & ~wm[n]; m[n[sel]] = True
    sel = ok & wm[n] & ~wm[o]; m[o[sel]] = True
    return m


def wall_mask(mesh, ncv):
    m = np.zeros(ncv, bool)
    with h5py.File(mesh,"r") as f:
        for k in WALL_IDS:
            if str(k) in f["BCONDS"]:
                ic = np.asarray(f["BCONDS"][str(k)]["iCells"]).ravel()
                m[ic[(ic>=0)&(ic<ncv)]] = True
    return m

def dist(a, b, m):
    """列ごとの相対 L2 [%] と L∞ [%]。運動量は共通分母。"""
    out = {}
    mom = np.sqrt(sum(a[k][m]**2 for k in ("roUx","roUy","roUz")))
    den_mom = np.sqrt(np.mean(mom**2))
    for k in list(a.keys()):
        if k not in b: continue
        d = b[k][m] - a[k][m]
        den = den_mom if k.startswith("roU") else np.sqrt(np.mean(a[k][m]**2))
        if den == 0: out[k] = (float("nan"), float("nan"), -1); continue
        l2 = 100.0*np.sqrt(np.mean(d**2))/den
        li = 100.0*np.abs(d).max()/den
        out[k] = (l2, li, int(np.argmax(np.abs(d))))
    return out

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--mesh", required=True)
    ap.add_argument("--pairs", nargs="+", required=True, help="LABEL=A.h5,B.h5 を並べる")
    a = ap.parse_args()
    ncv = None
    for spec in a.pairs:
        lbl, files = spec.split("=",1); fa, fb = files.split(",")
        A, B = load(fa), load(fb)
        if ncv is None:
            ncv = len(A["ro"]); wm = wall_mask(a.mesh, ncv); gm = np.ones(ncv, bool)
            am = adj_mask(a.mesh, ncv, wm)
            print(f"# CV {ncv}  壁ノード {int(wm.sum())}  第一内部ノード {int(am.sum())}")
        print(f"\n=== {lbl} ===")
        print(f"  {'列':8s} {'壁 L2 [%]':>12s} {'第一内部 L2':>12s} {'第一内部 L∞':>12s} {'全域 L2 [%]':>12s}")
        dw, da, dg = dist(A,B,wm), dist(A,B,am), dist(A,B,gm)
        for k in ("ro","roUx","roUy","roUz","roe","P"):
            if k not in dw: continue
            print(f"  {k:8s} {dw[k][0]:12.5f} {da[k][0]:12.5f} {da[k][1]:12.5f} {dg[k][0]:12.5f}")
