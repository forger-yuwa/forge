#!/usr/bin/env python3
"""実メッシュ (forge 変換済 h5) での LSQ 勾配の条件数センサスと精度比較。

目的 (plans/active/discretization-lsq-gradient.md の判断材料):
  1. **退化センサス**: 各ノード CV の LSQ 正規行列 M̂ = Σ_j d̂_j d̂_jᵀ (w=1/d² により
     単位方向のみの和) の固有値比 λ_min/λ_max の分布を測る。比が小さい = 近傍方向が
     共線/共面に近く「勾配を決める情報がデータに無い」真の退化 (演算精度では治らない)。
  2. **経路別の線形場勾配誤差**: 任意線形場 q = a + b·x に対し
       (a) 現行カーネル相当: M/b を double 蓄積 → float32 (flow_float) 格納 → double solve
           (`lsqGrad_accumInternal_d`/`lsqGrad_solve_d` の格納精度を忠実に模擬)
       (b) 全 double (理想)
       (c) 係数事前計算: double solve → 係数 c_ij を float32 格納 → float32 場に適用
           (提案方式: 実行時 double 不要)
     の 3 経路で max/percentile 誤差と near-wall 局在を比較する。

近傍セット (実装と同一): 内部双対面 (ip < nNormalPlanes) の相手 CV 中心 + 非 periodic
境界半割面の重心 (bvar 点)。periodic 境界は solver では DOF 合併で扱われるため本ツールでは
点を追加せず、periodic 隣接ノードは別枠で集計する (片側だけ見るぶんセンサスは悲観側)。

注意: 蓄積は numpy (double) で行うため、実カーネンの float32 atomicAdd 蓄積順ノイズ
(~1e-7, 全経路共通) は含まない。経路差 = solve に渡る格納精度の差のみを見る。

実行例:
  python3 check_lsq_gradient.py <mesh.h5> [--bcond bcondConfig.yaml] [--periodic 1,2,5,6]
  (--bcond 省略時はメッシュと同じディレクトリの bcondConfig.yaml を探す)
"""
import argparse
import os
import sys

import numpy as np
import h5py


def parse_plane_cells(strct, n_planes):
    """PLANES/STRUCT ([nn, nodes..., nc, cells...] 可変長) から各面のセル対を取り出す"""
    ic0 = np.empty(n_planes, dtype=np.int64)
    ic1 = np.full(n_planes, -1, dtype=np.int64)
    p = 0
    s = strct
    for ip in range(n_planes):
        nn = s[p]; p += 1 + nn
        nc = s[p]; p += 1
        ic0[ip] = s[p]
        if nc == 2:
            ic1[ip] = s[p + 1]
        p += nc
    return ic0, ic1


def load_periodic_ids(args, mesh_path):
    if args.periodic:
        return set(int(x) for x in args.periodic.split(","))
    bc_path = args.bcond or os.path.join(os.path.dirname(mesh_path), "bcondConfig.yaml")
    if os.path.exists(bc_path):
        import yaml
        with open(bc_path) as fp:
            bc = yaml.safe_load(fp)
        ids = {int(v["physID"]) for v in bc.values()
               if isinstance(v, dict) and v.get("kind") == "periodic"}
        print(f"periodic physID (from {bc_path}): {sorted(ids)}")
        return ids
    print("WARNING: bcondConfig.yaml が見つからず --periodic 未指定。全境界を bvar 点として扱う")
    return set()


def accum_sym3(idx, dx, dy, dz, w, n):
    """M += w d dᵀ を bincount で成分別に蓄積 (対称 6 成分)"""
    out = {}
    for name, a, c in (("xx", dx, dx), ("xy", dx, dy), ("xz", dx, dz),
                       ("yy", dy, dy), ("yz", dy, dz), ("zz", dz, dz)):
        out[name] = np.bincount(idx, weights=w*a*c, minlength=n)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mesh")
    ap.add_argument("--bcond", default=None)
    ap.add_argument("--periodic", default=None,
                    help="periodic の physID (カンマ区切り)。--bcond より優先")
    args = ap.parse_args()

    f = h5py.File(args.mesh, "r")
    nCells = int(f["MESH"].attrs["nCells"])
    nNormal = int(f["MESH"].attrs["nNormalPlanes"])
    nPlanes = int(f["MESH"].attrs["nPlanes"])
    cc32 = np.asarray(f["CELLS/centCoords"], dtype=np.float32).reshape(-1, 3)[:nCells]
    pc32 = np.asarray(f["PLANES/centCoords"], dtype=np.float32).reshape(-1, 3)
    wall_dist = (np.asarray(f["VALUE/wall_dist"], dtype=np.float64)
                 if "VALUE/wall_dist" in f else None)
    print(f"mesh: nCells(CV)={nCells}, nNormalPlanes={nNormal}, nPlanes={nPlanes}")

    print("PLANES/STRUCT を解析中...")
    ic0, ic1 = parse_plane_cells(np.asarray(f["PLANES/STRUCT"]), nPlanes)

    periodic_ids = load_periodic_ids(args, args.mesh)
    bnd_cells, bnd_planes, per_flag = [], [], np.zeros(nCells, dtype=bool)
    for bid in f["BCONDS"]:
        cells = np.asarray(f[f"BCONDS/{bid}/iCells"], dtype=np.int64)
        planes = np.asarray(f[f"BCONDS/{bid}/iPlanes"], dtype=np.int64)
        if int(bid) in periodic_ids:
            per_flag[cells[cells < nCells]] = True
        else:
            bnd_cells.append(cells); bnd_planes.append(planes)
    bnd_cells = np.concatenate(bnd_cells) if bnd_cells else np.empty(0, dtype=np.int64)
    bnd_planes = np.concatenate(bnd_planes) if bnd_planes else np.empty(0, dtype=np.int64)
    m = bnd_cells < nCells
    bnd_cells, bnd_planes = bnd_cells[m], bnd_planes[m]

    # ---- 近傍オフセット (座標は flow_float=float32 格納をそのまま double 化: カーネルと同一) ----
    i0, i1 = ic0[:nNormal], ic1[:nNormal]
    valid = (i0 >= 0) & (i1 >= 0) & (i0 < nCells) & (i1 < nCells)
    i0, i1 = i0[valid], i1[valid]
    d_in = (cc32[i1] - cc32[i0]).astype(np.float64)
    w_in = 1.0/np.maximum(np.sum(d_in*d_in, axis=1), 1e-300)
    d_bd = (pc32[bnd_planes] - cc32[bnd_cells]).astype(np.float64)
    w_bd = 1.0/np.maximum(np.sum(d_bd*d_bd, axis=1), 1e-300)

    # ---- M (double) を蓄積: 内部は両端点に同値、境界は片端 ----
    idx = np.concatenate([i0, i1, bnd_cells])
    DX = np.concatenate([d_in[:, 0], d_in[:, 0], d_bd[:, 0]])
    DY = np.concatenate([d_in[:, 1], d_in[:, 1], d_bd[:, 1]])
    DZ = np.concatenate([d_in[:, 2], d_in[:, 2], d_bd[:, 2]])
    W = np.concatenate([w_in, w_in, w_bd])
    Mc = accum_sym3(idx, DX, DY, DZ, W, nCells)
    M = np.empty((nCells, 3, 3))
    M[:, 0, 0] = Mc["xx"]; M[:, 0, 1] = M[:, 1, 0] = Mc["xy"]
    M[:, 0, 2] = M[:, 2, 0] = Mc["xz"]; M[:, 1, 1] = Mc["yy"]
    M[:, 1, 2] = M[:, 2, 1] = Mc["yz"]; M[:, 2, 2] = Mc["zz"]

    # ---- 2D 検出 (カーネル lsq_solve_sym3 と同一: mzz <= 1e-12·(mxx+myy) → xy 2×2) ----
    tol2d = 1.0e-12*(M[:, 0, 0] + M[:, 1, 1])
    is2d = M[:, 2, 2] <= tol2d
    n2d = int(is2d.sum())
    if n2d:
        print(f"2D ノード (z 近傍なし → 2×2 分岐): {n2d}/{nCells}")

    # ---- 1. 退化センサス: λ_min/λ_max (M は w=1/d² により Σ d̂d̂ᵀ そのもの) ----
    # 2D ノードは面内 2×2 の固有値比 (z の自明な退化はカーネルが分岐で処理するため除外)
    lam = np.linalg.eigvalsh(M)
    ratio = lam[:, 0]/np.maximum(lam[:, 2], 1e-300)
    if n2d:
        mxx, myy, mxy = M[is2d, 0, 0], M[is2d, 1, 1], M[is2d, 0, 1]
        t = 0.5*(mxx + myy)
        dsc = np.sqrt(np.maximum(0.25*(mxx - myy)**2 + mxy*mxy, 0.0))
        ratio[is2d] = (t - dsc)/np.maximum(t + dsc, 1e-300)
    inner = ~per_flag
    print("\n[1] 退化センサス λ_min/λ_max (1 に近いほど等方・小さいほど共線/共面退化)")
    for name, sel in (("periodic 非隣接", inner), ("periodic 隣接 (片側視点=悲観側)", per_flag)):
        if sel.sum() == 0:
            continue
        r = ratio[sel]
        q = np.percentile(r, [0, 0.1, 1, 50])
        print(f"  {name} ({sel.sum()} nodes): min={q[0]:.3e} p0.1={q[1]:.3e} "
              f"p1={q[2]:.3e} median={q[3]:.3e}")
        for th in (1e-2, 1e-3, 1e-4):
            nb = int((r < th).sum())
            print(f"    ratio < {th:.0e}: {nb} nodes ({nb/len(r)*100:.4f}%)")
    if wall_dist is not None:
        worst = np.argsort(ratio + per_flag*10.0)[:5]   # periodic 隣接を除いた worst
        print("  worst 5 (periodic 非隣接): " +
              ", ".join(f"ratio={ratio[i]:.2e}/wd={wall_dist[i]:.2e}" for i in worst))

    # ---- 2. 線形場勾配テスト ----
    rngb = np.array([0.7, -0.4, 0.0]) if n2d > nCells//2 else np.array([0.7, -0.4, 0.5])
    q32 = (0.3 + cc32.astype(np.float64) @ rngb).astype(np.float32)     # 場は float32 格納
    qb32 = (0.3 + pc32[bnd_planes].astype(np.float64) @ rngb).astype(np.float32)
    dq_in = (q32[i1].astype(np.float64) - q32[i0].astype(np.float64))
    dq_bd = (qb32.astype(np.float64) - q32[bnd_cells].astype(np.float64))
    WDQ = np.concatenate([w_in*dq_in, w_in*dq_in, w_bd*dq_bd])
    b = np.stack([np.bincount(idx, weights=WDQ*D, minlength=nCells)
                  for D in (DX, DY, DZ)], axis=-1)

    nb2 = np.linalg.norm(rngb)

    def report(name, g, sel=inner):
        err = np.linalg.norm(g - rngb, axis=1)/nb2
        e = err[sel]
        print(f"  {name}: max={e.max():.3e} p99.9={np.percentile(e, 99.9):.3e} "
              f"median={np.median(e):.3e}", end="")
        if wall_dist is not None:
            iw = np.argmax(err*sel)
            print(f"  (worst node wd={wall_dist[iw]:.2e})")
        else:
            print()
        return err

    def solve_masked(Mm, bm):
        """2D ノードは xy 2×2 で解き gz=0 (カーネル lsq_solve_sym3 と同一の分岐)"""
        g = np.zeros((nCells, 3))
        if n2d < nCells:
            g[~is2d] = np.linalg.solve(Mm[~is2d], bm[~is2d])
        if n2d:
            g[is2d, :2] = np.linalg.solve(Mm[is2d][:, :2, :2], bm[is2d, :2])
        return g

    def inv_masked(Mm):
        Minv = np.zeros_like(Mm)
        if n2d < nCells:
            Minv[~is2d] = np.linalg.inv(Mm[~is2d])
        if n2d:
            Minv[np.ix_(is2d, [0, 1], [0, 1])] = np.linalg.inv(Mm[is2d][:, :2, :2])
        return Minv

    print("\n[2] 線形場 q=a+b·x の勾配誤差 |g−b|/|b| (periodic 非隣接のみ集計)")
    # (b) 全 double
    g_dbl = solve_masked(M, b)
    report("(b) 全 double (理想)          ", g_dbl)
    # (a) 現行カーネル相当: M,b を float32 格納してから double solve
    Mf = M.astype(np.float32).astype(np.float64)
    bf = b.astype(np.float32).astype(np.float64)
    g_f32 = solve_masked(Mf, bf)
    report("(a) 現行相当 (float32 格納 M/b)", g_f32)
    # (c) 係数事前計算: c_ij = M⁻¹ w d (double) → float32 → float32 適用
    Minv = inv_masked(M)
    c_i0 = np.einsum("nij,nj->ni", Minv[i0], w_in[:, None]*d_in).astype(np.float32)
    c_i1 = np.einsum("nij,nj->ni", Minv[i1], w_in[:, None]*d_in).astype(np.float32)
    c_bd = np.einsum("nij,nj->ni", Minv[bnd_cells], w_bd[:, None]*d_bd).astype(np.float32)
    dq_in32 = (q32[i1] - q32[i0])
    dq_bd32 = (qb32 - q32[bnd_cells])
    g_pre = np.stack(
        [np.bincount(idx,
                     weights=np.concatenate([c_i0[:, k]*dq_in32, c_i1[:, k]*dq_in32,
                                             c_bd[:, k]*dq_bd32]).astype(np.float64),
                     minlength=nCells) for k in range(3)], axis=-1)
    report("(c) 事前計算係数 (float32 適用) ", g_pre)

    # (d) Green-Gauss (現行既定): 面値 fx 補間 + 双対面ベクトル和 / V (境界は厳密 bvar 値)
    # 軸対称 run は solver が A_planar/sx_planar (h5 に無い) で勾配を組むため本再現は不可 → skip
    sc_path = os.path.join(os.path.dirname(args.mesh), "solverConfig.yaml")
    if os.path.exists(sc_path):
        import yaml
        with open(sc_path) as fp:
            sc = yaml.safe_load(fp)
        # solver と同じ優先順 (mesh.isAxisymmetric、無ければ旧 physProp.isAxisymmetric)。以前はどの節でも 1 なら軸対称と
        # 読んでいて、mesh 0・physProp 1 (実効は平面) を軸対称にしていた
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        import check_solver_config as csc
        axi = csc.effective_axisymmetric(sc or {}) == 1
        if axi:
            print("  (d) Green-Gauss: 軸対称 run (planar 幾何が h5 に無い) のため N/A")
            if wall_dist is not None:
                near = inner & (wall_dist < np.percentile(wall_dist, 5))
                print(f"\n  near-wall (wall_dist 下位 5%, {near.sum()} nodes) の max 誤差:")
                for name, g in (("(b) double", g_dbl), ("(a) float32格納", g_f32),
                                ("(c) 事前計算", g_pre)):
                    e = (np.linalg.norm(g - rngb, axis=1)/nb2)[near]
                    print(f"    {name}: max={e.max():.3e} p99.9={np.percentile(e, 99.9):.3e}")
            return
    vol = np.asarray(f["CELLS/volume"], dtype=np.float64)[:nCells]
    sv = np.asarray(f["PLANES/surfVect"], dtype=np.float64).reshape(-1, 3)
    sv_in = sv[:nNormal][valid]
    # 向きを ic0 外向きに正規化 (dot(S, pc-cc0) > 0)
    pcd = pc32[:nNormal][valid].astype(np.float64)
    flip = np.sign(np.sum(sv_in*(pcd - cc32[i0].astype(np.float64)), axis=1))
    sv_in = sv_in*flip[:, None]
    dcc = cc32[i1].astype(np.float64) - cc32[i0].astype(np.float64)
    fx = np.clip(np.sum((cc32[i1].astype(np.float64) - pcd)*dcc, axis=1)
                 / np.maximum(np.sum(dcc*dcc, axis=1), 1e-300), 0.0, 1.0)
    qf = fx*q32[i0].astype(np.float64) + (1.0 - fx)*q32[i1].astype(np.float64)
    g_gg = np.stack(
        [np.bincount(np.concatenate([i0, i1]),
                     weights=np.concatenate([sv_in[:, k]*qf, -sv_in[:, k]*qf]),
                     minlength=nCells) for k in range(3)], axis=-1)
    sv_bd = sv[bnd_planes]
    flip_b = np.sign(np.sum(sv_bd*(pc32[bnd_planes].astype(np.float64)
                                   - cc32[bnd_cells].astype(np.float64)), axis=1))
    sv_bd = sv_bd*np.where(flip_b == 0, 1.0, flip_b)[:, None]
    for k in range(3):
        g_gg[:, k] += np.bincount(bnd_cells, weights=sv_bd[:, k]*qb32.astype(np.float64),
                                  minlength=nCells)
    g_gg /= vol[:, None]
    report("(d) Green-Gauss (現行既定)     ", g_gg)

    if wall_dist is not None:
        near = inner & (wall_dist < np.percentile(wall_dist, 5))
        print(f"\n  near-wall (wall_dist 下位 5%, {near.sum()} nodes) の max 誤差:")
        for name, g in (("(b) double", g_dbl), ("(a) float32格納", g_f32),
                        ("(c) 事前計算", g_pre), ("(d) GG", g_gg)):
            e = (np.linalg.norm(g - rngb, axis=1)/nb2)[near]
            print(f"    {name}: max={e.max():.3e} p99.9={np.percentile(e, 99.9):.3e}")


if __name__ == "__main__":
    main()
