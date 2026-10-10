#!/usr/bin/env python3
"""plan axisymmetric-freestream-hoop-gauge §4.6 の 2 の照合 (implementer の手元確認の台本を取り込んだもの)。
plan axisymmetric-freestream-hoop-gauge §4.1 #6 の手元確認 (変換器の出力)。

usage:
  verify_conv.py diff  OLD.h5 NEW.h5            # 新旧の幾何の差 (どのデータセットのどの面が変わったか)
  verify_conv.py check NEW.h5 [--rw]            # 双対の閉性 (Σ±S、--rw なら Σ±W と /PLANES/rSurfVect の独立な作り直しとの一致)
"""
import sys
import numpy as np
import h5py

EPS = np.finfo(np.float64).eps


def planes(f):
    ps = f['PLANES/STRUCT'][:]
    nP = f['PLANES/surfArea'].shape[0]
    nodes, cells = [], []
    p = 0
    for _ in range(nP):
        nn = ps[p]; nodes.append(ps[p + 1:p + 1 + nn]); p += 1 + nn
        nc = ps[p]; cells.append(ps[p + 1:p + 1 + nc]); p += 1 + nc
    pA = np.array([c[0] for c in cells], np.int64)
    pB = np.array([c[1] if len(c) > 1 else -1 for c in cells], np.int64)
    return nP, pA, pB


def cmd_diff(old, new):
    fo, fn = h5py.File(old, 'r'), h5py.File(new, 'r')
    nInt = int(fn['MESH'].attrs['nNormalPlanes'])
    nP = fn['PLANES/surfArea'].shape[0]
    print(f'nPlanes {nP} (interior {nInt}, boundary half {nP - nInt})')
    names = []
    fo.visititems(lambda n, o: names.append(n) if isinstance(o, h5py.Dataset) else None)
    namesN = []
    fn.visititems(lambda n, o: namesN.append(n) if isinstance(o, h5py.Dataset) else None)
    only_new = sorted(set(namesN) - set(names)); only_old = sorted(set(names) - set(namesN))
    print('datasets only in NEW:', only_new, ' only in OLD:', only_old)
    for n in sorted(set(names) & set(namesN)):
        a, b = fo[n][()], fn[n][()]
        if a.shape != b.shape or a.dtype != b.dtype:
            print(f'  {n}: shape/dtype differ {a.shape}/{a.dtype} vs {b.shape}/{b.dtype}'); continue
        if a.dtype.kind in 'fc':
            neq = a.view(np.uint8).reshape(a.shape + (-1,)) != b.view(np.uint8).reshape(b.shape + (-1,))
            neq = neq.any(axis=-1) if neq.ndim > 1 else neq
        else:
            neq = a != b
        cnt = int(np.count_nonzero(neq))
        if cnt == 0:
            continue
        msg = f'  {n}: {cnt} values differ'
        if n.startswith('PLANES/') and a.dtype.kind == 'f':
            per = a.size // nP
            idx = np.nonzero(neq)[0] // per
            msg += f' -> planes {idx.min()}..{idx.max()} (interior {np.count_nonzero(idx < nInt)}, boundary {np.count_nonzero(idx >= nInt)} values)'
            d = np.abs(b - a)
            msg += f'; max|Δ| {d.max():.3e}'
            if n == 'PLANES/surfVect':
                sa = np.repeat(fn['PLANES/surfArea'][:], 3)
                msg += f', max|Δ|/|S_f| {np.max(d / np.maximum(sa, 1e-300)):.3e}'
            elif n == 'PLANES/centCoords':
                # 面の大きさ (面積) と比べる
                sa = np.repeat(fn['PLANES/surfArea'][:], 3)
                msg += f', max|Δ|/|S_f| {np.max(d / np.maximum(sa, 1e-300)):.3e}, max|Δ|/|x| {np.max(d / np.maximum(np.abs(a), 1e-300)):.3e}'
            else:
                msg += f', max|Δ|/|old| {np.max(d / np.maximum(np.abs(a), 1e-300)):.3e}'
        elif a.dtype.kind == 'f':
            d = np.abs(b - a)
            msg += f'; max|Δ| {d.max():.3e}, max|Δ|/|old| {np.max(d[neq] / np.maximum(np.abs(a[neq]), 1e-300)):.3e}'
        print(msg)


def primal_cells(f):
    cw = f['VIZMESH/CONNE'][:]
    cells = []
    p = 0
    while p < len(cw):
        code = int(cw[p]); nn = {5: 4, 4: 3}[code]
        cells.append(cw[p + 1:p + 1 + nn]); p += 1 + nn
    return cells


def cmd_check(path, rw):
    f = h5py.File(path, 'r')
    X = f['MESH/COORD'][:].reshape(-1, 3); x, y = X[:, 0], X[:, 1]
    nN = X.shape[0]
    nP, pA, pB = planes(f)
    nInt = int(f['MESH'].attrs['nNormalPlanes'])
    SV = f['PLANES/surfVect'][:].reshape(-1, 3)
    PC = f['PLANES/centCoords'][:].reshape(-1, 3)
    Avol = f['CELLS/volume'][:]
    has_rw = 'PLANES/rSurfVect' in f
    print(f'{path}: nNodes {nN} nPlanes {nP} (interior {nInt}) /PLANES/rSurfVect present: {has_rw}'
          + (f' shape {f["PLANES/rSurfVect"].shape} dtype {f["PLANES/rSurfVect"].dtype}' if has_rw else ''))

    def node_sum(v):
        return (np.bincount(pA[:nInt], v[:nInt], nN) - np.bincount(pB[:nInt], v[:nInt], nN)
                + np.bincount(pA[nInt:], v[nInt:], nN))

    def node_abs(v):
        av = np.abs(v)
        return np.bincount(pA[:nInt], av[:nInt], nN) + np.bincount(pB[:nInt], av[:nInt], nN) + np.bincount(pA[nInt:], av[nInt:], nN)

    # 平面の閉性 Σ±S = 0
    for k, c in ((0, 'x'), (1, 'y')):
        s = node_sum(SV[:, k]); a = node_abs(SV[:, k])
        over = int(np.count_nonzero(np.abs(s) > 100 * EPS * a))
        print(f'  planar closure Σ±S_{c}: max|.|/Σ|S| {np.max(np.abs(s) / np.maximum(a, 1e-300)):.3e}, '
              f'max|.|/A {np.max(np.abs(s) / Avol):.3e}, CVs over 100·ε64·Σ|S|: {over}')
    if not rw:
        return
    if not has_rw:
        print('  ERROR: --rw but /PLANES/rSurfVect missing'); sys.exit(1)
    RW = f['PLANES/rSurfVect'][:].reshape(-1, 3)
    # Σ±W 閉性 (E_x, E_y)
    Sx = node_sum(RW[:, 0]); Sy = node_sum(RW[:, 1]); Ax = node_abs(RW[:, 0]); Ay = node_abs(RW[:, 1])
    Ex = np.abs(Sx) / Avol; Ey = np.abs(Sy - Avol) / Avol
    ox = int(np.count_nonzero(np.abs(Sx) > 100 * EPS * (Avol + Ax)))
    oy = int(np.count_nonzero(np.abs(Sy - Avol) > 100 * EPS * (Avol + Ay)))
    print(f'  r-weighted closure (rSurfVect): max E_x {Ex.max():.3e}  max E_y {Ey.max():.3e}  '
          f'CVs over 100·ε64·(A+Σ|W|): x {ox} y {oy}  (rSurfVect z max |.| {np.max(np.abs(RW[:, 2])):.1e})')
    # 比較: 今の r̄·S (同じ HDF5 の surfVect × max(pcy,1e-20))
    rf = np.maximum(PC[:, 1], 1e-20)
    Sx0 = node_sum(SV[:, 0] * rf); Sy0 = node_sum(SV[:, 1] * rf)
    print(f'  (reference) r̄·S on the same file: max E_x {np.max(np.abs(Sx0) / Avol):.3e}  max E_y {np.max(np.abs(Sy0 - Avol) / Avol):.3e}')
    # 向き: W と S の内積が正 (軸の上の面 W = 0 を除く)
    dot = np.einsum('ij,ij->i', RW, SV)
    nz = np.linalg.norm(RW, axis=1) > 0
    print(f'  orientation: faces with W·S <= 0 among W != 0: {int(np.count_nonzero(dot[nz] <= 0))}; '
          f'faces with W == 0: {int(np.count_nonzero(~nz))} (centroid r of those: max {PC[~nz, 1].max() if (~nz).any() else float("nan"):.1e})')

    # ---- 独立な作り直し (内部面: 区間 M̃→G、境界の半割面: N→M̃) ----
    cells = primal_cells(f)
    G = np.array([X[c].sum(axis=0) / len(c) if False else _vmean(X, c) for c in cells])
    key2ip = {}
    for ip in range(nInt):
        a, b = int(pA[ip]), int(pB[ip]); key2ip[(min(a, b), max(a, b))] = ip
    Wint = np.zeros((nInt, 2)); Sint = np.zeros((nInt, 2))
    for ic, c in enumerate(cells):
        m = len(c)
        for k in range(m):
            a, b = int(c[k]), int(c[(k + 1) % m])
            ip = key2ip.get((min(a, b), max(a, b)))
            if ip is None:
                continue  # 境界エッジ: 内部面は無い (区間 M̃→G は境界エッジの内部面として存在する場合のみ)
            A_, B_ = int(pA[ip]), int(pB[ip])
            Mx = 0.5 * (x[A_] + x[B_]); My = 0.5 * (y[A_] + y[B_])
            sx = G[ic, 0] - Mx; sy = G[ic, 1] - My
            nx, ny = sy, -sx
            if nx * (x[B_] - x[A_]) + ny * (y[B_] - y[A_]) < 0:
                nx, ny = -nx, -ny
            rk = My + 0.5 * sy
            Wint[ip, 0] += rk * nx; Wint[ip, 1] += rk * ny
            Sint[ip, 0] += nx; Sint[ip, 1] += ny
    # 境界エッジは内部面 (primal エッジ) としても存在する: 隣接セル 1 つの区間 M̃→G。key2ip にあれば上で入っている。
    dWi = np.linalg.norm(Wint - RW[:nInt, :2], axis=1) / np.maximum(np.linalg.norm(RW[:nInt, :2], axis=1), 1e-300)
    dSi = np.linalg.norm(Sint - SV[:nInt, :2], axis=1) / np.maximum(np.linalg.norm(SV[:nInt, :2], axis=1), 1e-300)
    print(f'  interior faces: rebuilt W vs rSurfVect max |ΔW|/|W| {dWi.max():.3e} (bit-equal {int(np.count_nonzero(dWi == 0))}/{nInt}); '
          f'rebuilt S vs surfVect max {dSi.max():.3e}')
    # 境界の半割面
    bc_of_plane = {}
    walls = []
    bk = {}
    for name in f['BCONDS']:
        g = f['BCONDS'][name]
        kind = g.attrs['bcondKind']; kind = kind.decode() if isinstance(kind, bytes) else str(kind)
        bk[name] = kind
        if kind in ('wall', 'wall_isothermal'):
            walls.append(name)
    hp_of = {}
    for name in f['BCONDS']:
        ipl = f['BCONDS'][name]['iPlanes'][:]
        for ip in ipl[ipl >= 0]:
            hp_of[(name, int(pA[ip]))] = int(ip)
    Wb = np.zeros((nP, 2)); Sb = np.zeros((nP, 2)); Cb = np.zeros((nP, 2)); Lb = np.zeros(nP)
    # 境界エッジの隣接 primal セル (外向き判定)
    edge2cell = {}
    for ic, c in enumerate(cells):
        m = len(c)
        for k in range(m):
            a, b = int(c[k]), int(c[(k + 1) % m])
            edge2cell.setdefault((min(a, b), max(a, b)), []).append(ic)
    nmiss = 0
    for name in f['BCONDS']:
        g = f['BCONDS'][name]
        if 'vizBfaceNodes' not in g:
            continue
        vb = g['vizBfaceNodes'][:]
        sz = g['vizBfaceSizes'][:] if 'vizBfaceSizes' in g else np.full(len(vb) // 2, 2)
        off = 0
        for s_ in sz:
            a, b = int(vb[off]), int(vb[off + 1]); off += s_
            q = edge2cell[(min(a, b), max(a, b))][0]
            s0, s1 = y[b] - y[a], -(x[b] - x[a])
            Mx = 0.5 * (x[a] + x[b]); My = 0.5 * (y[a] + y[b])
            if s0 * (Mx - G[q, 0]) + s1 * (My - G[q, 1]) < 0:
                s0, s1 = -s0, -s1
            for N in (a, b):
                h0, h1 = My - y[N], -(Mx - x[N])
                if h0 * s0 + h1 * s1 < 0:
                    h0, h1 = -h0, -h1
                r = 0.5 * (y[N] + My)
                ip = hp_of.get((name, N))
                if ip is None and bk[name].startswith('inlet_'):
                    for w in walls:
                        ip = hp_of.get((w, N))
                        if ip is not None:
                            break
                if ip is None:
                    nmiss += 1; continue
                L = np.hypot(h0, h1)
                Wb[ip] += (r * h0, r * h1); Sb[ip] += (h0, h1)
                Cb[ip] += (L * 0.5 * (x[N] + Mx), L * r); Lb[ip] += L
    hb = np.arange(nInt, nP)
    dWb = np.linalg.norm(Wb[hb] - RW[hb, :2], axis=1) / np.maximum(np.linalg.norm(RW[hb, :2], axis=1), 1e-300)
    dWb[np.linalg.norm(RW[hb, :2], axis=1) == 0] = np.linalg.norm(Wb[hb][np.linalg.norm(RW[hb, :2], axis=1) == 0], axis=1)
    dSb = np.linalg.norm(Sb[hb] - SV[hb, :2], axis=1) / np.maximum(np.linalg.norm(SV[hb, :2], axis=1), 1e-300)
    dCb = np.linalg.norm(Cb[hb] / Lb[hb, None] - PC[hb, :2], axis=1) / np.maximum(np.linalg.norm(SV[hb, :2], axis=1), 1e-300)
    print(f'  boundary half faces ({len(hb)}, unmapped {nmiss}): rebuilt W vs rSurfVect max |ΔW|/|W| {dWb.max():.3e}; '
          f'rebuilt S vs surfVect max |ΔS|/|S| {dSb.max():.3e}; rebuilt centroid vs centCoords max |Δc|/|S| {dCb.max():.3e}')


def _vmean(X, c):
    s = np.zeros(3)
    for n in c:
        s = s + X[n]
    return s / len(c)


if __name__ == '__main__':
    if sys.argv[1] == 'diff':
        cmd_diff(sys.argv[2], sys.argv[3])
    else:
        cmd_check(sys.argv[2], '--rw' in sys.argv)
